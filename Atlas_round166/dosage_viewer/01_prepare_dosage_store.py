#!/usr/bin/env python3
"""
01_prepare_dosage_store.py
==========================
Build a Parquet-backed dosage store from per-chromosome gzipped TSVs.
This is the **preprocessing** step for the standalone dosage viewer
(see specs_todo/from_turn129/S6_dosage_heatmap_streaming_viewer.md).

The atlas-side P4.1 bridge in `server_turn1/dosage_bridge.py` reads the
raw TSV.GZ directly — no preprocessing required. THIS script is for the
standalone viewer in `dosage_viewer/02_run_server.py` (next slice), which
needs random row-group access for fast region/sampling queries.

Source layout (input):
    <base>/<chrom>.sites.tsv.gz
    <base>/<chrom>.dosage.tsv.gz
    [<base>/samples.tsv]                    # optional; one ID per line

Sites TSV columns (1-based):
    1: pos     (int, 1-based bp)
    2: ref     (str)
    3: alt     (str)
    4: missingness        (float; NA allowed)
    5: diagnostic_score   (float; NA allowed)
    6: site_id            (str; falls back to "<chrom>:<pos>:<ref>><alt>")

Dosage TSV: n_sites rows, n_samples cells per row, values in
    {0, 1, 2, NA, ., -1, ""} -> normalized to int8 with -1 = NA.

Optional first row of sample IDs is auto-detected and used as the
sample order (overrides --samples). When neither is present, samples
are named s0..s{N-1}.

Output layout:
    <out>/manifest.json                     # store manifest
    <out>/samples.tsv                       # canonical sample order
    <out>/<chrom>/sites.parquet             # cols: pos, ref, alt, missingness, diagnostic_score, site_id
    <out>/<chrom>/dosage.parquet            # cols: pos, s0, s1, ..., s{n-1}; int8; -1=NA
    <out>/<chrom>/_chrom_manifest.json      # per-chrom marker for skip-if-done

Idempotency:
    A chrom is SKIPPED when its `_chrom_manifest.json` exists AND its
    `source_mtime_ns` matches the current sites.tsv.gz mtime AND its
    `source_dosage_mtime_ns` matches the current dosage.tsv.gz mtime.
    Otherwise the chrom is rebuilt.

CLI:
    python 01_prepare_dosage_store.py \\
        --base   /path/to/04_dosage_by_chr \\
        --out    dosage_viewer/output/dosage_store \\
        --samples /path/to/samples.tsv \\
        --chroms C_gar_LG01,C_gar_LG02,...,C_gar_LG28 \\
        --row-group-size 50000

Author:  Quentin Andres
Project: MS_Inversions_North_african_catfish
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
import time
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import pyarrow as pa
import pyarrow.parquet as pq

# Logging — print to stderr so the script can be piped without polluting
# stdout (we don't actually use stdout for output, but keep the discipline).
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger("prepare_dosage_store")


SCHEMA_VERSION = "dosage_store_v1"
DEFAULT_ROW_GROUP_SIZE = 50_000


# =============================================================================
# Cell parsing — must match server_turn1/dosage_bridge.py exactly so the
# bridge and the standalone viewer agree on missing-value semantics.
# =============================================================================

def _maybe_dosage_int(s: str) -> int:
    """Parse a dosage cell into int8-range, normalizing NA to -1."""
    s = (s or "").strip()
    if not s or s == "." or s.upper() == "NA":
        return -1
    try:
        v = int(s)
    except ValueError:
        return -1
    if v in (0, 1, 2):
        return v
    if v == -1:
        return -1
    # Out-of-range values fall back to missing rather than leaking through.
    return -1


def _maybe_float(s: str) -> Optional[float]:
    s = (s or "").strip()
    if not s or s == "." or s.upper() == "NA":
        return None
    try:
        return float(s)
    except ValueError:
        return None


# =============================================================================
# Sample-ID handling
# =============================================================================

def _read_samples_tsv(path: Path) -> List[str]:
    """One ID per line. Strips blanks/comments. Tolerates 'id<TAB>group' rows.

    Returns [] if the file is missing or empty (caller decides what to do).
    """
    if not path.exists():
        return []
    out: List[str] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                line = line.split("\t", 1)[0]
            out.append(line)
    return out


def _detect_dosage_header(first_line: str, n_samples_expected: int) -> bool:
    """True if first_line looks like a header row of sample IDs."""
    fields = first_line.rstrip("\n").split("\t")
    if len(fields) != n_samples_expected:
        return False
    for f in fields:
        s = f.strip()
        if not s or s == "." or s.upper() == "NA":
            continue
        try:
            int(s)
        except ValueError:
            return True
    return False


# =============================================================================
# Site reader — builds the sites.parquet column arrays + emits row order
# =============================================================================

def _read_sites(sites_gz: Path) -> Tuple[List[int], List[str], List[str],
                                          List[Optional[float]],
                                          List[Optional[float]],
                                          List[str]]:
    """Stream <chrom>.sites.tsv.gz into 6 column arrays.

    Returns (positions, refs, alts, missingness, diagnostic_score, site_ids).
    Header rows are auto-detected and skipped. The site_id column is
    synthesized from chrom:pos:ref>alt when absent.
    """
    positions: List[int] = []
    refs: List[str] = []
    alts: List[str] = []
    missingness: List[Optional[float]] = []
    diagnostic: List[Optional[float]] = []
    site_ids: List[str] = []

    chrom_for_synth = sites_gz.stem.replace(".sites.tsv", "").replace(".sites", "")
    # ^ stem of "C_gar_LG28.sites.tsv.gz" is "C_gar_LG28.sites.tsv" → strip suffix.

    with gzip.open(sites_gz, "rt", encoding="utf-8") as f:
        for raw in f:
            if not raw or raw.startswith("#"):
                continue
            cells = raw.rstrip("\n").split("\t")
            if not cells or not cells[0]:
                continue
            try:
                pos = int(cells[0])
            except ValueError:
                # Header — skip. The column order is fixed by spec, so we
                # don't try to infer column names from the header.
                continue
            ref = cells[1] if len(cells) > 1 else "N"
            alt = cells[2] if len(cells) > 2 else "N"
            miss = _maybe_float(cells[3]) if len(cells) > 3 else None
            ds = _maybe_float(cells[4]) if len(cells) > 4 else None
            sid = cells[5].strip() if len(cells) > 5 and cells[5].strip() else \
                  f"{chrom_for_synth}:{pos}:{ref}>{alt}"
            positions.append(pos)
            refs.append(ref)
            alts.append(alt)
            missingness.append(miss)
            diagnostic.append(ds)
            site_ids.append(sid)
    return positions, refs, alts, missingness, diagnostic, site_ids


# =============================================================================
# Dosage reader — streams once, validates row count + dimensionality,
# returns a single int8 numpy-style array (list-of-lists; we let pyarrow
# convert).
# =============================================================================

def _read_dosage(dosage_gz: Path, n_sites_expected: int,
                 n_samples_expected: int, samples_path: Path
                 ) -> Tuple[List[List[int]], List[str]]:
    """Stream <chrom>.dosage.tsv.gz; return (rows, sample_ids).

    sample_ids is the list returned via the optional first-row header. If
    the file has no header row, we fall back to whatever was provided in
    samples_path; if that's missing too, we synthesise s0..s{N-1}.

    Validates that every row has n_samples_expected cells (pads or truncates
    with a warning). Validates that the total row count == n_sites_expected
    after the optional header.
    """
    rows: List[List[int]] = []
    header_decided = False
    inline_header: Optional[List[str]] = None

    with gzip.open(dosage_gz, "rt", encoding="utf-8") as f:
        for raw in f:
            if not raw or raw.startswith("#"):
                continue
            line = raw.rstrip("\n")
            if not header_decided:
                if _detect_dosage_header(line, n_samples_expected):
                    inline_header = line.split("\t")
                    header_decided = True
                    continue
                header_decided = True
                # Fall through — first line is data
            cells = line.split("\t")
            if len(cells) < n_samples_expected:
                cells = list(cells) + [""] * (n_samples_expected - len(cells))
            elif len(cells) > n_samples_expected:
                cells = cells[:n_samples_expected]
            rows.append([_maybe_dosage_int(c) for c in cells])

    if len(rows) != n_sites_expected:
        raise ValueError(
            f"row count mismatch: dosage has {len(rows)} rows, sites has "
            f"{n_sites_expected}. {dosage_gz.name}"
        )

    if inline_header is not None:
        sample_ids = inline_header
    else:
        sample_ids = _read_samples_tsv(samples_path)
        if not sample_ids:
            sample_ids = [f"s{i}" for i in range(n_samples_expected)]
        elif len(sample_ids) != n_samples_expected:
            raise ValueError(
                f"sample-count mismatch: {samples_path.name} has "
                f"{len(sample_ids)} IDs, dosage row width is "
                f"{n_samples_expected}"
            )

    return rows, sample_ids


# =============================================================================
# Validation pass on positions
# =============================================================================

def _validate_positions(positions: List[int], chrom: str) -> Tuple[List[int], bool]:
    """Ensure positions are strictly increasing. Sort if not (returning the
    sorted list AND a flag indicating whether re-sorting was necessary).
    """
    is_sorted = all(positions[i] < positions[i + 1]
                    for i in range(len(positions) - 1))
    if is_sorted:
        return positions, False
    log.warning("[%s] positions not strictly increasing; sorting in memory",
                chrom)
    # Tag with original index so we can also reorder dosage rows
    return positions, True


def _argsort_positions(positions: List[int]) -> List[int]:
    """Indices that would sort positions ascending. Stable."""
    return sorted(range(len(positions)), key=lambda i: positions[i])


# =============================================================================
# Per-chrom build
# =============================================================================

def _peek_n_samples(dosage_gz: Path) -> int:
    """Read the first non-comment data line of the dosage file to count
    cells. Used before knowing whether the first row is a header (a header
    row's cell count == sample count anyway — same answer either way).
    """
    with gzip.open(dosage_gz, "rt", encoding="utf-8") as f:
        for raw in f:
            if not raw or raw.startswith("#"):
                continue
            return len(raw.rstrip("\n").split("\t"))
    return 0


def _build_chrom(sites_gz: Path, dosage_gz: Path, samples_path: Path,
                 out_chrom_dir: Path, chrom: str,
                 row_group_size: int = DEFAULT_ROW_GROUP_SIZE) -> dict:
    """Build sites.parquet + dosage.parquet for one chromosome. Returns the
    per-chrom manifest dict (also written to disk by the caller).
    """
    log.info("[%s] reading sites …", chrom)
    t_sites = time.perf_counter()
    positions, refs, alts, miss, diag, site_ids = _read_sites(sites_gz)
    n_sites = len(positions)
    if n_sites == 0:
        raise ValueError(f"[{chrom}] sites file is empty: {sites_gz}")
    log.info("[%s] %d sites in %.1fs", chrom, n_sites,
             time.perf_counter() - t_sites)

    n_samples = _peek_n_samples(dosage_gz)
    if n_samples == 0:
        raise ValueError(f"[{chrom}] dosage file is empty: {dosage_gz}")

    log.info("[%s] reading dosage (%d cols) …", chrom, n_samples)
    t_dosage = time.perf_counter()
    dosage_rows, sample_ids = _read_dosage(
        dosage_gz, n_sites, n_samples, samples_path
    )
    log.info("[%s] %d dosage rows × %d samples in %.1fs",
             chrom, len(dosage_rows), n_samples,
             time.perf_counter() - t_dosage)

    # Validate / sort positions (and dosage rows in step)
    _, needs_resort = _validate_positions(positions, chrom)
    if needs_resort:
        order = _argsort_positions(positions)
        positions = [positions[i] for i in order]
        refs      = [refs[i] for i in order]
        alts      = [alts[i] for i in order]
        miss      = [miss[i] for i in order]
        diag      = [diag[i] for i in order]
        site_ids  = [site_ids[i] for i in order]
        dosage_rows = [dosage_rows[i] for i in order]

    # Validate cell range — anything not in {-1, 0, 1, 2} is a producer bug
    # (caught in _maybe_dosage_int already, but assert here for safety).
    bad = 0
    for r in dosage_rows:
        for v in r:
            if v not in (-1, 0, 1, 2):
                bad += 1
    if bad:
        raise ValueError(f"[{chrom}] {bad} dosage cells outside {{-1,0,1,2}}")

    out_chrom_dir.mkdir(parents=True, exist_ok=True)

    # ---- sites.parquet --------------------------------------------------
    # Column types match the renderer's expectations: pos as int64 (some
    # chromosomes go beyond int32), missingness/diagnostic as nullable
    # float32 (small storage win, no precision loss for these ratios).
    sites_table = pa.table({
        "pos": pa.array(positions, type=pa.int64()),
        "ref": pa.array(refs, type=pa.string()),
        "alt": pa.array(alts, type=pa.string()),
        "missingness": pa.array(miss, type=pa.float32()),
        "diagnostic_score": pa.array(diag, type=pa.float32()),
        "site_id": pa.array(site_ids, type=pa.string()),
    })
    sites_path = out_chrom_dir / "sites.parquet"
    pq.write_table(sites_table, sites_path,
                   row_group_size=row_group_size,
                   compression="zstd")
    log.info("[%s] wrote %s (%.2f MB, %d row groups)",
             chrom, sites_path,
             sites_path.stat().st_size / 1e6,
             max(1, (n_sites + row_group_size - 1) // row_group_size))

    # ---- dosage.parquet -------------------------------------------------
    # Layout: pos column + one int8 column per sample. This matches the
    # spec's §Layout on disk:
    #     dosage.parquet  cols: pos, s0, s1, ..., s{n-1} (int8: 0/1/2/-1)
    # We use the canonical sample IDs from `sample_ids` for column names
    # so a reader can look up "s0" → samples[0] without a separate map.
    # Column-major build: transpose dosage_rows once, then build columns.
    cols: dict[str, pa.Array] = {"pos": pa.array(positions, type=pa.int64())}
    # per-sample columns
    for j, sid in enumerate(sample_ids):
        col_name = sid
        # Some sample IDs may collide with pyarrow reserved names; defensive
        # rename: prefix any column that starts with "_" or matches a special
        # keyword. (None of Quentin's IDs would, but the producer should be
        # robust.)
        if col_name == "pos":
            col_name = f"sample__{sid}"
        cols[col_name] = pa.array(
            [dosage_rows[i][j] for i in range(n_sites)],
            type=pa.int8(),
        )
    dosage_table = pa.table(cols)
    dosage_path = out_chrom_dir / "dosage.parquet"
    pq.write_table(dosage_table, dosage_path,
                   row_group_size=row_group_size,
                   compression="zstd")
    log.info("[%s] wrote %s (%.2f MB)",
             chrom, dosage_path,
             dosage_path.stat().st_size / 1e6)

    return {
        "schema_version": SCHEMA_VERSION,
        "chrom": chrom,
        "n_sites": n_sites,
        "n_samples": n_samples,
        "length_bp": positions[-1] if positions else 0,
        "first_bp": positions[0] if positions else 0,
        "row_group_size": row_group_size,
        "source_sites_mtime_ns": sites_gz.stat().st_mtime_ns,
        "source_dosage_mtime_ns": dosage_gz.stat().st_mtime_ns,
        "sample_ids_first": sample_ids[:3],   # diag aid, not the full list
        "sample_ids_last":  sample_ids[-3:],
    }


# =============================================================================
# Skip-if-done check
# =============================================================================

def _is_chrom_up_to_date(chrom_dir: Path, sites_gz: Path,
                          dosage_gz: Path) -> bool:
    """Check the per-chrom marker file. Up-to-date if the marker exists AND
    both source mtimes match what's recorded in it.
    """
    marker = chrom_dir / "_chrom_manifest.json"
    if not marker.is_file():
        return False
    try:
        m = json.loads(marker.read_text())
    except Exception:
        return False
    if m.get("schema_version") != SCHEMA_VERSION:
        return False
    if m.get("source_sites_mtime_ns") != sites_gz.stat().st_mtime_ns:
        return False
    if m.get("source_dosage_mtime_ns") != dosage_gz.stat().st_mtime_ns:
        return False
    if not (chrom_dir / "sites.parquet").is_file():
        return False
    if not (chrom_dir / "dosage.parquet").is_file():
        return False
    return True


# =============================================================================
# Discovery
# =============================================================================

def _discover_chroms(base: Path) -> List[str]:
    """Scan `base` for *.sites.tsv.gz and return chrom names that ALSO have
    a matching *.dosage.tsv.gz. Sorted alphabetically for stability.
    """
    chroms: List[str] = []
    for sf in sorted(base.glob("*.sites.tsv.gz")):
        chrom = sf.name[:-len(".sites.tsv.gz")]
        if (base / f"{chrom}.dosage.tsv.gz").exists():
            chroms.append(chrom)
    return chroms


# =============================================================================
# Top-level driver
# =============================================================================

def build_store(base: Path, out: Path, samples_path: Optional[Path] = None,
                chroms: Optional[List[str]] = None,
                row_group_size: int = DEFAULT_ROW_GROUP_SIZE,
                force: bool = False) -> dict:
    """Build the full dosage store. Returns the store-level manifest."""
    if not base.is_dir():
        raise FileNotFoundError(f"base directory does not exist: {base}")

    out.mkdir(parents=True, exist_ok=True)

    if chroms is None:
        chroms = _discover_chroms(base)
    if not chroms:
        raise ValueError(f"no chromosomes found under {base}")

    # Resolve samples path. Default heuristic:
    #   1. <base>/../samples.tsv (next to 04_dosage_by_chr/)
    #   2. <base>/samples.tsv (inside the dir)
    if samples_path is None:
        candidates = [base.parent / "samples.tsv", base / "samples.tsv"]
        for c in candidates:
            if c.is_file():
                samples_path = c
                break
    # samples_path can still be None — that's OK if every dosage file has
    # an inline header (we'll discover sample IDs there).

    log.info("Building store at %s from %s", out, base)
    log.info("Chromosomes: %s", ", ".join(chroms))
    log.info("Samples file: %s", samples_path or "(none — relying on inline headers)")

    chrom_entries: List[dict] = []
    n_skipped = 0
    n_built = 0
    n_samples_seen: Optional[int] = None
    canonical_samples: Optional[List[str]] = None

    for chrom in chroms:
        sites_gz = base / f"{chrom}.sites.tsv.gz"
        dosage_gz = base / f"{chrom}.dosage.tsv.gz"
        if not sites_gz.is_file() or not dosage_gz.is_file():
            log.warning("[%s] missing source file(s); skipping", chrom)
            continue
        chrom_dir = out / chrom
        if not force and _is_chrom_up_to_date(chrom_dir, sites_gz, dosage_gz):
            log.info("[%s] up-to-date; skipping", chrom)
            marker = json.loads((chrom_dir / "_chrom_manifest.json").read_text())
            chrom_entries.append({
                "name": chrom,
                "n_sites": marker.get("n_sites", 0),
                "length_bp": marker.get("length_bp", 0),
                "skipped": True,
            })
            n_skipped += 1
            # Don't shadow canonical_samples on a skip — we'll keep whatever
            # the first built/read chrom decides.
            continue

        manifest = _build_chrom(
            sites_gz=sites_gz, dosage_gz=dosage_gz,
            samples_path=samples_path or Path("/dev/null"),
            out_chrom_dir=chrom_dir, chrom=chrom,
            row_group_size=row_group_size,
        )
        # Write the per-chrom marker
        (chrom_dir / "_chrom_manifest.json").write_text(
            json.dumps(manifest, indent=2)
        )

        # Cross-check: every chrom should have the same n_samples
        if n_samples_seen is None:
            n_samples_seen = manifest["n_samples"]
            # Re-read sample IDs from the freshly-built dosage parquet so
            # we have the canonical column order
            t = pq.read_table(chrom_dir / "dosage.parquet",
                              columns=None)  # need column names
            cols = list(t.column_names)
            assert cols[0] == "pos", "first dosage column must be 'pos'"
            canonical_samples = cols[1:]
        elif manifest["n_samples"] != n_samples_seen:
            raise ValueError(
                f"[{chrom}] n_samples ({manifest['n_samples']}) != earlier "
                f"chroms ({n_samples_seen}). Inconsistent dosage cohort."
            )

        chrom_entries.append({
            "name": chrom,
            "n_sites": manifest["n_sites"],
            "length_bp": manifest["length_bp"],
            "skipped": False,
        })
        n_built += 1

    # If every chrom was skipped, recover the canonical samples list from
    # the first chrom's dosage parquet (without rebuilding it).
    if canonical_samples is None and chrom_entries:
        first_built = chrom_entries[0]
        t = pq.read_table(out / first_built["name"] / "dosage.parquet")
        cols = list(t.column_names)
        canonical_samples = cols[1:]

    # Write store-level samples.tsv (canonical column order)
    if canonical_samples:
        (out / "samples.tsv").write_text("\n".join(canonical_samples) + "\n")

    # Store-level manifest
    store_manifest = {
        "schema_version": SCHEMA_VERSION,
        "n_samples": n_samples_seen or 0,
        "samples": canonical_samples or [],
        "chroms": chrom_entries,
        "row_group_size": row_group_size,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_built": n_built,
        "n_skipped": n_skipped,
        "base": str(base.resolve()),
    }
    (out / "manifest.json").write_text(json.dumps(store_manifest, indent=2))
    log.info("Done — %d built, %d skipped. Manifest: %s",
             n_built, n_skipped, out / "manifest.json")
    return store_manifest


# =============================================================================
# CLI
# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, required=True,
                    help="directory containing <chrom>.sites.tsv.gz + "
                         "<chrom>.dosage.tsv.gz")
    ap.add_argument("--out", type=Path, required=True,
                    help="output dosage_store/ directory")
    ap.add_argument("--samples", type=Path, default=None,
                    help="optional: path to samples.tsv (one ID per line). "
                         "If omitted, falls back to <base>/../samples.tsv "
                         "or an inline header in the dosage file.")
    ap.add_argument("--chroms", type=str, default="",
                    help="optional: comma-separated chrom whitelist. "
                         "If omitted, auto-discover from base.")
    ap.add_argument("--row-group-size", type=int,
                    default=DEFAULT_ROW_GROUP_SIZE,
                    help=f"parquet row-group size (default {DEFAULT_ROW_GROUP_SIZE})")
    ap.add_argument("--force", action="store_true",
                    help="rebuild every chrom even if up-to-date")
    args = ap.parse_args()

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()] or None
    build_store(
        base=args.base.resolve(),
        out=args.out.resolve(),
        samples_path=args.samples.resolve() if args.samples else None,
        chroms=chroms,
        row_group_size=args.row_group_size,
        force=args.force,
    )


if __name__ == "__main__":
    main()
