"""
dosage_bridge.py — atlas-side bridge for the chromosome-wide dosage heatmap.

Endpoint: GET /api/dosage/chunk

Request (query params):
    chrom        str    e.g. "C_gar_LG28"   (required)
    start        int    1-based bp inclusive (required, >=1)
    end          int    1-based bp inclusive (required, >=start)
    cap          int    max sites returned (default 500, hard cap 20000)
    mode         str    "raw" (default) | "even" — controls site selection

Response (JSON, the renderer-compatible chunk shape):
    chrom
    start_bp
    end_bp
    samples       [str, ...]                       length = n_samples (string IDs)
    markers       [{pos_bp, missingness, diagnostic_score}, ...]
                                                    length = n_sites
    dosage        [[int, ...], ...]                 n_sites x n_samples; -1 = NA
    _meta         {ms_io, ms_select, n_sites_total, ...}

Why this exists
---------------
The atlas's `renderDosageHeatmap` (see Inversion_atlas.html ~line 15306) reads
chunks from `state.data.dosage_chunks.chunks[i].url`. Without a server, the
only way to populate that index is to bake static JSON chunks per chromosome
ahead of time. The bridge exposes one HTTP endpoint per chromosome that the
atlas treats as a "synthetic chunk", with templated URL parameters that get
substituted at fetch time.

Source data layout (from S6_dosage_heatmap_streaming_viewer.md §Source data):
    <base>/04_dosage_by_chr/
        <chrom>.sites.tsv.gz       (one row per site, sorted by pos)
        <chrom>.dosage.tsv.gz      (same row count as sites; one row per site)
    <base>/samples.tsv             (one ID per line, in dosage column order)

Sites TSV columns:
    1: pos     (int, 1-based bp)
    2: ref     (str)
    3: alt     (str)
    4: missingness        (float; "NA" allowed)
    5: diagnostic_score   (float; "NA" allowed)
    6: site_id            (str; falls back to "<chrom>:<pos>:<ref>><alt>")

Dosage TSV: n_sites rows, n_samples cells per row, values in
    {0, 1, 2, NA, ., -1, ""} → normalized to int8 with -1 = NA.

Sampling modes
--------------
- raw:  return all sites in [start, end]. If n_sites_total > cap, returns
        HTTP 400 with {"error": "raw_cap_exceeded", "n_sites_total": ...,
        "max_sites": ..., "suggestion": ...}. Per the spec, no silent
        auto-switch — keeps "raw means all" honest.
- even: pick `cap` sites evenly across [start, end] by genomic position.
        Use this when raw_cap_exceeded would otherwise fire.

Other modes (random, variance, hybrid, aggregate) belong to the standalone
viewer's richer endpoint, not the bridge. The bridge is deliberately minimal
because the atlas-side `selectTopMarkers` already does its own top-N work.

Renderer contract — DO NOT BREAK
--------------------------------
- dosage is sites × samples, NOT samples × sites
- missing is encoded as -1, NOT null (see Inversion_atlas.html line ~14602)
- chunk.samples[i] are STRING IDs, resolved against state.data.samples[].id
- chunk.dosage[i].length === chunk.samples.length for every i
- chunk.markers.length === chunk.dosage.length

These four facts come from reverse-engineering the live renderer in
`specs_todo/from_turn129/S6_dosage_heatmap_streaming_viewer.md`.
"""
from __future__ import annotations

import gzip
import io
import logging
import math
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from fastapi import HTTPException
from pydantic import BaseModel, Field, field_validator

log = logging.getLogger(__name__)


# =============================================================================
# Defaults & limits — defaults match the spec's §Hard caps
# =============================================================================

DEFAULT_MAX_SITES = 1000
DEFAULT_HARD_CAP = 20000
DEFAULT_MAX_REGION_BP = 50_000_000


# =============================================================================
# Request validation
# =============================================================================

class DosageChunkReq(BaseModel):
    """Validated query parameters for /api/dosage/chunk."""
    chrom: str = Field(..., min_length=1, max_length=128)
    start: int = Field(..., ge=1)
    end: int = Field(..., ge=1)
    cap: int = Field(DEFAULT_MAX_SITES, ge=1, le=DEFAULT_HARD_CAP)
    mode: str = Field("raw", pattern=r"^(raw|even)$")

    @field_validator("end")
    @classmethod
    def _end_ge_start(cls, v: int, info) -> int:
        # `info.data` carries already-validated fields. Pydantic validates
        # in declaration order so `start` is present here.
        start = info.data.get("start") if info and getattr(info, "data", None) else None
        if start is not None and v < start:
            raise ValueError("end must be >= start")
        return v

    @field_validator("chrom")
    @classmethod
    def _chrom_safe(cls, v: str) -> str:
        # Tight whitelist — chrom values are filenames, so we reject path
        # separators, parent-dir traversal, and shell metacharacters.
        bad = set("/\\.\0 \t\n\r\"'`;|&$<>(){}*?[]")
        if any(ch in bad for ch in v):
            raise ValueError(f"chrom contains forbidden chars: {v!r}")
        return v


# =============================================================================
# Helpers
# =============================================================================

def _maybe_float(s: str) -> Optional[float]:
    """Parse a TSV float cell, returning None for NA / empty / non-numeric.

    Matches the renderer convention: NA → None on the parser side, then
    propagated as null in JSON for `missingness` and `diagnostic_score`
    (these two fields ARE allowed to be null per the renderer contract;
    only the dosage matrix uses -1 for missing).
    """
    s = (s or "").strip()
    if not s or s == "." or s.upper() == "NA":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _maybe_dosage_int(s: str) -> int:
    """Parse a dosage cell into int8-range, normalizing NA to -1.

    Accepts: 0, 1, 2, -1, NA, ., '' → -1 for any missing/non-numeric cell.
    Out-of-range numeric values are also treated as missing (defensive;
    a real cohort will only ever have 0/1/2/NA).
    """
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
    # Any other value (e.g. 3 or 99) is treated as missing rather than
    # leaking through to the renderer where v < 0 is the missing test.
    return -1


def _open_gz_text(path: Path) -> io.TextIOWrapper:
    """Open a gzipped TSV as a text iterator."""
    return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")


def _detect_dosage_header(first_line: str, n_samples_expected: int) -> bool:
    """Return True if `first_line` looks like a dosage-header row of sample IDs.

    Rule (from §Source data layout): if line 1 has the same field count as
    samples.tsv AND parses to integers, it's data; otherwise it's a header.
    Since we only have the dosage cell-parser in this bridge, we use a
    simpler rule: if every field after stripping is in {'', '.', 'NA',
    integer}, it's data; if any field is a non-numeric string, treat as
    header.
    """
    fields = first_line.rstrip("\n").split("\t")
    if len(fields) != n_samples_expected:
        # Mismatched count: treat as data and let the loader complain
        # downstream — covering the case where samples.tsv has been
        # updated without re-emitting the dosage TSV.
        return False
    for f in fields:
        s = f.strip()
        if not s or s == "." or s.upper() == "NA":
            continue
        # If we can't parse to int, it's a header field
        try:
            int(s)
        except ValueError:
            return True
    return False


def _read_samples(path: Path) -> List[str]:
    """Read one-sample-ID-per-line file. Strips blanks and comments."""
    out: List[str] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Some pipelines emit "<sample>\t<group>" — take the first column.
            if "\t" in line:
                line = line.split("\t", 1)[0]
            out.append(line)
    return out


def _select_even(positions: List[int], cap: int) -> List[int]:
    """Pick `cap` indices from `positions` evenly across the bp range.

    Returns indices into `positions` (sorted ascending). When n <= cap
    returns range(n).

    "Even" here means equally spaced in genomic-bp space (not equally
    spaced in array-index space) — the spec's §Sampling modes table
    explicitly distinguishes this. We implement bp-even by binning the
    [start, end] range into `cap` bins and picking the closest site to
    each bin's midpoint. Ties broken toward the lower index.
    """
    n = len(positions)
    if n <= cap:
        return list(range(n))
    if cap <= 0:
        return []
    lo = positions[0]
    hi = positions[-1]
    if hi <= lo:
        # Degenerate: all sites at the same bp; return first cap by index
        return list(range(cap))
    span = hi - lo
    out: List[int] = []
    last_idx = -1
    for k in range(cap):
        # Midpoint of bin k
        target = lo + (k + 0.5) * span / cap
        # Binary search for the closest position
        # (positions is sorted ascending)
        idx = _closest_index(positions, target)
        # Avoid duplicates: if we picked the same index twice, advance.
        # This can happen when `cap` exceeds the unique-position density.
        if idx == last_idx:
            idx = idx + 1 if idx + 1 < n else idx
        out.append(idx)
        last_idx = idx
    # Dedup-and-sort (stable)
    seen = set()
    dedup: List[int] = []
    for i in out:
        if i not in seen and 0 <= i < n:
            seen.add(i)
            dedup.append(i)
    dedup.sort()
    return dedup


def _closest_index(sorted_vals: List[int], target: float) -> int:
    """Index of the value in `sorted_vals` closest to `target`."""
    n = len(sorted_vals)
    if n == 0:
        return 0
    # Binary search for insertion point
    lo, hi = 0, n
    while lo < hi:
        mid = (lo + hi) // 2
        if sorted_vals[mid] < target:
            lo = mid + 1
        else:
            hi = mid
    if lo == 0:
        return 0
    if lo == n:
        return n - 1
    # Compare lo-1 and lo
    a = abs(sorted_vals[lo - 1] - target)
    b = abs(sorted_vals[lo] - target)
    return (lo - 1) if a <= b else lo


# =============================================================================
# Data loaders
# =============================================================================

def _load_sites_in_window(sites_path: Path, start: int, end: int
                          ) -> Tuple[List[int], List[Optional[float]],
                                     List[Optional[float]], List[int]]:
    """Stream <chrom>.sites.tsv.gz, keeping rows where pos in [start, end].

    Returns (positions, missingness, diagnostic_score, line_indices).
    `line_indices` records the row index in the sites file (0-based, with
    any header skipped) so the dosage loader can fetch matching rows.

    Tolerates:
      - leading "#"-comment lines
      - an optional header row whose first cell is non-numeric
      - missingness/diagnostic_score columns absent (filled with None)
      - site_id column absent (we don't return it from this helper —
        the bridge renderer doesn't need it; selectTopMarkers reads
        pos_bp / missingness / diagnostic_score)
    """
    positions: List[int] = []
    missingness: List[Optional[float]] = []
    diagnostic: List[Optional[float]] = []
    line_indices: List[int] = []

    row_idx = -1   # Tracks logical row index (excluding comments + header)
    with _open_gz_text(sites_path) as f:
        for raw in f:
            if not raw or raw.startswith("#"):
                continue
            line = raw.rstrip("\n")
            cells = line.split("\t")
            if not cells or not cells[0]:
                continue
            # Detect a header row: first column non-numeric
            try:
                pos = int(cells[0])
            except ValueError:
                # Header — skip without bumping row_idx
                continue
            row_idx += 1
            if pos < start:
                continue
            if pos > end:
                # Sites are sorted by pos — we can stop early
                break
            miss = _maybe_float(cells[3]) if len(cells) > 3 else None
            ds = _maybe_float(cells[4]) if len(cells) > 4 else None
            positions.append(pos)
            missingness.append(miss)
            diagnostic.append(ds)
            line_indices.append(row_idx)
    return positions, missingness, diagnostic, line_indices


def _load_dosage_rows(dosage_path: Path, line_indices_sorted: List[int],
                      n_samples_expected: int) -> List[List[int]]:
    """Stream <chrom>.dosage.tsv.gz, returning rows at the requested indices.

    `line_indices_sorted` MUST be ascending — the loader streams forward
    and never seeks back. Each kept row is parsed into a list[int] of
    length n_samples_expected.

    Detects + skips an optional header row whose values are non-numeric.
    """
    if not line_indices_sorted:
        return []
    rows: Dict[int, List[int]] = {}
    wanted = set(line_indices_sorted)
    next_wanted_idx = 0
    next_wanted = line_indices_sorted[next_wanted_idx]

    row_idx = -1   # logical row index; -1 = haven't seen first data row yet
    header_decided = False
    with _open_gz_text(dosage_path) as f:
        for raw in f:
            if not raw or raw.startswith("#"):
                continue
            line = raw.rstrip("\n")
            if not header_decided:
                # First non-comment line — could be header or data
                if _detect_dosage_header(line, n_samples_expected):
                    header_decided = True
                    continue
                header_decided = True
                # Fall through to parse this line as data
            row_idx += 1
            if row_idx < next_wanted:
                continue
            if row_idx == next_wanted:
                cells = line.split("\t")
                # Don't fail hard if the row is short — pad with -1 (NA).
                if len(cells) < n_samples_expected:
                    cells = list(cells) + [""] * (n_samples_expected - len(cells))
                elif len(cells) > n_samples_expected:
                    cells = cells[:n_samples_expected]
                rows[row_idx] = [_maybe_dosage_int(c) for c in cells]
                next_wanted_idx += 1
                if next_wanted_idx >= len(line_indices_sorted):
                    break
                next_wanted = line_indices_sorted[next_wanted_idx]

    return [rows[i] for i in line_indices_sorted if i in rows]


# =============================================================================
# Top-level handler
# =============================================================================

def handle_dosage_chunk(req: DosageChunkReq, *,
                        dosage_dir: Path,
                        samples_path: Path,
                        max_region_bp: int = DEFAULT_MAX_REGION_BP) -> Dict[str, Any]:
    """Resolve one /api/dosage/chunk request and return the chunk payload.

    Pure function — no caching (the atlas already caches by URL+region in
    `_fetchAndCacheChunk`). Suitable for direct unit testing.
    """
    t0 = time.perf_counter()

    # Region-size guard (documented as MAX_REGION_BP in the spec). The bridge
    # only does raw + even, neither of which makes sense at chromosome scale.
    if req.end - req.start > max_region_bp:
        raise HTTPException(status_code=400, detail={
            "error": "region_too_large",
            "n_bp": req.end - req.start,
            "max_region_bp": max_region_bp,
            "suggestion": "Use the standalone viewer's mode=aggregate for "
                          "regions larger than max_region_bp.",
        })

    sites_path = dosage_dir / f"{req.chrom}.sites.tsv.gz"
    dos_path = dosage_dir / f"{req.chrom}.dosage.tsv.gz"
    if not sites_path.exists():
        raise HTTPException(status_code=404, detail={
            "error": "chrom_not_found",
            "chrom": req.chrom,
            "expected_file": str(sites_path),
        })
    if not dos_path.exists():
        raise HTTPException(status_code=404, detail={
            "error": "chrom_dosage_missing",
            "chrom": req.chrom,
            "expected_file": str(dos_path),
        })
    if not samples_path.exists():
        raise HTTPException(status_code=500, detail={
            "error": "samples_file_missing",
            "expected_file": str(samples_path),
        })

    samples = _read_samples(samples_path)
    if not samples:
        raise HTTPException(status_code=500, detail={
            "error": "samples_file_empty",
            "expected_file": str(samples_path),
        })

    # Sites in window
    positions, missingness, diagnostic, line_indices = _load_sites_in_window(
        sites_path, req.start, req.end
    )
    n_total = len(positions)
    t_io_sites = time.perf_counter() - t0

    # Selection
    t1 = time.perf_counter()
    if n_total == 0:
        sel_indices: List[int] = []
        downsampled = False
        warning: Optional[str] = "no sites in region"
    elif n_total <= req.cap:
        sel_indices = list(range(n_total))
        downsampled = False
        warning = None
    elif req.mode == "raw":
        # Raw cap policy: 400 with structured error, no auto-switch.
        raise HTTPException(status_code=400, detail={
            "error": "raw_cap_exceeded",
            "n_sites_total": n_total,
            "max_sites": req.cap,
            "suggestion": ("Use mode=even or raise cap (HARD_CAP={}). "
                           "Or use the standalone viewer for variance / "
                           "hybrid / aggregate modes.").format(DEFAULT_HARD_CAP),
        })
    else:  # mode == "even"
        sel_indices = _select_even(positions, req.cap)
        downsampled = True
        warning = ("even mode: {} positionally-uniform sites picked "
                   "from {} total").format(len(sel_indices), n_total)

    t_select = time.perf_counter() - t1

    # Dosage rows for selected sites
    t2 = time.perf_counter()
    sel_line_indices = [line_indices[i] for i in sel_indices]
    dosage_rows = _load_dosage_rows(dos_path, sel_line_indices, len(samples))
    t_io_dosage = time.perf_counter() - t2

    # If a dosage row didn't get loaded (e.g. truncated file), pad with all-NA
    while len(dosage_rows) < len(sel_indices):
        dosage_rows.append([-1] * len(samples))

    # Build response
    markers = []
    for i in sel_indices:
        markers.append({
            "pos_bp": positions[i],
            "missingness": missingness[i],
            "diagnostic_score": diagnostic[i],
        })

    return {
        "chrom": req.chrom,
        "start_bp": req.start,
        "end_bp": req.end,
        "samples": samples,
        "markers": markers,
        "dosage": dosage_rows,
        "_meta": {
            "ms_io_sites": int(t_io_sites * 1000),
            "ms_select": int(t_select * 1000),
            "ms_io_dosage": int(t_io_dosage * 1000),
            "n_sites_total": n_total,
            "n_sites_returned": len(sel_indices),
            "n_samples": len(samples),
            "downsampled": downsampled,
            "warning": warning,
            "mode": req.mode,
        },
    }
