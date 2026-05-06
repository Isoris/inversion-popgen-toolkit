"""
store_reader.py — read the Parquet dosage store built by 01_prepare_dosage_store.py.

Pure-function module. The standalone viewer's FastAPI server (02_run_server.py)
imports these helpers; tests import them directly without booting FastAPI.

The store layout (recap from 01_prepare_dosage_store.py):
    <store>/manifest.json
    <store>/samples.tsv              # canonical sample order (one ID per line)
    <store>/<chrom>/sites.parquet    # cols: pos, ref, alt, missingness, diagnostic_score, site_id
    <store>/<chrom>/dosage.parquet   # cols: pos, s0, s1, ..., s{n-1}; int8; -1=NA
    <store>/<chrom>/_chrom_manifest.json

Why two parquet files per chrom? Quoting the spec §Layout on disk:
    "site metadata is hot (every region query reads it); the dosage matrix
     is cold (only fetched after sites have been picked). Splitting them
     halves the bytes read for sample queries."

This reader honours that pattern:
    1. read_sites(chrom, start, end)        — small, hot
    2. read_dosage_for_positions(chrom, ps) — big, cold; only after selection

read_dosage_for_positions reads the FULL dosage table for the chrom (no
position filtering at parquet layer for now) and slices it in memory. A
future optimisation can use parquet predicate pushdown to skip row groups,
but for the cohort sizes Quentin is targeting (~226 samples × 5 M sites)
the full read is ~5 MB int8 — fine.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional, Tuple

import pyarrow as pa
import pyarrow.parquet as pq


SCHEMA_VERSION = "dosage_store_v1"


# =============================================================================
# Manifest + samples
# =============================================================================

def load_store_manifest(store_root: Path) -> dict:
    """Load <store>/manifest.json. Raises FileNotFoundError if missing."""
    p = store_root / "manifest.json"
    if not p.exists():
        raise FileNotFoundError(f"store manifest not found: {p}")
    with open(p, encoding="utf-8") as f:
        m = json.load(f)
    if m.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"manifest schema_version mismatch: expected {SCHEMA_VERSION!r}, "
            f"got {m.get('schema_version')!r}"
        )
    return m


def load_samples(store_root: Path) -> List[str]:
    """Load <store>/samples.tsv (one ID per line)."""
    p = store_root / "samples.tsv"
    if not p.exists():
        raise FileNotFoundError(f"samples.tsv not found: {p}")
    out: List[str] = []
    with open(p, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(line)
    return out


def list_chroms(store_root: Path) -> List[str]:
    """Return chroms found in the manifest, in manifest order."""
    m = load_store_manifest(store_root)
    return [c["name"] for c in m.get("chroms", [])]


# =============================================================================
# Site-level reads (hot path)
# =============================================================================

def read_sites_in_window(store_root: Path, chrom: str,
                         start_bp: int, end_bp: int
                         ) -> Tuple[List[int], List[Optional[float]],
                                    List[Optional[float]], List[str]]:
    """Read sites.parquet for `chrom`, filtered to [start_bp, end_bp].

    Returns (positions, missingness, diagnostic_score, site_ids), each a
    Python list. NaN floats round-trip as None so the JSON envelope can
    show them as null (matches the renderer contract).
    """
    p = store_root / chrom / "sites.parquet"
    if not p.exists():
        raise FileNotFoundError(f"sites.parquet not found for chrom={chrom}: {p}")
    # Predicate pushdown: pyarrow can filter by row-group min/max stats
    # if we hand it a pa.compute expression. For simplicity we read the
    # whole columns and slice — the sites file is small (60 bytes/site
    # at 6 cols × 10 bytes; 5M sites = 300 MB but most regions are < 5 Mb).
    # If perf becomes a problem, switch to:
    #   import pyarrow.compute as pc
    #   filt = (pc.field("pos") >= start_bp) & (pc.field("pos") <= end_bp)
    #   t = pq.read_table(p, filters=filt)
    t = pq.read_table(p, columns=["pos", "missingness",
                                  "diagnostic_score", "site_id"])
    pos_col = t.column("pos").to_pylist()
    miss_col = t.column("missingness").to_pylist()
    diag_col = t.column("diagnostic_score").to_pylist()
    sid_col = t.column("site_id").to_pylist()

    positions: List[int] = []
    missingness: List[Optional[float]] = []
    diagnostic: List[Optional[float]] = []
    site_ids: List[str] = []
    for i, p_bp in enumerate(pos_col):
        if p_bp < start_bp:
            continue
        if p_bp > end_bp:
            break    # sorted ascending; we can stop early
        positions.append(int(p_bp))
        # NaN → None
        m = miss_col[i]
        missingness.append(None if (m is None or _is_nan(m)) else float(m))
        d = diag_col[i]
        diagnostic.append(None if (d is None or _is_nan(d)) else float(d))
        site_ids.append(sid_col[i] or "")
    return positions, missingness, diagnostic, site_ids


def _is_nan(x) -> bool:
    """Tolerant NaN check that returns False for non-floats."""
    try:
        return x != x   # standard NaN check: NaN != NaN
    except Exception:
        return False


# =============================================================================
# Dosage reads (cold path)
# =============================================================================

def read_dosage_for_positions(store_root: Path, chrom: str,
                              positions: List[int]) -> List[List[int]]:
    """Read dosage.parquet rows whose `pos` is in `positions` (ascending).

    Returns a list of dosage rows, one per requested position, in the
    SAME ORDER as `positions`. Each row is a list of int (length = n_samples,
    -1 for NA).

    `positions` MUST be ascending. If a position isn't found in the
    parquet (e.g. caller passed a stale list), an all-NA row is returned
    for that slot (defensive — same convention as the bridge).
    """
    if not positions:
        return []
    p = store_root / chrom / "dosage.parquet"
    if not p.exists():
        raise FileNotFoundError(f"dosage.parquet not found for chrom={chrom}: {p}")

    t = pq.read_table(p)
    pos_col = t.column("pos").to_pylist()
    # Build a position -> row index map. positions is ascending in the
    # parquet (we sorted in 01_prepare). Linear scan is fine.
    pos_to_idx = {pos: i for i, pos in enumerate(pos_col)}

    sample_cols = [c for c in t.column_names if c != "pos"]
    n_samples = len(sample_cols)
    # Pull each sample column once; iterate over the requested positions
    # and gather the cells.
    sample_arrays: List[List[int]] = []
    for c in sample_cols:
        sample_arrays.append(t.column(c).to_pylist())

    rows: List[List[int]] = []
    for req_pos in positions:
        idx = pos_to_idx.get(req_pos)
        if idx is None:
            rows.append([-1] * n_samples)
            continue
        row = [int(sample_arrays[j][idx]) for j in range(n_samples)]
        rows.append(row)
    return rows


def read_dosage_window(store_root: Path, chrom: str,
                       start_bp: int, end_bp: int
                       ) -> Tuple[List[int], List[List[int]]]:
    """Read all dosage rows in [start_bp, end_bp], returning (positions, rows).

    Convenience wrapper used by aggregate mode and variance computation —
    they need every site in the window, not a pre-filtered subset.
    """
    p = store_root / chrom / "dosage.parquet"
    if not p.exists():
        raise FileNotFoundError(f"dosage.parquet not found for chrom={chrom}: {p}")
    t = pq.read_table(p)
    pos_col = t.column("pos").to_pylist()
    sample_cols = [c for c in t.column_names if c != "pos"]
    sample_arrays = [t.column(c).to_pylist() for c in sample_cols]
    positions: List[int] = []
    rows: List[List[int]] = []
    for i, p_bp in enumerate(pos_col):
        if p_bp < start_bp:
            continue
        if p_bp > end_bp:
            break
        positions.append(int(p_bp))
        rows.append([int(sample_arrays[j][i]) for j in range(len(sample_cols))])
    return positions, rows


# =============================================================================
# Helpers
# =============================================================================

def get_chrom_length_bp(store_root: Path, chrom: str) -> int:
    """Return the length_bp for `chrom` from its per-chrom manifest. 0 on error."""
    p = store_root / chrom / "_chrom_manifest.json"
    if not p.exists():
        return 0
    try:
        with open(p, encoding="utf-8") as f:
            m = json.load(f)
        return int(m.get("length_bp", 0))
    except Exception:
        return 0


def get_chrom_n_sites(store_root: Path, chrom: str) -> int:
    """Return n_sites for `chrom` from its per-chrom manifest. 0 on error."""
    p = store_root / chrom / "_chrom_manifest.json"
    if not p.exists():
        return 0
    try:
        with open(p, encoding="utf-8") as f:
            m = json.load(f)
        return int(m.get("n_sites", 0))
    except Exception:
        return 0
