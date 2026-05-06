"""
region_handler.py — orchestrates store_reader + sampling for the standalone
viewer's /api/region endpoint and its candidate/breakpoint wrappers.

Pure-function module. Tests import directly; the FastAPI server wraps with
HTTPException + JSONResponse.

Returns the rich envelope documented in the S6 spec §Endpoints / GET /api/region:

    {
        "chrom": ...,
        "start_bp": ..., "end_bp": ...,
        "sampling_mode": ...,
        "n_sites_total": ..., "n_sites_returned": ...,
        "downsampled": bool, "warning": str | null,
        "samples": [...],
        "site_ids": [...],
        "positions": [...],
        "missingness": [...], "diagnostic_score": [...],
        "matrix": [[...], ...],   # n_sites x n_samples, int8 with -1=NA
        "_meta": { ... }
    }

For aggregate mode the matrix is float (mean dosage per bin per sample, NaN
where empty), positions are bin centers, and site_ids/missingness are null
(no per-site identity at the aggregate level).
"""
from __future__ import annotations

import math
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from sampling import (
    select_raw, select_even, select_random,
    select_variance, select_hybrid, select_aggregate,
    per_site_variance,
)
from store_reader import (
    load_store_manifest, load_samples, list_chroms,
    read_sites_in_window, read_dosage_for_positions, read_dosage_window,
    get_chrom_length_bp, get_chrom_n_sites,
)

# Defaults — match dosage_bridge.py + spec §Hard caps
DEFAULT_MAX_SITES = 1000
DEFAULT_HARD_CAP = 20000
DEFAULT_MAX_REGION_BP = 50_000_000

ALL_MODES = ("raw", "even", "random", "variance", "hybrid", "aggregate")


# =============================================================================
# Region resolution — supports the standalone /api/region endpoint
# =============================================================================

class RegionError(Exception):
    """Raised by handle_region for cap_exceeded / region_too_large / unknown
    chrom. The server wraps these in HTTPException 400/404."""
    def __init__(self, code: str, **kw):
        self.code = code
        self.detail = {"error": code, **kw}
        super().__init__(code)


def handle_region(*, store_root: Path, chrom: str,
                  start_bp: int, end_bp: int,
                  max_sites: int = DEFAULT_MAX_SITES,
                  mode: str = "hybrid", seed: int = 1,
                  max_region_bp: int = DEFAULT_MAX_REGION_BP,
                  hard_cap: int = DEFAULT_HARD_CAP) -> Dict[str, Any]:
    """Resolve a region request and return the rich envelope.

    Modes raw/even/random/variance/hybrid return per-site indices into the
    store's positions; we then read those positions' dosage rows.

    Mode aggregate returns binned means and is structurally different —
    the matrix is float, positions are bin centers.

    Raises RegionError on:
        - chrom_not_found
        - region_too_large
        - bad_mode
        - cap_exceeded (only when mode='raw')
        - cap_above_hard_cap
    """
    t0 = time.perf_counter()

    if mode not in ALL_MODES:
        raise RegionError("bad_mode", mode=mode, allowed=list(ALL_MODES))
    if max_sites < 1:
        raise RegionError("bad_max_sites", max_sites=max_sites,
                          message="max_sites must be >= 1")
    if max_sites > hard_cap:
        raise RegionError("cap_above_hard_cap", max_sites=max_sites,
                          hard_cap=hard_cap)
    if start_bp < 1 or end_bp < start_bp:
        raise RegionError("bad_range", start_bp=start_bp, end_bp=end_bp)
    if end_bp - start_bp > max_region_bp:
        raise RegionError("region_too_large",
                          n_bp=end_bp - start_bp,
                          max_region_bp=max_region_bp,
                          suggestion=("Use mode=aggregate for chrom-scale "
                                      "or reduce the window."))

    # Validate chrom exists in store
    chroms = list_chroms(store_root)
    if chrom not in chroms:
        raise RegionError("chrom_not_found", chrom=chrom,
                          available=chroms)

    samples = load_samples(store_root)
    n_samples = len(samples)

    # ---------- Aggregate mode is structurally different ------------------
    if mode == "aggregate":
        return _handle_aggregate(
            store_root=store_root, chrom=chrom,
            start_bp=start_bp, end_bp=end_bp,
            max_sites=max_sites, samples=samples, t0=t0,
        )

    # ---------- Per-site modes (raw / even / random / variance / hybrid) -
    t_io = time.perf_counter()
    positions, missingness, diagnostic, site_ids = read_sites_in_window(
        store_root, chrom, start_bp, end_bp
    )
    n_total = len(positions)
    ms_io_sites = int((time.perf_counter() - t_io) * 1000)

    if n_total == 0:
        return _empty_envelope(chrom, start_bp, end_bp, mode, samples,
                                ms_io_sites, msg="no sites in region")

    # ---------- Selection ------------------------------------------------
    t_sel = time.perf_counter()
    sel_indices: List[int]
    downsampled = False
    warning: Optional[str] = None

    if mode == "raw":
        sel_indices, exceeded = select_raw(positions, max_sites)
        if exceeded:
            raise RegionError(
                "raw_cap_exceeded",
                n_sites_total=n_total,
                max_sites=max_sites,
                suggestion=("Use mode=even/random/hybrid/variance/aggregate, "
                            "or raise max_sites (HARD_CAP={}).").format(hard_cap),
            )
        # raw mode: all in window, no downsample
        downsampled = False
    elif mode == "even":
        sel_indices = select_even(positions, max_sites)
        downsampled = (len(sel_indices) < n_total)
        if downsampled:
            warning = (f"even mode: {len(sel_indices)} positionally-uniform "
                       f"sites picked from {n_total} total")
    elif mode == "random":
        sel_indices = select_random(positions, max_sites, seed=seed)
        downsampled = (len(sel_indices) < n_total)
        if downsampled:
            warning = (f"random mode (seed={seed}): {len(sel_indices)} sites "
                       f"picked from {n_total} total")
    elif mode == "variance":
        # Variance needs the dosage matrix; read once.
        _, dosage_full = read_dosage_window(store_root, chrom, start_bp, end_bp)
        variances = per_site_variance(dosage_full)
        sel_indices = select_variance(positions, variances, max_sites)
        downsampled = (len(sel_indices) < n_total)
        if downsampled:
            warning = (f"variance mode: top-{len(sel_indices)} by per-site "
                       f"dosage variance from {n_total} total — biased toward "
                       f"informative sites, not a neutral view")
    elif mode == "hybrid":
        # Hybrid needs variances too
        _, dosage_full = read_dosage_window(store_root, chrom, start_bp, end_bp)
        variances = per_site_variance(dosage_full)
        sel_indices = select_hybrid(positions, variances, max_sites, seed=seed)
        downsampled = (len(sel_indices) < n_total)
        if downsampled:
            warning = (f"hybrid mode: 70% even + 30% variance from "
                       f"{n_total} total — variance share biases toward "
                       f"informative sites")
    else:
        raise RegionError("bad_mode", mode=mode)

    ms_select = int((time.perf_counter() - t_sel) * 1000)

    # ---------- Read dosage for selected sites --------------------------
    t_dosage = time.perf_counter()
    sel_positions = [positions[i] for i in sel_indices]
    if mode in ("variance", "hybrid"):
        # We already read the full window's dosage; subset in-memory rather
        # than hit Parquet again.
        # `dosage_full` aligns row-for-row with `positions` (read_dosage_window
        # returns positions in ascending order, same order as read_sites_in_window).
        matrix = [dosage_full[i] for i in sel_indices]
    else:
        matrix = read_dosage_for_positions(store_root, chrom, sel_positions)
    ms_io_dosage = int((time.perf_counter() - t_dosage) * 1000)

    return {
        "chrom": chrom,
        "start_bp": start_bp,
        "end_bp": end_bp,
        "sampling_mode": mode,
        "n_sites_total": n_total,
        "n_sites_returned": len(sel_indices),
        "downsampled": downsampled,
        "warning": warning,
        "samples": samples,
        "site_ids": [site_ids[i] for i in sel_indices],
        "positions": sel_positions,
        "missingness": [missingness[i] for i in sel_indices],
        "diagnostic_score": [diagnostic[i] for i in sel_indices],
        "matrix": matrix,
        "_meta": {
            "ms_io_sites": ms_io_sites,
            "ms_select": ms_select,
            "ms_io_dosage": ms_io_dosage,
            "n_samples": n_samples,
            "seed": seed if mode in ("random", "hybrid") else None,
            "store_root": str(store_root),
        },
    }


def _empty_envelope(chrom: str, start_bp: int, end_bp: int, mode: str,
                    samples: List[str], ms_io_sites: int,
                    msg: str = "no sites") -> Dict[str, Any]:
    return {
        "chrom": chrom, "start_bp": start_bp, "end_bp": end_bp,
        "sampling_mode": mode,
        "n_sites_total": 0, "n_sites_returned": 0,
        "downsampled": False, "warning": msg,
        "samples": samples,
        "site_ids": [], "positions": [],
        "missingness": [], "diagnostic_score": [],
        "matrix": [],
        "_meta": {"ms_io_sites": ms_io_sites,
                  "ms_select": 0, "ms_io_dosage": 0,
                  "n_samples": len(samples), "seed": None},
    }


def _handle_aggregate(*, store_root: Path, chrom: str,
                       start_bp: int, end_bp: int, max_sites: int,
                       samples: List[str], t0: float) -> Dict[str, Any]:
    """Aggregate mode: bin into max_sites equal-bp bins, mean per bin per sample.

    Returns the same envelope shape but with `matrix` containing floats
    (NaN serialized as None for JSON-friendliness) and positions = bin
    centers. site_ids / missingness / diagnostic_score are null arrays.
    """
    t_io = time.perf_counter()
    positions, dosage_full = read_dosage_window(store_root, chrom, start_bp, end_bp)
    ms_io = int((time.perf_counter() - t_io) * 1000)
    n_total = len(positions)

    if n_total == 0:
        env = _empty_envelope(chrom, start_bp, end_bp, "aggregate", samples,
                              ms_io, msg="no sites in region (aggregate)")
        return env

    t_sel = time.perf_counter()
    centers, mean_matrix = select_aggregate(
        positions, dosage_full, n_samples=len(samples),
        n_bins=max_sites, bp_start=start_bp, bp_end=end_bp,
    )
    ms_sel = int((time.perf_counter() - t_sel) * 1000)

    # NaN -> None for JSON friendliness
    matrix_json: List[List[Optional[float]]] = []
    for row in mean_matrix:
        matrix_json.append([
            None if (v is None or (isinstance(v, float) and math.isnan(v)))
            else float(v)
            for v in row
        ])

    n_returned = len(centers)
    return {
        "chrom": chrom,
        "start_bp": start_bp,
        "end_bp": end_bp,
        "sampling_mode": "aggregate",
        "n_sites_total": n_total,
        "n_sites_returned": n_returned,
        "downsampled": True,
        "warning": (f"aggregate mode: {n_returned} bins from {n_total} sites — "
                    f"mean dosage per bin per sample, NA-aware. Bins are equal-bp."),
        "samples": samples,
        "site_ids": [None] * n_returned,
        "positions": centers,
        "missingness": [None] * n_returned,
        "diagnostic_score": [None] * n_returned,
        "matrix": matrix_json,
        "_meta": {
            "ms_io_dosage": ms_io,
            "ms_select": ms_sel,
            "n_samples": len(samples),
            "store_root": str(store_root),
            "is_aggregate": True,
        },
    }


# =============================================================================
# Manifest endpoint
# =============================================================================

def handle_manifest(store_root: Path) -> Dict[str, Any]:
    """Build the /api/manifest payload from store metadata.

    Returns the rich form documented in the spec §Endpoints / GET /api/manifest.
    """
    m = load_store_manifest(store_root)
    samples = load_samples(store_root)
    chroms_out: List[Dict[str, Any]] = []
    for c in m.get("chroms", []):
        chroms_out.append({
            "name": c["name"],
            "n_sites": int(c.get("n_sites", 0)),
            "length_bp": int(c.get("length_bp", 0)),
        })
    return {
        "schema_version": m.get("schema_version", "dosage_store_v1"),
        "n_samples": len(samples),
        "chroms": chroms_out,
        "samples": samples,
        "candidate_regions_loaded": False,   # set true by the server when CSV exists
        "sample_order_loaded": False,
        "limits": {
            "default_max_sites": DEFAULT_MAX_SITES,
            "hard_cap": DEFAULT_HARD_CAP,
            "max_region_bp": DEFAULT_MAX_REGION_BP,
        },
    }


# =============================================================================
# Candidate / breakpoint regions — config/candidate_regions.tsv
# =============================================================================

def load_candidate_regions(path: Path) -> Dict[str, Dict[str, Any]]:
    """Load candidate_regions.tsv into {candidate_id: {chrom, start, end, ...}}.

    Returns {} if the file is missing (graceful fallback so the server
    can run without candidate config).

    Schema (from spec §candidate_regions.tsv):
        candidate_id, chrom, start, end, [left_breakpoint, right_breakpoint, notes]
    """
    if not path.exists():
        return {}
    out: Dict[str, Dict[str, Any]] = {}
    with open(path, encoding="utf-8") as f:
        header: Optional[List[str]] = None
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cells = line.split("\t")
            if header is None:
                header = [c.strip() for c in cells]
                continue
            row = dict(zip(header, cells))
            cid = row.get("candidate_id", "").strip()
            if not cid:
                continue
            try:
                entry = {
                    "candidate_id": cid,
                    "chrom": row.get("chrom", "").strip(),
                    "start": int(row.get("start", "0")),
                    "end":   int(row.get("end",   "0")),
                }
                if row.get("left_breakpoint", "").strip():
                    entry["left_breakpoint"] = int(row["left_breakpoint"])
                if row.get("right_breakpoint", "").strip():
                    entry["right_breakpoint"] = int(row["right_breakpoint"])
                if row.get("notes"):
                    entry["notes"] = row["notes"]
            except (ValueError, KeyError):
                continue
            out[cid] = entry
    return out


def handle_candidate(*, store_root: Path, candidates: Dict[str, Dict[str, Any]],
                     candidate_id: str, max_sites: int = DEFAULT_MAX_SITES,
                     mode: str = "hybrid", seed: int = 1) -> Dict[str, Any]:
    """Wrap handle_region using a candidate's [chrom, start, end]."""
    cand = candidates.get(candidate_id)
    if not cand:
        raise RegionError("candidate_not_found", candidate_id=candidate_id,
                          available=list(candidates.keys())[:10])
    return handle_region(store_root=store_root, chrom=cand["chrom"],
                         start_bp=cand["start"], end_bp=cand["end"],
                         max_sites=max_sites, mode=mode, seed=seed)


def handle_breakpoint(*, store_root: Path, candidates: Dict[str, Dict[str, Any]],
                      candidate_id: str, side: str,
                      window_bp: int = 500_000,
                      max_sites: int = DEFAULT_MAX_SITES,
                      mode: str = "even", seed: int = 1) -> Dict[str, Any]:
    """Wrap handle_region with ± window_bp around a breakpoint."""
    if side not in ("left", "right"):
        raise RegionError("bad_side", side=side, allowed=["left", "right"])
    cand = candidates.get(candidate_id)
    if not cand:
        raise RegionError("candidate_not_found", candidate_id=candidate_id)
    bp_key = f"{side}_breakpoint"
    if bp_key not in cand:
        # Fallback: use start (left) / end (right) of the candidate.
        bp = cand["start"] if side == "left" else cand["end"]
    else:
        bp = cand[bp_key]
    start = max(1, bp - window_bp)
    end = bp + window_bp
    env = handle_region(store_root=store_root, chrom=cand["chrom"],
                        start_bp=start, end_bp=end,
                        max_sites=max_sites, mode=mode, seed=seed)
    env["_meta"]["breakpoint_bp"] = bp
    env["_meta"]["side"] = side
    env["_meta"]["window_bp"] = window_bp
    return env


def handle_tile(*, store_root: Path, chrom: str, width: int) -> Dict[str, Any]:
    """Whole-chromosome overview at width × n_samples resolution.

    Forces mode=aggregate (no other mode makes sense at chromosome scale).
    """
    if width < 2 or width > DEFAULT_HARD_CAP:
        raise RegionError("bad_tile_width", width=width,
                          allowed_min=2, allowed_max=DEFAULT_HARD_CAP)
    chroms = list_chroms(store_root)
    if chrom not in chroms:
        raise RegionError("chrom_not_found", chrom=chrom, available=chroms)
    length_bp = get_chrom_length_bp(store_root, chrom)
    if length_bp <= 0:
        raise RegionError("chrom_length_unknown", chrom=chrom)
    return handle_region(
        store_root=store_root, chrom=chrom,
        start_bp=1, end_bp=length_bp,
        max_sites=width, mode="aggregate",
        max_region_bp=length_bp + 1,   # tile spans the whole chrom by design
    )
