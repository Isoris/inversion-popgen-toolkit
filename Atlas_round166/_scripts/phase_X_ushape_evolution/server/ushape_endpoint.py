# =============================================================================
# server/ushape_endpoint.py
# =============================================================================
# Drop-in FastAPI module for popstats_server.py providing:
#
#     POST /api/ushape/candidate   compute U-shape evolution layer for ONE candidate
#
# How it integrates with the rest of popstats_server.py:
#
#     - Reuses _ensure_ready, _validate_groups_shape, SAMPLES, CACHE, ENGINES,
#       resolve_beagle, run_region_popstats from the host module.
#     - Calls region_popstats once over [start - flank, end + flank] with two
#       groups (HOMO_1 + HOMO_2). HET samples are dropped from the compute
#       (they don't enter dXY/FST between arrangements).
#     - Hands the per-window output to an Rscript subprocess (R/ushape_lib.R +
#       R/ushape_classify.R + R/ushape_io.R) which returns the candidate JSON
#       block on stdout.
#     - Caches the result keyed on (candidate_id, region, groups, params,
#       region_popstats engine hash, ushape_lib hash).
#
# Wiring in popstats_server.py — add at the bottom, near the other endpoints:
#
#     from ushape_endpoint import UshapeReq, ushape_candidate as _ushape_handler
#     @app.post("/api/ushape/candidate")
#     async def ushape_candidate(req: UshapeReq) -> Response:
#         return await _ushape_handler(
#             req,
#             ensure_ready=_ensure_ready, samples=SAMPLES, cache=CACHE,
#             engines=ENGINES, run_region_popstats=run_region_popstats,
#             resolve_beagle=resolve_beagle, beagle_dir=Path(CFG["beagle_dir"]),
#             ushape_r_dir=Path(CFG.get("ushape_r_dir", "R")),
#         )
#
# Atlas-side fetch (mirrors popstats_groupwise call):
#
#     POST {server}/api/ushape/candidate
#     Content-Type: application/json
#     {
#       "candidate_id": "LG28_INV_001",
#       "chrom": "C_gar_LG28",
#       "start_bp": 12345678, "end_bp": 15234567,
#       "left_breakpoint":  12345678, "right_breakpoint": 15234567,
#       "groups": {
#         "HOMO_1": ["sample001", ...],
#         "HOMO_2": ["sample002", ...],
#         "HET":    ["sample003", ...]   # optional, dropped from compute
#       },
#       "params": {  # all optional
#         "win_bp": 50000, "step_bp": 10000,
#         "n_windows_inside": 100, "max_flank_bp": 2000000,
#         "min_flank_bp": 100000, "min_window_bp": 1000,
#         "edge_fraction": 0.20, "center_fraction": 0.60,
#         "u_score_high": 1.5, "inside_flank_high": 1.5,
#         "internal_peak_high": 1.5, "asymmetry_log2_high": 1.0,
#         "fst_enrichment_high": 1.5, "oldness_min": 1.2,
#         "dxy_min_inside": 1e-4, "flatness_max_for_flat": 0.35
#       }
#     }
#
# Response: ushape_evolution_v1 candidate block as inline JSON
# (single-candidate; the atlas merges into its layer cache).
# =============================================================================

from __future__ import annotations

import asyncio
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Awaitable, Callable, Dict, List, Optional

from fastapi import HTTPException
from fastapi.responses import JSONResponse, Response
from pydantic import BaseModel, Field, field_validator


_CHROM_RE = re.compile(r"^[A-Za-z0-9_.\-]{1,80}$")


# -----------------------------------------------------------------------------
# request schema
# -----------------------------------------------------------------------------
class UshapeParams(BaseModel):
    win_bp: int = 50_000
    step_bp: int = 10_000
    n_windows_inside: int = 100
    max_flank_bp: int = 2_000_000
    min_flank_bp: int = 100_000
    min_window_bp: int = 1_000
    edge_fraction: float = 0.20
    center_fraction: float = 0.60
    u_score_high: float = 1.5
    inside_flank_high: float = 1.5
    internal_peak_high: float = 1.5
    asymmetry_log2_high: float = 1.0
    fst_enrichment_high: float = 1.5
    oldness_min: float = 1.2
    dxy_min_inside: float = 1e-4
    flatness_max_for_flat: float = 0.35
    min_inside_windows_for_shape: int = 8
    n_random_background: int = 1000
    matched_seed: int = 42


class UshapeReq(BaseModel):
    candidate_id: str
    chrom: str
    start_bp: int
    end_bp: int
    left_breakpoint: Optional[int] = None
    right_breakpoint: Optional[int] = None
    groups: Dict[str, List[str]]
    params: UshapeParams = Field(default_factory=UshapeParams)

    @field_validator("chrom")
    @classmethod
    def _chrom_ok(cls, v: str) -> str:
        if not _CHROM_RE.match(v):
            raise ValueError("chrom contains illegal characters")
        return v

    @field_validator("end_bp")
    @classmethod
    def _ordered(cls, v: int, info) -> int:
        s = info.data.get("start_bp")
        if s is not None and v <= s:
            raise ValueError("end_bp must be > start_bp")
        return v


# -----------------------------------------------------------------------------
# cache key
# -----------------------------------------------------------------------------
def _ushape_cache_key(*, candidate_id: str, chrom: str, start: int, end: int,
                      groups: Dict[str, List[str]], params: Dict[str, Any],
                      engines: Dict[str, str], r_lib_hash: str) -> str:
    payload = {
        "candidate_id": candidate_id,
        "chrom": chrom, "start": start, "end": end,
        "groups": {g: sorted(m) for g, m in groups.items()},
        "params": params,
        "engines": engines,
        "r_lib_hash": r_lib_hash,
    }
    blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def _hash_r_dir(r_dir: Path) -> str:
    h = hashlib.sha256()
    for f in sorted(r_dir.glob("*.R")):
        h.update(f.name.encode())
        h.update(f.read_bytes())
    return h.hexdigest()[:16]


# -----------------------------------------------------------------------------
# Rscript subprocess
# -----------------------------------------------------------------------------
_R_DRIVER = r"""
suppressPackageStartupMessages({
  library(jsonlite); library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
ushape_r_dir <- args[1]
in_path  <- args[2]
out_path <- args[3]

source(file.path(ushape_r_dir, "ushape_lib.R"))
source(file.path(ushape_r_dir, "ushape_classify.R"))
source(file.path(ushape_r_dir, "ushape_io.R"))

req <- fromJSON(in_path, simplifyVector = FALSE)
candidate <- req$candidate
params    <- req$params
popstats_payload <- req$popstats_payload
groups    <- req$groups

win <- compute_window_stats_from_popstats_json(popstats_payload, candidate, params)
if (nrow(win) == 0L) {
  write(toJSON(list(error = "empty_window_stats"), auto_unbox = TRUE), out_path)
  quit(status = 0)
}

raw <- candidate_raw_summary(win, candidate, groups)
sc  <- score_candidate(raw, win, params)
inside_n <- nrow(win[zone %in% c("left_edge","center","right_edge")])
cl  <- classify_candidate(sc, raw, params, inside_n)
classification <- list(
  primary_class = cl$primary_class, secondary_class = cl$secondary_class,
  confidence = cl$confidence, reason = cl$reason, flags = cl$flags
)

# matched-background z-scores from local flank resampling — same helper
# the offline U03 uses, so atlas-fetched and batch-exported numbers agree.
mbg <- flank_resample_zscores(win,
                              n_random = params$n_random_background %||% 1000L,
                              seed     = params$matched_seed         %||% 42L)

block <- ushape_candidate_block(win, raw, sc, classification, groups,
                                matched_bg = mbg)
write(toJSON(block, auto_unbox = TRUE, na = "null", digits = 6), out_path)
"""


async def _run_r_driver(*, ushape_r_dir: Path, payload: dict) -> dict:
    """Spawn Rscript with the in-memory driver script; return the parsed block."""
    with tempfile.TemporaryDirectory(prefix="ushape_") as td:
        tdp = Path(td)
        in_path = tdp / "in.json"
        out_path = tdp / "out.json"
        drv_path = tdp / "drive.R"
        in_path.write_text(json.dumps(payload, default=str))
        drv_path.write_text(_R_DRIVER)

        cmd = ["Rscript", "--vanilla", str(drv_path),
               str(ushape_r_dir), str(in_path), str(out_path)]
        proc = await asyncio.create_subprocess_exec(
            *cmd, stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        out, err = await proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f"Rscript failed: rc={proc.returncode}\n"
                               f"stderr:\n{err.decode(errors='replace')}")
        return json.loads(out_path.read_text())


# -----------------------------------------------------------------------------
# main handler — invoked from popstats_server.py
# -----------------------------------------------------------------------------
async def ushape_candidate(
    req: UshapeReq,
    *,
    ensure_ready: Callable[[], None],
    samples,            # SamplesRegistry
    cache,              # CacheLayer
    engines,            # EnginesRegistry
    run_region_popstats: Callable[..., Awaitable[Any]],
    resolve_beagle: Callable[[Path, str], Path],
    beagle_dir: Path,
    ushape_r_dir: Path,
    min_group_n: int = 3,
) -> Response:
    ensure_ready()

    # Drop HET from the compute — dXY/FST/π are arrangement-vs-arrangement.
    compute_groups = {
        "HOMO_1": list(req.groups.get("HOMO_1", [])),
        "HOMO_2": list(req.groups.get("HOMO_2", [])),
    }
    if (len(compute_groups["HOMO_1"]) < min_group_n or
        len(compute_groups["HOMO_2"]) < min_group_n):
        raise HTTPException(400,
            f"HOMO_1/HOMO_2 each need >= {min_group_n} samples")

    cleaned: Dict[str, List[str]] = {}
    for g, m in compute_groups.items():
        kept, dropped = samples.filter_known(m)
        if len(kept) < min_group_n:
            raise HTTPException(400,
                f"group '{g}' has only {len(kept)} known samples after "
                f"filtering against canonical sample_list (min={min_group_n}). "
                f"Dropped (unknown): {dropped[:5]}")
        cleaned[g] = kept

    # Region: candidate ± flank. Cap at chr_len when known (host server
    # already validates that for region_popstats).
    L = req.end_bp - req.start_bp + 1
    flank = max(req.params.min_flank_bp,
                min(req.params.max_flank_bp, L // 2))
    region_lo = max(1, req.start_bp - flank)
    region_hi = req.end_bp + flank

    # cache
    engine_hashes = engines.all_hashes()
    r_lib_hash = _hash_r_dir(ushape_r_dir)
    key = _ushape_cache_key(
        candidate_id=req.candidate_id, chrom=req.chrom,
        start=region_lo, end=region_hi,
        groups=cleaned, params=req.params.model_dump(),
        engines=engine_hashes, r_lib_hash=r_lib_hash,
    )
    cached = cache.get_json(key, "ushape")
    if cached is not None:
        cached["_cache"] = "hit"
        cached["_cache_key"] = key
        return JSONResponse(cached)

    # call region_popstats over the expanded region
    binary = engines.path("region_popstats")
    if not binary.exists():
        raise HTTPException(500, f"engine not compiled: {binary}")
    sample_list = Path(samples.path)
    try:
        beagle = resolve_beagle(beagle_dir, req.chrom)
    except FileNotFoundError as e:
        raise HTTPException(404, str(e))

    with tempfile.TemporaryDirectory(prefix="ushape_pop_",
                                     dir=str(cache.root)) as scratch:
        scratch_path = Path(scratch)
        try:
            outcome = await run_region_popstats(
                binary=binary, beagle=beagle, sample_list=sample_list,
                chrom=req.chrom,
                fixed_win=(req.params.win_bp, req.params.step_bp),
                win_type=2,
                downsample=1,
                ncores=4,
                groups=cleaned,
                region=(region_lo, region_hi),
                scratch=scratch_path,
            )
        except RuntimeError as e:
            raise HTTPException(500, f"region_popstats failed: {e}")

    popstats_payload = {
        "kind": "popstats_groupwise.v1",
        "chrom": req.chrom,
        "columns": outcome.columns,
        "windows": outcome.rows,
    }

    r_payload = {
        "candidate": {
            "candidate_id": req.candidate_id,
            "chrom": req.chrom,
            "start_bp": req.start_bp, "end_bp": req.end_bp,
        },
        "params": req.params.model_dump(),
        "groups": {g: cleaned.get(g, []) for g in ("HOMO_1", "HET", "HOMO_2")},
        "popstats_payload": popstats_payload,
    }
    try:
        block = await _run_r_driver(ushape_r_dir=ushape_r_dir, payload=r_payload)
    except RuntimeError as e:
        raise HTTPException(500, str(e))
    if isinstance(block, dict) and "error" in block:
        raise HTTPException(500, f"ushape compute returned: {block['error']}")

    payload = {
        "kind": "ushape_evolution_candidate.v1",
        "format_version": "ushape_evolution_v1",
        "candidate": block,
        "engine_hash": engine_hashes.get("region_popstats"),
        "ushape_r_lib_hash": r_lib_hash,
    }
    cache.put_json(key, "ushape", payload)
    payload["_cache"] = "miss"
    payload["_cache_key"] = key
    return JSONResponse(payload)
