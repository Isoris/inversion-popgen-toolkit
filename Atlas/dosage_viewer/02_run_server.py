#!/usr/bin/env python3
"""
02_run_server.py
================
Standalone FastAPI server for the dosage viewer.

Endpoints (all GET):
    /api/manifest                            store metadata + chrom list
    /api/region                              ?chrom=&start=&end=&max_sites=&mode=&seed=
    /api/dosage/chunk                        renderer-compatible bridge shape
                                             (alias of /api/region with mode=raw,
                                              suitable for atlas-side bridge use too)
    /api/candidate/{candidate_id}            wraps /api/region using candidate_regions.tsv
    /api/breakpoint/{candidate_id}/{side}    ± window around a breakpoint
    /api/tile                                ?chrom=&width=  (whole-chrom overview)
    /                                        static index.html (frontend)
    /static/*                                static assets

Bind: 127.0.0.1:8767 by default (different port from popstats's 8765 so they
can run side-by-side; spec §Decision points said 8767 was proposed).

Usage:
    python 02_run_server.py \\
        --store dosage_viewer/output/dosage_store \\
        --candidates dosage_viewer/config/candidate_regions.tsv \\
        --static dosage_viewer/static \\
        --port 8767

Author:  Quentin Andres
Project: MS_Inversions_North_african_catfish
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles

# Local sibling imports
sys.path.insert(0, str(Path(__file__).parent))
from region_handler import (
    RegionError,
    handle_region, handle_manifest,
    handle_candidate, handle_breakpoint, handle_tile,
    load_candidate_regions,
    DEFAULT_MAX_SITES, DEFAULT_HARD_CAP, DEFAULT_MAX_REGION_BP,
    ALL_MODES,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger("dosage_viewer")


# Module-level config — populated by main(). Keeping these globals (rather
# than threading via dependency injection) matches popstats_server.py's idiom.
STORE_ROOT: Optional[Path] = None
CANDIDATES: dict = {}
CANDIDATES_PATH: Optional[Path] = None
STATIC_DIR: Optional[Path] = None


def _ensure_store():
    if STORE_ROOT is None:
        raise HTTPException(status_code=500, detail={
            "error": "server_not_initialised",
            "message": "STORE_ROOT not set; run via __main__ entry point",
        })


def _wrap(coro_or_call):
    """Run a handler call; convert RegionError -> HTTPException."""
    try:
        return coro_or_call()
    except RegionError as e:
        # Map error codes to HTTP status. Most are 400 (client misuse);
        # chrom_not_found / candidate_not_found are 404; server errors none.
        status = 400
        if e.code in ("chrom_not_found", "candidate_not_found"):
            status = 404
        elif e.code in ("server_not_initialised",):
            status = 500
        raise HTTPException(status_code=status, detail=e.detail)
    except FileNotFoundError as e:
        raise HTTPException(status_code=500, detail={
            "error": "store_file_missing", "message": str(e),
        })


# =============================================================================
# App
# =============================================================================

app = FastAPI(title="dosage_viewer", version="1.0",
              description="Standalone dosage heatmap viewer (S6 spec).")

# CORS — local-only by default. Allow file:// (origin "null") and 127.0.0.1
# so a separately-served atlas can hit the viewer's /api endpoints.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:8767", "http://127.0.0.1:8767",
                   "http://localhost:8765", "http://127.0.0.1:8765",
                   "null"],
    allow_credentials=False,
    allow_methods=["GET"],
    allow_headers=["*"],
)


@app.get("/api/manifest")
async def get_manifest():
    _ensure_store()
    payload = _wrap(lambda: handle_manifest(STORE_ROOT))
    payload["candidate_regions_loaded"] = bool(CANDIDATES)
    if CANDIDATES_PATH is not None:
        payload["candidate_regions_path"] = str(CANDIDATES_PATH)
    return JSONResponse(payload)


@app.get("/api/region")
async def get_region(
    chrom: str,
    start: int = Query(..., ge=1),
    end: int = Query(..., ge=1),
    max_sites: int = DEFAULT_MAX_SITES,
    mode: str = "hybrid",
    seed: int = 1,
):
    _ensure_store()
    payload = _wrap(lambda: handle_region(
        store_root=STORE_ROOT, chrom=chrom,
        start_bp=start, end_bp=end,
        max_sites=max_sites, mode=mode, seed=seed,
    ))
    return JSONResponse(payload)


# Renderer-compatible shape used by the atlas-side bridge installer. We
# accept this here too so a single server (this one) can host both the
# bridge and the standalone viewer if the operator prefers to run only
# one process. Defaults to mode=raw for renderer compat (selectTopMarkers
# does its own selection on the atlas side).
@app.get("/api/dosage/chunk")
async def get_dosage_chunk(
    chrom: str,
    start: int = Query(..., ge=1),
    end: int = Query(..., ge=1),
    cap: int = DEFAULT_MAX_SITES,
    mode: str = "raw",
):
    """Return the renderer-compatible chunk shape (subset of /api/region's
    rich envelope, with the matrix at the top level + samples + markers
    structured the way the atlas's renderDosageHeatmap expects).
    """
    _ensure_store()
    rich = _wrap(lambda: handle_region(
        store_root=STORE_ROOT, chrom=chrom,
        start_bp=start, end_bp=end,
        max_sites=cap, mode=mode, seed=1,
    ))
    # Translate to the renderer-compatible shape (see specs_todo/from_turn129
    # /S6_dosage_heatmap_streaming_viewer.md §GET /api/dosage/chunk).
    chunk = {
        "chrom": rich["chrom"],
        "start_bp": rich["start_bp"],
        "end_bp": rich["end_bp"],
        "samples": rich["samples"],
        "markers": [
            {"pos_bp": p,
             "missingness": rich["missingness"][i] if i < len(rich["missingness"]) else None,
             "diagnostic_score": rich["diagnostic_score"][i] if i < len(rich["diagnostic_score"]) else None}
            for i, p in enumerate(rich["positions"])
        ],
        "dosage": rich["matrix"],
        "_meta": rich["_meta"],
    }
    return JSONResponse(chunk)


@app.get("/api/candidate/{candidate_id}")
async def get_candidate(candidate_id: str,
                         max_sites: int = DEFAULT_MAX_SITES,
                         mode: str = "hybrid", seed: int = 1):
    _ensure_store()
    payload = _wrap(lambda: handle_candidate(
        store_root=STORE_ROOT, candidates=CANDIDATES,
        candidate_id=candidate_id,
        max_sites=max_sites, mode=mode, seed=seed,
    ))
    return JSONResponse(payload)


@app.get("/api/breakpoint/{candidate_id}/{side}")
async def get_breakpoint(candidate_id: str, side: str,
                          window: int = 500_000,
                          max_sites: int = DEFAULT_MAX_SITES,
                          mode: str = "even", seed: int = 1):
    _ensure_store()
    payload = _wrap(lambda: handle_breakpoint(
        store_root=STORE_ROOT, candidates=CANDIDATES,
        candidate_id=candidate_id, side=side,
        window_bp=window, max_sites=max_sites,
        mode=mode, seed=seed,
    ))
    return JSONResponse(payload)


@app.get("/api/tile")
async def get_tile(chrom: str, width: int = 500):
    _ensure_store()
    payload = _wrap(lambda: handle_tile(
        store_root=STORE_ROOT, chrom=chrom, width=width,
    ))
    return JSONResponse(payload)


@app.get("/api/health")
async def get_health():
    return {
        "ok": True,
        "service": "dosage_viewer",
        "version": "1.0",
        "store_root": str(STORE_ROOT) if STORE_ROOT else None,
        "n_candidates": len(CANDIDATES),
        "candidates_path": str(CANDIDATES_PATH) if CANDIDATES_PATH else None,
    }


# =============================================================================
# Static frontend (mounted by main() if --static is provided)
# =============================================================================

def _mount_static(app: FastAPI, static_dir: Path) -> None:
    """Serve <static_dir>/index.html at /, plus <static_dir>/* under /static."""
    if not static_dir.is_dir():
        log.warning("static dir not found: %s — frontend disabled", static_dir)
        return
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")
    index_path = static_dir / "index.html"

    @app.get("/")
    async def get_index():
        if index_path.is_file():
            return FileResponse(str(index_path))
        return JSONResponse(
            {"error": "no_index", "expected_path": str(index_path)},
            status_code=404,
        )


# =============================================================================
# Entry point
# =============================================================================

def main():
    global STORE_ROOT, CANDIDATES, CANDIDATES_PATH, STATIC_DIR

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", type=Path, required=True,
                    help="dosage_store/ directory built by 01_prepare_dosage_store.py")
    ap.add_argument("--candidates", type=Path, default=None,
                    help="optional: candidate_regions.tsv (used by /api/candidate "
                         "and /api/breakpoint)")
    ap.add_argument("--static", type=Path, default=None,
                    help="optional: static frontend directory (containing index.html)")
    ap.add_argument("--host", default="127.0.0.1",
                    help="bind host (default 127.0.0.1; do not bind 0.0.0.0 "
                         "without auth)")
    ap.add_argument("--port", type=int, default=8767,
                    help="bind port (default 8767)")
    ap.add_argument("--reload", action="store_true",
                    help="dev-mode auto-reload")
    args = ap.parse_args()

    STORE_ROOT = args.store.resolve()
    if not STORE_ROOT.is_dir():
        raise SystemExit(f"store dir not found: {STORE_ROOT}")
    if not (STORE_ROOT / "manifest.json").is_file():
        raise SystemExit(f"store manifest missing: {STORE_ROOT}/manifest.json. "
                         "Run 01_prepare_dosage_store.py first.")

    if args.candidates:
        CANDIDATES_PATH = args.candidates.resolve()
        CANDIDATES = load_candidate_regions(CANDIDATES_PATH)
        log.info("loaded %d candidates from %s", len(CANDIDATES), CANDIDATES_PATH)

    if args.static:
        STATIC_DIR = args.static.resolve()
        _mount_static(app, STATIC_DIR)
        log.info("static frontend at /, served from %s", STATIC_DIR)

    log.info("dosage_viewer listening on http://%s:%d (store=%s)",
             args.host, args.port, STORE_ROOT)

    import uvicorn
    uvicorn.run(app, host=args.host, port=args.port, log_level="info",
                reload=args.reload)


if __name__ == "__main__":
    main()
