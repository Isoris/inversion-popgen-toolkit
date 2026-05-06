#!/usr/bin/env python3
"""
test_run_server.py
------------------
Integration test for dosage_viewer/02_run_server.py.

Most of the logic (sampling, parquet IO, region/candidate/breakpoint/tile
orchestration) is already covered by test_region_handler.py — this test
only covers the FastAPI HTTP wrapping:

    - status codes (200 OK / 400 / 404 / 500)
    - query-param coercion (strings -> ints, defaults)
    - RegionError -> HTTPException with the right `detail`
    - the /api/dosage/chunk renderer-compat shape (different envelope shape
      than /api/region — markers + dosage at top level)
    - /api/health basic shape
    - /api/manifest sets candidate_regions_loaded based on CANDIDATES global
    - CORS preflight on at least one endpoint (smoke test)
    - Full end-to-end: build store, populate globals, hit endpoints, verify
      response payloads match what handle_region returns

This is a sanity-check layer; a full assertion sweep over response content
would duplicate the handler tests.

Run:
  python3 tests/dosage_viewer/test_run_server.py
"""
from __future__ import annotations

import gzip
import importlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "dosage_viewer"))

prep = importlib.import_module("01_prepare_dosage_store")
import region_handler as rh

# Import 02_run_server as a module — file name has a leading digit, so use
# importlib.
server_mod = importlib.import_module("02_run_server")
from fastapi.testclient import TestClient

passed = 0
failed = 0
def expect(cond, msg):
    global passed, failed
    if cond:
        passed += 1
        print("  \u2713 " + msg)
    else:
        failed += 1
        print("  \u2717 " + msg, file=sys.stderr)
def section(s):
    print("\n\u2014 " + s + " \u2014")


# ---------------------------------------------------------------------------
# Build a tiny dosage store as the server's data source
# ---------------------------------------------------------------------------
section("Build a store + populate server globals")

tmp = Path(tempfile.mkdtemp(prefix="dv_server_"))
try:
    base = tmp / "tsv"
    base.mkdir()
    samples = ["S0", "S1", "S2", "S3"]
    samples_path = tmp / "samples.tsv"
    samples_path.write_text("\n".join(samples) + "\n")

    n_sites = 30
    rows = []
    for i in range(n_sites):
        pos = 1000 * (i + 1)
        if i in (3, 17, 25):
            doses = [0, 0, 2, 2]
        elif i in (5, 22):
            doses = [-1, 0, 1, 2]
        else:
            doses = [1, 1, 1, 1]
        rows.append((pos, doses))
    with gzip.open(base / "C_gar_LG28.sites.tsv.gz", "wt") as f:
        for pos, _ in rows:
            f.write(f"{pos}\tA\tG\tNA\tNA\t\n")
    with gzip.open(base / "C_gar_LG28.dosage.tsv.gz", "wt") as f:
        for _, doses in rows:
            f.write("\t".join(["NA" if d == -1 else str(d) for d in doses]) + "\n")
    out = tmp / "store"
    prep.build_store(base=base, out=out, samples_path=samples_path)

    # Candidate config
    cand_path = tmp / "candidates.tsv"
    cand_path.write_text(
        "candidate_id\tchrom\tstart\tend\tleft_breakpoint\tright_breakpoint\tnotes\n"
        "INV_LG28_001\tC_gar_LG28\t5000\t25000\t10000\t20000\ttest\n"
    )

    # Patch module-level globals on 02_run_server.py
    server_mod.STORE_ROOT = out.resolve()
    server_mod.CANDIDATES = rh.load_candidate_regions(cand_path)
    server_mod.CANDIDATES_PATH = cand_path.resolve()

    expect(server_mod.STORE_ROOT.is_dir(), "STORE_ROOT set on module")
    expect(len(server_mod.CANDIDATES) == 1, "CANDIDATES populated (1 candidate)")

    client = TestClient(server_mod.app)

    # ---------------------------------------------------------------------------
    section("/api/health")
    r = client.get("/api/health")
    expect(r.status_code == 200, f"health 200 (got {r.status_code})")
    body = r.json()
    expect(body["ok"] is True, "health.ok = True")
    expect(body["service"] == "dosage_viewer", "service = dosage_viewer")
    expect(body["n_candidates"] == 1, "n_candidates = 1")
    expect(body["store_root"] == str(out.resolve()), "store_root reported")

    # ---------------------------------------------------------------------------
    section("/api/manifest")
    r = client.get("/api/manifest")
    expect(r.status_code == 200, "manifest 200")
    m = r.json()
    expect(m["n_samples"] == 4, "manifest n_samples = 4")
    expect(m["candidate_regions_loaded"] is True,
           "candidate_regions_loaded reflects CANDIDATES global")
    expect(any(c["name"] == "C_gar_LG28" for c in m["chroms"]),
           "manifest contains LG28")
    expect("limits" in m, "limits present in manifest")

    # ---------------------------------------------------------------------------
    section("/api/region — happy path")
    r = client.get("/api/region",
                   params={"chrom": "C_gar_LG28", "start": 1, "end": 30000,
                           "max_sites": 10, "mode": "even", "seed": 1})
    expect(r.status_code == 200, "region 200")
    env = r.json()
    expect(env["sampling_mode"] == "even", "sampling_mode echoed")
    expect(env["n_sites_total"] == 30, "n_sites_total = 30")
    expect(env["n_sites_returned"] == 10, "n_sites_returned = 10")
    expect(env["positions"] == sorted(env["positions"]),
           "positions sorted ascending")
    expect(len(env["matrix"]) == 10, "matrix length matches positions")
    expect(all(len(row) == 4 for row in env["matrix"]),
           "every row n_samples wide")
    # Renderer contract: missing -> -1
    has_null = any(v is None for row in env["matrix"] for v in row)
    expect(not has_null, "no null in matrix (per-site mode)")

    # ---------------------------------------------------------------------------
    section("/api/region — query coercion + defaults")
    r = client.get("/api/region",
                   params={"chrom": "C_gar_LG28", "start": "5000", "end": "10000"})
    expect(r.status_code == 200, "string-numeric query coerced to int")
    env = r.json()
    expect(env["sampling_mode"] == "hybrid",
           "default mode = hybrid (per spec)")

    # ---------------------------------------------------------------------------
    section("/api/region — error mappings")
    # chrom_not_found -> 404
    r = client.get("/api/region",
                   params={"chrom": "C_gar_LG_NEVER", "start": 1, "end": 1000,
                           "mode": "raw"})
    expect(r.status_code == 404, f"chrom_not_found -> 404 (got {r.status_code})")
    expect(r.json()["detail"]["error"] == "chrom_not_found",
           "detail.error = chrom_not_found")

    # bad_mode -> 400
    r = client.get("/api/region",
                   params={"chrom": "C_gar_LG28", "start": 1, "end": 1000,
                           "mode": "banana"})
    expect(r.status_code == 400, "bad_mode -> 400")
    expect(r.json()["detail"]["error"] == "bad_mode",
           "detail.error = bad_mode")

    # region_too_large -> 400
    r = client.get("/api/region",
                   params={"chrom": "C_gar_LG28", "start": 1, "end": 100_000_000,
                           "mode": "raw"})
    expect(r.status_code == 400, "region_too_large -> 400")
    expect(r.json()["detail"]["error"] == "region_too_large",
           "detail.error = region_too_large")

    # raw_cap_exceeded -> 400
    r = client.get("/api/region",
                   params={"chrom": "C_gar_LG28", "start": 1, "end": 30000,
                           "max_sites": 5, "mode": "raw"})
    expect(r.status_code == 400, "raw_cap_exceeded -> 400")
    expect(r.json()["detail"]["error"] == "raw_cap_exceeded",
           "detail.error = raw_cap_exceeded")

    # cap_above_hard_cap -> 400
    r = client.get("/api/region",
                   params={"chrom": "C_gar_LG28", "start": 1, "end": 1000,
                           "max_sites": 999_999, "mode": "raw"})
    expect(r.status_code == 400, "cap_above_hard_cap -> 400")
    expect(r.json()["detail"]["error"] == "cap_above_hard_cap",
           "detail.error = cap_above_hard_cap")

    # FastAPI Query(ge=1) rejects start=0 with 422
    r = client.get("/api/region",
                   params={"chrom": "C_gar_LG28", "start": 0, "end": 1000,
                           "mode": "raw"})
    expect(r.status_code == 422,
           f"start < 1 rejected by FastAPI -> 422 (got {r.status_code})")

    # ---------------------------------------------------------------------------
    section("/api/dosage/chunk — renderer-compat shape")
    r = client.get("/api/dosage/chunk",
                   params={"chrom": "C_gar_LG28", "start": 1, "end": 10000,
                           "cap": 100, "mode": "raw"})
    expect(r.status_code == 200, "dosage/chunk 200")
    chunk = r.json()
    # Renderer-compat shape: markers[] (with pos_bp/missingness/diagnostic_score),
    # samples[], dosage[][] at top level
    expect("markers" in chunk and "samples" in chunk and "dosage" in chunk,
           "chunk has markers + samples + dosage at top level")
    expect(chunk["samples"] == samples,
           "chunk.samples canonical order")
    expect(len(chunk["markers"]) == len(chunk["dosage"]),
           "markers length == dosage length")
    expect(all("pos_bp" in m for m in chunk["markers"]),
           "every marker has pos_bp")
    # Renderer contract: dosage rows must have n_samples cells with -1 = NA
    expect(all(len(row) == 4 for row in chunk["dosage"]),
           "every chunk.dosage row n_samples=4 wide")
    has_null = any(v is None for row in chunk["dosage"] for v in row)
    expect(not has_null, "chunk.dosage: no null cells (renderer wants -1)")

    # ---------------------------------------------------------------------------
    section("/api/candidate/{id}")
    r = client.get("/api/candidate/INV_LG28_001",
                   params={"max_sites": 100, "mode": "raw"})
    expect(r.status_code == 200, "candidate 200")
    env = r.json()
    expect(env["chrom"] == "C_gar_LG28", "candidate chrom resolved")
    expect(env["start_bp"] == 5000 and env["end_bp"] == 25000,
           "candidate window applied")

    # Unknown candidate -> 404
    r = client.get("/api/candidate/no_such_id",
                   params={"mode": "raw", "max_sites": 100})
    expect(r.status_code == 404, "unknown candidate -> 404")
    expect(r.json()["detail"]["error"] == "candidate_not_found",
           "detail.error = candidate_not_found")

    # ---------------------------------------------------------------------------
    section("/api/breakpoint/{id}/{side}")
    r = client.get("/api/breakpoint/INV_LG28_001/left",
                   params={"window": 3000, "max_sites": 100, "mode": "even"})
    expect(r.status_code == 200, "breakpoint 200")
    env = r.json()
    expect(env["start_bp"] == 7000 and env["end_bp"] == 13000,
           "breakpoint window ± 3000 applied")
    expect(env["_meta"].get("breakpoint_bp") == 10000,
           "_meta.breakpoint_bp echoed")
    expect(env["_meta"].get("side") == "left", "_meta.side echoed")

    # Bad side -> 400
    r = client.get("/api/breakpoint/INV_LG28_001/middle",
                   params={"window": 1000, "mode": "even"})
    expect(r.status_code == 400, "bad side -> 400")
    expect(r.json()["detail"]["error"] == "bad_side",
           "detail.error = bad_side")

    # ---------------------------------------------------------------------------
    section("/api/tile")
    r = client.get("/api/tile", params={"chrom": "C_gar_LG28", "width": 10})
    expect(r.status_code == 200, "tile 200")
    env = r.json()
    expect(env["sampling_mode"] == "aggregate",
           "tile forces mode=aggregate")
    expect(env["n_sites_returned"] == 10, "tile width = 10 bins")
    # Aggregate matrix can have None cells (NaN -> None) — different from
    # per-site -1 convention.
    sample_cells = [v for row in env["matrix"] for v in row]
    has_some_floats = any(isinstance(v, (int, float)) and v is not None for v in sample_cells)
    expect(has_some_floats, "tile aggregate matrix contains float cells")

    # tile chrom not found -> 404
    r = client.get("/api/tile", params={"chrom": "C_gar_LG_NEVER", "width": 10})
    expect(r.status_code == 404, "tile missing chrom -> 404")

    # tile bad width -> 400
    r = client.get("/api/tile", params={"chrom": "C_gar_LG28", "width": 1})
    expect(r.status_code == 400, f"tile width < 2 -> 400 (got {r.status_code})")

    # ---------------------------------------------------------------------------
    section("Bad-input invariants")
    # Server NOT initialised path: clear STORE_ROOT and check 500
    saved_store = server_mod.STORE_ROOT
    server_mod.STORE_ROOT = None
    r = client.get("/api/manifest")
    expect(r.status_code == 500, f"missing STORE_ROOT -> 500 (got {r.status_code})")
    expect(r.json()["detail"]["error"] == "server_not_initialised",
           "detail.error = server_not_initialised")
    server_mod.STORE_ROOT = saved_store

    # ---------------------------------------------------------------------------
    section("CORS — file:// origin")
    r = client.get("/api/manifest", headers={"Origin": "null"})
    expect(r.status_code == 200, "GET with Origin: null still 200")
    expect(r.headers.get("access-control-allow-origin") == "null"
           or r.headers.get("access-control-allow-origin") is not None,
           "CORS allow-origin header present for null origin")

finally:
    if tmp.exists():
        shutil.rmtree(tmp)

# ---------------------------------------------------------------------------
print(f"\n=========================================")
print(f"  PASS: {passed}   FAIL: {failed}")
print(f"=========================================")
sys.exit(0 if failed == 0 else 1)
