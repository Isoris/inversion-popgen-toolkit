#!/usr/bin/env python3
"""
test_region_handler.py
----------------------
Tests for dosage_viewer/region_handler.py — orchestration of store_reader
+ sampling for the standalone viewer's /api/region endpoint.

Builds a small store and exercises:
  - handle_region for all 6 modes (raw / even / random / variance / hybrid / aggregate)
  - cap_exceeded path (raw mode, n_total > max_sites)
  - region_too_large guard
  - chrom_not_found
  - bad_mode / bad_max_sites / cap_above_hard_cap
  - bad_range
  - empty region returns clean envelope
  - random+variance+hybrid modes properly thread `seed`
  - handle_manifest envelope
  - load_candidate_regions: TSV parsing, missing file -> empty dict
  - handle_candidate: candidate -> region
  - handle_breakpoint: candidate + side + window
  - handle_tile: forces mode=aggregate, resolves whole chrom
  - Renderer-style invariants on the envelope: matrix len == positions len,
    every matrix row n_samples wide, missing as -1 (not null) for non-aggregate.

Run:
  python3 tests/dosage_viewer/test_region_handler.py
"""
from __future__ import annotations

import gzip
import math
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "dosage_viewer"))

import importlib
prep = importlib.import_module("01_prepare_dosage_store")
import region_handler as rh

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


tmp = Path(tempfile.mkdtemp(prefix="dv_handler_"))
try:
    # ---------------------------------------------------------------------------
    section("Build a real store with 30 sites x 4 samples")
    base = tmp / "tsv"
    base.mkdir()
    samples = ["S0", "S1", "S2", "S3"]
    samples_path = tmp / "samples.tsv"
    samples_path.write_text("\n".join(samples) + "\n")

    n_sites = 30
    rows = []
    for i in range(n_sites):
        pos = 1000 * (i + 1)         # 1000, 2000, ..., 30000
        # Variances: pick a few high-variance sites at known indices
        if i in (3, 17, 25):
            doses = [0, 0, 2, 2]    # high variance
        elif i in (5, 22):
            doses = [-1, 0, 1, 2]   # moderate, with NA
        else:
            doses = [1, 1, 1, 1]    # zero variance
        rows.append((pos, doses))
    with gzip.open(base / "C_gar_LG28.sites.tsv.gz", "wt") as f:
        for pos, _ in rows:
            f.write(f"{pos}\tA\tG\tNA\tNA\t\n")
    with gzip.open(base / "C_gar_LG28.dosage.tsv.gz", "wt") as f:
        for _, doses in rows:
            f.write("\t".join(["NA" if d == -1 else str(d) for d in doses]) + "\n")

    out = tmp / "store"
    prep.build_store(base=base, out=out, samples_path=samples_path)
    expect((out / "manifest.json").is_file(), "store built (30 sites x 4 samples)")

    # ---------------------------------------------------------------------------
    section("handle_manifest")
    m = rh.handle_manifest(out)
    expect(m["n_samples"] == 4, "manifest n_samples = 4")
    expect(m["samples"] == samples, "manifest samples canonical")
    expect(any(c["name"] == "C_gar_LG28" for c in m["chroms"]),
           "manifest contains LG28")
    expect(m["limits"]["hard_cap"] == 20000, "default hard_cap = 20000")
    expect(m["candidate_regions_loaded"] is False,
           "candidate_regions_loaded false by default")

    # ---------------------------------------------------------------------------
    section("handle_region — raw mode within cap")
    env = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                           start_bp=1000, end_bp=10_000,
                           max_sites=100, mode="raw")
    expect(env["sampling_mode"] == "raw", "mode = raw")
    expect(env["n_sites_total"] == 10, f"10 sites in [1000, 10000] (got {env['n_sites_total']})")
    expect(env["n_sites_returned"] == 10, "all 10 returned")
    expect(env["downsampled"] is False, "raw within cap: not downsampled")
    expect(len(env["matrix"]) == 10, "matrix has 10 rows")
    expect(all(len(row) == 4 for row in env["matrix"]),
           "every row n_samples=4 wide")
    # Renderer contract: missing must be -1, not null
    has_null = any(v is None for row in env["matrix"] for v in row)
    expect(not has_null, "no null in matrix — missing encoded as -1")

    # ---------------------------------------------------------------------------
    section("handle_region — raw cap exceeded")
    threw = False
    detail = None
    try:
        rh.handle_region(store_root=out, chrom="C_gar_LG28",
                         start_bp=1, end_bp=30_000,
                         max_sites=5, mode="raw")
    except rh.RegionError as e:
        threw = True
        detail = e.detail
    expect(threw, "raw + 30 sites > cap=5 -> RegionError")
    expect(detail["error"] == "raw_cap_exceeded",
           f"error code (got {detail['error']!r})")
    expect(detail["n_sites_total"] == 30, "error includes n_sites_total")

    # ---------------------------------------------------------------------------
    section("handle_region — even mode")
    env = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                           start_bp=1, end_bp=30_000,
                           max_sites=5, mode="even")
    expect(env["sampling_mode"] == "even", "mode = even")
    expect(env["n_sites_total"] == 30, "n_sites_total = 30")
    expect(env["n_sites_returned"] == 5, f"5 sites returned (got {env['n_sites_returned']})")
    expect(env["downsampled"] is True, "even is downsampled when n > cap")
    expect(env["positions"] == sorted(env["positions"]),
           f"even: positions sorted (got {env['positions']})")
    expect(len(env["matrix"]) == 5, "5 matrix rows")

    # ---------------------------------------------------------------------------
    section("handle_region — random mode (seed reproducibility)")
    env_a = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                             start_bp=1, end_bp=30_000,
                             max_sites=10, mode="random", seed=42)
    env_b = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                             start_bp=1, end_bp=30_000,
                             max_sites=10, mode="random", seed=42)
    expect(env_a["positions"] == env_b["positions"],
           "random seed=42 -> same positions on two calls")
    env_c = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                             start_bp=1, end_bp=30_000,
                             max_sites=10, mode="random", seed=999)
    expect(env_a["positions"] != env_c["positions"],
           "different seeds -> different selections")
    expect(env_a["positions"] == sorted(env_a["positions"]),
           "random: positions sorted")
    expect(env_a["_meta"]["seed"] == 42, "_meta.seed echoed")

    # ---------------------------------------------------------------------------
    section("handle_region — variance mode picks high-variance sites")
    env = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                           start_bp=1, end_bp=30_000,
                           max_sites=3, mode="variance")
    # Sites at indices 3, 17, 25 have variance > 0; positions are 4000, 18000, 26000
    expected = {4000, 18000, 26000}
    expect(set(env["positions"]) == expected,
           f"variance: top-3 = {expected} (got {set(env['positions'])})")
    expect(env["positions"] == sorted(env["positions"]),
           "variance: sorted by pos")

    # ---------------------------------------------------------------------------
    section("handle_region — hybrid mode")
    env = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                           start_bp=1, end_bp=30_000,
                           max_sites=10, mode="hybrid", seed=1)
    expect(env["sampling_mode"] == "hybrid", "mode = hybrid")
    expect(env["downsampled"] is True, "hybrid downsampled when n > cap")
    expect(env["positions"] == sorted(env["positions"]),
           "hybrid: sorted by pos")
    # The top-3 variance sites (4000, 18000, 26000) should appear in hybrid
    # since the 30% variance share = 3 picks.
    sel_pos = set(env["positions"])
    var_top = {4000, 18000, 26000}
    expect(len(sel_pos & var_top) >= 2,
           f"hybrid: at least 2 variance-top sites included (got {sel_pos & var_top})")

    # ---------------------------------------------------------------------------
    section("handle_region — aggregate mode")
    env = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                           start_bp=1, end_bp=30_000,
                           max_sites=5, mode="aggregate")
    expect(env["sampling_mode"] == "aggregate", "mode = aggregate")
    expect(env["n_sites_returned"] == 5, "5 bins (n_sites_returned)")
    expect(len(env["positions"]) == 5, "5 bin centers")
    expect(env["downsampled"] is True, "aggregate is always downsampled")
    expect(len(env["matrix"]) == 5, "5 matrix rows (one per bin)")
    expect(all(len(row) == 4 for row in env["matrix"]),
           "every bin row n_samples=4 wide")
    # Aggregate matrix is float (or None for empty bins) — different
    # from per-site -1 convention.
    sample_vals = env["matrix"][0]
    expect(all(v is None or isinstance(v, float) for v in sample_vals),
           "aggregate matrix cells are float or None (NaN -> None)")
    # site_ids are nulls
    expect(all(s is None for s in env["site_ids"]),
           "aggregate: site_ids null")
    expect(env["_meta"].get("is_aggregate") is True,
           "_meta.is_aggregate flag set")

    # ---------------------------------------------------------------------------
    section("handle_region — empty region")
    env = rh.handle_region(store_root=out, chrom="C_gar_LG28",
                           start_bp=900_000, end_bp=950_000,
                           max_sites=100, mode="hybrid")
    expect(env["n_sites_total"] == 0 and env["n_sites_returned"] == 0,
           "empty region: zeros throughout")
    expect(env["matrix"] == [], "empty region: empty matrix")
    expect("warning" in env and env["warning"], "empty region: warning set")

    # ---------------------------------------------------------------------------
    section("handle_region — error paths")
    threw = False
    try:
        rh.handle_region(store_root=out, chrom="C_gar_LG_NEVER",
                         start_bp=1, end_bp=100, max_sites=10, mode="raw")
    except rh.RegionError as e:
        threw = e.code == "chrom_not_found"
    expect(threw, "missing chrom -> chrom_not_found")

    threw = False
    try:
        rh.handle_region(store_root=out, chrom="C_gar_LG28",
                         start_bp=1, end_bp=100_000_000,
                         max_sites=10, mode="raw")
    except rh.RegionError as e:
        threw = e.code == "region_too_large"
    expect(threw, "region > 50 Mb -> region_too_large")

    threw = False
    try:
        rh.handle_region(store_root=out, chrom="C_gar_LG28",
                         start_bp=1, end_bp=1000,
                         max_sites=10, mode="banana")
    except rh.RegionError as e:
        threw = e.code == "bad_mode"
    expect(threw, "unknown mode -> bad_mode")

    threw = False
    try:
        rh.handle_region(store_root=out, chrom="C_gar_LG28",
                         start_bp=1000, end_bp=500,
                         max_sites=10, mode="raw")
    except rh.RegionError as e:
        threw = e.code == "bad_range"
    expect(threw, "end < start -> bad_range")

    threw = False
    try:
        rh.handle_region(store_root=out, chrom="C_gar_LG28",
                         start_bp=1, end_bp=100,
                         max_sites=999_999, mode="raw")
    except rh.RegionError as e:
        threw = e.code == "cap_above_hard_cap"
    expect(threw, "max_sites > hard_cap -> cap_above_hard_cap")

    # ---------------------------------------------------------------------------
    section("load_candidate_regions")
    cand_path = tmp / "candidates.tsv"
    cand_path.write_text(
        "candidate_id\tchrom\tstart\tend\tleft_breakpoint\tright_breakpoint\tnotes\n"
        "INV_LG28_001\tC_gar_LG28\t5000\t25000\t10000\t20000\ttest candidate\n"
        "INV_LG28_002\tC_gar_LG28\t8000\t12000\t\t\t\n"
    )
    cands = rh.load_candidate_regions(cand_path)
    expect(len(cands) == 2, f"2 candidates loaded (got {len(cands)})")
    c1 = cands["INV_LG28_001"]
    expect(c1["chrom"] == "C_gar_LG28", "candidate chrom")
    expect(c1["start"] == 5000 and c1["end"] == 25000, "candidate range")
    expect(c1.get("left_breakpoint") == 10000, "left breakpoint loaded")
    c2 = cands["INV_LG28_002"]
    expect("left_breakpoint" not in c2,
           "missing breakpoints field absent (graceful)")

    # Missing file -> empty dict (graceful)
    expect(rh.load_candidate_regions(tmp / "no_such_cands.tsv") == {},
           "missing TSV -> {}")

    # ---------------------------------------------------------------------------
    section("handle_candidate")
    env = rh.handle_candidate(store_root=out, candidates=cands,
                              candidate_id="INV_LG28_001",
                              max_sites=100, mode="raw")
    expect(env["chrom"] == "C_gar_LG28", "candidate chrom resolved")
    expect(env["start_bp"] == 5000 and env["end_bp"] == 25000,
           "candidate window applied")

    threw = False
    try:
        rh.handle_candidate(store_root=out, candidates=cands,
                            candidate_id="not_a_real_cand",
                            mode="raw")
    except rh.RegionError as e:
        threw = e.code == "candidate_not_found"
    expect(threw, "unknown candidate -> candidate_not_found")

    # ---------------------------------------------------------------------------
    section("handle_breakpoint")
    env = rh.handle_breakpoint(store_root=out, candidates=cands,
                                candidate_id="INV_LG28_001", side="left",
                                window_bp=3000, max_sites=100, mode="even")
    # Left breakpoint = 10000; window ±3000 -> [7000, 13000]
    expect(env["start_bp"] == 7000, f"breakpoint left start = 7000 (got {env['start_bp']})")
    expect(env["end_bp"] == 13000, f"breakpoint left end = 13000 (got {env['end_bp']})")
    expect(env["_meta"].get("breakpoint_bp") == 10000,
           "_meta.breakpoint_bp echoed")
    expect(env["_meta"].get("side") == "left", "_meta.side echoed")

    # Side fallback when breakpoint absent: uses start (left) / end (right)
    env = rh.handle_breakpoint(store_root=out, candidates=cands,
                                candidate_id="INV_LG28_002", side="right",
                                window_bp=2000, max_sites=100, mode="even")
    # No right_breakpoint in INV_LG28_002, falls back to end=12000; window ±2000 = [10000, 14000]
    expect(env["start_bp"] == 10000 and env["end_bp"] == 14000,
           f"breakpoint fallback uses end (got {env['start_bp']}-{env['end_bp']})")

    threw = False
    try:
        rh.handle_breakpoint(store_root=out, candidates=cands,
                              candidate_id="INV_LG28_001", side="middle",
                              window_bp=1000, mode="even")
    except rh.RegionError as e:
        threw = e.code == "bad_side"
    expect(threw, "bad side -> bad_side")

    # ---------------------------------------------------------------------------
    section("handle_tile")
    env = rh.handle_tile(store_root=out, chrom="C_gar_LG28", width=10)
    expect(env["sampling_mode"] == "aggregate",
           "tile forces mode=aggregate")
    expect(env["n_sites_returned"] == 10, "tile width = 10 bins")
    expect(env["start_bp"] == 1, "tile starts at 1")
    expect(env["end_bp"] >= 30_000, "tile spans whole chrom")

    threw = False
    try:
        rh.handle_tile(store_root=out, chrom="C_gar_LG28", width=999_999)
    except rh.RegionError as e:
        threw = e.code == "bad_tile_width"
    expect(threw, "tile width above hard cap -> bad_tile_width")

    threw = False
    try:
        rh.handle_tile(store_root=out, chrom="C_gar_LG_NEVER", width=10)
    except rh.RegionError as e:
        threw = e.code == "chrom_not_found"
    expect(threw, "tile missing chrom -> chrom_not_found")

finally:
    if tmp.exists():
        shutil.rmtree(tmp)

# ---------------------------------------------------------------------------
print(f"\n=========================================")
print(f"  PASS: {passed}   FAIL: {failed}")
print(f"=========================================")
sys.exit(0 if failed == 0 else 1)
