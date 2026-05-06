#!/usr/bin/env python3
"""
test_dosage_bridge.py
---------------------
End-to-end test for server_turn1/dosage_bridge.py — the atlas-side bridge
endpoint that hands the existing renderDosageHeatmap a chunk in the shape
it already expects.

What this exercises:
  - Synth a tiny <chrom>.sites.tsv.gz + <chrom>.dosage.tsv.gz + samples.tsv
    that match the §Source data layout in S6_dosage_heatmap_streaming_viewer.md.
  - Validate request-side (chrom whitelist, end >= start, cap caps).
  - Window filter (start/end) returns only sites in [start, end].
  - Dosage rows align row-for-row with markers.
  - Missing values are -1 (NOT null) on the wire — the renderer's v < 0
    test depends on this.
  - mode=raw + cap-exceeded -> raw_cap_exceeded (HTTP 400).
  - mode=even returns positionally-uniform sites, monotonic by pos.
  - Header detection on dosage TSV (with + without header row).
  - chrom_not_found -> HTTP 404.
  - chrom safety: path traversal etc. rejected at request validation.
  - region_too_large guard.

Run:
  python3 tests/dosage_bridge/test_dosage_bridge.py
"""
from __future__ import annotations

import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "server_turn1"))

from dosage_bridge import (
    DosageChunkReq, handle_dosage_chunk,
    DEFAULT_HARD_CAP, DEFAULT_MAX_REGION_BP,
    _maybe_dosage_int, _maybe_float, _select_even, _closest_index,
)
from fastapi import HTTPException

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
# Synth a tiny dosage store
# ---------------------------------------------------------------------------
section("Synthesizing a tiny store")

tmp = Path(tempfile.mkdtemp(prefix="dosage_bridge_"))
try:
    dose_dir = tmp / "04_dosage_by_chr"
    dose_dir.mkdir(parents=True)

    # 6 samples, 12 sites on a synthetic LG28 ranging 1_000_000 .. 1_011_000.
    samples = [f"CGA_{i:03d}" for i in range(6)]
    samples_path = tmp / "samples.tsv"
    samples_path.write_text("\n".join(samples) + "\n")

    # Site positions and content.
    # site_idx | pos       | dosages (per sample)
    #     0      1000000     0 0 0 1 1 2
    #     1      1001000     0 0 1 1 2 2
    #     2      1002000     0 1 1 2 2 2     (NA in col 0 to test missing)
    #     ...
    site_table = [
        (1_000_000, [ 0, 0, 0,  1, 1, 2], 0.00, None),
        (1_001_000, [ 0, 0, 1,  1, 2, 2], 0.05, 0.42),
        (1_002_000, [-1, 1, 1,  2, 2, 2], 0.17, None),  # col 0 NA
        (1_003_000, [ 0, 0, 1,  2, 2, 2], 0.03, 0.91),
        (1_004_000, [ 0, 1, 1,  2, 2, 2], 0.00, None),
        (1_005_000, [ 0, 0, 0,  0, 1, 2], 0.00, None),
        (1_006_000, [ 1, 1, 1,  1, 2, 2], 0.10, 0.55),
        (1_007_000, [ 0, 1, 1,  2, 2, 2], 0.04, None),
        (1_008_000, [ 0, 0, 0,  1, 2, 2], 0.00, 0.21),
        (1_009_000, [ 1, 1, 2,  2, 2, 2], 0.07, None),
        (1_010_000, [ 0, 0, 1,  1, 1, 2], 0.00, None),
        (1_011_000, [ 0, 0, 0,  0, 0, 1], 0.00, None),
    ]
    # Sites file: pos\tref\talt\tmissingness\tdiagnostic_score\tsite_id
    chrom = "C_gar_LG28"
    with gzip.open(dose_dir / f"{chrom}.sites.tsv.gz", "wt") as f:
        for pos, _ds, miss, diag in site_table:
            miss_s = "NA" if miss is None else f"{miss:.4f}"
            diag_s = "NA" if diag is None else f"{diag:.4f}"
            site_id = f"{chrom}:{pos}:A>G"
            f.write(f"{pos}\tA\tG\t{miss_s}\t{diag_s}\t{site_id}\n")
    # Dosage file: one row per site, 6 cells per row, NA encoded literally
    with gzip.open(dose_dir / f"{chrom}.dosage.tsv.gz", "wt") as f:
        for _pos, doses, _m, _d in site_table:
            cells = ["NA" if d == -1 else str(d) for d in doses]
            f.write("\t".join(cells) + "\n")

    expect((dose_dir / f"{chrom}.sites.tsv.gz").is_file(),
           "sites.tsv.gz written")
    expect((dose_dir / f"{chrom}.dosage.tsv.gz").is_file(),
           "dosage.tsv.gz written")
    expect(samples_path.is_file(), "samples.tsv written")

    # ---------------------------------------------------------------------------
    section("Helpers — _maybe_dosage_int")
    expect(_maybe_dosage_int("0") == 0, "'0' -> 0")
    expect(_maybe_dosage_int("1") == 1, "'1' -> 1")
    expect(_maybe_dosage_int("2") == 2, "'2' -> 2")
    expect(_maybe_dosage_int("NA") == -1, "'NA' -> -1")
    expect(_maybe_dosage_int(".") == -1, "'.' -> -1")
    expect(_maybe_dosage_int("") == -1, "'' -> -1")
    expect(_maybe_dosage_int("-1") == -1, "'-1' -> -1")
    expect(_maybe_dosage_int("3") == -1, "out-of-range '3' -> -1")
    expect(_maybe_dosage_int("xyz") == -1, "non-numeric -> -1")
    expect(_maybe_dosage_int(" 1 ") == 1, "whitespace stripped")

    section("Helpers — _maybe_float")
    expect(_maybe_float("0.05") == 0.05, "'0.05' -> 0.05")
    expect(_maybe_float("NA") is None, "'NA' -> None")
    expect(_maybe_float(".") is None, "'.' -> None")
    expect(_maybe_float("") is None, "'' -> None")
    expect(_maybe_float("xyz") is None, "non-numeric -> None")

    section("Helpers — _select_even / _closest_index")
    # Closest index works
    expect(_closest_index([10, 20, 30, 40], 14) == 0, "closest 14 -> 10")
    expect(_closest_index([10, 20, 30, 40], 16) == 1, "closest 16 -> 20")
    expect(_closest_index([10, 20, 30, 40], 30) == 2, "exact 30 -> idx 2")
    expect(_closest_index([10, 20, 30, 40], 100) == 3, "OOB high -> last")
    expect(_closest_index([10, 20, 30, 40], 0) == 0, "OOB low -> first")
    expect(_closest_index([], 5) == 0, "empty -> 0 (defensive)")

    # Even selection of 4 from 12 positions
    pos = [p for p, _, _, _ in site_table]
    sel = _select_even(pos, 4)
    expect(len(sel) == 4, f"even select 4 from 12 -> 4 indices (got {len(sel)})")
    expect(sel == sorted(sel), "selection sorted ascending")
    expect(all(0 <= i < len(pos) for i in sel), "all indices in range")
    expect(sel[0] == 0 or pos[sel[0]] - pos[0] < (pos[-1] - pos[0]) / 4,
           "first selection near start")
    expect(sel[-1] == len(pos) - 1 or pos[-1] - pos[sel[-1]] < (pos[-1] - pos[0]) / 4,
           "last selection near end")
    # When n <= cap, returns all indices
    sel_small = _select_even([1, 2, 3], 10)
    expect(sel_small == [0, 1, 2], "small input returns all indices")

    # ---------------------------------------------------------------------------
    section("Request validation")
    # Valid
    r = DosageChunkReq(chrom=chrom, start=1_000_000, end=1_011_000, cap=500, mode="raw")
    expect(r.chrom == chrom, "valid request constructs")
    # end < start rejected
    bad = False
    try:
        DosageChunkReq(chrom=chrom, start=1_010_000, end=1_000_000)
    except Exception:
        bad = True
    expect(bad, "end < start rejected")
    # path traversal in chrom rejected
    bad = False
    try:
        DosageChunkReq(chrom="../etc/passwd", start=1, end=2)
    except Exception:
        bad = True
    expect(bad, "chrom with '/' rejected")
    bad = False
    try:
        DosageChunkReq(chrom="C_gar_LG28; rm -rf /", start=1, end=2)
    except Exception:
        bad = True
    expect(bad, "chrom with shell metachars rejected")
    # cap > hard cap rejected
    bad = False
    try:
        DosageChunkReq(chrom=chrom, start=1, end=2, cap=DEFAULT_HARD_CAP + 1)
    except Exception:
        bad = True
    expect(bad, "cap above hard-cap rejected")
    # invalid mode rejected
    bad = False
    try:
        DosageChunkReq(chrom=chrom, start=1, end=2, mode="random")
    except Exception:
        bad = True
    expect(bad, "mode='random' rejected (only raw/even on bridge)")

    # ---------------------------------------------------------------------------
    section("Renderer-contract round-trip — raw mode, full window")
    req = DosageChunkReq(chrom=chrom, start=1_000_000, end=1_011_000, cap=500, mode="raw")
    out = handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    expect(out["chrom"] == chrom, "chrom round-trip")
    expect(out["start_bp"] == 1_000_000, "start_bp round-trip")
    expect(out["end_bp"] == 1_011_000, "end_bp round-trip")
    expect(out["samples"] == samples, "samples in canonical order")
    expect(len(out["markers"]) == 12, f"12 markers (got {len(out['markers'])})")
    expect(len(out["dosage"]) == 12, f"12 dosage rows (got {len(out['dosage'])})")
    # Renderer contract: dosage[i].length == samples.length for every i
    expect(all(len(row) == 6 for row in out["dosage"]),
           "every dosage row has 6 cells (= n_samples)")
    # Renderer contract: missing -> -1, NEVER null
    has_null = any(v is None for row in out["dosage"] for v in row)
    expect(not has_null,
           "no null values in dosage matrix (missing must be -1)")
    expect(out["dosage"][2][0] == -1,
           f"site 2 col 0 is -1 (synthesized NA) — got {out['dosage'][2][0]}")
    # First and last positions
    expect(out["markers"][0]["pos_bp"] == 1_000_000, "first marker pos = 1_000_000")
    expect(out["markers"][-1]["pos_bp"] == 1_011_000, "last marker pos = 1_011_000")
    # missingness / diagnostic_score CAN be null (renderer contract for these)
    expect(out["markers"][0]["missingness"] == 0.0,
           "marker 0 missingness round-trip")
    expect(out["markers"][1]["diagnostic_score"] == 0.42,
           "marker 1 diagnostic_score round-trip")
    expect(out["markers"][0]["diagnostic_score"] is None,
           "missing diagnostic_score round-trips as null")
    # Meta
    expect(out["_meta"]["n_sites_total"] == 12, "_meta.n_sites_total = 12")
    expect(out["_meta"]["n_sites_returned"] == 12, "_meta.n_sites_returned = 12")
    expect(out["_meta"]["downsampled"] is False, "_meta.downsampled false")
    expect(out["_meta"]["mode"] == "raw", "_meta.mode = raw")

    # ---------------------------------------------------------------------------
    section("Window filter")
    req = DosageChunkReq(chrom=chrom, start=1_003_000, end=1_007_000, cap=500, mode="raw")
    out = handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    pos_returned = [m["pos_bp"] for m in out["markers"]]
    expect(pos_returned == [1_003_000, 1_004_000, 1_005_000, 1_006_000, 1_007_000],
           f"window filter returns 5 sites in [1_003_000, 1_007_000] (got {pos_returned})")
    expect(len(out["dosage"]) == 5, "5 dosage rows for the windowed query")
    # Dosage at site 1_003_000 was [0,0,1,2,2,2]
    expect(out["dosage"][0] == [0, 0, 1, 2, 2, 2],
           f"site 1_003_000 dosage [0,0,1,2,2,2] (got {out['dosage'][0]})")

    # Empty region returns empty markers/dosage
    req = DosageChunkReq(chrom=chrom, start=2_000_000, end=2_001_000, cap=500, mode="raw")
    out = handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    expect(out["markers"] == [] and out["dosage"] == [],
           "empty region returns empty markers + dosage")
    expect(out["_meta"]["n_sites_total"] == 0, "_meta.n_sites_total = 0 for empty region")

    # ---------------------------------------------------------------------------
    section("raw cap policy")
    req = DosageChunkReq(chrom=chrom, start=1_000_000, end=1_011_000, cap=5, mode="raw")
    threw = False
    detail = None
    try:
        handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    except HTTPException as e:
        threw = True
        detail = e.detail
    expect(threw, "raw mode + cap exceeded -> HTTPException")
    expect(detail and detail.get("error") == "raw_cap_exceeded",
           f"raw_cap_exceeded error code (got {detail and detail.get('error')!r})")
    expect(detail and detail.get("n_sites_total") == 12,
           "error includes n_sites_total")
    expect(detail and detail.get("max_sites") == 5,
           "error includes max_sites")
    expect(detail and "suggestion" in detail,
           "error includes a suggestion string")

    # ---------------------------------------------------------------------------
    section("even mode")
    req = DosageChunkReq(chrom=chrom, start=1_000_000, end=1_011_000, cap=4, mode="even")
    out = handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    expect(out["_meta"]["downsampled"] is True, "even mode flags downsampled")
    expect(out["_meta"]["n_sites_returned"] == 4,
           f"even mode returns cap=4 sites (got {out['_meta']['n_sites_returned']})")
    pos_returned = [m["pos_bp"] for m in out["markers"]]
    expect(pos_returned == sorted(pos_returned),
           f"even-mode positions monotonic (got {pos_returned})")
    expect(len(out["dosage"]) == 4, "4 dosage rows for even mode")
    # Bp-even: gaps should be roughly equal (within ~1 step). Total span 11000
    # / 4 bins = 2750 bp per bin. Adjacent gaps should be approximately that.
    gaps = [pos_returned[i + 1] - pos_returned[i] for i in range(len(pos_returned) - 1)]
    expect(all(2000 <= g <= 4500 for g in gaps),
           f"even-mode gaps roughly uniform (got {gaps})")

    # When cap >= n, even mode also returns everything
    req = DosageChunkReq(chrom=chrom, start=1_000_000, end=1_011_000, cap=500, mode="even")
    out = handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    expect(out["_meta"]["n_sites_returned"] == 12,
           "even mode w/ cap >= n returns all sites")
    expect(out["_meta"]["downsampled"] is False,
           "even mode w/ cap >= n is not flagged downsampled")

    # ---------------------------------------------------------------------------
    section("Header-row tolerance on dosage TSV")
    # Build a second chrom with a sample-ID header row in the dosage file.
    chrom2 = "C_gar_LG29"
    with gzip.open(dose_dir / f"{chrom2}.sites.tsv.gz", "wt") as f:
        for pos in (2_000_000, 2_001_000, 2_002_000):
            f.write(f"{pos}\tA\tG\tNA\tNA\t{chrom2}:{pos}:A>G\n")
    with gzip.open(dose_dir / f"{chrom2}.dosage.tsv.gz", "wt") as f:
        # Header row first
        f.write("\t".join(samples) + "\n")
        f.write("0\t0\t1\t1\t2\t2\n")
        f.write("0\t1\t1\t2\t2\t2\n")
        f.write("0\t0\t0\t1\t1\t2\n")
    req = DosageChunkReq(chrom=chrom2, start=2_000_000, end=2_002_000, cap=500, mode="raw")
    out = handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    expect(len(out["markers"]) == 3, "header-row variant: 3 markers")
    expect(out["dosage"][0] == [0, 0, 1, 1, 2, 2],
           f"header-row variant: row 0 = [0,0,1,1,2,2] (got {out['dosage'][0]})")
    expect(out["dosage"][1] == [0, 1, 1, 2, 2, 2],
           "header-row variant: row 1 OK")

    # ---------------------------------------------------------------------------
    section("Error paths")
    req = DosageChunkReq(chrom="C_gar_LG99", start=1, end=1000, cap=500, mode="raw")
    threw = False
    detail = None
    try:
        handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    except HTTPException as e:
        threw = True
        detail = e.detail
    expect(threw, "missing chrom raises HTTPException")
    expect(detail and detail.get("error") == "chrom_not_found",
           "chrom_not_found error code")

    # Region too large
    req = DosageChunkReq(chrom=chrom, start=1, end=DEFAULT_MAX_REGION_BP + 2,
                         cap=500, mode="raw")
    threw = False
    detail = None
    try:
        handle_dosage_chunk(req, dosage_dir=dose_dir, samples_path=samples_path)
    except HTTPException as e:
        threw = True
        detail = e.detail
    expect(threw, "region > MAX_REGION_BP raises HTTPException")
    expect(detail and detail.get("error") == "region_too_large",
           "region_too_large error code")

    # Missing samples file
    req = DosageChunkReq(chrom=chrom, start=1_000_000, end=1_001_000, cap=500, mode="raw")
    threw = False
    try:
        handle_dosage_chunk(req, dosage_dir=dose_dir,
                            samples_path=tmp / "no_such_samples.tsv")
    except HTTPException as e:
        threw = True
    expect(threw, "missing samples file raises HTTPException")

finally:
    if tmp.exists():
        shutil.rmtree(tmp)

# ---------------------------------------------------------------------------
print(f"\n=========================================")
print(f"  PASS: {passed}   FAIL: {failed}")
print(f"=========================================")
sys.exit(0 if failed == 0 else 1)
