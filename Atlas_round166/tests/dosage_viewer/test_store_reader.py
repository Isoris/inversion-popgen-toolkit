#!/usr/bin/env python3
"""
test_store_reader.py
---------------------
End-to-end test for dosage_viewer/store_reader.py.

Strategy: build a tiny store with 01_prepare_dosage_store.py, then use
store_reader to round-trip everything. This is integration-shaped: a
real Parquet file goes through pyarrow on the way out, then back in.

What's covered
--------------
- load_store_manifest  — reads manifest.json, validates schema_version
- load_samples         — reads samples.tsv into a list
- list_chroms          — chrom names in manifest order
- read_sites_in_window — bp filter; NaN missingness/diag round-trip as None
- read_dosage_for_positions — pos lookup; NA cells preserved as -1
- read_dosage_window   — bulk read for a region (used by variance + agg)
- get_chrom_length_bp / get_chrom_n_sites — small accessors
- error paths: missing files, schema-version mismatch

Run:
  python3 tests/dosage_viewer/test_store_reader.py
"""
from __future__ import annotations

import gzip
import json
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "dosage_viewer"))

import importlib
prep = importlib.import_module("01_prepare_dosage_store")
import store_reader as sr

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
section("Build a store")
tmp = Path(tempfile.mkdtemp(prefix="dv_reader_"))
try:
    base = tmp / "tsv"
    base.mkdir()
    samples = ["S0", "S1", "S2"]
    samples_path = tmp / "samples.tsv"
    samples_path.write_text("\n".join(samples) + "\n")

    # Two chroms, deliberate-content fixtures
    lg = [
        # (pos, ref, alt, miss, diag, doses_for_3_samples)
        (1000, "A", "G", 0.0, None, [0, 1, 2]),
        (2000, "C", "T", 0.1, 0.42, [0, 0, 1]),
        (3000, "G", "A", None, None, [-1, 1, 2]),
        (4000, "T", "C", 0.0, 0.91, [2, 2, 2]),
        (5000, "A", "T", 0.05, None, [1, 1, -1]),
    ]
    with gzip.open(base / "C_gar_LG28.sites.tsv.gz", "wt") as f:
        for pos, r, a, m, d, _ in lg:
            ms = "NA" if m is None else f"{m:.4f}"
            ds = "NA" if d is None else f"{d:.4f}"
            f.write(f"{pos}\t{r}\t{a}\t{ms}\t{ds}\t\n")
    with gzip.open(base / "C_gar_LG28.dosage.tsv.gz", "wt") as f:
        for *_, doses in lg:
            f.write("\t".join(["NA" if d == -1 else str(d) for d in doses]) + "\n")

    # Second chrom with 2 sites
    lg12 = [
        (10_000, "A", "G", None, None, [1, 0, 0]),
        (20_000, "C", "T", 0.0, None, [2, 1, 1]),
    ]
    with gzip.open(base / "C_gar_LG12.sites.tsv.gz", "wt") as f:
        for pos, r, a, m, d, _ in lg12:
            ms = "NA" if m is None else f"{m:.4f}"
            ds = "NA" if d is None else f"{d:.4f}"
            f.write(f"{pos}\t{r}\t{a}\t{ms}\t{ds}\t\n")
    with gzip.open(base / "C_gar_LG12.dosage.tsv.gz", "wt") as f:
        for *_, doses in lg12:
            f.write("\t".join(str(d) for d in doses) + "\n")

    out = tmp / "store"
    prep.build_store(base=base, out=out, samples_path=samples_path)
    expect((out / "manifest.json").is_file(), "store built")

    # ---------------------------------------------------------------------------
    section("load_store_manifest")
    m = sr.load_store_manifest(out)
    expect(m["schema_version"] == "dosage_store_v1", "schema_version OK")
    expect(m["n_samples"] == 3, "n_samples = 3")
    chrom_names = [c["name"] for c in m["chroms"]]
    expect(set(chrom_names) == {"C_gar_LG12", "C_gar_LG28"},
           f"manifest has both chroms (got {chrom_names})")

    # ---------------------------------------------------------------------------
    section("load_samples + list_chroms")
    expect(sr.load_samples(out) == samples, "samples list round-trips")
    expect(set(sr.list_chroms(out)) == {"C_gar_LG12", "C_gar_LG28"},
           "list_chroms returns both")

    # ---------------------------------------------------------------------------
    section("read_sites_in_window — full window")
    pos, miss, diag, sids = sr.read_sites_in_window(out, "C_gar_LG28", 0, 1_000_000)
    expect(pos == [1000, 2000, 3000, 4000, 5000],
           f"all 5 positions returned (got {pos})")
    expect(miss[0] == 0.0 and miss[2] is None,
           "missingness: NaN -> None, finite values preserved")
    expect(diag[1] is not None and abs(diag[1] - 0.42) < 1e-5,
           f"diag[1] = 0.42 (got {diag[1]})")
    expect(diag[0] is None and diag[2] is None,
           "diag: NA round-trips as None")
    expect(all(s.startswith("C_gar_LG28:") for s in sids),
           "site_ids look canonical")

    # ---------------------------------------------------------------------------
    section("read_sites_in_window — narrow window")
    pos, miss, diag, sids = sr.read_sites_in_window(out, "C_gar_LG28", 2500, 4500)
    expect(pos == [3000, 4000],
           f"narrow window: 2 sites in [2500, 4500] (got {pos})")

    # Empty window
    pos, _, _, _ = sr.read_sites_in_window(out, "C_gar_LG28", 100_000, 200_000)
    expect(pos == [], "out-of-range window returns empty")

    # ---------------------------------------------------------------------------
    section("read_dosage_for_positions")
    pos, _, _, _ = sr.read_sites_in_window(out, "C_gar_LG28", 0, 1_000_000)
    rows = sr.read_dosage_for_positions(out, "C_gar_LG28", pos)
    expect(len(rows) == 5, f"5 dosage rows (got {len(rows)})")
    expect(all(len(r) == 3 for r in rows), "every row has n_samples=3 cells")
    expect(rows[0] == [0, 1, 2], f"row 0 dosage [0,1,2] (got {rows[0]})")
    expect(rows[2] == [-1, 1, 2], f"row 2 NA preserved as -1 (got {rows[2]})")
    expect(rows[4] == [1, 1, -1], "row 4 NA tail preserved")

    # Subset by position list
    rows = sr.read_dosage_for_positions(out, "C_gar_LG28", [2000, 4000])
    expect(rows == [[0, 0, 1], [2, 2, 2]],
           f"subset positions return only requested rows (got {rows})")

    # Position not in store -> all-NA row
    rows = sr.read_dosage_for_positions(out, "C_gar_LG28", [999_999])
    expect(rows == [[-1, -1, -1]],
           "position not in store -> all-NA row")

    # Empty list
    expect(sr.read_dosage_for_positions(out, "C_gar_LG28", []) == [],
           "empty position list returns []")

    # ---------------------------------------------------------------------------
    section("read_dosage_window")
    pos, rows = sr.read_dosage_window(out, "C_gar_LG28", 2500, 4500)
    expect(pos == [3000, 4000], "window positions match")
    expect(rows == [[-1, 1, 2], [2, 2, 2]],
           f"window dosage rows match (got {rows})")

    # Empty window
    pos, rows = sr.read_dosage_window(out, "C_gar_LG28", 100_000, 200_000)
    expect(pos == [] and rows == [], "empty window returns empty")

    # ---------------------------------------------------------------------------
    section("Per-chrom helpers")
    expect(sr.get_chrom_length_bp(out, "C_gar_LG28") == 5000,
           "length_bp for LG28 = 5000")
    expect(sr.get_chrom_length_bp(out, "C_gar_LG12") == 20_000,
           "length_bp for LG12 = 20000")
    expect(sr.get_chrom_n_sites(out, "C_gar_LG28") == 5,
           "n_sites for LG28 = 5")
    expect(sr.get_chrom_n_sites(out, "C_gar_LG_NEVER") == 0,
           "missing chrom -> 0")

    # ---------------------------------------------------------------------------
    section("Error paths")
    threw = False
    try:
        sr.load_store_manifest(tmp / "no_such_dir")
    except FileNotFoundError:
        threw = True
    expect(threw, "missing manifest -> FileNotFoundError")

    threw = False
    try:
        sr.read_sites_in_window(out, "C_gar_LG_NEVER", 0, 100)
    except FileNotFoundError:
        threw = True
    expect(threw, "missing chrom -> FileNotFoundError")

    threw = False
    try:
        sr.read_dosage_for_positions(out, "C_gar_LG_NEVER", [1000])
    except FileNotFoundError:
        threw = True
    expect(threw, "missing chrom dosage -> FileNotFoundError")

    # Schema-version mismatch
    bad_root = tmp / "store_bad"
    bad_root.mkdir()
    (bad_root / "manifest.json").write_text(json.dumps({
        "schema_version": "wrong_v9", "n_samples": 0, "chroms": []
    }))
    threw = False
    try:
        sr.load_store_manifest(bad_root)
    except ValueError as e:
        threw = "schema_version" in str(e)
    expect(threw, "schema_version mismatch -> ValueError")

finally:
    if tmp.exists():
        shutil.rmtree(tmp)

# ---------------------------------------------------------------------------
print(f"\n=========================================")
print(f"  PASS: {passed}   FAIL: {failed}")
print(f"=========================================")
sys.exit(0 if failed == 0 else 1)
