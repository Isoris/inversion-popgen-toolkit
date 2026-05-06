#!/usr/bin/env python3
"""
test_prepare_dosage_store.py
-----------------------------
End-to-end test for dosage_viewer/01_prepare_dosage_store.py.

What's exercised
----------------
- Synthesise tiny <chrom>.sites.tsv.gz + <chrom>.dosage.tsv.gz fixtures
  for two chromosomes with known content.
- Run build_store() against them.
- Verify:
    * dosage_store/<chrom>/sites.parquet & dosage.parquet exist
    * sites.parquet schema (pos int64, ref/alt strings, missingness +
      diagnostic_score float, site_id string)
    * dosage.parquet schema (pos + n_samples int8 columns)
    * data round-trips (positions, dosage cells, NA preserved as -1)
    * NA in missingness/diagnostic_score round-trips as null
    * dosage_store/samples.tsv contains the canonical sample IDs in order
    * dosage_store/manifest.json schema version + chroms + n_samples
    * second build with no source change -> chroms skipped (idempotency)
    * second build after touching one source -> only that chrom rebuilds
    * inline-header dosage TSV is auto-detected
    * unsorted positions are re-sorted in the parquet output
    * inconsistent n_samples across chroms -> ValueError
    * missing source files -> handled gracefully

Run:
  python3 tests/dosage_viewer/test_prepare_dosage_store.py
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "dosage_viewer"))

import importlib
prep = importlib.import_module("01_prepare_dosage_store")
import pyarrow as pa
import pyarrow.parquet as pq

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
# Synthesize a tiny store
# ---------------------------------------------------------------------------
section("Synthesizing fixtures")

tmp = Path(tempfile.mkdtemp(prefix="dv_prep_"))
try:
    base = tmp / "04_dosage_by_chr"
    base.mkdir()

    # 4 samples, 2 chroms (LG28 with 6 sites, LG12 with 4 sites)
    samples = ["CGA_A", "CGA_B", "CGA_C", "CGA_D"]
    samples_path = tmp / "samples.tsv"
    samples_path.write_text("\n".join(samples) + "\n")

    # LG28: 6 sites, mixture of dosages incl one NA cell
    lg28_sites = [
        # (pos, ref, alt, missingness, diag, [dosages])
        (1_000_000, "A", "G", 0.00, None, [0, 1, 2, -1]),
        (1_001_000, "C", "T", 0.05, 0.42, [0, 0, 1, 2]),
        (1_002_000, "G", "A", None, None, [1, 1, 2, 2]),
        (1_003_000, "A", "C", 0.01, 0.91, [0, 0, 0, 0]),
        (1_004_000, "T", "G", 0.00, None, [2, 2, 2, 2]),
        (1_005_000, "C", "A", 0.02, None, [0, 1, 2, -1]),
    ]
    with gzip.open(base / "C_gar_LG28.sites.tsv.gz", "wt") as f:
        for pos, ref, alt, miss, diag, _ in lg28_sites:
            ms = "NA" if miss is None else f"{miss:.4f}"
            ds = "NA" if diag is None else f"{diag:.4f}"
            f.write(f"{pos}\t{ref}\t{alt}\t{ms}\t{ds}\t\n")
    with gzip.open(base / "C_gar_LG28.dosage.tsv.gz", "wt") as f:
        # Headerless (canonical case)
        for _, _, _, _, _, doses in lg28_sites:
            cells = ["NA" if d == -1 else str(d) for d in doses]
            f.write("\t".join(cells) + "\n")

    # LG12: 4 sites with INLINE HEADER row in dosage TSV
    lg12_sites = [
        (2_000_000, "T", "A", 0.00, None, [0, 1, 1, 0]),
        (2_001_000, "G", "C", 0.10, 0.20, [1, 1, 2, 2]),
        (2_002_000, "A", "T", None, None, [0, 0, 0, 1]),
        (2_003_000, "C", "G", 0.00, None, [2, 2, 2, 2]),
    ]
    with gzip.open(base / "C_gar_LG12.sites.tsv.gz", "wt") as f:
        for pos, ref, alt, miss, diag, _ in lg12_sites:
            ms = "NA" if miss is None else f"{miss:.4f}"
            ds = "NA" if diag is None else f"{diag:.4f}"
            sid = f"LG12:{pos}:{ref}>{alt}_explicit"
            f.write(f"{pos}\t{ref}\t{alt}\t{ms}\t{ds}\t{sid}\n")
    with gzip.open(base / "C_gar_LG12.dosage.tsv.gz", "wt") as f:
        f.write("\t".join(samples) + "\n")    # inline header
        for _, _, _, _, _, doses in lg12_sites:
            cells = ["NA" if d == -1 else str(d) for d in doses]
            f.write("\t".join(cells) + "\n")

    expect((base / "C_gar_LG28.sites.tsv.gz").is_file(), "LG28 sites.tsv.gz written")
    expect((base / "C_gar_LG12.dosage.tsv.gz").is_file(), "LG12 dosage.tsv.gz written")

    # ---------------------------------------------------------------------------
    section("First build")
    out = tmp / "dosage_store"
    manifest = prep.build_store(base=base, out=out, samples_path=samples_path)

    expect(manifest["schema_version"] == "dosage_store_v1", "schema_version on store manifest")
    expect(manifest["n_samples"] == 4, "n_samples = 4 in store manifest")
    expect(manifest["samples"] == samples, "samples list canonical order")
    expect(manifest["n_built"] == 2 and manifest["n_skipped"] == 0,
           "first build: 2 built, 0 skipped")

    # Per-chrom outputs
    expect((out / "C_gar_LG28" / "sites.parquet").is_file(), "LG28 sites.parquet written")
    expect((out / "C_gar_LG28" / "dosage.parquet").is_file(), "LG28 dosage.parquet written")
    expect((out / "C_gar_LG28" / "_chrom_manifest.json").is_file(), "LG28 marker written")
    expect((out / "C_gar_LG12" / "sites.parquet").is_file(), "LG12 sites.parquet written")
    expect((out / "C_gar_LG12" / "dosage.parquet").is_file(), "LG12 dosage.parquet written")
    expect((out / "samples.tsv").is_file(), "store-level samples.tsv written")
    expect((out / "samples.tsv").read_text().strip().splitlines() == samples,
           "samples.tsv content matches canonical order")

    # ---------------------------------------------------------------------------
    section("sites.parquet round-trip — LG28")
    sites28 = pq.read_table(out / "C_gar_LG28" / "sites.parquet")
    cols = sites28.column_names
    expect(cols == ["pos", "ref", "alt", "missingness", "diagnostic_score", "site_id"],
           f"sites.parquet column order canonical (got {cols})")
    pos28 = sites28.column("pos").to_pylist()
    expect(pos28 == [1_000_000, 1_001_000, 1_002_000,
                     1_003_000, 1_004_000, 1_005_000],
           "LG28 positions round-trip")
    miss28 = sites28.column("missingness").to_pylist()
    expect(miss28[0] == 0.0,  "miss[0] = 0.0")
    expect(miss28[2] is None, "miss[2] NA -> None")
    diag28 = sites28.column("diagnostic_score").to_pylist()
    expect(diag28[1] is not None and abs(diag28[1] - 0.42) < 1e-5,
           "diag[1] = 0.42")
    expect(diag28[0] is None, "diag[0] NA -> None")
    site_ids28 = sites28.column("site_id").to_pylist()
    expect(site_ids28[0].startswith("C_gar_LG28:1000000:A>G"),
           f"synthesised site_id when source was empty (got {site_ids28[0]!r})")

    # LG12 used explicit site_id strings
    sites12 = pq.read_table(out / "C_gar_LG12" / "sites.parquet")
    sid12 = sites12.column("site_id").to_pylist()
    expect(sid12[0] == "LG12:2000000:T>A_explicit",
           "LG12 explicit site_id preserved")

    # ---------------------------------------------------------------------------
    section("dosage.parquet round-trip — LG28")
    dos28 = pq.read_table(out / "C_gar_LG28" / "dosage.parquet")
    expect(dos28.column_names[0] == "pos",
           "first dosage.parquet column is 'pos'")
    expect(dos28.column_names[1:] == samples,
           f"sample column names match canonical IDs (got {dos28.column_names[1:]})")
    # int8 storage check
    expect(dos28.schema.field("CGA_A").type == pa.int8(),
           "sample columns are int8")
    expect(dos28.schema.field("pos").type == pa.int64(),
           "pos column is int64")
    # Cell content
    rows = list(zip(*[dos28.column(c).to_pylist() for c in dos28.column_names]))
    # row 0: pos=1_000_000, doses=[0,1,2,-1]
    expect(rows[0] == (1_000_000, 0, 1, 2, -1),
           f"row 0 round-trips with NA -> -1 (got {rows[0]})")
    expect(rows[2] == (1_002_000, 1, 1, 2, 2),
           "row 2 round-trips")
    expect(rows[5] == (1_005_000, 0, 1, 2, -1),
           "row 5 round-trips with NA")

    # ---------------------------------------------------------------------------
    section("dosage.parquet — inline header detection (LG12)")
    dos12 = pq.read_table(out / "C_gar_LG12" / "dosage.parquet")
    expect(dos12.column_names[1:] == samples,
           "LG12: inline header used as sample IDs")
    rows12 = list(zip(*[dos12.column(c).to_pylist() for c in dos12.column_names]))
    expect(rows12[0] == (2_000_000, 0, 1, 1, 0),
           "LG12 row 0 round-trips (header skipped)")
    expect(len(rows12) == 4, "LG12 has 4 data rows (header skipped)")

    # ---------------------------------------------------------------------------
    section("Idempotency — second build is all-skip")
    manifest2 = prep.build_store(base=base, out=out, samples_path=samples_path)
    expect(manifest2["n_built"] == 0 and manifest2["n_skipped"] == 2,
           f"second build: 0 built, 2 skipped (got {manifest2['n_built']}/{manifest2['n_skipped']})")

    # ---------------------------------------------------------------------------
    section("Selective rebuild — touch one source, rebuild only that chrom")
    # Bump LG28 sites mtime
    new_t = time.time() + 5
    os.utime(base / "C_gar_LG28.sites.tsv.gz", (new_t, new_t))
    manifest3 = prep.build_store(base=base, out=out, samples_path=samples_path)
    expect(manifest3["n_built"] == 1 and manifest3["n_skipped"] == 1,
           f"selective rebuild: 1 built, 1 skipped (got {manifest3['n_built']}/{manifest3['n_skipped']})")

    # ---------------------------------------------------------------------------
    section("Force rebuild")
    manifest4 = prep.build_store(base=base, out=out,
                                  samples_path=samples_path, force=True)
    expect(manifest4["n_built"] == 2 and manifest4["n_skipped"] == 0,
           "force=True rebuilds all chroms")

    # ---------------------------------------------------------------------------
    section("Unsorted positions — re-sort happens in memory")
    base2 = tmp / "04_dosage_by_chr_unsorted"
    base2.mkdir()
    # 3 sites in REVERSE order
    rev_sites = [(3_000_000 - 1000 * i, "A", "G", None, None, [i, i, i, i])
                 for i in range(3)]
    # Becomes: pos=3_000_000 first (i=0, all 0s); pos=2_999_000 next (all 1s);
    # pos=2_998_000 last (all 2s). After sort: ascending, so doses become
    # [2,2,2,2] → [1,1,1,1] → [0,0,0,0].
    with gzip.open(base2 / "C_gar_LG07.sites.tsv.gz", "wt") as f:
        for pos, ref, alt, miss, diag, _ in rev_sites:
            f.write(f"{pos}\t{ref}\t{alt}\tNA\tNA\t\n")
    with gzip.open(base2 / "C_gar_LG07.dosage.tsv.gz", "wt") as f:
        for _, _, _, _, _, doses in rev_sites:
            f.write("\t".join(str(d) for d in doses) + "\n")
    out2 = tmp / "store_unsorted"
    prep.build_store(base=base2, out=out2, samples_path=samples_path)

    sites_re = pq.read_table(out2 / "C_gar_LG07" / "sites.parquet")
    pos_re = sites_re.column("pos").to_pylist()
    expect(pos_re == sorted(pos_re),
           f"re-sorted positions ascending (got {pos_re})")
    dos_re = pq.read_table(out2 / "C_gar_LG07" / "dosage.parquet")
    rows_re = list(zip(*[dos_re.column(c).to_pylist() for c in dos_re.column_names]))
    # After sort: 2_998_000 -> originally i=2, doses=[2,2,2,2]
    expect(rows_re[0] == (2_998_000, 2, 2, 2, 2),
           "first row of re-sorted dosage matches the LOWEST pos (originally last)")
    expect(rows_re[-1] == (3_000_000, 0, 0, 0, 0),
           "last row of re-sorted dosage matches the HIGHEST pos (originally first)")

    # ---------------------------------------------------------------------------
    section("Inconsistent n_samples across chroms — error")
    base3 = tmp / "04_dosage_by_chr_bad"
    base3.mkdir()
    # LG10: 4 samples; LG11: 5 samples
    with gzip.open(base3 / "C_gar_LG10.sites.tsv.gz", "wt") as f:
        f.write("100\tA\tG\tNA\tNA\t\n200\tC\tT\tNA\tNA\t\n")
    with gzip.open(base3 / "C_gar_LG10.dosage.tsv.gz", "wt") as f:
        f.write("0\t1\t2\t-1\n1\t1\t2\t2\n")
    with gzip.open(base3 / "C_gar_LG11.sites.tsv.gz", "wt") as f:
        f.write("100\tA\tG\tNA\tNA\t\n200\tC\tT\tNA\tNA\t\n")
    with gzip.open(base3 / "C_gar_LG11.dosage.tsv.gz", "wt") as f:
        f.write("0\t1\t2\t1\t2\n1\t1\t2\t2\t1\n")
    out3 = tmp / "store_bad"
    threw = False
    try:
        prep.build_store(base=base3, out=out3, samples_path=samples_path)
    except ValueError as e:
        msg = str(e).lower()
        # Accept either form of the error: "sample-count mismatch" (caught
        # at _read_dosage when samples_path width disagrees with the dosage
        # row width) or "n_samples ... != earlier chroms" (caught at the
        # cross-chrom check). Both indicate the same inconsistency.
        threw = ("sample" in msg and "mismatch" in msg) or "n_samples" in msg
    expect(threw, "inconsistent n_samples raises ValueError")

    # ---------------------------------------------------------------------------
    section("Missing source — graceful skip")
    base4 = tmp / "04_dosage_by_chr_partial"
    base4.mkdir()
    # Only sites file, no dosage file
    with gzip.open(base4 / "C_gar_LG20.sites.tsv.gz", "wt") as f:
        f.write("100\tA\tG\tNA\tNA\t\n")
    out4 = tmp / "store_partial"
    threw = False
    try:
        prep.build_store(base=base4, out=out4, samples_path=samples_path)
    except ValueError:
        threw = True
    expect(threw, "no usable chroms (sites without matching dosage) raises ValueError")

    # ---------------------------------------------------------------------------
    section("Empty sites file — error")
    base5 = tmp / "04_dosage_by_chr_empty"
    base5.mkdir()
    with gzip.open(base5 / "C_gar_LG30.sites.tsv.gz", "wt") as f:
        pass
    with gzip.open(base5 / "C_gar_LG30.dosage.tsv.gz", "wt") as f:
        pass
    out5 = tmp / "store_empty"
    threw = False
    try:
        prep.build_store(base=base5, out=out5, samples_path=samples_path)
    except ValueError:
        threw = True
    expect(threw, "empty sites file raises ValueError")

finally:
    if tmp.exists():
        shutil.rmtree(tmp)

# ---------------------------------------------------------------------------
print(f"\n=========================================")
print(f"  PASS: {passed}   FAIL: {failed}")
print(f"=========================================")
sys.exit(0 if failed == 0 else 1)
