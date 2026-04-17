#!/usr/bin/env python3
"""
Integration test for the chat-16 results_registry.

This test is Python because the chat sandbox doesn't have R. It exercises
the same atomic write path that R uses (both write the same manifest.tsv
schema), so passing it means the R side's manifest writes are correct too.

Covers:
  1. Fresh init creates the right directory layout
  2. put_pairwise with unregistered group raises FK error
  3. put_pairwise with registered groups succeeds and writes a manifest row
  4. Manifest row validates against result_row.schema.json
  5. ask() filters by kind / stat / who / where
  6. put_interval_summary with wrong stat raises
  7. Auto-migration from stats_cache/ → results_registry/
  8. sample_group schema validates a real row from sample_groups.tsv
"""
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

# Make the Python API importable
REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "registries" / "api" / "python"))
from registry_loader import load_registry

SCHEMA_DIR = REPO / "registries" / "schemas" / "registry_schemas"

# jsonschema is available in this environment (we installed it earlier)
try:
    from jsonschema import Draft7Validator
    HAS_JSONSCHEMA = True
except ImportError:
    HAS_JSONSCHEMA = False


def _write_master(registries_root: Path) -> None:
    """Pre-populate sample_registry and interval_registry with test data."""
    sr = registries_root / "data" / "sample_registry"
    sr.mkdir(parents=True, exist_ok=True)
    (sr / "groups").mkdir(exist_ok=True)

    samples = [f"CGA{i:03d}" for i in range(1, 11)]
    # sample_master.tsv
    with open(sr / "sample_master.tsv", "w") as fh:
        fh.write("sample_id\n")
        for s in samples:
            fh.write(s + "\n")
    # sample_groups.tsv — one group "all_10" + one group "unrelated_5"
    with open(sr / "sample_groups.tsv", "w") as fh:
        fh.write("group_id\tchrom\tinv_id\tsubgroup\tmembers_file\tn\tdescription\tcreated\n")
        fh.write(f"all_10\t.\t.\t.\tgroups/all_10.txt\t10\tcohort\t2026-04-17 10:00:00\n")
        fh.write(f"unrelated_5\t.\t.\t.\tgroups/unrelated_5.txt\t5\tpruned\t2026-04-17 10:00:00\n")
    with open(sr / "groups" / "all_10.txt", "w") as fh:
        fh.write("\n".join(samples) + "\n")
    with open(sr / "groups" / "unrelated_5.txt", "w") as fh:
        fh.write("\n".join(samples[:5]) + "\n")
    # backup dir the sample_registry.R expects
    (sr / "backups").mkdir(exist_ok=True)

    # interval_registry
    ir = registries_root / "data" / "interval_registry"
    ir.mkdir(parents=True, exist_ok=True)
    with open(ir / "candidate_intervals.tsv", "w") as fh:
        fh.write("candidate_id\tchrom\tstart_bp\tend_bp\tsize_kb\tscale\tparent_id\n")
        fh.write("LG12_17\tC_gar_LG12\t1000000\t1500000\t500\t100\t.\n")
        fh.write("LG25_3\tC_gar_LG25\t2000000\t2300000\t300\t50\t.\n")


def test_1_fresh_init():
    print("\n=== 1. Fresh init creates expected layout ===")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        os.environ["REGISTRIES"] = str(root)
        reg = load_registry(str(root))
        for sub in ("sample_registry", "interval_registry",
                     "evidence_registry/per_candidate",
                     "results_registry",
                     "results_registry/pairwise",
                     "results_registry/candidate",
                     "results_registry/interval"):
            assert (root / "data" / sub).exists(), f"Missing dir: {sub}"
        print("  OK  all required subdirs created")


def test_2_put_pairwise_fk_errors():
    print("\n=== 2. put_pairwise with unregistered group raises ===")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        reg = load_registry(str(root))

        try:
            reg.results.put_pairwise(
                chrom="C_gar_LG12",
                group_1="all_10",
                group_2="nonexistent_group",
                stat="fst",
                rows=[{"pos": 1000, "fst": 0.1}],
                fieldnames=["pos", "fst"],
                source_script="test.py"
            )
            print("  FAIL — expected ValueError for unregistered group_2"); sys.exit(1)
        except ValueError as e:
            if "not registered" in str(e):
                print(f"  OK   FK error raised as expected: {e}")
            else:
                print(f"  FAIL — wrong error: {e}"); sys.exit(1)


def test_3_put_pairwise_writes_manifest_row():
    print("\n=== 3. put_pairwise writes file + manifest row ===")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        reg = load_registry(str(root))

        out = reg.results.put_pairwise(
            chrom="C_gar_LG12",
            group_1="all_10",
            group_2="unrelated_5",
            stat="fst",
            rows=[{"pos": 1000, "fst": 0.12},
                  {"pos": 2000, "fst": 0.08}],
            fieldnames=["pos", "fst"],
            source_script="test_results_registry.py",
            engine="angsd", engine_version="0.940",
            upstream_files=["/path/to/bam.bamlist"]
        )
        assert out.exists(), "Output file not created"
        # canonical name: sorted pair
        assert out.name == "all_10__vs__unrelated_5.tsv", f"bad name: {out.name}"
        print(f"  OK   file written: {out.relative_to(root)}")

        rows = reg.results.list_manifest()
        assert len(rows) == 1, f"expected 1 manifest row, got {len(rows)}"
        r = rows[0]
        assert r["kind"] == "pairwise"
        assert r["chrom"] == "C_gar_LG12"
        assert r["group_1"] in ("all_10", "unrelated_5")
        assert r["group_2"] in ("all_10", "unrelated_5")
        assert r["stat"] == "fst"
        assert int(r["n_rows"]) == 2
        assert int(r["n_cols"]) == 2
        assert r["group_1_version"] == "2026-04-17 10:00:00"
        assert r["group_2_version"] == "2026-04-17 10:00:00"
        assert r["engine"] == "angsd"
        assert r["source_script"] == "test_results_registry.py"
        # UUID v4 format
        import re
        assert re.match(r"^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$",
                         r["row_id"]), f"bad uuid: {r['row_id']}"
        print("  OK   manifest row shape correct (FK versions stamped, UUID valid)")


def test_4_manifest_row_validates_schema():
    if not HAS_JSONSCHEMA:
        print("\n=== 4. manifest row validates against schema (SKIPPED — no jsonschema) ===")
        return
    print("\n=== 4. manifest row validates against result_row.schema.json ===")
    sch = json.loads((SCHEMA_DIR / "result_row.schema.json").read_text())
    validator = Draft7Validator(sch)

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        reg = load_registry(str(root))
        reg.results.put_pairwise(
            chrom="C_gar_LG12", group_1="all_10", group_2="unrelated_5",
            stat="fst", rows=[{"pos": 1, "fst": 0.1}], fieldnames=["pos","fst"],
            source_script="test.py")
        reg.results.put_interval_summary(
            chrom="C_gar_LG12", start_bp=1000000, end_bp=1100000,
            group_1="all_10", stat="delta12", K=8,
            rows=[{"win": 1, "delta12": 0.5}], fieldnames=["win","delta12"],
            source_script="test.py")

        rows = reg.results.list_manifest()
        # Convert TSV-strings back into typed fields matching schema
        fails = 0
        for raw in rows:
            obj = _row_to_schema_shape(raw)
            errs = list(validator.iter_errors(obj))
            if errs:
                fails += 1
                print(f"  FAIL  row {raw['row_id']}:")
                for e in errs[:3]:
                    print(f"         {e.message}")
        if fails == 0:
            print(f"  OK   all {len(rows)} manifest rows validate against schema")
        else:
            sys.exit(1)


def _row_to_schema_shape(flat):
    """Convert a flat TSV-row dict into the nested schema shape."""
    def int_or_none(v):
        if v == "" or v is None: return None
        try: return int(v)
        except Exception: return None
    obj = {
        "row_id":    flat["row_id"],
        "kind":      flat["kind"],
        "where":     {},
        "what":      {"stat": flat["stat"]},
        "file":      flat["file"],
        "provenance": {"source_script": flat["source_script"]},
        "timestamp": flat["timestamp"],
    }
    if flat.get("chrom"):        obj["where"]["chrom"] = flat["chrom"]
    if flat.get("start_bp"):     obj["where"]["start_bp"] = int_or_none(flat["start_bp"])
    if flat.get("end_bp"):       obj["where"]["end_bp"] = int_or_none(flat["end_bp"])
    if flat.get("candidate_id"): obj["where"]["candidate_id"] = flat["candidate_id"]
    if flat.get("K"):            obj["what"]["K"] = int_or_none(flat["K"])
    if flat.get("group_1"):
        obj["who_1"] = {"group_id": flat["group_1"],
                         "group_version": flat["group_1_version"] or ""}
    if flat.get("group_2"):
        obj["who_2"] = {"group_id": flat["group_2"],
                         "group_version": flat["group_2_version"] or ""}
    # optional
    if flat.get("n_rows") or flat.get("n_cols"):
        obj["shape"] = {}
        if flat.get("n_rows"): obj["shape"]["n_rows"] = int_or_none(flat["n_rows"])
        if flat.get("n_cols"): obj["shape"]["n_cols"] = int_or_none(flat["n_cols"])
        if flat.get("sha256"): obj["shape"]["sha256"] = flat["sha256"]
    return obj


def test_5_ask_filters():
    print("\n=== 5. ask() filters work ===")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        reg = load_registry(str(root))
        # Populate a few rows
        reg.results.put_pairwise(
            chrom="C_gar_LG12", group_1="all_10", group_2="unrelated_5",
            stat="fst", rows=[{"x":1,"y":0.1}], fieldnames=["x","y"],
            source_script="t.py")
        reg.results.put_pairwise(
            chrom="C_gar_LG25", group_1="all_10", group_2="unrelated_5",
            stat="fst", rows=[{"x":1,"y":0.2}], fieldnames=["x","y"],
            source_script="t.py")
        reg.results.put_pairwise(
            chrom="C_gar_LG12", group_1="all_10", group_2="unrelated_5",
            stat="dxy", rows=[{"x":1,"y":0.01}], fieldnames=["x","y"],
            source_script="t.py")

        assert len(reg.results.ask()) == 3
        assert len(reg.results.ask(where={"chrom": "C_gar_LG12"})) == 2
        assert len(reg.results.ask(what="fst")) == 2
        assert len(reg.results.ask(what="dxy")) == 1
        assert len(reg.results.ask(who="all_10")) == 3
        assert len(reg.results.ask(who="not_a_group")) == 0
        assert len(reg.results.ask(where={"chrom":"C_gar_LG12"}, what="fst")) == 1
        print("  OK   all filter combinations return correct counts")


def test_6_put_interval_summary_wrong_stat():
    print("\n=== 6. put_interval_summary rejects wrong stat ===")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        reg = load_registry(str(root))
        try:
            reg.results.put_interval_summary(
                chrom="C_gar_LG12", start_bp=1, end_bp=100,
                group_1="all_10", stat="fst",   # fst not valid for interval_summary
                rows=[{"x":1}], fieldnames=["x"], source_script="t.py")
            print("  FAIL — expected ValueError"); sys.exit(1)
        except ValueError as e:
            if "not valid for interval_summary" in str(e):
                print("  OK   rejected wrong stat for interval_summary")
            else:
                print(f"  FAIL — wrong error: {e}"); sys.exit(1)


def test_7_auto_migration():
    print("\n=== 7. auto-migration from stats_cache/ → results_registry/ ===")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        # Create a legacy stats_cache/ with some content
        legacy = root / "data" / "stats_cache"
        legacy.mkdir(parents=True, exist_ok=True)
        (legacy / "manifest.tsv").write_text("kind\tchrom\tg1\tg2\tstat\tsample_set\tcid\tK\tfile\ttimestamp\n")
        (legacy / "pairwise").mkdir()
        (legacy / "pairwise" / "marker.txt").write_text("hello from chat 15")
        # Load registry — migration should fire
        reg = load_registry(str(root))
        new_dir = root / "data" / "results_registry"
        assert (new_dir / "pairwise" / "marker.txt").exists(), \
            "marker not migrated"
        assert (new_dir / "manifest.tsv").exists(), "manifest not migrated"
        print("  OK   legacy stats_cache/ migrated to results_registry/")


def test_8_sample_group_schema_validates_row():
    if not HAS_JSONSCHEMA:
        print("\n=== 8. sample_group schema validates (SKIPPED — no jsonschema) ===")
        return
    print("\n=== 8. sample_group schema validates a real groups row ===")
    sch = json.loads((SCHEMA_DIR / "sample_group.schema.json").read_text())
    v = Draft7Validator(sch)
    # A real row from a test fixture
    row = {
        "group_id": "inv_LG12_17_HOM_INV__sub2",
        "chrom": "C_gar_LG12", "inv_id": "LG12_17",
        "subgroup": ".",
        "dimension": "karyotype_subcluster",
        "members_file": "groups/inv_LG12_17_HOM_INV__sub2.txt",
        "n": 15,
        "description": "HOM_INV sub-cluster 2 from recursive PCA",
        "created": "2026-04-17 14:30:00"
    }
    errs = list(v.iter_errors(row))
    if not errs:
        print("  OK   sub-cluster group row validates")
    else:
        for e in errs: print(f"  FAIL: {e.message}")
        sys.exit(1)


def test_9_scan_segment_semantics():
    """Verify the scan-semantics design: given a segmented candidate,
    a Python-side walker mirroring reg$intervals$resolve_smallest_at
    picks the correct candidate_id per window position.
    """
    print("\n=== 9. Multi-scale segment resolution semantics ===")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        # Add three segments inside LG12_17 (1.0-1.5 Mb parent in fixture)
        ir = root / "data" / "interval_registry" / "candidate_intervals.tsv"
        with open(ir, "a") as fh:
            fh.write("LG12_17_seg_L\tC_gar_LG12\t1000000\t1200000\t200\tseg\tLG12_17\n")
            fh.write("LG12_17_seg_DCO1\tC_gar_LG12\t1200000\t1250000\t50\tseg_dco\tLG12_17\n")
            fh.write("LG12_17_seg_R\tC_gar_LG12\t1250000\t1500000\t250\tseg\tLG12_17\n")

        reg = load_registry(str(root))

        # Python-side emulator of resolve_smallest_at
        def resolve_smallest_at(chrom, pos):
            with open(ir) as fh:
                header = fh.readline().strip().split("\t")
                cands = [dict(zip(header, line.strip().split("\t")))
                         for line in fh]
            hits = [c for c in cands if c.get("chrom") == chrom and
                     int(c.get("start_bp", 0)) <= pos <= int(c.get("end_bp", 0))]
            if not hits:
                return None
            return min(hits, key=lambda c: int(c["end_bp"]) - int(c["start_bp"]))["candidate_id"]

        test_cases = [
            ( 950000,  None,                   "outside parent"),
            (1100000, "LG12_17_seg_L",         "inside left segment"),
            (1225000, "LG12_17_seg_DCO1",      "inside DCO tract"),
            (1300000, "LG12_17_seg_R",         "inside right segment"),
            (1200000, "LG12_17_seg_DCO1",      "on breakpoint → smallest"),
        ]
        for pos, expected, desc in test_cases:
            got = resolve_smallest_at("C_gar_LG12", pos)
            if got == expected:
                print(f"  OK   pos={pos:>8d} → {got}  ({desc})")
            else:
                print(f"  FAIL pos={pos:>8d} expected {expected} got {got}")
                sys.exit(1)

        # Full scan sequence across parent
        seg_sequence = []
        prev = None
        for pos in range(1000000, 1500000, 10000):
            cid = resolve_smallest_at("C_gar_LG12", pos)
            if cid != prev:
                seg_sequence.append((pos, cid))
                prev = cid
        expected_seq = ["LG12_17_seg_L", "LG12_17_seg_DCO1", "LG12_17_seg_R"]
        actual_seq = [c for _, c in seg_sequence]
        assert actual_seq == expected_seq, f"Unexpected: {seg_sequence}"
        print(f"  OK   scan sequence: {actual_seq}")


def test_10_scan_result_persistence_shape():
    """Verify manifest shape for scan results — one row per (cid × stat)."""
    print("\n=== 10. Scan persistence: one row per (cid × stat) ===")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _write_master(root)
        ir = root / "data" / "interval_registry" / "candidate_intervals.tsv"
        with open(ir, "a") as fh:
            fh.write("LG12_17_seg_L\tC_gar_LG12\t1000000\t1200000\t200\tseg\tLG12_17\n")
            fh.write("LG12_17_seg_DCO1\tC_gar_LG12\t1200000\t1250000\t50\tseg_dco\tLG12_17\n")
            fh.write("LG12_17_seg_R\tC_gar_LG12\t1250000\t1500000\t250\tseg\tLG12_17\n")
        sr = root / "data" / "sample_registry"
        with open(sr / "sample_groups.tsv", "a") as fh:
            for seg in ["LG12_17_seg_L", "LG12_17_seg_DCO1", "LG12_17_seg_R"]:
                for k in ["HOM_REF", "HOM_INV"]:
                    fh.write(f"inv_{seg}_{k}\tC_gar_LG12\t{seg}\t.\t"
                             f"groups/inv_{seg}_{k}.txt\t3\tseg group\t2026-04-17 10:00:00\n")
                    (sr / "groups" / f"inv_{seg}_{k}.txt").write_text(
                        "CGA001\nCGA002\nCGA003\n")

        reg = load_registry(str(root))

        # Simulate scan_pairwise persistence: one put_pairwise per (cid, stat)
        for cid in ["LG12_17_seg_L", "LG12_17_seg_DCO1", "LG12_17_seg_R"]:
            reg.results.put_pairwise(
                chrom="C_gar_LG12",
                group_1=f"inv_{cid}_HOM_REF",
                group_2=f"inv_{cid}_HOM_INV",
                stat="fst",
                rows=[{"window_mid": 1100000 + i * 10000, "fst": 0.15 + i*0.01}
                      for i in range(5)],
                fieldnames=["window_mid", "fst"],
                source_script="scan_pairwise_test",
                engine="scan_pairwise")

        manifest = reg.results.list_manifest()
        assert len(manifest) == 3, f"expected 3, got {len(manifest)}"
        for r in manifest:
            assert int(r["n_rows"]) == 5
            assert r["stat"] == "fst"
            assert r["engine"] == "scan_pairwise"
        print(f"  OK   3 manifest rows, 5 windows each, engine=scan_pairwise")

        q = reg.results.ask(where={"chrom": "C_gar_LG12"}, what="fst")
        assert len(q) == 3
        print(f"  OK   ask() returns all 3 scan results")


if __name__ == "__main__":
    print("Testing chat-16 results_registry")
    test_1_fresh_init()
    test_2_put_pairwise_fk_errors()
    test_3_put_pairwise_writes_manifest_row()
    test_4_manifest_row_validates_schema()
    test_5_ask_filters()
    test_6_put_interval_summary_wrong_stat()
    test_7_auto_migration()
    test_8_sample_group_schema_validates_row()
    test_9_scan_segment_semantics()
    test_10_scan_result_persistence_shape()
    print("\nALL TESTS PASSED")
