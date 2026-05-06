#!/usr/bin/env python3
"""
test_step7_producers.py
-----------------------
End-to-end test for the step-7 producers (now covers all THREE per-candidate
layers — gt counts, evidence combinations, and the step-6 unlock support
matrix):

  1. Synthesize a tiny realistic input set (VCF + karyotype TSV +
     candidate JSON + evidence matrix + evidence_types JSON).
  2. Run STEP_SV_GT_AGG_aggregate_genotype_counts.py against it.
  3. Run STEP_SV_EVID_COMB_emit_combinations.py against it.
  4. Run STEP_SV_SUPPORT_emit_support_by_sample.py against it.
  5. Verify the per-candidate folder layout (all three layer files).
  6. Re-load the produced JSONs and check they pass the atlas-side
     validators (by spawning Node and calling the module's
     _validateCombinationsLayer + _validateSupportLayer + structural
     check on the gt layer).

Run:
  python3 tests/producers/test_step7_producers.py
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PRODUCERS = ROOT / "producers" / "sv_evidence"

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
# Synthesize inputs
# ---------------------------------------------------------------------------
section("Synthesizing inputs")

tmp = Path(tempfile.mkdtemp(prefix="sv_step7_"))
try:
    # 50 samples for a small realistic test (representative of how the
    # 226-sample cohort splits 25 / 16 / 9 across H1/H1, H1/H2, H2/H2):
    samples_h11 = [f"FL_{i:03d}" for i in range(0, 25)]
    samples_het = [f"FL_{i:03d}" for i in range(25, 41)]
    samples_h22 = [f"FL_{i:03d}" for i in range(41, 50)]
    all_samples = samples_h11 + samples_het + samples_h22

    # Karyotype TSV — page-21 lock format
    kary_path = tmp / "karyotype.tsv"
    with open(kary_path, "w") as f:
        f.write("sample_id\tlabel\n")
        for s in samples_h11: f.write(f"{s}\tHOMO_1\n")
        for s in samples_het: f.write(f"{s}\tHET\n")
        for s in samples_h22: f.write(f"{s}\tHOMO_2\n")

    # Candidate JSON
    cand_path = tmp / "candidate.json"
    cand = {
        "candidate_id":      "INV_LG28_test",
        "chrom":             "C_gar_LG28",
        "boundary_left_bp":  15_000_000,
        "boundary_right_bp": 18_000_000,
        "zone_definitions_bp": {
            "left_flank":     [14_000_000, 14_500_000],
            "left_boundary":  [14_500_000, 15_500_000],
            "inversion_body": [15_500_000, 17_500_000],
            "right_boundary": [17_500_000, 18_500_000],
            "right_flank":    [18_500_000, 19_000_000],
        },
    }
    cand_path.write_text(json.dumps(cand, indent=2))

    # SV calls — 4 distinct SVs across the candidate window
    # SV_A: at left boundary, strongly associated with H2/H2 carriers
    # SV_B: in body, associated with H1/H1 (opposite direction)
    # SV_C: at right boundary, het-specific
    # SV_D: in flank, uninformative noise
    sv_tsv = tmp / "calls.tsv"
    with open(sv_tsv, "w") as f:
        f.write("\t".join(["sv_id","chrom","position_bp","end_bp","sv_type",
                           "sample_id","GT","quality","callers"]) + "\n")
        # Helper to lay out per-sample GTs
        def emit(svid, pos, svtype, gtfn):
            for s in all_samples:
                gt = gtfn(s)
                f.write("\t".join([svid, "C_gar_LG28", str(pos), ".", svtype,
                                   s, gt, "PASS", "delly2"]) + "\n")
        # SV_A: BND at left boundary, almost all H2/H2 carry, almost no H1/H1 do
        emit("SV_A", 15_010_000, "BND",
             lambda s: ("1/1" if s in samples_h22 else
                        ("0/1" if s in samples_het else "0/0")))
        # SV_B: INV in body, H1/H1 carries, H2/H2 doesn't
        emit("SV_B", 16_500_000, "INV",
             lambda s: ("1/1" if s in samples_h11 else
                        ("0/1" if s in samples_het and int(s.split("_")[1]) % 2 == 0 else "0/0")))
        # SV_C: DEL at right boundary, het-specific (H1/H2 carry, both homs don't)
        emit("SV_C", 17_990_000, "DEL",
             lambda s: ("0/1" if s in samples_het else "0/0"))
        # SV_D: noise (uninformative)
        emit("SV_D", 18_700_000, "DUP",
             lambda s: ("0/1" if hash(s) % 5 == 0 else "0/0"))
    expect(sv_tsv.stat().st_size > 0, "calls TSV written")

    # Evidence matrix for the combinations producer
    ev_tsv = tmp / "evidence.tsv"
    ev_types = ["left_SA","right_SA","left_PE","right_PE",
                "Manta_INV_GT","DELLY_INV_GT","MAPQ0_left","MAPQ0_right"]
    with open(ev_tsv, "w") as f:
        f.write("sample_id\t" + "\t".join(ev_types) + "\n")
        # 4 distinct evidence patterns:
        #   - first 12 H2/H2 samples: full split-read+PE (the strongest pattern)
        #   - 5 of the H1/H1 samples: MAPQ0 only (the artefact pattern)
        #   - 4 het samples: caller-only
        #   - 2 het samples: left_SA only (single-sided)
        for s in all_samples:
            vec = ["0"] * len(ev_types)
            i = int(s.split("_")[1])
            if 41 <= i <= 49:                # H2/H2 — strongest
                vec[0] = vec[1] = vec[2] = vec[3] = "1"
            elif 0 <= i <= 4:                # 5 H1/H1 — MAPQ0 only
                vec[6] = vec[7] = "1"
            elif 25 <= i <= 28:              # 4 het — caller-only
                vec[4] = vec[5] = "1"
            elif 29 <= i <= 30:              # 2 het — left_SA only
                vec[0] = "1"
            f.write(s + "\t" + "\t".join(vec) + "\n")
    expect(ev_tsv.stat().st_size > 0, "evidence matrix TSV written")

    # Evidence types JSON
    et_path = tmp / "evidence_types.json"
    et_path.write_text(json.dumps([
        {"id":"left_SA",      "label":"Left split-read",  "side":"left",  "kind":"SA",     "tier":1},
        {"id":"right_SA",     "label":"Right split-read", "side":"right", "kind":"SA",     "tier":1},
        {"id":"left_PE",      "label":"Left PE",          "side":"left",  "kind":"PE",     "tier":2},
        {"id":"right_PE",     "label":"Right PE",         "side":"right", "kind":"PE",     "tier":2},
        {"id":"Manta_INV_GT", "label":"Manta INV",        "side":None,    "kind":"caller", "tier":3},
        {"id":"DELLY_INV_GT", "label":"DELLY INV",        "side":None,    "kind":"caller", "tier":3},
        {"id":"MAPQ0_left",   "label":"Left MAPQ0",       "side":"left",  "kind":"mapq0",  "tier":4},
        {"id":"MAPQ0_right",  "label":"Right MAPQ0",      "side":"right", "kind":"mapq0",  "tier":4},
    ], indent=2))

    out_root = tmp / "out"

    # ---------------------------------------------------------------------------
    section("Running STEP_SV_GT_AGG")
    r = subprocess.run([
        sys.executable,
        str(PRODUCERS / "STEP_SV_GT_AGG_aggregate_genotype_counts.py"),
        "--tsv", str(sv_tsv),
        "--candidate", str(cand_path),
        "--karyotype", str(kary_path),
        "--out-root", str(out_root),
        "--indent", "2",
    ], capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:\n" + r.stderr)
    expect(r.returncode == 0, "STEP_SV_GT_AGG exits 0")

    gt_path = out_root / "C_gar_LG28" / "candidates" / "INV_LG28_test" / "sv_genotype_counts.json"
    expect(gt_path.is_file(), f"gt counts file written at {gt_path}")
    gt = json.loads(gt_path.read_text())
    expect(gt["format_version"] == "sv_genotype_counts_v1", "format_version is correct")
    expect(gt["candidate_id"] == "INV_LG28_test", "candidate_id round-trips")
    expect(gt["boundary_left_bp"] == 15_000_000, "boundary_left preserved")
    expect(gt["boundary_right_bp"] == 18_000_000, "boundary_right preserved")
    expect(gt["groups_used"]["H1/H1"]["n"] == 25, "H1/H1 has 25 samples")
    expect(gt["groups_used"]["H1/H2"]["n"] == 16, "H1/H2 has 16 samples")
    expect(gt["groups_used"]["H2/H2"]["n"] == 9,  "H2/H2 has 9 samples")
    expect(len(gt["sv_calls"]) == 4, f"4 SV calls recorded (got {len(gt['sv_calls'])})")

    # SV_A should be strongly associated with H2/H2 → small FDR
    sv_a = next(s for s in gt["sv_calls"] if s["sv_id"] == "SV_A")
    expect(sv_a["fisher"]["fdr_bh"] is not None, "SV_A FDR computed")
    expect(sv_a["fisher"]["fdr_bh"] < 0.05, f"SV_A FDR < 0.05 (got {sv_a['fisher']['fdr_bh']})")
    expect(sv_a["zone"] == "left_boundary", f"SV_A zone classified (got {sv_a['zone']})")
    expect(sv_a["pattern_label"] in ("canonical_breakpoint_marker", "dominant_presence_marker"),
           f"SV_A is a strong marker (got {sv_a['pattern_label']})")

    # SV_C is het-specific
    sv_c = next(s for s in gt["sv_calls"] if s["sv_id"] == "SV_C")
    expect(sv_c["pattern_label"] == "het_specific_marker",
           f"SV_C classified as het_specific (got {sv_c['pattern_label']})")

    # boundary_summary should have non-zero counts
    expect(gt["boundary_summary"]["left"]["by_sv_type"]["BND"]["n_total"] >= 1,
           "boundary_summary left BND count >= 1")

    # ---------------------------------------------------------------------------
    section("Running STEP_SV_EVID_COMB")
    r = subprocess.run([
        sys.executable,
        str(PRODUCERS / "STEP_SV_EVID_COMB_emit_combinations.py"),
        "--evidence",        str(ev_tsv),
        "--evidence-types",  str(et_path),
        "--candidate-id",    "INV_LG28_test",
        "--chrom",           "C_gar_LG28",
        "--out-root",        str(out_root),
        "--top-n",           "20",
        "--indent",          "2",
    ], capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:\n" + r.stderr)
    expect(r.returncode == 0, "STEP_SV_EVID_COMB exits 0")

    cm_path = out_root / "C_gar_LG28" / "candidates" / "INV_LG28_test" / "sv_evidence_combinations.json"
    expect(cm_path.is_file(), f"combinations file written at {cm_path}")
    cm = json.loads(cm_path.read_text())
    expect(cm["format_version"] == "sv_evidence_combinations_v1", "format_version correct")
    expect(cm["candidate_id"] == "INV_LG28_test", "candidate_id round-trips")
    expect(len(cm["evidence_types"]) == 8, "8 evidence types preserved")

    # 4 distinct combinations from the synthesized data
    sizes = sorted([c["intersection_size"] for c in cm["combinations"]], reverse=True)
    expect(sizes[0] == 9,  f"top combination has 9 fish (the H2/H2 group, full evidence) (got {sizes[0]})")
    expect(sizes[1] == 5,  f"second combination has 5 fish (MAPQ0 only) (got {sizes[1]})")
    expect(sizes[2] == 4,  f"third combination has 4 fish (caller-only)")
    expect(sizes[3] == 2,  f"fourth combination has 2 fish (left_SA only)")

    # Check member sets
    bigm = next(c for c in cm["combinations"] if c["intersection_size"] == 9)
    expect(set(bigm["members"]) == {"left_SA","right_SA","left_PE","right_PE"},
           "top combination members are full split-read+PE")

    smallm = next(c for c in cm["combinations"] if c["intersection_size"] == 2)
    expect(set(smallm["members"]) == {"left_SA"},
           "smallest combination is left_SA only (single-sided)")

    # per_evidence_totals
    expect(cm["per_evidence_totals"]["left_SA"]["n_samples"] == 11,
           f"left_SA total = 11 (9 from full pattern + 2 from single-sided) (got {cm['per_evidence_totals']['left_SA']['n_samples']})")

    # ---------------------------------------------------------------------------
    section("Running STEP_SV_SUPPORT")
    r = subprocess.run([
        sys.executable,
        str(PRODUCERS / "STEP_SV_SUPPORT_emit_support_by_sample.py"),
        "--tsv",        str(sv_tsv),
        "--candidate",  str(cand_path),
        "--karyotype",  str(kary_path),
        "--out-root",   str(out_root),
        "--indent",     "2",
    ], capture_output=True, text=True)
    if r.returncode != 0:
        print("STDERR:\n" + r.stderr)
    expect(r.returncode == 0, "STEP_SV_SUPPORT exits 0")

    sp_path = out_root / "C_gar_LG28" / "candidates" / "INV_LG28_test" / "sv_support_by_sample.json"
    expect(sp_path.is_file(), f"support file written at {sp_path}")
    sp = json.loads(sp_path.read_text())
    expect(sp["format_version"] == "sv_support_by_sample_v1", "support format_version correct")
    expect(sp["candidate_id"] == "INV_LG28_test", "support candidate_id round-trips")
    expect(sp["encoding"] == "0=AA, 1=AB, 2=BB, .=miss", "support encoding string canonical")

    # 50 samples in the synthetic cohort (25 H1/H1 + 16 H1/H2 + 9 H2/H2)
    expect(len(sp["samples"]) == 50, f"50 samples in support layer (got {len(sp['samples'])})")
    expect(len(sp["dosage_compact"]) == 50, "dosage_compact has 50 rows")

    # 4 SVs in the candidate window: SV_A, SV_B, SV_C are inside the 14–19 Mb
    # window; SV_D is at 18.7 Mb which is in right_flank — so all 4 should be
    # included by the default 500 kb flank window.
    expect(len(sp["sv_ids"]) == 4, f"4 SVs in support layer (got {len(sp['sv_ids'])})")
    expect(set(sp["sv_ids"]) == {"SV_A", "SV_B", "SV_C", "SV_D"},
           "all 4 synthesized SVs included")

    # Each row's length == n_svs == 4
    expect(all(len(r) == 4 for r in sp["dosage_compact"]),
           "every row's compact string length matches n_svs=4")

    # row_groups should partition 0..49 cleanly: H1/H1 [0,24], H1/H2 [25,40], H2/H2 [41,49]
    expect(sp["row_groups"]["H1/H1"] == [0, 24], f"H1/H1 row range [0,24] (got {sp['row_groups']['H1/H1']})")
    expect(sp["row_groups"]["H1/H2"] == [25, 40], f"H1/H2 row range [25,40] (got {sp['row_groups']['H1/H2']})")
    expect(sp["row_groups"]["H2/H2"] == [41, 49], f"H2/H2 row range [41,49] (got {sp['row_groups']['H2/H2']})")

    # Sample ordering: first 25 are H1/H1
    expect(all(s.startswith("FL_") and int(s.split("_")[1]) < 25 for s in sp["samples"][:25]),
           "first 25 samples are H1/H1 (FL_000..FL_024)")
    expect(all(s.startswith("FL_") and 41 <= int(s.split("_")[1]) <= 49 for s in sp["samples"][41:50]),
           "samples 41..49 are H2/H2 (FL_041..FL_049)")

    # Dosage content sanity — SV_A is column index of "SV_A" in sv_ids; for H1/H1
    # samples the GT was 0/0 (AA -> '0'); for H2/H2 samples it was 1/1 (BB -> '2').
    sva_col = sp["sv_ids"].index("SV_A")
    h11_first = sp["dosage_compact"][0]    # FL_000
    h22_first = sp["dosage_compact"][41]   # FL_041
    expect(h11_first[sva_col] == "0", f"H1/H1 first row, SV_A col = '0' (AA) — got {h11_first[sva_col]!r}")
    expect(h22_first[sva_col] == "2", f"H2/H2 first row, SV_A col = '2' (BB) — got {h22_first[sva_col]!r}")

    # SV_C: het-specific in the synthesized data — H1/H2 samples carry 0/1 (AB -> '1')
    svc_col = sp["sv_ids"].index("SV_C")
    h12_first = sp["dosage_compact"][25]   # FL_025 (first H1/H2)
    expect(h12_first[svc_col] == "1", f"H1/H2 first row, SV_C col = '1' (AB het-specific) — got {h12_first[svc_col]!r}")
    expect(h11_first[svc_col] == "0", f"H1/H1 first row, SV_C col = '0' (AA non-carrier) — got {h11_first[svc_col]!r}")

    # ---------------------------------------------------------------------------
    section("Manifest aggregation")
    manifest_path = out_root / "C_gar_LG28" / "manifest.json"
    expect(manifest_path.is_file(), "chrom-level manifest.json written")
    mf = json.loads(manifest_path.read_text())
    expect(mf["format_version"] == "sv_evidence_chrom_manifest_v1",
           "manifest format_version correct")
    expect(any(c["candidate_id"] == "INV_LG28_test" for c in mf["candidates"]),
           "manifest lists INV_LG28_test")
    cand_entry = next(c for c in mf["candidates"] if c["candidate_id"] == "INV_LG28_test")
    expect("sv_genotype_counts_v1" in cand_entry["layers"], "manifest knows gt layer")
    expect("sv_evidence_combinations_v1" in cand_entry["layers"], "manifest knows combinations layer")
    expect("sv_support_by_sample_v1" in cand_entry["layers"], "manifest knows support layer")

    # ---------------------------------------------------------------------------
    section("Atlas-side validator round-trip")
    # Spawn Node and run the atlas validators on all three files
    node_script = f"""
const path = require('path');
process.chdir('{ROOT}');
const mod = require('./js/atlas_sv_evidence.js');
const fs = require('fs');
const gt = JSON.parse(fs.readFileSync('{gt_path}', 'utf-8'));
const cm = JSON.parse(fs.readFileSync('{cm_path}', 'utf-8'));
const sp = JSON.parse(fs.readFileSync('{sp_path}', 'utf-8'));
const v1 = mod._internals._validateCombinationsLayer(cm);
if (!v1.ok) {{ console.error('combinations validation failed:', v1.reason); process.exit(2); }}
const v2 = mod._internals._validateSupportLayer(sp);
if (!v2.ok) {{ console.error('support validation failed:', v2.reason); process.exit(6); }}
// Smoke check on gt: format_version, sv_calls is array, fisher object present
if (gt.format_version !== 'sv_genotype_counts_v1') {{ process.exit(3); }}
if (!Array.isArray(gt.sv_calls)) {{ process.exit(4); }}
if (!gt.sv_calls.every(s => s.fisher && typeof s.fisher === 'object')) {{ process.exit(5); }}
console.log('OK');
"""
    r = subprocess.run(["node", "-e", node_script], capture_output=True, text=True)
    if r.returncode != 0:
        print("STDOUT:", r.stdout)
        print("STDERR:", r.stderr)
    expect(r.returncode == 0, "all three layers validate against atlas-side validators")

finally:
    # Always clean up
    if tmp.exists():
        shutil.rmtree(tmp)

# ---------------------------------------------------------------------------
print(f"\n=========================================")
print(f"  PASS: {passed}   FAIL: {failed}")
print(f"=========================================")
sys.exit(0 if failed == 0 else 1)
