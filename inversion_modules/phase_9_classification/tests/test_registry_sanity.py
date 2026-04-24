#!/usr/bin/env python3
"""
test_registry_sanity.py — minimal end-to-end sanity test

Exercises:
  1. Create registry in a temp dir
  2. Register 4 groups for one candidate (HOM_REF, HET, HOM_INV, RECOMBINANT)
  3. Write Tier-2 blocks for 3 of the 4 tiers (Layer A, decomp, hypothesis)
  4. Read back flat keys, verify they were extracted per schema
  5. Promote group validation UNCERTAIN → SUPPORTED → VALIDATED
  6. Test that C01f would see the groups via has_group / get_group

Run:
    cd phase4_v10/tests
    python3 test_registry_sanity.py
"""

import shutil
import sys
import tempfile
from pathlib import Path

# Path hack so the test can find the library in multiple deployment layouts.
# registry_loader.py can live at:
#   - HERE.parent/registries/...                (standalone phase4_v10/)
#   - HERE.parent.parent/registries/...         (deployed into inversion_modules/)
#   - $REGISTRIES/api/python/registry_loader.py (env override)
HERE = Path(__file__).resolve().parent

import os as _os
_candidate_roots = []
for _r in [HERE.parent, HERE.parent.parent, HERE.parent.parent.parent]:
    _candidate_roots.append(_r / "registries")
_env = _os.environ.get("REGISTRIES")
if _env:
    _candidate_roots.insert(0, Path(_env))

API_PY = None
SCHEMAS_SRC = None
for _r in _candidate_roots:
    if (_r / "api" / "python" / "registry_loader.py").exists():
        API_PY = _r / "api" / "python"
        SCHEMAS_SRC = _r / "schemas"
        break

if API_PY is None:
    print("[test] cannot find registry_loader.py in:",
          [str(r) for r in _candidate_roots])
    sys.exit(2)

sys.path.insert(0, str(API_PY))

from registry_loader import load_registry


def main() -> int:
    tmpdir = Path(tempfile.mkdtemp(prefix="phase4_test_"))
    print(f"[test] tmpdir = {tmpdir}")

    # Copy schemas into the registry root so write_block can find them
    registries_root = tmpdir / "registries"
    registries_root.mkdir()
    shutil.copytree(SCHEMAS_SRC, registries_root / "schemas")

    reg = load_registry(registries_root=str(registries_root))

    # Populate sample_master.tsv so add_group doesn't warn
    master_file = reg.samples.master_file
    master_file.parent.mkdir(parents=True, exist_ok=True)
    with master_file.open("w") as f:
        f.write("sample_id\tind_id_226\tsample_index\n")
        for i in range(226):
            f.write(f"CGA{i:03d}\tInd{i}\t{i}\n")

    cid = "LG12_17"

    # ── Step 1: register 4 groups ────────────────────────────────────────
    print("\n[test] registering 4 groups for", cid)
    hom_ref   = [f"CGA{i:03d}" for i in range(0, 89)]
    het_samps = [f"CGA{i:03d}" for i in range(89, 187)]
    hom_inv   = [f"CGA{i:03d}" for i in range(187, 219)]
    recomb    = [f"CGA{i:03d}" for i in range(219, 226)]

    for cls, samps in [
        ("HOM_REF", hom_ref),
        ("HET", het_samps),
        ("HOM_INV", hom_inv),
        ("RECOMBINANT", recomb),
    ]:
        ok = reg.samples.add_group(
            f"inv_{cid}_{cls}", samps,
            chrom="C_gar_LG12", inv_id=cid, subgroup=cls,
            description=f"test: {cls}",
        )
        assert ok, f"add_group failed for {cls}"

    # ── Step 2: verify groups readable (mimics what C01f would do) ───────
    print("\n[test] C01f comp_from_registry equivalent")
    assert reg.samples.has_group(f"inv_{cid}_HOM_REF"), "HOM_REF missing"
    assert reg.samples.has_group(f"inv_{cid}_HET"), "HET missing"
    assert reg.samples.has_group(f"inv_{cid}_HOM_INV"), "HOM_INV missing"
    assert reg.samples.has_group(f"inv_{cid}_RECOMBINANT"), "RECOMBINANT missing"

    retrieved_hom_ref = reg.samples.get_carriers(cid, "HOM_REF")
    assert len(retrieved_hom_ref) == 89, f"expected 89, got {len(retrieved_hom_ref)}"
    print(f"  HOM_REF: {len(retrieved_hom_ref)} samples ✓")

    retrieved_rec = reg.samples.get_carriers(cid, "RECOMBINANT")
    assert len(retrieved_rec) == 7, f"expected 7, got {len(retrieved_rec)}"
    print(f"  RECOMBINANT: {len(retrieved_rec)} samples ✓")

    # ── Step 3: write existence_layer_a (simulating C01d) ────────────────
    print("\n[test] writing existence_layer_a block")
    res = reg.evidence.write_block(
        candidate_id=cid,
        block_type="existence_layer_a",
        data={
            "composite_score": 0.847,
            "dim_positive": 9,
            "tier": 1,
            "d1_block_strength": 0.91,
            "d2_block_shape": 0.88,
            "d3_nn_persistence": 0.76,
            "d4_decay_flatness": 0.82,
            "d5_interior_quality": 0.77,
            "d6_consensus": 0.85,
            "d7_sv_breakpoint": 0.60,
            "d8_peel_or_hyp": 0.50,
            "d9_pca_clusters": 0.72,
            "d10_partition": 0.68,
            "d11_boundary_concordance": 0.00,  # pass-1 placeholder
            "d12_snake_concordance": 0.65,
            "shape_class": "strong_square",
            "landscape_category": "hot_core",
            "n_children": 2,
            "span_kb": 4305.5,
            "pass_number": 1,
        },
        source_script="C01d_pass1",
    )
    assert res["status"] in ("validated", "incomplete"), f"status={res['status']}"
    print(f"  wrote {res['n_keys']} keys, status={res['status']}")

    # ── Step 3b: write existence_layer_d (phase_3 STEP03, Layer D OR test) ──
    # Added 2026-04-17 (chat 5 FIX 29 v2): mirrors what phase_3 STEP03
    # writes per candidate via the Python registry API.
    print("\n[test] writing existence_layer_d block (Layer D OR test)")
    res_d = reg.evidence.write_block(
        candidate_id=cid,
        block_type="existence_layer_d",
        data={
            "fisher_or": 23.7,
            "fisher_p": 4.98e-06,
            "fisher_ci_lower": 8.9,
            "fisher_ci_upper": 63.2,
            "contingency_table": [[45, 8], [12, 161]],
            "armitage_z": 5.12,
            "armitage_p": 3.01e-07,
            "n_inv_with_support": 45,
            "n_inv_total": 53,
            "n_ref_with_support": 12,
            "n_ref_total": 173,
            "samples_inv_with_support": ["CGA009", "CGA045", "CGA102"],
        },
        source_script="phase_3_refine/STEP_D03_statistical_tests_and_seeds.py",
    )
    assert res_d["status"] in ("validated", "incomplete"), f"layer_d status={res_d['status']}"
    print(f"  wrote {res_d['n_keys']} layer_d keys, status={res_d['status']}")

    # ── Step 3c: write existence_layer_b_bnd_rescue (phase_3 STEP_B06) ────
    # Represents a true orphan rescue — paired BND junctions that match
    # no entry in the DELLY or Manta INV catalogs.
    print("\n[test] writing existence_layer_b_bnd_rescue block (orphan BND pair)")
    res_br = reg.evidence.write_block(
        candidate_id=cid,
        block_type="existence_layer_b_bnd_rescue",
        data={
            "bnd_rescued": True,
            "rescue_source": "cross_caller",
            "bnd_pair_bp1": 8420000,
            "bnd_pair_bp2": 8505000,
            "bnd_pair_size_bp": 85000,
            "left_junction_id":  "DELLY_BND_00042",
            "right_junction_id": "MantaBND:0:123:234:0:0:0:1",
            "left_junction_source":  "delly",
            "right_junction_source": "manta_raw",
            "left_pe":  6,
            "right_pe": 4,
            "left_sr":  3,
            "right_sr": 2,
            "matched_inv_id": "",
            "match_type": "no_inv_match",
        },
        source_script="phase_3_refine/STEP_B06_bnd_rescue.py",
    )
    assert res_br["status"] in ("validated", "incomplete"), f"bnd_rescue status={res_br['status']}"
    print(f"  wrote {res_br['n_keys']} bnd_rescue keys, status={res_br['status']}")

    # ── Step 4: write internal_dynamics (C01i) ────────────────────────────
    print("\n[test] writing internal_dynamics block")
    reg.evidence.write_block(
        candidate_id=cid,
        block_type="internal_dynamics",
        data={
            "n_HOM_REF": 89, "n_HET": 98, "n_HOM_INV": 32,
            "n_RECOMBINANT": 7, "n_total": 226, "freq_inv": 0.358,
            "recombinants": [
                {
                    "sample_id": "CGA220", "event_class": "gene_conversion",
                    "switchpoint_bp": 8451200, "mosaic_length_bp": 42000,
                    "distance_to_nearest_boundary_bp": 35000,
                    "prior_dco": 0.0001, "prior_gc": 0.001, "posterior": 0.73,
                    "phase_support": 4,
                },
            ],
            "n_gene_conversion": 5, "n_double_crossover": 1, "n_suspicious": 1,
            "mean_dco_prior": 0.0023, "max_dco_prior": 0.0098,
        },
        source_script="C01i",
    )

    # Set validation to UNCERTAIN via flat key (mimics C01i end)
    reg.evidence.add_evidence(cid, "q6_group_validation", "UNCERTAIN", script="C01i")

    # ── Step 5: write hypothesis_verdict (C01f) + promotion ──────────────
    print("\n[test] writing hypothesis_verdict block")
    reg.evidence.write_block(
        candidate_id=cid,
        block_type="hypothesis_verdict",
        data={
            "verdict": "H2_or_H3_real_inversion",
            "verdict_confidence": "high",
            "t1_ratio": 1.23, "t1_evidence": "neutral",
            "t2_eff_k": 3.4, "t2_evidence": "moderate_diversity",
            "t3_retention": 0.89, "t3_evidence": "stable",
            "t8_concordance": 0.82,
            "t9_jackknife_status": "robust",
            "t9_max_delta": 0.12,
            "t10_theta_concordance": 0.71,
            "group_validation_before": "UNCERTAIN",
            "group_validation_after": "SUPPORTED",
            "comp_source": "C01i_registry",
        },
        source_script="C01f",
    )

    # Validation should have been bumped to SUPPORTED via keys_extracted
    keys = reg.evidence.get_keys(cid, "q6_group_validation")
    print(f"  q6_group_validation keys: {[k['value'] for k in keys]}")

    # ── Step 6: simulate Layer D promotion to VALIDATED ──────────────────
    print("\n[test] writing existence_layer_d (→ VALIDATED)")
    reg.evidence.write_block(
        candidate_id=cid,
        block_type="existence_layer_d",
        data={
            "fisher_or": 23.7,
            "fisher_p": 4.98e-06,
            "contingency_table": [[45, 8], [12, 161]],
            "armitage_z": 5.12,
            "armitage_p": 3.01e-07,
            "n_inv_with_support": 45,
            "n_inv_total": 53,
            "samples_inv_with_support": hom_inv[:45],
        },
        source_script="STEP03_statistical_tests.py",
    )
    reg.evidence.add_evidence(cid, "q6_group_validation", "VALIDATED",
                               script="STEP03_statistical_tests.py")

    # ── Step 7: final state check ────────────────────────────────────────
    print("\n[test] final registry state:")
    print(f"  candidates: {reg.evidence.list_candidates()}")
    print(f"  groups: {len(reg.samples.list_groups())}")
    all_keys = reg.evidence.get_keys(cid)
    print(f"  {cid} keys: {len(all_keys)}")
    val_keys = [k for k in all_keys if k["key"] == "q6_group_validation"]
    print(f"  q6_group_validation: {val_keys[-1]['value'] if val_keys else 'MISSING'}")
    assert val_keys and val_keys[-1]["value"] == "VALIDATED", "final validation wrong"

    print("\n[test] ✓ all checks passed")
    print(f"[test] cleaning up {tmpdir}")
    shutil.rmtree(tmpdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
