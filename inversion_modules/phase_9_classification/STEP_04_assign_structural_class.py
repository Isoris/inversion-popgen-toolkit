#!/usr/bin/env python3
"""
STEP_04_assign_structural_class.py — final structural-class synthesis

=============================================================================
PIPELINE POSITION
=============================================================================
  STEP 1  assembled-junction forensics       → mechanism_assembled.json
  STEP 2  single-sided BND support           → bnd_sided_support.json
  STEP 3A cross-species synteny bridge       → synteny_v6.json
  STEP 3B bp-pipeline bridge                 → fragment_distribution.json
→ STEP 4  structural-class assignment        ← THIS SCRIPT
                                              ↓
                                             final_label.json + final_catalog.tsv

Renamed from assign_structural_class_v7.py in pass 18 (2026-04-24). The
"v7" version tag now lives in the git history; the filename matches the
step number. v6→v7 logic changes (reading fragment_distribution instead
of interior_structure, adding complex_rearrangement_out_of_scope
terminal gate, adding interior-suffix modifiers) are preserved — only
the filename changed.

=============================================================================
ROLE
=============================================================================
Reads the four Tier-2 evidence blocks produced by STEPs 1/2/3A/3B plus
the per-candidate flat keys.tsv, and decides one of 16 structural class
labels per candidate. Writes per-candidate final_label.json and cohort
final_catalog.tsv.

This is the ONLY STEP that writes the terminal "q_overall_structural_class"
key. Everything upstream produces evidence; this step weighs it and
commits to a label.

=============================================================================
INPUTS
=============================================================================
--candidates           TSV with candidate_id column (the cohort list)
--keys_dir             dir containing per-candidate keys.tsv (flat keys
                        from registry: q1_*, q2_*, q3_*, q4_*, q5_*, q6_*,
                        q7_*). Expected layout:
                            <keys_dir>/<cid>/keys.tsv
                          or
                            <keys_dir>/<cid>.keys.tsv
--evidence_blocks_dir  dir containing per-candidate structured/ folders.
                        This is the same as STEP 1/2/3A/3B's --outdir.
                        Expected layout:
                            <evidence_blocks_dir>/<cid>/structured/
                                mechanism_assembled.json    ← STEP 1
                                bnd_sided_support.json      ← STEP 2
                                synteny_v6.json             ← STEP 3A
                                    (or synteny_dollo.json as legacy fallback)
                                fragment_distribution.json  ← STEP 3B
                                    (or interior_structure.json legacy)
--outdir               where to write final_label.json per candidate +
                        final_catalog.tsv cohort file

=============================================================================
DATA FLOW (how each decision input reaches this script)
=============================================================================

  FROM STEP 1 (mechanism_assembled.json):
    q4b_asm_precise_record_available → has_asm
    q4b_asm_junction_class           → asm_jc     (BLUNT/MH_short/MH_long/...)
    q4b_asm_source                   → asm_src    (DELLY or Manta)
      ↓ used by mechanism_hypothesis() to upgrade mechanism_hypothesis
      ↓ label with "_supported_by_assembly" vs "_hypothesis" suffix

  FROM STEP 2 (bnd_sided_support.json):
    q7b_bnd_left_support   → (currently loaded but not used in v7 rules;
    q7b_bnd_right_support     reserved for future tier-upgrade logic)

  FROM STEP 3A (synteny_v6.json):
    q5_conservation_class            → conservation
      ↓ used to fire supported_shared_between_species_inversion label
    q5_tree_polarization_direction   → justification text

  FROM STEP 3B (fragment_distribution.json):
    q2_interior_class                → interior_cls
      ↓ used by interior_suffix() to decide:
         complex_rearrangement_out_of_scope       → TERMINAL GATE
                                                    (overrides everything)
         deep_interior_recombinants               → _with_deep_recombinants
         bimodal_boundary_signal                  → _with_partial_interior_recombination
         clean_simple / edge_recombinants_only    → no suffix
         degenerate                               → no suffix (means
                                                    STEP 3B didn't have
                                                    phase_6 TSV inputs)

  FROM flat keys.tsv (written by all upstream phase writers):
    q7_layer_a_detected, q7_layer_b_detected,
    q7_layer_c_ghsl_detected, q7_layer_c_ghsl_quality,
    q7_layer_d_fisher_p, q7_layer_d_fisher_or
      ↓ existence_verdict() → (class, n_layers)

    q6_family_linkage, q1_family_likeness_mean,
    q7_t9_jackknife_status
      ↓ family_class() → clean / single_family_real /
                          family_confounded_fake / ambiguous

    q4_has_inverted_sd, q4_sd_identity_pct, q4_sd_length,
    q4_te_enrichment, q4_tr_enrichment
      ↓ substrate_verdict() → inverted_SD_strong / inverted_SD_weak /
                               TE_flanked / TR_rich / none / unknown

    q1_composite_flag, q1_n_children
      ↓ complex_unresolved_locus or supported_nested_inversion

    q6_freq_inv (for supported_fixed_inversion_SV_only decision)
    q7_verdict  (contradiction / artifact → terminal labels)

=============================================================================
CHANGES FROM v6 (preserved from old docstring)
=============================================================================
1. Reads fragment_distribution.json (from STEP_03B_bp_pipeline_bridge.py)
   and uses q2_interior_class instead of v6's interior_structure.json.
2. Adds terminal label `complex_rearrangement_out_of_scope` — fires
   when q2_interior_class = complex_rearrangement_out_of_scope. Overrides
   everything else; candidate excluded from mechanism/age/conservation.
3. Replaces `_with_interior_complexity` suffix with a graduated scale:
     edge_recombinants_only        → no suffix modifier (normal)
     bimodal_boundary_signal       → _with_partial_interior_recombination
     deep_interior_recombinants    → _with_deep_recombinants
     complex_rearrangement_out_of_scope → terminal label
4. Honors D.U.P-T.R.P/INV-DUP gating: those events land in
   complex_rearrangement_out_of_scope automatically because the fragment
   distribution won't look like a clean balanced inversion.

=============================================================================
LABELS (16 total, alphabetical)
=============================================================================
  candidate_weak_evidence
  complex_rearrangement_out_of_scope                   ← TERMINAL GATE
  complex_unresolved_locus
  contradicted_candidate
  diffuse_region_of_interest
  family_confounded_locus
  single_family_polymorphism
  supported_balanced_inversion_NAHR_like_supported_by_assembly
  supported_balanced_inversion_NAHR_like_hypothesis
  supported_balanced_inversion_NHEJ_like_supported_by_assembly
  supported_balanced_inversion_NHEJ_like_hypothesis
  supported_balanced_inversion_simple
  supported_balanced_inversion_with_substrate_mechanism_unresolved
  supported_balanced_inversion_with_deep_recombinants          ← NEW in v7
  supported_balanced_inversion_with_partial_interior_recombination ← NEW in v7
  supported_fixed_inversion_SV_only
  supported_locus_unreconstructed
  supported_nested_inversion
  supported_shared_between_species_inversion

Balanced-inversion labels can optionally receive an interior suffix from
interior_suffix() when q2_interior_class warrants it. The combinatorial
result is larger than 16 but the base labels above are what the
assign_label() function returns before interior-suffix decoration.

=============================================================================
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def load_keys_tsv(path):
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            k = row.get("key") or row.get("Key")
            v = row.get("value") or row.get("Value")
            if k:
                out[k] = v
    return out


def load_block(cand_dir, block_name):
    p = Path(cand_dir) / "structured" / f"{block_name}.json"
    if not p.exists():
        return {}
    try:
        j = json.loads(p.read_text())
        return j.get("data", {})
    except Exception:
        return {}


def to_bool(v, default=False):
    if v is None:
        return default
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    return s in ("true", "1", "yes")


def to_float(v, default=None):
    if v is None or v == "":
        return default
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


# -----------------------------------------------------------------------------
# Evidence aggregators (copy from v6)
# -----------------------------------------------------------------------------
def existence_verdict(keys):
    layers = 0
    if to_bool(keys.get("q7_layer_a_detected")): layers += 1
    if to_bool(keys.get("q7_layer_b_detected")): layers += 1
    if to_bool(keys.get("q7_layer_c_ghsl_detected")):
        q = (keys.get("q7_layer_c_ghsl_quality") or "").upper()
        if q in ("HIGH", "MODERATE"): layers += 1
    p = to_float(keys.get("q7_layer_d_fisher_p"), 1.0)
    or_val = to_float(keys.get("q7_layer_d_fisher_or"), 0.0)
    if p is not None and p < 0.05 and or_val > 5: layers += 1

    verdict = (keys.get("q7_verdict") or "").lower()
    if "artifact" in verdict or "family_structure" in verdict:
        return "artifact", layers
    if layers >= 3: return "confirmed", layers
    if layers == 2: return "likely", layers
    if layers == 1: return "candidate", layers
    return "insufficient", layers


def family_class(keys):
    linkage = (keys.get("q6_family_linkage") or "").lower()
    fk_lik = to_float(keys.get("q1_family_likeness_mean"), 0.0)
    jk = (keys.get("q7_t9_jackknife_status") or "").lower()
    if linkage == "pca_family_confounded":
        return "family_confounded_fake"
    if linkage == "single_family" or jk == "single_family_fragile":
        return "single_family_real"
    if linkage == "multi_family":
        return "clean"
    if fk_lik > 0.5:
        return "family_confounded_fake"
    return "ambiguous"


def substrate_verdict(keys):
    if to_bool(keys.get("q4_has_inverted_sd")):
        ident = to_float(keys.get("q4_sd_identity_pct"), 0)
        length = to_float(keys.get("q4_sd_length"), 0)
        if ident and length and ident > 90 and length > 1000:
            return "inverted_SD_strong"
        if ident and length and ident > 85 and length > 500:
            return "inverted_SD_weak"
    te_en = (keys.get("q4_te_enrichment") or "").upper()
    if te_en == "ENRICHED":
        return "TE_flanked"
    tr_en = (keys.get("q4_tr_enrichment") or "").upper()
    if tr_en == "ENRICHED":
        return "TR_rich"
    if keys.get("q4_has_inverted_sd") in ("FALSE", "False", False):
        return "none"
    return "unknown"


def assembled_junction_class(asm_data):
    if not asm_data:
        return False, None, None
    return (
        to_bool(asm_data.get("q4b_asm_precise_record_available")),
        asm_data.get("q4b_asm_junction_class"),
        asm_data.get("q4b_asm_source"),
    )


def mechanism_hypothesis(substrate, asm_junction_class, has_asm):
    if not substrate or substrate == "unknown":
        return None, None
    jc = asm_junction_class
    if substrate in ("inverted_SD_strong", "inverted_SD_weak"):
        if has_asm and jc == "MH_long":
            return "NAHR_like", "_supported_by_assembly"
        return "NAHR_like", "_hypothesis"
    if substrate == "none":
        if has_asm and jc == "blunt":
            return "NHEJ_like", "_supported_by_assembly"
        if has_asm and jc == "MH_short":
            return "MMEJ_like", "_supported_by_assembly"
        return "NHEJ_like", "_hypothesis"
    if substrate == "TE_flanked":
        return "TE_mediated_like", "_hypothesis"
    if substrate == "TR_rich":
        return "MMBIR_like", "_hypothesis"
    return None, None


def interior_suffix(interior_class):
    """
    Map q2_interior_class → suffix modifier for balanced-inversion labels.
    Returns (suffix, is_terminal_gate).
    """
    if interior_class == "complex_rearrangement_out_of_scope":
        return None, True
    if interior_class == "deep_interior_recombinants":
        return "_with_deep_recombinants", False
    if interior_class == "bimodal_boundary_signal":
        return "_with_partial_interior_recombination", False
    # clean_simple, edge_recombinants_only, degenerate → no suffix
    return "", False


# -----------------------------------------------------------------------------
# Master assignment (v7)
# -----------------------------------------------------------------------------
def assign_label(keys, asm_data, bnd_sided, synteny_data, fragment_data):
    """
    Returns (label, justifications[], weakest_component).
    """
    just = []
    ex_class, n_layers = existence_verdict(keys)
    fam = family_class(keys)
    subst = substrate_verdict(keys)
    has_asm, asm_jc, asm_src = assembled_junction_class(asm_data)

    # Read interior class from fragment_distribution block
    interior_cls = (fragment_data or {}).get("q2_interior_class", "")
    int_suffix, is_gate = interior_suffix(interior_cls)

    # -------- terminal labels first --------

    # 0. HARD GATE: complex rearrangement out of scope — overrides
    # everything. This is the DUP-TRP/INV-DUP-like case, or any
    # candidate where >70% of carriers are recombinant: cannot be
    # interpreted as a simple inversion with current data.
    if is_gate:
        return "complex_rearrangement_out_of_scope", [
            f"q2_interior_class={interior_cls}",
            f"q2_fragment_interior_recomb_fraction="
            f"{(fragment_data or {}).get('q2_fragment_interior_recomb_fraction')}",
            "excluded from mechanism/age/conservation analyses"
        ], "complex_interior"

    # 1. Contradicted
    verdict = (keys.get("q7_verdict") or "").lower()
    if "contradict" in verdict or "inconsistent" in verdict:
        return "contradicted_candidate", \
               ["q7_verdict contains contradiction"], "contradiction"

    # 2. Artifact / family fake
    if ex_class == "artifact":
        return "family_confounded_locus", [
            f"q7_verdict={verdict}", f"family_class={fam}"
        ], "existence_layer"

    if fam == "family_confounded_fake":
        return "family_confounded_locus", [
            "q6_family_linkage=pca_family_confounded",
        ], "family_confounding"

    # 3. Single-family real
    if fam == "single_family_real" and ex_class in ("confirmed", "likely"):
        return "single_family_polymorphism", [
            f"existence={ex_class} ({n_layers} layers)",
        ], "population_distribution"

    # 4. Q1 complex architecture (likely_composite)
    composite = (keys.get("q1_composite_flag") or "").lower()
    if composite == "likely_composite":
        return "complex_unresolved_locus", [
            f"q1_composite_flag={composite}",
        ], "q1_architecture"

    # 5. Insufficient
    if ex_class == "insufficient":
        return "candidate_weak_evidence", [
            f"only {n_layers}/4 existence layers pass"
        ], "existence_layer"

    # 6. Fixed SV-only
    if (not to_bool(keys.get("q7_layer_a_detected"))
            and to_bool(keys.get("q7_layer_b_detected"))):
        freq = to_float(keys.get("q6_freq_inv"), 0.0)
        if freq and freq > 0.85:
            return "supported_fixed_inversion_SV_only", [
                f"q6_freq_inv={freq:.2f}",
                "Layer_A(PCA) absent, Layer_B(SV) present"
            ], "PCA_contrast_missing"

    # 7. Shared between species
    conservation = (synteny_data or {}).get("q5_conservation_class", "")
    if conservation in ("shared_ancestral", "ASEAN_diagnostic") \
            and ex_class in ("confirmed", "likely"):
        lbl = "supported_shared_between_species_inversion"
        # Apply interior suffix if present
        if int_suffix:
            lbl = lbl + int_suffix
        return lbl, [
            f"q5_conservation_class={conservation}",
            f"q5_tree_polarization="
            f"{(synteny_data or {}).get('q5_tree_polarization_direction')}"
        ], "cross_species"

    # 8. Nested (Q1)
    n_children = to_float(keys.get("q1_n_children"), 0.0)
    if n_children == 1 and ex_class in ("confirmed", "likely"):
        return "supported_nested_inversion", [
            f"q1_n_children={keys.get('q1_n_children')}",
        ], "q1_architecture"

    # 9. Supported + mechanism paths (with interior-suffix decoration)
    mech_lbl, mech_suffix = mechanism_hypothesis(subst, asm_jc, has_asm)
    if ex_class in ("confirmed", "likely"):
        if mech_lbl and mech_suffix:
            lbl = f"supported_balanced_inversion_{mech_lbl}{mech_suffix}"
            if int_suffix:
                lbl = lbl + int_suffix
            just.append(f"substrate={subst}")
            if has_asm:
                just.append(f"assembled_junction_class={asm_jc} (source={asm_src})")
            if int_suffix:
                just.append(f"q2_interior_class={interior_cls}")
            weakest = "assembly_unavailable" if not has_asm \
                else ("interior" if int_suffix else "mechanism_hypothesis_only")
            return lbl, just, weakest

        if subst not in (None, "unknown", "none"):
            lbl = "supported_balanced_inversion_with_substrate_mechanism_unresolved"
            if int_suffix:
                lbl = lbl + int_suffix
            return lbl, [
                f"substrate={subst}",
                "no junction synthesis possible",
            ] + ([f"q2_interior_class={interior_cls}"] if int_suffix else []), \
                "junction_data"

        # Simple baseline
        lbl = "supported_balanced_inversion_simple"
        if int_suffix:
            lbl = lbl + int_suffix
        return lbl, [
            f"existence={ex_class} ({n_layers} layers)",
            "no substrate evidence, no mechanism claim",
        ] + ([f"q2_interior_class={interior_cls}"] if int_suffix else []), \
            "mechanism_unresolved"

    # 10. Candidate
    if ex_class == "candidate":
        return "candidate_weak_evidence", [
            f"only {n_layers}/4 existence layers pass"
        ], "existence_layer"

    # 11. Fallthrough
    return "supported_locus_unreconstructed", [
        f"existence={ex_class}"
    ], "fallthrough"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("=" * 77)[0])
    p.add_argument("--candidates", required=True)
    p.add_argument("--keys_dir", required=True)
    p.add_argument("--evidence_blocks_dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--dry_run", action="store_true")
    args = p.parse_args()

    cand_ids = []
    with open(args.candidates) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("candidate_id"):
                cand_ids.append(row["candidate_id"])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    keys_dir = Path(args.keys_dir)
    evidence_dir = Path(args.evidence_blocks_dir)

    summary = defaultdict(int)
    rows_out = []

    for cid in cand_ids:
        per_cand_keys = keys_dir / cid / "keys.tsv"
        if not per_cand_keys.exists():
            per_cand_keys = keys_dir / f"{cid}.keys.tsv"
        keys = load_keys_tsv(str(per_cand_keys)) if per_cand_keys.exists() else {}

        cand_evidence = evidence_dir / cid
        asm = load_block(cand_evidence, "mechanism_assembled")
        bnd = load_block(cand_evidence, "bnd_sided_support")
        syn = load_block(cand_evidence, "synteny_v6") or \
              load_block(cand_evidence, "synteny_dollo")
        fragment = load_block(cand_evidence, "fragment_distribution") or \
                   load_block(cand_evidence, "interior_structure")

        label, just, weakest = assign_label(keys, asm, bnd, syn, fragment)
        summary[label] += 1
        rows_out.append({
            "candidate_id": cid,
            "q_overall_structural_class": label,
            "weakest_component": weakest,
            "justification": "; ".join(just),
        })

        if not args.dry_run:
            cand_out = outdir / cid
            cand_out.mkdir(parents=True, exist_ok=True)
            (cand_out / "final_label.json").write_text(json.dumps({
                "candidate_id": cid,
                "q_overall_structural_class": label,
                "weakest_component": weakest,
                "justification": just,
                "source_script": "STEP_04_assign_structural_class.py",
            }, indent=2, default=str))

    catalog_path = outdir / "final_catalog.tsv"
    if not args.dry_run:
        with open(catalog_path, "w") as f:
            w = csv.DictWriter(
                f,
                fieldnames=["candidate_id", "q_overall_structural_class",
                            "weakest_component", "justification"],
                delimiter="\t")
            w.writeheader()
            for row in rows_out:
                w.writerow(row)
        print(f"[assign_v7] wrote {catalog_path}")

    print("[assign_v7] distribution:")
    for lbl, n in sorted(summary.items(), key=lambda kv: -kv[1]):
        print(f"  {n:5d}  {lbl}")


if __name__ == "__main__":
    main()
