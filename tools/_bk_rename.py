#!/usr/bin/env python3
"""
_bk_rename.py — apply the chat-15 BK key rename pass.

Reads RENAME_MAP, walks the three canonical schemas' keys_extracted
directives + compute_candidate_status.R q1/q2/q6 vectors, replaces
each old key with its new name, and verifies no unintended collisions.

Run after:  python3 tools/_schema_check.py   (structural validation)
Run before: python3 tools/_code_field_check.py (after rename we re-verify
            every new `from` path still maps to a source-code field)
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path

# (schema_filename, old_key, new_key)
# Q1 keys kept as-is: q1_composite_flag (v10.1 canonical).
# Q6 keys kept as-is: q6_n_HOM_REF_prelim etc (class counts are already scientific).
RENAME_MAP = [
    # ───────── internal_dynamics.schema.json ─────────
    ("internal_dynamics.schema.json", "q2_decomp_status",    "q2_pca_decomp_status"),
    ("internal_dynamics.schema.json", "q2_decomp_reason",    "q2_pca_decomp_skip_reason"),
    ("internal_dynamics.schema.json", "q2_seed_source",      "q2_pca_seed_source"),
    ("internal_dynamics.schema.json", "q2_n_seed_conflict",  "q2_n_seed_class_conflicts"),
    ("internal_dynamics.schema.json", "q2_silhouette_score", "q2_pca_cluster_silhouette"),
    ("internal_dynamics.schema.json", "q2_decomp_quality",   "q2_pca_cluster_separation_flag"),
    ("internal_dynamics.schema.json", "q2_phase_concordance","q2_het_phase_support_fraction"),
    ("internal_dynamics.schema.json", "q2_n_phase_supported","q2_n_het_phase_supported"),
    ("internal_dynamics.schema.json", "q2_flashlight_mode",  "q2_pca_seed_mode"),
    ("internal_dynamics.schema.json", "q2_n_discordant",     "q2_n_samples_sv_pca_discordant"),
    ("internal_dynamics.schema.json", "q2_n_cheat2_reclass", "q2_n_hetDEL_breakpoint_reclass"),
    ("internal_dynamics.schema.json", "q2_k_used",           "q2_pca_k_used"),
    # q2_bic_gap_k3_vs_k2 kept as-is (the name IS the measurement)
    # q2_n_seed_HOM_REF / HET / HOM_INV kept as-is (class counts are scientific)
    # q6_n_HOM_REF_prelim / HET_prelim / HOM_INV_prelim kept as-is

    # ───────── recombinant_map.schema.json ─────────
    ("recombinant_map.schema.json", "q2_n_recombinants",        "q2_n_recombinant_samples"),
    ("recombinant_map.schema.json", "q2_n_recombinant_gc",      "q2_n_recombinant_gc_samples"),
    ("recombinant_map.schema.json", "q2_n_recombinant_dco",     "q2_n_recombinant_dco_samples"),
    ("recombinant_map.schema.json", "q2_n_recomb_disputed",     "q2_n_regime_ghsl_disputed_samples"),
    ("recombinant_map.schema.json", "q2_n_recomb_ghsl_only",    "q2_n_ghsl_split_only_samples"),
    ("recombinant_map.schema.json", "q2_cheat24_version",       "q2_recomb_posterior_source"),
    ("recombinant_map.schema.json", "q2_combination_rule",      "q2_recomb_gate_rule_version"),
    ("recombinant_map.schema.json", "q2_recomb_min_dev_fraction","q2_recomb_min_regime_dev_fraction"),
    ("recombinant_map.schema.json", "q2_recomb_min_dev_bp",     "q2_recomb_min_regime_dev_bp"),
    ("recombinant_map.schema.json", "q2_recomb_min_dco_bp",     "q2_recomb_dco_threshold_bp"),

    # ───────── internal_ancestry_composition.schema.json ─────────
    # q1_composite_flag kept as-is (v10.1 canonical rename)
    ("internal_ancestry_composition.schema.json", "q1_ancestry_dominant_structure", "q1_ancestry_dominant_pattern"),
    ("internal_ancestry_composition.schema.json", "q2_pct_samples_two_block",    "q2_pct_samples_two_block_ancestry"),
    ("internal_ancestry_composition.schema.json", "q2_pct_samples_multi_block",  "q2_pct_samples_multi_block_ancestry"),
    ("internal_ancestry_composition.schema.json", "q2_pct_samples_homogeneous",  "q2_pct_samples_homogeneous_ancestry"),
    ("internal_ancestry_composition.schema.json", "q2_pct_samples_gradient",     "q2_pct_samples_gradient_ancestry"),
    ("internal_ancestry_composition.schema.json", "q2_pct_samples_diffuse",      "q2_pct_samples_diffuse_ancestry"),
    ("internal_ancestry_composition.schema.json", "q2_mean_fragmentation",       "q2_mean_ancestry_fragmentation"),
    ("internal_ancestry_composition.schema.json", "q2_mean_internal_entropy",    "q2_mean_ancestry_entropy_within_sample"),
    ("internal_ancestry_composition.schema.json", "q2_mean_switches",            "q2_mean_ancestry_switches_per_sample"),
    # q2_ancestry_K_used and q2_ancestry_n_samples_analyzed kept as-is
]

SCHEMA_DIR = Path("registries/schemas/structured_block_schemas")
SPEC_FILE = Path("inversion_modules/phase_4_postprocessing/4g_final_classification/compute_candidate_status.R")


def main() -> int:
    # 1. Update each schema file's keys_extracted in place
    per_schema: dict[str, list[tuple[str, str]]] = {}
    for fname, old, new in RENAME_MAP:
        per_schema.setdefault(fname, []).append((old, new))

    for fname, pairs in per_schema.items():
        p = SCHEMA_DIR / fname
        schema = json.loads(p.read_text())
        kex = schema.get("keys_extracted", [])
        renamed = 0
        for ke in kex:
            k = ke.get("key")
            for old, new in pairs:
                if k == old:
                    ke["key"] = new
                    renamed += 1
                    break
        # Preserve pretty formatting
        p.write_text(json.dumps(schema, indent=2) + "\n")
        print(f"  [{fname}] renamed {renamed}/{len(pairs)} keys")

    # 2. Update compute_candidate_status.R — replace quoted occurrences of each old name
    text = SPEC_FILE.read_text()
    total_in_r = 0
    for _, old, new in RENAME_MAP:
        old_quoted = f'"{old}"'
        new_quoted = f'"{new}"'
        if old_quoted in text:
            count = text.count(old_quoted)
            text = text.replace(old_quoted, new_quoted)
            total_in_r += count
    SPEC_FILE.write_text(text)
    print(f"  [compute_candidate_status.R] replaced {total_in_r} quoted occurrences")

    # 3. Verify no accidental collisions (any new name that was already in the target)
    all_new = {new for _, _, new in RENAME_MAP}
    for fname in set(f for f, _, _ in RENAME_MAP):
        schema = json.loads((SCHEMA_DIR / fname).read_text())
        seen: set[str] = set()
        dupes: list[str] = []
        for ke in schema.get("keys_extracted", []):
            if ke["key"] in seen:
                dupes.append(ke["key"])
            seen.add(ke["key"])
        if dupes:
            print(f"  ERROR [{fname}] duplicate keys after rename: {dupes}")
            return 1

    print("\nAll renames applied and verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
