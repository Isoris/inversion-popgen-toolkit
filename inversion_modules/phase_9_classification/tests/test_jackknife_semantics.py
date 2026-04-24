#!/usr/bin/env python3
"""
test_jackknife_semantics.py — unit tests for the richer compute_group_validation
with four-way jackknife classification.

Mirrors the R function compute_group_validation() in
patches/C01f_jackknife_semantics_patch.R as a Python pure function.
Tests distinguish between:
  - pca_family_confounded → SUSPECT (the OLD v10 "fragile" behavior)
  - single_family_fragile → NOT SUSPECT; tagged family_specific_polymorphism,
                             capped at SUPPORTED
  - few_family → NOT SUSPECT; tagged restricted_family_spread,
                             capped at SUPPORTED
  - multi_family_contributing → tagged partial_robustness, allowed to promote
  - robust_multi_family → no restrictions
"""
import math
import sys


def compute_group_validation(
    current_level="UNCERTAIN",
    t8_concordance=None,
    t9_jackknife_verdict=None,
    t9_max_delta=None,
    t10_theta_concordance=None,
    layer_d_fisher_p=None,
    layer_d_fisher_or=None,
    promotion_cap=None,
):
    """Python mirror of the R function in
    patches/C01f_jackknife_semantics_patch.R. Returns dict with
    level, quality_flags, family_linkage.
    """
    def finite(x):
        return x is not None and isinstance(x, (int, float)) and math.isfinite(x)

    if not current_level:
        current_level = "UNCERTAIN"

    quality_flags = []
    family_linkage = "unknown"
    internal_cap = None

    jk = (t9_jackknife_verdict or "").lower()

    if jk == "robust_multi_family":
        family_linkage = "multi_family"
    elif jk == "multi_family_contributing":
        family_linkage = "multi_family"
        quality_flags.append("partial_robustness")
    elif jk == "few_family":
        family_linkage = "few_family"
        quality_flags.append("restricted_family_spread")
        internal_cap = "SUPPORTED"
    elif jk == "single_family_fragile":
        family_linkage = "single_family"
        quality_flags.append("family_specific_polymorphism")
        internal_cap = "SUPPORTED"
    elif jk in ("fragile", "pca_family_driven", "pca_family_confounded"):
        family_linkage = "pca_family_confounded"
        quality_flags.append("pca_groups_family_confounded")
        return {
            "level": "SUSPECT",
            "quality_flags": quality_flags,
            "family_linkage": family_linkage,
        }

    if not jk and finite(t9_max_delta) and t9_max_delta > 0.3:
        quality_flags.append("high_jackknife_delta_no_verdict")
        internal_cap = "SUPPORTED"
        family_linkage = "uncertain"

    candidates = []
    if (finite(layer_d_fisher_p) and finite(layer_d_fisher_or)
            and layer_d_fisher_p < 0.05 and layer_d_fisher_or > 5):
        candidates.append("VALIDATED")
    if (finite(t8_concordance) and t8_concordance >= 0.70
            and jk == "robust_multi_family"):
        candidates.append("SUPPORTED")
    if (finite(t8_concordance) and t8_concordance >= 0.70
            and jk in ("single_family_fragile", "few_family")):
        candidates.append("SUPPORTED")

    rank = {"NONE": 0, "UNCERTAIN": 1, "SUPPORTED": 2, "VALIDATED": 3,
            "SUSPECT": -1}
    best = max([current_level] + candidates, key=lambda l: rank.get(l, 0))

    for cap in [promotion_cap, internal_cap]:
        if cap and cap in rank and rank[best] > rank[cap]:
            best = cap

    return {
        "level": best,
        "quality_flags": quality_flags if quality_flags else ["normal"],
        "family_linkage": family_linkage,
    }


def run_tests():
    failures = []

    def check(name, got, want):
        ok = all(got.get(k) == v for k, v in want.items())
        if ok:
            print(f"  ✓ {name}")
        else:
            failures.append((name, got, want))
            print(f"  ✗ {name}")
            print(f"      got:  {got}")
            print(f"      want: {want}")

    print("test_jackknife_semantics:")

    # ─── multi_family_contributing: the GOOD case ─────────────────────
    check("robust_multi_family + strong Layer D → VALIDATED, multi_family",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="robust_multi_family",
              layer_d_fisher_p=0.001, layer_d_fisher_or=20),
          {"level": "VALIDATED", "family_linkage": "multi_family",
           "quality_flags": ["normal"]})

    check("robust_multi_family + T8=0.85 → SUPPORTED",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="robust_multi_family",
              t8_concordance=0.85),
          {"level": "SUPPORTED", "family_linkage": "multi_family"})

    check("multi_family_contributing → tagged partial_robustness",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="multi_family_contributing"),
          {"level": "UNCERTAIN", "family_linkage": "multi_family",
           "quality_flags": ["partial_robustness"]})

    # ─── single_family_fragile: the REAL case that used to be wrongly SUSPECT
    check("single_family_fragile alone → UNCERTAIN + family tag",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="single_family_fragile"),
          {"level": "UNCERTAIN", "family_linkage": "single_family",
           "quality_flags": ["family_specific_polymorphism"]})

    check("single_family_fragile + T8=0.80 → SUPPORTED (NEW: allowed)",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="single_family_fragile",
              t8_concordance=0.80),
          {"level": "SUPPORTED", "family_linkage": "single_family"})

    check("single_family_fragile + strong Layer D → capped at SUPPORTED",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="single_family_fragile",
              layer_d_fisher_p=0.001, layer_d_fisher_or=20),
          {"level": "SUPPORTED", "family_linkage": "single_family"})

    check("single_family_fragile does NOT become SUSPECT",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="single_family_fragile"),
          {"family_linkage": "single_family"})
    # (explicit negative check to confirm single_family → not SUSPECT)

    # ─── few_family: similar treatment to single_family ────────────────
    check("few_family → capped at SUPPORTED",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="few_family",
              layer_d_fisher_p=0.001, layer_d_fisher_or=20),
          {"level": "SUPPORTED", "family_linkage": "few_family",
           "quality_flags": ["restricted_family_spread"]})

    check("few_family + T8=0.75 → SUPPORTED",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="few_family",
              t8_concordance=0.75),
          {"level": "SUPPORTED", "family_linkage": "few_family"})

    # ─── pca_family_confounded: the REAL SUSPECT case ─────────────────
    check("fragile (legacy) → SUSPECT",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="fragile"),
          {"level": "SUSPECT", "family_linkage": "pca_family_confounded"})

    check("pca_family_confounded → SUSPECT",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="pca_family_confounded"),
          {"level": "SUSPECT", "family_linkage": "pca_family_confounded"})

    check("SUSPECT beats any promotion evidence",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="pca_family_confounded",
              layer_d_fisher_p=0.001, layer_d_fisher_or=20,
              t8_concordance=0.99),
          {"level": "SUSPECT"})

    # ─── promotion_cap (external, from composite_flag) still works ─────
    check("robust_multi_family + external cap UNCERTAIN → UNCERTAIN",
          compute_group_validation("UNCERTAIN",
              t9_jackknife_verdict="robust_multi_family",
              layer_d_fisher_p=0.001, layer_d_fisher_or=20,
              promotion_cap="UNCERTAIN"),
          {"level": "UNCERTAIN", "family_linkage": "multi_family"})

    # ─── high delta without verdict → cap at SUPPORTED ────────────────
    check("high t9_max_delta no verdict → capped at SUPPORTED",
          compute_group_validation("UNCERTAIN",
              t9_max_delta=0.5,
              layer_d_fisher_p=0.001, layer_d_fisher_or=20),
          {"level": "SUPPORTED", "family_linkage": "uncertain",
           "quality_flags": ["high_jackknife_delta_no_verdict"]})

    # ─── unknown verdict → treat as unknown + standard promotion ─────
    check("no verdict, no evidence → UNCERTAIN unknown",
          compute_group_validation("UNCERTAIN"),
          {"level": "UNCERTAIN", "family_linkage": "unknown"})

    if failures:
        print(f"\n✗ {len(failures)} failure(s)")
        return 1
    print("\n✓ all jackknife-semantics tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(run_tests())
