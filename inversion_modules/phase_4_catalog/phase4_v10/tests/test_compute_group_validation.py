#!/usr/bin/env python3
"""
test_compute_group_validation.py — unit test for the pure-function version
of the C01f validation-promotion logic.

Mirrors the R function compute_group_validation() in
patches/C01f/C01f_comp_from_registry_patch.R. This Python port exists so
the logic can be tested without an R runtime on the test machine. Any
change to the R function must be mirrored here (and vice versa).

Test coverage:
  - starts at NONE with no evidence → stays NONE
  - UNCERTAIN → SUPPORTED on T8+T9
  - UNCERTAIN → VALIDATED on Layer D
  - SUPPORTED → VALIDATED on Layer D
  - ANY → SUSPECT on T9 fragile (demotion wins)
  - edge cases: missing fields, non-finite values, threshold boundaries
"""
import math
import sys


def compute_group_validation(
    current_level="UNCERTAIN",
    t8_concordance=None,
    t9_jackknife_status=None,
    t9_max_delta=None,
    t10_theta_concordance=None,
    layer_d_fisher_p=None,
    layer_d_fisher_or=None,
):
    """Python mirror of the R function. Must stay in sync."""

    def finite(x):
        return x is not None and isinstance(x, (int, float)) and math.isfinite(x)

    if not current_level:
        current_level = "UNCERTAIN"

    # 1. Demotion first
    is_fragile = (
        t9_jackknife_status == "fragile"
        or (finite(t9_max_delta) and t9_max_delta > 0.3)
    )
    if is_fragile:
        return "SUSPECT"

    # 2. Candidate promotions
    candidates = []

    if (finite(layer_d_fisher_p) and finite(layer_d_fisher_or)
            and layer_d_fisher_p < 0.05 and layer_d_fisher_or > 5):
        candidates.append("VALIDATED")

    if (finite(t8_concordance) and t8_concordance >= 0.70
            and t9_jackknife_status == "robust"):
        candidates.append("SUPPORTED")

    # 3. Take max of current + candidates
    rank = {"NONE": 0, "UNCERTAIN": 1, "SUPPORTED": 2, "VALIDATED": 3,
            "SUSPECT": -1}
    all_levels = [current_level] + candidates
    return max(all_levels, key=lambda lvl: rank.get(lvl, 0))


# =============================================================================
# Tests
# =============================================================================
def run_tests():
    failures = []

    def check(name, got, want):
        if got == want:
            print(f"  ✓ {name}")
        else:
            failures.append((name, got, want))
            print(f"  ✗ {name}: got {got!r}, want {want!r}")

    print("test_compute_group_validation:")

    # ── Baseline ──
    check("NONE + no evidence stays NONE",
          compute_group_validation("NONE"),
          "NONE")

    check("UNCERTAIN + no evidence stays UNCERTAIN",
          compute_group_validation("UNCERTAIN"),
          "UNCERTAIN")

    check("empty current defaults to UNCERTAIN",
          compute_group_validation(""),
          "UNCERTAIN")

    # ── SUPPORTED promotion ──
    check("T8=0.70 + robust → SUPPORTED",
          compute_group_validation(
              "UNCERTAIN", t8_concordance=0.70,
              t9_jackknife_status="robust"),
          "SUPPORTED")

    check("T8=0.85 + robust → SUPPORTED",
          compute_group_validation(
              "UNCERTAIN", t8_concordance=0.85,
              t9_jackknife_status="robust"),
          "SUPPORTED")

    check("T8=0.69 + robust does NOT promote (below threshold)",
          compute_group_validation(
              "UNCERTAIN", t8_concordance=0.69,
              t9_jackknife_status="robust"),
          "UNCERTAIN")

    check("T8=0.85 but jackknife=moderate does NOT promote",
          compute_group_validation(
              "UNCERTAIN", t8_concordance=0.85,
              t9_jackknife_status="moderate"),
          "UNCERTAIN")

    # ── VALIDATED promotion ──
    check("Fisher p=0.01, OR=10 → VALIDATED",
          compute_group_validation(
              "UNCERTAIN", layer_d_fisher_p=0.01,
              layer_d_fisher_or=10),
          "VALIDATED")

    check("Fisher p=0.049, OR=5.1 → VALIDATED (boundary)",
          compute_group_validation(
              "UNCERTAIN", layer_d_fisher_p=0.049,
              layer_d_fisher_or=5.1),
          "VALIDATED")

    check("Fisher p=0.051, OR=10 does NOT promote (p above threshold)",
          compute_group_validation(
              "UNCERTAIN", layer_d_fisher_p=0.051,
              layer_d_fisher_or=10),
          "UNCERTAIN")

    check("Fisher p=0.001, OR=5.0 does NOT promote (OR at threshold, need >5)",
          compute_group_validation(
              "UNCERTAIN", layer_d_fisher_p=0.001,
              layer_d_fisher_or=5.0),
          "UNCERTAIN")

    # ── Monotonicity: already VALIDATED stays VALIDATED ──
    check("SUPPORTED + new Layer D → VALIDATED",
          compute_group_validation(
              "SUPPORTED", layer_d_fisher_p=0.01,
              layer_d_fisher_or=10),
          "VALIDATED")

    check("VALIDATED + weaker evidence stays VALIDATED",
          compute_group_validation(
              "VALIDATED", t8_concordance=0.50),
          "VALIDATED")

    # ── SUSPECT demotion wins ──
    check("fragile jackknife demotes to SUSPECT",
          compute_group_validation(
              "UNCERTAIN", t9_jackknife_status="fragile"),
          "SUSPECT")

    check("T9 delta=0.35 demotes to SUSPECT (> 0.3)",
          compute_group_validation(
              "UNCERTAIN", t9_max_delta=0.35),
          "SUSPECT")

    check("T9 delta=0.30 does NOT demote (at threshold, need >0.3)",
          compute_group_validation(
              "UNCERTAIN", t9_max_delta=0.30),
          "UNCERTAIN")

    # Demotion beats promotion
    check("fragile + strong Layer D → SUSPECT (demotion wins)",
          compute_group_validation(
              "UNCERTAIN", t9_jackknife_status="fragile",
              layer_d_fisher_p=0.001, layer_d_fisher_or=50),
          "SUSPECT")

    check("VALIDATED + new fragile → SUSPECT",
          compute_group_validation(
              "VALIDATED", t9_jackknife_status="fragile"),
          "SUSPECT")

    # ── Missing/nonfinite data ──
    check("None p-value does not promote",
          compute_group_validation(
              "UNCERTAIN", layer_d_fisher_p=None,
              layer_d_fisher_or=10),
          "UNCERTAIN")

    check("inf p-value does not promote",
          compute_group_validation(
              "UNCERTAIN", layer_d_fisher_p=float("inf"),
              layer_d_fisher_or=10),
          "UNCERTAIN")

    check("nan p-value does not promote",
          compute_group_validation(
              "UNCERTAIN", layer_d_fisher_p=float("nan"),
              layer_d_fisher_or=10),
          "UNCERTAIN")

    # ── Combined: both promotions trigger, take max ──
    check("both SUPPORTED and VALIDATED conditions → VALIDATED",
          compute_group_validation(
              "UNCERTAIN",
              t8_concordance=0.85, t9_jackknife_status="robust",
              layer_d_fisher_p=0.001, layer_d_fisher_or=20),
          "VALIDATED")

    if failures:
        print(f"\n✗ {len(failures)} failure(s)")
        for name, got, want in failures:
            print(f"  {name}: got {got!r}, want {want!r}")
        return 1
    print("\n✓ all tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(run_tests())
