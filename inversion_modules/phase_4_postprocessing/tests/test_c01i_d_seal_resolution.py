#!/usr/bin/env python3
"""
test_c01i_d_seal_resolution.py — unit tests for the per-sample class
resolution rules in STEP_C01i_d_seal.R.

Mirrors the resolve_final_classes() logic as a pure Python function and
exercises all four rules + edge cases. Any change to the R logic must be
reflected here (and vice versa) — this is the rules-as-code spec.
"""
import sys


def resolve_sample(pca_class, recomb_info=None, structure_type=None):
    """Python mirror of R resolve_final_classes() per-sample logic.

    Args:
        pca_class: str, one of HOM_REF/HET/HOM_INV from C01i_decompose
        recomb_info: dict or None. If dict, has event_class key.
        structure_type: str or None. One of the six nested_composition
                         classes, or None if nested_composition had no data.

    Returns:
        dict with final_class, recomb_status, recomb_event_class,
        recomb_source, in_composite_region.
    """
    is_recomb = recomb_info is not None
    recomb_source = "signal_detection" if is_recomb else None
    event_class = recomb_info.get("event_class") if is_recomb else None

    # Rule 2: fragmented structure elevates to recombinant
    if not is_recomb and structure_type == "multi_block_fragmented":
        is_recomb = True
        recomb_source = "nested_comp_fragmented"
        event_class = "ambiguous"

    # Rule 3: two_block flag (doesn't change class, just marks it)
    in_composite = structure_type == "two_block_composite"

    final_class = "RECOMBINANT" if is_recomb else pca_class

    status_map = {
        "gene_conversion": "recomb_GC",
        "double_crossover": "recomb_DCO",
        "suspicious": "recomb_suspicious",
        "ambiguous": "recomb_ambiguous",
    }
    recomb_status = status_map.get(event_class, "recomb_ambiguous") if is_recomb else "NOT_RECOMB"

    return {
        "final_class": final_class,
        "recomb_status": recomb_status,
        "recomb_event_class": event_class,
        "recomb_source": recomb_source,
        "in_composite_region": in_composite,
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

    print("test_c01i_d_seal_resolution:")

    # ─── Rule 4 (default): just take pca_class ──────────────────────────
    check("HOM_REF with no extras stays HOM_REF",
          resolve_sample("HOM_REF"),
          {"final_class": "HOM_REF", "recomb_status": "NOT_RECOMB",
           "in_composite_region": False})

    check("HET with no extras stays HET",
          resolve_sample("HET"),
          {"final_class": "HET", "recomb_status": "NOT_RECOMB"})

    check("HOM_INV with no extras stays HOM_INV",
          resolve_sample("HOM_INV"),
          {"final_class": "HOM_INV", "recomb_status": "NOT_RECOMB"})

    check("homogeneous structure does not change class",
          resolve_sample("HOM_REF", structure_type="homogeneous"),
          {"final_class": "HOM_REF", "in_composite_region": False})

    check("dominant_plus_secondary does not change class",
          resolve_sample("HET", structure_type="dominant_plus_secondary"),
          {"final_class": "HET", "in_composite_region": False})

    # ─── Rule 1: recomb_info promotes to RECOMBINANT ───────────────────
    check("recomb_info with GC promotes HET → RECOMBINANT (recomb_GC)",
          resolve_sample("HET", recomb_info={"event_class": "gene_conversion"}),
          {"final_class": "RECOMBINANT", "recomb_status": "recomb_GC",
           "recomb_event_class": "gene_conversion",
           "recomb_source": "signal_detection"})

    check("recomb_info with DCO promotes HOM_REF → RECOMBINANT (recomb_DCO)",
          resolve_sample("HOM_REF", recomb_info={"event_class": "double_crossover"}),
          {"final_class": "RECOMBINANT", "recomb_status": "recomb_DCO",
           "recomb_event_class": "double_crossover"})

    check("recomb_info with suspicious event",
          resolve_sample("HOM_INV", recomb_info={"event_class": "suspicious"}),
          {"final_class": "RECOMBINANT", "recomb_status": "recomb_suspicious"})

    check("recomb overrides pca_class even for HOM_INV",
          resolve_sample("HOM_INV", recomb_info={"event_class": "gene_conversion"}),
          {"final_class": "RECOMBINANT", "recomb_status": "recomb_GC"})

    # ─── Rule 2: fragmented structure type when no recomb detected ──────
    check("multi_block_fragmented with no recomb → RECOMBINANT via rule 2",
          resolve_sample("HET", structure_type="multi_block_fragmented"),
          {"final_class": "RECOMBINANT", "recomb_status": "recomb_ambiguous",
           "recomb_event_class": "ambiguous",
           "recomb_source": "nested_comp_fragmented"})

    check("multi_block_fragmented on HOM_REF → RECOMBINANT",
          resolve_sample("HOM_REF", structure_type="multi_block_fragmented"),
          {"final_class": "RECOMBINANT", "recomb_source": "nested_comp_fragmented"})

    # Rule 1 beats Rule 2: explicit recomb signal wins over fragmentation inference
    check("recomb_info overrides rule 2 source",
          resolve_sample("HET",
                          recomb_info={"event_class": "gene_conversion"},
                          structure_type="multi_block_fragmented"),
          {"final_class": "RECOMBINANT", "recomb_status": "recomb_GC",
           "recomb_source": "signal_detection"})   # signal, not fragmented

    # ─── Rule 3: two_block_composite marks but keeps class ──────────────
    check("two_block_composite keeps HOM_REF, flags in_composite_region",
          resolve_sample("HOM_REF", structure_type="two_block_composite"),
          {"final_class": "HOM_REF", "in_composite_region": True,
           "recomb_status": "NOT_RECOMB"})

    check("two_block_composite keeps HET, flags composite",
          resolve_sample("HET", structure_type="two_block_composite"),
          {"final_class": "HET", "in_composite_region": True})

    check("two_block + recomb → RECOMBINANT but STILL in_composite",
          resolve_sample("HET",
                          recomb_info={"event_class": "gene_conversion"},
                          structure_type="two_block_composite"),
          {"final_class": "RECOMBINANT", "in_composite_region": True,
           "recomb_status": "recomb_GC"})

    # ─── Other structure types pass through ─────────────────────────────
    check("continuous_gradient keeps class, no composite flag",
          resolve_sample("HET", structure_type="continuous_gradient"),
          {"final_class": "HET", "in_composite_region": False})

    check("diffuse_mixed keeps class",
          resolve_sample("HOM_INV", structure_type="diffuse_mixed"),
          {"final_class": "HOM_INV", "in_composite_region": False})

    check("too_sparse structure_type is a no-op",
          resolve_sample("HOM_REF", structure_type="too_sparse"),
          {"final_class": "HOM_REF", "in_composite_region": False})

    # ─── Edge cases ─────────────────────────────────────────────────────
    check("None structure_type is fine (no nested data)",
          resolve_sample("HET", structure_type=None),
          {"final_class": "HET"})

    check("unknown event_class maps to recomb_ambiguous",
          resolve_sample("HET", recomb_info={"event_class": "wacky"}),
          {"final_class": "RECOMBINANT", "recomb_status": "recomb_ambiguous"})

    check("empty recomb_info dict (no event_class) is still recomb",
          resolve_sample("HET", recomb_info={}),
          {"final_class": "RECOMBINANT", "recomb_status": "recomb_ambiguous"})

    # ─── Summary ────────────────────────────────────────────────────────
    if failures:
        print(f"\n✗ {len(failures)} failure(s)")
        return 1
    print("\n✓ all resolution tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(run_tests())
