# SPEC — XP-EHH per-window track (popstats / ancestry pages)

**Status**: forward-looking spec. Not implemented. Created same turn as
`SPEC_focal_vs_background_widget.md`.

**Reading order**: this spec → `SPEC_focal_vs_background_widget.md` (related but
distinct: that's the focal-vs-background contrast widget; this is a per-window
selection-scan track) → page-8 popstats renderer (where this track will land).

**One-line summary**: add an XP-EHH (cross-population extended haplotype
homozygosity) per-window track to the popstats and ancestry pages. Blocked on
two prerequisites that don't exist yet.

---

## Purpose

XP-EHH (Sabeti et al. 2007) is a haplotype-based selection scan that detects
**recent positive selection** in one cohort relative to a reference cohort.
Inversions suppress recombination, which extends haplotype homozygosity inside
the inverted region for one arrangement. XP-EHH is therefore a sensitive
inversion-age and inversion-selection signal — high XP-EHH near inversion
boundaries suggests one arrangement is sweeping or has swept recently in one
cohort but not the other.

In the manuscript context: an XP-EHH track adjacent to F_ST and dXY on the
popstats page lets the user read **age + suppression + selection** in one
view.

---

## Why this matters

The atlas already shows F_ST and dXY tracks (or will, after the
focal-vs-background widget lands). These tell you:
- **F_ST** = arrangement frequency divergence + recombination suppression
- **dXY** = absolute divergence between arrangements (age proxy)

XP-EHH tells you something neither of those does:
- **XP-EHH** = haplotype-homozygosity asymmetry — *which* arrangement is
  carrying the long haplotypes (i.e. which one is younger / under selection)

Together: F_ST high + dXY low + XP-EHH asymmetric = **young, sweeping
inversion**. F_ST high + dXY high + XP-EHH symmetric = **old, balanced
inversion**. The three tracks separate cases the focal-vs-bg widget cannot.

---

## What this track DOES claim

- Per-window XP-EHH score, displayed alongside other popstats tracks
- Sign of XP-EHH indicates which cohort/arrangement carries longer
  haplotypes (positive = test cohort has longer EHH; negative = reference
  cohort does)
- A genome-wide XP-EHH percentile threshold for "outlier" calls (default
  top 1%, configurable)

## What this track does NOT claim

- Selection coefficient or sweep age (XP-EHH is qualitative for those, not
  quantitative)
- Causation between any specific variant and the haplotype signal
- That positive XP-EHH at an inversion boundary is selection — it could
  equally be founder effect, drift, or assortative mating in the test
  cohort

---

## Prerequisites (NOT YET MET)

XP-EHH cannot be computed from the current dosage matrix. Two upstream
prerequisites:

### 1. Phased haplotypes for the 226-sample cohort

Currently: BEAGLE-imputed dosage (allele dosage 0/1/2 per sample per site).
Required: phased VCF (haplotypes per sample). This is one BEAGLE step away —
BEAGLE 5.x outputs phased VCF natively when run in the right mode — but
hasn't been run yet for the 226-sample cohort.

Action: rerun BEAGLE 5.x with phased output, persist a `phased_226.vcf.gz`
on LANTA. ~24h compute; output is ~30 GB.

### 2. A reference cohort for the "X" in XP-EHH

XP-EHH is a **two-cohort** statistic. The 226-sample hatchery cohort alone
cannot produce it. Options:
- (a) **Wild *C. gariepinus* cohort** — if/when it exists. Currently zero
  samples on hand. This is the cleanest reference and the right manuscript
  framing.
- (b) **Karyotype class within the 226-sample cohort** — i.e. compare
  HOM_REF arrangement haplotypes vs HOM_INV arrangement haplotypes at each
  candidate inversion. This is a **within-cohort, between-arrangement**
  XP-EHH and is conceptually closer to the "inversion-age" reading. It
  requires phasing only (no second cohort), but interpretation is
  different — it's recombination-suppression-flavoured, not
  selection-flavoured.
- (c) **C. macrocephalus wild cohort** — Quentin's planned future cohort.
  Wrong reference for selection in the Gar lineage (different species),
  but right reference for *cross-species* comparative selection if the
  manuscript scope ever broadens.

**Day-1 path of least resistance**: option (b) — XP-EHH between karyotype
classes within the 226-sample cohort. Re-uses the karyotype assignments
already on page 7. Requires only the phasing step.

**Reviewer-grade path**: option (a) — a real wild reference cohort.
Several-month timeline.

This spec proceeds with option (b) as the day-1 framing and notes the
dependency on (a) for the eventual manuscript.

---

## Architecture

### Layer schema

```jsonc
// File: xpehh_per_window.json
{
  "schema_version": 1,
  "generated_at": "2026-MM-DD",
  "test_cohort_id":  "Gar_226_HOM_INV_at_LG28_inv",
  "ref_cohort_id":   "Gar_226_HOM_REF_at_LG28_inv",
  "n_test":          60,
  "n_ref":           60,
  "tool":            "selscan v2.0.x",
  "tool_args":       "--xpehh --vcf phased_226.vcf.gz --vcf-ref ... --map ...",
  "windows": {
    "<chrom>": {
      "window_start_bp": [...],
      "window_end_bp":   [...],
      "xpehh_mean":      [...],   // per-window mean of per-SNP XP-EHH
      "xpehh_max_abs":   [...],   // per-window max |XP-EHH|
      "n_snps":          [...],
      "norm_xpehh_mean": [...],   // genome-wide-normalized (z-score against autosomal background)
    }
  }
}
```

Per-SNP XP-EHH is what selscan emits; the atlas needs it aggregated to the
same per-window grid as the other popstats tracks (otherwise the tracks
don't visually align). Aggregation: mean of per-SNP XP-EHH within window,
plus max-absolute for the most-extreme-haplotype reading. Z-normalization
against autosomal background follows Sabeti's recommendation (selscan has
a `--norm` step).

### Atlas-side wiring

#### Popstats page (page 8)

Add `xpehh` to the popstats track list. Renders as a stacked panel, same
visual idiom as F_ST / dXY / Hobs:

- Y-axis: normalized XP-EHH (Z-scored)
- Two horizontal dashed lines at ±2 (top/bottom 2.5% under normal
  approximation) — outlier visual threshold
- Color: same accent as F_ST track
- Hover tooltip: window coords, raw XP-EHH, percentile rank, n_snps

Track order: F_ST → dXY → Hobs → **xpehh** → θπ. XP-EHH belongs near the
selection-flavoured tracks, not next to coverage.

#### Ancestry page (page 9)

XP-EHH is conceptually an ancestry-asymmetry signal. Adding it as a
sub-panel here is appropriate. Lower priority than the popstats mount —
ship popstats first.

#### Focal-vs-background widget integration

When the focal-vs-bg widget's metric pool gains `xpehh` (after this layer
lands), the widget's tail mode should default to `'two-sided'` — both
positive and negative XP-EHH peaks at boundaries are biologically
interesting (different sweeps).

### Page integration NOT in scope this spec

- The cross-species page 16 is **NOT** in scope here. XP-EHH is a
  within-species, two-cohort statistic. Cross-species XP-EHH between Gar
  and Mac is meaningless because there are no shared variants.
- Boundaries page 11 is **NOT** in scope. XP-EHH is a track-display
  thing (popstats), not a refinement thing. The focal-vs-bg widget
  consumes XP-EHH on page 11 once the layer is loaded.

---

## Module structure (suggested layout)

```
phase_Y_xpehh/
  README.md
  00_phasing/
    STEP_PHA_01_run_beagle5_phased.sh       # rerun BEAGLE for phased VCF
  01_xpehh/
    STEP_XPE_01_split_by_karyotype.py        # haplotype subsets per arrangement per candidate
    STEP_XPE_02_run_selscan_xpehh.sh         # selscan invocation, per chromosome
    STEP_XPE_03_normalize_xpehh.sh           # selscan --norm
    STEP_XPE_04_aggregate_per_window.R       # selscan output → per-window
    STEP_XPE_05_export_xpehh_json.py         # → atlas layer
  tests/
    test_aggregation.R
```

---

## Tests

- `xpehh_layer_loads_correctly`: layer JSON drops in, atlas state acquires
  it, popstats track renders.
- `xpehh_track_aligns_with_fst_track`: window grid matches F_ST exactly
  (same window_start_bp/window_end_bp arrays); visual alignment passes.
- `xpehh_outlier_threshold_lines_render_at_z_2`: dashed lines at ±2 visible.
- `xpehh_handles_missing_chrom`: chrom not in JSON → empty panel + hint.
- `xpehh_two_cohort_metadata_displayed_in_track_header`:
  test/ref cohort labels visible in the track title (so user knows which
  way the sign points).

---

## Open questions / deferred

- **Wild reference cohort** — when/if it exists, regenerate the layer and
  the manuscript framing changes from "between karyotype classes
  (recombination suppression flavour)" to "between hatchery and wild
  (selection flavour)." The atlas track does not need to change; only
  cohort labels in the JSON metadata.
- **iHS** (intra-cohort selection scan, single-cohort EHH-based) — also
  worth adding, doesn't require a second cohort, but requires phased
  haplotypes. Could ship simultaneously with the phasing prerequisite at
  near-zero extra cost. Flagged here, not specced separately.
- **n_eff sample size for XP-EHH stability**: Sabeti recommends n ≥ 30
  per cohort for stable XP-EHH. Within-karyotype-class sample sizes at
  Gar candidates: HOM_REF ≈ 60, HOM_INV ≈ 60, HET ≈ 106 (LG28 numbers
  from STATUS_REPORT). Comfortable for HOM_REF / HOM_INV pairs at
  high-frequency candidates; rare-allele candidates with HOM_INV count
  < 30 should NOT compute XP-EHH (gate at the layer-build step).

---

## Summary

XP-EHH track on popstats (and later ancestry) page. Blocked on phased VCF
+ choice of reference cohort. Day-1 path: phase 226-sample cohort, compute
XP-EHH between HOM_REF and HOM_INV haplotypes within each candidate.
Reviewer-grade path: wait for wild *C. gariepinus* cohort. Atlas-side
wiring is small once the layer JSON exists.

End of spec.
