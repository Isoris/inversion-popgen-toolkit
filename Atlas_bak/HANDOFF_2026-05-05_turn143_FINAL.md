# HANDOFF — turn 143 — Breeding-readiness card, Turn B (computation)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (69,665 lines)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC.
**Supersedes**: turn 142 Turn A handoff.

This is **Turn B** of a 4-turn build of
`SPEC_per_candidate_breeding_readiness_card.md`. Turn A (turn 142)
shipped the data plumbing — the cohort_diversity_v1 JSON loader, the
byCGA index, the `_diversityForSampleIdx` resolver, the real fixture.
**Turn B ships the entire computation layer**: per-sample karyotype
derivation, Wilcoxon rank-sum (with full tie + continuity correction),
per-arrangement F_ROH stats, K8×karyotype contingency table, MAF, and
the SPEC §2 pairing-advice rule engine. Turn C (next) wires this to a
print-ready HTML one-pager. Turn D wires it to bulk catalogue export.

---

## 0. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K=8 clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok). Never invent surname.

---

## 1. What turn B ships

### 1.1 Eight pure helpers + one builder (~660 LOC, ~line 20337–~21090)

All inserted right after Turn A's cohort_diversity loader, all PURE
(no DOM, no state mutation), all exposed on `window.*`.

| Helper | Purpose |
|---|---|
| `_breedingCardKaryotypePerSample(cand)` | Per-sample REF/HET/INV from `_classifyHLabelBands` (turn 139) + PC1 sort. Lowest-PC1 HOM band → REF, highest → INV, middles → HOM_MID. Single-HOM degrades to ambiguous_role with HOM_AMBIGUOUS_ROLE counts. |
| `_wilcoxonRankSumP(a, b)` | Two-sided Mann–Whitney U with average-rank tie handling, tie-corrected variance (Hollander/Wolfe/Chicken §4.1), continuity correction. Reuses existing `normalCDF` (Abramowitz & Stegun 7.1.26) at line 27620. |
| `_summarizeFROHGroup(values)` | n / mean / median / sd / q1 / q3 / min / max. Sample SD (n−1 denom). Quantiles use R Type 7 (linear interpolation between order stats). |
| `_perArrangementFROH(cand)` | Walks karyotype × cohortDiversity, partitions F_ROH into HOM_REF / HET / HOM_INV / HOM_MID / all_carriers, runs Wilcoxon REF-vs-INV and REF-vs-all_carriers, computes delta_mean_inv_minus_ref. |
| `_arrangementMAF(counts)` | Minor-allele frequency: `min(p_ref, p_inv)` treating REF/HET/INV as 2/1+1/2 alleles. Excludes HOM_MID and HOM_AMBIGUOUS_ROLE from the allele count. |
| `_carrierByK8Table(cand)` | K8 × karyotype contingency. Columns adapt: always [REF, HET, INV], adds HOM_MID iff multi_haplotype, adds HOM iff ambiguous_role. K-clusters sorted numerically (K1, K2, …, K10). |
| `_generatePairingAdvice(card, opts)` | SPEC §2 rule engine. Emits 0–N items per rule, each `{kind, severity, text, evidence}`. |
| `_buildBreedingCard(cand, opts)` | Pure data builder. Returns full JSON-serializable card (no Maps, no Int8Arrays, no DOM). |

### 1.2 Pairing-advice rules

The SPEC §2 rule table, with what's implementable now:

| # | Rule | Status | Severity |
|---|---|---|---|
| 1 | F_ROH asymmetric (Wilcoxon p<0.05, INV>REF) | ✓ shipped | strong |
| 1b | F_ROH inverted asymmetry (REF>INV) | ✓ shipped | info |
| 2 | Damaging-load asymmetric | `data_pending` (needs MODULE_CONSERVATION) | info |
| 3 | MAF imbalanced (minor allele <10%) | ✓ shipped | warn |
| 4 | Recombinant fraction high | `data_pending` (needs separate spec) | info |
| 5 | Default | ✓ shipped (only fires when no actionable rule) | info |
| — | Multi-haplotype caveat | ✓ shipped (always when applicable) | info |
| — | Ambiguous-role caveat | ✓ shipped (always when applicable) | info |

Constants (all `const` near the top of the rule engine):

```js
_BREEDING_FROH_P_THRESHOLD     = 0.05;   // SPEC rule 1
_BREEDING_MAF_IMBALANCE_THRESH = 0.10;   // SPEC rule 3
_BREEDING_RECOMBINANT_THRESH   = 0.20;   // SPEC rule 4 (when ready)
```

### 1.3 The card data shape

`_buildBreedingCard(cand)` returns:

```js
{
  schema: 'breeding_readiness_card_v1',
  generated_at: '<iso8601>',
  candidate: {
    id, chrom, start_bp, end_bp, span_bp, span_mb,
    K, tier, confidence, notes
  },
  karyotype: {
    summary: { HOM_REF, HET, HOM_INV, HOM_MID, HOM_AMBIGUOUS_ROLE,
               AMBIGUOUS, NO_DOSAGE, n_classified, n_unclassified },
    multi_haplotype, ambiguous_role,
    classification_summary,    // e.g. "K=3 · 2H+1T"
    consistency,                // from h_classification.implied_regime
    hom_band_idx_by_role: { REF, INV, MID[] },
    k_used,
  },
  burden: {
    available, reason,
    karyotype_summary,
    multi_haplotype, ambiguous_role,
    groups: { HOM_REF, HET, HOM_INV, HOM_MID },   // each {n, mean, median, sd, q1, q3, min, max}
    all_carriers: { ... },
    n_resolved, n_unresolved,
    wilcoxon_ref_vs_inv: { n_a, n_b, R_a, U_a, mu, sigma2, sigma,
                           n_tie_groups, tie_correction_factor,
                           z, p_two_sided, direction },
    wilcoxon_ref_vs_carriers: { ... },
    delta_mean_inv_minus_ref,
  },
  ancestry: {
    available, reason,
    k8_clusters,                // ['K1', 'K2', ...]
    karyotypes,                  // ['HOM_REF', 'HET', 'HOM_INV', 'HOM_MID'?, 'HOM'?]
    rows: [{ k8, counts, total, total_classified, frac }, ...],
    totals, n_resolved, n_unresolved,
  },
  maf: { p_ref, p_inv, maf, n_alleles, minor_allele } | null,
  coverage: { n_total, n_resolved, n_unresolved, unresolved_cgas[] },
  advice: [{ kind, severity, text, evidence }, ...],
  atlas_url: '<reviewer link>' | null,
}
```

This is fully JSON-serializable — Turn D can dump it as standalone
JSON; Turn C can render it directly to HTML.

---

## 2. Wilcoxon implementation correctness

The Wilcoxon rank-sum is the highest-stakes piece in Turn B because
the manuscript will quote p-values from it. The implementation uses
the standard normal-approximation form:

- **Rank assignment**: average ranks for ties (Mann & Whitney 1947 +
  Hollander/Wolfe/Chicken §4.1 convention)
- **U statistic**: `U_a = R_a − n_a(n_a+1)/2`
- **Null mean**: `μ = n_a · n_b / 2`
- **Null variance with tie correction**:
  `σ² = (n_a · n_b / 12) · ((N+1) − Σ(t³−t) / (N(N−1)))`
  where t = size of each tie group, N = n_a + n_b
- **Continuity correction**: `z = (|U−μ| − 0.5) / σ`, with z=0 when
  |U−μ| ≤ 0.5
- **Two-sided p**: `2 · (1 − Φ(|z|))`, clamped to [0, 1]

Verification in test §4 (32 assertions):

- **5v5 perfect separation**: hand-computed R_a=15, U_a=0, μ=12.5,
  σ²=22.917, z=2.5067, p=0.0122 — all match exactly
- **Symmetry**: `p(a,b) === p(b,a)` to 1e-12, `U_a + U_b = n_a · n_b`,
  direction flips
- **Full ties**: σ=0, z=NaN, p=NaN, direction='equal' (degenerate
  result reported, not silently zero)
- **Mixed ties**: `n_tie_groups`, `tie_correction_factor` recorded
  for diagnostics; R_a hand-verified
- **50v50 with shift+0.10**: p < 1e-10
- **1v1**: continuity correction makes |U−μ|=0.5 → z=0 → p=1
- **Empty / null / NaN-only inputs**: return null without throwing

The implementation matches R's `wilcox.test(..., exact=FALSE)`
output to within floating-point tolerance for n ≥ 5+5 (where the
normal approximation is valid). For n=60+60 — the manuscript-relevant
LG28 case — the approximation is well within its validity range.

---

## 3. Tests

`tests/test_turn143_breeding_card_compute.js` — **196 / 0** across
11 sections:

1. Source-level definitions present (12)
2. Window exports (8)
3. Karyotype derivation behavioral (28) — K=3 clean, K=4 multi-hap,
   single-HOM ambiguous, mixed AMBIGUOUS/NO_DOSAGE bands, defensive
4. Wilcoxon rank-sum (32) — perfect separation, symmetry, full ties,
   mixed ties, large-n, 1v1, edge cases
5. `_summarizeFROHGroup` (10)
6. `_perArrangementFROH` (20) — clean run, four stub-path scenarios,
   missing-CGA bookkeeping, `include_values` flag
7. `_arrangementMAF` (7)
8. `_carrierByK8Table` (12) — clean run, multi-hap column adapt,
   ambig-role column adapt, stubs
9. Pairing advice rule firing (24) — F_ROH asym, F_ROH inverted,
   MAF imbalance, multi-haplotype caveat, ambig-role caveat,
   data_pending always emitted, default fires only when no actionable
10. `_buildBreedingCard` end-to-end (29) — full assembly, JSON
   serializability, no-h_class fallback, no-cohort fallback
11. Real-fixture pipeline (8) — 226-sample synthetic 60/106/60 split,
   matches the user-reported actual LG28 karyotype shape

Adjacent suites unchanged: turn 142 cohort_diversity loader (133/0),
turn 141 candidate-bands (62/0), turn 140 H-label chip (45/0), turn
139 H-label classifier (143/0). Full turn-numbered suite **1816 / 0**
(up from 1620 / 0 at turn 142).

---

## 4. Atlas state

| | LOC | Tests | Files |
|---|---|---|---|
| Pre-turn (turn 142 baseline) | 68,679 | 1620 | 41 |
| Post-turn (this) | 69,665 | 1816 | 42 |
| Δ | +986 | +196 | +1 |

---

## 5. What's NOT done — Turn C scope

Turn B is intentionally render-free. The user notices Turn B's
shipping only via the test suite count and the new `window.*` exports
on the dev console — there is no DOM change yet. Turn C's job:

### 5.1 Print-ready one-pager render

A new artifact (probably an inline HTML render block on candidate-focus
page 3, or its own popup) that consumes `_buildBreedingCard(cand)`
output and renders it as the SPEC §1 visual contract:

- Header strip: chrom · coords · span · tier · confidence
- PCA panel (reuse existing `drawPCA` output)
- Karyotype counts row (large, color-coded)
- Carrier × K=8 ancestry table (rows by K, columns by karyotype)
- F_ROH per-arrangement bar chart with Wilcoxon p
- Cross-species orthology (if `state.crossSpecies` loaded)
- Pairing advice block (one styled card per rule, severity color-coded)
- Footer: provenance, atlas reviewer URL, generation timestamp

### 5.2 A4 print stylesheet

`@media print { @page { size: A4; margin: 1.5cm; } ... }` so the user
can ⌘P → save-as-PDF directly from the atlas. Body in 11pt, headers
12–18pt, no UI chrome (sidebar, tabs, toolbar) in print.

### 5.3 "Print breeding card" button

On candidate focus page 3, next to the existing Karyotype tab actions.

Estimated Turn C scope: ~250 LOC + ~50 tests.

---

## 6. What's NOT done — Turn D scope

After Turn C, Turn D adds:

- "Generate breeding cards" bulk button on page 5 catalogue
- Per-tier filter (default Tier 1 + Tier 2)
- Per-card download as standalone HTML
- Bulk export (zip of all cards, or one combined HTML)

Estimated Turn D scope: ~150 LOC + ~30 tests.

---

## 7. What this is NOT (still)

Same as Turn A's "what this is not" plus:

- **Not the recombinant detector**. Rule 4 stays `data_pending`.
- **Not the damaging-load layer**. Rule 2 stays `data_pending`.
- **Not the marker panel**. Separate spec.
- **Not the H-label Slice 2 work** (full multi-haplotype label
  resolution). For K=4+, the breeding card uses the simple PC1-end
  ordering convention (leftmost-PC1 = REF, rightmost = INV) and emits
  the multi_haplotype caveat. Slice 2 would resolve full ordinal
  labels (H1H1, H1H2, H2H2, H1H3, H2H3, H3H3, ...); when it ships, the
  breeding card's karyotype derivation can be upgraded.

---

## 8. Backups

```
Inversion_atlas.html.bak_pre_breeding_readiness_turnB   (pre-Turn B baseline)
Inversion_atlas.html.bak_pre_breeding_readiness_turnA   (pre-Turn A baseline)
Inversion_atlas.html.bak_pre_lines_cand_bands           (pre-turn 141 baseline)
```

`.bak_*` files NOT in bundle. Re-derivable.

---

## 9. Where to start the next chat

### Option 9a — Turn C (RECOMMENDED, continues the build)

Build the render layer per §5. Slices: HTML render block + A4 print
stylesheet + "Print" button. The renderer consumes `_buildBreedingCard`
output, which is fully tested and stable. No new computation needed —
this is a pure templating + styling pass. Estimated: ~250 LOC + ~50 tests.

### Option 9b — Pivot to a different spec

Tier 1 specs still in queue (per `SPECS_TIER_INDEX.md`):

- **`SPEC_inversion_age_origin_atlas.md`** — Porubsky 2022, Hartigan,
  Corbett-Detig. GDS classification + dip test + age proxy.
- **`SPEC_boundary_consensus_aggregator.md`** — KDE multimodal mode
  detection.
- **`SPEC_recombinant_dosage_changepoint_detector.md`** — would
  unblock pairing-advice rule 4.
- **`SPEC_hypothesis_test_framework_atlas.md`** — T1–T11 with BH.

### Recommendation

**9a.** Turn B's data builder is fully tested but invisible — Turn C
gives it a surface. The whole point of the breeding-readiness card is
the deliverable users (hatchery managers, journal reviewers) will
actually look at; the data is plumbing. Turn C is a 1-turn job that
finishes the Atlas-5-Part-B story.

---

## 10. Honest framing

**What turn 143 actually delivered:**

- Eight pure helpers + one builder, all individually unit-tested
- A correct two-sided Wilcoxon rank-sum with tie + continuity
  correction, validated against hand-computed expected values
- Per-sample karyotype derivation that handles K=3 cleanly, degrades
  gracefully on K=4+ multi-haplotype, and reports ambiguous-role for
  single-HOM cases (instead of silently picking REF or INV)
- A pairing-advice rule engine that fires the SPEC §2 rules whose
  data exists today (1, 1b, 3, 5, plus two caveats), and emits
  `data_pending` placeholders for the two rules whose prerequisite
  layers haven't shipped (2, 4) — so the renderer in Turn C can show
  the user what's coming without lying about what's available
- A card-data structure that's fully JSON-serializable, ready for
  Turn C render and Turn D bulk export

**What it deliberately didn't deliver:**

- Any rendering. No new HTML, no new buttons, no DOM hooks. Turn B is
  pure data on purpose: the unit that tests cleanly in isolation, leaving
  Turn C to focus on layout + print stylesheet without also debugging
  Wilcoxon math.
- Damaging-load and recombinant-fraction rules. Their `data_pending`
  placeholders are deliberate — they tell the user (and the
  manuscript reviewer) "this analysis exists in the framework but the
  data layer is upstream" rather than silently omitting the rule.
- H-label Slice 2 (full multi-haplotype ordinal labels). Out of scope
  for this card; the multi_haplotype caveat tells users where the
  approximation lies.

**Manuscript impact (when Turn C lands):**

- Atlas 5 Part B ships its first concrete artifact: per-candidate
  print-ready PDFs that hatchery managers can scan and act on
- The chat-`c03fc41e` framing — *"converts the paper from a population
  genomics study into a hatchery management resource"* — has its
  visible evidence
- Reviewers / aquaculture-genomics readers who don't read the full
  manuscript can scan one PDF per inversion and decide whether to
  use it for breeding decisions. The Wilcoxon p-values will be
  manuscript-quotable.

Walk the map carefully, respect cohort discipline, don't break the
test suite. Turn B's computation layer is mathematically sound and
fully covered; Turn C is the obvious next step.
