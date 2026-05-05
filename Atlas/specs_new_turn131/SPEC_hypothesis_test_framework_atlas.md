# SPEC — Hypothesis test framework atlas integration (T1–T11)

**Status**: drafted turn 130 final session. The hypothesis tests
T1–T11 already exist as `STEP_C01f_hypothesis_tests.R` on LANTA. This
spec defines how to **surface them in the atlas** as a per-candidate
verdict panel.

**Trigger** (from multiple chats, see `bf1625f9`, `c1a88579`):
> *"Hypothesis tests T1–T11: Within vs between band relatedness, ancestry
> diversity per band, kin-pruned signal retention, inner vs outer sample
> composition overlap, anchor stability across sub-regions, regime
> change at inner/outer boundary, within-carrier-only substructure,
> Clair3 indel genotype concordance (Cheat 9), ancestry jackknife
> (Cheat 6), theta het prior (Cheat 12), extended boundary suppression
> (Cheat 23)."*
>
> *"After Benjamini–Hochberg correction across the 10-test framework
> (q < 0.05), [N] putative inversions were confirmed."*

The atlas should display, per candidate, the verdicts of all 11 tests
plus the BH-corrected synthesis. Current atlas has tier-based
classification but doesn't surface the hypothesis-test trail.

---

## 1. The 11 tests (recap)

From chat `5d2ca953` (codebase consolidation):

| Test | Name | What it asks |
|---|---|---|
| **T1** | Within vs between band relatedness | Are samples within a band more related than between bands? (H1 family) |
| **T2** | Ancestry diversity per band | Do bands track Q components or are they ancestry-mixed? (H1 vs H2/H3) |
| **T3** | Kin-pruned signal retention | Does the signal survive when relatives are removed? (H1 vs H2) |
| **T4** | Inner vs outer sample composition | Do outer and inner halves have same band membership? (H2 simple vs H3 nested) |
| **T5** | Anchor stability across sub-regions | Do PCA anchors persist? (H2 hot core stability) |
| **T6** | Regime change at inner/outer boundary | Does the regime structure switch? (H3 nested) |
| **T7** | Within-carrier-only substructure | Do carriers form sub-clusters? (H4 sub-haplotype) |
| **T8** | Clair3 indel genotype concordance (Cheat 9) | Do PCA bands match indel genotype clusters? (independent) |
| **T9** | Ancestry jackknife (Cheat 6) | Does signal collapse without one Q-group? (founder linkage) |
| **T10** | Theta het prior (Cheat 12) | Do heterozygotes have elevated θ_P? (independent) |
| **T11** | Extended boundary suppression (Cheat 23) | Are there boundary-suppression effects? (H3) |

The 4 hypotheses:
- **H1**: Family/relatedness drives the pattern (artifact)
- **H2**: Broad inversion with internal hot core (real, simple)
- **H3**: Nested or composite inversion (real, complex)
- **H4**: Sub-haplotype/founder block within carriers (real, sub-structured)
- **H5**: Technical/informativeness artifact

## 2. Atlas surface — per-candidate verdict panel

A new section on the page-2 candidate-focus page (or a sub-tab):

```
Hypothesis test framework — Inversion {seq_num}
─────────────────────────────────────────────────────────────────
T1 Within/between band relatedness   [PASS] : multi-family bands
T2 Ancestry diversity per band       [PASS] : K=8 distributed across bands
T3 Kin-pruned signal retention       [PASS] : signal survives 81-pruned
T4 Inner/outer composition overlap   [PASS] : same composition (H2 simple)
T5 Anchor stability                   [PASS] : anchors stable across windows
T6 Regime change at boundary          [N/A]  : no inner boundary
T7 Within-carrier substructure        [WEAK] : sub-clusters present but small
T8 Clair3 indel concordance           [PASS] : OR=12.4, p=1e-8 (Cheat 9)
T9 Ancestry jackknife                 [ROBUST] : multi-family (max Δ=0.018)
T10 Theta het prior                   [PASS] : HET θ_P > HOM (Wilcoxon p=2e-5)
T11 Boundary suppression              [PASS] : Cheat 23 no extended effect

Decision tree verdict: H2 (simple inversion)
BH-corrected confidence: HIGH (q=4.3e-7)
Independent confirmations: 3/3 (T8, T9 robust, T10)
```

Each test row is **clickable** → opens a per-test detail popover
showing the actual statistic value and the test's specific output
(e.g., for T8: the contingency table; for T9: per-Q-group Δ values).

## 3. JSON layer schema

Atlas reads `hypothesis_tests_v1` JSON layer:

```json
{
  "schema": "hypothesis_tests_v1",
  "candidates": [
    {
      "id": "LG28_inv5",
      "tests": {
        "T1": { "verdict": "PASS", "statistic": 0.18, "p_value": 1.2e-3, "details": "..." },
        "T2": { "verdict": "PASS", "statistic": 6.8, "p_value": null, "details": "..." },
        ...
      },
      "decision_tree_verdict": "H2",
      "bh_corrected": {
        "q_value": 4.3e-7,
        "confidence": "HIGH",
        "n_independent_confirm": 3
      }
    }
  ]
}
```

Verdicts: `PASS` / `WEAK` / `FAIL` / `N/A` / `ROBUST` / `FRAGILE` / `UNKNOWN`.

## 4. Decision tree integration

The decision tree synthesis (per `c1a88579`) lives LANTA-side. Atlas
just renders the verdict. But atlas should surface the **decision
trail** — which steps led from T1 through T11 to the final verdict.

A small textual trail beneath the verdict:

```
Decision trail:
  STEP1 [T1]: Bands not driven by relatedness → continue
  STEP2 [T2]: Ancestry distributed across bands → continue
  STEP3 [T3]: Signal robust to kin pruning → not H1
  STEP4 [T4]: Inner = outer → not H3 nested
  STEP5 [T5]: Anchors stable → H2 simple confirmed
  → Verdict: H2 (simple inversion)
  → 3 independent confirmations (T8, T9, T10)
```

## 5. Implementation slices

### Slice 1 — JSON layer + verdict panel (~0.7 turn)
- Detect `hypothesis_tests_v1` layer
- Render the 11-test status grid on page 2
- Click → popover with test details

### Slice 2 — decision trail (~0.3 turn)
- Render the step-by-step decision trail
- Tooltip explanations for each hypothesis class

### Slice 3 — catalogue integration (~0.3 turn)
- Add `decision_tree_verdict` column to page-3 catalogue
- Filter by verdict (H1 / H2 / H3 / H4 / unresolved)

### Slice 4 — manuscript bundle inclusion (~0.2 turn)
- Per-candidate verdict in `results_boilerplate.md`
- Decision trail in supplementary

## 6. Open questions

1. **Single-test thresholds vs decision tree**: do we surface the raw
   test outputs or just the synthesized verdict? Both — verdict at the
   top, raw values in click-to-expand details.
2. **What about candidates missing some tests?** Some candidates can't
   run T6 if there's no inner boundary; T7 requires sufficient carriers.
   Render `[N/A]` cleanly without breaking the verdict synthesis.
3. **BH correction scope**: across all candidates per chromosome or
   genome-wide? LANTA-side decides; atlas just displays the output.

## 7. Tests

- Synthetic candidate with all 11 test outputs → all 11 rows render
  with correct verdict colors (green/yellow/red/grey for PASS/WEAK/FAIL/N/A).
- Decision trail renders steps in order without missing any.
- Catalogue filter by `verdict = H2` → only H2 candidates remain.

## 8. Dependencies

- LANTA-side `STEP_C01f_hypothesis_tests.R` (existing).
- BH correction logic (existing per `c1a88579`, `bf1625f9`).
- Engine B `region_stats_dispatcher.R` for T9 ancestry jackknife.

## 9. What this is NOT

- Not a re-implementation of the tests in the atlas. All compute
  remains LANTA-side. Atlas only renders.
- Not a hypothesis-discovery tool. The 4 hypotheses (H1–H4) are fixed;
  this surface displays which one applies per candidate.
- Not a substitute for human judgment on disputed cases. Where the
  decision tree returns "unresolved," the user reviews manually.

## 10. Cross-references

- Chat `c1a88579` — verdict synthesis design.
- Chat `bf1625f9` — manuscript text integration with test outputs.
- Chat `5d2ca953` — codebase consolidation, full T1-T11 description.
- `SPEC_manuscript_bundle_export.md` — verdicts feed boilerplate.
- `SPEC_recombinant_dosage_changepoint_detector.md` — T8 dependency
  (Clair3 dosage already used here).
