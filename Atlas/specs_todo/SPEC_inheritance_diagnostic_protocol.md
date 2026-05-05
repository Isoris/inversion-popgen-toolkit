# SPEC — Inheritance pipeline diagnostic protocol

**Status**: drafted turn 130. **This is a DIAGNOSTIC PROTOCOL, not a code
spec.** It tells future-Claude (or current-Quentin) what to check, in
what order, before committing to L2-sweep / sliding-window / unification
turns.

**Trigger**: turn 130 menu option 4 — *"first show me what the existing
candidate-based inheritance pipeline actually produces on real data
(now that the bridge is in place)."*

**Why a doc, not a code turn**: the diagnostic itself doesn't need new
code. It needs Quentin's eyes on the existing UI on a real chromosome.
The doc records what to look at and what would imply what next.

---

## 1. Setup (one-time, ~5 min)

1. Load the latest atlas with the turn-129 registry bridge.
2. Open a chromosome with at least 3 known candidates — LG28 (the
   primary validation chromosome) is the ideal.
3. Confirm the candidate strip shows your saved candidates. If it
   doesn't, the bridge isn't doing its job and this whole stack is
   moot.

## 2. Test 1 — does the inheritance pipeline now produce I·g pills?

Open the per-sample-lines panel (page 1, PC1 sub-panel). Look at the
**top of the plot** for small "I1·2g", "I2·3g" style pills above each
candidate's bp range.

| Observation | Interpretation | Next step |
|---|---|---|
| Pills visible, labels make sense | Bridge works, pipeline fires | Go to Test 2 |
| No pills at all, but candidates ARE in the strip | Pipeline ran but produced no output (n_groups < 1 or compute bailed) | Open browser console, check `state.inheritanceResult` — likely null |
| No pills AND no candidates in the strip | Bridge isn't bridging | Check `state.candidates` vs `state.candidateList` in console |
| Pills visible but labels look wrong (every candidate shows "1g") | Compute ran but Jaccard threshold too aggressive — everything collapsed to one group | Lower `_IGC_DEFAULT_COSINE_DIST_THRESHOLD` or change to a per-chromosome threshold |
| Pills visible but every candidate shows "Kg" (n_groups == K of candidate) | Threshold too lax — nothing merged | Raise threshold |

## 3. Test 2 — does the tooltip show useful group composition?

Hover an I·g pill. A tooltip should pop up with:
- Candidate id + bp span
- Per-band fish counts ("g0: 32 fish, g1: 78 fish, g2: 116 fish")
- Per-group composition ("group 0: spans bands g1 of I1 + g2 of I2")

| Observation | Interpretation |
|---|---|
| Tooltip useful | The data structure is there. Move on. |
| Tooltip empty / "no data" | `_inhTooltipBuildHtml` is reading the wrong field |
| Tooltip shows numbers but they don't match the strip | Stale-cache guard not firing — invalidate via `invalidateInheritanceCache()` and re-render |

## 4. Test 3 — does the candidate-page haplotype annotation panel show inh suggestions?

Open the candidate-focus page (page 2). For the active candidate, look at
the per-band haplotype annotation rows. Each row should show:
- Band swatch + index
- Vocab dropdown (REF/HET/INV or H1/H2/H3 etc.)
- Auto-classifier suggestion + confidence chip
- **`inh{id}` chip** linking the band to its inheritance group

| Observation | Interpretation |
|---|---|
| `inh{id}` chips visible | `_inheritanceSuggestionsForCandidate` works end-to-end. Done. |
| No chips, but pills work | The suggestions reader is reading from a different slot than the pills draw from. Check `state.inheritanceResult.cut.group_id_per_band` indexing. |

## 5. Test 4 — do groups make biological sense?

Pick 2 candidates on the chromosome. Look at the I·g pills:

- If both pills say "1g", the pipeline thinks both candidates segregate
  the same way (likely linked, or the threshold's wrong).
- If they say "2g" and "3g", they segregate differently.
- Use the matrix UI (`openInheritanceMatrix()` from the console) to see
  the Cramér's V between every candidate pair.

| Observation | Interpretation |
|---|---|
| High Cramér's V (≥ 0.7) between expected-linked candidates | Algorithm is finding biology |
| Low Cramér's V everywhere on a chromosome where you expect linkage | Either the K-means is producing noisy bands OR the candidates are genuinely independent OR the threshold needs work |
| High Cramér's V on candidates that biologically should be independent | The K-means is being driven by family structure, not karyotype — confounder |

## 6. What the diagnostic decides

Three branches:

### Branch A — pipeline produces useful output on LG28

- **Don't unify**, the candidate-based version works.
- Ship L2-sweep + sliding-window + fish-trajectory as **complementary**
  layers (they give different info than the candidate-based one does).
- Combinatorial selection + auto-promote on top.

### Branch B — pipeline produces output but it's noisy / wrong

- Audit `inheritanceGroupClustering` parameters
  (`_IGC_DEFAULT_COSINE_DIST_THRESHOLD`, `_IGC_MIN_BAND_SIZE`).
- Calibrate against LG28's known ~2.89 Mb sub-telomeric inversion.
- Then ship the four other specs.

### Branch C — pipeline produces nothing

- The bridge isn't actually fixing the empty-input problem in
  practice (despite the tests passing).
- Debug end-to-end: console.log `_gatherActiveCandidatesForInheritance()`,
  inspect items array, find why it's still empty.
- Fix that before any other work.

## 7. What to record

After running through Tests 1–4, write back:

- Which test failed first (or "all passed")
- Screenshot of the per-sample-lines panel with pills (or absence)
- Console log of `state.inheritanceResult` (full JSON)
- Console log of `state.candidateList.length` and
  `Object.keys(state.candidates).length` (these should match)

That data tells us whether to ship the four other specs as
complementary layers (Branch A), or fix the existing one first
(Branches B/C).

## 8. What this is NOT

- **Not a code change.** No new file is created by running this
  protocol.
- **Not a substitute for the four other specs.** The protocol just
  helps decide implementation order + threshold calibration.
