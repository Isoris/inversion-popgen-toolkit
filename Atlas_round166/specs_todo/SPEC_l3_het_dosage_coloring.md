# SPEC — Het (dosage) coloring on L3 contingency mini-PCAs

**Status**: Slice 1 **SHIPPED in turn 129**. Slice 2 (tracked-samples
always-on, per-sample-lines extension) deferred. Spec stays in
`specs_todo/` until Slice 2 lands or is formally cancelled.

**Slice 1 shipped contents** (turn 129, +330 lines in `Inversion_atlas.html`):
- `_hetRateColor(rate)` with `_HET_RAMP` constant (RdBu, anchored at 0.5)
- `_computeHetRateForL2(l2idx)` — sync, reads dosage chunk LRU via
  `window.popgenDosage.getCachedChunk`, memoized in `state.__hetRateCache`
- `_invalidateHetRateCache()`, `_getHetRateCache()`,
  `_refreshL3HetToggleAvailability()`
- `withAlpha` extended to handle `rgb(...)` inputs
- `drawMiniPCA` hook overriding non-tracked dot fill when toggle on
- Color-bar legend on focal mini-PCA, bottom-left
- L3 toolbar `#l3HetToggle` checkbox + localStorage persistence
- Cache invalidation + toggle re-availability on `dosage_chunks` layer load
- 58 tests in `tests/test_turn129_l3_het_coloring.js`

---

**Trigger**: Quentin's request (turn 128d):
> "Now that we have dosage, can we get Het (dosage) coloring implemented or
> all panels of l3 contingency? if activated? Add a button for it. but
> always active for tracked sample panel unless toggled off? Or how did we
> spec it last time?"

**Prior context (the "spec it last time" question)**: the existing
per-sample-lines color-mode picker (8 modes — kmeans, dosage, ghsl, het,
θπ, F_ROH, family, ⚠ confounder) was scaffolded in v3.99 turn 14e+. State
slot `state.linesColorMode` exists, the dropdown UI exists, layer-
availability validation exists. But the actual color resolvers
`_resolveSampleColorByMode` and `_resolveSampleScopeColor` are **stubs**
returning `null` for every mode except `family` (wired in v4 turn 108).
The dosage_chunks layer was the missing input — it now ships via the S6
bridge / standalone viewer.

---

## 1. The diagnostic question

When a candidate's K-means clustering produces three bands (g0/g1/g2),
the user wants to verify whether the middle band is genuinely
heterozygous in genotype-call terms, not just sitting between two
homozygous PC clouds for unrelated reasons (admixture, family structure,
confounder gradients).

Het coloring on the L3 mini-PCAs answers this in one glance:

- **Clean three-band karyotype**: g0 dots paint blue/cold (low het, ~0%
  expected for HOMO_REF), g1 dots paint warm-yellow (high het, ~50%
  expected for HET), g2 dots paint blue/cold (low het, ~0% for HOMO_INV).
  Strong PC1 gradient = strong inversion signal.

- **Hidden subdivision**: g0 dots paint heterogeneous colors → suggests
  the band is mixing two distinct populations, OR the K=3 split is at the
  wrong place. Worth investigating with K=6 / sub-clustering.

- **Misaligned bands**: g1 (the supposed HET band by PC1 ordering) doesn't
  paint warmer than g0/g2 → K-means assignment doesn't reflect biology.
  The candidate may be a false positive, or the K=3 partition is wrong.

This is the single most useful "is my karyotype call real?" diagnostic —
pure dosage truth without going through any clustering algorithm.

## 2. Where het coloring renders

Three rendering surfaces, three policies:

### 2.1 L3 contingency mini-PCAs (`drawMiniPCA`)
- **Default**: K-means cluster colors (existing behaviour)
- **When toggled**: every dot painted by per-sample het-rate at the L2's
  windows, scrubbed through a divergent ramp (cold→warm→cold). Tier 1
  surface — the most concrete payoff for the user.

### 2.2 Tracked-samples scatter (the right-aside list of dot+CGA labels)
- **Default**: tracked-color (current sample-identity palette)
- **Special**: **always-on** het halo for tracked samples even when the
  global toggle is off, with an explicit "off for tracked" override.
  Quentin's exact words: *"always active for tracked sample panel unless
  toggled off"*. Rationale: the tracked list is your watchlist — you
  always want the het signal on those specific dots.

### 2.3 Page-1 main `#pcaCanvas` (`drawPCA`)
- **Out of scope for this slice**. The main PCA already supports a
  `state.colorMode` switcher (`cluster | family | ancestry | manual |
  none | q_ancestry | cluster_dosage | cluster_theta_pi | cluster_ghsl`).
  Adding `'het'` to that list is straightforward but gets entangled with
  the lines-panel `state.linesColorMode` picker — the two state slots
  should converge or coexist deliberately. Defer.

### 2.4 Per-sample-lines panel (`drawLinesPanel`)
- **Already scaffolded** with `state.linesColorMode === 'het'` and the
  resolver stub. Wiring up the renderer body completes that scaffolding.
  Out of scope for Slice 1; included in Slice 2.

## 3. The het-rate computation

Per-sample, per-L2 (= scope of one mini-PCA):

```
hetRate(si, L2_idx) =
    count(dosage_chunks where windows ∈ L2 AND sample==si AND value==1)
  / count(dosage_chunks where windows ∈ L2 AND sample==si AND value != NA)
```

- `dosage_chunks` carries `markers` (per-window) × `samples` matrices
  with values in `{0, 1, 2, NA}`.
- A "window ∈ L2" maps via `state.windowToL2`. The dosage chunk's
  `markers` are window-indexed, so we filter to `mw` such that
  `state.windowToL2[mw] === l2idx`.
- Het rate is in `[0, 1]`. NA (-1) entries are excluded from numerator
  AND denominator.
- If denominator is 0 (no calls for this sample in this L2), color is
  `--ink-dimmer` (greyed) — the sample is genuinely uncalled here.

## 4. Color ramp

Divergent, anchored at 0.5 (the expected HET frequency under
Hardy-Weinberg equilibrium for a balanced diallelic locus):

| Het rate | Color | Semantic |
|---|---|---|
| 0.00 | `#2166AC` (cold blue) | All-homo (REF or INV) |
| 0.25 | midway blue→neutral | Sub-het |
| 0.50 | `#F7F7F7` (neutral) | Expected HET |
| 0.75 | midway neutral→red | Excess homo on the other side |
| 1.00 | `#B2182B` (warm red) | Pure-het (rare; data anomaly?) |

This is the **same FIG_C08 / ColorBrewer RdBu ramp** already used by the
dosage heatmap renderer, so the visual vocabulary is consistent across
the atlas.

Edge case: most populations have HET rates clustering around 0.5 within
the HET band and below 0.1 within the HOMO bands. The ramp's neutral
midpoint at 0.5 means HET-band samples appear pale (close to white),
HOMO-band samples appear blue (cold). The user reads "blue = homo, white
= het" — clean visual mapping.

## 5. UI controls

### 5.1 L3 toolbar het-toggle (Slice 1)

Add a `[ ] het` checkbox to the L3 toolbar (between the existing layout +
K-mode controls). State: `state.l3HetColoring` (boolean, default false).
Persists to `localStorage.pca_scrubber_v3.l3HetColoring`.

When ON:
- Every L3 mini-PCA paints dots by het-rate (instead of K-cluster).
- The L3 contingency tables stay K-cluster-colored (the tables' job is
  to show K-cluster correspondence, NOT het distribution).
- A small color-bar legend appears under the focal mini-PCA showing the
  ramp (~28×6px swatch with `0` and `1` end-tick labels).

When OFF (default):
- Existing K-cluster coloring (no behavioural change).

Tracked samples in the L3 mini-PCAs: their **halo** (the outer ring
around tracked dots) is K-cluster-colored regardless, so the user can
see "this sample is tracked AND it sits in band g0 BUT its het rate
suggests it might really be g1". The dot fill flips to het-rate when
the toggle is on.

### 5.2 Tracked-samples always-on subtoggle (Slice 2 — deferred)

A second checkbox in the tracked-samples aside header: `[x] het halos
on tracked` (default: ON). When the global L3 toggle is off, the halo
is K-cluster-colored as today. When ON, the halo gains an inner ring
in the het-rate color so the user sees the het signal even on a
K-cluster-colored mini-PCA.

State: `state.trackedSamplesHetHalo` (boolean, default true).

**Why deferred**: requires touching `drawMiniPCA`'s tracked-dot rendering
path, which has 3 visual modes today (shared / unique / dual rings).
The interaction with "dual" mode in particular needs design — there are
already two rings (focal + own-pane K-cluster) and we'd be adding a third.
Worth a separate turn with a screenshot from Quentin showing the desired
look.

## 6. Performance

Computing per-sample het rate per mini-PCA is `n_samples × n_windows_in_L2`
operations. Typical L2 = 5–50 windows; cohort = 226 samples. Per mini-
PCA: ~2k–11k array reads. The L3 panel renders up to 5 mini-PCAs at
once (focal + ±2 carousel). Total per redraw: ~10k–55k reads. Negligible.

Caching: memoize `hetRate(si, l2idx)` keyed by `(l2idx, dosageChunkVersion)`.
Invalidate when:
- A new dosage_chunks layer loads (via `state.layersPresent` add).
- The user uses the dosage-heatmap chunk-cache reset action.

Cache lives at `state.__hetRateCache`, a `Map<l2idx, Float32Array>`.

## 7. Layer requirements

- `dosage_chunks` (MUST be present). The het toggle is disabled with a
  tooltip "Needs dosage_chunks JSON layer" when missing.
- `state.windowToL2` (always present, built at JSON load).

No new JSON layer is needed. The same `dosage_chunks` layer used by the
dosage heatmap and S6 bridge is the source for het rates.

## 8. Implementation checklist (Slice 1, this turn)

- [x] Spec written.
- [x] Add `state.l3HetColoring: false` (boolean) + localStorage restore.
- [x] Add `[ ] het` checkbox to the L3 toolbar; wire to state slot.
- [x] New helper `_computeHetRateForL2(l2idx)` returns
      `Float32Array[n_samples]` of het rates ∈ [0,1] | NaN-for-missing.
      Memoizes via `state.__hetRateCache: Map`.
- [x] New helper `_hetRateColor(rate)` returns CSS color from the ramp.
      Returns `--ink-dimmer` when rate is NaN.
- [x] In `drawMiniPCA`: when `state.l3HetColoring`, override the dot fill
      with `_hetRateColor(hetRate[si])`. Keep the halo K-cluster-colored.
- [x] Color-bar legend under each mini-PCA when toggle is on.
- [x] Cache invalidation hook on `dosage_chunks` load.
- [x] Tests: pure helper for `_computeHetRateForL2` against a synthetic
      chunk. Source-level tests for the toggle UI + state restore.

## 9. Implementation checklist (Slice 2, deferred)

- [ ] `state.trackedSamplesHetHalo` (default true) + sub-checkbox in
      the tracked-samples aside header.
- [ ] Tracked-dot halo: when slot true AND `dosage_chunks` present,
      add an inner ring colored by het-rate even when the global L3
      toggle is off.
- [ ] Reconcile with `state.l3ColorMode === 'dual'` (already paints two
      rings) — design a tri-ring layout or grant het-halo priority.
- [ ] Wire `_resolveSampleScopeColor('het', si)` for callers that want
      het as a primary color.
- [ ] Wire the per-sample-lines `linesColorMode === 'het'` resolver
      (the existing scaffolded but stubbed mode).

## 10. Tests

Source-level + behavioural for Slice 1:

- `_computeHetRateForL2` returns an `n_samples`-length Float32Array.
- For a synthetic chunk where sample 0 is `[0,1,1,0]` across L2 windows,
  het rate is 0.5.
- For a sample with all NA, het rate is NaN.
- Cache hits avoid re-scan: second call with same `(l2idx, version)`
  returns the same Float32Array reference.
- `_hetRateColor(0)` ≈ cold-blue start; `_hetRateColor(0.5)` ≈ neutral;
  `_hetRateColor(1)` ≈ warm-red end.
- Toggle UI: source-level checks for the checkbox + handler + state
  slot + localStorage key.
- `drawMiniPCA`: source-level check that the het branch references
  `state.l3HetColoring` and `_hetRateColor`.

## 11. What this is NOT

- **Not a replacement for K-means clustering.** The K-cluster halos
  remain on tracked dots. The het coloring is an *overlay* diagnostic.
- **Not a global het track.** Page 6 (popstats) already has Hobs / Hexp
  tracks aggregated over the cohort. Het coloring here is per-sample,
  per-candidate — a different scope.
- **Not a substitute for dosage_chunks heatmap.** The heatmap shows
  per-marker × per-sample dosage in a 2D grid; het coloring averages
  across markers within the L2 to give one number per sample. Both
  views complement each other.
- **Not page-1 main PCA coloring.** That's `state.colorMode`; this is a
  separate `state.l3HetColoring`. The two slots are intentionally
  independent so the main PCA stays K-cluster while the L3 minis can
  switch to het without affecting the focal scatter.
