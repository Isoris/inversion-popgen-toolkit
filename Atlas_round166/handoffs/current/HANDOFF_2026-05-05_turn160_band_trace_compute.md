# HANDOFF — turn 160 — band-trace inheritance map (compute layer)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (74,417 lines, +313 LOC from 74,104)
**Working dir**: `/home/claude/work/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

**Closes** (compute half): the "fan-cosegregate-fan" intuition that's been
floating across multiple chats since `f74cf5d4-...` ("Integrating atlas
with phase 6 and 7"). Spec: `SPEC_distant_band_concordance_fish_trajectory.md`
Slice 4. Slice 1+2 already shipped turn 130 (lineage compute + per-sample-line
coloring). This turn ships the band-trace primitive — the layer that makes
"long-range guessing" of inversions possible by tracking a fish-set's band
occupancy across all L2s on a chromosome.

**Picked up from**: post-turn-159 working tree.

---

## 0. The intuition (paraphrased back)

Pick a fish-set — say, the union of bands {1, 2} at L2 #5 (some
reference position around 500 kb). Walk those same fish across every
other L2 on the chromosome. Three things can happen at each distant L2:

- **Co-segregating**: the fish-set is concentrated in 1–2 bands.
  Their founder haplotypes co-segregate at this distant locus too.
- **Partial**: spread across 2–3 of K bands but with a clear dominant.
- **Fanned**: spread roughly uniformly across all K bands — no
  inheritance signal, this distant locus is independent of the
  reference.

The pattern along genomic position is the **fan → co-seg → fan**
profile. Co-seg runs are inversion footprints. The boundaries of the
runs are inversion breakpoint candidates. Two distant co-seg runs
sharing the same fish-set → two inversions whose alleles
co-segregate in the founder pool → potential overlapping inversion
systems.

This is "long-range guessing" because the trace exposes inversions
the user never explicitly identified — the reference fish-set was
picked at one position, but the regime map reveals every chromosomal
region where those founder haplotypes happen to co-segregate.

---

## 1. What this turn ships

Three pure compute helpers, all on `window`:

### `_bandTraceShannonEntropy(fractions, K)` → number ∈ [0, 1]

Normalized Shannon entropy. 0 = degenerate (all in one band),
1 = uniform across K. Handles unnormalized counts (function divides
by total). K=1 returns 0 (no normalization possible). Empty input
returns 0.

### `_bandTraceForFishSet(fishSet, opts)` → trace object

Walks Hungarian-projected labels across every L2 in every chain on the
chromosome. At each L2, counts how many of the selected fish land in
each band, computes entropy, classifies into one of five regimes.

```js
{
  n_fish_selected, n_chains, n_total_L2, K,
  per_l2: [
    {
      l2_idx, chain_idx, chain_position,
      n_valid,                     // selected fish with non-(-1) label
      band_counts,                 // Int32Array[K]
      band_fractions,              // Float32Array[K], sums to 1
      entropy,                     // normalized Shannon, [0, 1]
      dominant_band,               // argmax(band_fractions), -1 if no_valid
      dominant_fraction,           // max(band_fractions)
      regime,                      // 'co_seg' | 'partial' | 'fanned' | 'sparse' | 'no_valid'
    },
    ...
  ],
}
```

Regime thresholds (defaults):
- `co_seg`:  entropy ≤ 0.40
- `fanned`:  entropy ≥ 0.85
- `partial`: in between
- `sparse`:  fewer than 3 valid fish at this L2 (overrides regime)
- `no_valid`: no fish in the set have a label here

Calibration on real LG28 data may shift these. All three are
overridable via `opts.coseg_max`, `opts.fanned_min`, `opts.min_valid`.

Reuses `_hungarianChainProjection` (turn 130). Does NOT touch
`state.lineageResult` — band-trace is a per-fish-SET trace, lineage
is a per-fish label. Orthogonal axes that share the same projection.

### `_bandTraceRegimeRuns(trace, opts)` → array of run objects

Walks `trace.per_l2` and emits runs of consecutive co-seg (and
optionally partial) L2s where the dominant band is stable. Runs do
not span chain breaks even when the regime is identical on both
sides — chain breaks mean Hungarian alignment couldn't stitch the
projection, so band IDs aren't comparable.

```js
[
  {
    start_l2_idx, end_l2_idx,
    start_chain_position, end_chain_position,
    chain_idx, n_L2,
    dominant_band,                 // mode of dominant_band across run
    mean_dominant_fraction,
    mean_entropy,
    n_co_seg, n_partial,
  }, ...
]
```

Filters:
- `min_run_length` (default 2): single-L2 blips dropped
- `allow_partial` (default true): include partial regimes in runs

---

## 2. What this turn does NOT do

Everything UI. No surface in the atlas. The layer is callable from
the console / programmatically via `window._bandTraceForFishSet`,
which is exactly what was in scope for "compute half ships first."

UI scope for next turn:
- Floating trace strip on the per-sample-lines panel showing one
  column per L2 with stacked bar of `band_fractions`, colored by
  regime
- Band-set picker: click bands at any L2 → populate
  `state._bandTraceSelection = { ref_l2, band_set }` → re-render strip
- Detected co-seg runs marked with thick top brackets =
  "candidate inversion footprint"
- Combinatorial enumeration: button that loops through all subsets
  of bands at the focal L2 (size 1..K-1, ≤ 2^K - 2 = 6 subsets at
  K=3) and ranks (subset, distant L2 region) pairs by run strength.
  Output is the long-range inversion candidate list.

---

## 3. Files touched

```
Inversion_atlas.html                                +313 LOC
  - _bandTraceShannonEntropy                        new function
  - _bandTraceForFishSet                            new function
  - _bandTraceRegimeRuns                            new function
  - 4 threshold constants + window exports

tests/test_turn160_band_trace_compute.js            new (58 assertions)
```

No existing code modified. No tests left in a broken state.

---

## 4. Test results

**Single test**: 58 / 0 across 9 sections:
1. Source-pattern checks (14 assertions)
2. `_bandTraceShannonEntropy` basics — degenerate, uniform, partial,
   empty, K=1, unnormalized counts (6)
3. `_bandTraceForFishSet` on synthetic two-inversion planted layout —
   30 fish, 10 L2s, K=3, inversion A at L2s 2-4 + inversion B at L2s
   6-8. Trace recovers both regimes with correct dominant bands (12)
4. Regime classification: 70/30 split → `partial`, sparse → `sparse`,
   all-invalid → `no_valid` (4)
5. Empty / null inputs → null (3)
6. `_bandTraceRegimeRuns` extracts both planted regimes with correct
   start/end L2s and dominant bands (8)
7. Chain breaks split runs even when regime is identical (3)
8. `min_run_length` and `allow_partial` filters work (4)
9. Regression checks for turns 130 / 156 / 159 (4)

**Full sweep**: **3094 / 0** (was 3036 / 0 at turn 159 close).
Zero regressions. JS syntax clean. HTML parser 0 errors.

---

## 5. What Quentin can do today (with the compute layer alone)

From the browser console after loading a chromosome:

```js
// Pick any fish-set: a candidate's band, manual group, lasso selection.
const fishSet = state.candidateList[0].locked_labels
  .map((b, si) => b === 0 ? si : -1)
  .filter(si => si >= 0);

const trace = _bandTraceForFishSet(fishSet);
// trace.per_l2[i].regime gives you the per-L2 regime call
// trace.per_l2[i].entropy is the continuous metric

const runs = _bandTraceRegimeRuns(trace);
// runs is the list of inversion-footprint candidates
console.table(runs);
```

The trace data is enough to build a CSV/TSV for spot-checking against
the candidate catalogue. Future-Claude can wrap that into a "📊 export
band-trace TSV" button when the UI lands.

---

## 6. What's NEXT

Next session continues SPEC_distant_band_concordance_fish_trajectory.md
Slice 4 with the UI half:

1. **Trace strip renderer** (~150 LOC) — Canvas above the per-sample-lines
   panel, one column per L2, stacked bar of band_fractions colored by
   regime + a top bracket per detected run.
2. **Band-set picker** (~100 LOC) — click any band at any L2 in the
   strip or main panel to add/remove from `state._bandTraceSelection`.
3. **Combinatorial enumeration button** (~80 LOC) — loops through all
   subsets of bands at a chosen reference L2 (capped at 2^K - 2),
   produces ranked list of (subset, run) pairs by `mean_dominant_fraction
   * n_L2`. This is the "auto-find overlapping inversions" button.

After UI lands, the spec moves from `specs_todo/` to `specs_done/`.

---

## 7. Honest framing

**What's solid:**
- Pure functions, no state mutation, fully sandbox-testable.
- Reuses existing `_hungarianChainProjection` — same chain-break
  semantics as the lineage compute, so projection behavior is
  consistent across both surfaces.
- Synthetic planted-inversion test recovers both regimes exactly,
  including the case where the same fish-set co-segregates in
  *different* dominant bands at the two inversions (band 0 at A,
  band 1 at B) — the regime call is band-independent (it's about
  concentration), the dominant_band annotates which band.
- Threshold constants are global and tunable; test confirms 70/30
  → partial, with the entropy in the documented [0.40, 0.85] window.

**What's risky:**
- Thresholds are guesses calibrated against synthetic data. Real
  LG28 with 226 fish and ~30 L2s will reveal whether `co_seg ≤ 0.40`
  is too strict (real data is noisier) or too lax (background fanned
  noise leaks into co_seg). First real-data run will probably need
  adjustment.
- Chain-break handling: if a chromosome has a Hungarian-fragile
  region in the middle, the projection breaks into multiple chains
  and the band-trace becomes per-chain. Runs cannot span chain breaks
  by design. This is correct (chain-break means band IDs aren't
  comparable) but it means "this inversion spans the chain break"
  is not detectable from this surface alone — needs the secondary
  marker-level evidence the spec calls out.
- The "sparse" regime is set at min_valid=3 hardcoded. If you
  routinely track 5-fish manual groups, that minimum may eat real
  signal. Adjust via `opts.min_valid` when calling, or change the
  global default.

**What's queued:**
- UI half of Slice 4 (next turn).
- Combinatorial enumeration button (next turn).
- After UI: spec moves to `specs_done/`.

End of handoff.
