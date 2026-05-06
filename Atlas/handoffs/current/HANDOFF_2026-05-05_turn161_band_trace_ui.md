# HANDOFF — turn 161 — band-trace UI strip

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (74,793 lines, +376 LOC from 74,417)
**Working dir**: `/home/claude/work/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

**Closes** the UI half of `SPEC_distant_band_concordance_fish_trajectory.md`
Slice 4. Compute layer shipped turn 160; this turn surfaces it as a
strip on the per-sample-lines panel.

**Picked up from**: post-turn-160 working tree.

---

## 0. Manuscript framing (the axiom)

Quentin's directive: 30 family hubs across 226 fish ⇒ effectively
random mixing of founder backgrounds. So when a fish-set picked at
position X co-segregates 15 Mb away into the same 3-of-6 bands, that's
not population structure — it's **physical co-inheritance** of a
haplotype block that didn't recombine. Whether it's one inversion,
two overlapping inversions, or some other recombination suppressor
is interpretation; the **observation** is real and stands.

The UI is built to that bar:
- **Strip is observation-only.** No "this is an inversion" markers.
- **Default OFF.** User opts in.
- **Stacked-bar visualization** of band-fraction at each L2, regime
  color (green=co_seg, amber=partial, grey=fanned) as a thin top
  stripe. The user reads the pattern.
- **Manuscript-friendly framing**: report "we observed N
  co-segregating haplotype blocks in a randomly-mixed hatchery
  cohort using long-range local-PCA signals." Don't push further.

---

## 1. What this turn ships

### State slots (4 new)

```js
state.bandTraceFishSet  : int[] | null   // sample indices to track
state.bandTraceOn       : bool            // strip visibility (default false)
state.bandTraceCache    : trace | null   // last-computed trace
state.bandTraceCacheKey : string | null  // chrom + fish-set fp + K + nE
```

All cleared on chrom switch (hooked into `applyData` alongside
existing `invalidateLineageCache()` call).

### Setter helpers + persistence

- `setBandTraceFishSet(arr)` — dedupes + integer-coerces + filters
  negatives + persists to `localStorage` (`inversion_atlas.bandTraceFishSet`)
  + invalidates cache + repaints. `[]` or `null` clears.
- `setBandTraceOn(bool)` — coerces + persists (`inversion_atlas.bandTraceOn`)
  + repaints.
- `_bandTraceClearCache()` — drops cache only, leaves fish-set.

### Auto-fill from focal candidate

`_bandTraceFromFocalCandidate(opts)` — picks fish from the
`state.candidate.locked_labels` array. By default takes the largest
band; `opts.bandIdx` overrides. Returns the fish-set array on
success, `null` if no candidate / no labels / no non-empty band.
Side effect: calls `setBandTraceFishSet` so localStorage + cache +
repaint all happen automatically.

### The strip drawer

`_drawBandTraceStrip(ctx, pad, plotW, plotH, mbMin, mbMax)` — Canvas
drawer. Skips silently when (toggle off) OR (no fish-set) OR (no
data). Otherwise per L2 in view range:
- Top 1 px stripe colored by regime (`co_seg`/`partial`/`fanned`/
  `sparse`/`no_valid` → 5 distinct colors)
- Stacked bars below showing each band's fraction, colored by the
  same band-color palette as the karyotype tab (`_gpKaryoColor`)

Strip dimensions: 7 px tall, sits at `pad.t - 21` (above the existing
lineage strip at `pad.t - 13`, height 5 — they coexist when both on).

### Cache layer

`_bandTraceCacheKey(chrom, fishSet, K, n_envelopes)` — deterministic
key. Order-insensitive in fish-set. Sensitive to chrom, K, and L2
envelope count. `_bandTraceGetOrCompute()` short-circuits when the
key matches; otherwise calls `_bandTraceForFishSet` (turn 160) and
stashes the result.

### Lines header UI

Two new elements inserted next to the existing `lineage` toggle:

```
[ ] band trace    🔍 trace
```

- Checkbox `linesBandTraceToggle`: default OFF, restores from
  `localStorage` if previously enabled. Toggle wires to `setBandTraceOn`.
- Button `linesBandTraceFromCandBtn` (🔍 trace): runs
  `_bandTraceFromFocalCandidate()`. On success: auto-enables the
  strip (sets toggle on if not already). On failure: alerts the user
  that they need to focus a candidate first.

### `drawLinesPanel` integration

Single-line addition right after the existing `_drawLineageStrip`
call: `_drawBandTraceStrip(ctx, pad, plotW, plotH, mbMin, mbMax)`,
guarded by `try/catch` so a drawer failure cannot break the rest of
the lines panel paint.

---

## 2. What this turn does NOT do

Deliberate scope discipline. None of the following are in:

- **Bracket overlays for detected runs.** `_bandTraceRegimeRuns`
  exists from turn 160 and the strip could draw a bracket per
  detected run. Skipped because that crosses into "this is an
  inversion" interpretation. Quentin's manuscript framing
  ("don't interpret too much") rules this out for now. Easy to
  add in a later turn if desired — one extra loop in the drawer.

- **Combinatorial enumeration button.** Looping through all
  subsets of bands at a focal L2 to rank long-range co-seg pairs.
  Listed as turn-160 followup; deferred to a separate turn so this
  one stays focused on "make the observation legible."

- **Hover tooltip.** No per-L2 tooltip yet (entropy, dominant
  band/fraction, regime). The data is in `state.bandTraceCache`
  for anyone inspecting from the console; UI tooltip is a small
  follow-up.

- **TSV export.** No "📊 export band-trace TSV" button. Same
  rationale — small follow-up.

- **Chain-break visualization in the strip.** When the Hungarian
  projection breaks into multiple chains, all L2s still get
  painted but chain identity isn't surfaced. Could add a thin tick
  between chain-break neighbors. Deferred.

---

## 3. Files touched

```
Inversion_atlas.html                                +376 LOC
  - 4 state slots (bandTraceFishSet, bandTraceOn, bandTraceCache,
                   bandTraceCacheKey)
  - _BTRACE_ON_LS_KEY, _BTRACE_FISH_SET_LS_KEY constants
  - _BTRACE_REGIME_COLOR palette
  - _bandTraceCacheKey(chrom, fishSet, K, n_envelopes)
  - _bandTraceGetOrCompute()
  - setBandTraceFishSet(arr)
  - setBandTraceOn(bool)
  - _bandTraceFromFocalCandidate(opts)
  - _drawBandTraceStrip(ctx, pad, plotW, plotH, mbMin, mbMax)
  - _bandTraceClearCache()
  - lines header: linesBandTraceToggle + linesBandTraceFromCandBtn
  - applyData: chrom-switch cache + fish-set clear
  - drawLinesPanel: invokes _drawBandTraceStrip
  - 6 window exports

tests/test_turn161_band_trace_ui.js                 new (77 assertions)
```

No existing functions modified. No tests left in a broken state.

---

## 4. Test results

**Single test**: 77 / 0 across 12 sections:
1. Source-pattern checks — state slots, exports, DOM ids, default-OFF (24)
2. Cache key determinism — order-insensitivity, sensitivity to chrom/K/nE (5)
3. `setBandTraceFishSet` — dedupe, persist, clear, invalidation (9)
4. `setBandTraceOn` — round-trip, boolean coercion (6)
5. `_bandTraceFromFocalCandidate` largest-band default (3)
6. `bandIdx` override + out-of-range fallback (3)
7. No-candidate / no-labels / all-invalid → null (3)
8. Drawer skip when toggle off + no fish-set (2)
9. Drawer paints expected primitives + cache populated (5)
10. Cache key invalidates on chrom switch (1)
11. `_bandTraceClearCache` (2)
12. Regression — turns 130/156/159/160 still wired (6)

**Full sweep**: **3171 / 0** (was 3094 / 0 at turn 160 close).
Zero regressions. JS syntax clean. HTML parser 0 errors.

The 10 environment-broken pre-turn-128 tests (missing fixture files)
remain unchanged.

---

## 5. What Quentin should exercise

The UX flow:

1. Load a chrom (e.g. LG28).
2. Promote/lock a candidate on page 1 (the lock-and-promote workflow).
3. In the per-sample-lines header, find the new `[ ] band trace 🔍 trace`
   pair, just to the right of the existing `[v] lineage` toggle.
4. Click **🔍 trace**. The strip auto-enables and paints. The fish
   from the focal candidate's largest band are now traced across
   every L2 on the chromosome.
5. Read the regime stripe (top 1 px of the strip):
   - **Green stretches** = those fish co-segregate at this position
   - **Amber stretches** = partial co-segregation
   - **Grey stretches** = fanned (no signal)
6. Read the stacked bars (the bottom 6 px): each color is one of K
   bands; the height proportions show how the traced fish distribute
   across bands at that L2.

Specific patterns to look for on real LG28:
- The strip should be **green at the candidate's own footprint**
  (sanity check — the fish that defined the band must co-segregate
  there).
- Every other green stretch = a long-range co-segregating haplotype
  block. **That's the observation.** It might be an inversion, or
  another suppressor, or co-inherited founder block — the strip
  doesn't claim which.
- If a clean K=3 candidate produces a strip that's almost entirely
  green: that fish-set is a deeply-conserved lineage (interesting,
  though not necessarily an inversion).
- If it's almost entirely fanned: the candidate's fish-set is
  randomly mixed at the genome scale (suggests the candidate's
  signal is local — maybe a small inversion or a local LD
  hotspot).

To trace a **different band** of the same candidate, drop into the
console: `_bandTraceFromFocalCandidate({bandIdx: 1})` (band 1 instead
of largest). Drop a fully custom fish-set: `setBandTraceFishSet([3, 7,
14, 22, ...])`.

---

## 6. What's NEXT (relative to the broader queue)

1. **Brackets for detected runs** (~30 LOC) — opt-in switch in the
   strip header. Draws thick top brackets over the L2 spans where
   `_bandTraceRegimeRuns` returned a run. Crosses into interpretation
   — only add if Quentin wants it after seeing real data.

2. **Hover tooltip on the strip** (~80 LOC) — per-L2 popup with
   entropy, dominant band/fraction, regime, n_valid. Mirrors the
   V-shape tooltip pattern (turn 157A).

3. **Combinatorial enumeration button** (~150 LOC) — loops through
   all subsets of bands at a focal L2 (capped at 2^K - 2), produces
   ranked list of (subset, run) pairs by signal strength. Output as
   a small dialog. This is the auto-finder.

4. **TSV export** (~60 LOC) — per-L2 row with l2_idx, n_valid,
   each band_fraction, entropy, regime, dominant_band. For
   manuscript supplementary data.

5. **Run the strip on real LG28** — calibrate the regime thresholds
   against actual data. The current `co_seg ≤ 0.40 / fanned ≥ 0.85`
   are guesses from synthetic tests; real founder-mixed cohort may
   need wider/narrower windows.

After 1+2+3+4 land, `SPEC_distant_band_concordance_fish_trajectory.md`
moves from `specs_todo/` to `specs_done/`.

---

## 7. Honest framing

**What's solid:**
- Compute + UI fully separated. The compute layer (turn 160) is
  testable in isolation; the UI layer (this turn) just plumbs.
- Every UI surface has a default that respects the manuscript's
  observational stance: toggle off, no interpretation overlays,
  observation-language labels.
- Cache invalidates correctly on (a) chrom switch, (b) fish-set
  change, (c) different K. Verified by test.
- The `🔍 trace` button is the one-click entry point Quentin asked
  for — focus a candidate, click, see the trace. No console
  required.

**What's risky:**
- Real-data thresholds are unknown. The synthetic test passes with
  `co_seg ≤ 0.40 / fanned ≥ 0.85` because the planted fixture is
  clean. Hatchery noise will probably need looser co_seg (~ 0.50)
  and looser fanned (~ 0.80). First real-LG28 run will tell.
- Strip + lineage strip together stack to 12 px above the plot.
  On small viewports this might crowd the panel header. No layout
  issue in current testing but worth eyeballing.
- The "trace from focal candidate" button picks the largest band by
  default. If the user wants to trace e.g. the HOM_INV band (which
  is often the smallest), they need the console call. A tiny
  band-picker dropdown next to the button would cover that — small
  follow-up.
- localStorage stores the fish-set as raw sample indices. If the
  user reloads a chrom that had a different cohort, the indices
  could refer to the wrong fish. Mitigated by the chrom-switch
  cache + fish-set clear, but a multi-page reload of the same chrom
  with different upstream data would still reuse stale indices.

**What's queued (in priority order):**
- Run on real LG28 → calibrate thresholds (Quentin)
- Run brackets + tooltip (small, mechanical)
- Combinatorial enumeration (medium, real "auto-find" surface)
- TSV export (small, manuscript)

End of handoff.
