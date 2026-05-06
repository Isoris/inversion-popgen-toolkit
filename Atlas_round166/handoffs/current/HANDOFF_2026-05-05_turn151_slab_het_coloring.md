# HANDOFF — turn 151 — Slab het-coloring parity shipped

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (72,199 lines, +109 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.
**Supersedes**: turn 150 handoff. This turn closes the het-coloring gap
that turn 129 (Slice 1 of `SPEC_l3_het_dosage_coloring`) left for slab
mode after turn 148–150 brought slab-mode L3 to structural parity with
L2 mode.

---

## 0. What Quentin asked for

> "We must first understand het dosage so try to connect it and we will
> see later. Last time the locked didn't work so we fixed it."

The earlier survey doc proposed a UV-rotation refactor using
locked-labels. Quentin redirected: **don't go there yet**, locked-labels
has a history of burning, **wire het dosage first**, look at it, then
decide.

Inspection found that the het-dosage path is ~80% scaffolded:
`state.l3HetColoring`, the `#l3HetToggle` checkbox in the L3 toolbar,
`_hetRateColor`, `_computeHetRateForL2`, the per-L2 cache, the
focal-pane gradient swatch, localStorage persistence, the
"disabled-until-dosage_chunks-loaded" gating — **all already exist
from turn 129**. The only missing piece for end-to-end was the slab
variant of the same wiring. This turn closes that.

After turn 151, when Quentin loads dosage_chunks and toggles
`#l3HetToggle`:
- L2-mode mini-PCAs paint by per-sample het rate (existing, turn 129)
- **Slab-mode mini-PCAs also paint by per-sample het rate (NEW)**
- Both paths render the same gradient swatch + 0/0.5/1 labels
- Both paths share the same `state.__hetRateCache` (different keys)

He can then look at what het actually shows on the contingency, and
either greenlight the UV-rotation work, request a different direction,
or ask for the FIG_C30 V-shape coherence plot next.

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K=8 clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse, pragmatic.
Multiple parallel chat sessions; tarball is the standard handoff format.

---

## 2. What shipped this turn

### 2.1 Helper refactor: range-keyed core + slab wrapper

`_computeHetRateForL2(l2idx)` was the only het-rate compute function
before. It inlined a chunk fetch + bp-range marker filter + per-sample
het tally. Refactored:

- **`_computeHetRateForRange(start_bp, end_bp, cacheKey)`** — extracted
  core. Takes a bp range and an arbitrary cache key. Cache key is
  `null`-safe (passing `null` skips caching, useful for tests).
- **`_computeHetRateForL2(l2idx)`** — now a thin wrapper. Pulls
  `env.start_bp` / `env.end_bp` from `state.data.l2_envelopes[l2idx]`
  and delegates with `cacheKey = l2idx`. Behaviour preserved
  byte-for-byte for callers (tests confirm).
- **`_computeHetRateForSlab(start_w, end_w)`** — NEW. Translates
  window indices to bp via `state.data.windows.start_bp[start_w]` and
  `state.data.windows.end_bp[end_w]`, delegates with
  `cacheKey = 'slab:<s>:<e>'`.

Validates window-index bounds (`< 0`, `>= n_windows`, `end < start`)
and falls through to all-NaN when the windows array isn't bp-annotated.

### 2.2 Cache contract

`state.__hetRateCache` is a single `Map`, shared between L2 and slab
paths. Keys are heterogeneous:
- **L2 path keys** are integers (`l2idx`)
- **Slab path keys** are strings (`'slab:<s>:<e>'`)

`Map` doesn't coerce keys, so the two namespaces don't collide. A test
explicitly verifies that L2(idx=0) and Slab(0..2) covering the same
bp range produce identical numbers but cache under different keys.

`_invalidateHetRateCache()` clears both kinds in one drop. No callsite
changes needed.

### 2.3 `drawSlabMiniPCA` honors the het toggle

The slab dot loop (line 45537–45544 pre-turn, now ~45550–45580) reads
`state.l3HetColoring`. When on, it calls `_computeHetRateForSlab` for
the slab range and overrides each dot's fill with `_hetRateColor(rate)`.
When off, behaviour unchanged: `groupColor(label)` per the K-means label.

Match parity with `drawMiniPCA`'s lines 47120–47144:
- Same `_hetMode` boolean
- Same `_hetRates` Float32Array
- Same fallback to label color when rate is NaN (via `_hetRateColor`)

### 2.4 Slab gradient swatch legend

Same 36×6 gradient with `_hetRateColor(0.00 / 0.25 / 0.50 / 0.75 / 1.00)`
color stops, same backdrop (`themeColor('bg')` at 0.92 alpha), same
"het" caption, same 0 / 0.5 / 1 tick labels under the swatch. Drawn
only when `state.l3HetColoring` is on.

The slab swatch uses the same anchor as the L2 swatch (bottom-left of
the plot area) so the two modes look identical at a glance. There is
no slab-specific paneOffset gate — slab mini-PCAs are typically a single
focal-style pane in the contingency, not a 5-pane carousel.

### 2.5 `state.__hetRateCache` reuse

Whatever invalidation paths existed for the L2 path
(`dosage_chunks` reload, manual reset) automatically apply to slab
keys too — the cache itself is shared. Verified in the test: after
populating both `0` (L2) and `'slab:0:2'`, calling
`_invalidateHetRateCache()` clears both.

---

## 3. What did NOT change

- **`drawMiniPCA`** — untouched. L2-mode het coloring already worked
  end-to-end from turn 129; nothing to fix.
- **`#l3HetToggle` checkbox + handler** — untouched. Already wires
  `state.l3HetColoring`, persists to localStorage, gates on
  `dosage_chunks` layer presence, force-off on layer disappearance.
- **Per-sample-lines color picker (`_resolveSampleColorByMode`)** — still
  a stub for `het` mode at line 32761. The picker UI is live but the
  per-window color override on the lines panel would require segment-
  per-window strokes which is a separate scope. NOT touched.
- **UV rotation** (`_getOrComputeUVRotation`) — untouched per Quentin's
  redirect.
- **Locked-label paths** — untouched per Quentin's redirect.

---

## 4. Test status

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 150)   | 72,090  | 2504            | 54    |
| Post-turn-151 (current)  | 72,199  | 2557            | 55    |
| **Δ session**            | +109    | +53             | +1    |

Full sweep at turn 151: **2557 / 0** across `tests/test_turn*.js`. JS
syntax check: clean. HTML parser: 0 errors.

`tests/test_turn151_slab_het_coloring.js` (52 tests) covers:

1. Source-pattern: helper declarations land + are exposed on window
2. L2 wrapper delegates to range core (no inlined chunk fetch)
3. Slab helper translates windows → bp, bounds-checked, NaN-safe
4. `drawSlabMiniPCA` calls `state.l3HetColoring`, `_computeHetRateForSlab`,
   `_hetRateColor`, gradient + ticks
5. Sandboxed unit tests on `_computeHetRateForRange` (no chunk, cache
   hit, invalid range, real synthetic chunk, marker filter)
6. Sandboxed unit tests on `_computeHetRateForSlab` (no chunk, oob,
   marker filter via window→bp mapping, cache key shape)
7. L2/slab parity: same numbers, different cache keys, single
   invalidation
8. L2 path regression preserved: missing env, no chunk, cache key shape

`tests/test_turn129_l3_het_coloring.js` patched (+1 test, +6 sandbox
load lines) so its in-sandbox runs of `_computeHetRateForL2` also load
the new range core. The patch was minimal: every
`vm.runInContext(fnComputeHet, sandbox)` is now preceded by
`vm.runInContext(fnComputeHetRange, sandbox)`. The L2 helper test
suite passes 59/59 with the patch.

---

## 5. What Quentin can do next session

The het toggle is now functional in both L2 and slab views. To exercise
it:

1. Load a chromosome (e.g. LG28).
2. Load dosage data via the existing heatmap toolbar / stripe quality
   button (these populate `popgenDosage.getCachedChunk`).
3. Click the L3 het toggle (`#l3HetToggle`).
4. Observe the contingency mini-PCAs paint by per-sample het rate:
   cold blue = ~0 het (looks homozygous for one or other arrangement),
   neutral = ~0.5 het (looks like a HET sample), warm red = >0.5 het
   (excess heterozygosity).
5. Switch to a slab view (slab cycler in the L2 panel). Toggle stays on.
   **NEW**: the slab mini-PCA also paints by het. The gradient swatch
   appears bottom-left.
6. Inspect: does the het coloring align with the K-cluster halos? Are
   there outlier samples (warm-red dots in a "homo" cluster, or
   cold-blue dots in a "HET" cluster)? These are the signal Quentin
   wanted to see.

Once Quentin has a sense of what the het pattern shows, he can pick a
direction:

- **Direction A (UV refactor)** — the original survey doc's plan.
  Use locked-groups for the rotation, DBSCAN per stripe. Now better
  informed by what het actually shows.
- **Direction B (FIG_C30 V-shape)** — port `(u, agreement_fraction)` as a
  diagnostic plot. Independent of UV refactor; uses the existing
  `computeStripeQuality` v3.94 port.
- **Direction C (per-sample-lines het coloring)** — fill in `_hetColor`
  in `_resolveSampleColorByMode`. Larger refactor: per-window colored
  segments on the line traces.
- **Direction D (something else entirely)** — Quentin pivots based on
  what het reveals.

---

## 6. Things I almost broke and fixed

- **Turn 129 test regression.** Refactoring `_computeHetRateForL2` to
  delegate broke turn 129's sandboxed runs that extracted that single
  function via `pullFunction()`. Fix: turn 129's sandbox setup now also
  loads `_computeHetRateForRange`. Six call sites patched; one new
  source-extraction assertion added.
- **Float32 precision in test assertions.** Initial tests used `===`
  comparisons against `0.42` and `0.6`. Float32Array stores rounded
  values (`0.41999998...`, `0.6000000238...`). Switched to
  `Math.abs(... - target) < 1e-6`.
- **Regex for combined `const swW = 36, swH = 6;` declarations.**
  Initial regex `const\s+swH\s*=\s*6` failed because `swH` is the
  second part of a comma-separated declaration with no `const`
  immediately preceding it. Switched to `swH\s*=\s*6`.

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 151 patched, 72,199 lines.
- `Inversion_atlas.html.original` — turn 150 backup, 72,090 lines.
- `tests/test_turn151_slab_het_coloring.js` — NEW, 52 tests.
- `tests/test_turn129_l3_het_coloring.js` — patched: +1 source assertion,
  +6 sandbox load lines.
- `HANDOFF_2026-05-05_turn151_slab_het_coloring.md` — this file.

---

## 8. Honest framing

**What's solid:**
- Slab mode now has full het-coloring parity with L2 mode.
- Test discipline maintained: 2557 / 0, every turn's loose ends closed.
- The refactor is a clean factoring (range-keyed core + two thin
  wrappers); future range-based callers can use the core directly.

**What's NOT done (and why that's right):**
- UV rotation refactor — Quentin redirected.
- Locked-labels integration — Quentin redirected.
- Per-sample-lines het coloring — separate scope (segment rendering
  is a larger commitment).
- FIG_C30 V-shape coherence plot — separate scope, awaits Quentin's
  direction after looking at het.
- An end-to-end DOM test against a real cohort fixture — test harness
  is structurally source-pattern + sandboxed-unit, no DOM emulation.
- Direct loading of dosage_chunks JSON — still requires the user to
  invoke the heatmap toolbar to populate `popgenDosage.getCachedChunk`.
  That's the existing gating; not changed this turn.

**What's queued:**
- Quentin reviews het mode visually, picks Direction A/B/C/D.
- All four Directions remain doable without conflicting with this
  turn's work.

End of handoff.
