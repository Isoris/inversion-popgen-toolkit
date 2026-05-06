# HANDOFF — turn 155 — threshold folded into inheritance cache key

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (72,917 lines, +56 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.

**Important context**: between turn 153 (mine) and this turn, a parallel
session shipped **turn 154** (G-panel inheritance compute UX hardening
— `inheritanceLastStatus`, auto-compute on mount, status pill, recompute
button feedback, distinct messages for 0-groups vs null-result). I picked
up from the post-154 working tree (2746 / 0 tests). Turn 155 is the
honesty fix flagged in turn 152 §6 and turn 154 §6.

---

## 0. What this turn fixes

The cache-key honesty hole that's been outstanding since turn 152:

- Pre-turn-155 `_inheritanceCacheKey(items, mode)` did NOT fold in the
  cosine-distance threshold.
- Two computes at different thresholds wrote the SAME cache key.
- The cache-hit branch in `runInheritanceCompute` would return stale
  clusters whenever the threshold changed without an intervening
  `force: true` or `invalidateInheritanceCache()` call.
- The workaround in `_gpInhSetThreshold` (calling
  `invalidateInheritanceCache` from the slider) was paper-thin: any
  consumer that bypassed the slider (e.g. a `runInheritanceCompute({
  threshold: x })` call at a non-default threshold) would write a
  wrongly-keyed entry, and any caller that built an "expected key"
  without reading `state.gPanelInheritanceThreshold` would mis-detect
  staleness.

The turn 152 handoff §6 and turn 154 §6 both flagged this as a "real fix
is ~10 LOC, separate turn". This is that turn.

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse, pragmatic.
Tarball is the standard handoff format.

---

## 2. What shipped this turn

### 2.1 `_inheritanceCacheKey` signature extended

```
function _inheritanceCacheKey(items, mode, threshold)
```

Threshold is OPTIONAL — backwards-compatible with old `(items, mode)`
callers. When omitted, falls back through:

1. `state.gPanelInheritanceThreshold` (the slider-driven canonical
   threshold)
2. `_IGC_DEFAULT_COSINE_DIST_THRESHOLD` (0.15) if state isn't set up
3. literal 0.15 if the constant somehow isn't loaded

This means **all six existing callers** (in `runInheritanceCompute`,
`_autoRegisterInheritanceOnCandidateChange`, `_drawInheritanceLabelsStrip`,
`runL2SweepInheritance`, the G-panel mount auto-compute, and the
cockpit/page-2 stale guard) keep working without code change. They each
now resolve the threshold from state by default, so they all agree by
construction. **One global source of truth** = `state.gPanelInheritanceThreshold`.

The threshold is rounded to 4 decimals (matching the slider's 0.01 step)
before going into the key, so floating-point jitter (e.g.
`0.150000000001`) doesn't fragment the cache.

Key shape:

```
{mode}@t{thresholdString}::{id1}@K{K1}@{bp1}@{fp1}|{id2}@K{K2}@{bp2}@{fp2}|...
```

### 2.2 `runInheritanceCompute` resolves threshold BEFORE cache key

Pre-155 the order was:
```
const cacheKey = _inheritanceCacheKey(items, mode);   // no threshold
const threshold = (opts && opts.threshold != null) ? +opts.threshold : _IGC_DEFAULT;
```

Post-155:
```
let threshold;
if (opts && opts.threshold != null) {
  threshold = +opts.threshold;
} else if (Number.isFinite(_state.gPanelInheritanceThreshold)) {
  threshold = _state.gPanelInheritanceThreshold;
} else {
  threshold = _IGC_DEFAULT_COSINE_DIST_THRESHOLD;
}
const cacheKey = _inheritanceCacheKey(items, mode, threshold);
```

Two consequences:
- **Cache key is honest** — folds in the threshold that compute will
  actually use.
- **Default threshold reads from state** — so `runInheritanceCompute()`
  with no opts now uses the user's slider value, not always 0.15. This
  closes the silent divergence where the lines-strip auto-fired
  compute always used 0.15 regardless of what the slider showed.

### 2.3 `inheritanceGroupClustering` receives resolved threshold

Pre-155 the call was `inheritanceGroupClustering(items, opts)`,
forwarding the original opts object. If `opts.threshold` was null/missing,
`inheritanceGroupClustering` would internally use `_IGC_DEFAULT` (0.15)
even if the cache key had been built at a different threshold. Now:

```
inheritanceGroupClustering(items, Object.assign({}, opts || {}, { threshold: threshold }))
```

The cluster compute and the cache key are guaranteed to agree.

### 2.4 `_gpInhSetThreshold`'s `invalidateInheritanceCache` retained

The slider's `invalidateInheritanceCache` call still runs, but it's now
**defense-in-depth** rather than the primary mechanism:

- After turn 155, a subsequent `runInheritanceCompute` would already
  cache-miss on threshold change (the cache key now differs).
- The invalidate still nulls `state.inheritanceResult` so any consumer
  that reads the result *directly* (without going through
  `runInheritanceCompute`) sees null instead of a stale result for the
  brief window until the next compute fires.

Comment updated to reflect the new role.

---

## 3. What did NOT change

- **`_hashLockedLabels`** — untouched.
- **`inheritanceGroupClustering`** — untouched. Just gets a more
  reliable opts object now.
- **`invalidateInheritanceCache`** — untouched.
- **All 6 callers of `_inheritanceCacheKey`** — untouched. They use the
  optional-arg default lookup, which resolves to
  `state.gPanelInheritanceThreshold`.
- **Turn 154 status tracking (`inheritanceLastStatus`)** — untouched.
- **Turn 153 auto-register hook** — untouched.
- **Turn 152 G-panel inheritance tab** — untouched.

---

## 4. Test status

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 154)   | 72,861  | 2746            | 58    |
| Post-turn-155 (current)  | 72,917  | 2776            | 59    |
| **Δ session**            | +56     | +30             | +1    |

Full sweep: **2776 / 0**. JS syntax: clean. HTML parser: 0 errors.

`tests/test_turn155_threshold_in_cache_key.js` (30 tests) covers:

1. Signature extended; defaults documented; rounding to 4 decimals;
   key string format
2. `runInheritanceCompute` resolves threshold first, builds key, passes
   resolved threshold to `inheritanceGroupClustering`
3. `_gpInhSetThreshold` retains `invalidateInheritanceCache` as
   defense-in-depth
4. Keys differ when threshold differs (sandboxed); same inputs yield
   identical keys; FP jitter coalesces; existing item / mode
   distinctions still hold
5. Optional threshold defaults to state then to `_IGC_DEFAULT`;
   2-arg backwards-compat works
6. Existing contracts (turn 129 / 152 / 153 / 154) all preserved

### Tests that needed inverting

- **`test_turn153_inheritance_auto_register.js`** — the assertion
  `_inheritanceCacheKey still exists` codified the 2-arg signature.
  Inverted to require the 3-arg signature with a comment explaining
  the turn 155 extension. Backwards-compat itself is exercised
  separately in turn 155's §5c.

### Tests that needed regex tightening

- **`test_turn154_compute_ux_hardening.js`** — the `wraps
  inheritanceGroupClustering in try/catch` assertion used a tight
  regex. Turn 155 inserted a comment between the `try {` and the call,
  so the regex now allows whitespace+content between them.

---

## 5. What Quentin can do next session

The cache contract is now honest. Practical implications:

1. **Slider drag at sub-100ms responsiveness.** Pre-155 the slider's
   debounced handler called `_gpInhSetThreshold(v)` →
   `invalidateInheritanceCache()` → `runInheritanceCompute({ force:
   true, threshold: v })`. Now the `force: true` is redundant — the
   key alone would mismatch. Behaviour is identical, just architecturally
   cleaner.

2. **Lines-strip auto-fire respects the slider.** Pre-155, opening the
   G-panel, dragging the slider to 0.20, then navigating to a page
   that draws the lines strip would auto-fire `runInheritanceCompute()`
   at threshold 0.15 (the function default), not 0.20. The pills would
   show clusters built at the wrong threshold. Post-155, the auto-fire
   uses `state.gPanelInheritanceThreshold = 0.20`. **The slider is now
   genuinely the single source of truth.**

3. **Cache hits where they should be.** A user who clicks ↻ recompute
   twice in a row at the same threshold gets a cache hit on the second
   click (turn 154's `cached` status), because the keys match. A user
   who clicks ↻ recompute, drags to a new threshold, and clicks ↻
   again, gets a fresh compute (turn 154's `computed` status), because
   the keys differ. Pre-155 these distinctions were unreliable.

To exercise:
- Open G-panel, click ↻ recompute → status pill shows `computed: N
  groups`
- Click ↻ again immediately → status pill shows `cache hit: N groups`
- Drag slider to new value → live update; auto-recompute → status pill
  shows `computed: N groups @ threshold 0.20`
- Drag slider back to 0.15 → cache hit on the original computation
  (status pill: `cache hit: N groups`). Pre-155 this would have been a
  silent stale-result scenario.

---

## 6. Things I almost broke and fixed

- **`runInheritanceCompute(opts)` forwarding stale opts to
  `inheritanceGroupClustering`.** Caught during the design phase —
  cache-key threshold and cluster-compute threshold would have
  silently disagreed when `opts.threshold` was null but
  `state.gPanelInheritanceThreshold` was non-default. Fixed by
  injecting the resolved threshold into a fresh opts object.
- **Lines-strip cache-key compare against state slider value.** Caught
  during call-site audit — without this, every redraw with a non-default
  slider value would have triggered an "infinite mismatch" loop where
  the auto-fire wrote keys at threshold=0.15 and the next draw expected
  threshold=0.20. Fixed by making `runInheritanceCompute()`'s default
  threshold ALSO read from `state.gPanelInheritanceThreshold`, which
  the optional-arg default in `_inheritanceCacheKey` already does.
- **Floating-point key fragmentation.** Slider `value` is a string that
  parses to (e.g.) `0.150000000001`. Without rounding, two consecutive
  identical drags could produce different cache keys. Fixed with
  `Math.round(t * 1e4) / 1e4` then `.toFixed(4)`.
- **Turn 153 test codified the 2-arg signature.** Per the convention,
  inverted in place rather than deleted. Backwards-compat is preserved
  *runtime* — only the signature documentation is updated.

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 155 patched, 72,917 lines.
- `tests/test_turn155_threshold_in_cache_key.js` — NEW, 30 tests.
- `tests/test_turn154_compute_ux_hardening.js` — patched (1 regex).
- `tests/test_turn153_inheritance_auto_register.js` — patched
  (1 inverted assertion per signature change).
- `HANDOFF_2026-05-05_turn155_threshold_in_cache_key.md` — this file.

Plus prior handoffs (carried for reference):
- `HANDOFF_2026-05-05_turn154_compute_ux_hardening.md`
- `HANDOFF_2026-05-05_turn153_inheritance_auto_register.md`
- `HANDOFF_2026-05-05_turn152_g_panel_inheritance_slice3.md`

---

## 8. Honest framing

**What's solid:**
- The honesty hole flagged across three handoffs is now closed. Single
  source of truth for the threshold; cache key includes it;
  `runInheritanceCompute` and `inheritanceGroupClustering` use the same
  resolved value.
- Backwards-compat preserved — old `(items, mode)` callers keep working
  via the optional-arg default that resolves from state.
- Float-jitter coalesced with 4-decimal rounding (matches slider
  step).
- Defense-in-depth retained: `_gpInhSetThreshold` still calls
  `invalidateInheritanceCache` so direct readers of
  `state.inheritanceResult` see null instead of stale.

**What's NOT done (and why that's right):**
- **B (FIG_C30 V-shape)** — was queued alongside E in this session's
  plan. Skipped this turn because turn 155 expanded slightly during the
  call-site audit (the lines-strip / cockpit threshold consistency
  issue I caught was real and needed fixing here, not deferred). B is
  still bounded ~150 LOC and the next obvious scientific deliverable.
- **No instrumentation for cache hit vs miss rate.** Turn 154 added
  the status pill which surfaces individual outcomes. A session-level
  aggregate ("you triggered N computes this session, M cache hits")
  isn't necessary unless Quentin reports the slider feels laggy.

**What's queued:**
- **B** — FIG_C30 V-shape `(u, agreement_fraction)` plot. ~150 LOC.
  Bounded scientific deliverable; no dependency on the inheritance
  stack.
- **A** — UV refactor on locked groups. Locked-labels caveat still
  applies.
- **C** — Per-sample-lines het coloring. Perf concerns at 226 × 30k.
- **D** — Pivot based on what het / inheritance reveals.
- **F** — Per-group "show fish" expand toggle in the G-panel
  inheritance card. ~50 LOC. Useful if Quentin wants to spot-check
  group membership before promoting to manual.

End of handoff.
