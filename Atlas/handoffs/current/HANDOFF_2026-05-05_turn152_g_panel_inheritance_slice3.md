# HANDOFF — turn 152 — G-panel Slice 3 (inheritance tab) shipped

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (72,592 lines, +393 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.
**Supersedes**: turn 151 handoff. This turn closes the final
G-panel slice (`SPEC_g_panel_unified_groups.md` Slice 3) — the
inheritance tab now renders the real cross-candidate Jaccard cluster
view instead of the placeholder.

---

## 0. What Quentin asked for

> "I would like to fix the grouping tab because when we have candidates
> and inheritance like the slices haven't been added and the G UI not
> yet merged"

The G-panel had three tabs (karyotype / inheritance / manual). Slice 1
(manual, turn 135) and Slice 2 (karyotype, turn 136) shipped real
content. Slice 3 (inheritance) was still a placeholder text body
saying "Slice 3 pending". This turn ships Slice 3.

Result: all three group flavours now live in one unified popup.

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse, pragmatic.
Tarball is the standard handoff format.

---

## 2. What shipped this turn

### 2.1 Architectural reuse — no new compute

The cross-candidate Jaccard pipeline (`runInheritanceCompute`,
`inheritanceGroupClustering`, `_gatherActiveCandidatesForInheritance`,
`crossCandidateMatrix`, `openInheritanceMatrix`) all already existed
from turn 115 and turn 122. Slice 3 is **purely a wiring layer**: it
exposes the existing compute through the G-panel UI without changing
any compute logic.

This is intentional. The compute is already used by the I·g pills on
the lines panel and by the Cramér's V matrix popup. Slice 3 adds a
third consumer (the G-panel tab) that shares the same cache key
(`state.inheritanceCacheKey`) and same result object
(`state.inheritanceResult`), so all three surfaces stay in sync
automatically.

### 2.2 New state slot

```
state.gPanelInheritanceThreshold: 0.15
```

Cosine-distance threshold for the inheritance compute. Default 0.15
(= cosine similarity ≥ 0.85 to merge). Bounds clamped at
[0.05, 0.40]. Persisted to localStorage under
`inversion_atlas.gPanelInheritanceThreshold`.

### 2.3 Eight new helpers

- **`_gpInhEnsureThreshold()`** — load + clamp + write back. Reads
  localStorage; falls back to state default if missing/invalid.
- **`_gpInhSetThreshold(t)`** — clamp to [MIN, MAX], persist, invalidate
  inheritance cache. Cache invalidation is essential because the
  cache key in `_inheritanceCacheKey` does NOT include the threshold —
  changing the threshold without invalidation would return stale
  clusters.
- **`_gpInhSiToCga(si)`** — sample-index → CGA name. Self-contained
  for testability; mirrors `_mgSiToCga`'s contract.
- **`_gpInhMembersForGroup(result, groupId)`** — returns `Int32Array`
  of unique sorted sample indices that belong to one inheritance
  group. Walks `result.fish_masks[b]` for every band `b` whose
  `cut.group_id_per_band[b]` matches the requested groupId.
- **`_gpInhMakeManualGroup(result, groupId)`** — wraps
  `addToManualGroup` with default name `"inh_g<gid>"`. Returns the new
  group, or null if `addToManualGroup` isn't available.
- **`_gpInhItemLabel(result, itemIdx)`** — formats `"I1"` / `"I2"` from
  `items_meta[itemIdx].seq_num`; falls back to id, or `"?"` if missing.
- **`_gpInhRenderEmpty(reason, items)`** — empty-state body when
  n_candidates < 2 or no result yet.
- **`_gpInhRenderResult(result, items, threshold)`** — renders the
  per-candidate summary line + per-group cards.

All eight exposed on `window.*` for tests.

### 2.4 Inheritance tab body

`_gPanelRenderTabInheritance()` no longer returns the "Slice 3 pending"
placeholder. It now:

1. Calls `_gatherActiveCandidatesForInheritance()` (the same gatherer
   the I·g pills use). This filters to candidates with `locked_labels`
   that aren't auto-pending — exactly the gate the compute itself uses.
2. Renders a header with `n eligible candidates`.
3. If < 2: renders the empty-state explaining the gate rules.
4. Else: renders the threshold slider (0.05–0.40 cosine distance),
   the recompute button, the "matrix view" shortcut, and (if a result
   exists in `state.inheritanceResult`) the result body.

The tab does NOT auto-trigger a compute on first open — the user
clicks ↻ recompute to start. This matches Slice 2's pattern (the
karyotype tab also requires the user to do the lock-and-promote
upstream first; both tabs are render-only of pre-computed state).

### 2.5 Modal interaction wiring

Added an `if (activeTab === 'inheritance')` branch to
`_renderGPanelModal()`'s post-render section:

- **Threshold slider** — 220 ms debounce. Live-updates the value label
  during drag (no jitter), debounce fires
  `_gpInhSetThreshold(v) → runInheritanceCompute({ force: true,
  threshold: v }) → _renderGPanelModal()`.
- **Recompute button** (`#gpInhComputeBtn`) — same path, immediate.
  Useful when candidates have been edited mid-session and the user
  wants to refresh.
- **Matrix view button** (`#gpInhMatrixBtn`) — opens the existing
  `openInheritanceMatrix()` modal. The matrix overlays the G-panel;
  closing it returns to the G-panel.
- **Per-group "+ make manual" buttons** — each per-group card has
  `<button data-gpinh-mkmg-gid="N">`. Click handler reads gid, calls
  `_gpInhMakeManualGroup`, gives a brief "✓ added" visual, then
  re-renders. The next tab switch to "manual" shows the new group in
  the list.

### 2.6 Cleaned stale UI copy

- Header subtitle changed from `"Slice 1 ships manual tab; karyotype
  (Slice 2) and inheritance (Slice 3) are placeholders"` to
  `"all three flavours wired; manual = chrom/cohort splits, karyotype =
  per-candidate band membership, inheritance = cross-candidate Jaccard
  groups"`.
- Tab strip no longer renders `"(slice N)"` after each tab label. The
  underlying `_GPANEL_TABS` const still carries `slice` metadata for
  any downstream reader, but it's no longer surfaced in the UI.

---

## 3. What did NOT change

- **`runInheritanceCompute` / `inheritanceGroupClustering` /
  `_gatherActiveCandidatesForInheritance`** — untouched. Slice 3 is a
  consumer, not a re-implementation.
- **`openInheritanceMatrix` popup** — untouched. The new matrix-view
  button just calls into it.
- **Manual tab (Slice 1)** — untouched.
- **Karyotype tab (Slice 2)** — untouched.
- **`addToManualGroup`** — untouched. `_gpInhMakeManualGroup` wraps it.
- **`_inheritanceCacheKey`** — still doesn't include threshold. The
  workaround is in `_gpInhSetThreshold`, which calls
  `invalidateInheritanceCache()` on every threshold change.

---

## 4. Test status

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 151)   | 72,199  | 2557            | 55    |
| Post-turn-152 (current)  | 72,592  | 2630            | 56    |
| **Δ session**            | +393    | +73             | +1    |

Full sweep at turn 152: **2630 / 0** across `tests/test_turn*.js`. JS
syntax check: clean. HTML parser: 0 errors.

`tests/test_turn152_g_panel_inheritance_slice3.js` (70 tests) covers:

1. State slot + helper declarations land + are exposed on window
2. Tab body — placeholder copy gone; new wiring (gatherer call,
   threshold helper, slider/buttons) in place
3. Result renderer — per-candidate summary, per-group cards,
   make-manual buttons, n_fish/n_bands surfaced, empty-group message
4. Modal wiring — slider debounce, recompute, matrix shortcut,
   per-group mkmg
5. Header / tab-strip cleaned of stale Slice-N copy
6. Sandboxed unit tests on all 8 helpers (defaults, clamping,
   localStorage roundtrip, cache invalidation, empty/missing inputs,
   union-and-sort behaviour)
7. Window exports
8. Existing flow preserved (compute, matrix, manual, karyotype tabs)

`tests/test_turn135_g_panel_slice1.js` patched (+8 assertions, two
existing assertions inverted per the "codified placeholder" rule from
turn 150's handoff): now verifies that the placeholder is gone and the
real Slice 3 wiring is rendered.

---

## 5. What Quentin can do next session

The G-panel popup is now feature-complete for all three group
flavours. To exercise:

1. Press `g` (or click `G ▾` in the L3 toolbar) → popup opens.
2. Switch to the **inheritance** tab.
3. If you have ≥ 2 candidates with locked labels, you'll see:
   - Header: `Inheritance groups · N eligible candidates`
   - Slider + threshold value + ↻ recompute + matrix view buttons
   - Empty-state until you click ↻ recompute
4. Click ↻ recompute → cards appear, one per inheritance group, with:
   - `g0` / `g1` / ... id chip
   - `n_fish=X · n_bands=Y` summary
   - `+ make manual` button per card
   - Member band list `I1·b0, I3·b2, ...`
5. Drag the threshold slider — clusters re-compute with debounce,
   tab re-renders.
6. Click any `+ make manual` → that inheritance group becomes a manual
   group named `inh_g<N>` in the manual tab and on the sidebar list.
7. Click `matrix view` → the existing Cramér's V matrix popup opens
   over the G-panel.

Quentin will see whether:
- the inheritance compute produces sensible group counts at default
  threshold,
- per-candidate group counts match the I·g pill labels on the lines
  panel,
- some inheritance groups make sense as manual groups (likely useful
  for "fish that share H1/H2 across all candidates").

---

## 6. Things I almost broke and fixed

- **Turn 135 test broke.** Two assertions codified the Slice 1
  placeholder behaviour (`"Slice 3 pending"` and the "Slice 2 / Slice 3
  placeholder framing" header copy). Per the turn 150 convention,
  pre-existing tests that codified flagged-but-not-fixed behaviour
  must be inverted when the underlying behaviour is corrected. Test
  6c rewritten with 8 new assertions exercising the real Slice 3
  body; the header-copy assertion inverted to check for the new
  "all three flavours wired" subtitle.
- **`_IGC_DEFAULT_COSINE_DIST_THRESHOLD` not in sandbox.** The new
  `_gpInhEnsureThreshold` references this constant as a fallback.
  The turn 135 sandbox didn't load it. Fixed by injecting the const
  directly into the sandbox before extracting the helpers.
- **Cache key doesn't include threshold.** `_inheritanceCacheKey`
  hashes locked-labels + bp ranges + mode but not threshold. Without
  invalidation, a slider change would return stale clusters. Fixed
  by calling `invalidateInheritanceCache()` from inside
  `_gpInhSetThreshold` so any slider movement guarantees a fresh
  compute.

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 152 patched, 72,592 lines.
- `tests/test_turn152_g_panel_inheritance_slice3.js` — NEW, 70 tests.
- `tests/test_turn135_g_panel_slice1.js` — patched: 8 new assertions,
  2 existing assertions inverted (codified placeholder behaviour
  removed per turn 150 convention).
- `HANDOFF_2026-05-05_turn152_g_panel_inheritance_slice3.md` — this file.

---

## 8. Honest framing

**What's solid:**
- All three G-panel slices ship real content. UI feature-complete for
  the unified groups workflow.
- The inheritance tab is a thin wiring layer over compute that already
  worked. Risk surface is small.
- 2630 / 0 tests passing. Test discipline maintained: pre-existing
  codified placeholder behaviours inverted, not just deleted.

**What's NOT done (and why that's right):**
- **Threshold inclusion in cache key.** Workaround via
  `invalidateInheritanceCache()` in the setter is correct but feels
  fragile. A proper fix would extend `_inheritanceCacheKey` to include
  the threshold; that's a separate ~10-LOC change touching turn-115
  code. Not in scope this turn.
- **Per-group fish-list expansion.** The card shows `n_fish=N` but
  not the actual CGA list (it'd be too much for 100+ fish per group).
  A "show fish" expand-toggle could come if Quentin wants it.
- **Per-group color preview.** Cards currently use a single accent
  color (`var(--good)`); coloring per-group from a palette would be
  visually nicer but isn't necessary.
- **Inheritance compute progress indicator.** Compute is synchronous;
  for 226 fish × ~10 candidates × 3-6 bands per candidate it's
  sub-100 ms in practice. If the user has many candidates it may feel
  laggy on slider drag. The 220-ms debounce mitigates but a real
  spinner would be nicer if it ever matters.

**What's queued:**
- The four directions from end of turn 151 are all still doable:
  - A — UV-rotation refactor on locked-labels
  - B — FIG_C30 V-shape (`(u, agreement_fraction)`) plot
  - C — per-sample-lines het coloring
  - D — Quentin pivots based on what het reveals

Slice 3 doesn't conflict with any of those.

End of handoff.
