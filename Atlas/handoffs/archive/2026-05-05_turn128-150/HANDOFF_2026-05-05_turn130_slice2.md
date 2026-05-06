# HANDOFF — turn 130 (Slice 2 follow-up)

**Date**: 2026-05-05
**Atlas line count**: 64,160 (was 63,901 → +259)
**Test status**: **1868 PASS / 0 FAIL** across all suites (was 1824 → +44)

---

## What shipped (Slice 2 of fish-trajectory)

The lineage compute from Slice 1 now surfaces in the UI. Two visible
features:

### 1. `lineage` colormode in the per-sample-lines picker

New `<option value="lineage">` in the lines header dropdown. When
selected, every fish line paints in its lineage's HSL color (golden-angle
rotation: lineage 0 = blue 210°, each next +137.508°). Always available
— no JSON layer prerequisite. Reads `state.lineageResult.lineage_id_per_sample`.

Resolver: `_lineageColor(si)` — handles missing result (auto-triggers
compute via `requestIdleCallback`, returns null until next paint),
out-of-range si, and lineage_id == -1 (no-call).

Wired into both `_resolveSampleColorByMode(si, gi, 'lineage')` and
`_resolveSampleScopeColor(si, 'lineage')` so lineage coloring works on
the per-sample-lines panel, the L3 mini-PCAs, and the tracked-samples
scatter.

### 2. Lineage strip above the per-sample-lines panel

New `_drawLineageStrip(ctx, pad, plotW, plotH, mbMin, mbMax)` — sibling
to the existing regime-breadth strip. For each L2 envelope:

- **In-chain L2** → colored bar with the dominant lineage of the
  largest band (same golden-angle palette as the line colors)
- **Chain-break L2** → faint diagonal-hatch grey ("Hungarian couldn't
  bridge here")

Toggle: `[x] lineage` in lines header (default ON), state slot
`state.linesLineageStripOn`, localStorage key
`inversion_atlas.linesLineageStripOn`. Setter `setLinesLineageStripOn(b)`.

### 3. Cache invalidation on chrom switch

`applyData` now calls `invalidateLineageCache()` and resets
`state._lineageComputeScheduled = false` when `state.data` swaps. The
next paint re-triggers compute on the new chromosome's L2 inventory.

### Tests

`tests/test_turn130_lineage_ui.js` — 44 tests:
- 22 source-level (mode registration, dropdown, resolver dispatch,
  strip drawer, state slot, setter, localStorage, toggle markup,
  applyData hook, fallback paths)
- 22 behavioural sandboxed (golden-angle hue math at lineage 0/1/2/3,
  edge cases — missing result auto-triggers, invalid si returns null,
  lineage_id -1 returns null, mode dispatch through both resolvers,
  fewer-than-3-L2s gates the auto-trigger off)

---

## File deltas

```
Inversion_atlas.html    63901 → 64160  (+259 lines)

NEW tests:
  tests/test_turn130_lineage_ui.js                   44

UPDATED specs:
  specs_todo/SPEC_distant_band_concordance_fish_trajectory.md
    Slice 1 + 2 marked shipped at top + checklist

Backups (drop before bundling):
  Inversion_atlas.html.bak_pre_lineage_ui
  Inversion_atlas.html.bak_post_lineage_ui
```

## What Quentin sees on next atlas reload

Open a chromosome with ≥3 L2 envelopes. The lineage compute fires
automatically on first paint via `requestIdleCallback`. Then:

1. **Top of per-sample-lines panel** → colored strip showing per-L2
   dominant lineage. Color changes mark where the underlying inheritance
   re-organizes along the chromosome. Grey-hatched gaps = chain breaks
   (Hungarian agreement < 0.5 between adjacent L2s).
2. **Lines header `color:` dropdown** → switch to `lineage` to recolor
   every fish line by its inheritance lineage.
3. **L3 mini-PCAs and tracked-samples scatter** → also honor lineage
   mode (the per-sample resolver dispatches there too).

The toggle `[x] lineage` in the header lets you turn the strip off if
it's noisy on a chromosome with chaotic K-means.

## Next-turn priorities

In rough order of value-per-turn:

1. **Slice 3** — co-membership matrix viewer modal. Click in the lines
   panel → see the 226×226 concordance heatmap reordered by hierarchical
   clustering. ~0.5 turn.

2. **L2-sweep Slice 1** — auto-promote candidates from sweep results.
   Removes the "user must promote candidates first" gate. ~1 turn.
   Spec: `SPEC_l2_sweep_inheritance.md`.

3. **Slice 4** — combinatorial multi-L2 lasso. Quentin's "select 2-3
   bands and see they visit different bands further along" UI.
   ~1 turn.

4. **Diagnostic protocol on real data** — Quentin opens the atlas with
   the new lineage strip + colormode and reports what the lineages look
   like on LG28. Decides whether thresholds need calibration.
