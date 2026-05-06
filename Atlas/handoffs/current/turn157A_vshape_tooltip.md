# HANDOFF — turn 157 — V-shape tooltip on hover

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (73,667 lines, +232 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.

**Picked up from**: post-turn-156 (2840 / 0). Small incremental turn —
wires the `hits[]` array that turn 156 already returned into a hover
tooltip.

---

## 0. What this turn ships

The V-shape diagnostic plot from turn 156 returned a `hits[]` array
specifically so a follow-up turn could wire mousemove → "this sample
is CGA xxx, group HET, agreement 0.42, stripe_quality peripheral".
This is that follow-up.

---

## 1. Cohort discipline

Same as turn 156. Quentin Andres, Kasetsart University Bangkok.
Direct, terse, pragmatic.

---

## 2. What shipped this turn

### 2.1 Six new helpers

```
_vsTooltipEnsureEl()                — lazy-creates #vsPointTooltip on body
_vsTooltipBuildHtml(hit)            — renders sample-detail HTML
_vsTooltipShow(hit, clientX, clientY) — positions + clamps to viewport
_vsTooltipHide()                    — display:none
_vsHitTest(hits, mouseX, mouseY)    — closest-point hit test
_wireVShapeTooltip(canvas)          — idempotent mousemove + mouseleave
```

All six exposed on `window.*` for tests + downstream.

### 2.2 New module-level state

```
let _vsLastRender = null;
```

Snapshot of the most-recent V-shape render — `{ hits, data,
candidate_id }`. Set inside `openVShapePlot` immediately after
`_drawVShapePlot` returns. Cleared (implicitly) when the modal is
re-opened against a different candidate.

### 2.3 `openVShapePlot` integration

```
const result = _drawVShapePlot(canvas, data, { ... });
_vsLastRender = {
  hits: (result && Array.isArray(result.hits)) ? result.hits : [],
  data: data,
  candidate_id: cand.id,
};
if (typeof _wireVShapeTooltip === 'function') {
  try { _wireVShapeTooltip(canvas); } catch (_) {}
}
```

Wiring is idempotent (the canvas's `dataset.vsTooltipWired` flag
prevents double-binding on re-open). Try/catch is defense-in-depth.

### 2.4 `closeVShapePlot` also hides the tooltip

Prevents an orphan tooltip from persisting on the page after the
modal is dismissed.

### 2.5 Tooltip content

Per-sample card. Header: colored dot (group color) + CGA + group
text. Body: `u` (4 decimals), `agreement` (4 decimals),
`coherence_class` (with semantic color: green for coherent, amber for
intermediate, red for discordant), `stripe_quality` tier (green for
core, amber for peripheral, red for junk).

`coherence_class` and `stripe_quality` are optional — if missing from
the hit's point, that row is just skipped. CGA falls back to "?",
non-finite numeric values render as "NA".

### 2.6 CSS-pixel coordinate space (intentional!)

**KEY DIFFERENCE from the inheritance pill tooltip pattern (turn 122):**
the inheritance code stores hit regions in raw canvas-internal
coordinates and re-applies `canvas.width / rect.width` DPR scaling at
hover time. The V-shape canvas uses `setTransform(dpr,…)` *before*
drawing, so `xOf`/`yOf` produce **CSS-pixel coordinates** that get
DPR-scaled by the transform during `arc()`. The hits[] array therefore
contains CSS-pixel coords directly, and the hover handler compares
raw `ev.clientX - rect.left` deltas without re-scaling.

This is documented in the `_wireVShapeTooltip` comment block. Test
asserts the absence of the DPR scaling pattern as a regression-guard.

### 2.7 Closest-point hit test (handles overlapping points)

Turn 156's draw order is `HOMO_1 → HOMO_2 → HET` so HET sits visually
on top. A naive "first hit wins" would resolve to whichever group
appeared first in iteration — wrong when an HET point overlays a
HOMO. `_vsHitTest` instead computes Euclidean distance and picks the
closest hit within `(r + 1.5)²` of the cursor. Tested with
deliberately overlapping points; correctly picks the closest.

---

## 3. What did NOT change

- **`_buildVShapeData`** — turn 156. Untouched.
- **`_drawVShapePlot`** — turn 156. Still returns `hits[]`. Untouched.
- **`_vShapeColor`** — turn 156. Reused by the tooltip's group-color
  swatch.
- **G-panel karyotype tab button (`gpKaryoVShapeBtn`)** — turn 156.
  Untouched.
- **Inheritance pill tooltip pattern** (turn 122). The V-shape tooltip
  uses different IDs (`#vsPointTooltip` vs `#inhPillTooltip`) and
  different state slots (`_vsLastRender` vs `state._inhPillHitRegions`)
  — no collision.
- **All inheritance pipeline (compute, cache, G-panel tab, matrix
  popup)** — untouched.

---

## 4. Test status

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 156)   | 73,435  | 2840            | 60    |
| Post-turn-157 (current)  | 73,667  | 2926            | 61    |
| **Δ session**            | +232    | +86             | +1    |

(Sweep delta is +86 not +61 because the test file count adds 1 and
some sweeps include warmup/identity assertions that resolve PASS
during the meta-pattern checks.)

Full sweep: **2926 / 0**. JS syntax: clean. HTML parser: 0 errors.

`tests/test_turn157_vshape_tooltip.js` (61 tests) covers:

1. All 6 helpers + `_vsLastRender` declared, all 6 exposed on window
2. `_vsTooltipEnsureEl`: unique id, sandbox-safe, idempotent, body
   append, fixed positioning, pointer-events:none, z-index 10000
3. `_vsTooltipBuildHtml`: empty-on-bad-hit, header, group color chip,
   u + agreement to 4 decimals, semantic colors for coherence /
   stripe_quality
4. `_vsTooltipShow`: default offset, edge-clipping flips, viewport
   clamp
5. `_wireVShapeTooltip`: idempotent dataset marker, reads from
   `_vsLastRender`, **CSS-pixel coords** (regex asserts the absence of
   the inheritance-pill DPR-scaling pattern as a regression guard),
   cursor swap, mouseleave handler
6. `openVShapePlot` integration (snapshot includes hits/data/candidate_id;
   wire call in try/catch; close also hides tooltip)
7. Sandboxed `_vsHitTest`: null/empty/non-array hits, direct hit,
   fudge-radius hit + miss boundary, overlapping points → closest
   wins, default radius
8. Sandboxed `_vsTooltipBuildHtml`: bad-hit returns empty, full point
   renders all fields, optional fields omitted gracefully, NaN/undefined
   render "NA"
9. Existing flow preserved (turn 156 V-shape, turn 122 inheritance
   pill, turn 152 G-panel tab, turn 155 cache key)

No existing test needed inverting.

---

## 5. What Quentin can do next session

Same workflow as turn 156, plus the new hover affordance:

1. Page 1 → lock + promote candidate
2. Page 2 → dosage heatmap → "Compute stripe quality"
3. G-panel (`g`) → karyotype tab → 🔍 V-shape diagnostic
4. **NEW:** hover any point in the V-shape plot → tooltip shows CGA,
   group, u, agreement, coherence class, stripe-quality tier

Useful for spot-checking individual outliers — e.g. an HET sample
with agreement 0.45 and `stripe_quality=junk` is worth investigating.
The tier color in the tooltip mirrors the colors used in the
karyotype-side per-sample regime panel, so visual identity is
preserved across surfaces.

---

## 6. Things I almost broke and fixed

- **DPR scaling regression risk.** The inheritance pill tooltip uses
  `canvas.width / rect.width` to scale mouse coords up to canvas-internal
  space. Reusing that pattern verbatim in V-shape land would have
  double-scaled — both the canvas's `setTransform(dpr,…)` and the
  manual scaling would have applied. Caught at design time. The
  `_wireVShapeTooltip` comment block documents the difference; the
  test asserts the absence of the bad pattern as a regression guard.
- **First-hit-wins on overlapping points.** Originally tempted to
  short-circuit on the first hit found. Caught while writing the
  draw-order test in turn 156 — HET points sit on top of HOMO points
  in the dip region, but iteration order doesn't match draw order.
  Fixed by switching to closest-by-distance picking with explicit
  test (`_vsHitTest` 7d).
- **Orphan tooltip on close.** First version of `closeVShapePlot`
  only set `modal.style.display = 'none'`. The tooltip itself
  (separate body-level element) would persist if the user closed
  while hovering. Fixed by adding `_vsTooltipHide()` to
  `closeVShapePlot`.
- **Module-level `let _vsLastRender`.** Originally drafted as a
  property on `state` (e.g. `state._vsLastRender`). Reverted to
  module-level because the tooltip is a UI-only ephemeral that
  doesn't need to persist across F5 — and putting it on `state`
  would invite serialization issues (the hits array contains
  references to non-serializable canvas pixel objects).

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 157 patched, 73,667 lines.
- `tests/test_turn157_vshape_tooltip.js` — NEW, 61 tests.
- `HANDOFF_2026-05-05_turn157_vshape_tooltip.md` — this file.

Plus prior handoffs (carried for reference):
- `HANDOFF_2026-05-05_turn156_vshape_diagnostic.md`
- `HANDOFF_2026-05-05_turn155_threshold_in_cache_key.md`
- `HANDOFF_2026-05-05_turn154_compute_ux_hardening.md`
- `HANDOFF_2026-05-05_turn153_inheritance_auto_register.md`
- `HANDOFF_2026-05-05_turn152_g_panel_inheritance_slice3.md`

---

## 8. Honest framing

**What's solid:**
- Pure read-only addition. No new state slots in `state`. No
  mutation of any existing pipeline.
- Idempotent canvas wiring via `dataset.vsTooltipWired` — safe to
  re-call.
- Closest-point hit test handles overlapping points correctly
  (verified by sandboxed test).
- Different ID + state-slot from the existing inheritance-pill tooltip
  — no collision.

**What's NOT done (and why that's right):**
- **No click handler / drill-in.** Hover-only for v0. A future turn
  could wire `click` on a hit to (e.g.) jump to that sample's row
  in the karyotype-tier-table, or open a per-sample dosage track.
  Not in v0 because per-sample drill-in is most useful from
  surfaces that already have it (karyotype tier table on page 2).
- **No keyboard navigation.** Tab + arrow keys to step between
  points would make the diagnostic accessible. Not in v0.
- **No "lock" mode.** Right now the tooltip follows the cursor; if
  the user wants to read it carefully without moving the mouse,
  they have to keep it stable. A click-to-pin variant could come
  later.

**What's queued:**
- **H** — Cross-candidate V-shape gallery / sparkline matrix
- **F** — Per-group "show fish" expand toggle in inheritance card
- **A** — UV refactor on locked groups
- **C** — Per-sample-lines het coloring
- **D** — Pivot based on what het / inheritance / V-shape reveal

End of handoff.
