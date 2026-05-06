# HANDOFF — turn 132 page-12 (θπ scrubber) frontend mirror, 7-slice build

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (65,507 lines, 856/0 tests passing)
**Working dir**: `/home/claude/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC (account `lt200308`).
**Supersedes**: `HANDOFF_2026-05-05_turn131_FINAL.md` (turn 131 spec-pass).

This handoff documents the full multi-turn build of page 12 (the
top-bar `2 local PCA θπ` tab, DOM `#page12`) from "empty-state scaffold"
to "6-of-8 panel renderers shipped." Each turn corresponds to one
slice; together they implement the page-1-frontend-mirror request from
the user.

---

## 0. Cohort discipline (NEVER conflate)

Same as turn 131 handoff:
1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok). Never invent surname.

---

## 1. The user request that triggered this build

> "Maybe the theta pi frontend (page 2) similarly to page 1 but with a
> slightly different hero panel as we discussed in another chat"

After clarification:
> "Its just that theta pi is another page index but in the top bar its
> the button 2."

So: **top-bar button "2"** → `data-page="page12"` → DOM id `page12` →
SCHEMA §22 θπ scrubber. NOT the candidate-focus deep-dive page.

**Hero panel** decision: per-carrier CUSUM strip + boundary-distribution
histogram (chat `487c7f04` turn 125). The slightly-different-from-page-1
panel — page 1 doesn't have CUSUM because dosage breakpoints are sharp
steps; θπ breakpoints are diffuse gradients, so per-carrier changepoint
distribution becomes the headline visualization.

---

## 2. The 7-slice build (turns A through G)

### Turn A — Slice 1: DOM scaffold

**Atlas changes**:
- Added `cusum_theta` row to `#thetaPiLayerStatus` grid (4th layer row,
  R-side `STEP_T05`)
- Inserted 9 mirrored panels inside `#page12` after the existing
  empty-state block, in correct top-to-bottom order:
  - `#thCtrlBar` (flex layout) — scrubber + play + viewmode
  - `#thCusumHeroPanel` — hero (per-carrier strip + cohort histogram +
    spread/concord badges) — *the slightly-different-from-page-1 panel*
  - `#thSimPanel` — sim_mat heatmap mirror
  - `#thZPanel` — diversity-|Z| waveform mirror
  - `#thLinesPanel` — per-sample θπ-PC1 lines mirror
  - `#thAnchorStripPanel` — envelope strip mirror
  - `#thPcaPanel` — θπ-PC1 vs PC2 K-means scatter mirror
  - `#thTrackedSamplesPanelCompact` — tracked-samples mirror (DOM only)
  - `#thL3Panel` — L3 contingency mirror (DOM only)
- Updated panel-layout description grid in the empty-state to document
  all 9 panels including the hero
- All panels start `display: none`; empty-state owns the page until
  layers load

**Tests**: `tests/test_turn132_page12_panel_mirror.js` — 78 / 0
**Backups**: `.bak_pre_page12_panel_mirror`, `.bak_post_page12_panel_mirror`

### Turn B — Slice 2: visibility wiring

**Atlas changes**:
- Added `_refreshThetaPiPanelVisibility()` adjacent to
  `_refreshThetaPiLayerStatus()` (~50 LOC)
- Layer→panel mapping:
  - `theta_pi_per_window` → `#thLinesPanel`, `#thTrackedSamplesPanelCompact`
  - `theta_pi_local_pca` → `#thSimPanel`, `#thZPanel`, `#thPcaPanel`, `#thL3Panel`
  - `theta_pi_envelopes` → `#thAnchorStripPanel`
  - `cusum_theta` → `#thCusumHeroPanel`
  - any θπ layer → `#thCtrlBar` shows (`flex`), `#thetaPiEmpty` hides
- Fires in onDataLoad refresh path right after `_refreshThetaPiLayerStatus()`
- Idempotent, reactive (layers removed → panels re-hide), defensive
  (no `state.layersPresent` / missing DOM / vm sandbox all guarded)

**Tests**: `tests/test_turn132_page12_visibility_wiring.js` — 50 / 0
**Backups**: `.bak_post_page12_visibility_wiring`

### Turn C — Slice 3: CUSUM hero renderer

**Atlas changes**:
- **JSON contract for `cusum_theta`** documented as a comment block in
  the atlas (since `SCHEMA_V2.md` isn't bundled). Carries only
  `persample[]` per-carrier observations — matches the walked-back
  `lib_persample_cusum.R` design from chat `487c7f04` ("observe first
  empirically, no assumed distributional shape"). Atlas computes
  display-only IQR-based spread classification descriptively.
- **`_drawThCusumHero()`** (~270 LOC). Paints two canvases + two
  badges:
  - **Lane 1 (`#thCusumStripCanvas`)** — per-carrier rug strip, sorted
    by `cp_bp` ascending, ticks colored by karyotype (HOM_REF blue /
    HET orange / HOM_INV purple / null grey), tick height encodes
    `strength`, draw direction encodes `asymmetry` sign, light Mb
    gridlines.
  - **Lane 2 (`#thCusumHistCanvas`)** — cohort changepoint histogram,
    ~√n bins clamped to [20,100], median line + IQR shading in amber.
  - **Spread badge** — IQR-based descriptive: `tight` <100 kb / 
    `intermediate` 100–500 kb / `ragged` >500 kb / `spread —` if <4
    carriers. Display threshold framing only — NOT a statistical claim.
  - **Concord badge** — reads `state.data.cusum_concordance`; shows
    `concord N/M · X%` when present, `concord —` until STEP_DC06 ships.

**Tests**: `tests/test_turn132_page12_cusum_hero_renderer.js` — 39 / 0
**Backups**: `.bak_post_page12_cusum_hero`

### Turn D — Slice 4: R-side `STEP_T05` spec

**No atlas code changes.** Wrote
`specs_new_turn131/SPEC_STEP_T05_theta_cusum_emitter.md` (357 lines).
Closes the loop on the hero: atlas-side renderer is shipped; R-side
now has a contract to build against.

10 sections: why this layer exists, inputs (theta_pi_per_window +
per-candidate karyotype TSV), algorithm (wraps the walked-back
`lib_persample_cusum.R`), output JSON shape, edge cases, repro pinning,
self-test fixture, companion specs to write later, atlas integration
already in place, build estimate.

Updated `README.md` and `SPECS_TIER_INDEX.md` in `specs_new_turn131/`.

### Turn E — Slice 5: per-sample θπ lines renderer

**Atlas changes**:
- **`_drawThLinesPanel()`** (~225 LOC). Minimum-viable mirror of page 1's
  lines panel. Honest scope decision documented inline:
  - Single source (θπ), single canvas — page 1 stacks PC1/PC2/GHSL/het
  - No lasso, no caching, no K-means coloring
  - No multi-source stacking
- **Dual-shape support**: Shape B (`data.theta_pi_per_window.{samples,
  windows, values}`) and Shape A (legacy `data.windows[w].theta`
  arrays) — matches existing detector tolerance.
- **Coloring cascade**: `state.data.cluster_labels_theta.labels` →
  `state.candidate.locked_labels` → neutral grey. Karyotype palette
  matches hero strip (HOM_REF blue / HET orange / HOM_INV purple).
  Two-pass draw so grouped lines paint over ungrouped greys.
- NaN-safe polylines (gaps break, don't bridge).
- Auto-creates canvas if `#thLinesCanvasContainer` empty; reuses
  existing canvas otherwise (idempotent — no DOM-node leaks).

**Tests**: `tests/test_turn132_page12_lines_renderer.js` — 38 / 0
**Backups**: `.bak_post_page12_lines_renderer`

### Turn F — Slice 6: paired sim_mat + |Z| renderers

**Atlas changes**:
- **`_drawThSimMatPanel()`** (~135 LOC). Reads
  `state.data.theta_pi_local_pca.sim_mat` (flat row-major) +
  `n_windows_thumb` (or √length fallback). Centered square geometry.
  Reuses page-1's `simColor` palette when globally available so the
  two pages stay visually consistent. Diagonal painted yellow for
  orientation. Persists `state._thSimGeom` for future click handlers.
- **`_drawThZPanel()`** (~110 LOC). Reads
  `state.data.theta_pi_local_pca.z`. Plots `|Z|` waveform vs Mb when
  `state.data.windows[].center_mb` available, falls back to raw window
  indices. Mb gridlines, dimmed `|Z|=2` reference line, NaN-safe
  polyline. Layout pad mirrors page-1 `drawZ`.
- Both vm-safe, idempotent, accept Float32Array typed inputs.
- **Documented contract** for `theta_pi_local_pca` inline in sim
  renderer comment block: required `z`, `sim_mat`, `n_windows_thumb`;
  optional `pc1`/`pc2`/`lam1`/`lam2`/`ratio`/`mds_coords` for future
  slices.
- Explicit deferral list in both renderer comments: L1/L2 overlays,
  click-to-jump, PDF-style triangle split, candidate strip, settings
  popover.

**Tests**: `tests/test_turn132_page12_simmat_z_renderers.js` — 50 / 0
**Backups**: `.bak_post_page12_sim_z_renderers`

### Turn G — Slice 7: paired anchor strip + PCA scatter renderers

**Atlas changes**:
- **`_drawThAnchorStripPanel()`** (~95 LOC). Reads
  `state.data.theta_pi_envelopes.{l1, l2}`. Paints L1 as blue
  rectangles on top half / L2 as green rectangles on bottom half.
  X-axis aligned to `state.data.windows[]` Mb range when available
  (so it stays visually-aligned with `#thLinesPanel`, `#thZPanel`),
  falls back to envelope extents. Skips malformed envelope rows.
  Side L1/L2 labels.
  - **Naming-clash documented inline**: page 1's `#anchorStripPanel`
    is anchor-concord; page 12's same-named panel is envelope geometry
    per the panel-layout grid. Naming kept for visual parity, semantics
    differ deliberately.
- **`_drawThPcaPanel()`** (~115 LOC). Reads
  `state.data.theta_pi_local_pca.{pc1, pc2}` (collapsed view) with
  optional per-window override via `pc1_by_window`/`pc2_by_window`
  when `state.cur` valid. Coloring cascade matches lines panel +
  hero strip. Two-pass draw (grey under, grouped on top). Frame,
  axis labels (rotated PC2 label using save/translate/rotate), 8%
  bound padding mirroring page-1.
- Both vm-safe, idempotent, accept Float32Array, NaN-safe.
- Explicit deferral list in PCA: window-trail, per-band legend
  overlay, sign-flip rules, lasso/hit-testing, per-cluster ellipses.

**Tests**: `tests/test_turn132_page12_anchor_pca_renderers.js` — 60 / 0
**Backups**: `.bak_post_page12_anchor_pca_renderers`

---

## 3. Cumulative test totals

| Turn | New tests | Cumulative |
|---|---|---|
| Baseline (turn 131) | — | 541 / 0 |
| A — DOM scaffold | +78 | 619 / 0 |
| B — visibility wiring | +50 | 669 / 0 |
| C — CUSUM hero | +39 | 708 / 0 |
| D — R-side spec | +0 | 708 / 0 |
| E — lines renderer | +38 | 746 / 0 |
| F — sim_mat + Z | +50 | 796 / 0 |
| G — anchor + PCA | +60 | 856 / 0 |

**Final: 856 / 0** across 53 working-format test files. Zero regressions.

Note: there are ~34 older tests in `tests/` with broken `require()`
calls from pre-turn-130 infrastructure (e.g. `test_atlas_dock_tabs.js`,
`test_atlas_export.js`). These have been broken since well before this
turn; the 856 count reflects working-format suites only. The handoff
from turn 131 reported "1900/0" — that was inferred from a different
runner pattern; the real working-format count was always 541, and is
now 856.

---

## 4. Atlas growth

| Turn | Atlas LOC | Δ |
|---|---|---|
| Baseline | 64,243 | — |
| A | 64,384 | +141 |
| B | 64,441 | +57 |
| C | 64,714 | +273 |
| D | 64,714 | 0 (specs only) |
| E | 64,954 | +240 |
| F | 65,246 | +292 |
| G | 65,507 | +261 |

**Net: +1,264 LOC** of page-12 frontend infrastructure across 7 slices.

---

## 5. State of page 12 after turn G

| Layer | Status | Source location |
|---|---|---|
| Empty-state docs + 4-row layer-status grid | ✅ | `<div id="page12">` |
| 9 mirrored panel DOM | ✅ | inside `#page12` |
| Visibility wiring | ✅ | `_refreshThetaPiPanelVisibility` |
| Hero CUSUM renderer | ✅ | `_drawThCusumHero` |
| R-side `STEP_T05` spec | ✅ | `specs_new_turn131/SPEC_STEP_T05_theta_cusum_emitter.md` |
| Per-sample θπ lines renderer | ✅ | `_drawThLinesPanel` |
| sim_mat heatmap renderer | ✅ | `_drawThSimMatPanel` |
| \|Z\| waveform renderer | ✅ | `_drawThZPanel` |
| L1/L2 envelope strip renderer | ✅ | `_drawThAnchorStripPanel` |
| θπ-PC1 vs PC2 scatter renderer | ✅ | `_drawThPcaPanel` |
| Tracked-samples renderer | ⚪ deferred | DOM only |
| L3 contingency renderer | ⚪ deferred | DOM only |
| Ctrl bar interactions | ⚪ deferred | DOM only |

**6 of 8 target renderers shipped.** The remaining two need a
state-machine slice (page-12 cur/tracked/anchor-window state), not
another renderer slice.

---

## 6. State / window slots added through turn G

```
state._thSimGeom               (Slice 6)  — sim_mat geometry for click handlers
                                            (matches page-1 state._simGeom convention)
```

```
window-exposed functions added (or implicitly defined at top level):
_refreshThetaPiPanelVisibility   (Slice 2)
_drawThCusumHero                 (Slice 3)
_drawThLinesPanel                (Slice 5)
_drawThSimMatPanel               (Slice 6)
_drawThZPanel                    (Slice 6)
_drawThAnchorStripPanel          (Slice 7)
_drawThPcaPanel                  (Slice 7)
```

All 7 fire in the onDataLoad refresh path right after
`_refreshThetaPiLayerStatus()`. None expose to `window.X` explicitly
(matches existing pattern — `typeof X === 'function'` guards).

---

## 7. Backups present

```
Inversion_atlas.html.bak_pre_page12_panel_mirror
Inversion_atlas.html.bak_post_page12_panel_mirror
Inversion_atlas.html.bak_post_page12_visibility_wiring
Inversion_atlas.html.bak_post_page12_cusum_hero
Inversion_atlas.html.bak_post_page12_lines_renderer
Inversion_atlas.html.bak_post_page12_sim_z_renderers
Inversion_atlas.html.bak_post_page12_anchor_pca_renderers
```

`.bak_*` files NOT in bundle to keep size down. Re-derivable from git
(or by re-running the slices).

---

## 8. JSON contracts pinned this turn

These are atlas-side contracts the renderers consume. Kept as inline
comment blocks in `Inversion_atlas.html` since `SCHEMA_V2.md` isn't
bundled. Make these authoritative if/when SCHEMA_V2 gets a §22 update.

### `cusum_theta`
See `specs_new_turn131/SPEC_STEP_T05_theta_cusum_emitter.md` §4 and
the atlas inline contract above `function _drawThCusumHero()`.
Per-carrier observations only (no aggregated modes).

### `theta_pi_local_pca`
Atlas inline contract above `function _drawThSimMatPanel()`:
```
{
  z:                [<n_windows>]      // diversity-|Z| per window
  sim_mat:          [<n×n flat>]       // window × window θπ-similarity
  n_windows_thumb:  <int>              // n; falls back to sqrt(sim_mat.length)
  pc1, pc2:         [<n_samples>]      // collapsed-view per-sample PCs (Slice 7)
  // optional:
  pc1_by_window, pc2_by_window:        // per-window override (window × n_samples)
  lam1, lam2, ratio, mds_coords        // not yet consumed
}
```

### `theta_pi_envelopes`
Atlas detector + Slice 7 anchor-strip consumer. Row shape matches
dosage envelope rows: `{ start_bp, end_bp }` minimum.

### `theta_pi_per_window`
Two valid shapes per the existing detector (dual-shape tolerant):
- **Shape B (current spec)**: top-level layer with `samples[]`,
  `windows[]`, `values[]` flat row-major
- **Shape A (legacy)**: `data.windows[w].theta = [<n_samples>]` per
  window (already in production via `export_precomp_to_json.R --theta`)

---

## 9. What's blocked on the R-side

Page 12 is fully scaffolded but largely empty until R-side ships:

| R-side script | Produces layer | Blocks atlas panels |
|---|---|---|
| `STEP_R39` (planned) | `theta_pi_per_window` | lines, tracked-samples |
| `STEP_R40` (planned) | `theta_pi_local_pca` | sim_mat, \|Z\|, PCA, L3 |
| `STEP_R41` (planned) | `theta_pi_envelopes` | anchor strip, ctrl bar viewmode |
| `STEP_T05` (spec'd this turn) | `cusum_theta` | hero panel |

When ANY of these land, the corresponding panel(s) reveal and paint
immediately with no further atlas changes. The hero is the most-polished
panel; once R-side ships `cusum_theta`, page 12 becomes useful even
before the other layers.

---

## 10. Where to start the next chat

### Option 10a — Finish page 12 (Slice 8: state machine + remaining renderers)

Build the page-12 state machine that the remaining 2 renderers need:
- `state.thCur` — current θπ-window scrubber position
- `state.thTracked` — tracked-samples set for page 12
- `state.thAnchorWin` — anchor window for concord
- `state.thViewMode` — genome / L1 zoom / L2 zoom

Then ship:
- Tracked-samples renderer for `#thTrackedSamplesPanelCompact`
- L3 contingency renderer for `#thL3Panel`
- Ctrl bar interactions (scrubber drives `state.thCur`, viewmode
  drives x-range)

Estimate: 2–3 turns. Expensive because it's not just rendering — it's
new control state.

### Option 10b — Pivot back to turn-131 priority queue

Per turn-131's `SPECS_TIER_INDEX.md` recommended build order:
1. **L2-sweep auto-promote Slice 1** — Tier 1 producer; biggest unlock
2. **G-panel scaffold Slice 1** — Tier 2 review surface, needs
   `_rebuildCandidateRegistries()` pre-fix first
3. **Trajectory matrix viewer Slice 3**
4. **Cross-chromosome lineages Slice 1**

Page 12 is "good enough" at 6/8 renderers — the remaining two need
control state most users don't touch from page 12 anyway. The L2-sweep
work has been on the critical path since turn 130.

### Option 10c — R-side build

The spec for `STEP_T05` is shipped. If Quentin wants to actually see
real cusum data in the page-12 hero, the R-side build is unblocked.
Spec estimate (from §10): ~1 LANTA session, gated only on `STEP_R39`
(`theta_pi_per_window`) shipping first.

### Recommendation

10b. Page 12 is in a usable state now (it'll paint as soon as R-side
data arrives), and the L2-sweep work has been the critical-path item
for two turns. Coming back to Slice 8 after L2-sweep ships is fine —
the page-12 state machine isn't blocking any other atlas work.

---

## 11. Honest framing

**What turn 132 actually delivered:**
- A faithful frontend mirror of page 1 on page 12, preserving the
  visual/structural pattern users already know
- One genuinely new visualization (the CUSUM hero) that page 1
  cannot produce because dosage breakpoints don't have the same
  geometry as θπ breakpoints
- A clean R-side contract that closes the loop on the hero

**What it deliberately didn't deliver:**
- Page-1's elaborate control state machinery (lasso, tracked samples,
  caching, multi-source stacking, K-means). Page 12 doesn't need most
  of these because the science driving the page is "paint the data and
  see," not "interactively navigate a candidate." When/if those
  interactions become valuable, they get their own slice.

**What this means for the manuscript:**
- The CUSUM hero is the figure that goes in the methods description
  for the boundary-distribution observations (chat `487c7f04` framing:
  *"the 5' breakpoint resolves to 14.01 Mb (IQR 70 kb across 60
  carriers), consistent with recent founder-locked transmission"*)
- The mirror panels (lines, sim_mat, |Z|, PCA, anchor strip) provide
  orthogonal-validation visuals that pair 1:1 with page-1 figures —
  use both pages side-by-side to demonstrate concordance / discordance

---

## 12. Bundle contents

`/mnt/user-data/outputs/Atlas_full_bundle_2026-05-05_turn132.tar.gz`

Contains:
- `Inversion_atlas.html` (current, 65,507 lines)
- `tests/` (all *.js, *.py)
- `specs_todo/` (active build queue, unchanged from turn 131)
- `specs_new_turn131/` (pending review queue, +1 spec from this turn:
  `SPEC_STEP_T05_theta_cusum_emitter.md`)
- All previous handoffs (kept for history)
- This handoff (`HANDOFF_2026-05-05_turn132_FINAL.md`)
- `OBSERVATIONS_TO_FIX.txt`

`.bak_*` files NOT in bundle.

Walk the map carefully, respect cohort discipline, don't break the
test suite. Page 12 is in a clean stopping state; pick up wherever
makes most sense.
