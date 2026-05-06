# HANDOFF — atlas spec-implementation session, 2026-05-03

## Where we are

Three turns shipped this session, working through `specs_todo/`. Turn 117
(focal-vs-background widget) is **partially shipped** — module fully built and
tested, atlas integration ~30% complete. Resume there.

| Turn | Spec | Status |
|---|---|---|
| 115 | `SPEC_cross_species_dotplot.md` | ✅ **SHIPPED** — full atlas integration, end-to-end tested |
| 116 | `SPEC_ncrna_density_layer.md` | ✅ **SHIPPED (MVP)** — storage layer + MVP renderer + R aggregator. Y-mode dropdown / view toggle / pan-zoom deferred to a polish turn |
| 117 | `SPEC_focal_vs_background_widget.md` | 🟡 **PARTIAL** — module written + 9-test suite passing + script tag + 2 mount divs. Renderer functions + radius binding + auto-discovery wiring NOT YET DONE |

Two specs untouched: `SPEC_ld_decay_overlay.md`, `SPEC_xpehh_track.md`. The
XP-EHH spec was moved to `specs_todo/later/` per user direction (blocked on
phased VCF + reference cohort decision; not to be pulled in next session).

## Atlas baseline tracking

| File | Lines | What changed |
|---|---|---|
| `Inversion_atlas.html` (turn 116 baseline) | 57,544 | +450 from turn 115 (ncRNA density layer) |
| `Inversion_atlas.html` (now, after partial t117) | 57,556 | +12 (two HTML mount divs + script tag for fvb module) |

The atlas remains in a clean, fail-soft state — even though the focal-vs-bg
renderer functions are not wired, the existing try/catch pattern around all
panel renderers ensures the missing wiring causes ZERO runtime errors. The
mount divs simply stay empty (`display: none` defaults). Safe to ship as-is.

`node --check` on the entire concatenated inline JS: clean.

## What's left in turn 117

The module (`atlas_focal_vs_bg.js`) is complete and tested. To finish
integration, the next session needs to add **6 small things** to the atlas:

### 1. Page-16 renderer — `_renderCrossSpeciesFocalVsBg()`

Sibling to `_renderCrossSpeciesDotplot()`. Mounts the widget into
`#csFocalVsBg`. Reads from `state.crossSpecies.breakpoints[]`, finds the
active breakpoint via `state._crossSpeciesUI.active_id`, builds focal_lo/hi
from the breakpoint's gar/mac coordinates ± `state._focalVsBg.csRadiusBp`
(default 2 000 000). Cohort_selector = null on day 1.

Insert next to `_renderCrossSpeciesDotplot` (currently at line ~20389 in the
turn-117 atlas). Wire the call into `_renderCrossSpeciesPage` right after
the dotplot try/catch block (currently at line ~18965).

### 2. Page-11 renderer — `_renderBndFocalVsBg()`

Sibling to `_renderNcRNADensityPanel`. Mounts the widget into `#bndFocalVsBg`.
Uses `_ensureBoundariesState().active_cand_id` + `_bndFindCandidate(id)` to
get the candidate, then `_boundaryScanRange(cand, scan_radius_bp, chromLen)`
to compute focal_lo/hi.

**Critical**: page 11 widget needs to **share** scan_radius_bp with the
parent boundaries page. Pattern:

```javascript
const widget = window.popgenFocalVsBg.makeFocalVsBgPanel({
  page_id: 'page11',
  ...
  radius_source: 'boundariesState',
  get_radius_bp: () => _ensureBoundariesState().scan_radius_bp,
  set_radius_bp: (r) => {
    const bs = _ensureBoundariesState();
    bs.scan_radius_bp = r;
    if (bs.active_cand_id != null) bs.cache.delete(bs.active_cand_id);
  },
  on_radius_change: (r) => _bndRefreshUI(),  // propagate to whole page
});
```

Wire the call into `_bndRefreshUI` next to the existing ncRNA renderer call.

### 3. Page-16 state init

```javascript
state._focalVsBg = state._focalVsBg || { csRadiusBp: 2000000 };
```

Anywhere top-level in the atlas init code is fine.

### 4. Cs-active-id change re-render

When the user clicks a different breakpoint in the cs catalogue, the
focal-vs-bg widget on page 16 needs to re-render. The catalogue click handler
already calls `_renderCrossSpeciesPage()` which calls `_renderCrossSpeciesFocus()`
— just hook `_renderCrossSpeciesFocalVsBg()` after the focus render.

### 5. Help-page table row

```html
<tr><td><b>17 focal-vs-bg widget</b></td><td>Spalax-style focal-vs-background
distribution comparison. Renders on page 11 (boundaries — focal = candidate
scan range) and page 16 (cross-species — focal = breakpoint ±2 Mb). Wilcoxon
rank-sum P value (live, in-browser, normal approximation w/ continuity
correction). Metric pool auto-discovers from loaded layers: |Z| (always),
F_ST/Hobs/θπ (when popstats live server hooked up), tRNA_all/rRNA_all/rRNA_5S/
ncRNA_all (when ncRNA density layer loaded), all_TE (when TE density layer
loaded). Z-transform toggle on by default. Bg scope: per-chrom (default) or
genome-wide-excl-focal. View mode: full or zoomed-to-scan-range.</td></tr>
```

### 6. Tarball delivery + CHANGELOG_turn117.md

Pattern: see `CHANGELOG_turn115.md` and `CHANGELOG_turn116.md`.

## Module file: `atlas_focal_vs_bg.js`

976 lines, UMD module. Public API:

```javascript
window.popgenFocalVsBg = {
  makeFocalVsBgPanel(opts) -> HTMLElement,
  wilcoxonRankSum(focal, bg) -> { W, U_focal, p_two_sided, logp_neg, z_stat, n_focal, n_bg },
  availableMetrics(state) -> [{id, label, family, getData}],
  computeFocalBgStats(opts) -> { focal_arr, bg_arr, n_focal, n_bg, focal_mean, bg_mean, ..., wilcoxon, scatter_xy },
};
```

**makeFocalVsBgPanel opts shape**:
- `page_id`: string, used as the localStorage prefs key — `'page11'` or `'page16'`
- `anchor_id`, `anchor_type`: identifiers (string + 'cs_breakpoint' | 'candidate_boundary')
- `chrom`: string
- `focal_lo_bp`, `focal_hi_bp`: integers
- `anchor_marks_mb`: array of Mb positions to draw as dashed red anchor lines
- `radius_source`: `'self'` | `'boundariesState'`
- `radius_default_bp`: number (used when `radius_source === 'self'`)
- `get_radius_bp`, `set_radius_bp`: callbacks (used when `radius_source === 'boundariesState'`)
- `on_radius_change`: callback fired when widget radius changes (host should re-render)
- `cohort_selector`: null on day 1 (forward-compat)
- `containerStyle`: optional inline CSS string

**Critical design refinement on the spec** discovered during build:

Only the `z` track in `POPSTATS_TRACKS` has a real `getData` — the others
(`fst_hom1_hom2`, `hobs_hexp`, `theta_invgt`) are placeholders depending on a
live popstats server. So the widget's `availableMetrics()` filters out
null-returning getters at render time. Day-1 metric pool = whatever's actually
loaded:
- `z` (always available from precomp)
- `fst_hom1_hom2`, `hobs_hexp_ratio`, `theta_pi_log2_ratio` — if popstatsLive present
- `ncrna_tRNA_all`, `ncrna_rRNA_all`, `ncrna_rRNA_5S`, `ncrna_ncRNA_all` — auto-discovered if turn-116 layer loaded
- `te_all_TE` — auto-discovered if repeat-density layer loaded

The ncRNA auto-discovery is the **direct payoff from turn 116** — six new
metrics in the focal-vs-bg pool with zero atlas-side wiring.

## What's already verified

### Turn 115 (cross-species dot plot) — full integration tests
- `node --check` on full inline JS: clean
- End-to-end render with cs_breakpoints v2 fixture (4 synteny blocks, one
  inversion): module loads, renderer runs, slot becomes visible, panel cache
  hits on identical re-render, cache invalidates on data change

### Turn 116 (ncRNA density) — full integration tests
- Detector recognizes `tool === 'ncrna_density_v1'`
- Store + persist round-trips through localStorage; restore picks it back up
- Renderer produces 4.9 KB SVG with chrom label, three family chip rows,
  correct Y-axis class label
- Class switching reactively re-renders with the new class label
- Empty-state hint shows when chrom isn't loaded
- `_classifyJSONKind` ordering preserved: cross_species → dotplot_mashmap →
  ncrna_density → repeat_density → chromosome

### Turn 117 partial — JS module unit tests (9 tests, all passing)
- Wilcoxon clean separation matches R's `wilcox.test`: `c(1,2,3) vs c(4,5,6)` → P = 0.0809 (R: 0.0809), U = 0
- Tied ranks handled correctly via average-rank assignment + tie correction in variance
- Identical samples → P = 1.0
- Spalax-like extreme separation (80 focal at 0.6, 4000 bg at -0.005) → P floored at 1e-16, z = 15.34
- NaN/null skip works correctly
- Empty input handled (returns NaN P)
- P-value floor at 1e-16 working
- `computeFocalBgStats` end-to-end: partitions correctly (3 focal, 7 bg from
  10 windows), computes means + medians, runs Wilcoxon on raw values
- `availableMetrics(state)` filters correctly: empty state → 0, with `state.data.windows` populated → `['z']`, with `state.ncRNADensity[chrom]` populated → adds `'ncrna_tRNA_all'`, `'ncrna_rRNA_all'` etc.
- End-to-end render with DOM shim: module produces 16,854 chars of HTML
  including 3 SVGs (focal histogram + bg histogram + scatter)

### Bug found and fixed during testing
The Z-transform branch divided by `bg_sd_raw` when the background was
constant-valued (e.g., `[0.1, 0.1, 0.1, ...]`). Floating-point noise gave
`sd ≈ 1.5e-17` instead of exactly 0, so the `bgSd > 0` check passed and the
focal mean displayed as `5e16`. Fix: epsilon guard `bgSd > 1e-12`. Fixed
in module before delivery.

## Working directory snapshot

```
/home/claude/work/
├── atlas/                        # Original baseline (read-only ref): 56,902 lines
│   └── Inversion_atlas.html
├── build/                        # Working dir; what next session pulls from
│   ├── Inversion_atlas.html       # 57,556 lines (turn 115 + turn 116 + partial t117)
│   ├── atlas_dotplot.js          # 844 lines — turn 115 module
│   ├── atlas_focal_vs_bg.js      # 976 lines — turn 117 module (COMPLETE)
│   ├── aggregate_ncrna_density.R  # 226 lines — turn 116 R aggregator
│   ├── sample_ncrna_density_LG28.json # synthetic fixture for testing
│   ├── CHANGELOG_turn115.md
│   └── CHANGELOG_turn116.md
├── fvb/                          # Turn 117 module work dir (mirror of build/atlas_focal_vs_bg.js)
│   └── atlas_focal_vs_bg.js
├── dotplot/                      # Turn 115 module work dir (with demo HTML)
│   ├── atlas_dotplot.js
│   ├── dotplot_demo.html         # standalone preview, no atlas needed
│   └── README.md
└── specs_todo/
    ├── README.md
    ├── SPEC_cross_species_dotplot.md  # IMPLEMENTED in t115 (move to specs_shipped/ if you want)
    ├── SPEC_focal_vs_background_widget.md  # IN PROGRESS in t117
    ├── SPEC_ncrna_density_layer.md  # IMPLEMENTED MVP in t116 (polish deferred)
    ├── SPEC_ld_decay_overlay.md  # untouched
    └── later/
        ├── README.md             # explains why
        └── SPEC_xpehh_track.md   # MOVED HERE per user direction
```

## How to resume in the next chat

1. Pull the latest atlas + module files (delivered as a single tarball)
2. **First action: re-read `ATLAS_REFERENCE_for_phase4b_synthesis.md` and this
   handoff** to re-anchor on the project state
3. Confirm the atlas baseline is 57,556 lines via `wc -l`
4. Continue turn 117 by adding the 6 small things listed above. Suggest doing
   them in a single turn — they're all small and tightly coupled
5. After turn 117 ships: next candidate is `SPEC_ld_decay_overlay.md`. It
   sits naturally next to focal-vs-bg on both pages (page 11 below the
   widget; page 16 below the dot plot)

## How NOT to break things in the next chat

- The atlas is currently in a fail-soft state. Adding more code is safe.
  Removing or modifying existing wiring is where bugs hide.
- `state._focalVsBg` is **not yet initialized** anywhere in the atlas. Any
  code that reads it without `state._focalVsBg = state._focalVsBg || {...}`
  will throw on page 16 boot. Initialize it the FIRST thing the new renderer
  does.
- The widget reads `window.state` directly via the closure. Make sure the
  page-16 renderer passes the right context — i.e. when iterating over
  multiple breakpoints in a future turn, each panel needs its own widget
  instance, not a shared one.
- Never assume `state.crossSpecies.breakpoints[].gar_chr` matches
  `state.data.chrom`. Page 16 mounts the widget for cross-species data
  (Cgar↔Cmac alignment), but the per-window `z` track lives on the cohort
  data (226 fish, single haplotype). The widget should make this clear via
  the metric pool — `z` on page 16 is "cohort divergence at this position",
  NOT "cross-species divergence". Surface this in the help-page row.

## What's deferred (low-priority, not blocking next turn)

- **Mode A** (precomputed JSON layer with locked focal+bg stats) — atlas-side
  reader scaffold is in the module, but no JSON shape is currently produced
  by any pipeline script. Day-1 path is Mode B (live computation) only.
  Mode A becomes interesting when sample sizes get large enough that
  in-browser Wilcoxon is slow.
- **ncRNA panel polish** (Y-mode dropdown, view-mode toggle, pan-zoom) —
  cloned from the canonical TE panel when needed
- **Page-16 cross-species ★ chips** for tRNA-cluster / rDNA-flank flags —
  small Python-side patch to `STEP_CS01_extract_breakpoints.py` plus a few
  lines of catalogue rendering. Now unblocked by turn 116
- **LD decay overlay** (`SPEC_ld_decay_overlay.md`) — natural next-spec
  after turn 117 ships

## Key files to bring into the next chat

The tarball includes these. From most-important to nice-to-have:

1. **HANDOFF_2026-05-03_session_summary.md** ← this file
2. **Inversion_atlas.html** — the working baseline (57,556 lines)
3. **atlas_dotplot.js** — turn 115 module
4. **atlas_focal_vs_bg.js** — turn 117 module (this is the file to integrate)
5. **aggregate_ncrna_density.R** — turn 116 R aggregator
6. **sample_ncrna_density_LG28.json** — fixture for testing turn 116 layer
7. **CHANGELOG_turn115.md**, **CHANGELOG_turn116.md** — what shipped each turn
8. **SPEC_focal_vs_background_widget.md** — the spec turn 117 implements
9. **SPEC_ld_decay_overlay.md** — next-up spec
10. **SPEC_cross_species_dotplot.md**, **SPEC_ncrna_density_layer.md** — already-shipped specs (reference)
11. **ATLAS_REFERENCE_for_phase4b_synthesis.md** — atlas conceptual frame
12. **specs_todo/later/** — the deferred XP-EHH spec + README explaining why

## Communication notes

User communication style: terse, direct. "Continue" = proceed to next agreed
item, the "next agreed item" being whatever Claude announced last. User works
across multiple parallel chats; don't expect them to remember exact module
filenames — refer to specs by topic ("the dot plot", "the focal-vs-bg widget").

User pushed cancel by mistake mid-`str_replace` during turn 117 integration —
the partial state is correctly captured; the cancellation only affected the
next-up patch (page-16 renderer function), not anything that had already
been written.
