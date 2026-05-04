# SPEC — Focal-vs-background widget (Spalax-style figure idiom)

**Status**: forward-looking spec. Not implemented. Created turn after turn 114c-partial.
Revised same-turn after a closer read of existing boundaries-page infrastructure
(see §0 below — many components already exist and the widget should reuse them
verbatim, not reinvent).

**Reading order**: this spec → §0 of this spec (settled infrastructure inventory)
→ `SPEC_OVERVIEW_multispecies_architecture.md` (multi-species framing) →
`SPEC_comparative_te_breakpoint_fragility.md` (sibling spec, same architectural
pattern: existing per-window data + optional aggregated JSON override) → atlas
`_renderRepeatDensityPanel` (line ~16223) and `_boundaryScanRange` (line
~14861) and `POPSTATS_TRACKS` (line ~52020) — the rendering primitives this
widget composes.

**Reference figure**: Spalax mole-rat fusion paper, Figure 4 panels C/D
(focal-vs-background histogram + Wilcoxon P) and E/F (along-chromosome scatter
of metric vs distance-from-breakpoint with focal zone shaded). Panels A/B
(rfmix admixture) are deliberately omitted — Quentin's cohort is single-species
(226 *C. gariepinus* hatchery broodstock) and rfmix requires two cohorts.
Once species 2 lands and the multi-cohort group-picker is wired (§5.5), an
rfmix-equivalent panel can be added in a follow-on spec.

---

## 0. Settled infrastructure — REUSE these, don't reinvent

The boundaries page (page 11) and popstats page (page 8) already encode most
of the range / interval / track / data-source logic this widget needs. The
widget is largely a new *combination* of existing primitives, with one new
small piece (Wilcoxon JS) and a new layout. Below is the inventory of what
exists. **Anything in this list is contractual: the widget reuses it and does
not duplicate.**

### 0.1 Range / interval logic — **reuse `_boundaryScanRange`**

`_boundaryScanRange(cand, scan_radius_bp, chromLen)` (atlas line ~14861)
already resolves a candidate → focal interval `{ start_bp, end_bp, win_lo,
win_hi }`. It handles:
- chromosome-end clamping
- huge-candidate auto-expansion (`CANDIDATE_HUGE_BP = 3 Mb` → `radius =
  max(radius, span × CANDIDATE_HUGE_RATIO)`)
- window-index resolution via `_bsearchWin`

The widget calls this directly for page 11 mounts. For page 16, the widget
calls a thin sibling `_focalRangeFromCsBreakpoint(bp, focal_radius_bp,
chromLen)` that uses the breakpoint position (rather than a candidate
[start,end] span) but produces the same `{ start_bp, end_bp, win_lo, win_hi }`
shape.

### 0.2 Scan-radius toolbar — **reuse the exact button layout from page 11**

The boundaries page (atlas lines ~21770–21793) already has the toolbar
pattern:

```
[1 kb] [5 kb] [10 kb] [25 kb] [100 kb] | [1 Mb] [1.5 Mb] [2 Mb] [5 Mb]
```

Visual separator divides kb-tier from Mb-tier. **The widget reuses this exact
toolbar HTML pattern** (extracted into a shared helper `_radiusToolbarHTML(idPrefix)`).
On page 11 the widget *inherits* the boundaries page's current radius selection
(reads from `_ensureBoundariesState().scan_radius_bp` rather than maintaining
its own). On page 16 the widget has its own copy of the toolbar but reuses
the same HTML helper, defaulting to 2 Mb (Spalax-matching).

### 0.3 Defaults — **reuse `BOUNDARY_DEFAULTS`**

```javascript
const BOUNDARY_DEFAULTS = {
  SCAN_RADIUS_BP: 1500000,        // 1.5 Mb each side
  ZONE_RADIUS_WINDOWS: 5,
  CANDIDATE_HUGE_BP: 3000000,
  CANDIDATE_HUGE_RATIO: 0.5,
  ...
};
```

The widget reads from this constant. If a future tweak changes the boundary
default radius, the widget tracks automatically. **No widget-specific
"FOCAL_DEFAULT_RADIUS" constant** — that's reinvention.

### 0.4 Per-window data — **reuse `POPSTATS_TRACKS` + `_galleryColumnGetData`**

The atlas already declares per-window tracks (atlas line ~52020):

```javascript
const POPSTATS_TRACKS = [
  { id: 'fst_hom1_hom2', label: 'Fst Hom1-Hom2', ... },
  { id: 'hobs_hexp',     label: 'Hobs/Hexp',     ... },
  { id: 'theta_invgt',   label: 'θπ by invgt',   ... },
  { id: 'z',             label: 'Z',             ... },
  ...
];
```

`_galleryColumnGetData(trackDef)` (line ~52177) returns
`{ mb, values, refLine }` from `popstatsLive.lastResponse` or from raw column
fallback. **The widget calls these existing getters** to obtain per-window
metric values. No new data path.

The metric pool the widget exposes:

| Widget metric ID | Widget label | Source `POPSTATS_TRACKS.id` | Notes |
|---|---|---|---|
| `fst` | F_ST (Hom1–Hom2) | `fst_hom1_hom2` | F_ST between karyotype classes — semantically what Spalax calls F_ST between cohorts |
| `theta_pi_ratio` | log₂ θπ ratio (HOM_INV / HOM_REF) | `theta_invgt` (multi-line; widget computes log2 ratio of the two relevant lines) | π asymmetry signal — already used to flag LG28 |
| `hobs_hexp` | Hobs / Hexp | `hobs_hexp` (multi-line; widget computes ratio per window) | low ratio inside an inversion = expected |
| `z` | robust \|Z\| | `z` | sanity track — should peak in focal zone if discovery worked |

**dXY is NOT in the day-1 metric pool** because the existing `POPSTATS_TRACKS`
list does not include a dxy track today. (This was an oversight in the v1 of
this spec.) When dxy lands as a popstats track, it slots into the metric pool
without needing widget changes — the widget auto-discovers any track whose
ID is added to `POPSTATS_TRACKS` (or auto-discovered via `state.data.tracks`)
that has scalar per-window values.

### 0.5 Empty-state pattern — **mirror `_renderRepeatDensityPanel` empty state**

Atlas line ~16238 already establishes the "data not loaded for this chrom"
hint with a "loaded for: chr1, chr5..." inline list. The widget reuses this
exact pattern when the selected metric isn't loaded for the active chrom.

### 0.6 View-mode toggle — **mirror `viewMode: 'full' | 'zoomed'`**

`_renderRepeatDensityPanel` already implements full-chrom vs zoomed-to-scan
view modes (atlas line ~16265, ~16298). The widget's along-chromosome scatter
(panel E equivalent) inherits this idiom: a `viewMode` toolbar chip toggles
between full chromosome (Spalax-style, E/F panels) and zoomed-to-scan (more
useful for boundary-scale work on page 11).

### 0.7 Cursor sync — **reuse `_navigateToCandidate` and the click-to-jump from turn 114a**

Turn 114a wired cs-breakpoint click-to-jump into the lines panel. The
widget's scatter inherits the same hit-testing pattern: clicking a window
in the scatter dispatches the same focal-window event, populating the
per-sample lines panel.

### 0.8 Sticky preferences — **mirror `_repeatDensityPrefs` storage pattern**

`_repeatDensityPrefs()` returns localStorage-backed prefs (yMode, viewMode).
The widget uses an analogous `_focalVsBgPrefs(pageId)` returning
`{ metric, bg_scope, z, viewMode }`, stored under
`focal_vs_bg::<pageId>`. Per-page (page 11 vs page 16 keep separate prefs).

---

## 1. Purpose

One reusable widget, two page mounts on day 1, that asks:

> Is the distribution of metric M inside a focal genomic zone different from
> the distribution of M outside that zone (the background)?

The widget produces:
1. Two histograms — focal and background — overlaid mean lines (panels C+D)
2. A Wilcoxon rank-sum P value comparing them
3. An along-chromosome scatter of M vs genomic position with the focal zone
   shaded (panel E)

The widget is **agnostic to what defines the focal zone**. Different atlas
pages anchor it differently:

| Page | Focal zone definition | Source of focal interval |
|---|---|---|
| 16 cross-species | ±W around a Gar–Mac karyotype breakpoint (default W = 2 Mb) | new helper `_focalRangeFromCsBreakpoint(bp, W, chromLen)` |
| 11 boundaries | scan range from `_boundaryScanRange(cand, scan_radius_bp, chromLen)` | **reuses existing function verbatim** |
| 8 popstats *(future)* | any user-drawn interval | TBD when the popstats interval-drawing UI exists |

The widget is **single-cohort on day 1** (the 226 *C. gariepinus* cohort).
Forward-compatible architecture mandatory for day-2 multi-cohort selection
(see §5.5).

---

## 2. Why this matters

The Spalax figure is clean visual logic for "is something special happening
near this structural feature?" The atlas already shows per-window metrics as
tracks (popstats page), but the *focal-vs-background statistical contrast* is
not surfaced today. Without it, the manuscript figure has to be built
externally and re-checked every time the candidate list shifts. With it, the
figure exists live, the Wilcoxon P updates as the user moves between
candidates / breakpoints, and the manuscript number lands on a fixed JSON
file when the analysis freezes.

This is also the natural home for **F_ST between karyotype classes**: the
same widget that draws Spalax-style F_ST near a fusion in Spalax draws F_ST
between HOM_REF and HOM_INV near an inversion boundary in Gar. The math is
identical; the biological interpretation is "recombination suppression
between arrangements" instead of "post-fusion divergence between species" —
but the plot reads the same way and the manuscript figure tells a cleaner
story because the inversion is the only structural variant in play.

The boundaries page (page 11) already drives F_ST as one of its 12 weighted
edge tracks (`fst_edge` weight 0.05). The widget surfaces the same data in
its other natural form — a focal-vs-background distribution test — without
duplicating the data path.

---

## 3. What the widget DOES / does NOT claim

### DOES

- The distribution of M in the focal zone differs from background, with a
  Wilcoxon rank-sum P value
- A visual sense of where in the chromosome the deviation concentrates
  (panel E)
- Per-chromosome OR genome-wide background, user-toggleable

### Does NOT

- Causation — "the inversion caused the F_ST elevation" — the widget cannot
  distinguish recombination suppression from selection from background
  heterogeneity
- Independence between metrics — F_ST and θπ are correlated; running the
  widget on each separately and counting hits is multiple testing the user
  has to think about
- That the focal zone is correctly placed — that's upstream (boundary
  refinement on page 11, breakpoint extraction on page 16)

---

## 4. Vocabulary discipline

| Allowed | Discouraged |
|---|---|
| "Z-F_ST is elevated in the focal zone (P = …)" | "F_ST shows the inversion is under selection" |
| "Focal-zone metric distribution differs from background" | "The breakpoint is selected" |
| "Recombination suppression between arrangements is consistent with elevated F_ST near boundaries" | "The inversion is causing differentiation" |
| "Wilcoxon P (rank-sum, two-sided) = …" | "highly significant" without P value |

P values displayed with explicit floor at `< 2.2e-16` (Spalax convention).

---

## 5. Architecture

### 5.1 Data model — three modes

The widget supports three data sources, in priority order. The first
available wins:

#### Mode A — precomputed aggregated JSON (preferred when present)

A new optional JSON layer per candidate / breakpoint, with focal+bg stats
locked at HPC-build time:

```jsonc
// File: focal_vs_bg_stats.json
{
  "schema_version": 1,
  "generated_at": "2026-MM-DD",
  "cohort_id": "Gar_226_hatchery",          // forward-compat (§5.5)
  "anchors": [
    {
      "anchor_id":      "cs_bp_LG28_001",   // matches cs_breakpoints_v1.id, OR
      "anchor_type":    "cs_breakpoint",    // "cs_breakpoint" | "candidate_boundary"
      "chrom":          "LG28",
      "focal_lo_bp":    14100000,
      "focal_hi_bp":    19270000,
      "focal_zone_def": "boundary_scan ±1.5Mb",   // free-form provenance string
      "background_scope": "per_chrom",      // "per_chrom" | "genome_wide_excl_focal"
      "metrics": {
        "fst":  { "n_focal": 78, "n_bg": 4202,
                  "focal_mean": 0.31,  "bg_mean": 0.04,
                  "focal_median": 0.28, "bg_median": 0.03,
                  "wilcoxon_W": 412934.0,
                  "wilcoxon_p": 5.84e-6,
                  "wilcoxon_logp_neg": 5.234,    // capped at 16
                  "z_focal_mean": 0.627,         // (focal_mean - bg_mean) / bg_sd
                  "z_focal_pct95": 1.92 },
        "theta_pi_ratio": { ... },
        "hobs_hexp": { ... }
      }
    }
  ]
}
```

When this file is loaded (drag-drop into the atlas, same loader pattern as
existing JSON layers): the widget displays its numbers as the canonical
Wilcoxon P, with a small lock icon indicating "manuscript-locked." Per-window
data is still pulled from `_galleryColumnGetData` for the visualization
(scatter + histogram bars) — Mode A overrides only the *summary statistics*,
not the visualization.

#### Mode B — live computation from existing per-window data

`_galleryColumnGetData(trackDef)` is called for the selected metric. The
widget partitions returned `{ mb, values }` into focal vs background using
the resolved focal interval (`_boundaryScanRange` or
`_focalRangeFromCsBreakpoint`). It computes:
- Histograms (Freedman–Diaconis bin width by default, fallback 30 bins)
- Means + medians
- Wilcoxon rank-sum two-sided (JS implementation, §5.2)
- Z-transform per window using background mean/sd

Live mode is fast at typical scale (~5000 windows per chrom × 4 metrics).
No async work required.

#### Mode C — empty-state

If neither precomputed nor live data is available for the selected metric on
the active chrom, the widget shows the same empty-state pattern as
`_renderRepeatDensityPanel` (atlas line ~16238): which chroms DO have data,
named JSONs to load.

### 5.2 Wilcoxon rank-sum in JS

Small standalone function in a new file `atlas_focal_vs_bg.js` (sibling to
`atlas_compare.js`, `atlas_ld.js`):

```javascript
// Returns { W, p_two_sided, logp_neg, z_stat }
// Normal approximation with continuity correction.
// Tied ranks averaged. NaN-skipping. Floor P at 1e-16.
function wilcoxon_rank_sum(focal_arr, bg_arr) { ... }
```

Tested against R's `wilcox.test(x, y, alternative='two.sided', exact=FALSE,
correct=TRUE)` on a small fixture (5–10 cases, mixed tie patterns and sample
sizes), committed to `tests/wilcoxon_fixture.json`. Tolerance: 1e-6 on
log10(P).

Exact P is not needed at the sample sizes this widget sees (n_focal ≥ 30 in
practice on a 1.5 Mb scan in 50 kb windows).

### 5.3 UI layout

```
┌──────────────────────────────────────────────────────────────────┐
│ [metric: F_ST ▾]  [bg: per-chrom ▾]  [Z ☑]  [view: full ▾]       │   ← row 1
│ scan radius:                                                     │   ← row 2
│ [1kb][5kb][10kb][25kb][100kb] | [1Mb][1.5Mb][2Mb][5Mb]    🔒 P=… │
├──────────────────────────────────────────────────────────────────┤
│  focal histogram                       │ background histogram     │
│  ▆▆▇▆▅▃▂▁                              │       ▁▃▆▇▆▃▁            │   ← panel C+D
│  mean = 0.627                          │ mean = -0.005            │
│  n = 78                                │ n = 4202                 │
├──────────────────────────────────────────────────────────────────┤
│ along-chromosome scatter, focal zone shaded                      │
│  · · ·  ▒▒▒▒▒  · · · · · · · · · · · · · · · · · · · · · · · · ·│   ← panel E
│         focal                                                     │
│ 0      |       breakpoint              |                  end Mb │
└──────────────────────────────────────────────────────────────────┘
```

**Row 1 (toolbar A) — picks**:
- `metric` dropdown: F_ST / θπ ratio / Hobs/Hexp / |Z|, gated by data
  availability via `_galleryColumnGetData(t).getData(state) != null`
- `bg` dropdown: `per-chrom` (default) | `genome-wide (excl focal)`
- `Z` checkbox: default ON, mirrors Spalax's Z-F_ST framing
- `view` dropdown: `full` (panel E Spalax-style) | `zoomed` (zoom to scan
  range, more useful at boundary scale on page 11)

**Row 2 (toolbar B) — focal radius**:
- Reuses the **exact** kb/Mb radius button layout from boundaries-page
  toolbar
- On page 11 widget mount: shares `_ensureBoundariesState().scan_radius_bp`
  with the parent page — clicking a radius button on the widget updates the
  whole boundaries page's scan radius, and vice versa. **One radius, one
  page.**
- On page 16 widget mount: maintains its own `state._focalVsBg.csRadiusBp`
  defaulting to 2 000 000 (Spalax). Independent from page 11's radius.
- P value display, right-aligned, with `(live)` (Mode B) or 🔒 (Mode A) icon

**Histograms (panels C+D)**:
- Same y-axis scale (`Math.max(maxFocalCount, maxBgCount)`); same x-axis
  range — visual comparison must be fair
- Mean line as dashed vertical, color `var(--accent)` for focal,
  `var(--ink-dim)` for bg. Mean labeled at top of line
- Bin width Freedman–Diaconis, capped at 60 bins
- Empty-state for either side: gray panel with "n=0 windows in focal" message

**Along-chromosome scatter (panel E)**:
- Reuses the SVG layout from `_renderRepeatDensityPanel` (margins, axes,
  toX/toY) — extracted into a small helper `_chromScatterSvg(opts)`
- X-axis: genomic position in Mb (full chromosome) or scan range (zoomed)
- Y-axis: metric value (or Z-metric if Z toggle on)
- Focal zone shaded `var(--accent)` α=0.12, dashed lines at focal_lo/hi
- Anchor marker:
  - page 16: dashed red line at breakpoint position (matches existing
    cs-breakpoint render style from turn 114a)
  - page 11: dashed marks at L_zone and R_zone positions
- Horizontal dashed line at y=0 (Spalax-style; meaningful for Z-transformed)
- Dot size 1.5 px, alpha 0.45, color `var(--accent)` inside focal,
  `var(--ink-dim)` outside
- Click on a window: dispatches focal-window event, populating per-sample
  lines panel — reuses turn-114a click-to-jump pattern

**No separate "control breakpoint" comparison panel** in v1 — Spalax's panel
F is shown by selecting a different breakpoint in the catalogue. Side-by-side
control comparator deferred until users ask.

### 5.4 Page integrations

#### Page 16 cross-species (`renderCrossSpeciesPage`)

Add a new section below `#csFlankCharts`, before `#csSynteny`:

```html
<div id="csFocalVsBg" style="display: none; margin-top: 14px;"></div>
```

In `_renderCrossSpeciesFocus`, after `flanks.innerHTML = _csBuildFlankCharts(bp);`:

```javascript
const fvb = document.getElementById('csFocalVsBg');
if (fvb) {
  fvb.style.display = 'block';
  const chromLen = _csIdeogramChromLength('gar', bp.gar_chr);
  const radius   = (state._focalVsBg && state._focalVsBg.csRadiusBp) || 2000000;
  const focal    = _focalRangeFromCsBreakpoint(bp, radius, chromLen);
  renderFocalVsBgWidget(fvb, {
    page_id:           'page16',
    anchor_id:         bp.id,
    anchor_type:       'cs_breakpoint',
    chrom:             bp.gar_chr,
    focal_lo_bp:       focal.start_bp,
    focal_hi_bp:       focal.end_bp,
    radius_source:     'self',                  // page 16 owns its radius
    radius_default_bp: 2000000,
    cohort_selector:   null,                    // day 1
  });
}
```

#### Page 11 boundaries (`renderBoundariesPage`)

Add a new section below `#bndSummary`:

```html
<div class="bnd-focal-vs-bg" id="bndFocalVsBg"></div>
```

In the function that re-renders when the user picks a candidate or changes
scan radius:

```javascript
const bs       = _ensureBoundariesState();
const cand     = _bndFindCandidate(bs.active_cand_id);
const chromLen = (state.data && state.data.n_bp) || null;
const scan     = _boundaryScanRange(cand, bs.scan_radius_bp, chromLen);   // EXISTING
const fvbEl    = document.getElementById('bndFocalVsBg');
if (fvbEl && cand && scan) {
  renderFocalVsBgWidget(fvbEl, {
    page_id:           'page11',
    anchor_id:         cand.id,
    anchor_type:       'candidate_boundary',
    chrom:             cand.chrom,
    focal_lo_bp:       scan.start_bp,
    focal_hi_bp:       scan.end_bp,
    radius_source:     'boundariesState',       // shared with parent page
    cohort_selector:   null,
  });
}
```

The widget's row-2 toolbar on page 11 binds clicks back into
`_ensureBoundariesState().scan_radius_bp` — so a radius change anywhere on
the page propagates to the existing scan-range visualizations (sim_mat,
repeat density, etc.), preserving consistency.

### 5.5 Forward-compatibility — multi-cohort group picker (day 2)

**This shapes the data model.** Quentin is adding species 2 (*C.
macrocephalus* wild cohort) and possibly more. The widget must accept
arbitrary cohort/group selections without rewriting.

#### Cohort registry (new top-level state key)

```javascript
state.cohort_registry = {
  cohorts: [
    { id: 'Gar_226_hatchery', label: 'C. gariepinus, 226 hatchery broodstock',
      species: 'C. gariepinus', n_samples: 226, sample_ids: [...],
      reference_genome: 'fClaHyb_Gar_LG.fa' },
    // Future:
    // { id: 'Mac_wild_NN', label: 'C. macrocephalus wild cohort', ... },
  ],
  groupings: [
    // Sub-cohort definitions (e.g. by karyotype class at a candidate)
    { id: 'LG28_inv_HOM_REF', cohort_id: 'Gar_226_hatchery', filter: {...} },
    { id: 'LG28_inv_HOM_INV', cohort_id: 'Gar_226_hatchery', filter: {...} },
  ],
};
```

#### Cohort-pair selector — UI

When the cohort registry has ≥2 cohorts (or ≥2 groupings), the widget
toolbar gains a `cohort_selector` chip:

```
[group A: HOM_REF ▾] vs [group B: HOM_INV ▾]
```

Both dropdowns list cohorts + groupings. The widget computes the metric
**between the two groups** (Spalax's two-cohort framing).

For day-1 single-cohort mode, the selector is hidden and the widget
implicitly compares karyotype classes (HOM_REF vs HOM_INV) — which is what
`fst_hom1_hom2` already means. **This is the right default; it's what the
cohort already gives.**

The `opts.cohort_selector` parameter is `null` on day 1 and becomes a real
selector spec on day 2:

```javascript
cohort_selector: {
  default_group_a: 'Gar_226_HOM_REF_at_anchor',
  default_group_b: 'Gar_226_HOM_INV_at_anchor',
  show_picker: true,
}
```

#### Migration path

1. Day 1: widget single-cohort, `cohort_selector: null`.
2. Day 2 (cohort registry + species 2): atlas adds the registry, widget
   grows the picker, no breaking change to existing pages — pages just
   start passing a non-null `cohort_selector`.
3. Aggregated JSON layer gains a `cohort_pair` field at schema_version 2.

**Key constraint**: `opts.cohort_selector` is in the widget API from day 1,
just always `null` for now. Don't bake "single-cohort" into the widget.

### 5.6 Per-page sticky preferences

```javascript
state._focalVsBgPrefs = {
  page16: { metric: 'fst', bg_scope: 'per_chrom', z: true, viewMode: 'full',
            csRadiusBp: 2000000 },
  page11: { metric: 'fst', bg_scope: 'per_chrom', z: true, viewMode: 'zoomed' },
  // page 11 does NOT store its own radius — uses parent state's
  // _ensureBoundariesState().scan_radius_bp
};
```

Survives candidate / breakpoint changes within the same page. localStorage
backed under `focal_vs_bg::<page_id>` per the existing
`_repeatDensityPrefs` pattern.

---

## 6. Module structure (suggested layout)

```
atlas/
  atlas_focal_vs_bg.js         # NEW — widget renderer + Wilcoxon JS +
                               # _focalRangeFromCsBreakpoint helper +
                               # _radiusToolbarHTML extraction (or moved here)
  Inversion_atlas.html         # MODIFIED — page 16 + page 11 widget mounts
  tests/
    test_focal_vs_bg.html      # NEW — widget rendering + interaction tests
    wilcoxon_fixture.json      # NEW — Wilcoxon JS vs R reference values

specs_todo/
  SPEC_focal_vs_background_widget.md   # this file

phase_X_focal_vs_bg/                   # OPTIONAL — server-side aggregation
  README.md
  STEP_FVB_01_extract_per_window_per_class.py
  STEP_FVB_02_aggregate_focal_bg_stats.R
  STEP_FVB_03_export_focal_vs_bg_stats_json.py
```

The server-side phase X is **optional**. The widget runs Mode B (live)
without it. Phase X exists only when the manuscript figure freezes.

---

## 7. Tests

### JS unit tests (`test_focal_vs_bg.html`)

- `wilcoxon_rank_sum_against_R_fixture`: each fixture row has `(focal, bg,
  expected_W, expected_p)` from R's `wilcox.test(x, y, alternative='two.sided',
  exact=FALSE, correct=TRUE)`. Tolerance 1e-6 on log10(P).
- `histogram_binning_fd_falls_back_at_low_n`: n=5 in focal still renders.
- `z_transform_uses_bg_only`: explicit check that mean/sd come from bg only.
- `bg_scope_genome_wide_excludes_focal`: scope = `genome_wide_excl_focal`
  removes focal windows before bg mean/sd.
- `widget_handles_null_cohort_selector_day_1`: opts with `cohort_selector:
  null` renders without errors; toolbar does not show cohort picker.
- `widget_uses_aggregated_json_when_present`: mocked aggregated JSON in
  state, widget shows the locked P value and 🔒 icon.
- `widget_falls_back_to_live_when_aggregated_missing`: widget shows
  `(live)` tag.
- `widget_reuses_boundary_scan_range`: page 11 mount with a candidate, the
  focal interval matches `_boundaryScanRange(cand, bs.scan_radius_bp, chromLen)`
  exactly — no off-by-one with the existing function.
- `widget_radius_change_propagates_to_boundaries_state`: clicking a radius
  button in the page 11 widget toolbar updates
  `_ensureBoundariesState().scan_radius_bp`, AND the existing repeat-density
  panel re-renders with the new radius (verifying no parallel-state bug).

### Page-integration tests

- `page16_widget_mounts_on_breakpoint_focus`: select a bp, `#csFocalVsBg`
  is non-empty.
- `page11_widget_mounts_on_candidate_select`: pick a candidate,
  `#bndFocalVsBg` is non-empty.
- `page16_focal_radius_change_recomputes_wilcoxon`: change radius from 2 Mb
  to 5 Mb, P value updates.
- `page11_widget_uses_existing_scan_radius`: `bs.scan_radius_bp = 1500000`
  → widget focal interval width = 3 Mb + candidate span (auto-expanded for
  huge candidates per `_boundaryScanRange`).

Target: 12–15 tests, all green before integration ships.

---

## 8. Open questions / explicitly deferred

- **Multiple-testing correction across candidates.** The widget reports a
  per-anchor Wilcoxon P. Across 200+ candidates, BH correction is needed
  for any genome-wide claim. Deferred to manuscript-time analysis. The
  catalogue page (page 5) could gain an FDR-q column derived from the
  Mode-A aggregated JSON. Spec deferred.
- **Paired-anchor (case/control) comparison panel.** Spalax's panels E vs
  F show a real breakpoint vs a control. The atlas does this by selecting
  between catalogue rows. Side-by-side comparator deferred.
- **Effect size beyond mean shift.** Cliff's δ or Mann-Whitney r — easy
  but cluttered. Deferred.
- **Asymmetric tail tests.** Hobs has a "low" tail of interest (LoH inside
  inversions). Widget should support `tail: 'high' | 'low' | 'two-sided'`,
  default `two-sided`, but page 11 mount might pass `'low'` for Hobs.
  5-line change; flagged so it doesn't get forgotten.
- **dXY** is not in the day-1 metric pool because there's no dxy track in
  `POPSTATS_TRACKS` today. When the dxy popstats track lands (sibling work
  on page 8), the widget auto-discovers it. Worth adding as a popstats
  track — distinguishing F_ST from dXY in the widget's view is the
  difference between "young sweeping" and "old balanced" inversions.
- **Auto-discovery from `state.data.tracks`.** The widget's metric pool
  should dynamically include any per-window scalar track auto-discovered
  via `state.data.tracks` (the GHSL+precomp tracks dict), not just the
  hard-coded `POPSTATS_TRACKS` list. Same auto-discovery pattern as
  `buildPopstatsTrackList()` (atlas line ~52146). This was implicit in
  v1; making explicit here.

---

## 9. Summary

One widget, two page mounts on day 1. Reuses `_boundaryScanRange`,
`POPSTATS_TRACKS`, `_galleryColumnGetData`, the page 11 scan-radius toolbar,
the `_renderRepeatDensityPanel` empty-state and viewMode patterns, the turn
114a click-to-jump, and the `_repeatDensityPrefs` localStorage pattern. New
code: Wilcoxon JS, `_focalRangeFromCsBreakpoint` (small), the histogram +
scatter renderers, and the cohort-selector forward-compat scaffold.

Three data modes: Mode A precomputed JSON (manuscript-locked, 🔒), Mode B
live JS (`(live)` tag), Mode C empty-state. Forward-compatible with
multi-cohort selection on day 2 (species 2). Mirrors Spalax Figure 4 panels
C–F. Panels A/B (rfmix) deferred.

End of spec.
