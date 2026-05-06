# SPEC — LD decay overlay (Spalax Fig. S62 idiom)

**Status**: forward-looking spec. Not implemented. Created same session as
`SPEC_focal_vs_background_widget.md` and `SPEC_xpehh_track.md`.

**Reading order**: this spec → `SPEC_focal_vs_background_widget.md`
(architectural sibling — same Mode A / Mode B / Mode C data-flow pattern,
same dual-page mount on page 16 + page 11) → atlas `atlas_ld.js` (the existing
LD module the renderer extends) → `server/server_turn11c_ld_fast/fast_ld_endpoint.py`
(the existing endpoint this spec consumes — no server changes required for
day 1) → `server/engine_fast_ld/` (reference: where decay-deciles + r² matrix
already get computed in C).

**Reference figure**: Spalax mole-rat fusion paper, Fig. S62 — LD decay
patterns for three fused chromosomes vs whole genome. Each panel: distance
on x (kb), r² on y, three traces (focal fusion zone / whole-fused-chromosome /
whole-genome reference). Focal LD elevation visible in panels A–C; absent in
panel D (control case).

**One-line summary**: fetch LD on three scopes (focal zone, focal chromosome,
whole genome) via the existing fast_ld endpoint; render the three (distance,
r²) traces on one axes with smoothed decay curves. Reuses the entire fast_ld
client + server infrastructure — no engine changes.

---

## 0. Settled infrastructure — REUSE these, don't reinvent

The atlas already has a complete LD pipeline. This spec composes existing
pieces; it does NOT add a new LD engine, server, or data path.

### 0.1 LD engine — **reuse `fast_ld.c`**

`fast_ld.c` (`server/engine_fast_ld/fast_ld.c`) already computes per-pair r²
across requested SNPs in C. Output per group:

- `pairs_b64`: uint8 r²·255, upper triangle, row-major
- `decay_deciles[10]`: per-distance-decile mean r²
- `median_r2_overall`, `median_r2_shelf`, `median_r2_flank`, `shelf_ratio`
- `pct_pairs_above_0_8`

The decay scatter renderer this spec adds **does not need a new server
endpoint**. It calls the existing `/fast_ld` endpoint three times (focal /
chrom / genome) and reconstructs (distance, r²) per pair from `pairs_b64` +
`sites.pos[]` client-side.

### 0.2 LD endpoint — **reuse `/fast_ld`**

`server/server_turn11c_ld_fast/fast_ld_endpoint.py` already accepts:

```python
class FastLDReq(BaseModel):
    chrom: str
    window_range: List[int]   # [lo, hi] inclusive window indices
    snp_cap: int = ...
    thin_to: int = ...
    shelf_bp: ...
    groups: ...
```

The decay overlay calls this endpoint with three different
`(chrom, window_range)` configurations per anchor — no schema change.

### 0.3 Atlas LD module — **extend `atlas_ld.js` (don't replace)**

`atlas/atlas_ld.js` already exposes (on `window.popgenLD`):
- `APEX_COLORS`, `apexColormap()`, `colorAt()`
- `decodePairsB64(b64, n)` → Uint8Array of pairs
- `pairAt(arr, n, i, j)` → uint8 r²·255
- `drawSplitHeatmap(canvas, payload, opts)`
- `makeLDSplitPanel(opts)` — current page-3 panel
- `ldSplitHeatmap(req, opts)` — request wrapper

The decay overlay adds two new exports to the same module (sibling to
`makeLDSplitPanel`):

- `makeLDDecayPanel(opts)` — the new SVG panel with the 3-trace overlay
- `ldDecayRequest(opts)` — dispatch helper that fires three `ldSplitHeatmap`
  requests (or pulls from the genome-wide cache)

**Do not create a new `atlas_ld_decay.js` file** — the LD code lives in one
module.

### 0.4 Page-3 LD heatmap stays as-is

The existing split-triangle LD heatmap on page 3 (catalogue) is a
**different visual idiom** (apex-colormap heatmap of all pairs in a focal
window). The decay overlay is a **complementary** idiom (distance-vs-r²
scatter with curves). They show different things; both are useful. No
deprecation, no merge.

### 0.5 fast_ld payload already has `sites.pos[]`

`fast_ld_endpoint.py` returns `sites: { idx[], pos[], maf_<g>[], ... }`. With
`pos[]` and `pairs_b64`, the client can reconstruct (distance, r²) for every
pair — `_recoverPairwiseDistanceR2(payload)` is a small client-side helper.

---

## 1. Purpose

Render Spalax-style LD decay overlay per anchor (cross-species breakpoint
or candidate inversion boundary):

```
  r²
  1.0 ┤ ░░
       ░░░░  focal fusion zone (light blue, raw scatter)
  0.8 ┤ ░░░░░░
       ░░░░░░░░  ← Hill curve fits overlaid
  0.6 ┤  ╲
       │   ╲___________________________
  0.4 ┤    ╲    ─────────── chrom curve (dark line)
       │     ╲___________ ──────────── genome curve (orange line)
  0.2 ┤        ╲────────────────────────────
       │
  0.0 ┼─────┬─────┬─────┬─────┬─────┬──── distance (kb)
       0   50   100  150  200  250  300
```

Three traces:
1. **Focal scatter** — every SNP-pair in the focal zone, plotted as raw
   `(distance_bp, r²)` dots (light blue, α≈0.25). Same density as Spalax.
2. **Focal-chromosome curve** — Hill-fit smoothed curve from the entire
   chromosome's pairs (dark blue solid line).
3. **Whole-genome curve** — Hill-fit smoothed curve from a precomputed
   genome-wide LD decay (orange solid line).

The visual logic mirrors Spalax exactly: dots above the curves =
fusion-zone LD elevation; dots tracking or below the curves = no
fusion-specific LD signal (panel-D-style control case).

---

## 2. Why this matters

LD elevation near structural breakpoints is a direct readout of
**recombination suppression**. Inversions suppress recombination across
their span, which inflates r² between markers inside (or near) the inversion
relative to genome-wide background. Spalax's fusion analogue is the same
phenomenon — fusion creates linkage between previously-unlinked chromosomes,
which inflates LD across the fusion point.

For Quentin's cohort, an LD-decay overlay at each candidate inversion
boundary shows:
- **Inverted region elevated** = expected canonical inversion signal
- **Inverted region tracks chromosome curve** = inversion not segregating
  in this cohort (no observed LD elevation despite the candidate call)
- **Inverted region BELOW chromosome curve** = unusual; could indicate the
  candidate is a population-specific recombination hotspot rather than a
  true inversion

The third case is the diagnostic Spalax made explicit with their panel D
control. It's the case the atlas currently can't show — page 3's
split-heatmap can show "this region has high r²" but cannot put that on
the same axes as a genome-wide expectation.

---

## 3. What this DOES / does NOT claim

### DOES

- A visual comparison of LD-vs-distance in the focal zone vs the
  chromosome vs the genome, on one axes
- A Hill-fit decay curve per scope, with fit parameters in the panel
  caption
- Detection of LD elevation or suppression at a glance

### Does NOT

- Causation — recombination suppression is consistent with inversions but
  also with selection, drift, or population structure
- A formal test against the genome curve (no Wilcoxon, no permutation in
  v1 — the focal-vs-bg widget does that for per-window F_ST, this widget
  is descriptive)
- That the focal zone is correctly placed (upstream concern)

---

## 4. Vocabulary discipline

| Allowed | Discouraged |
|---|---|
| "LD is elevated within the focal zone relative to chromosome and genome curves" | "LD proves the inversion is suppressing recombination" |
| "Decay shelf extends to ~N kb in the focal zone vs ~M kb genome-wide" | "Recombination is suppressed" (without distance specification) |
| "r² remains above 0.4 to 200 kb in the focal zone" | "very high LD" without numbers |
| "Consistent with recombination suppression across the inversion span" | "Confirms the inversion" |

---

## 5. Architecture

### 5.1 Data flow — three modes (mirrors focal-vs-bg spec §5.1)

#### Mode A — precomputed genome-wide decay JSON (lock the orange curve)

The orange "whole genome" curve is the same for every anchor on a given
cohort. Computing it on every request would be wasteful. **A precomputed
genome-wide decay JSON locks this curve once per cohort**:

```jsonc
// File: ld_decay_genome_wide.json
{
  "schema_version": 1,
  "cohort_id": "Gar_226_hatchery",
  "generated_at": "2026-MM-DD",
  "snp_cap_per_chrom": 5000,         // sub-sampling — reproducibility
  "thin_to": null,
  "min_maf": 0.05,
  "n_chroms_used": 28,
  "n_pairs_total": 3850000,
  "distance_bins_bp": [
    [0, 1000], [1000, 2000], ..., [299000, 300000]
  ],                                  // 100 bins, 3 kb each up to 300 kb
  "r2_mean_per_bin":   [0.51, 0.42, ..., 0.13],
  "r2_median_per_bin": [0.49, 0.40, ..., 0.12],
  "n_pairs_per_bin":   [128043, 117821, ..., 9123],
  "hill_fit": {                       // Hill et al. 1988 model, 4-param fit
    "C": 0.0089, "n_eff": 226, "rmse": 0.018
  }
}
```

When this file is loaded, it provides the orange whole-genome curve on
every panel without re-running fast_ld across all chromosomes. **Drag-drop
loader, same pattern as other JSON layers.**

When this file is absent, the panel renders without the orange trace and
shows a small `(genome curve unavailable — load ld_decay_genome_wide.json)`
hint inline. The blue and dark-blue traces still render.

#### Mode B — live, per-anchor

The dark-blue chromosome curve and light-blue focal scatter are computed
live, per-anchor:

1. Resolve focal interval from `_boundaryScanRange(cand,
   bs.scan_radius_bp, chromLen)` (page 11) or
   `_focalRangeFromCsBreakpoint(bp, radius_bp, chromLen)` (page 16).
2. Translate focal interval → `window_range` indices.
3. Fire `ldSplitHeatmap({chrom, window_range: [focal_lo_idx, focal_hi_idx],
   ...})` — the focal request.
4. Fire `ldSplitHeatmap({chrom, window_range: [0, max_idx], ...})` — the
   chromosome request.
5. For each, decode `pairs_b64`, build `(distance_bp, r²)` pairs from
   `sites.pos[]`, render scatter + binned mean curve.

LRU cache the chromosome-scope payload (it's the same for every anchor on
that chrom). `fast_ld_endpoint.py` already does cohort-keyed caching
server-side; the atlas adds anchor-→-payload caching client-side via
`state._ldDecayCache`.

#### Mode C — empty state

If the fast_ld server is unreachable or returns errors: same empty-state
pattern as `_renderRepeatDensityPanel` (atlas line ~16238) — named the
endpoint URL it tried, suggest `_loadFile` as the fallback for
diagnostics.

### 5.2 Hill curve fit (client-side JS)

The Hill et al. (1988) decay model:

```
E[r² | distance d, sample size n_eff] =
    [10 + C·d] / [(2 + C·d) · (11 + C·d)] ·
    [1 + (3 + C·d)(12 + 12·C·d + (C·d)²) / (n_eff · (2 + C·d) · (11 + C·d))]
```

- `C = 4·N_e·c` (population-scaled recombination rate per bp)
- `n_eff` = effective sample size

Fit `C` (and optionally `n_eff` if not supplied) via Levenberg-Marquardt or
simple bounded golden-section search on the residual sum of squares. JS
implementation in `atlas_ld.js` as `_hillFit(distances_bp, r2_values, n_eff)`,
returns `{ C, n_eff, rmse, predict: fn }`.

Test against R's `nls()` Hill fit on a small fixture (3–5 cases of varying
size and decay rate), committed to `tests/hill_fixture.json`. Tolerance
1e-3 on `C`.

The Hill fit is **only** used for the smooth dark-blue (chromosome) and
orange (genome) curves. The focal scatter is shown as raw dots, not fitted
— Spalax does the same; the focal sample is small enough that a fit would
overcommit.

### 5.3 UI — `makeLDDecayPanel(opts)`

```javascript
window.popgenLD.makeLDDecayPanel({
  page_id:           'page16',
  anchor_id:         bp.id,
  anchor_type:       'cs_breakpoint',
  chrom:             bp.gar_chr,
  focal_lo_bp:       focal.start_bp,
  focal_hi_bp:       focal.end_bp,
  max_distance_kb:   300,                 // x-axis cap
  show_genome_curve: true,                // false → omit orange trace
  cohort_selector:   null,                // forward-compat (§5.5)
});
```

Returns a DOM panel:

```
┌──────────────────────────────────────────────────────────────────┐
│ LD decay · LG28 · cs_bp_LG28_001 · focal 14.10–19.27 Mb          │   ← header
│ [scope: 0–300 kb ▾]      legend: ─ focal · ─ chrom · ─ genome     │
├──────────────────────────────────────────────────────────────────┤
│  r²                                                              │
│  1.0 ┤ ░░ ░ ░  ░                                                 │
│       ░░░░░░░░░░░░ ░ ░ ░  ░    ░    ░     ░       ░              │
│  0.5 ┤  ╲                                                        │
│       │   ╲_______                                               │
│       │           ──────────── chrom (Hill C = 0.012)            │
│  0.0 ┤              ──────────── genome (Hill C = 0.009, locked) │
│       └─────┴─────┴─────┴─────┴─────┴─────┴───── distance (kb)   │
│       0    50   100  150  200  250  300                          │
│                                                                  │
│ caption: focal-zone shelf-ratio = 1.85× chrom; focal pairs above │
│         r² = 0.4 = 73 / 412 (17.7%)                              │
└──────────────────────────────────────────────────────────────────┘
```

**Toolbar (header row 2)**:
- `scope` dropdown: 0–100 / 0–200 / 0–300 / 0–500 / 0–1000 kb x-axis cap.
  Default 300 kb (Spalax). Changing scope re-renders only — no new fetch.

**SVG panel**:
- X-axis: distance in kb (linear, 0 to scope cap)
- Y-axis: r², linear, [0, 1]
- Focal dots: light blue (`#9ec5e8`), radius 1.5 px, alpha 0.30
- Chrom curve: dark blue (`#1f4e79`), 1.5 px solid line
- Genome curve: orange (`#e8734a`), 1.5 px solid line
- Light dashed grid, font 10 px monospace axis labels
- Caption below: shelf-ratio (focal/chrom), pct-pairs-above-r²-0.4 in
  focal, Hill C for chrom and genome (with 🔒 icon next to genome C if
  Mode A locked)

**No legend overlapping the data area** — legend in header row, not in
plot. Spalax's inset legend works because their panels are huge; for an
embedded atlas widget, header legend reads better.

### 5.4 Page integrations

Same dual-mount pattern as the focal-vs-bg widget. The two widgets sit
**stacked** — focal-vs-bg first (histogram + scatter), then LD decay
beneath it — on both page 16 and page 11.

#### Page 16 cross-species

Add below `#csFocalVsBg`:

```html
<div id="csLDDecay" style="display: none; margin-top: 14px;"></div>
```

In `_renderCrossSpeciesFocus`:

```javascript
const ld = document.getElementById('csLDDecay');
if (ld) {
  ld.style.display = 'block';
  const chromLen = _csIdeogramChromLength('gar', bp.gar_chr);
  const radius   = (state._focalVsBg && state._focalVsBg.csRadiusBp) || 2000000;
  const focal    = _focalRangeFromCsBreakpoint(bp, radius, chromLen);
  ld.innerHTML = '';
  ld.appendChild(window.popgenLD.makeLDDecayPanel({
    page_id:           'page16',
    anchor_id:         bp.id,
    anchor_type:       'cs_breakpoint',
    chrom:             bp.gar_chr,
    focal_lo_bp:       focal.start_bp,
    focal_hi_bp:       focal.end_bp,
    show_genome_curve: !!state.ldDecayGenomeWide,
    cohort_selector:   null,
  }));
}
```

#### Page 11 boundaries

Add below `#bndFocalVsBg`:

```html
<div class="bnd-ld-decay" id="bndLDDecay"></div>
```

In the candidate-change re-render:

```javascript
const bs       = _ensureBoundariesState();
const cand     = _bndFindCandidate(bs.active_cand_id);
const chromLen = (state.data && state.data.n_bp) || null;
const scan     = _boundaryScanRange(cand, bs.scan_radius_bp, chromLen);
const ldEl     = document.getElementById('bndLDDecay');
if (ldEl && cand && scan) {
  ldEl.innerHTML = '';
  ldEl.appendChild(window.popgenLD.makeLDDecayPanel({
    page_id:           'page11',
    anchor_id:         cand.id,
    anchor_type:       'candidate_boundary',
    chrom:             cand.chrom,
    focal_lo_bp:       scan.start_bp,
    focal_hi_bp:       scan.end_bp,
    show_genome_curve: !!state.ldDecayGenomeWide,
    cohort_selector:   null,
  }));
}
```

The page 11 widget reuses the boundaries-page scan radius implicitly
through `_boundaryScanRange(cand, bs.scan_radius_bp, ...)` — when the user
changes the radius button, the LD decay panel re-renders along with
everything else.

### 5.5 Forward-compatibility — multi-cohort

**Identical pattern to focal-vs-bg spec §5.5.** The widget API accepts
`cohort_selector` from day 1, always `null` until species 2 lands. When a
cohort registry exists (cohort + grouping objects per
`SPEC_focal_vs_background_widget.md`), the LD-decay panel grows the same
"[group A ▾] vs [group B ▾]" picker — and computes LD **separately per
group** (two scatter clouds + two chrom curves), still on one axes.

The Mode A genome-wide JSON gains a `cohort_pair` field at schema_version
2 to handle multi-cohort genome curves.

### 5.6 Per-page sticky preferences

```javascript
state._ldDecayPrefs = {
  page16: { scope_kb: 300, show_genome_curve: true },
  page11: { scope_kb: 300, show_genome_curve: true },
};
```

localStorage backed under `ld_decay::<page_id>`.

### 5.7 Performance budget

For a typical anchor:
- Focal request: ~80 windows × ~5000 SNPs cap = 5000 SNPs → ~12.5 M pairs.
  fast_ld in C is ~0.5 s/M pairs → 6 s wallclock.
- Chrom request: full LG28 (~4300 windows) capped at 5000 SNPs total
  (existing endpoint behavior) → ~12.5 M pairs → 6 s. Cached after first
  call per chrom.
- Genome curve: precomputed (Mode A) — zero cost at render time.

First load per anchor: ~12 s (focal + chrom in parallel). Cached navigations
between anchors on the same chromosome: ~6 s (chrom payload reused).

If this is too slow in practice, the focal request can be sub-sampled
further (`snp_cap=2500`) and the chrom request can be sub-sampled to a
random window selection. Both are server-side options that already exist
on `FastLDReq`.

---

## 6. Module structure (suggested layout)

```
atlas/
  atlas_ld.js                  # MODIFIED — adds makeLDDecayPanel,
                               # ldDecayRequest, _hillFit,
                               # _recoverPairwiseDistanceR2
  Inversion_atlas.html         # MODIFIED — page 16 + page 11 widget mounts
  tests/
    test_ld_decay.html         # NEW — panel rendering + Hill fit tests
    hill_fixture.json          # NEW — Hill JS vs R reference values

specs_todo/
  SPEC_ld_decay_overlay.md     # this file

phase_X_genome_wide_ld/        # OPTIONAL but recommended for Mode A
  README.md
  STEP_LD_GW_01_sample_genome_pairs.py
  STEP_LD_GW_02_aggregate_decay_bins.R
  STEP_LD_GW_03_export_ld_decay_genome_wide_json.py
```

---

## 7. Tests

### JS unit tests (`test_ld_decay.html`)

- `recover_pairwise_distance_r2_correct`: given a payload with N=10 SNPs at
  known positions and a synthetic `pairs_b64`, the recovered
  `[distance, r²]` array has length N·(N-1)/2 with correct distance
  computation per pair.
- `hill_fit_against_R_fixture`: each fixture row has
  `(distances, r2, n_eff, expected_C, expected_rmse)`. Tolerance 1e-3 on C.
- `panel_renders_focal_only_when_chrom_request_fails`: mock chrom failure,
  panel still shows focal scatter and an inline "(chrom curve unavailable)"
  hint.
- `panel_omits_genome_trace_when_layer_absent`: no `state.ldDecayGenomeWide`
  → orange trace not rendered, header note shown.
- `panel_uses_genome_curve_when_layer_present`: mocked
  `ld_decay_genome_wide.json` → orange trace rendered with 🔒 icon next to
  Hill C.
- `scope_change_does_not_refetch`: change scope dropdown 300 → 500 kb,
  fetch counter unchanged.
- `chrom_payload_reused_between_anchors`: two anchors on same chrom, second
  render does not re-fetch chrom payload (one cache hit).
- `widget_handles_null_cohort_selector_day_1`: opts with `cohort_selector:
  null` renders without errors; toolbar does not show cohort picker.

### Page-integration tests

- `page16_ld_panel_mounts_on_breakpoint_focus`: select bp, `#csLDDecay`
  is non-empty.
- `page11_ld_panel_mounts_on_candidate_select`: pick candidate,
  `#bndLDDecay` is non-empty.
- `page11_ld_panel_uses_existing_scan_radius`: focal interval matches
  `_boundaryScanRange(cand, bs.scan_radius_bp, chromLen)` exactly.
- `radius_change_re_renders_ld_panel_on_page11`: radius 1.5 Mb → 5 Mb, LD
  panel re-fetches focal payload (chrom payload still cached).

Target: 12–14 tests, all green before integration ships.

---

## 8. Open questions / explicitly deferred

- **Statistical test of focal vs chrom decay** — the v1 panel is
  descriptive (visual + Hill C). A formal test (KS on r²-vs-distance
  distribution, or comparison of fitted C parameters with bootstrap CI)
  is deferred. The focal-vs-bg widget already provides Wilcoxon on
  per-window F_ST, which is the more interpretable per-anchor test.
- **Phased haplotypes** — fast_ld currently uses dosage. r² from dosage is
  a slight underestimate of phased r² but the bias is consistent across
  scopes and doesn't affect the comparison. When phased VCF lands (per
  `SPEC_xpehh_track.md` prerequisite), fast_ld can optionally consume
  phased data — flagged in `engine_fast_ld/README.md`, not specced here.
- **Per-karyotype-class decay curves** — split focal pairs by HOM_REF /
  HET / HOM_INV (using the existing `groups` parameter on `FastLDReq`).
  Would surface arrangement-specific recombination suppression. Deferred
  to a v2 of this spec because the visual gets cluttered with 5+ traces.
- **Genome-wide curve refresh trigger** — when the cohort changes (sample
  add/remove), `ld_decay_genome_wide.json` must be regenerated. No
  automated invalidation in v1; user reloads the file. Mode B shows the
  blue traces while waiting.
- **Whole-chromosome decay variation across LGs** — the dark-blue chrom
  curve is per-anchor's chromosome. For a clean manuscript figure, all
  panels should use the same y-scale; auto-scaling is the v1 default but
  should be locked to [0, 1] for export-mode (deferred polish).

---

## 9. Summary

One LD decay overlay panel, two page mounts on day 1 (page 16 + page 11),
sitting beneath the focal-vs-bg widget. Reuses the entire fast_ld stack
(C engine, FastAPI endpoint, atlas_ld.js client decoder) without
modification. New code: `makeLDDecayPanel` + `_hillFit` +
`_recoverPairwiseDistanceR2` in `atlas_ld.js`, an optional Mode-A
genome-wide-decay JSON layer, and the two page mounts.

Three traces per panel: focal scatter (light blue), chromosome curve
(dark blue), genome curve (orange, locked when JSON loaded). Mirrors
Spalax Fig. S62. Forward-compatible with multi-cohort selection on day 2
via the same `cohort_selector` parameter as the focal-vs-bg widget.

End of spec.
