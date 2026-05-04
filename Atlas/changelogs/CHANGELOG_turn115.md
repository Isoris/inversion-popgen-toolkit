# Turn 115 — Cross-species dot plot integration

**Spec implemented**: `SPEC_cross_species_dotplot.md`

## What shipped

The cross-species dot plot widget is now wired into page 16. It renders
below `#csSynteny`, sitting next to (not replacing) the existing
macro-synteny permutation panel.

### Two data sources, both optional

1. **wfmash synteny_blocks** — already inside `cs_breakpoints_v1.json`
   schema_version 2. Read directly from `state.crossSpecies`. No new
   file needed; the panel renders as soon as a v2 cs_breakpoints JSON
   is loaded. Default resolution = `wfmash`.

2. **mashmap multi-resolution** — a new optional layer
   `dotplot_mashmap_v1.json`. Generated once on the cluster
   (`STEP_DP01_run_mashmap_grid.sh` + `STEP_DP02_aggregate_mashmap_to_json.py`,
   per the spec). Drag-drop into the atlas; auto-detected; persists in
   `localStorage` under `inversion_atlas.dotplot_mashmap`. When loaded,
   adds a resolution dropdown to the panel header.

If neither source is loaded, the dotplot wrapper stays hidden — same
fail-soft pattern as `#csSynteny`.

### UX

- Mini panel inline (~300 px) in the page-16 focus column
- Hover-in: instant enlarge to a fixed-position popup (~720 px) with
  backdrop dim and color legend
- Click on enlarged popup: pins it (header gains a 📌 marker, backdrop
  darkens). Outside-click or × button closes; Esc also closes.
- Tooltip on hover in enlarged mode: chrom × chrom, bp ranges, strand,
  PI, block size

### Visual conventions

- 5-stop perceptual color ramp from slate-blue (PI 85%) to brick-red
  (PI 100%). Lives in the atlas accent family.
- Forward strand: solid α 0.85. Reverse strand: solid α 0.55 + 1.5 px
  dotted overlay. Strand encoded by texture, not color.
- LIS-flavour fattest-alignment diagonal layout — chroms reordered so
  heaviest synteny falls on the main diagonal. Reverse-flipped target
  chroms marked `*` in axis labels.

## Files modified

- `Inversion_atlas.html` (+192 lines)
  - `<script src="atlas_dotplot.js">` import (line ~56495)
  - `<div id="csDotplot">` mount under `#csSynteny` (page-16 markup)
  - `_renderCrossSpeciesDotplot()` renderer (next to `_renderCrossSpeciesSynteny`)
  - `dotplot_mashmap_v1.json` loader: `_isDotplotMashmapJSON`,
    `_storeDotplotMashmap`, `_persistDotplotMashmap`,
    `_restoreDotplotMashmap`, `_clearDotplotMashmap`
  - `_classifyJSONKind` updated to recognize `dotplot_mashmap`
  - Dispatch added to `loadMultipleJSONs`
  - Boot-time restore call after `_restoreCrossSpecies`
  - Help-page table row for `dotplot_mashmap_v1.json`
  - Page-16 render flow now calls `_renderCrossSpeciesDotplot()` after
    `_renderCrossSpeciesSynteny()`, fail-soft

- `atlas_dotplot.js` (NEW, 844 lines) — the widget renderer module.
  UMD; attaches to `window.popgenDotplot`. Exports
  `makeDotplotPanel`, `colormap`, `layoutChroms`.

## Cache behavior

Panel is cached on `state._csDotplotPanel` with a fingerprint key
`state._csDotplotPanelFp = "wf:<n_blocks>:<loaded_at>|mm:<n_resolutions>:<generated_at>"`.
Re-renders of page 16 (catalogue stepping, filter changes, etc.) hit
the cache and don't rebuild the panel — so the user's pinned-popup
state survives unrelated re-renders. Cache invalidates only when the
underlying data sources change.

## Tests run

End-to-end Node integration test with DOM shim + a synthetic
cs_breakpoints v2 fixture (4 synteny blocks including one inversion).
Verified:
- Module loads, public surface is `{makeDotplotPanel, colormap, layoutChroms}`
- Renderer runs without errors
- Slot becomes visible after first render
- Cache hit on identical re-render (panel object identity preserved)
- Cache invalidates on data change (mashmap layer added → fresh panel)

`node --check` on the entire concatenated inline JS of the patched
atlas: no syntax errors.

## Spec status

`SPEC_cross_species_dotplot.md`: **IMPLEMENTED** in turn 115.
Move to `specs_shipped/` per the README instructions.

### Open items deferred per spec §8

- Per-breakpoint highlight (yellow rectangle / red marker) — deferred
- Cross-link from dot plot to breakpoint catalogue — deferred
- Pan/zoom inside enlarged popup — deferred
- Cluster-side `STEP_DP01_run_mashmap_grid.sh` /
  `STEP_DP02_aggregate_mashmap_to_json.py` not in scope this turn —
  pipeline scripts are documented in the spec but not yet written.
  Day 1 path (wfmash synteny_blocks alone) works without them.

## Next turns

Suggest next: **`SPEC_focal_vs_background_widget.md`** — the keystone
that unlocks ncRNA + LD decay metrics for free. Or
**`SPEC_ncrna_density_layer.md`** if a smaller turn is preferred.
