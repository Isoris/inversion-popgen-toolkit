# Turn 116 — ncRNA density layer (rRNA / tRNA / ncRNA-other)

**Spec implemented**: `SPEC_ncrna_density_layer.md` (data-layer + MVP renderer)

## What shipped

A new sibling track to the existing repeat-density panel on the boundaries
page (page 11). Renders per-chromosome density (loci/Mb) of tRNA, rRNA,
and ncRNA-other features from the canonical GFF3 annotations Quentin
already has on disk.

### Three-family display with sub-class chips

- **tRNA**: `tRNA_all`, `tRNA_HC`, `tRNA_pseudo`, `tRNA_intronic`, `tRNA_Sec`
- **rRNA**: `rRNA_all`, `rRNA_5S`, `rRNA_5_8S`, `rRNA_18S`, `rRNA_28S`, `rRNA_partial`
- **ncRNA**: `ncRNA_all`, `ncRNA_snRNA`, `ncRNA_snoRNA`, `ncRNA_miRNA`, `ncRNA_other`

Click a chip to switch the panel's active class. Selection persists per chrom
in `localStorage` under `inversion_atlas.ncRNADensityActiveClass`.

### Three biology cases the panel makes visible

1. **tRNA cluster density** elevated near a candidate boundary — the Larkin
   2009 / Bermudez-Santana 2010 EBR-style breakpoint-hotspot signal
2. **rDNA array** (5S / 5.8S / 18S / 28S spikes) overlapping a candidate —
   independent fragility mechanism via tandem-repeat recombination
3. **3× 5S rRNA expansion in Mac vs Gar** (12,209 vs 4,272 loci genome-wide)
   — the per-chrom panels make this lineage difference immediately visible
   when both haplotypes are loaded

## Files modified / added

### Modified

- `Inversion_atlas.html` (+450 lines from turn 115 → 57,094 → 57,544)
  - **Storage layer** (~140 lines, sibling to the existing repeat-density block):
    `_isNcRNADensityJSON`, `_validateNcRNADensityChrom`, `_storeNcRNADensity`,
    `_persistNcRNADensity`, `_restoreNcRNADensity`, `_clearNcRNADensity`,
    `_getNcRNADensity`, `_ncRNADensityChromList`, `_resolveNcRNADensityClass`,
    `_setNcRNADensityActiveClass`, `_restoreNcRNADensityActiveClass`
  - **Renderer** (`_renderNcRNADensityPanel`, ~190 lines): SVG line track
    with class chips, scan-range highlight band (when a candidate is
    active), Y-axis grid, X-axis ticks. Fail-soft on missing data with a
    "drop a `<chrom>_ncrna_density.json` file" hint that lists currently-
    loaded chroms.
  - **Constants**: `NCRNA_DENSITY_TOOL`, `NCRNA_DENSITY_LS_PREFIX`,
    `NCRNA_DEFAULT_CLASS = 'tRNA_all'`, `NCRNA_CLASS_FAMILIES` (the
    three-family chip layout)
  - **Mount point**: `<div id="bndNcRNADensity">` added to
    `renderBoundariesPage`'s template literal, right after `#bndRepeatDensity`
  - **Refresh hook**: one line in `_bndRefreshUI` (next to the existing
    repeat-density call). Covers all 15 call sites of `_bndRefreshUI` for free.
  - **Loader dispatch**: `_classifyJSONKind` gains `ncrna_density` case;
    `loadMultipleJSONs` gains a dispatch block that stores + persists +
    re-renders the boundaries page when a matching JSON is dropped.
  - **Boot-time restore**: two `try { _restoreNcRNADensity }` /
    `try { _restoreNcRNADensityActiveClass }` calls alongside the existing
    cross-species and dotplot restores.
  - **Help-page row**: new table entry documenting `<chrom>_ncrna_density.json`.

### Added

- `aggregate_ncrna_density.R` (NEW, 226 lines) — sibling to
  `aggregate_repeat_density.R`. CLI:
  ```
  Rscript aggregate_ncrna_density.R \
    --trna-gff3 fClaHyb_Gar.tRNA.gff3 \
    --rrna-gff3 fClaHyb_Gar.rRNA.gff3 \
    --ncrna-gff3 fClaHyb_Gar.ncRNA.gff3 \
    --windows-tsv windows.tsv \
    --species "C. gariepinus" \
    --out-dir out/
  ```
  Reads tRNAscan-SE / barrnap / Rfam GFF3s; classifies into sub-categories
  using attribute heuristics (Pseudo flag, Name= patterns, Rfam type
  column); emits one `<chrom>_ncrna_density.json` per chromosome present
  in the windows grid.

- `sample_ncrna_density_LG28.json` — synthetic fixture for testing the
  loader without running the R aggregator. Contains a realistic rDNA
  cluster spike in window 3 to verify the panel surfaces structural
  signal correctly.

## Cache + reactive behavior

- Per-chrom storage on `state.ncRNADensity[chrom]` (not a global blob —
  parallels `state.repeatDensity[chrom]`)
- Active-class selection: `state.ncRNADensityActiveClass[chrom]`,
  localStorage-backed
- Class-switch click reactively re-renders just the panel (no full
  page-11 rebuild)
- Boot-time restore picks up previously-loaded chroms automatically

## Tests run

End-to-end Node integration test, simulating the full loader → render
flow with a 5-window fixture:
- Detector recognizes `tool === 'ncrna_density_v1'` ✓
- Store + persist round-trips through localStorage; restore picks it up ✓
- Renderer produces 4.9 KB of SVG with chrom label, three family chip rows,
  and the right Y-axis class label ✓
- Class switching reactively re-renders with new class label ✓
- Empty-state hint shows when chrom isn't loaded ✓
- `_classifyJSONKind` ordering preserved: cross_species → dotplot_mashmap
  → ncrna_density → repeat_density → chromosome ✓

`node --check` on the entire concatenated inline JS of the patched
atlas: no syntax errors.

## What's deferred to a polish turn

Per spec §0.4, the canonical TE panel has features the MVP renderer does
NOT yet replicate:
- Y-axis mode dropdown (linear / log / q99 / q95 / q90 / auto)
- Full-chrom vs zoomed-to-scan-range view toggle
- Pan/zoom inside zoomed view
- Double-click zoom

These are pure rendering polish and can be cloned in a follow-up if the
MVP isn't expressive enough in practice. None of the architecturally
important pieces (storage, persistence, mount, refresh, loader chain,
class switching) are missing.

## What this unblocks

- **Focal-vs-bg widget** (next turn candidate): once that lands, its
  metric pool will auto-discover `tRNA_all`, `rRNA_all`, `ncRNA_all` (and
  the sub-classes) as available metrics with **zero atlas changes** — same
  pattern the spec describes for `state.data.tracks` discovery.
- **Page-16 flanking chips** (`★ tRNA-cluster`, `★ rDNA-flank` per spec
  §5.4): now the data is loadable; the chip rendering is a small
  follow-up to `STEP_CS01_extract_breakpoints.py` (Python pipeline-side)
  + a tiny atlas patch.

## Spec status

`SPEC_ncrna_density_layer.md`: **PARTIALLY IMPLEMENTED in turn 116**.
Architecture + storage + MVP renderer + R aggregator shipped. Move to
`specs_in_progress/` (or keep in `specs_todo/` until polish + page-16
chip integration ships).

## Next turns

The keystone: **`SPEC_focal_vs_background_widget.md`**. With ncRNA
density now wired, focal-vs-bg gets four free metrics (F_ST, θπ-ratio,
Hobs/Hexp, |Z|) plus six new ones via auto-discovery (`tRNA_all`,
`rRNA_all`, `rRNA_5S`, `rRNA_18S`, `ncRNA_all`, `ncRNA_miRNA`) the
moment it's wired. That's the largest single payoff turn left in the
todo folder.
