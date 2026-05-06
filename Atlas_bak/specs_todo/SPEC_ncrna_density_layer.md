# SPEC — rRNA / tRNA / ncRNA density layer (boundaries + cross-species)

**Status**: forward-looking spec. Not implemented. Created same session as
`SPEC_focal_vs_background_widget.md`, `SPEC_xpehh_track.md`,
`SPEC_ld_decay_overlay.md`.

**Reading order**: this spec → atlas `_renderRepeatDensityPanel` (line ~16223,
the existing TE-density panel this spec mirrors) → atlas `state.repeatDensity`
storage pattern (line ~12410, the data shape this spec mirrors) →
`SPEC_comparative_te_breakpoint_fragility.md` (sibling spec with same
"per-chrom-density JSON layer" pattern).

**Reference biology**: tRNA gene clusters and rDNA arrays are independent
fragility signals from TE density. tRNA clusters are documented breakpoint
hotspots in vertebrate evolution (Larkin 2009 EBR work, Bermudez-Santana
2010); rDNA arrays mediate recombination through their tandem-repeat
architecture and are independently associated with structural rearrangements
in fish lineages. Adding these as a sibling track to the existing
TE-density panel gives the boundaries page a complete fragility picture,
not just the EDTA TE view.

**One-line summary**: data already exists on disk; this spec wires it into
the atlas as a density panel parallel to the existing repeat-density panel
on page 11, plus a flanking-density annotation column on page 16's
breakpoint catalogue.

---

## 0. What already exists on disk (confirmed via past chat search)

The user's annotation pipeline has already produced canonical files at:

```
fClaHyb_Gar/
├── fClaHyb_Gar.tRNA.gff3        tRNAscan-SE, 15,852 loci
├── fClaHyb_Gar.rRNA.gff3        barrnap, 4,286 loci (5S × 4,272, 5.8S × 4, 18S × 5, 28S × 5)
├── fClaHyb_Gar.ncRNA.gff3       Rfam/Infernal, ex-tRNA/rRNA: snRNA, snoRNA, miRNA, etc.
├── fClaHyb_Gar.miRNA.gff3       curated (FishmiRNA + Rfam supported)
└── analyses/ncRNA/              raw tRNAscan struct, barrnap raw, Rfam .tblout, intermediates

fClaHyb_Mac/
├── fClaHyb_Mac.tRNA.gff3        tRNAscan-SE, 15,699 loci
├── fClaHyb_Mac.rRNA.gff3        barrnap, 12,225 loci (5S × 12,209 — 3× expansion vs Gar!)
├── fClaHyb_Mac.ncRNA.gff3       Rfam/Infernal
└── ...
```

**Sub-categories worth surfacing in the atlas** (preserve the GFF3 `feature`
or `Name` attribute distinction):
- tRNA: `high_confidence` vs `pseudogene` vs `intronic`. tRNAscan-SE flags
  these natively; the GFF3 has them as attributes
- rRNA: `5S`, `5.8S`, `18S`, `28S`, `partial_fragment`. The 3× 5S expansion
  in Mac vs Gar is biologically interesting — visualize separately
- ncRNA: keep Rfam family ID (`RF00001`, `RF00005`, etc.) as the
  sub-category; group display by major class (snRNA, snoRNA, miRNA, etc.)

This spec does NOT generate the GFF3s — they're already shipped. It only
specifies (1) a small aggregator that converts GFF3 → per-window density
JSON parallel to TEfull, and (2) the atlas-side wiring.

---

## 1. Purpose

Surface rRNA and tRNA density as first-class atlas tracks, parallel to
the existing TE density:

1. **On page 11 boundaries**: a sibling density panel below the existing
   repeat-density panel, with class chips for tRNA / rRNA-5S /
   rRNA-large (5.8S+18S+28S) / ncRNA-other. Same panel shape as
   `_renderRepeatDensityPanel`, same scan-radius integration.

2. **On page 16 cross-species**: per-breakpoint flanking annotation
   showing tRNA cluster proximity and rDNA cluster proximity, displayed
   in the catalogue chips next to the existing all_TE chip. A breakpoint
   that overlaps a tRNA cluster gets a `★ tRNA-cluster` chip; one within
   ±100 kb of an rDNA array gets a `★ rDNA-flank` chip.

3. **Statistical hook**: the focal-vs-bg widget (per
   `SPEC_focal_vs_background_widget.md`) auto-discovers
   `trna_density` / `rrna_density` as additional metrics — same pipeline,
   no widget code changes.

---

## 2. Why this matters

Adding rRNA/tRNA density to the boundaries page changes what the atlas
can claim about a breakpoint:

| Evidence layer | What it supports |
|---|---|
| EDTA TE density (existing) | TE-mediated fragility — repeat-mediated recombination, the canonical Spalax mechanism |
| **tRNA cluster density (new)** | Vertebrate EBR-style breakpoint reuse — tRNA clusters are documented evolutionary breakpoint hotspots |
| **rDNA array density (new)** | rDNA-mediated recombination — independent mechanism, common in fish lineages, and Mac shows a 3× 5S expansion that's biologically interesting on its own |
| ncRNA density (new) | Functional cargo — does the inversion carry regulatory ncRNAs? |

Without these, the manuscript can only argue TE-mediated mechanism. With
them, the manuscript can present a layered fragility picture and triage
candidates by which fragility signal dominates.

---

## 3. What this layer DOES / does NOT claim

### DOES

- Per-window densities (loci/Mb and bp-occupancy/Mb) for tRNA, rRNA, and
  ncRNA classes
- Sub-class breakdown (rRNA-5S vs rRNA-large; tRNA-high-confidence vs
  tRNA-pseudogene)
- Flanking-density annotation per breakpoint on page 16, parallel to the
  existing TE flanking annotation

### Does NOT

- Causation — high tRNA density at a breakpoint ≠ tRNA cluster caused the
  rearrangement
- Phylogenetic age — these are present-day densities; age is a separate
  question (see `SPEC_phylogenetic_tree_integration.md`)
- Expression / functional impact — that's a downstream RNA-seq question

---

## 4. Vocabulary discipline

| Allowed | Discouraged |
|---|---|
| "tRNA cluster density elevated near breakpoint" | "tRNA cluster caused breakpoint" |
| "Inversion span overlaps rDNA array" | "rDNA recombination produced inversion" |
| "Mac shows 3× 5S rRNA expansion vs Gar — independent of breakpoint analysis" | "rDNA expansion drives chromosome rearrangement" |
| "Consistent with tRNA-cluster–mediated breakpoint reuse" | "tRNA cluster reuse confirmed" |

---

## 5. Architecture

### 5.1 Aggregator — convert GFF3 → per-window density JSON

A small R or Python script `aggregate_ncrna_density.R` (sibling to
`aggregate_repeat_density.R`) that:

1. Reads `<species>.tRNA.gff3`, `<species>.rRNA.gff3`,
   `<species>.ncRNA.gff3`
2. Bins loci into the same per-chrom windows used by the TE density
   layer (window size, step from chrom precomp)
3. Computes per-window:
   - `n_loci_per_window`: count of loci whose midpoint falls in the window
   - `bp_occupancy_per_window`: total bp covered by features in the window
   - `density_per_mb`: `n_loci_per_window / window_size_mb`
4. Stratifies by sub-category:
   - tRNA: `tRNA_all`, `tRNA_HC` (high-confidence), `tRNA_pseudo`,
     `tRNA_intronic`, `tRNA_Sec` (selenocysteine — there's exactly 1 in
     each haplotype, worth keeping visible)
   - rRNA: `rRNA_all`, `rRNA_5S`, `rRNA_5_8S`, `rRNA_18S`, `rRNA_28S`,
     `rRNA_partial`
   - ncRNA: `ncRNA_all`, `ncRNA_snRNA`, `ncRNA_snoRNA`, `ncRNA_miRNA`,
     `ncRNA_other`
5. Emits JSON with the same shape as the existing TEfull layer, just
   with a different category set under `by_class`:

```jsonc
// File: <chrom>_ncrna_density.json
{
  "schema_version": 1,
  "tool": "aggregate_ncrna_density.R v1",
  "generated_at": "2026-MM-DD",
  "input_files": [
    {"path": "fClaHyb_Gar.tRNA.gff3", "sha256": "...", "n_loci": 482},
    {"path": "fClaHyb_Gar.rRNA.gff3", "sha256": "...", "n_loci":  87},
    {"path": "fClaHyb_Gar.ncRNA.gff3", "sha256": "...", "n_loci": 1043}
  ],
  "chromosomes": [
    {
      "chrom": "LG28",
      "window_centers_mb": [...],
      "window_size_bp": 50000,
      "n_windows": 4302,
      "by_class": {
        "tRNA_all":     { "densities": [...], "n_loci_per_window": [...] },
        "tRNA_HC":      { "densities": [...], "n_loci_per_window": [...] },
        "tRNA_pseudo":  { "densities": [...], "n_loci_per_window": [...] },
        "rRNA_all":     { ... },
        "rRNA_5S":      { ... },
        "rRNA_18S":     { ... },
        "ncRNA_all":    { ... },
        "ncRNA_miRNA":  { ... }
      }
    }
  ]
}
```

The `densities[]` field is `loci_per_Mb` to match the TE layer's
unit-consistency expectation. Atlas readers don't need to translate.

### 5.2 Atlas storage — `state.ncRNADensity[chrom]`

Mirrors `state.repeatDensity[chrom]` exactly. Same loader pattern,
same drag-drop hook, same per-chrom indexing, same fail-soft empty
state when a chrom doesn't have a JSON loaded. **Use the same loader
function** and just dispatch by filename suffix (`*_repeat_density_TEfull.json`
→ `state.repeatDensity`; `*_ncrna_density.json` → `state.ncRNADensity`).

### 5.3 Page 11 — sibling panel below repeat density

```html
<div id="bndNcRNADensity" class="bnd-rd-wrap"
     style="margin-top: 14px;"></div>
```

`_renderNcRNADensityPanel()` is a near-clone of `_renderRepeatDensityPanel`:
- Same SVG layout (margins, axes, scan-range zoom)
- Same y-mode dropdown (linear / log / q99 / q95 / q90 / auto)
- Same view-mode toggle (full / zoomed)
- Same scan-radius integration (reads from
  `_ensureBoundariesState().scan_radius_bp`)
- Different class dropdown — populated from the ncRNA `by_class` keys
  instead of TE classes
- Different default class — `tRNA_all` (most informative breakpoint
  context)

The two panels stack vertically. Both update on radius change. Both
share the highlighted scan range.

### 5.4 Page 16 — flanking annotation in the catalogue

The existing breakpoint catalogue chips show event_type + Spalax all_TE
flag. Add two new chip slots:

- `★ tRNA-cluster` (orange) when ±100 kb flank around the breakpoint
  contains ≥ N tRNA loci (default N = 5, configurable via
  `state.csUI.trnaClusterMinLoci`)
- `★ rDNA-flank` (purple) when ±100 kb flank around the breakpoint
  contains ≥ 1 large-rRNA locus (18S/28S/5.8S — the structural rDNA
  array, NOT the dispersed 5S)

The flanking annotation is computed at JSON-load time, same as the
existing TE flanking annotation. Add a section to
`STEP_CS01_extract_breakpoints.py` (the existing pipeline script
the user pasted) that reads the ncRNA GFF3 alongside the TE JSON and
emits `flanking_ncrna_gar` / `flanking_ncrna_mac` per breakpoint, same
shape as `flanking_repeat_density_gar` / `_mac`.

**This is the only upstream-pipeline change in the spec.** Atlas reads
the new fields if present; gracefully degrades to "no ncRNA flanking
loaded" hint if absent.

### 5.5 Focal-vs-bg widget integration (free)

The widget auto-discovers any per-window scalar metric in
`state.data.tracks` or in a separate density layer. Adding tRNA / rRNA
density as widget metrics is **zero atlas code** once the JSON layer
loads — the widget's metric dropdown grows two new entries:
`trna_density (loci/Mb)` and `rrna_density (loci/Mb)`. Wilcoxon focal
vs bg works the same way as for F_ST.

This is the single biggest payoff of mirroring the existing TE-density
shape: every downstream consumer that already worked for TEs starts
working for tRNA/rRNA on day 1.

### 5.6 Help-page entry

Add to the help page glossary:
- "tRNA density" — tRNAscan-SE high-confidence + pseudogene loci per Mb
- "rRNA density" — barrnap detections, broken into 5S vs large-subunit
- "rDNA-flank chip" — present when breakpoint is within 100 kb of a
  large-subunit rRNA locus
- Citation row: tRNAscan-SE (Lowe & Eddy 1997, Chan & Lowe 2019);
  barrnap (Seemann); Rfam (Kalvari et al. 2021)

---

## 6. Module structure

```
phase_X_ncrna_density/
  README.md
  aggregate_ncrna_density.R     # NEW — GFF3 → per-window density JSON
  STEP_NC_01_aggregate.sh       # wrapper, one call per haplotype
  STEP_NC_02_validate_counts.R  # sanity check: total loci match GFF3 row count

atlas/
  Inversion_atlas.html          # MODIFIED — bndNcRNADensity panel mount,
                                # _renderNcRNADensityPanel function,
                                # csUI catalogue chip rendering update
  tests/
    test_ncrna_density.html     # NEW — panel renders + flanking chip logic

cross_species_breakpoints/
  STEP_CS01_extract_breakpoints.py   # MODIFIED — new --gar-ncrna-gff3 +
                                     # --mac-ncrna-gff3 args, emit
                                     # flanking_ncrna_gar / _mac
```

---

## 7. Tests

### Aggregator (`aggregate_ncrna_density.R`)

- `aggregate_trna_total_matches_gff3`: total `n_loci` summed across all
  chrom × all windows == total tRNA rows in the GFF3 (within tolerance
  for boundary-window double-counting on overlapping features)
- `aggregate_rrna_subclasses_match_canonical`: 5S + 5.8S + 18S + 28S
  counts in the JSON match the 4,272 + 4 + 5 + 5 = 4,286 from the user's
  pipeline
- `density_unit_is_loci_per_mb`: 5 loci in a 50 kb window → density = 100.0

### Atlas-side tests

- `ncrna_panel_mounts_below_te_panel`: page 11, ncRNA JSON loaded, both
  panels visible
- `ncrna_panel_radius_change_propagates`: change scan radius, ncRNA panel
  re-renders to match
- `cs_catalogue_trna_cluster_chip_renders`: breakpoint with ≥ 5 tRNA in
  flank → chip visible
- `cs_catalogue_rdna_flank_chip_renders`: breakpoint within 100 kb of an
  18S locus → chip visible
- `focal_vs_bg_widget_auto_discovers_trna`: load JSON, widget metric
  dropdown shows `trna_density`

Target: 8–10 tests.

---

## 8. Open questions / explicitly deferred

- **Centromere proximity track.** Centromeres are also breakpoint hotspots,
  and `fClaHyb_Gar.centromeres.bed` already exists. Same shape as this
  spec — could be added as a parallel panel (4th density: TE / tRNA /
  rRNA / centromere-distance). Deferred to its own one-line spec.
- **Telomere proximity** — `fClaHyb_Gar.telomeres.bed` exists. Same idea.
  Deferred.
- **Segmental duplications** — `fClaHyb_Gar.SD.biser.bed` exists. SD
  density at breakpoints is a *very* strong fragility signal in mammals
  (Bailey et al. 2002, evolutionary breakpoint regions). Worth adding as
  a sibling layer. Deferred to its own spec — this spec stays focused on
  ncRNAs as Quentin asked.
- **5S rRNA expansion in Mac.** The 3× expansion is interesting on its
  own. The atlas could add a small comparative panel (Gar vs Mac 5S
  density distribution per chrom). Deferred — out of scope for the
  boundary-page integration but flagged here so it doesn't get lost.
- **Per-isotype tRNA tracks** (Ala, Arg, ...). 20+ amino-acid tracks would
  clutter the panel. Day-1 default is `tRNA_all`; per-isotype is in
  `by_class` for users who want it via the dropdown.

---

## 9. Summary

Data already exists at canonical paths (`fClaHyb_*.tRNA.gff3`,
`fClaHyb_*.rRNA.gff3`, `fClaHyb_*.ncRNA.gff3`). This spec wires them
into the atlas: one new aggregator script, one parallel density panel
on page 11, two new flanking chips on page 16, free metrics in the
focal-vs-bg widget. Mirrors the existing repeat-density layer shape
exactly so all downstream consumers work without modification.

End of spec.
