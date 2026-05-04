# Specs to-do — atlas + pipeline forward-looking specs

This folder collects forward-looking specs that have NOT yet been implemented
on the atlas or the cluster-side pipelines. Each spec is self-contained,
written so a future Claude session (or Quentin) can pick it up cold.

## Adding a new spec

- One Markdown file, named `SPEC_<short_id>.md`
- Top of the file: `**Status**: forward-looking spec. Not implemented.`
- Top of the file: `**Reading order**: this spec → ...` (give the reader
  the prerequisite documents)
- Vocabulary discipline section if the spec touches manuscript-grade claims
- Tests section listing concrete test names that would gate "done"
- Open questions / deferred list at the bottom

The pattern matches the existing `SPEC_comparative_te_breakpoint_fragility.md`
and `SPEC_hypothesis_registry_and_multispecies.md` in the parent set.

## Current contents

| Spec | One-line | Blocking dependency |
|---|---|---|
| `SPEC_focal_vs_background_widget.md` | Spalax-style focal-vs-bg widget for page 16 + page 11 (Z-F_ST / θπ / Hobs); forward-compat with multi-cohort | none day-1 (live mode); aggregated JSON layer optional for manuscript-locked P |
| `SPEC_ld_decay_overlay.md` | Spalax Fig. S62 idiom — three-trace LD decay overlay (focal / chrom / genome) for page 16 + page 11; reuses entire fast_ld stack | none day-1 (live mode); precomputed `ld_decay_genome_wide.json` optional for the orange trace |
| `SPEC_ncrna_density_layer.md` | rRNA/tRNA/ncRNA density panel on page 11 + flanking chips on page 16; data already exists at `fClaHyb_*.{tRNA,rRNA,ncRNA}.gff3` | aggregator script `aggregate_ncrna_density.R` (small) |
| `SPEC_cross_species_dotplot.md` | Mummerplot-style dot plot on page 16 with hover-enlarge popup + click-to-pin; reads wfmash synteny_blocks (already in cs_breakpoints v2) + mashmap multi-resolution grid | mashmap pre-compute (~10s of seconds, one-shot) for the high-res tier |
| `SPEC_xpehh_track.md` | XP-EHH per-window track on popstats (and later ancestry) page | phased VCF for 226-sample cohort + choice of reference cohort |
| `SPEC_inversion_age_atlas_surface.md` | Two-task atlas surface for inversion age. **Task A** (between-inversion ranking → page 5 sortable column) and **Task B** (within-arrangement ranking, which group is older → page 3 candidate-focus panel). Calibrated with μ_year = 3×10⁻⁹/site/year (Liu et al. 2023 *M. electricus*, Siluriformes-anchored r8s estimate); generation time NOT assumed (sidesteps hatchery-vs-wild ambiguity). Consumes existing `STEP_C01f_c_burden_regression.R` v9.7 outputs via a small new emitter script. Soft-anchor framing throughout (μ uncertain, ranking robust). | one new HPC emitter script `STEP_C01f_d_emit_age_json.py` (~80 lines); atlas-side wiring is the bulk of the work |

## Companion specs (in the parent set, already exist)

- `SPEC_OVERVIEW_multispecies_architecture.md` — top-level multi-species framing
- `SPEC_hypothesis_registry_and_multispecies.md` — registry data model + classification UI
- `SPEC_busco_anchors.md` — BUSCO anchor layer
- `SPEC_phylogenetic_tree_integration.md` — tree/age polarization layer
- `SPEC_comparative_te_breakpoint_fragility.md` — multi-species TE density layer

## When a spec lands

When a spec is implemented and shipped, move it OUT of `specs_todo/` and
into a `specs_shipped/` folder (or wherever the project keeps post-ship
documentation). The `Status:` line at the top should change from
`forward-looking spec. Not implemented.` to `IMPLEMENTED in turn N (date)`.
