# SPECS_TIER_INDEX — Sorted by literature anchor

**Folder layout** (after turn 131):
- `specs_todo/` — active build queue (already reviewed by Quentin)
- `specs_new_turn131/` — pending review queue (this folder; 19 new specs from turn 131)

When a spec from `specs_new_turn131/` is approved for build, MOVE it to
`specs_todo/`. When fully shipped, move to `specs_todo/_archive/specs_done/`.

**Convention**:
- **Tier 1** = spec has direct literature reference (Porubsky 2022,
  Hudson 1992, Bhatia 2013, Hartigan 1985, Crane 2015, Corbett-Detig &
  Hartl 2012, Skotte 2013, etc.). The science is anchored.
- **Tier 2** = spec is atlas/UI engineering. No direct lit anchor;
  useful but build for ergonomics or pipeline structure.

Both tiers are needed. Tier 1 specs are higher priority because their
correctness is constrained by published methods — getting them wrong
means publishing wrong numbers. Tier 2 specs are about helping the
user navigate and review; mistakes are recoverable.

---

## Tier 1 — Literature-anchored

These specs implement methods from peer-reviewed literature. Build
quality matters because the output is publication numbers.

### Statistical / population-genomic methods
- **`SPEC_inversion_age_origin_atlas.md`** — Porubsky 2022 (Cell), 
  Hartigan 1985 (Annals of Statistics), Corbett-Detig & Hartl 2012
  (Genetics). GDS classification, dip test, age proxy.
- **`SPEC_recombinant_dosage_changepoint_detector.md`** — changepoint
  literature (PELT, binary segmentation). Per-carrier dosage step
  detection.
- **`SPEC_boundary_consensus_aggregator.md`** — KDE + multimodal mode
  detection. Handles the unimodal/bimodal carrier-distribution problem.
- **`SPEC_hypothesis_test_framework_atlas.md`** — T1-T11 hypothesis
  tests with BH correction. Mendelian transmission, ancestry jackknife
  (Cheat 6), Clair3 indel concordance (Cheat 9), theta-het prior
  (Cheat 12).
- **`SPEC_boundary_confirmation_5track.md`** — Bhatia 2013 (Hudson Fst),
  Crane 2015 (insulation). The "convicting" composite figure.
- **`SPEC_metric_overlay_priors.md`** — multi-metric convergence,
  cites the flashlight cheats inventory.
- **`SPEC_age_origin_panel.md`** (existing) — partly overlaps with
  `SPEC_inversion_age_origin_atlas.md`; consolidate later.

### Marker / breeding methods
- **`SPEC_marker_panel_design_atlas.md`** — private indel + dosage tier
  system (chat `1b4d8e12`). Cross-species F1 hybrid context (Cgar male
  × Cmac female) determines control structure.
- **`SPEC_per_candidate_breeding_readiness_card.md`** — Atlas 5 Part B
  from chat `c03fc41e`. Pairing implications from F_ROH + damaging-load
  asymmetry.

### R-side producer specs (atlas hero already consumes contract)
- **`SPEC_STEP_T05_theta_cusum_emitter.md`** — turn 132. Per-carrier
  CUSUM on θπ, produces `cusum_theta` layer that the atlas page-12
  hero (`_drawThCusumHero`, shipped turn 132) paints from. Wraps the
  walked-back `lib_persample_cusum.R` from chat `487c7f04`.

### Inversion-discovery methods
- **`SPEC_distant_band_concordance_fish_trajectory.md`** (existing,
  partly shipped) — band-trajectory clustering across candidates.
- **`SPEC_l2_sweep_inheritance.md`** (existing) — Cramér's V on L2
  envelopes; the auto-promote producer.
- **`SPEC_sliding_window_inheritance.md`** (existing) — fixed-tile
  Jaccard at fine resolution.
- **`SPEC_karyotype_per_interval_intersection.md`** — post-processing
  GHSL by triangle interval (chat `5b793a68`).

### Surface specs that consume Tier 1 producers
- **`SPEC_observable_allele_h_label_system.md`** — turn 137.
  Observable-allele framing for the H-label system (H1, H2, H3 are
  observed allele clusters, not asserted biological haplotypes) plus
  data-driven hom-vs-het band classification using dosage HET fraction
  as a cross-check. Replaces the fixed `_KARYO_DETAILED_LABELS` table
  with per-candidate labels that respect missing classes (e.g., we may
  never observe `H1/H1`). Atlas Slice 1 ~1.5 turns; companion R-side
  producer spec (`SPEC_STEP_HL01_h_label_classifier.md`) not yet drafted.
- **`SPEC_per_band_evidence_layer.md`** — per-band H, θπ, F_ROH,
  family. F_ROH-confounder + family-confounder flags from the
  hatchery genomics literature on confounded inversions.
- **`SPEC_multi_evidence_regime_framework.md`** — doctrine spec
  (chat `47cd29b9`). Treats K-means as geometric, requires multi-layer
  agreement for biological regime calls.

### Cross-species / synteny
- **`SPEC_comparative_te_breakpoint_fragility.md`** (existing) — TE
  enrichment at inversion breakpoints (NAHR mechanism literature).
- **`SPEC_cross_species_dotplot.md`** (existing) — wfmash-based
  alignment.
- **`SPEC_busco_anchors.md`** (existing) — synteny anchors from BUSCO.
- **`SPEC_phylogenetic_tree_integration.md`** (existing) — Dollo
  parsimony for inversion age across phylogeny.

---

## Tier 2 — Atlas/UI engineering

Build for ergonomics and pipeline structure. No direct lit anchor;
mistakes are recoverable.

### Master vision + orchestration
- **`SPEC_MASTER_full_automatic_pipeline.md`** — umbrella spec.
- **`SPEC_multichrom_load_orchestrator.md`** — bulk-load 28 JSONs.
- **`SPEC_genome_wide_ideogram.md`** — 28-chrom ideogram page.
- **`SPEC_cross_chromosome_lineages.md`** — Hungarian-aligned lineages
  across chromosomes.
- **`SPEC_manuscript_bundle_export.md`** — Markdown + TSV + JSON
  bundle export.

### Review surfaces
- **`SPEC_g_panel_unified_groups.md`** — popup with karyotype/inheritance/
  manual/auto/lineages tabs.
- **`SPEC_review_surfaces_auto_and_lineages.md`** — dashed-row UI for
  algorithm-proposed candidates.
- **`SPEC_l3_het_dosage_coloring.md`** (Slice 1 shipped turn 129) —
  het-rate coloring on L3 PCA dots.
- **`SPEC_lasso_inheritance_backgrounds.md`** — alpha-shaded intervals
  on per-sample-lines.
- **`SPEC_interval_collector.md`** — passive accumulator from tracked
  samples.
- **`SPEC_lines_panel_candidate_bands.md`** (existing).
- **`SPEC_l2_envelope_live_split.md`** (existing).

### Diagnostic / decision-gated
- **`SPEC_inheritance_diagnostic_protocol.md`** (existing).
- **`SPEC_inheritance_unification.md`** (existing, decision-gated).
- **`SPEC_focal_vs_background_widget.md`** (existing).

### Visual layers / overlays
- **`SPEC_gene_annotation_overlay.md`** — GFF3 integration.
- **`SPEC_ld_decay_overlay.md`** (existing).
- **`SPEC_ncrna_density_layer.md`** (existing).

### Evidence pages
- **`SPEC_sv_evidence_page.md`** (existing) + sub-specs
  (`per_candidate_folder.md`, `upset_redirect.md`).
- **`SPEC_hypothesis_registry_and_multispecies.md`** (existing) —
  hypothesis test results catalogue.

### Multi-species architecture
- **`SPEC_OVERVIEW_multispecies_architecture.md`** (existing).

---

## Build priority order

When the next chat decides what to build:

1. **L2-sweep auto-promote Slice 1** (Tier 1 producer; enables auto-review).
2. **G-panel scaffold Slice 1** (Tier 2 review surface; enables auto/lineages tabs).
3. **Trajectory matrix viewer Slice 3** (Tier 1 surface for already-shipped lineage compute).
4. **Cross-chromosome lineages** (Tier 2 producer; needs (1) and (3)).
5. **Multi-chrom load orchestrator** (Tier 2 infrastructure; enables genome-wide ideogram).
6. **Genome-wide ideogram** (Tier 2 navigation surface).
7. **Marker panel design atlas** (Tier 1 — has lit anchor, manuscript-ready).
8. **Per-candidate breeding card** (Tier 1 — converts paper to resource).
9. **Manuscript bundle export** (Tier 2; integrates everything).

Tier 1 producers above Tier 2 surfaces in the same priority slot.
Quentin can override if the manuscript deadline weights things differently.
