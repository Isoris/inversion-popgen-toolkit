# specs_new_turn131/ — 19 new specs from turn 131 spec-writing pass

These specs were written turn 131 (2026-05-05) in a single spec-writing
pass. They are kept SEPARATE from `specs_todo/` because that folder
already contains many specs that overlap with shipped atlas features
or shipped LANTA pipelines — putting new specs there pollutes it.

**Read `SPECS_TIER_INDEX.md` first** — it sorts everything by tier
(literature-anchored vs UI-engineering) and gives the recommended
build priority order.

## Inventory

### Master / orchestration (Tier 2)
- `SPEC_MASTER_full_automatic_pipeline.md` — umbrella vision spec
- `SPEC_multichrom_load_orchestrator.md` — bulk-load 28 JSONs
- `SPEC_genome_wide_ideogram.md` — 28-chrom ideogram page F1
- `SPEC_cross_chromosome_lineages.md` — Hungarian-aligned across chroms
- `SPEC_manuscript_bundle_export.md` — .zip bundle export

### Tier 1 (literature-anchored, highest correctness importance)
- `SPEC_recombinant_dosage_changepoint_detector.md` — replaces PCA-based recomb
- `SPEC_boundary_consensus_aggregator.md` — KDE stack-then-aggregate
- `SPEC_metric_overlay_priors.md` — 12+ secondary boundary detectors
- `SPEC_karyotype_per_interval_intersection.md` — 30-line post-processing
- `SPEC_hypothesis_test_framework_atlas.md` — T1-T11 with BH, H1-H4
- `SPEC_multi_evidence_regime_framework.md` — K-means provisional doctrine
- `SPEC_per_band_evidence_layer.md` — per-band H/θπ/F_ROH/family + flags
- `SPEC_inversion_age_origin_atlas.md` — Porubsky 2022 + Hartigan dip
- `SPEC_boundary_confirmation_5track.md` — F4 panel A composite
- `SPEC_marker_panel_design_atlas.md` — private indel tier system
- `SPEC_per_candidate_breeding_readiness_card.md` — Atlas 5 Part B
- `SPEC_gene_annotation_overlay.md` — GFF3 integration
- `SPEC_STEP_T05_theta_cusum_emitter.md` — R-side spec (turn 132) for the cusum_theta layer; atlas hero (`_drawThCusumHero`) already paints from this contract

### Tier 2 (UI/orchestration)
- `SPEC_lasso_inheritance_backgrounds.md` — alpha intervals on per-sample lines
- `SPEC_interval_collector.md` — passive accumulator from tracked samples

### Index
- `SPECS_TIER_INDEX.md` — global tier sort across BOTH folders + build priority

## Promotion rules

When a spec from this folder is reviewed and approved by Quentin for
build, MOVE it to `specs_todo/` (not copy). When it's fully shipped,
move it to `../specs_done/` (top-level, not a subdir of `specs_todo/`)
per the unified lifecycle defined in `specs_todo/README.md`.

This keeps `specs_todo/` as the active build queue, this folder as
the pending-review queue, and `specs_done/` as the shipped archive.

## Status of the 4+1 inheritance methods (re-confirmed turn 131)

All 4+1 inheritance methods are spec'd:

| Method | Spec | Status |
|---|---|---|
| 1. Candidate-only | (in atlas as cross-Cramér's V) | SHIPPED earlier |
| 2. Fish-trajectory | `specs_todo/SPEC_distant_band_concordance_fish_trajectory.md` | Slices 1+2 SHIPPED turn 130 |
| 3. L2-sweep | `specs_todo/SPEC_l2_sweep_inheritance.md` | spec only |
| 4. Sliding-window | `specs_todo/SPEC_sliding_window_inheritance.md` | spec only |
| 5. Cross-chromosome | `specs_new_turn131/SPEC_cross_chromosome_lineages.md` | spec only (NEW) |
| +1. Unification | `specs_todo/SPEC_inheritance_unification.md` | decision-gated |

So 2/5 shipped, 3/5 spec-only, +1 deferred. All accounted for, none missing.
