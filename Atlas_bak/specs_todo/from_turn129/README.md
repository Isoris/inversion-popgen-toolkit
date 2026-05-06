# Specs from turn 129 planning bundle

Dropped in from `atlas_handoff_v2_turn129.tar.gz`. Each is forward-looking;
none are implemented in `Inversion_atlas.html` yet.

Read the bundle's plans first (`../../handoffs_to_implement/turn129_bundle/plans/`)
for the priority context that produced these.

## Contents

| Spec | One-line | Blocked by |
|---|---|---|
| `S1_sv_evidence_tables.md` | 5-table SV schema (`sv_variant_catalog`, `sv_sample_genotypes`, `candidate_group_membership`, `sv_group_genotype_counts`, `sv_group_enrichment_results`) | SV calling outputs from MODULE_4 must be normalized into these tables first |
| `S2_indel_slope_burden.md` | Indel slope / burden layer placement; framed as "indel-burden support" not breakpoint proof | indel/burden script on HPC + JSON shim |
| `S3_breakpoint_scoring_bayes.md` | Bayesian beta-binomial breakpoint scoring framework; MAPQ0 use, evidence clustering | future module — needs phased-read evidence arrays first |
| `S4_double_crossover.md` | Recombinant / double-crossover detection extension | depends on S3 |
| `S5_sv_interpretation_rules.md` | Discipline doc: INV ≠ raw read, BND not a structural variant claim, evidence tiers | none — pure documentation |
| `S6_dosage_heatmap_streaming_viewer.md` | Dosage heatmap bridge (P4.1 minimum) + standalone streaming viewer with 6 sampling modes (raw / even / random / variance / hybrid / aggregate) | P1.5 stable server startup for the bridge slice |
| `S7_karyotype_breakpoint_internal_evidence.md` | Karyotype-stratified Bayesian scoring of boundary clusters AND paired internal-transition (recombinant) events. Builds on existing `STEP_03_per_read_evidence.py` ledger; consolidates user prompts on boundary scoring + double crossover. **Concrete extension** of S3 + S4 — defines architecture, vocabulary alignment (INV/HET/REF ↔ H1/H1, H1/H2, H2/H2), and required STEP_03 upgrades (MAPQ0 stream + cov_local emit) | LANTA validation of MAPQ0 + cov_local upgrades; one fully-junction-validated candidate for β prior calibration |

## Companion patches

The patches in `../../handoffs_to_implement/turn129_bundle/patches/` are the
*landed-soon* counterparts of these specs. In particular:

- `P4_1_dosage_chunk_endpoint.md` is the minimum-viable slice of `S6` — the
  bridge endpoint only, no preprocessing, no standalone viewer. Implement
  P4.1 first; grow into S6 later.
- `P5_1_boundaries_sv_skeleton.md` consumes the eventual outputs of `S1`
  (SV evidence tables) — the page skeleton stays empty-state until S1 lands.
- `P6_1_karyo_tier_catalogue_merge.md` is independent of these specs.

## Status discipline

When a spec lands as code, move it from `specs_todo/from_turn129/` into
`specs_done/` (create folder when first move happens) with a one-line
"shipped in turn N" header at the top. Don't delete; keep the spec as the
historical contract.
