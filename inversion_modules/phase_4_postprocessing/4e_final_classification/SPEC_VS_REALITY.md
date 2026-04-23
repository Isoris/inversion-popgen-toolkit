# Spec vs reality — what Phase 4e expects vs what gets written

- **Generated:** 2026-04-17 (chat 10)
- **Source of truth:** `compute_candidate_status.R::build_key_spec()` (FIX 42, chat 9)
- **v2 spec target:** 352 keys. Current spec list: **367 keys** (the extra 15 are
  intentional dual aliases for v9↔v10 migration — e.g. `q6_n_HOM_REF` +
  `q6_n_HOM_STD`, `q7_t9_jackknife_verdict` + `q7_t9_jackknife_status`,
  `q1_ancestry_composite_flag` + `q1_composite_flag`. A given value is
  written once under its canonical name; the alias exists so reader code
  that hasn't migrated keeps working.)

## Purpose

This document is the human-readable view of the `*_aspir` lists inside
`build_key_spec()`. It tells you, per Phase-4e spec key, whether a
current-tarball writer produces that key today or whether wiring the
writer is TODO work. The `*_aspir` lists in code are the
machine-readable version; this doc is the per-module view for
prioritising writer-wiring chats.

The completion denominator in `compute_completion()` excludes
aspirational keys, so "72% of Q7 complete" on a candidate today means
72% of *wired* Q7 keys are populated — not 72% of the aspirational v2
target. That's by design: aspirational keys need writers before they can
become meaningful completion signal.

## Summary

| Question | Spec total | Wired | Aspirational | % wired |
|---|---:|---:|---:|---:|
| Q1 | 49 | 49 | 0 | 100.0% |
| Q2 | 55 | 55 | 0 | 100.0% |
| Q3 | 74 | 65 | 9 | 87.8% |
| Q4 | 47 | 28 | 19 | 59.6% |
| Q5 | 39 | 28 | 11 | 71.8% |
| Q6 | 28 | 12 | 16 | 42.9% |
| Q7 | 75 | 34 | 41 | 45.3% |
| **Total** | **367** | **271** | **96** | **73.8%** |

**Takeaway:** Q1 + Q2 are completely wired. Q3 is near-done (only
RepeatMasker annotation missing). Q7 has the largest wiring backlog,
dominated by Q7B breakpoint-audit keys and Layer-C/D aspirationals.

## Aspirational keys grouped by producer module

The ordering below reflects the roadmap for wiring work. Items within a
module are the concrete flat-key writes that need adding — the
computation itself usually already exists.

### Phase 2 — upstream evidence modules

**GHSL v6 (MODULE_2E / STEP_C04_snake3_ghsl_v6 + STEP_C04b_snake3_ghsl_classify) — 7 keys**

The GHSL haplotype-contrast pipeline exists on HPC. Chat 14 (2026-04-18)
split it into the heavy-engine / light-classifier pair, added per-
chromosome annot / karyotypes / per-sample panel RDS shards, and wired
consumer reads (`run_all.R` 2d block scoring, `lib_ghsl_confirmation.R`
Tier-3, `STEP_C01d` D10). The per-candidate flat registry keys below are
still aspirational — they require adding a `keys_extracted` section to
a new `existence_layer_c.schema.json` Tier-2 block, or a direct
`add_evidence()` call from the GHSL summary script. Flagged under
Finding BK for chat 15.

- `q7_layer_c_ghsl_detected`
- `q7_layer_c_ghsl_contrast`
- `q7_layer_c_ghsl_n_pass`
- `q7_layer_c_ghsl_pct_pass`
- `q7_layer_c_ghsl_quality`
- `q7_layer_c_partition_stable`
- `q7_layer_c_ghsl_version`

### Phase 3 — breakpoint refinement

**STEP03 Cochran-Armitage (partial — Fisher wired, Armitage aspirational) — 3 keys**

STEP03 runs both Fisher exact and Cochran-Armitage on the 3×2 PE/SR
contingency table. Fisher outputs are registered
(`q7_layer_d_fisher_or`, `q7_layer_d_fisher_p`) and wired through to
`group_validation_gate.R` for the VALIDATED promotion gate. Armitage
outputs are computed in the same script but the flat-key writes haven't
landed in the chat-9 tarball.

- `q7_layer_d_armitage_z`
- `q7_layer_d_armitage_p`
- `q7_layer_d_concordance`

### Phase 4d — group-dependent cheats (module partially uploaded)

**cheat28 tandem-repeat context — 10 keys**

`cheat28_tandem_repeat_context.R` exists in `4d_group_dependent/` and
computes TRF + GMATA outputs at breakpoints. Wiring target: flat-key
extraction from the cheat28 output TSV to the registry.

- `q4_n_tr_at_bp_left`
- `q4_n_tr_at_bp_right`
- `q4_tr_density_ratio_left`
- `q4_tr_density_ratio_right`
- `q4_dominant_tr_class_left`
- `q4_dominant_tr_class_right`
- `q4_dominant_motif_left`
- `q4_dominant_motif_right`
- `q4_mmbir_signal_left`
- `q4_mmbir_signal_right`

**breakpoint_evidence_audit.py (Q7B audit) — 25 keys**

Pipeline partially exists; a subset of Q7B keys are written, others are
computed but not flat-extracted. This is the single largest aspirational
block — 25 keys across carrier-dropout, fragmentation, filter-pass, and
population-prior categories.

Dropout sub-block (7 keys):

- `q7b_delly_raw_carriers`
- `q7b_delly_strict_carriers`
- `q7b_delly_carrier_loss`
- `q7b_manta_raw_carriers`
- `q7b_manta_pass_carriers`
- `q7b_manta_carrier_loss`
- `q7b_expected_dropout_pct`

Fragmentation sub-block (8 keys):

- `q7b_delly_inv_present`
- `q7b_delly_bnd_3to3`
- `q7b_delly_bnd_5to5`
- `q7b_manta_inv_present`
- `q7b_manta_bnd_inv3`
- `q7b_manta_bnd_inv5`
- `q7b_bnd_rescued`
- `q7b_bnd_rescue_concordance`

Filter-pass sub-block (7 keys):

- `q7b_delly_site_in_raw`
- `q7b_delly_site_passes_strict`
- `q7b_delly_site_qual`
- `q7b_delly_site_pe`
- `q7b_manta_site_in_raw`
- `q7b_manta_site_passes_pass`
- `q7b_manta_site_qual`

Population-prior sub-block (3 keys):

- `q7b_pop_prior_applied`
- `q7b_pop_prior_freq_est`
- `q7b_pop_prior_n_rescued`

### Phase 5 / external annotation — requires RepeatMasker + synteny table

**Q3 repeat annotation (RepeatMasker) — 9 keys**

These require RepeatMasker to be run on the Gar subgenome reference. The pipeline
exists but hasn't been run on all chromosomes in the current tarball.
Wiring target: RepeatMasker run + TSV → registry key extraction in C01g.

- `q3_left_repeat_density`
- `q3_left_repeat_class`
- `q3_left_gc_content`
- `q3_right_repeat_density`
- `q3_right_repeat_class`
- `q3_right_gc_content`
- `q3_genome_wide_repeat_density`
- `q3_repeat_enrichment_left`
- `q3_repeat_enrichment_right`

**Q4 TE + GO enrichment — 9 keys**

Require RepeatMasker TE tracks + GO annotation on the Gar subgenome
reference (`fClaHyb_Gar_LG.fa`, the *C. gariepinus* haplotype extracted
from Section 1's haplotype-resolved F₁ hybrid assembly). These were
marked aspirational because TE family / GO annotation coverage on this
reference may be incomplete compared to long-established reference
species.

- `q4_te_family_left`
- `q4_te_family_right`
- `q4_te_enrichment`
- `q4_te_enrichment_fold`
- `q4_te_name_left`
- `q4_te_name_right`
- `q4_te_same_family`
- `q4_te_same_orientation`
- `q4_go_enrichment`

**Q5 Dollo reconstruction — 3 keys**

Require multi-species synteny table. Not available until the
macrosynteny-scan pipeline is wired through the registry.

- `q5_dollo_node`
- `q5_dollo_mya`
- `q5_n_species_sharing`

### Phase 3 / MODULE_3 — theta ladder per arrangement

**Q5 Tajima D + polymorphism counts — 8 keys**

The theta-ladder pipeline (MODULE_3) exists and produces
`thetaStat`-formatted output, but per-arrangement (HOM_REF / HOM_INV)
Tajima D values aren't yet registered as flat keys. Same for
segregating-sites / shared-polymorphism counts. Wiring target:
theta-ladder output → registry key extraction.

- `q5_tajima_D_std`
- `q5_tajima_D_inv`
- `q5_tajima_D_pooled`
- `q5_tajima_D_pooled_note`
- `q5_segregating_sites_std`
- `q5_segregating_sites_inv`
- `q5_fixed_differences`
- `q5_shared_polymorphisms`

**Q6 Tajima D note + selection pattern — 4 keys**

Same theta-ladder source, but under the Q6 namespace for v1
back-compat. Wiring is identical to Q5's; these are essentially
aliases.

- `q6_tajima_D_inside`
- `q6_tajima_D_flanking`
- `q6_tajima_D_note`
- `q6_selection_pattern`

### Phase 4b — frequency schema class-count extraction

**Frequency v2 schema nested class_counts — 6 keys**

`frequency.v2.schema.json` already extracts the top-level class-count
fields as a nested object (`class_counts: {HOM_REF: ..., HET: ...}`), but
the schema's `keys_extracted` directive doesn't flatten them to
individual top-level keys. Either extend the schema to extract them as
flat keys, or have seal explicitly write them.

- `q6_n_HOM_REF`
- `q6_n_HOM_STD`   *(v9 alias of q6_n_HOM_REF)*
- `q6_n_HET`
- `q6_n_HOM_INV`
- `q6_n_Recombinant`
- `q6_n_Unclassified`

**HWE block — 3 keys**

HWE is computed (schema writes `q6_hwe_p`, `q6_hwe_chi2`,
`q6_hwe_verdict`), but the expected-heterozygosity + observed/expected
ratio + deviation values are computed transiently and not emitted as
flat keys. Wiring target: schema extension to retain these three values.

- `q6_expected_het_hwe`
- `q6_observed_expected_ratio`
- `q6_hwe_deviation`

**Per-Q-group frequency scan — 3 keys**

A cross-Q-group frequency scan exists in principle but isn't currently
computed. Requires iterating Q-groups from `sample_registry` and
computing inversion-allele frequencies within each, then summarizing CV
across them.

- `q6_freq_per_qgroup`
- `q6_freq_cv_across_qgroups`
- `q6_jackknife_sensitive_qgroup`

### Chat-9 FIX 43 aspirationals — v9 back-compat placeholders

**T9 jackknife alias keys — 3 keys**

`q7_t9_jackknife_*` are aspirational v9-era placeholders kept for
reader back-compat. The canonical writer is `q6_family_linkage` (seal
initializes, C01f overwrites after jackknife). FIX 43 rewired all
readers to prefer `q6_family_linkage` and fall back to these only if
present. Wiring these is **low priority** — they exist for backward
compatibility with v9 registries, and the semantic info is already
available via `q6_family_linkage`.

- `q7_t9_jackknife_verdict`
- `q7_t9_jackknife_status`
- `q7_t9_max_delta`

**Other T-test concordance keys — 3 keys**

`q7_t8_clair3_concordance` and `q7_t10_theta_concordance` are computed
inside C01f but written as columns in `hypothesis_verdicts.tsv` rather
than as flat registry keys. `q7_fdr_qval` is a similar column-but-not-key
case. Wiring target: `register_C01f_keys` helper extension (Finding W,
chat 9) to cover these columns.

- `q7_t8_clair3_concordance`
- `q7_t10_theta_concordance`
- `q7_fdr_qval`

## Keys WITH writers (reference; NOT aspirational)

A complete list of wired keys would duplicate `build_key_spec()`. The
policy is: a key is wired iff it appears in `spec[[q]]` AND does NOT
appear in `spec[[paste0(q,"_aspir")]]`. Below are the key-producer
bindings for the largest wired blocks, for reader orientation:

**Phase 4a — C01d pass-1 / C01g / existence Layer B**

- Q1 D1–D12 keys (`q1_d01_block_strength` ... `q1_d12_snake_concordance`)
  are written by `STEP_C01d_candidate_scoring_wired_25_v934_registry.R`.
- Q1 composite keys (`q1_composite_score`, `q1_dim_positive`,
  `q1_shape_class`, etc.) are written by C01d and enriched by C01e.
- Q3 boundary keys (`q3_left_*`, `q3_right_*`, boundary concordance) are
  written by
  `STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R`.
- Layer B (`q7_layer_b_*`) keys from flashlight C00 output.

**Phase 4b — decompose / multi_recomb / nested_comp / seal**

- Q2 decomposition keys: all 55 of them, via
  `internal_dynamics.schema.json`, `recombinant_map.schema.json`, and
  `internal_ancestry_composition.schema.json` (FIX 40, chat 9
  consolidated `q1_composite_flag` to a single writer).
- Q1 composite-flag + ancestry-dominant-structure keys from
  `internal_ancestry_composition.schema.json`.
- Q6 `q6_group_validation`, `q6_validation_promotion_cap`,
  `q6_family_linkage` (initial value), `q6_polymorphism_class` from
  seal.
- Q6 `q6_freq_inv`, `q6_freq_class`, `q6_hwe_p`, `q6_hwe_chi2`,
  `q6_hwe_verdict`, `q6_genotype_balance` from
  `frequency.v2.schema.json`.

**Phase 4c — C01f hypothesis tests**

- Q6 `q6_family_linkage` (post-jackknife overwrite of seal's initial
  value) — FIX 43 + FIX 39 (chat 9) fixed the vocabulary mismatch so
  `few_family_contributing` flows through correctly.
- Q6 `q6_group_validation` promotion/demotion via
  `compute_group_validation`.
- Q6 `q6_jackknife_max_delta` from cheat6 output.
- Q7 `q7_verdict`, `q7_verdict_confidence`, `q7_composite_score`,
  `q7_tier`, `q7_dim_positive`, `q7_n_layers_tested`,
  `q7_n_layers_passed`, `q7_independence_class` via
  `register_C01f_keys` (Finding W — external helper, verify it handles
  the 5 new verdict columns).
- Q7 Layer D (`q7_layer_d_fisher_or`, `q7_layer_d_fisher_p`,
  `q7_layer_d_tested`) from Phase 3 STEP03 feeding into C01f.

**Phase 4d — cheat30 + Q5**

- Q5 GDS keys (`q5_gds_gap`, `q5_gds_gap_percentile`,
  `q5_gds_fst_spearman_rho`, etc.) from
  `cheat30_gds_by_genotype.R` (FIX 41, chat 9 fixed the alias scope
  bug so `compute_pairwise_ibs` is findable at module level).
- Q5 `q5_theta_pi_*`, `q5_dxy_*`, `q5_fst_b1b3` from Q5 age sub-block.
- Q5 `q5_origin_class`, `q5_origin_mechanism` from the origin
  classifier.

**Phase 4a pass-2 / 4e**

- Q1 `q1_d08_peel_or_hyp`, `q1_d11_boundary_concordance` (updated from
  pass-1 values) via C01d pass-2.
- Q7 `q7_tier` (final) via C01d pass-2.

## Wiring roadmap (chat 11+)

Priority order for closing out the 96 aspirational keys:

1. **Q7B breakpoint_evidence_audit flat-key extraction** (25 keys, one
   module, mechanical). Biggest single reduction.
2. **frequency.v2.schema `class_counts` flattening** (6 keys, schema-only
   change). Low risk.
3. **STEP03 Armitage flat-key writes** (3 keys, one function). Already
   computed, just needs the register call.
4. **cheat28 tandem-repeat flat-key extraction** (10 keys, one module).
5. **GHSL v5 Layer C registry wiring** (7 keys, needs new schema +
   pipeline integration). Higher-risk because it involves wiring a
   previously-standalone pipeline into the registry API.
6. **MODULE_3 theta ladder per-arrangement extraction** (12 keys
   between Q5 + Q6 — same source). Moderate; needs per-arrangement
   theta output decomposition.
7. **RepeatMasker + TE/GO annotation** (18 keys across Q3 + Q4).
   Requires running external annotation pipelines on the hybrid
   genome — work, but well-scoped.
8. **Q5 Dollo reconstruction** (3 keys). Requires synteny table
   wiring — lowest priority until multi-species comparison is a
   near-term goal.
9. **v9 back-compat placeholders** (3 keys). Only wire if a legacy
   reader needs them; `q6_family_linkage` covers the semantic case.

## Rules for maintaining this doc

1. When a writer gets wired, **move the key from `*_aspir` to the wired
   block** in `build_key_spec()` and remove its bullet here. The
   machine-readable source is code; this doc is a view of code.
2. New keys added to a Tier-2 schema's `keys_extracted` directive are
   automatically wired — they should be added to the `q<N>` list and NOT
   to `q<N>_aspir`.
3. To regenerate the summary table, run
   `Rscript /home/claude/walk_spec.R` (chat-10 helper) against
   `compute_candidate_status.R`. The script reads `build_key_spec()` and
   emits the current counts + aspirational per-Q listings.

## See also

- `4e_final_classification/compute_candidate_status.R::build_key_spec()`
  — the code source of truth.
- `4e_final_classification/characterize_candidate.R` — the per-Q
  convergence functions that consume these keys.
- `4e_final_classification/run_characterize.R` — chat-10 driver.
- `specs/INVERSION_REGISTRY_SPECIFICATION_v2.md` — the 352-key v2 target
  spec.
- `docs/PHASE4_ARCHITECTURE.md` — the registry-as-catalog design and
  per-question group-validation gates.
- Chat-9 `AUDIT_LOG_chat9_2026-04-17.md` FIX 42 — the source of the
  current spec list + aspirational markings.
