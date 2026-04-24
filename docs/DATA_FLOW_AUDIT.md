# Data Flow Audit — Chat A deliverable (2026-04-24)

**Scope:** map every structured block → schema → flat key → `build_key_spec()` registration → `characterize_candidate.R` consumer. Identify gaps in the wiring.

**Methodology:** static analysis of

- 37 JSON schemas in `registries/schemas/structured_block_schemas/`
- `build_key_spec()` in `inversion_modules/phase_9_classification/compute_candidate_status.R`
- Q-axis dispatchers in `inversion_modules/phase_9_classification/characterize_candidate.R`
- the two env-gated appenders (`_qc_shelf_reader.R`, `_axis5_final_label.R`)

Template substitution (`{side}` → `left`/`right`) is expanded. Aspirational keys declared in each `q<N>_aspir` list are distinguished from real phantoms. Schemas are keyed by filename because 7 schemas lack `block_type` and one `block_type` is shared between two schemas.

---

## Headline numbers

- **37 schemas** declaring **276 flat keys**
- **`build_key_spec()` declares 423 keys** across 9 axes (Q1–Q7, Q7B, Q_QC_SHELF), including **96 aspirational**
- **90 keys** actively read by `characterize_candidate.R`
- **149 orphan keys** (schemas produce, spec doesn't list) — **leaks**
- **208 real phantom keys** (spec declares non-aspirationally, nothing produces) — **missing writers**
- **88 aspirational phantoms** (spec declares *as aspirational*, tracked known-gap)
- **42 of 90 consumer reads (46%)** target real phantom keys — `characterize_candidate.R` is reading keys nothing writes

---

## 1. Structural schema issues

### 1.1 Schemas missing `block_type` field

**7 schemas** lack a `block_type` declaration. `registry_loader.R::load_schema(block_type)` looks up schemas by `block_type`, so a missing field means the writer can't find its schema at all (extraction silently produces 0 keys for that block).

- `boundary.schema.json` (question: Q3)
- `boundary_refined.schema.json` (question: Q3)
- `distance_concordance.schema.json` (question: None)
- `gene_conversion_tracts.schema.json` (question: None)
- `local_structure_segments.schema.json` (question: None)
- `regime_sample_dag.schema.json` (question: None)
- `regime_segments.schema.json` (question: None)

**Fix:** add `"block_type": "<stem>"` to each. The stem should match the schema filename without `.schema.json` extension, which is what `load_schema()` expects.

### 1.2 Duplicate `block_type` across schemas

Two schemas share the same `block_type`, which means `load_schema(bt)` is ambiguous:

- `block_type="frequency"` → ['frequency.schema.json', 'frequency.v2.schema.json']

**Fix:** either give the v2 schema its own block_type (e.g., `frequency_v2`) and update writers accordingly, or delete/archive the obsolete one. The current setup means whichever writer calls `write_block('frequency', ...)` picks up whichever schema `load_schema()` returns, which is filesystem-ordering dependent (currently returns `frequency.schema.json` because it sorts first).

## 2. Phantom keys (spec declares, nothing produces)

**208 non-aspirational phantom keys** — `build_key_spec()` lists them as real (non-aspirational) but no schema's `keys_extracted` produces them. These are the most concerning category: consumers read them expecting values, but the writers never fire.

### 2.1 Phantom distribution by axis

| Axis | Total spec keys | Real phantoms | Aspirational phantoms | % phantom |
|---|---|---|---|---|
| Q1 | 49 | 29 | 0 | 59% |
| Q2 | 91 | 37 | 0 | 40% |
| Q3 | 78 | 48 | 9 | 61% |
| Q4 | 47 | 27 | 18 | 57% |
| Q5 | 39 | 26 | 11 | 66% |
| Q6 | 31 | 1 | 16 | 3% |
| Q7 | 42 | 19 | 10 | 45% |
| Q7B | 33 | 8 | 24 | 24% |
| Q_QC_SHELF | 13 | 13 | 0 | 100% |

### 2.2 Phantoms that are actively consumed by `characterize_candidate.R`

These are the highest-priority fixes — `characterize_candidate.R` reads these keys and will get `NA` because nothing writes them:

#### `characterize_q1` — 4/8 consumer keys are real phantoms

- `q1_family_likeness_mean`
- `q1_inv_likeness_mean`
- `q1_nn_survives_40`
- `q1_squareness`

#### `characterize_q2` — 10/28 consumer keys are real phantoms

- `q2_class_entropy_mean`
- `q2_class_stability_pct`
- `q2_indel_class_concordance`
- `q2_n_double_crossover`
- `q2_n_recombinant`
- `q2_phase_block_ratio`
- `q2_recomb_center_n`
- `q2_recomb_left_n`
- `q2_recomb_right_n`
- `q2_switching_kw_pval`

#### `characterize_q3` — 1/9 consumer keys are real phantoms

- `q3_dropout_rate`

#### `characterize_q4` — 6/9 consumer keys are real phantoms

- `q4_biser2_concordance`
- `q4_has_inverted_sd`
- `q4_junction_type_left`
- `q4_mh_length_left`
- `q4_mh_length_right`
- `q4_sd_identity_pct`

#### `characterize_q5` — 9/11 consumer keys are real phantoms

- `q5_diversity_shape`
- `q5_dxy_ratio`
- `q5_fst_b1b3`
- `q5_gds_between`
- `q5_gds_fst_spearman_p`
- `q5_gds_gap_percentile`
- `q5_gds_het_pattern`
- `q5_gds_within_inv`
- `q5_gds_within_std`

#### `characterize_q6` — 1/6 consumer keys are real phantoms

- `q3_carrier_concordance`

#### `characterize_q7` — 4/12 consumer keys are real phantoms

- `q7_layer_a_detected`
- `q7_layer_b_detected`
- `q7_layer_d_tested`
- `q7b_observed_dropout_pct`

#### `characterize_q_qc_shelf` — 7/7 consumer keys are real phantoms

- `q_qc_shelf_coverage_ratio`
- `q_qc_shelf_flag`
- `q_qc_shelf_fst_enrichment_fold`
- `q_qc_shelf_hovere_het_inside`
- `q_qc_shelf_ran`
- `q_qc_shelf_snp_ratio`
- `q_qc_shelf_uncertain_ratio`

### 2.3 Full phantom list by axis (non-consumed)

Phantoms not currently read by `characterize_candidate.R` — spec declares them but nothing writes them AND nothing reads them. Likely candidates for removal from spec if they're no longer intended, or for wiring if they're planned features that were never marked aspirational.

**Q1** (25 non-consumed phantoms)

- `q1_ancestry_composite_flag`
- `q1_band_het_qgroup_fracs`
- `q1_band_inv_qgroup_fracs`
- `q1_band_std_qgroup_fracs`
- `q1_block_coherence`
- `q1_block_compactness`
- `q1_consensus_confidence`
- `q1_dip_pval_median`
- `q1_dominant_qgroup_inv`
- `q1_dominant_qgroup_std`
- `q1_dosage_het_cv_mean`
- `q1_flat_inv_score`
- `q1_fragmentation_score`
- `q1_inv_likeness_max`
- `q1_inv_likeness_sd`
- `q1_n_matrix_variants`
- `q1_n_windows`
- `q1_nn_birth`
- `q1_nn_survives_80`
- `q1_pve1_excess_mean`
- `q1_q_group_overlap`
- `q1_s_dip_mean`
- `q1_s_het_mean`
- `q1_s_pve_mean`
- `q1_spiky_inv_score`

**Q2** (27 non-consumed phantoms)

- `q2_accum_saturated_inv`
- `q2_accum_saturated_std`
- `q2_chao1_hom_inv`
- `q2_chao1_hom_std`
- `q2_decomp_quality_flags`
- `q2_diversity_gradient`
- `q2_diversity_gradient_r2`
- `q2_indel_ashman_D`
- `q2_mean_switch_rate_het`
- `q2_mean_switch_rate_inv`
- `q2_mean_switch_rate_rec`
- `q2_mean_switch_rate_std`
- `q2_n_founder_haplotypes_inv`
- `q2_n_founder_haplotypes_std`
- `q2_n_gene_conversion`
- `q2_n_indel_classes`
- `q2_n_windows_bimodal`
- `q2_n_windows_trimodal`
- `q2_n_windows_unimodal`
- `q2_phase_block_n50_flank`
- `q2_phase_block_n50_inside`
- `q2_phase_switch_rate_flank`
- `q2_phase_switch_rate_inside`
- `q2_profile_cor_observed`
- `q2_profile_cor_perm_p`
- `q2_recomb_family_spread`
- `q2_recomb_sample_ids`

**Q3** (46 non-consumed phantoms)

- `q3_concordance_pval_left`
- `q3_concordance_pval_right`
- `q3_left_cheats_list`
- `q3_left_concordance_kb`
- `q3_left_delta12_bp`
- `q3_left_fst_`
- `q3_left_fst_0kb`
- `q3_left_fst_bp`
- `q3_left_fst_m100kb`
- `q3_left_fst_m200kb`
- `q3_left_fst_m50kb`
- `q3_left_fst_p100kb`
- `q3_left_fst_p200kb`
- `q3_left_fst_p50kb`
- `q3_left_fst_step`
- `q3_left_hobs_`
- `q3_left_hobs_0kb`
- `q3_left_hobs_m100kb`
- `q3_left_hobs_p100kb`
- `q3_left_hobs_step`
- `q3_left_sv_bp`
- `q3_n_carriers_pca`
- `q3_n_carriers_sv`
- `q3_n_pca_carrier_sv_ref`
- `q3_n_rescued_by_prior`
- `q3_n_sv_carrier_pca_ref`
- `q3_population_prior_applied`
- `q3_right_cheats_list`
- `q3_right_concordance_kb`
- `q3_right_delta12_bp`
- `q3_right_fst_`
- `q3_right_fst_0kb`
- `q3_right_fst_bp`
- `q3_right_fst_m100kb`
- `q3_right_fst_m200kb`
- `q3_right_fst_m50kb`
- `q3_right_fst_p100kb`
- `q3_right_fst_p200kb`
- `q3_right_fst_p50kb`
- `q3_right_fst_step`
- `q3_right_hobs_`
- `q3_right_hobs_0kb`
- `q3_right_hobs_m100kb`
- `q3_right_hobs_p100kb`
- `q3_right_hobs_step`
- `q3_right_sv_bp`

**Q4** (21 non-consumed phantoms)

- `q4_decision_tree_path`
- `q4_gene_density_genome`
- `q4_gene_density_inside`
- `q4_gene_names_at_bp`
- `q4_junction_seq_left_50bp`
- `q4_junction_seq_right_50bp`
- `q4_junction_type_right`
- `q4_mechanism_support`
- `q4_mh_sequence`
- `q4_n_genes_inside`
- `q4_n_genes_spanning_bp`
- `q4_sd_alignment_length`
- `q4_sd_divergence_pct`
- `q4_sd_left_end`
- `q4_sd_left_start`
- `q4_sd_length`
- `q4_sd_n_mismatches`
- `q4_sd_orientation`
- `q4_sd_right_end`
- `q4_sd_right_start`
- `q4_tr_enrichment`

**Q5** (17 non-consumed phantoms)

- `q5_bimodality_dip_p`
- `q5_da_net_divergence`
- `q5_diversity_r2`
- `q5_dxy_std_inv`
- `q5_dxy_std_inv_flanking`
- `q5_gds_fst_spearman_rho`
- `q5_n_d_d_pairs`
- `q5_n_d_i_pairs`
- `q5_n_i_i_pairs`
- `q5_origin_gds_dip_p`
- `q5_origin_mechanism`
- `q5_theta_pi_flanking`
- `q5_theta_pi_het`
- `q5_theta_pi_inside`
- `q5_theta_pi_inv`
- `q5_theta_pi_std`
- `q5_theta_ratio`

**Q6** (1 non-consumed phantoms)

- `q6_validation_promotion_cap`

**Q7** (16 non-consumed phantoms)

- `q7_independence_class`
- `q7_layer_a_beta_pval`
- `q7_layer_a_beta_qval`
- `q7_layer_a_core_family`
- `q7_layer_a_inv_likeness`
- `q7_layer_a_pa_pattern`
- `q7_layer_b_bnd_triang`
- `q7_layer_b_ciend`
- `q7_layer_b_cipos`
- `q7_layer_b_delly`
- `q7_layer_b_manta`
- `q7_layer_b_n_carriers`
- `q7_layer_b_pe_support`
- `q7_layer_b_sr_support`
- `q7_n_layers_passed`
- `q7_n_layers_tested`

**Q7B** (7 non-consumed phantoms)

- `q7b_dropout_suspected_fraction`
- `q7b_n_sv_sites_audited`
- `q7b_pca_carrier_sv_support_pct`
- `q7b_pca_carriers_no_sv`
- `q7b_pca_carriers_strong_sv`
- `q7b_pca_carriers_weak_sv`
- `q7b_true_absence_fraction`

**Q_QC_SHELF** (6 non-consumed phantoms)

- `q_qc_shelf_coverage_cv_ratio`
- `q_qc_shelf_fst_hom1_hom2_inside`
- `q_qc_shelf_fst_hom1_hom2_outside`
- `q_qc_shelf_hovere_hom1_inside`
- `q_qc_shelf_hovere_hom2_inside`
- `q_qc_shelf_z_flatness`

## 3. Orphan keys (schema produces, spec doesn't list)

**149 orphan keys** — produced by schemas (so they land in keys.tsv) but never registered in `build_key_spec()`. These keys exist on disk per-candidate but don't contribute to any completion %, aren't exposed to the consumer layer cleanly, and aren't documented as part of the public key surface.

Could be either: (a) real features the writer added without updating spec, (b) debug/diagnostic keys that shouldn't be registered, or (c) keys that used to be in spec and got dropped without removing the writer.

### 3.1 Orphans by schema file

**`age_evidence.schema.json`** — block_type=`age_evidence`, question=Q5 — 9 orphans

- `q5_age_class`
- `q5_age_confidence`
- `q5_dxy_inside`
- `q5_fst_decay_rate`
- `q5_fst_inside`
- `q5_gds_dip_p`
- `q5_gds_percentile_global`
- `q5_theta_pi_hom_inv`
- `q5_theta_pi_hom_ref`

**`ancestral_fragments_summary.schema.json`** — block_type=`ancestral_fragments_summary`, question=Q3 — 6 orphans

- `q3_frag_left_bp_mode`
- `q3_frag_left_ci_width_kb`
- `q3_frag_n_carriers`
- `q3_frag_right_bp_mode`
- `q3_frag_right_ci_width_kb`
- `q3_frag_status`

**`band_composition.schema.json`** — block_type=`band_composition`, question=Q1 — 2 orphans

- `q1_ancestry_div_hom_ref_vs_hom_inv`
- `q1_k_used`

**`block_detect.schema.json`** — block_type=`block_detect`, question=Q1, Q3 — 6 orphans

- `q1_block_confidence`
- `q1_contrast_ratio`
- `q3_inner_left_bp`
- `q3_inner_right_bp`
- `q3_outer_left_bp`
- `q3_outer_right_bp`

**`bnd_sided_support.schema.json`** — block_type=`bnd_sided_support`, question=Q7 — 5 orphans

- `q7b_bnd_left_ct_count`
- `q7b_bnd_left_support`
- `q7b_bnd_right_ct_count`
- `q7b_bnd_right_support`
- `q7b_bnd_sided_support_class`

**`boundary_refined.schema.json`** — block_type=`(missing)`, question=Q3 — 12 orphans

- `q3_refined_left_bp`
- `q3_refined_left_ci_width_kb`
- `q3_refined_left_n_methods`
- `q3_refined_left_n_methods_agreeing`
- `q3_refined_left_primary_source`
- `q3_refined_left_shift_kb`
- `q3_refined_right_bp`
- `q3_refined_right_ci_width_kb`
- `q3_refined_right_n_methods`
- `q3_refined_right_n_methods_agreeing`
- `q3_refined_right_primary_source`
- `q3_refined_right_shift_kb`

**`breakpoints_per_method.schema.json`** — block_type=`breakpoints_per_method`, question=Q3 — 1 orphans

- `q3_n_methods_total`

**`burden.schema.json`** — block_type=`burden`, question=Q6 — 6 orphans

- `q6_burden_beta`
- `q6_burden_ci_lower`
- `q6_burden_ci_upper`
- `q6_burden_p`
- `q6_n_del_sites`
- `q6_selection_signal`

**`carrier_reconciliation.schema.json`** — block_type=`carrier_reconciliation`, question=Q3, Q6, Q7 — 3 orphans

- `q7_n_disagreeing_samples`
- `q7_n_reconciliation_sources`
- `q7_reconciliation_agreement`

**`dosage_blocks.schema.json`** — block_type=`dosage_blocks`, question=Q3 — 6 orphans

- `q3_dosage_core_n_markers`
- `q3_dosage_ext_left_bp`
- `q3_dosage_ext_right_bp`
- `q3_dosage_shift_left_kb`
- `q3_dosage_shift_right_kb`
- `q3_dosage_status`

**`encoding_robustness.schema.json`** — block_type=`encoding_robustness`, question=Q2 — 3 orphans

- `q2_encoding_mean_ari`
- `q2_encoding_min_ari`
- `q2_encoding_verdict`

**`existence_layer_a.schema.json`** — block_type=`existence_layer_a`, question=Q1, Q7 — 2 orphans

- `q7_cheat25_status`
- `q7_fdr_q_value`

**`existence_layer_b.schema.json`** — block_type=`existence_layer_b`, question=Q7 — 6 orphans

- `q7_bp1_concordance_kb`
- `q7_bp2_concordance_kb`
- `q7_n_delly_calls`
- `q7_n_manta_calls`
- `q7_n_sv_calls`
- `q7_sv_concordance_class`

**`existence_layer_b_bnd_rescue.schema.json`** — block_type=`existence_layer_b_bnd_rescue`, question=Q7 — 10 orphans

- `q7b_bnd_left_pe`
- `q7b_bnd_left_sr`
- `q7b_bnd_match_type`
- `q7b_bnd_matched_inv_id`
- `q7b_bnd_pair_bp1`
- `q7b_bnd_pair_bp2`
- `q7b_bnd_pair_size_bp`
- `q7b_bnd_rescue_source`
- `q7b_bnd_right_pe`
- `q7b_bnd_right_sr`

**`existence_layer_c.schema.json`** — block_type=`existence_layer_c`, question=Q7 — 7 orphans

- `q7_ghsl_concordance_class`
- `q7_ghsl_contrast_p`
- `q7_ghsl_contrast_score`
- `q7_ghsl_vs_pca_rand`
- `q7_n_clair3_snps`
- `q7_triangle_insulation_left`
- `q7_triangle_insulation_right`

**`existence_layer_d.schema.json`** — block_type=`existence_layer_d`, question=Q7 — 4 orphans

- `q7_layer_d_fisher_ci_lower`
- `q7_layer_d_fisher_ci_upper`
- `q7_layer_d_n_inv_total`
- `q7_layer_d_n_inv_with_support`

**`flank_coherence.schema.json`** — block_type=`flank_coherence`, question=Q3 — 3 orphans

- `q3_flank_coherence_class`
- `q3_left_flank_coherence_score`
- `q3_right_flank_coherence_score`

**`fragment_distribution.schema.json`** — block_type=`fragment_distribution`, question=Q2+Q3 — 10 orphans

- `q2_fragment_interior_recomb_fraction`
- `q2_interior_class`
- `q2_n_fragment_carriers`
- `q3_breakpoint_precision_class`
- `q3_final_left_bp`
- `q3_final_right_bp`
- `q3_left_ci_width_kb`
- `q3_n_methods_agreeing_left`
- `q3_n_methods_agreeing_right`
- `q3_right_ci_width_kb`

**`gene_conversion_tracts.schema.json`** — block_type=`(missing)`, question=None — 3 orphans

- `q2_gc_n_snps_pass_qc_diagnostic`
- `q2_gc_total_samples_with_gc`
- `q2_gc_total_tracts`

**`hypothesis_verdict.schema.json`** — block_type=`hypothesis_verdict`, question=Q7 — 3 orphans

- `q7_t1_ratio`
- `q7_t2_eff_k`
- `q7_t3_retention`

**`inv_detect_blocks.schema.json`** — block_type=`inv_detect_blocks`, question=Q1, Q3 — 6 orphans

- `q1_consensus_strength`
- `q1_n_scales_detected`
- `q1_staircase_rank`
- `q1_staircase_score`
- `q3_regime_change_left`
- `q3_regime_change_right`

**`mechanism.schema.json`** — block_type=`mechanism`, question=Q4 — 8 orphans

- `q4_junction_orientation`
- `q4_mechanism_class`
- `q4_microhomology_length_bp`
- `q4_sd_concordance`
- `q4_sd_identity_left`
- `q4_sd_identity_right`
- `q4_sd_length_left_bp`
- `q4_sd_length_right_bp`

**`mechanism_assembled.schema.json`** — block_type=`mechanism_assembled`, question=Q4 — 6 orphans

- `q4b_asm_consensus_available`
- `q4b_asm_homlen`
- `q4b_asm_junction_class`
- `q4b_asm_precise_record_available`
- `q4b_asm_source`
- `q4b_asm_vs_ref_concordance`

**`morphology.schema.json`** — block_type=`morphology`, question=Q1 — 8 orphans

- `q1_architecture`
- `q1_aspect_ratio`
- `q1_interior_homogeneity`
- `q1_patchiness_score`
- `q1_sample_grouping`
- `q1_size_class`
- `q1_spatial_consistency`
- `q1_stripe_count`

**`peel_diagnostic.schema.json`** — block_type=`peel_diagnostic`, question=Q7 — 4 orphans

- `q7_peel_combined`
- `q7_peel_l1b_score`
- `q7_peel_l2_score`
- `q7_peel_survives`

**`synteny_dollo.schema.json`** — block_type=`synteny_dollo`, question=Q5 — 5 orphans

- `q5_conservation_class`
- `q5_dollo_state`
- `q5_outgroup_state`
- `q5_present_c_gariepinus`
- `q5_present_c_macrocephalus`

**`synteny_v6.schema.json`** — block_type=`synteny_v6`, question=Q5 — 6 orphans

- `q5_bs_event_overlap`
- `q5_bs_event_type`
- `q5_conservation_class`
- `q5_dollo_vs_tree_concordance`
- `q5_tree_polarization_confidence`
- `q5_tree_polarization_direction`

## 4. Per-schema data flow map

One row per schema. `spec coverage` = fraction of the schema's declared keys that `build_key_spec()` recognizes. Where the ratio is below 100%, the orphan keys listed in §3 identify what's missing.

| Schema | block_type | question | keys declared | spec coverage | source_script |
|---|---|---|---|---|---|
| `age_evidence.schema.json` | `age_evidence` | Q5 | 11 | 2/11 | cheat30_gds_by_genotype + unified_ancestry/region_popstats |
| `ancestral_fragments_summary.schema.json` | `ancestral_fragments_summary` | Q3 | 6 | 0/6 | 02_ancestral_fragments.R |
| `band_composition.schema.json` | `band_composition` | Q1 | 2 | 0/2 | C01i + unified_ancestry/instant_q |
| `block_detect.schema.json` | `block_detect` | Q1, Q3 | 6 | 0/6 | PHASE_01C block_detect |
| `bnd_sided_support.schema.json` | `bnd_sided_support` | Q7 | 5 | 0/5 | bnd_sided_support.py |
| `boundary.schema.json` | `**MISSING**` | Q3 | 18 | 18/18 | STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_regist |
| `boundary_refined.schema.json` | `**MISSING**` | Q3 | 12 | 0/12 | 03_consensus_merge.R |
| `breakpoints_per_method.schema.json` | `breakpoints_per_method` | Q3 | 1 | 0/1 | 03_consensus_merge.R |
| `burden.schema.json` | `burden` | Q6 | 6 | 0/6 | burden_regression + STEP_C01f_c_burden_regression.R |
| `carrier_reconciliation.schema.json` | `carrier_reconciliation` | Q3, Q6, Q7 | 3 | 0/3 | audit_cores_vs_sv (new) + carrier_audit.R |
| `distance_concordance.schema.json` | `**MISSING**` | — | 3 | 3/3 | — |
| `dosage_blocks.schema.json` | `dosage_blocks` | Q3 | 6 | 0/6 | 01_dosage_signal.R |
| `encoding_robustness.schema.json` | `encoding_robustness` | Q2 | 3 | 0/3 | 05_four_encoding_diagnostic.R |
| `existence_layer_a.schema.json` | `existence_layer_a` | Q1, Q7 | 23 | 21/23 | C01a + C01b + C01d |
| `existence_layer_b.schema.json` | `existence_layer_b` | Q7 | 6 | 0/6 | C00 build_flashlight + MODULE_4H_ALL_Manta + MODULE_4_DELLY_ |
| `existence_layer_b_bnd_rescue.schema.json` | `existence_layer_b_bnd_rescue` | Q7 | 11 | 1/11 | phase_3_refine/STEP_B06_bnd_rescue.py |
| `existence_layer_c.schema.json` | `existence_layer_c` | Q7 | 7 | 0/7 | C04 GHSL v5 (within-sample haplotype divergence) + STEP_C01c |
| `existence_layer_d.schema.json` | `existence_layer_d` | Q7 | 8 | 4/8 | STEP03_statistical_tests.py |
| `flank_coherence.schema.json` | `flank_coherence` | Q3 | 3 | 0/3 | cross_species_bridge.py |
| `fragment_distribution.schema.json` | `fragment_distribution` | Q2+Q3 | 10 | 0/10 | bp_pipeline_bridge.py (reads breakpoint_pipeline/01+02+03 ou |
| `frequency.schema.json` | `frequency` | Q6 | 8 | 8/8 | C01i + C01f |
| `frequency.v2.schema.json` | `frequency` | Q6 | 10 | 10/10 | C01i_d_seal + C01f |
| `gene_conversion_tracts.schema.json` | `**MISSING**` | — | 3 | 0/3 | — |
| `hypothesis_verdict.schema.json` | `hypothesis_verdict` | Q7 | 14 | 11/14 | STEP_C01f_hypothesis_tests_wired_6_9_12_23_v934_registry.R |
| `internal_ancestry_composition.schema.json` | `internal_ancestry_composition` | Q1, Q2 | 12 | 12/12 | STEP_C01i_c_nested_composition.py |
| `internal_dynamics.schema.json` | `internal_dynamics` | Q1, Q2 | 19 | 19/19 | STEP_C01i_decompose.R |
| `inv_detect_blocks.schema.json` | `inv_detect_blocks` | Q1, Q3 | 6 | 0/6 | inv_detect_v9.3 + 12_bridge_to_codebase.R |
| `local_structure_segments.schema.json` | `**MISSING**` | — | 6 | 6/6 | — |
| `mechanism.schema.json` | `mechanism` | Q4 | 10 | 2/10 | cheat14 + cheat21 + cheat27 + cheat28 + cheat29 |
| `mechanism_assembled.schema.json` | `mechanism_assembled` | Q4 | 6 | 0/6 | cheat29b_assembled_junction.py |
| `morphology.schema.json` | `morphology` | Q1 | 10 | 2/10 | C01a precompute + C01d scoring |
| `peel_diagnostic.schema.json` | `peel_diagnostic` | Q7 | 4 | 0/4 | C01n peel / inv_detect Phase 9 |
| `recombinant_map.schema.json` | `recombinant_map` | Q2, Q4 | 10 | 10/10 | STEP_C01i_b_multi_recomb.R |
| `regime_sample_dag.schema.json` | `**MISSING**` | — | 5 | 5/5 | — |
| `regime_segments.schema.json` | `**MISSING**` | — | 4 | 4/4 | — |
| `synteny_dollo.schema.json` | `synteny_dollo` | Q5 | 5 | 0/5 | synteny_inversion_detection.sh + dollo_parsimony_inversions. |
| `synteny_v6.schema.json` | `synteny_v6` | Q5 | 6 | 0/6 | cross_species_bridge_v6.py |

## 5. Env-gated appenders

Two scripts append keys into the in-memory `keys` list after the main spec is built:

- `_qc_shelf_reader.R` — gated by `QC_SHELF_EVIDENCE_DIR` (pass 13 public API)
- `_axis5_final_label.R` — gated by `V7_FINAL_DIR` (pass 8 axis-5 wiring)

### 5.1 `_qc_shelf_reader.R`

Keys touched (detected by static analysis): 2

- `q_qc_shelf_flag` — **in spec**
- `q_qc_shelf_ran` — **in spec**

### 5.2 `_axis5_final_label.R`

_No keys detected by static analysis — may use dynamic key names or the detector missed them. Manual inspection recommended._

## 6. Recommended chat B fix priority

Ordered by impact, not by difficulty:

**P1 — schemas missing `block_type` (§1.1).** These are silent failures: `write_block()` can't locate the schema, extraction produces zero keys for the block, and there is no warning. 7 schemas affected. Pure mechanical fix (one JSON field per schema). Highest-leverage item in the audit.

**P2 — duplicate `block_type` on frequency schemas (§1.2).** Ambiguous schema lookup, resolution depends on filesystem ordering. Fix by renaming `frequency.v2.schema.json`'s `block_type` to `frequency_v2` or archiving whichever is obsolete. Quentin's call.

**P3 — consumer-visible phantoms (§2.2).** `characterize_candidate.R` actively reads these keys; they always arrive as NA. Fixes here are per-axis:

- Check whether each phantom is genuinely unmet (needs a writer) or should be marked aspirational (if the writer is scheduled but not yet built).
- For Q3's `q3_left_hardness`, `q3_left_verdict`, etc. — `boundary.schema.json` uses `{side}` templating; verify `registry_loader.R` actually performs the substitution (the schema's note says it does). If expansion is runtime-only and produces keys not visible to static analysis, the audit false-positives here — add a code inspection step.
- For Q5 and Q_QC_SHELF (9/11 and 7/7 phantoms respectively), most likely the writers exist but key names drifted. Audit the actual writer outputs.

**P4 — non-consumed phantoms (§2.3).** Dead weight in the spec: declared but neither produced nor consumed. Best fix is usually removal from spec rather than building a writer, unless the biological intent is still live.

**P5 — orphan keys (§3).** Schemas produce keys the spec doesn't register. Lowest priority because the keys do land in keys.tsv and are reachable — just not formally part of the API. Either register them (add to relevant `q<N>` vector) or decide they're diagnostic and leave them.

---

## 7. What this audit does NOT cover

Honest list of static-analysis limits so chat B doesn't over-trust the numbers:

- **Runtime key expansion.** The `{side}` template is expanded here; any other runtime substitutions (none detected but possible) would produce false-positive phantoms.
- **Writers bypassing `write_block()`.** If a script writes directly to keys.tsv without going through `registry_loader.R::write_block()`, my analysis won't see it. Worth grepping `keys.tsv` writers.
- **Consumer keys via dynamic lookup.** `characterize_candidate.R` uses static `keys[["literal"]]` patterns that I parse; if any dispatcher builds key names with `paste0()` at runtime, they're not captured.
- **Env-gated appender parsing.** The two appenders are scanned by simple regex. Recommend manual spot-check in chat B.
- **Aspirational-vs-real classification.** Based on presence in `q<N>_aspir` lists. If a key is meant to be aspirational but wasn't added to the list, it shows up as a real phantom here.

---

_Generated 2026-04-24 by `audit/generate_audit.py` against the post-pass-15 tree. Re-run after any schema or `build_key_spec()` change to refresh._