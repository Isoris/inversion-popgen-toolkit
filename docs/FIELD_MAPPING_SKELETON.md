# Field-mapping skeleton for HELPER_MISSING schemas

For each of the 6 schemas that are waiting on helpers in `utils/registry_key_helpers.R`, this document lists: the schema's declared `keys_extracted` contract, the writer script(s) that should supply the data, and the actual column/variable names in each writer.

Match status legend:

- `MATCH` — schema's `from` path appears as-is in the writer
- `DRIFT` — similar name exists (e.g., `sharpness` vs `cheat4_sharpness`)
- `MISSING` — no matching or similar name in the writer
- `NA_TODO` — writer initializes the column to NA and never populates it (feature archived or not wired)

**This document records mismatches for Quentin to resolve. Claude did not pick any mappings.**

## `boundary`

- **Helper function:** `register_C01g_boundary(bd, cid, side, outdir)`
- **Writer(s):** `inversion_modules/phase_4_catalog/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R`
- **Notes:** Called per-side (left/right) inside a loop in C01g around L1430. Helper needs to build a block_data list and call write_block_safe(reg, cid, paste0('boundary_', side), data, ...).

### Writer: `inversion_modules/phase_4_catalog/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q3_{side}_bp` | `boundary_bp` | **MATCH** | L220,237,267 |
| `q3_{side}_hardness` | `hardness` | **MATCH** | L520,559 |
| `q3_{side}_sharpness` | `sharpness` | **MATCH** | L756 |
| `q3_{side}_type` | `boundary_type` | **MATCH** | L224,238,268 |
| `q3_{side}_verdict` | `boundary_verdict` | **MATCH** | L1354,1385,1386 |
| `q3_{side}_n_cheats` | `n_cheats_supporting` | **MATCH** | L1333,1356,1357 |
| `q3_{side}_clip_count` | `clip_count` | **DRIFT** | `cheat11_clip_bimodal`, `cheat11_clip_enrichment`, `cheat11_clip_score`, `clip_frac`, `n_clip` |
| `q3_{side}_depth_anomaly` | `depth_anomaly` | **DRIFT** | `cheat10_depth_dip`, `cheat10_depth_ratio`, `depth` |
| `q3_{side}_is_fossil` | `is_fossil` | **DRIFT** | `n_fossil` |

## `existence_layer_a`

- **Helper function:** `register_C01d_keys(cd, cid, outdir)`
- **Writer(s):** `inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
- **Notes:** Called per-candidate inside a loop in C01d around L1043. cd is one row of cand_dt (a data.table). Also references a store_C01d_results() that writes diagnostic files.

### Writer: `inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q1_composite_score` | `composite_score` | **MISSING** | — |
| `q1_dim_positive` | `dim_positive` | **MATCH** | L367,437,836 |
| `q1_d01_block_strength` | `d1_block_strength` | **MATCH** | L422 |
| `q1_d02_block_shape` | `d2_block_shape` | **MATCH** | L423 |
| `q1_d03_nn_persistence` | `d3_nn_persistence` | **MATCH** | L424 |
| `q1_d04_decay_flatness` | `d4_decay_flatness` | **MATCH** | L425 |
| `q1_d05_interior_quality` | `d5_interior_quality` | **MATCH** | L426 |
| `q1_d06_consensus` | `d6_consensus` | **MATCH** | L427 |
| `q1_d07_sv_breakpoint` | `d7_sv_breakpoint` | **MATCH** | L428,645 |
| `q1_d08_peel_or_hyp` | `d8_peel_or_hyp` | **MATCH** | L429,647 |
| `q1_d09_pca_clusters` | `d9_pca_clusters` | **MATCH** | L430 |
| `q1_d10_partition` | `d10_partition` | **MATCH** | L431 |
| `q1_d11_boundary_concordance` | `d11_boundary_concordance` | **MATCH** | L433,618,646 |
| `q1_d12_snake_concordance` | `d12_snake_concordance` | **MATCH** | L434,530,543 |
| `q1_shape_class` | `shape_class` | **MATCH** | L211,391,441 |
| `q1_landscape_category` | `landscape_category` | **MATCH** | L392,442 |
| `q1_n_children` | `n_children` | **MATCH** | L390 |
| `q1_span_kb` | `span_kb` | **DRIFT** | `span_mb` |
| `q7_tier` | `tier` | **MATCH** | L382,387,438 |
| `q7_composite_score` | `composite_score` | **MISSING** | — |
| `q7_dim_positive` | `dim_positive` | **MATCH** | L367,437,836 |
| `q7_cheat25_status` | `cheat25_status` | **MATCH** | L640,669 |
| `q7_fdr_q_value` | `fdr_q_value` | **MISSING** | — |

## `frequency`

- **Helper function:** `N/A (no helper declared — writer currently uses add_evidence directly)`
- **Writer(s):** `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R`, `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`
- **Notes:** C01i_d_seal writes 5 flat keys via add_evidence() (q6_family_linkage, q6_polymorphism_class, q6_group_validation, q6_validation_promotion_cap, q2_decomp_quality_flags). Full frequency block never written as write_block. Decision: keep add_evidence or promote to write_block + frequency schema? frequency.schema.json and frequency.v2.schema.json collide on block_type.

### Writer: `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q6_freq_inv` | `freq_inv` | **MATCH** | L369,429 |
| `q6_freq_class` | `freq_class` | **DRIFT** | `freq_inv` |
| `q6_n_total` | `n_total` | **MATCH** | L429 |
| `q6_hwe_chi2` | `hwe_chi2` | **MISSING** | — |
| `q6_hwe_p` | `hwe_p` | **MISSING** | — |
| `q6_hwe_verdict` | `hwe_verdict` | **MISSING** | — |
| `q6_genotype_balance` | `genotype_balance` | **MISSING** | — |
| `q6_family_linkage` | `family_linkage` | **MATCH** | L433 |

### Writer: `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q6_freq_inv` | `freq_inv` | **MISSING** | — |
| `q6_freq_class` | `freq_class` | **MISSING** | — |
| `q6_n_total` | `n_total` | **MISSING** | — |
| `q6_hwe_chi2` | `hwe_chi2` | **MISSING** | — |
| `q6_hwe_p` | `hwe_p` | **MISSING** | — |
| `q6_hwe_verdict` | `hwe_verdict` | **DRIFT** | `jackknife_verdict`, `t9_jackknife_verdict`, `verdict`, `verdict_role` |
| `q6_genotype_balance` | `genotype_balance` | **MISSING** | — |
| `q6_family_linkage` | `family_linkage` | **MATCH** | L513,519,521 |

## `frequency.v2`

- **Helper function:** `N/A`
- **Writer(s):** `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R`, `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`
- **Notes:** See frequency. Superset schema (adds q6_jackknife_max_delta, q6_polymorphism_class). Choose this or archive v1.

### Writer: `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q6_freq_inv` | `freq_inv` | **MATCH** | L369,429 |
| `q6_freq_class` | `freq_class` | **DRIFT** | `freq_inv` |
| `q6_n_total` | `n_total` | **MATCH** | L429 |
| `q6_hwe_chi2` | `hwe_chi2` | **MISSING** | — |
| `q6_hwe_p` | `hwe_p` | **MISSING** | — |
| `q6_hwe_verdict` | `hwe_verdict` | **MISSING** | — |
| `q6_genotype_balance` | `genotype_balance` | **MISSING** | — |
| `q6_family_linkage` | `family_linkage` | **MATCH** | L433 |
| `q6_jackknife_max_delta` | `jackknife_max_delta` | **MISSING** | — |
| `q6_polymorphism_class` | `polymorphism_class` | **MATCH** | L434 |

### Writer: `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q6_freq_inv` | `freq_inv` | **MISSING** | — |
| `q6_freq_class` | `freq_class` | **MISSING** | — |
| `q6_n_total` | `n_total` | **MISSING** | — |
| `q6_hwe_chi2` | `hwe_chi2` | **MISSING** | — |
| `q6_hwe_p` | `hwe_p` | **MISSING** | — |
| `q6_hwe_verdict` | `hwe_verdict` | **DRIFT** | `jackknife_verdict`, `t9_jackknife_verdict`, `verdict`, `verdict_role` |
| `q6_genotype_balance` | `genotype_balance` | **MISSING** | — |
| `q6_family_linkage` | `family_linkage` | **MATCH** | L513,519,521 |
| `q6_jackknife_max_delta` | `jackknife_max_delta` | **DRIFT** | `delta`, `jackknife`, `jackknife_verdict`, `max_delta`, `step9_jackknife` |
| `q6_polymorphism_class` | `polymorphism_class` | **MISSING** | — |

## `hypothesis_verdict`

- **Helper function:** `register_C01f_keys(vd, cid, outdir)`
- **Writer(s):** `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`
- **Notes:** Called per-candidate inside a loop in C01f L2535. vd is the per-candidate verdict data frame.

### Writer: `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q7_verdict` | `verdict` | **MATCH** | L1398,1404,1405 |
| `q7_verdict_confidence` | `verdict_confidence` | **DRIFT** | `confidence`, `jackknife_verdict`, `t9_jackknife_verdict`, `verdict`, `verdict_role` |
| `q7_t1_ratio` | `t1_ratio` | **MATCH** | L2416 |
| `q7_t2_eff_k` | `t2_eff_k` | **MATCH** | L2417 |
| `q7_t3_retention` | `t3_retention` | **MATCH** | L2418 |
| `q7_t8_clair3_concordance` | `t8_concordance` | **MATCH** | L503,2381,2426 |
| `q7_t9_jackknife_status` | `t9_jackknife_status` | **DRIFT** | `jackknife`, `jackknife_verdict`, `step9_jackknife`, `t9_jackknife_verdict` |
| `q7_t9_max_delta` | `t9_max_delta` | **MATCH** | L505,2383,2433 |
| `q7_t10_theta_concordance` | `t10_theta_concordance` | **MATCH** | L506,2384,2439 |
| `q3_extended_suppression` | `t11_has_extended` | **MATCH** | L2445 |
| `q3_suppression_extent_kb` | `t11_extent_kb` | **MATCH** | L2446 |
| `q3_fst_decay_rate` | `t11_fst_decay_rate` | **DRIFT** | `fst_decay_rate`, `moderate`, `rate` |
| `q6_group_validation` | `group_validation_after` | **MATCH** | L2398 |
| `q6_family_linkage` | `t9_jackknife_status` | **DRIFT** | `jackknife`, `jackknife_verdict`, `step9_jackknife`, `t9_jackknife_verdict` |

## `morphology`

- **Helper function:** `(no helper declared yet — precomputed and consumed, never written as block)`
- **Writer(s):** `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R`, `inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
- **Notes:** Schema says produced by 'C01a precompute + C01d scoring' but neither writes a morphology block today. The values (shape_class, architecture, patchiness_score, etc.) flow as columns of cand_dt but don't go through write_block. Possibly a candidate for auto-promotion via register_C01d_keys.

### Writer: `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q1_span_kb` | `span_kb` | **DRIFT** | `chrom_span_mb`, `span` |
| `q1_aspect_ratio` | `aspect_ratio` | **DRIFT** | `family_fst_ratio`, `inv_eigen_ratio`, `test05_family_fst_ratio` |
| `q1_interior_homogeneity` | `interior_homogeneity` | **MISSING** | — |
| `q1_patchiness_score` | `patchiness_score` | **MISSING** | — |
| `q1_stripe_count` | `stripe_count` | **MISSING** | — |
| `q1_n_children` | `n_children` | **MISSING** | — |
| `q1_architecture` | `architecture` | **MISSING** | — |
| `q1_size_class` | `size_class` | **DRIFT** | `size` |
| `q1_sample_grouping` | `sample_grouping` | **DRIFT** | `sample_inv_states` |
| `q1_spatial_consistency` | `spatial_consistency` | **MISSING** | — |

### Writer: `inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q1_span_kb` | `span_kb` | **DRIFT** | `span_mb` |
| `q1_aspect_ratio` | `aspect_ratio` | **DRIFT** | `cheat5_family_fst_ratio`, `family_fst_ratio`, `far_near_ratio` |
| `q1_interior_homogeneity` | `interior_homogeneity` | **DRIFT** | `d5_interior_quality`, `homogeneity` |
| `q1_patchiness_score` | `patchiness_score` | **DRIFT** | `patchiness` |
| `q1_stripe_count` | `stripe_count` | **MATCH** | L239 |
| `q1_n_children` | `n_children` | **MATCH** | L390 |
| `q1_architecture` | `architecture` | **MISSING** | — |
| `q1_size_class` | `size_class` | **DRIFT** | `base_size`, `size` |
| `q1_sample_grouping` | `sample_grouping` | **MISSING** | — |
| `q1_spatial_consistency` | `spatial_consistency` | **MISSING** | — |
