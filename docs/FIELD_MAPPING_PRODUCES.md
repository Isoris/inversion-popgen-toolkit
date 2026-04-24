# Field-mapping analysis — PRODUCES_BUT_NOT_WIRED schemas

For the 6 schemas where a writer script exists and produces output files but never calls a registry function. Same MATCH/DRIFT/MISSING legend as `FIELD_MAPPING_SKELETON.md`.

**Difference from HELPER_MISSING:** these scripts don't need a helper file — they can call `write_block()` or `write_block_safe()` directly. The decision is where in the script to add the call and which fields to include.

## `band_composition`

- **Writer(s):** `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R`
- **Note:** C01i_decompose already writes an internal_dynamics block via write_block_safe. band_composition was declared separately but it looks like the data it wants is computed as part of decompose. Decision: either (a) promote band_composition keys into internal_dynamics via a schema update, or (b) add a second write_block_safe call for band_composition.

### Writer: `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q1_k_used` | `k_used` | **MATCH** | L523 |
| `q1_ancestry_div_hom_ref_vs_hom_inv` | `ancestry_divergence_hom_ref_vs_hom_inv` | **MISSING** | — |

## `block_detect`

- **Writer(s):** `inversion_modules/phase_2_discovery/2c_precomp/PHASE_01C_block_detect.R`
- **Note:** PHASE_01C is the block-detection stage — produces RDS / TSV diagnostic outputs but never registers keys. Expected to populate q1/q3 block-structure keys.

### Writer: `inversion_modules/phase_2_discovery/2c_precomp/PHASE_01C_block_detect.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q1_block_confidence` | `confidence_score` | **MISSING** | — |
| `q3_outer_left_bp` | `outer_boundary_left_bp` | **DRIFT** | `at_contig_boundary`, `boundary_type`, `inner_soft_boundary_candidate`, `n_inner_soft_boundary_candidate` |
| `q3_outer_right_bp` | `outer_boundary_right_bp` | **DRIFT** | `at_contig_boundary`, `boundary_type`, `inner_soft_boundary_candidate`, `n_inner_soft_boundary_candidate` |
| `q3_inner_left_bp` | `inner_boundary_left_bp` | **DRIFT** | `at_contig_boundary`, `boundary_type`, `inner`, `inner_ambiguous`, `inner_hard_assembly` |
| `q3_inner_right_bp` | `inner_boundary_right_bp` | **DRIFT** | `at_contig_boundary`, `boundary_type`, `inner`, `inner_ambiguous`, `inner_hard_assembly` |
| `q1_contrast_ratio` | `contrast_ratio` | **MISSING** | — |

## `existence_layer_b`

- **Writer(s):** `inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R`
- **Note:** Layer B = SV prior (Manta + DELLY INV). C00 builds the prior catalog as a TSV but doesn't register per-candidate layer_b keys. Downstream Layer B block is also populated by STEP_B06_bnd_rescue (confirmed OK in TRIAGE_LEDGER), so existence_layer_b may be partially wired already — this schema is the upstream 'signal at seed time' variant.

### Writer: `inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q7_n_sv_calls` | `n_sv_calls` | **DRIFT** | `calls`, `inv_calls`, `n_inv_calls` |
| `q7_n_delly_calls` | `n_delly_inv_calls` | **DRIFT** | `calls`, `delly`, `inv_calls`, `n_inv_calls`, `n_inv_delly` |
| `q7_n_manta_calls` | `n_manta_inv_calls` | **DRIFT** | `calls`, `inv_calls`, `manta`, `n_inv_calls`, `n_inv_manta` |
| `q7_sv_concordance_class` | `sv_concordance_class` | **DRIFT** | `concordance_frac` |
| `q7_bp1_concordance_kb` | `bp1_concordance_kb` | **DRIFT** | `concordance_frac` |
| `q7_bp2_concordance_kb` | `bp2_concordance_kb` | **DRIFT** | `concordance_frac` |

## `existence_layer_c`

- **Writer(s):** `inversion_modules/phase_2_discovery/2e_ghsl/STEP_C04_snake3_ghsl_v6.R`, `inversion_modules/phase_2_discovery/2e_ghsl/STEP_C04b_snake3_ghsl_classify.R`
- **Note:** Layer C = GHSL. Runs, produces RDS/TSV output per chromosome. Per-candidate wiring never added. The schema's source_script mentions STEP_C01c triangle insulation too — check if that's current (was there a v5 -> v6 transition that dropped it?).

### Writer: `inversion_modules/phase_2_discovery/2e_ghsl/STEP_C04_snake3_ghsl_v6.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q7_ghsl_contrast_score` | `ghsl_contrast_score` | **DRIFT** | `ghsl`, `ghsl_div` |
| `q7_ghsl_contrast_p` | `ghsl_contrast_p` | **DRIFT** | `ghsl`, `ghsl_div` |
| `q7_ghsl_vs_pca_rand` | `ghsl_vs_pca_rand_index` | **DRIFT** | `ghsl`, `ghsl_div`, `index` |
| `q7_ghsl_concordance_class` | `ghsl_concordance_class` | **DRIFT** | `ghsl`, `ghsl_div` |
| `q7_n_clair3_snps` | `n_clair3_snps_in_region` | **MISSING** | — |
| `q7_triangle_insulation_left` | `triangle_insulation_score` | **MISSING** | — |
| `q7_triangle_insulation_right` | `triangle_insulation_right` | **MISSING** | — |

### Writer: `inversion_modules/phase_2_discovery/2e_ghsl/STEP_C04b_snake3_ghsl_classify.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q7_ghsl_contrast_score` | `ghsl_contrast_score` | **DRIFT** | `contrast_mad`, `contrast_med`, `div_contrast`, `div_contrast_z`, `ghsl_` |
| `q7_ghsl_contrast_p` | `ghsl_contrast_p` | **DRIFT** | `contrast_mad`, `contrast_med`, `div_contrast`, `div_contrast_z`, `ghsl_` |
| `q7_ghsl_vs_pca_rand` | `ghsl_vs_pca_rand_index` | **DRIFT** | `ghsl_`, `ghsl_v6_score`, `ghsl_v6_status` |
| `q7_ghsl_concordance_class` | `ghsl_concordance_class` | **DRIFT** | `ghsl_`, `ghsl_v6_score`, `ghsl_v6_status` |
| `q7_n_clair3_snps` | `n_clair3_snps_in_region` | **MISSING** | — |
| `q7_triangle_insulation_left` | `triangle_insulation_score` | **MISSING** | — |
| `q7_triangle_insulation_right` | `triangle_insulation_right` | **MISSING** | — |

## `flank_coherence`

- **Writer(s):** `inversion_modules/phase_8_evidence_biology/q5_age_and_origin/cross_species_bridge_v6.py`
- **Note:** cross_species_bridge_v6.py already writes a synteny_v6 block (OK in TRIAGE_LEDGER). The flank_coherence schema is a different aspect — the old cross_species_bridge.py (v5?) may have produced this. Check if flank_coherence keys belong in synteny_v6 now.

### Writer: `inversion_modules/phase_8_evidence_biology/q5_age_and_origin/cross_species_bridge_v6.py`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q3_left_flank_coherence_score` | `q3_left_flank_coherence_score` | **DRIFT** | `flank_class`, `flank_cls`, `flank_idx` |
| `q3_right_flank_coherence_score` | `q3_right_flank_coherence_score` | **DRIFT** | `flank_class`, `flank_cls`, `flank_idx` |
| `q3_flank_coherence_class` | `q3_flank_coherence_class` | **DRIFT** | `flank_class`, `flank_cls`, `flank_idx` |
