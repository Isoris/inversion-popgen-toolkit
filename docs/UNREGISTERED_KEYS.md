# Unregistered keys ledger

Keys that writers produce (or would produce once wired) but `build_key_spec()` in `phase_9_classification/compute_candidate_status.R` doesn't list in any q-axis vector. These keys still land in `keys.tsv` via schema `keys_extracted`, but consumers of `characterize_candidate.R` can't read them via the fast-path `keys[["..."]]` lookup — they'd have to use `read_block()` for the full JSON.

**Two sections:**

1. **Wired schemas** — writer is currently running, data IS flowing to `keys.tsv`. Consumer reads will miss these 70 keys unless they're registered. Each key needs a decision: register in spec (fast-path) or leave as block-only (diagnostic).
2. **Not-yet-wired schemas** — writer isn't running (HELPER_MISSING, SCRIPT_IS_LIB, etc.). These 89 keys are moot until the writer is wired, but listing them helps with the downstream decision.

Plus 6 truly-unowned keys (figure paths + registry markers) that don't need spec registration.

**Totals:** wired-but-unregistered: 70, not-wired-unregistered: 89, orphan: 6

## Section 1: Wired schemas — data flowing, not registered

These keys are in `keys.tsv` right now. Each needs a per-key decision.

### `boundary_refined` — 12 unregistered

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

### `fragment_distribution` — 10 unregistered

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

### `existence_layer_b_bnd_rescue` — 10 unregistered

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

### `dosage_blocks` — 6 unregistered

- `q3_dosage_core_n_markers`
- `q3_dosage_ext_left_bp`
- `q3_dosage_ext_right_bp`
- `q3_dosage_shift_left_kb`
- `q3_dosage_shift_right_kb`
- `q3_dosage_status`

### `ancestral_fragments_summary` — 6 unregistered

- `q3_frag_left_bp_mode`
- `q3_frag_left_ci_width_kb`
- `q3_frag_n_carriers`
- `q3_frag_right_bp_mode`
- `q3_frag_right_ci_width_kb`
- `q3_frag_status`

### `mechanism_assembled` — 6 unregistered

- `q4b_asm_consensus_available`
- `q4b_asm_homlen`
- `q4b_asm_junction_class`
- `q4b_asm_precise_record_available`
- `q4b_asm_source`
- `q4b_asm_vs_ref_concordance`

### `synteny_v6` — 6 unregistered

- `q5_bs_event_overlap`
- `q5_bs_event_type`
- `q5_conservation_class`
- `q5_dollo_vs_tree_concordance`
- `q5_tree_polarization_confidence`
- `q5_tree_polarization_direction`

### `bnd_sided_support` — 5 unregistered

- `q7b_bnd_left_ct_count`
- `q7b_bnd_left_support`
- `q7b_bnd_right_ct_count`
- `q7b_bnd_right_support`
- `q7b_bnd_sided_support_class`

### `existence_layer_d` — 4 unregistered

- `q7_layer_d_fisher_ci_lower`
- `q7_layer_d_fisher_ci_upper`
- `q7_layer_d_n_inv_total`
- `q7_layer_d_n_inv_with_support`

### `encoding_robustness` — 3 unregistered

- `q2_encoding_mean_ari`
- `q2_encoding_min_ari`
- `q2_encoding_verdict`

### `breakpoints_per_method` — 1 unregistered

- `q3_n_methods_total`

### `synteny_dollo` — 1 unregistered

- `q5_conservation_class`

## Section 2: Not-yet-wired schemas — unregistered (moot until wired)

These writers aren't currently running. When they're wired (via helpers or runners), each of these keys will face the same register-or-block-only decision as Section 1.

### `boundary` — 9 unregistered

- `q3_right_bp`
- `q3_right_clip_count`
- `q3_right_depth_anomaly`
- `q3_right_hardness`
- `q3_right_is_fossil`
- `q3_right_n_cheats`
- `q3_right_sharpness`
- `q3_right_type`
- `q3_right_verdict`

### `age_evidence` — 9 unregistered

- `q5_age_class`
- `q5_age_confidence`
- `q5_dxy_inside`
- `q5_fst_decay_rate`
- `q5_fst_inside`
- `q5_gds_dip_p`
- `q5_gds_percentile_global`
- `q5_theta_pi_hom_inv`
- `q5_theta_pi_hom_ref`

### `morphology` — 8 unregistered

- `q1_architecture`
- `q1_aspect_ratio`
- `q1_interior_homogeneity`
- `q1_patchiness_score`
- `q1_sample_grouping`
- `q1_size_class`
- `q1_spatial_consistency`
- `q1_stripe_count`

### `mechanism` — 8 unregistered

- `q4_junction_orientation`
- `q4_mechanism_class`
- `q4_microhomology_length_bp`
- `q4_sd_concordance`
- `q4_sd_identity_left`
- `q4_sd_identity_right`
- `q4_sd_length_left_bp`
- `q4_sd_length_right_bp`

### `existence_layer_c` — 7 unregistered

- `q7_ghsl_concordance_class`
- `q7_ghsl_contrast_p`
- `q7_ghsl_contrast_score`
- `q7_ghsl_vs_pca_rand`
- `q7_n_clair3_snps`
- `q7_triangle_insulation_left`
- `q7_triangle_insulation_right`

### `block_detect` — 6 unregistered

- `q1_block_confidence`
- `q1_contrast_ratio`
- `q3_inner_left_bp`
- `q3_inner_right_bp`
- `q3_outer_left_bp`
- `q3_outer_right_bp`

### `inv_detect_blocks` — 6 unregistered

- `q1_consensus_strength`
- `q1_n_scales_detected`
- `q1_staircase_rank`
- `q1_staircase_score`
- `q3_regime_change_left`
- `q3_regime_change_right`

### `burden` — 6 unregistered

- `q6_burden_beta`
- `q6_burden_ci_lower`
- `q6_burden_ci_upper`
- `q6_burden_p`
- `q6_n_del_sites`
- `q6_selection_signal`

### `existence_layer_b` — 6 unregistered

- `q7_bp1_concordance_kb`
- `q7_bp2_concordance_kb`
- `q7_n_delly_calls`
- `q7_n_manta_calls`
- `q7_n_sv_calls`
- `q7_sv_concordance_class`

### `synteny_dollo` — 4 unregistered

- `q5_dollo_state`
- `q5_outgroup_state`
- `q5_present_c_gariepinus`
- `q5_present_c_macrocephalus`

### `peel_diagnostic` — 4 unregistered

- `q7_peel_combined`
- `q7_peel_l1b_score`
- `q7_peel_l2_score`
- `q7_peel_survives`

### `gene_conversion_tracts` — 3 unregistered

- `q2_gc_n_snps_pass_qc_diagnostic`
- `q2_gc_total_samples_with_gc`
- `q2_gc_total_tracts`

### `flank_coherence` — 3 unregistered

- `q3_flank_coherence_class`
- `q3_left_flank_coherence_score`
- `q3_right_flank_coherence_score`

### `carrier_reconciliation` — 3 unregistered

- `q7_n_disagreeing_samples`
- `q7_n_reconciliation_sources`
- `q7_reconciliation_agreement`

### `hypothesis_verdict` — 3 unregistered

- `q7_t1_ratio`
- `q7_t2_eff_k`
- `q7_t3_retention`

### `band_composition` — 2 unregistered

- `q1_ancestry_div_hom_ref_vs_hom_inv`
- `q1_k_used`

### `existence_layer_a` — 2 unregistered

- `q7_cheat25_status`
- `q7_fdr_q_value`

## Section 3: Orphan unregistered keys (no schema owner)

Written via `add_evidence()` or dynamic paths. Not registered in spec because they're not analysis keys — they're annotations.

- `figure_breakpoint_diagnostic_pdf_path`
- `figure_breakpoint_diagnostic_png_path`
- `figure_stream_graph_bands_tsv_path`
- `figure_stream_graph_pdf_path`
- `figure_stream_graph_png_path`
- `scored`

These don't need spec registration — they're fine as-is (file-path annotations + `scored` self-marker).
