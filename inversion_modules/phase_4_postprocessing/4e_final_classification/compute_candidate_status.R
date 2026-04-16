#!/usr/bin/env Rscript
# =============================================================================
# compute_candidate_status.R — Three-axis tiering + completion calculator
# =============================================================================
#
# THE THREE AXES:
#   AXIS 1: CONFIDENCE TIER (is it real?)
#     → Pathway-based, from Q7 independence layers
#     → Tier 1 / 2 / 3 / 4 / SV-only
#
#   AXIS 2: COMPLETION (how much do we know?)
#     → 317 keys across 7 questions, percentage resolved
#     → Per-question breakdown: Q1–Q7 each has its own %
#
#   AXIS 3: EVOLUTIONARY CLASS (what kind of inversion?)
#     → From Q4 mechanism + Q5 age + Q6 frequency + Q3 boundaries
#     → Only assigned when enough keys are resolved
#
# READS:
#   - Registry key-value TSVs from each pipeline stage
#   - candidate_scores.tsv.gz (from C01d)
#   - hypothesis_verdicts.tsv.gz (from C01f)
#
# OUTPUTS:
#   - candidate_status.tsv           — one row per candidate, 3 axes
#   - candidate_completion.tsv       — per-question breakdown
#   - candidate_status_summary.txt   — human-readable report
#   - completion_heatmap.png         — visual summary
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

# =============================================================================
# REGISTRY KEY SPECIFICATION (317 keys across 7 questions)
# =============================================================================

build_key_spec <- function() {
  # Each key: name, question, can_be_NA (not_applicable in some contexts)
  keys <- list()

  # ── Q1: WHAT IS IT? (49 keys) ──
  q1 <- c(
    "q1_n_windows", "q1_span_kb", "q1_block_compactness", "q1_block_coherence",
    "q1_squareness", "q1_shape_class", "q1_n_children", "q1_landscape_category",
    "q1_nn_birth", "q1_nn_survives_40", "q1_nn_survives_80",
    "q1_n_matrix_variants", "q1_consensus_confidence",
    "q1_flat_inv_score", "q1_spiky_inv_score", "q1_fragmentation_score",
    "q1_s_het_mean", "q1_s_pve_mean", "q1_s_dip_mean",
    "q1_inv_likeness_mean", "q1_inv_likeness_max", "q1_inv_likeness_sd",
    "q1_family_likeness_mean", "q1_dosage_het_cv_mean",
    "q1_pve1_excess_mean", "q1_dip_pval_median",
    "q1_d01_block_strength", "q1_d02_block_shape", "q1_d03_nn_persistence",
    "q1_d04_decay_flatness", "q1_d05_interior_quality", "q1_d06_consensus",
    "q1_d07_sv_breakpoint", "q1_d08_peel_or_hyp", "q1_d09_pca_clusters",
    "q1_d10_partition", "q1_d11_boundary_concordance", "q1_d12_snake_concordance",
    "q1_composite_score", "q1_dim_positive",
    "q1_band_std_qgroup_fracs", "q1_band_het_qgroup_fracs", "q1_band_inv_qgroup_fracs",
    "q1_dominant_qgroup_std", "q1_dominant_qgroup_inv", "q1_q_group_overlap"
  )

  # ── Q2: WHAT'S HAPPENING INSIDE? (40 keys) ──
  q2 <- c(
    "q2_n_recombinant", "q2_n_gene_conversion", "q2_n_double_crossover",
    "q2_recomb_sample_ids",
    "q2_mean_switch_rate_std", "q2_mean_switch_rate_het",
    "q2_mean_switch_rate_inv", "q2_mean_switch_rate_rec",
    "q2_switching_kw_pval", "q2_recomb_family_spread",
    "q2_recomb_left_n", "q2_recomb_center_n", "q2_recomb_right_n",
    "q2_diversity_gradient_r2", "q2_diversity_gradient",
    "q2_profile_cor_observed", "q2_profile_cor_perm_p",
    "q2_class_entropy_mean", "q2_n_windows_trimodal",
    "q2_n_windows_bimodal", "q2_n_windows_unimodal",
    "q2_class_stability_pct",
    "q2_phase_block_n50_inside", "q2_phase_block_n50_flank",
    "q2_phase_block_ratio", "q2_phase_switch_rate_inside", "q2_phase_switch_rate_flank",
    "q2_n_indel_classes", "q2_indel_class_concordance", "q2_indel_ashman_D",
    "q2_chao1_hom_std", "q2_chao1_hom_inv",
    "q2_accum_saturated_std", "q2_accum_saturated_inv",
    "q2_n_founder_haplotypes_std", "q2_n_founder_haplotypes_inv"
  )

  # ── Q3: WHAT ARE THE BOUNDARIES DOING? (73 keys) ──
  q3_left <- c(
    "q3_left_bp", "q3_left_sv_bp", "q3_left_delta12_bp", "q3_left_fst_bp",
    "q3_left_hardness", "q3_left_sharpness", "q3_left_type",
    "q3_left_n_cheats", "q3_left_cheats_list", "q3_left_verdict",
    "q3_left_clip_count", "q3_left_depth_anomaly", "q3_left_is_fossil"
  )
  q3_right <- gsub("left", "right", q3_left)
  q3_concordance <- c(
    "q3_left_concordance_kb", "q3_right_concordance_kb",
    "q3_concordance_pval_left", "q3_concordance_pval_right",
    "q3_extended_suppression", "q3_suppression_extent_kb", "q3_fst_decay_rate"
  )
  q3_fst_profile <- c(
    paste0("q3_left_fst_", c("m200kb","m100kb","m50kb","0kb","p50kb","p100kb","p200kb","step")),
    paste0("q3_right_fst_", c("m200kb","m100kb","m50kb","0kb","p50kb","p100kb","p200kb","step"))
  )
  q3_hobs <- c(
    paste0("q3_left_hobs_", c("m100kb","0kb","p100kb","step")),
    paste0("q3_right_hobs_", c("m100kb","0kb","p100kb","step"))
  )
  q3_repeat <- c(
    "q3_left_repeat_density", "q3_left_repeat_class", "q3_left_gc_content",
    "q3_right_repeat_density", "q3_right_repeat_class", "q3_right_gc_content",
    "q3_genome_wide_repeat_density",
    "q3_repeat_enrichment_left", "q3_repeat_enrichment_right"
  )
  # NEW v2: Carrier reconciliation (PCA vs SV)
  q3_carrier <- c(
    "q3_n_carriers_pca", "q3_n_carriers_sv", "q3_carrier_concordance",
    "q3_n_pca_carrier_sv_ref", "q3_n_sv_carrier_pca_ref",
    "q3_dropout_rate", "q3_population_prior_applied", "q3_n_rescued_by_prior"
  )
  q3 <- c(q3_left, q3_right, q3_concordance, q3_fst_profile, q3_hobs, q3_repeat, q3_carrier)

  # ── Q4: HOW DID IT FORM? (47 keys) ──
  q4 <- c(
    "q4_has_inverted_sd", "q4_sd_left_start", "q4_sd_left_end",
    "q4_sd_right_start", "q4_sd_right_end", "q4_sd_length",
    "q4_sd_identity_pct", "q4_sd_orientation",
    "q4_biser2_concordance",
    "q4_junction_type_left", "q4_junction_type_right",
    "q4_mh_length_left", "q4_mh_length_right",
    "q4_te_family_left", "q4_te_family_right",
    "q4_te_enrichment", "q4_te_enrichment_fold",
    "q4_tr_enrichment",
    "q4_mechanism_support", "q4_mechanism_confidence", "q4_decision_tree_path",
    "q4_junction_seq_left_50bp", "q4_junction_seq_right_50bp",
    "q4_mh_sequence", "q4_te_name_left", "q4_te_name_right",
    "q4_te_same_family", "q4_te_same_orientation",
    "q4_sd_alignment_length", "q4_sd_n_mismatches", "q4_sd_divergence_pct",
    "q4_n_genes_inside", "q4_n_genes_spanning_bp", "q4_gene_names_at_bp",
    "q4_gene_density_inside", "q4_gene_density_genome", "q4_go_enrichment"
  )

  # ── Q5: HOW OLD IS IT? (39 keys) ──
  q5 <- c(
    "q5_gds_gap", "q5_gds_gap_percentile", "q5_fst_b1b3",
    "q5_gds_fst_spearman_rho", "q5_gds_fst_spearman_p",
    "q5_diversity_r2", "q5_diversity_shape",
    "q5_theta_pi_inside", "q5_theta_pi_flanking", "q5_theta_ratio",
    "q5_dollo_node", "q5_dollo_mya", "q5_n_species_sharing",
    "q5_theta_pi_std", "q5_theta_pi_inv", "q5_theta_pi_het",
    "q5_dxy_std_inv", "q5_dxy_std_inv_flanking", "q5_dxy_ratio",
    "q5_da_net_divergence",
    "q5_tajima_D_std", "q5_tajima_D_inv", "q5_tajima_D_pooled", "q5_tajima_D_pooled_note",
    "q5_segregating_sites_std", "q5_segregating_sites_inv",
    "q5_fixed_differences", "q5_shared_polymorphisms",
    "q5_origin_gds_dip_p", "q5_origin_class", "q5_origin_mechanism",
    "q5_gds_within_std", "q5_gds_within_inv", "q5_gds_between", "q5_gds_het_pattern"
  )

  # ── Q6: HOW COMMON IS IT? (28 keys) ──
  q6 <- c(
    "q6_freq_inv", "q6_n_HOM_STD", "q6_n_HET", "q6_n_HOM_INV",
    "q6_n_Recombinant", "q6_n_Unclassified", "q6_n_total",
    "q6_expected_het_hwe", "q6_observed_expected_ratio",
    "q6_hwe_deviation", "q6_hwe_chisq_p",
    "q6_freq_per_qgroup", "q6_freq_cv_across_qgroups",
    "q6_family_linkage", "q6_jackknife_max_delta", "q6_jackknife_sensitive_qgroup",
    "q6_tajima_D_inside", "q6_tajima_D_flanking", "q6_tajima_D_note",
    "q6_selection_pattern"
  )

  # ── Q7: IS IT REAL? (68 keys — updated for 4-layer framework) ──
  q7_layer_a <- c(
    "q7_layer_a_detected", "q7_layer_a_inv_likeness", "q7_layer_a_beta_pval",
    "q7_layer_a_beta_qval", "q7_layer_a_core_family", "q7_layer_a_pa_pattern"
  )
  q7_layer_b <- c(
    "q7_layer_b_detected", "q7_layer_b_delly", "q7_layer_b_manta",
    "q7_layer_b_bnd_triang", "q7_layer_b_n_carriers",
    "q7_layer_b_pe_support", "q7_layer_b_sr_support",
    "q7_layer_b_cipos", "q7_layer_b_ciend"
  )
  # NEW: Layer C = GHSL haplotype contrast (independent from PCA)
  q7_layer_c <- c(
    "q7_layer_c_ghsl_detected", "q7_layer_c_ghsl_contrast",
    "q7_layer_c_ghsl_n_pass", "q7_layer_c_ghsl_pct_pass",
    "q7_layer_c_ghsl_quality", "q7_layer_c_partition_stable",
    "q7_layer_c_ghsl_version"
  )
  # RENAMED: Layer D = genotype-breakpoint association (was Layer C in v1)
  q7_layer_d <- c(
    "q7_layer_d_tested", "q7_layer_d_fisher_or", "q7_layer_d_fisher_p",
    "q7_layer_d_armitage_z", "q7_layer_d_armitage_p", "q7_layer_d_concordance"
  )
  q7_independence <- c(
    "q7_n_layers_tested", "q7_n_layers_passed", "q7_independence_class"
  )
  q7_hypothesis <- c(
    "q7_t8_clair3_concordance", "q7_t9_jackknife_status", "q7_t9_max_delta",
    "q7_t10_theta_concordance"
  )
  q7_overall <- c(
    "q7_composite_score", "q7_tier", "q7_dim_positive",
    "q7_fdr_qval", "q7_verdict", "q7_verdict_confidence"
  )
  # NEW: SV caller audit keys (failure mode tracking)
  q7b_dropout <- c(
    "q7b_delly_raw_carriers", "q7b_delly_strict_carriers", "q7b_delly_carrier_loss",
    "q7b_manta_raw_carriers", "q7b_manta_pass_carriers", "q7b_manta_carrier_loss",
    "q7b_expected_dropout_pct", "q7b_observed_dropout_pct"
  )
  q7b_fragmentation <- c(
    "q7b_delly_inv_present", "q7b_delly_bnd_3to3", "q7b_delly_bnd_5to5",
    "q7b_manta_inv_present", "q7b_manta_bnd_inv3", "q7b_manta_bnd_inv5",
    "q7b_bnd_rescued", "q7b_bnd_rescue_concordance"
  )
  q7b_filter <- c(
    "q7b_delly_site_in_raw", "q7b_delly_site_passes_strict",
    "q7b_delly_site_qual", "q7b_delly_site_pe",
    "q7b_manta_site_in_raw", "q7b_manta_site_passes_pass", "q7b_manta_site_qual"
  )
  q7b_prior <- c(
    "q7b_pop_prior_applied", "q7b_pop_prior_freq_est", "q7b_pop_prior_n_rescued"
  )
  q7 <- c(q7_layer_a, q7_layer_b, q7_layer_c, q7_layer_d,
           q7_independence, q7_hypothesis, q7_overall,
           q7b_dropout, q7b_fragmentation, q7b_filter, q7b_prior)

  list(
    Q1 = q1, Q2 = q2, Q3 = q3, Q4 = q4, Q5 = q5, Q6 = q6, Q7 = q7,
    all_keys = c(q1, q2, q3, q4, q5, q6, q7),
    counts = c(Q1 = length(q1), Q2 = length(q2), Q3 = length(q3),
               Q4 = length(q4), Q5 = length(q5), Q6 = length(q6),
               Q7 = length(q7))
  )
}


# =============================================================================
# AXIS 1: CONFIDENCE TIER (pathway-based from Q7)
# =============================================================================
#
# NOT count-based anymore. Uses EVIDENCE PATHWAYS:
#
# Pathway A (PCA + SV convergence):
#   Layer A (PCA) strong + Layer B (SV) present + Layer C (association) tested
#   → Tier 1: VALIDATED
#
# Pathway B (PCA-only, strong):
#   Layer A strong + NN persistent + snake concordant + peel stable
#   Layer B absent (no SV) + Layer C not testable
#   → Tier 2: STRONG CANDIDATE (breakpoints unresolved)
#
# Pathway C (mixed/moderate):
#   Some evidence from multiple sources but not convergent
#   → Tier 3: SUPPORTED (needs resolution)
#
# Pathway D (weak/dead):
#   Cheat 25 DEAD, or noise, or family artifact
#   → Tier 4: WEAK
#
# Pathway E (SV-only):
#   No PCA core, but SV callers detect INV
#   → SV-1: both callers / SV-2: single caller
#
# =============================================================================

assign_confidence_tier <- function(keys) {
  # Extract Q7 layer results — 4 LAYERS (A/B/C/D)
  # A = Local PCA, B = SV callers, C = GHSL haplotype contrast, D = association
  layer_a <- as.logical(keys[["q7_layer_a_detected"]] %||% FALSE)
  layer_b <- as.logical(keys[["q7_layer_b_detected"]] %||% FALSE)
  layer_c <- as.logical(keys[["q7_layer_c_ghsl_detected"]] %||% FALSE)
  layer_d <- as.logical(keys[["q7_layer_d_tested"]] %||% FALSE)
  fisher_p <- as.numeric(keys[["q7_layer_d_fisher_p"]] %||% 1)
  n_layers <- as.integer(keys[["q7_n_layers_passed"]] %||% 0)
  ghsl_quality <- keys[["q7_layer_c_ghsl_quality"]] %||% "ABSENT"
  layer_c_strong <- layer_c && ghsl_quality %in% c("HIGH", "MODERATE")

  # PCA strength indicators
  d1 <- as.numeric(keys[["q1_d01_block_strength"]] %||% 0)
  d3 <- as.numeric(keys[["q1_d03_nn_persistence"]] %||% 0)
  d5 <- as.numeric(keys[["q1_d05_interior_quality"]] %||% 0)
  d8 <- as.numeric(keys[["q1_d08_peel_or_hyp"]] %||% 0)
  d12 <- as.numeric(keys[["q1_d12_snake_concordance"]] %||% 0)
  dim_pos <- as.integer(keys[["q1_dim_positive"]] %||% 0)

  # SV detail
  delly <- as.logical(keys[["q7_layer_b_delly"]] %||% FALSE)
  manta <- as.logical(keys[["q7_layer_b_manta"]] %||% FALSE)

  # Family linkage
  family <- keys[["q6_family_linkage"]] %||% "unknown"
  jackknife <- keys[["q7_t9_jackknife_status"]] %||% "unknown"

  # Verdict from C01f
  verdict <- keys[["q7_verdict"]] %||% "not_tested"

  # ── Pathway A: Full convergence (PCA + SV + GHSL + association) ──
  # 4/4 layers: the gold standard
  if (layer_a && layer_b && layer_c_strong && layer_d && fisher_p < 0.05) {
    if (family == "multi_family" || jackknife == "robust_multi_family") {
      return(list(tier = 1L, pathway = "A_full_4layer_convergence",
                  reason = "4/4 layers (PCA+SV+GHSL+Fisher) + multi-family"))
    } else {
      return(list(tier = 1L, pathway = "A_full_4layer_single_family",
                  reason = "4/4 layers converge (single family carrier)"))
    }
  }

  # 3/3 tested layers all pass (GHSL not available or not tested)
  if (layer_a && layer_b && layer_d && fisher_p < 0.05 && !layer_c) {
    n_tested <- 3L  # A, B, D tested; C not available
    if (family == "multi_family" || jackknife == "robust_multi_family") {
      return(list(tier = 1L, pathway = "A_strong_3layer_convergence",
                  reason = "3/3 tested layers (PCA+SV+Fisher) + multi-family; GHSL not tested"))
    } else {
      return(list(tier = 1L, pathway = "A_strong_3layer_single_family",
                  reason = "3/3 tested layers converge; GHSL not tested"))
    }
  }

  # ── Pathway B: Two independent methods agree, third missing ──

  # B1: PCA + GHSL convergence, no SV breakpoints
  if (layer_a && layer_c_strong && !layer_b && pca_strong) {
    return(list(tier = 2L, pathway = "B1_pca_ghsl_convergence",
                reason = "PCA + GHSL agree (independent data: dosage vs phased GT), no SV breakpoints"))
  }

  # B2: PCA + SV, no GHSL, association marginal
  if (layer_a && layer_b && !layer_c && d1 >= 0.5) {
    return(list(tier = 2L, pathway = "B2_pca_sv_partial",
                reason = "PCA + SV present, GHSL not tested, association weak/untested"))
  }

  # B3: Strong PCA only (no SV, no GHSL yet but strong structural evidence)
  if (layer_a && !layer_b && !layer_c && pca_strong && pca_supported) {
    return(list(tier = 2L, pathway = "B3_pca_only_strong",
                reason = "Strong PCA (D1+D3+D5) + peel/snake, no SV or GHSL yet"))
  }

  # ── Pathway C: Mixed evidence ──
  if (layer_a && dim_pos >= 4) {
    return(list(tier = 3L, pathway = "C_mixed",
                reason = paste0(dim_pos, "/12 dimensions, needs resolution")))
  }

  # GHSL-only (no PCA core, but GHSL found signal — unusual)
  if (layer_c_strong && !layer_a) {
    return(list(tier = 3L, pathway = "C_ghsl_only",
                reason = "GHSL signal without PCA core — investigate"))
  }

  # ── Pathway E: SV-only ──
  if (layer_b && !layer_a && !layer_c) {
    sv_tier <- if (delly && manta) "SV-1" else "SV-2"
    return(list(tier = sv_tier, pathway = "E_sv_only",
                reason = paste0("SV callers detect INV, no PCA/GHSL (",
                                if(delly) "DELLY" else "", if(delly&&manta) "+" else "",
                                if(manta) "Manta" else "", ")")))
  }

  # ── Pathway D: Weak/dead ──
  if (grepl("artifact|family_structure|H1_family", verdict)) {
    return(list(tier = 4L, pathway = "D_artifact",
                reason = paste0("Verdict: ", verdict)))
  }

  return(list(tier = 4L, pathway = "D_insufficient",
              reason = paste0("Only ", dim_pos, " dimensions positive")))
}


# =============================================================================
# AXIS 2: COMPLETION (per-question % of 317 keys resolved)
# =============================================================================

compute_completion <- function(keys, spec) {
  result <- list()
  total_resolved <- 0L
  total_applicable <- 0L

  for (q in paste0("Q", 1:7)) {
    q_keys <- spec[[q]]
    n_total <- length(q_keys)
    n_resolved <- 0L
    n_na <- 0L
    missing <- character()

    for (k in q_keys) {
      val <- keys[[k]]
      if (is.null(val) || (length(val) == 1 && is.na(val))) {
        missing <- c(missing, k)
      } else if (identical(val, "not_tested") || identical(val, "not_assessed")) {
        n_na <- n_na + 1L
      } else {
        n_resolved <- n_resolved + 1L
      }
    }

    applicable <- n_total - n_na
    pct <- if (applicable > 0) round(n_resolved / applicable * 100, 1) else 0
    result[[q]] <- list(
      resolved = n_resolved,
      applicable = applicable,
      total = n_total,
      pct = pct,
      n_missing = length(missing),
      top_missing = head(missing, 5)
    )
    total_resolved <- total_resolved + n_resolved
    total_applicable <- total_applicable + applicable
  }

  overall_pct <- if (total_applicable > 0) round(total_resolved / total_applicable * 100, 1) else 0

  list(
    overall = overall_pct,
    total_resolved = total_resolved,
    total_applicable = total_applicable,
    per_question = result
  )
}


# =============================================================================
# AXIS 3: EVOLUTIONARY CLASS (requires sufficient Q3+Q4+Q5+Q6 keys)
# =============================================================================

assign_evolutionary_class <- function(keys, completion) {
  # Minimum requirements: Q4 ≥ 30% AND Q5 ≥ 30% to assign class
  q4_pct <- completion$per_question$Q4$pct
  q5_pct <- completion$per_question$Q5$pct
  q6_pct <- completion$per_question$Q6$pct

  if (q4_pct < 30 || q5_pct < 30) {
    return(list(class = "unresolved",
                reason = paste0("Insufficient data (Q4=", q4_pct, "%, Q5=", q5_pct, "%)")))
  }

  # Read age, mechanism, frequency, boundary
  age <- keys[["q5_origin_class"]] %||% keys[["age_class"]] %||% "unranked"
  mechanism <- keys[["q4_mechanism_confidence"]] %||% "unclassified"
  freq <- as.numeric(keys[["q6_freq_inv"]] %||% NA)
  boundary <- keys[["boundary_status"]] %||% "not_assessed"
  gds_pct <- as.numeric(keys[["q5_gds_gap_percentile"]] %||% 50)

  # Classification
  if (!is.na(freq) && freq > 0.95) {
    return(list(class = "fixed_species_diagnostic",
                reason = "freq > 0.95, likely species-diagnostic"))
  }

  if (gds_pct >= 67 && grepl("eroded|soft", boundary)) {
    return(list(class = "old_polymorphic",
                reason = paste0("High GDS gap (", gds_pct, " pctl), eroded boundaries")))
  }

  if (gds_pct >= 67) {
    return(list(class = "old_polymorphic",
                reason = paste0("High GDS gap (", gds_pct, " pctl)")))
  }

  if (gds_pct <= 33 && grepl("hard_hard", boundary)) {
    return(list(class = "young_polymorphic",
                reason = paste0("Low GDS gap (", gds_pct, " pctl), sharp boundaries")))
  }

  if (gds_pct <= 33) {
    return(list(class = "young_polymorphic",
                reason = paste0("Low GDS gap (", gds_pct, " pctl)")))
  }

  n_children <- as.integer(keys[["q1_n_children"]] %||% 0)
  landscape <- keys[["q1_landscape_category"]] %||% "standard"
  if (n_children >= 2 || landscape == "complex_system") {
    return(list(class = "complex_nested",
                reason = paste0(n_children, " nested children, ", landscape)))
  }

  return(list(class = "intermediate_polymorphic",
              reason = paste0("GDS gap percentile ", gds_pct)))
}


# =============================================================================
# MAIN: Read registry, compute all three axes
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript compute_candidate_status.R <registry_dir> <outdir>\n")
  cat("  registry_dir: contains per-candidate key files (*.keys.tsv)\n")
  cat("  outdir: output directory\n")
  quit(save = "no", status = 1)
}

registry_dir <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Build key specification
spec <- build_key_spec()
cat("[STATUS] Registry specification: ", length(spec$all_keys), " keys across 7 questions\n")
cat("[STATUS]   Q1:", spec$counts["Q1"], "  Q2:", spec$counts["Q2"],
    "  Q3:", spec$counts["Q3"], "  Q4:", spec$counts["Q4"],
    "  Q5:", spec$counts["Q5"], "  Q6:", spec$counts["Q6"],
    "  Q7:", spec$counts["Q7"], "\n")

# Find all candidate key files
key_files <- list.files(registry_dir, pattern = "\\.keys\\.tsv$",
                         full.names = TRUE, recursive = TRUE)

# Also try reading from a single merged registry TSV
merged_registry <- file.path(registry_dir, "evidence_registry.tsv")
if (length(key_files) == 0 && file.exists(merged_registry)) {
  cat("[STATUS] Reading merged registry: ", merged_registry, "\n")
  reg_dt <- fread(merged_registry)
  # Convert to per-candidate key lists
  candidates <- unique(reg_dt$candidate_id)
} else if (length(key_files) > 0) {
  cat("[STATUS] Found ", length(key_files), " per-candidate key files\n")
  candidates <- gsub("\\.keys\\.tsv$", "", basename(key_files))
} else {
  # Try reading from candidate_scores.tsv.gz + other pipeline outputs
  cat("[STATUS] No key files found. Attempting to reconstruct from pipeline outputs...\n")

  # Look for standard pipeline output files
  score_file <- NULL
  for (f in c(
    file.path(registry_dir, "candidate_scores.tsv.gz"),
    file.path(registry_dir, "..", "candidate_scores", "candidate_scores.tsv.gz"),
    Sys.glob(file.path(registry_dir, "..", "*", "candidate_scores.tsv.gz"))
  )) {
    if (file.exists(f)) { score_file <- f; break }
  }

  if (is.null(score_file)) {
    stop("[STATUS] No registry data found in ", registry_dir)
  }

  cat("[STATUS] Reconstructing from: ", score_file, "\n")
  sdt <- fread(score_file)
  if (!"candidate_id" %in% names(sdt)) {
    sdt[, candidate_id := paste0(chrom, "_", interval_id)]
  }
  candidates <- sdt$candidate_id
}

cat("[STATUS] Processing ", length(candidates), " candidates...\n")

# ── Process each candidate ──

status_rows <- list()
completion_rows <- list()
report_lines <- character()

for (ci in seq_along(candidates)) {
  cid <- candidates[ci]

  # Load keys for this candidate
  keys <- list()

  # Method 1: per-candidate file
  kf <- file.path(registry_dir, paste0(cid, ".keys.tsv"))
  if (file.exists(kf)) {
    kdt <- fread(kf)
    for (ki in seq_len(nrow(kdt))) {
      keys[[kdt$key[ki]]] <- kdt$value[ki]
    }
  }

  # Method 2: from merged registry
  if (exists("reg_dt") && nrow(reg_dt) > 0) {
    cand_rows <- reg_dt[candidate_id == cid]
    for (ki in seq_len(nrow(cand_rows))) {
      keys[[cand_rows$key[ki]]] <- cand_rows$value[ki]
    }
  }

  # Method 3: from scores data.table (map column names to q1_ keys)
  if (exists("sdt") && cid %in% sdt$candidate_id) {
    row <- sdt[candidate_id == cid][1]
    col_to_key <- c(
      d1_block_strength = "q1_d01_block_strength",
      d2_block_shape = "q1_d02_block_shape",
      d3_nn_persistence = "q1_d03_nn_persistence",
      d4_decay_flatness = "q1_d04_decay_flatness",
      d5_interior_quality = "q1_d05_interior_quality",
      d6_consensus = "q1_d06_consensus",
      d7_sv_breakpoint = "q1_d07_sv_breakpoint",
      d8_peel_or_hyp = "q1_d08_peel_or_hyp",
      d9_pca_clusters = "q1_d09_pca_clusters",
      d10_partition = "q1_d10_partition",
      d11_boundary_concordance = "q1_d11_boundary_concordance",
      d12_snake_concordance = "q1_d12_snake_concordance",
      final_score = "q1_composite_score",
      dim_positive = "q1_dim_positive",
      shape_class = "q1_shape_class",
      span_mb = "q1_span_kb",  # needs conversion
      tier = "q7_tier"
    )
    for (col in names(col_to_key)) {
      if (col %in% names(row) && !is.na(row[[col]])) {
        val <- row[[col]]
        if (col == "span_mb") val <- as.numeric(val) * 1000  # Mb → kb
        keys[[col_to_key[col]]] <- as.character(val)
      }
    }
  }

  # ── Compute three axes ──
  conf <- assign_confidence_tier(keys)
  comp <- compute_completion(keys, spec)
  evol <- assign_evolutionary_class(keys, comp)

  # ── Store results ──
  status_rows[[ci]] <- data.table(
    candidate_id = cid,
    confidence_tier = as.character(conf$tier),
    confidence_pathway = conf$pathway,
    confidence_reason = conf$reason,
    completion_pct = comp$overall,
    completion_resolved = comp$total_resolved,
    completion_applicable = comp$total_applicable,
    evolutionary_class = evol$class,
    evolutionary_reason = evol$reason,
    # Per-question completion
    Q1_pct = comp$per_question$Q1$pct,
    Q2_pct = comp$per_question$Q2$pct,
    Q3_pct = comp$per_question$Q3$pct,
    Q4_pct = comp$per_question$Q4$pct,
    Q5_pct = comp$per_question$Q5$pct,
    Q6_pct = comp$per_question$Q6$pct,
    Q7_pct = comp$per_question$Q7$pct
  )

  # Per-question detail
  for (q in paste0("Q", 1:7)) {
    qd <- comp$per_question[[q]]
    completion_rows[[length(completion_rows) + 1]] <- data.table(
      candidate_id = cid, question = q,
      resolved = qd$resolved, applicable = qd$applicable,
      total = qd$total, pct = qd$pct,
      top_missing = paste(qd$top_missing, collapse = "; ")
    )
  }

  # Report line
  report_lines <- c(report_lines, sprintf(
    "%-25s  Tier %-4s  %5.1f%% complete  %-25s  [%s]",
    cid, conf$tier, comp$overall, evol$class, conf$pathway
  ))

  if (ci %% 10 == 0) cat("  ", ci, "/", length(candidates), "\n")
}

# ── Write outputs ──
status_dt <- rbindlist(status_rows, fill = TRUE)
compl_dt <- rbindlist(completion_rows, fill = TRUE)

fwrite(status_dt, file.path(outdir, "candidate_status.tsv"), sep = "\t")
fwrite(compl_dt, file.path(outdir, "candidate_completion.tsv"), sep = "\t")

# Human-readable report
report_header <- c(
  paste(rep("=", 100), collapse = ""),
  sprintf("CANDIDATE STATUS REPORT  (%d candidates, %d keys per candidate)",
          nrow(status_dt), length(spec$all_keys)),
  paste(rep("=", 100), collapse = ""),
  "",
  sprintf("%-25s  %-6s  %7s  %-25s  %s",
          "CANDIDATE", "TIER", "DONE", "EVOLUTION", "PATHWAY"),
  paste(rep("-", 100), collapse = ""),
  report_lines,
  "",
  paste(rep("=", 100), collapse = ""),
  "SUMMARY:",
  sprintf("  Tier 1 (validated):     %d", sum(status_dt$confidence_tier == "1")),
  sprintf("  Tier 2 (strong):        %d", sum(status_dt$confidence_tier == "2")),
  sprintf("  Tier 3 (supported):     %d", sum(status_dt$confidence_tier == "3")),
  sprintf("  Tier 4 (weak):          %d", sum(status_dt$confidence_tier == "4")),
  sprintf("  SV-only (SV-1/SV-2):    %d", sum(grepl("SV", status_dt$confidence_tier))),
  "",
  sprintf("  Median completion:      %.1f%%", median(status_dt$completion_pct)),
  sprintf("  Mean completion:        %.1f%%", mean(status_dt$completion_pct)),
  "",
  "  Per-question median completion:",
  sprintf("    Q1 (what is it):          %.1f%%", median(status_dt$Q1_pct)),
  sprintf("    Q2 (internal dynamics):   %.1f%%", median(status_dt$Q2_pct)),
  sprintf("    Q3 (boundaries):          %.1f%%", median(status_dt$Q3_pct)),
  sprintf("    Q4 (mechanism):           %.1f%%", median(status_dt$Q4_pct)),
  sprintf("    Q5 (age):                 %.1f%%", median(status_dt$Q5_pct)),
  sprintf("    Q6 (frequency):           %.1f%%", median(status_dt$Q6_pct)),
  sprintf("    Q7 (is it real):          %.1f%%", median(status_dt$Q7_pct)),
  "",
  "  Evolutionary classes:",
  paste0("    ", names(table(status_dt$evolutionary_class)), ": ",
         table(status_dt$evolutionary_class)),
  paste(rep("=", 100), collapse = "")
)

writeLines(report_header, file.path(outdir, "candidate_status_summary.txt"))

cat("\n")
cat(paste(report_header, collapse = "\n"))
cat("\n\n[DONE] → ", outdir, "\n")
