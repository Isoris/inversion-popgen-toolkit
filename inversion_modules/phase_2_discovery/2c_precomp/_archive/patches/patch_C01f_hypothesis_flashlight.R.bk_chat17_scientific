#!/usr/bin/env Rscript
# =============================================================================
# PATCH: Flashlight integration for STEP_C01f_hypothesis_tests.R
#
# Adds test T9: SV GENOTYPE CONCORDANCE
#
# Do DELLY/Manta per-sample genotypes predict PCA cluster membership?
#
# Method:
#   1. Load flashlight anchors for the candidate region
#   2. Get PCA band assignments from precomp (k=3 on PC1)
#   3. Build 2×3 contingency table:
#        rows = SV genotype (HOM_REF, HET, HOM_INV)
#        cols = PCA band (1, 2, 3)
#   4. Fisher exact test: is SV genotype associated with PCA band?
#   5. If significant → PCA clusters = real inversion genotypes (supports H2/H3)
#   6. If not significant → PCA clusters ≠ SV genotypes (supports H1/H5)
#
# Additional sub-test: Cheat 2 verification
#   For samples with het-DEL at breakpoint: are they enriched in band 1/2
#   (reference/het) vs band 3 (HOM_INV)?
#   If yes → breakpoint evidence independently confirms the inversion model.
#
# INSERT: After existing tests (T1-T8), before verdict computation.
# =============================================================================

# ─── Source flashlight ───────────────────────────────────────────────

fl_loader <- Sys.getenv("FLASHLIGHT_LOADER", "")
if (!nzchar(fl_loader)) {
  for (p in c(
    file.path(dirname(outdir), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path),
    file.path(dirname(dirname(outdir)), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path)
  )) if (file.exists(p)) { fl_loader <- p; break }
}
.hyptest_has_flashlight <- FALSE
if (nzchar(fl_loader) && file.exists(fl_loader)) {
  tryCatch({
    source(fl_loader)
    .hyptest_has_flashlight <- TRUE
    message("[C01f] Flashlight sourced — T9 SV concordance test available")
  }, error = function(e) message("[C01f] Flashlight: ", e$message))
}

# ─── T9: SV GENOTYPE CONCORDANCE TEST ───────────────────────────────

test_sv_concordance <- function(chr, start_bp, end_bp, pc, sample_names,
                                 real_names = NULL) {
  result <- list(
    test = "T9_sv_concordance",
    p_value = NA_real_,
    odds_ratio = NA_real_,
    n_anchors = 0L,
    n_concordant = 0L,
    concordance_rate = NA_real_,
    orientation = NA_character_,
    cheat2_p = NA_real_,
    cheat2_n_del = 0L,
    evidence = "no_flashlight",
    supports = "neutral"
  )

  if (!.hyptest_has_flashlight) return(result)

  fl <- load_flashlight(chr)
  if (is.null(fl) || nrow(fl$inv_calls) == 0) {
    result$evidence <- "no_inv_calls"
    return(result)
  }

  # Get anchor samples
  anchors <- get_sv_anchors(chr, start_bp, end_bp, min_confidence = "LOW")
  if (nrow(anchors) < 10) {
    result$evidence <- paste0("too_few_anchors_", nrow(anchors))
    return(result)
  }

  # Get PCA bands from precomp
  dt <- pc$dt
  win_idx <- which(dt$start_bp >= start_bp & dt$end_bp <= end_bp)
  if (length(win_idx) < 5) {
    result$evidence <- "too_few_windows"
    return(result)
  }

  # Use precomp sample names for PC columns
  pc_names <- sample_names
  pc1_cols <- paste0("PC_1_", pc_names)
  available <- intersect(pc1_cols, names(dt))
  if (length(available) < 20) {
    result$evidence <- "too_few_pc_columns"
    return(result)
  }

  # k=3 bands from average PC1 loadings in the candidate region
  mat <- as.matrix(dt[win_idx, ..available])
  avg <- colMeans(mat, na.rm = TRUE)
  valid <- is.finite(avg)
  if (sum(valid) < 20) {
    result$evidence <- "insufficient_valid_pc1"
    return(result)
  }
  vals <- avg[valid]
  snames <- sub("^PC_1_", "", names(vals))

  km <- tryCatch(kmeans(vals, centers = 3, nstart = 10), error = function(e) NULL)
  if (is.null(km)) {
    result$evidence <- "kmeans_failed"
    return(result)
  }

  co <- order(km$centers[, 1])
  bands <- integer(length(vals))
  for (b in 1:3) bands[km$cluster == co[b]] <- b
  names(bands) <- snames

  # Map precomp names to CGA names for matching with anchors
  if (!is.null(real_names) && length(real_names) == length(sample_names)) {
    name_map <- setNames(real_names, sample_names)
    band_cga <- setNames(bands, name_map[names(bands)])
  } else {
    band_cga <- bands
  }

  # Build contingency table: SV genotype × PCA band
  anchors_unique <- anchors[!duplicated(sample_id)]
  matched <- anchors_unique[sample_id %in% names(band_cga)]
  if (nrow(matched) < 10) {
    result$evidence <- paste0("few_matched_", nrow(matched))
    return(result)
  }

  matched[, pca_band := band_cga[sample_id]]
  matched <- matched[!is.na(pca_band)]

  # Try both orientations: REF→1,HET→2,INV→3 and REF→3,HET→2,INV→1
  expected_fwd <- c(HOM_REF = 1L, HET = 2L, HOM_INV = 3L)
  expected_rev <- c(HOM_REF = 3L, HET = 2L, HOM_INV = 1L)

  matched[, exp_fwd := expected_fwd[sv_genotype]]
  matched[, exp_rev := expected_rev[sv_genotype]]
  n_fwd <- sum(matched$pca_band == matched$exp_fwd, na.rm = TRUE)
  n_rev <- sum(matched$pca_band == matched$exp_rev, na.rm = TRUE)
  n_best <- max(n_fwd, n_rev)
  orientation <- if (n_fwd >= n_rev) "forward" else "reversed"
  concordance_rate <- n_best / nrow(matched)

  # Fisher exact test on 3×3 contingency table
  # Simplify to 2×2 for power: (REF vs INV) × (band1 vs band3)
  ref_samples <- matched[sv_genotype == "HOM_REF"]
  inv_samples <- matched[sv_genotype == "HOM_INV"]

  if (nrow(ref_samples) >= 3 && nrow(inv_samples) >= 3) {
    if (orientation == "forward") {
      a <- sum(ref_samples$pca_band == 1)  # REF in band 1 (expected)
      b <- sum(ref_samples$pca_band == 3)  # REF in band 3 (unexpected)
      c <- sum(inv_samples$pca_band == 1)  # INV in band 1 (unexpected)
      d <- sum(inv_samples$pca_band == 3)  # INV in band 3 (expected)
    } else {
      a <- sum(ref_samples$pca_band == 3)
      b <- sum(ref_samples$pca_band == 1)
      c <- sum(inv_samples$pca_band == 3)
      d <- sum(inv_samples$pca_band == 1)
    }

    ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2)),
                   error = function(e) NULL)
    if (!is.null(ft)) {
      result$p_value <- ft$p.value
      result$odds_ratio <- if (is.finite(ft$estimate)) ft$estimate else NA_real_
    }
  }

  result$n_anchors <- nrow(matched)
  result$n_concordant <- n_best
  result$concordance_rate <- round(concordance_rate, 3)
  result$orientation <- orientation

  # ── Cheat 2 sub-test: het-DEL breakpoint enrichment ──
  bp_dels_left  <- get_breakpoint_dels(chr, start_bp, window = 50000L)
  bp_dels_right <- get_breakpoint_dels(chr, end_bp, window = 50000L)
  all_del_carriers <- unique(c(
    unlist(bp_dels_left$het_carriers),
    unlist(bp_dels_right$het_carriers)
  ))
  del_carriers_in_region <- intersect(all_del_carriers, names(band_cga))

  if (length(del_carriers_in_region) >= 3) {
    del_bands <- band_cga[del_carriers_in_region]
    # Cheat 2 predicts: het-DEL carriers should be in band 1 or 2 (not 3)
    # if orientation is forward, or band 3 or 2 if reversed
    if (orientation == "forward") {
      n_expected <- sum(del_bands %in% c(1, 2))  # REF or HET
      n_unexpected <- sum(del_bands == 3)          # HOM_INV (unexpected)
    } else {
      n_expected <- sum(del_bands %in% c(2, 3))
      n_unexpected <- sum(del_bands == 1)
    }

    # Binomial test: is fraction in expected bands higher than chance?
    bt <- tryCatch(binom.test(n_expected, n_expected + n_unexpected, p = 2/3,
                               alternative = "greater"),
                   error = function(e) NULL)
    if (!is.null(bt)) result$cheat2_p <- bt$p.value
    result$cheat2_n_del <- length(del_carriers_in_region)
  }

  # Interpretation
  if (!is.na(result$p_value) && result$p_value < 0.001 && concordance_rate > 0.70) {
    result$evidence <- "strong_concordance"
    result$supports <- "H2_or_H3"  # real inversion
  } else if (!is.na(result$p_value) && result$p_value < 0.05 && concordance_rate > 0.50) {
    result$evidence <- "moderate_concordance"
    result$supports <- "H2_or_H3"
  } else if (!is.na(result$concordance_rate) && concordance_rate < 0.35) {
    result$evidence <- "no_concordance"
    result$supports <- "H1"  # family structure
  } else {
    result$evidence <- "inconclusive"
    result$supports <- "neutral"
  }

  result
}

# ─── USAGE: Inside the main test loop ───────────────────────────────
# After existing tests T1-T8:
#
#   # T9: SV genotype concordance (flashlight)
#   t9 <- test_sv_concordance(chr, start_bp, end_bp, pc, precomp_sample_names,
#                              real_names)
#   all_results[[length(all_results) + 1]] <- data.table(
#     chrom = chr, interval_id = iid, test = "T9_sv_concordance",
#     statistic = t9$concordance_rate, p_value = t9$p_value,
#     effect_size = t9$odds_ratio, n_samples = t9$n_anchors,
#     evidence = t9$evidence, supports = t9$supports,
#     detail = paste0("n_anchors=", t9$n_anchors,
#                     " conc=", t9$concordance_rate,
#                     " orient=", t9$orientation,
#                     " cheat2_p=", round(t9$cheat2_p, 4),
#                     " cheat2_n=", t9$cheat2_n_del)
#   )
#
# In verdict computation, add T9 to the evidence aggregation:
#
#   # T9 affects verdict: strong SV concordance → supports real inversion
#   if (t9$supports == "H2_or_H3") {
#     inv_evidence <- inv_evidence + 2  # strong boost
#   }
#   if (!is.na(t9$cheat2_p) && t9$cheat2_p < 0.05) {
#     inv_evidence <- inv_evidence + 1  # additional Cheat 2 boost
#   }
