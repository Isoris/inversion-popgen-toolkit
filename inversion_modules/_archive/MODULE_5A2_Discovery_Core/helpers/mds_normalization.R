#!/usr/bin/env Rscript

# =============================================================================
# mds_normalization.R  (v8.0)
#
# ROBUST MDS SCALING UTILITIES for the inversion discovery framework.
#
# PROBLEM: Raw MDS coordinates and derived similarity scores are
# chromosome-specific and NOT comparable across chromosomes because:
#   1. make_sim_mat normalizes by within-chromosome max/quantile
#   2. z-scores use per-chromosome mean/sd
#   3. MDS spread depends on number and diversity of windows
#
# This module provides:
#   A. Explicit documentation of what is/isn't comparable
#   B. Robust normalization functions (quantile, IQR, MAD)
#   C. Cross-chromosome percentile-based threshold helpers
#   D. Diagnostic summaries of scale consistency
#
# WHAT IS CHROMOSOME-RELATIVE (not globally comparable):
#   - MDS coordinates (raw cmdscale output)
#   - z-scores (per-chromosome mean/sd)
#   - sim_mat values (normalized by per-chr distance quantile)
#   - continuity scores (depend on sim_mat)
#   - Snake 1 seed threshold crossing
#
# WHAT IS GLOBALLY COMPARABLE:
#   - window_id (from master registry)
#   - start_bp / end_bp (genomic coordinates)
#   - Snake 2 NN preservation (0-1 scale, comparable)
#   - Snake 3 GHSL rates (per-kb, comparable)
#   - Consensus support signatures (categorical)
#   - Rescue labels (categorical)
#
# WHAT THRESHOLDS ASSUME COMPARABILITY:
#   - S1_ACCEPT_THRESH (0.55) — applied to sim_mat-derived scores (chr-relative)
#   - SEED_MIN_Z (2.5) — applied to per-chr z-scores (chr-relative)
#   - NN_PASS_THRESH (0.60) — applied to NN preservation (globally comparable)
#   - GHSL_PASS_THRESH (2.0) — applied to per-kb rates (globally comparable)
#
# Usage:
#   source("mds_normalization.R")
#   # Then call functions directly
# =============================================================================

# =============================================================================
# A. ROBUST DISTANCE NORMALIZATION
# =============================================================================

#' Normalize a distance matrix using robust quantile scaling
#'
#' @param dmat square distance matrix
#' @param quantile_cap percentile for max distance (default 0.95)
#' @return list(sim = similarity matrix, dmax = normalization constant, method = string)
robust_sim_mat <- function(dmat, quantile_cap = 0.95) {
  finite_vals <- dmat[is.finite(dmat) & dmat > 0]
  if (length(finite_vals) == 0) {
    sim <- matrix(1, nrow(dmat), ncol(dmat))
    diag(sim) <- 1
    return(list(sim = sim, dmax = 1, method = "fallback_all_zero"))
  }

  dmax <- quantile(finite_vals, quantile_cap, na.rm = TRUE)
  if (!is.finite(dmax) || dmax <= 0) dmax <- max(finite_vals, na.rm = TRUE)
  if (!is.finite(dmax) || dmax <= 0) dmax <- 1

  sim <- 1 - pmin(dmat / dmax, 1)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1

  list(sim = sim, dmax = dmax, method = paste0("quantile_", quantile_cap))
}

#' Normalize MDS coordinates using IQR-based robust spread
#'
#' @param mds_mat N x K matrix of MDS coordinates
#' @return list(scaled = scaled matrix, spreads = per-axis spread, method = string)
robust_mds_spread <- function(mds_mat) {
  spreads <- apply(mds_mat, 2, function(x) {
    iq <- IQR(x, na.rm = TRUE)
    if (is.finite(iq) && iq > 0) {
      iq * 2.5  # ~covers 99% of normal distribution
    } else {
      rng <- diff(range(x, na.rm = TRUE))
      if (is.finite(rng) && rng > 0) rng else 1
    }
  })
  list(spreads = spreads, method = "IQR_2.5x")
}

#' Compute MAD-based z-scores (more robust than mean/sd)
#'
#' @param x numeric vector
#' @return z-scores using median and MAD
mad_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  m <- mad(x, na.rm = TRUE)
  if (!is.finite(m) || m == 0) m <- sd(x, na.rm = TRUE)
  if (!is.finite(m) || m == 0) return(rep(0, length(x)))
  (x - med) / m
}

# =============================================================================
# B. CROSS-CHROMOSOME PERCENTILE THRESHOLDS
# =============================================================================

#' Compute genome-wide percentile-based threshold from per-chromosome scores
#'
#' Instead of using a fixed threshold like 0.55, compute the threshold
#' as the Pth percentile of all scores across chromosomes.
#'
#' @param score_list named list of numeric vectors (one per chromosome)
#' @param percentile target percentile (e.g. 0.95 for top 5%)
#' @return list(threshold, per_chr_thresholds, global_distribution_summary)
compute_percentile_threshold <- function(score_list, percentile = 0.95) {
  all_scores <- unlist(score_list)
  all_scores <- all_scores[is.finite(all_scores)]

  global_thresh <- quantile(all_scores, percentile, na.rm = TRUE)

  per_chr <- vapply(score_list, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    quantile(x, percentile, na.rm = TRUE)
  }, numeric(1))

  list(
    global_threshold = as.numeric(global_thresh),
    per_chr_thresholds = per_chr,
    n_total = length(all_scores),
    global_mean = mean(all_scores),
    global_sd = sd(all_scores),
    global_median = median(all_scores),
    global_mad = mad(all_scores)
  )
}

# =============================================================================
# C. SCALE CONSISTENCY DIAGNOSTICS
# =============================================================================

#' Diagnose scale consistency across chromosomes
#'
#' Reports how much normalization constants vary, which helps assess
#' whether fixed thresholds are appropriate.
#'
#' @param per_chr_results list from STEP10 (per_chr structure)
#' @return data.table with per-chromosome scale diagnostics
diagnose_scale_consistency <- function(per_chr_results) {
  rows <- list()
  for (chr in names(per_chr_results)) {
    r <- per_chr_results[[chr]]
    if (is.null(r)) next

    dmat <- r$dmat
    mds_points <- if (!is.null(r$mds$points)) r$mds$points else
                  if (!is.null(r$mds)) as.matrix(r$mds) else NULL

    # Distance matrix scale
    finite_d <- dmat[is.finite(dmat) & dmat > 0]
    d_max <- if (length(finite_d) > 0) max(finite_d) else NA
    d_q95 <- if (length(finite_d) > 0) quantile(finite_d, 0.95) else NA
    d_median <- if (length(finite_d) > 0) median(finite_d) else NA

    # MDS spread
    mds_spread <- NA_real_
    if (!is.null(mds_points) && nrow(mds_points) > 1) {
      sp <- apply(mds_points, 2, function(x) IQR(x, na.rm = TRUE))
      mds_spread <- max(sp, na.rm = TRUE)
    }

    n_win <- nrow(dmat)

    rows[[chr]] <- data.table(
      chrom = chr,
      n_windows = n_win,
      dmat_max = round(as.numeric(d_max), 4),
      dmat_q95 = round(as.numeric(d_q95), 4),
      dmat_median = round(as.numeric(d_median), 4),
      dmat_q95_to_max_ratio = round(as.numeric(d_q95 / d_max), 4),
      mds_max_iqr_spread = round(mds_spread, 4)
    )
  }

  if (length(rows) == 0) return(data.table())
  dt <- rbindlist(rows)

  # Flag chromosomes with unusual scale
  if (nrow(dt) > 1) {
    med_q95 <- median(dt$dmat_q95, na.rm = TRUE)
    mad_q95 <- mad(dt$dmat_q95, na.rm = TRUE)
    if (is.finite(mad_q95) && mad_q95 > 0) {
      dt[, scale_z := round((dmat_q95 - med_q95) / mad_q95, 2)]
      dt[, scale_flag := fifelse(abs(scale_z) > 2, "UNUSUAL", "OK")]
    } else {
      dt[, scale_z := 0]
      dt[, scale_flag := "OK"]
    }
  }

  dt
}

message("[mds_normalization] Loaded robust MDS scaling utilities (v8.0)")
