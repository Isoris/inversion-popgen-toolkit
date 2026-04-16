#!/usr/bin/env Rscript
# =============================================================================
# cheat5_qaxis_fst_scan.R — Q-axis Fst scan
#
# For each window: use instant_q → local Q, then for each K component,
# split samples HIGH vs LOW and compute Fst. Compare to PC1-band Fst.
#
# family_fst_ratio ≈ 1.0 → Q grouping as good as PC1 → family-driven
# family_fst_ratio ≈ 0.0 → Q can't capture it → INVERSION
#
# REQUIRES: load_bridge.R sourced (get_Q, get_region_stats)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

Q_SPLIT_QUANTILE <- 0.30  # top/bottom 30% for HIGH/LOW split

#' Compute Q-axis Fst for one window
#' @param chr Chromosome
#' @param start_bp, end_bp Window bounds
#' @param sample_names CGA names
#' @param pc1_bands Named integer vector (sample → band 1/2/3)
#' @param K Number of ancestry components
#' @return data.table row
compute_qaxis_fst_window <- function(chr, start_bp, end_bp, sample_names,
                                      pc1_bands, K = 8L) {
  result <- data.table(
    chrom = chr, start_bp = start_bp, end_bp = end_bp,
    fst_pc1 = NA_real_, fst_q_best = NA_real_, fst_q_best_k = NA_integer_,
    family_fst_ratio = NA_real_, n_q_components_tested = 0L
  )

  # Get local Q for this window
  local_q <- tryCatch(get_Q(chr, start_bp, end_bp), error = function(e) NULL)
  if (is.null(local_q) || nrow(local_q) < 20) return(result)

  # PC1-band Fst (the standard grouping)
  band1 <- names(pc1_bands)[pc1_bands == 1]
  band3 <- names(pc1_bands)[pc1_bands == 3]
  if (length(band1) < 5 || length(band3) < 5) return(result)

  fst_pc1 <- tryCatch({
    s <- get_region_stats(chr, start_bp, end_bp, what = "Fst",
                           groups = list(band1 = band1, band3 = band3))
    s$Fst$band1_vs_band3
  }, error = function(e) NA_real_)

  result$fst_pc1 <- fst_pc1

  # Per-Q-component Fst
  q_fsts <- numeric(K)
  n_tested <- 0L

  for (k in seq_len(K)) {
    q_col <- paste0("Q", k)
    if (!q_col %in% names(local_q)) next

    q_vals <- local_q[[q_col]]
    names(q_vals) <- local_q$sample_id

    # Split: top 30% vs bottom 30%
    q_lo <- quantile(q_vals, Q_SPLIT_QUANTILE, na.rm = TRUE)
    q_hi <- quantile(q_vals, 1 - Q_SPLIT_QUANTILE, na.rm = TRUE)
    high_q <- names(q_vals)[q_vals >= q_hi]
    low_q  <- names(q_vals)[q_vals <= q_lo]

    if (length(high_q) < 5 || length(low_q) < 5) { q_fsts[k] <- NA; next }

    fst_k <- tryCatch({
      s <- get_region_stats(chr, start_bp, end_bp, what = "Fst",
                             groups = list(high = high_q, low = low_q))
      s$Fst$high_vs_low
    }, error = function(e) NA_real_)

    q_fsts[k] <- fst_k
    n_tested <- n_tested + 1L
  }

  if (n_tested == 0 || !is.finite(fst_pc1) || fst_pc1 == 0) return(result)

  fst_q_best <- max(q_fsts, na.rm = TRUE)
  fst_q_best_k <- which.max(q_fsts)

  ratio <- fst_q_best / fst_pc1
  # Clamp to [0, 2] for sanity
  ratio <- pmin(2, pmax(0, ratio))

  result$fst_q_best <- round(fst_q_best, 4)
  result$fst_q_best_k <- fst_q_best_k
  result$family_fst_ratio <- round(ratio, 4)
  result$n_q_components_tested <- n_tested

  result
}

#' Scan across windows for a chromosome, computing family_fst_ratio
#' Returns a data.table to merge into inv_like_dt
scan_qaxis_fst_chromosome <- function(chr, dt, sample_names, pc1_bands, K = 8L,
                                       step = 5L) {
  # Only scan a subset of windows for speed (every step-th window)
  scan_idx <- seq(1, nrow(dt), by = step)
  rows <- vector("list", length(scan_idx))

  for (ii in seq_along(scan_idx)) {
    wi <- scan_idx[ii]
    rows[[ii]] <- compute_qaxis_fst_window(
      chr, dt$start_bp[wi], dt$end_bp[wi], sample_names, pc1_bands, K
    )
    if (ii %% 200 == 0) message("[cheat5] ", chr, " window ", ii, "/", length(scan_idx))
  }

  out <- rbindlist(rows, fill = TRUE)
  # Interpolate to all windows (nearest-neighbor)
  if (nrow(out) > 0 && nrow(out) < nrow(dt)) {
    out[, mid_bp := (start_bp + end_bp) / 2]
    dt[, mid_bp := (start_bp + end_bp) / 2]
    setkey(out, mid_bp)
    setkey(dt, mid_bp)
    merged <- out[dt, roll = "nearest", on = "mid_bp"]
    merged[, mid_bp := NULL]
    return(merged[, .(chrom, start_bp = i.start_bp, end_bp = i.end_bp,
                       fst_pc1, fst_q_best, fst_q_best_k,
                       family_fst_ratio, n_q_components_tested)])
  }
  out
}
