#!/usr/bin/env Rscript

# =============================================================================
# STEP32_candidate_multiscale_windows.R
#
# Multiscale re-analysis within fixed candidate intervals.
# Reruns local PCA at 250, 500, 1000 SNP windows to test whether the
# inversion-like signal remains coherent at larger window scales.
#
# NOT for genome-wide re-discovery. Same candidate coordinates, different
# internal SNP-window granularity.
#
# For each candidate and each scale:
#   1. Split candidate-region dosage into windows of size W
#   2. Compute local PCA per window (eigenvectors + eigenvalues)
#   3. Run k-means K=3 per window and classify topology
#   4. Compute window agreement matrix
#   5. Summarize scale-level coherence and compare across scales
#
# Inputs:
#   - Dosage matrix from DOSAGE_DIR
#   - Regional sites from DOSAGE_DIR
#   - Candidate table
#
# Outputs per candidate:
#   - candidate_multiscale_summary.tsv         one row per scale
#   - candidate_multiscale_windows_W<N>.tsv    per-window detail per scale
#   - candidate_multiscale_comparison.tsv      cross-scale comparison
#
# Usage:
#   Rscript STEP32_candidate_multiscale_windows.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(RSpectra)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

WINDOW_SIZES <- c(100L, 250L, 500L, 1000L)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

# ── Local PCA for one window ─────────────────────────────────────────────────
local_pca_window <- function(X_win, npc = 2) {
  # X_win: SNPs × samples
  n_snps <- nrow(X_win)
  n_samp <- ncol(X_win)
  if (n_snps < 5 || n_samp < 5) return(NULL)

  rm <- rowMeans(X_win, na.rm = TRUE)
  Xc <- X_win - rm
  Xc[!is.finite(Xc)] <- 0
  covmat <- cov(Xc, use = "pairwise.complete.obs")
  if (anyNA(covmat)) return(NULL)

  k <- min(npc, n_samp - 1)
  ee <- tryCatch(
    RSpectra::eigs_sym(covmat, k = k, which = "LM"),
    error = function(e) {
      ee2 <- eigen(covmat, symmetric = TRUE)
      list(values = ee2$values[seq_len(k)], vectors = ee2$vectors[, seq_len(k), drop = FALSE])
    }
  )

  list(eigval = ee$values, eigvec = ee$vectors)
}

# ── Classify one window ─────────────────────────────────────────────────────
classify_one_window <- function(eigvec, n_samples) {
  if (is.null(eigvec) || nrow(eigvec) < 5) return("insufficient")
  pc1 <- eigvec[, 1]

  set.seed(42)
  km <- tryCatch(
    kmeans(matrix(pc1, ncol = 1), centers = min(3, n_samples), nstart = 20),
    error = function(e) NULL
  )
  if (is.null(km)) return("clustering_failed")

  centers <- sort(km$centers[, 1])
  if (length(centers) < 3) return("weak")

  gap12 <- centers[2] - centers[1]
  gap23 <- centers[3] - centers[2]
  grp <- match(km$cluster, order(km$centers[, 1]))
  within_sd <- mean(sqrt(tapply(pc1, grp, var, na.rm = TRUE)), na.rm = TRUE)

  separation <- if (within_sd > 0) min(gap12, gap23) / within_sd else 0

  if (separation >= 2) return("clean_3band")
  if (separation >= 1.5) return("tilted")
  if (separation >= 1) return("split_het")
  if (separation >= 0.5) return("weak")
  return("ambiguous")
}

# ── Main processing loop ────────────────────────────────────────────────────
for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(cand_dir)

  # Load dosage
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) next

  dos <- fread(dos_file)
  sites <- fread(sites_file)

  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 100) next

  sites_reg <- sites[keep]
  dos_reg <- dos[keep]

  sample_cols <- setdiff(names(dos_reg), "marker")
  # Try Ind-style mapping
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (file.exists(rot_file)) {
    rot <- fread(rot_file)
    if (all(grepl("^Ind", sample_cols)) && length(sample_cols) == nrow(rot)) {
      setnames(dos_reg, old = sample_cols, new = rot$sample)
      sample_cols <- rot$sample
    } else {
      sample_cols <- intersect(rot$sample, sample_cols)
    }
  }
  if (length(sample_cols) < 10) next

  X <- as.matrix(dos_reg[, ..sample_cols])
  storage.mode(X) <- "double"
  n_snps_total <- nrow(X)
  n_samples <- ncol(X)

  message("[INFO] Candidate ", cid, " (", chr, "): ", n_snps_total,
          " SNPs × ", n_samples, " samples — multiscale analysis")

  scale_summaries <- list()

  for (W in WINDOW_SIZES) {
    n_win <- floor(n_snps_total / W)
    if (n_win < 2) {
      message("[INFO]   W=", W, ": only ", n_win, " windows — skipping")
      next
    }

    win_results <- list()
    win_groups <- list()

    for (w in seq_len(n_win)) {
      i1 <- (w - 1L) * W + 1L
      i2 <- w * W
      X_win <- X[i1:i2, , drop = FALSE]
      pos_start <- sites_reg$pos[i1]
      pos_end <- sites_reg$pos[i2]

      pca <- local_pca_window(X_win, npc = 2)
      topo <- classify_one_window(pca$eigvec, n_samples)

      # k-means assignments for agreement
      grp_assign <- NULL
      if (!is.null(pca$eigvec) && nrow(pca$eigvec) >= 5) {
        pc1 <- pca$eigvec[, 1]
        set.seed(42)
        km <- tryCatch(
          kmeans(matrix(pc1, ncol = 1), centers = min(3, n_samples), nstart = 20),
          error = function(e) NULL
        )
        if (!is.null(km)) {
          grp_assign <- match(km$cluster, order(km$centers[, 1]))
        }
      }

      win_results[[w]] <- data.table(
        candidate_id = cid, window_size = W, window_idx = w,
        pos_start = pos_start, pos_end = pos_end,
        n_snps = W, topology = topo
      )
      win_groups[[w]] <- grp_assign
    }

    win_dt <- rbindlist(win_results)

    # Compute mean pairwise agreement (ARI) across windows
    ari_vals <- c()
    for (i in seq_len(n_win)) {
      if (i < n_win) {
        for (j in (i + 1):n_win) {
          if (!is.null(win_groups[[i]]) && !is.null(win_groups[[j]]) &&
              length(win_groups[[i]]) == length(win_groups[[j]])) {
            tab <- table(win_groups[[i]], win_groups[[j]])
            n <- sum(tab)
            sc <- sum(choose(tab, 2))
            sa <- sum(choose(rowSums(tab), 2))
            sb <- sum(choose(colSums(tab), 2))
            ex <- sa * sb / choose(n, 2)
            mx <- 0.5 * (sa + sb)
            ari <- if (mx == ex) 1 else (sc - ex) / (mx - ex)
            ari_vals <- c(ari_vals, ari)
          }
        }
      }
    }

    n_clean <- sum(win_dt$topology == "clean_3band")
    frac_clean <- n_clean / n_win
    mean_ari <- if (length(ari_vals) > 0) round(mean(ari_vals, na.rm = TRUE), 4) else NA_real_

    scale_summaries[[as.character(W)]] <- data.table(
      candidate_id = cid, window_size = W,
      n_windows = n_win, n_clean_3band = n_clean,
      frac_clean = round(frac_clean, 4),
      mean_pairwise_ari = mean_ari,
      n_topologies = uniqueN(win_dt$topology),
      dominant_topology = win_dt[, .N, by = topology][which.max(N), topology]
    )

    fwrite(win_dt, file.path(cand_dir, paste0("candidate_multiscale_windows_W", W, ".tsv")),
           sep = "\t")

    message("[INFO]   W=", W, ": ", n_win, " windows, clean=",
            n_clean, " (", round(frac_clean * 100, 1), "%), mean_ARI=",
            ifelse(is.na(mean_ari), "NA", mean_ari))
  }

  if (length(scale_summaries) > 0) {
    ms_dt <- rbindlist(scale_summaries)

    # Cross-scale comparison: does signal persist?
    if (nrow(ms_dt) >= 2) {
      ms_dt[, signal_persistent := frac_clean >= 0.3]
      ms_dt[, scale_improves := c(NA, diff(frac_clean) > 0)]
    }

    fwrite(ms_dt, file.path(cand_dir, "candidate_multiscale_summary.tsv"), sep = "\t")
  }
}

message("[DONE] STEP32 complete")
