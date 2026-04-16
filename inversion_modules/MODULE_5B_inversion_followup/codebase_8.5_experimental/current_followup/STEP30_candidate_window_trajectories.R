#!/usr/bin/env Rscript

# =============================================================================
# STEP30_candidate_window_trajectories.R
#
# Per-sample window-state trajectory layer.
#
# For each candidate, across ordered 100-SNP windows:
#   1. Assign per-sample per-window simplified state (AA/AB/BB/AMB)
#   2. Build sample × window state matrix
#   3. Count state switches per sample
#   4. Identify stable vs unstable samples
#   5. Cluster trajectories
#   6. Export regime tracks
#
# Inputs:
#   - STEP09 window_pca RDS files
#   - STEP20 window summary (overlapping windows)
#   - STEP21 candidate_pca_rotated.tsv (for coarse group labels)
#
# Outputs per candidate:
#   - candidate_window_state_matrix.tsv.gz   (sample × window matrix)
#   - candidate_trajectory_summary.tsv       (per-sample switch counts)
#   - candidate_trajectory_clusters.tsv      (trajectory-based grouping)
#
# Usage:
#   Rscript STEP30_candidate_window_trajectories.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

# Load STEP09 RDS files
step09_rds <- list.files(c(STEP09_DIR, INV_ROOT),
                         pattern = "\\.window_pca\\.rds$",
                         full.names = TRUE, recursive = TRUE)
step09_by_chr <- list()
for (f in step09_rds) {
  obj <- readRDS(f)
  step09_by_chr[[obj$chrom]] <- obj
}

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  if (!dir.exists(cand_dir)) next

  # Load window summary
  win_file <- file.path(cand_dir, "candidate_window_summary.tsv")
  if (!file.exists(win_file)) next
  win <- fread(win_file)
  if (nrow(win) < 3) next

  # Load coarse groups
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (!file.exists(rot_file)) next
  rot <- fread(rot_file)

  # Get STEP09 data
  step09_obj <- step09_by_chr[[chr]]
  if (is.null(step09_obj)) next

  sample_names <- step09_obj$sample_names
  pca_df <- step09_obj$pca
  n_samples <- length(sample_names)

  message("[INFO] Candidate ", cid, ": trajectory analysis (",
          nrow(win), " windows × ", n_samples, " samples)")

  # Build sample × window state matrix
  n_win <- nrow(win)
  state_mat <- matrix("AMB", nrow = n_samples, ncol = n_win)
  rownames(state_mat) <- sample_names
  colnames(state_mat) <- paste0("W", win$window_id)

  for (wi in seq_len(n_win)) {
    wid <- win$window_id[wi]
    row_idx <- which(pca_df$window_id == wid)
    if (length(row_idx) != 1) next

    # Extract PC1 loadings for this window
    pc1_cols <- paste0("PC_1_", sample_names)
    avail <- intersect(pc1_cols, names(pca_df))
    if (length(avail) != n_samples) next

    pc1 <- as.numeric(pca_df[row_idx, avail, with = FALSE])

    # k-means K=3 on PC1
    set.seed(42)
    km <- tryCatch(
      kmeans(matrix(pc1, ncol = 1), centers = min(3, n_samples), nstart = 20),
      error = function(e) NULL
    )
    if (is.null(km)) next

    # Order groups by centroid
    centers <- km$centers[, 1]
    ord <- order(centers)
    grp <- match(km$cluster, ord)

    # Map to AA/AB/BB
    labels <- c("AA", "AB", "BB")
    if (length(unique(grp)) == 3) {
      state_mat[, wi] <- labels[grp]
    } else if (length(unique(grp)) == 2) {
      state_mat[, wi] <- ifelse(grp == 1, "AA", "BB")
    }
  }

  # ── Per-sample trajectory summary ──────────────────────────────────────
  traj_list <- list()
  for (si in seq_len(n_samples)) {
    states <- state_mat[si, ]
    non_amb <- states[states != "AMB"]

    # Count switches
    switches <- if (length(non_amb) > 1) sum(non_amb[-1] != non_amb[-length(non_amb)]) else 0L

    # Dominant state
    state_tab <- table(non_amb)
    dominant <- if (length(state_tab) > 0) names(state_tab)[which.max(state_tab)] else "AMB"
    dominant_frac <- if (length(state_tab) > 0) max(state_tab) / sum(state_tab) else 0

    # Stability: 1 - switch_rate
    switch_rate <- if (length(non_amb) > 1) switches / (length(non_amb) - 1) else 0
    stability <- 1 - switch_rate

    # Fraction in each state
    frac_aa <- mean(states == "AA")
    frac_ab <- mean(states == "AB")
    frac_bb <- mean(states == "BB")
    frac_amb <- mean(states == "AMB")

    traj_list[[si]] <- data.table(
      candidate_id = cid,
      sample = sample_names[si],
      n_windows = n_win,
      n_non_ambiguous = length(non_amb),
      dominant_state = dominant,
      dominant_fraction = round(dominant_frac, 4),
      switch_count = switches,
      switch_rate = round(switch_rate, 4),
      stability = round(stability, 4),
      frac_AA = round(frac_aa, 4),
      frac_AB = round(frac_ab, 4),
      frac_BB = round(frac_bb, 4),
      frac_AMB = round(frac_amb, 4)
    )
  }

  traj_dt <- rbindlist(traj_list)

  # Merge with coarse group
  traj_dt <- merge(traj_dt,
                   rot[, .(sample, coarse_group_refined)],
                   by = "sample", all.x = TRUE)

  # Trajectory stability class
  traj_dt[, trajectory_class := fifelse(
    stability >= 0.8 & dominant_fraction >= 0.7, "stable",
    fifelse(stability >= 0.5, "moderate", "unstable")
  )]

  # ── Trajectory-based clustering ────────────────────────────────────────
  # Encode states as numeric for clustering: AA=0, AB=1, BB=2, AMB=NA
  state_numeric <- matrix(NA_real_, nrow = n_samples, ncol = n_win)
  state_numeric[state_mat == "AA"] <- 0
  state_numeric[state_mat == "AB"] <- 1
  state_numeric[state_mat == "BB"] <- 2

  # Hamming distance between trajectories
  if (n_samples >= 5 && n_win >= 3) {
    valid_cols <- colSums(is.finite(state_numeric)) >= n_samples * 0.5
    if (sum(valid_cols) >= 3) {
      sn_valid <- state_numeric[, valid_cols, drop = FALSE]
      sn_valid[is.na(sn_valid)] <- 1  # impute ambiguous as AB

      dmat <- dist(sn_valid, method = "manhattan")
      hc <- hclust(dmat, method = "ward.D2")
      k_traj <- min(5, n_samples - 1)
      traj_cluster <- cutree(hc, k = k_traj)
      traj_dt[, trajectory_cluster := traj_cluster[match(sample, sample_names)]]
    }
  }

  # ── Save outputs ───────────────────────────────────────────────────────
  # State matrix
  state_out <- data.table(sample = sample_names)
  for (wi in seq_len(n_win)) {
    state_out[[paste0("W", win$window_id[wi])]] <- state_mat[, wi]
  }
  fwrite(state_out, file.path(cand_dir, "candidate_window_state_matrix.tsv.gz"), sep = "\t")

  # Trajectory summary
  fwrite(traj_dt, file.path(cand_dir, "candidate_trajectory_summary.tsv"), sep = "\t")

  # Summary
  stable_n <- sum(traj_dt$trajectory_class == "stable")
  moderate_n <- sum(traj_dt$trajectory_class == "moderate")
  unstable_n <- sum(traj_dt$trajectory_class == "unstable")
  message("[INFO]   Trajectories: stable=", stable_n,
          " moderate=", moderate_n, " unstable=", unstable_n)
}

message("[DONE] STEP30 complete")
