#!/usr/bin/env Rscript

# =============================================================================
# STEP20b_anchor_geometry_module.R
#
# Anchor-based (u,v) transformation + diagnostic split experiments.
#
# For each candidate:
#   1. Load regional PCA from STEP12 (PC1, PC2)
#   2. Identify dense core points for left/middle/right anchors
#   3. Compute centroids → define (u,v) coordinate system
#   4. Project ALL samples into anchor space
#   5. Compute per-sample geometric features:
#      u, v, center_distance, broad_axis_position, left_affinity,
#      right_affinity, branch_balance, ambiguity_score
#   6. Run diagnostic split experiments:
#      VIEW A: all samples (reference)
#      VIEW B: homo-only (hidden structure in extremes)
#      VIEW C: het-only (hidden structure in bridge)
#
# Usage:
#   Rscript STEP20b_anchor_geometry_module.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

has_ggplot <- suppressWarnings(require(ggplot2, quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]
message("[INFO] Processing ", nrow(cand), " candidates")

# ── Helper: identify dense core points ──────────────────────────────────────
# Uses density estimation on PC1 to find the core of each group
find_core_points <- function(pc1, group, frac = 0.5) {
  # For each group, keep the central `frac` of points by PC1 density
  core_flags <- rep(FALSE, length(pc1))
  for (g in unique(group)) {
    idx <- which(group == g)
    if (length(idx) < 5) { core_flags[idx] <- TRUE; next }
    vals <- pc1[idx]
    q_lo <- quantile(vals, (1 - frac) / 2)
    q_hi <- quantile(vals, 1 - (1 - frac) / 2)
    core_flags[idx] <- vals >= q_lo & vals <= q_hi
  }
  core_flags
}

# ── Helper: compute anchor projection features ─────────────────────────────
compute_anchor_features <- function(pc1, pc2, group) {
  groups <- unique(group)
  n_groups <- length(groups)

  # Group centroids ordered by PC1
  centroids <- data.table(group = group, pc1 = pc1, pc2 = pc2)[
    , .(c1 = mean(pc1), c2 = mean(pc2), n = .N), by = group
  ][order(c1)]

  if (n_groups < 2) {
    return(data.table(
      u = pc1, v = pc2,
      center_distance = NA_real_, broad_axis_position = pc1,
      left_affinity = NA_real_, right_affinity = NA_real_,
      branch_balance = NA_real_, ambiguity_score = NA_real_
    ))
  }

  # Identify left / middle / right
  if (n_groups >= 3) {
    C_L <- c(centroids$c1[1], centroids$c2[1])
    C_M <- c(centroids$c1[2], centroids$c2[2])  # middle centroid
    C_R <- c(centroids$c1[n_groups], centroids$c2[n_groups])
  } else {
    # k=2: use midpoint as pseudo-middle
    C_L <- c(centroids$c1[1], centroids$c2[1])
    C_R <- c(centroids$c1[2], centroids$c2[2])
    C_M <- (C_L + C_R) / 2
  }

  # Stripe axis: from C_L to C_R
  axis_vec <- C_R - C_L
  axis_len <- sqrt(sum(axis_vec^2))
  if (axis_len < 1e-10) axis_len <- 1
  axis_unit <- axis_vec / axis_len

  # Perpendicular direction
  perp_unit <- c(-axis_unit[2], axis_unit[1])

  # Project all samples relative to C_M
  centered_pc1 <- pc1 - C_M[1]
  centered_pc2 <- pc2 - C_M[2]

  # u = projection onto stripe axis (broad_axis_position)
  u <- centered_pc1 * axis_unit[1] + centered_pc2 * axis_unit[2]
  # v = projection onto perpendicular (deviation from axis)
  v <- centered_pc1 * perp_unit[1] + centered_pc2 * perp_unit[2]

  # Center distance
  center_distance <- sqrt(centered_pc1^2 + centered_pc2^2)

  # Left and right affinity (distance to left/right cores)
  dist_to_left <- sqrt((pc1 - C_L[1])^2 + (pc2 - C_L[2])^2)
  dist_to_right <- sqrt((pc1 - C_R[1])^2 + (pc2 - C_R[2])^2)

  # Normalize affinities (0 = far, 1 = at centroid)
  max_dist <- max(c(dist_to_left, dist_to_right), na.rm = TRUE)
  if (max_dist < 1e-10) max_dist <- 1
  left_affinity <- 1 - dist_to_left / max_dist
  right_affinity <- 1 - dist_to_right / max_dist

  # Branch balance: positive = right-leaning, negative = left-leaning
  branch_balance <- right_affinity - left_affinity

  # Ambiguity: high when both affinities are similar (sample is between groups)
  extreme_affinity <- pmax(left_affinity, right_affinity)
  ambiguity_score <- 1 - abs(branch_balance) / pmax(extreme_affinity, 0.01)
  ambiguity_score <- pmax(pmin(ambiguity_score, 1), 0)

  data.table(
    u = round(u, 6),
    v = round(v, 6),
    center_distance = round(center_distance, 6),
    broad_axis_position = round(u, 6),  # same as u
    left_affinity = round(left_affinity, 6),
    right_affinity = round(right_affinity, 6),
    branch_balance = round(branch_balance, 6),
    ambiguity_score = round(ambiguity_score, 6),
    dist_to_left_core = round(dist_to_left, 6),
    dist_to_right_core = round(dist_to_right, 6)
  )
}

# ── Helper: run split re-analysis ───────────────────────────────────────────
split_reanalysis <- function(X, samples_keep, label) {
  if (length(samples_keep) < 5) {
    return(data.table(
      n_samples = length(samples_keep), label = label,
      axis1_var = NA_real_, axis2_var = NA_real_,
      internal_cluster_count = NA_integer_,
      compactness = NA_real_, fragmentation_flag = FALSE, notes = "too_few_samples"
    ))
  }

  X_sub <- X[, samples_keep, drop = FALSE]
  pc <- tryCatch(prcomp(t(X_sub), center = TRUE, scale. = FALSE),
                 error = function(e) NULL)
  if (is.null(pc)) {
    return(data.table(
      n_samples = length(samples_keep), label = label,
      axis1_var = NA_real_, axis2_var = NA_real_,
      internal_cluster_count = NA_integer_,
      compactness = NA_real_, fragmentation_flag = FALSE, notes = "pca_failed"
    ))
  }

  var_expl <- pc$sdev^2 / sum(pc$sdev^2)
  pc1 <- pc$x[, 1]

  # Try k=2,3 and pick by WSS ratio
  km2 <- tryCatch(kmeans(matrix(pc1, ncol = 1), centers = 2, nstart = 20),
                   error = function(e) NULL)
  km3 <- tryCatch(kmeans(matrix(pc1, ncol = 1), centers = 3, nstart = 20),
                   error = function(e) NULL)

  best_k <- 1L
  if (!is.null(km2) && !is.null(km3)) {
    wss_ratio <- km3$tot.withinss / km2$tot.withinss
    best_k <- if (wss_ratio < 0.7) 3L else 2L
  } else if (!is.null(km2)) {
    best_k <- 2L
  }

  separation <- if (best_k >= 2 && !is.null(km2)) {
    centers <- sort(km2$centers[, 1])
    if (length(centers) >= 2) {
      gap <- abs(diff(centers))
      pooled_sd <- sqrt(km2$tot.withinss / max(length(pc1) - 2, 1))
      if (pooled_sd > 0) gap / pooled_sd else Inf
    } else NA_real_
  } else NA_real_

  data.table(
    n_samples = length(samples_keep),
    label = label,
    axis1_var = round(var_expl[1], 4),
    axis2_var = round(if (length(var_expl) >= 2) var_expl[2] else 0, 4),
    internal_cluster_count = best_k,
    compactness = round(sd(pc1), 4),
    fragmentation_flag = best_k >= 3,
    separation_score = round(separation, 4),
    notes = ""
  )
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════════

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  plot_dir <- file.path(PLOTS_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(cand_dir)
  ensure_dir(plot_dir)

  message("\n[INFO] Candidate ", cid, " — anchor geometry + split experiments")

  # ── Load STEP12 regional PCA ─────────────────────────────────────────────
  pca_files <- list.files(STEP12_DIR,
                          pattern = paste0("candidate_", cid, "\\.regional_pca_samples"),
                          full.names = TRUE)
  if (length(pca_files) == 0) { message("[SKIP] No STEP12 PCA for candidate ", cid); next }
  step12 <- fread(pca_files[1])
  if (nrow(step12) < 10) { message("[SKIP] Too few samples"); next }

  pc1 <- step12$PC1
  pc2 <- step12$PC2
  samples <- step12$sample

  # Coarse group
  grp_col <- intersect(c("group_label", "ordered_group"), names(step12))
  if (length(grp_col) == 0) { message("[SKIP] No group column"); next }
  group <- as.character(step12[[grp_col[1]]])

  # ── Find core points ──────────────────────────────────────────────────────
  core_flag <- find_core_points(pc1, group, frac = 0.6)
  message("[INFO] Core points: ", sum(core_flag), " / ", length(core_flag))

  # ── Compute anchor features using core-defined centroids ─────────────────
  # Build centroids from core points only
  core_centroids <- data.table(
    group = group[core_flag], pc1 = pc1[core_flag], pc2 = pc2[core_flag]
  )[, .(c1 = mean(pc1), c2 = mean(pc2), n = .N), by = group][order(c1)]

  # But project ALL samples
  features <- compute_anchor_features(pc1, pc2, group)

  # Build output table
  anchor_dt <- data.table(
    sample = samples,
    candidate_id = cid,
    PC1 = round(pc1, 6),
    PC2 = round(pc2, 6),
    coarse_group = group,
    is_core = core_flag
  )
  anchor_dt <- cbind(anchor_dt, features)

  # Add provisional group from anchor features
  anchor_dt[, provisional_anchor_group := fifelse(
    broad_axis_position < quantile(broad_axis_position, 0.25), "LEFT",
    fifelse(broad_axis_position > quantile(broad_axis_position, 0.75), "RIGHT", "MIDDLE")
  )]

  fwrite(anchor_dt,
         file.path(cand_dir, "candidate_anchor_projection.tsv"), sep = "\t")

  # ── Anchor summary ───────────────────────────────────────────────────────
  anchor_summary <- data.table(
    candidate_id = cid,
    n_left_core = core_centroids[1, n],
    n_mid_core = if (nrow(core_centroids) >= 3) core_centroids[2, n] else NA_integer_,
    n_right_core = core_centroids[nrow(core_centroids), n],
    left_mid_distance = if (nrow(core_centroids) >= 3)
      round(sqrt((core_centroids$c1[1] - core_centroids$c1[2])^2 +
                 (core_centroids$c2[1] - core_centroids$c2[2])^2), 4) else NA_real_,
    right_mid_distance = if (nrow(core_centroids) >= 3)
      round(sqrt((core_centroids$c1[nrow(core_centroids)] - core_centroids$c1[2])^2 +
                 (core_centroids$c2[nrow(core_centroids)] - core_centroids$c2[2])^2), 4) else NA_real_,
    left_right_distance = round(sqrt(
      (core_centroids$c1[1] - core_centroids$c1[nrow(core_centroids)])^2 +
      (core_centroids$c2[1] - core_centroids$c2[nrow(core_centroids)])^2), 4),
    anchor_quality = if (nrow(core_centroids) >= 3 &&
                          all(core_centroids$n >= 5)) "good" else "limited"
  )

  fwrite(anchor_summary,
         file.path(cand_dir, "candidate_anchor_summary.tsv"), sep = "\t")

  # ═════════════════════════════════════════════════════════════════════════
  # DIAGNOSTIC SPLIT EXPERIMENTS
  # ═════════════════════════════════════════════════════════════════════════

  # Load dosage for split experiments
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file)) { message("[SKIP] No dosage for splits"); next }

  dos <- fread(dos_file)
  sites <- fread(sites_file)
  dos_sample_cols <- setdiff(names(dos), "marker")

  # Map Ind-style if needed
  if (all(grepl("^Ind", dos_sample_cols)) && length(dos_sample_cols) == length(samples)) {
    setnames(dos, old = dos_sample_cols, new = samples)
    dos_sample_cols <- samples
  }

  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 20) next
  X <- as.matrix(dos[keep, ..dos_sample_cols])
  storage.mode(X) <- "double"
  colnames(X) <- dos_sample_cols

  # Identify homo and het samples
  het_labels <- c("Het", "HET", "G2")  # possible het labels
  homo_labels <- setdiff(unique(group), c(het_labels, "AMBIGUOUS", "ambiguous"))

  homo_idx <- which(group %in% homo_labels)
  het_idx <- which(group %in% het_labels)

  # VIEW A: all samples (reference — just record)
  view_a <- split_reanalysis(X, seq_len(ncol(X)), "all_samples")

  # VIEW B: homo-only
  if (length(homo_idx) >= 5) {
    homo_cols <- match(samples[homo_idx], dos_sample_cols)
    homo_cols <- homo_cols[!is.na(homo_cols)]
    view_b <- split_reanalysis(X, homo_cols, "homo_only")
    view_b[, homo1_n := sum(group[homo_idx] %in% c("Homo_1", "HOMO_1", "G1"))]
    view_b[, homo2_n := sum(group[homo_idx] %in% c("Homo_2", "HOMO_2", "G3", "G2"))]
  } else {
    view_b <- data.table(n_samples = length(homo_idx), label = "homo_only",
                          notes = "too_few_homos")
  }

  # VIEW C: het-only
  if (length(het_idx) >= 5) {
    het_cols <- match(samples[het_idx], dos_sample_cols)
    het_cols <- het_cols[!is.na(het_cols)]
    view_c <- split_reanalysis(X, het_cols, "het_only")
    view_c[, het_branch_hint := fifelse(
      !is.na(internal_cluster_count) & internal_cluster_count >= 2,
      "possible_branches", "single_cloud"
    )]
  } else {
    view_c <- data.table(n_samples = length(het_idx), label = "het_only",
                          notes = "too_few_hets")
  }

  split_summary <- rbindlist(list(view_a, view_b, view_c), fill = TRUE)
  split_summary[, candidate_id := cid]

  fwrite(split_summary,
         file.path(cand_dir, "candidate_split_reanalysis_summary.tsv"), sep = "\t")

  # ── Diagnostic plots ─────────────────────────────────────────────────────
  if (has_ggplot) {
    # PCA colored by branch_balance
    p1 <- ggplot(anchor_dt, aes(x = PC1, y = PC2, color = branch_balance)) +
      geom_point(size = 2, alpha = 0.7) +
      geom_point(data = anchor_dt[is_core == TRUE],
                 aes(x = PC1, y = PC2), shape = 1, size = 4, stroke = 0.8, color = "black") +
      scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#1B7837", midpoint = 0) +
      theme_bw(base_size = 11) +
      labs(title = paste0("Candidate ", cid, " — branch balance"),
           subtitle = "Circles = core anchor points",
           color = "Branch\nbalance")
    ggsave(file.path(plot_dir, "anchor_pca_branch_balance.png"),
           p1, width = 7, height = 5.5, dpi = 300)

    # u vs v
    p2 <- ggplot(anchor_dt, aes(x = u, y = v, color = coarse_group)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = GROUP_COLORS, na.value = "grey60") +
      theme_bw(base_size = 11) +
      labs(title = paste0("Candidate ", cid, " — anchor (u,v) space"),
           x = "u (broad stripe axis)", y = "v (perpendicular)")
    ggsave(file.path(plot_dir, "anchor_uv_space.png"),
           p2, width = 7, height = 5.5, dpi = 300)

    # center_distance vs branch_balance
    p3 <- ggplot(anchor_dt, aes(x = branch_balance, y = center_distance, color = coarse_group)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = GROUP_COLORS, na.value = "grey60") +
      theme_bw(base_size = 11) +
      labs(title = paste0("Candidate ", cid, " — branch balance vs center distance"),
           x = "Branch balance", y = "Center distance")
    ggsave(file.path(plot_dir, "anchor_balance_vs_center.png"),
           p3, width = 7, height = 5.5, dpi = 300)

    # PCA colored by ambiguity_score
    p4 <- ggplot(anchor_dt, aes(x = PC1, y = PC2, color = ambiguity_score)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_viridis_c(option = "inferno", direction = -1) +
      theme_bw(base_size = 11) +
      labs(title = paste0("Candidate ", cid, " — ambiguity score"),
           color = "Ambiguity")
    ggsave(file.path(plot_dir, "anchor_pca_ambiguity.png"),
           p4, width = 7, height = 5.5, dpi = 300)
  }

  message("[INFO] Candidate ", cid, " — anchor geometry + splits complete")
}

message("\n[DONE] STEP20b anchor geometry module complete")
