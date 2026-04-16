#!/usr/bin/env Rscript

# =============================================================================
# STEP22_candidate_subclustering.R
#
# Within-stripe density clustering + grid-based PCA shape encoding.
#
# For each candidate:
#   1. Load rotated PCA from STEP21
#   2. For each coarse stripe: run DBSCAN on v (+PC3 if available)
#   3. Classify stripe geometry: compact / split_discrete /
#      continuous_gradient / curved / diffuse
#   4. Overlay 8×8 grid on PCA space for shape quantification
#   5. Compute global + per-stripe grid features
#   6. Build candidate motif feature vector
#
# Outputs per candidate:
#   - candidate_subclusters.tsv
#   - candidate_stripe_geometry.tsv
#   - candidate_grid_features.tsv
#   - candidate_motif_vector.tsv
#
# Usage:
#   Rscript STEP22_candidate_subclustering.R <config.R> [candidate_id=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dbscan)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

# ── Read candidate table ─────────────────────────────────────────────────────
cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

# ── Stripe geometry classifier ───────────────────────────────────────────────
classify_stripe_geometry <- function(v_vals, pc3_vals = NULL) {
  n <- length(v_vals)
  if (n < 3) return("insufficient")

  v_range <- diff(range(v_vals, na.rm = TRUE))
  v_sd <- sd(v_vals, na.rm = TRUE)
  v_iqr <- IQR(v_vals, na.rm = TRUE)

  # Silverman's test for multimodality: check density peaks
  if (n >= 8) {
    dens <- density(v_vals, bw = "SJ", n = 512)
    # Count sign changes in derivative
    dy <- diff(dens$y)
    peaks <- sum(diff(sign(dy)) == -2)
  } else {
    peaks <- 1
  }

  # Curvature: if PC3 available, check for arc pattern
  curvature_score <- 0
  if (!is.null(pc3_vals) && sum(is.finite(pc3_vals)) >= 5) {
    r <- suppressWarnings(cor(v_vals, pc3_vals^2, use = "complete.obs"))
    if (is.finite(r)) curvature_score <- abs(r)
  }

  # Gradient: correlation of v with rank (sorted by something)
  gradient_score <- 0
  if (n >= 5) {
    # Check if v is spread smoothly
    sorted_v <- sort(v_vals)
    expected <- seq_along(sorted_v) / length(sorted_v)
    ks <- suppressWarnings(ks.test(v_vals, "punif",
                                   min(v_vals), max(v_vals))$statistic)
    gradient_score <- 1 - ks
  }

  # Compactness: coefficient of variation
  cv <- if (mean(abs(v_vals)) > 0) v_sd / mean(abs(v_vals)) else 0

  # Decision tree
  if (peaks >= 2 && v_range > 2 * v_iqr) {
    return("split_discrete")
  }
  if (curvature_score > 0.5) {
    return("curved")
  }
  if (gradient_score > 0.8 && peaks == 1 && v_range > v_sd * 2) {
    return("continuous_gradient")
  }
  if (cv < 0.5 && peaks <= 1) {
    return("compact")
  }

  return("diffuse")
}

# ── Grid-based PCA shape encoding ────────────────────────────────────────────
compute_grid_features <- function(x, y, groups, grid_size = 8) {
  # Overlay grid on PCA space
  x_range <- range(x, na.rm = TRUE)
  y_range <- range(y, na.rm = TRUE)

  # Pad ranges slightly
  x_pad <- diff(x_range) * 0.05
  y_pad <- diff(y_range) * 0.05
  x_breaks <- seq(x_range[1] - x_pad, x_range[2] + x_pad, length.out = grid_size + 1)
  y_breaks <- seq(y_range[1] - y_pad, y_range[2] + y_pad, length.out = grid_size + 1)

  # Assign cells
  x_cell <- findInterval(x, x_breaks, rightmost.closed = TRUE)
  y_cell <- findInterval(y, y_breaks, rightmost.closed = TRUE)
  x_cell <- pmin(pmax(x_cell, 1), grid_size)
  y_cell <- pmin(pmax(y_cell, 1), grid_size)

  cell_ids <- paste0(x_cell, "_", y_cell)

  # Global features
  occupied <- length(unique(cell_ids))
  total_cells <- grid_size^2
  frac_occupied <- occupied / total_cells

  # Occupancy counts per cell
  cell_counts <- table(cell_ids)
  n_total <- sum(cell_counts)

  # Entropy
  probs <- as.numeric(cell_counts) / n_total
  entropy <- -sum(probs * log2(probs + 1e-12))
  max_entropy <- log2(occupied)
  norm_entropy <- if (max_entropy > 0) entropy / max_entropy else 0

  # Concentration
  sorted_counts <- sort(as.numeric(cell_counts), decreasing = TRUE)
  top1_conc <- sorted_counts[1] / n_total
  top3_conc <- sum(sorted_counts[1:min(3, length(sorted_counts))]) / n_total

  # Connected components (simple 4-connected)
  occ_mat <- matrix(FALSE, grid_size, grid_size)
  for (cc in names(cell_counts)) {
    ij <- as.integer(strsplit(cc, "_")[[1]])
    occ_mat[ij[1], ij[2]] <- TRUE
  }
  n_components <- count_connected(occ_mat)

  # Per-stripe features
  stripe_features <- list()
  for (g in sort(unique(groups))) {
    idx <- groups == g
    if (sum(idx) < 2) next

    s_cells <- cell_ids[idx]
    s_counts <- table(s_cells)
    s_n <- sum(s_counts)
    s_occupied <- length(s_counts)

    # Bounding box in cell coordinates
    s_x <- x_cell[idx]
    s_y <- y_cell[idx]
    bb_width <- diff(range(s_x))
    bb_height <- diff(range(s_y))

    # Centroid in PCA space
    cx <- mean(x[idx])
    cy <- mean(y[idx])

    # Dispersion
    disp <- sqrt(mean((x[idx] - cx)^2 + (y[idx] - cy)^2))

    # Stripe entropy
    s_probs <- as.numeric(s_counts) / s_n
    s_entropy <- -sum(s_probs * log2(s_probs + 1e-12))

    stripe_features[[as.character(g)]] <- data.table(
      group = g,
      n_samples = sum(idx),
      occupied_cells = s_occupied,
      occupied_rows = length(unique(s_y)),
      occupied_cols = length(unique(s_x)),
      frac_occupied = round(s_occupied / total_cells, 4),
      bb_width = bb_width,
      bb_height = bb_height,
      centroid_x = round(cx, 4),
      centroid_y = round(cy, 4),
      dispersion = round(disp, 4),
      entropy = round(s_entropy, 4),
      top_cell_conc = round(max(s_counts) / s_n, 4)
    )
  }

  list(
    global = data.table(
      grid_size = grid_size,
      total_occupied = occupied,
      frac_occupied = round(frac_occupied, 4),
      occupancy_entropy = round(entropy, 4),
      norm_entropy = round(norm_entropy, 4),
      top1_concentration = round(top1_conc, 4),
      top3_concentration = round(top3_conc, 4),
      n_connected_components = n_components
    ),
    per_stripe = if (length(stripe_features) > 0) rbindlist(stripe_features) else NULL,
    cell_assignments = data.table(x_cell = x_cell, y_cell = y_cell, cell_id = cell_ids)
  )
}

# Simple connected-component count for 4-connected grid
count_connected <- function(mat) {
  visited <- matrix(FALSE, nrow(mat), ncol(mat))
  nr <- nrow(mat); nc <- ncol(mat)
  count <- 0L

  bfs <- function(si, sj) {
    queue <- list(c(si, sj))
    visited[si, sj] <<- TRUE
    while (length(queue) > 0) {
      pos <- queue[[1]]; queue <- queue[-1]
      i <- pos[1]; j <- pos[2]
      for (d in list(c(-1,0), c(1,0), c(0,-1), c(0,1))) {
        ni <- i + d[1]; nj <- j + d[2]
        if (ni >= 1 && ni <= nr && nj >= 1 && nj <= nc &&
            !visited[ni, nj] && mat[ni, nj]) {
          visited[ni, nj] <<- TRUE
          queue <- c(queue, list(c(ni, nj)))
        }
      }
    }
  }

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (mat[i, j] && !visited[i, j]) {
        bfs(i, j)
        count <- count + 1L
      }
    }
  }
  count
}

# ── Pairwise stripe features ────────────────────────────────────────────────
compute_pairwise_stripe <- function(u, v, groups) {
  g_levels <- sort(unique(groups))
  if (length(g_levels) < 2) return(NULL)

  pairs <- list()
  for (i in seq_along(g_levels)) {
    for (j in seq(i + 1, length(g_levels))) {
      if (j > length(g_levels)) break
      g1 <- g_levels[i]; g2 <- g_levels[j]
      idx1 <- groups == g1; idx2 <- groups == g2

      c1 <- c(mean(u[idx1]), mean(v[idx1]))
      c2 <- c(mean(u[idx2]), mean(v[idx2]))
      cdist <- sqrt(sum((c1 - c2)^2))

      # Gap along u
      u_gap <- min(u[idx2]) - max(u[idx1])

      pairs[[length(pairs) + 1]] <- data.table(
        group_1 = g1, group_2 = g2,
        centroid_distance = round(cdist, 4),
        u_gap = round(u_gap, 4),
        n_group_1 = sum(idx1),
        n_group_2 = sum(idx2)
      )
    }
  }
  rbindlist(pairs)
}

# ── Main processing loop ────────────────────────────────────────────────────
all_subclusters <- list()
all_geometries <- list()
all_grids <- list()
all_motifs <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom

  # Load STEP21 rotated PCA
  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (!file.exists(rot_file)) {
    message("[SKIP] No STEP21 output for candidate ", cid)
    next
  }

  rot <- fread(rot_file)
  n <- nrow(rot)
  if (n < 5) next

  message("[INFO] Candidate ", cid, " (", chr, "): n=", n)

  # ── Within-stripe DBSCAN ──────────────────────────────────────────────
  sub_results <- list()
  geom_results <- list()

  for (g in sort(unique(rot$coarse_group_refined))) {
    idx <- which(rot$coarse_group_refined == g)
    ng <- length(idx)

    if (ng < 3) {
      sub_results[[g]] <- data.table(
        candidate_id = cid, sample = rot$sample[idx],
        coarse_group = g, subcluster_label = paste0(g, "_sub1"),
        subcluster_method = "none", cluster_confidence = NA_real_,
        is_noise = FALSE, u = rot$u[idx], v = rot$v[idx]
      )
      geom_results[[g]] <- data.table(
        candidate_id = cid, coarse_group = g,
        geometry_label = "insufficient", n_samples = ng,
        n_subclusters = 1
      )
      next
    }

    # Build feature matrix for DBSCAN: v + PC3
    v_vals <- rot$v[idx]
    feats <- matrix(v_vals, ncol = 1)
    if ("PC3" %in% names(rot) && sum(is.finite(rot$PC3[idx])) >= ng * 0.8) {
      feats <- cbind(feats, rot$PC3[idx])
      # Scale features
      feats <- scale(feats)
      feats[!is.finite(feats)] <- 0
    } else {
      feats <- scale(feats)
      feats[!is.finite(feats)] <- 0
    }

# Auto-tune epsilon from k-nearest-neighbor distances
k_nn <- min(5, ng - 1)
knn_dists <- dbscan::kNNdist(feats, k = k_nn)

# kNNdist() may return either:
# - a matrix (multi-dimensional features)
# - a vector (1D feature case)
knn_last <- if (is.matrix(knn_dists) || is.data.frame(knn_dists)) {
  knn_dists[, k_nn]
} else {
  knn_dists
}

knn_last <- knn_last[is.finite(knn_last)]
if (length(knn_last) == 0) {
  eps_val <- 0.5
} else {
  eps_val <- median(knn_last) * 0.8
  if (!is.finite(eps_val) || eps_val <= 0) eps_val <- 0.5
}


    # Run DBSCAN
    db <- dbscan::dbscan(feats, eps = eps_val, minPts = max(3, floor(ng * 0.08)))

    # Label subclusters
    n_clust <- max(db$cluster)
    if (n_clust == 0) n_clust <- 1  # all noise → single cluster

    sub_labels <- ifelse(
      db$cluster == 0,
      paste0(g, "_noise"),
      paste0(g, "_sub", db$cluster)
    )

knn_first <- if (is.matrix(knn_dists) || is.data.frame(knn_dists)) {
  knn_dists[, 1]
} else {
  knn_dists
}

knn_first <- ifelse(is.finite(knn_first), knn_first, eps_val)

cluster_conf <- ifelse(db$cluster == 0, 0, 1 - pmin(knn_first / eps_val, 1))


    sub_results[[g]] <- data.table(
      candidate_id = cid,
      sample = rot$sample[idx],
      coarse_group = g,
      subcluster_label = sub_labels,
      subcluster_method = "DBSCAN",
      cluster_confidence = round(cluster_conf, 4),
      is_noise = db$cluster == 0,
      u = round(rot$u[idx], 4),
      v = round(rot$v[idx], 4)
    )

    # Stripe geometry
    pc3_vals <- if ("PC3" %in% names(rot)) rot$PC3[idx] else NULL
    geom_label <- classify_stripe_geometry(v_vals, pc3_vals)
    geom_results[[g]] <- data.table(
      candidate_id = cid,
      coarse_group = g,
      geometry_label = geom_label,
      n_samples = ng,
      n_subclusters = max(1, n_clust),
      eps_used = round(eps_val, 4),
      minPts_used = max(3, floor(ng * 0.08))
    )
  }

  sub_dt <- rbindlist(sub_results, fill = TRUE)
  geom_dt <- rbindlist(geom_results, fill = TRUE)

  # ── Grid-based PCA shape encoding ──────────────────────────────────────
  grid <- compute_grid_features(rot$u, rot$v, rot$coarse_group_refined, grid_size = 8)
  grid_global <- grid$global
  grid_global[, candidate_id := cid]

  # Pairwise stripe features
  pw <- compute_pairwise_stripe(rot$u, rot$v, rot$coarse_group_refined)
  if (!is.null(pw)) pw[, candidate_id := cid]

  # ── Motif vector ───────────────────────────────────────────────────────
#  subgroup_counts <- sub_dt[!is_noise, .(n_sub = uniqueN(subcluster_label)), by = coarse_group]
subgroup_counts <- sub_dt[!is.na(is_noise) & is_noise == FALSE,
                          .(n_sub = uniqueN(subcluster_label)),
                          by = coarse_group]

  homo1_sub <- subgroup_counts[coarse_group == "HOMO_1", n_sub]
  het_sub <- subgroup_counts[coarse_group == "HET", n_sub]
  homo2_sub <- subgroup_counts[coarse_group == "HOMO_2", n_sub]

  motif_dt <- data.table(
    candidate_id = cid,
    chrom = chr,
    HOMO_1_subgroups = if (length(homo1_sub) > 0) homo1_sub else 1L,
    HET_subgroups = if (length(het_sub) > 0) het_sub else 1L,
    HOMO_2_subgroups = if (length(homo2_sub) > 0) homo2_sub else 1L,
    grid_occupied = grid_global$total_occupied,
    grid_frac_occupied = grid_global$frac_occupied,
    grid_entropy = grid_global$occupancy_entropy,
    grid_top1_conc = grid_global$top1_concentration,
    grid_components = grid_global$n_connected_components,
    dominant_geometry = paste(geom_dt$geometry_label, collapse = "/")
  )

  # ── Save per-candidate ─────────────────────────────────────────────────
  fwrite(sub_dt, file.path(cand_dir, "candidate_subclusters.tsv"), sep = "\t")
  fwrite(geom_dt, file.path(cand_dir, "candidate_stripe_geometry.tsv"), sep = "\t")
  fwrite(grid_global, file.path(cand_dir, "candidate_grid_features.tsv"), sep = "\t")
  if (!is.null(grid$per_stripe))
    fwrite(grid$per_stripe, file.path(cand_dir, "candidate_grid_per_stripe.tsv"), sep = "\t")
  if (!is.null(pw))
    fwrite(pw, file.path(cand_dir, "candidate_pairwise_stripe.tsv"), sep = "\t")
  fwrite(motif_dt, file.path(cand_dir, "candidate_motif_vector.tsv"), sep = "\t")

  all_subclusters[[length(all_subclusters) + 1]] <- sub_dt
  all_geometries[[length(all_geometries) + 1]] <- geom_dt
  all_grids[[length(all_grids) + 1]] <- grid_global
  all_motifs[[length(all_motifs) + 1]] <- motif_dt

  message("[INFO]   Subclusters: ",
          paste(geom_dt$coarse_group, geom_dt$n_subclusters, geom_dt$geometry_label,
                sep = "=", collapse = " | "))
}

# ── Global outputs ───────────────────────────────────────────────────────────
if (length(all_subclusters) > 0) {
  fwrite(rbindlist(all_subclusters, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_subclusters.tsv.gz"), sep = "\t")
  fwrite(rbindlist(all_geometries, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_stripe_geometry.tsv.gz"), sep = "\t")
  fwrite(rbindlist(all_grids, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_grid_features.tsv.gz"), sep = "\t")
  fwrite(rbindlist(all_motifs, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_motif_vectors.tsv.gz"), sep = "\t")
}

message("[DONE] STEP22 complete")
