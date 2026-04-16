#!/usr/bin/env Rscript

# =============================================================================
# STEP21_candidate_state_assignment.R
#
# Stripe-aware coarse state assignment + rotation refinement.
#
# For each candidate:
#   A. Coarse K=3 grouping on PC1 (GMM preferred, k-means fallback)
#   B. Order groups by PC1 centroid → HOMO_1 / HET / HOMO_2
#   C. Fit broad stripe axis through centroids → rotation to u / v
#   D. Optionally refine grouping on u instead of PC1
#   E. Assign confidence and AMBIGUOUS flags
#
# Inputs:
#   - STEP12 regional_pca_samples tables
#   - Candidate table
#
# Outputs per candidate (in FOLLOWUP_DIR/<chrom>.candidate_<id>/):
#   - candidate_pca_rotated.tsv
#
# Global output:
#   - all_candidates_pca_rotated.tsv.gz
#
# Usage:
#   Rscript STEP21_candidate_state_assignment.R <config.R> [candidate_id=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

# ── Try loading mclust for GMM; fall back to k-means ────────────────────────
use_gmm <- suppressWarnings(require(mclust, quietly = TRUE))
if (use_gmm) {
  message("[INFO] Using mclust GMM for coarse grouping")
} else {
  message("[INFO] mclust not available — using k-means for coarse grouping")
}

# ── Read candidate table ─────────────────────────────────────────────────────
cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]
message("[INFO] Processing ", nrow(cand), " candidates")

# ── Find STEP12 regional PCA files ──────────────────────────────────────────
pca_files <- list.files(STEP12_DIR,
                        pattern = "\\.regional_pca_samples\\.tsv\\.gz$",
                        full.names = TRUE)
message("[INFO] Found ", length(pca_files), " STEP12 PCA files")

# Index by candidate_id
pca_index <- data.table(
  file = pca_files,
  basename = basename(pca_files)
)
pca_index[, candidate_id := as.integer(
  sub(".*candidate_([0-9]+)\\.regional.*", "\\1", basename)
)]

# ── Rotation helper ──────────────────────────────────────────────────────────
rotate_to_stripe_axis <- function(pc1, pc2, group) {
  # Compute centroids per group
  centroids <- data.table(group = group, pc1 = pc1, pc2 = pc2)[
    , .(c1 = mean(pc1), c2 = mean(pc2)), by = group
  ][order(c1)]

  if (nrow(centroids) < 2) {
    return(list(u = pc1, v = pc2, angle = 0))
  }

  # Fit line through centroids (endpoints or all)
  if (nrow(centroids) >= 3) {
    # Use extreme centroids
    dx <- centroids$c1[nrow(centroids)] - centroids$c1[1]
    dy <- centroids$c2[nrow(centroids)] - centroids$c2[1]
  } else {
    dx <- centroids$c1[2] - centroids$c1[1]
    dy <- centroids$c2[2] - centroids$c2[1]
  }

  angle <- atan2(dy, dx)

  # Rotation matrix
  cos_a <- cos(-angle)
  sin_a <- sin(-angle)
  u <-  pc1 * cos_a + pc2 * sin_a
  v <- -pc1 * sin_a + pc2 * cos_a

  list(u = u, v = v, angle = angle * 180 / pi)
}

# ── Confidence scorer ────────────────────────────────────────────────────────
compute_group_confidence <- function(pc1, group, method = "distance") {
  centroids <- tapply(pc1, group, mean, na.rm = TRUE)
  centroids <- sort(centroids)

  if (length(centroids) < 3) return(rep(0.5, length(pc1)))

  # For each sample: distance to assigned centroid / distance to nearest other
  conf <- numeric(length(pc1))
  for (i in seq_along(pc1)) {
    g <- group[i]
    own_center <- centroids[as.character(g)]
    other_centers <- centroids[names(centroids) != as.character(g)]
    d_own <- abs(pc1[i] - own_center)
    d_nearest <- min(abs(pc1[i] - other_centers))
    conf[i] <- if ((d_own + d_nearest) > 0) d_nearest / (d_own + d_nearest) else 0.5
  }
  conf
}

# ── Main processing loop ────────────────────────────────────────────────────
all_rotated <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  # Find STEP12 file
  pca_row <- pca_index[candidate_id == cid]
  if (nrow(pca_row) == 0) {
    message("[SKIP] No STEP12 file for candidate ", cid)
    next
  }

  pcs <- fread(pca_row$file[1])
  if (!all(c("sample", "PC1", "PC2") %in% names(pcs))) {
    message("[SKIP] Missing PC columns for candidate ", cid)
    next
  }

  n <- nrow(pcs)
  if (n < 5) {
    message("[SKIP] Too few samples (", n, ") for candidate ", cid)
    next
  }

  # ── A. Coarse grouping on PC1 ──────────────────────────────────────────
  pc1_vec <- pcs$PC1
  n_groups <- 3L

  if (use_gmm && n >= 10) {
    mc <- tryCatch(
      mclust::Mclust(pc1_vec, G = 3, modelNames = "V", verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(mc)) {
      raw_group <- mc$classification
      posteriors <- apply(mc$z, 1, max)
    } else {
      set.seed(42)
      km <- kmeans(matrix(pc1_vec, ncol = 1), centers = 3, nstart = 50)
      raw_group <- km$cluster
      posteriors <- rep(NA_real_, n)
    }
  } else {
    set.seed(42)
    km <- kmeans(matrix(pc1_vec, ncol = 1), centers = min(3, n), nstart = 50)
    raw_group <- km$cluster
    posteriors <- rep(NA_real_, n)
    n_groups <- length(unique(km$cluster))
  }

  # Order by PC1 centroid
  grp_means <- tapply(pc1_vec, raw_group, mean)
  ord <- order(grp_means)
  ordered_group <- match(raw_group, ord)

  # Labels
  if (n_groups == 3) {
    coarse_labels <- c("HOMO_1", "HET", "HOMO_2")
  } else {
    coarse_labels <- paste0("G", seq_len(n_groups))
  }
  coarse_group <- coarse_labels[ordered_group]

  # ── B. Rotation to u / v ───────────────────────────────────────────────
  rot <- rotate_to_stripe_axis(pcs$PC1, pcs$PC2, ordered_group)

  # ── C. Optional refinement on u ────────────────────────────────────────
  # Re-run k-means on u to see if grouping improves
  set.seed(42)
  km_u <- tryCatch(
    kmeans(matrix(rot$u, ncol = 1), centers = min(3, n), nstart = 50),
    error = function(e) NULL
  )

  if (!is.null(km_u)) {
    grp_means_u <- tapply(rot$u, km_u$cluster, mean)
    ord_u <- order(grp_means_u)
    refined_group_idx <- match(km_u$cluster, ord_u)
    refined_group <- coarse_labels[refined_group_idx]
  } else {
    refined_group <- coarse_group
  }

  # ── D. Confidence and ambiguity ────────────────────────────────────────
  confidence <- compute_group_confidence(rot$u, refined_group)
  ambiguous <- confidence < 0.4

  # ── E. Add PC3 if available ────────────────────────────────────────────
  pc3 <- if ("PC3" %in% names(pcs)) pcs$PC3 else rep(NA_real_, n)

  # ── F. Preserve existing regional_het if present ───────────────────────
  reg_het <- if ("regional_het" %in% names(pcs)) pcs$regional_het else rep(NA_real_, n)
  hap_score <- if ("regional_hap_score" %in% names(pcs)) pcs$regional_hap_score else rep(NA_real_, n)

  # ── Build output table ─────────────────────────────────────────────────
  out_dt <- data.table(
    candidate_id       = cid,
    chrom              = chr,
    start_bp           = c_start,
    end_bp             = c_end,
    sample             = pcs$sample,
    PC1                = round(pcs$PC1, 4),
    PC2                = round(pcs$PC2, 4),
    PC3                = round(pc3, 4),
    u                  = round(rot$u, 4),
    v                  = round(rot$v, 4),
    rotation_angle_deg = round(rot$angle, 2),
    coarse_group_initial = coarse_group,
    coarse_group_refined = refined_group,
    coarse_confidence  = round(confidence, 4),
    ambiguous_flag     = ambiguous,
    regional_het       = round(reg_het, 6),
    regional_hap_score = round(hap_score, 6)
  )

  # Save per-candidate
  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(cand_dir)
  fwrite(out_dt, file.path(cand_dir, "candidate_pca_rotated.tsv"), sep = "\t")

  all_rotated[[length(all_rotated) + 1]] <- out_dt

  # Summary
  grp_tab <- table(refined_group)
  grp_str <- paste(names(grp_tab), grp_tab, sep = "=", collapse = " ")
  message("[INFO] Candidate ", cid, " (", chr, "): ", grp_str,
          " | rotation=", round(rot$angle, 1), "° | ambiguous=", sum(ambiguous))
}

# ── Global output ────────────────────────────────────────────────────────────
if (length(all_rotated) > 0) {
  global_dt <- rbindlist(all_rotated, fill = TRUE)
  fwrite(global_dt,
         file.path(FOLLOWUP_DIR, "all_candidates_pca_rotated.tsv.gz"),
         sep = "\t")
  message("[DONE] Global rotated PCA table: ", nrow(global_dt), " records across ",
          uniqueN(global_dt$candidate_id), " candidates")
}

message("[DONE] STEP21 complete")
