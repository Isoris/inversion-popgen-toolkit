#!/usr/bin/env Rscript

# =============================================================================
# STEP33_candidate_within_stripe_analysis.R
#
# Dedicated within-stripe analyses: stripe-3-only, stripe-2-only, and
# stripe-2 vs stripe-3 comparison.
#
# For each stripe:
#   1. PCA on stripe samples only
#   2. Geometry metrics (centroid, covariance, elongation, curvature)
#   3. Quality-aware subgroup detection using coherence tiers
#   4. Cross-stripe comparison
#   5. AB asymmetry analysis (stripe 2 vs stripe 1/3 separately)
#
# Inputs:
#   - STEP21 candidate_pca_rotated.tsv
#   - STEP29 candidate_sample_coherence.tsv
#   - Dosage matrix
#
# Outputs per candidate:
#   - candidate_stripe3_pca.tsv
#   - candidate_stripe2_pca.tsv
#   - candidate_stripe_geometry_detailed.tsv
#   - candidate_stripe_comparison.tsv
#   - candidate_ab_asymmetry.tsv
#
# Usage:
#   Rscript STEP33_candidate_within_stripe_analysis.R <config.R> [cid=all]
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

# ── Stripe geometry metrics ──────────────────────────────────────────────────
compute_stripe_geometry <- function(pc1, pc2) {
  n <- length(pc1)
  if (n < 3) return(data.table(n = n, computable = FALSE))

  cx <- mean(pc1, na.rm = TRUE)
  cy <- mean(pc2, na.rm = TRUE)

  dists <- sqrt((pc1 - cx)^2 + (pc2 - cy)^2)
  med_dist <- median(dists)
  mad_dist <- mad(dists)
  mean_dist <- mean(dists)

  # Covariance eigenvalues for elongation
  if (n >= 5) {
    covmat <- cov(cbind(pc1, pc2), use = "pairwise.complete.obs")
    evals <- eigen(covmat, symmetric = TRUE)$values
    elongation <- if (evals[2] > 0) evals[1] / evals[2] else Inf
    principal_angle <- atan2(eigen(covmat)$vectors[2, 1],
                              eigen(covmat)$vectors[1, 1]) * 180 / pi
  } else {
    evals <- c(NA, NA)
    elongation <- NA_real_
    principal_angle <- NA_real_
  }

  # Density concentration: fraction within 1 MAD of centroid
  conc_1mad <- mean(dists <= med_dist + mad_dist)

  # Multimodality: Hartigan's dip test on PC1
  dip_p <- tryCatch({
    dt <- diptest::dip.test(pc1)
    dt$p.value
  }, error = function(e) NA_real_)

  data.table(
    n = n,
    computable = TRUE,
    centroid_pc1 = round(cx, 4),
    centroid_pc2 = round(cy, 4),
    median_dist_to_centroid = round(med_dist, 4),
    mad_dist_to_centroid = round(mad_dist, 4),
    mean_dist_to_centroid = round(mean_dist, 4),
    eigenvalue_1 = round(evals[1], 4),
    eigenvalue_2 = round(evals[2], 4),
    elongation_ratio = round(elongation, 4),
    principal_angle_deg = round(principal_angle, 2),
    concentration_1mad = round(conc_1mad, 4),
    dip_test_pvalue = round(dip_p, 6)
  )
}

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (!file.exists(rot_file)) next

  rot <- fread(rot_file)
  if (nrow(rot) < 10) next

  # Load coherence if available
  coh <- tryCatch(fread(file.path(cand_dir, "candidate_sample_coherence.tsv")),
                  error = function(e) NULL)
  if (!is.null(coh)) {
    rot <- merge(rot, coh[, .(sample, coherence_class, stripe_quality)],
                 by = "sample", all.x = TRUE)
  }

  message("[INFO] Candidate ", cid, " (", chr, "): within-stripe analysis")

  # Load dosage for stripe-specific PCA
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  has_dosage <- file.exists(dos_file) && file.exists(sites_file)

  dos_reg <- NULL; X <- NULL
  if (has_dosage) {
    dos <- fread(dos_file)
    sites <- fread(sites_file)
    keep <- which(sites$pos >= c_start & sites$pos <= c_end)
    if (length(keep) >= 20) {
      dos_reg <- dos[keep]
      # Handle Ind-style columns
      sc <- setdiff(names(dos_reg), "marker")
      if (all(grepl("^Ind", sc)) && length(sc) == nrow(rot)) {
        setnames(dos_reg, old = sc, new = rot$sample)
      }
      common <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))
      if (length(common) >= 10) {
        X <- as.matrix(dos_reg[, ..common])
        storage.mode(X) <- "double"
      }
    }
  }

  # ── Per-stripe analysis ────────────────────────────────────────────────
  geom_results <- list()

  for (g in c("HOMO_1", "HET", "HOMO_2")) {
    stripe <- rot[coarse_group_refined == g]
    ng <- nrow(stripe)
    if (ng < 3) next

    # Geometry on existing PCA coords
    geom <- compute_stripe_geometry(stripe$PC1, stripe$PC2)
    geom[, stripe := g]
    geom[, candidate_id := cid]

    # Also compute on u/v
    geom_uv <- compute_stripe_geometry(stripe$u, stripe$v)
    setnames(geom_uv, names(geom_uv),
             paste0(names(geom_uv), "_uv"))
    geom <- cbind(geom, geom_uv)

    # Quality tier distribution
    if ("stripe_quality" %in% names(stripe)) {
      tier_tab <- table(stripe$stripe_quality)
      geom[, n_core := as.integer(tier_tab["core"])]
      geom[, n_peripheral := as.integer(tier_tab["peripheral"])]
      geom[, n_junk := as.integer(tier_tab["junk"])]
      for (tc in c("n_core", "n_peripheral", "n_junk")) {
        if (is.na(geom[[tc]])) geom[[tc]] <- 0L
      }
    }

    geom_results[[g]] <- geom

    # Stripe-specific PCA from dosage
    if (!is.null(X) && ng >= 5) {
      stripe_cols <- intersect(stripe$sample, colnames(X))
      if (length(stripe_cols) >= 5) {
        X_stripe <- X[, stripe_cols, drop = FALSE]
        # Select variable markers
        v_var <- apply(X_stripe, 1, var, na.rm = TRUE)
        top_m <- order(v_var, decreasing = TRUE)[1:min(nrow(X_stripe), 3000)]

        pca_s <- tryCatch(
          prcomp(t(X_stripe[top_m, ]), center = TRUE, scale. = FALSE),
          error = function(e) NULL
        )

        if (!is.null(pca_s)) {
          stripe_pca_dt <- data.table(
            candidate_id = cid,
            sample = stripe_cols,
            stripe = g,
            stripe_PC1 = round(pca_s$x[, 1], 4),
            stripe_PC2 = round(pca_s$x[, 2], 4)
          )
          if ("stripe_quality" %in% names(rot)) {
            stripe_pca_dt <- merge(stripe_pca_dt,
                                   rot[, .(sample, stripe_quality)],
                                   by = "sample", all.x = TRUE)
          }

          outname <- paste0("candidate_stripe",
                           ifelse(g == "HOMO_1", "1",
                                  ifelse(g == "HET", "2", "3")),
                           "_pca.tsv")
          fwrite(stripe_pca_dt, file.path(cand_dir, outname), sep = "\t")
        }
      }
    }
  }

  if (length(geom_results) > 0) {
    geom_dt <- rbindlist(geom_results, fill = TRUE)
    fwrite(geom_dt, file.path(cand_dir, "candidate_stripe_geometry_detailed.tsv"),
           sep = "\t")
  }

  # ── AB asymmetry analysis ──────────────────────────────────────────────
  # For each HET sample, compute similarity to HOMO_1 vs HOMO_2
  if (!is.null(X)) {
    het_samples <- rot[coarse_group_refined == "HET", sample]
    h1_samples <- rot[coarse_group_refined == "HOMO_1", sample]
    h2_samples <- rot[coarse_group_refined == "HOMO_2", sample]

    het_cols <- intersect(het_samples, colnames(X))
    h1_cols <- intersect(h1_samples, colnames(X))
    h2_cols <- intersect(h2_samples, colnames(X))

    if (length(het_cols) >= 3 && length(h1_cols) >= 3 && length(h2_cols) >= 3) {
      # Use informative markers
      h1_mean <- rowMeans(X[, h1_cols, drop = FALSE], na.rm = TRUE)
      h2_mean <- rowMeans(X[, h2_cols, drop = FALSE], na.rm = TRUE)
      delta <- abs(h2_mean - h1_mean)
      info_idx <- which(delta >= quantile(delta, 0.75, na.rm = TRUE))

      asym_list <- list()
      for (sid in het_cols) {
        x_s <- X[info_idx, sid]
        d_h1 <- mean(abs(x_s - h1_mean[info_idx]), na.rm = TRUE)
        d_h2 <- mean(abs(x_s - h2_mean[info_idx]), na.rm = TRUE)

        # AA-side similarity (lower distance = more similar)
        aa_sim <- if ((d_h1 + d_h2) > 0) 1 - d_h1 / (d_h1 + d_h2) else 0.5
        bb_sim <- if ((d_h1 + d_h2) > 0) 1 - d_h2 / (d_h1 + d_h2) else 0.5

        asym_list[[length(asym_list) + 1]] <- data.table(
          candidate_id = cid,
          sample = sid,
          dist_to_HOMO_1 = round(d_h1, 4),
          dist_to_HOMO_2 = round(d_h2, 4),
          aa_side_similarity = round(aa_sim, 4),
          bb_side_similarity = round(bb_sim, 4),
          ab_balance = round(aa_sim - bb_sim, 4)  # positive = more AA-like
        )
      }

      asym_dt <- rbindlist(asym_list)
      fwrite(asym_dt, file.path(cand_dir, "candidate_ab_asymmetry.tsv"), sep = "\t")

      msg_bal <- round(mean(abs(asym_dt$ab_balance)), 3)
      message("[INFO]   AB asymmetry: mean |balance| = ", msg_bal,
              " (0 = symmetric, 1 = completely one-sided)")
    }
  }

  # ── Cross-stripe comparison ────────────────────────────────────────────
  if (length(geom_results) >= 2) {
    pairs <- list()
    gs <- names(geom_results)
    for (i in seq_along(gs)) {
      for (j in seq(i + 1, length(gs))) {
        if (j > length(gs)) break
        g1 <- geom_results[[gs[i]]]
        g2 <- geom_results[[gs[j]]]
        if (!g1$computable || !g2$computable) next

        cdist <- sqrt((g1$centroid_pc1 - g2$centroid_pc1)^2 +
                       (g1$centroid_pc2 - g2$centroid_pc2)^2)
        pairs[[length(pairs) + 1]] <- data.table(
          candidate_id = cid,
          stripe_1 = gs[i], stripe_2 = gs[j],
          centroid_distance = round(cdist, 4),
          n_stripe_1 = g1$n, n_stripe_2 = g2$n,
          elongation_1 = g1$elongation_ratio,
          elongation_2 = g2$elongation_ratio
        )
      }
    }
    if (length(pairs) > 0) {
      comp_dt <- rbindlist(pairs)
      fwrite(comp_dt, file.path(cand_dir, "candidate_stripe_comparison.tsv"), sep = "\t")
    }
  }

  message("[INFO]   Within-stripe analysis done")
}

message("[DONE] STEP33 complete")
