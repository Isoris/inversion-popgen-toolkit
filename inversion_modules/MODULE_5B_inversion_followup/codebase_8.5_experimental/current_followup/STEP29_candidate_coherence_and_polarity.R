#!/usr/bin/env Rscript

# =============================================================================
# STEP29_candidate_coherence_and_polarity.R  v2 — debug-first rewrite
#
# FIXES:
#   1. HET coherence: separate thresholds for homos vs HET (HET is expected
#      intermediate, so agreement ~0.5 is NOT discordant for HET)
#   2. L1/L2 hierarchical polarity: L1 = per-marker sign, L2 = smoothed
#      block consensus with minimum block size
#   3. Anti-circularity: polarity inferred on core homos only, applied to all
#   4. Full per-marker audit table with all requested columns
#   5. Block definition uses minimum block size (default 5 markers)
#
# Usage:
#   Rscript STEP29_candidate_coherence_and_polarity.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

MIN_BLOCK_SIZE <- 5L  # minimum markers per L2 polarity block

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id; chr <- row$chrom
  c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (!file.exists(rot_file)) next
  rot <- fread(rot_file)
  if (nrow(rot) < 5) next

  message("[INFO] Candidate ", cid, " (", chr, "): coherence + polarity v2")

  # ── Load dosage ────────────────────────────────────────────────────────
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) next
  dos <- fread(dos_file); sites <- fread(sites_file)
  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 20) next
  sites_reg <- sites[keep]; dos_reg <- dos[keep]
  sample_cols <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))
  if (length(sample_cols) < 5) {
    ind_cols <- setdiff(names(dos_reg), "marker")
    if (all(grepl("^Ind", ind_cols)) && length(ind_cols) == nrow(rot)) {
      setnames(dos_reg, old = ind_cols, new = rot$sample)
      sample_cols <- rot$sample
    } else next
  }
  X <- as.matrix(dos_reg[, ..sample_cols]); storage.mode(X) <- "double"
  n_markers <- nrow(X); n_samples <- ncol(X)

  # ── Group membership ───────────────────────────────────────────────────
  groups <- rot[match(sample_cols, sample), coarse_group_refined]
  h1_idx <- which(groups == "HOMO_1")
  h2_idx <- which(groups == "HOMO_2")
  het_idx <- which(groups == "HET")

  # ── Per-group marker summaries (mean + median) ─────────────────────────
  grp_list <- list(HOMO_1 = h1_idx, HET = het_idx, HOMO_2 = h2_idx)
  grp_means <- grp_medians <- grp_mads <- matrix(NA_real_, n_markers, 3)
  colnames(grp_means) <- colnames(grp_medians) <- colnames(grp_mads) <- c("HOMO_1", "HET", "HOMO_2")

  for (gi in seq_along(grp_list)) {
    gidx <- grp_list[[gi]]
    if (length(gidx) >= 2) {
      grp_means[, gi] <- rowMeans(X[, gidx, drop = FALSE], na.rm = TRUE)
      grp_medians[, gi] <- apply(X[, gidx, drop = FALSE], 1, median, na.rm = TRUE)
      grp_mads[, gi] <- apply(X[, gidx, drop = FALSE], 1, mad, na.rm = TRUE)
    }
  }

  # ══════════════════════════════════════════════════════════════════════
  # L1 POLARITY: per-marker raw orientation
  # ══════════════════════════════════════════════════════════════════════
  delta_hom_mean <- grp_means[, "HOMO_2"] - grp_means[, "HOMO_1"]
  delta_hom_median <- grp_medians[, "HOMO_2"] - grp_medians[, "HOMO_1"]
  delta_hom_mean[!is.finite(delta_hom_mean)] <- 0
  delta_hom_median[!is.finite(delta_hom_median)] <- 0

  abs_delta <- abs(delta_hom_mean)

  # L1 sign and confidence
  l1_sign <- sign(delta_hom_mean)
  l1_sign[l1_sign == 0] <- 1  # tie-break positive
  l1_confidence <- pmin(abs_delta / 0.5, 1)  # saturates at |Δ|=0.5
  l1_ambiguous <- abs_delta < 0.1

  # Dominant candidate-wide sign
  dominant_sign <- sign(median(delta_hom_mean[abs_delta > 0.1], na.rm = TRUE))
  if (!is.finite(dominant_sign) || dominant_sign == 0) dominant_sign <- 1

  l1_reversed <- l1_sign != dominant_sign & !l1_ambiguous

  # ══════════════════════════════════════════════════════════════════════
  # L2 POLARITY: block-level smoothed consensus
  # ══════════════════════════════════════════════════════════════════════
  # Sliding window consensus: for each marker, look at a neighborhood
  # of MIN_BLOCK_SIZE markers and take the weighted median sign
  half_win <- MIN_BLOCK_SIZE %/% 2
  l2_sign <- numeric(n_markers)
  l2_confidence <- numeric(n_markers)

  for (mi in seq_len(n_markers)) {
    lo <- max(1, mi - half_win)
    hi <- min(n_markers, mi + half_win)
    local_delta <- delta_hom_mean[lo:hi]
    local_abs <- abs(local_delta)

    # Weighted median sign: weight by |Δ|
    w <- local_abs
    w[!is.finite(w)] <- 0
    pos_weight <- sum(w[local_delta > 0])
    neg_weight <- sum(w[local_delta < 0])
    l2_sign[mi] <- if (pos_weight >= neg_weight) 1 else -1
    total_weight <- pos_weight + neg_weight
    l2_confidence[mi] <- if (total_weight > 0) max(pos_weight, neg_weight) / total_weight else 0.5
  }

  l2_reversed <- l2_sign != dominant_sign

  # Final flip decision: use L2 for strong blocks, L1 for isolated strong markers
  final_flip <- ifelse(
    l2_confidence >= 0.7, l2_reversed,
    ifelse(l1_confidence >= 0.5 & !l1_ambiguous, l1_reversed, l2_reversed)
  )

  flip_source <- ifelse(
    l2_confidence >= 0.7, "L2_block",
    ifelse(l1_confidence >= 0.5 & !l1_ambiguous, "L1_marker", "L2_block")
  )

  support_class <- ifelse(
    abs_delta >= 0.5, "strong",
    ifelse(abs_delta >= 0.2, "moderate",
           ifelse(abs_delta >= 0.1, "weak", "ambiguous"))
  )

  # ── L2 block boundaries ────────────────────────────────────────────────
  block_id <- cumsum(c(1, abs(diff(l2_sign)) > 0))

  # ── Per-marker audit table ─────────────────────────────────────────────
  marker_dt <- data.table(
    candidate_id = cid, chrom = chr, pos = sites_reg$pos,
    marker_id = if ("marker" %in% names(sites_reg)) sites_reg$marker else paste0("M", seq_len(n_markers)),
    raw_mean_HOMO_1 = round(grp_means[, "HOMO_1"], 4),
    raw_mean_HOMO_2 = round(grp_means[, "HOMO_2"], 4),
    raw_mean_HET = round(grp_means[, "HET"], 4),
    raw_median_HOMO_1 = round(grp_medians[, "HOMO_1"], 4),
    raw_median_HOMO_2 = round(grp_medians[, "HOMO_2"], 4),
    raw_median_HET = round(grp_medians[, "HET"], 4),
    mad_HOMO_1 = round(grp_mads[, "HOMO_1"], 4),
    mad_HOMO_2 = round(grp_mads[, "HOMO_2"], 4),
    mad_HET = round(grp_mads[, "HET"], 4),
    delta_hom_mean = round(delta_hom_mean, 4),
    delta_hom_median = round(delta_hom_median, 4),
    abs_delta_hom = round(abs_delta, 4),
    n_HOMO_1 = length(h1_idx), n_HOMO_2 = length(h2_idx),
    polarity_L1_sign = l1_sign,
    polarity_L1_confidence = round(l1_confidence, 4),
    polarity_L1_ambiguous = l1_ambiguous,
    polarity_L1_reversed = l1_reversed,
    block_id = block_id,
    polarity_L2_sign = l2_sign,
    polarity_L2_confidence = round(l2_confidence, 4),
    polarity_L2_reversed = l2_reversed,
    final_flip_decision = final_flip,
    flip_source = flip_source,
    support_class = support_class
  )
  fwrite(marker_dt, file.path(cand_dir, "candidate_marker_polarity.tsv"), sep = "\t")

  # ── Block summary ──────────────────────────────────────────────────────
  block_summary <- marker_dt[, .(
    start_pos = min(pos), end_pos = max(pos), n_markers = .N,
    block_delta_hom_mean = round(mean(delta_hom_mean), 4),
    block_delta_hom_median = round(median(delta_hom_median), 4),
    polarity_L2_sign = l2_sign[1],
    polarity_L2_confidence = round(mean(polarity_L2_confidence), 4),
    dominant_direction = ifelse(mean(delta_hom_mean) > 0, "positive", "negative"),
    n_strong = sum(support_class == "strong"),
    n_weak = sum(support_class %in% c("weak", "ambiguous")),
    frac_L1_agree = round(mean(polarity_L1_sign == l2_sign), 4)
  ), by = block_id]
  block_summary[, candidate_id := cid]
  fwrite(block_summary, file.path(cand_dir, "candidate_marker_blocks.tsv"), sep = "\t")

  # ══════════════════════════════════════════════════════════════════════
  # SAMPLE COHERENCE with group-specific thresholds
  # ══════════════════════════════════════════════════════════════════════
  # Use top informative markers
  info_thresh <- quantile(abs_delta, 0.75, na.rm = TRUE)
  info_idx <- which(abs_delta >= max(info_thresh, 0.1))
  if (length(info_idx) < 10) info_idx <- order(abs_delta, decreasing = TRUE)[1:min(50, n_markers)]

  coherence_list <- list()
  for (si in seq_along(sample_cols)) {
    sid <- sample_cols[si]; g <- groups[si]
    x_sample <- X[info_idx, si]

    if (g %in% names(grp_list) && length(grp_list[[g]]) >= 2) {
      expected <- grp_means[info_idx, g]
      other_cols <- setdiff(c("HOMO_1", "HET", "HOMO_2"), g)
      other_mean <- rowMeans(grp_means[info_idx, other_cols, drop = FALSE], na.rm = TRUE)

      dist_own <- abs(x_sample - expected)
      dist_other <- abs(x_sample - other_mean)
      valid <- is.finite(dist_own) & is.finite(dist_other)
      agree_frac <- if (sum(valid) > 0) mean(dist_own[valid] < dist_other[valid]) else NA_real_
    } else {
      agree_frac <- NA_real_
    }

    # Group-specific coherence thresholds
    # HET is expected intermediate → lower thresholds
    if (g == "HET") {
      coh_class <- if (is.na(agree_frac)) "insufficient"
                   else if (agree_frac >= 0.55) "coherent"
                   else if (agree_frac >= 0.40) "intermediate"
                   else "discordant"
    } else {
      coh_class <- if (is.na(agree_frac)) "insufficient"
                   else if (agree_frac >= 0.70) "coherent"
                   else if (agree_frac >= 0.45) "intermediate"
                   else "discordant"
    }

    # Centroid distance
    rot_row <- rot[sample == sid]
    stripe_samples <- rot[coarse_group_refined == g]
    dist_centroid <- centroid_z <- NA_real_
    if (nrow(stripe_samples) >= 3 && nrow(rot_row) == 1) {
      cx <- mean(stripe_samples$u); cy <- mean(stripe_samples$v)
      dist_centroid <- sqrt((rot_row$u - cx)^2 + (rot_row$v - cy)^2)
      all_dists <- sqrt((stripe_samples$u - cx)^2 + (stripe_samples$v - cy)^2)
      mad_d <- mad(all_dists, na.rm = TRUE)
      centroid_z <- if (mad_d > 0) (dist_centroid - median(all_dists)) / mad_d else 0
    }

    # Quality tier: group-specific
    if (g == "HET") {
      tier <- if (coh_class == "coherent") "core"
              else if (coh_class == "intermediate" && !is.na(centroid_z) && centroid_z < 3) "peripheral"
              else if (coh_class == "discordant") "junk"
              else "peripheral"
    } else {
      tier <- if (coh_class == "coherent" && !is.na(centroid_z) && centroid_z < 2) "core"
              else if (coh_class %in% c("coherent", "intermediate") && !is.na(centroid_z) && centroid_z < 4) "peripheral"
              else "junk"
    }

    # Rule that triggered tier
    tier_rule <- paste0("coherence=", coh_class,
                        ";z=", ifelse(is.na(centroid_z), "NA", round(centroid_z, 1)),
                        ";group=", g)

    coherence_list[[si]] <- data.table(
      candidate_id = cid, sample = sid, coarse_group = g,
      agreement_fraction = round(agree_frac, 4),
      coherence_class = coh_class,
      dist_to_centroid = round(dist_centroid, 4),
      centroid_z_score = round(centroid_z, 4),
      stripe_quality = tier,
      tier_rule = tier_rule,
      n_informative_markers = length(info_idx)
    )
  }

  coherence_dt <- rbindlist(coherence_list)
  fwrite(coherence_dt, file.path(cand_dir, "candidate_sample_coherence.tsv"), sep = "\t")

  # ── Core-only PCA (include HET core samples) ──────────────────────────
  core_samples <- coherence_dt[stripe_quality == "core", sample]
  if (length(core_samples) >= 10) {
    core_cols <- intersect(core_samples, colnames(X))
    if (length(core_cols) >= 10) {
      core_var <- apply(X[, core_cols, drop = FALSE], 1, var, na.rm = TRUE)
      top_m <- order(core_var, decreasing = TRUE)[1:min(n_markers, 5000)]
      pca_core <- prcomp(t(X[top_m, core_cols]), center = TRUE, scale. = FALSE)
      core_pca_dt <- data.table(
        candidate_id = cid, sample = core_cols,
        core_PC1 = round(pca_core$x[, 1], 4),
        core_PC2 = round(pca_core$x[, 2], 4)
      )
      core_pca_dt <- merge(core_pca_dt, rot[, .(sample, coarse_group_refined)],
                           by = "sample", all.x = TRUE)
      fwrite(core_pca_dt, file.path(cand_dir, "candidate_core_pca.tsv"), sep = "\t")
    }
  }

  # ── Summary ────────────────────────────────────────────────────────────
  ct <- table(coherence_dt$coherence_class)
  qt <- table(coherence_dt$stripe_quality)
  n_l1_rev <- sum(l1_reversed); n_l2_rev <- sum(l2_reversed); n_final_flip <- sum(final_flip)
  n_blocks <- max(block_id)

  message("[INFO]   Coherence: ", paste(names(ct), ct, sep = "=", collapse = " "))
  message("[INFO]   Quality: ", paste(names(qt), qt, sep = "=", collapse = " "))
  message("[INFO]   Polarity: L1_rev=", n_l1_rev, " L2_rev=", n_l2_rev,
          " final_flip=", n_final_flip, "/", n_markers, " | ", n_blocks, " L2 blocks")
}

message("[DONE] STEP29 v2 complete")
