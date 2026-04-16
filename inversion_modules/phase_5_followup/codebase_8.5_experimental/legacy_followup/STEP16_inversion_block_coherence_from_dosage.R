#!/usr/bin/env Rscript

# =============================================================================
# STEP16_inversion_block_coherence_from_dosage.R
#
# Per-sample haploblock coherence metrics from dosage matrices.
# For each candidate region, uses the k-means groups from STEP12 to:
#   1. Identify informative markers (strong dosage difference between groups)
#   2. Compute per-sample support for each arrangement
#   3. Measure block coherence (switch rate, long-range agreement)
#   4. Assign provisional simple labels
#
# Conceptual basis:
#   Inversions suppress recombination → markers inherited as a block →
#   coherent samples show consistent support for one arrangement across
#   many markers with low switching.
#
# Usage:
#   Rscript STEP16_inversion_block_coherence_from_dosage.R \
#     <candidate_table> <step12_dir> <dosage_dir> <outdir> \
#     [min_dosage_diff=0.5] [min_informative=20]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript STEP16_... <candidate_table> <step12_dir> <dosage_dir> <outdir> [min_dosage_diff=0.5] [min_informative=20]")
}

candidate_file   <- args[1]
step12_dir       <- args[2]
dosage_dir       <- args[3]
outdir           <- args[4]
min_dosage_diff  <- if (length(args) >= 5) as.numeric(args[5]) else 0.5
min_informative  <- if (length(args) >= 6) as.integer(args[6]) else 20L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Read candidates ────────────────────────────────────────────────────────
cand <- fread(candidate_file)
stopifnot(all(c("candidate_id", "chrom", "start_bp", "end_bp") %in% names(cand)))
message("[INFO] Candidates: ", nrow(cand))

all_results <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- as.character(row$candidate_id)
  chr <- row$chrom
  inv_start <- as.numeric(row$start_bp)
  inv_end   <- as.numeric(row$end_bp)

  # ── Find STEP12 PCA file ──────────────────────────────────────────────
  pca_pattern <- paste0("STEP12_", chr, ".candidate_", cid, ".regional_pca_samples.tsv.gz")
  pca_file <- file.path(step12_dir, pca_pattern)
  if (!file.exists(pca_file)) {
    message("[WARN] Missing STEP12 file for candidate ", cid, ": ", pca_file)
    next
  }

  pcs <- fread(pca_file)
  if (!all(c("sample", "ordered_group") %in% names(pcs))) {
    message("[WARN] STEP12 file missing required columns for candidate ", cid)
    next
  }

  # ── Read dosage matrix for this chromosome ────────────────────────────
  dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(dosage_dir, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[WARN] Missing dosage files for ", chr)
    next
  }

  dos <- fread(dos_file)
  sites <- fread(sites_file)

  if (!identical(sites$marker, dos$marker)) {
    message("[WARN] Sites/dosage marker mismatch for ", chr)
    next
  }

  # ── Handle Ind-style columns ──────────────────────────────────────────
  sample_cols <- setdiff(names(dos), "marker")
  if (all(grepl("^Ind[0-9]+$", sample_cols)) && "sample" %in% names(pcs)) {
    # Map from pcs sample names if possible
    # This requires a sample order file — skip if we can't map
    if (length(sample_cols) == nrow(pcs)) {
      # Assume same order as bamlist
      setnames(dos, old = sample_cols, new = pcs$sample)
      sample_cols <- pcs$sample
    } else {
      message("[WARN] Cannot map Ind columns for candidate ", cid, " (length mismatch)")
      next
    }
  }

  # ── Filter to candidate region SNPs ──────────────────────────────────
  keep <- which(sites$chrom == chr & sites$pos >= inv_start & sites$pos <= inv_end)
  if (length(keep) < min_informative) {
    message("[WARN] Too few SNPs (", length(keep), ") in candidate ", cid)
    next
  }

  dos_reg <- dos[keep]
  sites_reg <- sites[keep]

  # Build matrix: SNPs × samples (only samples in pcs)
  common_samples <- intersect(sample_cols, pcs$sample)
  if (length(common_samples) < 10) {
    message("[WARN] Too few common samples for candidate ", cid)
    next
  }

  X <- as.matrix(dos_reg[, ..common_samples])
  storage.mode(X) <- "double"

  # ── Identify groups ──────────────────────────────────────────────────
  pcs_sub <- pcs[sample %in% common_samples]
  grp_a <- pcs_sub[ordered_group == 1, sample]
  grp_b <- pcs_sub[ordered_group == max(ordered_group), sample]

  if (length(grp_a) < 3 || length(grp_b) < 3) {
    message("[WARN] Extreme groups too small for candidate ", cid)
    next
  }

  # ── Compute mean dosage per group per SNP ────────────────────────────
  mean_a <- rowMeans(X[, common_samples %in% grp_a, drop = FALSE], na.rm = TRUE)
  mean_b <- rowMeans(X[, common_samples %in% grp_b, drop = FALSE], na.rm = TRUE)
  dosage_diff <- abs(mean_a - mean_b)

  # Informative markers: strong dosage difference between extreme groups
  informative <- which(dosage_diff >= min_dosage_diff)
  n_info <- length(informative)

  if (n_info < min_informative) {
    message("[INFO] Candidate ", cid, ": only ", n_info, " informative markers (threshold=", min_informative, ")")
  }

  # ── Per-sample metrics ──────────────────────────────────────────────
  # For each sample: classify each informative marker as A-like or B-like
  # based on which group mean it's closer to
  results_list <- list()

  for (sid in common_samples) {
    x_sample <- X[informative, sid]

    # Distance to group A vs group B mean at each informative marker
    dist_a <- abs(x_sample - mean_a[informative])
    dist_b <- abs(x_sample - mean_b[informative])

    a_support <- sum(dist_a < dist_b, na.rm = TRUE)
    b_support <- sum(dist_b < dist_a, na.rm = TRUE)
    total_classified <- a_support + b_support

    mean_a_support <- if (total_classified > 0) a_support / total_classified else NA_real_
    mean_b_support <- if (total_classified > 0) b_support / total_classified else NA_real_
    support_margin <- abs(mean_a_support - mean_b_support)

    # Switch rate: how often does the marker assignment flip between consecutive markers
    closer_to <- ifelse(dist_a < dist_b, "A", ifelse(dist_b < dist_a, "B", "tie"))
    closer_to <- closer_to[closer_to != "tie"]
    switches <- if (length(closer_to) > 1) sum(closer_to[-1] != closer_to[-length(closer_to)]) else 0L
    switch_rate <- if (length(closer_to) > 1) switches / (length(closer_to) - 1) else NA_real_
    block_coherence <- if (!is.na(switch_rate)) 1 - switch_rate else NA_real_

    # Long-range agreement: split markers into 4 quartiles by position,
    # check if majority assignment is the same across quartiles
    if (n_info >= 8) {
      pos_info <- sites_reg$pos[informative]
      quartile <- as.integer(cut(pos_info, breaks = 4, labels = FALSE))
      q_majority <- tapply(closer_to[seq_along(closer_to)],
                           quartile[!closer_to %in% "tie"][seq_along(closer_to)],
                           function(v) {
                             tt <- table(v)
                             names(tt)[which.max(tt)]
                           })
      long_range_agreement <- if (length(q_majority) >= 3) {
        sum(q_majority == q_majority[1]) / length(q_majority)
      } else NA_real_
    } else {
      long_range_agreement <- NA_real_
    }

    # Provisional label
    prov_label <- if (is.na(support_margin) || n_info < min_informative) {
      "low_info"
    } else if (support_margin >= 0.6 && block_coherence >= 0.7) {
      if (mean_a_support > mean_b_support) "coherent_A" else "coherent_B"
    } else if (support_margin < 0.3) {
      "mixed_or_het"
    } else {
      "ambiguous"
    }

    results_list[[length(results_list) + 1]] <- data.table(
      candidate_id           = as.integer(cid),
      chrom                  = chr,
      sample                 = sid,
      n_informative_markers  = n_info,
      mean_A_support         = round(mean_a_support, 4),
      mean_B_support         = round(mean_b_support, 4),
      support_margin         = round(support_margin, 4),
      switch_rate            = round(switch_rate, 4),
      block_coherence        = round(block_coherence, 4),
      long_range_agreement   = round(long_range_agreement, 4),
      provisional_label      = prov_label
    )
  }

  if (length(results_list) > 0) {
    all_results[[length(all_results) + 1]] <- rbindlist(results_list)
  }

  message("[INFO] Candidate ", cid, " (", chr, "): ", n_info, " informative markers, ",
          length(common_samples), " samples processed")
}

if (length(all_results) == 0) {
  message("[WARN] No candidates processed successfully")
  quit("no", status = 0)
}

result_dt <- rbindlist(all_results, fill = TRUE)

# ── Write outputs ────────────────────────────────────────────────────────
fwrite(result_dt,
       file.path(outdir, "inversion_block_coherence_per_sample.tsv.gz"),
       sep = "\t")

# Summary per candidate
summary_dt <- result_dt[, .(
  n_samples              = .N,
  n_informative_markers  = n_informative_markers[1],
  mean_block_coherence   = round(mean(block_coherence, na.rm = TRUE), 4),
  median_block_coherence = round(median(block_coherence, na.rm = TRUE), 4),
  mean_support_margin    = round(mean(support_margin, na.rm = TRUE), 4),
  n_coherent_A           = sum(provisional_label == "coherent_A"),
  n_coherent_B           = sum(provisional_label == "coherent_B"),
  n_mixed_or_het         = sum(provisional_label == "mixed_or_het"),
  n_ambiguous            = sum(provisional_label == "ambiguous"),
  n_low_info             = sum(provisional_label == "low_info")
), by = .(candidate_id, chrom)]

fwrite(summary_dt,
       file.path(outdir, "inversion_block_coherence_summary.tsv"),
       sep = "\t")

message("[DONE] Wrote:")
message("  ", file.path(outdir, "inversion_block_coherence_per_sample.tsv.gz"))
message("  ", file.path(outdir, "inversion_block_coherence_summary.tsv"))
