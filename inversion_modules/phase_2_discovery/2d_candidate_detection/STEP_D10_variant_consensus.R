#!/usr/bin/env Rscript
# ============================================================================
# STEP_D10_variant_consensus.R — Cross-variant block consensus
# ============================================================================
#
# Matches candidates detected on EACH matrix variant (raw, distcorr,
# localnorm, denoised, resid_bg) by coordinate overlap.
#
# Classifies relationships: same, merged, split, nested, lost, novel
# Final consensus: blocks appearing in ≥ min_variants = high confidence
#
# Usage:
#   source("00_config.R")
#   source("STEP_D10_variant_consensus.R")
#   consensus <- build_consensus(blocks_by_variant, scores_all)
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================================
# MAIN: Build consensus from blocks detected on each variant
# ============================================================================

build_consensus <- function(blocks_by_variant,
                             scores_all = NULL,
                             overlap_thresh = CFG$CONSENSUS_OVERLAP,
                             min_variants = CFG$CONSENSUS_MIN_VARIANTS) {

  variant_names <- names(blocks_by_variant)
  n_variants <- length(variant_names)
  cat("Consensus builder | variants:", paste(variant_names, collapse = ", "), "\n")
  cat("  overlap threshold:", overlap_thresh,
      "| min variants:", min_variants, "\n")

  # ---- Step 1: Build a master block list from all variants ----
  # Each unique genomic interval gets one entry
  master <- data.table()
  master_id <- 0L

  for (vname in variant_names) {
    blocks <- blocks_by_variant[[vname]]
    if (is.null(blocks) || nrow(blocks) == 0) next

    for (r in seq_len(nrow(blocks))) {
      bs <- blocks$start[r]
      be <- blocks$end[r]

      # Check if this overlaps an existing master entry
      matched <- FALSE
      if (nrow(master) > 0) {
        for (m in seq_len(nrow(master))) {
          ov <- compute_recip_overlap(bs, be, master$start[m], master$end[m])
          if (ov >= overlap_thresh) {
            # Add this variant to the existing master entry
            old_variants <- master$variants[[m]]
            master$variants[[m]] <- c(old_variants, vname)
            master$n_variants[m] <- length(master$variants[[m]])

            # Update coordinates to union
            master$start[m] <- min(master$start[m], bs)
            master$end[m]   <- max(master$end[m], be)
            master$width[m] <- master$end[m] - master$start[m] + 1L
            matched <- TRUE
            break
          }
        }
      }

      if (!matched) {
        master_id <- master_id + 1L
        master <- rbind(master, data.table(
          consensus_id = master_id,
          start        = bs,
          end          = be,
          width        = be - bs + 1L,
          start_mb     = bin_to_mb(bs),
          end_mb       = bin_to_mb(be),
          variants     = list(vname),
          n_variants   = 1L
        ), fill = TRUE)
      }
    }
  }

  if (nrow(master) == 0) {
    cat("  No blocks found in any variant.\n")
    return(data.table())
  }

  # ---- Step 2: Flatten variants list to string ----
  master[, variant_list := sapply(variants, paste, collapse = ",")]
  master[, variants := NULL]

  # ---- Step 3: Classify each consensus entry ----
  master[, confidence := ifelse(n_variants >= min_variants, "HIGH",
                           ifelse(n_variants >= 2, "MODERATE", "LOW"))]

  # ---- Step 4: Classify variant relationships ----
  master[, relationship := classify_variant_relationship(
    variant_list, n_variants, n_variants_total = n_variants
  ), by = consensus_id]

  # ---- Step 5: Merge with scores if available ----
  if (!is.null(scores_all) && nrow(scores_all) > 0) {
    # For each consensus entry, pull the best score across variants
    master <- merge_best_scores(master, scores_all, blocks_by_variant)
  }

  # Sort by confidence then position
  setorder(master, -n_variants, start)

  cat("\n=== Consensus Results ===\n")
  cat("  Total entries:", nrow(master), "\n")
  cat("  HIGH confidence (≥", min_variants, "variants):",
      sum(master$confidence == "HIGH"), "\n")
  cat("  MODERATE (2 variants):", sum(master$confidence == "MODERATE"), "\n")
  cat("  LOW (1 variant only):", sum(master$confidence == "LOW"), "\n")

  return(master)
}


# ============================================================================
# VARIANT RELATIONSHIP CLASSIFIER
# ============================================================================

classify_variant_relationship <- function(variant_list, n_variants,
                                           n_variants_total) {
  if (n_variants == n_variants_total) return("universal")
  if (n_variants >= n_variants_total - 1) return("near_universal")

  variants <- strsplit(variant_list, ",")[[1]]

  # Check which variants found it
  has_raw      <- "raw" %in% variants
  has_distcorr <- "distcorr" %in% variants
  has_locnorm  <- "localnorm" %in% variants
  has_denoise  <- "denoised" %in% variants
  has_resid    <- "resid_bg" %in% variants

  if (!has_raw && (has_distcorr || has_resid)) {
    return("hidden_under_haze")  # only visible after distance/background correction
  }
  if (!has_raw && has_denoise) {
    return("too_patchy_for_raw")  # only visible after denoising
  }
  if (has_raw && n_variants == 1) {
    return("raw_only")  # might be noise
  }

  return("partial")
}


# ============================================================================
# MERGE BEST SCORES
# ============================================================================

merge_best_scores <- function(master, scores_all, blocks_by_variant) {
  # For each consensus entry, find the matching score with best contrast
  master$best_contrast  <- NA_real_
  master$best_squareness <- NA_real_
  master$best_occupancy  <- NA_real_
  master$best_sharpness  <- NA_real_
  master$best_shape      <- NA_character_
  master$mean_patchiness <- NA_real_

  for (m in seq_len(nrow(master))) {
    ms <- master$start[m]
    me <- master$end[m]

    # Find matching block_ids across all variants
    matching_scores <- scores_all[FALSE, ]

    for (vname in names(blocks_by_variant)) {
      blocks <- blocks_by_variant[[vname]]
      if (is.null(blocks) || nrow(blocks) == 0) next

      for (r in seq_len(nrow(blocks))) {
        ov <- compute_recip_overlap(ms, me, blocks$start[r], blocks$end[r])
        if (ov >= 0.5) {
          # Find this block's score
          sc <- scores_all[block_id == blocks$block_id[r] & variant == vname]
          if (nrow(sc) > 0) matching_scores <- rbind(matching_scores, sc)
        }
      }
    }

    if (nrow(matching_scores) > 0) {
      master$best_contrast[m]   <- max(matching_scores$contrast, na.rm = TRUE)
      master$best_squareness[m] <- max(matching_scores$squareness, na.rm = TRUE)
      master$best_occupancy[m]  <- max(matching_scores$occupancy, na.rm = TRUE)
      master$best_sharpness[m]  <- max(matching_scores$sharpness, na.rm = TRUE)
      master$mean_patchiness[m] <- mean(matching_scores$patchiness, na.rm = TRUE)

      # Best shape = from highest-contrast variant
      best_row <- which.max(matching_scores$contrast)
      master$best_shape[m] <- matching_scores$shape_class[best_row]
    }
  }

  return(master)
}


# ============================================================================
# HELPER: Reciprocal overlap (same as in 09)
# ============================================================================

compute_recip_overlap <- function(s1, e1, s2, e2) {
  w1 <- e1 - s1 + 1
  w2 <- e2 - s2 + 1
  if (w1 <= 0 || w2 <= 0) return(0)
  ov_start <- max(s1, s2)
  ov_end   <- min(e1, e2)
  ov <- max(0, ov_end - ov_start + 1)
  return(min(ov / w1, ov / w2))
}


# ============================================================================
# PRINT CONSENSUS TABLE
# ============================================================================

print_consensus <- function(consensus, top_n = 30) {
  cat("\n=== Top", min(top_n, nrow(consensus)), "Consensus Candidates ===\n")
  show <- head(consensus, top_n)
  for (r in seq_len(nrow(show))) {
    cat(sprintf("  %s #%d: %d-%d (%.1f-%.1f Mb, %d bins) | %d/%s variants | %s | %s\n",
                show$confidence[r], show$consensus_id[r],
                show$start[r], show$end[r],
                show$start_mb[r], show$end_mb[r], show$width[r],
                show$n_variants[r], show$variant_list[r],
                show$relationship[r],
                ifelse(is.na(show$best_shape[r]), "unscored", show$best_shape[r])))
  }
}
