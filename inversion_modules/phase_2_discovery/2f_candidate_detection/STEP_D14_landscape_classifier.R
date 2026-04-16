#!/usr/bin/env Rscript
# ============================================================================
# STEP_D14_landscape_classifier.R — Full-chromosome block classification + confidence
# ============================================================================
#
# Not just "find blocks" — label EVERYTHING on the chromosome.
# Every region gets a structural category, a confidence score,
# and a descriptive label for the annotated plot.
#
# CATEGORIES:
#
#   BLOCK-LEVEL (from staircase blocks):
#     strong_inversion     — high squareness, sharp edges, persists at nn40+
#     diffuse_inversion    — moderate squareness, some patchiness, still square
#     complex_system       — parent block with 2+ children = nested inversions
#     nested_fixed         — inner block at very high frequency (near-fixed)
#     nested_rare          — inner block at low frequency (rare sub-haplotype)
#
#   DIAGONAL-LEVEL (from decay shape):
#     family_ld_band       — smooth decay, high bg_level, no sharp edges
#     diffuse_diagonal     — elevated diagonal but with many small blocks inside
#     extended_suppression — broad region of mild elevation (old/weak inversion?)
#
#   ARTIFACT-LEVEL:
#     artifact_stripe      — blue cross = single window with low sim to neighbors
#     recombinant_dip      — narrow dip inside a block (double crossover?)
#     assembly_gap         — cluster of artifact stripes (assembly problem zone)
#
#   INTER-BLOCK:
#     linked_blocks        — two separate blocks with elevated cross-similarity
#     transition_zone      — gradual change between block and background
#     clean_background     — no structure, near-zero similarity
#
# CONFIDENCE SCORE (0-100%):
#   Computed from: squareness, contrast, NN persistence, vote count,
#   peel stability, step sharpness, internal homogeneity.
#   Each metric contributes 0-1, weighted, averaged, scaled to %.
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# MAIN
# ============================================================================

classify_landscape <- function(blocks, smat, vote_profile = NULL,
                                step_table = NULL, nn_evidence = NULL,
                                peel_evidence = NULL, scores = NULL,
                                window_size_bp = 50000L) {

  N <- nrow(smat)
  if (nrow(blocks) == 0) return(data.table())

  # ---- Compute nesting depth ----
  blocks <- compute_depths(blocks)

  # ---- Background statistics ----
  bg_vals <- smat[upper.tri(smat)]
  bg_vals <- bg_vals[is.finite(bg_vals)]
  bg_median <- median(bg_vals)
  bg_mad <- mad(bg_vals, constant = 1.4826)

  # ---- Per-block classification ----
  results <- list()

  for (r in seq_len(nrow(blocks))) {
    b <- blocks[r]
    bi <- b$start; be <- b$end; bw <- be - bi + 1

    # ---- Gather all available evidence ----

    # Bloc scores (squareness, contrast, occupancy, patchiness)
    sc <- NULL
    if (!is.null(scores) && nrow(scores) > 0) {
      sc <- scores[block_id == b$block_id]
      if ("variant" %in% names(sc)) sc <- sc[variant == "raw"]
      if (nrow(sc) > 0) sc <- sc[1] else sc <- NULL
    }

    squareness <- if (!is.null(sc) && "squareness" %in% names(sc)) sc$squareness else NA
    contrast   <- if (!is.null(sc) && "contrast" %in% names(sc)) sc$contrast else {
      if (!is.na(b$height)) b$height - bg_median else NA
    }
    occupancy  <- if (!is.null(sc) && "occupancy" %in% names(sc)) sc$occupancy else NA
    patchiness <- if (!is.null(sc) && "patchiness" %in% names(sc)) sc$patchiness else NA
    shape      <- if (!is.null(sc) && "shape_class" %in% names(sc)) sc$shape_class else "unknown"

    # NN persistence
    survives_nn40 <- NA; survives_nn80 <- NA
    if (!is.null(nn_evidence) && nrow(nn_evidence) > 0) {
      nn_row <- nn_evidence[candidate_id == b$block_id]
      if (nrow(nn_row) > 0) {
        survives_nn40 <- if ("survives_nn40" %in% names(nn_row)) nn_row$survives_nn40[1] else NA
        survives_nn80 <- if ("survives_nn80" %in% names(nn_row)) nn_row$survives_nn80[1] else NA
      }
    }

    # Peel diagnostic
    peel_l1b <- NA_character_; peel_l2 <- NA_character_
    if (!is.null(peel_evidence) && nrow(peel_evidence) > 0) {
      pe_l1b <- peel_evidence[block_id == b$block_id & peel_mode == "L1b_chrlocal_kin"]
      pe_l2  <- peel_evidence[block_id == b$block_id & peel_mode == "L2_local_coseg"]
      if (nrow(pe_l1b) > 0) peel_l1b <- pe_l1b$effect_class[1]
      if (nrow(pe_l2) > 0) peel_l2 <- pe_l2$effect_class[1]
    }

    # Vote count at boundaries
    vote_left <- 0; vote_right <- 0
    if (!is.null(vote_profile) && nrow(vote_profile) > 0) {
      if (bi >= 1 && bi <= nrow(vote_profile))
        vote_left <- vote_profile$vote_smooth[bi]
      if (be >= 1 && be <= nrow(vote_profile))
        vote_right <- vote_profile$vote_smooth[be]
    }

    # Step sharpness (from block table)
    step_left  <- if ("step_down_left" %in% names(b)) b$step_down_left else NA
    step_right <- if ("step_down_right" %in% names(b)) b$step_down_right else NA

    # n_children
    n_children <- sum(blocks$parent_id == b$block_id, na.rm = TRUE)

    # Internal artifacts
    n_artifacts <- if ("n_artifacts" %in% names(b)) b$n_artifacts else 0L

    # ---- Structural category ----
    category <- classify_block(
      squareness = squareness, contrast = contrast, occupancy = occupancy,
      patchiness = patchiness, depth = b$depth, n_children = n_children,
      survives_nn40 = survives_nn40, survives_nn80 = survives_nn80,
      peel_l1b = peel_l1b, peel_l2 = peel_l2,
      width = bw, n_artifacts = n_artifacts,
      height = b$height, bg_median = bg_median
    )

    # ---- Confidence score (0-100%) ----
    confidence <- compute_confidence(
      squareness = squareness, contrast = contrast, occupancy = occupancy,
      survives_nn40 = survives_nn40, survives_nn80 = survives_nn80,
      peel_l1b = peel_l1b, peel_l2 = peel_l2,
      step_left = step_left, step_right = step_right,
      vote_left = vote_left, vote_right = vote_right,
      patchiness = patchiness, width = bw
    )

    # ---- Descriptive label ----
    label <- build_rich_label(
      block_id = b$block_id, category = category, confidence = confidence,
      start_mb = b$start_mb, end_mb = b$end_mb, depth = b$depth,
      n_children = n_children, n_artifacts = n_artifacts,
      squareness = squareness, height = b$height,
      survives_nn40 = survives_nn40, peel_l1b = peel_l1b
    )

    results[[r]] <- data.table(
      block_id       = b$block_id,
      start          = bi, end = be, width = bw,
      start_mb       = b$start_mb, end_mb = b$end_mb,
      depth          = b$depth,
      parent_id      = b$parent_id,
      n_children     = n_children,
      category       = category,
      confidence_pct = confidence,
      label          = label,
      # Raw evidence
      height         = b$height,
      contrast       = round(contrast %||% NA_real_, 4),
      squareness     = round(squareness %||% NA_real_, 4),
      occupancy      = round(occupancy %||% NA_real_, 4),
      patchiness     = round(patchiness %||% NA_real_, 4),
      survives_nn40  = survives_nn40,
      survives_nn80  = survives_nn80,
      peel_l1b       = peel_l1b,
      peel_l2        = peel_l2,
      vote_left      = vote_left,
      vote_right     = vote_right,
      step_left      = round(step_left %||% NA_real_, 4),
      step_right     = round(step_right %||% NA_real_, 4),
      n_artifacts    = n_artifacts
    )
  }

  classified <- rbindlist(results, fill = TRUE)

  # ---- Label inter-block regions ----
  gaps <- classify_gaps(classified, smat, N, bg_median, bg_mad, window_size_bp)
  if (nrow(gaps) > 0) {
    classified <- rbind(classified, gaps, fill = TRUE)
    setorder(classified, start)
  }

  return(classified)
}


# ============================================================================
# BLOCK CATEGORY CLASSIFIER
# ============================================================================

classify_block <- function(squareness, contrast, occupancy, patchiness,
                            depth, n_children, survives_nn40, survives_nn80,
                            peel_l1b, peel_l2, width, n_artifacts,
                            height, bg_median) {

  sq <- squareness %||% NA; cn <- contrast %||% NA
  oc <- occupancy %||% NA; pa <- patchiness %||% NA
  nn40 <- survives_nn40 %||% NA; nn80 <- survives_nn80 %||% NA

  # Thresholds from config (with fallback defaults)
  th_strong <- if (exists("CFG")) CFG$LANDSCAPE_STRONG_SQ %||% 0.65 else 0.65
  th_diffuse <- if (exists("CFG")) CFG$LANDSCAPE_DIFFUSE_SQ %||% 0.40 else 0.40
  th_family <- if (exists("CFG")) CFG$LANDSCAPE_FAMILY_SQ %||% 0.35 else 0.35
  th_fixed_h <- if (exists("CFG")) CFG$LANDSCAPE_FIXED_HEIGHT %||% 0.85 else 0.85

  if (n_children >= 2) return("complex_system")

  if (depth >= 1) {
    if (!is.na(height) && height > th_fixed_h) return("nested_fixed")
    if (!is.na(oc) && oc < 0.3)                return("nested_rare")
    return("nested_sub_block")
  }

  if ((!is.na(sq) && sq > th_strong) &&
      (!is.na(cn) && cn > 0.05) &&
      (is.na(nn40) || isTRUE(nn40))) {
    return("strong_inversion")
  }

  if ((!is.na(sq) && sq > th_diffuse) &&
      (!is.na(cn) && cn > 0.03)) {
    return("diffuse_inversion")
  }

  if ((!is.na(sq) && sq < th_family) &&
      (!is.na(pa) && pa > 0.25)) {
    if (!is.na(peel_l1b) && peel_l1b == "disappeared") return("family_ld_band")
    return("diffuse_diagonal")
  }

  if (n_artifacts > 3 && !is.na(sq) && sq < 0.5) {
    return("diffuse_high_density")
  }

  if (!is.na(cn) && cn < 0.02) return("weak_signal")

  return("unclassified")
}


# ============================================================================
# CONFIDENCE SCORE (0-100%)
# ============================================================================

compute_confidence <- function(squareness, contrast, occupancy,
                                survives_nn40, survives_nn80,
                                peel_l1b, peel_l2,
                                step_left, step_right,
                                vote_left, vote_right,
                                patchiness, width) {

  scores <- numeric(0)
  weights <- numeric(0)

  # Scaling constants from config
  contrast_scale <- if (exists("CFG")) CFG$LANDSCAPE_CONTRAST_SCALE %||% 0.15 else 0.15
  vote_scale     <- if (exists("CFG")) CFG$LANDSCAPE_VOTE_SCALE %||% 50 else 50
  step_scale     <- if (exists("CFG")) CFG$LANDSCAPE_STEP_SCALE %||% 0.10 else 0.10

  if (!is.na(squareness)) {
    scores <- c(scores, min(1, squareness))
    weights <- c(weights, 3)
  }
  if (!is.na(contrast)) {
    scores <- c(scores, min(1, contrast / contrast_scale))
    weights <- c(weights, 2)
  }
  if (!is.na(occupancy)) {
    scores <- c(scores, occupancy)
    weights <- c(weights, 1)
  }
  if (!is.na(survives_nn40)) {
    scores <- c(scores, if (isTRUE(survives_nn40)) 1 else 0)
    weights <- c(weights, 3)
  }
  if (!is.na(survives_nn80)) {
    scores <- c(scores, if (isTRUE(survives_nn80)) 1 else 0)
    weights <- c(weights, 2)
  }
  if (!is.na(peel_l1b)) {
    peel_score <- switch(peel_l1b,
      stable = 1.0, weakened = 0.5, disappeared = 0.0,
      revealed_child = 0.7, 0.3)
    scores <- c(scores, peel_score)
    weights <- c(weights, 3)
  }
  if (!is.na(peel_l2)) {
    peel_score <- switch(peel_l2,
      stable = 1.0, weakened = 0.6, disappeared = 0.1,
      revealed_child = 0.8, 0.3)
    scores <- c(scores, peel_score)
    weights <- c(weights, 2)
  }
  mean_step <- mean(c(step_left %||% NA, step_right %||% NA), na.rm = TRUE)
  if (is.finite(mean_step)) {
    scores <- c(scores, min(1, mean_step / step_scale))
    weights <- c(weights, 1)
  }
  mean_vote <- mean(c(vote_left, vote_right), na.rm = TRUE)
  if (is.finite(mean_vote) && mean_vote > 0) {
    scores <- c(scores, min(1, mean_vote / vote_scale))
    weights <- c(weights, 2)
  }
  if (!is.na(patchiness)) {
    scores <- c(scores, max(0, 1 - patchiness * 3))
    weights <- c(weights, 1)
  }

  if (length(scores) == 0) return(NA_real_)

  # Weighted average → scale to 0-100
  conf <- sum(scores * weights) / sum(weights) * 100
  return(round(min(100, max(0, conf)), 0))
}


# ============================================================================
# RICH LABEL BUILDER
# ============================================================================

build_rich_label <- function(block_id, category, confidence,
                              start_mb, end_mb, depth, n_children,
                              n_artifacts, squareness, height,
                              survives_nn40, peel_l1b) {

  span <- round(end_mb - start_mb, 1)

  # Category short name
  cat_label <- switch(category,
    strong_inversion   = "Strong Inv",
    diffuse_inversion  = "Diffuse Inv",
    complex_system     = paste0("Complex (", n_children, " sub)"),
    nested_fixed       = "Nested Fixed",
    nested_rare        = "Nested Rare",
    nested_sub_block   = "Nested Sub",
    family_ld_band     = "Family LD",
    diffuse_diagonal   = "Diffuse Diag",
    diffuse_high_density = "Dense Scatter",
    extended_suppression = "Ext Suppress",
    weak_signal        = "Weak",
    "??"
  )

  # Confidence badge
  conf_badge <- if (!is.na(confidence)) {
    paste0(confidence, "%")
  } else "?%"

  # NN badge
  nn_badge <- if (!is.na(survives_nn40)) {
    if (isTRUE(survives_nn40)) "nn40+" else "nn40-"
  } else ""

  # Peel badge
  peel_badge <- if (!is.na(peel_l1b)) {
    switch(peel_l1b,
      stable = "peel:OK", weakened = "peel:weak",
      disappeared = "peel:GONE", revealed_child = "peel:REVEAL", "")
  } else ""

  # Depth prefix
  depth_prefix <- if (depth == 0) "Outer" else paste0("L", depth)

  # Assemble
  line1 <- paste0("B", block_id, " ", depth_prefix, " ", cat_label,
                  " [", conf_badge, "]")
  line2 <- paste0(round(start_mb, 1), "-", round(end_mb, 1),
                  " Mb (", span, " Mb)")
  line3 <- paste(c(nn_badge, peel_badge), collapse = " ")
  line3 <- trimws(line3)

  paste(c(line1, line2, if (nchar(line3) > 0) line3), collapse = "\n")
}


# ============================================================================
# CLASSIFY GAPS BETWEEN BLOCKS
# ============================================================================

classify_gaps <- function(classified, smat, N, bg_median, bg_mad,
                           window_size_bp) {
  if (nrow(classified) == 0) return(data.table())

  setorder(classified, start)
  gaps <- data.table()

  # Check regions before first block, between blocks, after last block
  all_starts <- c(1L, classified$end + 1L)
  all_ends   <- c(classified$start - 1L, N)

  for (g in seq_along(all_starts)) {
    gs <- all_starts[g]; ge <- all_ends[g]
    if (ge - gs < 5) next

    # Sample similarity in this gap region
    sub <- smat[gs:ge, gs:ge]
    ut <- sub[upper.tri(sub)]
    ut <- ut[is.finite(ut)]
    if (length(ut) < 10) next

    gap_median <- median(ut)
    gap_elevation <- gap_median - bg_median

    # Off-diagonal check: similarity to distant windows
    far_range <- max(1, gs - 500):min(N, ge + 500)
    far_range <- far_range[far_range < gs | far_range > ge]
    if (length(far_range) > 100) far_range <- sample(far_range, 100)
    far_sims <- c()
    for (fi in far_range) {
      for (gi in seq(gs, ge, length.out = min(10, ge - gs + 1))) {
        v <- smat[as.integer(gi), fi]
        if (is.finite(v)) far_sims <- c(far_sims, v)
      }
    }
    far_median <- if (length(far_sims) > 0) median(far_sims) else bg_median

    # Classify
    if (gap_elevation < bg_mad) {
      cat_gap <- "clean_background"
    } else if (far_median > bg_median + bg_mad) {
      cat_gap <- "family_ld_band"
    } else if (gap_elevation > 3 * bg_mad) {
      cat_gap <- "extended_suppression"
    } else {
      cat_gap <- "transition_zone"
    }

    gaps <- rbind(gaps, data.table(
      block_id    = NA_integer_,
      start       = gs, end = ge, width = ge - gs + 1,
      start_mb    = (gs - 1) * window_size_bp / 1e6,
      end_mb      = ge * window_size_bp / 1e6,
      depth       = NA_integer_,
      parent_id   = NA_integer_,
      n_children  = 0L,
      category    = cat_gap,
      confidence_pct = NA_real_,
      label       = paste0(cat_gap, " [",
                           round((gs-1)*window_size_bp/1e6, 1), "-",
                           round(ge*window_size_bp/1e6, 1), " Mb]"),
      height      = round(gap_median, 4),
      contrast    = round(gap_elevation, 4)
    ), fill = TRUE)
  }

  return(gaps)
}


# ============================================================================
# HELPER
# ============================================================================

compute_depths <- function(blocks) {
  if (nrow(blocks) == 0) return(blocks)
  blocks$depth <- 0L
  for (r in seq_len(nrow(blocks))) {
    pid <- blocks$parent_id[r]
    d <- 0L
    while (!is.na(pid) && d < 10) {
      d <- d + 1L
      pr <- which(blocks$block_id == pid)
      if (length(pr) == 0) break
      pid <- blocks$parent_id[pr[1]]
    }
    blocks$depth[r] <- d
  }
  return(blocks)
}


# ============================================================================
# CLI
# ============================================================================

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  precomp_file <- args[1]
  blocks_file  <- args[2]
  outdir       <- if (length(args) > 2) args[3] else "."

  pc <- readRDS(precomp_file)
  blocks <- fread(blocks_file)
  chr <- pc$chrom

  # Try loading optional evidence
  scores_file <- sub("blocks_", "block_scores_", blocks_file)
  nn_file     <- sub("blocks_", "nn_persistence_", blocks_file)
  peel_file   <- sub("blocks_", "peel_diagnostic_", blocks_file)
  vote_file   <- Sys.glob(file.path(dirname(blocks_file), paste0("*", chr, "*votes*")))

  scores <- if (file.exists(scores_file)) fread(scores_file) else NULL
  nn_ev  <- if (file.exists(nn_file)) fread(nn_file) else NULL
  peel   <- if (file.exists(peel_file)) fread(peel_file) else NULL
  votes  <- if (length(vote_file) > 0 && file.exists(vote_file[1])) fread(vote_file[1]) else NULL

  result <- classify_landscape(blocks, pc$sim_mat,
                                vote_profile = votes,
                                nn_evidence = nn_ev,
                                peel_evidence = peel,
                                scores = scores)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(result, file.path(outdir, paste0("landscape_", chr, ".tsv")), sep = "\t")

  # Print summary
  cat("\n=== LANDSCAPE CLASSIFICATION:", chr, "===\n")
  tab <- table(result$category)
  for (cat_name in sort(names(tab))) {
    cat("  ", cat_name, ":", tab[cat_name], "\n")
  }

  # Print top candidates
  strong <- result[category %in% c("strong_inversion", "complex_system") &
                    !is.na(confidence_pct)]
  if (nrow(strong) > 0) {
    setorder(strong, -confidence_pct)
    cat("\nTop candidates:\n")
    for (r in seq_len(min(10, nrow(strong)))) {
      cat("  ", strong$label[r], "\n")
    }
  }

  cat("\n[DONE] -> ", outdir, "\n")
}
