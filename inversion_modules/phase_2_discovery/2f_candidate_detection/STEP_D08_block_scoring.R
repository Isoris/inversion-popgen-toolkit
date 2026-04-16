#!/usr/bin/env Rscript
# ============================================================================
# STEP_D08_block_scoring.R — Per-block quality metrics
# ============================================================================
#
# For each detected block (from staircase), compute:
#   - inside_mean, outside_mean, contrast
#   - squareness (far/near ratio — does sim decay with distance?)
#   - sharpness (edge vs interior contrast)
#   - occupancy (fraction of cells above background)
#   - patchiness (internal heterogeneity / CV)
#   - shape_class: strong_square, diffuse_square, diagonal_band, noise
#
# Usage:
#   source("00_config.R")
#   source("STEP_D08_block_scoring.R")
#   scored <- score_blocks(blocks, smat)
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# MAIN: Score all blocks on a single matrix variant
# ============================================================================

score_blocks <- function(blocks, smat, variant_name = "raw",
                          near_frac = CFG$BLOCK_NEAR_FRAC,
                          far_frac  = CFG$BLOCK_FAR_FRAC,
                          edge_depth = CFG$BLOCK_EDGE_DEPTH) {

  if (nrow(blocks) == 0) return(blocks)

  N <- nrow(smat)
  cat("Block scoring |", nrow(blocks), "blocks |", variant_name, "\n")

  # Pre-compute whole-matrix background
  # Sample off-diagonal values far from any block
  bg_vals <- sample_background(smat, blocks, n_sample = 5000)
  bg_median <- median(bg_vals, na.rm = TRUE)
  bg_mad    <- mad(bg_vals, constant = 1.4826, na.rm = TRUE)
  cat("  Background: median =", round(bg_median, 4),
      "MAD =", round(bg_mad, 4), "\n")

  # Score each block
  scores <- list()

  for (r in seq_len(nrow(blocks))) {
    bi <- blocks$start[r]
    be <- blocks$end[r]
    bw <- be - bi + 1
    bid <- blocks$block_id[r]

    if (bw < 4 || bi < 1 || be > N) {
      scores[[r]] <- make_empty_score(bid, variant_name)
      next
    }

    sub <- smat[bi:be, bi:be]

    # ---- inside_mean / outside_mean / contrast ----
    ut_vals <- sub[upper.tri(sub)]
    ut_vals <- ut_vals[is.finite(ut_vals)]
    inside_mean <- if (length(ut_vals) > 0) mean(ut_vals) else NA_real_
    inside_med  <- if (length(ut_vals) > 0) median(ut_vals) else NA_real_

    # Outside: flanking regions of same width
    flank_w <- min(bw, 50L)
    left_start  <- max(1L, bi - flank_w)
    right_end   <- min(N, be + flank_w)

    # Cross-block similarity (block × left flank)
    outside_vals <- c()
    if (bi - left_start > 0) {
      cross_left <- smat[bi:be, left_start:(bi - 1)]
      outside_vals <- c(outside_vals, cross_left[is.finite(cross_left)])
    }
    if (right_end - be > 0) {
      cross_right <- smat[bi:be, (be + 1):right_end]
      outside_vals <- c(outside_vals, cross_right[is.finite(cross_right)])
    }
    outside_mean <- if (length(outside_vals) > 0) mean(outside_vals) else bg_median

    contrast <- inside_mean - outside_mean

    # ---- squareness (far/near ratio) ----
    # Near = close to diagonal (small |i-j|), Far = far from diagonal (large |i-j|)
    near_limit <- max(1L, round(bw * near_frac))
    far_start  <- max(1L, round(bw * (1 - far_frac)))

    near_vals <- c()
    far_vals  <- c()
    for (ii in 1:bw) {
      for (jj in (ii + 1):bw) {
        d <- jj - ii
        v <- sub[ii, jj]
        if (!is.finite(v)) next
        if (d <= near_limit) near_vals <- c(near_vals, v)
        if (d >= far_start)  far_vals  <- c(far_vals, v)
      }
    }

    near_med <- if (length(near_vals) > 0) median(near_vals) else NA_real_
    far_med  <- if (length(far_vals) > 0)  median(far_vals)  else NA_real_
    squareness <- if (!is.na(near_med) && near_med > 0) far_med / near_med else NA_real_

    # ---- sharpness (edge contrast) ----
    # Compare bins at the block boundary (inside edge vs just outside)
    edge_inside_left  <- smat[bi:min(bi + edge_depth - 1, be),
                              bi:min(bi + edge_depth - 1, be)]
    edge_outside_left <- if (bi > edge_depth) {
      smat[(bi - edge_depth):(bi - 1), bi:min(bi + edge_depth - 1, be)]
    } else NULL

    edge_inside_right  <- smat[max(be - edge_depth + 1, bi):be,
                               max(be - edge_depth + 1, bi):be]
    edge_outside_right <- if (be + edge_depth <= N) {
      smat[(be + 1):(be + edge_depth), max(be - edge_depth + 1, bi):be]
    } else NULL

    edge_in_vals <- c(edge_inside_left[is.finite(edge_inside_left)],
                      edge_inside_right[is.finite(edge_inside_right)])
    edge_out_vals <- c()
    if (!is.null(edge_outside_left))
      edge_out_vals <- c(edge_out_vals, edge_outside_left[is.finite(edge_outside_left)])
    if (!is.null(edge_outside_right))
      edge_out_vals <- c(edge_out_vals, edge_outside_right[is.finite(edge_outside_right)])

    edge_in_med  <- if (length(edge_in_vals) > 0)  median(edge_in_vals)  else NA_real_
    edge_out_med <- if (length(edge_out_vals) > 0) median(edge_out_vals) else NA_real_
    sharpness <- if (!is.na(edge_in_med) && !is.na(edge_out_med)) {
      edge_in_med - edge_out_med
    } else NA_real_

    # ---- occupancy (fraction above background) ----
    occ_threshold <- bg_median + bg_mad
    occupancy <- if (length(ut_vals) > 0) {
      sum(ut_vals > occ_threshold) / length(ut_vals)
    } else NA_real_

    # ---- patchiness (CV of interior) ----
    interior_cv <- if (length(ut_vals) > 2 && inside_mean > 0) {
      sd(ut_vals, na.rm = TRUE) / inside_mean
    } else NA_real_
    patchiness <- interior_cv  # higher = more patchy

    # ---- shape classification ----
    shape_class <- classify_shape(squareness, contrast, occupancy,
                                   patchiness, sharpness, bw)

    scores[[r]] <- data.table(
      block_id      = bid,
      variant       = variant_name,
      inside_mean   = round(inside_mean, 4),
      inside_med    = round(inside_med, 4),
      outside_mean  = round(outside_mean, 4),
      contrast      = round(contrast, 4),
      squareness    = round(squareness, 4),
      near_med      = round(near_med, 4),
      far_med       = round(far_med, 4),
      sharpness     = round(sharpness, 4),
      occupancy     = round(occupancy, 4),
      patchiness    = round(patchiness, 4),
      shape_class   = shape_class,
      bg_median     = round(bg_median, 4),
      bg_mad        = round(bg_mad, 4)
    )
  }

  return(rbindlist(scores))
}


# ============================================================================
# SHAPE CLASSIFIER
# ============================================================================

classify_shape <- function(squareness, contrast, occupancy,
                            patchiness, sharpness, width) {
  # Conservative rules — prefer "unknown" over false classification
  if (is.na(squareness) || is.na(contrast)) return("unknown")

  # Strong square: high squareness + high contrast + low patchiness + sharp edges
  if (squareness > 0.7 && contrast > 0.05 && !is.na(occupancy) &&
      occupancy > 0.5 && !is.na(patchiness) && patchiness < 0.3) {
    return("strong_square")
  }

  # Diffuse square: moderate squareness but some patchiness
  if (squareness > 0.5 && contrast > 0.03 && !is.na(occupancy) &&
      occupancy > 0.3) {
    return("diffuse_square")
  }

  # Diagonal band: low squareness (sim decays with distance)
  if (squareness < 0.4 && contrast > 0.02) {
    return("diagonal_band")
  }

  # Noise: low contrast, low occupancy
  if (contrast < 0.02 || (!is.na(occupancy) && occupancy < 0.2)) {
    return("noise")
  }

  return("ambiguous")
}


# ============================================================================
# BACKGROUND SAMPLING
# ============================================================================

sample_background <- function(smat, blocks, n_sample = 5000) {
  N <- nrow(smat)

  # Mark block regions
  in_block <- rep(FALSE, N)
  for (r in seq_len(nrow(blocks))) {
    bi <- max(1L, blocks$start[r])
    be <- min(N, blocks$end[r])
    in_block[bi:be] <- TRUE
  }

  # Sample from non-block regions
  non_block <- which(!in_block)
  if (length(non_block) < 10) {
    # Fallback: use corners / edges of matrix
    non_block <- c(1:min(20, N), max(1, N - 20):N)
  }

  vals <- numeric(0)
  attempts <- 0L
  while (length(vals) < n_sample && attempts < n_sample * 5) {
    i <- sample(non_block, 1)
    j <- sample(non_block, 1)
    if (abs(i - j) > 10) {
      v <- smat[min(i, j), max(i, j)]
      if (is.finite(v)) vals <- c(vals, v)
    }
    attempts <- attempts + 1L
  }

  return(vals)
}


# ============================================================================
# EMPTY SCORE ROW
# ============================================================================

make_empty_score <- function(bid, variant_name) {
  data.table(
    block_id = bid, variant = variant_name,
    inside_mean = NA_real_, inside_med = NA_real_,
    outside_mean = NA_real_, contrast = NA_real_,
    squareness = NA_real_, near_med = NA_real_, far_med = NA_real_,
    sharpness = NA_real_, occupancy = NA_real_,
    patchiness = NA_real_, shape_class = "unknown",
    bg_median = NA_real_, bg_mad = NA_real_
  )
}


# ============================================================================
# SCORE BLOCKS ACROSS ALL MATRIX VARIANTS
# ============================================================================

score_blocks_all_variants <- function(blocks, variants) {
  all_scores <- list()
  for (vname in names(variants)) {
    cat("\n--- Scoring on variant:", vname, "---\n")
    sc <- score_blocks(blocks, variants[[vname]], variant_name = vname)
    all_scores[[vname]] <- sc
  }
  return(rbindlist(all_scores))
}
