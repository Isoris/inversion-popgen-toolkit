# ============================================================================
# STEP_D04_interior_cv.R — Block interior homogeneity (CV, patchiness, stripes)
# ============================================================================
# For each candidate, compute the coefficient of variation of similarity
# values within the block.
#
# Low CV = uniform signal = real inversion
# High CV = patchy signal = family block or noise
#
# Also computes stripe count: number of dark bands (local sim minima)
# crossing the interior of the block.
#
# BUGFIX 2026-04-17 (FIX 17, SILENT): the "patchiness" column here was a
# tail-fraction (fraction of cells < median/2) with semantics completely
# different from D08's CV-based `patchiness`. Both landed on tab during
# Phase 8 merges, producing `patchiness.x` (D04 tail fraction) and
# `patchiness.y` (D08 CV), so C01d's `iv$patchiness` at line 222 was NULL
# and D5_interior_quality silently skewed (`1 - 0*3 = 1` saturated
# patch_score). Renaming to `below_median_frac` preserves both quantities
# without collision. C01d reads D08's `patchiness` which is what the
# scoring logic at L222 `1 - patchiness * 3` actually assumes (CV range).
#
# Returns data.frame with columns:
#   candidate_id, interior_cv, homogeneity, stripe_count, below_median_frac
# ============================================================================

compute_interior_cv <- function(candidates, smat) {
  results <- list()

  for (i in seq_len(nrow(candidates))) {
    si <- candidates$start[i]
    ei <- candidates$end[i]
    width <- ei - si + 1
    cid <- candidates$candidate_id[i]

    if (width < 4) {
      results[[i]] <- data.frame(
        candidate_id = cid,
        interior_cv  = NA_real_,
        homogeneity  = NA_real_,
        stripe_count = NA_integer_,
        below_median_frac = NA_real_,  # BUGFIX 2026-04-17 (FIX 17)
        stringsAsFactors = FALSE
      )
      next
    }

    # Extract sub-matrix
    sub <- smat[si:ei, si:ei]
    ut_vals <- sub[upper.tri(sub)]
    ut_vals <- ut_vals[is.finite(ut_vals)]

    if (length(ut_vals) < 10) {
      results[[i]] <- data.frame(
        candidate_id = cid,
        interior_cv  = NA_real_,
        homogeneity  = NA_real_,
        stripe_count = NA_integer_,
        below_median_frac = NA_real_,  # BUGFIX 2026-04-17 (FIX 17)
        stringsAsFactors = FALSE
      )
      next
    }

    # CV
    mu <- mean(ut_vals, na.rm=TRUE)
    sd_val <- sd(ut_vals, na.rm=TRUE)
    cv <- if (mu > 0) sd_val / mu else NA
    homogeneity <- if (!is.na(cv)) max(0, 1 - cv) else NA

    # Stripe count: look for dark vertical/horizontal lines
    # Compute column means within the sub-matrix
    col_means <- colMeans(sub, na.rm=TRUE)

    # Smooth slightly to avoid single-bin noise
    if (length(col_means) > 5) {
      kernel <- rep(1/3, 3)
      col_smoothed <- stats::filter(col_means, kernel, sides=2)
      col_smoothed[is.na(col_smoothed)] <- col_means[is.na(col_smoothed)]
    } else {
      col_smoothed <- col_means
    }

    # Count local minima that are significantly below the block mean
    stripe_threshold <- mu * 0.7  # minima below 70% of mean
    stripes <- count_local_minima(col_smoothed, stripe_threshold)

    # Below-median fraction: fraction of within-block values below half the
    # block median. Previously called `patchiness` but that name clashed
    # with D08's CV-based patchiness on merge. See FIX 17 header comment.
    block_median <- median(ut_vals, na.rm=TRUE)
    below_median <- sum(ut_vals < block_median * 0.5) / length(ut_vals)

    results[[i]] <- data.frame(
      candidate_id = cid,
      interior_cv  = round(cv, 4),
      homogeneity  = round(homogeneity, 4),
      stripe_count = stripes,
      below_median_frac = round(below_median, 4),  # BUGFIX 2026-04-17 (FIX 17)
      stringsAsFactors = FALSE
    )
  }

  return(do.call(rbind, results))
}


# ---- Helper: Count local minima below threshold ----
count_local_minima <- function(x, threshold) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3) return(0L)

  count <- 0L
  for (i in 2:(n-1)) {
    if (!is.na(x[i]) && !is.na(x[i-1]) && !is.na(x[i+1])) {
      if (x[i] < x[i-1] && x[i] < x[i+1] && x[i] < threshold) {
        count <- count + 1L
      }
    }
  }
  return(count)
}
