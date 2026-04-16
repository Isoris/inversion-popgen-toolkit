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
# Returns data.frame with columns:
#   candidate_id, interior_cv, homogeneity, stripe_count, patchiness
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
        patchiness   = NA_real_,
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
        patchiness   = NA_real_,
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

    # Patchiness: fraction of within-block values below the block median
    block_median <- median(ut_vals, na.rm=TRUE)
    below_median <- sum(ut_vals < block_median * 0.5) / length(ut_vals)

    results[[i]] <- data.frame(
      candidate_id = cid,
      interior_cv  = round(cv, 4),
      homogeneity  = round(homogeneity, 4),
      stripe_count = stripes,
      patchiness   = round(below_median, 4),
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
