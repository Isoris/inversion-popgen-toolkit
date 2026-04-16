#!/usr/bin/env Rscript
# =============================================================================
# cheat7_beta_adaptive_threshold.R — Beta distribution adaptive thresholding
#
# Fit beta distribution to per-chromosome inv_likeness (or any 0-1 score).
# Background windows follow Beta(α, β) with mass near 0. Inversion windows
# sit in the extreme tail.
#
# Per-chromosome: a noisy chromosome has right-shifted beta → higher threshold.
# A clean chromosome has left-shifted beta → catches weaker inversions.
#
# Two levels:
#   1. Beta prior per variable: reduces chromosome-specific noise
#   2. Joint classifier: PVE1 × het_CV joint distribution
#
# REQUIRES: MASS::fitdistr (for beta MLE)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

#' Fit beta distribution and compute per-window p-values
#' @param scores Numeric vector of 0-1 scores (e.g., inv_likeness)
#' @param alpha_prior Prior alpha for regularization (default 0.5)
#' @return data.table with original scores + beta_pval + adaptive_seed
beta_adaptive_pvalues <- function(scores, alpha_prior = 0.5) {
  valid <- scores[is.finite(scores) & scores > 0.001 & scores < 0.999]
  if (length(valid) < 50) {
    return(data.table(score = scores, beta_pval = NA_real_,
                       adaptive_seed = FALSE))
  }

  # MLE fit of Beta distribution
  fit <- tryCatch({
    if (requireNamespace("MASS", quietly = TRUE)) {
      MASS::fitdistr(valid, "beta", start = list(shape1 = 1, shape2 = 5))
    } else {
      # Method of moments fallback
      m <- mean(valid); v <- var(valid)
      if (v >= m * (1 - m)) v <- m * (1 - m) * 0.99  # regularize
      alpha <- m * (m * (1 - m) / v - 1)
      beta_param  <- (1 - m) * (m * (1 - m) / v - 1)
      list(estimate = c(shape1 = max(0.1, alpha), shape2 = max(0.1, beta_param)))
    }
  }, error = function(e) {
    # Fallback: assume shape1=1, shape2=5 (right-skewed background)
    list(estimate = c(shape1 = 1, shape2 = 5))
  })

  a <- fit$estimate["shape1"]
  b <- fit$estimate["shape2"]

  # P-value: probability of seeing this score or higher under the fitted beta
  pvals <- 1 - pbeta(scores, a, b)
  pvals[!is.finite(scores)] <- NA_real_

  data.table(
    score = scores,
    beta_alpha = round(a, 3),
    beta_beta = round(b, 3),
    beta_pval = pvals,
    adaptive_seed = pvals < 0.01 & is.finite(pvals)
  )
}

#' Apply beta adaptive thresholding to a chromosome's inv_likeness
#' @param dt Precomp data.table (must have inv_likeness column)
#' @return dt with added columns: beta_pval, adaptive_seed, beta_alpha, beta_beta
apply_beta_threshold_chromosome <- function(dt) {
  if (!"inv_likeness" %in% names(dt) || nrow(dt) < 50) {
    dt[, `:=`(beta_pval = NA_real_, adaptive_seed = FALSE,
              beta_alpha = NA_real_, beta_beta = NA_real_)]
    return(dt)
  }

  bt <- beta_adaptive_pvalues(dt$inv_likeness)
  dt[, `:=`(
    beta_pval      = bt$beta_pval,
    adaptive_seed  = bt$adaptive_seed,
    beta_alpha     = bt$beta_alpha,
    beta_beta      = bt$beta_beta
  )]

  n_seeds <- sum(dt$adaptive_seed, na.rm = TRUE)
  message("[cheat7] Beta adaptive: α=", round(bt$beta_alpha[1], 2),
          " β=", round(bt$beta_beta[1], 2),
          " → ", n_seeds, " adaptive seeds (p<0.01)")

  dt
}

#' Joint classifier: PVE1 × het_contrast bivariate beta
#' Fits marginal betas to both, then flags windows where BOTH are extreme
#' @param dt Precomp data.table (must have inv_pve1, inv_het_contrast)
#' @return dt with joint_beta_pval column
apply_joint_beta_classifier <- function(dt) {
  if (!all(c("inv_pve1", "inv_het_contrast") %in% names(dt))) {
    dt[, joint_beta_pval := NA_real_]
    return(dt)
  }

  # Normalize both to 0-1 for beta fitting
  pve1_norm <- pmin(1, pmax(0.001, dt$inv_pve1))
  het_norm  <- pmin(1, pmax(0.001, dt$inv_het_contrast / max(dt$inv_het_contrast, na.rm = TRUE)))

  bt_pve1 <- beta_adaptive_pvalues(pve1_norm)
  bt_het  <- beta_adaptive_pvalues(het_norm)

  # Joint p-value: Fisher's method (combine two independent p-values)
  p1 <- bt_pve1$beta_pval
  p2 <- bt_het$beta_pval
  # -2 * (log(p1) + log(p2)) ~ chi-sq(4) under null
  joint_stat <- -2 * (log(pmax(p1, 1e-300)) + log(pmax(p2, 1e-300)))
  joint_pval <- pchisq(joint_stat, df = 4, lower.tail = FALSE)
  joint_pval[!is.finite(joint_stat)] <- NA_real_

  dt[, joint_beta_pval := joint_pval]

  n_joint <- sum(joint_pval < 0.01, na.rm = TRUE)
  message("[cheat7] Joint beta classifier: ", n_joint, " windows with p<0.01")

  dt
}
