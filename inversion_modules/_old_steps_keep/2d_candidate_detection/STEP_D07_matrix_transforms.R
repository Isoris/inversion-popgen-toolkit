#!/usr/bin/env Rscript
# ============================================================================
# STEP_D07_matrix_transforms.R — Downstream transforms on precomputed sim_mat
# ============================================================================
#
# Each transform produces a NEW matrix variant.  Original is never modified.
# All transforms are DETERMINISTIC.  Parameters exposed in 00_config.R.
#
# Transforms:
#   A1: Distance / background correction  → sim_distcorr
#   A2: Local contrast normalization       → sim_localnorm
#   A3: 2D denoising (fused lasso proxy)   → sim_denoised
#   A4: Broad background subtraction       → sim_resid_bg
#   A5: Edge / boundary maps               → sim_edge
#   A6: Thresholded support maps           → sim_support
#
# Usage:
#   source("00_config.R")
#   source("STEP_D07_matrix_transforms.R")
#   variants <- generate_all_variants(smat)
#   # variants is a named list: $raw, $distcorr, $localnorm, ...
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================================
# MASTER FUNCTION — generate all enabled variants
# ============================================================================

generate_all_variants <- function(smat, cfg = CFG) {
  N <- nrow(smat)
  cat("Matrix transforms v9.2 | N =", N, "\n")

  variants <- list(raw = smat)

  if (cfg$MTX_DISTCORR_ENABLED) {
    cat("  A1: Distance correction...\n")
    variants$distcorr <- mtx_distance_correction(smat)
  }

  if (cfg$MTX_LOCALNORM_ENABLED) {
    cat("  A2: Local contrast normalization...\n")
    variants$localnorm <- mtx_local_contrast(smat, window = cfg$MTX_LOCALNORM_WINDOW)
  }

  if (cfg$MTX_DENOISED_ENABLED) {
    cat("  A3: 2D denoising (median proxy)...\n")
    variants$denoised <- mtx_denoise_proxy(smat, lambda = cfg$MTX_DENOISED_LAMBDA)
  }

  if (cfg$MTX_RESIDBG_ENABLED) {
    cat("  A4: Background residualization...\n")
    variants$resid_bg <- mtx_background_residual(smat, rank_k = cfg$MTX_RESIDBG_RANK)
  }

  if (cfg$MTX_EDGE_ENABLED) {
    cat("  A5: Edge map...\n")
    variants$edge <- mtx_edge_map(smat)
  }

  if (cfg$MTX_SUPPORT_ENABLED) {
    cat("  A6: Support map...\n")
    variants$support <- mtx_support_map(smat, quantile_thresh = cfg$MTX_SUPPORT_QUANTILE)
  }

  cat("  Generated", length(variants), "matrix variants.\n")
  return(variants)
}


# ============================================================================
# A1: DISTANCE / BACKGROUND CORRECTION
# ============================================================================
# For each off-diagonal distance d, compute median sim at that distance,
# then subtract it.  Result: only ABOVE-expected signal remains.
# Family LD (smooth decay) → flat near zero.
# Inversions (elevated above expectation) → positive blocks.

mtx_distance_correction <- function(smat) {
  N <- nrow(smat)
  cat("    A1: vectorized distance correction (N=", N, ")...\n")

  # Compute median similarity at each distance d = 1..N-1
  # VECTORIZED: extract each diagonal with a single index operation
  max_d <- N - 1L
  dist_medians <- numeric(max_d)

  for (d in seq_len(max_d)) {
    # Extract d-th diagonal: indices (1,1+d), (2,2+d), ..., (N-d,N)
    diag_idx <- cbind(seq_len(N - d), seq_len(N - d) + d)
    diag_vals <- smat[diag_idx]
    diag_vals <- diag_vals[is.finite(diag_vals)]
    dist_medians[d] <- if (length(diag_vals) > 0) median(diag_vals) else 0
  }

  # Subtract expected — VECTORIZED: build correction matrix from diag offsets
  # The correction for cell (i,j) is dist_medians[|i-j|]
  # Build a Toeplitz-like matrix
  correction <- matrix(0, nrow = N, ncol = N)
  for (d in seq_len(max_d)) {
    # Set both diagonals at offset d
    idx_upper <- cbind(seq_len(N - d), seq_len(N - d) + d)
    idx_lower <- cbind(seq_len(N - d) + d, seq_len(N - d))
    correction[idx_upper] <- dist_medians[d]
    correction[idx_lower] <- dist_medians[d]
  }

  out <- smat - correction
  diag(out) <- 0
  out[!is.finite(out)] <- 0

  return(out)
}


# ============================================================================
# A2: LOCAL CONTRAST NORMALIZATION
# ============================================================================
# FAST VERSION: Instead of extracting a patch around every cell (O(N²×W²)),
# compute per-row running median and MAD, then per-column, and average.
# This is O(N²×W) — much faster for large matrices.

mtx_local_contrast <- function(smat, window = 50L) {
  N <- nrow(smat)
  cat("    A2: fast local contrast (N=", N, ", window=", window, ")...\n")

  # Step 1: Per-row running median and MAD
  row_med <- matrix(0, N, N)
  row_mad <- matrix(1, N, N)

  for (i in seq_len(N)) {
    row_vals <- smat[i, ]
    # Running median across columns
    for (j in seq_len(N)) {
      j_lo <- max(1L, j - window)
      j_hi <- min(N, j + window)
      patch <- row_vals[j_lo:j_hi]
      patch <- patch[is.finite(patch)]
      if (length(patch) >= 5) {
        row_med[i, j] <- median(patch)
        row_mad[i, j] <- max(mad(patch, constant = 1.4826), 1e-10)
      }
    }
    if (i %% 1000 == 0) cat("    row", i, "/", N, "\n")
  }

  # Step 2: Per-column running median and MAD
  col_med <- matrix(0, N, N)
  col_mad <- matrix(1, N, N)

  for (j in seq_len(N)) {
    col_vals <- smat[, j]
    for (i in seq_len(N)) {
      i_lo <- max(1L, i - window)
      i_hi <- min(N, i + window)
      patch <- col_vals[i_lo:i_hi]
      patch <- patch[is.finite(patch)]
      if (length(patch) >= 5) {
        col_med[i, j] <- median(patch)
        col_mad[i, j] <- max(mad(patch, constant = 1.4826), 1e-10)
      }
    }
  }

  # Step 3: Average the two directions, compute z-score
  local_med <- (row_med + col_med) / 2
  local_mad <- (row_mad + col_mad) / 2

  out <- (smat - local_med) / local_mad
  out[!is.finite(out)] <- 0

  # Symmetrize
  out <- (out + t(out)) / 2

  return(out)
}


# ============================================================================
# A3: 2D DENOISING — Practical proxy (iterated median filter)
# ============================================================================
# True fused lasso needs gfl/flsa — may not be on HPC.
# Proxy: iterated 2D median filter.  Preserves edges while smoothing noise.
# For real fused lasso, see mtx_denoise_fused_lasso() below.

mtx_denoise_proxy <- function(smat, lambda = 1.0) {
  N <- nrow(smat)
  out <- smat

  # Number of iterations proportional to lambda
  n_iter <- max(1L, round(lambda * 3))
  # Kernel size proportional to lambda
  kernel_half <- max(1L, round(lambda * 2))

  for (iter in seq_len(n_iter)) {
    prev <- out
    for (i in seq_len(N)) {
      for (j in i:N) {
        i_lo <- max(1L, i - kernel_half)
        i_hi <- min(N, i + kernel_half)
        j_lo <- max(1L, j - kernel_half)
        j_hi <- min(N, j + kernel_half)

        patch <- prev[i_lo:i_hi, j_lo:j_hi]
        patch_vals <- patch[is.finite(patch)]
        if (length(patch_vals) > 0) {
          med_val <- median(patch_vals)
          out[i, j] <- med_val
          out[j, i] <- med_val
        }
      }
    }
    cat("    iteration", iter, "/", n_iter, "\n")
  }

  return(out)
}


# ---- True fused lasso (if flsa or gfl available) ----
# Call this instead of the proxy if you have the package
mtx_denoise_fused_lasso <- function(smat, lambda = 1.0) {
  if (!requireNamespace("flsa", quietly = TRUE)) {
    cat("  flsa not available. Use mtx_denoise_proxy() instead.\n")
    cat("  Or install: install.packages('flsa')\n")
    return(mtx_denoise_proxy(smat, lambda))
  }

  N <- nrow(smat)

  # flsa operates on vectors — reshape matrix, apply, reshape back
  # For symmetric matrix, operate on upper triangle
  vec <- smat[upper.tri(smat, diag = TRUE)]
  n_elem <- length(vec)

  # Build adjacency for 2D grid (upper triangle mapped to 1D)
  # This is complex for arbitrary triangles — use the simpler row-wise approach
  # Apply 1D fused lasso to each row, then average with column-wise result

  out <- smat
  for (i in seq_len(N)) {
    row_vals <- smat[i, ]
    # 1D fused lasso on each row
    result <- tryCatch(
      flsa::flsa(row_vals, lambda2 = lambda),
      error = function(e) row_vals
    )
    out[i, ] <- result
  }

  # Symmetrize
  out <- (out + t(out)) / 2

  return(out)
}


# ============================================================================
# A4: BROAD BACKGROUND SUBTRACTION (low-rank residual)
# ============================================================================
# SVD decomposition: top-k singular values capture broad smooth structure.
# Residual = sim_raw - low_rank → contains local block structure.

mtx_background_residual <- function(smat, rank_k = 3L) {
  N <- nrow(smat)

  # Handle NAs
  smat_clean <- smat
  smat_clean[!is.finite(smat_clean)] <- median(smat, na.rm = TRUE)

  # SVD
  sv <- svd(smat_clean, nu = rank_k, nv = rank_k)

  # Low-rank reconstruction
  low_rank <- sv$u[, 1:rank_k, drop = FALSE] %*%
              diag(sv$d[1:rank_k], nrow = rank_k) %*%
              t(sv$v[, 1:rank_k, drop = FALSE])

  # Residual
  residual <- smat_clean - low_rank

  # Report variance explained
  total_var <- sum(sv$d^2)
  explained <- sum(sv$d[1:rank_k]^2) / total_var * 100
  cat("    SVD rank", rank_k, ": explains", round(explained, 1), "% of variance\n")
  cat("    Residual contains local block structure\n")

  return(residual)
}


# ---- Alternative: large-kernel median subtraction ----
mtx_background_smooth <- function(smat, kernel = 200L) {
  N <- nrow(smat)
  half <- kernel %/% 2

  # Compute smoothed background
  bg <- matrix(0, nrow = N, ncol = N)
  for (i in seq_len(N)) {
    for (j in i:N) {
      i_lo <- max(1L, i - half)
      i_hi <- min(N, i + half)
      j_lo <- max(1L, j - half)
      j_hi <- min(N, j + half)

      patch <- smat[i_lo:i_hi, j_lo:j_hi]
      vals <- patch[is.finite(patch)]
      med <- if (length(vals) > 0) median(vals) else 0
      bg[i, j] <- med
      bg[j, i] <- med
    }
    if (i %% 500 == 0) cat("    row", i, "/", N, "\n")
  }

  return(smat - bg)
}


# ============================================================================
# A5: EDGE / BOUNDARY MAP
# ============================================================================
# Gradient magnitude at each cell.
# High values = block boundaries.  Low values = interiors and background.

mtx_edge_map <- function(smat) {
  N <- nrow(smat)
  edge <- matrix(0, nrow = N, ncol = N)

  for (i in 2:(N - 1)) {
    for (j in 2:(N - 1)) {
      # Horizontal gradient (Sobel-like)
      gx <- (smat[i, j + 1] - smat[i, j - 1]) / 2

      # Vertical gradient
      gy <- (smat[i + 1, j] - smat[i - 1, j]) / 2

      # Handle NAs
      if (!is.finite(gx)) gx <- 0
      if (!is.finite(gy)) gy <- 0

      edge[i, j] <- sqrt(gx^2 + gy^2)
    }
  }

  # Symmetrize (gradient magnitude should be symmetric for symmetric matrix)
  edge <- (edge + t(edge)) / 2

  return(edge)
}


# ============================================================================
# A6: THRESHOLDED SUPPORT MAP
# ============================================================================
# Adaptive threshold using quantile.  Returns soft occupancy map.

mtx_support_map <- function(smat, quantile_thresh = 0.75) {
  N <- nrow(smat)

  # Compute threshold
  all_vals <- smat[upper.tri(smat)]
  all_vals <- all_vals[is.finite(all_vals)]
  threshold <- quantile(all_vals, probs = quantile_thresh, na.rm = TRUE)

  cat("    Support threshold (q", quantile_thresh, "):", round(threshold, 4), "\n")

  # Soft support: max(0, (sim - threshold)) / (1 - threshold)
  out <- pmax(0, smat - threshold) / max(1 - threshold, 1e-10)
  out[!is.finite(out)] <- 0

  return(out)
}


# ============================================================================
# SAVE ALL VARIANTS
# ============================================================================

save_matrix_variants <- function(variants, chr, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  for (name in names(variants)) {
    if (name == "raw") next  # don't re-save the original
    outfile <- file.path(outdir, sprintf("%s_sim_%s.rds", chr, name))
    saveRDS(variants[[name]], outfile)
    cat("  Saved:", basename(outfile), "\n")
  }
}
