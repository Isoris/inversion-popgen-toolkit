# ============================================================================
# STEP_D08b_local_pca.R — Local PCA / ICA per candidate
# ============================================================================
# For each candidate, extract the sub-matrix and run PCA (and optionally
# ICA) to check whether samples separate into clean genotype clusters.
#
# Inversion → PC1 shows 3 discrete clusters (INV/INV, INV/non, non/non)
#           → high silhouette score, trimodal distribution
# Family LD → PC1 shows a gradient, no discrete clusters
#           → low silhouette score, unimodal distribution
#
# ICA (if fastICA available): looks for the maximally non-Gaussian
# component. A trimodal IC is strong inversion evidence even when
# family structure dominates the PCA.
#
# Returns data.frame with columns:
#   candidate_id, pc1_silhouette, pc1_gap_k, pc1_trimodality,
#   ica_best_kurtosis, ica_trimodality
# ============================================================================

compute_local_pca <- function(candidates, smat, do_ica=TRUE) {
  # Check if fastICA is available
  has_ica <- do_ica && requireNamespace("fastICA", quietly=TRUE)
  if (do_ica && !has_ica) {
    cat("  fastICA not available. Skipping ICA. Install with:\n")
    cat("  install.packages('fastICA')\n")
  }

  results <- list()

  for (i in seq_len(nrow(candidates))) {
    si <- candidates$start[i]
    ei <- candidates$end[i]
    width <- ei - si + 1
    cid <- candidates$candidate_id[i]

    if (width < 10) {
      results[[i]] <- make_empty_pca_row(cid, has_ica)
      next
    }

    # Extract sub-matrix: windows as rows, windows as columns
    # Each column = one window's similarity profile to all other windows
    sub <- smat[si:ei, si:ei]
    n <- nrow(sub)

    # Replace NAs with column medians
    for (j in seq_len(ncol(sub))) {
      na_idx <- is.na(sub[, j]) | !is.finite(sub[, j])
      if (any(na_idx)) sub[na_idx, j] <- median(sub[!na_idx, j], na.rm=TRUE)
    }

    # If still NAs, replace with 0
    sub[is.na(sub)] <- 0

    # ---- PCA on the sub-matrix ----
    # Transpose: we want to cluster COLUMNS (windows) based on their
    # similarity profiles. Each column becomes a "sample" in PCA space.
    # But actually — what we really want is to cluster the SAMPLES
    # based on their behavior in this region. The sim_mat doesn't
    # directly give us per-sample scores.
    #
    # However, the sim_mat encodes sample information implicitly:
    # two windows have high similarity if the SAME samples drive
    # the local PCA in similar ways. So PCA on the window×window
    # sim sub-matrix captures the dominant mode of variation among
    # windows — which IS the genotype structure.
    #
    # PC1 of the sub-matrix separates windows into groups that have
    # different similarity patterns → different genotype compositions.

    pc <- tryCatch({
      prcomp(sub, center=TRUE, scale.=FALSE, rank.=5)
    }, error=function(e) NULL)

    if (is.null(pc) || ncol(pc$x) < 2) {
      results[[i]] <- make_empty_pca_row(cid, has_ica)
      next
    }

    pc1 <- pc$x[, 1]
    pc2 <- pc$x[, 2]
    var_explained_1 <- summary(pc)$importance[2, 1]  # proportion
    var_explained_2 <- summary(pc)$importance[2, 2]

    # ---- Silhouette score for k=3 clustering on PC1 ----
    sil_k3 <- compute_silhouette_k(pc1, k=3)
    sil_k2 <- compute_silhouette_k(pc1, k=2)

    # Best k by gap-like heuristic
    best_k <- if (sil_k3 > sil_k2 && sil_k3 > 0.3) 3L
              else if (sil_k2 > 0.3) 2L
              else 1L

    # ---- Trimodality of PC1 ----
    # Hartigan's dip test approximation: check if the distribution
    # of PC1 values has multiple modes
    trimod <- estimate_trimodality(pc1)

    # ---- ICA (if available) ----
    ica_kurtosis  <- NA_real_
    ica_trimod    <- NA_real_

    if (has_ica && n >= 10) {
      ica_result <- tryCatch({
        fastICA::fastICA(sub, n.comp=min(3, n-1), alg.typ="parallel",
                         fun="logcosh", verbose=FALSE)
      }, error=function(e) NULL)

      if (!is.null(ica_result) && !is.null(ica_result$S)) {
        # Find the component with highest kurtosis (most non-Gaussian)
        n_comp <- ncol(ica_result$S)
        kurtoses <- numeric(n_comp)
        trimods  <- numeric(n_comp)
        for (j in seq_len(n_comp)) {
          ic <- ica_result$S[, j]
          kurtoses[j] <- compute_kurtosis(ic)
          trimods[j]  <- estimate_trimodality(ic)
        }
        best_ic <- which.max(abs(kurtoses))
        ica_kurtosis <- kurtoses[best_ic]
        ica_trimod   <- trimods[best_ic]
      }
    }

    results[[i]] <- data.frame(
      candidate_id    = cid,
      pc1_var_expl    = round(var_explained_1, 4),
      pc2_var_expl    = round(var_explained_2, 4),
      pc1_silhouette  = round(max(sil_k2, sil_k3), 4),
      pc1_best_k      = best_k,
      pc1_trimodality = round(trimod, 4),
      ica_kurtosis    = round(ica_kurtosis, 4),
      ica_trimodality = round(ica_trimod, 4),
      stringsAsFactors = FALSE
    )
  }

  return(do.call(rbind, results))
}


# ---- Empty row ----
make_empty_pca_row <- function(cid, has_ica) {
  data.frame(
    candidate_id    = cid,
    pc1_var_expl    = NA_real_,
    pc2_var_expl    = NA_real_,
    pc1_silhouette  = NA_real_,
    pc1_best_k      = NA_integer_,
    pc1_trimodality = NA_real_,
    ica_kurtosis    = NA_real_,
    ica_trimodality = NA_real_,
    stringsAsFactors = FALSE
  )
}


# ---- Silhouette score for k-means with k clusters ----
compute_silhouette_k <- function(x, k=3) {
  n <- length(x)
  if (n < k + 1) return(0)

  # k-means on 1D data
  km <- tryCatch({
    kmeans(x, centers=k, nstart=10)
  }, error=function(e) NULL)

  if (is.null(km)) return(0)

  # Simplified silhouette: for each point, compute
  # (b - a) / max(a, b) where a = mean dist to own cluster,
  # b = mean dist to nearest other cluster
  labels <- km$cluster
  sil_vals <- numeric(n)

  for (i in seq_len(n)) {
    own   <- labels[i]
    own_members <- which(labels == own)
    own_members <- setdiff(own_members, i)

    if (length(own_members) == 0) { sil_vals[i] <- 0; next }

    a <- mean(abs(x[i] - x[own_members]))

    # Distance to nearest other cluster
    other_clusters <- setdiff(unique(labels), own)
    if (length(other_clusters) == 0) { sil_vals[i] <- 0; next }

    b_vals <- numeric(length(other_clusters))
    for (j in seq_along(other_clusters)) {
      members <- which(labels == other_clusters[j])
      b_vals[j] <- mean(abs(x[i] - x[members]))
    }
    b <- min(b_vals)

    sil_vals[i] <- (b - a) / max(a, b, 1e-10)
  }

  return(mean(sil_vals))
}


# ---- Trimodality estimate ----
# Simple heuristic: fit a density, count peaks
estimate_trimodality <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 10) return(0)

  # Kernel density estimation
  d <- density(x, n=256)

  # Count local maxima in the density
  peaks <- 0
  for (i in 2:(length(d$y) - 1)) {
    if (d$y[i] > d$y[i-1] && d$y[i] > d$y[i+1]) {
      peaks <- peaks + 1
    }
  }

  # Normalize: 3+ peaks → high trimodality score
  # Use proportion of peaks relative to expected for uniform
  score <- min(1, (peaks - 1) / 2)  # 1 peak → 0, 3 peaks → 1
  return(max(0, score))
}


# ---- Excess kurtosis ----
compute_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4) return(0)
  mu <- mean(x)
  s  <- sd(x)
  if (s == 0) return(0)
  m4 <- mean((x - mu)^4)
  return(m4 / s^4 - 3)  # excess kurtosis
}
