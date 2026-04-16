# ============================================================================
# STEP_D05_ghsl_stability.R — GHSL partition stability from sim_mat
# ============================================================================
# For each candidate, measure whether the sample partition (which samples
# are in each genotype cluster) is consistent across windows within the
# interval.
#
# Note 2026-04-17: Only compute_ghsl_from_simmat() is actually wired
# into the pipeline — run_all.R calls it directly (Phase 7c). The
# file-based compute_ghsl_stability() below is an unused stub from the
# design phase; it was never completed because the sim_mat fallback
# turned out to produce sufficient signal for D10 and the direct
# GHSL-file reader was not needed. Kept here for potential future use
# if per-window band assignments from phase_2/2e_ghsl/STEP_C04_snake3_ghsl_v5.R
# become worth reading directly.
#
# compute_ghsl_from_simmat() returns a data.frame with columns:
#   candidate_id, partition_stability, n_consistent_windows,
#   dominant_partition_frac, partition_entropy
# These are consumed by C01d at the D10 dimension via the scoring_table
# merge in run_all.R Phase 8.
# ============================================================================

# ---- UNUSED STUB — see header note ----
#
# File-based GHSL reader. Never finished because the sim_mat fallback
# (compute_ghsl_from_simmat below) turned out to be sufficient. Do not
# wire this into run_all.R as-is — the band-file format-parsing is a
# placeholder and returns all-NA rows. If you revive it, the format
# should match what STEP_C04_snake3_ghsl_v5.R writes.
compute_ghsl_stability <- function(candidates, ghsl_dir, chr) {
  # Look for GHSL band assignment files
  band_file <- file.path(ghsl_dir, sprintf("%s_band_assignments.tsv", chr))

  if (!file.exists(band_file)) {
    # Try alternative naming patterns
    candidates_files <- list.files(ghsl_dir, pattern=chr, full.names=TRUE)
    band_files <- grep("band|assign|ghsl", candidates_files, value=TRUE)
    if (length(band_files) > 0) {
      band_file <- band_files[1]
    } else {
      cat("  No GHSL band files found. Falling back to sim_mat method.\n")
      return(NULL)  # Caller will use sim_mat fallback
    }
  }

  cat("  Reading GHSL bands from:", basename(band_file), "\n")
  bands <- read.delim(band_file, stringsAsFactors=FALSE)

  # Expected format: rows = samples, columns = windows (or vice versa)
  # Adapt to actual format...
  # For now, return placeholder — actual GHSL format parsing depends on
  # the exact output of STEP_C04_snake3_ghsl_v5.R

  results <- list()
  for (i in seq_len(nrow(candidates))) {
    cid <- candidates$candidate_id[i]
    si <- candidates$start[i]
    ei <- candidates$end[i]

    # Extract band assignments for windows in this interval
    # TODO: map window indices to GHSL output indices
    results[[i]] <- data.frame(
      candidate_id           = cid,
      partition_stability    = NA_real_,
      n_consistent_windows   = NA_integer_,
      dominant_partition_frac = NA_real_,
      partition_entropy      = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  return(do.call(rbind, results))
}


# ---- Fallback: compute partition stability from sim_mat directly ----
compute_ghsl_from_simmat <- function(candidates, smat) {
  N <- nrow(smat)
  results <- list()

  for (i in seq_len(nrow(candidates))) {
    si <- candidates$start[i]
    ei <- candidates$end[i]
    width <- ei - si + 1
    cid <- candidates$candidate_id[i]

    if (width < 6) {
      results[[i]] <- data.frame(
        candidate_id           = cid,
        partition_stability    = NA_real_,
        n_consistent_windows   = NA_integer_,
        dominant_partition_frac = NA_real_,
        partition_entropy      = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Strategy: for each window (column) in the block, the column of the
    # sim_mat gives a "similarity profile" — how similar this window is
    # to every other window. Windows within the same genotype class should
    # have correlated profiles.
    #
    # Cluster the columns of the sub-matrix and measure consistency.

    sub <- smat[si:ei, si:ei]

    # Compute pairwise correlation between columns
    # Each column = one window's similarity profile to all other windows in block
    n_cols <- ncol(sub)

    if (n_cols > 200) {
      # Subsample for speed
      col_idx <- sort(sample(n_cols, 200))
      sub <- sub[col_idx, col_idx]
      n_cols <- 200
    }

    # Correlation matrix of columns
    col_cors <- cor(sub, use="pairwise.complete.obs")
    col_cors[!is.finite(col_cors)] <- 0

    # Simple k=3 clustering on the correlation matrix
    # Using hierarchical clustering for robustness
    if (n_cols < 4) {
      results[[i]] <- data.frame(
        candidate_id           = cid,
        partition_stability    = NA_real_,
        n_consistent_windows   = NA_integer_,
        dominant_partition_frac = NA_real_,
        partition_entropy      = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    d <- as.dist(1 - col_cors)
    hc <- hclust(d, method="ward.D2")
    k3 <- cutree(hc, k=min(3, n_cols))

    # Partition stability: for each pair of adjacent windows,
    # do they stay in the same cluster? High stability = consistent partition
    n_same <- 0
    n_pairs <- 0
    for (j in 1:(n_cols - 1)) {
      n_pairs <- n_pairs + 1
      if (k3[j] == k3[j + 1]) n_same <- n_same + 1
    }
    stability <- if (n_pairs > 0) n_same / n_pairs else NA

    # Dominant partition: fraction of windows in the largest cluster
    tab <- table(k3)
    dominant_frac <- max(tab) / sum(tab)

    # Entropy of partition
    probs <- tab / sum(tab)
    entropy <- -sum(probs * log2(probs + 1e-10))

    # How many windows are "consistent" — in a run of same-cluster
    runs <- rle(k3)
    max_run <- max(runs$lengths)
    n_consistent <- max_run

    results[[i]] <- data.frame(
      candidate_id           = cid,
      partition_stability    = round(stability, 4),
      n_consistent_windows   = n_consistent,
      dominant_partition_frac = round(dominant_frac, 4),
      partition_entropy      = round(entropy, 4),
      stringsAsFactors = FALSE
    )
  }

  return(do.call(rbind, results))
}
