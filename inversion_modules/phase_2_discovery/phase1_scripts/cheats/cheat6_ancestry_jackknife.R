#!/usr/bin/env Rscript
# =============================================================================
# cheat6_ancestry_jackknife.R — Leave-one-family-out Fst robustness
#
# For each boundary (or window), systematically remove one ancestry group
# at a time and measure how the Fst changes.
#
# max_delta_Q LOW  + n_families HIGH  → INVERSION (robust, multi-family)
# max_delta_Q HIGH + n_families = 1   → FAMILY LD (fragile, single-family)
#
# Side product: per-family inversion carrier probability.
# Cross-check: if removing a family makes boundary SHARPER → recombinants.
#
# REQUIRES: load_bridge.R sourced (get_region_stats, reg)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

DELTA_THRESH <- 0.02  # Fst drop must exceed this to count as "contributing"

#' Run ancestry jackknife for one region
#' @param chr Chromosome
#' @param start_bp, end_bp Region bounds
#' @param pc1_bands Named integer vector (sample → band 1/2/3)
#' @param ancestry_groups Named list of character vectors (Q1→samples, Q2→...)
#' @return data.table with per-group jackknife results
ancestry_jackknife <- function(chr, start_bp, end_bp, pc1_bands, ancestry_groups) {
  band1 <- names(pc1_bands)[pc1_bands == 1]
  band3 <- names(pc1_bands)[pc1_bands == 3]

  if (length(band1) < 5 || length(band3) < 5) {
    return(data.table(group_removed = "none", fst = NA_real_,
                       delta_fst = NA_real_, direction = "none"))
  }

  # Full Fst (all samples)
  fst_full <- tryCatch({
    s <- get_region_stats(chr, start_bp, end_bp, what = "Fst",
                           groups = list(b1 = band1, b3 = band3))
    s$Fst$b1_vs_b3
  }, error = function(e) NA_real_)

  if (!is.finite(fst_full)) {
    return(data.table(group_removed = "none", fst = fst_full,
                       delta_fst = NA_real_, direction = "none"))
  }

  rows <- list()
  rows[[1]] <- data.table(group_removed = "none", fst = round(fst_full, 4),
                            delta_fst = 0, direction = "baseline", n_removed = 0L)

  for (gname in names(ancestry_groups)) {
    remove_samples <- ancestry_groups[[gname]]
    b1_kept <- setdiff(band1, remove_samples)
    b3_kept <- setdiff(band3, remove_samples)

    if (length(b1_kept) < 3 || length(b3_kept) < 3) {
      rows[[length(rows) + 1]] <- data.table(
        group_removed = gname, fst = NA_real_, delta_fst = NA_real_,
        direction = "insufficient", n_removed = length(remove_samples))
      next
    }

    fst_without <- tryCatch({
      s <- get_region_stats(chr, start_bp, end_bp, what = "Fst",
                             groups = list(b1 = b1_kept, b3 = b3_kept))
      s$Fst$b1_vs_b3
    }, error = function(e) NA_real_)

    delta <- if (is.finite(fst_without)) fst_full - fst_without else NA_real_
    direction <- if (is.finite(delta)) {
      if (delta > DELTA_THRESH) "drops"      # removing this family reduces Fst
      else if (delta < -DELTA_THRESH) "rises" # removing sharpens → recombinants
      else "stable"
    } else "unknown"

    rows[[length(rows) + 1]] <- data.table(
      group_removed = gname, fst = round(fst_without, 4),
      delta_fst = round(delta, 4), direction = direction,
      n_removed = length(intersect(remove_samples, c(band1, band3))))
  }

  out <- rbindlist(rows)

  # Summary metrics
  real_deltas <- out[direction %in% c("drops", "rises", "stable") & group_removed != "none"]
  if (nrow(real_deltas) > 0) {
    out[, max_delta := max(abs(real_deltas$delta_fst), na.rm = TRUE)]
    out[, n_contributing := sum(real_deltas$direction == "drops")]
    out[, n_sharpening  := sum(real_deltas$direction == "rises")]

    # Interpretation
    md <- max(abs(real_deltas$delta_fst), na.rm = TRUE)
    nc <- sum(real_deltas$direction == "drops")
    out[, jackknife_verdict := if (md < DELTA_THRESH) "robust_multi_family"
                               else if (nc == 1) "single_family_fragile"
                               else if (nc <= 2) "few_family"
                               else "multi_family_contributing"]
  }

  out
}

#' Get ancestry groups from registry (K=8 NGSadmix groups)
get_ancestry_groups_from_registry <- function() {
  if (is.null(reg)) return(list())

  groups <- list()
  for (k in 1:8) {
    gid <- paste0("ancestry_K8_Q", k)
    if (reg$has_group(gid)) {
      members <- reg$get_group(gid)
      if (length(members) > 0) groups[[paste0("Q", k)]] <- members
    }
  }
  groups
}
