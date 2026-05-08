#!/usr/bin/env Rscript
# =============================================================================
# cheat26_kin_pruned_retention.R — Kin-pruned signal retention
#
# BIOLOGY:
#   With Ne~20, family structure drives many local PCA signals. If a block
#   disappears when we restrict to 81 unrelated samples → FAMILY_DRIVEN.
#   If it persists → STRUCTURAL (real inversion).
#
# INPUT:  pruned_samples.txt (81 IDs), BEAGLE dosage, precomp
# OUTPUT: retention_ratio, classification per core
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
STRUCTURAL_THRESHOLD <- 0.80
MIXED_THRESHOLD      <- 0.50
MIN_PRUNED_SAMPLES   <- 20L
MIN_WINDOWS          <- 5L

# ── Compute inv_likeness for a sample subset ──────────────────────────

compute_subset_inv_likeness <- function(dosage_matrix, sample_ids,
                                         window_starts, window_ends,
                                         window_size = 100L) {
  avail <- intersect(sample_ids, colnames(dosage_matrix))
  if (length(avail) < MIN_PRUNED_SAMPLES) return(rep(NA_real_, length(window_starts)))

  sub_dosage <- dosage_matrix[, avail, drop = FALSE]
  n_windows <- length(window_starts)
  il_vals <- numeric(n_windows)

  for (wi in seq_len(n_windows)) {
    # Get SNPs in this window
    snp_idx <- which(seq_len(nrow(sub_dosage)) >= window_starts[wi] &
                      seq_len(nrow(sub_dosage)) <= window_ends[wi])
    if (length(snp_idx) < 10) { il_vals[wi] <- NA_real_; next }

    # PCA on subset
    chunk <- t(sub_dosage[snp_idx, , drop = FALSE])
    chunk <- chunk[, apply(chunk, 2, var, na.rm = TRUE) > 0, drop = FALSE]
    if (ncol(chunk) < 5) { il_vals[wi] <- NA_real_; next }

    pca <- tryCatch(prcomp(chunk, center = TRUE, scale. = FALSE),
                     error = function(e) NULL)
    if (is.null(pca)) { il_vals[wi] <- NA_real_; next }

    pve1 <- pca$sdev[1]^2 / sum(pca$sdev^2)
    pc1  <- pca$x[, 1]

    # Simple inv_likeness proxy: PVE1 × bimodality indicator
    # Dip test for trimodality (if diptest available)
    dip_p <- tryCatch({
      if (requireNamespace("diptest", quietly = TRUE))
        diptest::dip.test(pc1)$p.value
      else 0.5
    }, error = function(e) 0.5)

    # inv_likeness = PVE1 if dip_p < 0.1, else PVE1 * 0.5
    il_vals[wi] <- pve1 * (if (dip_p < 0.1) 1.0 else 0.5)
  }
  il_vals
}

# ── Compute retention ratio ───────────────────────────────────────────

compute_retention_ratio <- function(il_full, il_pruned) {
  valid <- is.finite(il_full) & is.finite(il_pruned)
  if (sum(valid) < MIN_WINDOWS) return(NA_real_)

  mean_full   <- mean(il_full[valid], na.rm = TRUE)
  mean_pruned <- mean(il_pruned[valid], na.rm = TRUE)

  if (mean_full <= 0) return(NA_real_)
  mean_pruned / mean_full
}

# ── Classify retention ────────────────────────────────────────────────

classify_retention <- function(retention_ratio) {
  if (is.na(retention_ratio)) return("UNTESTED")
  if (retention_ratio > STRUCTURAL_THRESHOLD) return("STRUCTURAL")
  if (retention_ratio > MIXED_THRESHOLD) return("MIXED")
  "FAMILY_DRIVEN"
}

# ── Search mode ────────────────────────────────────────────────────────

search_kin_pruned <- function(chr, zone_start, zone_end, dt = NULL,
                               pruned_ids = NULL, ...) {
  empty <- data.table(method = "kin_pruned_retention", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(dt) || nrow(dt) == 0 || is.null(pruned_ids)) return(empty)

  # Quick check: does inv_likeness exist for pruned subset?
  region_dt <- dt[start_bp >= zone_start & end_bp <= zone_end]
  if (nrow(region_dt) < MIN_WINDOWS) return(empty)

  if (!"inv_likeness" %in% names(region_dt)) return(empty)

  full_il <- mean(region_dt$inv_likeness, na.rm = TRUE)

  # If we have pruned inv_likeness stored
  if ("inv_likeness_pruned" %in% names(region_dt)) {
    pruned_il <- mean(region_dt$inv_likeness_pruned, na.rm = TRUE)
    ratio <- if (full_il > 0) pruned_il / full_il else NA
    sc <- if (!is.na(ratio)) pmin(1, ratio) else 0
    return(data.table(method = "kin_pruned_retention",
                       best_bp = as.integer((zone_start + zone_end) / 2),
                       score = round(sc, 3), is_precise = FALSE,
                       detail = paste0("retention=", round(ratio, 3))))
  }

  data.table(method = "kin_pruned_retention",
             best_bp = as.integer((zone_start + zone_end) / 2),
             score = 0.5, is_precise = FALSE,
             detail = "pruned_il_not_precomputed")
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat26 <- function(chr, candidate_start, candidate_end,
                         dt_full, pruned_ids = NULL,
                         dosage_matrix = NULL) {
  message("[cheat26] ", chr, ":", round(candidate_start/1e6,1), "-",
          round(candidate_end/1e6,1), " Mb")

  # Load pruned sample IDs
  if (is.null(pruned_ids)) {
    pruned_file <- Sys.getenv("PRUNED_SAMPLES", "")
    if (nzchar(pruned_file) && file.exists(pruned_file)) {
      pruned_ids <- readLines(pruned_file)
      message("[cheat26] Loaded ", length(pruned_ids), " pruned samples")
    } else {
      message("[cheat26] No pruned samples file")
      return(list(retention_ratio = NA_real_, status = "UNTESTED",
                  search_result = search_kin_pruned(chr, candidate_start,
                    candidate_end)))
    }
  }

  # Get full inv_likeness from precomp
  region_dt <- dt_full[start_bp >= candidate_start & end_bp <= candidate_end]
  if (nrow(region_dt) < MIN_WINDOWS) {
    message("[cheat26] Too few windows")
    return(list(retention_ratio = NA_real_, status = "UNTESTED",
                search_result = search_kin_pruned(chr, candidate_start,
                  candidate_end)))
  }

  il_full <- if ("inv_likeness" %in% names(region_dt))
    region_dt$inv_likeness else rep(NA_real_, nrow(region_dt))

  # Compute pruned inv_likeness
  if (!is.null(dosage_matrix)) {
    # Use dosage matrix to recompute PCA on pruned subset
    il_pruned <- compute_subset_inv_likeness(
      dosage_matrix, pruned_ids,
      seq_len(nrow(region_dt)),  # simplified window indices
      seq_len(nrow(region_dt)))
  } else if ("inv_likeness_pruned" %in% names(region_dt)) {
    il_pruned <- region_dt$inv_likeness_pruned
  } else {
    message("[cheat26] Cannot compute pruned inv_likeness without dosage or precomp")
    return(list(retention_ratio = NA_real_, status = "UNTESTED",
                search_result = search_kin_pruned(chr, candidate_start,
                  candidate_end, dt_full, pruned_ids)))
  }

  ratio <- compute_retention_ratio(il_full, il_pruned)
  status <- classify_retention(ratio)

  message("[cheat26] Retention: ", round(ratio, 3), " → ", status)
  message("[cheat26] Full mean IL: ", round(mean(il_full, na.rm = TRUE), 4),
          " | Pruned mean IL: ", round(mean(il_pruned, na.rm = TRUE), 4))

  list(retention_ratio = round(ratio, 4), status = status,
       full_mean_il = round(mean(il_full, na.rm = TRUE), 4),
       pruned_mean_il = round(mean(il_pruned, na.rm = TRUE), 4),
       n_pruned_samples = length(pruned_ids),
       search_result = search_kin_pruned(chr, candidate_start, candidate_end,
                                          dt_full, pruned_ids))
}
