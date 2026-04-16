#!/usr/bin/env Rscript
# =============================================================================
# cheat12_theta_het_prior.R — Per-sample θ_P heterozygosity prior
#
# WHAT:
#   Uses per-sample windowed θ_P (nucleotide diversity from ANGSD saf2theta)
#   as an independent genotype classifier. HET samples carry two divergent
#   haplotypes inside an inversion → elevated θ_P. HOM samples (either
#   arrangement) have low θ_P because both haplotypes are the same.
#
# WHY:
#   Completely independent of local PCA eigenvectors, SV callers, Fst,
#   phasing, and all other cheats. Computed from ANGSD SAF→SFS→theta,
#   a fundamentally different statistical framework (site frequency spectrum
#   vs genotype likelihood PCA). If PCA says "HET" and θ_P says "high"
#   for the same sample → strong concordance.
#
# INPUT:
#   Per-sample thetaStat output from MODULE_3 (ANGSD windowed theta):
#     Window size: 50,000 bp, step: 10,000 bp
#     File pattern: {theta_dir}/{sample_id}.thetaStat.pestPG.gz
#     Columns: chr, wincenter, tW, tP, tF, tH, tL, Tajima, fuf, fud, fayH, zeng
#
#   Per-chromosome PC1 bands from precomp:
#     Fixed k=3 bands: band1 (HOM_REF-like), band2 (HET-like), band3 (HOM_INV-like)
#
# OUTPUT per candidate region:
#   Per-sample: theta_class (HET_like / HOM_like / AMBIGUOUS), confidence (0-1)
#   Per-boundary: band separation score (gap between HET and HOM θ_P distributions)
#
# RESOLUTION:
#   θ_P windows are 50kb/10kb step = coarser than PCA windows (100 SNPs, ~2.5kb).
#   Adequate for genotype classification within 1-4 Mb candidate regions
#   (100-400 θ_P windows per candidate). NOT useful for fine boundary detection.
#   Confidence is proportional to band separation — narrow gap → low prior,
#   wide gap → strong prior.
#
# REQUIRES: load_bridge.R (for smap), precomp RDS (for band assignments)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# =============================================================================
# PARAMETERS
# =============================================================================

THETA_WIN_SIZE   <- 50000L   # ANGSD thetaStat window size
THETA_STEP_SIZE  <- 10000L   # ANGSD thetaStat step size
MIN_WINDOWS      <- 5L       # minimum θ_P windows in candidate region
MIN_SEPARATION   <- 0.5      # minimum z-score gap for confident classification
HIGH_SEPARATION  <- 1.5      # z-score gap for high-confidence classification

# =============================================================================
# LOAD PER-SAMPLE THETA FILES
# =============================================================================

#' Load one sample's thetaStat pestPG output for a chromosome
#' @param theta_dir Directory containing {sample_id}.thetaStat.pestPG.gz
#' @param sample_id CGA sample ID
#' @param chr Chromosome name
#' @return data.table with wincenter, tP columns (filtered to chr)
load_sample_theta <- function(theta_dir, sample_id, chr) {
  # Try multiple naming patterns
  patterns <- c(
    file.path(theta_dir, paste0(sample_id, ".thetaStat.pestPG.gz")),
    file.path(theta_dir, paste0(sample_id, ".", chr, ".thetaStat.pestPG.gz")),
    file.path(theta_dir, chr, paste0(sample_id, ".thetaStat.pestPG.gz"))
  )

  ff <- NULL
  for (p in patterns) {
    if (file.exists(p)) { ff <- p; break }
  }
  if (is.null(ff)) return(data.table())

  dt <- tryCatch(
    fread(ff, header = TRUE, select = c("Chr", "WinCenter", "tP")),
    error = function(e) data.table()
  )
  if (nrow(dt) == 0) return(data.table())

  # Rename and filter to target chromosome
  setnames(dt, c("Chr", "WinCenter", "tP"), c("chr", "wincenter", "tP"),
           skip_absent = TRUE)
  dt <- dt[chr == ..chr]
  dt[, sample_id := sample_id]
  dt
}

#' Load theta for all samples in a region
#' @param theta_dir Base directory
#' @param sample_ids Character vector of CGA sample IDs
#' @param chr Chromosome name
#' @param start_bp, end_bp Region bounds
#' @return data.table: sample_id, wincenter, tP (region-filtered)
load_region_theta <- function(theta_dir, sample_ids, chr, start_bp, end_bp) {
  all_dt <- list()

  for (sid in sample_ids) {
    dt <- load_sample_theta(theta_dir, sid, chr)
    if (nrow(dt) == 0) next

    # Filter to region (with some padding for edge windows)
    pad <- THETA_WIN_SIZE / 2
    dt <- dt[wincenter >= (start_bp - pad) & wincenter <= (end_bp + pad)]
    if (nrow(dt) < MIN_WINDOWS) next

    all_dt[[length(all_dt) + 1]] <- dt
  }

  if (length(all_dt) > 0) rbindlist(all_dt) else data.table()
}

# =============================================================================
# Z-SCORE NORMALIZATION (within-chromosome)
# =============================================================================

#' Compute per-sample z-scored θ_P within a chromosome
#' Z-score relative to the full chromosome distribution for that sample
#' @param theta_dir Base directory
#' @param sample_ids Character vector
#' @param chr Chromosome name
#' @param start_bp, end_bp Candidate region bounds
#' @return data.table: sample_id, wincenter, tP, tP_z
compute_theta_zscores <- function(theta_dir, sample_ids, chr, start_bp, end_bp) {
  # Load full-chromosome theta for z-scoring, then subset
  all_dt <- list()

  for (sid in sample_ids) {
    dt <- load_sample_theta(theta_dir, sid, chr)
    if (nrow(dt) < 20) next

    # Z-score within whole chromosome
    chr_mean <- mean(dt$tP, na.rm = TRUE)
    chr_sd   <- sd(dt$tP, na.rm = TRUE)
    if (!is.finite(chr_sd) || chr_sd == 0) next

    dt[, tP_z := (tP - chr_mean) / chr_sd]

    # Filter to candidate region
    pad <- THETA_WIN_SIZE / 2
    dt <- dt[wincenter >= (start_bp - pad) & wincenter <= (end_bp + pad)]
    if (nrow(dt) < MIN_WINDOWS) next

    all_dt[[length(all_dt) + 1]] <- dt
  }

  if (length(all_dt) > 0) rbindlist(all_dt) else data.table()
}

# =============================================================================
# BAND-STRATIFIED THETA ANALYSIS
# =============================================================================

#' Compute per-band mean θ_P z-score within a candidate region
#' @param theta_z data.table from compute_theta_zscores
#' @param band_assignments Named vector: sample_id → band (1, 2, 3)
#' @return data.table: band, mean_tP_z, sd_tP_z, n_samples, n_windows
compute_band_theta <- function(theta_z, band_assignments) {
  if (nrow(theta_z) == 0) return(data.table())

  theta_z[, band := band_assignments[sample_id]]
  theta_z <- theta_z[!is.na(band)]

  if (nrow(theta_z) == 0) return(data.table())

  # Per-sample mean z-score across all windows in region
  per_sample <- theta_z[, .(
    mean_tP_z = mean(tP_z, na.rm = TRUE),
    sd_tP_z = sd(tP_z, na.rm = TRUE),
    n_windows = .N
  ), by = .(sample_id, band)]

  # Per-band summary
  per_band <- per_sample[, .(
    band_mean_tP_z = mean(mean_tP_z, na.rm = TRUE),
    band_sd_tP_z = sd(mean_tP_z, na.rm = TRUE),
    n_samples = .N,
    total_windows = sum(n_windows)
  ), by = band]

  list(per_sample = per_sample, per_band = per_band)
}

# =============================================================================
# CLASSIFICATION ENGINE
# =============================================================================

#' Classify samples using θ_P heterozygosity prior
#' @param theta_z data.table from compute_theta_zscores
#' @param band_assignments Named vector: sample_id → band (1, 2, 3)
#' @return data.table: sample_id, theta_class, confidence, prior_weight,
#'   mean_tP_z, assigned_band, band_separation
classify_by_theta <- function(theta_z, band_assignments) {
  if (nrow(theta_z) == 0) {
    return(data.table(
      sample_id = character(0), theta_class = character(0),
      confidence = numeric(0), prior_weight = numeric(0)
    ))
  }

  band_result <- compute_band_theta(theta_z, band_assignments)
  per_sample <- band_result$per_sample
  per_band <- band_result$per_band

  if (nrow(per_band) < 2) {
    return(data.table(
      sample_id = per_sample$sample_id,
      theta_class = "AMBIGUOUS",
      confidence = 0, prior_weight = 0,
      mean_tP_z = per_sample$mean_tP_z,
      assigned_band = per_sample$band,
      band_separation = 0
    ))
  }

  # Find the HET-like band (highest mean θ_P z) and HOM-like bands (lowest)
  per_band <- per_band[order(band_mean_tP_z)]
  het_band <- per_band[which.max(band_mean_tP_z)]
  hom_band <- per_band[which.min(band_mean_tP_z)]

  # Band separation = gap between HET and HOM band means
  band_sep <- het_band$band_mean_tP_z - hom_band$band_mean_tP_z

  # For each sample: how close is its θ_P to the HET vs HOM band mean?
  het_mean <- het_band$band_mean_tP_z
  hom_mean <- hom_band$band_mean_tP_z
  midpoint <- (het_mean + hom_mean) / 2

  per_sample[, `:=`(
    dist_to_het = abs(mean_tP_z - het_mean),
    dist_to_hom = abs(mean_tP_z - hom_mean),
    above_midpoint = mean_tP_z > midpoint
  )]

  # Classify
  per_sample[, theta_class := fifelse(
    band_sep < MIN_SEPARATION, "AMBIGUOUS",
    fifelse(above_midpoint, "HET_like", "HOM_like")
  )]

  # Confidence: based on how far from midpoint AND how wide the band gap is
  per_sample[, raw_confidence := fifelse(
    band_sep < MIN_SEPARATION, 0,
    abs(mean_tP_z - midpoint) / (band_sep / 2)  # 1.0 = at band center
  )]
  per_sample[, confidence := pmin(1, pmax(0, raw_confidence))]

  # Prior weight: scales with band separation
  #   Separation < 0.5 → weight = 0 (bands overlap, θ_P uninformative)
  #   Separation 0.5-1.5 → weight increases linearly
  #   Separation > 1.5 → weight = 1 (strong prior)
  prior_weight <- pmin(1, pmax(0, (band_sep - MIN_SEPARATION) /
                                     (HIGH_SEPARATION - MIN_SEPARATION)))

  per_sample[, `:=`(
    prior_weight = prior_weight,
    band_separation = band_sep
  )]

  # Clean up
  per_sample[, c("dist_to_het", "dist_to_hom", "above_midpoint",
                  "raw_confidence") := NULL]

  per_sample
}

# =============================================================================
# SEARCH MODE: Scan within a boundary zone for θ_P transition
# =============================================================================

#' Find θ_P transition point within a boundary zone
#' The position where band-stratified θ_P lines converge or cross
#' @param theta_z data.table from compute_theta_zscores
#' @param band_assignments Named vector
#' @param zone_start, zone_end Zone bounds in bp
#' @return data.table(method, best_bp, score, is_precise, detail)
search_theta_transition <- function(theta_z, band_assignments, zone_start, zone_end) {
  if (nrow(theta_z) == 0) {
    return(data.table(method = "theta_transition", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data"))
  }

  theta_z[, band := band_assignments[sample_id]]
  theta_z <- theta_z[!is.na(band)]

  # Get unique window positions within zone
  zone_wins <- sort(unique(theta_z[wincenter >= zone_start &
                                     wincenter <= zone_end]$wincenter))

  if (length(zone_wins) < 3) {
    return(data.table(method = "theta_transition", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "too_few_windows"))
  }

  # At each window position, compute band separation
  separations <- numeric(length(zone_wins))
  for (wi in seq_along(zone_wins)) {
    wc <- zone_wins[wi]
    # Use ±2 windows for smoothing (±20 kb at 10kb step)
    nearby <- theta_z[abs(wincenter - wc) <= 2 * THETA_STEP_SIZE]
    if (nrow(nearby) < 10) { separations[wi] <- NA; next }

    band_means <- nearby[, .(m = mean(tP_z, na.rm = TRUE)), by = band]
    if (nrow(band_means) < 2) { separations[wi] <- NA; next }

    separations[wi] <- max(band_means$m) - min(band_means$m)
  }

  valid <- is.finite(separations)
  if (sum(valid) < 3) {
    return(data.table(method = "theta_transition", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "insufficient_valid"))
  }

  # Transition = position where separation changes most (steepest gradient)
  # = position where bands converge (boundary of inversion)
  diffs <- abs(diff(separations[valid]))
  if (length(diffs) == 0) {
    return(data.table(method = "theta_transition", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_gradient"))
  }

  best_idx <- which.max(diffs)
  best_bp <- as.integer(zone_wins[valid][best_idx])

  # Score based on how steep the transition is
  max_diff <- max(diffs, na.rm = TRUE)
  score <- pmin(1, max_diff / 1.0)  # saturates at diff=1.0 z-score units

  data.table(
    method = "theta_transition",
    best_bp = best_bp,
    score = round(score, 3),
    is_precise = FALSE,  # θ_P resolution is 10 kb at best
    detail = paste0("max_gradient=", round(max_diff, 3),
                    ",n_windows=", sum(valid))
  )
}

# =============================================================================
# CONVENIENCE: Full Cheat 12 for one candidate
# =============================================================================

#' Run Cheat 12 for a candidate inversion region
#' @param chr Chromosome
#' @param start_bp, end_bp Candidate region bounds
#' @param theta_dir Path to per-sample thetaStat output
#' @param sample_ids All sample IDs
#' @param band_assignments Named vector: sample_id → band (1/2/3)
#' @return list(classification, band_summary, search_result)
run_cheat12 <- function(chr, start_bp, end_bp, theta_dir,
                         sample_ids, band_assignments) {
  message("[cheat12] ", chr, ":", round(start_bp/1e6, 1), "-",
          round(end_bp/1e6, 1), " Mb (", length(sample_ids), " samples)")

  # Compute z-scored theta in the region
  theta_z <- compute_theta_zscores(theta_dir, sample_ids, chr, start_bp, end_bp)

  if (nrow(theta_z) == 0) {
    message("[cheat12] No theta data available")
    return(list(
      classification = data.table(),
      band_summary = data.table(),
      search_result = data.table(method = "theta_transition",
                                  best_bp = NA_integer_, score = 0,
                                  is_precise = FALSE, detail = "no_data")
    ))
  }

  n_samples_with_data <- length(unique(theta_z$sample_id))
  n_windows <- length(unique(theta_z$wincenter))
  message("[cheat12] Loaded: ", n_samples_with_data, " samples × ",
          n_windows, " windows")

  # Classify
  classification <- classify_by_theta(theta_z, band_assignments)

  # Band summary
  band_result <- compute_band_theta(theta_z, band_assignments)

  # Report
  if (nrow(classification) > 0) {
    n_het <- sum(classification$theta_class == "HET_like")
    n_hom <- sum(classification$theta_class == "HOM_like")
    n_amb <- sum(classification$theta_class == "AMBIGUOUS")
    sep <- if (nrow(classification) > 0) classification$band_separation[1] else 0
    pw <- if (nrow(classification) > 0) classification$prior_weight[1] else 0

    message("[cheat12] Classification: ", n_het, " HET_like, ",
            n_hom, " HOM_like, ", n_amb, " AMBIGUOUS")
    message("[cheat12] Band separation: ", round(sep, 2),
            " | Prior weight: ", round(pw, 2))
  }

  list(
    classification = classification,
    band_summary = band_result$per_band,
    search_result = search_theta_transition(theta_z, band_assignments,
                                             start_bp, end_bp)
  )
}

# =============================================================================
# CONCORDANCE WITH OTHER CHEATS
# =============================================================================

#' Check concordance between Cheat 12 (θ_P) and PCA band assignment
#' @param classification data.table from classify_by_theta
#' @return data.table with concordance summary
check_theta_pca_concordance <- function(classification) {
  if (nrow(classification) == 0) return(data.table())

  # Band 2 should be HET_like, Bands 1 and 3 should be HOM_like
  classification[, expected := fifelse(band == 2, "HET_like", "HOM_like")]
  classification[, concordant := (theta_class == expected) |
                                   theta_class == "AMBIGUOUS"]

  n_total <- nrow(classification[theta_class != "AMBIGUOUS"])
  n_concordant <- sum(classification$concordant[classification$theta_class != "AMBIGUOUS"])
  concordance_rate <- if (n_total > 0) n_concordant / n_total else NA_real_

  # Fisher test: is θ_P classification associated with PCA band?
  if (n_total >= 10) {
    tab <- table(
      classification[theta_class != "AMBIGUOUS"]$theta_class,
      classification[theta_class != "AMBIGUOUS"]$band
    )
    fisher_p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
  } else {
    fisher_p <- NA_real_
  }

  data.table(
    n_classified = n_total,
    n_concordant = n_concordant,
    concordance_rate = round(concordance_rate, 3),
    fisher_p = fisher_p,
    significant = !is.na(fisher_p) & fisher_p < 0.001
  )
}
