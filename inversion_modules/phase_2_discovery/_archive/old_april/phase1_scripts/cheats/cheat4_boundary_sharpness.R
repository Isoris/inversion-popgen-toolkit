#!/usr/bin/env Rscript
# =============================================================================
# cheat4_boundary_sharpness.R — Fst/dXY/Hobs step-function at boundaries
#
# For each candidate boundary position, computes population genetic statistics
# between LEFT and RIGHT sides using Engine B (region_stats_dispatcher).
#
# Sharp step = real structural boundary.
# Diffuse transition = family LD or noise.
#
# REQUIRES: load_bridge.R sourced (provides get_region_stats, smap, reg)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# Guard: need dispatcher
if (!exists("get_region_stats", mode = "function")) {
  message("[cheat4] get_region_stats not available — skipping boundary sharpness")
}

HALF_FLANK <- 20L  # windows each side

#' Compute boundary sharpness for a single boundary position
#' @param chr Chromosome name
#' @param boundary_bp Boundary position in bp
#' @param dt Precomp data.table (with start_bp, end_bp, PC_1_* columns)
#' @param sample_names CGA sample names
#' @return data.table row with sharpness metrics
compute_boundary_sharpness <- function(chr, boundary_bp, dt, sample_names) {
  # Find the window index closest to boundary
  mid_bps <- (dt$start_bp + dt$end_bp) / 2
  bnd_idx <- which.min(abs(mid_bps - boundary_bp))

  left_idx  <- max(1, bnd_idx - HALF_FLANK):(bnd_idx - 1)
  right_idx <- (bnd_idx + 1):min(nrow(dt), bnd_idx + HALF_FLANK)
  if (length(left_idx) < 5 || length(right_idx) < 5) {
    return(data.table(boundary_bp = boundary_bp, sharpness = NA_real_))
  }

  left_start  <- dt$start_bp[min(left_idx)]
  left_end    <- dt$end_bp[max(left_idx)]
  right_start <- dt$start_bp[min(right_idx)]
  right_end   <- dt$end_bp[max(right_idx)]

  # Get k=3 band assignments from PC1 in the broader region
  all_idx <- c(left_idx, right_idx)
  pc1_cols <- paste0("PC_1_", sample_names)
  pc1_cols <- intersect(pc1_cols, names(dt))
  if (length(pc1_cols) < 20) return(data.table(boundary_bp = boundary_bp, sharpness = NA_real_))

  pc1_mat <- as.matrix(dt[all_idx, ..pc1_cols])
  avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
  valid <- is.finite(avg_pc1)
  if (sum(valid) < 20) return(data.table(boundary_bp = boundary_bp, sharpness = NA_real_))

  km <- tryCatch(kmeans(avg_pc1[valid], centers = 3, nstart = 10), error = function(e) NULL)
  if (is.null(km)) return(data.table(boundary_bp = boundary_bp, sharpness = NA_real_))

  co <- order(km$centers[, 1])
  bands <- integer(sum(valid))
  for (b in 1:3) bands[km$cluster == co[b]] <- b
  names(bands) <- sub("^PC_1_", "", names(avg_pc1)[valid])

  band1_samples <- names(bands)[bands == 1]
  band3_samples <- names(bands)[bands == 3]
  het_samples   <- names(bands)[bands == 2]

  if (length(band1_samples) < 5 || length(band3_samples) < 5) {
    return(data.table(boundary_bp = boundary_bp, sharpness = NA_real_))
  }

  # Use Engine B for Fst between band1 vs band3
  fst_left <- fst_right <- hobs_left <- hobs_right <- NA_real_

  tryCatch({
    groups_arg <- list(band1 = band1_samples, band3 = band3_samples)

    stats_left <- get_region_stats(chr, left_start, left_end,
                                    what = c("Fst", "Hobs"), groups = groups_arg)
    stats_right <- get_region_stats(chr, right_start, right_end,
                                     what = c("Fst", "Hobs"), groups = groups_arg)

    fst_left  <- stats_left$Fst$band1_vs_band3
    fst_right <- stats_right$Fst$band1_vs_band3

    # Hobs for HET samples specifically
    if (length(het_samples) >= 5) {
      het_groups <- list(het = het_samples)
      hobs_l <- get_region_stats(chr, left_start, left_end,
                                  what = "Hobs", groups = het_groups)
      hobs_r <- get_region_stats(chr, right_start, right_end,
                                  what = "Hobs", groups = het_groups)
      hobs_left  <- hobs_l$Hobs$het
      hobs_right <- hobs_r$Hobs$het
    }
  }, error = function(e) {
    message("[cheat4] Engine B error: ", e$message)
  })

  # Compute step metrics
  fst_step <- if (is.finite(fst_left) && is.finite(fst_right)) {
    abs(fst_left - fst_right)
  } else NA_real_

  hobs_contrast <- if (is.finite(hobs_left) && is.finite(hobs_right)) {
    abs(hobs_left - hobs_right)
  } else NA_real_

  # Combined sharpness: weighted
  sharpness <- 0
  n_contrib <- 0
  if (is.finite(fst_step)) {
    sharpness <- sharpness + pmin(1, fst_step / 0.15) * 0.5
    n_contrib <- n_contrib + 1
  }
  if (is.finite(hobs_contrast)) {
    sharpness <- sharpness + pmin(1, hobs_contrast / 0.05) * 0.5
    n_contrib <- n_contrib + 1
  }
  if (n_contrib > 0) sharpness <- sharpness / (n_contrib * 0.5)

  data.table(
    boundary_bp = boundary_bp,
    fst_left = round(fst_left, 4), fst_right = round(fst_right, 4),
    fst_step = round(fst_step, 4),
    hobs_left = round(hobs_left, 4), hobs_right = round(hobs_right, 4),
    hobs_contrast = round(hobs_contrast, 4),
    sharpness = round(sharpness, 4),
    sharp_class = if (is.finite(sharpness)) {
      if (sharpness >= 0.7) "sharp" else if (sharpness >= 0.3) "moderate" else "diffuse"
    } else "unknown"
  )
}

#' Scan all boundaries for a candidate region
#' @param chr Chromosome
#' @param start_bp, end_bp Region bounds
#' @param dt Precomp data.table
#' @param sample_names CGA names
#' @return data.table with sharpness per boundary (left + right of region)
scan_candidate_boundaries <- function(chr, start_bp, end_bp, dt, sample_names) {
  boundaries <- c(start_bp, end_bp)
  rbindlist(lapply(boundaries, function(bp) {
    compute_boundary_sharpness(chr, bp, dt, sample_names)
  }), fill = TRUE)
}
