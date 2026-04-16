#!/usr/bin/env Rscript
# =============================================================================
# cheat10_depth_boundary.R — Read depth anomaly at sim_mat boundaries
#
# At real inversion breakpoints, read depth often shows:
#   - DIPS: reads that span the breakpoint map poorly → local coverage drop
#   - SPIKES: duplicated/rearranged segments → coverage excess
#   - VARIANCE: HET samples have mixed-mapping reads → higher per-sample CV
#
# Multi-scale analysis at each boundary:
#   TIGHT:  ±80 windows  (~200 kb at 2.5 kb/window) — breakpoint-proximal
#   MEDIUM: ±160 windows (~400 kb) — captures broader rearrangement effects
#   WIDE:   ±320 windows (~800 kb) — context for normalisation
#
# Boundary types (from merge/landscape):
#   HARD: sharp sim_mat drop → expect depth anomaly → structural breakpoint
#   SOFT: gradual transition → expect NO anomaly → recombination zone
#   If hard boundary has NO depth anomaly → maybe family LD, not structural
#   If soft boundary HAS depth anomaly → hidden breakpoint, investigate
#
# DATA SOURCE: mosdepth per-sample BED files (already computed for PAV pipeline)
#   Path: delly_sv/00_markdup/{sample}.markdup.mosdepth.regions.bed.gz
#   OR:   MODULE_4D_PAV/mosdepth_output/{sample}.regions.bed.gz
#   Format: chrom, start, end, depth (per ~500bp bin)
#
# REQUIRES: load_bridge.R sourced (get_sample_ids)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# Scales: number of windows (not bp — actual bp depends on window size)
SCALES <- list(
  tight  = 80L,
  medium = 160L,
  wide   = 320L
)

#' Find mosdepth BED file for a sample
#' @param sample_id CGA sample name
#' @param base_dir BASE path
#' @return Path to mosdepth regions BED, or NULL
find_mosdepth_bed <- function(sample_id, base_dir = NULL) {
  if (is.null(base_dir)) base_dir <- Sys.getenv("BASE",
    "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04")

  candidates <- c(
    file.path(base_dir, "delly_sv/00_markdup",
              paste0(sample_id, ".markdup.mosdepth.regions.bed.gz")),
    file.path(base_dir, "MODULE_4D_PAV/mosdepth_output",
              paste0(sample_id, ".regions.bed.gz")),
    file.path(base_dir, "delly_sv/00_markdup",
              paste0(sample_id, ".mosdepth.regions.bed.gz"))
  )
  for (p in candidates) if (file.exists(p)) return(p)
  NULL
}

#' Load mosdepth depth for a region across all samples (or subset)
#' @param chr Chromosome
#' @param start_bp, end_bp Region bounds
#' @param sample_ids CGA sample names to query
#' @param base_dir BASE path
#' @return data.table: sample_id, bin_start, bin_end, depth
load_depth_region <- function(chr, start_bp, end_bp, sample_ids, base_dir = NULL) {
  rows <- list()
  for (sid in sample_ids) {
    bed <- find_mosdepth_bed(sid, base_dir)
    if (is.null(bed)) next

    # Use tabix if indexed, otherwise grep
    cmd <- if (file.exists(paste0(bed, ".tbi"))) {
      sprintf("tabix %s %s:%d-%d 2>/dev/null", bed, chr, start_bp, end_bp)
    } else {
      sprintf("zcat %s | awk '$1==\"%s\" && $2>=%d && $3<=%d'",
              bed, chr, start_bp, end_bp)
    }

    dt <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t"),
                    error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next
    setnames(dt, c("chrom", "bin_start", "bin_end", "depth"))
    dt[, sample_id := sid]
    rows[[length(rows) + 1]] <- dt
  }
  if (length(rows) > 0) rbindlist(rows) else data.table()
}

#' Compute depth metrics for one side of a boundary
#' @param depth_dt data.table from load_depth_region
#' @param sample_ids Samples to analyze
#' @return list(mean_depth, cv_depth, per_sample_means)
summarize_depth_side <- function(depth_dt, sample_ids) {
  if (nrow(depth_dt) == 0) return(list(mean_depth = NA_real_, cv_depth = NA_real_,
                                         per_sample_means = setNames(rep(NA_real_, length(sample_ids)), sample_ids)))

  per_sample <- depth_dt[, .(mean_depth = mean(depth, na.rm = TRUE)), by = sample_id]
  mean_all <- mean(per_sample$mean_depth, na.rm = TRUE)
  sd_all <- sd(per_sample$mean_depth, na.rm = TRUE)
  cv <- if (mean_all > 0) sd_all / mean_all else NA_real_

  psm <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
  psm[per_sample$sample_id] <- per_sample$mean_depth

  list(mean_depth = mean_all, cv_depth = cv, per_sample_means = psm)
}

#' Compute depth boundary score at one position, multiple scales
#' @param chr Chromosome
#' @param boundary_bp Boundary position
#' @param dt Precomp data.table (for window grid)
#' @param sample_ids CGA names
#' @param pc1_bands Named integer vector (sample → band) for per-group analysis
#' @param boundary_type "hard" or "soft"
#' @return data.table row per scale
compute_depth_boundary <- function(chr, boundary_bp, dt, sample_ids,
                                    pc1_bands = NULL, boundary_type = "hard") {
  mid_bps <- (dt$start_bp + dt$end_bp) / 2
  bnd_idx <- which.min(abs(mid_bps - boundary_bp))

  rows <- list()
  for (scale_name in names(SCALES)) {
    half <- SCALES[[scale_name]]
    left_start_idx  <- max(1, bnd_idx - half)
    left_end_idx    <- max(1, bnd_idx - 1)
    right_start_idx <- min(nrow(dt), bnd_idx + 1)
    right_end_idx   <- min(nrow(dt), bnd_idx + half)

    if (left_end_idx <= left_start_idx || right_end_idx <= right_start_idx) {
      rows[[length(rows) + 1]] <- data.table(
        boundary_bp = boundary_bp, scale = scale_name,
        boundary_type = boundary_type,
        depth_left = NA_real_, depth_right = NA_real_,
        depth_ratio = NA_real_, cv_left = NA_real_, cv_right = NA_real_,
        depth_dip = NA_real_, depth_score = NA_real_)
      next
    }

    left_start_bp  <- dt$start_bp[left_start_idx]
    left_end_bp    <- dt$end_bp[left_end_idx]
    right_start_bp <- dt$start_bp[right_start_idx]
    right_end_bp   <- dt$end_bp[right_end_idx]

    # Sample subset: use ≤30 samples for speed (stratified by band if available)
    sub_samples <- sample_ids
    if (length(sample_ids) > 30) {
      if (!is.null(pc1_bands)) {
        sub_samples <- c(
          sample(names(pc1_bands)[pc1_bands == 1], min(10, sum(pc1_bands == 1))),
          sample(names(pc1_bands)[pc1_bands == 2], min(10, sum(pc1_bands == 2))),
          sample(names(pc1_bands)[pc1_bands == 3], min(10, sum(pc1_bands == 3)))
        )
      } else {
        sub_samples <- sample(sample_ids, 30)
      }
    }

    depth_left  <- load_depth_region(chr, left_start_bp, left_end_bp, sub_samples)
    depth_right <- load_depth_region(chr, right_start_bp, right_end_bp, sub_samples)

    sl <- summarize_depth_side(depth_left, sub_samples)
    sr <- summarize_depth_side(depth_right, sub_samples)

    # Depth ratio: sharp change = structural breakpoint
    depth_ratio <- if (is.finite(sl$mean_depth) && is.finite(sr$mean_depth) &&
                        max(sl$mean_depth, sr$mean_depth) > 0) {
      min(sl$mean_depth, sr$mean_depth) / max(sl$mean_depth, sr$mean_depth)
    } else NA_real_

    # Depth dip: at the boundary itself (narrow ±5 kb)
    dip_depth <- load_depth_region(chr, boundary_bp - 5000, boundary_bp + 5000, sub_samples)
    dip_mean <- if (nrow(dip_depth) > 0) mean(dip_depth$depth, na.rm = TRUE) else NA_real_
    context_mean <- mean(c(sl$mean_depth, sr$mean_depth), na.rm = TRUE)
    depth_dip <- if (is.finite(dip_mean) && is.finite(context_mean) && context_mean > 0) {
      1 - dip_mean / context_mean  # positive = dip, negative = spike
    } else NA_real_

    # Per-band CV comparison (HET samples should have higher CV at breakpoints)
    cv_het <- cv_hom <- NA_real_
    if (!is.null(pc1_bands)) {
      het_ids <- intersect(names(pc1_bands)[pc1_bands == 2], sub_samples)
      hom_ids <- intersect(names(pc1_bands)[pc1_bands %in% c(1, 3)], sub_samples)
      if (length(het_ids) >= 5) {
        het_depths <- c(
          sl$per_sample_means[het_ids],
          sr$per_sample_means[het_ids]
        )
        het_depths <- het_depths[is.finite(het_depths)]
        cv_het <- if (length(het_depths) >= 5 && mean(het_depths) > 0)
          sd(het_depths) / mean(het_depths) else NA_real_
      }
      if (length(hom_ids) >= 5) {
        hom_depths <- c(
          sl$per_sample_means[hom_ids],
          sr$per_sample_means[hom_ids]
        )
        hom_depths <- hom_depths[is.finite(hom_depths)]
        cv_hom <- if (length(hom_depths) >= 5 && mean(hom_depths) > 0)
          sd(hom_depths) / mean(hom_depths) else NA_real_
      }
    }

    # Composite depth score
    score <- 0; n_c <- 0
    if (is.finite(depth_dip) && abs(depth_dip) > 0.05) {
      score <- score + pmin(1, abs(depth_dip) / 0.3) * 0.4; n_c <- n_c + 1
    }
    if (is.finite(depth_ratio) && depth_ratio < 0.9) {
      score <- score + (1 - depth_ratio) * 0.3; n_c <- n_c + 1
    }
    if (is.finite(cv_het) && is.finite(cv_hom) && cv_het > cv_hom * 1.2) {
      score <- score + pmin(1, (cv_het - cv_hom) / 0.1) * 0.3; n_c <- n_c + 1
    }
    if (n_c > 0) score <- score / (n_c * max(0.3, 0.4))
    score <- pmin(1, pmax(0, score))

    rows[[length(rows) + 1]] <- data.table(
      boundary_bp = boundary_bp, scale = scale_name,
      boundary_type = boundary_type,
      depth_left = round(sl$mean_depth, 2),
      depth_right = round(sr$mean_depth, 2),
      depth_ratio = round(depth_ratio, 4),
      cv_left = round(sl$cv_depth, 4),
      cv_right = round(sr$cv_depth, 4),
      cv_het = round(cv_het, 4),
      cv_hom = round(cv_hom, 4),
      depth_dip = round(depth_dip, 4),
      depth_score = round(score, 4)
    )
  }
  rbindlist(rows, fill = TRUE)
}

#' Scan all boundaries of a candidate at all scales
scan_depth_boundaries <- function(chr, boundaries_dt, precomp_dt, sample_ids,
                                   pc1_bands = NULL) {
  # boundaries_dt must have: boundary_bp, boundary_type
  if (nrow(boundaries_dt) == 0) return(data.table())

  results <- lapply(seq_len(nrow(boundaries_dt)), function(bi) {
    compute_depth_boundary(
      chr, boundaries_dt$boundary_bp[bi], precomp_dt, sample_ids,
      pc1_bands, boundaries_dt$boundary_type[bi]
    )
  })
  rbindlist(results, fill = TRUE)
}
