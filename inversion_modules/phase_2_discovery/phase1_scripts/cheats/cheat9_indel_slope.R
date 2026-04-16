#!/usr/bin/env Rscript
# =============================================================================
# cheat9_indel_slope.R — Cumulative indel slope heterokaryotype detector
#
# HET samples carry two divergent haplotypes → MORE het indels per kb.
# HOM_INV samples have arrangement-specific indels fixed → MORE hom indels.
#
# het_slope/hom_slope ratio separates three genotype classes independently
# from PCA, SV callers, Fst, and phasing. Pure counting.
#
# Data source: Clair3 postprocess `all_variants_with_phase.tsv` per sample.
#   Columns: IS_INDEL (col 14), GT_CLASS (col 20), ABS_INDEL_LEN (col 18)
#
# REQUIRES: load_bridge.R sourced (get_sample_ids), Clair3 output available
# =============================================================================

suppressPackageStartupMessages(library(data.table))

WINDOW_KB <- 50L  # sliding window for slope computation

#' Compute cumulative indel slopes for one sample in one region
#' @param phase_file Path to all_variants_with_phase.tsv
#' @param start_bp, end_bp Region bounds
#' @return list(het_slope, hom_slope, n_het_indels, n_hom_indels)
compute_indel_slopes <- function(phase_file, start_bp, end_bp) {
  result <- list(het_slope = NA_real_, hom_slope = NA_real_,
                 n_het_indels = 0L, n_hom_indels = 0L,
                 het_slope_ratio = NA_real_)

  if (!file.exists(phase_file)) return(result)

  # Read relevant columns: POS (col 2), IS_INDEL (col 14), GT_CLASS (col 20)
  dt <- tryCatch({
    fread(phase_file, select = c(2, 14, 20), header = FALSE,
          colClasses = list(integer = c(1, 2), character = 3))
  }, error = function(e) NULL)

  if (is.null(dt) || nrow(dt) == 0) return(result)
  setnames(dt, c("POS", "IS_INDEL", "GT_CLASS"))

  # Filter to region and indels only
  dt <- dt[POS >= start_bp & POS <= end_bp & IS_INDEL == 1]
  if (nrow(dt) < 5) return(result)

  # Split by genotype class
  het_indels <- dt[GT_CLASS == "het"]
  hom_indels <- dt[GT_CLASS %in% c("hom_alt", "hom_var")]

  result$n_het_indels <- nrow(het_indels)
  result$n_hom_indels <- nrow(hom_indels)

  region_kb <- (end_bp - start_bp) / 1000

  # Slopes: indels per kb
  result$het_slope <- if (nrow(het_indels) >= 2) {
    nrow(het_indels) / region_kb
  } else 0

  result$hom_slope <- if (nrow(hom_indels) >= 2) {
    nrow(hom_indels) / region_kb
  } else 0

  # Ratio: high → likely HET carrier (many het indels relative to hom)
  if (result$hom_slope > 0) {
    result$het_slope_ratio <- result$het_slope / result$hom_slope
  } else if (result$het_slope > 0) {
    result$het_slope_ratio <- Inf  # no hom indels at all → strong HET signal
  }

  result
}

#' Scan all samples for indel slopes in a candidate region
#' @param chr Chromosome
#' @param start_bp, end_bp Region bounds
#' @param sample_names CGA sample names
#' @param clair3_dir Path to Clair3 postprocess_results
#' @return data.table with per-sample slopes + predicted genotype
scan_indel_slopes <- function(chr, start_bp, end_bp, sample_names, clair3_dir) {
  rows <- list()

  for (sid in sample_names) {
    phase_file <- file.path(clair3_dir, chr, sid, "all_variants_with_phase.tsv")
    slopes <- compute_indel_slopes(phase_file, start_bp, end_bp)

    rows[[length(rows) + 1]] <- data.table(
      sample_id = sid,
      het_slope = round(slopes$het_slope, 4),
      hom_slope = round(slopes$hom_slope, 4),
      het_slope_ratio = round(slopes$het_slope_ratio, 4),
      n_het_indels = slopes$n_het_indels,
      n_hom_indels = slopes$n_hom_indels
    )
  }

  out <- rbindlist(rows)

  if (nrow(out) < 20) return(out)

  # Classify using bimodal detection on het_slope
  valid <- out[is.finite(het_slope) & het_slope > 0]
  if (nrow(valid) < 20) {
    out[, indel_genotype := "UNKNOWN"]
    return(out)
  }

  # k-means k=2 on het_slope → bimodal split
  km2 <- tryCatch(kmeans(valid$het_slope, centers = 2, nstart = 10),
                   error = function(e) NULL)

  if (!is.null(km2)) {
    co <- order(km2$centers[, 1])
    # Lower het_slope cluster = HOM carriers (similar haplotypes)
    # Higher het_slope cluster = HET carriers (divergent haplotypes)
    valid[, indel_class := ifelse(km2$cluster == co[1], "HOM_like", "HET_like")]

    # Try k=3 to separate HOM_REF from HOM_INV using hom_slope
    hom_like <- valid[indel_class == "HOM_like" & is.finite(hom_slope)]
    if (nrow(hom_like) >= 10) {
      km_hom <- tryCatch(kmeans(hom_like$hom_slope, centers = 2, nstart = 10),
                          error = function(e) NULL)
      if (!is.null(km_hom)) {
        co_h <- order(km_hom$centers[, 1])
        hom_like[, indel_subclass := ifelse(km_hom$cluster == co_h[1],
                                             "HOM_REF_like", "HOM_INV_like")]
        valid[hom_like, on = "sample_id", indel_class := i.indel_subclass]
      }
    }

    out[valid, on = "sample_id", indel_genotype := i.indel_class]
    out[is.na(indel_genotype), indel_genotype := "UNKNOWN"]

    # Report
    tab <- table(out$indel_genotype)
    message("[cheat9] Indel slopes (", chr, ":", start_bp, "-", end_bp, "): ",
            paste(names(tab), tab, sep = "=", collapse = " "))
  } else {
    out[, indel_genotype := "UNKNOWN"]
  }

  out
}
