#!/usr/bin/env Rscript
# =============================================================================
# cheat17_fossil_breakpoints.R — Fossil breakpoint detection
#
# BIOLOGY:
#   Some boundaries in the sim_mat look structural but have no active
#   inversion behind them. These may be relics of old inversions that
#   broke apart through subsequent rearrangement. A "fossil breakpoint"
#   shows the sim_mat color transition but no trimodal PCA, no SV calls,
#   and no carrier/non-carrier distinction.
#
# INPUT:  boundary catalog, candidate regions, precomp RDS
# OUTPUT: fossil classification for orphaned boundaries
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
INV_LIKENESS_THRESHOLD <- 0.3   # elevated if above this
DIP_P_THRESHOLD        <- 0.05  # significant trimodality
MIN_BAND_SIZE          <- 10    # minimum samples per band for real inversion
FST_STEP_THRESHOLD     <- 0.05  # meaningful Fst step

# ── Identify orphaned boundaries ──────────────────────────────────────

identify_orphaned_boundaries <- function(boundary_catalog, candidate_regions) {
  if (nrow(boundary_catalog) == 0) return(data.table())
  if (nrow(candidate_regions) == 0) return(boundary_catalog)

  # A boundary is orphaned if it is NOT inside any candidate region
  boundary_catalog[, orphaned := TRUE]
  for (i in seq_len(nrow(candidate_regions))) {
    cr <- candidate_regions[i]
    boundary_catalog[
      chr == cr$chr &
      boundary_bp >= cr$start_bp &
      boundary_bp <= cr$end_bp,
      orphaned := FALSE
    ]
  }
  # Filter to structural types only
  structural_types <- c("STRUCTURAL_SHARP", "STRUCTURAL_MODERATE",
                         "STRUCTURAL_DIFFUSE")
  if ("boundary_type" %in% names(boundary_catalog)) {
    boundary_catalog[orphaned == TRUE &
                      boundary_type %in% structural_types]
  } else {
    boundary_catalog[orphaned == TRUE]
  }
}

# ── Test fossil signal at an orphaned boundary ────────────────────────

test_fossil_signal <- function(orphaned_bp, chr, precomp_dt,
                                sample_names = NULL, fl = NULL) {
  if (nrow(precomp_dt) == 0)
    return(list(classification = "NO_DATA", inv_likeness_elevated = NA,
                trimodal = NA, sv_overlap = NA, fst_step = NA))

  # Window nearest to boundary
  bp_pos <- orphaned_bp
  precomp_dt[, dist_to_bp := abs((start_bp + end_bp)/2 - bp_pos)]
  nearest <- precomp_dt[order(dist_to_bp)][1:min(5, nrow(precomp_dt))]

  # 1. inv_likeness elevated?
  il <- if ("inv_likeness" %in% names(nearest))
    mean(nearest$inv_likeness, na.rm = TRUE) else NA_real_
  il_elevated <- !is.na(il) && il > INV_LIKENESS_THRESHOLD

  # 2. Trimodality (dip test if diptest available)
  trimodal <- FALSE
  if ("PC1_values" %in% names(attributes(nearest)) ||
      "pve1" %in% names(nearest)) {
    # Use PVE1 as proxy: high PVE1 suggests structured
    pve <- if ("pve1" %in% names(nearest)) mean(nearest$pve1, na.rm = TRUE)
           else NA_real_
    trimodal <- !is.na(pve) && pve > 0.15
  }

  # 3. SV overlap
  sv_overlap <- FALSE
  if (!is.null(fl)) {
    sv_call <- tryCatch(flashlight_summary(chr, bp_pos - 50000, bp_pos + 50000),
                         error = function(e) list(n_inv_overlapping = 0))
    sv_overlap <- sv_call$n_inv_overlapping > 0
  }

  # 4. Fst step (left vs right of boundary)
  left_dt  <- precomp_dt[end_bp < bp_pos][order(-end_bp)][1:min(5, sum(precomp_dt$end_bp < bp_pos))]
  right_dt <- precomp_dt[start_bp > bp_pos][order(start_bp)][1:min(5, sum(precomp_dt$start_bp > bp_pos))]
  fst_step <- NA_real_
  if (nrow(left_dt) > 0 && nrow(right_dt) > 0 &&
      "inv_likeness" %in% names(left_dt)) {
    fst_step <- abs(mean(right_dt$inv_likeness, na.rm = TRUE) -
                     mean(left_dt$inv_likeness, na.rm = TRUE))
  }
  has_fst_step <- !is.na(fst_step) && fst_step > FST_STEP_THRESHOLD

  # Classification decision tree
  classification <- if (il_elevated && trimodal && sv_overlap) {
    "MISSED_INVERSION"
  } else if (il_elevated && !trimodal && !sv_overlap) {
    "FOSSIL_CANDIDATE"
  } else if (!il_elevated) {
    "FALSE_BOUNDARY"
  } else if (has_fst_step && !trimodal) {
    "POSSIBLE_FOSSIL"
  } else {
    "AMBIGUOUS"
  }

  precomp_dt[, dist_to_bp := NULL]

  list(classification = classification,
       inv_likeness_elevated = il_elevated,
       inv_likeness_value = round(il, 4),
       trimodal = trimodal,
       sv_overlap = sv_overlap,
       fst_step = round(fst_step, 4),
       fst_step_significant = has_fst_step)
}

# ── Search mode ────────────────────────────────────────────────────────

search_fossil_breakpoint <- function(chr, zone_start, zone_end,
                                      precomp_dt = NULL, ...) {
  empty <- data.table(method = "fossil_breakpoint", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(precomp_dt) || nrow(precomp_dt) == 0) return(empty)

  mid <- as.integer((zone_start + zone_end) / 2)
  result <- test_fossil_signal(mid, chr, precomp_dt)
  sc <- if (result$classification == "FOSSIL_CANDIDATE") 0.8
        else if (result$classification == "POSSIBLE_FOSSIL") 0.5
        else if (result$classification == "MISSED_INVERSION") 1.0
        else 0.1

  data.table(method = "fossil_breakpoint", best_bp = mid,
             score = round(sc, 3), is_precise = FALSE,
             detail = result$classification)
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat17 <- function(chr, boundary_catalog, candidate_regions,
                         precomp_dt, fl = NULL) {
  message("[cheat17] ", chr, ": scanning for fossil breakpoints")

  orphans <- identify_orphaned_boundaries(boundary_catalog, candidate_regions)
  message("[cheat17] Found ", nrow(orphans), " orphaned structural boundaries")

  if (nrow(orphans) == 0)
    return(list(fossil_boundaries = data.table(),
                search_result = data.table(method = "fossil_breakpoint",
                  best_bp = NA_integer_, score = 0, is_precise = FALSE,
                  detail = "no_orphans")))

  results <- list()
  for (i in seq_len(nrow(orphans))) {
    bp <- orphans$boundary_bp[i]
    fossil <- test_fossil_signal(bp, chr, precomp_dt, fl = fl)
    results[[i]] <- data.table(
      chr = chr, boundary_bp = bp,
      classification = fossil$classification,
      inv_likeness = fossil$inv_likeness_value,
      trimodal = fossil$trimodal,
      sv_overlap = fossil$sv_overlap,
      fst_step = fossil$fst_step)
    message("[cheat17]   bp=", bp, " → ", fossil$classification)
  }
  fossil_dt <- rbindlist(results)

  tab <- table(fossil_dt$classification)
  message("[cheat17] Summary: ", paste(names(tab), tab, sep = "=", collapse = ", "))

  list(fossil_boundaries = fossil_dt,
       search_result = search_fossil_breakpoint(chr, min(orphans$boundary_bp),
                        max(orphans$boundary_bp), precomp_dt))
}
