#!/usr/bin/env Rscript
# =============================================================================
# cheat13_breakpoint_proximity_weight.R — Breakpoint-proximal weighting
#
# BIOLOGY:
#   Windows near inversion breakpoints retain cleaner historical signal
#   than central windows. Near breakpoints, variation comes mainly from
#   new mutations (double crossovers negligible, only gene conversion).
#   Central windows accumulate exchange variation from both arrangements.
#   Breakpoint-proximal windows therefore carry MORE weight in scoring.
#
# INPUT:  boundary positions, precomp dt (per-window metrics)
# OUTPUT: weight vector, weighted_fst, breakpoint_signal_ratio
# REQUIRES: precomp RDS
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
BP_WEIGHT_FLOOR    <- 0.05
BP_DECAY_DIVISOR   <- 10     # decay_kb = span_kb / this

# ── Core: compute per-window proximity weights ─────────────────────────

compute_proximity_weights <- function(candidate_start, candidate_end,
                                       left_boundary_bp, right_boundary_bp,
                                       window_positions, decay_kb = NULL) {
  if (length(window_positions) == 0) return(numeric(0))
  if (is.na(left_boundary_bp) || is.na(right_boundary_bp))
    return(rep(1, length(window_positions)))

  span_bp <- right_boundary_bp - left_boundary_bp
  if (span_bp <= 0) return(rep(1, length(window_positions)))

  if (is.null(decay_kb)) {
    decay_kb <- max(20, (span_bp / 1000) / BP_DECAY_DIVISOR)
  }
  decay_bp <- decay_kb * 1000

  dist_nearest <- pmin(abs(window_positions - left_boundary_bp),
                        abs(window_positions - right_boundary_bp))
  weights <- 1 / (1 + dist_nearest / decay_bp)
  pmax(weights, BP_WEIGHT_FLOOR)
}

# ── Apply weights to a score vector ────────────────────────────────────

apply_proximity_weights <- function(score_vector, weights) {
  if (length(score_vector) == 0 || length(weights) == 0) return(NA_real_)
  ok <- is.finite(score_vector) & is.finite(weights)
  if (sum(ok) == 0) return(NA_real_)
  sum(score_vector[ok] * weights[ok]) / sum(weights[ok])
}

# ── Weighted Fst at candidate ──────────────────────────────────────────

weighted_fst_at_candidate <- function(dt, sample_names = NULL,
                                       band_assignments = NULL,
                                       left_bp, right_bp,
                                       decay_kb = NULL,
                                       fst_col = "fst_band13") {
  if (nrow(dt) == 0)
    return(list(weighted_fst = NA_real_, unweighted_fst = NA_real_,
                breakpoint_signal_ratio = NA_real_, weights = numeric(0)))

  win_pos <- (dt$start_bp + dt$end_bp) / 2
  wts <- compute_proximity_weights(min(dt$start_bp), max(dt$end_bp),
                                    left_bp, right_bp, win_pos, decay_kb)
  fst <- if (fst_col %in% names(dt)) dt[[fst_col]]
         else if ("inv_likeness" %in% names(dt)) dt$inv_likeness
         else rep(NA_real_, nrow(dt))

  w_fst <- apply_proximity_weights(fst, wts)
  u_fst <- mean(fst, na.rm = TRUE)
  bsr   <- if (is.finite(w_fst) && is.finite(u_fst) && u_fst > 0)
             round(w_fst / u_fst, 4) else NA_real_

  list(weighted_fst = w_fst, unweighted_fst = u_fst,
       breakpoint_signal_ratio = bsr, weights = wts)
}

# ── Search mode ────────────────────────────────────────────────────────

search_proximity_weight <- function(chr, zone_start, zone_end, dt = NULL,
                                     candidate_start = NULL,
                                     candidate_end = NULL) {
  empty <- data.table(method = "proximity_weight", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(dt) || nrow(dt) == 0) return(empty)
  if (is.null(candidate_start)) candidate_start <- min(dt$start_bp)
  if (is.null(candidate_end))   candidate_end   <- max(dt$end_bp)

  zone_dt <- dt[start_bp <= zone_end & end_bp >= zone_start]
  if (nrow(zone_dt) < 3) return(empty)
  sig_col <- if ("inv_likeness" %in% names(zone_dt)) "inv_likeness" else return(empty)
  zone_pos <- (zone_dt$start_bp + zone_dt$end_bp) / 2
  zone_sig <- zone_dt[[sig_col]]

  other_bp <- if (abs(zone_start - candidate_start) < abs(zone_end - candidate_end))
    candidate_end else candidate_start
  test_bps <- seq(zone_start, zone_end, length.out = min(20, nrow(zone_dt)))
  best_sc <- 0; best_bp <- NA_integer_
  for (tbp in test_bps) {
    wts <- compute_proximity_weights(candidate_start, candidate_end,
                                      tbp, other_bp, zone_pos)
    ratio <- apply_proximity_weights(zone_sig, wts) /
             mean(zone_sig, na.rm = TRUE)
    sc <- abs(ratio - 1)
    if (is.finite(sc) && sc > best_sc) { best_sc <- sc; best_bp <- as.integer(round(tbp)) }
  }
  data.table(method = "proximity_weight", best_bp = best_bp,
             score = round(pmin(1, best_sc), 3), is_precise = FALSE,
             detail = paste0("bsr_dev=", round(best_sc, 4)))
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat13 <- function(chr, start_bp, end_bp, left_bp, right_bp,
                         dt, sample_names = NULL, band_assignments = NULL,
                         decay_kb = NULL) {
  message("[cheat13] ", chr, ":", round(start_bp/1e6,1), "-",
          round(end_bp/1e6,1), " Mb | span ",
          round((right_bp - left_bp)/1e3,0), " kb")

  region_dt <- dt[start_bp >= (..start_bp) & end_bp <= (..end_bp)]
  if (nrow(region_dt) == 0) {
    message("[cheat13] No windows in region")
    return(list(weights = numeric(0),
                weighted_fst_result = list(weighted_fst = NA_real_,
                  unweighted_fst = NA_real_, breakpoint_signal_ratio = NA_real_),
                search_result = data.table(method = "proximity_weight",
                  best_bp = NA_integer_, score = 0, is_precise = FALSE,
                  detail = "no_windows")))
  }

  wfr <- weighted_fst_at_candidate(region_dt, sample_names, band_assignments,
                                    left_bp, right_bp, decay_kb)
  bsr <- wfr$breakpoint_signal_ratio
  interp <- if (!is.finite(bsr)) "NA"
    else if (bsr > 1.3) "bp-sharp (old inversion)"
    else if (bsr > 0.9) "uniform (young/short)"
    else "center-sharp (unusual)"
  message("[cheat13] wFst=", round(wfr$weighted_fst,4),
          " uFst=", round(wfr$unweighted_fst,4),
          " BSR=", bsr, " → ", interp)

  list(weights = wfr$weights, weighted_fst_result = wfr,
       search_result = search_proximity_weight(chr, start_bp, end_bp,
                                                region_dt, start_bp, end_bp))
}
