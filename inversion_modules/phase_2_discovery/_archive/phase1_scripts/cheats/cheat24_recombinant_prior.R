#!/usr/bin/env Rscript
# =============================================================================
# cheat24_recombinant_prior.R — Position-aware recombinant prior
#
# BIOLOGY:
#   Recombination inside inversions is not uniform. Near breakpoints:
#   double crossovers ~10^-4, only gene conversion. Center of large
#   inversions: rare double crossovers ~10^-3. This gradient should
#   serve as a Bayesian prior when interpreting recombinant signals.
#
# INPUT:  boundary positions, candidate regions, C01h recombinant flags
# OUTPUT: per-position prior, event classification, recombinant map
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
PRIOR_GENE_CONVERSION <- 0.001   # constant at all positions
PRIOR_MAX_DCO         <- 0.01    # maximum double crossover rate at center
MOSAIC_SHORT_BP       <- 50000L  # short mosaic tract (gene conversion)
MOSAIC_LONG_BP        <- 200000L # long mosaic tract (double crossover)

# ── Logistic function for position-dependent prior ────────────────────

logistic_prior <- function(x, midpoint, steepness) {
  1 / (1 + exp(-steepness * (x - midpoint)))
}

# ── Compute recombinant prior at a position ───────────────────────────

compute_recombinant_prior <- function(position_bp, left_bp, right_bp) {
  span <- right_bp - left_bp
  if (span <= 0) return(list(prior_dco = 0, prior_gc = PRIOR_GENE_CONVERSION,
                              prior_total = PRIOR_GENE_CONVERSION, rel_pos = NA))

  dist_to_nearest <- min(abs(position_bp - left_bp),
                          abs(position_bp - right_bp))
  rel_pos <- (position_bp - left_bp) / span

  # Double crossover prior: logistic increase toward center
  midpoint <- span / 4
  steepness <- 8 / span
  prior_dco <- PRIOR_MAX_DCO * logistic_prior(dist_to_nearest, midpoint, steepness)

  # Gene conversion: constant
  prior_gc <- PRIOR_GENE_CONVERSION

  list(prior_dco = round(prior_dco, 6),
       prior_gc = round(prior_gc, 6),
       prior_total = round(prior_dco + prior_gc, 6),
       rel_pos = round(rel_pos, 4),
       dist_to_nearest_bp = dist_to_nearest)
}

# ── Compute prior across a vector of positions ───────────────────────

compute_prior_vector <- function(positions_bp, left_bp, right_bp) {
  results <- list()
  for (i in seq_along(positions_bp)) {
    p <- compute_recombinant_prior(positions_bp[i], left_bp, right_bp)
    results[[i]] <- data.table(position_bp = positions_bp[i],
                                rel_pos = p$rel_pos,
                                prior_dco = p$prior_dco,
                                prior_gc = p$prior_gc,
                                prior_total = p$prior_total,
                                dist_nearest_bp = p$dist_to_nearest_bp)
  }
  rbindlist(results)
}

# ── Classify recombinant events ───────────────────────────────────────

classify_recombinant_event <- function(recomb_position, left_bp, right_bp,
                                        mosaic_length = NA_integer_) {
  prior <- compute_recombinant_prior(recomb_position, left_bp, right_bp)
  span <- right_bp - left_bp

  dist_frac <- prior$dist_to_nearest_bp / span
  near_bp <- dist_frac < 0.15  # within 15% of breakpoint
  at_center <- dist_frac > 0.35

  short_mosaic <- !is.na(mosaic_length) && mosaic_length < MOSAIC_SHORT_BP
  long_mosaic  <- !is.na(mosaic_length) && mosaic_length >= MOSAIC_LONG_BP

  # Classification matrix
  event_class <- if (near_bp && (short_mosaic || is.na(mosaic_length))) {
    "gene_conversion"
  } else if (at_center && long_mosaic) {
    "double_crossover"
  } else if (near_bp && long_mosaic) {
    "suspicious"  # long mosaic near breakpoint is unexpected
  } else if (at_center && short_mosaic) {
    "ambiguous"   # short mosaic at center could be either
  } else {
    "ambiguous"
  }

  # Posterior = P(real | observed) ∝ P(observed | real) × prior
  # Simplified: use prior as weight
  posterior <- prior$prior_total * (if (event_class == "suspicious") 0.1
    else if (event_class == "gene_conversion") 0.5
    else if (event_class == "double_crossover") 0.8
    else 0.3)

  list(event_class = event_class,
       posterior_probability = round(posterior, 6),
       position_prior = prior$prior_total,
       dist_fraction = round(dist_frac, 3),
       near_breakpoint = near_bp)
}

# ── Build recombinant map for a candidate ─────────────────────────────

recombinant_map <- function(recombinant_events, left_bp, right_bp) {
  if (is.null(recombinant_events) || nrow(recombinant_events) == 0)
    return(data.table())

  results <- list()
  for (i in seq_len(nrow(recombinant_events))) {
    pos <- recombinant_events$position_bp[i]
    mosaic_len <- if ("mosaic_length" %in% names(recombinant_events))
      recombinant_events$mosaic_length[i] else NA_integer_

    classif <- classify_recombinant_event(pos, left_bp, right_bp, mosaic_len)
    results[[i]] <- data.table(
      position_bp = pos,
      event_class = classif$event_class,
      posterior = classif$posterior_probability,
      prior = classif$position_prior,
      dist_fraction = classif$dist_fraction,
      near_bp = classif$near_breakpoint)
  }
  rbindlist(results)
}

# ── Search mode ────────────────────────────────────────────────────────

search_recombinant_prior <- function(chr, zone_start, zone_end,
                                      left_bp = NULL, right_bp = NULL,
                                      ...) {
  empty <- data.table(method = "recombinant_prior", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(left_bp) || is.null(right_bp)) return(empty)

  mid <- as.integer((zone_start + zone_end) / 2)
  prior <- compute_recombinant_prior(mid, left_bp, right_bp)
  sc <- round(prior$prior_total * 100, 3)  # scale for scoring

  data.table(method = "recombinant_prior", best_bp = mid,
             score = pmin(1, sc), is_precise = FALSE,
             detail = paste0("prior=", prior$prior_total,
                              ",rel_pos=", prior$rel_pos))
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat24 <- function(chr, candidate_start, candidate_end,
                         left_bp, right_bp,
                         recombinant_events = NULL) {
  span_kb <- (right_bp - left_bp) / 1000
  message("[cheat24] ", chr, ":", round(candidate_start/1e6,1), "-",
          round(candidate_end/1e6,1), " Mb | span ", round(span_kb,0), " kb")

  # Compute prior across candidate
  positions <- seq(left_bp, right_bp, length.out = min(100, span_kb))
  prior_dt <- compute_prior_vector(positions, left_bp, right_bp)
  message("[cheat24] Prior range: ", min(prior_dt$prior_total), " - ",
          max(prior_dt$prior_total))

  # Classify events if provided
  event_map <- if (!is.null(recombinant_events) && nrow(recombinant_events) > 0) {
    rmap <- recombinant_map(recombinant_events, left_bp, right_bp)
    if (nrow(rmap) > 0) {
      tab <- table(rmap$event_class)
      message("[cheat24] Events: ",
              paste(names(tab), tab, sep = "=", collapse = ", "))
    }
    rmap
  } else {
    message("[cheat24] No recombinant events provided")
    data.table()
  }

  list(prior_table = prior_dt, event_map = event_map,
       search_result = search_recombinant_prior(chr, candidate_start,
                        candidate_end, left_bp, right_bp))
}
