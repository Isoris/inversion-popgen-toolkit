# =============================================================================
# R/ushape_classify.R
# =============================================================================
# Pure scoring + classification functions for the U-shape evolution module.
# Sourceable from CLI scripts and from the popstats-server R subprocess.
#
# Classification model (the spec; deviates slightly from the original prompt
# for the reasons documented in the handoff conversation):
#
#   Step A — AGE GATE
#       compute oldness_score from depth + private SNP signal.
#       if oldness_score < threshold OR dxy_inside_mean < threshold:
#           label = young_weak_divergence
#           STOP — shape statistics are not informative on young inversions.
#       else continue.
#
#   Step B — SHAPE VERDICT (only on candidates that pass the gate)
#       The candidate enters one of:
#           neutral_like_U_shape
#           locally_adapted_internal_peak_like
#           locally_adapted_breakpoint_like
#           flat_high_deep_structural_haplotype
#           asymmetric_edge
#           complex_mixed
#       chosen by a small ordered rule set on the U/peak/asymmetry/flatness scores.
#
# The rule logic intentionally matches Berdan/Guerrero Fig 2: the U-shape only
# appears once T_in_Ne is large, so we gate on "old enough" first. This avoids
# the failure mode where a young candidate with noise in dxy_u_score gets
# called "neutral_like_U_shape".
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })


# -----------------------------------------------------------------------------
# tiny helpers
# -----------------------------------------------------------------------------
.safe_div <- function(num, den, eps = 1e-9) {
  if (is.null(num) || is.null(den)) return(NA_real_)
  if (!is.finite(num) || !is.finite(den)) return(NA_real_)
  if (abs(den) < eps) return(NA_real_)
  num / den
}

.safe_log2 <- function(x) {
  if (is.null(x) || !is.finite(x) || x <= 0) return(NA_real_)
  log2(x)
}


# -----------------------------------------------------------------------------
# shape scores
# -----------------------------------------------------------------------------
# Input: a one-row data.table (candidate_raw_summary) PLUS the per-window
# data.table for that candidate (so we can compute peak/flatness scores
# on the actual inside-windows).
score_candidate <- function(raw_row, window_dt, params) {
  d_in   <- raw_row$dxy_inside_mean
  d_fl   <- raw_row$dxy_flank_mean
  d_le   <- raw_row$dxy_left_edge_mean
  d_re   <- raw_row$dxy_right_edge_mean
  d_ce   <- raw_row$dxy_center_mean
  f_in   <- raw_row$fst_inside_mean
  f_fl   <- raw_row$fst_flank_mean
  f_le   <- raw_row$fst_left_edge_mean
  f_re   <- raw_row$fst_right_edge_mean
  f_ce   <- raw_row$fst_center_mean
  p1_in  <- raw_row$pi_homo1_inside_mean
  p2_in  <- raw_row$pi_homo2_inside_mean
  afd_in <- raw_row$allele_freq_delta_inside_mean

  # inside-only window subset for peak / flatness
  inside <- window_dt[zone %in% c("left_edge", "center", "right_edge")]
  edge   <- window_dt[zone %in% c("left_edge", "right_edge")]
  centr  <- window_dt[zone == "center"]

  # ratios
  dxy_inside_flank_ratio <- .safe_div(d_in, d_fl)
  dxy_u_score            <- .safe_div(mean(c(d_le, d_re), na.rm = TRUE), d_ce)
  dxy_asymmetry_score    <- .safe_div(d_le, d_re)
  abs_log2_asymmetry     <- abs(.safe_log2(dxy_asymmetry_score))
  dxy_internal_peak_score <- .safe_div(
    suppressWarnings(max(centr$dxy_homo1_homo2, na.rm = TRUE)),
    mean(edge$dxy_homo1_homo2, na.rm = TRUE)
  )
  dxy_flatness_score <- .safe_div(
    sd(inside$dxy_homo1_homo2, na.rm = TRUE),
    mean(inside$dxy_homo1_homo2, na.rm = TRUE)
  )
  fst_inside_flank_ratio <- .safe_div(f_in, f_fl)
  fst_u_score            <- .safe_div(mean(c(f_le, f_re), na.rm = TRUE), f_ce)
  fst_internal_peak_score <- .safe_div(
    suppressWarnings(max(centr$fst_homo1_homo2, na.rm = TRUE)),
    mean(edge$fst_homo1_homo2, na.rm = TRUE)
  )

  pi_ratio_homo1_homo2 <- .safe_div(p1_in, p2_in)
  pi_imbalance_score   <- abs(.safe_log2(pi_ratio_homo1_homo2))

  # composite scores in [0, ~3] range
  arrangement_contrast_score <- mean(c(
    pmin(3, dxy_inside_flank_ratio %||% NA),
    pmin(3, fst_inside_flank_ratio %||% NA),
    pmin(3, (afd_in %||% 0) * 5)
  ), na.rm = TRUE)

  # oldness: deep dXY relative to within-arrangement diversity AND high private
  # SNP load (HOMO_1 and HOMO_2 each carry their own variants since the split).
  pi_avg <- mean(c(p1_in, p2_in), na.rm = TRUE)
  oldness_score <- mean(c(
    pmin(3, dxy_inside_flank_ratio %||% NA),
    pmin(3, .safe_div(d_in, pi_avg) %||% NA),
    pmin(3, fst_inside_flank_ratio %||% NA),
    if (is.finite(raw_row$private_snp_homo1_count) &&
        is.finite(raw_row$private_snp_homo2_count) &&
        is.finite(raw_row$shared_snp_count) &&
        raw_row$shared_snp_count > 0)
      pmin(3, (raw_row$private_snp_homo1_count + raw_row$private_snp_homo2_count) /
              max(1L, raw_row$shared_snp_count)) else NA_real_
  ), na.rm = TRUE)

  # youngness is the inverse axis; we keep it as a continuous field for plotting
  youngness_score <- if (is.finite(oldness_score)) max(0, 2 - oldness_score) else NA_real_

  # shape composites
  neutral_u_shape_score <- mean(c(
    pmin(3, dxy_u_score %||% NA),
    pmin(3, fst_u_score %||% NA),
    -pmin(3, dxy_internal_peak_score %||% NA)
  ), na.rm = TRUE)

  local_adaptation_internal_score <- mean(c(
    pmin(3, dxy_internal_peak_score %||% NA),
    pmin(3, fst_internal_peak_score %||% NA)
  ), na.rm = TRUE)

  breakpoint_adaptation_score <- mean(c(
    pmin(3, .safe_div(mean(c(d_le, d_re), na.rm = TRUE), d_fl) %||% NA),
    pmin(3, .safe_div(mean(c(f_le, f_re), na.rm = TRUE), f_fl) %||% NA),
    pmin(3, dxy_u_score %||% NA)
  ), na.rm = TRUE)

  data.table(
    candidate_id = raw_row$candidate_id,
    log10_length_bp = log10(pmax(1L, raw_row$length_bp)),
    dxy_inside_flank_ratio = dxy_inside_flank_ratio,
    dxy_u_score = dxy_u_score,
    dxy_asymmetry_score = dxy_asymmetry_score,
    abs_log2_asymmetry = abs_log2_asymmetry,
    dxy_internal_peak_score = dxy_internal_peak_score,
    dxy_flatness_score = dxy_flatness_score,
    fst_inside_flank_ratio = fst_inside_flank_ratio,
    fst_u_score = fst_u_score,
    fst_internal_peak_score = fst_internal_peak_score,
    pi_ratio_homo1_homo2 = pi_ratio_homo1_homo2,
    pi_imbalance_score = pi_imbalance_score,
    arrangement_contrast_score = arrangement_contrast_score,
    oldness_score = oldness_score,
    youngness_score = youngness_score,
    neutral_u_shape_score = neutral_u_shape_score,
    local_adaptation_internal_score = local_adaptation_internal_score,
    breakpoint_adaptation_score = breakpoint_adaptation_score
  )
}


# -----------------------------------------------------------------------------
# classification (age-gated)
# -----------------------------------------------------------------------------
# Returns list(primary_class, secondary_class, confidence, reason, flags).
classify_candidate <- function(scores, raw_row, params, n_inside_windows) {
  thr <- list(
    u_high          = params$u_score_high           %||% 1.5,
    inflank_high    = params$inside_flank_high       %||% 1.5,
    peak_high       = params$internal_peak_high      %||% 1.5,
    asym_high       = params$asymmetry_log2_high     %||% 1.0,
    fst_high        = params$fst_enrichment_high     %||% 1.5,
    oldness_min     = params$oldness_min             %||% 1.2,
    dxy_min_inside  = params$dxy_min_inside          %||% 1e-4,
    min_windows     = params$min_inside_windows_for_shape %||% 8L,
    flatness_max    = params$flatness_max_for_flat   %||% 0.35
  )
  flags <- character(0)

  # --- insufficient data --------------------------------------------------
  if (!is.finite(scores$dxy_inside_flank_ratio) ||
      n_inside_windows < thr$min_windows) {
    return(list(
      primary_class = "insufficient_data",
      secondary_class = NA_character_,
      confidence = "low",
      reason = sprintf("only %d inside windows or non-finite ratios",
                       n_inside_windows),
      flags = flags
    ))
  }

  # --- AGE GATE -----------------------------------------------------------
  if (!is.finite(scores$oldness_score) ||
      scores$oldness_score < thr$oldness_min ||
      raw_row$dxy_inside_mean < thr$dxy_min_inside) {
    flags <- c(flags, "below_age_gate")
    return(list(
      primary_class = "young_weak_divergence",
      secondary_class = NA_character_,
      confidence = if (is.finite(scores$dxy_inside_flank_ratio) &&
                       scores$dxy_inside_flank_ratio > 1.0) "medium" else "low",
      reason = sprintf("oldness=%.2f<%.2f or dxy_inside=%.4f<%.4f — shape uninformative",
                       scores$oldness_score %||% NA, thr$oldness_min,
                       raw_row$dxy_inside_mean %||% NA, thr$dxy_min_inside),
      flags = flags
    ))
  }

  # --- SHAPE VERDICT ------------------------------------------------------
  # ordered rules — first match wins; record runner-up as secondary_class.
  signals <- list(
    asymmetric_edge = is.finite(scores$abs_log2_asymmetry) &&
                      scores$abs_log2_asymmetry >= thr$asym_high,
    locally_adapted_internal_peak_like =
                      (is.finite(scores$dxy_internal_peak_score) &&
                       scores$dxy_internal_peak_score >= thr$peak_high) ||
                      (is.finite(scores$fst_internal_peak_score) &&
                       scores$fst_internal_peak_score >= thr$peak_high),
    neutral_like_U_shape =
                      is.finite(scores$dxy_u_score) &&
                      scores$dxy_u_score >= thr$u_high &&
                      (is.finite(scores$dxy_internal_peak_score) &&
                       scores$dxy_internal_peak_score < thr$peak_high),
    locally_adapted_breakpoint_like =
                      is.finite(scores$breakpoint_adaptation_score) &&
                      scores$breakpoint_adaptation_score >= thr$inflank_high &&
                      (!is.finite(scores$dxy_u_score) ||
                       scores$dxy_u_score <  thr$u_high),
    flat_high_deep_structural_haplotype =
                      is.finite(scores$dxy_inside_flank_ratio) &&
                      scores$dxy_inside_flank_ratio >= thr$inflank_high &&
                      is.finite(scores$dxy_flatness_score) &&
                      scores$dxy_flatness_score   <= thr$flatness_max
  )
  hits <- names(signals)[unlist(signals)]

  if (length(hits) == 0L) {
    primary <- "complex_mixed"
    secondary <- NA_character_
    reason <- "no rule fired above threshold"
  } else if (length(hits) == 1L) {
    primary <- hits[1]
    secondary <- NA_character_
    reason <- sprintf("one rule fired: %s", primary)
  } else {
    # multiple — pick by a fixed priority but keep the runner-up
    priority <- c("asymmetric_edge",
                  "locally_adapted_internal_peak_like",
                  "neutral_like_U_shape",
                  "locally_adapted_breakpoint_like",
                  "flat_high_deep_structural_haplotype")
    ranked <- intersect(priority, hits)
    primary <- ranked[1]
    secondary <- ranked[2]
    if (length(hits) >= 3L) flags <- c(flags, "complex_multi_signal")
    reason <- sprintf("multiple signals fired: %s", paste(hits, collapse = "+"))
    if (length(hits) >= 3L) primary <- "complex_mixed"
  }

  # confidence
  conf <- "medium"
  if (primary == "complex_mixed") conf <- "low"
  if (n_inside_windows < 2L * thr$min_windows) conf <- "low"
  if (primary != "complex_mixed" && n_inside_windows >= 4L * thr$min_windows) {
    # margin check: how far above threshold is the deciding score?
    margin <- switch(
      primary,
      neutral_like_U_shape = scores$dxy_u_score - thr$u_high,
      locally_adapted_internal_peak_like = max(scores$dxy_internal_peak_score,
                                               scores$fst_internal_peak_score, na.rm = TRUE)
                                           - thr$peak_high,
      locally_adapted_breakpoint_like = scores$breakpoint_adaptation_score - thr$inflank_high,
      flat_high_deep_structural_haplotype = scores$dxy_inside_flank_ratio - thr$inflank_high,
      asymmetric_edge = scores$abs_log2_asymmetry - thr$asym_high,
      0
    )
    if (is.finite(margin) && margin >= 0.5) conf <- "high"
  }

  list(primary_class = primary,
       secondary_class = secondary,
       confidence = conf,
       reason = reason,
       flags = flags)
}


`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
