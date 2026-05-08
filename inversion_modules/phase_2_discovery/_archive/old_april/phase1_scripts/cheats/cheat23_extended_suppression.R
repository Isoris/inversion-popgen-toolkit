#!/usr/bin/env Rscript
# =============================================================================
# cheat23_extended_suppression.R — Extended boundary suppression
#
# BIOLOGY:
#   Inversion effects spill beyond the formal breakpoint coordinates.
#   Recombination suppression extends into flanking regions, creating an
#   "effective inversion boundary" wider than the structural breakpoint.
#   Markers just outside the inversion may show elevated Fst.
#
# INPUT:  precomp RDS, decomposition classes, boundary positions, Engine B
# OUTPUT: four-bin statistics, suppression extent, decay rate
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
DEFAULT_EXT_FLANK <- 200000L   # 200 kb for external flank bins
SCAN_STEP_KB      <- 20L       # scan outward in 20 kb steps
FST_BG_THRESHOLD  <- 0.02      # background Fst level
MAX_SCAN_KB       <- 500L      # max scan distance

# ── Four-bin statistics ───────────────────────────────────────────────

compute_four_bin_stats <- function(dt, sample_names, decomp_classes,
                                    left_bp, right_bp,
                                    external_flank = DEFAULT_EXT_FLANK) {
  if (nrow(dt) == 0 || length(decomp_classes) == 0)
    return(data.table())

  span <- right_bp - left_bp
  if (span <= 0) return(data.table())

  # Define four bins
  bins <- data.table(
    bin = c("deep_external", "near_external", "near_internal", "center"),
    start = c(left_bp - external_flank * 2, left_bp - external_flank,
              left_bp, left_bp + min(external_flank, span/4)),
    end   = c(left_bp - external_flank, left_bp,
              left_bp + min(external_flank, span/4), right_bp - min(external_flank, span/4))
  )
  # Extend: also consider right-side bins (averaged)
  bins_r <- data.table(
    bin = c("deep_external", "near_external", "near_internal", "center"),
    start = c(right_bp + external_flank, right_bp,
              right_bp - min(external_flank, span/4), left_bp + min(external_flank, span/4)),
    end   = c(right_bp + external_flank * 2, right_bp + external_flank,
              right_bp, right_bp - min(external_flank, span/4))
  )

  # Get sample groups
  ref_ids <- names(decomp_classes[decomp_classes == "HOM_REF"])
  inv_ids <- names(decomp_classes[decomp_classes == "HOM_INV"])

  results <- list()
  for (i in seq_len(nrow(bins))) {
    b <- bins[i]
    br <- bins_r[i]

    # Windows in left-side bin
    win_l <- dt[start_bp >= b$start & end_bp <= b$end]
    win_r <- dt[start_bp >= br$start & end_bp <= br$end]
    win_all <- rbind(win_l, win_r)

    if (nrow(win_all) == 0) {
      results[[i]] <- data.table(bin = b$bin, fst = NA_real_,
                                  theta_inv = NA_real_, theta_ref = NA_real_,
                                  het_contrast = NA_real_, n_windows = 0L)
      next
    }

    # Fst proxy: inv_likeness or direct per-window metric
    fst_val <- if ("inv_likeness" %in% names(win_all))
      mean(win_all$inv_likeness, na.rm = TRUE) else NA_real_

    # Het contrast (if available)
    hc <- if ("het_contrast" %in% names(win_all))
      mean(win_all$het_contrast, na.rm = TRUE) else NA_real_

    results[[i]] <- data.table(bin = b$bin, fst = round(fst_val, 4),
                                theta_inv = NA_real_, theta_ref = NA_real_,
                                het_contrast = round(hc, 4),
                                n_windows = nrow(win_all))
  }
  rbindlist(results)
}

# ── Test extended suppression ─────────────────────────────────────────

test_extended_suppression <- function(four_bin_stats) {
  if (nrow(four_bin_stats) == 0)
    return(list(has_extended_suppression = FALSE,
                suppression_extent_kb = 0,
                fst_decay_rate = NA_real_))

  deep_ext <- four_bin_stats[bin == "deep_external"]$fst
  near_ext <- four_bin_stats[bin == "near_external"]$fst

  if (length(deep_ext) == 0 || length(near_ext) == 0 ||
      is.na(deep_ext) || is.na(near_ext))
    return(list(has_extended_suppression = FALSE,
                suppression_extent_kb = 0,
                fst_decay_rate = NA_real_))

  # Extended suppression: near_external > deep_external
  has_ext <- near_ext > deep_ext * 1.5 && near_ext > FST_BG_THRESHOLD

  # Decay rate
  decay_rate <- if (has_ext && deep_ext > 0)
    round((near_ext - deep_ext) / deep_ext, 3) else NA_real_

  list(has_extended_suppression = has_ext,
       fst_near_external = round(near_ext, 4),
       fst_deep_external = round(deep_ext, 4),
       suppression_extent_kb = if (has_ext) DEFAULT_EXT_FLANK / 1000 else 0,
       fst_decay_rate = decay_rate)
}

# ── Outward scan for suppression extent ───────────────────────────────

scan_suppression_extent <- function(dt, left_bp, right_bp,
                                     step_kb = SCAN_STEP_KB,
                                     max_kb = MAX_SCAN_KB) {
  step_bp <- step_kb * 1000
  results_left  <- list()
  results_right <- list()

  # Scan leftward from left breakpoint
  for (d in seq(step_bp, max_kb * 1000, by = step_bp)) {
    win <- dt[end_bp <= left_bp - d + step_bp & start_bp >= left_bp - d]
    fst <- if (nrow(win) > 0 && "inv_likeness" %in% names(win))
      mean(win$inv_likeness, na.rm = TRUE) else NA_real_
    results_left[[length(results_left)+1]] <-
      data.table(side = "left", distance_kb = d / 1000, fst = fst)
  }

  # Scan rightward from right breakpoint
  for (d in seq(step_bp, max_kb * 1000, by = step_bp)) {
    win <- dt[start_bp >= right_bp + d - step_bp & end_bp <= right_bp + d]
    fst <- if (nrow(win) > 0 && "inv_likeness" %in% names(win))
      mean(win$inv_likeness, na.rm = TRUE) else NA_real_
    results_right[[length(results_right)+1]] <-
      data.table(side = "right", distance_kb = d / 1000, fst = fst)
  }

  scan_dt <- rbindlist(c(results_left, results_right))
  if (nrow(scan_dt) == 0) return(list(extent_kb = 0, scan_data = data.table()))

  # Find where Fst drops to background
  for (side_name in c("left", "right")) {
    side_dt <- scan_dt[side == side_name & !is.na(fst)][order(distance_kb)]
    if (nrow(side_dt) < 2) next
    bg_idx <- which(side_dt$fst <= FST_BG_THRESHOLD)[1]
    if (is.na(bg_idx)) bg_idx <- nrow(side_dt)
    scan_dt[side == side_name,
            above_bg := distance_kb <= side_dt$distance_kb[bg_idx]]
  }

  extent_kb <- max(scan_dt[above_bg == TRUE]$distance_kb, 0, na.rm = TRUE)
  list(extent_kb = extent_kb, scan_data = scan_dt)
}

# ── Diversity-differentiation decoupling ──────────────────────────────

diversity_differentiation_decoupling <- function(four_bin_stats) {
  if (nrow(four_bin_stats) < 3)
    return(list(decoupling_score = NA_real_, correlation = NA_real_))

  fst_vals <- four_bin_stats$fst
  hc_vals  <- four_bin_stats$het_contrast

  valid <- is.finite(fst_vals) & is.finite(hc_vals)
  if (sum(valid) < 3)
    return(list(decoupling_score = NA_real_, correlation = NA_real_))

  cor_val <- cor(fst_vals[valid], hc_vals[valid], method = "spearman")
  # Decoupling: low correlation = Fst and diversity are independent signals
  decoupling <- 1 - abs(cor_val)

  list(decoupling_score = round(decoupling, 3),
       correlation = round(cor_val, 3),
       interpretation = if (decoupling > 0.5) "decoupled_real_signal"
         else "coupled_simple_pattern")
}

# ── Search mode ────────────────────────────────────────────────────────

search_extended_suppression <- function(chr, zone_start, zone_end,
                                         dt = NULL, decomp_classes = NULL,
                                         ...) {
  empty <- data.table(method = "extended_suppression", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(dt) || nrow(dt) == 0) return(empty)

  mid <- as.integer((zone_start + zone_end) / 2)
  # Quick check: is inv_likeness elevated just outside the zone?
  outside_l <- dt[end_bp < zone_start & end_bp >= zone_start - 100000]
  outside_r <- dt[start_bp > zone_end & start_bp <= zone_end + 100000]
  outside_fst <- mean(c(
    if (nrow(outside_l) > 0) mean(outside_l$inv_likeness, na.rm = TRUE) else NA,
    if (nrow(outside_r) > 0) mean(outside_r$inv_likeness, na.rm = TRUE) else NA
  ), na.rm = TRUE)

  sc <- if (is.finite(outside_fst)) pmin(1, outside_fst * 5) else 0
  data.table(method = "extended_suppression", best_bp = mid,
             score = round(sc, 3), is_precise = FALSE,
             detail = paste0("outside_fst=", round(outside_fst, 4)))
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat23 <- function(chr, candidate_start, candidate_end,
                         left_bp, right_bp, dt, sample_names,
                         decomp_classes) {
  message("[cheat23] ", chr, ":", round(candidate_start/1e6,1), "-",
          round(candidate_end/1e6,1), " Mb")

  four_bin <- compute_four_bin_stats(dt, sample_names, decomp_classes,
                                      left_bp, right_bp)
  message("[cheat23] Four-bin Fst:\n",
          paste(capture.output(print(four_bin)), collapse = "\n"))

  suppression <- test_extended_suppression(four_bin)
  message("[cheat23] Extended suppression: ", suppression$has_extended_suppression,
          " | Extent: ", suppression$suppression_extent_kb, " kb")

  scan <- scan_suppression_extent(dt, left_bp, right_bp)
  message("[cheat23] Scan extent: ", scan$extent_kb, " kb")

  decoupling <- diversity_differentiation_decoupling(four_bin)
  message("[cheat23] Decoupling: ", decoupling$decoupling_score,
          " (", decoupling$interpretation, ")")

  list(four_bin = four_bin, suppression = suppression,
       scan = scan, decoupling = decoupling,
       search_result = search_extended_suppression(chr, candidate_start,
                        candidate_end, dt, decomp_classes))
}
