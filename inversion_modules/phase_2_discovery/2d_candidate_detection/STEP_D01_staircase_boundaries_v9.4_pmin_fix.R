#!/usr/bin/env Rscript
# ============================================================================
# STEP_D01_staircase_boundaries.R — Vote-based boundary detection (v9.3)
# ============================================================================
#
# SIMPLE IDEA:
#   1. For each window i, read sim_mat[i, ] and find WHERE it steps down.
#      Don't care about the exact height — just the POSITIONS of the steps.
#   2. Count votes: position X is a boundary if MANY windows see a step at X.
#      100 windows agree = real. 3 windows agree = noise.
#   3. Blocks = intervals between high-vote boundaries.
#
# WHY THIS WORKS:
#   Inside a block, all windows see the SAME steps at the SAME positions.
#   Window 500 sees a step at 540. Window 501 sees a step at 540. Window 502
#   sees a step at 540. That's 40+ votes for a boundary at 540.
#   A noise dip at position 347 is seen only by windows 345-349. That's 5 votes.
#   The vote count IS the noise filter. No MAD, no adaptive threshold.
#
# WHY RANKS NOT HEIGHTS:
#   The sim_mat values depend on many factors (sample composition, window
#   content, local diversity). Height 0.82 on LG01 might mean the same thing
#   as height 0.65 on LG10. But a step-down is a step-down regardless of
#   the absolute level. Converting to ranks (step 1, step 2, step 3) makes
#   the detection robust to the overall similarity scale.
#
# STEP WIDTH MATTERS:
#   A wide plateau (100 bins) means 100 windows vote for the same boundary
#   at its edge. A narrow plateau (5 bins) means only 5 votes. Wide plateaus
#   = strong blocks. Narrow plateaus = weak or noisy. The vote count
#   naturally encodes plateau width.
#
# OUTPUTS:
#   $step_table      — every step seen by every window (raw)
#   $vote_profile    — votes per position (how many windows see a step here)
#   $boundaries      — positions with high vote counts (real boundaries)
#   $blocks          — intervals between boundaries, scored
#   $rank_table      — per-window rank structure (for comparison)
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

detect_blocks_staircase <- function(smat,
                                     min_block_width = CFG$STAIR_MIN_BLOCK_WIDTH,
                                     min_drop        = CFG$STAIR_MIN_DROP,
                                     persist_n       = CFG$STAIR_PERSIST_N,
                                     smooth_span     = CFG$STAIR_SMOOTH_SPAN,
                                     ema_alpha       = CFG$STAIR_EMA_ALPHA,
                                     max_steps       = 8L,
                                     vote_blur       = 2L,
                                     chrom           = NA_character_,
                                     global_window_ids = NULL) {

  N <- nrow(smat)

  # Identity args: chrom + global_window_ids are optional. When provided they
  # get stamped onto every output data.table so downstream modules can join
  # the staircase results back to the precompute window grid without needing
  # to match on bp coordinates. When NULL, the rows still carry the local
  # integer window_id (1..N) as before.
  if (!is.null(global_window_ids)) {
    if (length(global_window_ids) != N) {
      warning("[staircase] global_window_ids length (", length(global_window_ids),
              ") != smat rows (", N, "); ignoring")
      global_window_ids <- NULL
    }
  }

  cat("Staircase detector v9.3 | N =", N, "\n")
  if (!is.na(chrom)) cat("  chrom =", chrom, "\n")
  cat("  min_drop =", min_drop, "| persist_n =", persist_n,
      "| vote_blur =", vote_blur, "\n")

  # ==== STEP 1: For each window, find step positions ====
  # Read row i of sim_mat, smooth it, find where sustained drops happen.
  # Record the POSITION of each step, its height before/after, and the
  # final background level (similarity to distant windows = family LD signal).
  cat("Step 1: Finding step positions per window...\n")

  step_list   <- vector("list", N)
  bg_level    <- numeric(N)  # per-window background level (far-distance similarity)

  for (i in seq_len(N)) {
    row_vals <- smat[i, ]
    row_vals[!is.finite(row_vals)] <- NA_real_

    right_info <- find_step_positions(row_vals, center = i, direction = "right",
                                       min_drop = min_drop, persist_n = persist_n,
                                       smooth_span = smooth_span, max_steps = max_steps,
                                       ema_alpha = ema_alpha)

    left_info <- find_step_positions(row_vals, center = i, direction = "left",
                                      min_drop = min_drop, persist_n = persist_n,
                                      smooth_span = smooth_span, max_steps = max_steps,
                                      ema_alpha = ema_alpha)

    step_list[[i]] <- list(right = right_info, left = left_info)

    # Background = average of both directions' background
    bg_r <- right_info$background
    bg_l <- left_info$background
    bg_vals <- c(bg_r, bg_l)
    bg_vals <- bg_vals[is.finite(bg_vals)]
    bg_level[i] <- if (length(bg_vals) > 0) mean(bg_vals) else NA_real_

    if (i %% 2000 == 0) cat("  window", i, "/", N, "\n")
  }

  # Flag family LD windows: background level is elevated compared to chromosome median
  chr_bg_median <- median(bg_level, na.rm = TRUE)
  chr_bg_mad    <- mad(bg_level, constant = 1.4826, na.rm = TRUE)
  family_ld_flag <- bg_level > chr_bg_median + 2 * chr_bg_mad

  cat("  Done. Family LD windows:", sum(family_ld_flag, na.rm = TRUE),
      "/ ", N, "\n")
  cat("  Chr background: median=", round(chr_bg_median, 3),
      " MAD=", round(chr_bg_mad, 3), "\n")

  # ==== STEP 2: Vote — count how many windows see a step at each position ====
  cat("Step 2: Counting boundary votes...\n")

  vote_right <- numeric(N)
  vote_left  <- numeric(N)

  for (i in seq_len(N)) {
    for (sp in step_list[[i]]$right$positions) {
      lo <- max(1L, sp - vote_blur)
      hi <- min(N, sp + vote_blur)
      vote_right[lo:hi] <- vote_right[lo:hi] + 1
    }
    for (sp in step_list[[i]]$left$positions) {
      lo <- max(1L, sp - vote_blur)
      hi <- min(N, sp + vote_blur)
      vote_left[lo:hi] <- vote_left[lo:hi] + 1
    }
  }

  # Combined vote: total evidence at each position from either direction.
  #
  # BUGFIX (2026-04-18, CRITICAL): the old v9.3 code used
  #   vote_combined <- pmin(vote_right, vote_left)
  # with the stated intent "a boundary should be seen from BOTH directions".
  # The intent is wrong: at a block's LEFT edge (say pos L), only windows
  # INSIDE the block vote for L as a LEFT-step; windows OUTSIDE the block
  # don't see L as any kind of right-step (their right-steps land further
  # right). So vote_left[L] is high and vote_right[L] is zero — pmin = 0.
  # Same symmetric problem at the RIGHT edge. pmin silently zeros out
  # every real boundary found on both synthetic and real data.
  #
  # The correct aggregation is total evidence at each position:
  #   vote_combined <- vote_right + vote_left
  # This accumulates contributions from both walking directions without
  # requiring the same position to be seen as both a right-step and a
  # left-step — which it almost never is for a real block boundary.
  #
  # Max votes at a strong boundary is approximately:
  #   width_of_block windows (looking outward from inside the block)
  # which is consistent with the "votes_threshold = min_block_width"
  # heuristic in the comment above Step 3.
  vote_combined <- vote_right + vote_left

  # Smooth votes to merge nearby peaks
  vote_smooth <- running_med_vec(vote_combined, span = smooth_span)

  cat("  Max votes:", max(vote_smooth), "at position",
      which.max(vote_smooth), "\n")

  # ==== STEP 3: Find boundary positions (peaks in the vote profile) ====
  # Real boundaries = local maxima in the vote profile with enough votes.
  # "Enough" = the vote count from a block of width W is roughly W
  # (each window inside the block votes for the boundary at its edge).
  # So min_block_width IS the vote threshold.
  cat("Step 3: Identifying boundaries...\n")

  vote_threshold <- max(min_block_width, 3L)
  boundaries <- find_vote_peaks(vote_smooth, threshold = vote_threshold,
                                 min_gap = min_block_width)

  cat("  Boundaries found:", length(boundaries),
      "(threshold:", vote_threshold, "votes)\n")

  # ==== STEP 4: Build blocks from intervals between boundaries ====
  cat("Step 4: Building blocks...\n")

  blocks <- build_blocks_from_boundaries(boundaries, smat, N,
                                          min_block_width, min_drop,
                                          vote_smooth = vote_smooth)

  cat("  Blocks:", nrow(blocks), "\n")

  # ==== STEP 5: Build step table (for inspection — keeps EVERYTHING) ====
  step_table <- build_step_table(step_list, N, blocks = blocks,
                                  chrom = chrom,
                                  global_window_ids = global_window_ids)

  # ==== STEP 6: Build rank table (per-window rank structure + background) ====
  rank_table <- build_rank_table(step_list, N, bg_level,
                                  chrom = chrom,
                                  global_window_ids = global_window_ids)

  # ==== STEP 7: Parent/child nesting ====
  if (nrow(blocks) > 1) {
    blocks <- assign_parents(blocks)
  }

  cat("  Final blocks:", nrow(blocks), "\n")

  # ---- Stamp identity fields on blocks + vote_profile ----
  # blocks: a block spans windows start..end, so stamp paired gwids. chrom
  # is a scalar column (same value for every row).
  if (nrow(blocks) > 0) {
    blocks[, chrom := chrom]
    if (!is.null(global_window_ids)) {
      # Guard each lookup — start/end are 1-based local positions on smat.
      safe_gwid <- function(i) if (is.finite(i) && i >= 1 && i <= length(global_window_ids))
                                 global_window_ids[i] else NA_integer_
      blocks[, start_global_window_id := vapply(start, safe_gwid, integer(1))]
      blocks[, end_global_window_id   := vapply(end,   safe_gwid, integer(1))]
    } else {
      blocks[, start_global_window_id := NA_integer_]
      blocks[, end_global_window_id   := NA_integer_]
    }
    # Reorder so chrom + identity cols come first (readable columns leading)
    id_cols  <- intersect(c("chrom", "block_id", "start_global_window_id", "end_global_window_id"), names(blocks))
    rest     <- setdiff(names(blocks), id_cols)
    setcolorder(blocks, c(id_cols, rest))
  }

  vote_profile <- data.table(
    chrom         = chrom,
    position      = seq_len(N),
    global_window_id = if (!is.null(global_window_ids)) global_window_ids else NA_integer_,
    vote_right    = vote_right,
    vote_left     = vote_left,
    vote_combined = vote_combined,
    vote_smooth   = vote_smooth,
    bg_level      = round(bg_level, 4),
    family_ld     = family_ld_flag
  )

  # Boundaries as a proper data.table (was a bare integer vector). Back-compat
  # note: if an older downstream caller expected an integer vector, it should
  # read `$position` instead. We keep the same name because it's clearer.
  boundaries_dt <- if (length(boundaries) > 0) {
    data.table(
      chrom = chrom,
      position = boundaries,
      global_window_id = if (!is.null(global_window_ids)) global_window_ids[boundaries] else NA_integer_,
      votes = vote_smooth[boundaries]
    )
  } else {
    data.table(chrom = character(), position = integer(),
               global_window_id = integer(), votes = numeric())
  }

  # ---- Summary statistics (non-destructive; new table in the return list) ----
  # Per-chromosome stats so a log scan / QC pass has meaningful numbers.
  n_steps_per_window <- vapply(step_list, function(sl) {
    length(sl$right$positions) + length(sl$left$positions)
  }, integer(1))

  step_drops <- if (nrow(step_table) > 0) step_table$drop_size[is.finite(step_table$drop_size)] else numeric(0)
  block_widths <- if (nrow(blocks) > 0) blocks$width else integer(0)
  block_heights <- if (nrow(blocks) > 0) blocks$height else numeric(0)
  block_contrasts <- if (nrow(blocks) > 0) blocks$contrast else numeric(0)
  vote_peaks <- vote_smooth[vote_smooth > 0]

  summary_stats <- data.table(
    chrom = chrom,
    n_windows = N,

    # Step-detection stats
    n_steps_total         = nrow(step_table),
    n_steps_per_win_mean  = round(mean(n_steps_per_window), 3),
    n_steps_per_win_median = stats::median(n_steps_per_window),
    n_steps_per_win_max   = max(n_steps_per_window),
    n_windows_with_steps  = sum(n_steps_per_window > 0),
    n_windows_no_steps    = sum(n_steps_per_window == 0),
    frac_windows_with_steps = if (N > 0) round(sum(n_steps_per_window > 0) / N, 4) else NA_real_,

    # Drop-size distribution (per step)
    drop_size_mean   = if (length(step_drops) > 0) round(mean(step_drops), 4) else NA_real_,
    drop_size_median = if (length(step_drops) > 0) round(stats::median(step_drops), 4) else NA_real_,
    drop_size_p95    = if (length(step_drops) > 0) round(stats::quantile(step_drops, 0.95, names = FALSE), 4) else NA_real_,
    drop_size_max    = if (length(step_drops) > 0) round(max(step_drops), 4) else NA_real_,

    # Vote-profile stats
    vote_smooth_max     = if (length(vote_peaks) > 0) round(max(vote_peaks), 2) else 0,
    vote_smooth_mean    = if (length(vote_peaks) > 0) round(mean(vote_peaks), 2) else 0,
    vote_smooth_p95     = if (length(vote_peaks) > 0) round(stats::quantile(vote_peaks, 0.95, names = FALSE), 2) else 0,
    n_positions_voted   = length(vote_peaks),

    # Boundary stats
    n_boundaries       = length(boundaries),
    boundary_votes_mean = if (length(boundaries) > 0) round(mean(vote_smooth[boundaries]), 2) else NA_real_,
    boundary_votes_min  = if (length(boundaries) > 0) round(min(vote_smooth[boundaries]), 2) else NA_real_,

    # Block-shape distribution
    n_blocks            = nrow(blocks),
    block_width_mean    = if (length(block_widths) > 0) round(mean(block_widths), 1) else NA_real_,
    block_width_median  = if (length(block_widths) > 0) stats::median(block_widths) else NA_real_,
    block_width_max     = if (length(block_widths) > 0) as.integer(max(block_widths)) else NA_integer_,
    block_width_min     = if (length(block_widths) > 0) as.integer(min(block_widths)) else NA_integer_,
    block_height_mean   = if (length(block_heights) > 0) round(mean(block_heights), 4) else NA_real_,
    block_contrast_mean = if (length(block_contrasts) > 0) round(mean(block_contrasts), 4) else NA_real_,
    n_nested_blocks     = if ("parent_id" %in% names(blocks)) sum(!is.na(blocks$parent_id)) else 0L,

    # Background / family-LD stats
    chr_bg_median      = round(chr_bg_median, 4),
    chr_bg_mad         = round(chr_bg_mad, 4),
    n_family_ld_windows = sum(family_ld_flag, na.rm = TRUE),
    frac_family_ld      = if (N > 0) round(sum(family_ld_flag, na.rm = TRUE) / N, 4) else NA_real_
  )

  # Per-chrom log line so tail -f log.txt shows the highlights
  cat(sprintf("[staircase-summary] %s  N=%d  blocks=%d (mean_width=%s max=%s)  steps/win=%.2f  n_bdy=%d (votes_mean=%s)  family_LD=%d (%.1f%%)\n",
    if (is.na(chrom)) "(no_chrom)" else chrom,
    N, nrow(blocks),
    if (length(block_widths) > 0) as.character(as.integer(round(mean(block_widths)))) else "NA",
    if (length(block_widths) > 0) as.character(as.integer(max(block_widths))) else "NA",
    mean(n_steps_per_window),
    length(boundaries),
    if (length(boundaries) > 0) sprintf("%.1f", mean(vote_smooth[boundaries])) else "NA",
    sum(family_ld_flag, na.rm = TRUE),
    100 * if (N > 0) sum(family_ld_flag, na.rm = TRUE) / N else NA_real_
  ))

  return(list(
    blocks        = blocks,
    step_table    = step_table,
    rank_table    = rank_table,
    vote_profile  = vote_profile,
    boundaries    = boundaries_dt,
    summary_stats = summary_stats,
    artifacts     = data.table(),
    chrom         = chrom
  ))
}


# ============================================================================
# FIND STEP POSITIONS — for one window, one direction
# ============================================================================
# Read the similarity values from center outward, find where sustained
# drops happen. Return POSITIONS AND HEIGHTS of each step, plus the
# final background level.
#
# Returns list with:
#   $positions — integer vector of step positions (absolute)
#   $heights   — numeric vector: similarity level BEFORE each step
#   $drops     — numeric vector: similarity level AFTER each step
#   $background — the final level after all steps (far from center)

find_step_positions <- function(row_vals, center, direction,
                                 min_drop, persist_n, smooth_span,
                                 max_steps, ema_alpha = 0.05) {
  N <- length(row_vals)
  result <- list(positions = integer(0), heights = numeric(0),
                 drops = numeric(0), background = NA_real_)

  if (direction == "right") {
    idx <- (center + 1):N
  } else {
    idx <- (center - 1):1
  }

  if (length(idx) < persist_n + 1) return(result)

  # Get values in traversal order
  vals <- row_vals[idx]
  vals[!is.finite(vals)] <- NA_real_

  # Smooth
  if (length(vals) > smooth_span * 2) {
    vals <- running_med_vec(vals, span = smooth_span)
  }

  n <- length(vals)
  current_level <- median(vals[seq_len(min(5L, n))], na.rm = TRUE)
  if (!is.finite(current_level) || current_level <= 0) return(result)

  positions <- integer(0)
  heights   <- numeric(0)
  drops     <- numeric(0)

  # ---- Fixed step finder: detect EACH intermediate level ----
  # Problem with old code: a cascade like 0.85 → 0.60 → 0.20 got merged
  # into one step (0.85 → 0.20) because below_count kept incrementing
  # across BOTH drops. The 0.60 shelf was invisible.
  #
  # Fix: when we're in a "below" streak and the value drops FURTHER
  # (to a new level significantly below the FIRST drop level), record
  # the first step and start tracking the second drop separately.
  #
  # This correctly separates:
  #   0.85 → 0.60 (step 1 at position of first drop)
  #   0.60 → 0.20 (step 2 at position of second drop)
  # even if the 0.60 shelf is only 1-2 bins wide.

  n_steps <- 0L
  first_drop_k <- 0L        # position where current drop streak started
  first_drop_level <- NA_real_  # level at start of current drop streak

  for (k in seq_len(n)) {
    v <- vals[k]
    if (is.na(v)) next

    if (v < current_level - min_drop) {
      # We're below the current plateau

      if (first_drop_k == 0L) {
        # Start of a new drop streak
        first_drop_k <- k
        first_drop_level <- v
      } else {
        # Already in a drop streak. Check: did we drop to a NEW level?
        # (significantly below the first drop level)
        if (v < first_drop_level - min_drop) {
          # YES — the value dropped FURTHER. Record the first step
          # (from current_level to first_drop_level) at first_drop_k.
          # Then start tracking the deeper drop.
          n_steps <- n_steps + 1L
          step_abs_pos <- idx[min(first_drop_k, length(idx))]
          positions <- c(positions, step_abs_pos)
          heights   <- c(heights, current_level)
          drops     <- c(drops, first_drop_level)

          if (n_steps >= max_steps) break

          # The intermediate level becomes the new "current"
          current_level <- first_drop_level
          first_drop_k <- k
          first_drop_level <- v
        }
        # Otherwise: still in the same drop, just accumulating
      }

      # Check if drop has persisted enough
      bins_in_drop <- k - first_drop_k + 1L
      if (bins_in_drop >= persist_n) {
        # Confirmed drop. Record the step.
        n_steps <- n_steps + 1L
        step_abs_pos <- idx[min(first_drop_k, length(idx))]
        positions <- c(positions, step_abs_pos)
        heights   <- c(heights, current_level)
        drops     <- c(drops, v)

        if (n_steps >= max_steps) break

        current_level <- v
        first_drop_k <- 0L
        first_drop_level <- NA_real_
      }

    } else {
      # Value is back within plateau range
      if (first_drop_k > 0L) {
        # We were in a drop streak that RECOVERED (transient dip)
        # Don't record as step. Reset.
        first_drop_k <- 0L
        first_drop_level <- NA_real_
      }
      # Gently track current level
      current_level <- current_level * (1 - ema_alpha) + v * ema_alpha
    }
  }

  # Background = the level after the last step (median of last 10% of values)
  tail_n <- max(5L, n %/% 10)
  tail_vals <- vals[seq(max(1, n - tail_n + 1), n)]
  tail_vals <- tail_vals[is.finite(tail_vals)]
  background <- if (length(tail_vals) > 0) median(tail_vals) else current_level

  list(positions = positions, heights = heights,
       drops = drops, background = background)
}


# ============================================================================
# FIND PEAKS IN VOTE PROFILE
# ============================================================================
# Local maxima above threshold, with minimum gap between peaks.

find_vote_peaks <- function(votes, threshold, min_gap) {
  n <- length(votes)
  if (n < 3) return(integer(0))

  peaks <- integer(0)

  for (i in 2:(n - 1)) {
    if (votes[i] >= threshold &&
        votes[i] >= votes[i - 1] &&
        votes[i] >= votes[i + 1]) {
      # Check minimum gap from last peak
      if (length(peaks) == 0 || (i - peaks[length(peaks)]) >= min_gap) {
        peaks <- c(peaks, i)
      } else if (votes[i] > votes[peaks[length(peaks)]]) {
        # Replace last peak if this one is higher
        peaks[length(peaks)] <- i
      }
    }
  }

  return(peaks)
}


# ============================================================================
# BUILD BLOCKS FROM BOUNDARY POSITIONS
# ============================================================================
# TWO-PASS LOGIC:
#   Pass 1: Standard blocks (min_block_width = 5+). Outer and large blocks.
#   Pass 2: Tiny blocks (min_block_width = 2) but ONLY if:
#           (a) they sit inside an already-detected block from Pass 1
#           (b) their boundary has high votes (≥ 2× standard threshold)
#   This catches 3-window rare inversions nested inside large common ones.

build_blocks_from_boundaries <- function(boundaries, smat, N,
                                          min_block_width, min_drop,
                                          vote_smooth = NULL) {
  all_bounds <- sort(unique(c(1L, boundaries, N)))
  bg_sim <- median(smat[upper.tri(smat)], na.rm = TRUE)

  # ---- Pass 1: Standard width filter ----
  blocks <- data.table()
  block_id <- 0L

  for (k in seq_len(length(all_bounds) - 1)) {
    bi <- all_bounds[k]
    be <- all_bounds[k + 1]
    bw <- be - bi + 1L

    if (bw < min_block_width) next

    sub <- smat[bi:be, bi:be]
    ut <- sub[upper.tri(sub)]
    ut <- ut[is.finite(ut)]
    if (length(ut) == 0) next

    block_height <- median(ut)
    if (block_height < bg_sim + min_drop) next

    block_id <- block_id + 1L
    step_left  <- compute_boundary_drop(smat, bi, "left", 5)
    step_right <- compute_boundary_drop(smat, be, "right", 5)

    blocks <- rbind(blocks, data.table(
      block_id        = block_id,
      start           = bi,
      end             = be,
      width           = bw,
      start_bp        = bin_to_bp(bi),
      end_bp          = bin_to_bp(be),
      start_mb        = bin_to_mb(bi),
      end_mb          = bin_to_mb(be),
      height          = round(block_height, 4),
      bg_sim          = round(bg_sim, 4),
      contrast        = round(block_height - bg_sim, 4),
      step_down_left  = round(step_left, 4),
      step_down_right = round(step_right, 4),
      parent_id       = NA_integer_,
      n_artifacts     = 0L
    ))
  }

  # ---- Pass 2: Tiny blocks inside existing blocks ----
  # Catch 2-4 window inversions nested inside larger ones.
  # Require: (a) inside a Pass 1 block, (b) high vote support if available.
  min_tiny <- 2L
  high_vote_thresh <- max(min_block_width * 2, 10L)

  for (k in seq_len(length(all_bounds) - 1)) {
    bi <- all_bounds[k]
    be <- all_bounds[k + 1]
    bw <- be - bi + 1L

    # Only consider intervals that were TOO SMALL for Pass 1
    if (bw >= min_block_width || bw < min_tiny) next

    # Must be inside an existing block
    inside_parent <- FALSE
    for (r in seq_len(nrow(blocks))) {
      if (bi >= blocks$start[r] && be <= blocks$end[r]) {
        inside_parent <- TRUE
        break
      }
    }
    if (!inside_parent) next

    # If we have vote data, require high votes at BOTH boundaries
    if (!is.null(vote_smooth) && length(vote_smooth) >= N) {
      left_votes  <- if (bi >= 1 && bi <= length(vote_smooth)) vote_smooth[bi] else 0
      right_votes <- if (be >= 1 && be <= length(vote_smooth)) vote_smooth[be] else 0
      if (min(left_votes, right_votes) < high_vote_thresh) next
    }

    # Check similarity: must be HIGHER than parent's interior
    sub <- smat[bi:be, bi:be]
    ut <- sub[upper.tri(sub)]
    ut <- ut[is.finite(ut)]
    if (length(ut) == 0) next

    block_height <- median(ut)
    # Must be elevated: brighter than the surrounding parent block
    parent_height <- blocks$height[r]  # r from the loop above
    if (block_height < parent_height + min_drop * 0.5) next

    block_id <- block_id + 1L
    step_left  <- compute_boundary_drop(smat, bi, "left", min(3, bw))
    step_right <- compute_boundary_drop(smat, be, "right", min(3, bw))

    blocks <- rbind(blocks, data.table(
      block_id        = block_id,
      start           = bi,
      end             = be,
      width           = bw,
      start_bp        = bin_to_bp(bi),
      end_bp          = bin_to_bp(be),
      start_mb        = bin_to_mb(bi),
      end_mb          = bin_to_mb(be),
      height          = round(block_height, 4),
      bg_sim          = round(bg_sim, 4),
      contrast        = round(block_height - bg_sim, 4),
      step_down_left  = round(step_left, 4),
      step_down_right = round(step_right, 4),
      parent_id       = NA_integer_,
      n_artifacts     = 0L
    ))
  }

  return(blocks)
}


# ============================================================================
# BUILD STEP TABLE (for inspection — EVERY step from EVERY window)
# ============================================================================
# NOTHING is removed. Low-vote steps inside blocks are flagged as
# "internal_feature" not deleted. Could be recombinant double crossover
# or assembly error — downstream modules decide.

build_step_table <- function(step_list, N, blocks = NULL,
                              chrom = NA_character_,
                              global_window_ids = NULL) {
  rows <- list()
  gwid_of <- function(i) if (!is.null(global_window_ids) && i >= 1 && i <= length(global_window_ids))
                           global_window_ids[i] else NA_integer_
  for (i in seq_len(N)) {
    sl <- step_list[[i]]
    for (j in seq_along(sl$right$positions)) {
      rows[[length(rows) + 1]] <- data.table(
        chrom         = chrom,
        window_id     = i,
        global_window_id = gwid_of(i),
        direction     = "right",
        step_rank     = j,
        step_position = sl$right$positions[j],
        step_global_window_id = gwid_of(sl$right$positions[j]),
        height_before = round(sl$right$heights[j], 4),
        height_after  = round(sl$right$drops[j], 4),
        drop_size     = round(sl$right$heights[j] - sl$right$drops[j], 4)
      )
    }
    for (j in seq_along(sl$left$positions)) {
      rows[[length(rows) + 1]] <- data.table(
        chrom         = chrom,
        window_id     = i,
        global_window_id = gwid_of(i),
        direction     = "left",
        step_rank     = j,
        step_position = sl$left$positions[j],
        step_global_window_id = gwid_of(sl$left$positions[j]),
        height_before = round(sl$left$heights[j], 4),
        height_after  = round(sl$left$drops[j], 4),
        drop_size     = round(sl$left$heights[j] - sl$left$drops[j], 4)
      )
    }
  }
  st <- if (length(rows) > 0) rbindlist(rows) else data.table()

  # Flag steps that are INSIDE detected blocks
  # These are internal features — keep in table, don't use for block boundaries
  if (nrow(st) > 0 && !is.null(blocks) && nrow(blocks) > 0) {
    st[, inside_block := FALSE]
    st[, containing_block_id := NA_integer_]
    for (bi in seq_len(nrow(blocks))) {
      b_start <- blocks$start[bi]
      b_end   <- blocks$end[bi]
      b_id    <- blocks$block_id[bi]
      st[step_position >= b_start & step_position <= b_end,
         `:=`(inside_block = TRUE, containing_block_id = b_id)]
    }
  }

  return(st)
}


# ============================================================================
# BUILD RANK TABLE (per-window rank signature + background level)
# ============================================================================
# For each window, its "rank signature" = the sorted list of step positions.
# Windows with the same rank signature are in the same structural context.
# Background level = family LD indicator.

build_rank_table <- function(step_list, N, bg_level,
                              chrom = NA_character_,
                              global_window_ids = NULL) {
  rows <- list()
  for (i in seq_len(N)) {
    sl <- step_list[[i]]
    all_steps <- sort(unique(c(sl$right$positions, sl$left$positions)))
    n_steps <- length(all_steps)

    sig <- if (n_steps > 0) paste(all_steps, collapse = "_") else "none"
    gwid <- if (!is.null(global_window_ids) && i >= 1 && i <= length(global_window_ids))
              global_window_ids[i] else NA_integer_

    rows[[length(rows) + 1]] <- data.table(
      chrom       = chrom,
      window_id   = i,
      global_window_id = gwid,
      n_steps     = n_steps,
      signature   = sig,
      steps_right = paste(sl$right$positions, collapse = ","),
      steps_left  = paste(sl$left$positions, collapse = ","),
      bg_level    = round(bg_level[i], 4)
    )
  }
  if (length(rows) > 0) rbindlist(rows) else data.table()
}


# ============================================================================
# BOUNDARY DROP
# ============================================================================

compute_boundary_drop <- function(smat, pos, direction, depth = 5) {
  N <- nrow(smat)
  if (direction == "left") {
    inside  <- max(1, pos)     : min(N, pos + depth - 1)
    outside <- max(1, pos - depth) : max(1, pos - 1)
  } else {
    inside  <- max(1, pos - depth + 1) : min(N, pos)
    outside <- min(N, pos + 1) : min(N, pos + depth)
  }
  if (length(inside) < 1 || length(outside) < 1) return(0)

  cross <- smat[inside, outside]
  cross <- cross[is.finite(cross)]
  if (length(inside) > 1) {
    in_sub <- smat[inside, inside]
    in_vals <- in_sub[upper.tri(in_sub)]
    in_vals <- in_vals[is.finite(in_vals)]
  } else {
    in_vals <- 1
  }

  in_med <- median(in_vals, na.rm = TRUE)
  cr_med <- median(cross, na.rm = TRUE)
  if (is.na(in_med) || is.na(cr_med)) return(0)
  in_med - cr_med
}


# ============================================================================
# PARENT/CHILD NESTING
# ============================================================================

assign_parents <- function(blocks) {
  setorder(blocks, -width)
  n <- nrow(blocks)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      if (blocks$start[i] >= blocks$start[j] &&
          blocks$end[i]   <= blocks$end[j] &&
          blocks$width[i] <  blocks$width[j]) {
        if (is.na(blocks$parent_id[i]) ||
            blocks$width[blocks$block_id == blocks$parent_id[i]] > blocks$width[j]) {
          blocks$parent_id[i] <- blocks$block_id[j]
        }
      }
    }
  }
  return(blocks)
}


# ============================================================================
# HELPERS
# ============================================================================

running_med_vec <- function(x, span = 3) {
  n <- length(x)
  if (n <= span) return(x)
  half <- span %/% 2
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - half)
    hi <- min(n, i + half)
    out[i] <- median(x[lo:hi], na.rm = TRUE)
  }
  return(out)
}


# ============================================================================
# CLI ENTRY POINT
# ============================================================================

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  sim_mat_file <- args[1]
  outdir <- if (length(args) > 1) args[2] else "."

  if (!exists("CFG")) {
    CFG <- list(STAIR_MIN_BLOCK_WIDTH = 5L, STAIR_MIN_DROP = 0.03,
                STAIR_PERSIST_N = 3L, STAIR_SMOOTH_SPAN = 3L,
                STAIR_EMA_ALPHA = 0.05, STAIR_VOTE_BLUR = 2L,
                WINDOW_SIZE_BP = 50000L)
    bin_to_bp <- function(b, w = 50000L) (b - 1) * w
    bin_to_mb <- function(b, w = 50000L) (b - 1) * w / 1e6
  }

  # Check for --compact flag
  compact_mode <- "--compact" %in% args

  cat("Loading sim_mat:", sim_mat_file, "\n")
  obj <- readRDS(sim_mat_file)

  # Pull identity fields from the RDS when available. Standard precomp object
  # contains $sim_mat, $chrom, $dt (with global_window_id). Bare matrix input
  # still works — we just don't get the identity stamps.
  smat <- if (is.list(obj) && "sim_mat" %in% names(obj)) obj$sim_mat else obj
  chrom_val <- if (is.list(obj) && "chrom" %in% names(obj)) as.character(obj$chrom) else NA_character_
  gwids <- NULL
  if (is.list(obj) && !is.null(obj$dt) && "global_window_id" %in% names(obj$dt)) {
    gwids <- obj$dt$global_window_id
    if (length(gwids) != nrow(smat)) {
      cat("  WARNING: obj$dt rows (", length(gwids), ") != smat rows (",
          nrow(smat), ") — ignoring gwids\n")
      gwids <- NULL
    }
  }

  result <- detect_blocks_staircase(smat, chrom = chrom_val,
                                    global_window_ids = gwids)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  bn <- gsub("\\.rds$", "", basename(sim_mat_file))

  fwrite(result$blocks, file.path(outdir, paste0(bn, "_blocks.tsv")), sep = "\t")

  # Step table: compact mode only saves steps at high-vote positions
  if (compact_mode && nrow(result$step_table) > 0) {
    # Only keep steps at positions that got ≥ vote_threshold votes
    vote_thresh <- max(CFG$STAIR_MIN_BLOCK_WIDTH, 3L)
    high_vote_pos <- result$vote_profile[vote_smooth >= vote_thresh]$position
    compact_steps <- result$step_table[step_position %in% high_vote_pos]
    fwrite(compact_steps, file.path(outdir, paste0(bn, "_steps.tsv.gz")), sep = "\t")
    cat("  Compact steps:", nrow(compact_steps), "of", nrow(result$step_table), "\n")
  } else {
    fwrite(result$step_table, file.path(outdir, paste0(bn, "_steps.tsv.gz")), sep = "\t")
  }

  fwrite(result$rank_table, file.path(outdir, paste0(bn, "_ranks.tsv.gz")), sep = "\t")
  fwrite(result$vote_profile, file.path(outdir, paste0(bn, "_votes.tsv.gz")), sep = "\t")
  fwrite(result$boundaries, file.path(outdir, paste0(bn, "_boundaries.tsv")), sep = "\t")
  fwrite(result$summary_stats, file.path(outdir, paste0(bn, "_summary.tsv")), sep = "\t")
  saveRDS(result, file.path(outdir, paste0(bn, "_staircase.rds")))
  cat("Done. Summary TSV:", file.path(outdir, paste0(bn, "_summary.tsv")), "\n")
}
