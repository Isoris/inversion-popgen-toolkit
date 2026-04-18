#!/usr/bin/env Rscript
# =============================================================================
# PATCH: Flashlight integration for STEP_C01b_2_merge.R
#
# TWO modifications to the merge engine:
#
# 1. SV BREAKPOINT MERGE BLOCKER (Cheat 1+2)
#    Inside run_merge_fuzzy(), before accepting a merge, check if a
#    HIGH/VERY_HIGH confidence SV breakpoint falls in the gap between
#    cores. If yes → block the merge (the SV says these are different
#    structural features separated by a real breakpoint).
#
# 2. BOUNDARY SNAP TO SV BREAKPOINT
#    After merging, if a merged region boundary is within 50 kb of a
#    precise SV breakpoint, snap the boundary to the SV position.
#    This upgrades "fuzzy PCA boundary" to "bp-resolution SV boundary".
#
# INSERT POINTS:
#   - Source sv_prior_loader.R at script top (after library() calls)
#   - Add sv_prior gate inside run_merge_fuzzy() while loop
#   - Add boundary refinement after each chromosome's merge loop
#
# Cheat 2 usage: het-DEL pileup at breakpoint positions provides
# additional evidence that a boundary is real (not an artifact).
# =============================================================================

# ─── Source sv_prior (add near top of script, after library calls) ──

# <<< INSERT after: `%||%` <- function(a, b) if (is.null(a)) b else a

fl_loader <- Sys.getenv("SV_PRIOR_LOADER", "")
if (!nzchar(fl_loader)) {
  for (p in c(
    file.path(dirname(dirname(outdir)), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path),
    file.path(dirname(outdir), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path)
  )) if (file.exists(p)) { fl_loader <- p; break }
}
.merge_has_sv_prior <- FALSE
if (nzchar(fl_loader) && file.exists(fl_loader)) {
  tryCatch({
    source(fl_loader)
    .merge_has_sv_prior <- TRUE
    message("[merge] Flashlight loader sourced: ", fl_loader)
  }, error = function(e) message("[merge] Flashlight load failed: ", e$message))
} else {
  message("[merge] No sv_prior — merge runs without SV breakpoint gates")
}

# ─── GATE 4: SV BREAKPOINT BLOCKER ──────────────────────────────────
# INSERT this block inside run_merge_fuzzy(), in the while loop,
# AFTER Gate 2 (fuzzy score) and BEFORE Gate 3 (landscape adjustment).
#
# Location: inside run_merge_fuzzy(), after:
#   ms <- compute_merge_score(dt, ref_idxs, nxt$idxs, ...)
# and before:
#   la <- landscape_adjust(ms$score, ...)

sv_prior_merge_gate <- function(chr, core_a_end_bp, core_b_start_bp) {
  # Returns: list(blocked = T/F, reason = "...")
  if (!.merge_has_sv_prior) return(list(blocked = FALSE, reason = "no_sv_prior"))

  fl <- load_sv_prior(chr)
  if (is.null(fl) || nrow(fl$inv_calls) == 0) {
    return(list(blocked = FALSE, reason = "no_inv_calls"))
  }

  # Check for HIGH/VERY_HIGH confidence SV breakpoints in the gap
  gap_start <- core_a_end_bp
  gap_end   <- core_b_start_bp

  # bp1 = left breakpoint, bp2 = right breakpoint of INV calls
  bp1_in_gap <- fl$inv_calls[
    bp1 >= gap_start & bp1 <= gap_end &
    confidence_level %in% c("HIGH", "VERY_HIGH")
  ]
  bp2_in_gap <- fl$inv_calls[
    bp2 >= gap_start & bp2 <= gap_end &
    confidence_level %in% c("HIGH", "VERY_HIGH")
  ]

  sv_bps <- rbind(
    bp1_in_gap[, .(inv_id, bp_pos = bp1, confidence_level, caller)],
    bp2_in_gap[, .(inv_id, bp_pos = bp2, confidence_level, caller)]
  )

  if (nrow(sv_bps) == 0) {
    return(list(blocked = FALSE, reason = "no_sv_bp_in_gap"))
  }

  # Additional evidence from Cheat 2: het-DELs at this breakpoint
  has_del_support <- FALSE
  if (nrow(fl$breakpoint_dels) > 0) {
    for (bp in sv_bps$bp_pos) {
      bp_dels <- fl$breakpoint_dels[abs(bp_pos - bp) <= 30000]
      if (nrow(bp_dels) > 0) {
        has_del_support <- TRUE
        break
      }
    }
  }

  best <- sv_bps[1]  # Already HIGH/VERY_HIGH
  reason <- paste0("SV_bp_blocker:", best$inv_id,
                   "@", best$bp_pos, "(", best$confidence_level, ")",
                   if (has_del_support) "+het_DEL_support" else "")

  list(blocked = TRUE, reason = reason,
       inv_id = best$inv_id, bp_pos = best$bp_pos,
       confidence = best$confidence_level,
       del_support = has_del_support)
}

# ─── USAGE: Insert in the while loop of run_merge_fuzzy() ───────────
# After computing merge score (ms) and BEFORE landscape adjustment:
#
#   # Gate 4: Flashlight SV breakpoint blocker
#   sv_gate <- sv_prior_merge_gate(chr, dt$end_bp[cur_end], dt$start_bp[nxt_start])
#   if (sv_gate$blocked) {
#     score_log[[length(score_log) + 1]] <- data.table(
#       chrom = chr, merge_family = params$name,
#       cur_end_bp = dt$end_bp[cur_end],
#       nxt_start_bp = dt$start_bp[nxt_start],
#       gap_windows = gap_size,
#       fuzzy_score = round(ms$score, 4),
#       membership_score = round(ms$membership_score, 4),
#       geometric_score = round(ms$geometric_score, 4),
#       landscape_adj = 0,
#       landscape_reason = sv_gate$reason,
#       adjusted_score = 0,
#       threshold = params$accept_threshold,
#       combination = combination,
#       decision = "blocked_sv_breakpoint"
#     )
#     break
#   }


# ─── BOUNDARY SNAP TO SV BREAKPOINTS ────────────────────────────────
# After merge regions are finalized, snap boundaries to nearby precise
# SV breakpoints. This is a POST-PROCESSING step, not a gate.
#
# INSERT after the merge QC loop, before writing output.

snap_boundaries_to_sv <- function(region_rows, chroms_with_regions) {
  if (!.merge_has_sv_prior || length(region_rows) == 0) {
    return(list(region_rows = region_rows, n_snapped = 0L, snap_log = data.table()))
  }

  SNAP_WINDOW <- 50000L  # ±50 kb
  n_snapped <- 0L
  snap_log_rows <- list()

  for (ri in seq_along(region_rows)) {
    r <- region_rows[[ri]]
    if (!is.data.table(r) || nrow(r) != 1) next

    chr <- r$chrom
    fl <- load_sv_prior(chr)
    if (is.null(fl) || nrow(fl$inv_calls) == 0) next

    orig_start <- r$start_bp
    orig_end   <- r$end_bp

    # Check left boundary
    near_left <- fl$inv_calls[
      abs(bp1 - orig_start) <= SNAP_WINDOW |
      abs(bp2 - orig_start) <= SNAP_WINDOW
    ]
    if (nrow(near_left) > 0) {
      # Find closest breakpoint to our left boundary
      all_bps_left <- c(near_left$bp1, near_left$bp2)
      dists_left <- abs(all_bps_left - orig_start)
      best_idx <- which.min(dists_left)
      if (dists_left[best_idx] <= SNAP_WINDOW) {
        new_start <- all_bps_left[best_idx]
        snap_log_rows[[length(snap_log_rows) + 1]] <- data.table(
          chrom = chr, region_id = r$region_id,
          boundary = "left", orig_bp = orig_start, snapped_bp = new_start,
          snap_dist = new_start - orig_start,
          inv_id = near_left$inv_id[ceiling(best_idx / 2)]
        )
        region_rows[[ri]]$start_bp <- new_start
        n_snapped <- n_snapped + 1L
      }
    }

    # Check right boundary
    near_right <- fl$inv_calls[
      abs(bp1 - orig_end) <= SNAP_WINDOW |
      abs(bp2 - orig_end) <= SNAP_WINDOW
    ]
    if (nrow(near_right) > 0) {
      all_bps_right <- c(near_right$bp1, near_right$bp2)
      dists_right <- abs(all_bps_right - orig_end)
      best_idx <- which.min(dists_right)
      if (dists_right[best_idx] <= SNAP_WINDOW) {
        new_end <- all_bps_right[best_idx]
        snap_log_rows[[length(snap_log_rows) + 1]] <- data.table(
          chrom = chr, region_id = r$region_id,
          boundary = "right", orig_bp = orig_end, snapped_bp = new_end,
          snap_dist = new_end - orig_end,
          inv_id = near_right$inv_id[ceiling(best_idx / 2)]
        )
        region_rows[[ri]]$end_bp <- new_end
        n_snapped <- n_snapped + 1L
      }
    }
  }

  snap_log <- if (length(snap_log_rows) > 0) rbindlist(snap_log_rows) else data.table()

  list(region_rows = region_rows, n_snapped = n_snapped, snap_log = snap_log)
}

# ─── USAGE: After the main chromosome merge loop completes ──────────
#
# # Snap boundaries to SV breakpoints (post-processing)
# snap_result <- snap_boundaries_to_sv(all_region_rows, chroms)
# all_region_rows <- snap_result$region_rows
# if (snap_result$n_snapped > 0) {
#   message("[merge] SV boundary snap: ", snap_result$n_snapped, " boundaries refined")
#   fwrite(snap_result$snap_log,
#          file.path(outdir, "sv_prior_boundary_snaps.tsv"), sep = "\t")
# }
