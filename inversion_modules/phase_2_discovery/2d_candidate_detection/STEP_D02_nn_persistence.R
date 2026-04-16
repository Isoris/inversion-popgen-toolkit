#!/usr/bin/env Rscript
# ============================================================================
# STEP_D02_nn_persistence.R — NN-scale persistence via independent staircases
# ============================================================================
#
# SESSION 7 REWRITE:
#   Instead of "check if nn0 block survives in nn40 matrix", this now:
#   1. Runs staircase on EACH NN scale independently
#   2. Matches block registries across scales by coordinate overlap
#   3. Reports survival, merge, split, disappear for each candidate
#
# This is the LIGHTWEIGHT version of 09_nn_sweep_tree.R — uses only the
# standard 4 NN scales (0, 20, 40, 80) without the full sweep.
#
# Usage:
#   source("00_config.R")
#   source("STEP_D01_staircase_boundaries.R")
#   source("STEP_D02_nn_persistence.R")
#   nn_ev <- compute_nn_persistence_v2(sim_mats_by_nn)
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


compute_nn_persistence_v2 <- function(sim_mats_by_nn,
                                       staircase_fn = detect_blocks_staircase,
                                       overlap_thresh = CFG$NN_OVERLAP_THRESH) {

  nn_scales <- sort(as.integer(names(sim_mats_by_nn)))
  cat("NN persistence v2 | scales:", paste(nn_scales, collapse = ", "), "\n")

  # ---- Run staircase on each NN independently ----
  all_blocks <- list()
  for (nn in nn_scales) {
    nn_str <- as.character(nn)
    cat("  Running staircase on nn", nn, "...\n")
    result <- staircase_fn(sim_mats_by_nn[[nn_str]])
    all_blocks[[nn_str]] <- result$blocks
    cat("    Blocks:", nrow(result$blocks), "\n")
  }

  # ---- Use finest scale (nn0 or lowest) as reference ----
  ref_nn <- as.character(nn_scales[1])
  ref_blocks <- all_blocks[[ref_nn]]

  if (nrow(ref_blocks) == 0) {
    cat("  No blocks at reference scale nn", nn_scales[1], "\n")
    return(data.table())
  }

  # ---- Match each ref block to blocks at other scales ----
  results <- list()

  for (r in seq_len(nrow(ref_blocks))) {
    rs <- ref_blocks$start[r]
    re <- ref_blocks$end[r]
    bid <- ref_blocks$block_id[r]

    row <- data.table(candidate_id = bid, ref_start = rs, ref_end = re,
                      ref_width = re - rs + 1L)

    for (nn in nn_scales) {
      nn_str <- as.character(nn)
      nn_blocks <- all_blocks[[nn_str]]
      col_sim <- paste0("nn", nn, "_match")
      col_status <- paste0("nn", nn, "_status")

      if (nrow(nn_blocks) == 0) {
        row[[col_sim]] <- NA_character_
        row[[col_status]] <- "absent"
        next
      }

      # Find best overlapping blocks at this scale
      best_ids <- c()
      for (bi in seq_len(nrow(nn_blocks))) {
        ov <- compute_recip_ov(rs, re, nn_blocks$start[bi], nn_blocks$end[bi])
        if (ov >= overlap_thresh * 0.5) {
          best_ids <- c(best_ids, bi)
        }
      }

      if (length(best_ids) == 0) {
        row[[col_sim]] <- "-"
        row[[col_status]] <- "disappears"
      } else if (length(best_ids) == 1) {
        row[[col_sim]] <- as.character(best_ids)
        row[[col_status]] <- "stable"
      } else {
        # Multiple matches — check if ref block maps to multiple blocks (split)
        # or multiple ref blocks map to one (merge)
        row[[col_sim]] <- paste(best_ids, collapse = "+")
        row[[col_status]] <- "splits"
      }
    }

    # ---- Determine overall survival ----
    status_cols <- grep("_status$", names(row), value = TRUE)
    statuses <- as.character(row[1, ..status_cols])

    # Highest NN where it's still stable or present
    max_survive_nn <- 0L
    for (nn in rev(nn_scales)) {
      col <- paste0("nn", nn, "_status")
      if (col %in% names(row) && row[[col]] %in% c("stable", "splits")) {
        max_survive_nn <- nn
        break
      }
    }

    row$max_survive_nn <- max_survive_nn
    row$survives_nn40  <- max_survive_nn >= 40
    row$survives_nn80  <- max_survive_nn >= 80

    # Overall topology
    if (all(statuses %in% c("stable", ""))) {
      row$nn_topology <- "stable"
    } else if (any(statuses == "disappears")) {
      first_disappear <- nn_scales[which(statuses == "disappears")[1]]
      row$nn_topology <- paste0("disappears_nn", first_disappear)
    } else if (any(statuses == "splits")) {
      row$nn_topology <- "splits"
    } else {
      row$nn_topology <- "complex"
    }

    results[[r]] <- row
  }

  return(rbindlist(results, fill = TRUE))
}


# ---- Helper ----
compute_recip_ov <- function(s1, e1, s2, e2) {
  w1 <- e1 - s1 + 1; w2 <- e2 - s2 + 1
  if (w1 <= 0 || w2 <= 0) return(0)
  ov <- max(0, min(e1, e2) - max(s1, s2) + 1)
  return(min(ov / w1, ov / w2))
}
