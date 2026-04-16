#!/usr/bin/env Rscript
# ============================================================================
# STEP_D09_nn_sweep_tree.R — Interval tree of blocks across NN-scale sweep
# ============================================================================
#
# THE PRIMARY CLASSIFIER.  No separate "evidence module" needed.
#
# Architecture:
#   1. Sweep NN scales from large (coarse) to small (fine)
#   2. At each NN scale, run the per-window staircase
#   3. As NN decreases, blocks SPLIT → build tree from containment
#   4. Each tree node gets nn_birth (coarsest NN where it appears separately)
#   5. nn_birth IS the classification:
#      - High nn_birth (>200)  = strong structure = INVERSION
#      - Medium nn_birth (80-200) = moderate = CANDIDATE
#      - Low nn_birth (<80)    = weak = FAMILY_LD or NOISE
#
# Usage:
#   source("00_config.R")
#   source("STEP_D01_staircase_boundaries.R")
#   source("STEP_D09_nn_sweep_tree.R")
#   tree <- build_nn_sweep_tree(sim_mats_by_nn)  # named list: "0", "20", ...
#   # OR: provide a function to compute sim_mat on the fly
#   tree <- build_nn_sweep_tree_lazy(compute_simmat_fn, nn_scales)
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================================
# MAIN: Build the interval tree from pre-loaded sim_mats
# ============================================================================

build_nn_sweep_tree <- function(sim_mats_by_nn, staircase_fn = detect_blocks_staircase) {

  nn_scales <- sort(as.integer(names(sim_mats_by_nn)), decreasing = TRUE)
  cat("NN sweep tree | scales:", paste(nn_scales, collapse = ", "), "\n")

  # ---- Run staircase at each NN scale ----
  all_blocks <- list()
  for (nn in nn_scales) {
    nn_str <- as.character(nn)
    cat("\n=== NN =", nn, "===\n")
    smat <- sim_mats_by_nn[[nn_str]]
    result <- staircase_fn(smat)
    blocks <- result$blocks
    if (nrow(blocks) > 0) {
      blocks[, nn_scale := nn]
    }
    all_blocks[[nn_str]] <- blocks
    cat("  Blocks at nn", nn, ":", nrow(blocks), "\n")
  }

  # ---- Build tree by matching across scales ----
  cat("\nBuilding interval tree...\n")
  tree <- match_across_scales(all_blocks, nn_scales)

  cat("Tree built:", nrow(tree), "nodes\n")
  return(tree)
}


# ============================================================================
# MATCH BLOCKS ACROSS NN SCALES
# ============================================================================

match_across_scales <- function(all_blocks, nn_scales,
                                 overlap_thresh = CFG$NN_OVERLAP_THRESH) {

  # Start from the coarsest (highest) NN scale
  # Each block at the coarsest scale is a root node
  # Walk down to finer scales, matching children to parents

  # Initialize tree nodes from the coarsest scale
  coarsest_nn <- nn_scales[1]
  coarsest_blocks <- all_blocks[[as.character(coarsest_nn)]]

  if (nrow(coarsest_blocks) == 0) {
    # No blocks at coarsest → start from next
    for (nn in nn_scales) {
      b <- all_blocks[[as.character(nn)]]
      if (nrow(b) > 0) {
        coarsest_nn <- nn
        coarsest_blocks <- b
        break
      }
    }
  }

  if (nrow(coarsest_blocks) == 0) {
    return(data.table(
      node_id = integer(0), start = integer(0), end = integer(0),
      width = integer(0), nn_birth = integer(0), nn_death = character(0),
      parent_node = integer(0), children = character(0), height = numeric(0),
      topology = character(0)
    ))
  }

  # Initialize nodes
  node_id_counter <- 0L
  tree_nodes <- list()

  for (r in seq_len(nrow(coarsest_blocks))) {
    node_id_counter <- node_id_counter + 1L
    tree_nodes[[node_id_counter]] <- list(
      node_id    = node_id_counter,
      start      = coarsest_blocks$start[r],
      end        = coarsest_blocks$end[r],
      width      = coarsest_blocks$width[r],
      nn_birth   = coarsest_nn,
      nn_death   = "leaf",
      parent_node = NA_integer_,
      children   = integer(0),
      height     = coarsest_blocks$height[r],
      start_bp   = coarsest_blocks$start_bp[r],
      end_bp     = coarsest_blocks$end_bp[r],
      start_mb   = coarsest_blocks$start_mb[r],
      end_mb     = coarsest_blocks$end_mb[r],
      topology   = "root"
    )
  }

  # Walk from coarse to fine
  prev_nn <- coarsest_nn
  prev_active <- seq_len(node_id_counter)  # active leaf node IDs

  for (scale_idx in 2:length(nn_scales)) {
    nn <- nn_scales[scale_idx]
    current_blocks <- all_blocks[[as.character(nn)]]

    if (nrow(current_blocks) == 0) next

    # Match current blocks to active parent leaves
    matched_parents <- integer(nrow(current_blocks))  # 0 = unmatched
    parent_match_count <- integer(length(prev_active))  # how many children

    for (cb in seq_len(nrow(current_blocks))) {
      cb_start <- current_blocks$start[cb]
      cb_end   <- current_blocks$end[cb]

      best_parent <- 0L
      best_overlap <- 0

      for (pi in seq_along(prev_active)) {
        pid <- prev_active[pi]
        p_start <- tree_nodes[[pid]]$start
        p_end   <- tree_nodes[[pid]]$end

        # Reciprocal overlap
        ov <- compute_reciprocal_overlap(cb_start, cb_end, p_start, p_end)
        if (ov > best_overlap) {
          best_overlap <- ov
          best_parent <- pi
        }
      }

      if (best_overlap >= overlap_thresh) {
        matched_parents[cb] <- best_parent
        parent_match_count[best_parent] <- parent_match_count[best_parent] + 1L
      }
    }

    # Process matches
    new_active <- integer(0)

    for (pi in seq_along(prev_active)) {
      pid <- prev_active[pi]
      n_children <- parent_match_count[pi]

      if (n_children == 0) {
        # DISAPPEARS at this NN — block was an artifact of oversmoothing
        tree_nodes[[pid]]$nn_death <- paste0("disappears_nn", nn)
        tree_nodes[[pid]]$topology <- "disappears"
        # Keep as leaf but mark it
        new_active <- c(new_active, pid)

      } else if (n_children == 1) {
        # STABLE — same block, maybe slight boundary shift
        child_idx <- which(matched_parents == pi)
        # Update boundaries of existing node (slight refinement)
        tree_nodes[[pid]]$start <- current_blocks$start[child_idx]
        tree_nodes[[pid]]$end   <- current_blocks$end[child_idx]
        tree_nodes[[pid]]$width <- current_blocks$width[child_idx]
        tree_nodes[[pid]]$height <- current_blocks$height[child_idx]
        new_active <- c(new_active, pid)

      } else {
        # SPLIT — parent splits into multiple children
        tree_nodes[[pid]]$nn_death <- paste0("splits_nn", nn)
        tree_nodes[[pid]]$topology <- "splits"

        child_indices <- which(matched_parents == pi)
        child_ids <- integer(0)

        for (ci in child_indices) {
          node_id_counter <- node_id_counter + 1L
          new_id <- node_id_counter

          tree_nodes[[new_id]] <- list(
            node_id    = new_id,
            start      = current_blocks$start[ci],
            end        = current_blocks$end[ci],
            width      = current_blocks$width[ci],
            nn_birth   = nn,
            nn_death   = "leaf",
            parent_node = pid,
            children   = integer(0),
            height     = current_blocks$height[ci],
            start_bp   = current_blocks$start_bp[ci],
            end_bp     = current_blocks$end_bp[ci],
            start_mb   = current_blocks$start_mb[ci],
            end_mb     = current_blocks$end_mb[ci],
            topology   = "child"
          )
          child_ids <- c(child_ids, new_id)
          new_active <- c(new_active, new_id)
        }

        tree_nodes[[pid]]$children <- child_ids
      }
    }

    # NEW blocks that appeared at this NN (no parent match)
    unmatched <- which(matched_parents == 0L)
    for (cb in unmatched) {
      # Check if it's inside an existing parent
      cb_start <- current_blocks$start[cb]
      cb_end   <- current_blocks$end[cb]

      parent_found <- NA_integer_
      for (pid in prev_active) {
        if (cb_start >= tree_nodes[[pid]]$start &&
            cb_end   <= tree_nodes[[pid]]$end) {
          parent_found <- pid
          break
        }
      }

      node_id_counter <- node_id_counter + 1L
      tree_nodes[[node_id_counter]] <- list(
        node_id    = node_id_counter,
        start      = cb_start,
        end        = cb_end,
        width      = current_blocks$width[cb],
        nn_birth   = nn,
        nn_death   = "leaf",
        parent_node = if (!is.na(parent_found)) parent_found else NA_integer_,
        children   = integer(0),
        height     = current_blocks$height[cb],
        start_bp   = current_blocks$start_bp[cb],
        end_bp     = current_blocks$end_bp[cb],
        start_mb   = current_blocks$start_mb[cb],
        end_mb     = current_blocks$end_mb[cb],
        topology   = "novel"
      )
      new_active <- c(new_active, node_id_counter)
    }

    prev_active <- unique(new_active)
    cat("  nn", nn, "→", length(new_active), "active nodes,",
        length(unmatched), "novel\n")
  }

  # Convert to data.table
  tree_dt <- rbindlist(lapply(tree_nodes, function(nd) {
    data.table(
      node_id     = nd$node_id,
      start       = nd$start,
      end         = nd$end,
      width       = nd$width,
      start_mb    = nd$start_mb,
      end_mb      = nd$end_mb,
      nn_birth    = nd$nn_birth,
      nn_death    = nd$nn_death,
      parent_node = nd$parent_node,
      children    = paste(nd$children, collapse = ","),
      height      = nd$height,
      topology    = nd$topology
    )
  }))

  # Classify based on nn_birth
  tree_dt[, classification := classify_nn_birth(nn_birth)]

  return(tree_dt)
}


# ============================================================================
# NN BIRTH CLASSIFIER
# ============================================================================

classify_nn_birth <- function(nn_birth) {
  ifelse(nn_birth >= 200, "INVERSION",
    ifelse(nn_birth >= 80, "CANDIDATE",
      ifelse(nn_birth >= 40, "WEAK_CANDIDATE",
        "FAMILY_LD")))
}


# ============================================================================
# RECIPROCAL OVERLAP
# ============================================================================

compute_reciprocal_overlap <- function(s1, e1, s2, e2) {
  w1 <- e1 - s1 + 1
  w2 <- e2 - s2 + 1
  if (w1 <= 0 || w2 <= 0) return(0)

  ov_start <- max(s1, s2)
  ov_end   <- min(e1, e2)
  ov <- max(0, ov_end - ov_start + 1)

  # Reciprocal: min of (overlap/w1, overlap/w2)
  return(min(ov / w1, ov / w2))
}


# ============================================================================
# TRACKING TABLE — one row per inversion system
# ============================================================================

build_tracking_table <- function(all_blocks, nn_scales,
                                  overlap_thresh = CFG$NN_OVERLAP_THRESH) {
  # Match blocks across adjacent scales and build tracking rows
  # Each row = one "system" tracked across all NN scales

  # Start from finest (nn0) — find all blocks
  finest <- nn_scales[length(nn_scales)]
  finest_blocks <- all_blocks[[as.character(finest)]]

  if (nrow(finest_blocks) == 0) return(data.table())

  tracking <- data.table(
    label  = paste0("L", sprintf("%03d", seq_len(nrow(finest_blocks)))),
    status = rep("UNKNOWN", nrow(finest_blocks))
  )

  # For each NN scale, find which finest block matches
  for (nn in nn_scales) {
    nn_str <- as.character(nn)
    nn_blocks <- all_blocks[[nn_str]]
    col_name <- paste0("nn", nn, "_block")
    tracking[[col_name]] <- NA_character_

    if (nrow(nn_blocks) == 0) next

    for (fi in seq_len(nrow(finest_blocks))) {
      fs <- finest_blocks$start[fi]
      fe <- finest_blocks$end[fi]

      matches <- c()
      for (bi in seq_len(nrow(nn_blocks))) {
        ov <- compute_reciprocal_overlap(fs, fe,
                                          nn_blocks$start[bi], nn_blocks$end[bi])
        if (ov >= overlap_thresh * 0.5) {
          matches <- c(matches, bi)
        }
      }

      if (length(matches) > 0) {
        tracking[[col_name]][fi] <- paste(matches, collapse = "+")
      } else {
        tracking[[col_name]][fi] <- "-"
      }
    }
  }

  # Determine status
  nn_cols <- grep("^nn.*_block$", names(tracking), value = TRUE)
  for (fi in seq_len(nrow(tracking))) {
    vals <- as.character(tracking[fi, ..nn_cols])

    # Count how many scales have a match
    n_present <- sum(vals != "-" & !is.na(vals))

    if (n_present == length(nn_cols)) {
      tracking$status[fi] <- "PERSISTS"
    } else if (n_present >= length(nn_cols) * 0.5) {
      tracking$status[fi] <- "WEAKENS"
    } else {
      tracking$status[fi] <- "DISAPPEARS"
    }

    # Check for merges: if nn40 block matches multiple nn0 blocks
    for (col in nn_cols) {
      v <- as.character(tracking[[col]][fi])
      if (grepl("\\+", v)) {
        tracking$status[fi] <- paste0("MERGES_at_", gsub("_block", "", col))
        break
      }
    }
  }

  return(tracking)
}


# ============================================================================
# SUMMARY PRINTER
# ============================================================================

print_tree_summary <- function(tree_dt) {
  cat("\n=== NN Sweep Tree Summary ===\n")
  cat("Total nodes:", nrow(tree_dt), "\n")
  cat("Root nodes:", sum(is.na(tree_dt$parent_node)), "\n")
  cat("Leaf nodes:", sum(tree_dt$nn_death == "leaf"), "\n")

  cat("\nClassification breakdown:\n")
  tab <- table(tree_dt$classification)
  for (cls in names(tab)) {
    cat("  ", cls, ":", tab[cls], "\n")
  }

  cat("\nTopology breakdown:\n")
  tab2 <- table(tree_dt$topology)
  for (tp in names(tab2)) {
    cat("  ", tp, ":", tab2[tp], "\n")
  }

  # Show strong candidates
  strong <- tree_dt[classification %in% c("INVERSION", "CANDIDATE")]
  if (nrow(strong) > 0) {
    cat("\nStrong candidates (nn_birth ≥ 80):\n")
    setorder(strong, -nn_birth)
    for (r in seq_len(min(20, nrow(strong)))) {
      cat(sprintf("  node %d: %d-%d (%.1f-%.1f Mb) | nn_birth=%d | %s\n",
                  strong$node_id[r], strong$start[r], strong$end[r],
                  strong$start_mb[r], strong$end_mb[r],
                  strong$nn_birth[r], strong$classification[r]))
    }
  }
}
