#!/usr/bin/env Rscript
# ============================================================================
# STEP_D09_nn_sweep_tree.R — Merged NN persistence + sweep tree (v1.4)
# ============================================================================
#
# v1.4 CHANGELOG
# --------------
# Span-aware classification. Previously a tree node was classified purely
# by `nn_birth` (persistence across NN scales). A 25-kb block that happened
# to survive to nn320 was still labeled INVERSION — biologically
# implausible for low-coverage WGS. v1.4 adds `span_mb = end_mb - start_mb`
# as a secondary criterion:
#
#   INVERSION     : nn_birth >= 200 AND span_mb >= inv_min_span_mb  (default 0.5 Mb)
#   CANDIDATE     : nn_birth >= 80  AND span_mb >= cand_min_span_mb (default 0.1 Mb)
#                   (or INVERSION-class node demoted for being too small)
#   WEAK_CANDIDATE: nn_birth >= 40
#                   (or CANDIDATE demoted for being too small)
#   FAMILY_LD     : everything else
#
# Thresholds are tunable via new args `inv_min_span_mb`, `cand_min_span_mb`.
# Pass them to build_nn_sweep_tree() or keep the defaults.
#
# v1.3 (historical): precomputed_blocks parameter.
# v1.2 (historical): dt is mandatory.
# ============================================================================
#
# Replaces the old split of:
#   STEP_D02_nn_persistence.R  (flat per-block survival flags, "fast")
#   STEP_D09_nn_sweep_tree.R   (full persistence tree, "deep")
#
# ARCHITECTURE
# ------------
# One function, two modes. Deep is the default because:
#   - Inputs are identical (list of sim_mats by NN scale).
#   - First step (staircase at each scale) is shared.
#   - Tree construction is cheap compared to the staircase work.
#   - Every column D02 used to emit is either directly present or
#     derivable from the tree.
#
# The function returns a single result list:
#   $summary          — flat per-block table (always produced)
#   $tree             — tree node table (deep mode only; NULL in fast)
#   $blocks_by_scale  — raw staircase output per scale
#   $mode             — "deep" or "fast" (after auto-downgrade, if any)
#   $nn_scales        — sorted integer vector of scales actually used
#
# SUMMARY COLUMNS (always present)
#   block_id          — ID of the block from finest-scale staircase (reference)
#   ref_start,        — block coordinates at the finest (nn0) scale
#   ref_end,
#   ref_width
#   nn<k>_status      — one per scale in ladder:
#                         "stable"     match with reciprocal overlap >= thresh
#                         "splits"     multiple matches at this scale
#                         "disappears" no matches at this scale
#                         "absent"     sim_mat_nn<k> not available
#   nn<k>_match       — matching block IDs at scale k ("3" or "3+4" on split)
#   survives_nn<k>    — logical, one per scale: max_survive_nn >= k
#   max_survive_nn    — coarsest scale at which this block still has a match
#   nn_topology       — stable / disappears_nn<k> / splits_at_nn<k> / complex
#
# DEEP-MODE EXTRAS (attached to summary, also exposed as separate tree table)
#   nn_birth          — coarsest scale at which this block first appears
#                       as its own distinct unit (tree-derived)
#   nn_death          — "leaf" if persists to finest scale, else "splits_nn<k>"
#                       or "disappears_nn<k>"
#   parent_node       — node_id of parent in tree, NA if root
#   classification    — INVERSION (nn_birth >= 200),
#                       CANDIDATE (80-199),
#                       WEAK_CANDIDATE (40-79),
#                       FAMILY_LD (<40)
#
# USAGE
#   source("STEP_D01_staircase_boundaries_v9.3.R")   # provides detect_blocks_staircase
#   source("STEP_D09_nn_sweep_tree.R")
#
#   # sim_mats is a named list: names are NN scales as strings
#   #   "0" = raw, "20" = nn20, "40" = nn40, ..., "320" = nn320
#   result <- build_nn_sweep_tree(sim_mats)   # deep by default
#
#   fwrite(result$summary, "nn_summary_<chr>.tsv", sep = "\t")
#   if (!is.null(result$tree)) {
#     fwrite(result$tree, "nn_tree_<chr>.tsv", sep = "\t")
#   }
#
# DEPENDENCIES
#   - data.table
#   - detect_blocks_staircase() from STEP_D01 (or any function with the
#     same contract: takes a sim_mat, returns list with $blocks data.table
#     containing block_id, start, end, width, start_bp, end_bp,
#     start_mb, end_mb, height)
#   - CFG$NN_OVERLAP_THRESH (default 0.70 if not set)
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

build_nn_sweep_tree <- function(sim_mats_by_nn,
                                 dt,
                                 staircase_fn = NULL,
                                 mode = c("deep", "fast"),
                                 overlap_thresh = NULL,
                                 nn_birth_inversion_threshold = 200L,
                                 nn_birth_candidate_threshold = 80L,
                                 nn_birth_weak_threshold = 40L,
                                 inv_min_span_mb  = 0.5,
                                 cand_min_span_mb = 0.1,
                                 precomputed_blocks = NULL,
                                 verbose = TRUE) {

  mode <- match.arg(mode)

  # v9.6: dt is mandatory. D01 needs it to look up real window coordinates.
  if (missing(dt) || is.null(dt)) {
    stop("[build_nn_sweep_tree] dt is required. Pass pc$dt from the ",
         "precomp RDS.")
  }
  if (!all(c("start_bp", "end_bp") %in% names(dt))) {
    stop("[build_nn_sweep_tree] dt must contain start_bp and end_bp columns")
  }

  if (is.null(staircase_fn)) {
    if (!exists("detect_blocks_staircase", mode = "function")) {
      stop("staircase_fn not provided and detect_blocks_staircase() not found. ",
           "Source STEP_D01 first, or pass staircase_fn explicitly.")
    }
    staircase_fn <- detect_blocks_staircase
  }

  if (is.null(overlap_thresh)) {
    overlap_thresh <- if (exists("CFG") && !is.null(CFG$NN_OVERLAP_THRESH))
      CFG$NN_OVERLAP_THRESH else 0.70
  }

  if (length(sim_mats_by_nn) == 0) {
    stop("build_nn_sweep_tree: sim_mats_by_nn is empty")
  }

  nn_scales <- sort(as.integer(names(sim_mats_by_nn)))
  if (any(is.na(nn_scales))) {
    stop("build_nn_sweep_tree: sim_mats_by_nn names must be integer NN scales ",
         "as strings (e.g. '0', '20', '40', '80'). Got: ",
         paste(names(sim_mats_by_nn), collapse = ", "))
  }

  # Safety rail: tree needs >=3 scales
  if (mode == "deep" && length(nn_scales) < 3) {
    warning("build_nn_sweep_tree: only ", length(nn_scales),
            " scales available (need >=3 for tree). Falling back to fast mode.")
    mode <- "fast"
  }

  # Also need >=2 scales for anything to do at all
  if (length(nn_scales) < 2) {
    warning("build_nn_sweep_tree: only 1 scale provided, no persistence to compute. ",
            "Returning empty result.")
    return(list(
      summary = data.table(),
      tree = NULL,
      blocks_by_scale = list(),
      mode = "fast",
      nn_scales = nn_scales
    ))
  }

  if (verbose) {
    cat("build_nn_sweep_tree | mode =", mode,
        "| scales:", paste(nn_scales, collapse = ", "),
        "| overlap_thresh =", overlap_thresh, "\n")
  }

  # ==== Step 1: per-scale block detection ====
  # By default this runs D01 staircase on every sim_mat (slow: ~3 min each).
  # If precomputed_blocks is supplied (a named list "0","20",... of block
  # data.tables, typically loaded from disk via fread of the nn*.tsv files
  # D01 already wrote), we skip step 1 entirely and go straight to tree
  # construction. Block counts and coordinates are validated against the
  # nn_scales before use.
  if (!is.null(precomputed_blocks)) {
    if (!is.list(precomputed_blocks)) {
      stop("[build_nn_sweep_tree] precomputed_blocks must be a named list")
    }
    missing_scales <- setdiff(as.character(nn_scales), names(precomputed_blocks))
    if (length(missing_scales) > 0) {
      stop("[build_nn_sweep_tree] precomputed_blocks missing scales: ",
           paste(missing_scales, collapse = ", "))
    }
    # Validate each table has required columns + stamp nn_scale
    all_blocks <- list()
    for (nn in nn_scales) {
      nn_str <- as.character(nn)
      bl <- as.data.table(precomputed_blocks[[nn_str]])
      if (nrow(bl) > 0) {
        req_cols <- c("start", "end", "width", "start_bp", "end_bp")
        missing_cols <- setdiff(req_cols, names(bl))
        if (length(missing_cols) > 0) {
          stop("[build_nn_sweep_tree] precomputed_blocks[['", nn_str,
               "']] missing cols: ", paste(missing_cols, collapse = ", "))
        }
        bl[, nn_scale := nn]
      }
      all_blocks[[nn_str]] <- bl
    }
    if (verbose) {
      cat("Using precomputed blocks (skipping staircase rerun):\n")
      for (nn in nn_scales) {
        cat(sprintf("  nn%d: %d blocks\n", nn, nrow(all_blocks[[as.character(nn)]])))
      }
    }
  } else {
    all_blocks <- run_staircase_at_each_scale(sim_mats_by_nn, staircase_fn,
                                               verbose, dt = dt)
  }

  # ==== Step 2: build flat per-block summary (always) ====
  if (verbose) cat("Building per-block summary...\n")
  summary <- build_per_block_summary(all_blocks, nn_scales, overlap_thresh)

  if (verbose) {
    cat("  Summary rows:", nrow(summary), "\n")
    if (nrow(summary) > 0) {
      cat("  Topology breakdown:\n")
      print(table(summary$nn_topology))
    }
  }

  # ==== Step 3: build persistence tree (deep mode only) ====
  tree <- NULL
  if (mode == "deep") {
    if (verbose) cat("Building persistence tree...\n")
    tree <- match_across_scales_build_tree(all_blocks, nn_scales, overlap_thresh)

    if (!is.null(tree) && nrow(tree) > 0) {
      # v1.4: compute per-node span in Mb and pass to classifier
      tree[, span_mb := end_mb - start_mb]

      # Classify tree nodes (persistence + span)
      tree[, classification := classify_nn_birth(
        nn_birth,
        span_mb          = span_mb,
        inv_thresh       = nn_birth_inversion_threshold,
        cand_thresh      = nn_birth_candidate_threshold,
        weak_thresh      = nn_birth_weak_threshold,
        inv_min_span_mb  = inv_min_span_mb,
        cand_min_span_mb = cand_min_span_mb
      )]

      if (verbose) {
        cat("  Tree nodes:", nrow(tree), "\n")
        cat("  Classification:\n")
        print(table(tree$classification))
      }

      # ==== Step 4: attach tree fields back to summary ====
      summary <- attach_tree_fields_to_summary(summary, tree)
    } else {
      if (verbose) cat("  Tree is empty (no blocks at any scale)\n")
    }
  }

  list(
    summary         = summary,
    tree            = tree,
    blocks_by_scale = all_blocks,
    mode            = mode,
    nn_scales       = nn_scales
  )
}


# ============================================================================
# STEP 1: run staircase at each NN scale
# ============================================================================

run_staircase_at_each_scale <- function(sim_mats_by_nn, staircase_fn, verbose,
                                         dt) {
  if (missing(dt) || is.null(dt)) {
    stop("[run_staircase_at_each_scale] dt is required (v9.6: no fallback)")
  }
  nn_scales <- sort(as.integer(names(sim_mats_by_nn)))
  all_blocks <- list()

  for (nn in nn_scales) {
    nn_str <- as.character(nn)
    if (verbose) cat("  Running staircase on nn", nn, "... ", sep = "")

    t0 <- proc.time()
    res <- tryCatch(
      staircase_fn(sim_mats_by_nn[[nn_str]], dt = dt),
      error = function(e) {
        warning("staircase failed on nn", nn, ": ", conditionMessage(e))
        list(blocks = data.table())
      }
    )
    blocks <- if (!is.null(res$blocks)) as.data.table(res$blocks) else data.table()
    if (nrow(blocks) > 0) {
      blocks[, nn_scale := nn]
    }
    all_blocks[[nn_str]] <- blocks

    dt_s <- round((proc.time() - t0)[3], 1)
    if (verbose) cat("blocks =", nrow(blocks), " (", dt_s, "s)\n", sep = "")
  }

  all_blocks
}


# ============================================================================
# STEP 2: per-block summary (D02-equivalent)
# ============================================================================
#
# Reference = finest scale (nn0 if available, else smallest scale).
# Each reference block gets one row with matches + status at every scale.

build_per_block_summary <- function(all_blocks, nn_scales, overlap_thresh) {
  ref_nn <- nn_scales[1L]
  ref_blocks <- all_blocks[[as.character(ref_nn)]]

  if (is.null(ref_blocks) || nrow(ref_blocks) == 0) {
    return(data.table())
  }

  # Build one row per reference block, with a column per scale
  rows <- vector("list", nrow(ref_blocks))

  for (r in seq_len(nrow(ref_blocks))) {
    rs <- ref_blocks$start[r]
    re <- ref_blocks$end[r]
    bid <- ref_blocks$block_id[r]

    row <- data.table(
      block_id   = bid,
      ref_start  = rs,
      ref_end    = re,
      ref_width  = re - rs + 1L
    )

    max_survive <- ref_nn  # starts at finest scale

    for (nn in nn_scales) {
      nn_str <- as.character(nn)
      nn_blocks <- all_blocks[[nn_str]]
      match_col  <- paste0("nn", nn, "_match")
      status_col <- paste0("nn", nn, "_status")

      if (is.null(nn_blocks) || nrow(nn_blocks) == 0) {
        row[[match_col]]  <- "-"
        row[[status_col]] <- "absent"
        next
      }

      # Find all blocks at scale `nn` with reciprocal overlap >= threshold
      matches <- integer(0)
      for (bi in seq_len(nrow(nn_blocks))) {
        ov <- compute_reciprocal_overlap(rs, re,
                                          nn_blocks$start[bi],
                                          nn_blocks$end[bi])
        if (ov >= overlap_thresh) {
          matches <- c(matches, nn_blocks$block_id[bi])
        }
      }

      if (length(matches) == 0) {
        # Relax slightly to catch near-misses
        for (bi in seq_len(nrow(nn_blocks))) {
          ov <- compute_reciprocal_overlap(rs, re,
                                            nn_blocks$start[bi],
                                            nn_blocks$end[bi])
          if (ov >= overlap_thresh * 0.5) {
            matches <- c(matches, nn_blocks$block_id[bi])
          }
        }
        if (length(matches) == 0) {
          row[[match_col]]  <- "-"
          row[[status_col]] <- "disappears"
          next
        } else {
          # Loose match only — still not a "stable" label
          row[[match_col]]  <- paste(matches, collapse = "+")
          row[[status_col]] <- if (length(matches) > 1) "splits" else "weak_match"
          if (nn > max_survive) max_survive <- nn
          next
        }
      }

      row[[match_col]]  <- paste(matches, collapse = "+")
      row[[status_col]] <- if (length(matches) > 1) "splits" else "stable"
      if (nn > max_survive) max_survive <- nn
    }

    row[, max_survive_nn := max_survive]
    row[, nn_topology   := summarise_topology(row, nn_scales)]

    rows[[r]] <- row
  }

  summary <- rbindlist(rows, fill = TRUE)

  # Add survives_nn<k> logical columns for every scale in ladder
  for (nn in nn_scales) {
    col <- paste0("survives_nn", nn)
    summary[[col]] <- summary$max_survive_nn >= nn
  }

  summary
}


# Topology summariser: looks at per-scale status columns and produces
# a one-word label describing what happens to this block across scales.
summarise_topology <- function(row, nn_scales) {
  statuses <- vapply(nn_scales, function(nn) {
    v <- row[[paste0("nn", nn, "_status")]]
    if (is.null(v) || length(v) == 0) NA_character_ else as.character(v[[1]])
  }, character(1))

  # Drop absent (scale unavailable) for classification purposes
  valid <- statuses[statuses != "absent"]

  if (length(valid) == 0) return("no_data")

  if (all(valid == "stable")) return("stable")

  # First disappearance tells us the death scale
  idx_dis <- which(statuses == "disappears")
  if (length(idx_dis) > 0) {
    first_dis_scale <- nn_scales[idx_dis[1]]
    # But if something splits before disappearing, call it complex
    idx_split <- which(statuses == "splits")
    if (length(idx_split) > 0 && idx_split[1] < idx_dis[1]) return("complex")
    return(paste0("disappears_nn", first_dis_scale))
  }

  idx_split <- which(statuses == "splits")
  if (length(idx_split) > 0) {
    return(paste0("splits_at_nn", nn_scales[idx_split[1]]))
  }

  "complex"
}


# ============================================================================
# STEP 3: persistence tree (D09-equivalent)
# ============================================================================
#
# Walks from coarsest scale to finest. At each step, match blocks at
# the current scale to blocks at the previous (coarser) scale.
#   - 1-to-1 match        -> node stays stable, boundaries refined
#   - many-to-1 match     -> parent splits; each new block gets a new
#                            node with parent = old node
#   - 0 matches           -> node disappears at current scale
#   - novel block         -> either attached to an existing parent by
#                            coordinate containment, or new root
#
# Each node carries:
#   node_id, start, end, width, start_bp, end_bp, start_mb, end_mb, height
#   nn_birth    coarsest scale at which it first appeared distinctly
#   nn_death    "leaf" if still alive at finest scale, else "splits_nn<k>"
#               or "disappears_nn<k>"
#   parent_node node_id of parent, or NA for root
#   children    comma-separated list of child node_ids
#   topology    "root" | "stable" | "splits" | "disappears" | "novel"

match_across_scales_build_tree <- function(all_blocks, nn_scales, overlap_thresh) {
  # Walk coarse -> fine
  nn_sorted_coarse_first <- sort(nn_scales, decreasing = TRUE)

  tree_nodes <- list()
  node_id_counter <- 0L

  # Initialise tree with the coarsest scale — every block is a root
  coarsest <- nn_sorted_coarse_first[1]
  coarsest_blocks <- all_blocks[[as.character(coarsest)]]

  if (is.null(coarsest_blocks) || nrow(coarsest_blocks) == 0) {
    # Try next scales until we find one with blocks
    for (nn in nn_sorted_coarse_first[-1]) {
      cb <- all_blocks[[as.character(nn)]]
      if (!is.null(cb) && nrow(cb) > 0) {
        coarsest <- nn
        coarsest_blocks <- cb
        break
      }
    }
  }

  if (is.null(coarsest_blocks) || nrow(coarsest_blocks) == 0) {
    return(data.table())   # no blocks anywhere
  }

  # Active node IDs = nodes still alive at previous iteration
  prev_active <- integer(0)

  for (b in seq_len(nrow(coarsest_blocks))) {
    node_id_counter <- node_id_counter + 1L
    tree_nodes[[node_id_counter]] <- make_node(
      id = node_id_counter,
      block_row = coarsest_blocks[b],
      nn_birth = coarsest,
      parent = NA_integer_,
      topology = "root"
    )
    prev_active <- c(prev_active, node_id_counter)
  }

  # Walk finer scales
  remaining <- nn_sorted_coarse_first[nn_sorted_coarse_first < coarsest]

  for (nn in remaining) {
    current_blocks <- all_blocks[[as.character(nn)]]
    if (is.null(current_blocks) || nrow(current_blocks) == 0) {
      # Nothing at this scale — every active node "disappears" here
      for (pid in prev_active) {
        if (is.na(tree_nodes[[pid]]$nn_death) ||
            tree_nodes[[pid]]$nn_death == "leaf") {
          tree_nodes[[pid]]$nn_death <- paste0("disappears_nn", nn)
          tree_nodes[[pid]]$topology <- "disappears"
        }
      }
      next
    }

    # For each current-scale block, find the best-matching previous node
    matched_parent <- integer(nrow(current_blocks))
    for (cb in seq_len(nrow(current_blocks))) {
      best_parent <- 0L
      best_ov <- 0
      for (pid in prev_active) {
        nd <- tree_nodes[[pid]]
        ov <- compute_reciprocal_overlap(
          nd$start, nd$end,
          current_blocks$start[cb], current_blocks$end[cb]
        )
        if (ov >= overlap_thresh && ov > best_ov) {
          best_parent <- pid
          best_ov <- ov
        }
      }
      matched_parent[cb] <- best_parent
    }

    # Count how many children each previous node got
    parent_match_count <- tabulate(matched_parent[matched_parent > 0],
                                     nbins = length(tree_nodes))

    new_active <- integer(0)

    # Update each previous node based on its match count
    for (pid in prev_active) {
      n_children <- if (pid <= length(parent_match_count))
        parent_match_count[pid] else 0L

      if (n_children == 0) {
        # Disappears at this scale
        tree_nodes[[pid]]$nn_death <- paste0("disappears_nn", nn)
        tree_nodes[[pid]]$topology <- "disappears"
        # Don't carry it forward

      } else if (n_children == 1) {
        # Stable — update boundaries to the refined match
        child_idx <- which(matched_parent == pid)[1]
        cb_row <- current_blocks[child_idx]
        tree_nodes[[pid]]$start    <- cb_row$start
        tree_nodes[[pid]]$end      <- cb_row$end
        tree_nodes[[pid]]$width    <- cb_row$width
        tree_nodes[[pid]]$height   <- cb_row$height %||% tree_nodes[[pid]]$height
        tree_nodes[[pid]]$start_bp <- cb_row$start_bp %||% tree_nodes[[pid]]$start_bp
        tree_nodes[[pid]]$end_bp   <- cb_row$end_bp   %||% tree_nodes[[pid]]$end_bp
        tree_nodes[[pid]]$start_mb <- cb_row$start_mb %||% tree_nodes[[pid]]$start_mb
        tree_nodes[[pid]]$end_mb   <- cb_row$end_mb   %||% tree_nodes[[pid]]$end_mb
        if (tree_nodes[[pid]]$topology == "root") {
          # Root staying stable — keep topology as "root" for roots,
          # or "stable" for non-roots
        } else {
          tree_nodes[[pid]]$topology <- "stable"
        }
        new_active <- c(new_active, pid)

      } else {
        # Splits — record the death, then create child nodes
        tree_nodes[[pid]]$nn_death <- paste0("splits_nn", nn)
        tree_nodes[[pid]]$topology <- "splits"

        child_indices <- which(matched_parent == pid)
        child_ids <- integer(0)

        for (ci in child_indices) {
          node_id_counter <- node_id_counter + 1L
          tree_nodes[[node_id_counter]] <- make_node(
            id = node_id_counter,
            block_row = current_blocks[ci],
            nn_birth = nn,
            parent = pid,
            topology = "child"
          )
          child_ids <- c(child_ids, node_id_counter)
          new_active <- c(new_active, node_id_counter)
        }
        tree_nodes[[pid]]$children <- paste(child_ids, collapse = ",")
      }
    }

    # Novel blocks at this scale — not matched to any active parent
    unmatched <- which(matched_parent == 0L)
    for (cb in unmatched) {
      cb_row <- current_blocks[cb]
      # Look for a containment parent among all nodes (not just active)
      contain_parent <- NA_integer_
      for (nd_id in seq_along(tree_nodes)) {
        nd <- tree_nodes[[nd_id]]
        if (!is.null(nd) &&
            cb_row$start >= nd$start &&
            cb_row$end   <= nd$end &&
            cb_row$width <  nd$width) {
          if (is.na(contain_parent) ||
              nd$width < tree_nodes[[contain_parent]]$width) {
            contain_parent <- nd_id
          }
        }
      }

      node_id_counter <- node_id_counter + 1L
      tree_nodes[[node_id_counter]] <- make_node(
        id = node_id_counter,
        block_row = cb_row,
        nn_birth = nn,
        parent = contain_parent,
        topology = "novel"
      )
      new_active <- c(new_active, node_id_counter)
    }

    prev_active <- unique(new_active)
  }

  # Anything still active at the finest scale is a "leaf"
  for (pid in prev_active) {
    if (is.na(tree_nodes[[pid]]$nn_death)) {
      tree_nodes[[pid]]$nn_death <- "leaf"
    }
  }

  # Convert to data.table
  as_tree_dt(tree_nodes)
}


make_node <- function(id, block_row, nn_birth, parent, topology) {
  list(
    node_id     = id,
    start       = as.integer(block_row$start),
    end         = as.integer(block_row$end),
    width       = as.integer(block_row$width),
    height      = as.numeric(block_row$height %||% NA_real_),
    start_bp    = as.numeric(block_row$start_bp %||% NA_real_),
    end_bp      = as.numeric(block_row$end_bp   %||% NA_real_),
    start_mb    = as.numeric(block_row$start_mb %||% NA_real_),
    end_mb      = as.numeric(block_row$end_mb   %||% NA_real_),
    nn_birth    = as.integer(nn_birth),
    nn_death    = NA_character_,
    parent_node = as.integer(parent),
    children    = "",
    topology    = as.character(topology)
  )
}


as_tree_dt <- function(tree_nodes) {
  tree_nodes <- tree_nodes[!vapply(tree_nodes, is.null, logical(1))]
  if (length(tree_nodes) == 0) return(data.table())

  rbindlist(lapply(tree_nodes, function(nd) {
    data.table(
      node_id     = nd$node_id,
      start       = nd$start,
      end         = nd$end,
      width       = nd$width,
      start_mb    = nd$start_mb,
      end_mb      = nd$end_mb,
      start_bp    = nd$start_bp,
      end_bp      = nd$end_bp,
      height      = nd$height,
      nn_birth    = nd$nn_birth,
      nn_death    = nd$nn_death,
      parent_node = nd$parent_node,
      children    = nd$children,
      topology    = nd$topology
    )
  }))
}


# ============================================================================
# STEP 4: attach tree fields back to the flat summary
# ============================================================================
#
# For each block in the summary, find the tree node whose interval has
# the highest reciprocal overlap with the block. Take that node's
# nn_birth, nn_death, parent_node, classification. If no node overlaps,
# leave as NA.

attach_tree_fields_to_summary <- function(summary, tree) {
  if (is.null(tree) || nrow(tree) == 0 || nrow(summary) == 0) {
    summary[, nn_birth       := NA_integer_]
    summary[, nn_death       := NA_character_]
    summary[, parent_node    := NA_integer_]
    summary[, classification := NA_character_]
    return(summary)
  }

  nn_birth_per_block       <- integer(nrow(summary))
  nn_death_per_block       <- character(nrow(summary))
  parent_node_per_block    <- integer(nrow(summary))
  classification_per_block <- character(nrow(summary))

  for (r in seq_len(nrow(summary))) {
    bs <- summary$ref_start[r]
    be <- summary$ref_end[r]

    best_node <- NA_integer_
    best_ov <- 0

    for (ti in seq_len(nrow(tree))) {
      ov <- compute_reciprocal_overlap(bs, be, tree$start[ti], tree$end[ti])
      if (ov > best_ov) {
        best_ov <- ov
        best_node <- ti
      }
    }

    if (!is.na(best_node) && best_ov > 0) {
      nn_birth_per_block[r]       <- tree$nn_birth[best_node]
      nn_death_per_block[r]       <- tree$nn_death[best_node]
      parent_node_per_block[r]    <- tree$parent_node[best_node]
      classification_per_block[r] <- tree$classification[best_node]
    } else {
      nn_birth_per_block[r]       <- NA_integer_
      nn_death_per_block[r]       <- NA_character_
      parent_node_per_block[r]    <- NA_integer_
      classification_per_block[r] <- NA_character_
    }
  }

  summary[, nn_birth       := nn_birth_per_block]
  summary[, nn_death       := nn_death_per_block]
  summary[, parent_node    := parent_node_per_block]
  summary[, classification := classification_per_block]

  summary
}


# ============================================================================
# CLASSIFICATION from nn_birth
# ============================================================================
#
# Thresholds are exposed as function args of build_nn_sweep_tree; the
# defaults (200/80/40) match the audit-era convention and what C01d's D3
# score saturates at.

classify_nn_birth <- function(nn_birth,
                               span_mb = NULL,
                               inv_thresh        = 200L,
                               cand_thresh       = 80L,
                               weak_thresh       = 40L,
                               inv_min_span_mb   = 0.5,
                               cand_min_span_mb  = 0.1) {
  # Persistence-based classification
  out <- character(length(nn_birth))
  out[] <- "FAMILY_LD"
  out[nn_birth >= weak_thresh] <- "WEAK_CANDIDATE"
  out[nn_birth >= cand_thresh] <- "CANDIDATE"
  out[nn_birth >= inv_thresh]  <- "INVERSION"

  # v1.4: span-based demotion. A node can survive to coarse NN scales
  # (nn_birth >= 200) and still be too small to be a real inversion —
  # these get demoted one tier. Sub-100-kb "CANDIDATE" blocks are LD
  # fragments, not inversion candidates, so they demote to WEAK.
  if (!is.null(span_mb) && length(span_mb) == length(nn_birth)) {
    too_small_inv  <- out == "INVERSION" & is.finite(span_mb) &
                       span_mb < inv_min_span_mb
    out[too_small_inv] <- "CANDIDATE"

    too_small_cand <- out == "CANDIDATE" & is.finite(span_mb) &
                       span_mb < cand_min_span_mb
    out[too_small_cand] <- "WEAK_CANDIDATE"
  }

  out[is.na(nn_birth)] <- NA_character_
  out
}


# ============================================================================
# RECIPROCAL OVERLAP HELPER
# ============================================================================
#
# min(overlap / width_a, overlap / width_b)
# Returns 0 if no overlap, 1 if identical intervals.

compute_reciprocal_overlap <- function(s1, e1, s2, e2) {
  w1 <- e1 - s1 + 1
  w2 <- e2 - s2 + 1
  if (w1 <= 0 || w2 <= 0) return(0)
  ov_start <- max(s1, s2)
  ov_end   <- min(e1, e2)
  ov <- max(0, ov_end - ov_start + 1)
  min(ov / w1, ov / w2)
}


# Tiny fallback operator in case the caller hasn't defined it
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}


# ============================================================================
# BACKWARDS-COMPAT WRAPPER (matches the old D02 signature)
# ============================================================================
#
# Existing run_all.R code that called compute_nn_persistence_v2() can
# keep working without change. It just gets the fast-mode summary.

compute_nn_persistence_v2 <- function(sim_mats_by_nn,
                                       dt,
                                       staircase_fn = NULL,
                                       overlap_thresh = NULL) {
  result <- build_nn_sweep_tree(
    sim_mats_by_nn  = sim_mats_by_nn,
    dt              = dt,
    staircase_fn    = staircase_fn,
    mode            = "fast",
    overlap_thresh  = overlap_thresh,
    verbose         = FALSE
  )
  result$summary
}


# ============================================================================
# PRINT HELPERS (for interactive use)
# ============================================================================

print_tree_summary <- function(tree_dt, top_n = 20L) {
  if (is.null(tree_dt) || nrow(tree_dt) == 0) {
    cat("Tree is empty.\n")
    return(invisible(NULL))
  }

  cat("\n=== NN Sweep Tree ===\n")
  cat("Total nodes:", nrow(tree_dt), "\n")
  cat("Roots:",       sum(is.na(tree_dt$parent_node)), "\n")
  cat("Leaves:",      sum(tree_dt$nn_death == "leaf", na.rm = TRUE), "\n")

  cat("\nTopology:\n")
  print(table(tree_dt$topology, useNA = "ifany"))

  cat("\nClassification:\n")
  print(table(tree_dt$classification, useNA = "ifany"))

  strong <- tree_dt[classification %in% c("INVERSION", "CANDIDATE")]
  if (nrow(strong) > 0) {
    setorder(strong, -nn_birth)
    cat(sprintf("\nTop %d strong nodes:\n", min(top_n, nrow(strong))))
    print(head(strong[, .(node_id, start, end,
                          start_mb = round(start_mb, 2),
                          end_mb   = round(end_mb, 2),
                          nn_birth, nn_death,
                          parent_node, classification)], top_n))
  }

  invisible(tree_dt)
}


print_summary_overview <- function(summary_dt) {
  if (is.null(summary_dt) || nrow(summary_dt) == 0) {
    cat("Summary is empty.\n")
    return(invisible(NULL))
  }

  cat("\n=== Per-block NN Summary ===\n")
  cat("Blocks:", nrow(summary_dt), "\n")
  cat("\nTopology:\n")
  print(table(summary_dt$nn_topology, useNA = "ifany"))

  if ("classification" %in% names(summary_dt)) {
    cat("\nClassification:\n")
    print(table(summary_dt$classification, useNA = "ifany"))
  }

  survives_cols <- grep("^survives_nn\\d+$", names(summary_dt), value = TRUE)
  if (length(survives_cols) > 0) {
    cat("\nSurvival counts (# blocks surviving at each scale):\n")
    surv_counts <- sapply(survives_cols, function(col) sum(summary_dt[[col]], na.rm = TRUE))
    print(surv_counts)
  }

  invisible(summary_dt)
}


# ============================================================================
# CLI ENTRY POINT
# ============================================================================
#
# Usage:
#   Rscript STEP_D09_nn_sweep_tree.R <chr> <sim_mat_dir> <out_dir>
#       [--mode deep|fast] [--staircase-src /path/to/STEP_D01_staircase_boundaries.R]
#
#   <sim_mat_dir> should contain files named <chr>.sim_mat_nn<k>.rds
#   (or <chr>_sim_mat_nn<k>.rds — either pattern is accepted).

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript STEP_D09_nn_sweep_tree.R <chr> <sim_mat_dir> <out_dir> ",
         "[--mode deep|fast] [--staircase-src <path>]")
  }

  chr          <- args[1]
  sim_mat_dir  <- args[2]
  out_dir      <- args[3]
  cli_mode     <- "deep"
  staircase_src <- NULL

  i <- 4L
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--mode" && i < length(args)) {
      cli_mode <- args[i + 1]; i <- i + 2L
    } else if (a == "--staircase-src" && i < length(args)) {
      staircase_src <- args[i + 1]; i <- i + 2L
    } else {
      i <- i + 1L
    }
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Source the staircase
  if (!is.null(staircase_src)) {
    source(staircase_src)
  } else if (!exists("detect_blocks_staircase", mode = "function")) {
    # Try same dir as this script
    this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
    candidates <- list.files(this_dir, pattern = "STEP_D01.*staircase.*\\.R$",
                              full.names = TRUE)
    if (length(candidates) > 0) {
      message("[D09] sourcing staircase from ", candidates[1])
      source(candidates[1])
    } else {
      stop("detect_blocks_staircase() not found. Pass --staircase-src /path.")
    }
  }

  # Load sim_mats for this chromosome
  patterns <- c(
    paste0("^", chr, "\\.sim_mat_nn(\\d+)\\.rds$"),
    paste0("^", chr, "_sim_mat_nn(\\d+)\\.rds$"),
    paste0("^sim_mat_", chr, "_nn(\\d+)\\.rds$")
  )
  all_files <- list.files(sim_mat_dir, full.names = FALSE)

  sim_mats <- list()
  for (pat in patterns) {
    hits <- grep(pat, all_files, value = TRUE)
    for (f in hits) {
      nn <- as.integer(sub(pat, "\\1", f))
      if (!as.character(nn) %in% names(sim_mats)) {
        sim_mats[[as.character(nn)]] <- readRDS(file.path(sim_mat_dir, f))
      }
    }
  }

  # Also try to load nn0 from a precomp RDS if present and no nn0 found yet
  if (!("0" %in% names(sim_mats))) {
    precomp_rds <- file.path(sim_mat_dir, paste0(chr, ".precomp.rds"))
    if (file.exists(precomp_rds)) {
      pc <- readRDS(precomp_rds)
      if (!is.null(pc$sim_mat)) sim_mats[["0"]] <- pc$sim_mat
    }
  }

  if (length(sim_mats) == 0) {
    stop("[D09] no sim_mats found for chrom '", chr, "' in ", sim_mat_dir)
  }

  cat("[D09] loaded", length(sim_mats), "scales for", chr, ":",
      paste(sort(as.integer(names(sim_mats))), collapse = ", "), "\n")

  result <- build_nn_sweep_tree(sim_mats, mode = cli_mode)

  fwrite(result$summary, file.path(out_dir, paste0("nn_summary_", chr, ".tsv")),
         sep = "\t")
  if (!is.null(result$tree)) {
    fwrite(result$tree, file.path(out_dir, paste0("nn_tree_", chr, ".tsv")),
           sep = "\t")
  }

  print_summary_overview(result$summary)
  if (!is.null(result$tree)) print_tree_summary(result$tree)

  cat("\n[D09] done. Mode =", result$mode,
      "| scales:", paste(result$nn_scales, collapse = ", "), "\n")
}
