#!/usr/bin/env Rscript
# ============================================================================
# STEP_D13_plot_annotated_simmat.R — v2 (sandbox debug edition)
# ============================================================================
#
# CHROMOSOME-WIDE ANNOTATED SIMILARITY MATRIX with full debug overlay stack.
#
# Backwards compatible with v1: when only --precomp and --blocks are passed,
# behavior is identical (minus two v1 bugs fixed below). New --blocks_01c,
# --tree, --sv_prior, --scales_summary flags enable additional debug layers
# that compare detectors and validate block edges against known breakpoints.
#
# LAYERS (enabled by the presence of the corresponding input)
# -----------------------------------------------------------
#   1. Sim_mat heatmap (base)                                     --precomp
#   2. Mini-block density zones (orange semi-transparent)         --blocks
#   3. D01 block outlines, colored by inversion-system root       --blocks
#   4. 01C block outlines, dashed teal (comparator detector)      --blocks_01c
#   5. D09 tree classification ring (overrides system colors)     --tree
#   6. Per-block labels with Mb range, type, classification       --blocks
#   7. Internal-feature X marks                                   --steps
#   8. SV breakpoint ticks (DELLY/Manta INVs)                     --sv_prior
#   9. NN persistence strip (right margin, per-block)             --scales_summary
#  10. Vote profile as separate track plot                        --votes
#
# INPUT SCHEMA NOTES
# ------------------
#   --blocks            D01 block table; needs block_id, start_mb, end_mb,
#                       parent_id, width, height; optional: shape_class, depth
#   --blocks_01c        01C block_registry TSV; needs start_bp, end_bp OR
#                       start_mb, end_mb, and a block label column
#   --tree              D09 nn_tree TSV; needs start_mb, end_mb, classification
#                       (INVERSION / CANDIDATE / WEAK_CANDIDATE / FAMILY_LD)
#   --sv_prior          sv_prior_<chr>.rds from STEP_C00; reads $inv_calls with
#                       bp1, bp2, svlen, caller
#   --scales_summary    D09 nn_summary TSV; needs block_id, ref_start, ref_end,
#                       survives_nn{20,40,80,120,160,200,240,320} logicals
#
# USAGE
# -----
#   Rscript STEP_D13_plot_annotated_simmat.R \
#     --precomp         precomp/C_gar_LG28.precomp.rds \
#     --blocks          d01/blocks_C_gar_LG28.tsv \
#     --votes           d01/votes_C_gar_LG28.tsv.gz \
#     --steps           d01/steps_C_gar_LG28.tsv.gz \
#     --blocks_01c      landscape/block_registry_C_gar_LG28.tsv.gz \
#     --tree            d09/nn_tree_C_gar_LG28.tsv \
#     --scales_summary  d09/nn_summary_C_gar_LG28.tsv \
#     --sv_prior        sv_prior/sv_prior_C_gar_LG28.rds \
#     --outdir          plots/ \
#     --chr             C_gar_LG28
#
# FIXES vs v1
# -----------
#   - DEPTH_COLORS typo in vote-profile section replaced with SYSTEM_PALETTE
#     (v1 crashed whenever the vote-profile track tried to render)
#   - sys.frame(1)$ofile guarded with tryCatch for non-sourced contexts
#   - classification override logic routes through a single compute_outline_color()
#     so D01 system palette and D09 tree classification don't fight
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a


# ============================================================================
# CLI
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
precomp_file    <- NULL
blocks_file     <- NULL
votes_file      <- NULL
steps_file      <- NULL
blocks_01c_file <- NULL
tree_file       <- NULL
sv_prior_file   <- NULL
scales_summary_file <- NULL
outdir          <- "plots"
chr_label       <- NULL
max_display     <- 600L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  take_next <- function() { v <- args[i + 1L]; i <<- i + 2L; v }
  if (a == "--precomp"         && i < length(args)) { precomp_file    <- take_next() }
  else if (a == "--blocks"     && i < length(args)) { blocks_file     <- take_next() }
  else if (a == "--votes"      && i < length(args)) { votes_file      <- take_next() }
  else if (a == "--steps"      && i < length(args)) { steps_file      <- take_next() }
  else if (a == "--blocks_01c" && i < length(args)) { blocks_01c_file <- take_next() }
  else if (a == "--tree"       && i < length(args)) { tree_file       <- take_next() }
  else if (a == "--sv_prior"   && i < length(args)) { sv_prior_file   <- take_next() }
  else if (a == "--scales_summary" && i < length(args)) { scales_summary_file <- take_next() }
  else if (a == "--outdir"     && i < length(args)) { outdir          <- take_next() }
  else if (a == "--chr"        && i < length(args)) { chr_label       <- take_next() }
  else if (a == "--max_display"&& i < length(args)) { max_display     <- as.integer(take_next()) }
  else { i <- i + 1L }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# STYLE CONFIG
# ============================================================================

# System palette — one color per ROOT block (independent inversion system).
# All nested children inherit the root's color.
SYSTEM_PALETTE <- c(
  "#1a1a1a",  # 1: near-black
  "#1e40af",  # 2: dark blue
  "#991b1b",  # 3: dark red
  "#166534",  # 4: dark green
  "#7c2d12",  # 5: dark brown
  "#6b21a8",  # 6: purple
  "#0e7490",  # 7: dark cyan
  "#b45309",  # 8: amber
  "#4338ca",  # 9: indigo
  "#be123c",  # 10: rose
  "#0f766e",  # 11: teal
  "#a16207",  # 12: yellow-brown
  "#7e22ce",  # 13: violet
  "#0369a1",  # 14: sky blue
  "#9f1239"   # 15: deep pink
)

# D09 tree-classification palette — overrides system colors when --tree given
TREE_CLASS_COLORS <- c(
  INVERSION       = "#dc2626",  # red — strong persistence (nn_birth >= 200)
  CANDIDATE       = "#ea580c",  # orange — moderate (80-199)
  WEAK_CANDIDATE  = "#facc15",  # yellow — weak (40-79)
  FAMILY_LD       = "#6b7280"   # grey — noise / family LD (<40)
)

# 01C comparator outline — single distinctive style so it reads as "the other
# detector said this too"
STYLE_01C <- list(
  color    = "#0d9488",   # teal
  linetype = "longdash",
  lwd      = 0.7
)

# SV breakpoint styling
STYLE_SV <- list(
  delly_color = "#1e40af",   # blue for DELLY
  manta_color = "#ea580c",   # orange for Manta
  both_color  = "#7c2d12",   # dark brown when both callers agree
  alpha       = 0.55,
  lwd         = 0.45
)

# Gap / non-block regions (kept from v1 for legend / future use)
GAP_COLORS <- c(
  family_ld_band       = "#9ca3af",
  diffuse_diagonal     = "#78716c",
  extended_suppression = "#0d9488",
  clean_background     = "#d4d4d8",
  transition_zone      = "#a3a3a3",
  mini_block_zone      = "#ea580c"
)

# Line width / linetype / label size by nesting depth
DEPTH_LWD       <- c(1.3, 0.95, 0.7, 0.55, 0.4, 0.35, 0.3)
DEPTH_LINETYPE  <- c("solid", "solid", "dashed", "dashed", "dotted", "dotted", "dotted")
DEPTH_LABEL_SIZE<- c(2.8, 2.3, 2.0, 1.7, 1.5, 1.3, 1.2)

# Sim_mat color scale
SIM_COLORS <- c("#0C1E3C", "#1E3A5F", "#4A7FB5", "#C4A035", "#D4712A", "#B8282E")


# ============================================================================
# LOAD DATA
# ============================================================================

load_precomp <- function(f) {
  stopifnot(!is.null(f), file.exists(f))
  pc <- readRDS(f)
  list(smat = pc$sim_mat, dt = pc$dt,
       chrom = pc$chrom, n_windows = nrow(pc$sim_mat))
}

safe_fread <- function(f) {
  if (is.null(f) || !file.exists(f)) return(data.table())
  tryCatch(fread(f), error = function(e) {
    message("[PLOT] failed to read ", f, ": ", conditionMessage(e))
    data.table()
  })
}

safe_readRDS <- function(f) {
  if (is.null(f) || !file.exists(f)) return(NULL)
  tryCatch(readRDS(f), error = function(e) {
    message("[PLOT] failed to read ", f, ": ", conditionMessage(e))
    NULL
  })
}

pc        <- load_precomp(precomp_file)
smat      <- pc$smat
dt        <- pc$dt
chr_label <- chr_label %||% pc$chrom
N         <- pc$n_windows

blocks          <- safe_fread(blocks_file)
votes           <- safe_fread(votes_file)
steps           <- safe_fread(steps_file)
blocks_01c      <- safe_fread(blocks_01c_file)
tree            <- safe_fread(tree_file)
scales_summary  <- safe_fread(scales_summary_file)
sv_prior        <- safe_readRDS(sv_prior_file)

message("[PLOT] ", chr_label, " | ", N, " windows")
message("  D01 blocks:       ", nrow(blocks))
message("  01C blocks:       ", nrow(blocks_01c))
message("  D09 tree nodes:   ", nrow(tree))
message("  NN summary rows:  ", nrow(scales_summary))
message("  SV prior:         ",
        if (is.null(sv_prior)) "none" else
        paste0(nrow(sv_prior$inv_calls %||% data.table()), " INV calls"))

# Optional auto-load of landscape classification if --blocks dir has it
landscape <- data.table()
if (!is.null(blocks_file)) {
  landscape_file <- sub("blocks_", "landscape_", blocks_file)
  if (file.exists(landscape_file)) {
    landscape <- fread(landscape_file)
    message("  Landscape classification: ", nrow(landscape), " entries")
  }
}


# ============================================================================
# NORMALISE INCOMING BLOCK TABLES TO A COMMON COLUMN SET
# ============================================================================

normalise_blocks <- function(bt, window_size_bp = 50000L) {
  if (nrow(bt) == 0) return(bt)

  # Ensure Mb columns exist
  if (!"start_mb" %in% names(bt)) {
    if ("start_bp" %in% names(bt)) {
      bt[, start_mb := start_bp / 1e6]
    } else if ("start" %in% names(bt)) {
      bt[, start_mb := (start - 1L) * window_size_bp / 1e6]
    }
  }
  if (!"end_mb" %in% names(bt)) {
    if ("end_bp" %in% names(bt)) {
      bt[, end_mb := end_bp / 1e6]
    } else if ("end" %in% names(bt)) {
      bt[, end_mb := (end - 1L) * window_size_bp / 1e6]
    }
  }

  # parent_id
  if (!"parent_id" %in% names(bt)) bt[, parent_id := NA_integer_]

  # block_id
  if (!"block_id" %in% names(bt)) {
    bt[, block_id := seq_len(.N)]
  }

  # width
  if (!"width" %in% names(bt)) {
    if ("start" %in% names(bt) && "end" %in% names(bt)) {
      bt[, width := end - start + 1L]
    } else {
      bt[, width := as.integer(round((end_mb - start_mb) * 1e6 / window_size_bp))]
    }
  }

  # height
  if (!"height" %in% names(bt)) bt[, height := NA_real_]

  bt
}

blocks     <- normalise_blocks(blocks)
blocks_01c <- normalise_blocks(blocks_01c)


# ============================================================================
# COMPUTE NESTING DEPTH (v1 logic, lightly cleaned)
# ============================================================================

compute_depths <- function(bt) {
  if (nrow(bt) == 0) return(bt)
  bt[, depth := 0L]
  for (r in seq_len(nrow(bt))) {
    pid <- bt$parent_id[r]
    d <- 0L
    while (!is.na(pid) && d < 10) {
      d <- d + 1L
      parent_row <- which(bt$block_id == pid)
      if (length(parent_row) == 0) break
      pid <- bt$parent_id[parent_row[1]]
    }
    bt$depth[r] <- d
  }
  bt
}

if (nrow(blocks) > 0) blocks <- compute_depths(blocks)


# ============================================================================
# BUILD HEATMAP DATA (subsampled)
# ============================================================================

step_size <- max(1L, N %/% max_display)
idx <- seq(1L, N, by = step_size)
ns <- length(idx)

pos_mb <- if ("start_bp" %in% names(dt)) {
  (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6
} else {
  (idx - 1L) * 50000 / 1e6
}

message("[PLOT] Building heatmap data (", ns, "x", ns, ")...")

# Vectorised: build a full grid, then look up values
grid <- CJ(ii = seq_len(ns), jj = seq_len(ns))
grid[, x := pos_mb[ii]]
grid[, y := pos_mb[jj]]
grid[, value := smat[cbind(idx[ii], idx[jj])]]
grid <- grid[is.finite(value), .(x, y, value)]

hm_dt <- grid
sim_med <- median(hm_dt$value, na.rm = TRUE)

sim_scale <- scale_fill_gradientn(
  colours = SIM_COLORS,
  values  = scales::rescale(c(0, sim_med * 0.5, sim_med * 0.85,
                               sim_med * 1.1, sim_med * 1.3, 1)),
  name    = "Similarity",
  limits  = c(0, 1)
)


# ============================================================================
# SYSTEM COLORS — root-based (v1 logic)
# ============================================================================

assign_system_colors <- function(bt) {
  if (nrow(bt) == 0) return(character(0))

  block_root <- integer(nrow(bt))
  for (r in seq_len(nrow(bt))) {
    bid <- bt$block_id[r]
    pid <- bt$parent_id[r]
    root <- bid
    k <- 0L
    while (!is.na(pid) && k < 20) {
      root <- pid
      pr <- which(bt$block_id == pid)
      if (length(pr) == 0) break
      pid <- bt$parent_id[pr[1]]
      k <- k + 1L
    }
    block_root[r] <- root
  }

  unique_roots <- unique(block_root)
  root_color_map <- setNames(
    SYSTEM_PALETTE[(seq_along(unique_roots) - 1L) %% length(SYSTEM_PALETTE) + 1L],
    as.character(unique_roots)
  )

  colors <- root_color_map[as.character(block_root)]
  names(colors) <- as.character(bt$block_id)
  colors
}


# ============================================================================
# TREE CLASSIFICATION LOOKUP — attach classification to each D01 block
# by reciprocal-overlap matching to tree nodes
# ============================================================================

attach_tree_classification <- function(bt, tr) {
  if (nrow(bt) == 0 || nrow(tr) == 0) return(character(0))
  if (!"classification" %in% names(tr)) return(character(0))

  # Tree may have start/end in bins or Mb — use whichever we have
  tree_start_mb <- if ("start_mb" %in% names(tr)) tr$start_mb else
    if ("start" %in% names(tr)) (tr$start - 1L) * 50000 / 1e6 else numeric(0)
  tree_end_mb <- if ("end_mb" %in% names(tr)) tr$end_mb else
    if ("end" %in% names(tr)) (tr$end - 1L) * 50000 / 1e6 else numeric(0)

  if (length(tree_start_mb) == 0 || length(tree_end_mb) == 0) return(character(0))

  out <- character(nrow(bt))
  out[] <- NA_character_

  for (r in seq_len(nrow(bt))) {
    bs <- bt$start_mb[r]; be <- bt$end_mb[r]
    bw <- be - bs
    if (bw <= 0) next

    best_ov <- 0
    best_cls <- NA_character_
    for (tr_i in seq_len(nrow(tr))) {
      ts <- tree_start_mb[tr_i]; te <- tree_end_mb[tr_i]
      tw <- te - ts
      if (tw <= 0) next
      ov <- max(0, min(be, te) - max(bs, ts))
      rec <- min(ov / bw, ov / tw)
      if (rec > best_ov) {
        best_ov  <- rec
        best_cls <- tr$classification[tr_i]
      }
    }
    if (best_ov > 0.3) out[r] <- best_cls
  }

  names(out) <- as.character(bt$block_id)
  out
}


# ============================================================================
# OUTLINE COLOR RESOLVER
# Tree classification takes priority; falls back to system-root color.
# ============================================================================

compute_outline_color <- function(sys_color, tree_class) {
  if (!is.na(tree_class) && tree_class %in% names(TREE_CLASS_COLORS)) {
    TREE_CLASS_COLORS[[tree_class]]
  } else {
    sys_color
  }
}


# ============================================================================
# NN SURVIVAL STRIP — per-block badge showing which scales it survives at
# ============================================================================

build_nn_strip_data <- function(bt, scales_summary) {
  if (nrow(bt) == 0 || nrow(scales_summary) == 0) return(data.table())

  survives_cols <- grep("^survives_nn\\d+$", names(scales_summary), value = TRUE)
  if (length(survives_cols) == 0) return(data.table())

  nn_scales <- as.integer(sub("survives_nn", "", survives_cols))

  # Join on block_id
  key_col <- if ("block_id" %in% names(scales_summary)) "block_id" else
             if ("candidate_id" %in% names(scales_summary)) "candidate_id" else NA_character_
  if (is.na(key_col)) return(data.table())

  merged <- merge(
    bt[, .(block_id, end_mb)],
    scales_summary[, c(key_col, survives_cols), with = FALSE],
    by.x = "block_id", by.y = key_col, all.x = TRUE
  )

  # Long format for plotting
  long <- list()
  for (ci in seq_along(survives_cols)) {
    col <- survives_cols[ci]
    nn  <- nn_scales[ci]
    long[[ci]] <- data.table(
      block_id = merged$block_id,
      y        = merged$end_mb,
      scale    = nn,
      survives = as.logical(merged[[col]])
    )
  }
  rbindlist(long)
}


# ============================================================================
# BUILD BLOCK OUTLINE RECTANGLES (D01) — extended for tree classification
# ============================================================================

build_block_outlines <- function(bt, landscape = NULL, tree = NULL) {
  if (nrow(bt) == 0) return(data.table())

  sys_colors <- assign_system_colors(bt)

  tree_class <- if (!is.null(tree) && nrow(tree) > 0) {
    attach_tree_classification(bt, tree)
  } else setNames(rep(NA_character_, nrow(bt)), as.character(bt$block_id))

  outlines <- list()
  for (r in seq_len(nrow(bt))) {
    b <- bt[r]
    d <- min((b$depth %||% 0L) + 1L, length(DEPTH_LWD))

    x_lo <- b$start_mb
    x_hi <- b$end_mb

    bid_str <- as.character(b$block_id)
    sys_col <- if (bid_str %in% names(sys_colors)) sys_colors[bid_str] else "#1a1a1a"
    t_cls   <- if (bid_str %in% names(tree_class)) tree_class[bid_str] else NA_character_
    color   <- compute_outline_color(sys_col, t_cls)

    # Rich label
    rich_label <- NULL
    confidence <- NA_real_
    cat_name <- "unclassified"

    if (!is.null(landscape) && nrow(landscape) > 0) {
      ls_row <- landscape[block_id == b$block_id]
      if (nrow(ls_row) > 0) {
        cat_name   <- ls_row$category[1]
        rich_label <- ls_row$label[1]
        confidence <- ls_row$confidence_pct[1]
      }
    }
    if (is.null(rich_label) || is.na(rich_label)) rich_label <- build_block_label(b, t_cls)

    outlines[[r]] <- data.table(
      block_id       = b$block_id,
      xmin           = x_lo,
      xmax           = x_hi,
      ymin           = x_lo,
      ymax           = x_hi,
      depth          = b$depth %||% 0L,
      category       = cat_name,
      confidence     = confidence,
      tree_class     = t_cls,
      color          = as.character(color),
      system_color   = as.character(sys_col),
      lwd            = DEPTH_LWD[d],
      linetype       = DEPTH_LINETYPE[d],
      label_size     = DEPTH_LABEL_SIZE[d],
      label          = rich_label,
      label_x        = x_lo + (x_hi - x_lo) * 0.02,
      label_y        = x_hi - (x_hi - x_lo) * 0.05
    )
  }
  rbindlist(outlines)
}


# ============================================================================
# 01C BLOCK OUTLINES — single-style dashed teal rectangles
# ============================================================================

build_01c_outlines <- function(bt) {
  if (nrow(bt) == 0) return(data.table())

  rows <- list()
  for (r in seq_len(nrow(bt))) {
    b <- bt[r]
    x_lo <- b$start_mb
    x_hi <- b$end_mb

    block_name <- if ("block" %in% names(b)) as.character(b$block) else
                   if ("label" %in% names(b)) as.character(b$label) else
                   paste0("C", b$block_id)

    rows[[r]] <- data.table(
      xmin       = x_lo,
      xmax       = x_hi,
      ymin       = x_lo,
      ymax       = x_hi,
      color      = STYLE_01C$color,
      linetype   = STYLE_01C$linetype,
      lwd        = STYLE_01C$lwd,
      label      = block_name,
      label_x    = x_hi - (x_hi - x_lo) * 0.02,
      label_y    = x_lo + (x_hi - x_lo) * 0.05
    )
  }
  rbindlist(rows)
}


# ============================================================================
# MINI-BLOCK DENSITY ZONES (v1 logic preserved)
# ============================================================================

build_density_zones <- function(bt, window_size_bp = 50000L,
                                 zone_window_mb = 2.0,
                                 min_blocks_per_zone = 3L) {
  if (nrow(bt) == 0) return(data.table())

  small <- bt[width < 50]
  if (nrow(small) < min_blocks_per_zone) return(data.table())

  max_mb <- max(small$end_mb, na.rm = TRUE)
  zones <- list()

  zone_start <- 0
  while (zone_start < max_mb) {
    zone_end <- zone_start + zone_window_mb
    zone_blocks <- small[start_mb >= zone_start & end_mb <= zone_end]
    n_in_zone <- nrow(zone_blocks)

    if (n_in_zone >= min_blocks_per_zone) {
      zones[[length(zones) + 1L]] <- data.table(
        zone_start_mb = zone_start,
        zone_end_mb   = zone_end,
        n_blocks      = n_in_zone,
        density       = round(n_in_zone / zone_window_mb, 1),
        mean_width    = round(mean(zone_blocks$width), 1),
        label         = paste0("Dense: ", n_in_zone, " mini-blocks / ",
                               zone_window_mb, " Mb")
      )
    }
    zone_start <- zone_start + zone_window_mb * 0.5
  }

  if (length(zones) > 0) rbindlist(zones) else data.table()
}


build_block_label <- function(b, tree_class = NA_character_) {
  type_short <- if ("shape_class" %in% names(b) && !is.na(b$shape_class)) {
    switch(as.character(b$shape_class),
      strong_square  = "Strong",
      diffuse_square = "Diffuse",
      diagonal_band  = "Band",
      noise          = "Noise",
      "")
  } else ""

  parent_info <- if (!is.na(b$parent_id)) paste0("child of B", b$parent_id) else "outer"
  span <- round(b$end_mb - b$start_mb, 1)

  cls_line <- if (!is.na(tree_class)) paste0("\n[", tree_class, "]") else ""

  paste0("B", b$block_id, " [", round(b$start_mb, 1), "-",
         round(b$end_mb, 1), " Mb] ", span, "Mb ", type_short, "\n",
         parent_info,
         if (!is.na(b$height)) paste0(" h=", round(b$height, 3)) else "",
         cls_line)
}


# ============================================================================
# SV PRIOR — extract INV breakpoint Mb positions
# ============================================================================

build_sv_breakpoints <- function(sv_prior, chr_label) {
  if (is.null(sv_prior)) return(data.table())
  if (!is.list(sv_prior) || is.null(sv_prior$inv_calls)) return(data.table())
  inv_calls <- as.data.table(sv_prior$inv_calls)
  if (nrow(inv_calls) == 0) return(data.table())

  # Filter to this chromosome if a chrom column is present
  if ("chrom" %in% names(inv_calls)) {
    inv_calls <- inv_calls[chrom == chr_label]
    if (nrow(inv_calls) == 0) return(data.table())
  }

  caller_col <- if ("caller" %in% names(inv_calls)) inv_calls$caller else
                  rep("unknown", nrow(inv_calls))

  # Both breakpoints become ticks
  rows <- list()
  for (i in seq_len(nrow(inv_calls))) {
    rows[[length(rows) + 1L]] <- data.table(
      inv_id = inv_calls$inv_id[i] %||% paste0("INV", i),
      bp_mb  = inv_calls$bp1[i] / 1e6,
      side   = "bp1",
      caller = caller_col[i]
    )
    rows[[length(rows) + 1L]] <- data.table(
      inv_id = inv_calls$inv_id[i] %||% paste0("INV", i),
      bp_mb  = inv_calls$bp2[i] / 1e6,
      side   = "bp2",
      caller = caller_col[i]
    )
  }
  bps <- rbindlist(rows)
  bps[, tick_color := fifelse(caller == "delly", STYLE_SV$delly_color,
                       fifelse(caller == "manta", STYLE_SV$manta_color,
                               STYLE_SV$both_color))]
  bps
}


# ============================================================================
# INTERNAL FEATURE MARKS (v1 logic with window-size awareness)
# ============================================================================

build_internal_marks <- function(steps, bt, window_size_bp = 50000L) {
  if (nrow(steps) == 0 || nrow(bt) == 0) return(data.table())
  if (!"inside_block" %in% names(steps)) return(data.table())

  internal <- steps[inside_block == TRUE]
  if (nrow(internal) == 0) return(data.table())

  unique_pos <- sort(unique(internal$step_position))
  data.table(
    position_mb = (unique_pos - 1L) * window_size_bp / 1e6,
    type        = "internal_feature"
  )
}


# ============================================================================
# VOTE TRACK (v1 logic, unchanged)
# ============================================================================

build_vote_track <- function(votes, ns) {
  if (nrow(votes) == 0) return(data.table())
  votes_sub <- votes[seq(1L, nrow(votes), length.out = min(ns, nrow(votes))), ]
  votes_sub[, pos_mb := (position - 1L) * 50000 / 1e6]
  votes_sub
}


# ============================================================================
# ASSEMBLE THE MAIN PLOT
# ============================================================================

outline_dt     <- build_block_outlines(blocks, landscape, tree)
outline_01c_dt <- build_01c_outlines(blocks_01c)
density_zones  <- build_density_zones(blocks)
marks_dt       <- build_internal_marks(steps, blocks)
sv_bp_dt       <- build_sv_breakpoints(sv_prior, chr_label)
nn_strip_dt    <- build_nn_strip_data(blocks, scales_summary)
vote_track     <- build_vote_track(votes, ns)

if (nrow(density_zones) > 0) {
  message("[PLOT] Mini-block density zones: ", nrow(density_zones))
}

message("[PLOT] Assembling annotated heatmap...")

p <- ggplot()

# Layer 1: sim_mat heatmap
p <- p + geom_raster(data = hm_dt, aes(x = x, y = y, fill = value),
                      interpolate = TRUE) +
  sim_scale

# Layer 1.5: mini-block density zones
if (nrow(density_zones) > 0) {
  for (dz in seq_len(nrow(density_zones))) {
    dzr <- density_zones[dz]
    p <- p + annotate(
      "rect",
      xmin = dzr$zone_start_mb, xmax = dzr$zone_end_mb,
      ymin = dzr$zone_start_mb, ymax = dzr$zone_end_mb,
      fill = "#ea580c", alpha = 0.06, color = "#ea580c",
      linewidth = 0.35, linetype = "dotted"
    ) + annotate(
      "text",
      x = dzr$zone_start_mb, y = dzr$zone_end_mb,
      label = dzr$label,
      size = 1.8, color = "#ea580c",
      hjust = 0, vjust = 1, fontface = "italic"
    )
  }
}

# Layer 2: 01C block outlines (teal dashed) — drawn BEFORE D01 so D01 sits on top
if (nrow(outline_01c_dt) > 0) {
  for (r in seq_len(nrow(outline_01c_dt))) {
    od <- outline_01c_dt[r]
    p <- p + geom_rect(
      data = od,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA,
      color = od$color,
      linewidth = od$lwd,
      linetype = od$linetype
    ) + annotate(
      "text", x = od$label_x, y = od$label_y,
      label = paste0("(01C) ", od$label),
      size = 1.9, color = STYLE_01C$color,
      hjust = 1, vjust = 0, fontface = "italic", alpha = 0.8
    )
  }
}

# Layer 3: D01 block outlines — colored by tree class (if --tree) else system root
if (nrow(outline_dt) > 0) {
  for (r in seq_len(nrow(outline_dt))) {
    od <- outline_dt[r]
    p <- p + geom_rect(
      data = od,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA,
      color = od$color,
      linewidth = od$lwd,
      linetype = od$linetype
    )
  }

  # Layer 4: D01 block labels
  for (r in seq_len(nrow(outline_dt))) {
    od <- outline_dt[r]
    p <- p + annotate(
      "text",
      x = od$label_x, y = od$label_y,
      label = od$label,
      size = od$label_size,
      color = od$color,
      hjust = 0, vjust = 1,
      fontface = "bold",
      alpha = 0.9
    )
  }
}

# Layer 5: internal-feature X marks
if (nrow(marks_dt) > 0) {
  p <- p + geom_point(
    data = marks_dt,
    aes(x = position_mb, y = position_mb),
    shape = 4, size = 0.8,
    color = "#dc2626",
    alpha = 0.6
  )
}

# Layer 6: SV breakpoint ticks at the top margin.
# Short vertical segments just above the diagonal maximum.
if (nrow(sv_bp_dt) > 0) {
  y_max <- max(hm_dt$y, na.rm = TRUE)
  tick_top    <- y_max + 0.015 * (y_max - min(hm_dt$y, na.rm = TRUE))
  tick_bottom <- y_max - 0.005 * (y_max - min(hm_dt$y, na.rm = TRUE))
  for (s in seq_len(nrow(sv_bp_dt))) {
    bs <- sv_bp_dt[s]
    p <- p + annotate(
      "segment",
      x = bs$bp_mb, xend = bs$bp_mb,
      y = tick_bottom, yend = tick_top,
      color = bs$tick_color,
      linewidth = STYLE_SV$lwd,
      alpha = STYLE_SV$alpha
    )
  }
}

# Layer 7: NN survival strip at the right margin.
# Horizontal band of cells to the right of the heatmap, one per scale.
if (nrow(nn_strip_dt) > 0) {
  scales_present <- sort(unique(nn_strip_dt$scale))
  x_base <- max(hm_dt$x, na.rm = TRUE)
  x_span <- x_base - min(hm_dt$x, na.rm = TRUE)
  cell_w <- 0.012 * x_span

  for (bi in seq_len(nrow(blocks))) {
    b <- blocks[bi]
    # Centre the strip vertically on the block's mid-line
    y_cell <- (b$start_mb + b$end_mb) / 2
    for (si in seq_along(scales_present)) {
      nn <- scales_present[si]
      row <- nn_strip_dt[block_id == b$block_id & scale == nn]
      if (nrow(row) == 0) next
      survives <- isTRUE(row$survives[1])
      cell_x <- x_base + (si - 0.5) * cell_w
      cell_fill <- if (survives) "#16a34a" else "#e5e7eb"
      p <- p + annotate(
        "tile",
        x = cell_x, y = y_cell,
        width = cell_w * 0.9,
        height = max(0.4, (b$end_mb - b$start_mb) * 0.9),
        fill = cell_fill, alpha = 0.8
      )
      # Scale label, only at the top of the plot
      if (bi == 1) {
        p <- p + annotate(
          "text",
          x = cell_x,
          y = max(hm_dt$y, na.rm = TRUE),
          label = paste0("nn", nn),
          size = 1.7, color = "grey30",
          angle = 90, hjust = -0.1, vjust = 0.5
        )
      }
    }
  }
}

# Cosmetics
caption_lines <- c(
  "D01 blocks: same color = same inversion system (root + all children)",
  "Line style: thick solid = outer; thin dashed = inner nested; dotted = deep"
)
if (nrow(outline_01c_dt) > 0) caption_lines <- c(caption_lines,
  "01C comparator: teal long-dash outlines")
if (!is.null(tree) && nrow(tree) > 0) caption_lines <- c(caption_lines,
  "Tree class: red = INVERSION; orange = CANDIDATE; yellow = WEAK; grey = FAMILY_LD")
if (nrow(sv_bp_dt) > 0) caption_lines <- c(caption_lines,
  "SV ticks (top): blue = DELLY; orange = Manta; brown = both")
if (nrow(nn_strip_dt) > 0) caption_lines <- c(caption_lines,
  "Right strip: per-block NN survival (green = survives at scale)")
caption_lines <- c(caption_lines,
  "Red x = internal feature step inside block")

p <- p +
  coord_fixed() +
  labs(
    x = paste0(chr_label, " (Mb)"),
    y = paste0(chr_label, " (Mb)"),
    title = paste0(chr_label, " -- Annotated Similarity Matrix (D13 v2)"),
    subtitle = paste0(
      N, " windows (", ns, "x", ns, " display) | ",
      nrow(blocks), " D01 blocks | ",
      if (nrow(blocks_01c) > 0) paste0(nrow(blocks_01c), " 01C blocks | ") else "",
      if (!is.null(tree) && nrow(tree) > 0) paste0(nrow(tree), " tree nodes | ") else "",
      "median sim = ", round(sim_med, 3)
    ),
    caption = paste(caption_lines, collapse = "\n")
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    plot.caption     = element_text(size = 6, hjust = 0, color = "#8A8A8A"),
    legend.key.size  = unit(0.4, "cm")
  )


# ============================================================================
# SAVE MAIN PLOT
# ============================================================================

out_base <- paste0(chr_label, "_annotated_simmat_v2")

for (ext in c("png", "pdf")) {
  out_file <- file.path(outdir, paste0(out_base, ".", ext))
  tryCatch({
    ggsave(out_file, p, width = 14, height = 13, dpi = 350)
    message("[SAVED] ", out_file)
  }, error = function(e) message("[FAIL] ", ext, ": ", conditionMessage(e)))
}


# ============================================================================
# VOTE PROFILE (separate plot) — BUGFIX: DEPTH_COLORS was undefined in v1;
# use system palette or outline color derived from blocks
# ============================================================================

if (nrow(vote_track) > 0 && "vote_smooth" %in% names(vote_track)) {
  p_vote <- ggplot(vote_track, aes(x = pos_mb, y = vote_smooth)) +
    geom_area(fill = "#4A7FB5", alpha = 0.3) +
    geom_line(color = "#1e40af", linewidth = 0.4)

  if (nrow(blocks) > 0) {
    # Reuse the D01 outline color per block for consistency
    for (r in seq_len(nrow(outline_dt))) {
      od <- outline_dt[r]
      p_vote <- p_vote +
        geom_vline(
          xintercept = od$xmin,
          color = od$color,
          linewidth = od$lwd * 0.5,
          linetype = od$linetype
        ) +
        geom_vline(
          xintercept = od$xmax,
          color = od$color,
          linewidth = od$lwd * 0.5,
          linetype = od$linetype
        )
    }
  }

  if ("bg_level" %in% names(vote_track)) {
    vm <- max(vote_track$vote_smooth, na.rm = TRUE)
    if (is.finite(vm) && vm > 0) {
      p_vote <- p_vote +
        geom_line(aes(y = bg_level * vm),
                  color = "#dc2626", linewidth = 0.3, alpha = 0.5)
    }
  }

  # Mark SV breakpoints on the vote profile as short down-ticks on the x-axis
  if (nrow(sv_bp_dt) > 0) {
    for (s in seq_len(nrow(sv_bp_dt))) {
      bs <- sv_bp_dt[s]
      p_vote <- p_vote +
        annotate(
          "segment",
          x = bs$bp_mb, xend = bs$bp_mb,
          y = 0, yend = -0.04 * max(vote_track$vote_smooth, na.rm = TRUE),
          color = bs$tick_color,
          linewidth = STYLE_SV$lwd,
          alpha = STYLE_SV$alpha
        )
    }
  }

  p_vote <- p_vote +
    labs(
      x = paste0(chr_label, " (Mb)"),
      y = "Boundary votes",
      title = paste0(chr_label, " -- Vote Profile"),
      subtitle = "Blue area = vote_smooth | Red line = bg_level (family LD) | Vertical lines = D01 block boundaries | Bottom ticks = SV breakpoints"
    ) +
    theme_minimal(base_size = 9) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  for (ext in c("png", "pdf")) {
    out_file <- file.path(outdir, paste0(chr_label, "_vote_profile_v2.", ext))
    tryCatch({
      ggsave(out_file, p_vote, width = 14, height = 4, dpi = 350)
      message("[SAVED] ", out_file)
    }, error = function(e) message("[FAIL] ", ext, ": ", conditionMessage(e)))
  }
}

message("\n[DONE] Annotated plots -> ", outdir)
