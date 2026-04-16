#!/usr/bin/env Rscript
# ============================================================================
# STEP_D13_plot_annotated_simmat.R — Chromosome-wide annotated sim_mat heatmap
# ============================================================================
#
# Reads detector v9.3 output (blocks, vote profile, step table) and the
# precomp RDS. Plots the sim_mat heatmap with:
#
#   - Nested block outlines: outer blocks in thick solid lines,
#     inner blocks in thinner lines, deeper nesting = thinner + dashed
#   - Color per nesting level: level 0 = black, level 1 = dark blue,
#     level 2 = dark red, level 3 = dark green, etc.
#   - Text labels directly on the heatmap: block ID, Mb range, type
#   - Vote profile as a marginal track
#   - Background level track (family LD indicator)
#   - Internal features marked with small crosses
#
# Usage:
#   Rscript 13_plot_annotated_simmat.R \
#     --precomp precomp/C_gar_LG01.precomp.rds \
#     --blocks inv_detect_out_v9.3/blocks_C_gar_LG01.tsv \
#     --votes inv_detect_out_v9.3/*_votes.tsv.gz \
#     --steps inv_detect_out_v9.3/*_steps.tsv.gz \
#     --outdir plots/
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- CLI ----
args <- commandArgs(trailingOnly = TRUE)
precomp_file <- NULL; blocks_file <- NULL; votes_file <- NULL
steps_file <- NULL; outdir <- "plots"; chr_label <- NULL
max_display <- 600L  # subsample to max_display × max_display for plotting

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--precomp" && i < length(args))  { precomp_file <- args[i+1]; i <- i+2 }
  else if (a == "--blocks" && i < length(args)) { blocks_file <- args[i+1]; i <- i+2 }
  else if (a == "--votes" && i < length(args))  { votes_file <- args[i+1]; i <- i+2 }
  else if (a == "--steps" && i < length(args))  { steps_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args)) { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--chr" && i < length(args))    { chr_label <- args[i+1]; i <- i+2 }
  else if (a == "--max_display" && i < length(args)) { max_display <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# COLOR AND STYLE CONFIG
# ============================================================================

# System palette: each ROOT block (independent inversion system) gets one color.
# All children/grandchildren inherit the same color.
# Color identifies the SYSTEM. Thickness identifies the DEPTH.
SYSTEM_PALETTE <- c(
  "#1a1a1a",  # system 1: near-black
  "#1e40af",  # system 2: dark blue
  "#991b1b",  # system 3: dark red
  "#166534",  # system 4: dark green
  "#7c2d12",  # system 5: dark brown
  "#6b21a8",  # system 6: purple
  "#0e7490",  # system 7: dark cyan
  "#b45309",  # system 8: amber
  "#4338ca",  # system 9: indigo
  "#be123c",  # system 10: rose
  "#0f766e",  # system 11: teal
  "#a16207",  # system 12: yellow-brown
  "#7e22ce",  # system 13: violet
  "#0369a1",  # system 14: sky blue
  "#9f1239"   # system 15: deep pink
)

# For gap/non-block regions
GAP_COLORS <- c(
  family_ld_band       = "#9ca3af",
  diffuse_diagonal     = "#78716c",
  extended_suppression = "#0d9488",
  clean_background     = "#d4d4d8",
  transition_zone      = "#a3a3a3",
  mini_block_zone      = "#ea580c"
)

# Line width by nesting depth (decreasing) — same color, thinner lines deeper
DEPTH_LWD <- c(1.3, 0.95, 0.7, 0.55, 0.4, 0.35, 0.3)

# Linetype by nesting depth — solid then dashed then dotted
DEPTH_LINETYPE <- c("solid", "solid", "dashed", "dashed", "dotted", "dotted", "dotted")

# Label size by depth — smaller for deeper nesting
DEPTH_LABEL_SIZE <- c(2.8, 2.3, 2.0, 1.7, 1.5, 1.3, 1.2)

# Sim_mat color scale
SIM_COLORS <- c("#0C1E3C", "#1E3A5F", "#4A7FB5", "#C4A035", "#D4712A", "#B8282E")


# ============================================================================
# LOAD DATA
# ============================================================================

# Load precomp
if (!is.null(precomp_file) && file.exists(precomp_file)) {
  pc <- readRDS(precomp_file)
  smat <- pc$sim_mat
  dt <- pc$dt
  chr_label <- chr_label %||% pc$chrom
  N <- nrow(smat)
  message("[PLOT] ", chr_label, " | ", N, " windows")
} else {
  stop("--precomp required")
}

# Load blocks
blocks <- if (!is.null(blocks_file) && file.exists(blocks_file)) {
  fread(blocks_file)
} else {
  data.table()
}

# Load votes
votes <- if (!is.null(votes_file) && file.exists(votes_file)) {
  fread(votes_file)
} else {
  data.table()
}

# Load steps (for internal feature marks)
steps <- if (!is.null(steps_file) && file.exists(steps_file)) {
  fread(steps_file)
} else {
  data.table()
}

# Load or compute landscape classification
landscape <- data.table()
landscape_file <- sub("blocks_", "landscape_", blocks_file %||% "")
if (file.exists(landscape_file)) {
  landscape <- fread(landscape_file)
  message("[PLOT] Loaded landscape classification: ", nrow(landscape), " entries")
} else if (nrow(blocks) > 0) {
  # Try to compute on the fly if the classifier is available
  classifier_file <- file.path(dirname(sys.frame(1)$ofile %||% "."), "14_landscape_classifier.R")
  if (file.exists(classifier_file)) {
    source(classifier_file)
    landscape <- classify_landscape(blocks, smat)
    message("[PLOT] Computed landscape classification: ", nrow(landscape), " entries")
  }
}


# ============================================================================
# COMPUTE NESTING DEPTH PER BLOCK
# ============================================================================

compute_depths <- function(blocks) {
  if (nrow(blocks) == 0) return(blocks)
  blocks$depth <- 0L
  for (r in seq_len(nrow(blocks))) {
    pid <- blocks$parent_id[r]
    d <- 0L
    while (!is.na(pid) && d < 10) {
      d <- d + 1L
      parent_row <- which(blocks$block_id == pid)
      if (length(parent_row) == 0) break
      pid <- blocks$parent_id[parent_row[1]]
    }
    blocks$depth[r] <- d
  }
  return(blocks)
}

if (nrow(blocks) > 0) blocks <- compute_depths(blocks)


# ============================================================================
# BUILD HEATMAP DATA (subsampled for display)
# ============================================================================

step_size <- max(1L, N %/% max_display)
idx <- seq(1L, N, by = step_size)
ns <- length(idx)

# Mb positions
pos_mb <- if ("start_bp" %in% names(dt)) {
  (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6
} else {
  idx * 50000 / 1e6  # fallback
}

# Build full matrix data
message("[PLOT] Building heatmap data (", ns, "×", ns, ")...")
hm_rows <- list()
k <- 0L
for (ii in seq_len(ns)) {
  for (jj in seq_len(ns)) {
    v <- smat[idx[ii], idx[jj]]
    if (is.finite(v)) {
      k <- k + 1L
      hm_rows[[k]] <- list(x = pos_mb[ii], y = pos_mb[jj], value = v)
    }
  }
}
hm_dt <- rbindlist(hm_rows)

sim_med <- median(hm_dt$value, na.rm = TRUE)

# Sim_mat color scale
sim_scale <- scale_fill_gradientn(
  colours = SIM_COLORS,
  values = scales::rescale(c(0, sim_med * 0.5, sim_med * 0.85,
                              sim_med * 1.1, sim_med * 1.3, 1)),
  name = "Similarity", limits = c(0, 1)
)


# ============================================================================
# ASSIGN SYSTEM COLORS — each root block gets a unique color,
# all its children/grandchildren inherit the same color
# ============================================================================

assign_system_colors <- function(blocks) {
  if (nrow(blocks) == 0) return(character(0))

  # Find root for each block (walk up parent chain)
  block_root <- integer(nrow(blocks))
  for (r in seq_len(nrow(blocks))) {
    bid <- blocks$block_id[r]
    pid <- blocks$parent_id[r]
    root <- bid
    steps <- 0L
    while (!is.na(pid) && steps < 20) {
      root <- pid
      pr <- which(blocks$block_id == pid)
      if (length(pr) == 0) break
      pid <- blocks$parent_id[pr[1]]
      steps <- steps + 1L
    }
    block_root[r] <- root
  }

  # Assign colors: one per unique root
  unique_roots <- unique(block_root)
  root_color_map <- setNames(
    SYSTEM_PALETTE[seq_along(unique_roots) %% length(SYSTEM_PALETTE) + 1],
    as.character(unique_roots)
  )

  # Map each block to its root's color
  colors <- root_color_map[as.character(block_root)]
  names(colors) <- as.character(blocks$block_id)
  return(colors)
}


# ============================================================================
# BUILD BLOCK OUTLINE RECTANGLES
# ============================================================================

build_block_outlines <- function(blocks, landscape = NULL) {
  if (nrow(blocks) == 0) return(data.table())

  # Assign system colors (root-based)
  sys_colors <- assign_system_colors(blocks)

  outlines <- list()
  for (r in seq_len(nrow(blocks))) {
    b <- blocks[r]
    d <- min((b$depth %||% 0L) + 1L, length(DEPTH_LWD))

    x_lo <- b$start_mb
    x_hi <- b$end_mb

    # System color from root assignment
    bid_str <- as.character(b$block_id)
    color <- if (bid_str %in% names(sys_colors)) sys_colors[bid_str] else "#1a1a1a"

    # Get rich label from landscape if available
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

    if (is.null(rich_label)) rich_label <- build_block_label(b)

    outlines[[r]] <- data.table(
      block_id    = b$block_id,
      xmin        = x_lo,
      xmax        = x_hi,
      ymin        = x_lo,
      ymax        = x_hi,
      depth       = b$depth %||% 0L,
      category    = cat_name,
      confidence  = confidence,
      color       = as.character(color),
      lwd         = DEPTH_LWD[d],
      linetype    = DEPTH_LINETYPE[d],
      label_size  = DEPTH_LABEL_SIZE[d],
      label       = rich_label,
      label_x     = x_lo + (x_hi - x_lo) * 0.02,
      label_y     = x_hi - (x_hi - x_lo) * 0.05
    )
  }
  rbindlist(outlines)
}


# ============================================================================
# MINI-BLOCK DENSITY ZONES — detect regions with many tiny squares
# ============================================================================
# In the diffuse diagonal, there are often many small blocks packed together.
# Instead of labeling each one individually (unreadable), detect ZONES of
# high block density and label the zone.

build_density_zones <- function(blocks, window_size_bp = 50000L,
                                 zone_window_mb = 2.0, min_blocks_per_zone = 3L) {
  if (nrow(blocks) == 0) return(data.table())

  # Only consider small blocks (width < 50 bins = ~2.5 Mb)
  small <- blocks[width < 50]
  if (nrow(small) < min_blocks_per_zone) return(data.table())

  # Bin small blocks by position (zone_window_mb sliding window)
  max_mb <- max(small$end_mb, na.rm = TRUE)
  zones <- list()

  zone_start <- 0
  while (zone_start < max_mb) {
    zone_end <- zone_start + zone_window_mb
    zone_blocks <- small[start_mb >= zone_start & end_mb <= zone_end]
    n_in_zone <- nrow(zone_blocks)

    if (n_in_zone >= min_blocks_per_zone) {
      zones[[length(zones) + 1]] <- data.table(
        zone_start_mb = zone_start,
        zone_end_mb   = zone_end,
        n_blocks      = n_in_zone,
        density       = round(n_in_zone / zone_window_mb, 1),
        mean_width    = round(mean(zone_blocks$width), 1),
        label         = paste0("Dense zone: ", n_in_zone, " mini-blocks / ",
                              zone_window_mb, " Mb\n",
                              round(zone_start, 1), "-", round(zone_end, 1), " Mb")
      )
    }
    zone_start <- zone_start + zone_window_mb * 0.5  # 50% overlap
  }

  if (length(zones) > 0) rbindlist(zones) else data.table()
}

build_block_label <- function(b) {
  # Build a concise label for the block
  type_short <- if ("shape_class" %in% names(b) && !is.na(b$shape_class)) {
    switch(as.character(b$shape_class),
      strong_square = "Strong",
      diffuse_square = "Diffuse",
      diagonal_band = "Band",
      noise = "Noise",
      "")
  } else ""

  parent_info <- if (!is.na(b$parent_id)) {
    paste0("child of B", b$parent_id)
  } else "outer"

  span <- round(b$end_mb - b$start_mb, 1)

  paste0("B", b$block_id, " [", round(b$start_mb, 1), "-",
         round(b$end_mb, 1), " Mb] ",
         span, "Mb ", type_short, "\n",
         parent_info,
         if (!is.na(b$height)) paste0(" h=", b$height) else "")
}

outline_dt <- build_block_outlines(blocks, landscape)
density_zones <- build_density_zones(blocks)

if (nrow(density_zones) > 0) {
  message("[PLOT] Mini-block density zones: ", nrow(density_zones))
}


# ============================================================================
# BUILD INTERNAL FEATURE MARKS
# ============================================================================

build_internal_marks <- function(steps, blocks) {
  if (nrow(steps) == 0 || nrow(blocks) == 0) return(data.table())
  if (!"inside_block" %in% names(steps)) return(data.table())

  internal <- steps[inside_block == TRUE]
  if (nrow(internal) == 0) return(data.table())

  # Unique positions of internal steps
  unique_pos <- sort(unique(internal$step_position))

  # Convert to Mb (approximate: bin * 50kb)
  marks <- data.table(
    position_mb = unique_pos * 50000 / 1e6,
    type = "internal_feature"
  )
  return(marks)
}

marks_dt <- build_internal_marks(steps, blocks)


# ============================================================================
# BUILD VOTE TRACK (marginal)
# ============================================================================

build_vote_track <- function(votes) {
  if (nrow(votes) == 0) return(data.table())
  # Subsample to match display
  votes_sub <- votes[seq(1L, nrow(votes), length.out = min(ns, nrow(votes))), ]
  votes_sub[, pos_mb := position * 50000 / 1e6]
  return(votes_sub)
}

vote_track <- build_vote_track(votes)


# ============================================================================
# ASSEMBLE THE PLOT
# ============================================================================

message("[PLOT] Assembling annotated heatmap...")

p <- ggplot()

# Layer 1: Sim_mat heatmap
p <- p + geom_raster(data = hm_dt, aes(x = x, y = y, fill = value),
                      interpolate = TRUE) +
  sim_scale

# Layer 1.5: Mini-block density zones (semi-transparent orange rectangles)
if (nrow(density_zones) > 0) {
  for (dz in seq_len(nrow(density_zones))) {
    dzr <- density_zones[dz]
    p <- p + annotate(
      "rect",
      xmin = dzr$zone_start_mb, xmax = dzr$zone_end_mb,
      ymin = dzr$zone_start_mb, ymax = dzr$zone_end_mb,
      fill = "#ea580c", alpha = 0.08, color = "#ea580c",
      linewidth = 0.4, linetype = "dotted"
    ) + annotate(
      "text",
      x = dzr$zone_start_mb, y = dzr$zone_end_mb,
      label = dzr$label,
      size = 1.8, color = "#ea580c",
      hjust = 0, vjust = 1, fontface = "italic"
    )
  }
}

# Layer 2: Block outlines — each block gets its own color from system root
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

  # Layer 3: Block labels
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
      alpha = 0.85
    )
  }
}

# Layer 4: Internal feature marks (small X marks on the diagonal)
if (nrow(marks_dt) > 0) {
  p <- p + geom_point(
    data = marks_dt,
    aes(x = position_mb, y = position_mb),
    shape = 4,  # X mark
    size = 0.8,
    color = "#dc2626",
    alpha = 0.6
  )
}

# Cosmetics
p <- p +
  coord_fixed() +
  labs(
    x = paste0(chr_label, " (Mb)"),
    y = paste0(chr_label, " (Mb)"),
    title = paste0(chr_label, " — Annotated Similarity Matrix"),
    subtitle = paste0(N, " windows (", ns, "×", ns, " display) | ",
                     nrow(blocks), " blocks detected | ",
                     "median sim = ", round(sim_med, 3)),
    caption = paste0(
      "Same color = same inversion system (root + all children)\n",
      "Thick solid = outer block | Thin dashed = inner nested | Dotted = deep\n",
      "Orange dotted zone = high density of mini-blocks\n",
      "Red x = internal feature | [%] = confidence score"
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    plot.caption = element_text(size = 6, hjust = 0, color = "#8A8A8A"),
    legend.key.size = unit(0.4, "cm")
  )


# ============================================================================
# SAVE
# ============================================================================

out_base <- paste0(chr_label, "_annotated_simmat")

for (ext in c("png", "pdf")) {
  out_file <- file.path(outdir, paste0(out_base, ".", ext))
  tryCatch({
    ggsave(out_file, p, width = 14, height = 13, dpi = 350)
    message("[SAVED] ", out_file)
  }, error = function(e) message("[FAIL] ", ext, ": ", e$message))
}


# ============================================================================
# ALSO SAVE: Vote profile as separate track plot
# ============================================================================

if (nrow(vote_track) > 0 && "vote_smooth" %in% names(vote_track)) {
  p_vote <- ggplot(vote_track, aes(x = pos_mb, y = vote_smooth)) +
    geom_area(fill = "#4A7FB5", alpha = 0.3) +
    geom_line(color = "#1e40af", linewidth = 0.4)

  # Add boundary positions as vertical lines
  if (nrow(blocks) > 0) {
    for (r in seq_len(nrow(blocks))) {
      b <- blocks[r]
      d <- min(b$depth + 1L, length(DEPTH_COLORS))
      p_vote <- p_vote +
        geom_vline(xintercept = b$start_mb, color = DEPTH_COLORS[d],
                   linewidth = DEPTH_LWD[d] * 0.5, linetype = DEPTH_LINETYPE[d]) +
        geom_vline(xintercept = b$end_mb, color = DEPTH_COLORS[d],
                   linewidth = DEPTH_LWD[d] * 0.5, linetype = DEPTH_LINETYPE[d])
    }
  }

  # Background level track
  if ("bg_level" %in% names(vote_track)) {
    p_vote <- p_vote +
      geom_line(aes(y = bg_level * max(vote_smooth, na.rm = TRUE)),
                color = "#dc2626", linewidth = 0.3, alpha = 0.5)
  }

  p_vote <- p_vote +
    labs(x = paste0(chr_label, " (Mb)"), y = "Boundary votes",
         title = paste0(chr_label, " — Vote Profile"),
         subtitle = "Blue = boundary agreement | Red = background level (family LD indicator)") +
    theme_minimal(base_size = 9) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  for (ext in c("png", "pdf")) {
    out_file <- file.path(outdir, paste0(chr_label, "_vote_profile.", ext))
    tryCatch({
      ggsave(out_file, p_vote, width = 14, height = 4, dpi = 350)
      message("[SAVED] ", out_file)
    }, error = function(e) message("[FAIL] ", ext, ": ", e$message))
  }
}

message("\n[DONE] Annotated plots -> ", outdir)
