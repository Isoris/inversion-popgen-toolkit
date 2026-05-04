#!/usr/bin/env Rscript
# ============================================================================
# STEP_D15_plot_zoomed_regions.R — Saturation-driven zoomed per-region plots
# ============================================================================
#
# IDEA: Order blocks by position. Accumulate how many unique windows
# they cover. The SHAPE of that curve tells you where to zoom:
#
#   - PLATEAU in accumulation = many blocks packing into the same space
#     (overlapping, nested, dense mini-inversions). One zoomed plot.
#   - STEEP RISE = new territory. Gap between zoom regions.
#
# Cut at plateaus → each plateau = one plotting region.
# Adapts to the actual structure. No fixed thresholds.
#
# Output:
#   saturation_curve_<chr>.png    — the accumulation plot itself
#   zoom_region_<N>_<chr>.png/pdf — one zoomed heatmap per region
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
steps_file <- NULL; landscape_file <- NULL; outdir <- "plots/zoomed"
chr_label <- NULL; max_zoom_bins <- 800L; gap_frac <- 0.05

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--precomp" && i < length(args))      { precomp_file <- args[i+1]; i <- i+2 }
  else if (a == "--blocks" && i < length(args))   { blocks_file <- args[i+1]; i <- i+2 }
  else if (a == "--votes" && i < length(args))    { votes_file <- args[i+1]; i <- i+2 }
  else if (a == "--steps" && i < length(args))    { steps_file <- args[i+1]; i <- i+2 }
  else if (a == "--landscape" && i < length(args)) { landscape_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))   { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--chr" && i < length(args))      { chr_label <- args[i+1]; i <- i+2 }
  else if (a == "--gap_frac" && i < length(args)) { gap_frac <- as.numeric(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
pc <- readRDS(precomp_file)
smat <- pc$sim_mat; dt <- pc$dt; N <- nrow(smat)
chr_label <- chr_label %||% pc$chrom
blocks <- fread(blocks_file)
votes <- if (!is.null(votes_file) && file.exists(votes_file)) fread(votes_file) else data.table()
steps <- if (!is.null(steps_file) && file.exists(steps_file)) fread(steps_file) else data.table()
landscape <- if (!is.null(landscape_file) && file.exists(landscape_file)) fread(landscape_file) else data.table()

window_bp <- 50000L
message("[ZOOM] ", chr_label, " | ", N, " windows | ", nrow(blocks), " blocks")

if (nrow(blocks) == 0) { message("No blocks. Done."); quit(save = "no") }

# ---- Compute depths if missing ----
if (!"depth" %in% names(blocks)) {
  blocks$depth <- 0L
  for (r in seq_len(nrow(blocks))) {
    pid <- blocks$parent_id[r]; d <- 0L
    while (!is.na(pid) && d < 10) {
      d <- d + 1L; pr <- which(blocks$block_id == pid)
      if (length(pr) == 0) break; pid <- blocks$parent_id[pr[1]]
    }
    blocks$depth[r] <- d
  }
}

# ---- System colors ----
SYSTEM_PALETTE <- c("#1a1a1a", "#1e40af", "#991b1b", "#166534", "#7c2d12",
                     "#6b21a8", "#0e7490", "#b45309", "#4338ca", "#be123c",
                     "#0f766e", "#a16207", "#7e22ce", "#0369a1", "#9f1239")
DEPTH_LWD <- c(1.3, 0.95, 0.7, 0.55, 0.4, 0.35, 0.3)
DEPTH_LINETYPE <- c("solid", "solid", "dashed", "dashed", "dotted", "dotted", "dotted")
DEPTH_LABEL_SIZE <- c(3.2, 2.6, 2.2, 1.9, 1.6, 1.4, 1.3)
SIM_COLORS <- c("#0C1E3C", "#1E3A5F", "#4A7FB5", "#C4A035", "#D4712A", "#B8282E")

get_root <- function(bid, blocks) {
  pid <- blocks$parent_id[blocks$block_id == bid]
  if (length(pid) == 0 || is.na(pid)) return(bid)
  get_root(pid, blocks)
}
sys_colors <- character(0)
if (nrow(blocks) > 0) {
  roots <- sapply(blocks$block_id, function(bid) get_root(bid, blocks))
  uq <- unique(roots)
  cmap <- setNames(SYSTEM_PALETTE[seq_along(uq) %% length(SYSTEM_PALETTE) + 1], as.character(uq))
  sys_colors <- setNames(cmap[as.character(roots)], as.character(blocks$block_id))
}


# ============================================================================
# STEP 1: Build the cumulative saturation curve
# ============================================================================

setorder(blocks, start)

# Mark which windows are covered by at least one block
covered <- rep(FALSE, N)
cum_new <- integer(nrow(blocks))  # how many NEW windows this block adds
cum_total <- integer(nrow(blocks))  # running total of covered windows

for (r in seq_len(nrow(blocks))) {
  bi <- blocks$start[r]; be <- blocks$end[r]
  new_windows <- sum(!covered[bi:be])
  covered[bi:be] <- TRUE
  cum_new[r] <- new_windows
  cum_total[r] <- sum(covered)
}

blocks$cum_new <- cum_new
blocks$cum_total <- cum_total
blocks$cum_frac <- cum_total / N

message("[SATURATION] Total coverage: ", sum(covered), "/", N,
        " (", round(100 * sum(covered) / N, 1), "%)")


# ============================================================================
# STEP 2: Plot the saturation curve
# ============================================================================

sat_dt <- data.table(
  block_order = seq_len(nrow(blocks)),
  block_id    = blocks$block_id,
  start_mb    = blocks$start_mb,
  cum_new     = cum_new,
  cum_total   = cum_total,
  cum_frac    = cum_total / N
)

p_sat <- ggplot(sat_dt, aes(x = block_order, y = cum_frac)) +
  geom_line(color = "#1e40af", linewidth = 0.6) +
  geom_point(aes(size = cum_new), color = "#1e40af", alpha = 0.5) +
  scale_size_continuous(range = c(0.3, 3), name = "New windows") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dotted", color = "#9ca3af") +
  labs(x = "Blocks (ordered by position)",
       y = "Cumulative fraction of chromosome covered",
       title = paste0(chr_label, " — Block Saturation Curve"),
       subtitle = paste0(nrow(blocks), " blocks cover ", round(100 * sum(covered)/N, 1),
                         "% of ", N, " windows"),
       caption = paste0("Plateaus = many blocks in same space (dense zone)\n",
                        "Steep rises = new territory | Flat top = saturation")) +
  theme_minimal(base_size = 9) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(face = "bold"))

for (ext in c("png", "pdf")) {
  f <- file.path(outdir, paste0(chr_label, "_saturation_curve.", ext))
  tryCatch({ ggsave(f, p_sat, width = 10, height = 5, dpi = 350)
    message("[SAVED] ", basename(f)) }, error = function(e) message("[FAIL] ", e$message))
}


# ============================================================================
# STEP 3: Cut the chromosome into zoom regions at saturation gaps
# ============================================================================
# A "gap" in the saturation = a block that adds many NEW windows
# (steep rise) = it's in uncovered territory = it starts a new zoom region.
#
# Threshold: if a block adds more than gap_frac * N new windows at once,
# it's entering new territory → cut here.

gap_threshold <- max(10L, round(N * gap_frac))

# Find region boundaries: where cum_new spikes above threshold
region_cuts <- c(1L)  # always start at block 1
for (r in 2:nrow(blocks)) {
  if (cum_new[r] > gap_threshold) {
    # This block jumps to new territory — start a new region
    region_cuts <- c(region_cuts, r)
  }
}
region_cuts <- c(region_cuts, nrow(blocks) + 1L)  # sentinel

n_regions <- length(region_cuts) - 1
message("[ZOOM] Saturation regions: ", n_regions, " (gap threshold: ",
        gap_threshold, " new windows)")

# ---- Bug 9 fix: if only 1 region covers >70% of chromosome, split it ----
# at the highest internal vote boundaries.
if (n_regions == 1 && sum(covered) > N * 0.7 && nrow(votes) > 0 &&
    "vote_smooth" %in% names(votes)) {
  message("[ZOOM] Single region covers >70% — splitting at highest vote peaks")

  reg_blocks <- blocks[region_cuts[1]:(region_cuts[2] - 1)]
  reg_start <- min(reg_blocks$start)
  reg_end   <- max(reg_blocks$end)

  # Find vote peaks within this region
  v_sub <- votes[position >= reg_start & position <= reg_end]
  if (nrow(v_sub) > 0) {
    # Take top 3-5 vote peaks as split points
    setorder(v_sub, -vote_smooth)
    n_splits <- min(4L, max(1L, nrow(v_sub) %/% 100))
    split_positions <- sort(v_sub$position[seq_len(n_splits)])

    # Convert positions to block indices
    new_cuts <- c(1L)
    for (sp in split_positions) {
      bi <- which.min(abs(blocks$start - sp))
      if (bi > 1 && bi < nrow(blocks)) new_cuts <- c(new_cuts, bi)
    }
    new_cuts <- sort(unique(c(new_cuts, nrow(blocks) + 1L)))

    if (length(new_cuts) > 2) {
      region_cuts <- new_cuts
      n_regions <- length(region_cuts) - 1
      message("[ZOOM] Split into ", n_regions, " sub-regions")
    }
  }
}


# ============================================================================
# STEP 4: Plot each zoom region
# ============================================================================

for (reg in seq_len(n_regions)) {
  reg_start_idx <- region_cuts[reg]
  reg_end_idx   <- region_cuts[reg + 1] - 1
  if (reg_end_idx < reg_start_idx) next

  reg_blocks <- blocks[reg_start_idx:reg_end_idx]
  n_blocks_in_region <- nrow(reg_blocks)

  # Genomic range of this region
  zoom_start <- min(reg_blocks$start)
  zoom_end   <- max(reg_blocks$end)
  zoom_width <- zoom_end - zoom_start + 1

  # Skip tiny regions (< 10 bins)
  if (zoom_width < 10) next

  # Add 10% padding
  pad <- max(10L, zoom_width %/% 10)
  zoom_start <- max(1L, zoom_start - pad)
  zoom_end   <- min(N, zoom_end + pad)
  zoom_width <- zoom_end - zoom_start + 1

  start_mb <- round((zoom_start - 1) * window_bp / 1e6, 1)
  end_mb   <- round(zoom_end * window_bp / 1e6, 1)
  span_mb  <- round(end_mb - start_mb, 1)

  message("  Region ", reg, ": bins ", zoom_start, "-", zoom_end,
          " (", start_mb, "-", end_mb, " Mb, ", n_blocks_in_region, " blocks)")

  # Subsample if needed
  step_s <- max(1L, zoom_width %/% max_zoom_bins)
  local_idx <- seq(zoom_start, zoom_end, by = step_s)
  ns <- length(local_idx)
  pos_mb <- (local_idx - 1) * window_bp / 1e6

  # Build heatmap
  hm_rows <- list(); k <- 0L
  for (ii in seq_len(ns)) {
    for (jj in seq_len(ns)) {
      v <- smat[local_idx[ii], local_idx[jj]]
      if (is.finite(v)) { k <- k + 1L; hm_rows[[k]] <- list(x = pos_mb[ii], y = pos_mb[jj], value = v) }
    }
  }
  hm_dt <- rbindlist(hm_rows)
  if (nrow(hm_dt) == 0) next

  sim_med <- median(hm_dt$value, na.rm = TRUE)
  sim_scale <- scale_fill_gradientn(
    colours = SIM_COLORS,
    values = scales::rescale(c(0, sim_med*0.5, sim_med*0.85, sim_med*1.1, sim_med*1.3, 1)),
    name = "Sim", limits = c(0, 1))

  p <- ggplot() +
    geom_raster(data = hm_dt, aes(x = x, y = y, fill = value), interpolate = TRUE) +
    sim_scale

  # All blocks that overlap this zoom range
  all_zoom_blocks <- blocks[start <= zoom_end & end >= zoom_start]

  if (nrow(all_zoom_blocks) > 0) {
    for (br in seq_len(nrow(all_zoom_blocks))) {
      b <- all_zoom_blocks[br]
      d <- min((b$depth %||% 0L) + 1L, length(DEPTH_LWD))
      bid_str <- as.character(b$block_id)
      color <- if (bid_str %in% names(sys_colors)) sys_colors[bid_str] else "#1a1a1a"
      x_lo <- b$start_mb; x_hi <- b$end_mb

      lbl <- paste0("B", b$block_id)
      if (nrow(landscape) > 0) {
        ls <- landscape[block_id == b$block_id]
        if (nrow(ls) > 0) lbl <- ls$label[1]
      }

      p <- p +
        annotate("rect", xmin = x_lo, xmax = x_hi, ymin = x_lo, ymax = x_hi,
                 fill = NA, color = color,
                 linewidth = DEPTH_LWD[d], linetype = DEPTH_LINETYPE[d]) +
        annotate("text", x = x_lo + (x_hi - x_lo) * 0.02,
                 y = x_hi - (x_hi - x_lo) * 0.04,
                 label = lbl, size = DEPTH_LABEL_SIZE[d],
                 color = color, hjust = 0, vjust = 1, fontface = "bold", alpha = 0.9)
    }
  }

  # Internal features
  if (nrow(steps) > 0 && "inside_block" %in% names(steps)) {
    int_s <- steps[inside_block == TRUE & step_position >= zoom_start & step_position <= zoom_end]
    if (nrow(int_s) > 0) {
      up <- sort(unique(int_s$step_position))
      mm <- (up - 1) * window_bp / 1e6
      p <- p + annotate("point", x = mm, y = mm, shape = 4, size = 1.2,
                         color = "#dc2626", alpha = 0.7)
    }
  }

  p <- p + coord_fixed() +
    labs(x = paste0(chr_label, " (Mb)"), y = paste0(chr_label, " (Mb)"),
         title = paste0(chr_label, " — Region ", reg, " of ", n_regions,
                        " [", start_mb, "-", end_mb, " Mb]"),
         subtitle = paste0(n_blocks_in_region, " blocks (", nrow(all_zoom_blocks),
                          " overlap) | ", span_mb, " Mb | ",
                          ns, "×", ns, " display")) +
    theme_minimal(base_size = 9) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, color = "grey40"))

  fname <- paste0(chr_label, "_zoom_region_", reg)
  for (ext in c("png", "pdf")) {
    f <- file.path(outdir, paste0(fname, ".", ext))
    tryCatch({ ggsave(f, p, width = 12, height = 11, dpi = 350)
      message("  [SAVED] ", basename(f)) }, error = function(e) message("  [FAIL] ", e$message))
  }
}

message("\n[DONE] ", n_regions, " zoom regions + saturation curve -> ", outdir)
