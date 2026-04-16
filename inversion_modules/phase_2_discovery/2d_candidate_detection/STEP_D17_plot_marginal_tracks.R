#!/usr/bin/env Rscript
# ============================================================================
# STEP_D17_plot_marginal_tracks.R — Marginal-track chromosome plot with geometric upper/lower triangles
# ============================================================================
#
# Produces a multi-panel figure with:
#   TOP:    Nesting depth track (TH-score style from Hi-C TAD papers)
#   LEFT:   Same nesting depth track (vertical, mirrored)
#   CENTER: Split heatmap (upper=raw, lower=distance-corrected)
#   BOTTOM: Vote profile + cumulative block coverage curve
#
# The nesting depth track:
#   For each window position, count how many block outlines it sits inside.
#   Position in background = 0. Inside one block = 1. Inside child of parent = 2.
#   Weighted version: deeper nesting = higher score (like TH score).
#   Plot as a colored bar: white=0, light=1, medium=2, dark=3+.
#
# The split heatmap:
#   Upper triangle = raw sim_mat (orange/yellow scale)
#   Lower triangle = distance-corrected sim_mat (red-white-blue diverging)
#   This directly shows what the distance correction reveals.
#
# Usage:
#   Rscript 17_plot_marginal_tracks.R \
#     --precomp precomp/C_gar_LG01.precomp.rds \
#     --blocks blocks_C_gar_LG01.tsv \
#     --votes *_votes.tsv.gz \
#     [--distcorr distcorr_matrix.rds]
#     --outdir plots/
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
precomp_file <- NULL; blocks_file <- NULL; votes_file <- NULL
distcorr_file <- NULL; outdir <- "plots"; chr_label <- NULL
max_display <- 600L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--precomp" && i < length(args))   { precomp_file <- args[i+1]; i <- i+2 }
  else if (a == "--blocks" && i < length(args)) { blocks_file <- args[i+1]; i <- i+2 }
  else if (a == "--votes" && i < length(args))  { votes_file <- args[i+1]; i <- i+2 }
  else if (a == "--distcorr" && i < length(args)) { distcorr_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args)) { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--chr" && i < length(args))    { chr_label <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load ----
pc <- readRDS(precomp_file)
smat <- pc$sim_mat; dt <- pc$dt; N <- nrow(smat)
chr_label <- chr_label %||% pc$chrom
blocks <- if (!is.null(blocks_file) && file.exists(blocks_file)) fread(blocks_file) else data.table()
votes <- if (!is.null(votes_file) && file.exists(votes_file)) fread(votes_file) else data.table()

window_bp <- 50000L
message("[TRACKS] ", chr_label, " | N=", N, " | ", nrow(blocks), " blocks")

# ---- Compute nesting depths if missing ----
if (nrow(blocks) > 0 && !"depth" %in% names(blocks)) {
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


# ============================================================================
# 1. NESTING DEPTH TRACK (TH-score style)
# ============================================================================
# For each window, count how many blocks contain it.
# Simple version: count = number of containing blocks.
# Weighted version (TH-score): sum(1/rank_of_block) or just max depth.

nesting_count <- integer(N)     # how many blocks contain this window
nesting_max_depth <- integer(N)  # deepest nesting level at this window

if (nrow(blocks) > 0) {
  for (r in seq_len(nrow(blocks))) {
    bi <- blocks$start[r]; be <- blocks$end[r]
    bi <- max(1L, bi); be <- min(N, be)
    nesting_count[bi:be] <- nesting_count[bi:be] + 1L
    d <- blocks$depth[r] + 1L  # depth 0 = level 1
    nesting_max_depth[bi:be] <- pmax(nesting_max_depth[bi:be], d)
  }
}

# Subsample for display
step_s <- max(1L, N %/% max_display)
idx <- seq(1L, N, by = step_s)
ns <- length(idx)
pos_mb <- (idx - 1) * window_bp / 1e6

depth_dt <- data.table(
  pos_mb = pos_mb,
  nesting = nesting_count[idx],
  max_depth = nesting_max_depth[idx]
)

# Nesting depth bar plot
DEPTH_PALETTE <- c("0" = "#f5f5f5", "1" = "#fee2e2", "2" = "#fca5a5",
                    "3" = "#ef4444", "4" = "#b91c1c", "5" = "#7f1d1d")

depth_dt[, depth_cat := as.character(pmin(max_depth, 5))]

p_depth <- ggplot(depth_dt, aes(x = pos_mb, y = 1, fill = depth_cat)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = DEPTH_PALETTE, name = "Nesting\ndepth") +
  labs(x = NULL, y = NULL,
       title = paste0(chr_label, " — Nesting Depth (TH-score style)")) +
  theme_minimal(base_size = 8) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA))

# Also as a line plot showing nesting count
p_nesting_line <- ggplot(depth_dt, aes(x = pos_mb, y = nesting)) +
  geom_area(fill = "#ef4444", alpha = 0.3) +
  geom_line(color = "#b91c1c", linewidth = 0.4) +
  labs(x = paste0(chr_label, " (Mb)"), y = "# blocks\ncontaining",
       subtitle = "Higher = deeper nesting = complex inversion system") +
  theme_minimal(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", color = NA))


# ============================================================================
# 2. CUMULATIVE BLOCK COVERAGE CURVE
# ============================================================================

if (nrow(blocks) > 0) {
  setorder(blocks, start)
  covered <- rep(FALSE, N)
  cum_frac <- numeric(nrow(blocks))
  for (r in seq_len(nrow(blocks))) {
    bi <- max(1L, blocks$start[r]); be <- min(N, blocks$end[r])
    covered[bi:be] <- TRUE
    cum_frac[r] <- sum(covered) / N
  }

  cum_dt <- data.table(
    block_order = seq_len(nrow(blocks)),
    pos_mb = blocks$start_mb,
    cum_frac = cum_frac,
    width = blocks$width
  )

  p_cumul <- ggplot(cum_dt, aes(x = pos_mb, y = cum_frac)) +
    geom_step(color = "#1e40af", linewidth = 0.5) +
    geom_point(aes(size = width), color = "#1e40af", alpha = 0.4) +
    scale_size_continuous(range = c(0.3, 2.5), guide = "none") +
    geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dotted", color = "#9ca3af") +
    labs(x = paste0(chr_label, " (Mb)"), y = "Cumul.\ncoverage",
         subtitle = "Plateaus = dense blocks | Steps = new territory") +
    theme_minimal(base_size = 8) +
    theme(plot.background = element_rect(fill = "white", color = NA))
} else {
  p_cumul <- ggplot() + theme_void()
}


# ============================================================================
# 2b. LOCAL CONTRAST RATIO — their intra_max / inter formula
# ============================================================================
# At each position p, with depth d bins:
#   L = mean(smat[p-d:p, p-d:p])   — left triangle (upstream of p)
#   R = mean(smat[p:p+d, p:p+d])   — right triangle (downstream of p)
#   X = mean(smat[p-d:p, p:p+d])   — cross region (between L and R)
#   intra_max = max(L, R)
#   ratio = intra_max / X           — high at boundaries, ~1 in background
#
# This is the direct analogue of the Hi-C boundary scores from image 7-8.

d_contrast <- max(10L, N %/% 100)  # depth for L/R/X triangles
cat("[TRACKS] Computing local contrast ratio (d=", d_contrast, ")...\n")

contrast_ratio <- numeric(N)
intra_max_vec  <- numeric(N)
inter_vec      <- numeric(N)

for (p in seq_len(N)) {
  l_lo <- max(1L, p - d_contrast); l_hi <- p
  r_lo <- p; r_hi <- min(N, p + d_contrast)

  # L = left triangle
  if (l_hi - l_lo >= 2) {
    L_sub <- smat[l_lo:l_hi, l_lo:l_hi]
    L_vals <- L_sub[upper.tri(L_sub)]
    L_vals <- L_vals[is.finite(L_vals)]
    L_mean <- if (length(L_vals) > 0) mean(L_vals) else NA
  } else L_mean <- NA

  # R = right triangle
  if (r_hi - r_lo >= 2) {
    R_sub <- smat[r_lo:r_hi, r_lo:r_hi]
    R_vals <- R_sub[upper.tri(R_sub)]
    R_vals <- R_vals[is.finite(R_vals)]
    R_mean <- if (length(R_vals) > 0) mean(R_vals) else NA
  } else R_mean <- NA

  # X = cross region
  if (l_hi > l_lo && r_hi > r_lo) {
    X_sub <- smat[l_lo:l_hi, r_lo:r_hi]
    X_vals <- X_sub[is.finite(X_sub)]
    X_mean <- if (length(X_vals) > 0) mean(X_vals) else NA
  } else X_mean <- NA

  if (!is.na(L_mean) && !is.na(R_mean) && !is.na(X_mean) && X_mean > 0.001) {
    intra_max_vec[p] <- max(L_mean, R_mean)
    inter_vec[p] <- X_mean
    contrast_ratio[p] <- max(L_mean, R_mean) / X_mean
  }
}

# Subsample
cr_dt <- data.table(
  pos_mb = pos_mb,
  ratio = contrast_ratio[idx],
  intra_max = intra_max_vec[idx],
  inter = inter_vec[idx]
)

p_contrast <- ggplot(cr_dt, aes(x = pos_mb, y = ratio)) +
  geom_line(color = "#7c2d12", linewidth = 0.4) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "#9ca3af") +
  labs(x = paste0(chr_label, " (Mb)"), y = "intra/inter\nratio",
       subtitle = "Peaks = block boundaries (high L or R vs low cross)") +
  theme_minimal(base_size = 8) +
  theme(plot.background = element_rect(fill = "white", color = NA))

if (nrow(votes) > 0 && "vote_smooth" %in% names(votes)) {
  vote_sub <- votes[seq(1L, nrow(votes), length.out = min(ns, nrow(votes)))]
  vote_sub[, pos_mb := (position - 1) * window_bp / 1e6]

  p_votes <- ggplot(vote_sub, aes(x = pos_mb, y = vote_smooth)) +
    geom_area(fill = "#4A7FB5", alpha = 0.3) +
    geom_line(color = "#1e40af", linewidth = 0.3) +
    labs(x = paste0(chr_label, " (Mb)"), y = "Boundary\nvotes") +
    theme_minimal(base_size = 8) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  # Background level track
  if ("bg_level" %in% names(vote_sub)) {
    p_bg <- ggplot(vote_sub, aes(x = pos_mb, y = bg_level)) +
      geom_line(color = "#dc2626", linewidth = 0.3) +
      geom_hline(yintercept = median(vote_sub$bg_level, na.rm = TRUE),
                 linetype = "dashed", color = "#9ca3af", linewidth = 0.3) +
      labs(x = paste0(chr_label, " (Mb)"), y = "Background\nlevel",
           subtitle = "High = family LD (elevated far-distance similarity)") +
      theme_minimal(base_size = 8) +
      theme(plot.background = element_rect(fill = "white", color = NA))
  } else {
    p_bg <- NULL
  }
} else {
  p_votes <- ggplot() + theme_void()
  p_bg <- NULL
}


# ============================================================================
# 4. SPLIT HEATMAP (upper=raw, lower=distance-corrected)
# ============================================================================

SIM_COLORS <- c("#0C1E3C", "#1E3A5F", "#4A7FB5", "#C4A035", "#D4712A", "#B8282E")
DIV_COLORS <- c("#2166ac", "#4393c3", "#92c5de", "#f7f7f7", "#f4a582", "#d6604d", "#b2182b")

# Build upper triangle (raw)
upper_rows <- list(); k <- 0L
for (ii in seq_len(ns)) {
  for (jj in ii:ns) {
    v <- smat[idx[ii], idx[jj]]
    if (is.finite(v)) {
      k <- k + 1L
      upper_rows[[k]] <- list(x = pos_mb[ii], y = pos_mb[jj], value = v)
    }
  }
}
upper_dt <- rbindlist(upper_rows)
sim_med <- median(upper_dt$value, na.rm = TRUE)

# Try to load or compute distcorr for lower triangle
distcorr_mat <- NULL
if (!is.null(distcorr_file) && file.exists(distcorr_file)) {
  distcorr_mat <- readRDS(distcorr_file)
  message("[TRACKS] Loaded distance-corrected matrix")
} else {
  # Try computing A1 on the fly (fast with vectorized version)
  mtx_file <- file.path(dirname(sys.frame(1)$ofile %||% "."), "STEP_D07_matrix_transforms.R")
  if (file.exists(mtx_file)) {
    # Only source if we have config
    if (!exists("CFG")) {
      CFG <- list(MTX_DISTCORR_ENABLED = TRUE, MTX_LOCALNORM_ENABLED = FALSE,
                  MTX_DENOISED_ENABLED = FALSE, MTX_RESIDBG_ENABLED = FALSE,
                  MTX_EDGE_ENABLED = FALSE, MTX_SUPPORT_ENABLED = FALSE)
    } else {
      # Temporarily disable all except distcorr
      saved_cfg <- CFG
      CFG$MTX_LOCALNORM_ENABLED <- FALSE; CFG$MTX_DENOISED_ENABLED <- FALSE
      CFG$MTX_RESIDBG_ENABLED <- FALSE; CFG$MTX_EDGE_ENABLED <- FALSE
      CFG$MTX_SUPPORT_ENABLED <- FALSE
    }
    tryCatch({
      source(mtx_file)
      variants <- generate_all_variants(smat)
      if ("distcorr" %in% names(variants)) {
        distcorr_mat <- variants$distcorr
        message("[TRACKS] Computed distance correction on the fly")
      }
    }, error = function(e) message("[TRACKS] Could not compute distcorr: ", e$message))
  }
}

has_split <- !is.null(distcorr_mat)

if (has_split) {
  # Build lower triangle (distance-corrected)
  lower_rows <- list(); k <- 0L
  for (ii in seq_len(ns)) {
    for (jj in 1:ii) {
      v <- distcorr_mat[idx[ii], idx[jj]]
      if (is.finite(v)) {
        k <- k + 1L
        lower_rows[[k]] <- list(x = pos_mb[jj], y = pos_mb[ii], value = v)
      }
    }
  }
  lower_dt <- rbindlist(lower_rows)

  # Diverging scale for corrected (centered at 0)
  max_abs <- quantile(abs(lower_dt$value), 0.95, na.rm = TRUE)

  p_split <- ggplot() +
    # Upper triangle: raw (warm scale)
    geom_raster(data = upper_dt, aes(x = x, y = y, fill = value), interpolate = TRUE) +
    scale_fill_gradientn(
      colours = SIM_COLORS,
      values = scales::rescale(c(0, sim_med*0.5, sim_med*0.85, sim_med*1.1, sim_med*1.3, 1)),
      name = "Raw sim", limits = c(0, 1))

  # Need ggnewscale for second fill scale
  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    p_split <- p_split +
      ggnewscale::new_scale_fill() +
      geom_raster(data = lower_dt, aes(x = x, y = y, fill = value), interpolate = TRUE) +
      scale_fill_gradientn(
        colours = DIV_COLORS,
        values = scales::rescale(c(-max_abs, -max_abs*0.5, -max_abs*0.1,
                                    0, max_abs*0.1, max_abs*0.5, max_abs)),
        name = "Corrected", limits = c(-max_abs, max_abs))
  } else {
    # Fallback: just plot raw both triangles
    lower_dt_raw <- copy(upper_dt)
    lower_dt_raw[, c("x","y") := .(y, x)]
    p_split <- ggplot(rbind(upper_dt, lower_dt_raw[x != y]),
                       aes(x = x, y = y, fill = value)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_gradientn(
        colours = SIM_COLORS,
        values = scales::rescale(c(0, sim_med*0.5, sim_med*0.85, sim_med*1.1, sim_med*1.3, 1)),
        name = "Similarity", limits = c(0, 1))
    message("[TRACKS] ggnewscale not available — using raw for both triangles")
  }

  p_split <- p_split +
    coord_fixed() +
    labs(x = paste0(chr_label, " (Mb)"), y = paste0(chr_label, " (Mb)"),
         title = paste0(chr_label, " — Split: Raw (upper) vs Distance-corrected (lower)"),
         subtitle = paste0(ns, "×", ns, " | Red/blue = above/below expected at each distance"),
         caption = "Upper: raw similarity | Lower: residual after subtracting distance-expected\nBlue = below expected (boundary) | Red = above expected (block signal)") +
    theme_minimal(base_size = 9) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.caption = element_text(size = 6, hjust = 0, color = "#8A8A8A"))
}


# ============================================================================
# 5. SAVE ALL TRACKS
# ============================================================================

# Individual track plots
for (ext in c("png", "pdf")) {
  tryCatch({
    ggsave(file.path(outdir, paste0(chr_label, "_nesting_depth.", ext)),
           p_depth, width = 14, height = 1.5, dpi = 350)
    ggsave(file.path(outdir, paste0(chr_label, "_nesting_line.", ext)),
           p_nesting_line, width = 14, height = 3, dpi = 350)
    ggsave(file.path(outdir, paste0(chr_label, "_contrast_ratio.", ext)),
           p_contrast, width = 14, height = 3, dpi = 350)
    if (nrow(blocks) > 0)
      ggsave(file.path(outdir, paste0(chr_label, "_cumulative_coverage.", ext)),
             p_cumul, width = 14, height = 3, dpi = 350)
    if (nrow(votes) > 0)
      ggsave(file.path(outdir, paste0(chr_label, "_vote_track.", ext)),
             p_votes, width = 14, height = 3, dpi = 350)
    if (!is.null(p_bg))
      ggsave(file.path(outdir, paste0(chr_label, "_bg_level.", ext)),
             p_bg, width = 14, height = 3, dpi = 350)
    if (has_split)
      ggsave(file.path(outdir, paste0(chr_label, "_split_raw_corrected.", ext)),
             p_split, width = 14, height = 13, dpi = 350)
  }, error = function(e) message("[FAIL] ", ext, ": ", e$message))
}

# Combined multi-panel (stacked tracks)
if (requireNamespace("patchwork", quietly = TRUE) ||
    requireNamespace("gridExtra", quietly = TRUE)) {

  track_list <- list(p_depth, p_nesting_line, p_contrast)
  if (nrow(votes) > 0) track_list <- c(track_list, list(p_votes))
  if (!is.null(p_bg)) track_list <- c(track_list, list(p_bg))
  if (nrow(blocks) > 0) track_list <- c(track_list, list(p_cumul))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    p_combined <- Reduce(`/`, track_list)
    heights <- rep(1, length(track_list))
    heights[1] <- 0.5  # depth bar is short
  } else {
    library(gridExtra)
    p_combined <- do.call(gridExtra::grid.arrange,
                          c(track_list, list(ncol = 1)))
  }

  for (ext in c("png", "pdf")) {
    tryCatch({
      ggsave(file.path(outdir, paste0(chr_label, "_all_tracks.", ext)),
             p_combined, width = 14, height = 2.5 * length(track_list), dpi = 350)
      message("[SAVED] ", chr_label, "_all_tracks.", ext)
    }, error = function(e) message("[FAIL] combined: ", e$message))
  }
}

message("\n[DONE] Marginal tracks -> ", outdir)
