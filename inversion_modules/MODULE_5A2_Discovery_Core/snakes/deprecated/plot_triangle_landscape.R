#!/usr/bin/env Rscript

# =============================================================================
# plot_triangle_landscape.R  (v1)
#
# UNIFIED PLOT SCRIPT for triangle insulation results.
# Reads outputs, produces all figures. No detection logic.
#
# KEY INNOVATION: Detects triangles on NN-smoothed sim_mats (nn20/40/80),
# then projects those boundaries back onto the RAW sim_mat (nn0) for visual
# validation. If a triangle found at nn40 looks like a real block on nn0,
# it's real. If it's diffuse noise on nn0, it's an artifact.
#
# INPUTS:
#   --precomp-dir   : precomp RDS directory (for raw sim_mat + window grid)
#   --triangle-dir  : multiscale triangle output (with nn0/, nn20/, nn40/, nn80/)
#   --bloc-dir      : (optional) PHASE_01C bloc detect output (landscape/)
#   --outdir        : where to write plots
#   --chrom         : single chromosome (default: all)
#
# OUTPUTS (per chromosome):
#   A1: Whole-chr rotated triangle heatmap (nn0) with NN-detected intervals
#   A2: Same but zoomed into top intervals (one per zoom)
#   B1: Multi-scale comparison strip (nn0 vs nn20 vs nn40 vs nn80 side by side)
#   C1: Interval persistence table (which intervals survive across scales)
#   D1: Combined panel (heatmap + strip + insulation + scores)
#
# Usage:
#   Rscript plot_triangle_landscape.R \
#     --precomp-dir precomp/ \
#     --triangle-dir triangles_v3/ \
#     --bloc-dir landscape/ \
#     --outdir plots_landscape/ \
#     --chrom C_gar_LG01
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# PARSE ARGS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
precomp_dir <- NULL; triangle_dir <- NULL; bloc_dir <- NULL
outdir <- "plots_landscape"; chrom_filter <- NULL
ZOOM_TOP <- 10L; BIN_SIZE <- 2L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--precomp-dir" && i < length(args))  { precomp_dir <- args[i+1]; i <- i+2L }
  else if (a == "--triangle-dir" && i < length(args)) { triangle_dir <- args[i+1]; i <- i+2L }
  else if (a == "--bloc-dir" && i < length(args)) { bloc_dir <- args[i+1]; i <- i+2L }
  else if (a == "--outdir" && i < length(args))   { outdir <- args[i+1]; i <- i+2L }
  else if (a == "--chrom" && i < length(args))    { chrom_filter <- args[i+1]; i <- i+2L }
  else if (a == "--zoom-top" && i < length(args)) { ZOOM_TOP <- as.integer(args[i+1]); i <- i+2L }
  else if (a == "--bin" && i < length(args))      { BIN_SIZE <- as.integer(args[i+1]); i <- i+2L }
  else { i <- i+1L }
}

if (is.null(precomp_dir)) stop("Need --precomp-dir")
if (is.null(triangle_dir)) stop("Need --triangle-dir")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "zooms"), recursive = TRUE, showWarnings = FALSE)

# Theme
THEME_FILE <- file.path(Sys.getenv("CODEBASE", ""), "utils", "theme_systems_plate.R")
if (!file.exists(THEME_FILE)) {
  THEME_FILE <- file.path(
    Sys.getenv("BASE", "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"),
    "inversion_codebase_v8.5", "utils", "theme_systems_plate.R")
}
if (file.exists(THEME_FILE)) { source(THEME_FILE); THEME <- theme_plate(base_size = 9)
} else { THEME <- theme_minimal(base_size = 9) }

DPI <- 400

# Scale colors for NN overlay
NN_COLORS <- c(
  "nn0"  = "#78909C",   # grey — raw (for reference only)
  "nn20" = "#FF1744",   # red
  "nn40" = "#FF9100",   # orange
  "nn80" = "#FFEA00"    # yellow
)

# Interval type line styles
TYPE_LTY <- c(
  "strong_triangle" = "solid",
  "moderate_triangle" = "solid",
  "sharp_but_not_square" = "dashed",
  "patchy_signal" = "dotted",
  "diffuse_zone" = "dotted",
  "weak_zone" = "dotted",
  "transition" = "dotted"
)

message("[PLOT] Triangle Landscape Plotter")
message("[PLOT] Precomp:  ", precomp_dir)
message("[PLOT] Triangle: ", triangle_dir)
message("[PLOT] Bloc:     ", bloc_dir %||% "(none)")
message("[PLOT] Output:   ", outdir)

# =============================================================================
# LOAD PRECOMP
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No precomp RDS files")

precomp_list <- list()
for (f in rds_files) { obj <- readRDS(f); precomp_list[[obj$chrom]] <- obj }
chroms <- names(precomp_list)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)
message("[PLOT] ", length(chroms), " chromosomes")

# =============================================================================
# LOAD TRIANGLE RESULTS (all NN scales)
# =============================================================================

nn_dirs <- list.dirs(triangle_dir, recursive = FALSE, full.names = TRUE)
nn_dirs <- nn_dirs[grepl("nn[0-9]+$", basename(nn_dirs))]
nn_dirs <- nn_dirs[order(as.integer(sub(".*nn", "", basename(nn_dirs))))]

all_intervals <- list()
for (nd in nn_dirs) {
  scale_tag <- basename(nd)
  # Find interval files: try chr-prefixed first, then unprefixed
  iv_files <- list.files(nd, pattern = "triangle_intervals\\.tsv\\.gz$", full.names = TRUE)
  if (length(iv_files) == 0) next
  for (ivf in iv_files) {
    iv <- fread(ivf)
    iv[, nn_scale := scale_tag]
    all_intervals[[paste0(scale_tag, "_", basename(ivf))]] <- iv
  }
  message("[PLOT] Loaded ", scale_tag, ": ", sum(sapply(iv_files, function(f) nrow(fread(f, nrows = 0)))), " files, ",
          nrow(iv), " intervals (last file)")
}

if (length(all_intervals) == 0) stop("No triangle interval files found in ", triangle_dir)
all_iv <- rbindlist(all_intervals, fill = TRUE)

# =============================================================================
# LOAD BLOC DETECT BOUNDARIES (optional)
# =============================================================================

bloc_boundaries <- data.table()
bloc_blocks <- data.table()
if (!is.null(bloc_dir) && dir.exists(bloc_dir)) {
  for (chr in chroms) {
    bf <- file.path(bloc_dir, paste0("boundary_catalog_", chr, ".tsv.gz"))
    if (file.exists(bf)) bloc_boundaries <- rbind(bloc_boundaries, fread(bf), fill = TRUE)
    blf <- file.path(bloc_dir, paste0("block_registry_", chr, ".tsv.gz"))
    if (file.exists(blf)) bloc_blocks <- rbind(bloc_blocks, fread(blf), fill = TRUE)
  }
  if (nrow(bloc_boundaries) > 0) message("[PLOT] Bloc boundaries: ", nrow(bloc_boundaries))
  if (nrow(bloc_blocks) > 0) message("[PLOT] Bloc blocks: ", nrow(bloc_blocks))
}

# =============================================================================
# HELPER: Build rotated triangle heatmap data from sim_mat
# =============================================================================

build_triangle_data <- function(sim_mat, dt, bin_size = 2L, max_bins = 500L) {
  n <- nrow(sim_mat)

  # Direct subsampling: pick every step-th window, no binning loop needed
  total_step <- max(1L, n %/% max_bins) * bin_size
  idx <- seq(1, n, by = total_step)
  ns <- length(idx)

  # Subsample sim_mat directly (instant — just matrix indexing)
  sub_mat <- sim_mat[idx, idx]

  # Bin positions
  bin_pos <- (dt$start_bp[idx] + dt$end_bp[pmin(idx + total_step - 1L, n)]) / 2e6
  dx <- median(diff(bin_pos))

  # Build triangle: upper triangle only, cap distance at half chromosome
  max_dist_idx <- ns %/% 2
  rows <- vector("list", ns * min(max_dist_idx, ns))
  ri <- 0L
  for (ii in seq_len(ns)) {
    jj_max <- min(ns, ii + max_dist_idx)
    for (jj in ii:jj_max) {
      ri <- ri + 1L
      rows[[ri]] <- list(
        x_mb = (bin_pos[ii] + bin_pos[jj]) / 2,
        y_dist = bin_pos[jj] - bin_pos[ii],
        sim = sub_mat[ii, jj]
      )
    }
  }
  tri_dt <- rbindlist(rows[seq_len(ri)])

  message("[PLOT]   Triangle data: ", ns, " points (step=", total_step,
          "), ", nrow(tri_dt), " tiles, dx=", round(dx, 3), " Mb")

  list(tri_dt = tri_dt, bin_pos = bin_pos, sub_mat = sub_mat, ns = ns,
       step = total_step, idx = idx, dx = dx)
}

# =============================================================================
# HELPER: Build triangle outlines for interval set
# =============================================================================

build_interval_outlines <- function(intervals, scale_tag) {
  if (nrow(intervals) == 0) return(data.table())
  outlines <- list()
  for (r in seq_len(nrow(intervals))) {
    Lm <- as.numeric(intervals$start_mb[r])
    Rm <- as.numeric(intervals$end_mb[r])
    mid <- (Lm + Rm) / 2
    h <- (Rm - Lm) / 2  # height in Mb
    sq <- as.numeric(intervals$squareness[r])
    itype <- intervals$interval_type[r]
    
    outlines[[r]] <- data.table(
      x = c(Lm, mid, Rm, Lm),
      y = c(0, h, 0, 0),
      grp = paste0(scale_tag, "_I", intervals$interval_id[r]),
      scale = scale_tag,
      interval_id = as.integer(intervals$interval_id[r]),
      squareness = sq,
      interval_type = itype,
      lw = 0.5 + 2.0 * min(sq, 1),
      alpha = 0.3 + 0.6 * min(sq, 1)
    )
  }
  rbindlist(outlines)
}

# =============================================================================
# PLOT A1: Whole-chromosome nn0 heatmap with NN-detected triangles overlaid
# =============================================================================

for (chr in chroms) {
  pc <- precomp_list[[chr]]
  if (is.null(pc) || pc$n_windows < 50) next
  dt <- pc$dt; sim_mat <- pc$sim_mat
  
  message("\n[PLOT] === ", chr, " ===")
  
  # Build triangle heatmap from raw (nn0) sim_mat
  tri_data <- build_triangle_data(sim_mat, dt, BIN_SIZE)
  tri_dt <- tri_data$tri_dt
  
  # Collect interval outlines from ALL NN scales (except nn0)
  chr_iv <- all_iv[chrom == chr]
  nn_scales_present <- sort(unique(chr_iv$nn_scale))
  nn_scales_overlay <- nn_scales_present[nn_scales_present != "nn0"]
  
  # Also keep nn0 for reference (thinner, greyer)
  all_outlines <- list()
  for (stag in nn_scales_present) {
    iv_sub <- chr_iv[nn_scale == stag]
    # For overlay on nn0 heatmap, only show strong/moderate from NN scales
    if (stag != "nn0") {
      iv_sub <- iv_sub[interval_type %in% c("strong_triangle", "moderate_triangle",
                                              "sharp_but_not_square")]
    }
    ol <- build_interval_outlines(iv_sub, stag)
    if (nrow(ol) > 0) all_outlines[[stag]] <- ol
  }
  outlines_dt <- if (length(all_outlines) > 0) rbindlist(all_outlines, fill = TRUE) else data.table()
  
  # Bloc detect boundaries
  chr_bloc_bnd <- bloc_boundaries[chrom == chr & boundary_type %in% c("clear_hard", "clear_soft")]
  
  # ── PLOT A1: whole chromosome ──
  y_max <- max(tri_dt$y_dist, na.rm = TRUE)
  vmax <- quantile(tri_dt$sim[is.finite(tri_dt$sim)], 0.995, na.rm = TRUE)
  dx <- tri_data$dx

  pA1 <- ggplot(tri_dt, aes(x = x_mb, y = y_dist, fill = sim)) +
    geom_tile(width = dx * 1.05, height = dx * 1.05) +
    scale_fill_gradientn(
      colours = c("#080808", "#1a0a2e", "#3d1261", "#6a1b8a",
                  "#9c2e7a", "#cc4466", "#e8734a", "#f5a623"),
      name = "Similarity", na.value = "white",
      limits = c(0, vmax)) +
    scale_y_reverse(limits = c(y_max * 1.02, -y_max * 0.015))
  
  # Add NN triangle outlines
  if (nrow(outlines_dt) > 0) {
    for (stag in unique(outlines_dt$scale)) {
      ol_sub <- outlines_dt[scale == stag]
      col <- NN_COLORS[stag] %||% "white"
      lw_base <- if (stag == "nn0") 0.3 else 0.8
      alpha_base <- if (stag == "nn0") 0.3 else 0.7
      
      pA1 <- pA1 + geom_path(
        data = ol_sub, inherit.aes = FALSE,
        aes(x = x, y = y, group = grp),
        color = col, linewidth = lw_base, alpha = alpha_base)
    }
  }
  
  # Add bloc detect boundaries as vertical dashes
  if (nrow(chr_bloc_bnd) > 0) {
    pA1 <- pA1 + geom_vline(
      xintercept = chr_bloc_bnd$pos_mb,
      color = "#00E5FF", linewidth = 0.3, linetype = "dashed", alpha = 0.5)
  }
  
  # Legend for NN scales
  nn_legend_items <- list()
  for (stag in nn_scales_present) {
    col <- NN_COLORS[stag] %||% "grey"
    nn_legend_items[[stag]] <- annotate("segment", x = -Inf, xend = -Inf,
                                         y = -Inf, yend = -Inf,
                                         color = col, linewidth = 1.5)
  }
  
  # Manual legend via dummy data
  legend_dt <- data.table(
    x = rep(min(tri_dt$x_mb), length(nn_scales_present)),
    y = rep(-1, length(nn_scales_present)),
    scale = nn_scales_present
  )
  
  pA1 <- pA1 +
    geom_path(data = legend_dt, inherit.aes = FALSE,
              aes(x = x, y = y, color = scale), linewidth = 0, alpha = 0) +
    scale_color_manual(values = NN_COLORS, name = "Detection scale") +
    labs(title = paste0(chr, " -- Triangle Landscape (NN-detected on nn0 heatmap)"),
         subtitle = paste0(
           "nn0: ", sum(chr_iv$nn_scale == "nn0"), " intervals | ",
           paste(sapply(nn_scales_overlay, function(s) {
             paste0(s, ": ", sum(chr_iv$nn_scale == s &
                                   chr_iv$interval_type %in% c("strong_triangle", "moderate_triangle", "sharp_but_not_square")))
           }), collapse = " | "),
           if (nrow(chr_bloc_bnd) > 0) paste0(" | bloc boundaries: ", nrow(chr_bloc_bnd)) else ""
         ),
         x = paste0(chr, " position (Mb)"), y = "Pairwise distance (Mb)") +
    THEME +
    theme(plot.title = element_text(size = 11, face = "bold"),
          legend.position = "right")
  
  f_A1 <- file.path(outdir, paste0(chr, "_A1_landscape_wholechr.png"))
  tryCatch(ggsave(f_A1, pA1, width = 16, height = 8, dpi = DPI),
           error = function(e) message("  [PLOT] ", e$message))
  message("[PLOT] ", f_A1)
  
  # ── PLOT A2: zoomed views for top NN-detected intervals ──
  # Take top intervals by squareness from nn20+ (not nn0)
  nn_iv <- chr_iv[nn_scale != "nn0" &
                    interval_type %in% c("strong_triangle", "moderate_triangle")]
  if (nrow(nn_iv) > 0) {
    nn_iv[, squareness := as.numeric(squareness)]
    nn_iv <- nn_iv[order(-squareness)]
    # Deduplicate overlapping intervals across scales (keep highest scale)
    nn_iv[, start_mb := as.numeric(start_mb)]
    nn_iv[, end_mb := as.numeric(end_mb)]
    
    zoom_count <- 0L
    zoomed_regions <- list()
    
    for (zi in seq_len(nrow(nn_iv))) {
      if (zoom_count >= ZOOM_TOP) break
      iv <- nn_iv[zi]
      
      # Check overlap with already-zoomed regions
      overlaps <- FALSE
      for (zr in zoomed_regions) {
        if (iv$start_mb < zr[2] && iv$end_mb > zr[1]) { overlaps <- TRUE; break }
      }
      if (overlaps) next
      
      pad <- max(1, (iv$end_mb - iv$start_mb) * 0.5)
      r_start <- iv$start_mb - pad
      r_end <- iv$end_mb + pad
      
      # Subset triangle data for zoom region
      zoom_dt <- tri_dt[x_mb >= r_start & x_mb <= r_end & y_dist <= (r_end - r_start)]
      
      if (nrow(zoom_dt) < 10) next
      
      # Subset outlines
      zoom_ol <- outlines_dt[
        sapply(seq_len(nrow(outlines_dt)), function(r) {
          outlines_dt$x[r] >= r_start * 0.9 & outlines_dt$x[r] <= r_end * 1.1
        })]
      
      vmax_z <- quantile(zoom_dt$sim[is.finite(zoom_dt$sim)], 0.995, na.rm = TRUE)
      y_max_z <- max(zoom_dt$y_dist, na.rm = TRUE)
      
      pZ <- ggplot(zoom_dt, aes(x = x_mb, y = y_dist, fill = sim)) +
        geom_tile(width = dx * 1.05, height = dx * 1.05) +
        scale_fill_gradientn(
          colours = c("#080808", "#1a0a2e", "#3d1261", "#6a1b8a",
                      "#9c2e7a", "#cc4466", "#e8734a", "#f5a623"),
          name = "Similarity", limits = c(0, vmax_z)) +
        scale_y_reverse(limits = c(y_max_z * 1.02, -y_max_z * 0.02))
      
      # Overlays
      if (nrow(zoom_ol) > 0) {
        for (stag in unique(zoom_ol$scale)) {
          ol_s <- zoom_ol[scale == stag]
          col <- NN_COLORS[stag] %||% "white"
          lw <- if (stag == "nn0") 0.3 else 1.0
          al <- if (stag == "nn0") 0.3 else 0.8
          pZ <- pZ + geom_path(data = ol_s, inherit.aes = FALSE,
                                aes(x = x, y = y, group = grp),
                                color = col, linewidth = lw, alpha = al)
        }
      }
      
      # Bloc boundaries
      zoom_bloc <- chr_bloc_bnd[pos_mb >= r_start & pos_mb <= r_end]
      if (nrow(zoom_bloc) > 0) {
        pZ <- pZ + geom_vline(xintercept = zoom_bloc$pos_mb,
                               color = "#00E5FF", linewidth = 0.4, linetype = "dashed", alpha = 0.6)
      }
      
      pZ <- pZ +
        labs(title = paste0("I", iv$interval_id, " -- ", iv$interval_type,
                            " (sq=", round(iv$squareness, 2), ") -- ", chr),
             subtitle = paste0("Detected at ", iv$nn_scale, " | ",
                               round(iv$start_mb, 2), "-", round(iv$end_mb, 2), " Mb"),
             x = paste0(chr, " (Mb)"), y = "Pairwise distance (Mb)") +
        THEME + theme(plot.title = element_text(size = 11, face = "bold"))
      
      f_zoom <- file.path(outdir, "zooms",
                           paste0(chr, "_zoom_", iv$nn_scale, "_I", iv$interval_id,
                                  "_sq", round(iv$squareness, 2), ".png"))
      tryCatch(ggsave(f_zoom, pZ, width = 12, height = 6, dpi = DPI),
               error = function(e) message("  [PLOT] ", e$message))
      
      zoom_count <- zoom_count + 1L
      zoomed_regions[[zoom_count]] <- c(iv$start_mb, iv$end_mb)
    }
    message("[PLOT] ", zoom_count, " zoomed views")
  }
  
  # ── PLOT C1: Scale persistence table ──
  # Which intervals persist across NN scales?
  if (length(nn_scales_overlay) >= 2) {
    # Build overlap matrix: for each nn0 interval, check if it overlaps with
    # an interval at each higher NN scale
    nn0_iv <- chr_iv[nn_scale == "nn0"]
    if (nrow(nn0_iv) > 0) {
      persist_rows <- list()
      for (r in seq_len(nrow(nn0_iv))) {
        iv0 <- nn0_iv[r]
        s0 <- as.numeric(iv0$start_mb); e0 <- as.numeric(iv0$end_mb)
        row_data <- list(
          interval_id = as.integer(iv0$interval_id),
          type_nn0 = iv0$interval_type,
          sq_nn0 = round(as.numeric(iv0$squareness), 3),
          start_mb = s0, end_mb = e0,
          size_mb = round(e0 - s0, 3)
        )
        
        for (stag in nn_scales_overlay) {
          nn_sub <- chr_iv[nn_scale == stag]
          # Find best overlapping interval
          overlap_frac <- 0; best_sq <- 0; best_type <- "none"
          for (ri in seq_len(nrow(nn_sub))) {
            sn <- as.numeric(nn_sub$start_mb[ri]); en <- as.numeric(nn_sub$end_mb[ri])
            ol <- max(0, min(e0, en) - max(s0, sn))
            span <- max(e0 - s0, en - sn)
            frac <- if (span > 0) ol / span else 0
            if (frac > overlap_frac) {
              overlap_frac <- frac
              best_sq <- as.numeric(nn_sub$squareness[ri])
              best_type <- nn_sub$interval_type[ri]
            }
          }
          row_data[[paste0("overlap_", stag)]] <- round(overlap_frac, 2)
          row_data[[paste0("sq_", stag)]] <- round(best_sq, 3)
          row_data[[paste0("type_", stag)]] <- best_type
        }
        
        # Persistence score: how many NN scales have >50% overlap?
        n_persist <- sum(sapply(nn_scales_overlay, function(s) {
          row_data[[paste0("overlap_", s)]] > 0.3
        }))
        row_data$n_scales_persist <- n_persist
        
        persist_rows[[r]] <- as.data.table(row_data)
      }
      
      persist_dt <- rbindlist(persist_rows, fill = TRUE)
      persist_dt <- persist_dt[order(-n_scales_persist, -sq_nn0)]
      
      f_persist <- file.path(outdir, paste0(chr, "_C1_scale_persistence.tsv"))
      fwrite(persist_dt, f_persist, sep = "\t")
      message("[PLOT] Persistence table: ", f_persist)
      
      # Summary
      message("[PLOT] Scale persistence for ", chr, ":")
      message("  Persist 0 scales: ", sum(persist_dt$n_scales_persist == 0),
              " (noise at nn0)")
      for (ns in seq_along(nn_scales_overlay)) {
        message("  Persist ", ns, " scales: ",
                sum(persist_dt$n_scales_persist == ns))
      }
      message("  Persist ALL (", length(nn_scales_overlay), "): ",
              sum(persist_dt$n_scales_persist == length(nn_scales_overlay)),
              " (most confident)")
    }
  }
}

message("\n[PLOT] Done. Output: ", outdir)
