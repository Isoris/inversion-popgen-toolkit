#!/usr/bin/env Rscript
# =============================================================================
# STEP_D17c_overlay_plot.R
#
# Heatmap with D17b catalogue overlay. Per chromosome, produces a tiled PDF
# (5 tiles by default) showing:
#   - sim_mat at the chosen NN smoothing scale (default nn80)
#   - D17b candidate boundaries drawn as outlines, with line width reduced
#     by 2x compared to the previous D13d-style overlay, and color coded by
#     recursion_depth (0 = atomic D17, 1 = first sub-block split, 2 = nested
#     sub-block of a sub-block)
#   - ggrepel labels with candidate_id and Mb range
#   - status-aware styling: PASS solid, PASS_CHILD solid, PASS_PARENT dashed,
#     PASS_COMPLEX dotted, FAIL/FRAGMENTED greyed out
#
# Inputs:
#   --precomp     <slim precomp .rds for window coordinates>
#   --sim_mat     <smoothed sim_mat .rds, default tries sim_mat_nn80.rds>
#   --catalogue   <D17b catalogue TSV>
#   --chr         <chromosome label>
#   --outdir      <output directory>
#   --n_tiles     <number of tiles to split the chromosome into; default 5>
#
# Output:
#   <outdir>/<chr>_d17c_overlay.pdf
#
# Author: Claude (Anthropic)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)   # for layered fill scales (sim + Z)
  library(scales)
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ---- CLI --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

precomp_f    <- get_arg("--precomp")
sim_mat_f    <- get_arg("--sim_mat")
catalogue_f  <- get_arg("--catalogue")
chr_label    <- get_arg("--chr", "chr")
outdir       <- get_arg("--outdir", ".")
n_tiles      <- as.integer(get_arg("--n_tiles", "5"))
nn_label     <- get_arg("--nn_label", "nn80")
page_w       <- as.numeric(get_arg("--page_width_in", "12"))
page_h       <- as.numeric(get_arg("--page_height_in", "12"))

# Outline line width — reduced by 2x compared to typical 0.6
outline_lw   <- as.numeric(get_arg("--outline_lw", "0.30"))
label_size   <- as.numeric(get_arg("--label_size", "2.6"))

# Split-heatmap Z mode: lower triangle = raw similarity, upper triangle = Z.
# Z mode controls how Z is computed per pixel:
#   "diagonal" (default): Z relative to pixels at same distance from diagonal.
#                          Best for separating real inversions from proximity-LD.
#   "global"            : Z relative to all pixels in the chromosome.
#   "row"               : Z relative to other pixels in the same row.
#   "off"               : disable Z; render full heatmap as raw similarity.
z_mode       <- get_arg("--z_mode", "diagonal")
z_clip       <- as.numeric(get_arg("--z_clip", "5"))   # clip Z to +/- this
stopifnot(z_mode %in% c("diagonal", "global", "row", "off"))

if (is.na(precomp_f)   || !file.exists(precomp_f))   stop("[D17c] --precomp required")
if (is.na(catalogue_f) || !file.exists(catalogue_f)) stop("[D17c] --catalogue required")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[D17c] precomp:    ", precomp_f, "\n")
cat("[D17c] sim_mat:    ", sim_mat_f %||% "(auto)", "\n")
cat("[D17c] catalogue:  ", catalogue_f, "\n")
cat("[D17c] chr:        ", chr_label, "\n")
cat("[D17c] outdir:     ", outdir, "\n")
cat("[D17c] n_tiles:    ", n_tiles, "\n")
cat("[D17c] outline_lw: ", outline_lw, "\n")

# ---- Load ------------------------------------------------------------------

cat("[D17c] loading precomp\n")
pc_obj <- readRDS(precomp_f)
if (!is.null(pc_obj$dt)) {
  pc <- pc_obj
} else if (!is.null(pc_obj$pc)) {
  pc <- pc_obj$pc
} else stop("[D17c] cannot find pc$dt in precomp file")
dt_pc <- as.data.table(pc$dt)
n_windows_total <- nrow(dt_pc)
window_start_bp <- dt_pc$start
window_end_bp   <- dt_pc$end
cat("[D17c] N windows: ", n_windows_total, "\n")

if (is.na(sim_mat_f)) {
  candidates_path <- c(
    file.path(dirname(precomp_f), "sim_mat_nn80.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn160.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn40.rds")
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0) stop("[D17c] could not auto-find sim_mat_nn*.rds")
  sim_mat_f <- hit[1]
}
cat("[D17c] loading sim_mat: ", sim_mat_f, "\n")
sm_obj <- readRDS(sim_mat_f)
if (is.matrix(sm_obj)) {
  sim_mat <- sm_obj
} else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat)) {
  sim_mat <- sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) == 1) {
  sim_mat <- sm_obj[[1]]
} else stop("[D17c] sim_mat structure not recognized")
storage.mode(sim_mat) <- "double"

cat("[D17c] loading catalogue\n")
cat_dt <- fread(catalogue_f)
cat("[D17c] catalogue rows: ", nrow(cat_dt), "\n", sep = "")

# Restrict to PASS/PASS_CHILD/PASS_PARENT/PASS_COMPLEX/RESCUED for display.
# (FAIL/FRAGMENTED are greyed out optionally.)
display_status <- c("PASS", "PASS_CHILD", "PASS_PARENT", "PASS_COMPLEX",
                    "RESCUED", "FRAGMENTED", "FAIL")
cat_dt <- cat_dt[status %in% display_status]
cat_dt <- cat_dt[!is.na(start_w) & !is.na(end_w)]

# ---- Z-score matrix (precompute once, used across all tiles) ----------------
# When z_mode != "off", the upper triangle of each tile renders Z instead of
# raw similarity. The Z matrix is computed from the FULL chromosome sim_mat
# so the values are comparable across tiles.

compute_z_matrix <- function(M, mode) {
  N <- nrow(M)
  if (mode == "off") {
    return(NULL)
  } else if (mode == "global") {
    mu <- mean(M, na.rm = TRUE)
    sg <- sd(as.numeric(M), na.rm = TRUE)
    if (!is.finite(sg) || sg < 1e-9) sg <- 1
    return((M - mu) / sg)
  } else if (mode == "row") {
    row_mu <- rowMeans(M, na.rm = TRUE)
    row_sd <- apply(M, 1, sd, na.rm = TRUE)
    row_sd[!is.finite(row_sd) | row_sd < 1e-9] <- 1
    Z <- sweep(M, 1, row_mu, "-")
    Z <- sweep(Z, 1, row_sd, "/")
    return(Z)
  } else if (mode == "diagonal") {
    # For each off-diagonal distance d, compute mu and sd over all pixels at
    # |i - j| == d, then Z each such pixel against (mu_d, sd_d).
    # Use direct integer indexing rather than a full logical mask to keep
    # memory allocation small.
    cat("[D17c] computing diagonal-distance Z\n")
    Z <- matrix(0, nrow = N, ncol = N)
    for (d in 0:(N - 1L)) {
      idx_i <- seq_len(N - d)
      idx_j <- idx_i + d
      idx_upper <- cbind(idx_i, idx_j)
      idx_lower <- cbind(idx_j, idx_i)
      vals_upper <- M[idx_upper]
      vals_lower <- M[idx_lower]
      vals <- if (d == 0L) vals_upper else c(vals_upper, vals_lower)
      mu_d <- mean(vals, na.rm = TRUE)
      sd_d <- sd(vals, na.rm = TRUE)
      if (!is.finite(sd_d) || sd_d < 1e-9) sd_d <- 1
      Z[idx_upper] <- (vals_upper - mu_d) / sd_d
      if (d > 0L) Z[idx_lower] <- (vals_lower - mu_d) / sd_d
      if (d > 0L && d %% 1000L == 0L) {
        cat(sprintf("[D17c]   diag d=%d / %d\n", d, N - 1L))
      }
    }
    return(Z)
  }
  stop("[D17c] unknown z_mode: ", mode)
}

z_mat <- NULL
if (z_mode != "off") {
  cat("[D17c] computing Z matrix in mode '", z_mode, "'\n", sep = "")
  z_mat <- compute_z_matrix(sim_mat, z_mode)
  zr <- range(z_mat, na.rm = TRUE)
  cat(sprintf("[D17c] Z range: [%.2f, %.2f]; clipping at +/- %.1f\n",
              zr[1], zr[2], z_clip))
  # Clip extreme values for stable color rendering
  z_mat[z_mat >  z_clip] <-  z_clip
  z_mat[z_mat < -z_clip] <- -z_clip
}

# ---- Color & line styling --------------------------------------------------

# Color by discovery_pass: 1 large=deep blue, 2 medium=orange, 3 small=purple.
# Children of clustering keep their parent's pass color (inherited).
depth_palette <- c(
  "1" = "#1F3A6E",   # deep blue (pass 1: large W={800,400})
  "2" = "#E07B3F",   # orange    (pass 2: medium W={200,100})
  "3" = "#7A4FBF"    # purple    (pass 3: small W={50})
)

# Line type by status
status_linetype <- c(
  "PASS"         = "solid",
  "PASS_CHILD"   = "solid",
  "PASS_PARENT"  = "dashed",
  "PASS_COMPLEX" = "dotted",
  "RESCUED"      = "solid",
  "FRAGMENTED"   = "11",       # tight dots
  "FAIL"         = "11"
)

# Alpha by status (greyed-out for fail/fragmented)
status_alpha <- function(s) {
  ifelse(s %in% c("FRAGMENTED", "FAIL"), 0.35, 0.95)
}

# ---- Tile setup ------------------------------------------------------------

windows_per_tile <- ceiling(n_windows_total / n_tiles)
tiles <- lapply(seq_len(n_tiles), function(t) {
  s <- (t - 1L) * windows_per_tile + 1L
  e <- min(t * windows_per_tile, n_windows_total)
  list(tile_idx = t, w_start = s, w_end = e)
})

# ---- Build long-form heatmap data per tile and plot ------------------------

build_heatmap_long <- function(w_start, w_end) {
  W <- w_end - w_start + 1L
  sub_sim <- sim_mat[w_start:w_end, w_start:w_end]
  ii <- rep(seq_len(W), times = W)
  jj <- rep(seq_len(W), each = W)
  i_global <- ii + w_start - 1L
  j_global <- jj + w_start - 1L

  if (z_mode == "off" || is.null(z_mat)) {
    # No Z mode — render whole matrix as raw similarity in lower-triangle scale
    return(data.table(
      i = i_global,
      j = j_global,
      tri = "sim",
      sim_val = as.numeric(sub_sim),
      z_val = NA_real_
    ))
  }

  sub_z <- z_mat[w_start:w_end, w_start:w_end]
  # Lower triangle (j <= i): raw sim
  # Upper triangle (j >  i): Z
  lower_mask <- jj <= ii
  upper_mask <- jj >  ii

  out <- data.table(
    i = i_global,
    j = j_global,
    tri = ifelse(lower_mask, "sim", "z"),
    sim_val = ifelse(lower_mask, as.numeric(sub_sim), NA_real_),
    z_val   = ifelse(upper_mask, as.numeric(sub_z),   NA_real_)
  )
  out
}

candidates_in_tile <- function(w_start, w_end) {
  # Keep candidates whose window range overlaps the tile
  dt <- cat_dt[!(end_w < w_start | start_w > w_end)]
  if (nrow(dt) == 0L) return(dt)
  # Clip to tile
  dt[, draw_start := pmax(start_w, w_start)]
  dt[, draw_end   := pmin(end_w,   w_end)]
  dt[, mid_w := (draw_start + draw_end) / 2]
  dt[, mid_mb := (window_start_bp[round(mid_w)] +
                  window_end_bp[round(mid_w)]) / 2 / 1e6]
  dt[is.na(pattern_class), pattern_class := "NA"]
  dt[, label_text := sprintf("%s [%s] %.2f-%.2f Mb",
                             candidate_id,
                             pattern_class,
                             start_bp / 1e6, end_bp / 1e6)]
  # Color discriminator: discovery_pass (1/2/3). Children inherit parent's pass.
  if ("discovery_pass" %in% names(dt)) {
    dt[, depth_chr := as.character(discovery_pass)]
  } else if ("recursion_depth" %in% names(dt)) {
    # Backward compat: if running on old D17b output with no discovery_pass,
    # fall back to recursion_depth coloring (mapped 0->1, 1->2, 2->3).
    dt[, depth_chr := as.character(pmin(3L, as.integer(recursion_depth) + 1L))]
  } else {
    dt[, depth_chr := "1"]
  }
  # NA-safe
  dt[is.na(depth_chr) | depth_chr == "NA", depth_chr := "1"]
  dt
}

cat("[D17c] building plots for ", n_tiles, " tiles\n")

plots <- list()
for (tile in tiles) {
  w_start <- tile$w_start
  w_end   <- tile$w_end
  W <- w_end - w_start + 1L
  cat(sprintf("[D17c]   tile %d/%d: w=%d-%d (%.2f-%.2f Mb)\n",
              tile$tile_idx, n_tiles, w_start, w_end,
              window_start_bp[w_start] / 1e6,
              window_end_bp[w_end] / 1e6))

  hm <- build_heatmap_long(w_start, w_end)
  cands <- candidates_in_tile(w_start, w_end)

  # Compute per-tile color scale quantiles for robust contrast (sim layer)
  sim_layer <- hm[tri == "sim"]
  q_lo <- as.numeric(quantile(sim_layer$sim_val, 0.05, na.rm = TRUE))
  q_hi <- as.numeric(quantile(sim_layer$sim_val, 0.98, na.rm = TRUE))

  if (z_mode == "off" || is.null(z_mat)) {
    # Single-layer rendering: full matrix as similarity
    p <- ggplot(hm, aes(x = i, y = j, fill = sim_val)) +
      geom_raster() +
      scale_fill_gradientn(
        colours = c("#F8F8F8", "#A8DBC2", "#F2DC78", "#E08838", "#7E1F1F"),
        values = scales::rescale(c(0, q_lo, (q_lo + q_hi)/2, q_hi, 1)),
        limits = c(0, 1),
        oob = scales::squish,
        name = "Similarity"
      )
  } else {
    # Two-layer rendering: lower triangle = sim, upper triangle = Z
    z_layer <- hm[tri == "z"]
    # Z color: diverging blue-white-red, symmetric range
    z_max <- max(abs(z_layer$z_val), na.rm = TRUE)
    z_max <- min(z_max, z_clip)
    if (!is.finite(z_max) || z_max < 1e-3) z_max <- 1

    p <- ggplot() +
      # Lower triangle: similarity
      geom_raster(
        data = sim_layer,
        aes(x = i, y = j, fill = sim_val)
      ) +
      scale_fill_gradientn(
        colours = c("#F8F8F8", "#A8DBC2", "#F2DC78", "#E08838", "#7E1F1F"),
        values = scales::rescale(c(0, q_lo, (q_lo + q_hi)/2, q_hi, 1)),
        limits = c(0, 1),
        oob = scales::squish,
        name = "Similarity\n(lower △)"
      ) +
      # Open a new fill scale for the Z layer
      ggnewscale::new_scale_fill() +
      geom_raster(
        data = z_layer,
        aes(x = i, y = j, fill = z_val)
      ) +
      scale_fill_gradient2(
        low = "#2C5AA0", mid = "#FAFAFA", high = "#B22222",
        midpoint = 0,
        limits = c(-z_max, z_max),
        oob = scales::squish,
        name = sprintf("Z (%s)\n(upper △)", z_mode)
      )
  }

  # Safe secondary Mb axis:
  # ggplot2 sec_axis() requires a strictly monotonic transform.
  # So we keep the transform as identity (~ .) and only convert labels to Mb.
  axis_breaks <- pretty(c(w_start, w_end), n = 6)
  axis_breaks <- unique(as.integer(round(axis_breaks)))
  axis_breaks <- axis_breaks[axis_breaks >= w_start & axis_breaks <= w_end]
  axis_labels <- sprintf("%.2f",
    window_start_bp[pmin(pmax(axis_breaks, 1L), n_windows_total)] / 1e6
  )

  p <- p +
    coord_equal() +
    scale_x_continuous(
      expand = c(0, 0),
      name = "window index",
      breaks = axis_breaks,
      sec.axis = sec_axis(
        ~ .,
        breaks = axis_breaks,
        labels = axis_labels,
        name = "Mb"
      )
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      name = "window index",
      breaks = axis_breaks,
      sec.axis = sec_axis(
        ~ .,
        breaks = axis_breaks,
        labels = axis_labels,
        name = "Mb"
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.2, "cm"),
      plot.title = element_text(size = 11, face = "bold")
    ) +
    labs(title = sprintf("%s | %s | tile %d/%d | w %d-%d (%.2f-%.2f Mb)",
                         chr_label, nn_label, tile$tile_idx, n_tiles,
                         w_start, w_end,
                         window_start_bp[w_start] / 1e6,
                         window_end_bp[w_end] / 1e6))

  if (nrow(cands) > 0L) {
    # Outline rectangles, color by depth, line type by status, thinner
    rect_dt <- copy(cands)
    rect_dt[, alpha_v := status_alpha(status)]

    p <- p +
      geom_rect(
        data = rect_dt,
        aes(xmin = draw_start - 0.5, xmax = draw_end + 0.5,
            ymin = draw_start - 0.5, ymax = draw_end + 0.5,
            colour = depth_chr, linetype = status, alpha = alpha_v),
        fill = NA,
        linewidth = outline_lw,
        inherit.aes = FALSE
      ) +
      scale_colour_manual(
        values = depth_palette,
        breaks = c("1", "2", "3"),
        labels = c("Pass 1 (large W=800/400)",
                   "Pass 2 (medium W=200/100)",
                   "Pass 3 (small W=50)"),
        name = "Discovery pass",
        drop = FALSE
      ) +
      scale_linetype_manual(
        values = status_linetype,
        breaks = names(status_linetype),
        name = "Status",
        drop = FALSE
      ) +
      scale_alpha_identity() +
      guides(
        colour = guide_legend(override.aes = list(linewidth = outline_lw * 2)),
        linetype = guide_legend(override.aes = list(linewidth = outline_lw * 2))
      )

    # Labels near the top-right corner of each rectangle, ggrepel
    label_dt <- copy(cands)
    label_dt[, lab_x := (draw_start + draw_end) / 2]
    label_dt[, lab_y := draw_end + (W * 0.02)]
    p <- p +
      ggrepel::geom_text_repel(
        data = label_dt,
        aes(x = lab_x, y = lab_y, label = label_text, colour = depth_chr),
        size = label_size,
        fontface = "bold",
        bg.colour = "white",
        bg.r = 0.10,
        segment.size = outline_lw,
        segment.alpha = 0.5,
        min.segment.length = 0,
        seed = 1,
        max.overlaps = Inf,
        force = 1.2,
        box.padding = 0.18,
        point.padding = 0.05,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

  plots[[tile$tile_idx]] <- p
}

# ---- Write multi-page PDF --------------------------------------------------

outpath <- file.path(outdir, paste0(chr_label, "_d17c_overlay.pdf"))
cat("[D17c] writing PDF: ", outpath, "\n")

pdf(outpath, width = page_w, height = page_h)
for (p in plots) print(p)
dev.off()

cat("[D17c] done\n")
