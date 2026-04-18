#!/usr/bin/env Rscript
# =============================================================================
# C01a_diag_matrices.R — Similarity matrix heatmaps and transforms
# (plots 06, 14, 15, 16, 17 + window-index + distance-normalized variants)
#
# v8.5.2 — Plot overhaul:
#   G1: sys.frame → env var sourcing
#   PC2/PC3: 3 variants (Mb, window-index, distance-normalized)
#   PC6: multi-scale sim_mat (+80, +160)
#   PC11: sim_mat with marginal density
#   Higher DPI, increased color saturation
#
# Usage: Rscript C01a_diag_matrices.R <precomp_dir> <outdir> [chrom]
# =============================================================================

SCRIPT_DIR <- Sys.getenv("SNAKES_DIR",
  file.path(Sys.getenv("BASE", "."),
            "inversion_codebase_v8.5/MODULE_5A2_Discovery_Core/snakes"))
source(file.path(SCRIPT_DIR, "STEP_C01a_diag_common.R"))

cli <- parse_diag_args()
data <- load_diag_data(cli$precomp_dir, cli$chrom_filter)
precomp_list <- data$precomp_list
chroms <- data$chroms
inv_like_dt <- data$inv_like_dt
dir.create(cli$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPER: Build sim_mat tile data
# =============================================================================

# axis_mode: "mb" (default), "window_index", "distance_normalized"
build_simmat_tiledata <- function(sim_mat, dt, step, idx, axis_mode = "mb") {
  n_idx <- length(idx)
  if (axis_mode == "mb") {
    coords <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6
  } else if (axis_mode == "window_index") {
    coords <- seq_along(idx)
  } else {
    coords <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6
  }

  pdt <- data.table(
    x = rep(coords, each = n_idx),
    y = rep(coords, times = n_idx),
    sim = as.vector(sim_mat[idx, idx])
  )

  if (axis_mode == "distance_normalized") {
    pdt[, dist := abs(x - y)]
    # Normalize by genomic distance: subtract expected sim at that distance
    pdt[, dist_bin := cut(dist, breaks = 50)]
    dist_means <- pdt[is.finite(sim), .(mean_sim = mean(sim, na.rm = TRUE)), by = dist_bin]
    pdt <- merge(pdt, dist_means, by = "dist_bin", all.x = TRUE)
    pdt[, sim := sim - mean_sim]
    pdt[, c("dist", "dist_bin", "mean_sim") := NULL]
  }

  pdt[is.finite(sim)]
}

# =============================================================================
# MATRIX BUILDERS
# =============================================================================

# -- 06: Full similarity matrix heatmap (PC11: marginal density option) ------
build_sim_heatmap <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 10) return(NULL)

  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  pdt <- build_simmat_tiledata(sim_mat, dt, step, idx, "mb")
  sim_med <- median(pdt$sim, na.rm = TRUE)

  p <- ggplot(pdt, aes(x = x, y = y, fill = sim)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = SIM_COLS,
      values = scales::rescale(c(0, sim_med * 0.6, sim_med, sim_med * 1.1, sim_med * 1.3, 1)),
      name = "Similarity") +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " \u2014 Similarity Matrix (", n, " windows)"),
         subtitle = paste0(length(idx), "\u00d7", length(idx), " subsampled | median sim = ", round(sim_med, 3)),
         caption = "Red blocks on diagonal = shared PCA structure (inversion candidates)\nBlue = low similarity (background or boundary)") +
    coord_fixed() + THEME_BASE
  p
}

# -- 06b: Sim_mat by window index (no gaps) ----------------------------------
build_sim_heatmap_idx <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 10) return(NULL)

  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  pdt <- build_simmat_tiledata(sim_mat, dt, step, idx, "window_index")
  sim_med <- median(pdt$sim, na.rm = TRUE)

  p <- ggplot(pdt, aes(x = x, y = y, fill = sim)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = SIM_COLS,
      values = scales::rescale(c(0, sim_med * 0.6, sim_med, sim_med * 1.1, sim_med * 1.3, 1)),
      name = "Similarity") +
    labs(x = "Window index", y = "Window index",
         title = paste0(chr, " \u2014 Similarity Matrix (window index, no gaps)"),
         subtitle = paste0(length(idx), "\u00d7", length(idx), " | removes SNP density distortion"),
         caption = "Window index view eliminates genomic-distance distortion") +
    coord_fixed() + THEME_BASE
  p
}

# -- 14: Local similarity heatmap (band around diagonal, PC6: multi-scale) ----
build_local_sim_heatmap <- function(chr, bandwidth = 80L) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  if (is.null(sim_mat) || nrow(sim_mat) < 20) return(NULL)

  n_w <- nrow(sim_mat)
  step <- max(1L, n_w %/% 1000L)
  idx <- seq(1, n_w, by = step)

  rows <- list()
  for (i in idx) {
    j_min <- max(1L, i - bandwidth)
    j_max <- min(n_w, i + bandwidth)
    j_range <- seq(j_min, j_max, by = max(1L, step %/% 2L))
    for (j in j_range) {
      s <- sim_mat[i, j]
      if (is.finite(s)) {
        rows[[length(rows) + 1]] <- data.table(
          x_mb = (dt$start_bp[i] + dt$end_bp[i]) / 2e6,
          y_mb = (dt$start_bp[j] + dt$end_bp[j]) / 2e6,
          sim = s
        )
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  pdt <- rbindlist(rows)
  med_sim <- median(pdt$sim, na.rm = TRUE)

  # Increased saturation
  p <- ggplot(pdt, aes(x = x_mb, y = y_mb, fill = sim)) +
    geom_tile(width = step * (dt$end_bp[2] - dt$start_bp[1]) / 1e6,
              height = step * (dt$end_bp[2] - dt$start_bp[1]) / 1e6) +
    scale_fill_gradientn(
      colours = c("#0C1E3C", "#2563EB", "#F8F9FA", "#E8913A", "#C0392B"),
      values = scales::rescale(c(0, med_sim * 0.7, med_sim, med_sim * 1.15, 1)),
      name = "Similarity"
    ) +
    coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " \u2014 Local Similarity (\u00b1", bandwidth, " windows)"),
         subtitle = paste0("Median sim = ", round(med_sim, 3)),
         caption = "White = median | Red = high (block) | Blue = low (gap)") +
    THEME_BASE
  p
}

# PC6: +160 variant
build_local_sim_heatmap_160 <- function(chr) {
  build_local_sim_heatmap(chr, bandwidth = 160L)
}

# -- 16: Contrast sim_mat (row Z-score) — PC3: 3 variants -------------------
build_simmat_zscore_core <- function(chr, axis_mode = "mb") {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 50) return(NULL)

  z_mat <- sim_mat
  for (ri in seq_len(n)) {
    row_vals <- sim_mat[ri, ]
    rm <- mean(row_vals, na.rm = TRUE)
    rsd <- sd(row_vals, na.rm = TRUE)
    if (is.finite(rsd) && rsd > 1e-8) z_mat[ri, ] <- (row_vals - rm) / rsd
    else z_mat[ri, ] <- 0
  }
  z_mat <- (z_mat + t(z_mat)) / 2

  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  sub <- z_mat[idx, idx]

  if (axis_mode == "mb") {
    coords <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6
    x_lab <- paste0(chr, " (Mb)"); y_lab <- x_lab
  } else {
    coords <- seq_along(idx)
    x_lab <- "Window index"; y_lab <- x_lab
  }

  pdt <- data.table(
    x = rep(coords, each = length(idx)),
    y = rep(coords, times = length(idx)),
    z = as.vector(sub)
  )
  pdt <- pdt[is.finite(z)]
  pdt[, z_clipped := pmin(pmax(z, -3), 3)]

  mode_label <- if (axis_mode == "mb") "Mb position" else "window index"

  p <- ggplot(pdt, aes(x = x, y = y, fill = z_clipped)) +
    geom_tile() +
    scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH,
                         midpoint = 0, limits = c(-3, 3),
                         name = "Row Z-score") +
    coord_fixed() +
    labs(x = x_lab, y = y_lab,
         title = paste0(chr, " \u2014 Contrast Sim_Mat (Row Z-score, ", mode_label, ")"),
         subtitle = paste0(n, " windows (", length(idx), "\u00d7", length(idx), " subsampled) | Symmetrised"),
         caption = "Row Z removes family baseline | Red = higher-than-expected | Blue = lower | Clipped \u00b13 SD") +
    THEME_BASE
  p
}

build_simmat_zscore <- function(chr) build_simmat_zscore_core(chr, "mb")
build_simmat_zscore_idx <- function(chr) build_simmat_zscore_core(chr, "window_index")

# -- 17: Contrast sim_mat (row quantile-normalized) — PC2: 3 variants --------
build_simmat_quantile_core <- function(chr, axis_mode = "mb") {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 50) return(NULL)

  q_mat <- sim_mat
  for (ri in seq_len(n)) {
    row_vals <- sim_mat[ri, ]
    valid <- is.finite(row_vals)
    if (sum(valid) < 10) { q_mat[ri, ] <- 0.5; next }
    q_mat[ri, valid] <- frank(row_vals[valid]) / sum(valid)
    q_mat[ri, !valid] <- 0.5
  }
  q_mat <- (q_mat + t(q_mat)) / 2

  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  sub <- q_mat[idx, idx]

  if (axis_mode == "mb") {
    coords <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6
    x_lab <- paste0(chr, " (Mb)"); y_lab <- x_lab
  } else {
    coords <- seq_along(idx)
    x_lab <- "Window index"; y_lab <- x_lab
  }

  pdt <- data.table(
    x = rep(coords, each = length(idx)),
    y = rep(coords, times = length(idx)),
    q = as.vector(sub)
  )
  pdt <- pdt[is.finite(q)]
  mode_label <- if (axis_mode == "mb") "Mb position" else "window index"

  p <- ggplot(pdt, aes(x = x, y = y, fill = q)) +
    geom_tile() +
    scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Row quantile") +
    coord_fixed() +
    labs(x = x_lab, y = y_lab,
         title = paste0(chr, " \u2014 Contrast Sim_Mat (Row Quantile, ", mode_label, ")"),
         subtitle = paste0(n, " windows (", length(idx), "\u00d7", length(idx), " subsampled) | Symmetrised"),
         caption = "Row-quantile removes per-row baseline\nRed = top quantile | Blue = bottom quantile") +
    THEME_BASE
  p
}

build_simmat_quantile <- function(chr) build_simmat_quantile_core(chr, "mb")
build_simmat_quantile_idx <- function(chr) build_simmat_quantile_core(chr, "window_index")

# -- 17c: Distance-normalized contrast sim_mat --------------------------------
build_simmat_distnorm <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 50) return(NULL)

  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  pdt <- build_simmat_tiledata(sim_mat, dt, step, idx, "distance_normalized")
  pdt[, sim_clipped := pmin(pmax(sim, -0.3), 0.3)]

  p <- ggplot(pdt, aes(x = x, y = y, fill = sim_clipped)) +
    geom_tile() +
    scale_fill_gradient2(low = DIV_LOW, mid = DIV_MID, high = DIV_HIGH,
                         midpoint = 0, limits = c(-0.3, 0.3),
                         name = "Dist-norm sim") +
    coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " \u2014 Distance-Normalized Similarity"),
         subtitle = paste0(n, " windows | Expected similarity at each distance subtracted"),
         caption = "Removes diagonal proximity effect | Red = excess similarity | Blue = deficit") +
    THEME_BASE
  p
}

# -- 15: Composite per-chr overview ------------------------------------------
build_composite_page <- function(chr) {
  # Collect individual plots
  p_sim  <- tryCatch(build_sim_heatmap(chr), error = function(e) NULL)
  p_zscore <- tryCatch(build_simmat_zscore(chr), error = function(e) NULL)
  p_quant  <- tryCatch(build_simmat_quantile(chr), error = function(e) NULL)

  plots <- Filter(Negate(is.null), list(p_sim, p_zscore, p_quant))
  if (length(plots) < 2) return(NULL)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    combined <- wrap_plots(plots, ncol = min(3, length(plots))) +
      plot_annotation(
        title = paste0(chr, " \u2014 Matrix Overview"),
        theme = theme(plot.title = element_text(size = 14, face = "bold")))
    return(combined)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    return(gridExtra::arrangeGrob(grobs = plots, ncol = min(3, length(plots)),
      top = grid::textGrob(paste0(chr, " \u2014 Matrix Overview"),
                           gp = grid::gpar(fontsize = 14, fontface = "bold"))))
  }
  NULL
}

# =============================================================================
# DISPATCH
# =============================================================================

plot_builders <- list(
  "06_sim_mat_heatmap"              = build_sim_heatmap,
  "06b_sim_mat_window_idx"          = build_sim_heatmap_idx,
  "14_local_sim_heatmap_80"         = build_local_sim_heatmap,
  "14b_local_sim_heatmap_160"       = build_local_sim_heatmap_160,
  "15_composite_overview"           = build_composite_page,
  "16_simmat_zscore_contrast"       = build_simmat_zscore,
  "16b_simmat_zscore_window_idx"    = build_simmat_zscore_idx,
  "17_simmat_quantile_contrast"     = build_simmat_quantile,
  "17b_simmat_quantile_window_idx"  = build_simmat_quantile_idx,
  "17c_simmat_distance_normalized"  = build_simmat_distnorm
)

plot_dims <- list(
  "06_sim_mat_heatmap"              = c(w = 9, h = 8),
  "06b_sim_mat_window_idx"          = c(w = 9, h = 8),
  "14_local_sim_heatmap_80"         = c(w = 9, h = 8),
  "14b_local_sim_heatmap_160"       = c(w = 9, h = 8),
  "15_composite_overview"           = c(w = 24, h = 8),
  "16_simmat_zscore_contrast"       = c(w = 9, h = 8),
  "16b_simmat_zscore_window_idx"    = c(w = 9, h = 8),
  "17_simmat_quantile_contrast"     = c(w = 9, h = 8),
  "17b_simmat_quantile_window_idx"  = c(w = 9, h = 8),
  "17c_simmat_distance_normalized"  = c(w = 9, h = 8)
)

render_plots(plot_builders, plot_dims, chroms, precomp_list, cli$outdir)
message("[DONE] Matrices -> ", cli$outdir)
