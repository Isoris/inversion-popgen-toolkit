#!/usr/bin/env Rscript
# =============================================================================
# C01a_diag_matrices.R — Similarity matrix heatmaps and transforms (plots 06, 14, 16, 17)
# Computationally heavy (matrix operations on subsampled sim_mat).
# Usage: Rscript C01a_diag_matrices.R <precomp_dir> <outdir> [chrom]
# =============================================================================

SCRIPT_DIR <- dirname(sys.frame(1)$ofile %||% ".")
source(file.path(SCRIPT_DIR, "C01a_diag_common.R"))

cli <- parse_diag_args()
data <- load_diag_data(cli$precomp_dir, cli$chrom_filter)
precomp_list <- data$precomp_list
chroms <- data$chroms
inv_like_dt <- data$inv_like_dt
dir.create(cli$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# MATRIX BUILDERS
# =============================================================================

# -- 06: Similarity matrix heatmap -----------------------------------
build_sim_heatmap <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  n <- nrow(sim_mat)
  if (n < 10) return(NULL)

  # Subsample for large chromosomes (max 500x500)
  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  sub <- sim_mat[idx, idx]
  pos_mb <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6

  pdt <- data.table(
    x = rep(pos_mb, each = length(idx)),
    y = rep(pos_mb, times = length(idx)),
    sim = as.vector(sub)
  )
  pdt <- pdt[is.finite(sim)]

  # Center color on chr median for better contrast
  sim_med <- median(pdt$sim, na.rm = TRUE)

  p <- ggplot(pdt, aes(x = x, y = y, fill = sim)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("navy", "blue3", "steelblue", "gold", "orange", "red3"),
      values = scales::rescale(c(0, sim_med * 0.6, sim_med, sim_med * 1.1, sim_med * 1.3, 1)),
      name = "Similarity") +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " -- Similarity Matrix (", n, " windows)"),
         subtitle = paste0(length(idx), "x", length(idx), " subsampled | ",
                          "median sim = ", round(sim_med, 3)),
         caption = paste0("Source: MDS-space correlation between local PCA loading vectors\n",
                         "Red blocks on diagonal = windows sharing PCA structure (inversion candidates)\n",
                         "Off-diagonal blocks = long-range similarity (shared inversion across distant windows)\n",
                         "Blue = low similarity (background or boundary between systems)")) +
    coord_fixed() + THEME_BASE
  p
}

# -- 14: Local similarity heatmap (band around diagonal) --------------
build_local_sim_heatmap <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  if (is.null(sim_mat) || nrow(sim_mat) < 20) return(NULL)

  n_w <- nrow(sim_mat)
  bandwidth <- 80L

  # Denser subsampling (max 1000 points along diagonal)
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

  # Center color scale on median similarity for better contrast
  med_sim <- median(pdt$sim, na.rm = TRUE)

  p <- ggplot(pdt, aes(x = x_mb, y = y_mb, fill = sim)) +
    geom_tile(width = step * (dt$end_bp[2] - dt$start_bp[1]) / 1e6,
              height = step * (dt$end_bp[2] - dt$start_bp[1]) / 1e6) +
    scale_fill_gradientn(
      colours = c("navy", "steelblue", "white", "orange", "red3"),
      values = scales::rescale(c(0, med_sim * 0.7, med_sim, med_sim * 1.15, 1)),
      name = "Similarity"
    ) +
    coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " -- Local Similarity (+-", bandwidth, " windows)"),
         subtitle = paste0("Snake extension context | median sim = ", round(med_sim, 3)),
         caption = "White = median similarity | Red = high (block) | Blue = low (gap)") +
    THEME_BASE
  p
}


# -- 15: Composite per-chr overview (all key plots on one page) -------
build_composite_page <- function(chr) {
  # Collect individual plots
  p_eigen <- tryCatch(build_eigen_profile(chr), error = function(e) NULL)
  p_inv   <- tryCatch(build_inv_likeness_profile(chr), error = function(e) NULL)
  p_z     <- tryCatch(build_z_profile(chr), error = function(e) NULL)
  p_nn    <- tryCatch(build_nn_profile(chr), error = function(e) NULL)
  p_mds   <- tryCatch(build_mds_raw(chr), error = function(e) NULL)
  p_knn   <- tryCatch(build_knn_structure(chr), error = function(e) NULL)

  # Filter NULLs
  plots <- Filter(Negate(is.null), list(p_eigen, p_inv, p_z, p_nn, p_mds, p_knn))
  if (length(plots) < 3) return(NULL)

  # Use patchwork if available, otherwise gridExtra
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    # Layout: top row = 3 wide panels, bottom row = 3 panels
    n_plots <- length(plots)
    if (n_plots >= 6) {
      composite <- (plots[[1]] | plots[[2]] | plots[[3]]) /
                   (plots[[4]] | plots[[5]] | plots[[6]]) +
        plot_annotation(
          title = paste0(chr, " -- Complete Diagnostic Overview"),
          subtitle = paste0(precomp_list[[chr]]$n_windows, " windows"),
          theme = theme(plot.title = element_text(size = 14, face = "bold"))
        )
    } else {
      composite <- wrap_plots(plots, ncol = min(3, n_plots)) +
        plot_annotation(
          title = paste0(chr, " -- Complete Diagnostic Overview"),
          theme = theme(plot.title = element_text(size = 14, face = "bold"))
        )
    }
    return(composite)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    n_plots <- length(plots)
    ncol_g <- min(3L, n_plots)
    nrow_g <- ceiling(n_plots / ncol_g)
    composite <- gridExtra::arrangeGrob(
      grobs = plots, ncol = ncol_g, nrow = nrow_g,
      top = grid::textGrob(paste0(chr, " -- Complete Diagnostic Overview"),
                           gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
    return(composite)
  } else {
    message("  [SKIP] composite: install patchwork or gridExtra")
    return(NULL)
  }
}

# =============================================================================
# RENDER PER-CHROMOSOME PLOTS
# =============================================================================


# -- 16: Contrast-enhanced sim_mat (row-wise Z-score) ----------------
# Each row of the sim_mat is Z-scored: subtract row mean, divide by row SD.
# This removes the global baseline inflation from family structure. Windows
# that are "more similar than expected" relative to their own row stand out.
# NOT a replacement for raw sim_mat — an additional diagnostic view.
build_simmat_zscore <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 50) return(NULL)

  # Row-wise Z-score
  z_mat <- sim_mat
  for (ri in seq_len(n)) {
    row_vals <- sim_mat[ri, ]
    rm <- mean(row_vals, na.rm = TRUE)
    rsd <- sd(row_vals, na.rm = TRUE)
    if (is.finite(rsd) && rsd > 1e-8) {
      z_mat[ri, ] <- (row_vals - rm) / rsd
    } else {
      z_mat[ri, ] <- 0
    }
  }

  # Symmetrise: average with transpose for cleaner visual
  z_mat <- (z_mat + t(z_mat)) / 2

  # Subsample (max 500x500)
  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  sub <- z_mat[idx, idx]
  pos_mb <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6

  pdt <- data.table(
    x = rep(pos_mb, each = length(idx)),
    y = rep(pos_mb, times = length(idx)),
    z = as.vector(sub)
  )
  pdt <- pdt[is.finite(z)]
  # Clip extremes for display
  pdt[, z_clipped := pmin(pmax(z, -3), 3)]

  p <- ggplot(pdt, aes(x = x, y = y, fill = z_clipped)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "grey95", high = "red3",
                         midpoint = 0, limits = c(-3, 3),
                         name = "Row Z-score") +
    coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " -- Contrast Sim_Mat (Row-wise Z-score)"),
         subtitle = paste0(n, " windows (", length(idx), "x", length(idx),
                          " subsampled) | Symmetrised"),
         caption = paste0("Row Z removes family baseline inflation\n",
                         "Red = higher-than-expected similarity (inversion block)\n",
                         "Blue = lower-than-expected (gap/boundary)\n",
                         "Clipped at +-3 SD")) +
    THEME_BASE
  p
}


# -- 17: Contrast-enhanced sim_mat (row quantile-normalized) ----------
# Per-row rank transform: replace each value with its within-row quantile.
# Produces a [0,1] matrix where blocks stand out as regions where nearby
# windows have unusually high rank relative to far-away windows.
build_simmat_quantile <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 50) return(NULL)

  # Row-wise quantile normalization
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
  pos_mb <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6

  pdt <- data.table(
    x = rep(pos_mb, each = length(idx)),
    y = rep(pos_mb, times = length(idx)),
    q = as.vector(sub)
  )
  pdt <- pdt[is.finite(q)]

  p <- ggplot(pdt, aes(x = x, y = y, fill = q)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red3",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Row quantile") +
    coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " -- Contrast Sim_Mat (Row Quantile-Normalized)"),
         subtitle = paste0(n, " windows (", length(idx), "x", length(idx),
                          " subsampled) | Symmetrised"),
         caption = paste0("Row-quantile removes per-row baseline\n",
                         "Red = similarity in top quantile for that row\n",
                         "Blue = similarity in bottom quantile\n",
                         "Blocks visible where nearby windows outrank distant")) +
    THEME_BASE
  p
}


# =============================================================================
# DISPATCH
# =============================================================================

plot_builders <- list(
  "06_sim_mat_heatmap" = build_sim_heatmap,
  "14_local_sim_heatmap" = build_local_sim_heatmap,
  "15_composite_overview" = build_composite_page,
  "16_simmat_zscore_contrast" = build_simmat_zscore,
  "17_simmat_quantile_contrast" = build_simmat_quantile
)

plot_dims <- list(
  "06_sim_mat_heatmap" = c(w = 8, h = 7),
  "14_local_sim_heatmap" = c(w = 8, h = 7),
  "15_composite_overview" = c(w = 20, h = 14),
  "16_simmat_zscore_contrast" = c(w = 8, h = 7),
  "17_simmat_quantile_contrast" = c(w = 8, h = 7)
)

render_plots(plot_builders, plot_dims, chroms, precomp_list, cli$outdir)
message("[DONE] Matrices -> ", cli$outdir)
