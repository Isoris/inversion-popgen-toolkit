#!/usr/bin/env Rscript
###############################################################################
# plot_evaladmix_diagnostic.R — Per-K evalAdmix residual diagnostics
#
# For each K: upper-triangle heatmap of residual correlations, per-sample
# burden barplot. Produces one PDF per K and a faceted grid overview.
#
# Usage:
#   Rscript plot_evaladmix_diagnostic.R \
#     --eval-dir structure_results/evaladmix/ \
#     --file-prefix wholegenome_thin500_all226 \
#     --sample-file samples.txt \
#     --theme-file utils/theme_systems_plate.R \
#     --out-dir structure_results/figures/ \
#     --k-min 2 --k-max 20
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(optparse)
})

option_list <- list(
  make_option("--eval-dir", type = "character"),
  make_option("--file-prefix", type = "character"),
  make_option("--sample-file", type = "character"),
  make_option("--theme-file", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--k-min", type = "integer", default = 2),
  make_option("--k-max", type = "integer", default = 20),
  make_option("--grid-only", action = "store_true", default = FALSE,
              help = "Only produce grid overview, skip per-K PDFs")
)
opt <- parse_args(OptionParser(option_list = option_list))

source(opt[["theme-file"]])

dir.create(opt[["out-dir"]], showWarnings = FALSE, recursive = TRUE)

samples <- trimws(readLines(opt[["sample-file"]], warn = FALSE))
samples <- samples[nzchar(samples)]
n <- length(samples)

# =============================================================================
# LOAD AND PLOT EACH K
# =============================================================================

grid_panels <- list()
grid_K <- integer()

for (K in seq(opt[["k-min"]], opt[["k-max"]])) {
  stem <- sprintf("%s_K%02d_best", opt[["file-prefix"]], K)
  cor_file <- file.path(opt[["eval-dir"]], paste0(stem, ".corres.txt"))

  if (!file.exists(cor_file)) {
    cat("[SKIP] K=", K, " — missing\n")
    next
  }

  M <- tryCatch(
    as.matrix(read.table(cor_file, header = FALSE)),
    error = function(e) NULL
  )
  if (is.null(M) || nrow(M) != n) {
    cat("[SKIP] K=", K, " — bad matrix\n")
    next
  }

  diag(M) <- NA

  # Per-sample burden
  burden <- rowMeans(abs(M), na.rm = TRUE)
  burden_dt <- data.table(
    sample_idx = seq_len(n),
    sample = samples,
    burden = burden
  )

  # Melt matrix for heatmap
  rownames(M) <- colnames(M) <- seq_len(n)
  melt_dt <- data.table(
    row = rep(seq_len(n), n),
    col = rep(seq_len(n), each = n),
    value = as.vector(M)
  )

  # ---- Per-K figure (heatmap + burden bar) ----
  if (!opt[["grid-only"]]) {
    p_heat <- ggplot(melt_dt, aes(x = col, y = n - row + 1, fill = value)) +
      geom_tile(color = NA) +
      scale_fill_gradientn(
        colors = c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
                   "#FDDBC7", "#EF8A62", "#B2182B"),
        limits = c(-0.3, 0.3), oob = scales::squish,
        name = "Residual\ncorrelation"
      ) +
      theme_plate(grid = "none") +
      theme(
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), axis.line = element_blank()
      ) +
      coord_fixed() +
      labs(title = paste0("K=", K, " evalAdmix residuals"))

    p_bar <- ggplot(burden_dt, aes(x = sample_idx, y = burden)) +
      geom_col(width = 1, fill = "#E34A33", alpha = 0.7) +
      theme_plate_compact(grid = "y", panel_border = TRUE) +
      labs(x = "Sample", y = "Mean |resid|") +
      scale_x_continuous(expand = c(0, 0))

    p_k <- p_heat / p_bar + plot_layout(heights = c(4, 1))

    out_k <- file.path(opt[["out-dir"]],
                       paste0(opt[["file-prefix"]], "_evaladmix_K", sprintf("%02d", K), ".pdf"))
    ggsave(out_k, p_k, width = 10, height = 11, device = cairo_pdf, bg = PLATE_BG)
    cat("[OK] K=", K, " →", out_k, "\n")
  }

  # ---- Grid panel (burden only) ----
  p_grid <- ggplot(burden_dt, aes(x = sample_idx, y = burden)) +
    geom_col(width = 1, fill = "#E34A33", alpha = 0.7) +
    theme_plate_compact(grid = "none", panel_border = TRUE) +
    labs(title = paste0("K=", K), x = NULL, y = NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(
      plot.title = element_text(size = 7, face = "bold"),
      axis.text = element_text(size = 4),
      plot.margin = margin(2, 2, 2, 2)
    )

  grid_panels[[length(grid_panels) + 1]] <- p_grid
  grid_K <- c(grid_K, K)
}

# =============================================================================
# GRID OVERVIEW
# =============================================================================

if (length(grid_panels) > 0) {
  n_panels <- length(grid_panels)
  ncols <- min(5, n_panels)
  nrows <- ceiling(n_panels / ncols)

  grid_fig <- wrap_plots(grid_panels, ncol = ncols) +
    plot_annotation(
      title = "evalAdmix per-sample burden across K values",
      subtitle = paste0(opt[["file-prefix"]], " — per-sample mean |residual correlation|"),
      theme = theme(
        plot.title = element_text(size = SIZE_MAIN_TITLE, face = "bold", color = TEXT_TITLE),
        plot.subtitle = element_text(size = SIZE_SUBTITLE, color = TEXT_SUBTITLE),
        plot.background = element_rect(fill = PLATE_BG, color = NA)
      )
    )

  out_grid <- file.path(opt[["out-dir"]],
                        paste0(opt[["file-prefix"]], "_evaladmix_grid.pdf"))
  ggsave(out_grid, grid_fig, width = 16, height = 3 * nrows + 1,
         device = cairo_pdf, bg = PLATE_BG)
  ggsave(sub("\\.pdf$", ".png", out_grid), grid_fig,
         width = 16, height = 3 * nrows + 1, dpi = 300, bg = PLATE_BG)
  cat("[DONE] Grid:", out_grid, "\n")
}
