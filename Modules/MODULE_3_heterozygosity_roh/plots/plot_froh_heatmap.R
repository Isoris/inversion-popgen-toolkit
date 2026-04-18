#!/usr/bin/env Rscript
# =============================================================================
# plot_froh_heatmap.R  (v3)
# =============================================================================
# Generates one F_ROH heatmap per palette, side by side, so you can pick the
# best one visually. By default runs all palettes listed in PALETTES_TO_RUN.
#
# Palettes produced:
#   warm    — pale cream → orange → red → dark red-brown (biology paper look)
#   rocket  — viridisLite rocket (purple → red → yellow, warm & colorblind-safe)
#   mako    — viridisLite mako   (navy → teal → green, cool)
#   turbo   — viridisLite turbo  (rainbow, punchy, not perceptually uniform)
#   ylorbr  — RColorBrewer-style YlOrBr-like (warm sequential)
#   firered — pale yellow → orange → red → black (strong high-end contrast)
#
# Design baked in for all palettes:
#   - Rows ordered by total F_ROH (highest at top)
#   - NA cells grey
#   - F_ROH capped at 1; overshoot cells flagged with cyan ×
#   - 600 dpi PNGs for print; PDF is vector
#   - Bottom horizontal colour bar
#
# USAGE:
#   Rscript plot_froh_heatmap.R per_chr_roh_summary.tsv out_dir \\
#       [sample_list.txt] [palette1,palette2,...]
#
# Example — all palettes, pruned81:
#   Rscript plot_froh_heatmap.R per_chr_roh_summary.tsv ./heatmaps/ \\
#       pruned81_samples.tsv
#
# Example — just two palettes for all226:
#   Rscript plot_froh_heatmap.R per_chr_roh_summary.tsv ./heatmaps/ \\
#       NONE warm,rocket
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridisLite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: plot_froh_heatmap.R per_chr_roh_summary.tsv out_dir [sample_list.txt|NONE] [palettes_csv]")
}
infile       <- args[1]
out_dir      <- args[2]
samples_file <- if (length(args) >= 3 && args[3] != "NONE") args[3] else NULL
palettes_arg <- if (length(args) >= 4) args[4] else "warm,rocket,mako,turbo,ylorbr,firered"
PALETTES_TO_RUN <- trimws(strsplit(palettes_arg, ",")[[1]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load + prep data once ──────────────────────────────────────────────────
dt <- fread(infile)
stopifnot(all(c("sample", "chrom", "froh_chr") %in% names(dt)))

if (!is.null(samples_file) && file.exists(samples_file)) {
  keep <- fread(samples_file, header = FALSE)$V1
  dt <- dt[sample %in% keep]
  cat(sprintf("Restricted to %d samples from %s\n",
              length(unique(dt$sample)), basename(samples_file)))
}

dt[, froh_plot := pmin(froh_chr, 1.0)]
dt[, overshoot := froh_chr > 1.0]
n_overshoot <- sum(dt$overshoot, na.rm = TRUE)

chrom_levels <- unique(dt$chrom)
chrom_levels <- chrom_levels[order(as.numeric(sub(".*LG0?", "", chrom_levels)))]
dt[, chrom := factor(chrom, levels = chrom_levels)]

sample_order <- dt[, .(total_froh = sum(froh_chr, na.rm = TRUE)), by = sample]
setorder(sample_order, -total_froh)
dt[, sample := factor(sample, levels = sample_order$sample)]

n_samples <- length(levels(dt$sample))
show_labels <- n_samples <= 100
subset_tag <- if (!is.null(samples_file)) sub("\\.tsv$", "", basename(samples_file)) else "all"

cat(sprintf("Data: %d samples × %d chromosomes. %d cells with F_ROH > 1.\n",
            n_samples, length(chrom_levels), n_overshoot))

# ── Palette dispatcher ─────────────────────────────────────────────────────
get_scale <- function(pal) {
  switch(pal,
    "warm" = scale_fill_gradientn(
      colours = c("#FFFBE6", "#FDD49E", "#FDAE61", "#F46D43", "#D73027", "#7F0000"),
      values = scales::rescale(c(0, 0.15, 0.30, 0.50, 0.75, 1.0)),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      na.value = "grey88", name = expression(F[ROH]),
      guide = guide_colorbar(barwidth = 12, barheight = 0.7,
                             title.position = "top", title.hjust = 0.5)
    ),
    "rocket" = scale_fill_viridis_c(
      option = "rocket", direction = -1, limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      na.value = "grey88", name = expression(F[ROH]),
      guide = guide_colorbar(barwidth = 12, barheight = 0.7,
                             title.position = "top", title.hjust = 0.5)
    ),
    "mako" = scale_fill_viridis_c(
      option = "mako", direction = -1, limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      na.value = "grey88", name = expression(F[ROH]),
      guide = guide_colorbar(barwidth = 12, barheight = 0.7,
                             title.position = "top", title.hjust = 0.5)
    ),
    "turbo" = scale_fill_viridis_c(
      option = "turbo", limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      na.value = "grey88", name = expression(F[ROH]),
      guide = guide_colorbar(barwidth = 12, barheight = 0.7,
                             title.position = "top", title.hjust = 0.5)
    ),
    "ylorbr" = scale_fill_gradientn(
      colours = c("#FFFFE5", "#FFF7BC", "#FEE391", "#FEC44F", "#FE9929",
                  "#EC7014", "#CC4C02", "#993404", "#662506"),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      na.value = "grey88", name = expression(F[ROH]),
      guide = guide_colorbar(barwidth = 12, barheight = 0.7,
                             title.position = "top", title.hjust = 0.5)
    ),
    "firered" = scale_fill_gradientn(
      colours = c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404", "#000000"),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      na.value = "grey88", name = expression(F[ROH]),
      guide = guide_colorbar(barwidth = 12, barheight = 0.7,
                             title.position = "top", title.hjust = 0.5)
    ),
    stop("Unknown palette: ", pal)
  )
}

# ── Loop and render ────────────────────────────────────────────────────────
for (pal in PALETTES_TO_RUN) {
  cat(sprintf("  rendering palette: %s\n", pal))
  p <- ggplot(dt, aes(x = chrom, y = sample, fill = froh_plot)) +
    geom_tile(colour = NA) +
    get_scale(pal) +
    geom_point(data = dt[overshoot == TRUE],
               aes(x = chrom, y = sample),
               shape = 4, colour = "cyan3", size = 0.6, inherit.aes = FALSE) +
    scale_x_discrete(labels = function(x) sub(".*LG", "LG", x), expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Chromosome", y = "Sample",
         title = sprintf("Sample × chromosome F_ROH (%d samples, palette: %s)",
                         n_samples, pal),
         subtitle = "Rows ordered by total F_ROH (highest at top). Cyan × = F_ROH > 1. Grey = no data.",
         caption = sprintf("Cells with F_ROH > 1: %d (%.2f%%).",
                           n_overshoot, 100 * n_overshoot / nrow(dt))) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey40", linewidth = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = if (show_labels) element_text(size = 5) else element_blank(),
      axis.ticks.y = if (show_labels) element_line() else element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(colour = "grey40", size = 8),
      plot.caption = element_text(colour = "grey50", size = 7, hjust = 0)
    )

  h <- max(4, min(20, n_samples * 0.12 + 2.5))
  outbase <- file.path(out_dir,
                       sprintf("heatmap_froh_%s_%s", subset_tag, pal))
  ggsave(paste0(outbase, ".pdf"), p, width = 11, height = h, useDingbats = FALSE)
  ggsave(paste0(outbase, ".png"), p, width = 11, height = h, dpi = 600)
}

cat(sprintf("\nDone. Wrote %d palette variants to %s/\n",
            length(PALETTES_TO_RUN), out_dir))
