#!/usr/bin/env Rscript
# =============================================================================
# plot_theta_pi_multiscale.R
# =============================================================================
# Publication-quality θ_π landscape plots at multiple scales.
#
# Generates three types of plots per run:
#   1. theta_pi_per_chromosome_<scale>.{pdf,png}
#      Per-chromosome facet grid, one panel per LG
#   2. theta_pi_genome_wide_<scale>.{pdf,png}
#      Single-track genome-wide view
#   3. theta_pi_landscape_combined.{pdf,png}
#      Stacked genome-wide view at TWO scales (coarse + fine) for the
#      supplementary figure. This is the main deliverable.
#
# Design choices baked in:
#   - coord_cartesian(ylim=...) NOT ylim() — keeps outliers in data, clips view
#   - Y-axis capped at 99th percentile; outliers shown as red triangles at cap
#   - Dashed line at genome-wide median for reference
#   - LOESS smoother (span tuned to window density) overlays per-window points
#   - Log-y optional for fine scales (θ_π spans 2-3 decades at 5kb)
#   - Per-chr facets share y-axis for comparability
#
# USAGE:
#   Rscript plot_theta_pi_multiscale.R <out_dir> \\
#       <coarse_window_tsv> <coarse_label> \\
#       [<fine_window_tsv> <fine_label>]
#
# The "window tsv" files are ST_theta_pi_per_window.tsv outputs from
# build_theta_supp_tables.py.
#
# Example:
#   Rscript plot_theta_pi_multiscale.R ./theta_figs \\
#       13_theta_pi_supp/win500000_step500000/ST_theta_pi_per_window.tsv "500 kb" \\
#       13_theta_pi_supp/win50000_step10000/ST_theta_pi_per_window.tsv "50 kb"
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plot_theta_pi_multiscale.R <out_dir> <coarse.tsv> <coarse_label> [<fine.tsv> <fine_label>]")
}
out_dir       <- args[1]
coarse_file   <- args[2]
coarse_label  <- args[3]
fine_file     <- if (length(args) >= 5) args[4] else NULL
fine_label    <- if (length(args) >= 5) args[5] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "grey40", size = 8),
        axis.text.x = element_text(size = 7))

# ── Helper: load + aggregate per-window ─────────────────────────────────────
load_windows <- function(path) {
  d <- fread(path)
  # Standardise column names (handle UTF-8 theta column)
  setnames(d,
           old = c("Chromosome", "Window centre (bp)", "n sites",
                   "θ_π (per site)", "θ_W (per site)", "Tajima D"),
           new = c("chrom", "pos", "n_sites",
                   "theta_pi", "theta_W", "tajima"),
           skip_absent = TRUE)
  # Collapse samples → per-window mean (one value per genomic window)
  agg <- d[, .(theta_pi        = mean(theta_pi, na.rm = TRUE),
               theta_pi_median = median(theta_pi, na.rm = TRUE),
               n_samples       = .N),
           by = .(chrom, pos)]
  # Natural chromosome order
  chrom_levels <- unique(agg$chrom)
  chrom_levels <- chrom_levels[order(as.numeric(sub(".*LG0?", "", chrom_levels)))]
  agg[, chrom := factor(chrom, levels = chrom_levels)]
  # Cumulative position for genome-wide view
  agg[, pos_mb := pos / 1e6]
  agg[]
}

# ── Helper: one plot ────────────────────────────────────────────────────────
make_genomewide <- function(dt, label, title_suffix = "") {
  ymax_99 <- quantile(dt$theta_pi, 0.99, na.rm = TRUE)
  med     <- median(dt$theta_pi, na.rm = TRUE)
  mu      <- mean(dt$theta_pi, na.rm = TRUE)
  dt[, is_outlier := theta_pi > ymax_99]
  dt[, plot_y := pmin(theta_pi, ymax_99 * 1.02)]
  n_out <- sum(dt$is_outlier)

  # Build cumulative positions so chromosomes plot end-to-end
  chrom_sizes <- dt[, .(chrom_max = max(pos_mb, na.rm = TRUE)), by = chrom]
  setorder(chrom_sizes, chrom)
  chrom_sizes[, cum_start := cumsum(c(0, chrom_max[-.N]))]
  chrom_sizes[, chrom_mid := cum_start + chrom_max / 2]
  dt <- merge(dt, chrom_sizes[, .(chrom, cum_start)], by = "chrom", all.x = TRUE)
  dt[, abs_mb := pos_mb + cum_start]

  # LOESS span tuned to window count
  span_val <- if (nrow(dt) > 50000) 0.01
              else if (nrow(dt) > 10000) 0.03
              else 0.08

  p <- ggplot(dt, aes(x = abs_mb)) +
    # Chromosome boundaries
    geom_vline(data = chrom_sizes,
               aes(xintercept = cum_start + chrom_max),
               colour = "grey85", linewidth = 0.2) +
    geom_hline(yintercept = med, colour = "grey50",
               linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(y = plot_y), colour = "darkgreen",
               size = 0.25, alpha = 0.35) +
    geom_smooth(aes(y = plot_y), method = "loess", span = span_val,
                se = FALSE, colour = "red", linewidth = 0.6) +
    geom_point(data = dt[is_outlier == TRUE],
               aes(y = ymax_99 * 1.02), shape = 17,
               colour = "darkred", size = 1.0) +
    coord_cartesian(ylim = c(0, ymax_99 * 1.05)) +
    scale_x_continuous(breaks = chrom_sizes$chrom_mid,
                       labels = sub(".*LG", "", chrom_sizes$chrom),
                       expand = c(0, 0)) +
    scale_y_continuous(labels = scientific_format(digits = 2)) +
    labs(title = paste0("θ_π landscape — ", label, " windows", title_suffix),
         subtitle = sprintf(
           "mean = %.2e  |  median = %.2e  |  y-axis capped at 99th percentile (%.2e)  |  %d outlier windows (red ▲)",
           mu, med, ymax_99, n_out
         ),
         x = "Linkage group", y = expression(theta[pi]~"per site")) +
    theme_pub

  attr(p, "label") <- label
  p
}

make_per_chr <- function(dt, label) {
  ymax_99 <- quantile(dt$theta_pi, 0.99, na.rm = TRUE)
  med     <- median(dt$theta_pi, na.rm = TRUE)
  mu      <- mean(dt$theta_pi, na.rm = TRUE)
  dt[, is_outlier := theta_pi > ymax_99]
  dt[, plot_y := pmin(theta_pi, ymax_99 * 1.02)]

  ggplot(dt, aes(x = pos_mb)) +
    geom_hline(yintercept = med, colour = "grey60",
               linetype = "dashed", linewidth = 0.3) +
    geom_line(aes(y = plot_y), colour = "red",
              linewidth = 0.35, alpha = 0.85) +
    geom_point(data = dt[is_outlier == TRUE],
               aes(y = ymax_99 * 1.02), shape = 17,
               colour = "darkred", size = 1.0) +
    facet_wrap(~ chrom, ncol = 4, scales = "free_x") +
    coord_cartesian(ylim = c(0, ymax_99 * 1.05)) +
    scale_y_continuous(labels = scientific_format(digits = 2)) +
    scale_x_continuous(breaks = pretty_breaks(n = 4)) +
    labs(title = paste0("θ_π per chromosome — ", label, " windows"),
         subtitle = sprintf("mean = %.2e  |  median = %.2e  |  outliers as ▲ at y cap (99th pct = %.2e)",
                            mu, med, ymax_99),
         x = "Position (Mb)", y = expression(theta[pi]~"per site")) +
    theme_pub +
    theme(panel.spacing = unit(0.3, "lines"),
          strip.text = element_text(size = 7))
}

# ── Main ────────────────────────────────────────────────────────────────────
cat("Loading coarse:", coarse_file, "\n")
coarse_dt <- load_windows(coarse_file)
cat(sprintf("  %d windows, %d chromosomes\n", nrow(coarse_dt), uniqueN(coarse_dt$chrom)))

p_gw_coarse <- make_genomewide(coarse_dt, coarse_label)
p_pc_coarse <- make_per_chr(coarse_dt, coarse_label)

ggsave(file.path(out_dir, sprintf("theta_pi_genome_wide_%s.pdf",
                                   gsub("[^A-Za-z0-9]", "_", coarse_label))),
       p_gw_coarse, width = 13, height = 3.5, useDingbats = FALSE)
ggsave(file.path(out_dir, sprintf("theta_pi_genome_wide_%s.png",
                                   gsub("[^A-Za-z0-9]", "_", coarse_label))),
       p_gw_coarse, width = 13, height = 3.5, dpi = 150)

ggsave(file.path(out_dir, sprintf("theta_pi_per_chromosome_%s.pdf",
                                   gsub("[^A-Za-z0-9]", "_", coarse_label))),
       p_pc_coarse, width = 13, height = 11, useDingbats = FALSE)
ggsave(file.path(out_dir, sprintf("theta_pi_per_chromosome_%s.png",
                                   gsub("[^A-Za-z0-9]", "_", coarse_label))),
       p_pc_coarse, width = 13, height = 11, dpi = 150)

cat("  wrote coarse-scale plots.\n")

if (!is.null(fine_file)) {
  cat("Loading fine:", fine_file, "\n")
  fine_dt <- load_windows(fine_file)
  cat(sprintf("  %d windows, %d chromosomes\n", nrow(fine_dt), uniqueN(fine_dt$chrom)))

  p_gw_fine <- make_genomewide(fine_dt, fine_label)
  p_pc_fine <- make_per_chr(fine_dt, fine_label)

  ggsave(file.path(out_dir, sprintf("theta_pi_genome_wide_%s.pdf",
                                     gsub("[^A-Za-z0-9]", "_", fine_label))),
         p_gw_fine, width = 13, height = 3.5, useDingbats = FALSE)
  ggsave(file.path(out_dir, sprintf("theta_pi_genome_wide_%s.png",
                                     gsub("[^A-Za-z0-9]", "_", fine_label))),
         p_gw_fine, width = 13, height = 3.5, dpi = 150)

  ggsave(file.path(out_dir, sprintf("theta_pi_per_chromosome_%s.pdf",
                                     gsub("[^A-Za-z0-9]", "_", fine_label))),
         p_pc_fine, width = 13, height = 11, useDingbats = FALSE)
  ggsave(file.path(out_dir, sprintf("theta_pi_per_chromosome_%s.png",
                                     gsub("[^A-Za-z0-9]", "_", fine_label))),
         p_pc_fine, width = 13, height = 11, dpi = 150)

  cat("  wrote fine-scale plots.\n")

  # ── Combined supplementary figure: both scales stacked ───────────────────
  combined <- p_gw_coarse / p_gw_fine +
    plot_annotation(
      title = "Nucleotide diversity (θ_π) landscape at two scales",
      subtitle = sprintf(
        "Top: %s windows (chromosome-scale view, used for statistical summaries). Bottom: %s windows (fine-scale view, for inspection).",
        coarse_label, fine_label
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 12),
                    plot.subtitle = element_text(colour = "grey40", size = 9))
    )

  ggsave(file.path(out_dir, "theta_pi_landscape_combined.pdf"),
         combined, width = 13, height = 8, useDingbats = FALSE)
  ggsave(file.path(out_dir, "theta_pi_landscape_combined.png"),
         combined, width = 13, height = 8, dpi = 150)

  cat("  wrote combined supplementary figure.\n")
}

cat("\nDone. Output:", out_dir, "\n")
