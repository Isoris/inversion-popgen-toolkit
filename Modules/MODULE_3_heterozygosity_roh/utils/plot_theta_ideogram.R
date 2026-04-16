#!/usr/bin/env Rscript
# =============================================================================
# plot_theta_ideogram.R — Local theta/diversity ideogram-like tracks
# =============================================================================
# Adapted from plot_theta_ideogram_simple.R
# Produces genome-wide ideogram view of local theta (diversity proxy).
#
# NOTE: These are LOCAL DIVERSITY TRACKS, not literal per-site Hobs.
#
# Usage:
#   Rscript plot_theta_ideogram.R \
#     <theta_dir> <chrom_sizes.tsv> <out_dir> [sample_list.txt] [max_samples=10]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plot_theta_ideogram.R <theta_dir> <chrom_sizes.tsv> <out_dir> [sample_list] [max_samples]")
}

theta_dir  <- args[1]
sizes_file <- args[2]
out_dir    <- args[3]
sample_file <- if (length(args) >= 4 && file.exists(args[4])) args[4] else NULL
max_samples <- if (length(args) >= 5) as.integer(args[5]) else 10

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))

# ── Chromosome sizes ─────────────────────────────────────────────────────
sizes <- fread(sizes_file, header = FALSE, col.names = c("chrom", "len"))
chrom_order <- sizes$chrom
sizes[, cum_start := cumsum(c(0, len[-.N]))]

# ── Find theta window files ──────────────────────────────────────────────
theta_files <- list.files(theta_dir, pattern = "\\.pestPG$", full.names = TRUE)
if (length(theta_files) == 0) {
  cat("No .pestPG files found in", theta_dir, "\n")
  quit(save = "no")
}

# Extract sample names from filenames
get_sample <- function(f) {
  bn <- basename(f)
  sub("\\.win.*$", "", bn)
}

samples <- unique(sapply(theta_files, get_sample))

# Optionally filter to sample list
if (!is.null(sample_file)) {
  keep <- fread(sample_file, header = FALSE)$V1
  samples <- intersect(samples, keep)
}
if (length(samples) > max_samples) {
  cat("Plotting first", max_samples, "of", length(samples), "samples\n")
  samples <- samples[1:max_samples]
}

# ── Read and combine theta windows ───────────────────────────────────────
all_theta <- list()
for (samp in samples) {
  f <- theta_files[grepl(paste0("^", samp, "\\."), basename(theta_files))]
  if (length(f) == 0) next
  f <- f[1]

  t <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(t) || nrow(t) == 0) next

  # Auto-detect columns (thetaStat output varies)
  cn <- names(t)
  chr_col <- cn[cn %in% c("Chr", "chr", "Chromo", "chrom")][1]
  center_col <- cn[cn %in% c("WinCenter", "midPos", "wincenter")][1]
  tp_col <- cn[cn %in% c("tP", "thetaPi", "ThetaPi", "pi")][1]
  nsites_col <- cn[cn %in% c("nSites", "nsite", "Nsites")][1]

  if (any(is.na(c(chr_col, center_col, tp_col)))) next

  sub_t <- t[, .SD, .SDcols = c(chr_col, center_col, tp_col,
                                  if (!is.na(nsites_col)) nsites_col)]
  setnames(sub_t, c(chr_col, center_col, tp_col), c("chrom", "center", "tP"))
  if (!is.na(nsites_col)) setnames(sub_t, nsites_col, "nSites")

  sub_t[, sample := samp]
  # Normalize: theta per site
  if ("nSites" %in% names(sub_t)) {
    sub_t[, theta_per_site := tP / nSites]
  } else {
    sub_t[, theta_per_site := tP]
  }
  all_theta[[samp]] <- sub_t
}

if (length(all_theta) == 0) {
  cat("No theta data loaded.\n")
  quit(save = "no")
}

dt <- rbindlist(all_theta, fill = TRUE)

# Compute genome-wide position
dt <- merge(dt, sizes[, .(chrom, cum_start)], by = "chrom")
dt[, gw_pos := cum_start + center]
dt[, chrom := factor(chrom, levels = chrom_order)]

# ── Per-sample ideogram ─────────────────────────────────────────────────
for (samp in unique(dt$sample)) {
  sub <- dt[sample == samp]
  if (nrow(sub) < 10) next

  # Genome-wide view
  p <- ggplot(sub, aes(x = gw_pos / 1e6, y = theta_per_site)) +
    geom_point(size = 0.3, alpha = 0.4, color = "steelblue") +
    geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "red", linewidth = 0.5) +
    labs(
      title = paste0(samp, " — Local theta diversity proxy (NOT literal Hobs)"),
      x = "Genomic position (Mb, genome-wide)",
      y = "Theta/site (diversity proxy)"
    ) +
    theme_pub

  ggsave(file.path(out_dir, paste0(samp, "_theta_genomewide.png")), p, width = 14, height = 4, dpi = 300)

  # Per-chromosome faceted view
  p2 <- ggplot(sub, aes(x = center / 1e6, y = theta_per_site)) +
    geom_point(size = 0.2, alpha = 0.3, color = "steelblue") +
    facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
    labs(
      title = paste0(samp, " — Per-chromosome local diversity proxy"),
      x = "Position (Mb)", y = "Theta/site"
    ) +
    theme_pub +
    theme(axis.text.x = element_text(size = 6))

  ggsave(file.path(out_dir, paste0(samp, "_theta_per_chr.png")), p2, width = 14, height = 10, dpi = 300)
}

# ── Multi-sample overlay (mean across samples per window) ────────────────
mean_dt <- dt[, .(
  mean_theta = mean(theta_per_site, na.rm = TRUE),
  sd_theta = sd(theta_per_site, na.rm = TRUE),
  n_samples = .N
), by = .(chrom, center)]

mean_dt <- merge(mean_dt, sizes[, .(chrom, cum_start)], by = "chrom")
mean_dt[, gw_pos := cum_start + center]
mean_dt[, chrom := factor(chrom, levels = chrom_order)]

p_mean <- ggplot(mean_dt, aes(x = gw_pos / 1e6, y = mean_theta)) +
  geom_point(size = 0.3, alpha = 0.4, color = "darkgreen") +
  geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "red", linewidth = 0.5) +
  labs(
    title = "Mean local theta diversity proxy across all samples",
    subtitle = "(NOT literal per-site Hobs)",
    x = "Genomic position (Mb)", y = "Mean theta/site"
  ) +
  theme_pub

ggsave(file.path(out_dir, "mean_theta_genomewide.png"), p_mean, width = 14, height = 4, dpi = 300)
ggsave(file.path(out_dir, "mean_theta_genomewide.pdf"), p_mean, width = 14, height = 4)

cat("Theta ideogram plots written to:", out_dir, "\n")
