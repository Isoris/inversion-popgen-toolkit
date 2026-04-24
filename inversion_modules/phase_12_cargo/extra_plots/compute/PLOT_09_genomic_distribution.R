#!/usr/bin/env Rscript

# =============================================================================
# PLOT_09_genomic_distribution.R
#
# Stacked windowed-density Manhattan of private vs shared variants along the
# genome. Replaces the circos plot in image 1F with something that scales to
# 28 chromosomes and doesn't require circlize.
#
# Two facets:
#   Top:    private variant density (per 1Mb window)
#   Bottom: shared variant density (per 1Mb window)
# X axis: chromosome offset.
#
# A variant is "private" if its carrier count is 1.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

VARIANT_MASTER <- Sys.getenv("VARIANT_MASTER")
REF_FAI        <- Sys.getenv("REF_FAI")
EXTRAS_FIG_DIR <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR <- Sys.getenv("EXTRAS_TBL_DIR")
WIN <- 1e6  # 1Mb windows

if (!file.exists(VARIANT_MASTER)) { message("[PLOT_09] [skip]"); quit(status = 0) }
if (!file.exists(REF_FAI))        { message("[PLOT_09] [skip] missing FAI"); quit(status = 0) }

v <- fread(VARIANT_MASTER)
nc_col <- intersect(c("n_carriers", "carrier_count", "AC_samples"), names(v))
if (length(nc_col) > 0) {
  v[, n_carriers := get(nc_col[1])]
} else if ("samples_with_alt" %in% names(v)) {
  v[, n_carriers := lengths(strsplit(samples_with_alt, ";"))]
} else {
  message("[PLOT_09] [skip] no carrier counts"); quit(status = 0)
}
chr_col <- if ("chr" %in% names(v)) "chr" else "chrom"
v[, status := ifelse(n_carriers == 1, "private", "shared")]

fai <- fread(REF_FAI, header = FALSE,
             col.names = c("chrom", "len", "off", "x1", "x2"))
chr_order <- fai$chrom
chr_off   <- setNames(c(0, cumsum(fai$len[-nrow(fai)])), chr_order)

# Window assignment
v[, win_start := (get(chr_col) %in% chr_order) * (floor(pos / WIN) * WIN)]
v[, win_id := paste0(get(chr_col), ":", win_start)]
agg <- v[get(chr_col) %in% chr_order,
         .N, by = .(chrom = get(chr_col), win_start, status)]
agg[, abs_x := chr_off[chrom] + win_start + WIN / 2]

fwrite(agg, file.path(EXTRAS_TBL_DIR, "PLOT_09_window_density.tsv"), sep = "\t")

# Chromosome midpoint labels for x-axis
chr_mid <- data.table(
  chrom = chr_order,
  abs_mid = chr_off + fai$len / 2,
  abs_end = chr_off + fai$len)
agg[, chrom := factor(chrom, levels = chr_order)]
agg[, chr_band := match(chrom, chr_order) %% 2]

p <- ggplot(agg, aes(x = abs_x, y = N, color = factor(chr_band))) +
  geom_point(size = 0.6, alpha = 0.7) +
  facet_wrap(~ status, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"),
                       guide = "none") +
  scale_x_continuous(breaks = chr_mid$abs_mid, labels = chr_mid$chrom,
                       expand = c(0.005, 0.005)) +
  labs(title = "Genomic distribution of private vs shared variants",
       subtitle = paste0(format(WIN/1e6, big.mark=","), " Mb windows"),
       x = NULL, y = "Variants per window") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_09_genomic_distribution.pdf"),
       p, width = 16, height = 7, device = cairo_pdf)
message("[PLOT_09] Done")
