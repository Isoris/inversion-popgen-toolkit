#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

option_list <- list(
  make_option("--summary_dir", type="character"),
  make_option("--plot_dir", type="character"),
  make_option("--ref_fai", type="character"),
  make_option("--samples_unrelated", type="character", default="")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$plot_dir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank())

col_main <- "#2E8B8B"
col_light <- "#7ECACA"

master <- fread(file.path(opt$summary_dir, "master_DUP_annotation_226.tsv"))
persamp <- fread(file.path(opt$summary_dir, "per_sample_DUP_summary.tsv"))
freq_summ <- fread(file.path(opt$summary_dir, "frequency_class_summary.tsv"))
size_summ <- fread(file.path(opt$summary_dir, "size_class_summary.tsv"))
chr_summ <- fread(file.path(opt$summary_dir, "per_chromosome_summary.tsv"))

unrel <- character(0)
if (file.exists(opt$samples_unrelated)) {
  unrel <- scan(opt$samples_unrelated, what = "character", quiet = TRUE)
}

# Frequency
freq_order <- c("singleton","rare","low_frequency","common","widespread","near_fixed")
freq_summ[, frequency_class := factor(frequency_class, levels = freq_order)]
p_freq <- ggplot(freq_summ, aes(x = frequency_class, y = count)) +
  geom_col(fill = col_main) +
  geom_text(aes(label = count), vjust = -0.3, size = 3) +
  labs(x = "Frequency class", y = "Count", title = "DUP frequency spectrum") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(opt$plot_dir, "freq_class_barplot.pdf"), p_freq, width = 8, height = 5, dpi = 400)

# Size class
size_summ[, size_class := factor(size_class, levels = c("sub_SV","small","medium","large"))]
p_size <- ggplot(size_summ, aes(x = size_class, y = count)) +
  geom_col(fill = col_main) +
  geom_text(aes(label = count), vjust = -0.3, size = 3) +
  labs(x = "Size class", y = "Count", title = "DUP count by size class") +
  theme_pub
ggsave(file.path(opt$plot_dir, "size_class_barplot.pdf"), p_size, width = 6, height = 5, dpi = 400)

# Functional class
func_order <- c("intergenic","intronic","exon_overlap","CDS_overlap")
master[, func_class := factor(func_class, levels = func_order)]
p_func <- ggplot(master, aes(x = func_class)) +
  geom_bar(fill = col_main) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3, size = 3) +
  labs(x = "Functional class", y = "Count", title = "DUP count by functional class") +
  theme_pub
ggsave(file.path(opt$plot_dir, "func_class_barplot.pdf"), p_func, width = 7, height = 5, dpi = 400)

# Per-sample burden
persamp[, subset_fill := ifelse(sample %in% unrel, "Unrelated (81)", "Related")]
persamp <- persamp[order(-total_DUP)]
persamp[, sample := factor(sample, levels = sample)]

p_burden <- ggplot(persamp, aes(x = sample, y = total_DUP, fill = subset_fill)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("Unrelated (81)" = col_main, "Related" = col_light), name = "Subset") +
  labs(x = "Sample", y = "DUP count", title = "Per-sample DUP burden") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 2.5))
ggsave(file.path(opt$plot_dir, "per_sample_DUP_burden.pdf"), p_burden, width = 16, height = 5, dpi = 400)

# PCA
p_pca <- ggplot(persamp, aes(x = PC1, y = PC2, color = subset_fill)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Unrelated (81)" = col_main, "Related" = col_light), name = "Subset") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "PCA on DUP presence/absence") +
  theme_pub
ggsave(file.path(opt$plot_dir, "PCA_DUP.pdf"), p_pca, width = 8, height = 6, dpi = 400)

# Chromosome density
chr_summ[, chrom_label := gsub("C_gar_", "", chrom)]
chr_summ[, chrom_label := factor(chrom_label, levels = chrom_label)]
p_chr <- ggplot(chr_summ, aes(x = chrom_label, y = DUP_per_Mb)) +
  geom_col(fill = col_main) +
  labs(x = "Chromosome", y = "DUPs per Mb", title = "DUP density by chromosome") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(opt$plot_dir, "DUP_density_per_Mb.pdf"), p_chr, width = 10, height = 5, dpi = 400)

# Composite
comp <- plot_grid(
  p_freq, p_size, p_func, p_chr,
  labels = c("A","B","C","D"), ncol = 2
)
ggsave(file.path(opt$plot_dir, "composite_DUP_summary.pdf"), comp, width = 14, height = 10, dpi = 400)
