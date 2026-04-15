#!/usr/bin/env Rscript
# =============================================================================
# 06_plot_inv_results.R
# Publication-style INV catalog summary plots
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(stats)
})

option_list <- list(
  make_option("--summary_dir", type="character"),
  make_option("--final_dir",   type="character"),
  make_option("--plot_dir",    type="character"),
  make_option("--ref_fai",     type="character"),
  make_option("--samples_unrelated", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$plot_dir, showWarnings=FALSE, recursive=TRUE)

theme_pub <- theme_classic(base_size=10) +
  theme(plot.title = element_text(face="bold", size=12),
        axis.text = element_text(color="black"),
        legend.position = "right", panel.grid = element_blank())

col_main  <- "#2E8B8B"
col_light <- "#7ECACA"
col_dark  <- "#1A5252"

cat("=== Loading data ===\n")

per_sample   <- fread(file.path(opt$summary_dir, "per_sample_INV_counts.tsv"))
per_chr      <- fread(file.path(opt$summary_dir, "per_chromosome_INV_counts.tsv"))
pairwise     <- fread(file.path(opt$summary_dir, "pairwise_shared_INV.tsv"))
binary_gt    <- fread(file.path(opt$summary_dir, "INV_binary_genotype_matrix.tsv"))
fai          <- fread(opt$ref_fai, header=FALSE,
                      col.names=c("chrom","length","offset","line_bases","line_bytes"))

unrelated_samples <- character(0)
if (file.exists(opt$samples_unrelated))
  unrelated_samples <- scan(opt$samples_unrelated, what="character", quiet=TRUE)

# Standardize count column name
count_col <- names(per_sample)[2]
setnames(per_sample, count_col, "n_SVs")

cat("  Samples:", nrow(per_sample), "\n")
cat("  INV sites:", nrow(binary_gt), "\n")

sample_order <- per_sample[order(-n_SVs), sample]

# =============================================================================
# PANEL B: Per-sample INV count barplot
# =============================================================================
cat("=== Panel B: Per-sample barplot ===\n")

per_sample[, source := ifelse(sample %in% unrelated_samples, "Unrelated (81)", "Related")]
per_sample[, sample := factor(sample, levels=sample_order)]

p_b <- ggplot(per_sample, aes(x=sample, y=n_SVs, fill=source)) +
  geom_col(width=0.8) +
  scale_fill_manual(values=c("Unrelated (81)"=col_main, "Related"=col_light), name="Subset") +
  labs(x="Sample", y="INV count", title="B") +
  theme_pub +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=3),
        plot.title = element_text(face="bold", size=14))

ggsave(file.path(opt$plot_dir, "panel_B_per_sample_INV.pdf"), p_b, width=16, height=5, dpi=400)
ggsave(file.path(opt$plot_dir, "panel_B_per_sample_INV.png"), p_b, width=16, height=5, dpi=400)
cat("  Panel B saved.\n")

# =============================================================================
# PANEL C: PCA on binary INV genotype matrix
# =============================================================================
cat("=== Panel C: PCA ===\n")

gt_for_pca <- as.matrix(binary_gt[, -1, with=FALSE])
pca_input <- t(gt_for_pca)
rownames(pca_input) <- colnames(binary_gt)[-1]
col_var <- apply(pca_input, 2, var)
pca_input <- pca_input[, col_var > 0, drop=FALSE]

cat("  PCA input:", nrow(pca_input), "samples x", ncol(pca_input), "INVs\n")

if (ncol(pca_input) > 1) {
  pca_res <- prcomp(pca_input, center=TRUE, scale.=FALSE)
  pca_df <- data.table(
    sample = rownames(pca_input),
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2]
  )
  pca_df[, source := ifelse(sample %in% unrelated_samples, "Unrelated (81)", "Related")]
  var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  p_c <- ggplot(pca_df, aes(x=PC1, y=PC2, color=source)) +
    geom_point(size=2.5, alpha=0.8) +
    scale_color_manual(values=c("Unrelated (81)"=col_main, "Related"=col_light), name="Subset") +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") +
    geom_vline(xintercept=0, linetype="dashed", color="grey50") +
    labs(x=paste0("PC1 (", var_pct[1], "%)"),
         y=paste0("PC2 (", var_pct[2], "%)"),
         title="C") +
    theme_pub + theme(plot.title = element_text(face="bold", size=14))
  
  ggsave(file.path(opt$plot_dir, "panel_C_PCA_INV.pdf"), p_c, width=7, height=6, dpi=400)
  ggsave(file.path(opt$plot_dir, "panel_C_PCA_INV.png"), p_c, width=7, height=6, dpi=400)
  cat("  Panel C saved.\n")
} else {
  cat("  Skipping PCA: insufficient variant sites.\n")
}

# =============================================================================
# PANEL D: Pairwise shared INV heatmap
# =============================================================================
cat("=== Panel D: Pairwise shared heatmap ===\n")

pw_wide <- dcast(pairwise, sample_i ~ sample_j, value.var="n_shared", fill=0)
pw_samples <- as.character(pw_wide$sample_i)
pw_mat <- as.matrix(pw_wide[, -1, with=FALSE])
rownames(pw_mat) <- pw_samples
for (i in seq_along(pw_samples))
  for (j in seq_along(pw_samples))
    if (pw_mat[i,j] == 0 && pw_mat[j,i] > 0) pw_mat[i,j] <- pw_mat[j,i]

hc <- hclust(as.dist(max(pw_mat) - pw_mat), method="ward.D2")
col_fun_d <- colorRamp2(
  quantile(pw_mat[upper.tri(pw_mat)], c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE),
  c("#440154", "#31688e", "#21918c", "#5ec962", "#fde725")
)

source_colors <- c("Unrelated (81)"=col_main, "Related"=col_light)
sample_source <- ifelse(pw_samples %in% unrelated_samples, "Unrelated (81)", "Related")
row_anno <- rowAnnotation(Source = sample_source, col = list(Source = source_colors),
                          show_annotation_name = FALSE)

pdf(file.path(opt$plot_dir, "panel_D_pairwise_shared.pdf"), width=12, height=11)
ht_d <- Heatmap(pw_mat, name="# Shared\nINV", col=col_fun_d,
  cluster_rows=hc, cluster_columns=hc,
  show_row_names=(nrow(pw_mat) <= 80),
  show_column_names=(ncol(pw_mat) <= 80),
  row_names_gp=gpar(fontsize=4), column_names_gp=gpar(fontsize=4),
  left_annotation=row_anno,
  column_title="D  Pairwise shared INVs",
  column_title_gp=gpar(fontsize=12, fontface="bold"),
  use_raster=TRUE, raster_quality=3)
draw(ht_d)
dev.off()
cat("  Panel D saved.\n")

# =============================================================================
# PANEL E: INV size distribution
# =============================================================================
cat("=== Panel E: INV size distribution ===\n")

svlen_file <- file.path(opt$summary_dir, "INV_svlen_distribution.tsv")
if (file.exists(svlen_file)) {
  svlen_dist <- fread(svlen_file)
  svlen_dist[, svlen_kb := svlen_bp / 1000]
  
  p_e <- ggplot(svlen_dist, aes(x=svlen_bp)) +
    geom_histogram(bins=80, fill=col_main, color="white", linewidth=0.2) +
    scale_x_log10(labels=scales::comma) +
    annotation_logticks(sides="b") +
    labs(x="INV size (bp, log scale)", y="Count", title="E") +
    theme_pub +
    theme(plot.title = element_text(face="bold", size=14))
  
  ggsave(file.path(opt$plot_dir, "panel_E_INV_size_distribution.pdf"), p_e,
         width=7, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "panel_E_INV_size_distribution.png"), p_e,
         width=7, height=5, dpi=400)
  cat("  Panel E saved.\n")
}

# =============================================================================
# PANEL F: Per-chromosome INV counts
# =============================================================================
cat("=== Panel F: Per-chromosome ===\n")

per_chr[, chrom_num := as.integer(gsub("\\D+", "", chrom))]
per_chr <- per_chr[order(chrom_num)]
per_chr[, chrom_label := gsub("C_gar_", "", chrom)]
per_chr[, chrom_label := factor(chrom_label, levels=chrom_label)]

count_col_chr <- names(per_chr)[2]
setnames(per_chr, count_col_chr, "n_SVs")

p_f <- ggplot(per_chr, aes(x=chrom_label, y=n_SVs)) +
  geom_col(fill=col_main, width=0.7) +
  labs(x="Chromosome", y="INV count", title="F") +
  theme_pub +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(face="bold", size=14))

ggsave(file.path(opt$plot_dir, "panel_F_per_chromosome.pdf"), p_f, width=10, height=5, dpi=400)
ggsave(file.path(opt$plot_dir, "panel_F_per_chromosome.png"), p_f, width=10, height=5, dpi=400)
cat("  Panel F saved.\n")

# =============================================================================
# COMPOSITE
# =============================================================================
cat("=== Composite figure ===\n")
plots_list <- list(p_b + theme(legend.position="none"), p_f)
if (exists("p_c")) plots_list <- c(list(p_b + theme(legend.position="none")), list(p_c + theme(legend.position="none")), list(p_f))
if (exists("p_e")) plots_list <- c(plots_list, list(p_e))

composite <- plot_grid(plotlist=plots_list, ncol=2)
ggsave(file.path(opt$plot_dir, "composite_INV.pdf"), composite, width=16, height=8, dpi=400)
ggsave(file.path(opt$plot_dir, "composite_INV.png"), composite, width=16, height=8, dpi=400)
cat("  Composite saved.\n")

cat("\nAll plots written to:", opt$plot_dir, "\n")
cat("Done.\n")
