#!/usr/bin/env Rscript
# =============================================================================
# 06_plot_bnd_results.R
# Publication-style BND catalog summary plots
#
# Panels:
#   A — Genome-wide BND density heatmap (1-Mb windows, per sample)
#   B — Per-sample BND count barplot
#   C — PCA on binary BND genotype matrix
#   D — Pairwise shared BND heatmap
#   F — Per-chromosome BND count barplot
#
# Notes:
#   - A and D are written as both PDF and PNG
#   - Composite is explicitly BCF only
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
  library(grid)
})

option_list <- list(
  make_option("--summary_dir", type = "character"),
  make_option("--final_dir",   type = "character"),
  make_option("--plot_dir",    type = "character"),
  make_option("--ref_fai",     type = "character"),
  make_option("--samples_unrelated", type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$plot_dir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    panel.grid = element_blank()
  )

col_main  <- "#2E8B8B"
col_light <- "#7ECACA"
col_dark  <- "#1A5252"

cat("=== Loading data ===\n")

per_sample <- fread(file.path(opt$summary_dir, "per_sample_BND_counts.tsv"))
per_chr    <- fread(file.path(opt$summary_dir, "per_chromosome_BND_counts.tsv"))
pairwise   <- fread(file.path(opt$summary_dir, "pairwise_shared_BND.tsv"))
binary_gt  <- fread(file.path(opt$summary_dir, "BND_binary_genotype_matrix.tsv"))
fai        <- fread(opt$ref_fai, header = FALSE,
                    col.names = c("chrom","length","offset","line_bases","line_bytes"))

unrelated_samples <- character(0)
if (file.exists(opt$samples_unrelated)) {
  unrelated_samples <- scan(opt$samples_unrelated, what = "character", quiet = TRUE)
}

setnames(per_sample, names(per_sample)[2], "n_SVs")
setnames(per_chr, names(per_chr)[2], "n_SVs")

cat("  Samples:", nrow(per_sample), "\n")
cat("  BND sites:", nrow(binary_gt), "\n")

sample_order <- per_sample[order(-n_SVs), sample]

# =============================================================================
# PANEL A: Genome-wide BND density heatmap per sample
# =============================================================================
cat("=== Panel A: Genome-wide heatmap ===\n")

gt_matrix_file <- file.path(opt$final_dir, "catalog_226.BND.GT_matrix.tsv")
has_panel_a <- file.exists(gt_matrix_file)

if (has_panel_a) {
  gt_matrix <- fread(gt_matrix_file)

  if (nrow(gt_matrix) > 0 && ncol(gt_matrix) >= 6) {
    samples <- colnames(gt_matrix)[6:ncol(gt_matrix)]
    gt_matrix[, window_1Mb := floor(POS / 1e6) * 1e6]

    fai_dt <- data.table(chrom = fai$chrom, length = fai$length)
    fai_dt[, cum_offset := cumsum(c(0, length[-.N]))]
    fai_dt[, cum_end := cum_offset + length]

    window_sample_list <- list()
    for (s in samples) {
      carried <- gt_matrix[get(s) %in% c("0/1","1/1","0|1","1|0","1|1"),
                           .(CHROM, window_1Mb)]
      if (nrow(carried) > 0) {
        cts <- carried[, .N, by = .(CHROM, window_1Mb)]
        cts[, sample := s]
        window_sample_list[[s]] <- cts
      }
    }

    if (length(window_sample_list) > 0) {
      window_sample_dt <- rbindlist(window_sample_list, fill = TRUE)
      setnames(window_sample_dt, "N", "n_SVs")

      window_sample_dt <- merge(
        window_sample_dt,
        fai_dt[, .(chrom, cum_offset)],
        by.x = "CHROM",
        by.y = "chrom",
        all.x = TRUE
      )
      window_sample_dt[, cum_pos := cum_offset + window_1Mb + 5e5]

      all_windows <- fai_dt[, .(window_1Mb = seq(0, length - 1, by = 1e6)), by = chrom]
      all_windows <- merge(all_windows, fai_dt[, .(chrom, cum_offset)], by = "chrom")
      all_windows[, window_id := paste0(chrom, ":", window_1Mb)]

      wide_dt <- dcast(
        window_sample_dt,
        sample ~ paste0(CHROM, ":", window_1Mb),
        value.var = "n_SVs",
        fill = 0
      )
      wide_dt <- wide_dt[match(sample_order, sample)]

      mat_a <- as.matrix(wide_dt[, -1, with = FALSE])
      rownames(mat_a) <- wide_dt$sample
      col_order <- all_windows[order(cum_offset, window_1Mb), window_id]
      col_order <- col_order[col_order %in% colnames(mat_a)]
      mat_a <- mat_a[, col_order, drop = FALSE]

      cum_counts <- colSums(mat_a)

      col_fun_a <- colorRamp2(c(0, 1, 2, 4, 8), c("#440154", "#31688e", "#35b779", "#90d743", "#fde725"))

      top_anno <- HeatmapAnnotation(
        cumulative = anno_barplot(
          cum_counts,
          height = unit(2, "cm"),
          gp = gpar(fill = col_dark, col = NA),
          axis_param = list(gp = gpar(fontsize = 6))
        ),
        annotation_name_side = "left"
      )

      right_anno <- rowAnnotation(
        total = anno_barplot(
          per_sample[match(sample_order, sample), n_SVs],
          width = unit(2, "cm"),
          gp = gpar(fill = col_main, col = NA),
          axis_param = list(gp = gpar(fontsize = 6))
        ),
        annotation_name_side = "top"
      )

      make_ht_a <- function() {
        Heatmap(
          mat_a,
          name = "# BND",
          col = col_fun_a,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = (nrow(mat_a) <= 80),
          show_column_names = FALSE,
          row_names_gp = gpar(fontsize = 5),
          top_annotation = top_anno,
          right_annotation = right_anno,
          column_title = "A  Genome-wide BND density (1-Mb windows)",
          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          row_title = "Samples",
          use_raster = TRUE,
          raster_quality = 3,
          heatmap_legend_param = list(title = "# BND", title_gp = gpar(fontsize = 8))
        )
      }

      pdf(file.path(opt$plot_dir, "panel_A_genome_heatmap.pdf"), width = 18, height = 10)
      draw(make_ht_a())
      dev.off()

      png(file.path(opt$plot_dir, "panel_A_genome_heatmap.png"), width = 18, height = 10,
          units = "in", res = 400)
      draw(make_ht_a())
      dev.off()

      cat("  Panel A saved.\n")
    } else {
      cat("  Skipping Panel A: no carried BND calls in GT matrix.\n")
      has_panel_a <- FALSE
    }
  } else {
    cat("  Skipping Panel A: BND GT matrix is empty.\n")
    has_panel_a <- FALSE
  }
} else {
  cat("  Skipping Panel A: catalog_226.BND.GT_matrix.tsv not found.\n")
}

# =============================================================================
# PANEL B: Per-sample BND count barplot
# =============================================================================
cat("=== Panel B: Per-sample barplot ===\n")

per_sample[, source := ifelse(sample %in% unrelated_samples, "Unrelated (81)", "Related")]
per_sample[, sample := factor(sample, levels = sample_order)]

p_b <- ggplot(per_sample, aes(x = sample, y = n_SVs, fill = source)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("Unrelated (81)" = col_main, "Related" = col_light), name = "Subset") +
  labs(x = "Sample", y = "BND count", title = "B") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave(file.path(opt$plot_dir, "panel_B_per_sample_BND.pdf"), p_b, width = 16, height = 5, dpi = 400)
ggsave(file.path(opt$plot_dir, "panel_B_per_sample_BND.png"), p_b, width = 16, height = 5, dpi = 400)
cat("  Panel B saved.\n")

# =============================================================================
# PANEL C: PCA on binary BND genotype matrix
# =============================================================================
cat("=== Panel C: PCA ===\n")

gt_for_pca <- as.matrix(binary_gt[, -1, with = FALSE])
pca_input <- t(gt_for_pca)
rownames(pca_input) <- colnames(binary_gt)[-1]
col_var <- apply(pca_input, 2, var)
pca_input <- pca_input[, col_var > 0, drop = FALSE]

cat("  PCA input:", nrow(pca_input), "samples x", ncol(pca_input), "BNDs\n")

has_panel_c <- FALSE
if (ncol(pca_input) > 1) {
  pca_res <- prcomp(pca_input, center = TRUE, scale. = FALSE)
  pca_df <- data.table(
    sample = rownames(pca_input),
    PC1 = pca_res$x[,1],
    PC2 = pca_res$x[,2]
  )
  pca_df[, source := ifelse(sample %in% unrelated_samples, "Unrelated (81)", "Related")]
  var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  p_c <- ggplot(pca_df, aes(x = PC1, y = PC2, color = source)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = c("Unrelated (81)" = col_main, "Related" = col_light), name = "Subset") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = paste0("PC1 (", var_pct[1], "%)"),
         y = paste0("PC2 (", var_pct[2], "%)"),
         title = "C") +
    theme_pub +
    theme(plot.title = element_text(face = "bold", size = 14))

  ggsave(file.path(opt$plot_dir, "panel_C_PCA_BND.pdf"), p_c, width = 7, height = 6, dpi = 400)
  ggsave(file.path(opt$plot_dir, "panel_C_PCA_BND.png"), p_c, width = 7, height = 6, dpi = 400)
  cat("  Panel C saved.\n")
  has_panel_c <- TRUE
} else {
  cat("  Skipping PCA: insufficient variant sites.\n")
}

# =============================================================================
# PANEL D: Pairwise shared BND heatmap
# =============================================================================
cat("=== Panel D: Pairwise shared heatmap ===\n")

pw_wide <- dcast(pairwise, sample_i ~ sample_j, value.var = "n_shared", fill = 0)
pw_samples <- as.character(pw_wide$sample_i)
pw_mat <- as.matrix(pw_wide[, -1, with = FALSE])
rownames(pw_mat) <- pw_samples

for (i in seq_along(pw_samples)) {
  for (j in seq_along(pw_samples)) {
    if (pw_mat[i,j] == 0 && pw_mat[j,i] > 0) pw_mat[i,j] <- pw_mat[j,i]
  }
}

hc <- hclust(as.dist(max(pw_mat) - pw_mat), method = "ward.D2")
col_fun_d <- colorRamp2(
  quantile(pw_mat[upper.tri(pw_mat)], c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
  c("#440154", "#31688e", "#21918c", "#5ec962", "#fde725")
)

source_colors <- c("Unrelated (81)" = col_main, "Related" = col_light)
sample_source <- ifelse(pw_samples %in% unrelated_samples, "Unrelated (81)", "Related")
row_anno <- rowAnnotation(Source = sample_source, col = list(Source = source_colors),
                          show_annotation_name = FALSE)

make_ht_d <- function() {
  Heatmap(
    pw_mat,
    name = "# Shared\nBND",
    col = col_fun_d,
    cluster_rows = hc, cluster_columns = hc,
    show_row_names = (nrow(pw_mat) <= 80),
    show_column_names = (ncol(pw_mat) <= 80),
    row_names_gp = gpar(fontsize = 4),
    column_names_gp = gpar(fontsize = 4),
    left_annotation = row_anno,
    column_title = "D  Pairwise shared BNDs",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    use_raster = TRUE, raster_quality = 3
  )
}

pdf(file.path(opt$plot_dir, "panel_D_pairwise_shared.pdf"), width = 12, height = 11)
draw(make_ht_d())
dev.off()

png(file.path(opt$plot_dir, "panel_D_pairwise_shared.png"), width = 12, height = 11,
    units = "in", res = 400)
draw(make_ht_d())
dev.off()

cat("  Panel D saved.\n")

# =============================================================================
# PANEL F: Per-chromosome BND counts
# =============================================================================
cat("=== Panel F: Per-chromosome ===\n")

per_chr[, chrom_num := suppressWarnings(as.integer(gsub("\\D+", "", chrom)))]
per_chr <- per_chr[order(chrom_num, chrom)]
per_chr[, chrom_label := gsub("C_gar_", "", chrom)]
per_chr[, chrom_label := factor(chrom_label, levels = chrom_label)]

p_f <- ggplot(per_chr, aes(x = chrom_label, y = n_SVs)) +
  geom_col(fill = col_main, width = 0.7) +
  labs(x = "Chromosome", y = "BND count", title = "F") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave(file.path(opt$plot_dir, "panel_F_per_chromosome.pdf"), p_f, width = 10, height = 5, dpi = 400)
ggsave(file.path(opt$plot_dir, "panel_F_per_chromosome.png"), p_f, width = 10, height = 5, dpi = 400)
cat("  Panel F saved.\n")

# =============================================================================
# COMPOSITE (BCF only)
# =============================================================================
cat("=== Composite figure (BCF only) ===\n")

plots_list <- list(p_b + theme(legend.position = "none"))
if (has_panel_c) plots_list <- c(plots_list, list(p_c + theme(legend.position = "none")))
plots_list <- c(plots_list, list(p_f))

composite <- plot_grid(plotlist = plots_list, ncol = 2)
ggsave(file.path(opt$plot_dir, "composite_BND_BCF.pdf"), composite, width = 16, height = 8, dpi = 400)
ggsave(file.path(opt$plot_dir, "composite_BND_BCF.png"), composite, width = 16, height = 8, dpi = 400)
cat("  Composite BCF saved.\n")

legend_file <- file.path(opt$plot_dir, "figure_legend.txt")
legend_lines <- c(
  "Figure X. Breakend variants in catfish hatchery cohort.",
  "",
  paste0("(A) Genome-wide distribution of breakend variants across ", nrow(per_sample),
         " samples. Heatmap shows BND count per non-overlapping 1-Mb window."),
  paste0("(B) Per-sample BND counts, sorted by total. ",
         "Teal = unrelated subset (n=", length(unrelated_samples), ")."),
  if (has_panel_c) "(C) PCA based on binary BND genotype matrix." else "(C) PCA omitted due to insufficient variable sites.",
  "(D) Pairwise shared BND heatmap with hierarchical clustering (Ward's method).",
  "(F) BND counts per chromosome/linkage group."
)
writeLines(legend_lines, legend_file)

cat("\nAll plots written to:", opt$plot_dir, "\n")
cat("Done.\n")
