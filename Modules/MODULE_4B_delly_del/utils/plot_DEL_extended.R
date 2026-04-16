#!/usr/bin/env Rscript
# =============================================================================
# 11_plot_DEL_extended.R — Extended DEL analysis plotting suite
#
# Extends the basic 06_plot_delly_results.R with:
#   - Size-stratified views (small/medium/large)
#   - Ancestry-ordered heatmaps and colored PCA/burden
#   - Repeat-aware views
#   - Gene/functional class views
#   - QC relationship plots
#   - 81-only versions
#   - Frequency class charts
#   - Per-chromosome functional breakdown
#
# Usage:
#   Rscript 11_plot_DEL_extended.R \
#     --summary_dir DIR --final_dir DIR --plot_dir DIR \
#     --ref_fai PATH --samples_unrelated PATH \
#     [--qopt PATH --qopt_samples PATH]
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
  make_option("--samples_unrelated", type="character", default=""),
  make_option("--qopt", type="character", default=""),
  make_option("--qopt_samples", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$plot_dir, showWarnings=FALSE, recursive=TRUE)

# ── Theme ───────────────────────────────────────────────────────────────────
theme_pub <- theme_classic(base_size=10) +
  theme(plot.title=element_text(face="bold", size=12),
        axis.text=element_text(color="black"),
        legend.position="right", panel.grid=element_blank())

# Ancestry palette (K up to 12)
ancestry_pal <- c(
  Q1="#E41A1C", Q2="#377EB8", Q3="#4DAF4A", Q4="#984EA3",
  Q5="#FF7F00", Q6="#FFFF33", Q7="#A65628", Q8="#F781BF",
  Q9="#66C2A5", Q10="#FC8D62", Q11="#8DA0CB", Q12="#E78AC3"
)

col_main  <- "#2E8B8B"
col_light <- "#7ECACA"
col_dark  <- "#1A5252"

cat("=== Loading data ===\n")

# ── Load data ───────────────────────────────────────────────────────────────
master  <- fread(file.path(opt$summary_dir, "master_DEL_annotation_226.tsv"))
persamp <- fread(file.path(opt$summary_dir, "per_sample_DEL_summary.tsv"))
fai     <- fread(opt$ref_fai, header=FALSE,
                 col.names=c("chrom","length","offset","line_bases","line_bytes"))

# Frequency and size summaries
freq_summ <- tryCatch(fread(file.path(opt$summary_dir, "frequency_class_summary.tsv")), error=function(e) NULL)
size_summ <- tryCatch(fread(file.path(opt$summary_dir, "size_class_summary.tsv")), error=function(e) NULL)

# Unrelated sample list
unrelated_samples <- character(0)
if (file.exists(opt$samples_unrelated))
  unrelated_samples <- scan(opt$samples_unrelated, what="character", quiet=TRUE)

# Ancestry
has_ancestry <- FALSE
if (nchar(opt$qopt) > 0 && file.exists(opt$qopt) && nchar(opt$qopt_samples) > 0 && file.exists(opt$qopt_samples)) {
  qsamp <- scan(opt$qopt_samples, what="character", quiet=TRUE)
  qmat <- as.matrix(fread(opt$qopt, header=FALSE))
  if (nrow(qmat) == length(qsamp)) {
    K <- ncol(qmat)
    anc_dt <- data.table(sample=qsamp)
    for (k in 1:K) anc_dt[, paste0("Q", k) := qmat[, k]]
    anc_dt[, cluster := paste0("Q", apply(qmat, 1, which.max))]
    anc_dt[, max_Q := apply(qmat, 1, max)]
    has_ancestry <- TRUE
    cat("  Ancestry loaded: K=", K, ", samples=", nrow(anc_dt), "\n")
  }
}

cat("  Master annotation:", nrow(master), "DELs\n")
cat("  Per-sample summary:", nrow(persamp), "samples\n")

# =============================================================================
# 1. FREQUENCY CLASS CHART
# =============================================================================
cat("=== Frequency class chart ===\n")
freq_order <- c("singleton","rare","low_frequency","common","widespread","near_fixed")
if (!is.null(freq_summ)) {
  freq_summ[, frequency_class := factor(frequency_class, levels=freq_order)]
  p_freq <- ggplot(freq_summ, aes(x=frequency_class, y=count)) +
    geom_col(fill=col_main, width=0.7) +
    geom_text(aes(label=count), vjust=-0.3, size=3) +
    labs(x="Frequency class", y="DEL count", title="DEL frequency spectrum") +
    theme_pub + theme(axis.text.x=element_text(angle=30, hjust=1))
  ggsave(file.path(opt$plot_dir, "freq_class_barplot.pdf"), p_freq, width=8, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "freq_class_barplot.png"), p_freq, width=8, height=5, dpi=400)
}

# =============================================================================
# 2. SIZE-STRATIFIED BURDEN PLOTS
# =============================================================================
cat("=== Size-stratified views ===\n")
if ("burden_small" %in% names(persamp)) {
  ps_long <- melt(persamp, id.vars=c("sample","subset"),
                  measure.vars=c("burden_small","burden_medium","burden_large"),
                  variable.name="size_class", value.name="count")
  ps_long[, size_class := gsub("burden_", "", size_class)]
  ps_long[, size_class := factor(size_class, levels=c("small","medium","large"))]

  p_size_burden <- ggplot(ps_long, aes(x=reorder(sample, -count), y=count, fill=size_class)) +
    geom_col(position="stack", width=0.8) +
    scale_fill_manual(values=c(small="#7ECACA", medium="#2E8B8B", large="#1A5252"), name="Size class") +
    labs(x="Sample", y="DEL count", title="Per-sample DEL burden by size class") +
    theme_pub + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=2.5))
  ggsave(file.path(opt$plot_dir, "size_stratified_burden.pdf"), p_size_burden, width=16, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "size_stratified_burden.png"), p_size_burden, width=16, height=5, dpi=400)

  # Size class histogram from master
  p_size_hist <- ggplot(master, aes(x=factor(size_class, levels=c("sub_SV","small","medium","large")))) +
    geom_bar(fill=col_main) +
    geom_text(stat="count", aes(label=after_stat(count)), vjust=-0.3, size=3) +
    labs(x="Size class", y="Count", title="DEL count by size class") +
    theme_pub
  ggsave(file.path(opt$plot_dir, "size_class_barplot.pdf"), p_size_hist, width=6, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "size_class_barplot.png"), p_size_hist, width=6, height=5, dpi=400)
}

# =============================================================================
# 3. FUNCTIONAL CLASS VIEWS
# =============================================================================
cat("=== Functional class views ===\n")
func_order <- c("intergenic","intronic","exon_overlap","CDS_overlap")
if ("func_class" %in% names(master)) {
  master[, func_class_f := factor(func_class, levels=func_order)]

  p_func <- ggplot(master, aes(x=func_class_f)) +
    geom_bar(fill=col_main) +
    geom_text(stat="count", aes(label=after_stat(count)), vjust=-0.3, size=3) +
    labs(x="Functional class", y="Count", title="DEL count by functional class") +
    theme_pub
  ggsave(file.path(opt$plot_dir, "func_class_barplot.pdf"), p_func, width=7, height=5, dpi=400)

  # Functional class by chromosome
  func_chr <- master[, .N, by=.(chrom, func_class)]
  func_chr[, chrom_num := as.integer(gsub("\\D+", "", chrom))]
  func_chr <- func_chr[order(chrom_num)]
  func_chr[, chrom_label := gsub("C_gar_", "", chrom)]
  func_chr[, chrom_label := factor(chrom_label, levels=unique(chrom_label))]
  func_chr[, func_class := factor(func_class, levels=func_order)]

  p_func_chr <- ggplot(func_chr, aes(x=chrom_label, y=N, fill=func_class)) +
    geom_col(position="stack") +
    scale_fill_manual(values=c(intergenic="#B0BEC5", intronic="#7ECACA",
                               exon_overlap="#2E8B8B", CDS_overlap="#E69F00"),
                      name="Functional class") +
    labs(x="Chromosome", y="DEL count", title="Functional class by chromosome") +
    theme_pub + theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(file.path(opt$plot_dir, "func_class_by_chromosome.pdf"), p_func_chr, width=12, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "func_class_by_chromosome.png"), p_func_chr, width=12, height=5, dpi=400)

  # Functional class by size class
  if ("size_class" %in% names(master)) {
    func_size <- master[, .N, by=.(size_class, func_class)]
    func_size[, size_class := factor(size_class, levels=c("sub_SV","small","medium","large"))]
    func_size[, func_class := factor(func_class, levels=func_order)]

    p_func_size <- ggplot(func_size, aes(x=size_class, y=N, fill=func_class)) +
      geom_col(position="stack") +
      scale_fill_manual(values=c(intergenic="#B0BEC5", intronic="#7ECACA",
                                 exon_overlap="#2E8B8B", CDS_overlap="#E69F00"),
                        name="Functional class") +
      labs(x="Size class", y="Count", title="Functional class by size") +
      theme_pub
    ggsave(file.path(opt$plot_dir, "func_class_by_size.pdf"), p_func_size, width=7, height=5, dpi=400)
  }
}

# =============================================================================
# 4. REPEAT-AWARE VIEWS
# =============================================================================
cat("=== Repeat-aware views ===\n")
if ("repeat_50pct" %in% names(master)) {
  master[, repeat_status := ifelse(repeat_50pct == 1, ">=50% repeat", "non-repeat")]

  p_rep <- ggplot(master, aes(x=repeat_status)) +
    geom_bar(fill=col_main) +
    geom_text(stat="count", aes(label=after_stat(count)), vjust=-0.3, size=3) +
    labs(x="Repeat status", y="Count", title="DEL repeat overlap") +
    theme_pub
  ggsave(file.path(opt$plot_dir, "repeat_status_barplot.pdf"), p_rep, width=5, height=5, dpi=400)

  # Repeat by size class
  rep_size <- master[, .N, by=.(size_class, repeat_status)]
  rep_size[, size_class := factor(size_class, levels=c("sub_SV","small","medium","large"))]

  p_rep_size <- ggplot(rep_size, aes(x=size_class, y=N, fill=repeat_status)) +
    geom_col(position="dodge") +
    scale_fill_manual(values=c(">=50% repeat"="#E69F00", "non-repeat"=col_main)) +
    labs(x="Size class", y="Count", title="Repeat overlap by size") +
    theme_pub
  ggsave(file.path(opt$plot_dir, "repeat_by_size.pdf"), p_rep_size, width=7, height=5, dpi=400)
}

# =============================================================================
# 5. ANCESTRY-COLORED PCA
# =============================================================================
cat("=== Ancestry-colored PCA ===\n")
if (has_ancestry && "PC1" %in% names(persamp)) {
  pca_dt <- merge(persamp[, .(sample, PC1, PC2, subset)],
                  anc_dt[, .(sample, cluster, max_Q)], by="sample", all.x=TRUE)
  pca_dt[is.na(cluster), cluster := "unknown"]

  p_anc_pca <- ggplot(pca_dt, aes(x=PC1, y=PC2, color=cluster)) +
    geom_point(size=2.5, alpha=0.8) +
    scale_color_manual(values=ancestry_pal, name="Ancestry") +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") +
    geom_vline(xintercept=0, linetype="dashed", color="grey50") +
    labs(title="PCA colored by ancestry cluster") +
    theme_pub
  ggsave(file.path(opt$plot_dir, "PCA_ancestry_colored.pdf"), p_anc_pca, width=8, height=6, dpi=400)
  ggsave(file.path(opt$plot_dir, "PCA_ancestry_colored.png"), p_anc_pca, width=8, height=6, dpi=400)

  # Ancestry-colored burden plot
  burden_anc <- merge(persamp[, .(sample, total_DEL)],
                      anc_dt[, .(sample, cluster, max_Q)], by="sample", all.x=TRUE)
  burden_anc[is.na(cluster), cluster := "unknown"]
  burden_anc <- burden_anc[order(cluster, -max_Q)]
  burden_anc[, sample := factor(sample, levels=sample)]

  p_anc_burden <- ggplot(burden_anc, aes(x=sample, y=total_DEL, fill=cluster)) +
    geom_col(width=0.8) +
    scale_fill_manual(values=ancestry_pal, name="Ancestry") +
    labs(x="Sample (ancestry-ordered)", y="DEL count",
         title="Per-sample DEL burden ordered by ancestry") +
    theme_pub + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=2.5))
  ggsave(file.path(opt$plot_dir, "burden_ancestry_ordered.pdf"), p_anc_burden, width=16, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "burden_ancestry_ordered.png"), p_anc_burden, width=16, height=5, dpi=400)
}

# =============================================================================
# 6. ANCESTRY-ORDERED GENOME HEATMAP
# =============================================================================
cat("=== Ancestry-ordered genome heatmap ===\n")
if (has_ancestry) {
  gt_matrix <- fread(file.path(opt$final_dir, "catalog_226.DEL.GT_matrix.tsv"))
  samples <- colnames(gt_matrix)[6:ncol(gt_matrix)]
  gt_matrix[, window_1Mb := floor(POS / 1e6) * 1e6]

  # Ancestry sample order
  anc_order <- anc_dt[order(cluster, -max_Q), sample]
  anc_order <- anc_order[anc_order %in% samples]
  # Add any samples not in ancestry at the end
  remaining <- samples[!samples %in% anc_order]
  sample_order <- c(anc_order, remaining)

  # Per-sample per-window matrix
  wsl <- list()
  for (s in samples) {
    carried <- gt_matrix[get(s) %in% c("0/1","1/1","0|1","1|0","1|1"), .(CHROM, window_1Mb)]
    if (nrow(carried) > 0) {
      cts <- carried[, .N, by=.(CHROM, window_1Mb)]
      cts[, sample := s]
      wsl[[s]] <- cts
    }
  }
  wsd <- rbindlist(wsl)
  setnames(wsd, "N", "n_DELs")

  fai_dt <- data.table(chrom=fai$chrom, length=fai$length)
  all_windows <- fai_dt[, .(window_1Mb = seq(0, length - 1, by=1e6)), by=chrom]
  all_windows[, window_id := paste0(chrom, ":", window_1Mb)]
  col_order <- all_windows[, window_id]

  wide <- dcast(wsd, sample ~ paste0(CHROM, ":", window_1Mb), value.var="n_DELs", fill=0)
  wide <- wide[match(sample_order, sample)]
  mat <- as.matrix(wide[, -1, with=FALSE])
  rownames(mat) <- wide$sample
  avail_cols <- col_order[col_order %in% colnames(mat)]
  mat <- mat[, avail_cols, drop=FALSE]

  # Ancestry color bar
  sample_clusters <- anc_dt[match(sample_order, sample), cluster]
  sample_clusters[is.na(sample_clusters)] <- "unknown"
  anc_colors <- ancestry_pal[1:length(unique(sample_clusters))]
  names(anc_colors) <- sort(unique(sample_clusters))

  col_fun <- colorRamp2(c(0, 5, 10, 15), c("#440154", "#31688e", "#35b779", "#fde725"))

  left_anno <- rowAnnotation(
    Ancestry = sample_clusters,
    col = list(Ancestry = anc_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize=8)
  )

  pdf(file.path(opt$plot_dir, "genome_heatmap_ancestry_ordered.pdf"), width=18, height=12)
  ht <- Heatmap(mat,
    name="# DEL", col=col_fun,
    cluster_rows=FALSE, cluster_columns=FALSE,
    show_row_names=(nrow(mat) <= 80),
    show_column_names=FALSE,
    row_names_gp=gpar(fontsize=4),
    left_annotation=left_anno,
    column_title="Genome-wide DEL density (ancestry-ordered samples)",
    use_raster=TRUE, raster_quality=3
  )
  draw(ht)
  dev.off()
  cat("  Ancestry-ordered heatmap saved.\n")
}

# =============================================================================
# 7. QC RELATIONSHIP PLOTS
# =============================================================================
cat("=== QC relationship plots ===\n")
# These require per-sample missingness from the summary
if ("missingness" %in% names(persamp)) {
  p_miss <- ggplot(persamp, aes(x=missingness, y=total_DEL)) +
    geom_point(color=col_main, alpha=0.7, size=2) +
    geom_smooth(method="lm", se=TRUE, color="black", linewidth=0.5) +
    labs(x="Genotype missingness", y="DEL count",
         title="DEL burden vs missingness") +
    theme_pub
  ggsave(file.path(opt$plot_dir, "qc_burden_vs_missingness.pdf"), p_miss, width=6, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "qc_burden_vs_missingness.png"), p_miss, width=6, height=5, dpi=400)
}

# =============================================================================
# 8. 81-ONLY VIEWS
# =============================================================================
cat("=== 81-only subset views ===\n")
if (length(unrelated_samples) > 0) {
  ps81 <- persamp[sample %in% unrelated_samples]

  if (nrow(ps81) > 0) {
    # Burden plot 81
    p_b81 <- ggplot(ps81[order(-total_DEL)], aes(x=reorder(sample, -total_DEL), y=total_DEL)) +
      geom_col(fill=col_main, width=0.8) +
      labs(x="Sample", y="DEL count", title="Per-sample DEL burden (81 unrelated)") +
      theme_pub + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=4))
    ggsave(file.path(opt$plot_dir, "81_burden.pdf"), p_b81, width=14, height=5, dpi=400)
    ggsave(file.path(opt$plot_dir, "81_burden.png"), p_b81, width=14, height=5, dpi=400)

    # PCA 81
    if ("PC1" %in% names(ps81)) {
      p_pca81 <- ggplot(ps81, aes(x=PC1, y=PC2))
      if (has_ancestry) {
        pca81m <- merge(ps81, anc_dt[, .(sample, cluster)], by="sample", all.x=TRUE)
        pca81m[is.na(cluster), cluster := "unknown"]
        p_pca81 <- ggplot(pca81m, aes(x=PC1, y=PC2, color=cluster)) +
          scale_color_manual(values=ancestry_pal, name="Ancestry")
      }
      p_pca81 <- p_pca81 + geom_point(size=3, alpha=0.8) +
        geom_hline(yintercept=0, linetype="dashed", color="grey50") +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        labs(title="PCA (81 unrelated)") + theme_pub
      ggsave(file.path(opt$plot_dir, "81_PCA.pdf"), p_pca81, width=8, height=6, dpi=400)
      ggsave(file.path(opt$plot_dir, "81_PCA.png"), p_pca81, width=8, height=6, dpi=400)
    }
  }
}

# =============================================================================
# 9. MARKER SUMMARY PLOTS
# =============================================================================
cat("=== Marker summary plots ===\n")
marker_files <- c(
  tier1 = file.path(opt$summary_dir, "markers_tier1_strict.tsv"),
  tier2 = file.path(opt$summary_dir, "markers_tier2_gene.tsv"),
  tier3 = file.path(opt$summary_dir, "markers_tier3_group_informative.tsv")
)
marker_counts <- data.table(
  tier = c("Tier 1: Strict", "Tier 2: Gene", "Tier 3: Group"),
  count = sapply(marker_files, function(f) {
    if (file.exists(f)) nrow(fread(f)) - 0 else 0  # fread auto-skips header
  })
)

if (any(marker_counts$count > 0)) {
  p_markers <- ggplot(marker_counts, aes(x=tier, y=count)) +
    geom_col(fill=c("#2E8B8B", "#E69F00", "#984EA3"), width=0.6) +
    geom_text(aes(label=count), vjust=-0.3, size=4) +
    labs(x="Marker tier", y="Count", title="DEL marker selection summary") +
    theme_pub
  ggsave(file.path(opt$plot_dir, "marker_tier_summary.pdf"), p_markers, width=6, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "marker_tier_summary.png"), p_markers, width=6, height=5, dpi=400)

  # Tier 1 size distribution
  if (file.exists(marker_files["tier1"])) {
    t1 <- fread(marker_files["tier1"])
    if (nrow(t1) > 0 && "svlen" %in% names(t1)) {
      p_t1_size <- ggplot(t1, aes(x=svlen)) +
        geom_histogram(bins=50, fill=col_main, color="white", linewidth=0.2) +
        scale_x_log10(labels=scales::comma) +
        annotation_logticks(sides="b") +
        labs(x="DEL size (bp, log)", y="Count",
             title=paste0("Tier 1 marker size distribution (n=", nrow(t1), ")")) +
        theme_pub
      ggsave(file.path(opt$plot_dir, "marker_tier1_size.pdf"), p_t1_size, width=7, height=5, dpi=400)
    }
  }
}

# =============================================================================
# 10. PER-CHROMOSOME FUNCTIONAL SUMMARY (from 10_gene_summary_tables output)
# =============================================================================
cat("=== Per-chromosome functional summary ===\n")
chr_summ_file <- file.path(opt$summary_dir, "per_chromosome_summary.tsv")
if (file.exists(chr_summ_file)) {
  chr_summ <- fread(chr_summ_file)
  chr_summ[, chrom_label := gsub("C_gar_", "", chrom)]
  chr_summ[, chrom_label := factor(chrom_label, levels=chrom_label)]

  p_del_mb <- ggplot(chr_summ, aes(x=chrom_label, y=DEL_per_Mb)) +
    geom_col(fill=col_main, width=0.7) +
    labs(x="Chromosome", y="DELs per Mb", title="DEL density by chromosome") +
    theme_pub + theme(axis.text.x=element_text(angle=45, hjust=1))
  ggsave(file.path(opt$plot_dir, "DEL_density_per_Mb.pdf"), p_del_mb, width=10, height=5, dpi=400)
  ggsave(file.path(opt$plot_dir, "DEL_density_per_Mb.png"), p_del_mb, width=10, height=5, dpi=400)
}

# =============================================================================
# COMPOSITE FIGURE 2: Size + Frequency + Repeat + Function
# =============================================================================
cat("=== Composite figure 2 ===\n")
plots_fig2 <- list()
if (exists("p_size_hist"))  plots_fig2[["A"]] <- p_size_hist + ggtitle("A")
if (exists("p_freq") && !is.null(freq_summ)) plots_fig2[["B"]] <- p_freq + ggtitle("B")
if (exists("p_rep"))        plots_fig2[["C"]] <- p_rep + ggtitle("C")
if (exists("p_func"))       plots_fig2[["D"]] <- p_func + ggtitle("D")

if (length(plots_fig2) == 4) {
  comp2 <- plot_grid(plotlist=plots_fig2, ncol=2, nrow=2)
  ggsave(file.path(opt$plot_dir, "composite_fig2_size_freq_repeat_func.pdf"),
         comp2, width=14, height=10, dpi=400)
  ggsave(file.path(opt$plot_dir, "composite_fig2_size_freq_repeat_func.png"),
         comp2, width=14, height=10, dpi=400)
}

cat("\n=== Extended plotting complete ===\n")
cat("Output:", opt$plot_dir, "\n")
