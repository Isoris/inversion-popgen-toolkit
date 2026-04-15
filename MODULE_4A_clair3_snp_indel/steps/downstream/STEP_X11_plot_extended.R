#!/usr/bin/env Rscript
# =============================================================================
# 11_plot_extended.R — Extended Clair3 variant analysis plots
#
# Extends 10_plot_main.R with:
#   - Frequency class barplots
#   - Ancestry-colored PCA and burden
#   - Burden by frequency class (stacked)
#   - QC: burden vs missingness
#   - 81-only subset views
#   - Marker tier summary
#   - Rare sharing network summary
#   - Per-sample window density heatmap (ComplexHeatmap)
#
# Usage:
#   Rscript 11_plot_extended.R \
#     --ds_chrom_dir DIR --chrom CHROM --ref_fai PATH \
#     --samples_unrelated PATH \
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
  make_option("--ds_chrom_dir", type="character"),
  make_option("--chrom",        type="character"),
  make_option("--ref_fai",      type="character"),
  make_option("--samples_unrelated", type="character", default=""),
  make_option("--qopt",         type="character", default=""),
  make_option("--qopt_samples", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option_list))

plot_dir <- file.path(opt$ds_chrom_dir, "figures")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

theme_pub <- theme_classic(base_size=10) +
  theme(plot.title=element_text(face="bold", size=12),
        axis.text=element_text(color="black"),
        legend.position="right", panel.grid=element_blank())

ancestry_pal <- c(Q1="#E41A1C", Q2="#377EB8", Q3="#4DAF4A", Q4="#984EA3",
                   Q5="#FF7F00", Q6="#FFFF33", Q7="#A65628", Q8="#F781BF")
col_snp    <- "#2E8B8B"
col_indel  <- "#E69F00"
col_combined <- "#5B4A8A"
col_light  <- "#7ECACA"
col_dark   <- "#1A5252"
type_colors <- c(SNP=col_snp, INDEL=col_indel, COMBINED=col_combined)

freq_class_colors <- c(monomorphic="#BDBDBD", singleton="#FDD835",
                        rare="#FB8C00", low_freq="#E53935",
                        common="#2E8B8B", high_freq="#1A5252")

# ── Load ancestry ──
unrelated_samples <- character(0)
if (file.exists(opt$samples_unrelated))
  unrelated_samples <- scan(opt$samples_unrelated, what="character", quiet=TRUE)

has_ancestry <- FALSE
anc_dt <- NULL
if (nchar(opt$qopt) > 0 && file.exists(opt$qopt) &&
    nchar(opt$qopt_samples) > 0 && file.exists(opt$qopt_samples)) {
  qsamp <- scan(opt$qopt_samples, what="character", quiet=TRUE)
  qmat <- as.matrix(fread(opt$qopt, header=FALSE))
  if (nrow(qmat) == length(qsamp)) {
    K <- ncol(qmat)
    anc_dt <- data.table(sample=qsamp)
    anc_dt[, cluster := paste0("Q", apply(qmat, 1, which.max))]
    anc_dt[, max_Q := apply(qmat, 1, max)]
    has_ancestry <- TRUE
    cat("Ancestry loaded: K=", K, "\n")
  }
}

cat("=== Extended plots ===\n")

for (vtype in c("SNP", "INDEL", "COMBINED")) {
  cat("\n=== ", vtype, " ===\n")
  vcol <- type_colors[vtype]
  prefix <- tolower(vtype)

  ps_file  <- file.path(opt$ds_chrom_dir, "per_sample", paste0("per_sample_", vtype, "_summary.tsv"))
  ann_file <- file.path(opt$ds_chrom_dir, "annotation", paste0("master_", vtype, "_annotation.tsv"))
  fq_file  <- file.path(opt$ds_chrom_dir, "annotation", paste0("frequency_class_summary.", vtype, ".tsv"))
  wd_file  <- file.path(opt$ds_chrom_dir, "distances", paste0("window_density.", vtype, ".tsv"))
  rn_file  <- file.path(opt$ds_chrom_dir, "distances", paste0("rare_sharing_nodes.", vtype, ".tsv"))
  mk_file  <- file.path(opt$ds_chrom_dir, "markers", paste0("markers_tier1_strict.", vtype, ".tsv"))

  per_samp <- tryCatch(fread(ps_file), error=function(e) NULL)
  master   <- tryCatch(fread(ann_file), error=function(e) NULL)
  freq_sum <- tryCatch(fread(fq_file), error=function(e) NULL)

  # ═════════════════════════════════════════════════════════════════════════
  # 1. FREQUENCY CLASS BARPLOT
  # ═════════════════════════════════════════════════════════════════════════
  if (!is.null(freq_sum) && nrow(freq_sum) > 0) {
    cat("  Frequency class barplot\n")
    fc_order <- c("singleton","rare","low_freq","common","high_freq")
    freq_sum <- freq_sum[frequency_class %in% fc_order]
    freq_sum[, frequency_class := factor(frequency_class, levels=fc_order)]

    p_fc <- ggplot(freq_sum, aes(x=frequency_class, y=count)) +
      geom_col(fill=vcol, width=0.7) +
      geom_text(aes(label=count), vjust=-0.3, size=3) +
      labs(x="Frequency class", y="Count",
           title=paste0(vtype, " frequency class spectrum")) +
      theme_pub + theme(axis.text.x=element_text(angle=30, hjust=1))

    ggsave(file.path(plot_dir, paste0(prefix, "_freq_class.pdf")), p_fc, width=7, height=5, dpi=400)
    ggsave(file.path(plot_dir, paste0(prefix, "_freq_class.png")), p_fc, width=7, height=5, dpi=400)
  }

  # ═════════════════════════════════════════════════════════════════════════
  # 2. ANCESTRY-COLORED PCA
  # ═════════════════════════════════════════════════════════════════════════
  if (!is.null(per_samp) && has_ancestry && "PC1" %in% names(per_samp)) {
    cat("  Ancestry PCA\n")
    pca_m <- merge(per_samp, anc_dt[, .(sample, cluster)],
                   by.x="SAMPLE", by.y="sample", all.x=TRUE)
    pca_m[is.na(cluster), cluster := "unknown"]

    p_anc_pca <- ggplot(pca_m, aes(x=as.numeric(PC1), y=as.numeric(PC2), color=cluster)) +
      geom_point(size=2.5, alpha=0.8) +
      scale_color_manual(values=ancestry_pal, name="Ancestry") +
      geom_hline(yintercept=0, linetype="dashed", color="grey50") +
      geom_vline(xintercept=0, linetype="dashed", color="grey50") +
      labs(x="PC1", y="PC2", title=paste0("PCA — ", vtype, " (ancestry colored)")) +
      theme_pub

    ggsave(file.path(plot_dir, paste0(prefix, "_PCA_ancestry.pdf")), p_anc_pca, width=8, height=6, dpi=400)
    ggsave(file.path(plot_dir, paste0(prefix, "_PCA_ancestry.png")), p_anc_pca, width=8, height=6, dpi=400)
  }

  # ═════════════════════════════════════════════════════════════════════════
  # 3. BURDEN BY FREQUENCY CLASS (stacked)
  # ═════════════════════════════════════════════════════════════════════════
  if (!is.null(per_samp) && "N_SINGLETON" %in% names(per_samp)) {
    cat("  Burden by frequency class\n")
    burden_cols <- c("N_SINGLETON","N_RARE","N_LOW_FREQ","N_COMMON","N_HIGH_FREQ")
    avail <- burden_cols[burden_cols %in% names(per_samp)]
    if (length(avail) > 0) {
      ps_long <- melt(per_samp, id.vars="SAMPLE", measure.vars=avail,
                       variable.name="FREQ_CLASS", value.name="COUNT")
      ps_long[, FREQ_CLASS := gsub("N_", "", tolower(FREQ_CLASS))]
      ps_long[, COUNT := as.numeric(COUNT)]
      total_order <- per_samp[, .(total=as.numeric(TOTAL_VARIANTS)), by=SAMPLE][order(-total), SAMPLE]
      ps_long[, SAMPLE := factor(SAMPLE, levels=total_order)]

      p_burden <- ggplot(ps_long, aes(x=SAMPLE, y=COUNT, fill=FREQ_CLASS)) +
        geom_col(position="stack", width=0.8) +
        scale_fill_manual(values=freq_class_colors, name="Frequency class") +
        labs(x="Sample", y=paste0(vtype, " count"),
             title=paste0("Per-sample ", vtype, " burden by frequency class")) +
        theme_pub +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=2))

      ggsave(file.path(plot_dir, paste0(prefix, "_burden_by_freq.pdf")), p_burden, width=16, height=5, dpi=400)
      ggsave(file.path(plot_dir, paste0(prefix, "_burden_by_freq.png")), p_burden, width=16, height=5, dpi=400)
    }
  }

  # ═════════════════════════════════════════════════════════════════════════
  # 4. QC: BURDEN vs MISSINGNESS
  # ═════════════════════════════════════════════════════════════════════════
  if (!is.null(per_samp) && "MISSINGNESS" %in% names(per_samp)) {
    cat("  QC: burden vs missingness\n")
    per_samp[, MISS := as.numeric(MISSINGNESS)]
    per_samp[, TOTAL := as.numeric(TOTAL_VARIANTS)]

    p_qc <- ggplot(per_samp, aes(x=MISS, y=TOTAL)) +
      geom_point(color=vcol, alpha=0.7, size=2) +
      geom_smooth(method="lm", se=TRUE, color="black", linewidth=0.5) +
      labs(x="Genotype missingness", y=paste0(vtype, " count"),
           title=paste0("QC: ", vtype, " burden vs missingness")) +
      theme_pub

    ggsave(file.path(plot_dir, paste0(prefix, "_qc_burden_vs_miss.pdf")), p_qc, width=6, height=5, dpi=400)
    ggsave(file.path(plot_dir, paste0(prefix, "_qc_burden_vs_miss.png")), p_qc, width=6, height=5, dpi=400)
  }

  # ═════════════════════════════════════════════════════════════════════════
  # 5. 81-ONLY PCA
  # ═════════════════════════════════════════════════════════════════════════
  if (!is.null(per_samp) && length(unrelated_samples) > 0 && "PC1" %in% names(per_samp)) {
    cat("  81-only PCA\n")
    ps81 <- per_samp[SAMPLE %in% unrelated_samples]
    if (nrow(ps81) > 0) {
      p81 <- ggplot(ps81, aes(x=as.numeric(PC1), y=as.numeric(PC2)))
      if (has_ancestry) {
        ps81m <- merge(ps81, anc_dt[, .(sample, cluster)],
                       by.x="SAMPLE", by.y="sample", all.x=TRUE)
        ps81m[is.na(cluster), cluster := "unknown"]
        p81 <- ggplot(ps81m, aes(x=as.numeric(PC1), y=as.numeric(PC2), color=cluster)) +
          scale_color_manual(values=ancestry_pal, name="Ancestry")
      }
      p81 <- p81 + geom_point(size=3, alpha=0.8) +
        geom_hline(yintercept=0, linetype="dashed", color="grey50") +
        geom_vline(xintercept=0, linetype="dashed", color="grey50") +
        labs(title=paste0("PCA — ", vtype, " (81 unrelated)")) + theme_pub

      ggsave(file.path(plot_dir, paste0(prefix, "_81_PCA.pdf")), p81, width=8, height=6, dpi=400)
      ggsave(file.path(plot_dir, paste0(prefix, "_81_PCA.png")), p81, width=8, height=6, dpi=400)
    }
  }

  # ═════════════════════════════════════════════════════════════════════════
  # 6. RARE SHARING NETWORK NODE DEGREE DISTRIBUTION
  # ═════════════════════════════════════════════════════════════════════════
  if (file.exists(rn_file)) {
    cat("  Rare sharing network summary\n")
    rn <- fread(rn_file)
    if (nrow(rn) > 0 && "n_rare_partners" %in% names(rn)) {
      p_rn <- ggplot(rn, aes(x=n_rare_partners, y=total_rare_shared)) +
        geom_point(color=vcol, alpha=0.6, size=2) +
        labs(x="Number of rare-sharing partners", y="Total rare variants shared",
             title=paste0("Rare variant sharing — ", vtype)) +
        theme_pub

      ggsave(file.path(plot_dir, paste0(prefix, "_rare_sharing_scatter.pdf")), p_rn, width=7, height=5, dpi=400)
      ggsave(file.path(plot_dir, paste0(prefix, "_rare_sharing_scatter.png")), p_rn, width=7, height=5, dpi=400)
    }
  }

  # ═════════════════════════════════════════════════════════════════════════
  # 7. PER-SAMPLE WINDOW DENSITY HEATMAP (ComplexHeatmap)
  # ═════════════════════════════════════════════════════════════════════════
  if (file.exists(wd_file)) {
    cat("  Window density heatmap\n")
    wd <- fread(wd_file)
    if (nrow(wd) > 0 && ncol(wd) > 3) {
      wd_samples <- colnames(wd)[4:ncol(wd)]
      mat_wd <- as.matrix(wd[, 4:ncol(wd), with=FALSE])
      mat_wd <- t(mat_wd)  # samples × windows
      rownames(mat_wd) <- wd_samples
      colnames(mat_wd) <- paste0(wd$WINDOW_START / 1e6, "Mb")

      # Sort samples by total
      row_totals <- rowSums(mat_wd)
      mat_wd <- mat_wd[order(-row_totals), , drop=FALSE]

      col_fun <- colorRamp2(c(0, quantile(mat_wd[mat_wd > 0], c(0.5, 0.9, 1))),
                             c("white", "#7ECACA", "#2E8B8B", "#1A5252"))

      pdf(file.path(plot_dir, paste0(prefix, "_window_heatmap.pdf")), width=16, height=10)
      ht <- Heatmap(mat_wd, name=paste0("# ", vtype), col=col_fun,
        cluster_rows=FALSE, cluster_columns=FALSE,
        show_row_names=(nrow(mat_wd) <= 80),
        show_column_names=FALSE,
        row_names_gp=gpar(fontsize=4),
        column_title=paste0(vtype, " density across ", opt$chrom, " (1-Mb windows)"),
        use_raster=TRUE, raster_quality=3)
      draw(ht)
      dev.off()
    }
  }

  # ═════════════════════════════════════════════════════════════════════════
  # 8. MARKER TIER SUMMARY
  # ═════════════════════════════════════════════════════════════════════════
  mk_sum <- file.path(opt$ds_chrom_dir, "markers", "marker_selection_summary.tsv")
  if (file.exists(mk_sum)) {
    cat("  Marker summary\n")
    mks <- fread(mk_sum)
    if (nrow(mks) > 0) {
      mks_long <- melt(mks, id.vars="TYPE", variable.name="TIER", value.name="COUNT")
      mks_long[, TIER := gsub("_", " ", TIER)]
      p_mk <- ggplot(mks_long, aes(x=TYPE, y=COUNT, fill=TIER)) +
        geom_col(position="dodge", width=0.6) +
        scale_fill_manual(values=c("TIER1 STRICT"=col_snp, "TIER3 GROUP INFORMATIVE"=col_indel), name="Tier") +
        geom_text(aes(label=COUNT), position=position_dodge(0.6), vjust=-0.3, size=3) +
        labs(x="Variant type", y="Count", title="Marker selection summary") +
        theme_pub

      ggsave(file.path(plot_dir, paste0("marker_summary.pdf")), p_mk, width=8, height=5, dpi=400)
      ggsave(file.path(plot_dir, paste0("marker_summary.png")), p_mk, width=8, height=5, dpi=400)
    }
    break  # marker summary is the same for all types
  }

  cat("  ", vtype, " extended plots done.\n")
}

cat("\n=== Extended plots complete ===\n")
cat("Output:", plot_dir, "\n")
