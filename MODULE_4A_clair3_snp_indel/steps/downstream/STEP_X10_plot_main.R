#!/usr/bin/env Rscript
# =============================================================================
# 10_plot_main.R — Main Clair3 variant catalog publication plots
#
# Adapted from DELLY 06_plot_delly_results.R for Clair3 SNP + INDEL results.
# Generates per-type (SNP, INDEL, COMBINED) panels:
#   A — Variant density across chromosome (100-kb windows)
#   B — Per-sample variant count barplot (sorted)
#   C — PCA on binary genotype matrix
#   D — Pairwise shared variant heatmap
#   E — Site frequency spectrum (carrier count histogram)
#   F — INDEL size distribution (log-scale, INDEL only)
#
# Usage:
#   Rscript 10_plot_main.R \
#     --ds_chrom_dir <downstream_results/CHROM> \
#     --chrom C_gar_LG01 \
#     --ref_fai PATH \
#     --samples_unrelated PATH
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

# ── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--ds_chrom_dir", type="character"),
  make_option("--chrom",        type="character"),
  make_option("--ref_fai",      type="character"),
  make_option("--samples_unrelated", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option_list))

plot_dir <- file.path(opt$ds_chrom_dir, "figures")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

# ── Theme ───────────────────────────────────────────────────────────────────
theme_pub <- theme_classic(base_size=10) +
  theme(plot.title=element_text(face="bold", size=12),
        axis.text=element_text(color="black"),
        legend.position="right", panel.grid=element_blank())

# ── Colors ──────────────────────────────────────────────────────────────────
col_snp    <- "#2E8B8B"   # teal for SNP
col_indel  <- "#E69F00"   # orange for INDEL
col_main   <- "#2E8B8B"
col_light  <- "#7ECACA"
col_dark   <- "#1A5252"
col_combined <- "#5B4A8A" # purple for combined

type_colors <- c(SNP=col_snp, INDEL=col_indel, COMBINED=col_combined)

# ── Load data ───────────────────────────────────────────────────────────────
cat("=== Loading data ===\n")

fai <- fread(opt$ref_fai, header=FALSE,
             col.names=c("chrom","length","offset","line_bases","line_bytes"))
chr_len <- fai[chrom == opt$chrom, length]

unrelated_samples <- character(0)
if (file.exists(opt$samples_unrelated))
  unrelated_samples <- scan(opt$samples_unrelated, what="character", quiet=TRUE)

# ── Loop per type ───────────────────────────────────────────────────────────
for (vtype in c("SNP", "INDEL", "COMBINED")) {

  cat("\n=== ", vtype, " ===\n")
  vcol <- type_colors[vtype]

  # Paths
  cat_file  <- file.path(opt$ds_chrom_dir, "catalogs", paste0(vtype, "_catalog.tsv"))
  freq_file <- file.path(opt$ds_chrom_dir, "catalogs", paste0(vtype, "_frequency_spectrum.tsv"))
  bin_file  <- file.path(opt$ds_chrom_dir, "matrices", paste0(vtype, "_binary_genotype_matrix.tsv"))  # fixed name
  bin_file2 <- file.path(opt$ds_chrom_dir, "matrices", paste0("binary_genotype_matrix.", vtype, ".tsv"))
  pw_file   <- file.path(opt$ds_chrom_dir, "distances", paste0("pairwise_shared.", vtype, ".tsv"))
  ps_file   <- file.path(opt$ds_chrom_dir, "per_sample", paste0("per_sample_", vtype, "_summary.tsv"))
  wd_file   <- file.path(opt$ds_chrom_dir, "gene_tables", paste0("variant_density_windows.", vtype, ".tsv"))

  bf <- if (file.exists(bin_file)) bin_file else bin_file2
  if (!file.exists(bf)) { cat("  No binary matrix, skipping\n"); next }

  # ── Load ──
  catalog   <- tryCatch(fread(cat_file), error=function(e) NULL)
  freq_spec <- tryCatch(fread(freq_file), error=function(e) NULL)
  binary_gt <- fread(bf)
  pairwise  <- tryCatch(fread(pw_file), error=function(e) NULL)
  per_samp  <- tryCatch(fread(ps_file), error=function(e) NULL)
  win_dens  <- tryCatch(fread(wd_file), error=function(e) NULL)

  n_vars <- nrow(binary_gt)
  samples <- colnames(binary_gt)[-1]
  n_samp <- length(samples)
  cat("  Variants:", n_vars, " Samples:", n_samp, "\n")

  prefix <- tolower(vtype)

  # ═══════════════════════════════════════════════════════════════════════════
  # PANEL A: Variant density along chromosome
  # ═══════════════════════════════════════════════════════════════════════════
  if (!is.null(win_dens) && nrow(win_dens) > 0) {
    cat("  Panel A: density plot\n")
    win_dens[, MID := (WINDOW_START + WINDOW_END) / 2 / 1e6]

    p_a <- ggplot(win_dens, aes(x=MID, y=N_TOTAL)) +
      geom_area(fill=vcol, alpha=0.6) +
      geom_line(color=vcol, linewidth=0.4) +
      labs(x=paste0(opt$chrom, " position (Mb)"),
           y=paste0("# ", vtype, " per 100-kb window"),
           title=paste0("A  ", vtype, " density — ", opt$chrom)) +
      theme_pub

    # If combined, show SNP + INDEL stacked
    if (vtype == "COMBINED" && "N_SNP" %in% names(win_dens)) {
      wd_long <- melt(win_dens, id.vars=c("CHROM","WINDOW_START","WINDOW_END","MID"),
                       measure.vars=c("N_SNP","N_INDEL"),
                       variable.name="TYPE", value.name="COUNT")
      wd_long[, TYPE := gsub("N_", "", TYPE)]
      p_a <- ggplot(wd_long, aes(x=MID, y=COUNT, fill=TYPE)) +
        geom_area(alpha=0.7, position="stack") +
        scale_fill_manual(values=c(SNP=col_snp, INDEL=col_indel), name="Type") +
        labs(x=paste0(opt$chrom, " position (Mb)"),
             y="# variants per 100-kb window",
             title=paste0("A  Variant density — ", opt$chrom)) +
        theme_pub
    }

    ggsave(file.path(plot_dir, paste0(prefix, "_panel_A_density.pdf")), p_a, width=12, height=4, dpi=400)
    ggsave(file.path(plot_dir, paste0(prefix, "_panel_A_density.png")), p_a, width=12, height=4, dpi=400)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PANEL B: Per-sample variant count barplot
  # ═══════════════════════════════════════════════════════════════════════════
  if (!is.null(per_samp) && nrow(per_samp) > 0) {
    cat("  Panel B: per-sample barplot\n")
    per_samp[, source := ifelse(SAMPLE %in% unrelated_samples, "Unrelated (81)", "Related")]
    per_samp[, SAMPLE := factor(SAMPLE, levels=per_samp[order(-TOTAL_VARIANTS), SAMPLE])]

    p_b <- ggplot(per_samp, aes(x=SAMPLE, y=TOTAL_VARIANTS, fill=source)) +
      geom_col(width=0.8) +
      scale_fill_manual(values=c("Unrelated (81)"=vcol, "Related"=col_light), name="Subset") +
      labs(x="Sample", y=paste0(vtype, " count"), title=paste0("B  Per-sample ", vtype, " count")) +
      theme_pub +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=2.5))

    ggsave(file.path(plot_dir, paste0(prefix, "_panel_B_per_sample.pdf")), p_b, width=16, height=5, dpi=400)
    ggsave(file.path(plot_dir, paste0(prefix, "_panel_B_per_sample.png")), p_b, width=16, height=5, dpi=400)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PANEL C: PCA
  # ═══════════════════════════════════════════════════════════════════════════
  cat("  Panel C: PCA\n")
  gt_mat <- as.matrix(binary_gt[, -1, with=FALSE])
  pca_input <- t(gt_mat)
  rownames(pca_input) <- samples
  col_var <- apply(pca_input, 2, var)
  pca_input <- pca_input[, col_var > 0, drop=FALSE]

  if (ncol(pca_input) > 2) {
    pca_res <- prcomp(pca_input, center=TRUE, scale.=FALSE)
    pca_df <- data.table(
      sample=rownames(pca_input),
      PC1=pca_res$x[,1], PC2=pca_res$x[,2]
    )
    pca_df[, source := ifelse(sample %in% unrelated_samples, "Unrelated (81)", "Related")]
    var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

    p_c <- ggplot(pca_df, aes(x=PC1, y=PC2, color=source)) +
      geom_point(size=2.5, alpha=0.8) +
      scale_color_manual(values=c("Unrelated (81)"=vcol, "Related"=col_light), name="Subset") +
      geom_hline(yintercept=0, linetype="dashed", color="grey50") +
      geom_vline(xintercept=0, linetype="dashed", color="grey50") +
      labs(x=paste0("PC1 (", var_pct[1], "%)"),
           y=paste0("PC2 (", var_pct[2], "%)"),
           title=paste0("C  PCA — ", vtype)) +
      theme_pub

    ggsave(file.path(plot_dir, paste0(prefix, "_panel_C_PCA.pdf")), p_c, width=7, height=6, dpi=400)
    ggsave(file.path(plot_dir, paste0(prefix, "_panel_C_PCA.png")), p_c, width=7, height=6, dpi=400)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PANEL D: Pairwise shared heatmap
  # ═══════════════════════════════════════════════════════════════════════════
  if (!is.null(pairwise) && nrow(pairwise) > 0) {
    cat("  Panel D: pairwise heatmap\n")
    pw_wide <- dcast(pairwise, sample_i ~ sample_j, value.var="n_shared", fill=0)
    pw_samples <- as.character(pw_wide$sample_i)
    pw_mat <- as.matrix(pw_wide[, -1, with=FALSE])
    rownames(pw_mat) <- pw_samples
    # Fill symmetric
    for (i in seq_along(pw_samples))
      for (j in seq_along(pw_samples))
        if (pw_mat[i,j] == 0 && j <= ncol(pw_mat) && pw_mat[j,i] > 0) pw_mat[i,j] <- pw_mat[j,i]

    hc <- hclust(as.dist(max(pw_mat) - pw_mat), method="ward.D2")
    col_fun_d <- colorRamp2(
      quantile(pw_mat[upper.tri(pw_mat)], c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE),
      c("#440154", "#31688e", "#21918c", "#5ec962", "#fde725")
    )

    sample_source <- ifelse(pw_samples %in% unrelated_samples, "Unrelated", "Related")
    row_anno <- rowAnnotation(
      Source=sample_source,
      col=list(Source=c(Unrelated=vcol, Related=col_light)),
      show_annotation_name=FALSE
    )

    pdf(file.path(plot_dir, paste0(prefix, "_panel_D_pairwise.pdf")), width=12, height=11)
    ht <- Heatmap(pw_mat, name=paste0("# Shared\n", vtype),
      col=col_fun_d, cluster_rows=hc, cluster_columns=hc,
      show_row_names=(nrow(pw_mat) <= 80),
      show_column_names=(ncol(pw_mat) <= 80),
      row_names_gp=gpar(fontsize=4), column_names_gp=gpar(fontsize=4),
      left_annotation=row_anno,
      column_title=paste0("D  Pairwise shared ", vtype),
      column_title_gp=gpar(fontsize=12, fontface="bold"),
      use_raster=TRUE, raster_quality=3)
    draw(ht)
    dev.off()
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PANEL E: Site frequency spectrum
  # ═══════════════════════════════════════════════════════════════════════════
  if (!is.null(freq_spec) && nrow(freq_spec) > 0) {
    cat("  Panel E: frequency spectrum\n")
    # Exclude monomorphic (0 carriers)
    fs <- freq_spec[n_carriers > 0]
    p_e <- ggplot(fs, aes(x=n_carriers, y=count)) +
      geom_col(fill=vcol, width=0.8) +
      labs(x="Number of carrier samples", y=paste0("# ", vtype, " sites"),
           title=paste0("E  Site frequency spectrum — ", vtype)) +
      theme_pub

    # If many categories, use log y
    if (max(fs$count) > 100 * min(fs$count[fs$count > 0]))
      p_e <- p_e + scale_y_log10()

    ggsave(file.path(plot_dir, paste0(prefix, "_panel_E_SFS.pdf")), p_e, width=10, height=5, dpi=400)
    ggsave(file.path(plot_dir, paste0(prefix, "_panel_E_SFS.png")), p_e, width=10, height=5, dpi=400)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # PANEL F: INDEL size distribution (INDEL/COMBINED only)
  # ═══════════════════════════════════════════════════════════════════════════
  if (vtype %in% c("INDEL", "COMBINED") && !is.null(catalog)) {
    cat("  Panel F: INDEL size distribution\n")
    indels <- catalog[VAR_TYPE != "SNP" & !is.na(ABS_INDEL_LEN)]
    indels[, abs_len := as.numeric(ABS_INDEL_LEN)]
    indels <- indels[abs_len > 0]

    if (nrow(indels) > 0) {
      p_f <- ggplot(indels, aes(x=abs_len)) +
        geom_histogram(bins=80, fill=col_indel, color="white", linewidth=0.2) +
        scale_x_log10(labels=scales::comma) +
        annotation_logticks(sides="b") +
        labs(x="INDEL size (bp, log scale)", y="Count",
             title=paste0("F  INDEL size distribution")) +
        theme_pub

      ggsave(file.path(plot_dir, paste0(prefix, "_panel_F_indel_size.pdf")), p_f, width=7, height=5, dpi=400)
      ggsave(file.path(plot_dir, paste0(prefix, "_panel_F_indel_size.png")), p_f, width=7, height=5, dpi=400)
    }
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # COMPOSITE: B + C + E (+ F if INDEL)
  # ═══════════════════════════════════════════════════════════════════════════
  cat("  Composite figure\n")
  panel_list <- list()
  if (exists("p_b")) panel_list[["B"]] <- p_b + theme(legend.position="none")
  if (exists("p_c")) panel_list[["C"]] <- p_c + theme(legend.position="none")
  if (exists("p_e")) panel_list[["E"]] <- p_e

  if (length(panel_list) >= 3) {
    comp <- plot_grid(plotlist=panel_list, ncol=2, nrow=2)
    ggsave(file.path(plot_dir, paste0(prefix, "_composite_main.pdf")), comp, width=18, height=10, dpi=400)
    ggsave(file.path(plot_dir, paste0(prefix, "_composite_main.png")), comp, width=18, height=10, dpi=400)
  }

  cat("  ", vtype, " plots done.\n")

  # Clean up per-iteration objects
  rm(list=intersect(ls(), c("p_a","p_b","p_c","p_e","p_f","catalog","freq_spec",
                             "binary_gt","pairwise","per_samp","win_dens",
                             "pca_res","pca_df","pw_wide","pw_mat")))
}

cat("\n=== All main plots complete ===\n")
cat("Output:", plot_dir, "\n")
