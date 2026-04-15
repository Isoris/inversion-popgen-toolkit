#!/usr/bin/env Rscript
# =============================================================================
# 12_plot_phase_signatures.R — Phase block and inversion marker plots
#
# Visualizes:
#   1. Phase block statistics per sample (TIER1 vs TIER2)
#   2. Phase block size/span distributions
#   3. Phase signature distance MDS
#   4. Inversion marker candidate density along chromosome
#   5. Phase block coverage along chromosome (fraction of samples phased per window)
#
# Usage:
#   Rscript 12_plot_phase_signatures.R \
#     --ds_chrom_dir DIR --chrom CHROM --ref_fai PATH \
#     [--samples_unrelated PATH] [--qopt PATH --qopt_samples PATH]
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(viridis)
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

col_t1 <- "#2E8B8B"
col_t2 <- "#E69F00"

ancestry_pal <- c(Q1="#E41A1C", Q2="#377EB8", Q3="#4DAF4A", Q4="#984EA3",
                   Q5="#FF7F00", Q6="#FFFF33", Q7="#A65628", Q8="#F781BF")

# Load ancestry
has_ancestry <- FALSE; anc_dt <- NULL
if (nchar(opt$qopt) > 0 && file.exists(opt$qopt) &&
    nchar(opt$qopt_samples) > 0 && file.exists(opt$qopt_samples)) {
  qsamp <- scan(opt$qopt_samples, what="character", quiet=TRUE)
  qmat <- as.matrix(fread(opt$qopt, header=FALSE))
  if (nrow(qmat) == length(qsamp)) {
    anc_dt <- data.table(sample=qsamp, cluster=paste0("Q", apply(qmat, 1, which.max)))
    has_ancestry <- TRUE
  }
}

unrelated_samples <- character(0)
if (file.exists(opt$samples_unrelated))
  unrelated_samples <- scan(opt$samples_unrelated, what="character", quiet=TRUE)

# ── Load phase data ──
phase_dir <- file.path(opt$ds_chrom_dir, "phase_prep")
blocks_file <- file.path(phase_dir, "phase_block_catalog.tsv")
ps_sum_file <- file.path(phase_dir, "phase_block_summary_per_sample.tsv")
inv_file    <- file.path(phase_dir, "inversion_marker_candidates.tsv")
dist_file   <- file.path(opt$ds_chrom_dir, "distances", "signature_distance_matrix.tsv")

cat("=== Phase signature plots ===\n")

# ═════════════════════════════════════════════════════════════════════════════
# 1. PER-SAMPLE PHASE STATS
# ═════════════════════════════════════════════════════════════════════════════
if (file.exists(ps_sum_file)) {
  cat("  Per-sample phase stats\n")
  ps <- fread(ps_sum_file)

  # Stacked bar: TIER1 vs TIER2 vs UNPHASED
  ps_long <- melt(ps, id.vars="SAMPLE",
                   measure.vars=c("N_TIER1_VARIANTS","N_TIER2_VARIANTS","N_UNPHASED"),
                   variable.name="TIER", value.name="COUNT")
  ps_long[, TIER := gsub("N_TIER1_VARIANTS", "WhatsHap", TIER)]
  ps_long[, TIER := gsub("N_TIER2_VARIANTS", "Read-pair", TIER)]
  ps_long[, TIER := gsub("N_UNPHASED", "Unphased", TIER)]
  ps_long[, TIER := factor(TIER, levels=c("Unphased","Read-pair","WhatsHap"))]
  total_order <- ps[, .(total=N_TIER1_VARIANTS+N_TIER2_VARIANTS+N_UNPHASED), by=SAMPLE][order(-total), SAMPLE]
  ps_long[, SAMPLE := factor(SAMPLE, levels=total_order)]

  p1 <- ggplot(ps_long, aes(x=SAMPLE, y=COUNT, fill=TIER)) +
    geom_col(position="stack", width=0.8) +
    scale_fill_manual(values=c(WhatsHap=col_t1, `Read-pair`=col_t2, Unphased="#BDBDBD"), name="Phase source") +
    labs(x="Sample", y="Variant count", title=paste0("Phase composition per sample — ", opt$chrom)) +
    theme_pub + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=2))

  ggsave(file.path(plot_dir, "phase_per_sample_stacked.pdf"), p1, width=16, height=5, dpi=400)
  ggsave(file.path(plot_dir, "phase_per_sample_stacked.png"), p1, width=16, height=5, dpi=400)

  # Block count scatter: T1 vs T2
  p1b <- ggplot(ps, aes(x=N_TIER1_BLOCKS, y=N_TIER2_BLOCKS)) +
    geom_point(color=col_t1, alpha=0.7, size=2) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
    labs(x="WhatsHap blocks", y="Read-pair blocks",
         title=paste0("Phase blocks per sample — ", opt$chrom)) +
    theme_pub

  ggsave(file.path(plot_dir, "phase_blocks_T1_vs_T2.pdf"), p1b, width=6, height=5, dpi=400)
  ggsave(file.path(plot_dir, "phase_blocks_T1_vs_T2.png"), p1b, width=6, height=5, dpi=400)
}

# ═════════════════════════════════════════════════════════════════════════════
# 2. PHASE BLOCK SIZE/SPAN DISTRIBUTIONS
# ═════════════════════════════════════════════════════════════════════════════
if (file.exists(blocks_file)) {
  cat("  Phase block distributions\n")
  blocks <- fread(blocks_file)

  if (nrow(blocks) > 0) {
    # Block span distribution by tier
    p2a <- ggplot(blocks, aes(x=BLOCK_SPAN + 1, fill=PHASE_TIER)) +
      geom_histogram(bins=60, alpha=0.7, position="identity") +
      scale_x_log10(labels=scales::comma) +
      scale_fill_manual(values=c(TIER_1_WHATSHAP=col_t1, TIER_2_READPAIR=col_t2), name="Phase tier") +
      annotation_logticks(sides="b") +
      labs(x="Block span (bp, log)", y="Count",
           title=paste0("Phase block span distribution — ", opt$chrom)) +
      theme_pub

    ggsave(file.path(plot_dir, "phase_block_span_dist.pdf"), p2a, width=8, height=5, dpi=400)
    ggsave(file.path(plot_dir, "phase_block_span_dist.png"), p2a, width=8, height=5, dpi=400)

    # Variants per block
    p2b <- ggplot(blocks, aes(x=N_VARIANTS, fill=PHASE_TIER)) +
      geom_histogram(binwidth=1, alpha=0.7, position="identity") +
      scale_fill_manual(values=c(TIER_1_WHATSHAP=col_t1, TIER_2_READPAIR=col_t2), name="Phase tier") +
      coord_cartesian(xlim=c(1, 20)) +
      labs(x="Variants per block", y="Count",
           title=paste0("Phase block size distribution — ", opt$chrom)) +
      theme_pub

    ggsave(file.path(plot_dir, "phase_block_size_dist.pdf"), p2b, width=8, height=5, dpi=400)
    ggsave(file.path(plot_dir, "phase_block_size_dist.png"), p2b, width=8, height=5, dpi=400)
  }
}

# ═════════════════════════════════════════════════════════════════════════════
# 3. SIGNATURE DISTANCE MDS
# ═════════════════════════════════════════════════════════════════════════════
if (file.exists(dist_file)) {
  cat("  Signature distance MDS\n")
  dist_dt <- fread(dist_file)
  dist_samples <- dist_dt[[1]]
  dist_mat <- as.matrix(dist_dt[, -1, with=FALSE])
  rownames(dist_mat) <- dist_samples
  # Ensure symmetric
  dist_mat[is.na(dist_mat)] <- 1
  dist_mat <- (dist_mat + t(dist_mat)) / 2
  diag(dist_mat) <- 0

  if (nrow(dist_mat) > 3) {
    mds_res <- cmdscale(as.dist(dist_mat), k=2)
    mds_df <- data.table(sample=dist_samples, MDS1=mds_res[,1], MDS2=mds_res[,2])
    mds_df[, source := ifelse(sample %in% unrelated_samples, "Unrelated", "Related")]

    p3 <- ggplot(mds_df, aes(x=MDS1, y=MDS2))

    if (has_ancestry) {
      mds_df <- merge(mds_df, anc_dt, by="sample", all.x=TRUE)
      mds_df[is.na(cluster), cluster := "unknown"]
      p3 <- ggplot(mds_df, aes(x=MDS1, y=MDS2, color=cluster)) +
        scale_color_manual(values=ancestry_pal, name="Ancestry")
    } else {
      p3 <- ggplot(mds_df, aes(x=MDS1, y=MDS2, color=source)) +
        scale_color_manual(values=c(Unrelated=col_t1, Related="#7ECACA"), name="Subset")
    }

    p3 <- p3 + geom_point(size=2.5, alpha=0.8) +
      labs(x="MDS dimension 1", y="MDS dimension 2",
           title=paste0("Phase signature MDS — ", opt$chrom)) +
      theme_pub

    ggsave(file.path(plot_dir, "phase_signature_MDS.pdf"), p3, width=8, height=6, dpi=400)
    ggsave(file.path(plot_dir, "phase_signature_MDS.png"), p3, width=8, height=6, dpi=400)
  }
}

# ═════════════════════════════════════════════════════════════════════════════
# 4. INVERSION MARKER CANDIDATE DENSITY
# ═════════════════════════════════════════════════════════════════════════════
if (file.exists(inv_file)) {
  cat("  Inversion marker density\n")
  inv <- fread(inv_file)

  if (nrow(inv) > 0) {
    # Per 100-kb window: number of inversion marker candidates
    inv[, POS := as.numeric(POS)]
    inv[, window := floor(POS / 1e5) * 1e5]

    # Count unique positions per window (across all samples)
    win_counts <- inv[, .(n_candidates=uniqueN(POS), n_samples=uniqueN(SAMPLE)), by=window]
    win_counts[, mid_Mb := (window + 50000) / 1e6]

    p4 <- ggplot(win_counts, aes(x=mid_Mb, y=n_candidates)) +
      geom_col(fill=col_t1, width=0.08) +
      labs(x=paste0(opt$chrom, " position (Mb)"),
           y="Phased het SNP positions (per 100-kb)",
           title=paste0("Inversion marker candidate density — ", opt$chrom)) +
      theme_pub

    ggsave(file.path(plot_dir, "inversion_marker_density.pdf"), p4, width=12, height=4, dpi=400)
    ggsave(file.path(plot_dir, "inversion_marker_density.png"), p4, width=12, height=4, dpi=400)

    # Phase GT pattern: 0|1 vs 1|0 ratio per window
    inv_gt <- inv[PHASE_GT %in% c("0|1","1|0")]
    if (nrow(inv_gt) > 0) {
      gt_window <- inv_gt[, .N, by=.(window, PHASE_GT)]
      gt_window[, mid_Mb := (window + 50000) / 1e6]
      gt_wide <- dcast(gt_window, window + mid_Mb ~ PHASE_GT, value.var="N", fill=0)
      if ("0|1" %in% names(gt_wide) && "1|0" %in% names(gt_wide)) {
        gt_wide[, total := `0|1` + `1|0`]
        gt_wide[, frac_01 := `0|1` / total]
        gt_wide <- gt_wide[total >= 5]  # only plot windows with enough data

        if (nrow(gt_wide) > 0) {
          p4b <- ggplot(gt_wide, aes(x=mid_Mb, y=frac_01)) +
            geom_point(size=1, alpha=0.5, color=col_t1) +
            geom_hline(yintercept=0.5, linetype="dashed", color="red") +
            labs(x=paste0(opt$chrom, " position (Mb)"),
                 y="Fraction 0|1 (vs 1|0)",
                 title=paste0("Phase GT balance along ", opt$chrom,
                              " (deviation from 0.5 may indicate inversions)")) +
            ylim(0, 1) + theme_pub

          ggsave(file.path(plot_dir, "inversion_phase_balance.pdf"), p4b, width=12, height=4, dpi=400)
          ggsave(file.path(plot_dir, "inversion_phase_balance.png"), p4b, width=12, height=4, dpi=400)
        }
      }
    }
  }
}

# ═════════════════════════════════════════════════════════════════════════════
# 5. PHASE BLOCK COVERAGE ALONG CHROMOSOME
# ═════════════════════════════════════════════════════════════════════════════
if (file.exists(blocks_file)) {
  cat("  Phase block coverage\n")
  blocks <- fread(blocks_file)
  if (nrow(blocks) > 0) {
    # For each 100-kb window, count how many samples have at least one block
    blocks[, window := floor(BLOCK_START / 1e5) * 1e5]
    cov_dt <- blocks[, .(n_samples=uniqueN(SAMPLE)), by=.(window, PHASE_TIER)]
    cov_dt[, mid_Mb := (window + 50000) / 1e6]

    p5 <- ggplot(cov_dt, aes(x=mid_Mb, y=n_samples, fill=PHASE_TIER)) +
      geom_col(position="dodge", width=0.08) +
      scale_fill_manual(values=c(TIER_1_WHATSHAP=col_t1, TIER_2_READPAIR=col_t2), name="Phase tier") +
      labs(x=paste0(opt$chrom, " position (Mb)"),
           y="# samples with phase block",
           title=paste0("Phase block coverage — ", opt$chrom)) +
      theme_pub

    ggsave(file.path(plot_dir, "phase_block_coverage.pdf"), p5, width=12, height=4, dpi=400)
    ggsave(file.path(plot_dir, "phase_block_coverage.png"), p5, width=12, height=4, dpi=400)
  }
}

# ═════════════════════════════════════════════════════════════════════════════
# COMPOSITE: Phase summary figure
# ═════════════════════════════════════════════════════════════════════════════
cat("  Composite phase figure\n")
panels <- list()
if (exists("p1"))  panels[["A"]] <- p1 + ggtitle("A") + theme(legend.position="bottom")
if (exists("p2a")) panels[["B"]] <- p2a + ggtitle("B")
if (exists("p3"))  panels[["C"]] <- p3 + ggtitle("C")
if (exists("p4"))  panels[["D"]] <- p4 + ggtitle("D")

if (length(panels) >= 3) {
  comp <- plot_grid(plotlist=panels, ncol=2, nrow=2)
  ggsave(file.path(plot_dir, "composite_phase_summary.pdf"), comp, width=18, height=10, dpi=400)
  ggsave(file.path(plot_dir, "composite_phase_summary.png"), comp, width=18, height=10, dpi=400)
}

cat("\n=== Phase signature plots complete ===\n")
cat("Output:", plot_dir, "\n")
