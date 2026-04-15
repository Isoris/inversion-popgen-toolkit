#!/usr/bin/env Rscript
# =============================================================================
# 06_plot_manta_results.R — v3: publication panels with tier composition
# =============================================================================
#
# Panels:
#   A. Tier composition waterfall (stacked bars: PASS → excluded → tier1 → tier2 → tier3)
#   B. Genome-wide SV density heatmap (DEL, 1-Mb windows)
#   C. Per-sample SV burden — PASS level (stacked by type)
#   D. Per-sample SV burden — tier2 level (stacked by type)
#   E. Size distribution per SV type (violin + boxplot, log scale)
#   F. Private vs shared spectrum per type
#   G. Per-chromosome SV counts (grouped bars)
#   H. Evidence distribution (sum_total_alt) per SV type
#   I. Precision composition per type (PRECISE vs IMPRECISE)
#   J. Pairwise shared DEL heatmap (226 cohort)
#   K. Tier composition per type — 81 cohort (separate panel)
#
# Usage:
#   Rscript 06_plot_manta_results.R \
#     --summary_dir DIR_SUMMARY \
#     --final_dir DIR_FINAL \
#     --filter_dir 10_filtered_catalogs \
#     --plot_dir DIR_PLOTS \
#     --ref_fai REF.fai
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(scales)
  library(optparse)
})

opts <- list(
  make_option("--summary_dir", type="character"),
  make_option("--final_dir", type="character"),
  make_option("--filter_dir", type="character"),
  make_option("--plot_dir", type="character"),
  make_option("--ref_fai", type="character")
)
args <- parse_args(OptionParser(option_list=opts))

summary_dir <- args$summary_dir
final_dir   <- args$final_dir
filter_dir  <- args$filter_dir
plot_dir    <- args$plot_dir
ref_fai     <- args$ref_fai

dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)
cat("Plot output:", plot_dir, "\n")

# ── Reference ──────────────────────────────────────────────────────────────
fai <- fread(ref_fai, header=FALSE, col.names=c("chrom","length","offset","bases","bytes"))
chrom_order <- fai$chrom

# ── Colors ─────────────────────────────────────────────────────────────────
sv_colors <- c(
  DEL="#E41A1C", DUP="#377EB8", INV="#4DAF4A",
  INS_small="#FF7F00", INS_large="#FFCC00", BND="#984EA3"
)
tier_colors <- c(
  excluded="#CCCCCC", tier1_only="#FDB462", tier2_only="#80B1D3", tier3="#4DAF4A"
)
precision_colors <- c(PRECISE="#2166AC", IMPRECISE="#F4A582")
sv_types <- c("DEL","DUP","INV","INS_small","INS_large","BND")

save_plot <- function(name, p, w=12, h=6) {
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), p, width=w, height=h)
  ggsave(file.path(plot_dir, paste0(name, ".png")), p, width=w, height=h, dpi=200)
  cat("  Saved:", name, "\n")
}

# =============================================================================
# PANEL A: Tier composition waterfall — both cohorts
# =============================================================================
cat("Panel A: Tier composition...\n")

for (cohort_n in c("226","81")) {
  tc_f <- file.path(summary_dir, paste0("tier_composition_", cohort_n, ".tsv"))
  if (!file.exists(tc_f)) next

  tc <- fread(tc_f)
  tc_long <- melt(tc[, .(svtype, excluded_by_tier1, tier1_only, tier2_only, tier3)],
                  id.vars="svtype", variable.name="segment", value.name="count")
  tc_long[, segment := factor(segment,
    levels=c("excluded_by_tier1","tier1_only","tier2_only","tier3"),
    labels=c("excluded","tier1_only","tier2_only","tier3"))]
  tc_long[, svtype := factor(svtype, levels=sv_types)]

  pA <- ggplot(tc_long, aes(x=svtype, y=count, fill=segment)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=tier_colors, name="Filter tier") +
    labs(x="SV type", y="Variant count",
         title=paste0("Tier composition (", cohort_n, " cohort, PASS → tier3)")) +
    theme_cowplot(font_size=11) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

  save_plot(paste0("A_tier_composition_", cohort_n), pA, w=10, h=6)
}

# =============================================================================
# PANEL B: Genome-wide DEL density heatmap
# =============================================================================
cat("Panel B: Genome heatmap...\n")
wf <- file.path(summary_dir, "DEL_windows_1Mb_226.tsv")
if (file.exists(wf)) {
  win <- fread(wf, header=FALSE, col.names=c("chrom","window_start","window_end","n_SVs"))
  win[, chrom := factor(chrom, levels=chrom_order)]
  win[, pos_Mb := window_start / 1e6]

  pB <- ggplot(win, aes(x=pos_Mb, y=chrom, fill=n_SVs)) +
    geom_tile(width=1, height=0.8) +
    scale_fill_viridis(name="DEL\ncount", option="inferno", trans="sqrt") +
    labs(x="Position (Mb)", y="", title="DEL density (1-Mb windows, 226 PASS)") +
    theme_minimal(base_size=10) +
    theme(panel.grid=element_blank())

  save_plot("B_genome_DEL_heatmap", pB, w=14, h=6)
}

# =============================================================================
# PANEL C+D: Per-sample burden — PASS and tier2
# =============================================================================
cat("Panel C/D: Per-sample burden...\n")

for (tier_label in c("PASS","tier2")) {
  for (cohort_n in c("226","81")) {
    psf <- file.path(summary_dir, paste0("per_sample_", tier_label, "_", cohort_n, ".tsv"))
    if (!file.exists(psf)) next

    ps <- fread(psf)
    # Melt to long
    type_cols <- intersect(sv_types, names(ps))
    if (length(type_cols) == 0) next

    ps_long <- melt(ps, id.vars=c("sample","total"),
                    measure.vars=type_cols, variable.name="svtype", value.name="count")
    ps_long[, svtype := factor(svtype, levels=sv_types)]

    sample_order <- ps[order(-total), sample]
    ps_long[, sample := factor(sample, levels=sample_order)]

    p <- ggplot(ps_long, aes(x=sample, y=count, fill=svtype)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=sv_colors, name="SV type") +
      labs(x=paste0("Sample (n=", length(sample_order), ", sorted by total)"),
           y="SV count",
           title=paste0("Per-sample burden — ", tier_label, " (", cohort_n, ")")) +
      theme_cowplot(font_size=10) +
      theme(axis.text.x=element_blank(), panel.grid.major.x=element_blank())

    panel_letter <- ifelse(tier_label == "PASS", "C", "D")
    save_plot(paste0(panel_letter, "_burden_", tier_label, "_", cohort_n), p, w=14, h=5)
  }
}

# =============================================================================
# PANEL E: Size distribution per SV type
# =============================================================================
cat("Panel E: Size distributions...\n")

for (cohort_n in c("226","81")) {
  size_list <- list()
  for (svt in sv_types[sv_types != "BND"]) {
    sf <- file.path(summary_dir, paste0(svt, "_svlen_", cohort_n, ".tsv"))
    if (!file.exists(sf)) next
    dt <- fread(sf, header=FALSE, col.names="svlen_bp")
    if (nrow(dt) > 0) {
      dt[, svtype := svt]
      size_list[[svt]] <- dt
    }
  }
  if (length(size_list) == 0) next

  sizes <- rbindlist(size_list)
  sizes[, svtype := factor(svtype, levels=sv_types)]
  sizes <- sizes[svlen_bp > 0]

  pE <- ggplot(sizes, aes(x=svtype, y=svlen_bp, fill=svtype)) +
    geom_violin(alpha=0.6, scale="width") +
    geom_boxplot(width=0.15, outlier.size=0.3, alpha=0.8) +
    scale_y_log10(labels=comma, breaks=c(50,500,5e3,5e4,5e5,5e6)) +
    scale_fill_manual(values=sv_colors, guide="none") +
    labs(x="SV type", y="Size (bp, log scale)",
         title=paste0("SV size distribution (", cohort_n, " PASS)")) +
    theme_cowplot(font_size=11) +
    annotation_logticks(sides="l")

  save_plot(paste0("E_size_dist_", cohort_n), pE, w=9, h=5)
}

# =============================================================================
# PANEL F: Sharing spectrum
# =============================================================================
cat("Panel F: Sharing spectrum...\n")

for (cohort_n in c("226","81")) {
  share_list <- list()
  for (svt in sv_types[sv_types != "BND"]) {
    sf <- file.path(summary_dir, paste0(svt, "_sharing_", cohort_n, ".tsv"))
    if (!file.exists(sf)) next
    dt <- fread(sf)
    if (nrow(dt) > 0) {
      dt[, svtype := svt]
      share_list[[svt]] <- dt
    }
  }
  if (length(share_list) == 0) next

  shares <- rbindlist(share_list, fill=TRUE)
  shares[, svtype := factor(svtype, levels=sv_types)]
  shares[, category := factor(category,
    levels=c("private","shared_2_5","shared_6_20","shared_21plus","fixed"))]

  pF <- ggplot(shares, aes(x=category, y=count, fill=svtype)) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=sv_colors, name="SV type") +
    labs(x="Sharing category", y="Count",
         title=paste0("Allele sharing spectrum (", cohort_n, ")")) +
    theme_cowplot(font_size=11) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

  save_plot(paste0("F_sharing_", cohort_n), pF, w=10, h=5)
}

# =============================================================================
# PANEL G: Per-chromosome SV counts
# =============================================================================
cat("Panel G: Per-chromosome...\n")

for (cohort_n in c("226","81")) {
  chr_list <- list()
  for (svt in sv_types[sv_types != "BND"]) {
    cf <- file.path(summary_dir, paste0(svt, "_per_chr_", cohort_n, ".tsv"))
    if (!file.exists(cf)) next
    dt <- fread(cf, header=FALSE, col.names=c("chrom","n_SVs"))
    if (nrow(dt) > 0) {
      dt[, svtype := svt]
      chr_list[[svt]] <- dt
    }
  }
  if (length(chr_list) == 0) next

  chrs <- rbindlist(chr_list, fill=TRUE)
  chrs[, svtype := factor(svtype, levels=sv_types)]
  chrs[, chrom := factor(chrom, levels=chrom_order)]

  pG <- ggplot(chrs, aes(x=chrom, y=n_SVs, fill=svtype)) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=sv_colors, name="SV type") +
    labs(x="Chromosome", y="SV count",
         title=paste0("Per-chromosome SV counts (", cohort_n, " PASS)")) +
    theme_cowplot(font_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=7))

  save_plot(paste0("G_per_chr_", cohort_n), pG, w=14, h=5)
}

# =============================================================================
# PANEL H: Evidence distribution (sum_total_alt) per type
# =============================================================================
cat("Panel H: Evidence distribution...\n")

for (cohort_n in c("226","81")) {
  ev_list <- list()
  for (svt in sv_types) {
    summary_f <- file.path(filter_dir, paste0(cohort_n, ".", svt, ".summary.tsv"))
    if (!file.exists(summary_f)) next
    dt <- fread(summary_f)
    if (nrow(dt) > 0 && "sum_total_alt" %in% names(dt)) {
      ev_list[[svt]] <- data.table(svtype=svt, total_alt=dt$sum_total_alt)
    }
  }
  if (length(ev_list) == 0) next

  ev <- rbindlist(ev_list)
  ev[, svtype := factor(svtype, levels=sv_types)]
  ev[, total_alt_capped := pmin(total_alt, 50)]  # cap at 50 for visibility

  pH <- ggplot(ev, aes(x=svtype, y=total_alt_capped, fill=svtype)) +
    geom_violin(alpha=0.6, scale="width") +
    geom_boxplot(width=0.15, outlier.size=0.3, alpha=0.8) +
    geom_hline(yintercept=c(4, 6, 10), linetype="dashed", color="grey40", linewidth=0.4) +
    annotate("text", x=0.5, y=c(4,6,10), label=c("tier1","tier2","tier3"),
             hjust=0, vjust=-0.3, size=3, color="grey30") +
    scale_fill_manual(values=sv_colors, guide="none") +
    labs(x="SV type", y="Total alt evidence (PR+SR, capped at 50)",
         title=paste0("Evidence support per SV (", cohort_n, " PASS)")) +
    theme_cowplot(font_size=11)

  save_plot(paste0("H_evidence_", cohort_n), pH, w=9, h=5)
}

# =============================================================================
# PANEL I: Precision composition per type (PRECISE vs IMPRECISE)
# =============================================================================
cat("Panel I: Precision composition...\n")

for (cohort_n in c("226","81")) {
  prec_list <- list()
  for (svt in sv_types) {
    summary_f <- file.path(filter_dir, paste0(cohort_n, ".", svt, ".summary.tsv"))
    if (!file.exists(summary_f)) next
    dt <- fread(summary_f)
    if (nrow(dt) > 0 && "is_imprecise" %in% names(dt)) {
      counts <- dt[, .N, by=is_imprecise]
      counts[, svtype := svt]
      prec_list[[svt]] <- counts
    }
  }
  if (length(prec_list) == 0) next

  prec <- rbindlist(prec_list)
  prec[, svtype := factor(svtype, levels=sv_types)]

  pI <- ggplot(prec, aes(x=svtype, y=N, fill=is_imprecise)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=precision_colors, name="Breakpoint") +
    labs(x="SV type", y="Variant count",
         title=paste0("Precision composition (", cohort_n, " PASS)")) +
    theme_cowplot(font_size=11) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

  save_plot(paste0("I_precision_", cohort_n), pI, w=9, h=5)
}

# =============================================================================
# PANEL J: Pairwise shared DEL heatmap (226)
# =============================================================================
cat("Panel J: Pairwise DEL heatmap...\n")

pw_f <- file.path(summary_dir, "pairwise_shared_DEL_226.tsv")
if (file.exists(pw_f)) {
  pw <- fread(pw_f)
  if (nrow(pw) > 0) {
    all_samples <- sort(unique(c(pw$sample_i, pw$sample_j)))
    ns <- length(all_samples)
    mat <- matrix(0, nrow=ns, ncol=ns, dimnames=list(all_samples, all_samples))
    for (i in seq_len(nrow(pw))) {
      si <- pw$sample_i[i]; sj <- pw$sample_j[i]; n <- pw$n_shared[i]
      mat[si, sj] <- n; mat[sj, si] <- n
    }

    d <- as.dist(max(mat) - mat)
    hc <- hclust(d, method="ward.D2")
    ord <- hc$order
    mat_ord <- mat[ord, ord]

    melted <- data.table(
      row=rep(seq_len(ns), ns),
      col=rep(seq_len(ns), each=ns),
      value=as.vector(mat_ord)
    )

    pJ <- ggplot(melted, aes(x=col, y=row, fill=value)) +
      geom_tile() +
      scale_fill_viridis(name="Shared\nDELs", option="magma") +
      labs(x="", y="", title="Pairwise shared DELs (226 PASS)") +
      theme_minimal(base_size=8) +
      theme(axis.text=element_blank(), panel.grid=element_blank())

    save_plot("J_pairwise_DEL_226", pJ, w=8, h=7)
  }
}

# =============================================================================
cat("\n=== All panels complete ===\n")
cat("Output:", plot_dir, "\n")
flist <- list.files(plot_dir, pattern="\\.(pdf|png)$")
cat(paste(flist, collapse="\n"), "\n")
