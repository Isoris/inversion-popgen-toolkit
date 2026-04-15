#!/usr/bin/env Rscript
# =============================================================================
# run_stats.R — Statistics for het/ROH module (sample + chromosome level)
# =============================================================================
# Usage:
#   Rscript run_stats.R <master_summary.tsv> <stats_dir> [ancestry_labels.tsv] \
#     [per_chr_roh_summary.tsv]
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: run_stats.R <master.tsv> <stats_dir> [ancestry] [per_chr]")

master_file <- args[1]
stats_dir   <- args[2]
anc_file    <- if (length(args) >= 3 && nchar(args[3]) > 0 && file.exists(args[3])) args[3] else NULL
chr_file    <- if (length(args) >= 4 && nchar(args[4]) > 0 && file.exists(args[4])) args[4] else NULL

dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

d <- fread(master_file)
for (col in c("het_genomewide", "froh", "roh_total_bp", "longest_roh")) {
  if (col %in% names(d)) d[, (col) := as.numeric(get(col))]
}

# ── Merge ancestry ───────────────────────────────────────────────────────
if (!is.null(anc_file)) {
  anc <- fread(anc_file)
  if (ncol(anc) >= 2) {
    setnames(anc, 1:2, c("sample", "ancestry"))
    d <- merge(d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
  }
}
if (!"ancestry" %in% names(d)) d[, ancestry := "all"]

# ═══════════════════════════════════════════════════════════════════════════
# SAMPLE-LEVEL STATS
# ═══════════════════════════════════════════════════════════════════════════

# 1. Spearman correlations
cor_results <- list()
for (pair in list(c("het_genomewide","froh"), c("het_genomewide","roh_total_bp"),
                  c("het_genomewide","longest_roh"))) {
  v1 <- pair[1]; v2 <- pair[2]
  if (!all(c(v1, v2) %in% names(d))) next
  x <- d[[v1]]; y <- d[[v2]]
  ok <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  if (sum(ok) < 5) next
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  cor_results[[length(cor_results) + 1]] <- data.table(
    var1 = v1, var2 = v2, n = sum(ok), rho = ct$estimate,
    statistic = ct$statistic, p_value = ct$p.value,
    direction = ifelse(ct$estimate > 0, "positive", "negative"))
}
if (length(cor_results) > 0) {
  fwrite(rbindlist(cor_results), file.path(stats_dir, "spearman_correlations.tsv"), sep = "\t")
  cat("Wrote spearman_correlations.tsv\n")
}

# 2. Group comparisons
groups <- unique(d$ancestry[!is.na(d$ancestry)])
if (length(groups) > 1) {
  group_results <- list()
  for (v in c("het_genomewide", "froh", "roh_total_bp", "longest_roh", "n_tracts")) {
    if (!v %in% names(d)) next
    vals <- as.numeric(d[[v]])
    ok <- !is.na(vals) & is.finite(vals) & !is.na(d$ancestry)
    if (sum(ok) < 5) next
    dsub <- d[ok]; ng <- length(unique(dsub$ancestry))
    if (ng == 2) {
      wt <- wilcox.test(as.formula(paste(v, "~ ancestry")), data = dsub)
      group_results[[length(group_results) + 1]] <- data.table(
        variable = v, grouping = "ancestry", test = "Wilcoxon_rank_sum",
        n = nrow(dsub), n_groups = ng, statistic = wt$statistic,
        p_value = wt$p.value, adj_p_value = NA_real_, direction = "")
    } else if (ng > 2) {
      kt <- kruskal.test(as.formula(paste(v, "~ ancestry")), data = dsub)
      group_results[[length(group_results) + 1]] <- data.table(
        variable = v, grouping = "ancestry", test = "Kruskal_Wallis",
        n = nrow(dsub), n_groups = ng, statistic = kt$statistic,
        p_value = kt$p.value, adj_p_value = NA_real_, direction = "")
      pw <- pairwise.wilcox.test(dsub[[v]], dsub$ancestry, p.adjust.method = "BH")
      pm <- pw$p.value
      for (i in seq_len(nrow(pm))) for (j in seq_len(ncol(pm)))
        if (!is.na(pm[i, j]))
          group_results[[length(group_results) + 1]] <- data.table(
            variable = v, grouping = paste0(rownames(pm)[i], "_vs_", colnames(pm)[j]),
            test = "pairwise_Wilcoxon_BH", n = nrow(dsub), n_groups = 2,
            statistic = NA_real_, p_value = NA_real_, adj_p_value = pm[i, j], direction = "")
    }
  }
  if (length(group_results) > 0) {
    fwrite(rbindlist(group_results, fill = TRUE), file.path(stats_dir, "group_comparisons.tsv"), sep = "\t")
    cat("Wrote group_comparisons.tsv\n")
  }
}

# 3. Descriptive summary
desc_list <- list()
for (v in c("het_genomewide", "froh", "roh_total_bp", "longest_roh", "n_tracts")) {
  if (!v %in% names(d)) next
  x <- as.numeric(d[[v]]); x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) next
  desc_list[[length(desc_list) + 1]] <- data.table(
    variable = v, n = length(x), mean = mean(x), median = median(x),
    sd = sd(x), min = min(x), max = max(x),
    q25 = quantile(x, 0.25), q75 = quantile(x, 0.75))
}
if (length(desc_list) > 0) {
  fwrite(rbindlist(desc_list), file.path(stats_dir, "descriptive_summary.tsv"), sep = "\t")
  cat("Wrote descriptive_summary.tsv\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# CHROMOSOME-LEVEL STATS (if per_chr_roh_summary.tsv provided)
# ═══════════════════════════════════════════════════════════════════════════
if (!is.null(chr_file)) {
  chr_d <- fread(chr_file)
  for (col in c("roh_bp", "froh_chr", "n_tracts", "longest_tract"))
    if (col %in% names(chr_d)) chr_d[, (col) := as.numeric(get(col))]
  if (!is.null(anc_file))
    chr_d <- merge(chr_d, unique(d[, .(sample, ancestry)]), by = "sample", all.x = TRUE)
  if (!"ancestry" %in% names(chr_d)) chr_d[, ancestry := "all"]

  # Per-chr descriptive
  chr_desc <- list()
  for (v in c("froh_chr", "roh_bp", "n_tracts", "longest_tract")) {
    if (!v %in% names(chr_d)) next
    tmp <- chr_d[!is.na(get(v)), .(
      n = .N, mean = mean(get(v), na.rm = TRUE), median = median(get(v), na.rm = TRUE),
      sd = sd(get(v), na.rm = TRUE), min = min(get(v), na.rm = TRUE),
      max = max(get(v), na.rm = TRUE)
    ), by = chrom]
    tmp[, variable := v]
    chr_desc[[length(chr_desc) + 1]] <- tmp
  }
  if (length(chr_desc) > 0) {
    fwrite(rbindlist(chr_desc, fill = TRUE),
           file.path(stats_dir, "per_chr_descriptive_summary.tsv"), sep = "\t")
    cat("Wrote per_chr_descriptive_summary.tsv\n")
  }

  # Per-chr group comparisons
  if (length(groups) > 1 && "froh_chr" %in% names(chr_d)) {
    chr_grp <- list()
    for (chrom_val in unique(chr_d$chrom)) {
      dsub <- chr_d[chrom == chrom_val & !is.na(froh_chr) & !is.na(ancestry)]
      ng <- length(unique(dsub$ancestry))
      if (nrow(dsub) < 5 || ng < 2) next
      if (ng == 2) {
        wt <- tryCatch(wilcox.test(froh_chr ~ ancestry, data = dsub), error = function(e) NULL)
        if (!is.null(wt))
          chr_grp[[length(chr_grp) + 1]] <- data.table(
            chrom = chrom_val, variable = "froh_chr", test = "Wilcoxon",
            n = nrow(dsub), n_groups = ng, statistic = wt$statistic, p_value = wt$p.value)
      } else {
        kt <- tryCatch(kruskal.test(froh_chr ~ ancestry, data = dsub), error = function(e) NULL)
        if (!is.null(kt))
          chr_grp[[length(chr_grp) + 1]] <- data.table(
            chrom = chrom_val, variable = "froh_chr", test = "Kruskal_Wallis",
            n = nrow(dsub), n_groups = ng, statistic = kt$statistic, p_value = kt$p.value)
      }
    }
    if (length(chr_grp) > 0) {
      fwrite(rbindlist(chr_grp, fill = TRUE),
             file.path(stats_dir, "per_chr_group_comparisons.tsv"), sep = "\t")
      cat("Wrote per_chr_group_comparisons.tsv\n")
    }
  }
}

cat("Stats complete:", stats_dir, "\n")
