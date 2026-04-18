#!/usr/bin/env Rscript
# =============================================================================
# C01a_diag_summaries.R — Genome-wide summary plots (S01-S08)
# Runs fast, produces PNG summaries across all chromosomes.
#
# v8.5.2 — Plot overhaul:
#   G1: sys.frame → env var sourcing
#   Fixed: outdir variable bug (was using undefined 'outdir' in some places)
#   Themed: all palettes from diag_common
#   S07: remove fixed 0.70 annotation, make faint
#   S05: use adaptive thresholds in description
#
# Usage: Rscript C01a_diag_summaries.R <precomp_dir> <outdir> [chrom]
# =============================================================================

SCRIPT_DIR <- Sys.getenv("SNAKES_DIR",
  file.path(Sys.getenv("BASE", "."),
            "inversion_codebase_v8.5/MODULE_5A2_Discovery_Core/snakes"))
source(file.path(SCRIPT_DIR, "STEP_C01a_diag_common.R"))

cli <- parse_diag_args()
data <- load_diag_data(cli$precomp_dir, cli$chrom_filter)
precomp_list <- data$precomp_list
chroms <- data$chroms
inv_like_dt <- data$inv_like_dt

summ_dir <- file.path(cli$outdir, "summaries")
dir.create(summ_dir, recursive = TRUE, showWarnings = FALSE)

message("[diag] Generating genome-wide summaries...")

# =============================================================================
# S01: Inv-likeness violin per chromosome
# =============================================================================
if (nrow(inv_like_dt) > 0) {
  il_plot <- inv_like_dt[chrom %in% chroms & is.finite(inv_likeness)]
  il_plot[, chrom := factor(chrom, levels = chroms)]

  # Per-chromosome adaptive q75 (shown instead of fixed 0.90)
  chr_q75 <- il_plot[, .(q75 = quantile(inv_likeness, 0.75, na.rm = TRUE)), by = chrom]

  p <- ggplot(il_plot, aes(x = chrom, y = inv_likeness)) +
    geom_violin(fill = "#DBEAFE", color = "#3B82F6", linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "#93C5FD") +
    geom_point(data = chr_q75, aes(x = chrom, y = q75),
               shape = 18, size = 2, color = "#DC2626") +
    labs(x = NULL, y = "Inv-likeness",
         title = "Inversion-Likeness Distribution per Chromosome",
         subtitle = paste0(nrow(il_plot), " windows across ", length(chroms), " chromosomes"),
         caption = "Red diamonds = per-chr adaptive q75 (seed gate)\nViolin width scaled to window count") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S01_inv_likeness_genome.png"),
         p, width = 14, height = 6, dpi = DPI)
  message("[diag] S01_inv_likeness_genome.png")
}

# =============================================================================
# S02: Z-score distribution per chromosome
# =============================================================================
z_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  if ("max_abs_z" %in% names(dt))
    z_rows[[chr]] <- data.table(chrom = chr, z = abs(dt$max_abs_z))
}
if (length(z_rows) > 0) {
  z_all <- rbindlist(z_rows)[is.finite(z)]
  z_all[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(z_all, aes(x = z)) +
    geom_histogram(bins = 100, fill = "#DBEAFE", color = "#3B82F6", linewidth = 0.2) +
    geom_vline(xintercept = 1.2, linetype = "dashed", color = COL_Q85) +
    geom_vline(xintercept = 1.8, linetype = "dashed", color = COL_Q90) +
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "#2563EB") +
    facet_wrap(~ chrom, scales = "free_y", ncol = 7) +
    labs(x = "|Robust z|", y = "Count",
         title = "Robust Z-Score Distribution per Chromosome",
         subtitle = paste0(nrow(z_all), " windows | Thresholds: S1L=1.2 S1M=1.8 S1S=2.5")) +
    THEME_BASE +
    theme(strip.text = element_text(size = 6))

  ggsave(file.path(summ_dir, "S02_z_score_distribution.png"),
         p, width = 16, height = 10, dpi = DPI)
  message("[diag] S02_z_score_distribution.png")
}

# =============================================================================
# S03: Background baseline comparison
# =============================================================================
bg_rows <- list()
for (chr in chroms) {
  bg_q <- precomp_list[[chr]]$bg_continuity_quantiles
  if (!is.null(bg_q))
    bg_rows[[chr]] <- data.table(chrom = chr, quantile = names(bg_q), value = as.numeric(bg_q))
}
if (length(bg_rows) > 0) {
  bg_dt <- rbindlist(bg_rows)
  bg_dt[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(bg_dt, aes(x = chrom, y = value, color = quantile, group = quantile)) +
    geom_point(size = 2) + geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("50%" = "#6B7280", "75%" = "#93C5FD",
                                   "80%" = "#3B82F6", "85%" = COL_Q85,
                                   "90%" = COL_Q90, "95%" = COL_Q95)) +
    labs(x = NULL, y = "Background continuity",
         title = "Adaptive Threshold Baselines per Chromosome",
         subtitle = "S1S uses q95 (red), S1M uses q90 (teal), S1L uses q85 (amber)") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S03_bg_baseline_comparison.png"),
         p, width = 14, height = 6, dpi = DPI)
  message("[diag] S03_bg_baseline_comparison.png")
}

# =============================================================================
# S04: Window count barplot
# =============================================================================
wc_dt <- data.table(
  chrom = chroms,
  n_windows = vapply(chroms, function(chr) precomp_list[[chr]]$n_windows, integer(1))
)
wc_dt[, chrom := factor(chrom, levels = chroms)]

p <- ggplot(wc_dt, aes(x = chrom, y = n_windows)) +
  geom_col(fill = "#3B82F6") +
  geom_text(aes(label = n_windows), vjust = -0.3, size = 2.5) +
  labs(x = NULL, y = "Windows",
       title = "Window Count per Chromosome",
       subtitle = paste0("Total: ", sum(wc_dt$n_windows), " windows across ",
                        nrow(wc_dt), " chromosomes"),
       caption = "100-SNP windows with step-20 dense registry") +
  THEME_BASE +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(summ_dir, "S04_window_count_barplot.png"),
       p, width = 12, height = 5, dpi = DPI)
message("[diag] S04_window_count_barplot.png")

# =============================================================================
# S05: Seed eligibility per family per chr
# =============================================================================
seed_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  n <- nrow(dt)
  has_z <- "max_abs_z" %in% names(dt)
  has_nn <- "seed_nn_dist" %in% names(dt)

  # Use adaptive NN thresholds
  nn_vals <- if (has_nn) dt$seed_nn_dist[is.finite(dt$seed_nn_dist)] else numeric(0)
  nn_qq <- if (length(nn_vals) > 10) quantile(nn_vals, c(0.25, 0.50, 0.75), na.rm = TRUE)
           else c("25%" = 0.70, "50%" = 0.80, "75%" = 0.90)

  for (fam_info in list(
    list(name = "S1S", z_min = 2.5, nn_max = nn_qq["25%"]),
    list(name = "S1M", z_min = 1.8, nn_max = nn_qq["50%"]),
    list(name = "S1L", z_min = 1.2, nn_max = nn_qq["75%"])
  )) {
    n_eligible <- n
    if (has_z) n_eligible <- sum(abs(dt$max_abs_z) >= fam_info$z_min, na.rm = TRUE)
    if (has_nn) n_eligible <- min(n_eligible, sum(dt$seed_nn_dist < fam_info$nn_max, na.rm = TRUE))

    seed_rows[[length(seed_rows) + 1]] <- data.table(
      chrom = chr, family = fam_info$name,
      n_eligible = n_eligible, n_total = n,
      pct = round(100 * n_eligible / n, 1)
    )
  }
}
if (length(seed_rows) > 0) {
  seed_dt <- rbindlist(seed_rows)
  seed_dt[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(seed_dt, aes(x = chrom, y = pct, fill = family)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = PAL_FAMILY) +
    labs(x = NULL, y = "Eligible seeds (%)",
         title = "Seed Eligibility per Family (Adaptive NN Thresholds)",
         subtitle = "Based on z-score and per-chr adaptive NN quantiles",
         caption = "S1S: z\u22652.5 + NN<q25 | S1M: z\u22651.8 + NN<q50 | S1L: z\u22651.2 + NN<q75") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S05_seed_eligibility.png"),
         p, width = 14, height = 6, dpi = DPI)
  message("[diag] S05_seed_eligibility.png")
}

# =============================================================================
# S06: Eigenvalue ratio distribution
# =============================================================================
er_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (!is.null(l1_col) && !is.null(l2_col)) {
    ratio <- dt[[l1_col]] / dt[[l2_col]]
    er_rows[[chr]] <- data.table(chrom = chr, eigen_ratio = ratio[is.finite(ratio)])
  }
}
if (length(er_rows) > 0) {
  er_all <- rbindlist(er_rows)
  er_all[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(er_all, aes(x = chrom, y = pmin(eigen_ratio, 20))) +
    geom_violin(fill = "#FEF3C7", color = "#D97706", linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "#F59E0B") +
    geom_hline(yintercept = 3, linetype = "dashed", color = "#DC2626", linewidth = 0.3) +
    labs(x = NULL, y = "\u03BB1/\u03BB2 (capped at 20)",
         title = "Eigenvalue Ratio Distribution per Chromosome",
         subtitle = "High ratio = strong PC1 dominance = inversion signal",
         caption = "Inversions typically \u03BB1/\u03BB2 > 3 (red dashed)") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S06_eigenvalue_ratio_genome.png"),
         p, width = 14, height = 6, dpi = DPI)
  message("[diag] S06_eigenvalue_ratio_genome.png")
}

# =============================================================================
# S07: NN adaptive quantiles across chromosomes
# =============================================================================
nn_q_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  if ("seed_nn_dist" %in% names(dt)) {
    nn <- dt$seed_nn_dist[is.finite(dt$seed_nn_dist)]
    if (length(nn) > 10) {
      qq <- quantile(nn, c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
      for (qi in seq_along(qq))
        nn_q_rows[[length(nn_q_rows) + 1]] <- data.table(
          chrom = chr, quantile = names(qq)[qi], value = as.numeric(qq[qi]))
    }
  }
}
if (length(nn_q_rows) > 0) {
  nn_q_dt <- rbindlist(nn_q_rows)
  nn_q_dt[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(nn_q_dt, aes(x = chrom, y = value, color = quantile, group = quantile)) +
    geom_point(size = 2) + geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("10%" = "#D1D5DB", "25%" = "#2563EB",
                                   "50%" = "#059669", "75%" = "#D97706",
                                   "90%" = "#DC2626")) +
    # S07: fixed 0.70 very faint (alpha=0.15) instead of prominent
    geom_hline(yintercept = 0.70, linetype = "dotted", color = "#DC2626",
               linewidth = 0.3, alpha = 0.15) +
    labs(x = NULL, y = "NN Distance",
         title = "Adaptive NN Quantiles Across Chromosomes",
         subtitle = "Proposed: S1S<q25 (blue) S1M<q50 (teal) S1L<q75 (amber)",
         caption = "Adaptive quantiles calibrate to per-chromosome data range") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  tryCatch({
    ggsave(file.path(summ_dir, "S07_nn_adaptive_quantiles.png"),
           p, width = 14, height = 6, dpi = DPI)
    message("[diag] S07_nn_adaptive_quantiles.png")
  }, error = function(e) message("[FAIL] S07: ", e$message))
}

# =============================================================================
# S08: Seed eligibility ADAPTIVE vs FIXED comparison
# =============================================================================
seed_compare_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  n <- nrow(dt)
  if (!all(c("max_abs_z", "seed_nn_dist") %in% names(dt))) next
  nn <- dt$seed_nn_dist[is.finite(dt$seed_nn_dist)]
  if (length(nn) < 10) next
  qq <- quantile(nn, c(0.25, 0.50, 0.75), na.rm = TRUE)

  for (fam_info in list(
    list(name = "S1S", z_min = 2.5, nn_fixed = 0.70, nn_adapt = qq["25%"]),
    list(name = "S1M", z_min = 1.8, nn_fixed = 0.80, nn_adapt = qq["50%"]),
    list(name = "S1L", z_min = 1.2, nn_fixed = 0.90, nn_adapt = qq["75%"])
  )) {
    z_pass <- abs(dt$max_abs_z) >= fam_info$z_min
    seed_compare_rows[[length(seed_compare_rows) + 1]] <- data.table(
      chrom = chr, family = fam_info$name,
      fixed_eligible = sum(z_pass & dt$seed_nn_dist < fam_info$nn_fixed, na.rm = TRUE),
      adaptive_eligible = sum(z_pass & dt$seed_nn_dist < fam_info$nn_adapt, na.rm = TRUE),
      n_total = n
    )
  }
}
if (length(seed_compare_rows) > 0) {
  sc_dt <- rbindlist(seed_compare_rows)
  sc_dt[, chrom := factor(chrom, levels = chroms)]

  sc_long <- melt(sc_dt, id.vars = c("chrom", "family", "n_total"),
                  measure.vars = c("fixed_eligible", "adaptive_eligible"),
                  variable.name = "method", value.name = "n_eligible")
  sc_long[, pct := round(100 * n_eligible / n_total, 1)]
  sc_long[, method := fifelse(method == "fixed_eligible", "Fixed (0.70/0.80/0.90)",
                              "Adaptive (q25/q50/q75)")]

  p <- ggplot(sc_long, aes(x = chrom, y = pct, fill = method)) +
    geom_col(position = "dodge") +
    facet_wrap(~ family, ncol = 1) +
    scale_fill_manual(values = PAL_METHOD) +
    labs(x = NULL, y = "Eligible seeds (%)",
         title = "Seed Eligibility: Fixed vs Adaptive NN Thresholds",
         subtitle = "Fixed thresholds pass everything \u2014 adaptive actually filters",
         caption = "Blue = adaptive (per-chr quantile) | Grey = fixed") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(face = "bold"))

  tryCatch({
    ggsave(file.path(summ_dir, "S08_seed_fixed_vs_adaptive.png"),
           p, width = 14, height = 10, dpi = DPI)
    message("[diag] S08_seed_fixed_vs_adaptive.png")
  }, error = function(e) message("[FAIL] S08: ", e$message))
}

message("[DONE] Summaries -> ", summ_dir)
