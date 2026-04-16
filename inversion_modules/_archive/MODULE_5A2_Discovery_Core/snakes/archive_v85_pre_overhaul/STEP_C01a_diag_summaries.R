#!/usr/bin/env Rscript
# =============================================================================
# C01a_diag_summaries.R — Genome-wide summary plots (S01-S06)
# Runs fast, produces PNG summaries across all chromosomes.
# Usage: Rscript C01a_diag_summaries.R <precomp_dir> <outdir> [chrom]
# =============================================================================

SCRIPT_DIR <- dirname(sys.frame(1)$ofile %||% ".")
source(file.path(SCRIPT_DIR, "C01a_diag_common.R"))

cli <- parse_diag_args()
data <- load_diag_data(cli$precomp_dir, cli$chrom_filter)
precomp_list <- data$precomp_list
chroms <- data$chroms
inv_like_dt <- data$inv_like_dt

summ_dir <- file.path(cli$outdir, "summaries")
dir.create(summ_dir, recursive = TRUE, showWarnings = FALSE)

message("[diag] Generating genome-wide summaries...")

# =============================================================================
# GENOME-WIDE SUMMARY PLOTS
# =============================================================================

message("[C01a_diag] Generating genome-wide summaries...")
summ_dir <- file.path(outdir, "summaries")
dir.create(summ_dir, recursive = TRUE, showWarnings = FALSE)

# -- S01: Inv-likeness violin per chromosome --------------------------
if (nrow(inv_like_dt) > 0) {
  il_plot <- inv_like_dt[chrom %in% chroms & is.finite(inv_likeness)]
  il_plot[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(il_plot, aes(x = chrom, y = inv_likeness)) +
    geom_violin(fill = "aliceblue", color = "blue4", linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "lightblue") +
    geom_hline(yintercept = 0.90, color = "red3", linetype = "dashed") +
    labs(x = NULL, y = "Inv-likeness",
         title = "Inversion-Likeness Distribution per Chromosome",
         subtitle = paste0(nrow(il_plot), " windows across ", length(chroms), " chromosomes"),
         caption = "Red dashed = seed gate (0.90)\nViolin width scaled to window count") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S01_inv_likeness_genome.png"),
         p, width = 14, height = 6, dpi = DPI)
  message("[C01a_diag] S01_inv_likeness_genome.png")
}

# -- S02: Z-score distribution ----------------------------------------
z_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  if ("max_abs_z" %in% names(dt)) {
    z_rows[[chr]] <- data.table(chrom = chr, z = abs(dt$max_abs_z))
  }
}
if (length(z_rows) > 0) {
  z_all <- rbindlist(z_rows)
  z_all <- z_all[is.finite(z)]
  z_all[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(z_all, aes(x = z)) +
    geom_histogram(bins = 100, fill = "lightblue", color = "blue4", linewidth = 0.2) +
    geom_vline(xintercept = 1.2, linetype = "dashed", color = "darkorange") +
    geom_vline(xintercept = 1.8, linetype = "dashed", color = "green4") +
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "blue4") +
    facet_wrap(~ chrom, scales = "free_y", ncol = 7) +
    labs(x = "|Robust z|", y = "Count",
         title = "Robust Z-Score Distribution per Chromosome",
         subtitle = paste0(nrow(z_all), " windows | Thresholds: S1L=1.2 S1M=1.8 S1S=2.5"),
         caption = "Note: y-axis free per chromosome") +
    THEME_BASE +
    theme(strip.text = element_text(size = 6))

  ggsave(file.path(summ_dir, "S02_z_score_distribution.png"),
         p, width = 16, height = 10, dpi = DPI)
  message("[C01a_diag] S02_z_score_distribution.png")
}

# -- S03: Background baseline comparison ------------------------------
bg_rows <- list()
for (chr in chroms) {
  bg_q <- precomp_list[[chr]]$bg_continuity_quantiles
  if (!is.null(bg_q)) {
    bg_rows[[chr]] <- data.table(
      chrom = chr,
      quantile = names(bg_q),
      value = as.numeric(bg_q)
    )
  }
}
if (length(bg_rows) > 0) {
  bg_dt <- rbindlist(bg_rows)
  bg_dt[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(bg_dt, aes(x = chrom, y = value, color = quantile, group = quantile)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("50%" = "gray40", "75%" = "lightblue",
                                   "80%" = "cornflowerblue", "85%" = "darkorange",
                                   "90%" = "green4", "95%" = "red3")) +
    labs(x = NULL, y = "Background continuity",
         title = "Adaptive Threshold Baselines per Chromosome",
         subtitle = "S1S uses q95 (red), S1M uses q90 (green), S1L uses q85 (orange)",
         caption = "Higher baseline = stricter adaptive threshold\nVariation reflects chromosome-specific PCA structure") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S03_bg_baseline_comparison.png"),
         p, width = 14, height = 6, dpi = DPI)
  message("[C01a_diag] S03_bg_baseline_comparison.png")
}

# -- S04: Window count barplot ----------------------------------------
wc_dt <- data.table(
  chrom = chroms,
  n_windows = vapply(chroms, function(chr) precomp_list[[chr]]$n_windows, integer(1))
)
wc_dt[, chrom := factor(chrom, levels = chroms)]

p <- ggplot(wc_dt, aes(x = chrom, y = n_windows)) +
  geom_col(fill = "royalblue") +
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
message("[C01a_diag] S04_window_count_barplot.png")

# -- S05: Seed eligibility per family per chr -------------------------
seed_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  n <- nrow(dt)
  has_z <- "max_abs_z" %in% names(dt)
  has_nn <- "seed_nn_dist" %in% names(dt)
  has_inv <- "inv_likeness" %in% names(dt)

  for (fam_info in list(
    list(name = "S1S", z_min = 2.5, nn_max = 0.70),
    list(name = "S1M", z_min = 1.8, nn_max = 0.80),
    list(name = "S1L", z_min = 1.2, nn_max = 0.90)
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
    scale_fill_manual(values = c("S1S" = "blue4", "S1M" = "green4", "S1L" = "darkorange")) +
    labs(x = NULL, y = "Eligible seeds (%)",
         title = "Seed Eligibility per Family per Chromosome",
         subtitle = "Based on z-score and NN distance thresholds (before inv-likeness gate)",
         caption = "S1S: z>=2.5 + NN<0.70 | S1M: z>=1.8 + NN<0.80 | S1L: z>=1.2 + NN<0.90") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S05_seed_eligibility.png"),
         p, width = 14, height = 6, dpi = DPI)
  message("[C01a_diag] S05_seed_eligibility.png")
}

# -- S06: Eigenvalue ratio distribution -------------------------------
er_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (!is.null(l1_col) && !is.null(l2_col)) {
    ratio <- dt[[l1_col]] / dt[[l2_col]]
    ratio <- ratio[is.finite(ratio)]
    er_rows[[chr]] <- data.table(chrom = chr, eigen_ratio = ratio)
  }
}
if (length(er_rows) > 0) {
  er_all <- rbindlist(er_rows)
  er_all[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(er_all, aes(x = chrom, y = pmin(eigen_ratio, 20))) +
    geom_violin(fill = "lemonchiffon", color = "darkorange", linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "gold") +
    labs(x = NULL, y = "lambda1/lambda2 (capped at 20)",
         title = "Eigenvalue Ratio Distribution per Chromosome",
         subtitle = "High ratio = strong PC1 dominance = inversion signal",
         caption = "Inversions typically produce lambda1/lambda2 > 3") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(summ_dir, "S06_eigenvalue_ratio_genome.png"),
         p, width = 14, height = 6, dpi = DPI)
  message("[C01a_diag] S06_eigenvalue_ratio_genome.png")
}

message("\n[DONE] STEP_C01a_diag complete -> ", outdir)

# -- S07: NN adaptive quantiles across chromosomes --------------------
nn_q_rows <- list()
for (chr in chroms) {
  dt <- precomp_list[[chr]]$dt
  if ("seed_nn_dist" %in% names(dt)) {
    nn <- dt$seed_nn_dist[is.finite(dt$seed_nn_dist)]
    if (length(nn) > 10) {
      qq <- quantile(nn, c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
      for (qi in seq_along(qq)) {
        nn_q_rows[[length(nn_q_rows) + 1]] <- data.table(
          chrom = chr, quantile = names(qq)[qi], value = as.numeric(qq[qi])
        )
      }
    }
  }
}
if (length(nn_q_rows) > 0) {
  nn_q_dt <- rbindlist(nn_q_rows)
  nn_q_dt[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(nn_q_dt, aes(x = chrom, y = value, color = quantile, group = quantile)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("10%" = "grey50", "25%" = "blue4",
                                   "50%" = "green4", "75%" = "darkorange",
                                   "90%" = "red3")) +
    geom_hline(yintercept = 0.70, linetype = "dotted", color = "red3", linewidth = 0.4) +
    annotate("text", x = 1, y = 0.72, label = "current fixed S1S=0.70",
             size = 2.5, color = "red3", hjust = 0) +
    labs(x = NULL, y = "NN Distance",
         title = "Adaptive NN Quantiles Across Chromosomes",
         subtitle = "Proposed: S1S<q25 (blue) S1M<q50 (green) S1L<q75 (orange)",
         caption = "Current fixed thresholds (0.70/0.80/0.90) are above ALL quantiles\nAdaptive quantiles would actually filter windows") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  tryCatch({
    ggsave(file.path(summ_dir, "S07_nn_adaptive_quantiles.png"),
           p, width = 14, height = 6, dpi = DPI)
    message("[C01a_diag] S07_nn_adaptive_quantiles.png")
  }, error = function(e) message("[FAIL] S07: ", e$message))
}

# -- S08: Seed eligibility ADAPTIVE vs FIXED comparison ----------------
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
    nn_pass_fixed <- dt$seed_nn_dist < fam_info$nn_fixed
    nn_pass_adapt <- dt$seed_nn_dist < fam_info$nn_adapt

    seed_compare_rows[[length(seed_compare_rows) + 1]] <- data.table(
      chrom = chr, family = fam_info$name,
      fixed_eligible = sum(z_pass & nn_pass_fixed, na.rm = TRUE),
      adaptive_eligible = sum(z_pass & nn_pass_adapt, na.rm = TRUE),
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
    scale_fill_manual(values = c("Fixed (0.70/0.80/0.90)" = "grey70",
                                  "Adaptive (q25/q50/q75)" = "steelblue")) +
    labs(x = NULL, y = "Eligible seeds (%)",
         title = "Seed Eligibility: Fixed vs Adaptive NN Thresholds",
         subtitle = "Fixed thresholds pass everything -- adaptive actually filters",
         caption = "Blue = adaptive (per-chr quantile), Grey = fixed\nAdaptive reduces eligible seeds by selecting truly clustered windows") +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(face = "bold"))

  tryCatch({
    ggsave(file.path(summ_dir, "S08_seed_fixed_vs_adaptive.png"),
           p, width = 14, height = 10, dpi = DPI)
    message("[C01a_diag] S08_seed_fixed_vs_adaptive.png")
  }, error = function(e) message("[FAIL] S08: ", e$message))
}

message("\n[DONE] STEP_C01a_diag (all plots) complete -> ", outdir)

message("[DONE] Summaries -> ", summ_dir)
