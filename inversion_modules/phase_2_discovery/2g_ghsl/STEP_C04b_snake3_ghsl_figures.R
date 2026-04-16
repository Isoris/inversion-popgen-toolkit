#!/usr/bin/env Rscript

# =============================================================================
# STEP_C04b_snake3_ghsl_figures.R
#
# Diagnostic + publication figures for Snake 3 v3 GHSL.
#
# FIGURES:
#   G1: GHSL contrast profile (per chr)
#   G2: GHSL concordance heatmap (per candidate region)
#   G3: Het band resolution scatter (per candidate)
#   G4: Phase coverage heatmap (per chr)
#   G5: GHSL vs eigenvalue comparison (per chr)
#   G6: Pairwise karyotype divergence (per candidate)
#   G7: Per-candidate composite (publication)
#   G8: Genome-wide GHSL summary strip
#
# Usage:
#   Rscript STEP_C04b_snake3_ghsl_figures.R <snake3_outdir> <precomp_dir> \
#     <figure_outdir> [--candidates <file>] [--chrom <chr>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: see header")

snake3_dir  <- args[1]
precomp_dir <- args[2]
fig_outdir  <- args[3]

cand_file    <- NULL
chrom_filter <- NULL
i <- 4L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--candidates" && i < length(args)) { cand_file <- args[i+1]; i <- i+2L }
  else if (a == "--chrom" && i < length(args))  { chrom_filter <- args[i+1]; i <- i+2L }
  else { i <- i+1L }
}

dir.create(fig_outdir, recursive = TRUE, showWarnings = FALSE)
DPI <- 300

THEME_BASE <- theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "gray35"),
        plot.caption = element_text(size = 6, color = "gray55", hjust = 0),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

# =============================================================================
# LOAD DATA
# =============================================================================

message("[G-figs] Loading Snake 3 v3 outputs...")

track_dt <- fread(file.path(snake3_dir, "snake3v3_window_track.tsv.gz"))
scores_dt <- tryCatch(fread(file.path(snake3_dir, "snake3v3_sample_scores.tsv.gz")),
                       error = function(e) data.table())
het_dt <- tryCatch(fread(file.path(snake3_dir, "snake3v3_het_resolution.tsv.gz")),
                    error = function(e) data.table())

message("[G-figs] Track: ", nrow(track_dt), " windows")
message("[G-figs] Scores: ", nrow(scores_dt), " sample-window records")
message("[G-figs] Het res: ", nrow(het_dt), " records")

# Load precomp for PC1 data
rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
precomp_list <- list()
for (f in rds_files) { obj <- readRDS(f); precomp_list[[obj$chrom]] <- obj }
chroms <- names(precomp_list)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)

# Load candidates if available
cand_dt <- data.table()
if (!is.null(cand_file) && file.exists(cand_file)) {
  cand_dt <- fread(cand_file)
  message("[G-figs] Candidates: ", nrow(cand_dt))
}

# Inv-likeness for comparison
inv_like_file <- file.path(dirname(precomp_dir), "snake_inv_likeness.tsv.gz")
inv_like_dt <- if (file.exists(inv_like_file)) fread(inv_like_file) else data.table()

# =============================================================================
# G1: GHSL Contrast Profile (per chromosome)
# =============================================================================

message("[G-figs] G1: GHSL contrast profiles...")
g1_dir <- file.path(fig_outdir, "G1_contrast_profile")
dir.create(g1_dir, recursive = TRUE, showWarnings = FALSE)

for (chr in chroms) {
  chr_track <- track_dt[chrom == chr & is.finite(ghsl_contrast)]
  if (nrow(chr_track) == 0) next

  chr_track[, status_col := factor(ghsl_status, levels = c("PASS", "WEAK", "FAIL"))]

  p <- ggplot(chr_track, aes(x = pos_mb, y = ghsl_contrast, color = status_col)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("PASS" = "green4", "WEAK" = "darkorange", "FAIL" = "grey70"),
                        name = "GHSL status") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "green4", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "GHSL contrast (within - between)",
         title = paste0(chr, " -- GHSL Haplotype Contrast Profile"),
         subtitle = paste0(nrow(chr_track), " windows | PASS:", sum(chr_track$ghsl_status == "PASS"),
                          " WEAK:", sum(chr_track$ghsl_status == "WEAK")),
         caption = "Positive = groups have distinct haplotypes (inversion signal)\nGreen dotted = PASS threshold (0.05)") +
    THEME_BASE

  # Overlay candidate regions if available
  if (nrow(cand_dt) > 0) {
    chr_cand <- cand_dt[chrom == chr]
    if (nrow(chr_cand) > 0) {
      p <- p + geom_rect(data = chr_cand, inherit.aes = FALSE,
                          aes(xmin = start_bp / 1e6, xmax = end_bp / 1e6,
                              ymin = -Inf, ymax = Inf),
                          fill = "steelblue", alpha = 0.08)
    }
  }

  tryCatch(ggsave(file.path(g1_dir, paste0(chr, "_G1_ghsl_contrast.png")),
                   p, width = 14, height = 5, dpi = DPI),
           error = function(e) message("  [FAIL] G1 ", chr, ": ", e$message))
}

# =============================================================================
# G2: GHSL Concordance Heatmap (per candidate region)
# =============================================================================

message("[G-figs] G2: GHSL concordance heatmaps...")
g2_dir <- file.path(fig_outdir, "G2_concordance_heatmap")
dir.create(g2_dir, recursive = TRUE, showWarnings = FALSE)

# Find concordance matrix files
mat_dir <- file.path(snake3_dir, "concordance_matrices")
mat_files <- if (dir.exists(mat_dir)) list.files(mat_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE) else character(0)

if (length(mat_files) > 0 && nrow(scores_dt) > 0) {
  # Group matrices by chromosome and merge into region-level heatmaps
  for (chr in chroms) {
    chr_mats <- mat_files[grepl(paste0("^", chr, "_"), basename(mat_files))]
    if (length(chr_mats) == 0) next

    # Load and merge all concordance data for this chr
    all_conc <- rbindlist(lapply(chr_mats, function(f) tryCatch(fread(f), error = function(e) data.table())),
                           fill = TRUE)
    if (nrow(all_conc) == 0) next

    # Average concordance per sample pair across all PASS windows
    pair_avg <- all_conc[, .(mean_conc = mean(concordance, na.rm = TRUE),
                              n_windows = .N),
                          by = .(sample_i, sample_j)]

    # Make symmetric
    pair_sym <- rbind(pair_avg[, .(s1 = sample_i, s2 = sample_j, conc = mean_conc)],
                       pair_avg[, .(s1 = sample_j, s2 = sample_i, conc = mean_conc)])

    # Get band assignments for ordering
    chr_scores_sub <- scores_dt[chrom == chr]
    if (nrow(chr_scores_sub) > 0) {
      band_assign <- chr_scores_sub[, .(band = as.integer(names(which.max(table(pc1_band))))),
                                     by = sample_id]
      setkey(band_assign, sample_id)

      # Order samples by band
      all_samps <- sort(unique(c(pair_sym$s1, pair_sym$s2)))
      samp_order <- merge(data.table(sample_id = all_samps), band_assign,
                           by = "sample_id", all.x = TRUE)
      samp_order[is.na(band), band := 0L]
      samp_order <- samp_order[order(band, sample_id)]

      pair_sym[, s1 := factor(s1, levels = samp_order$sample_id)]
      pair_sym[, s2 := factor(s2, levels = samp_order$sample_id)]

      # Add self-concordance = 1
      self_dt <- data.table(s1 = samp_order$sample_id,
                             s2 = samp_order$sample_id, conc = 1.0)
      self_dt[, s1 := factor(s1, levels = samp_order$sample_id)]
      self_dt[, s2 := factor(s2, levels = samp_order$sample_id)]
      pair_sym <- rbind(pair_sym, self_dt)

      p <- ggplot(pair_sym, aes(x = s1, y = s2, fill = conc)) +
        geom_tile() +
        scale_fill_gradient2(low = "navy", mid = "white", high = "red3",
                             midpoint = 0.5, limits = c(0, 1), name = "Concordance") +
        # Band separator lines
        {
          band_breaks <- samp_order[, .N, by = band][order(band)]
          band_breaks[, cumN := cumsum(N)]
          if (nrow(band_breaks) > 1) {
            cuts <- band_breaks$cumN[-nrow(band_breaks)] + 0.5
            list(
              geom_hline(yintercept = cuts, color = "black", linewidth = 0.5),
              geom_vline(xintercept = cuts, color = "black", linewidth = 0.5)
            )
          }
        } +
        labs(title = paste0(chr, " -- GHSL Concordance Heatmap"),
             subtitle = paste0(length(all_samps), " samples | ",
                              length(chr_mats), " PASS windows averaged"),
             x = NULL, y = NULL,
             caption = "Ordered by PC1 band | Red = high concordance | Blue = low\nBlock structure = real inversion (distinct haplotype groups)") +
        coord_fixed() + THEME_BASE +
        theme(axis.text = element_blank(), axis.ticks = element_blank())

      tryCatch(ggsave(file.path(g2_dir, paste0(chr, "_G2_concordance_heatmap.png")),
                       p, width = 10, height = 9, dpi = DPI),
               error = function(e) message("  [FAIL] G2 ", chr, ": ", e$message))
    }
  }
} else {
  message("[G-figs] G2: no concordance matrices found (use --store_matrices in main engine)")
}

# =============================================================================
# G3: Het Band Resolution Scatter (per chromosome)
# =============================================================================

message("[G-figs] G3: Het band resolution...")
g3_dir <- file.path(fig_outdir, "G3_het_resolution")
dir.create(g3_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(het_dt) > 0) {
  for (chr in unique(het_dt$chrom)) {
    chr_het <- het_dt[chrom == chr]
    if (nrow(chr_het) < 5) next

    # Average per sample across windows
    het_avg <- chr_het[, .(mean_conc_b1 = mean(concordance_with_band1, na.rm = TRUE),
                            mean_conc_b3 = mean(concordance_with_band3, na.rm = TRUE),
                            mean_het_score = mean(het_score, na.rm = TRUE),
                            dominant_genotype = names(which.max(table(inferred_genotype)))),
                        by = sample_id]
    het_avg <- het_avg[is.finite(mean_conc_b1) & is.finite(mean_conc_b3)]
    if (nrow(het_avg) < 3) next

    p <- ggplot(het_avg, aes(x = mean_conc_b1, y = mean_conc_b3, color = dominant_genotype)) +
      geom_point(size = 2, alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("INV_HET" = "green4", "INV_LIKE_BAND1" = "blue4",
                                     "INV_LIKE_BAND3" = "red3", "UNKNOWN" = "grey60"),
                          name = "Inferred genotype") +
      labs(x = "Mean concordance with Band 1 (hom group A)",
           y = "Mean concordance with Band 3 (hom group B)",
           title = paste0(chr, " -- Het Band Resolution (Band 2 samples)"),
           subtitle = paste0(nrow(het_avg), " band-2 samples | averaged across GHSL windows"),
           caption = "Diagonal = equal concordance with both groups (true het)\nOff-diagonal = misassigned (belongs in band 1 or 3)") +
      THEME_BASE

    tryCatch(ggsave(file.path(g3_dir, paste0(chr, "_G3_het_resolution.png")),
                     p, width = 8, height = 7, dpi = DPI),
             error = function(e) message("  [FAIL] G3 ", chr, ": ", e$message))
  }
}

# =============================================================================
# G4: Phase Coverage Heatmap (per chromosome)
# =============================================================================

message("[G-figs] G4: Phase coverage heatmaps...")
g4_dir <- file.path(fig_outdir, "G4_phase_coverage")
dir.create(g4_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(scores_dt) > 0) {
  for (chr in chroms) {
    chr_sc <- scores_dt[chrom == chr]
    if (nrow(chr_sc) == 0) next

    # Get window positions
    chr_track_sub <- track_dt[chrom == chr, .(global_window_id, pos_mb)]

    # Merge
    cov_dt <- merge(chr_sc[, .(sample_id, global_window_id, n_het_phased, n_hom)],
                     chr_track_sub, by = "global_window_id")
    cov_dt[, total_sites := n_het_phased + n_hom]

    # Subsample windows for display (max 200)
    wids <- sort(unique(cov_dt$global_window_id))
    step_w <- max(1L, length(wids) %/% 200L)
    wids_sub <- wids[seq(1, length(wids), by = step_w)]
    cov_sub <- cov_dt[global_window_id %in% wids_sub]

    if (nrow(cov_sub) < 100) next

    p <- ggplot(cov_sub, aes(x = pos_mb, y = sample_id, fill = pmin(n_het_phased, 20))) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "navy", name = "Phased hets\n(capped 20)") +
      labs(x = paste0(chr, " (Mb)"), y = NULL,
           title = paste0(chr, " -- Phase Coverage Heatmap"),
           subtitle = paste0(length(unique(cov_sub$sample_id)), " samples × ",
                            length(wids_sub), " windows (subsampled)"),
           caption = "Dark = many phased het SNPs (good GHSL signal)\nWhite = no phasing (GHSL blind here)") +
      THEME_BASE +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

    tryCatch(ggsave(file.path(g4_dir, paste0(chr, "_G4_phase_coverage.png")),
                     p, width = 14, height = 8, dpi = DPI),
             error = function(e) message("  [FAIL] G4 ", chr, ": ", e$message))
  }
}

# =============================================================================
# G5: GHSL vs Eigenvalue Comparison (per chromosome)
# =============================================================================

message("[G-figs] G5: GHSL vs eigenvalue comparison...")
g5_dir <- file.path(fig_outdir, "G5_ghsl_vs_eigenvalue")
dir.create(g5_dir, recursive = TRUE, showWarnings = FALSE)

for (chr in chroms) {
  chr_track_sub <- track_dt[chrom == chr & is.finite(ghsl_contrast)]
  if (nrow(chr_track_sub) == 0) next

  # Get inv-likeness and het_contrast from precomp
  pc <- precomp_list[[chr]]
  dt <- pc$dt

  panels <- list()

  # Panel A: inv-likeness
  if ("inv_likeness" %in% names(dt)) {
    panels[["A: Inv-likeness"]] <- data.table(
      pos_mb = (dt$start_bp + dt$end_bp) / 2e6,
      value = dt$inv_likeness, panel = "A: Inv-likeness"
    )
  }

  # Panel B: GHSL contrast
  panels[["B: GHSL contrast"]] <- data.table(
    pos_mb = chr_track_sub$pos_mb,
    value = chr_track_sub$ghsl_contrast, panel = "B: GHSL contrast"
  )

  # Panel C: het_contrast (if available)
  het_col <- if ("inv_het_contrast" %in% names(dt)) "inv_het_contrast"
             else if ("het_contrast" %in% names(dt)) "het_contrast" else NULL
  if (!is.null(het_col)) {
    panels[["C: Het-contrast"]] <- data.table(
      pos_mb = (dt$start_bp + dt$end_bp) / 2e6,
      value = dt[[het_col]], panel = "C: Het-contrast"
    )
  }

  # Panel D: phase coverage
  panels[["D: Phase coverage"]] <- data.table(
    pos_mb = chr_track_sub$pos_mb,
    value = chr_track_sub$phase_coverage, panel = "D: Phase coverage"
  )

  pdt <- rbindlist(panels)
  pdt <- pdt[is.finite(value)]
  pdt[, panel := factor(panel, levels = sort(unique(panel)))]

  p <- ggplot(pdt, aes(x = pos_mb, y = value)) +
    geom_point(size = 0.4, alpha = 0.5, color = "steelblue") +
    geom_smooth(method = "loess", span = 0.1, se = FALSE, color = "red3",
                linewidth = 0.5, na.rm = TRUE) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " -- GHSL vs Eigenvalue Signals"),
         subtitle = "Red line = LOESS smooth | Coincident peaks = high-confidence inversion",
         caption = "A: eigenvalue-based | B: haplotype-based | C: het variance-based | D: data quality") +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold", size = 9),
          panel.spacing = unit(0.8, "lines"))

  tryCatch(ggsave(file.path(g5_dir, paste0(chr, "_G5_ghsl_vs_eigenvalue.png")),
                   p, width = 14, height = 12, dpi = DPI),
           error = function(e) message("  [FAIL] G5 ", chr, ": ", e$message))
}

# =============================================================================
# G6: Pairwise Karyotype Divergence (requires concordance matrices)
# =============================================================================

message("[G-figs] G6: Pairwise karyotype divergence...")
g6_dir <- file.path(fig_outdir, "G6_pairwise_divergence")
dir.create(g6_dir, recursive = TRUE, showWarnings = FALSE)

if (length(mat_files) > 0) {
  for (chr in chroms) {
    chr_mats <- mat_files[grepl(paste0("^", chr, "_"), basename(mat_files))]
    if (length(chr_mats) < 10) next

    # Load all concordance data
    all_conc <- rbindlist(lapply(chr_mats, function(f) tryCatch(fread(f), error = function(e) data.table())),
                           fill = TRUE)
    if (nrow(all_conc) == 0) next

    # Get window positions
    wid_pos <- track_dt[chrom == chr, .(global_window_id, pos_mb)]
    all_conc <- merge(all_conc, wid_pos, by = "global_window_id")

    # Find top divergent pairs (low concordance = different karyotypes)
    pair_means <- all_conc[, .(mean_conc = mean(concordance)), by = .(sample_i, sample_j)]
    divergent <- pair_means[order(mean_conc)][1:min(10, nrow(pair_means))]

    if (nrow(divergent) == 0) next

    # Plot concordance along chromosome for divergent pairs
    div_conc <- merge(all_conc, divergent[, .(sample_i, sample_j)],
                       by = c("sample_i", "sample_j"))
    div_conc[, pair_label := paste0(sample_i, " vs ", sample_j)]

    p <- ggplot(div_conc, aes(x = pos_mb, y = concordance, group = pair_label)) +
      geom_line(alpha = 0.3, linewidth = 0.3, color = "steelblue") +
      stat_summary(aes(group = 1), fun = mean, geom = "line",
                   color = "red3", linewidth = 1, na.rm = TRUE) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
      labs(x = paste0(chr, " (Mb)"), y = "Pairwise concordance",
           title = paste0(chr, " -- Karyotype Divergence Profile"),
           subtitle = paste0("Top ", nrow(divergent), " most divergent sample pairs"),
           caption = "Blue = individual pairs | Red = mean\nDrops in concordance mark inversion breakpoints") +
      THEME_BASE

    tryCatch(ggsave(file.path(g6_dir, paste0(chr, "_G6_pairwise_divergence.png")),
                     p, width = 14, height = 5, dpi = DPI),
             error = function(e) message("  [FAIL] G6 ", chr, ": ", e$message))
  }
}

# =============================================================================
# G8: Genome-wide GHSL Summary Strip
# =============================================================================

message("[G-figs] G8: Genome-wide summary strip...")

if (nrow(track_dt) > 0) {
  # Smooth GHSL contrast per chr for strip display
  strip_rows <- list()
  for (chr in chroms) {
    chr_t <- track_dt[chrom == chr & is.finite(ghsl_contrast)]
    if (nrow(chr_t) == 0) next
    chr_t[, chrom_f := factor(chr, levels = chroms)]
    strip_rows[[length(strip_rows) + 1L]] <- chr_t[, .(chrom = chr, pos_mb, ghsl_contrast, ghsl_status)]
  }

  if (length(strip_rows) > 0) {
    strip_dt <- rbindlist(strip_rows)
    chr_order <- chroms[order(as.integer(gsub("\\D", "", chroms)))]
    strip_dt[, chrom := factor(chrom, levels = chr_order)]

    p <- ggplot(strip_dt, aes(x = pos_mb, y = chrom, fill = pmin(pmax(ghsl_contrast, -0.1), 0.3))) +
      geom_tile(width = 0.1) +
      scale_fill_gradient2(low = "navy", mid = "grey95", high = "red3",
                           midpoint = 0, name = "GHSL contrast") +
      labs(x = "Position (Mb)", y = NULL,
           title = "Genome-wide GHSL Haplotype Contrast",
           subtitle = paste0(sum(track_dt$ghsl_status == "PASS", na.rm = TRUE), " PASS windows across ",
                            length(chroms), " chromosomes"),
           caption = "Red = high contrast (inversion signal) | Blue = negative (noise)\nGrey = no signal") +
      THEME_BASE +
      theme(axis.text.y = element_text(size = 7))

    tryCatch(ggsave(file.path(fig_outdir, "G8_genome_wide_ghsl_strip.png"),
                     p, width = 16, height = 8, dpi = DPI),
             error = function(e) message("  [FAIL] G8: ", e$message))
  }
}

message("\n[DONE] Snake 3 v3 figures -> ", fig_outdir)
