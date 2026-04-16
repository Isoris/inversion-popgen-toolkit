#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01a_diag.R  (v8.5)
#
# DIAGNOSTIC PLOTS for C01a precomputed data.
# Standalone -- reads .precomp.rds + snake_inv_likeness.tsv.gz
#
# Per-chromosome (PDF with one page per chr):
#   01_eigenvalue_spectrum.pdf       -- lambda1, lambda2, PVE1 across windows
#   02_inv_likeness_profile.pdf      -- inv_likeness score + components
#   03_max_abs_z_profile.pdf          -- robust z-score (median/MAD) + seed threshold
#   04_bg_continuity_baseline.pdf    -- background continuity distribution + quantiles
#   05_nn_distance_profile.pdf       -- nearest-neighbor distance vs position
#   06_sim_mat_heatmap.pdf           -- similarity matrix (subsampled for large chr)
#   07_mds_scatter_raw.pdf           -- MDS1 vs MDS2, no snake coloring
#   16_simmat_zscore_contrast.pdf    -- row-wise Z-scored sim_mat (enhanced contrast)
#   17_simmat_quantile_contrast.pdf  -- row-wise quantile-normalized sim_mat
#   18_mds_by_het_contrast.pdf       -- MDS colored by het_contrast (v8.5)
#   19_mds_by_lambda2.pdf            -- MDS colored by lambda2 (low=inversion)
#   20_mds_by_pve1_excess.pdf        -- MDS colored by PVE1 excess over chr median
#   21_mds_by_nn_distance.pdf        -- MDS colored by NN distance (low=clustered)
#   22_mds_by_eigenvalue_ratio.pdf   -- MDS colored by lambda1/lambda2
#   23_mds_faceted_all_metrics.pdf   -- 6-panel faceted MDS (all metrics at once)
#
# Genome-wide summary (PNG):
#   S01_inv_likeness_genome.png      -- all chr inv_likeness violin/boxplot
#   S02_z_score_distribution.png     -- genome-wide z histogram + per-chr
#   S03_bg_baseline_comparison.png   -- bg_q50/q85/q90/q95 across chromosomes
#   S04_window_count_barplot.png     -- windows per chromosome
#   S05_seed_eligibility.png         -- seeds eligible per family per chr
#   S06_eigenvalue_ratio_genome.png  -- lambda1/lambda2 ratio distribution
#
# Usage:
#   Rscript STEP_C01a_diag.R <precomp_dir> <outdir> [chrom]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript STEP_C01a_diag.R <precomp_dir> <outdir> [chrom]")

precomp_dir <- args[1]
outdir      <- args[2]
chrom_filter <- if (length(args) >= 3 && args[3] != "all") args[3] else NULL

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DPI <- 350
THEME_BASE <- theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "gray35"),
        plot.caption = element_text(size = 6, color = "gray55", hjust = 0),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

# Color palettes
PAL_INV <- c("PASS" = "green4", "WEAK" = "orange", "FAIL" = "red3")
PAL_Z   <- c("z>=5" = "navy", "z>=4" = "blue4", "z>=3" = "royalblue",
             "z>=2" = "cornflowerblue", "z>=1.5" = "lightblue", "z<1.5" = "gray90")

# =============================================================================
# LOAD DATA
# =============================================================================

message("[C01a_diag] Loading precomputed data...")
rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

precomp_list <- list()
for (f in rds_files) {
  obj <- readRDS(f)
  precomp_list[[obj$chrom]] <- obj
}

chroms <- names(precomp_list)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)
message("[C01a_diag] ", length(chroms), " chromosomes")

# Inv-likeness table
inv_like_file <- file.path(dirname(precomp_dir), "snake_inv_likeness.tsv.gz")
inv_like_dt <- if (file.exists(inv_like_file)) fread(inv_like_file) else data.table()
if (nrow(inv_like_dt) > 0) message("[C01a_diag] Inv-likeness: ", nrow(inv_like_dt), " windows")

# =============================================================================
# PER-CHROMOSOME PLOTS (PDF, one page per chr)
# =============================================================================

message("[C01a_diag] Generating per-chromosome plots...")

# -- 01: Eigenvalue spectrum ------------------------------------------
build_eigen_profile <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  # Column names: lam_1/lam_2 from STEP10_v2, or lambda1/lambda2
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (is.null(l1_col) || is.null(l2_col)) return(NULL)

  pdt <- data.table(
    pos_mb = (dt$start_bp + dt$end_bp) / 2e6,
    lambda1 = dt[[l1_col]],
    lambda2 = dt[[l2_col]]
  )
  pdt[, pve1 := lambda1 / (lambda1 + lambda2)]
  pdt <- pdt[is.finite(lambda1) & is.finite(lambda2) & lambda1 > 0]
  if (nrow(pdt) < 10) return(NULL)

  lam_max <- max(pdt$lambda1, na.rm = TRUE)
  if (!is.finite(lam_max) || lam_max == 0) lam_max <- 1

  p <- ggplot(pdt) +
    geom_line(aes(x = pos_mb, y = lambda1), color = "blue4", linewidth = 0.3, alpha = 0.7) +
    geom_line(aes(x = pos_mb, y = lambda2), color = "red3", linewidth = 0.3, alpha = 0.5) +
    geom_line(aes(x = pos_mb, y = pve1 * lam_max),
              color = "green4", linewidth = 0.4, alpha = 0.6) +
    scale_y_continuous(sec.axis = sec_axis(~ . / lam_max, name = "PVE1 (green)")) +
    labs(x = paste0(chr, " (Mb)"), y = "Eigenvalue",
         title = paste0(chr, " -- Eigenvalue Spectrum"),
         subtitle = paste0("Blue: lambda1, Red: lambda2, Green: PVE1\n",
                          nrow(pdt), " windows"),
         caption = paste0("Median PVE1: ", round(median(pdt$pve1, na.rm = TRUE), 3),
                         " | lambda1/lambda2 median: ",
                         round(median(pdt$lambda1/pdt$lambda2, na.rm = TRUE), 2))) +
    THEME_BASE
  p
}

# -- 02: Inv-likeness profile ----------------------------------------
build_inv_likeness_profile <- function(chr) {
  if (nrow(inv_like_dt) == 0) return(NULL)
  il_chr <- inv_like_dt[chrom == chr]
  if (nrow(il_chr) == 0) return(NULL)

  # inv_like_dt may not have start_bp if generated by older C01a -- join if needed
  if (!"start_bp" %in% names(il_chr)) {
    dt <- precomp_list[[chr]]$dt
    if ("global_window_id" %in% names(dt) && "global_window_id" %in% names(il_chr)) {
      il_chr <- merge(il_chr, dt[, .(global_window_id, start_bp, end_bp)],
                       by = "global_window_id", all.x = TRUE)
    } else {
      return(NULL)
    }
  }

  il_chr <- il_chr[is.finite(start_bp) & is.finite(inv_likeness)]
  if (nrow(il_chr) == 0) return(NULL)

  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]
  il_chr[, status := fifelse(inv_likeness >= 0.5, "PASS",
                     fifelse(inv_likeness >= 0.25, "WEAK", "FAIL"))]

  chr_med <- median(il_chr$inv_likeness, na.rm = TRUE)
  chr_q75 <- quantile(il_chr$inv_likeness, 0.75, na.rm = TRUE)

  p <- ggplot(il_chr, aes(x = pos_mb, y = inv_likeness, color = status)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_color_manual(values = PAL_INV) +
    geom_hline(yintercept = chr_med, linetype = "dashed", color = "gray40", linewidth = 0.3) +
    geom_hline(yintercept = chr_q75, linetype = "dotted", color = "gray20", linewidth = 0.3) +
    geom_hline(yintercept = 0.90, linetype = "solid", color = "red3", linewidth = 0.3, alpha = 0.5) +
    annotate("text", x = min(il_chr$pos_mb), y = chr_med + 0.02,
             label = paste0("median=", round(chr_med, 3)), size = 2.5, hjust = 0) +
    annotate("text", x = min(il_chr$pos_mb), y = chr_q75 + 0.02,
             label = paste0("q75=", round(chr_q75, 3)), size = 2.5, hjust = 0) +
    labs(x = paste0(chr, " (Mb)"), y = "Inv-likeness",
         title = paste0(chr, " -- Inversion-Likeness Score"),
         subtitle = paste0(nrow(il_chr), " windows | ",
                          sum(il_chr$status == "PASS"), " PASS, ",
                          sum(il_chr$status == "WEAK"), " WEAK, ",
                          sum(il_chr$status == "FAIL"), " FAIL"),
         caption = paste0("Components: PVE1_excess (0.25) + eigen_ratio (0.20) + trimodality (0.25) + het_contrast (0.30)\n",
                         "Seed gate at 0.90 (red line) | chr median (dashed) | q75 (dotted)")) +
    THEME_BASE
  p
}

# -- 03: Robust z-score profile --------------------------------------
build_z_profile <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  if (!"max_abs_z" %in% names(dt)) return(NULL)

  pdt <- dt[, .(pos_mb = (start_bp + end_bp) / 2e6, z = abs(max_abs_z))]
  pdt <- pdt[is.finite(z)]

  pdt[, z_class := fifelse(z >= 5, "z>=5",
                   fifelse(z >= 4, "z>=4",
                   fifelse(z >= 3, "z>=3",
                   fifelse(z >= 2, "z>=2",
                   fifelse(z >= 1.5, "z>=1.5", "z<1.5")))))]

  p <- ggplot(pdt, aes(x = pos_mb, y = z, color = z_class)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_color_manual(values = PAL_Z) +
    geom_hline(yintercept = 1.2, linetype = "dashed", color = "darkorange", linewidth = 0.3) +
    geom_hline(yintercept = 1.8, linetype = "dashed", color = "green4", linewidth = 0.3) +
    geom_hline(yintercept = 2.5, linetype = "dashed", color = "blue4", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "|Robust z|",
         title = paste0(chr, " -- Robust Z-Score Profile"),
         subtitle = paste0(nrow(pdt), " windows | z>=2.5: ", sum(pdt$z >= 2.5),
                          " | z>=1.8: ", sum(pdt$z >= 1.8),
                          " | z>=1.2: ", sum(pdt$z >= 1.2)),
         caption = "Thresholds: S1S=2.5 (blue) | S1M=1.8 (green) | S1L=1.2 (orange)") +
    THEME_BASE
  p
}

# -- 04: Background continuity baseline ------------------------------
build_bg_baseline <- function(chr) {
  pc <- precomp_list[[chr]]
  bg_q <- pc$bg_continuity_quantiles
  if (is.null(bg_q) || length(bg_q) < 4) return(NULL)

  # Ensure names are accessible (might be "50%" or just numeric indices)
  if (is.null(names(bg_q))) names(bg_q) <- paste0(c(50, 75, 80, 85, 90, 95), "%")[seq_along(bg_q)]

  # Adjacent-pair sim_mat values for histogram
  sim_mat <- pc$sim_mat
  n_w <- nrow(sim_mat)
  if (n_w < 5) return(NULL)
  adj_sims <- vapply(seq_len(n_w - 1), function(i) {
    s <- sim_mat[i, i + 1L]
    if (is.finite(s)) s else NA_real_
  }, numeric(1))
  adj_sims <- adj_sims[is.finite(adj_sims)]
  if (length(adj_sims) < 10) return(NULL)

  pdt <- data.table(sim = adj_sims)

  p <- ggplot(pdt, aes(x = sim)) +
    geom_histogram(bins = 80, fill = "lightblue", color = "blue4", linewidth = 0.2) +
    geom_vline(xintercept = as.numeric(bg_q["50%"]), color = "gray40", linetype = "dashed") +
    geom_vline(xintercept = as.numeric(bg_q["85%"]), color = "darkorange", linetype = "dashed") +
    geom_vline(xintercept = as.numeric(bg_q["90%"]), color = "green4", linetype = "dashed") +
    geom_vline(xintercept = as.numeric(bg_q["95%"]), color = "red3", linetype = "dashed") +
    annotate("text", x = as.numeric(bg_q["85%"]), y = Inf, vjust = 1.5,
             label = paste0("q85=", round(bg_q["85%"], 3)), size = 2.5, color = "darkorange") +
    annotate("text", x = as.numeric(bg_q["90%"]), y = Inf, vjust = 3,
             label = paste0("q90=", round(bg_q["90%"], 3)), size = 2.5, color = "green4") +
    annotate("text", x = as.numeric(bg_q["95%"]), y = Inf, vjust = 4.5,
             label = paste0("q95=", round(bg_q["95%"], 3)), size = 2.5, color = "red3") +
    labs(x = "Adjacent-window similarity", y = "Count",
         title = paste0(chr, " -- Background Continuity Baseline"),
         subtitle = paste0("From adjacent below-median-inv window pairs\n",
                          "Adaptive thresholds: S1S=q95 S1M=q90 S1L=q85"),
         caption = paste0("q50=", round(bg_q["50%"], 3),
                         " q75=", round(bg_q["75%"], 3),
                         " q85=", round(bg_q["85%"], 3),
                         " q90=", round(bg_q["90%"], 3),
                         " q95=", round(bg_q["95%"], 3))) +
    THEME_BASE
  p
}

# -- 05: NN distance profile (adaptive + smoothed regime bands) -------
build_nn_profile <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  if (!"seed_nn_dist" %in% names(dt)) return(NULL)

  pdt <- dt[, .(pos_mb = (start_bp + end_bp) / 2e6, nn = seed_nn_dist)]
  pdt <- pdt[is.finite(nn)]
  if (nrow(pdt) < 30) return(NULL)

  # Adaptive quantiles
  qq <- quantile(pdt$nn, c(0.25, 0.50, 0.75), na.rm = TRUE)
  y_cap <- quantile(pdt$nn, 0.99, na.rm = TRUE) * 1.1  # cap y-axis at q99 + 10%

  # Smoothed NN for regime detection
  smooth_w <- min(31L, nrow(pdt) %/% 10L)
  if (smooth_w %% 2 == 0) smooth_w <- smooth_w + 1L
  smooth_w <- max(5L, smooth_w)
  pdt[, nn_smooth := frollmedian(nn, n = smooth_w, align = "center", na.rm = TRUE)]

  # Regime segmentation on smoothed track
  pdt[, nn_state := fifelse(nn_smooth <= qq["25%"], "core_low",
                   fifelse(nn_smooth <= qq["50%"], "mid", "background"))]
  pdt[, nn_state := nafill(nn_state, type = "locf")]
  pdt[is.na(nn_state), nn_state := "background"]
  pdt[, run_id := rleid(nn_state)]

  runs <- pdt[, .(state = first(nn_state), start_mb = min(pos_mb),
                   end_mb = max(pos_mb), n_win = .N), by = run_id]
  runs <- runs[n_win >= 10]  # drop tiny runs

  p <- ggplot(pdt, aes(x = pos_mb, y = nn)) +
    # Regime background bands
    {if (nrow(runs) > 0) geom_rect(
      data = runs, inherit.aes = FALSE,
      aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf, fill = state),
      alpha = 0.12
    )} +
    scale_fill_manual(values = c("core_low" = "blue3", "mid" = "darkseagreen3",
                                  "background" = "grey80"), name = "Regime") +
    # Raw points
    geom_point(size = 0.3, alpha = 0.35, color = "slateblue") +
    # Smoothed line
    geom_line(aes(y = nn_smooth), color = "black", linewidth = 0.5, na.rm = TRUE) +
    # Adaptive thresholds
    geom_hline(yintercept = qq["25%"], linetype = "dashed", color = "blue4", linewidth = 0.4) +
    geom_hline(yintercept = qq["50%"], linetype = "dashed", color = "green4", linewidth = 0.4) +
    geom_hline(yintercept = qq["75%"], linetype = "dashed", color = "darkorange", linewidth = 0.4) +
    labs(x = paste0(chr, " (Mb)"), y = "NN Distance (MDS-space)",
         title = paste0(chr, " -- Nearest-Neighbor Distance Along Chromosome"),
         subtitle = paste0(nrow(pdt), " windows (100-SNP, step-20) | ",
                          "Adaptive quantiles: q25=", round(qq["25%"], 3),
                          " q50=", round(qq["50%"], 3), " q75=", round(qq["75%"], 3)),
         caption = paste0("Source: MDS embedding of local PCA similarity matrix (chunked_2x, 20 dims)\n",
                         "Low NN = tightly clustered windows = consistent PCA signal (inversion core)\n",
                         "Blue band = core_low (<q25) | Green = mid (q25-q50) | Grey = background\n",
                         "Black line = rolling median (w=", smooth_w, ")")) +
    coord_cartesian(ylim = c(0, y_cap)) +
    THEME_BASE
  p
}

# -- 06: Similarity matrix heatmap -----------------------------------
build_sim_heatmap <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  n <- nrow(sim_mat)
  if (n < 10) return(NULL)

  # Subsample for large chromosomes (max 500x500)
  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  sub <- sim_mat[idx, idx]
  pos_mb <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6

  pdt <- data.table(
    x = rep(pos_mb, each = length(idx)),
    y = rep(pos_mb, times = length(idx)),
    sim = as.vector(sub)
  )
  pdt <- pdt[is.finite(sim)]

  # Center color on chr median for better contrast
  sim_med <- median(pdt$sim, na.rm = TRUE)

  p <- ggplot(pdt, aes(x = x, y = y, fill = sim)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("navy", "blue3", "steelblue", "gold", "orange", "red3"),
      values = scales::rescale(c(0, sim_med * 0.6, sim_med, sim_med * 1.1, sim_med * 1.3, 1)),
      name = "Similarity") +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " -- Similarity Matrix (", n, " windows)"),
         subtitle = paste0(length(idx), "x", length(idx), " subsampled | ",
                          "median sim = ", round(sim_med, 3)),
         caption = paste0("Source: MDS-space correlation between local PCA loading vectors\n",
                         "Red blocks on diagonal = windows sharing PCA structure (inversion candidates)\n",
                         "Off-diagonal blocks = long-range similarity (shared inversion across distant windows)\n",
                         "Blue = low similarity (background or boundary between systems)")) +
    coord_fixed() + THEME_BASE
  p
}

# -- 07: MDS scatter (raw, no snake coloring) ------------------------
build_mds_raw <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     pos_mb = (dt$start_bp + dt$end_bp) / 2e6)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2)]

  p <- ggplot(pdt, aes(x = MDS1, y = MDS2, color = pos_mb)) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_viridis_c(name = "Position (Mb)") +
    labs(title = paste0(chr, " -- MDS1 vs MDS2 (Raw)"),
         subtitle = paste0(nrow(pdt), " windows, colored by genomic position"),
         caption = "Clusters = candidate inversion regions\nOutliers = strong PCA signal windows") +
    THEME_BASE
  p
}

# -- 07b: Lostruct-style MDS1 + MDS2 vs position (3 panels) ----------
build_mds_vs_position <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)

  pos_mb <- (dt$start_bp + dt$end_bp) / 2e6
  pdt <- data.table(
    pos_mb = rep(pos_mb, 3),
    value = c(mds_mat[, 1], mds_mat[, 2], pos_mb),
    panel = rep(c("MDS1 vs position", "MDS2 vs position", "MDS1 vs MDS2"),
                each = nrow(dt)),
    x_val = c(pos_mb, pos_mb, mds_mat[, 1]),
    color_pos = rep(pos_mb, 3)
  )
  # For the MDS1vsMDS2 panel, y = MDS2
  pdt[panel == "MDS1 vs MDS2", value := mds_mat[, 2]]
  pdt <- pdt[is.finite(value) & is.finite(x_val)]

  # Use 3 colors via kmeans on MDS1 for cluster coloring (like lostruct paper)
  km <- tryCatch(kmeans(mds_mat[, 1:2], centers = 3, nstart = 10), error = function(e) NULL)
  if (!is.null(km)) {
    cluster_col <- rep(km$cluster, 3)
    pdt[, cluster := as.factor(cluster_col[seq_len(.N)])]
    pal_km <- c("1" = "black", "2" = "green4", "3" = "darkorange")

    p <- ggplot(pdt, aes(x = x_val, y = value, color = cluster)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_manual(values = pal_km, guide = "none") +
      facet_wrap(~ panel, ncol = 3, scales = "free") +
      labs(title = paste0(chr, " -- MDS Coordinates vs Position (Lostruct-style)"),
           subtitle = paste0(nrow(dt), " windows | k=3 clustering on MDS1+MDS2"),
           caption = "Left/center: MDS coords along chromosome show inversion breakpoints\nRight: MDS scatter shows cluster structure") +
      THEME_BASE +
      theme(strip.text = element_text(face = "bold", size = 8))
  } else {
    p <- ggplot(pdt, aes(x = x_val, y = value, color = color_pos)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_viridis_c(name = "Position (Mb)") +
      facet_wrap(~ panel, ncol = 3, scales = "free") +
      labs(title = paste0(chr, " -- MDS Coordinates vs Position"),
           subtitle = paste0(nrow(dt), " windows")) +
      THEME_BASE
  }
  p
}

# -- 08: NN distance histogram with adaptive quantiles ----------------
build_nn_histogram <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  if (!"seed_nn_dist" %in% names(dt)) return(NULL)
  nn <- dt$seed_nn_dist[is.finite(dt$seed_nn_dist)]
  if (length(nn) < 20) return(NULL)

  qq <- quantile(nn, c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
  x_cap <- quantile(nn, 0.99, na.rm = TRUE) * 1.1
  pdt <- data.table(nn = nn[nn <= x_cap])  # trim outliers for display

  p <- ggplot(pdt, aes(x = nn)) +
    geom_histogram(bins = 80, fill = "lavender", color = "slateblue", linewidth = 0.2) +
    geom_vline(xintercept = qq["25%"], color = "blue4", linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = qq["50%"], color = "green4", linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = qq["75%"], color = "darkorange", linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = 0.70, color = "red3", linetype = "dotted", linewidth = 0.3) +
    annotate("text", x = qq["25%"], y = Inf, vjust = 1.5,
             label = paste0("q25=", round(qq["25%"], 3)), size = 2.5, color = "blue4") +
    annotate("text", x = qq["50%"], y = Inf, vjust = 3,
             label = paste0("q50=", round(qq["50%"], 3)), size = 2.5, color = "green4") +
    annotate("text", x = qq["75%"], y = Inf, vjust = 4.5,
             label = paste0("q75=", round(qq["75%"], 3)), size = 2.5, color = "darkorange") +
    labs(x = "NN Distance (MDS-space)", y = "Count",
         title = paste0(chr, " -- NN Distance Distribution"),
         subtitle = paste0(length(nn), " windows | trimmed at q99 for display\n",
                          "Proposed adaptive seed thresholds: S1S<q25 S1M<q50 S1L<q75"),
         caption = paste0("Source: pairwise distance in MDS embedding between adjacent windows\n",
                         "q10=", round(qq["10%"], 3),
                         " q25=", round(qq["25%"], 3),
                         " q50=", round(qq["50%"], 3),
                         " q75=", round(qq["75%"], 3),
                         " q90=", round(qq["90%"], 3),
                         "\nRed dotted = current fixed threshold (0.70) -- above all data")) +
    THEME_BASE
  p
}

# -- 09: MDS colored by inv_likeness ----------------------------------
build_mds_inv <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  if (!"inv_likeness" %in% names(dt)) {
    # Try joining from inv_like_dt
    if (nrow(inv_like_dt) > 0 && "global_window_id" %in% names(dt)) {
      dt <- merge(dt, inv_like_dt[chrom == chr, .(global_window_id, inv_likeness)],
                   by = "global_window_id", all.x = TRUE)
    } else return(NULL)
  }

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     inv = dt$inv_likeness)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(inv)]
  if (nrow(pdt) < 10) return(NULL)

  # Center color on per-chr median for contrast
  inv_med <- median(pdt$inv, na.rm = TRUE)

  p <- ggplot(pdt[order(inv)], aes(x = MDS1, y = MDS2, color = inv)) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_gradientn(
      colours = c("navy", "steelblue", "grey80", "orange", "red3"),
      values = scales::rescale(c(0, inv_med * 0.7, inv_med, inv_med * 1.2, 1)),
      name = "Inv-likeness"
    ) +
    labs(title = paste0(chr, " -- MDS colored by Inv-Likeness"),
         subtitle = paste0(nrow(pdt), " windows"),
         caption = "Red = high inv-likeness (inversion signal)\nGrey = low (background)") +
    THEME_BASE
  p
}

# -- 10: MDS colored by z-score ---------------------------------------
build_mds_z <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  if (!"max_abs_z" %in% names(dt)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     z = dt$max_abs_z)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(z)]
  if (nrow(pdt) < 10) return(NULL)

  p <- ggplot(pdt[order(z)], aes(x = MDS1, y = MDS2, color = pmin(z, 5))) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_gradient2(low = "grey80", mid = "steelblue", high = "navy",
                          midpoint = 2.5, name = "|z| (capped 5)") +
    labs(title = paste0(chr, " -- MDS colored by Z-Score"),
         subtitle = paste0(nrow(pdt), " windows"),
         caption = "Dark blue = high z (strong MDS outlier)\nGrey = low z (background)") +
    THEME_BASE
  p
}

# -- 11: NN threshold calibration plot --------------------------------
build_nn_calibration <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  if (!"seed_nn_dist" %in% names(dt)) return(NULL)

  pdt <- dt[is.finite(seed_nn_dist), .(pos_mb = (start_bp + end_bp) / 2e6,
                                        nn = seed_nn_dist)]
  if (nrow(pdt) < 20) return(NULL)

  qq <- quantile(pdt$nn, c(0.10, 0.25, 0.50, 0.75), na.rm = TRUE)
  y_cap <- quantile(pdt$nn, 0.99, na.rm = TRUE) * 1.1

  # Color by which adaptive threshold each window would pass
  pdt[, adaptive_class := fifelse(nn <= qq["25%"], "S1S (q25)",
                          fifelse(nn <= qq["50%"], "S1M (q50)",
                          fifelse(nn <= qq["75%"], "S1L (q75)", "above q75")))]

  pal_adapt <- c("S1S (q25)" = "blue4", "S1M (q50)" = "green4",
                 "S1L (q75)" = "darkorange", "above q75" = "grey70")

  p <- ggplot(pdt, aes(x = pos_mb, y = nn, color = adaptive_class)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_color_manual(values = pal_adapt) +
    geom_hline(yintercept = qq["25%"], linetype = "dashed", color = "blue4", linewidth = 0.3) +
    geom_hline(yintercept = qq["50%"], linetype = "dashed", color = "green4", linewidth = 0.3) +
    geom_hline(yintercept = qq["75%"], linetype = "dashed", color = "darkorange", linewidth = 0.3) +
    geom_hline(yintercept = 0.70, linetype = "dotted", color = "red3", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "NN Distance (MDS-space)",
         title = paste0(chr, " -- Seed Threshold Calibration (NN Distance)"),
         subtitle = paste0("Adaptive quantiles: q25=", round(qq["25%"], 3),
                          " q50=", round(qq["50%"], 3),
                          " q75=", round(qq["75%"], 3),
                          " | y-axis capped at q99"),
         caption = paste0("Source: per-window nearest-neighbor distance in MDS embedding\n",
                         "Blue = S1S candidates (<q25, most clustered) | ",
                         "Green = S1M (<q50) | Orange = S1L (<q75)\n",
                         "Red dotted = current fixed 0.70 (above all data, non-discriminative)")) +
    coord_cartesian(ylim = c(0, y_cap)) +
    THEME_BASE
  p
}

# -- 12: Combined 4-panel profile ------------------------------------
build_combined_profile <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  n <- nrow(dt)
  pos_mb <- (dt$start_bp + dt$end_bp) / 2e6

  # Build panel data
  panels <- list()

  # Panel A: eigenvalue
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (!is.null(l1_col) && !is.null(l2_col)) {
    pve1 <- dt[[l1_col]] / (dt[[l1_col]] + dt[[l2_col]])
    panels[["A_PVE1"]] <- data.table(pos_mb = pos_mb, value = pve1, panel = "A: PVE1")
  }

  # Panel B: z-score
  if ("max_abs_z" %in% names(dt)) {
    panels[["B_zscore"]] <- data.table(pos_mb = pos_mb, value = dt$max_abs_z,
                                        panel = "B: |z-score|")
  }

  # Panel C: inv-likeness
  if ("inv_likeness" %in% names(dt)) {
    panels[["C_inv"]] <- data.table(pos_mb = pos_mb, value = dt$inv_likeness,
                                     panel = "C: Inv-likeness")
  } else if (nrow(inv_like_dt) > 0) {
    il_chr <- merge(data.table(global_window_id = dt$global_window_id, pos_mb = pos_mb),
                     inv_like_dt[chrom == chr, .(global_window_id, inv_likeness)],
                     by = "global_window_id", all.x = TRUE)
    panels[["C_inv"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$inv_likeness,
                                     panel = "C: Inv-likeness")
  }

  # Panel D: NN distance
  if ("seed_nn_dist" %in% names(dt)) {
    panels[["D_nn"]] <- data.table(pos_mb = pos_mb, value = dt$seed_nn_dist,
                                    panel = "D: NN distance")
  }

  if (length(panels) < 2) return(NULL)

  pdt <- rbindlist(panels)
  pdt <- pdt[is.finite(value)]
  pdt[, panel := factor(panel, levels = c("A: PVE1", "B: |z-score|",
                                           "C: Inv-likeness", "D: NN distance"))]

  # Color by signal strength within each panel
  pdt[, value_q := frank(value, na.last = "keep") / .N, by = panel]

  p <- ggplot(pdt, aes(x = pos_mb, y = value, color = value_q)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_gradientn(colours = c("grey80", "steelblue", "gold", "red3"),
                          name = "Quantile") +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " -- Combined Signal Profile"),
         subtitle = paste0(n, " windows -- red = high signal, grey = low"),
         caption = "Coincident red peaks across panels A+B+C with blue D = strong inversion") +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "grey90", color = NA),
          panel.spacing = unit(0.8, "lines"))
  p
}

# -- 13: k-NN graph structure from sim_mat ----------------------------
# Two metrics per window:
# (a) mean_nn_sim: average similarity to k nearest neighbors (high = in a block)
# (b) block_score: ratio of k-NN sim to background sim (high = distinct block)
build_knn_structure <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  if (is.null(sim_mat) || nrow(sim_mat) < 20) return(NULL)

  n_w <- nrow(sim_mat)
  k <- 5L

  # Background: mean similarity for each window to ALL other windows
  bg_sim <- rowMeans(sim_mat, na.rm = TRUE)

  knn_rows <- list()
  for (wi in seq_len(n_w)) {
    sims <- sim_mat[wi, ]
    sims[wi] <- -Inf
    top_k <- order(sims, decreasing = TRUE)[seq_len(k)]
    nn_sim <- mean(sims[top_k], na.rm = TRUE)
    # Block score: how much better are NN than background?
    block_sc <- if (bg_sim[wi] > 0) nn_sim / bg_sim[wi] else 1
    knn_rows[[wi]] <- data.table(
      pos_mb = (dt$start_bp[wi] + dt$end_bp[wi]) / 2e6,
      mean_nn_sim = nn_sim,
      bg_sim = bg_sim[wi],
      block_score = block_sc
    )
  }
  pdt <- rbindlist(knn_rows)

  # Color by block_score quantile
  bs_q75 <- quantile(pdt$block_score, 0.75, na.rm = TRUE)
  bs_q90 <- quantile(pdt$block_score, 0.90, na.rm = TRUE)
  pdt[, signal := fifelse(block_score >= bs_q90, "strong_block",
                  fifelse(block_score >= bs_q75, "moderate_block", "background"))]

  pal_block <- c("strong_block" = "red3", "moderate_block" = "darkorange", "background" = "grey70")

  p <- ggplot(pdt, aes(x = pos_mb)) +
    geom_point(aes(y = mean_nn_sim, color = signal), size = 0.6, alpha = 0.6) +
    geom_line(aes(y = bg_sim), color = "steelblue", linewidth = 0.3, alpha = 0.4) +
    scale_color_manual(values = pal_block, name = "Block signal") +
    labs(x = paste0(chr, " (Mb)"),
         y = "Mean k-NN similarity",
         title = paste0(chr, " -- k-NN Block Structure (k=5)"),
         subtitle = paste0(n_w, " windows | Blue line = background (mean sim to all)\n",
                          "Points = mean sim to 5 nearest neighbors"),
         caption = paste0("Red = NN sim >> background (distinct block = inversion candidate)\n",
                         "Orange = moderate block | Grey = background level\n",
                         "q75=", round(bs_q75, 2), " q90=", round(bs_q90, 2))) +
    THEME_BASE
  p
}

# -- 14: Local similarity heatmap (band around diagonal) --------------
build_local_sim_heatmap <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  if (is.null(sim_mat) || nrow(sim_mat) < 20) return(NULL)

  n_w <- nrow(sim_mat)
  bandwidth <- 80L

  # Denser subsampling (max 1000 points along diagonal)
  step <- max(1L, n_w %/% 1000L)
  idx <- seq(1, n_w, by = step)

  rows <- list()
  for (i in idx) {
    j_min <- max(1L, i - bandwidth)
    j_max <- min(n_w, i + bandwidth)
    j_range <- seq(j_min, j_max, by = max(1L, step %/% 2L))
    for (j in j_range) {
      s <- sim_mat[i, j]
      if (is.finite(s)) {
        rows[[length(rows) + 1]] <- data.table(
          x_mb = (dt$start_bp[i] + dt$end_bp[i]) / 2e6,
          y_mb = (dt$start_bp[j] + dt$end_bp[j]) / 2e6,
          sim = s
        )
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  pdt <- rbindlist(rows)

  # Center color scale on median similarity for better contrast
  med_sim <- median(pdt$sim, na.rm = TRUE)

  p <- ggplot(pdt, aes(x = x_mb, y = y_mb, fill = sim)) +
    geom_tile(width = step * (dt$end_bp[2] - dt$start_bp[1]) / 1e6,
              height = step * (dt$end_bp[2] - dt$start_bp[1]) / 1e6) +
    scale_fill_gradientn(
      colours = c("navy", "steelblue", "white", "orange", "red3"),
      values = scales::rescale(c(0, med_sim * 0.7, med_sim, med_sim * 1.15, 1)),
      name = "Similarity"
    ) +
    coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " -- Local Similarity (+-", bandwidth, " windows)"),
         subtitle = paste0("Snake extension context | median sim = ", round(med_sim, 3)),
         caption = "White = median similarity | Red = high (block) | Blue = low (gap)") +
    THEME_BASE
  p
}

# -- 15: Composite per-chr overview (all key plots on one page) -------
build_composite_page <- function(chr) {
  # Collect individual plots
  p_eigen <- tryCatch(build_eigen_profile(chr), error = function(e) NULL)
  p_inv   <- tryCatch(build_inv_likeness_profile(chr), error = function(e) NULL)
  p_z     <- tryCatch(build_z_profile(chr), error = function(e) NULL)
  p_nn    <- tryCatch(build_nn_profile(chr), error = function(e) NULL)
  p_mds   <- tryCatch(build_mds_raw(chr), error = function(e) NULL)
  p_knn   <- tryCatch(build_knn_structure(chr), error = function(e) NULL)

  # Filter NULLs
  plots <- Filter(Negate(is.null), list(p_eigen, p_inv, p_z, p_nn, p_mds, p_knn))
  if (length(plots) < 3) return(NULL)

  # Use patchwork if available, otherwise gridExtra
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    # Layout: top row = 3 wide panels, bottom row = 3 panels
    n_plots <- length(plots)
    if (n_plots >= 6) {
      composite <- (plots[[1]] | plots[[2]] | plots[[3]]) /
                   (plots[[4]] | plots[[5]] | plots[[6]]) +
        plot_annotation(
          title = paste0(chr, " -- Complete Diagnostic Overview"),
          subtitle = paste0(precomp_list[[chr]]$n_windows, " windows"),
          theme = theme(plot.title = element_text(size = 14, face = "bold"))
        )
    } else {
      composite <- wrap_plots(plots, ncol = min(3, n_plots)) +
        plot_annotation(
          title = paste0(chr, " -- Complete Diagnostic Overview"),
          theme = theme(plot.title = element_text(size = 14, face = "bold"))
        )
    }
    return(composite)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    n_plots <- length(plots)
    ncol_g <- min(3L, n_plots)
    nrow_g <- ceiling(n_plots / ncol_g)
    composite <- gridExtra::arrangeGrob(
      grobs = plots, ncol = ncol_g, nrow = nrow_g,
      top = grid::textGrob(paste0(chr, " -- Complete Diagnostic Overview"),
                           gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
    return(composite)
  } else {
    message("  [SKIP] composite: install patchwork or gridExtra")
    return(NULL)
  }
}

# =============================================================================
# RENDER PER-CHROMOSOME PLOTS
# =============================================================================

# -- 16: Contrast-enhanced sim_mat (row-wise Z-score) ----------------
# Each row of the sim_mat is Z-scored: subtract row mean, divide by row SD.
# This removes the global baseline inflation from family structure. Windows
# that are "more similar than expected" relative to their own row stand out.
# NOT a replacement for raw sim_mat — an additional diagnostic view.
build_simmat_zscore <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 50) return(NULL)

  # Row-wise Z-score
  z_mat <- sim_mat
  for (ri in seq_len(n)) {
    row_vals <- sim_mat[ri, ]
    rm <- mean(row_vals, na.rm = TRUE)
    rsd <- sd(row_vals, na.rm = TRUE)
    if (is.finite(rsd) && rsd > 1e-8) {
      z_mat[ri, ] <- (row_vals - rm) / rsd
    } else {
      z_mat[ri, ] <- 0
    }
  }

  # Symmetrise: average with transpose for cleaner visual
  z_mat <- (z_mat + t(z_mat)) / 2

  # Subsample (max 500x500)
  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  sub <- z_mat[idx, idx]
  pos_mb <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6

  pdt <- data.table(
    x = rep(pos_mb, each = length(idx)),
    y = rep(pos_mb, times = length(idx)),
    z = as.vector(sub)
  )
  pdt <- pdt[is.finite(z)]
  # Clip extremes for display
  pdt[, z_clipped := pmin(pmax(z, -3), 3)]

  p <- ggplot(pdt, aes(x = x, y = y, fill = z_clipped)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "grey95", high = "red3",
                         midpoint = 0, limits = c(-3, 3),
                         name = "Row Z-score") +
    coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " -- Contrast Sim_Mat (Row-wise Z-score)"),
         subtitle = paste0(n, " windows (", length(idx), "x", length(idx),
                          " subsampled) | Symmetrised"),
         caption = paste0("Row Z removes family baseline inflation\n",
                         "Red = higher-than-expected similarity (inversion block)\n",
                         "Blue = lower-than-expected (gap/boundary)\n",
                         "Clipped at +-3 SD")) +
    THEME_BASE
  p
}

# -- 17: Contrast-enhanced sim_mat (row quantile-normalized) ----------
# Per-row rank transform: replace each value with its within-row quantile.
# Produces a [0,1] matrix where blocks stand out as regions where nearby
# windows have unusually high rank relative to far-away windows.
build_simmat_quantile <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  n <- nrow(sim_mat)
  if (n < 50) return(NULL)

  # Row-wise quantile normalization
  q_mat <- sim_mat
  for (ri in seq_len(n)) {
    row_vals <- sim_mat[ri, ]
    valid <- is.finite(row_vals)
    if (sum(valid) < 10) { q_mat[ri, ] <- 0.5; next }
    q_mat[ri, valid] <- frank(row_vals[valid]) / sum(valid)
    q_mat[ri, !valid] <- 0.5
  }
  q_mat <- (q_mat + t(q_mat)) / 2

  step <- max(1L, n %/% 500L)
  idx <- seq(1, n, by = step)
  sub <- q_mat[idx, idx]
  pos_mb <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6

  pdt <- data.table(
    x = rep(pos_mb, each = length(idx)),
    y = rep(pos_mb, times = length(idx)),
    q = as.vector(sub)
  )
  pdt <- pdt[is.finite(q)]

  p <- ggplot(pdt, aes(x = x, y = y, fill = q)) +
    geom_tile() +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red3",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Row quantile") +
    coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " -- Contrast Sim_Mat (Row Quantile-Normalized)"),
         subtitle = paste0(n, " windows (", length(idx), "x", length(idx),
                          " subsampled) | Symmetrised"),
         caption = paste0("Row-quantile removes per-row baseline\n",
                         "Red = similarity in top quantile for that row\n",
                         "Blue = similarity in bottom quantile\n",
                         "Blocks visible where nearby windows outrank distant")) +
    THEME_BASE
  p
}

# -- 18: MDS colored by het_contrast (v8.5 dimension) -----------------
build_mds_het <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  het_col <- if ("inv_het_contrast" %in% names(dt)) "inv_het_contrast"
             else if ("het_contrast" %in% names(dt)) "het_contrast"
             else NULL
  if (is.null(het_col)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     het = dt[[het_col]])
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(het)]
  if (nrow(pdt) < 10) return(NULL)

  het_med <- median(pdt$het, na.rm = TRUE)
  p <- ggplot(pdt[order(het)], aes(x = MDS1, y = MDS2, color = het)) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_gradientn(
      colours = c("grey85", "lightblue", "steelblue", "darkorange", "red3"),
      values = scales::rescale(c(0, het_med * 0.5, het_med, het_med * 1.5,
                                  max(pdt$het, na.rm = TRUE))),
      name = "Het contrast"
    ) +
    labs(title = paste0(chr, " -- MDS colored by Het-Contrast"),
         subtitle = paste0(nrow(pdt), " windows | median=", round(het_med, 3)),
         caption = "High het_contrast = clear between/within variance separation (k=3)\nRed = strong inversion signal, grey = background") +
    THEME_BASE
  p
}

# -- 19: MDS colored by lambda2 (low = clean single-axis inversion) ---
build_mds_lambda2 <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2"
            else if ("lambda2" %in% names(dt)) "lambda2"
            else NULL
  if (is.null(l2_col)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     l2 = dt[[l2_col]])
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(l2) & l2 > 0]
  if (nrow(pdt) < 10) return(NULL)

  l2_med <- median(pdt$l2, na.rm = TRUE)
  # Invert: low lambda2 = good inversion, so color low=red, high=grey
  p <- ggplot(pdt[order(-l2)], aes(x = MDS1, y = MDS2, color = l2)) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_gradientn(
      colours = c("red3", "darkorange", "gold", "grey80"),
      name = "Lambda2"
    ) +
    labs(title = paste0(chr, " -- MDS colored by Lambda2"),
         subtitle = paste0(nrow(pdt), " windows | median=", round(l2_med, 3)),
         caption = "Low lambda2 (red) = clean single-axis inversion\nHigh lambda2 (grey) = complex multi-axis or background") +
    THEME_BASE
  p
}

# -- 20: MDS colored by PVE1 excess over chr median (v8.5) ------------
build_mds_pve1_excess <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  if ("inv_pve1_excess" %in% names(dt)) {
    excess_col <- "inv_pve1_excess"
  } else {
    # Compute on the fly
    l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
    l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
    if (is.null(l1_col) || is.null(l2_col)) return(NULL)
    dt[, .pve1_tmp := get(l1_col) / (get(l1_col) + get(l2_col))]
    chr_med <- median(dt$.pve1_tmp, na.rm = TRUE)
    dt[, .pve1_excess := .pve1_tmp - chr_med]
    excess_col <- ".pve1_excess"
  }

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     excess = dt[[excess_col]])
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(excess)]
  if (nrow(pdt) < 10) return(NULL)

  p <- ggplot(pdt[order(excess)], aes(x = MDS1, y = MDS2, color = excess)) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_gradient2(low = "navy", mid = "grey80", high = "red3",
                          midpoint = 0, name = "PVE1 excess") +
    labs(title = paste0(chr, " -- MDS colored by PVE1 Excess"),
         subtitle = paste0(nrow(pdt), " windows | excess = PVE1 - chr_median"),
         caption = "Red = PVE1 above chromosome median (PC1-dominant)\nBlue = below median") +
    THEME_BASE
  p
}

# -- 21: MDS colored by NN distance -----------------------------------
build_mds_nn <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  if (!"seed_nn_dist" %in% names(dt)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     nn = dt$seed_nn_dist)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(nn)]
  if (nrow(pdt) < 10) return(NULL)

  nn_q25 <- quantile(pdt$nn, 0.25, na.rm = TRUE)
  # Low NN = tightly clustered = inversion signal → red
  p <- ggplot(pdt[order(-nn)], aes(x = MDS1, y = MDS2, color = nn)) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_gradientn(
      colours = c("red3", "darkorange", "gold", "grey80"),
      name = "NN distance"
    ) +
    labs(title = paste0(chr, " -- MDS colored by NN Distance"),
         subtitle = paste0(nrow(pdt), " windows | q25=", round(nn_q25, 3)),
         caption = "Low NN (red) = tightly clustered in MDS = inversion core\nHigh NN (grey) = isolated/background") +
    THEME_BASE
  p
}

# -- 22: MDS colored by eigenvalue ratio (lambda1/lambda2) ------------
build_mds_eigratio <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (is.null(l1_col) || is.null(l2_col)) return(NULL)

  ratio <- dt[[l1_col]] / dt[[l2_col]]
  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     ratio = ratio)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(ratio)]
  if (nrow(pdt) < 10) return(NULL)

  # Cap for display
  pdt[, ratio_cap := pmin(ratio, 15)]

  p <- ggplot(pdt[order(ratio_cap)], aes(x = MDS1, y = MDS2, color = ratio_cap)) +
    geom_point(size = 1.0, alpha = 0.6) +
    scale_color_gradientn(
      colours = c("grey85", "lightblue", "steelblue", "navy"),
      name = "L1/L2 (cap 15)"
    ) +
    labs(title = paste0(chr, " -- MDS colored by Eigenvalue Ratio"),
         subtitle = paste0(nrow(pdt), " windows"),
         caption = "High ratio (dark) = strong PC1 dominance\nInversions typically L1/L2 > 3") +
    THEME_BASE
  p
}

# -- 23: Faceted MDS — all colorings on one page ----------------------
# 6 small panels: position, inv-likeness, z-score, het_contrast, lambda2, NN
build_mds_faceted_all <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)

  n <- nrow(dt)
  pos_mb <- (dt$start_bp + dt$end_bp) / 2e6

  # Gather available metrics, rescale all to [0,1] for unified color
  rescale01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (r[2] == r[1]) return(rep(0.5, length(x)))
    (x - r[1]) / (r[2] - r[1])
  }

  panels <- list()

  # Position
  panels[["A: Position"]] <- data.table(
    MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
    value = rescale01(pos_mb), panel = "A: Position"
  )

  # Inv-likeness
  inv_vals <- if ("inv_likeness" %in% names(dt)) dt$inv_likeness
              else if (nrow(inv_like_dt) > 0 && "global_window_id" %in% names(dt)) {
                m <- merge(data.table(global_window_id = dt$global_window_id),
                           inv_like_dt[chrom == chr, .(global_window_id, inv_likeness)],
                           by = "global_window_id", all.x = TRUE)
                m$inv_likeness
              } else NULL
  if (!is.null(inv_vals)) {
    panels[["B: Inv-likeness"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = rescale01(inv_vals), panel = "B: Inv-likeness"
    )
  }

  # Z-score
  if ("max_abs_z" %in% names(dt)) {
    panels[["C: |z-score|"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = rescale01(pmin(abs(dt$max_abs_z), 5)), panel = "C: |z-score|"
    )
  }

  # Het contrast
  het_col <- if ("inv_het_contrast" %in% names(dt)) "inv_het_contrast"
             else if ("het_contrast" %in% names(dt)) "het_contrast"
             else NULL
  if (!is.null(het_col)) {
    panels[["D: Het-contrast"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = rescale01(dt[[het_col]]), panel = "D: Het-contrast"
    )
  }

  # Lambda2 (inverted: low = high signal)
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2"
            else if ("lambda2" %in% names(dt)) "lambda2"
            else NULL
  if (!is.null(l2_col)) {
    l2_vals <- dt[[l2_col]]
    l2_valid <- is.finite(l2_vals) & l2_vals > 0
    l2_rescaled <- rep(0.5, n)
    l2_rescaled[l2_valid] <- 1 - rescale01(l2_vals[l2_valid])  # invert
    panels[["E: Low-lambda2"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = l2_rescaled, panel = "E: Low-lambda2"
    )
  }

  # NN distance (inverted: low = high signal)
  if ("seed_nn_dist" %in% names(dt)) {
    nn_vals <- dt$seed_nn_dist
    nn_valid <- is.finite(nn_vals)
    nn_rescaled <- rep(0.5, n)
    nn_rescaled[nn_valid] <- 1 - rescale01(nn_vals[nn_valid])  # invert
    panels[["F: Low-NN"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = nn_rescaled, panel = "F: Low-NN"
    )
  }

  if (length(panels) < 3) return(NULL)

  pdt <- rbindlist(panels)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(value)]
  # Panels ordered by label
  pdt[, panel := factor(panel, levels = sort(unique(panel)))]

  n_panels <- length(unique(pdt$panel))
  ncol_f <- min(3L, n_panels)

  p <- ggplot(pdt, aes(x = MDS1, y = MDS2, color = value)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_gradientn(
      colours = c("grey85", "lightblue", "steelblue", "darkorange", "red3"),
      name = "Signal (0-1)"
    ) +
    facet_wrap(~ panel, ncol = ncol_f, scales = "free") +
    labs(title = paste0(chr, " -- MDS Faceted by All Metrics"),
         subtitle = paste0(nrow(dt), " windows | ", n_panels, " panels | all rescaled [0,1]"),
         caption = paste0("Red = high signal for each metric\n",
                         "E/F are inverted (low lambda2 / low NN = high signal)\n",
                         "Coincident red across panels = strongest inversion evidence")) +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold", size = 8),
          strip.background = element_rect(fill = "grey92", color = NA),
          panel.spacing = unit(0.5, "lines"),
          legend.key.width = unit(0.8, "cm"))
  p
}

# -- 24: Local entropy profile along chromosome ----------------------
build_entropy_profile <- function(chr) {
  il_chr <- if (nrow(inv_like_dt) > 0) inv_like_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0 || !"local_entropy" %in% names(il_chr)) return(NULL)
  il_chr <- il_chr[is.finite(local_entropy) & is.finite(start_bp)]
  if (nrow(il_chr) < 20) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  p <- ggplot(il_chr, aes(x = pos_mb)) +
    geom_point(aes(y = local_entropy), size = 0.4, alpha = 0.5, color = "darkorange") +
    geom_hline(yintercept = log(3), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "Local Shannon entropy (H)",
         title = paste0(chr, " -- Local Structure Entropy Profile"),
         subtitle = paste0(nrow(il_chr), " windows | chr-wide fixed k=3 bands, per-window soft membership"),
         caption = paste0("Source: per-window PC1 loadings measured against chromosome-wide fixed k=3 bands\n",
                         "Soft membership via inverse distance to fixed band means\n",
                         "H = -\u03A3 P_i \u00B7 ln(P_i) averaged across 226 samples\n",
                         "Low H = clean single group dominates | High H = mixed | Dashed = ln(3) = max for k=3")) +
    THEME_BASE
  p
}

# -- 25: ENA profile along chromosome --------------------------------
build_ena_profile <- function(chr) {
  il_chr <- if (nrow(inv_like_dt) > 0) inv_like_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0 || !"local_ena" %in% names(il_chr)) return(NULL)
  il_chr <- il_chr[is.finite(local_ena) & is.finite(start_bp)]
  if (nrow(il_chr) < 20) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  p <- ggplot(il_chr, aes(x = pos_mb, y = local_ena)) +
    geom_point(size = 0.4, alpha = 0.5, color = "red3") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "green4", linewidth = 0.3) +
    geom_hline(yintercept = 3, linetype = "dotted", color = "grey50", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "Effective Number of Backgrounds (ENA)",
         title = paste0(chr, " -- ENA Profile"),
         subtitle = paste0(nrow(il_chr), " windows | ENA = exp(H), chr-wide fixed bands"),
         caption = paste0("Source: exp(local Shannon entropy) from chr-wide fixed k=3 band memberships\n",
                         "ENA \u2248 1 = single clean background (green dashed)\n",
                         "ENA \u2248 3 = max mixing for k=3 (grey dotted)\n",
                         "No smoothing applied: raw per-window values")) +
    THEME_BASE
  p
}

# -- 26: Multi-source het assessment (6-panel combined) ---------------
build_het_assessment <- function(chr) {
  il_chr <- if (nrow(inv_like_dt) > 0) inv_like_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0) return(NULL)
  if (!any(c("pc1_bimodality", "het_pc1_gap", "het_mid_fraction",
             "local_delta12", "local_entropy") %in% names(il_chr))) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  panels <- list()

  if ("inv_het_contrast" %in% names(il_chr) && any(is.finite(il_chr$inv_het_contrast)))
    panels[["A: Het-contrast\n(between/within var)"]] <- data.table(
      pos_mb = il_chr$pos_mb, value = il_chr$inv_het_contrast,
      panel = "A: Het-contrast\n(between/within var)")

  if ("pc1_bimodality" %in% names(il_chr) && any(is.finite(il_chr$pc1_bimodality)))
    panels[["B: PC1 bimodality\n(Ashman D)"]] <- data.table(
      pos_mb = il_chr$pos_mb, value = il_chr$pc1_bimodality,
      panel = "B: PC1 bimodality\n(Ashman D)")

  if ("het_pc1_gap" %in% names(il_chr) && any(is.finite(il_chr$het_pc1_gap)))
    panels[["C: PC1 gap\n(min center dist)"]] <- data.table(
      pos_mb = il_chr$pos_mb, value = il_chr$het_pc1_gap,
      panel = "C: PC1 gap\n(min center dist)")

  if ("het_mid_fraction" %in% names(il_chr) && any(is.finite(il_chr$het_mid_fraction)))
    panels[["D: Middle band fraction"]] <- data.table(
      pos_mb = il_chr$pos_mb, value = il_chr$het_mid_fraction,
      panel = "D: Middle band fraction")

  if ("local_delta12" %in% names(il_chr) && any(is.finite(il_chr$local_delta12)))
    panels[["E: Local \u039412\n(dominance)"]] <- data.table(
      pos_mb = il_chr$pos_mb, value = il_chr$local_delta12,
      panel = "E: Local \u039412\n(dominance)")

  if ("local_entropy" %in% names(il_chr) && any(is.finite(il_chr$local_entropy)))
    panels[["F: Local entropy (H)"]] <- data.table(
      pos_mb = il_chr$pos_mb, value = il_chr$local_entropy,
      panel = "F: Local entropy (H)")

  if (length(panels) < 3) return(NULL)
  pdt <- rbindlist(panels)
  pdt <- pdt[is.finite(value)]
  pdt[, panel := factor(panel, levels = sort(unique(panel)))]

  p <- ggplot(pdt, aes(x = pos_mb, y = value)) +
    geom_point(size = 0.3, alpha = 0.4, color = "steelblue") +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " -- Multi-Source Het Assessment"),
         subtitle = paste0(nrow(il_chr), " windows | chr-wide fixed k=3 bands, 6 metrics"),
         caption = paste0("A: k=3 between/within variance ratio (het_contrast, per-window k-means)\n",
                         "B: Ashman D between outer fixed bands (bimodality of hom groups)\n",
                         "C: minimum gap between fixed band means at this window (group separation)\n",
                         "D: fraction of samples in middle fixed band (expected ~50% for common inv)\n",
                         "E: mean \u039412 = max(P) - second(P) from soft fixed-band memberships\n",
                         "F: mean Shannon entropy from soft fixed-band memberships\n",
                         "No smoothing: raw values | Coincident peaks A-E with trough F = inversion")) +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold", size = 7),
          panel.spacing = unit(0.6, "lines"))
  p
}

plot_builders <- list(
  "01_eigenvalue_spectrum" = build_eigen_profile,
  "02_inv_likeness_profile" = build_inv_likeness_profile,
  "03_max_abs_z_profile" = build_z_profile,
  "04_bg_continuity_baseline" = build_bg_baseline,
  "05_nn_distance_profile" = build_nn_profile,
  "06_sim_mat_heatmap" = build_sim_heatmap,
  "07_mds_scatter_raw" = build_mds_raw,
  "07b_mds_vs_position" = build_mds_vs_position,
  "08_nn_histogram_adaptive" = build_nn_histogram,
  "09_mds_by_inv_likeness" = build_mds_inv,
  "10_mds_by_zscore" = build_mds_z,
  "11_nn_calibration" = build_nn_calibration,
  "12_combined_profile" = build_combined_profile,
  "13_knn_graph_structure" = build_knn_structure,
  "14_local_sim_heatmap" = build_local_sim_heatmap,
  "16_simmat_zscore_contrast" = build_simmat_zscore,
  "17_simmat_quantile_contrast" = build_simmat_quantile,
  "18_mds_by_het_contrast" = build_mds_het,
  "19_mds_by_lambda2" = build_mds_lambda2,
  "20_mds_by_pve1_excess" = build_mds_pve1_excess,
  "21_mds_by_nn_distance" = build_mds_nn,
  "22_mds_by_eigenvalue_ratio" = build_mds_eigratio,
  "23_mds_faceted_all_metrics" = build_mds_faceted_all,
  "24_local_entropy_profile" = build_entropy_profile,
  "25_ena_profile" = build_ena_profile,
  "26_het_assessment_multi" = build_het_assessment
)

# Page dimensions per plot type
plot_dims <- list(
  "01_eigenvalue_spectrum" = c(w = 14, h = 5),
  "02_inv_likeness_profile" = c(w = 14, h = 5),
  "03_max_abs_z_profile" = c(w = 14, h = 5),
  "04_bg_continuity_baseline" = c(w = 8, h = 5),
  "05_nn_distance_profile" = c(w = 14, h = 5),
  "06_sim_mat_heatmap" = c(w = 8, h = 7),
  "07_mds_scatter_raw" = c(w = 8, h = 7),
  "07b_mds_vs_position" = c(w = 18, h = 6),
  "08_nn_histogram_adaptive" = c(w = 8, h = 5),
  "09_mds_by_inv_likeness" = c(w = 8, h = 7),
  "10_mds_by_zscore" = c(w = 8, h = 7),
  "11_nn_calibration" = c(w = 14, h = 5),
  "12_combined_profile" = c(w = 16, h = 16),
  "13_knn_graph_structure" = c(w = 14, h = 5),
  "14_local_sim_heatmap" = c(w = 8, h = 7),
  "16_simmat_zscore_contrast" = c(w = 8, h = 7),
  "17_simmat_quantile_contrast" = c(w = 8, h = 7),
  "18_mds_by_het_contrast" = c(w = 8, h = 7),
  "19_mds_by_lambda2" = c(w = 8, h = 7),
  "20_mds_by_pve1_excess" = c(w = 8, h = 7),
  "21_mds_by_nn_distance" = c(w = 8, h = 7),
  "22_mds_by_eigenvalue_ratio" = c(w = 8, h = 7),
  "23_mds_faceted_all_metrics" = c(w = 16, h = 12),
  "24_local_entropy_profile" = c(w = 14, h = 5),
  "25_ena_profile" = c(w = 14, h = 5),
  "26_het_assessment_multi" = c(w = 16, h = 14)
)

for (pname in names(plot_builders)) {
  message("[C01a_diag] Building: ", pname, " ...")
  builder <- plot_builders[[pname]]
  dims <- plot_dims[[pname]]

  plots <- list()
  for (chr in chroms) {
    p <- tryCatch(builder(chr), error = function(e) {
      message("  [WARN] BUILD ", pname, " / ", chr, ": ", e$message); NULL
    })
    if (!is.null(p)) {
      message("  [OK] ", pname, " / ", chr)
      plots[[chr]] <- p
    }
  }

  if (length(plots) > 0) {
    # PDF output (use cairo_pdf for robust rendering on headless HPC)
    f_pdf <- file.path(outdir, paste0(pname, ".pdf"))
    message("[C01a_diag] Rendering PDF: ", f_pdf, " (", length(plots), " pages)")
    use_cairo <- capabilities("cairo")
    if (use_cairo) {
      cairo_pdf(f_pdf, width = dims["w"], height = dims["h"], onefile = TRUE)
    } else {
      pdf(f_pdf, width = dims["w"], height = dims["h"])
    }
    for (chr in names(plots)) {
      tryCatch({
        print(plots[[chr]])
      }, error = function(e) {
        message("  [RENDER FAIL] ", pname, " / ", chr, ": ", e$message)
        plot.new()
        text(0.5, 0.5, paste0("RENDER ERROR: ", chr, "\n", e$message), cex = 0.8)
      })
    }
    dev.off()

    # PNG output
    png_dir <- file.path(outdir, "png", pname)
    dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
    for (chr in names(plots)) {
      f_png <- file.path(png_dir, paste0(chr, ".png"))
      tryCatch(
        ggsave(f_png, plots[[chr]], width = dims["w"], height = dims["h"],
               dpi = DPI, device = if (use_cairo) "png" else "png"),
        error = function(e) message("  [PNG FAIL] ", chr, ": ", e$message)
      )
    }
    message("[C01a_diag] ", pname, ": DONE")
  } else {
    message("[C01a_diag] ", pname, ": no plots generated (all chr returned NULL)")
  }
}

# =============================================================================
# COMPOSITE PER-CHROMOSOME PAGES (one page = all plots for one chr)
# =============================================================================

message("[C01a_diag] Building composite overview pages...")
composite_plots <- list()
for (chr in chroms) {
  cp <- tryCatch(build_composite_page(chr), error = function(e) {
    message("  [WARN] composite / ", chr, ": ", e$message); NULL
  })
  if (!is.null(cp)) composite_plots[[chr]] <- cp
}

if (length(composite_plots) > 0) {
  f_comp <- file.path(outdir, "15_composite_overview.pdf")
  use_cairo <- capabilities("cairo")
  if (use_cairo) {
    cairo_pdf(f_comp, width = 20, height = 14, onefile = TRUE)
  } else {
    pdf(f_comp, width = 20, height = 14)
  }
  for (chr in names(composite_plots)) {
    tryCatch({
      if (inherits(composite_plots[[chr]], "patchwork")) {
        print(composite_plots[[chr]])
      } else {
        grid::grid.draw(composite_plots[[chr]])
        grid::grid.newpage()
      }
    }, error = function(e) {
      message("  [COMPOSITE FAIL] ", chr, ": ", e$message)
      plot.new(); text(0.5, 0.5, paste0("ERROR: ", chr))
    })
  }
  dev.off()
  message("[C01a_diag] Composite: ", length(composite_plots), " pages -> ", f_comp)
}

# =============================================================================
# GENOME-WIDE SUMMARY PLOTS (PNG)
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
