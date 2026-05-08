#!/usr/bin/env Rscript
# =============================================================================
# C01a_diag_mds.R — MDS scatter plots and sample-space diagnostics
# (plots 07-11, 13, 18-23, 26)
#
# v8.5.2 — Plot overhaul:
#   G1: sys.frame → env var sourcing
#   PC1: MDS faceted: larger panels, viridis, better strip labels
#   PC5: MDS vs position: taller panels, smoothed trend
#   PC7: k-NN: clarify blue line, use inv_likeness for classification
#   PC9: NN calibration: remove fixed 0.70 red line
#   PC10: NN histogram: remove fixed 0.70, composite distribution
#
# Usage: Rscript C01a_diag_mds.R <precomp_dir> <outdir> [chrom]
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
dir.create(cli$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# MDS BUILDERS
# =============================================================================

# -- 07: MDS scatter (raw, viridis) ------------------------------------------
build_mds_raw <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     pos_mb = (dt$start_bp + dt$end_bp) / 2e6)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2)]

  p <- ggplot(pdt, aes(x = MDS1, y = MDS2, color = pos_mb)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_viridis_c(name = "Position (Mb)", option = "C") +
    labs(title = paste0(chr, " \u2014 MDS1 vs MDS2"),
         subtitle = paste0(nrow(pdt), " windows, colored by genomic position"),
         caption = "Clusters = candidate inversion regions") +
    THEME_BASE
  p
}

# -- 07b: Lostruct-style MDS vs position (PC5: taller panels, smoothed) ------
build_mds_vs_position <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)

  pos_mb <- (dt$start_bp + dt$end_bp) / 2e6
  pdt <- data.table(
    pos_mb = rep(pos_mb, 3),
    value = c(mds_mat[, 1], mds_mat[, 2], pos_mb),
    panel = rep(c("MDS1 vs position", "MDS2 vs position", "MDS1 vs MDS2"), each = nrow(dt)),
    x_val = c(pos_mb, pos_mb, mds_mat[, 1]),
    color_pos = rep(pos_mb, 3)
  )
  pdt[panel == "MDS1 vs MDS2", value := mds_mat[, 2]]
  pdt <- pdt[is.finite(value) & is.finite(x_val)]

  # Themed k=3 palette
  km <- tryCatch(kmeans(mds_mat[, 1:2], centers = 3, nstart = 10), error = function(e) NULL)
  if (!is.null(km)) {
    cluster_col <- rep(km$cluster, 3)
    pdt[, cluster := as.factor(cluster_col[seq_len(.N)])]
    pal_km <- c("1" = "#2563EB", "2" = "#059669", "3" = "#D97706")

    p <- ggplot(pdt, aes(x = x_val, y = value, color = cluster)) +
      geom_point(size = 0.6, alpha = 0.5) +
      # PC5: smoothed trend per cluster
      geom_smooth(method = "loess", span = 0.2, se = FALSE, linewidth = 0.5) +
      scale_color_manual(values = pal_km, guide = "none") +
      facet_wrap(~ panel, ncol = 3, scales = "free") +
      labs(title = paste0(chr, " \u2014 MDS Coordinates vs Position (Lostruct-style)"),
           subtitle = paste0(nrow(dt), " windows | k=3 clustering on MDS1+MDS2"),
           caption = "Smoothed trends show breakpoint transitions\nLeft/center: MDS along chr | Right: MDS scatter") +
      THEME_BASE +
      theme(strip.text = element_text(face = "bold", size = 9))
  } else {
    p <- ggplot(pdt, aes(x = x_val, y = value, color = color_pos)) +
      geom_point(size = 0.6, alpha = 0.5) +
      scale_color_viridis_c(name = "Position (Mb)", option = "C") +
      facet_wrap(~ panel, ncol = 3, scales = "free") +
      labs(title = paste0(chr, " \u2014 MDS Coordinates vs Position"),
           subtitle = paste0(nrow(dt), " windows")) +
      THEME_BASE
  }
  p
}

# -- 08: NN distance histogram (PC10: remove 0.70, themed) ------------------
build_nn_histogram <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  if (!"seed_nn_dist" %in% names(dt)) return(NULL)
  nn <- dt$seed_nn_dist[is.finite(dt$seed_nn_dist)]
  if (length(nn) < 20) return(NULL)

  qq <- quantile(nn, c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
  x_cap <- quantile(nn, 0.99, na.rm = TRUE) * 1.1
  pdt <- data.table(nn = nn[nn <= x_cap])

  p <- ggplot(pdt, aes(x = nn)) +
    geom_histogram(bins = 80, fill = "#DBEAFE", color = "#3B82F6", linewidth = 0.2) +
    geom_vline(xintercept = qq["25%"], color = "#2563EB", linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = qq["50%"], color = "#059669", linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = qq["75%"], color = "#D97706", linetype = "dashed", linewidth = 0.5) +
    # PC10: NO fixed 0.70 red line
    annotate("text", x = qq["25%"], y = Inf, vjust = 1.5,
             label = paste0("q25=", round(qq["25%"], 3)), size = 2.5, color = "#2563EB") +
    annotate("text", x = qq["50%"], y = Inf, vjust = 3,
             label = paste0("q50=", round(qq["50%"], 3)), size = 2.5, color = "#059669") +
    annotate("text", x = qq["75%"], y = Inf, vjust = 4.5,
             label = paste0("q75=", round(qq["75%"], 3)), size = 2.5, color = "#D97706") +
    labs(x = "NN Distance (MDS-space)", y = "Count",
         title = paste0(chr, " \u2014 NN Distance Distribution"),
         subtitle = paste0(length(nn), " windows | trimmed at q99\n",
                          "Adaptive seed thresholds: S1S<q25 S1M<q50 S1L<q75"),
         caption = paste0("q10=", round(qq["10%"], 3),
                         " q25=", round(qq["25%"], 3),
                         " q50=", round(qq["50%"], 3),
                         " q75=", round(qq["75%"], 3),
                         " q90=", round(qq["90%"], 3))) +
    THEME_BASE
  p
}

# -- 09: MDS colored by inv_likeness -----------------------------------------
build_mds_inv <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  if (!"inv_likeness" %in% names(dt)) {
    if (nrow(inv_like_dt) > 0 && "global_window_id" %in% names(dt))
      dt <- merge(dt, inv_like_dt[chrom == chr, .(global_window_id, inv_likeness)],
                   by = "global_window_id", all.x = TRUE)
    else return(NULL)
  }

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2], inv = dt$inv_likeness)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(inv)]
  if (nrow(pdt) < 10) return(NULL)
  inv_med <- median(pdt$inv, na.rm = TRUE)

  p <- ggplot(pdt[order(inv)], aes(x = MDS1, y = MDS2, color = inv)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_gradientn(
      colours = c("#E5E7EB", "#93C5FD", "#3B82F6", "#D97706", "#DC2626"),
      values = scales::rescale(c(0, inv_med * 0.7, inv_med, inv_med * 1.2, 1)),
      name = "Inv-likeness") +
    labs(title = paste0(chr, " \u2014 MDS colored by Inv-Likeness"),
         subtitle = paste0(nrow(pdt), " windows"),
         caption = "Red = high inv-likeness | Grey = background") +
    THEME_BASE
  p
}

# -- 10: MDS colored by z-score -----------------------------------------------
build_mds_z <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  if (!"max_abs_z" %in% names(dt)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2], z = dt$max_abs_z)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(z)]
  if (nrow(pdt) < 10) return(NULL)

  p <- ggplot(pdt[order(z)], aes(x = MDS1, y = MDS2, color = pmin(z, 5))) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_gradient2(low = "#E5E7EB", mid = "#3B82F6", high = "#0C1E3C",
                          midpoint = 2.5, name = "|z| (cap 5)") +
    labs(title = paste0(chr, " \u2014 MDS colored by Z-Score"),
         subtitle = paste0(nrow(pdt), " windows"),
         caption = "Dark = high z (strong MDS outlier) | Grey = background") +
    THEME_BASE
  p
}

# -- 11: NN calibration (PC9: remove fixed 0.70) -----------------------------
build_nn_calibration <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  if (!"seed_nn_dist" %in% names(dt)) return(NULL)

  pdt <- dt[is.finite(seed_nn_dist), .(pos_mb = (start_bp + end_bp) / 2e6, nn = seed_nn_dist)]
  if (nrow(pdt) < 20) return(NULL)

  qq <- quantile(pdt$nn, c(0.10, 0.25, 0.50, 0.75), na.rm = TRUE)
  y_cap <- quantile(pdt$nn, 0.99, na.rm = TRUE) * 1.1

  pdt[, adaptive_class := fifelse(nn <= qq["25%"], "S1S (q25)",
                          fifelse(nn <= qq["50%"], "S1M (q50)",
                          fifelse(nn <= qq["75%"], "S1L (q75)", "above q75")))]

  p <- ggplot(pdt, aes(x = pos_mb, y = nn, color = adaptive_class)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = PAL_ADAPT) +
    geom_hline(yintercept = qq["25%"], linetype = "dashed", color = "#2563EB", linewidth = 0.3) +
    geom_hline(yintercept = qq["50%"], linetype = "dashed", color = "#059669", linewidth = 0.3) +
    geom_hline(yintercept = qq["75%"], linetype = "dashed", color = "#D97706", linewidth = 0.3) +
    # PC9: fixed 0.70 removed entirely
    labs(x = paste0(chr, " (Mb)"), y = "NN Distance (MDS-space)",
         title = paste0(chr, " \u2014 Seed Threshold Calibration (NN Distance)"),
         subtitle = paste0("Adaptive quantiles: q25=", round(qq["25%"], 3),
                          " q50=", round(qq["50%"], 3),
                          " q75=", round(qq["75%"], 3),
                          " | y-axis capped at q99"),
         caption = paste0("Blue = S1S (<q25, most clustered) | Teal = S1M (<q50) | Amber = S1L (<q75)\n",
                         "Adaptive quantiles calibrate to per-chromosome data range")) +
    coord_cartesian(ylim = c(0, y_cap)) +
    THEME_BASE
  p
}

# -- 13: k-NN block structure (PC7: clarify, use inv_likeness) ---------------
build_knn_structure <- function(chr) {
  sim_mat <- precomp_list[[chr]]$sim_mat
  dt <- precomp_list[[chr]]$dt
  if (is.null(sim_mat) || nrow(sim_mat) < 20) return(NULL)

  n_w <- nrow(sim_mat)
  k <- 5L
  bg_sim <- rowMeans(sim_mat, na.rm = TRUE)

  knn_rows <- list()
  for (wi in seq_len(n_w)) {
    sims <- sim_mat[wi, ]
    sims[wi] <- -Inf
    top_k <- order(sims, decreasing = TRUE)[seq_len(k)]
    nn_sim <- mean(sims[top_k], na.rm = TRUE)
    block_sc <- if (bg_sim[wi] > 0) nn_sim / bg_sim[wi] else 1
    knn_rows[[wi]] <- data.table(
      pos_mb = (dt$start_bp[wi] + dt$end_bp[wi]) / 2e6,
      mean_nn_sim = nn_sim, bg_sim = bg_sim[wi], block_score = block_sc
    )
  }
  pdt <- rbindlist(knn_rows)

  # PC7: Use inv_likeness to separate inversion from family blocks
  has_inv <- "inv_likeness" %in% names(dt)
  if (has_inv) {
    inv_vals <- dt$inv_likeness
    inv_med <- median(inv_vals, na.rm = TRUE)
    bs_q75 <- quantile(pdt$block_score, 0.75, na.rm = TRUE)
    pdt[, signal := fifelse(block_score >= bs_q75 & inv_vals >= inv_med, "inversion_block",
                   fifelse(block_score >= bs_q75, "family_block", "background"))]
    pal_sig <- c("inversion_block" = "#DC2626", "family_block" = "#D97706", "background" = "#D1D5DB")
  } else {
    bs_q75 <- quantile(pdt$block_score, 0.75, na.rm = TRUE)
    bs_q90 <- quantile(pdt$block_score, 0.90, na.rm = TRUE)
    pdt[, signal := fifelse(block_score >= bs_q90, "strong_block",
                   fifelse(block_score >= bs_q75, "moderate_block", "background"))]
    pal_sig <- PAL_BLOCK
  }

  p <- ggplot(pdt, aes(x = pos_mb)) +
    geom_point(aes(y = mean_nn_sim, color = signal), size = 0.6, alpha = 0.6) +
    # PC7: clarify blue line meaning
    geom_line(aes(y = bg_sim), color = "#3B82F6", linewidth = 0.3, alpha = 0.5) +
    scale_color_manual(values = pal_sig, name = "Block signal") +
    labs(x = paste0(chr, " (Mb)"),
         y = "Mean k-NN similarity",
         title = paste0(chr, " \u2014 k-NN Block Structure (k=5)"),
         subtitle = paste0(n_w, " windows | Blue line = background (mean sim to ALL windows)\n",
                          "Points = mean sim to 5 nearest neighbors"),
         caption = paste0(if (has_inv) "Red = inversion block (high NN + high inv_likeness)\nAmber = family block (high NN + low inv_likeness)\n"
                          else "Red = strong block (NN >> background) | Amber = moderate\n",
                         "Grey = background level | q75=", round(bs_q75, 2))) +
    THEME_BASE
  p
}

# -- 18: MDS by het_contrast --------------------------------------------------
build_mds_het <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  het_col <- if ("inv_het_contrast" %in% names(dt)) "inv_het_contrast"
             else if ("het_contrast" %in% names(dt)) "het_contrast"
             else NULL
  if (is.null(het_col)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2], het = dt[[het_col]])
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(het)]
  if (nrow(pdt) < 10) return(NULL)
  het_med <- median(pdt$het, na.rm = TRUE)

  p <- ggplot(pdt[order(het)], aes(x = MDS1, y = MDS2, color = het)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_gradientn(
      colours = c("#E5E7EB", "#93C5FD", "#3B82F6", "#D97706", "#DC2626"),
      values = scales::rescale(c(0, het_med * 0.5, het_med, het_med * 1.5,
                                  max(pdt$het, na.rm = TRUE))),
      name = "Het contrast") +
    labs(title = paste0(chr, " \u2014 MDS by Het-Contrast"),
         subtitle = paste0(nrow(pdt), " windows | median=", round(het_med, 3)),
         caption = "High het_contrast = clear between/within separation") +
    THEME_BASE
  p
}

# -- 19: MDS by lambda2 -------------------------------------------------------
build_mds_lambda2 <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (is.null(l2_col)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2], l2 = dt[[l2_col]])
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(l2) & l2 > 0]
  if (nrow(pdt) < 10) return(NULL)

  p <- ggplot(pdt[order(-l2)], aes(x = MDS1, y = MDS2, color = l2)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_gradientn(colours = c("#DC2626", "#D97706", "#F59E0B", "#D1D5DB"),
                          name = "Lambda2") +
    labs(title = paste0(chr, " \u2014 MDS by Lambda2"),
         subtitle = paste0(nrow(pdt), " windows"),
         caption = "Low \u03BB2 (red) = clean single-axis inversion | High (grey) = complex/background") +
    THEME_BASE
  p
}

# -- 20: MDS by PVE1 excess ---------------------------------------------------
build_mds_pve1_excess <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  if ("inv_pve1_excess" %in% names(dt)) {
    excess_col <- "inv_pve1_excess"
  } else {
    l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
    l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
    if (is.null(l1_col) || is.null(l2_col)) return(NULL)
    dt[, .pve1_tmp := get(l1_col) / (get(l1_col) + get(l2_col))]
    dt[, .pve1_excess := .pve1_tmp - median(.pve1_tmp, na.rm = TRUE)]
    excess_col <- ".pve1_excess"
  }

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2], excess = dt[[excess_col]])
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(excess)]
  if (nrow(pdt) < 10) return(NULL)

  p <- ggplot(pdt[order(excess)], aes(x = MDS1, y = MDS2, color = excess)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_gradient2(low = DIV_LOW, mid = "#E5E7EB", high = DIV_HIGH,
                          midpoint = 0, name = "PVE1 excess") +
    labs(title = paste0(chr, " \u2014 MDS by PVE1 Excess"),
         subtitle = paste0(nrow(pdt), " windows | excess = PVE1 - chr_median"),
         caption = "Red = PVE1 above median (PC1-dominant) | Blue = below") +
    THEME_BASE
  p
}

# -- 21: MDS by NN distance ---------------------------------------------------
build_mds_nn <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  if (!"seed_nn_dist" %in% names(dt)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2], nn = dt$seed_nn_dist)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(nn)]
  if (nrow(pdt) < 10) return(NULL)

  p <- ggplot(pdt[order(-nn)], aes(x = MDS1, y = MDS2, color = nn)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_gradientn(colours = c("#DC2626", "#D97706", "#F59E0B", "#D1D5DB"),
                          name = "NN distance") +
    labs(title = paste0(chr, " \u2014 MDS by NN Distance"),
         subtitle = paste0(nrow(pdt), " windows"),
         caption = "Low NN (red) = tightly clustered = inversion core | High (grey) = background") +
    THEME_BASE
  p
}

# -- 22: MDS by eigenvalue ratio -----------------------------------------------
build_mds_eigratio <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (is.null(l1_col) || is.null(l2_col)) return(NULL)

  pdt <- data.table(MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
                     ratio = dt[[l1_col]] / dt[[l2_col]])
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(ratio)]
  if (nrow(pdt) < 10) return(NULL)
  pdt[, ratio_cap := pmin(ratio, 15)]

  p <- ggplot(pdt[order(ratio_cap)], aes(x = MDS1, y = MDS2, color = ratio_cap)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_color_gradientn(colours = c("#E5E7EB", "#93C5FD", "#3B82F6", "#0C1E3C"),
                          name = "L1/L2 (cap 15)") +
    labs(title = paste0(chr, " \u2014 MDS by Eigenvalue Ratio"),
         subtitle = paste0(nrow(pdt), " windows"),
         caption = "High ratio (dark) = strong PC1 dominance | Inversions typically L1/L2 > 3") +
    THEME_BASE
  p
}

# -- 23: Faceted MDS — all metrics (PC1: larger, viridis, better labels) -----
build_mds_faceted_all <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  mds_mat <- precomp_list[[chr]]$mds_mat
  if (is.null(mds_mat) || ncol(mds_mat) < 2) return(NULL)

  n <- nrow(dt)
  pos_mb <- (dt$start_bp + dt$end_bp) / 2e6

  panels <- list()

  # PC1: improved strip text
  panels[["A: Chromosome position (Mb)"]] <- data.table(
    MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
    value = rescale01(pos_mb), panel = "A: Chromosome position (Mb)"
  )

  inv_vals <- if ("inv_likeness" %in% names(dt)) dt$inv_likeness
              else if (nrow(inv_like_dt) > 0 && "global_window_id" %in% names(dt)) {
                m <- merge(data.table(global_window_id = dt$global_window_id),
                           inv_like_dt[chrom == chr, .(global_window_id, inv_likeness)],
                           by = "global_window_id", all.x = TRUE)
                m$inv_likeness
              } else NULL
  if (!is.null(inv_vals))
    panels[["B: Inv-likeness score"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = rescale01(inv_vals), panel = "B: Inv-likeness score")

  if ("max_abs_z" %in% names(dt))
    panels[["C: |Robust z-score|"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = rescale01(pmin(abs(dt$max_abs_z), 5)), panel = "C: |Robust z-score|")

  het_col <- if ("inv_het_contrast" %in% names(dt)) "inv_het_contrast"
             else if ("het_contrast" %in% names(dt)) "het_contrast"
             else NULL
  if (!is.null(het_col))
    panels[["D: Het-contrast"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = rescale01(dt[[het_col]]), panel = "D: Het-contrast")

  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (!is.null(l2_col)) {
    l2_vals <- dt[[l2_col]]
    l2_valid <- is.finite(l2_vals) & l2_vals > 0
    l2_rescaled <- rep(0.5, n)
    l2_rescaled[l2_valid] <- 1 - rescale01(l2_vals[l2_valid])
    panels[["E: Low-lambda2 (inverted)"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = l2_rescaled, panel = "E: Low-lambda2 (inverted)")
  }

  if ("seed_nn_dist" %in% names(dt)) {
    nn_vals <- dt$seed_nn_dist
    nn_valid <- is.finite(nn_vals)
    nn_rescaled <- rep(0.5, n)
    nn_rescaled[nn_valid] <- 1 - rescale01(nn_vals[nn_valid])
    panels[["F: Low-NN distance (inverted)"]] <- data.table(
      MDS1 = mds_mat[, 1], MDS2 = mds_mat[, 2],
      value = nn_rescaled, panel = "F: Low-NN distance (inverted)")
  }

  if (length(panels) < 3) return(NULL)
  pdt <- rbindlist(panels)
  pdt <- pdt[is.finite(MDS1) & is.finite(MDS2) & is.finite(value)]
  pdt[, panel := factor(panel, levels = sort(unique(panel)))]

  n_panels <- length(unique(pdt$panel))
  # PC1: 2 columns for more vertical space per panel
  ncol_f <- min(2L, n_panels)

  p <- ggplot(pdt, aes(x = MDS1, y = MDS2, color = value)) +
    geom_point(size = 0.7, alpha = 0.55) +
    # PC1: viridis instead of default orange-blue
    scale_color_viridis_c(option = "C", name = "Signal (0-1)") +
    facet_wrap(~ panel, ncol = ncol_f, scales = "free") +
    labs(title = paste0(chr, " \u2014 MDS Faceted by All Metrics"),
         subtitle = paste0(nrow(dt), " windows | ", n_panels, " panels | all rescaled [0,1]"),
         caption = paste0("Red/yellow = high signal | Purple = low\n",
                         "E/F inverted (low \u03BB2 / low NN = high signal)\n",
                         "Coincident bright across panels = strongest inversion evidence")) +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold", size = 8),
          strip.background = element_rect(fill = "#F0F0F0", color = NA),
          panel.spacing = unit(0.6, "lines"),
          legend.key.width = unit(0.8, "cm"))
  p
}

# -- 26: Multi-source het assessment ------------------------------------------
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
    geom_point(size = 0.3, alpha = 0.4, color = "#3B82F6") +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " \u2014 Multi-Source Het Assessment"),
         subtitle = paste0(nrow(il_chr), " windows | chr-wide fixed k=3 bands, 6 metrics"),
         caption = "Coincident peaks A-E with trough F = inversion signal") +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold", size = 7),
          panel.spacing = unit(0.6, "lines"))
  p
}

# =============================================================================
# DISPATCH
# =============================================================================

plot_builders <- list(
  "07_mds_scatter_raw"            = build_mds_raw,
  "07b_mds_vs_position"           = build_mds_vs_position,
  "08_nn_histogram_adaptive"      = build_nn_histogram,
  "09_mds_by_inv_likeness"        = build_mds_inv,
  "10_mds_by_zscore"              = build_mds_z,
  "11_nn_calibration"             = build_nn_calibration,
  "13_knn_graph_structure"        = build_knn_structure,
  "18_mds_by_het_contrast"        = build_mds_het,
  "19_mds_by_lambda2"             = build_mds_lambda2,
  "20_mds_by_pve1_excess"         = build_mds_pve1_excess,
  "21_mds_by_nn_distance"         = build_mds_nn,
  "22_mds_by_eigenvalue_ratio"    = build_mds_eigratio,
  "23_mds_faceted_all_metrics"    = build_mds_faceted_all,
  "26_het_assessment_multi"       = build_het_assessment
)

plot_dims <- list(
  "07_mds_scatter_raw"            = c(w = 8, h = 7),
  # PC5: reduce width, increase height for squarer panels
  "07b_mds_vs_position"           = c(w = 16, h = 7),
  "08_nn_histogram_adaptive"      = c(w = 8, h = 5),
  "09_mds_by_inv_likeness"        = c(w = 8, h = 7),
  "10_mds_by_zscore"              = c(w = 8, h = 7),
  "11_nn_calibration"             = c(w = 14, h = 5),
  "13_knn_graph_structure"        = c(w = 14, h = 5),
  "18_mds_by_het_contrast"        = c(w = 8, h = 7),
  "19_mds_by_lambda2"             = c(w = 8, h = 7),
  "20_mds_by_pve1_excess"         = c(w = 8, h = 7),
  "21_mds_by_nn_distance"         = c(w = 8, h = 7),
  "22_mds_by_eigenvalue_ratio"    = c(w = 8, h = 7),
  # PC1: larger figure, 2-col layout
  "23_mds_faceted_all_metrics"    = c(w = 14, h = 18),
  "26_het_assessment_multi"       = c(w = 16, h = 14)
)

render_plots(plot_builders, plot_dims, chroms, precomp_list, cli$outdir)
message("[DONE] MDS plots -> ", cli$outdir)
