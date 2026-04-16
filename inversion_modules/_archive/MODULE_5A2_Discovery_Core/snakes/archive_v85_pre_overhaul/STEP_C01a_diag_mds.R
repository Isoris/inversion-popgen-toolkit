#!/usr/bin/env Rscript
# =============================================================================
# C01a_diag_mds.R — MDS scatter plots and sample-space diagnostics
# (plots 07-11, 13, 18-23, 26)
# Usage: Rscript C01a_diag_mds.R <precomp_dir> <outdir> [chrom]
# =============================================================================

SCRIPT_DIR <- dirname(sys.frame(1)$ofile %||% ".")
source(file.path(SCRIPT_DIR, "C01a_diag_common.R"))

cli <- parse_diag_args()
data <- load_diag_data(cli$precomp_dir, cli$chrom_filter)
precomp_list <- data$precomp_list
chroms <- data$chroms
inv_like_dt <- data$inv_like_dt
dir.create(cli$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# MDS BUILDERS
# =============================================================================

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

# =============================================================================
# DISPATCH
# =============================================================================

plot_builders <- list(
  "07_mds_scatter_raw" = build_mds_raw,
  "07b_mds_vs_position" = build_mds_vs_position,
  "08_nn_histogram_adaptive" = build_nn_histogram,
  "09_mds_by_inv_likeness" = build_mds_inv,
  "10_mds_by_zscore" = build_mds_z,
  "11_nn_calibration" = build_nn_calibration,
  "13_knn_graph_structure" = build_knn_structure,
  "18_mds_by_het_contrast" = build_mds_het,
  "19_mds_by_lambda2" = build_mds_lambda2,
  "20_mds_by_pve1_excess" = build_mds_pve1_excess,
  "21_mds_by_nn_distance" = build_mds_nn,
  "22_mds_by_eigenvalue_ratio" = build_mds_eigratio,
  "23_mds_faceted_all_metrics" = build_mds_faceted_all,
  "26_het_assessment_multi" = build_het_assessment
)

plot_dims <- list(
  "07_mds_scatter_raw" = c(w = 8, h = 7),
  "07b_mds_vs_position" = c(w = 18, h = 6),
  "08_nn_histogram_adaptive" = c(w = 8, h = 5),
  "09_mds_by_inv_likeness" = c(w = 8, h = 7),
  "10_mds_by_zscore" = c(w = 8, h = 7),
  "11_nn_calibration" = c(w = 14, h = 5),
  "13_knn_graph_structure" = c(w = 14, h = 5),
  "18_mds_by_het_contrast" = c(w = 8, h = 7),
  "19_mds_by_lambda2" = c(w = 8, h = 7),
  "20_mds_by_pve1_excess" = c(w = 8, h = 7),
  "21_mds_by_nn_distance" = c(w = 8, h = 7),
  "22_mds_by_eigenvalue_ratio" = c(w = 8, h = 7),
  "23_mds_faceted_all_metrics" = c(w = 16, h = 12),
  "26_het_assessment_multi" = c(w = 16, h = 14)
)

render_plots(plot_builders, plot_dims, chroms, precomp_list, cli$outdir)
message("[DONE] MDS plots -> ", cli$outdir)
