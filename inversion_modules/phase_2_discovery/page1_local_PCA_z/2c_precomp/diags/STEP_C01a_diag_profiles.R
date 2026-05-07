#!/usr/bin/env Rscript
# =============================================================================
# C01a_diag_profiles.R — Profile plots along chromosome
# (plots 01-05, 12, 24, 25, 27-32, N1_dashboard)
#
# v8.5.2 — Plot overhaul:
#   G1: sys.frame → env var sourcing
#   PC8:  Combined profile: taller, core boundaries, 5th panel
#   PC12: NN distance: ylim zoom, adaptive color, no fixed threshold
#   PC13: Background baseline: improved subtitle
#   PC14: Inv-likeness: SEED semantics, adaptive threshold, continuous color
#   PC15: Eigenvalue: smoothed, softer colors
#   N1:   Composite window summary dashboard
#
# Usage: Rscript C01a_diag_profiles.R <precomp_dir> <outdir> [chrom]
# =============================================================================

# G1: No sys.frame — resolve from env var
SCRIPT_DIR <- Sys.getenv("SNAKES_DIR",
  file.path(Sys.getenv("BASE", "."),
            "inversion_codebase_v8.5/MODULE_5A2_Discovery_Core/snakes"))
source(file.path(SCRIPT_DIR, "STEP_C01a_diag_common.R"))

cli <- parse_diag_args()
data <- load_diag_data(cli$precomp_dir, cli$chrom_filter)
precomp_list <- data$precomp_list
chroms <- data$chroms
window_dt <- data$window_dt
dir.create(cli$outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PROFILE BUILDERS
# =============================================================================

# -- 01: Eigenvalue spectrum (PC15: smoothed, softer colors) -----------------
build_eigen_profile <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (is.null(l1_col) || is.null(l2_col)) return(NULL)

  pdt <- data.table(
    pos_mb  = (dt$start_bp + dt$end_bp) / 2e6,
    lambda1 = dt[[l1_col]],
    lambda2 = dt[[l2_col]]
  )
  pdt[, pve1 := lambda1 / (lambda1 + lambda2)]
  pdt <- pdt[is.finite(lambda1) & is.finite(lambda2) & lambda1 > 0]
  if (nrow(pdt) < 10) return(NULL)

  smooth_w <- min(51L, nrow(pdt) %/% 5L)
  if (smooth_w %% 2 == 0) smooth_w <- smooth_w + 1L
  smooth_w <- max(5L, smooth_w)
  pdt[, lam1_smooth := frollmean(lambda1, n = smooth_w, align = "center", na.rm = TRUE)]
  pdt[, lam2_smooth := frollmean(lambda2, n = smooth_w, align = "center", na.rm = TRUE)]
  pdt[, pve1_smooth := frollmean(pve1, n = smooth_w, align = "center", na.rm = TRUE)]

  lam_max <- max(pdt$lambda1, na.rm = TRUE)
  if (!is.finite(lam_max) || lam_max == 0) lam_max <- 1

  p <- ggplot(pdt) +
    geom_point(aes(x = pos_mb, y = lambda1), color = COL_LAMBDA1, size = 0.15, alpha = 0.15) +
    geom_point(aes(x = pos_mb, y = lambda2), color = COL_LAMBDA2, size = 0.15, alpha = 0.15) +
    geom_line(aes(x = pos_mb, y = lam1_smooth), color = COL_LAMBDA1, linewidth = 0.6, na.rm = TRUE) +
    geom_line(aes(x = pos_mb, y = lam2_smooth), color = COL_LAMBDA2, linewidth = 0.6, na.rm = TRUE) +
    geom_line(aes(x = pos_mb, y = pve1_smooth * lam_max),
              color = COL_PVE1, linewidth = 0.7, na.rm = TRUE) +
    geom_hline(yintercept = 0.5 * lam_max, linetype = "dotted", color = "#9CA3AF", linewidth = 0.3) +
    scale_y_continuous(sec.axis = sec_axis(~ . / lam_max, name = "PVE1 (teal)")) +
    labs(x = paste0(chr, " (Mb)"), y = "Eigenvalue",
         title = paste0(chr, " \u2014 Eigenvalue Spectrum"),
         subtitle = paste0("Steel blue: \u03BB1 | Coral: \u03BB2 | Teal: PVE1 | ",
                          nrow(pdt), " windows, smoothed (w=", smooth_w, ")"),
         caption = paste0("Median PVE1: ", round(median(pdt$pve1, na.rm = TRUE), 3),
                         " | \u03BB1/\u03BB2 median: ",
                         round(median(pdt$lambda1/pdt$lambda2, na.rm = TRUE), 2),
                         " | Dotted = PVE1=0.5")) +
    THEME_BASE
  p
}

# -- 02: Inv-likeness profile (PC14: SEED semantics, continuous color) -------
build_inv_likeness_profile <- function(chr) {
  if (nrow(window_dt) == 0) return(NULL)
  il_chr <- window_dt[chrom == chr]
  if (nrow(il_chr) == 0) return(NULL)

  if (!"start_bp" %in% names(il_chr)) {
    dt <- precomp_list[[chr]]$dt
    if ("global_window_id" %in% names(dt) && "global_window_id" %in% names(il_chr)) {
      il_chr <- merge(il_chr, dt[, .(global_window_id, start_bp, end_bp)],
                       by = "global_window_id", all.x = TRUE)
    } else return(NULL)
  }

  il_chr <- il_chr[is.finite(start_bp) & is.finite(inv_likeness)]
  if (nrow(il_chr) == 0) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  chr_med <- median(il_chr$inv_likeness, na.rm = TRUE)
  chr_q75 <- quantile(il_chr$inv_likeness, 0.75, na.rm = TRUE)

  il_chr[, status := fifelse(inv_likeness >= chr_q75, "SEED_ELIGIBLE",
                     fifelse(inv_likeness >= chr_med, "SEED_MARGINAL", "SEED_INELIGIBLE"))]

  p <- ggplot(il_chr, aes(x = pos_mb, y = inv_likeness, color = inv_likeness)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_gradientn(
      colours = c("#E5E7EB", "#93C5FD", "#2563EB", "#D97706", "#DC2626"),
      values = scales::rescale(c(0, chr_med * 0.5, chr_med, chr_q75, 1)),
      name = "Inv-likeness"
    ) +
    geom_hline(yintercept = chr_med, linetype = "dashed", color = "#6B7280", linewidth = 0.3) +
    geom_hline(yintercept = chr_q75, linetype = "dashed", color = "#2563EB", linewidth = 0.3) +
    annotate("text", x = min(il_chr$pos_mb), y = chr_med + 0.02,
             label = paste0("median=", round(chr_med, 3)), size = 2.5, hjust = 0, color = "#6B7280") +
    annotate("text", x = min(il_chr$pos_mb), y = chr_q75 + 0.02,
             label = paste0("q75=", round(chr_q75, 3)), size = 2.5, hjust = 0, color = "#2563EB") +
    labs(x = paste0(chr, " (Mb)"), y = "Inv-likeness",
         title = paste0(chr, " \u2014 Inversion-Likeness Score"),
         subtitle = paste0(nrow(il_chr), " windows | ",
                          sum(il_chr$status == "SEED_ELIGIBLE"), " eligible, ",
                          sum(il_chr$status == "SEED_MARGINAL"), " marginal, ",
                          sum(il_chr$status == "SEED_INELIGIBLE"), " ineligible"),
         caption = paste0("inv_likeness = PVE1 (0.35) + eigen_ratio (0.25) + trimodality (0.40)\n",
                         "Seed gate selects windows eligible for core extension\n",
                         "Blue dashed = per-chr adaptive q75 | Grey dashed = median")) +
    THEME_BASE
  p
}

# -- 03: Robust z-score profile -----------------------------------------------
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
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = PAL_Z) +
    geom_hline(yintercept = 1.2, linetype = "dashed", color = COL_Q85, linewidth = 0.3) +
    geom_hline(yintercept = 1.8, linetype = "dashed", color = COL_Q90, linewidth = 0.3) +
    geom_hline(yintercept = 2.5, linetype = "dashed", color = "#2563EB", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "|Robust z|",
         title = paste0(chr, " \u2014 Robust Z-Score Profile"),
         subtitle = paste0(nrow(pdt), " windows | z\u22652.5: ", sum(pdt$z >= 2.5),
                          " | z\u22651.8: ", sum(pdt$z >= 1.8),
                          " | z\u22651.2: ", sum(pdt$z >= 1.2)),
         caption = "Thresholds: S1S=2.5 (blue) | S1M=1.8 (teal) | S1L=1.2 (amber)") +
    THEME_BASE
  p
}

# -- 04: Background continuity baseline (PC13: clear explanation) -------------
build_bg_baseline <- function(chr) {
  pc <- precomp_list[[chr]]
  bg_q <- pc$bg_continuity_quantiles
  if (is.null(bg_q) || length(bg_q) < 4) return(NULL)
  if (is.null(names(bg_q))) names(bg_q) <- paste0(c(50, 75, 80, 85, 90, 95), "%")[seq_along(bg_q)]

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
    geom_histogram(bins = 80, fill = "#DBEAFE", color = "#3B82F6", linewidth = 0.2) +
    geom_vline(xintercept = as.numeric(bg_q["50%"]), color = "#6B7280", linetype = "dashed") +
    geom_vline(xintercept = as.numeric(bg_q["85%"]), color = COL_Q85, linetype = "dashed") +
    geom_vline(xintercept = as.numeric(bg_q["90%"]), color = COL_Q90, linetype = "dashed") +
    geom_vline(xintercept = as.numeric(bg_q["95%"]), color = COL_Q95, linetype = "dashed") +
    annotate("text", x = as.numeric(bg_q["85%"]), y = Inf, vjust = 1.5,
             label = paste0("q85=", round(bg_q["85%"], 3)), size = 2.5, color = COL_Q85) +
    annotate("text", x = as.numeric(bg_q["90%"]), y = Inf, vjust = 3,
             label = paste0("q90=", round(bg_q["90%"], 3)), size = 2.5, color = COL_Q90) +
    annotate("text", x = as.numeric(bg_q["95%"]), y = Inf, vjust = 4.5,
             label = paste0("q95=", round(bg_q["95%"], 3)), size = 2.5, color = COL_Q95) +
    labs(x = "Adjacent-window similarity", y = "Count",
         title = paste0(chr, " \u2014 Background Continuity Baseline"),
         subtitle = paste0("Adjacent-window similarity for background (non-inversion) windows\n",
                          "Adaptive thresholds define min similarity for core membership:\n",
                          "S1S = q95 (strictest) | S1M = q90 (moderate) | S1L = q85 (generous)"),
         caption = paste0("q50=", round(bg_q["50%"], 3),
                         " q75=", round(bg_q["75%"], 3),
                         " q85=", round(bg_q["85%"], 3),
                         " q90=", round(bg_q["90%"], 3),
                         " q95=", round(bg_q["95%"], 3))) +
    THEME_BASE
  p
}

# -- 05: NN distance profile (PC12: zoom, adaptive color, no fixed 0.70) -----
build_nn_profile <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  if (!"seed_nn_dist" %in% names(dt)) return(NULL)

  pdt <- dt[, .(pos_mb = (start_bp + end_bp) / 2e6, nn = seed_nn_dist)]
  pdt <- pdt[is.finite(nn)]
  if (nrow(pdt) < 30) return(NULL)

  qq <- quantile(pdt$nn, c(0.25, 0.50, 0.75), na.rm = TRUE)
  y_cap <- quantile(pdt$nn, 0.99, na.rm = TRUE) * 1.1

  pdt[, nn_class := fifelse(nn <= qq["25%"], "S1S (q25)",
                   fifelse(nn <= qq["50%"], "S1M (q50)",
                   fifelse(nn <= qq["75%"], "S1L (q75)", "above q75")))]

  smooth_w <- min(31L, nrow(pdt) %/% 10L)
  if (smooth_w %% 2 == 0) smooth_w <- smooth_w + 1L
  smooth_w <- max(5L, smooth_w)
  pdt[, nn_smooth := frollmedian(nn, n = smooth_w, align = "center", na.rm = TRUE)]

  pdt[, nn_state := fifelse(nn_smooth <= qq["25%"], "core_low",
                   fifelse(nn_smooth <= qq["50%"], "mid", "background"))]
  pdt[, nn_state := nafill(nn_state, type = "locf")]
  pdt[is.na(nn_state), nn_state := "background"]
  pdt[, run_id := rleid(nn_state)]

  runs <- pdt[, .(state = first(nn_state), start_mb = min(pos_mb),
                   end_mb = max(pos_mb), n_win = .N), by = run_id]
  runs <- runs[n_win >= 10]

  p <- ggplot(pdt, aes(x = pos_mb, y = nn)) +
    {if (nrow(runs) > 0) geom_rect(
      data = runs, inherit.aes = FALSE,
      aes(xmin = start_mb, xmax = end_mb, ymin = -Inf, ymax = Inf, fill = state),
      alpha = 0.12
    )} +
    scale_fill_manual(values = PAL_REGIME, name = "Regime") +
    geom_point(aes(color = nn_class), size = 0.4, alpha = 0.5) +
    scale_color_manual(values = PAL_ADAPT, name = "Adaptive class") +
    geom_line(aes(y = nn_smooth), color = "#1A1A1A", linewidth = 0.5, na.rm = TRUE) +
    geom_hline(yintercept = qq["25%"], linetype = "dashed", color = "#2563EB", linewidth = 0.3) +
    geom_hline(yintercept = qq["50%"], linetype = "dashed", color = "#059669", linewidth = 0.3) +
    geom_hline(yintercept = qq["75%"], linetype = "dashed", color = "#D97706", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "NN Distance (MDS-space)",
         title = paste0(chr, " \u2014 Nearest-Neighbor Distance Along Chromosome"),
         subtitle = paste0(nrow(pdt), " windows | Adaptive: q25=", round(qq["25%"], 3),
                          " q50=", round(qq["50%"], 3), " q75=", round(qq["75%"], 3)),
         caption = paste0("Low NN = tightly clustered = inversion core\n",
                         "Blue band = core_low (<q25) | Green = mid | Grey = background\n",
                         "Black line = rolling median (w=", smooth_w, ")")) +
    coord_cartesian(ylim = c(0, y_cap)) +
    THEME_BASE
  p
}

# -- 12: Combined 5-panel profile (PC8: taller, 5th panel) -------------------
build_combined_profile <- function(chr) {
  dt <- precomp_list[[chr]]$dt
  n <- nrow(dt)
  pos_mb <- (dt$start_bp + dt$end_bp) / 2e6

  panels <- list()
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (!is.null(l1_col) && !is.null(l2_col)) {
    pve1 <- dt[[l1_col]] / (dt[[l1_col]] + dt[[l2_col]])
    panels[["A_PVE1"]] <- data.table(pos_mb = pos_mb, value = pve1, panel = "A: PVE1")
  }
  if ("max_abs_z" %in% names(dt))
    panels[["B_zscore"]] <- data.table(pos_mb = pos_mb, value = dt$max_abs_z, panel = "B: |z-score|")
  if ("inv_likeness" %in% names(dt)) {
    panels[["C_inv"]] <- data.table(pos_mb = pos_mb, value = dt$inv_likeness, panel = "C: Inv-likeness")
  } else if (nrow(window_dt) > 0) {
    il_chr <- merge(data.table(global_window_id = dt$global_window_id, pos_mb = pos_mb),
                     window_dt[chrom == chr, .(global_window_id, inv_likeness)],
                     by = "global_window_id", all.x = TRUE)
    panels[["C_inv"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$inv_likeness, panel = "C: Inv-likeness")
  }
  if ("seed_nn_dist" %in% names(dt))
    panels[["D_nn"]] <- data.table(pos_mb = pos_mb, value = dt$seed_nn_dist, panel = "D: NN distance")

  # PC8: 5th panel — diagonal sim summary
  sim_mat <- precomp_list[[chr]]$sim_mat
  if (!is.null(sim_mat) && nrow(sim_mat) >= 20) {
    bw <- 20L
    diag_mean <- vapply(seq_len(n), function(i) {
      j_lo <- max(1L, i - bw); j_hi <- min(n, i + bw)
      mean(sim_mat[i, j_lo:j_hi], na.rm = TRUE)
    }, numeric(1))
    panels[["E_diag"]] <- data.table(pos_mb = pos_mb, value = diag_mean,
                                      panel = "E: Diagonal sim (\u00b120 win)")
  }

  if (length(panels) < 2) return(NULL)
  pdt <- rbindlist(panels)
  pdt <- pdt[is.finite(value)]
  pdt[, panel := factor(panel, levels = c("A: PVE1", "B: |z-score|",
                                           "C: Inv-likeness", "D: NN distance",
                                           "E: Diagonal sim (\u00b120 win)"))]
  pdt[, value_q := frank(value, na.last = "keep") / .N, by = panel]

  p <- ggplot(pdt, aes(x = pos_mb, y = value, color = value_q)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_gradientn(colours = c("#E5E7EB", "#3B82F6", "#F59E0B", "#DC2626"),
                          name = "Quantile") +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " \u2014 Combined Signal Profile"),
         subtitle = paste0(n, " windows \u2014 red = high signal, grey = low"),
         caption = "Coincident red A+B+C with blue D = strong inversion\nE: mean local similarity within \u00b120 window band") +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "#F0F0F0", color = NA),
          panel.spacing = unit(1.0, "lines"))
  p
}


# -- 24: Local entropy profile ------------------------------------------------
build_entropy_profile <- function(chr) {
  il_chr <- if (nrow(window_dt) > 0) window_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0 || !"local_entropy" %in% names(il_chr)) return(NULL)
  il_chr <- il_chr[is.finite(local_entropy) & is.finite(start_bp)]
  if (nrow(il_chr) < 20) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  p <- ggplot(il_chr, aes(x = pos_mb)) +
    geom_point(aes(y = local_entropy), size = 0.4, alpha = 0.5, color = "#D97706") +
    geom_hline(yintercept = log(3), linetype = "dashed", color = "#9CA3AF", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "Local Shannon entropy (H)",
         title = paste0(chr, " \u2014 Local Structure Entropy"),
         subtitle = paste0(nrow(il_chr), " windows | chr-wide fixed k=3 bands"),
         caption = "Low H = single group | High H = mixed | Dashed = ln(3) = max k=3") +
    THEME_BASE
  p
}

# -- 25: ENA profile -----------------------------------------------------------
build_ena_profile <- function(chr) {
  il_chr <- if (nrow(window_dt) > 0) window_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0 || !"local_ena" %in% names(il_chr)) return(NULL)
  il_chr <- il_chr[is.finite(local_ena) & is.finite(start_bp)]
  if (nrow(il_chr) < 20) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  p <- ggplot(il_chr, aes(x = pos_mb, y = local_ena)) +
    geom_point(size = 0.4, alpha = 0.5, color = "#DC2626") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#059669", linewidth = 0.3) +
    geom_hline(yintercept = 3, linetype = "dotted", color = "#9CA3AF", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "ENA",
         title = paste0(chr, " \u2014 ENA Profile"),
         subtitle = paste0(nrow(il_chr), " windows | ENA = exp(H)"),
         caption = "ENA \u2248 1 = single background | ENA \u2248 3 = max mixing k=3") +
    THEME_BASE
  p
}

# -- 27: Three-track decomposition -------------------------------------------
build_three_track_profile <- function(chr) {
  il_chr <- if (nrow(window_dt) > 0) window_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0) return(NULL)
  has_str <- "structure_likeness" %in% names(il_chr) && any(is.finite(il_chr$structure_likeness))
  if (!has_str) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  tracks <- list()
  tracks[["Structure"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$structure_likeness,
                                       track = "A: Structure-likeness")
  tracks[["Inversion"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$inv_likeness,
                                       track = "B: Inv-likeness")
  if ("family_likeness" %in% names(il_chr) && any(is.finite(il_chr$family_likeness)))
    tracks[["Family"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$family_likeness,
                                      track = "C: Family-likeness")

  pdt <- rbindlist(tracks)
  pdt <- pdt[is.finite(value)]
  pdt[, track := factor(track, levels = unique(track))]

  p <- ggplot(pdt, aes(x = pos_mb, y = value)) +
    geom_point(size = 0.3, alpha = 0.4, color = "#3B82F6") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "#9CA3AF", linewidth = 0.3) +
    facet_wrap(~ track, ncol = 1, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = "Score [0-1]",
         title = paste0(chr, " \u2014 Three-Track Decomposition"),
         subtitle = "A: all structure | B: inversion-specific | C: family-driven",
         caption = "Ideal inversion: A high + B high + C low | Family block: A+C high") +
    THEME_BASE +
    theme(strip.text = element_text(face = "bold", size = 8), panel.spacing = unit(0.5, "lines"))
  p
}

# -- 28: Family-likeness profile -----------------------------------------------
build_family_profile <- function(chr) {
  il_chr <- if (nrow(window_dt) > 0) window_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0 || !"family_likeness" %in% names(il_chr)) return(NULL)
  il_chr <- il_chr[is.finite(family_likeness)]
  if (nrow(il_chr) < 20) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  p <- ggplot(il_chr, aes(x = pos_mb, y = family_likeness)) +
    geom_point(size = 0.4, alpha = 0.5, color = "#7C3AED") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "#DC2626", linewidth = 0.3) +
    labs(x = paste0(chr, " (Mb)"), y = "Family-likeness",
         title = paste0(chr, " \u2014 Family-Likeness Profile"),
         subtitle = paste0(nrow(il_chr), " windows | Q variance by k=3 PC1 bands"),
         caption = "High = family-driven grouping | Low = structural signal (inversion)") +
    THEME_BASE
  p
}

# -- 29: Dosage het rate profile ------------------------------------------------
build_dosage_het_profile <- function(chr) {
  il_chr <- if (nrow(window_dt) > 0) window_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0) return(NULL)
  has_mean <- "dosage_het_rate_median" %in% names(il_chr) && any(is.finite(il_chr$dosage_het_rate_median))
  has_cv <- "dosage_het_rate_cv" %in% names(il_chr) && any(is.finite(il_chr$dosage_het_rate_cv))
  if (!has_mean && !has_cv) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  tracks <- list()
  if (has_mean)
    tracks[["median"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$dosage_het_rate_median,
      track = "A: Median het rate (226 samples)")
  if (has_cv)
    tracks[["cv"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$dosage_het_rate_cv,
      track = "B: Het rate CV (bimodality signal)")
  if ("dosage_het_rate_sd" %in% names(il_chr) && any(is.finite(il_chr$dosage_het_rate_sd)))
    tracks[["sd"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$dosage_het_rate_sd,
      track = "C: Het rate SD")

  pdt <- rbindlist(tracks)
  pdt <- pdt[is.finite(value)]
  pdt[, track := factor(track, levels = unique(track))]

  p <- ggplot(pdt, aes(x = pos_mb, y = value)) +
    geom_point(size = 0.3, alpha = 0.4, color = "#059669") +
    facet_wrap(~ track, ncol = 1, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " \u2014 Dosage Het Rate Profile"),
         subtitle = "From raw dosage matrix (independent from PCA)",
         caption = "A: MEDIAN het | B: CV = sd/median (high CV = inversion) | C: raw SD") +
    THEME_BASE + theme(strip.text = element_text(face = "bold", size = 8))
  p
}

# -- 30: Band fraction profile -------------------------------------------------
build_band_fraction_profile <- function(chr) {
  il_chr <- if (nrow(window_dt) > 0) window_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0 || !"band1_frac" %in% names(il_chr)) return(NULL)
  il_chr <- il_chr[is.finite(band1_frac)]
  if (nrow(il_chr) < 20) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  pdt <- melt(il_chr[, .(pos_mb, band1_frac, band2_frac, band3_frac)],
              id.vars = "pos_mb", variable.name = "band", value.name = "fraction")
  pdt[, band := factor(band, levels = c("band1_frac", "band2_frac", "band3_frac"),
                         labels = c("Band 1 (low PC1)", "Band 2 (mid/het)", "Band 3 (high PC1)"))]

  p <- ggplot(pdt[is.finite(fraction)], aes(x = pos_mb, y = fraction, color = band)) +
    geom_point(size = 0.3, alpha = 0.3) +
    scale_color_manual(values = PAL_BAND) +
    geom_hline(yintercept = 1/3, linetype = "dotted", color = "#9CA3AF") +
    labs(x = paste0(chr, " (Mb)"), y = "Fraction of samples",
         title = paste0(chr, " \u2014 Band Fraction Profile (k=3)"),
         subtitle = "Fraction of 226 samples in each chr-wide k=3 band per window",
         caption = "Balanced ~33% = common inversion | Unbalanced = rare/weak | Dotted = 1/3",
         color = NULL) +
    THEME_BASE
  p
}

# -- 31: Local Q profile -------------------------------------------------------
build_localQ_profile <- function(chr) {
  il_chr <- if (nrow(window_dt) > 0) window_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0) return(NULL)
  has_lq <- "localQ_delta12" %in% names(il_chr) && any(is.finite(il_chr$localQ_delta12))
  if (!has_lq) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  tracks <- list()
  tracks[["d12"]] <- data.table(pos_mb = il_chr$pos_mb,
    value = il_chr$localQ_delta12, track = "A: Local Q \u039412 (dominance)")
  if ("localQ_entropy" %in% names(il_chr))
    tracks[["H"]] <- data.table(pos_mb = il_chr$pos_mb,
      value = il_chr$localQ_entropy, track = "B: Local Q entropy (H)")
  if ("localQ_ena" %in% names(il_chr))
    tracks[["ena"]] <- data.table(pos_mb = il_chr$pos_mb,
      value = il_chr$localQ_ena, track = "C: Local Q ENA")

  pdt <- rbindlist(tracks)
  pdt <- pdt[is.finite(value)]
  pdt[, track := factor(track, levels = unique(track))]

  p <- ggplot(pdt, aes(x = pos_mb, y = value)) +
    geom_point(size = 0.3, alpha = 0.4, color = "#D97706") +
    facet_wrap(~ track, ncol = 1, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " \u2014 Local Q Profile (ancestry_bridge)"),
         subtitle = "Real admixture proportions (independent from PCA)",
         caption = "A: \u039412 per window | B: Shannon entropy | C: ENA = exp(H)") +
    THEME_BASE + theme(strip.text = element_text(face = "bold", size = 8))
  p
}

# -- 32: PC1 het metrics profile ------------------------------------------------
build_het_metrics_profile <- function(chr) {
  il_chr <- if (nrow(window_dt) > 0) window_dt[chrom == chr] else data.table()
  if (nrow(il_chr) == 0) return(NULL)
  il_chr[, pos_mb := (start_bp + end_bp) / 2e6]

  tracks <- list()
  if ("het_pc1_gap" %in% names(il_chr) && any(is.finite(il_chr$het_pc1_gap)))
    tracks[["gap"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$het_pc1_gap,
      track = "A: PC1 gap (min dist between band centers)")
  if ("het_mid_fraction" %in% names(il_chr) && any(is.finite(il_chr$het_mid_fraction)))
    tracks[["mid"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$het_mid_fraction,
      track = "B: Middle band fraction (expected ~0.5 for inv)")
  if ("pc1_bimodality" %in% names(il_chr) && any(is.finite(il_chr$pc1_bimodality)))
    tracks[["bimod"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$pc1_bimodality,
      track = "C: PC1 bimodality (Ashman D)")
  if ("local_delta12" %in% names(il_chr) && any(is.finite(il_chr$local_delta12)))
    tracks[["d12"]] <- data.table(pos_mb = il_chr$pos_mb, value = il_chr$local_delta12,
      track = "D: \u039412 from fixed-band soft membership")

  if (length(tracks) < 2) return(NULL)
  pdt <- rbindlist(tracks)
  pdt <- pdt[is.finite(value)]
  pdt[, track := factor(track, levels = unique(track))]

  p <- ggplot(pdt, aes(x = pos_mb, y = value)) +
    geom_point(size = 0.3, alpha = 0.4, color = "#3B82F6") +
    facet_wrap(~ track, ncol = 1, scales = "free_y") +
    labs(x = paste0(chr, " (Mb)"), y = NULL,
         title = paste0(chr, " \u2014 Het Metrics Profile (PC1-derived)"),
         subtitle = "Multiple het assessment channels from k=3 band analysis",
         caption = "A: gap between centers | B: middle band fraction | C: Ashman D | D: \u039412") +
    THEME_BASE + theme(strip.text = element_text(face = "bold", size = 7))
  p
}

# =============================================================================
# N1: Composite Window Summary Dashboard (per chromosome)
# =============================================================================
build_dashboard <- function(chr) {
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)
  library(patchwork)

  dt <- precomp_list[[chr]]$dt
  n <- nrow(dt)
  pos_mb <- (dt$start_bp + dt$end_bp) / 2e6

  # Panel B: inv-likeness + smooth
  p_inv <- NULL
  inv_vals <- if ("inv_likeness" %in% names(dt)) dt$inv_likeness else {
    if (nrow(window_dt) > 0 && "global_window_id" %in% names(dt)) {
      m <- merge(data.table(global_window_id = dt$global_window_id),
                 window_dt[chrom == chr, .(global_window_id, inv_likeness)],
                 by = "global_window_id", all.x = TRUE)
      m$inv_likeness
    } else NULL
  }
  if (!is.null(inv_vals)) {
    inv_dt <- data.table(pos_mb = pos_mb, inv = inv_vals)[is.finite(inv)]
    if (nrow(inv_dt) > 20) {
      p_inv <- ggplot(inv_dt, aes(x = pos_mb, y = inv)) +
        geom_point(size = 0.3, alpha = 0.3, color = "#2563EB") +
        geom_smooth(method = "loess", span = 0.1, se = FALSE,
                    color = "#1E3A5F", linewidth = 0.6) +
        labs(y = "Inv-likeness") + THEME_BASE +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
  }

  # Panel E: NN distance
  p_nn <- NULL
  if ("seed_nn_dist" %in% names(dt)) {
    nn_dt <- data.table(pos_mb = pos_mb, nn = dt$seed_nn_dist)[is.finite(nn)]
    if (nrow(nn_dt) > 20) {
      y_cap <- quantile(nn_dt$nn, 0.99, na.rm = TRUE) * 1.1
      p_nn <- ggplot(nn_dt, aes(x = pos_mb, y = nn)) +
        geom_point(size = 0.3, alpha = 0.3, color = "#7C3AED") +
        coord_cartesian(ylim = c(0, y_cap)) +
        labs(y = "NN dist") + THEME_BASE +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    }
  }

  # Panel F: eigenvalue ratio
  p_eigr <- NULL
  l1_col <- if ("lam_1" %in% names(dt)) "lam_1" else if ("lambda1" %in% names(dt)) "lambda1" else NULL
  l2_col <- if ("lam_2" %in% names(dt)) "lam_2" else if ("lambda2" %in% names(dt)) "lambda2" else NULL
  if (!is.null(l1_col) && !is.null(l2_col)) {
    eig_dt <- data.table(pos_mb = pos_mb, ratio = dt[[l1_col]] / dt[[l2_col]])[is.finite(ratio)]
    if (nrow(eig_dt) > 20) {
      p_eigr <- ggplot(eig_dt, aes(x = pos_mb, y = pmin(ratio, 15))) +
        geom_point(size = 0.3, alpha = 0.3, color = "#059669") +
        geom_hline(yintercept = 3, linetype = "dashed", color = "#9CA3AF", linewidth = 0.3) +
        labs(x = paste0(chr, " (Mb)"), y = "\u03BB1/\u03BB2") + THEME_BASE
    }
  }

  plots <- Filter(Negate(is.null), list(p_inv, p_nn, p_eigr))
  if (length(plots) < 2) return(NULL)

  wrap_plots(plots, ncol = 1) +
    plot_annotation(
      title = paste0(chr, " \u2014 Window Summary Dashboard"),
      subtitle = paste0(n, " windows | All panels share x-axis (Mb)"),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", color = "#1A1A1A"),
        plot.subtitle = element_text(size = 10, color = "#5A5A5A"))
    )
}

# =============================================================================
# DISPATCH
# =============================================================================

plot_builders <- list(
  "01_eigenvalue_spectrum"       = build_eigen_profile,
  "02_inv_likeness_profile"      = build_inv_likeness_profile,
  "03_max_abs_z_profile"         = build_z_profile,
  "04_bg_continuity_baseline"    = build_bg_baseline,
  "05_nn_distance_profile"       = build_nn_profile,
  "12_combined_profile"          = build_combined_profile,
  "24_local_entropy_profile"     = build_entropy_profile,
  "25_ena_profile"               = build_ena_profile,
  "27_three_track_profile"       = build_three_track_profile,
  "28_family_likeness_profile"   = build_family_profile,
  "29_dosage_het_rate_profile"   = build_dosage_het_profile,
  "30_band_fraction_profile"     = build_band_fraction_profile,
  "31_local_Q_profile"           = build_localQ_profile,
  "32_het_metrics_profile"       = build_het_metrics_profile,
  "N1_window_summary_dashboard"  = build_dashboard
)

plot_dims <- list(
  "01_eigenvalue_spectrum"       = c(w = 14, h = 5),
  "02_inv_likeness_profile"      = c(w = 14, h = 5),
  "03_max_abs_z_profile"         = c(w = 14, h = 5),
  "04_bg_continuity_baseline"    = c(w = 8, h = 6),
  "05_nn_distance_profile"       = c(w = 14, h = 6),
  "12_combined_profile"          = c(w = 16, h = 20),
  "24_local_entropy_profile"     = c(w = 14, h = 5),
  "25_ena_profile"               = c(w = 14, h = 5),
  "27_three_track_profile"       = c(w = 14, h = 10),
  "28_family_likeness_profile"   = c(w = 14, h = 5),
  "29_dosage_het_rate_profile"   = c(w = 14, h = 10),
  "30_band_fraction_profile"     = c(w = 14, h = 5),
  "31_local_Q_profile"           = c(w = 14, h = 10),
  "32_het_metrics_profile"       = c(w = 14, h = 12),
  "N1_window_summary_dashboard"  = c(w = 14, h = 14)
)

render_plots(plot_builders, plot_dims, chroms, precomp_list, cli$outdir)
message("[DONE] Profiles -> ", cli$outdir)
