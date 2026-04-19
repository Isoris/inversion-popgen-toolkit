#!/usr/bin/env Rscript
# =============================================================================
# q04_compose_plot.R  (v3 — polished with chromosome-browser styling)
#
# Composes a multi-track per-chromosome diagnostic PDF with genome-browser look:
#
#   Top strip 1: ideogram of chromosome with shelf highlighted in red
#   Top strip 2: sim_mat 1D collapse (row-means near-diagonal similarity)
#
#   Track 1: Robust Z
#   Track 2: SNPs per local-PCA window
#   Track 3: BEAGLE posterior uncertainty
#   Track 4: Mean coverage + CV
#   Track 5: # samples low-cov (bars)
#   Track 6: theta_pi mean + CV across samples
#   Track 7: ancestry delta12 + regime-class strip above
#   Track 7b: multi-scale delta12 (1x/5x/10x)
#   Track 8: theta_pi by inv-genotype (Hom1/Het/Hom2)
#   Track 9: Hom1-vs-Hom2 Fst
#
# Style touches applied uniformly:
#   - double-stroke lines (thick light halo + thin dark core) on smoothed tracks
#   - red dashed breakpoint vertical lines across all tracks if specified
#   - amber shelf shading across all tracks
#   - right-edge interpretation labels on Δ12/H/ENA-style tracks
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
PRECOMP    <- get_arg("--precomp")
CHROM      <- get_arg("--chrom")
SNP_FILE   <- get_arg("--snp_track")
UNC_FILE   <- get_arg("--unc_track")
COV_FILE   <- get_arg("--cov_track")
THETA_FILE <- get_arg("--theta_track")
ANC_FILE   <- get_arg("--ancestry_track")
ANC_MULTI  <- get_arg("--ancestry_multiscale", NULL)
POPSTATS_INVGT_FILE <- get_arg("--popstats_invgt", NULL)
INVGT_ASSIGN_FILE   <- get_arg("--invgt_assignments", NULL)
SIMMAT_FILE <- get_arg("--sim_mat", NULL)
OUT        <- get_arg("--out")
SHELF_A    <- as.numeric(get_arg("--shelf_start_mb", NA))
SHELF_B    <- as.numeric(get_arg("--shelf_end_mb",   NA))
BP1_MB     <- as.numeric(get_arg("--breakpoint1_mb", NA))
BP2_MB     <- as.numeric(get_arg("--breakpoint2_mb", NA))
SMOOTH_WIN <- as.integer(get_arg("--smooth_win", "1"))

stopifnot(!is.null(PRECOMP), file.exists(PRECOMP))

rolling_median <- function(x, k) {
  if (is.na(k) || k <= 1) return(x)
  if (k %% 2 == 0) k <- k + 1
  n <- length(x); half <- (k - 1L) %/% 2L
  out <- x
  for (i in seq_len(n)) {
    lo <- max(1L, i - half); hi <- min(n, i + half)
    out[i] <- stats::median(x[lo:hi], na.rm = TRUE)
  }
  out
}

# =============================================================================
# Styling helpers
# =============================================================================

# Breakpoint dashed verticals: infer from shelf if BP not explicit
resolve_bps <- function() {
  bps <- c()
  if (is.finite(BP1_MB)) bps <- c(bps, BP1_MB)
  if (is.finite(BP2_MB)) bps <- c(bps, BP2_MB)
  if (length(bps) == 0 && is.finite(SHELF_A) && is.finite(SHELF_B)) {
    bps <- c(SHELF_A, SHELF_B)
  }
  bps
}
BPS <- resolve_bps()

# Shelf shading + breakpoint lines applied to every panel
decorate <- function(p) {
  if (is.finite(SHELF_A) && is.finite(SHELF_B)) {
    p <- p + annotate("rect", xmin = SHELF_A, xmax = SHELF_B,
                      ymin = -Inf, ymax = Inf,
                      fill = "#f5a524", alpha = 0.10)
  }
  if (length(BPS) > 0) {
    p <- p + geom_vline(xintercept = BPS, linetype = "dashed",
                        color = "#e0555c", linewidth = 0.35, alpha = 0.8)
  }
  p
}

# Double-stroke line: thick halo underneath, thin core on top
double_stroke <- function(halo_color = "#ffffff", core_color = "#1f4e79",
                          halo_width = 1.4, core_width = 0.45,
                          halo_alpha = 0.75) {
  list(
    geom_line(color = halo_color, linewidth = halo_width, alpha = halo_alpha),
    geom_line(color = core_color, linewidth = core_width)
  )
}

# Right-edge interpretation labels (top/bottom)
edge_labels <- function(top_label, bot_label, top_color = "#2c3e50",
                        bot_color = "#2c3e50") {
  # Use annotation to draw; relies on x-axis being set correctly later.
  # Returns a list of grobs to add.
  list(
    annotate("text", x = Inf, y = Inf, label = top_label,
             hjust = -0.05, vjust = 1.4, size = 2.2,
             color = top_color, family = "mono"),
    annotate("text", x = Inf, y = -Inf, label = bot_label,
             hjust = -0.05, vjust = -0.5, size = 2.2,
             color = bot_color, family = "mono"),
    coord_cartesian(clip = "off")
  )
}

# Theme: track body vs track last-in-stack (has x-axis)
base_theme <- theme_classic(base_size = 9) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(2, 46, 2, 10),   # extra right padding for edge labels
        panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.3))

bottom_theme <- theme_classic(base_size = 9) +
  theme(plot.margin = margin(2, 46, 5, 10),
        panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.3))

placeholder <- function(label) {
  ggplot() + labs(y = label) + theme_void() +
    theme(plot.margin = margin(2, 46, 2, 10))
}

# =============================================================================
# Load Z from precomp
# =============================================================================
pc <- readRDS(PRECOMP)
dt <- as.data.table(pc$dt)
dt[, mb := (start_bp + end_bp) / 2 / 1e6]
z_col <- NULL
for (cand in c("robust_z", "z_robust", "z", "mds_z_robust", "mds_z1_robust")) {
  if (cand %in% names(dt)) { z_col <- cand; break }
}
if (is.null(z_col)) stop("No Z column in precomp$dt")
dt[, z := get(z_col)]

chrom_start_mb <- min(dt$mb, na.rm = TRUE)
chrom_end_mb   <- max(dt$mb, na.rm = TRUE)
common_x <- scale_x_continuous(limits = c(chrom_start_mb, chrom_end_mb),
                               expand = expansion(mult = 0))

# =============================================================================
# TOP STRIP 1: Ideogram
# =============================================================================
ideogram <- {
  ideo_df <- data.frame(x0 = chrom_start_mb, x1 = chrom_end_mb,
                        y0 = 0, y1 = 1)
  p <- ggplot() +
    # Base chromosome body (rounded look via tall rect with grey)
    geom_rect(data = ideo_df, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
              fill = "#e8ebef", color = "#555e69", linewidth = 0.3) +
    # Faint ticks at every 2 Mb
    geom_segment(data = data.frame(x = seq(ceiling(chrom_start_mb),
                                           floor(chrom_end_mb), by = 2)),
                 aes(x = x, xend = x, y = 0.02, yend = 0.98),
                 color = "#c8cdd2", linewidth = 0.25)
  if (is.finite(SHELF_A) && is.finite(SHELF_B)) {
    p <- p + geom_rect(aes(xmin = SHELF_A, xmax = SHELF_B, ymin = 0, ymax = 1),
                       fill = "#e0555c", alpha = 0.55, color = "#8a2b30",
                       linewidth = 0.3)
  }
  if (length(BPS) > 0) {
    p <- p + geom_segment(data = data.frame(x = BPS),
                          aes(x = x, xend = x, y = -0.2, yend = 1.2),
                          color = "#e0555c", linewidth = 0.5)
  }
  p + common_x +
    scale_y_continuous(limits = c(-0.3, 1.3), expand = c(0,0)) +
    labs(title = paste0(CHROM, " — shelf diagnostic stack"),
         subtitle = sprintf("%d local-PCA windows%s%s", nrow(dt),
           if (is.finite(SHELF_A)) sprintf("  ·  shelf %g–%g Mb", SHELF_A, SHELF_B) else "",
           if (SMOOTH_WIN > 1) sprintf("  ·  smooth_win=%d", SMOOTH_WIN) else "")) +
    theme_void() +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, color = "#555e69"),
          plot.margin = margin(6, 46, 2, 10))
}

# =============================================================================
# TOP STRIP 2: sim_mat 1D collapse
# =============================================================================
# Near-diagonal similarity at each window = row mean of sim_mat values within
# +/- k off-diagonal. If sim_mat absent, skip.
sim_strip <- placeholder("")
sim_source <- pc$sim_mat
if (is.null(sim_source) && !is.null(SIMMAT_FILE) && file.exists(SIMMAT_FILE)) {
  sim_source <- tryCatch(readRDS(SIMMAT_FILE), error = function(e) NULL)
}
if (!is.null(sim_source) && is.matrix(sim_source) &&
    nrow(sim_source) == nrow(dt)) {
  N <- nrow(sim_source)
  # off-diagonal band width for collapse (in windows)
  K <- max(5L, floor(N * 0.01))
  collapse_val <- numeric(N)
  for (i in seq_len(N)) {
    lo <- max(1L, i - K); hi <- min(N, i + K)
    collapse_val[i] <- mean(sim_source[i, lo:hi], na.rm = TRUE)
  }
  sim_dt <- data.table(mb = dt$mb, v = collapse_val)
  sim_strip <- decorate(
    ggplot(sim_dt, aes(mb, v)) +
      geom_area(fill = "#cfa973", alpha = 0.6) +
      geom_line(color = "#8a6b3b", linewidth = 0.3)
  ) + common_x +
    labs(y = "sim_mat") +
    theme_classic(base_size = 8) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(1, 46, 1, 10),
          panel.grid = element_blank())
}

# =============================================================================
# Track 1: Z
# =============================================================================
p1 <- decorate(
  ggplot(dt, aes(mb, z)) +
    geom_hline(yintercept = c(-2.5, 0, 2.5), linetype = "dashed",
               color = c("#2c3e50", "#999999", "#2c3e50"), linewidth = 0.25) +
    geom_point(size = 0.22, color = "#1f4e79", alpha = 0.55)
) + common_x + labs(y = "Robust Z") + base_theme +
  edge_labels("outlier", "typical", "#2c3e50", "#999999")

# =============================================================================
# Track 2: n_snps
# =============================================================================
p2 <- decorate(
  ggplot(dt, aes(mb, n_snps)) +
    geom_point(size = 0.22, color = "#2c7a39", alpha = 0.55) +
    geom_hline(yintercept = median(dt$n_snps, na.rm = TRUE),
               linetype = "dashed", color = "#2c3e50", linewidth = 0.25)
) + common_x + labs(y = "SNPs/win") + base_theme +
  edge_labels("dense", "sparse")

# =============================================================================
# Track 3: BEAGLE uncertainty
# =============================================================================
p3 <- placeholder("uncertainty (no Q02)")
if (!is.null(UNC_FILE) && file.exists(UNC_FILE)) {
  unc <- fread(UNC_FILE)
  p3 <- decorate(
    ggplot(unc, aes(bin_mid_mb, mean_uncertain_frac)) +
      double_stroke("#ffffff", "#7b3294", 1.2, 0.45) +
      geom_hline(yintercept = median(unc$mean_uncertain_frac, na.rm = TRUE),
                 linetype = "dashed", color = "#2c3e50", linewidth = 0.25)
  ) + common_x + labs(y = "uncert frac") + base_theme +
    edge_labels("low-conf", "confident")
}

# =============================================================================
# Track 4: coverage mean + CV color
# =============================================================================
p4 <- placeholder("coverage (no Q03)")
if (!is.null(COV_FILE) && file.exists(COV_FILE)) {
  cov <- fread(COV_FILE)
  p4 <- decorate(
    ggplot(cov, aes(bin_mid_mb, mean_cov)) +
      geom_line(color = "#ffffff", linewidth = 1.2, alpha = 0.75) +
      geom_line(aes(color = cv_across_samples), linewidth = 0.4) +
      scale_color_gradient(low = "#1f4e79", high = "#c0504d", name = "CV") +
      geom_hline(yintercept = median(cov$mean_cov, na.rm = TRUE),
                 linetype = "dashed", color = "#2c3e50", linewidth = 0.25)
  ) + common_x + labs(y = "mean cov") + base_theme +
    theme(legend.position = "right",
          legend.key.width = unit(0.25, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))
}

# =============================================================================
# Track 5: # low-cov samples
# =============================================================================
p5 <- placeholder("# low-cov samples")
if (!is.null(COV_FILE) && file.exists(COV_FILE)) {
  cov2 <- fread(COV_FILE)
  p5 <- decorate(
    ggplot(cov2, aes(bin_mid_mb, n_samples_low_cov)) +
      geom_col(fill = "#c0504d", alpha = 0.7, width = 0.04)
  ) + common_x + labs(y = "# low cov") + base_theme +
    edge_labels("many low", "all OK")
}

# =============================================================================
# Track 6: theta mean + CV
# =============================================================================
p6 <- placeholder("theta mean / CV (no Q05)")
theta_summary <- NULL
if (!is.null(THETA_FILE) && file.exists(THETA_FILE)) {
  theta_raw <- fread(THETA_FILE)
  theta_raw <- theta_raw[chrom == CHROM]
  if (nrow(theta_raw) > 0) {
    theta_summary <- theta_raw[, .(
      mean_theta = mean(theta_pi_persite, na.rm = TRUE),
      cv_theta   = {
        x <- theta_pi_persite; m <- mean(x, na.rm = TRUE)
        if (is.na(m) || m == 0) NA_real_ else stats::sd(x, na.rm = TRUE) / m
      },
      n_samples  = .N
    ), by = win_center_bp]
    theta_summary[, win_mid_mb := win_center_bp / 1e6]
    setorder(theta_summary, win_mid_mb)
    theta_summary[, mean_theta_sm := rolling_median(mean_theta, SMOOTH_WIN)]
    theta_summary[, cv_theta_sm   := rolling_median(cv_theta,   SMOOTH_WIN)]
    rng_m <- range(theta_summary$mean_theta_sm, na.rm = TRUE)
    rng_c <- range(theta_summary$cv_theta_sm,   na.rm = TRUE)
    if (is.finite(rng_c[1]) && is.finite(rng_c[2]) && diff(rng_c) > 0) {
      theta_summary[, cv_overlay := rng_m[1] + (cv_theta_sm - rng_c[1]) /
                     diff(rng_c) * diff(rng_m)]
    } else {
      theta_summary[, cv_overlay := NA_real_]
    }
    p6 <- decorate(
      ggplot(theta_summary, aes(win_mid_mb)) +
        geom_line(aes(y = mean_theta_sm),
                  color = "#ffffff", linewidth = 1.2, alpha = 0.75) +
        geom_line(aes(y = mean_theta_sm), color = "#c0504d", linewidth = 0.45) +
        geom_line(aes(y = cv_overlay), color = "#7b3294", linewidth = 0.35,
                  linetype = "22")
    ) + common_x + labs(y = "θπ mean / CV") + base_theme +
      edge_labels("diverse", "low π")
  }
}

# =============================================================================
# Track 7: ancestry delta12 WITH regime-class strip above
# =============================================================================
p7 <- placeholder("ancestry Δ12 (no Q06)")
anc_dt <- NULL
if (!is.null(ANC_FILE) && file.exists(ANC_FILE)) {
  anc_dt <- fread(ANC_FILE)
  anc_dt <- anc_dt[chrom == CHROM]
  if (nrow(anc_dt) > 0) {
    setorder(anc_dt, window_mid_mb)
    anc_dt[, delta12_sm := rolling_median(delta12, SMOOTH_WIN)]

    # Regime class strip above the delta12 track
    # Color by dominant ancestry label at each window
    anc_dt[, maxQ_factor := factor(maxQ_label,
                                   levels = sort(unique(maxQ_label %||% "K?")))]
    # Palette
    lvls <- levels(anc_dt$maxQ_factor)
    strip_palette <- setNames(
      c("#1f4e79", "#f5a524", "#3cc08a", "#e0555c", "#b07cf7",
        "#f0d56a", "#7ad3db", "#ff8c6e")[seq_along(lvls)],
      lvls
    )
    # Top 8% of panel height is the regime strip
    y_max <- max(anc_dt$delta12_sm, na.rm = TRUE)
    y_min <- min(anc_dt$delta12_sm, na.rm = TRUE)
    pad <- (y_max - y_min) * 0.14
    strip_y_bot <- y_max + pad * 0.10
    strip_y_top <- y_max + pad

    p7 <- decorate(
      ggplot(anc_dt, aes(window_mid_mb, delta12_sm)) +
        # regime strip
        geom_rect(aes(xmin = window_start_bp / 1e6,
                      xmax = window_end_bp / 1e6,
                      ymin = strip_y_bot, ymax = strip_y_top,
                      fill = maxQ_factor),
                  color = NA) +
        scale_fill_manual(values = strip_palette, name = "dom. Q",
                          na.value = "#cccccc") +
        # delta12 double-stroke
        geom_line(color = "#ffffff", linewidth = 1.2, alpha = 0.75) +
        geom_line(color = "#2c7a39", linewidth = 0.45) +
        geom_hline(yintercept = c(0.5, 0.8), linetype = "dashed",
                   color = "#999999", linewidth = 0.25) +
        scale_y_continuous(limits = c(y_min, strip_y_top),
                           breaks = c(0, 0.5, 1.0))
    ) + common_x + labs(y = "Δ12") + base_theme +
      theme(legend.position = "right",
            legend.key.width = unit(0.25, "cm"),
            legend.key.height = unit(0.35, "cm"),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 7)) +
      edge_labels("clear", "ambiguous", "#2c7a39", "#2c3e50")
  }
}

# =============================================================================
# Track 7b: multi-scale Δ12
# =============================================================================
p7b <- NULL
if (!is.null(ANC_MULTI) && file.exists(ANC_MULTI)) {
  ms <- fread(ANC_MULTI)
  ms <- ms[chrom == CHROM]
  if (nrow(ms) > 0 && "scale" %in% names(ms)) {
    ms[, scale := factor(scale, levels = c("1x", "5x", "10x",
                                             setdiff(unique(scale),
                                                     c("1x", "5x", "10x"))))]
    setorder(ms, scale, window_mid_mb)
    if (SMOOTH_WIN > 1) {
      ms[, delta12_sm := rolling_median(delta12, SMOOTH_WIN), by = scale]
    } else {
      ms[, delta12_sm := delta12]
    }
    p7b <- decorate(
      ggplot(ms, aes(window_mid_mb, delta12_sm, color = scale)) +
        geom_line(linewidth = 0.4, alpha = 0.9) +
        scale_color_manual(values = c("1x" = "#2c7a39",
                                      "5x" = "#f5a524",
                                      "10x" = "#7b3294"),
                           name = "scale")
    ) + common_x + labs(y = "Δ12 (multi)") + base_theme +
      theme(legend.position = "right",
            legend.key.width = unit(0.3, "cm"),
            legend.key.height = unit(0.35, "cm"),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 7)) +
      edge_labels("scale-stable", "scale-dep")
  }
}

# =============================================================================
# Tracks 8 + 9: per-invgt θ_π and Hom1-vs-Hom2 Fst (from Q07)
# =============================================================================
p8 <- NULL; p9 <- NULL
if (!is.null(POPSTATS_INVGT_FILE) && file.exists(POPSTATS_INVGT_FILE)) {
  ps <- fread(POPSTATS_INVGT_FILE)
  setnames(ps, tolower(names(ps)))
  start_col <- intersect(c("window_start", "start_bp", "start", "ws"), names(ps))[1]
  end_col   <- intersect(c("window_end",   "end_bp",   "end",   "we"), names(ps))[1]
  if (!is.null(start_col) && !is.null(end_col)) {
    ps[, mb := (get(start_col) + get(end_col)) / 2 / 1e6]
  }
  theta_cols <- setdiff(grep("^theta_pi_", names(ps), value = TRUE),
                        "theta_pi_all")
  if (length(theta_cols) > 0 && "mb" %in% names(ps)) {
    long <- melt(ps[, c("mb", theta_cols), with = FALSE],
                 id.vars = "mb", variable.name = "group",
                 value.name = "theta_pi")
    long[, group := sub("^theta_pi_", "", as.character(group))]
    grp_order <- intersect(c("Hom1", "Het", "Hom2"), unique(long$group))
    grp_order <- c(grp_order, setdiff(unique(long$group), grp_order))
    long[, group := factor(group, levels = grp_order)]
    if (SMOOTH_WIN > 1) {
      long[, theta_pi_sm := rolling_median(theta_pi, SMOOTH_WIN), by = group]
    } else {
      long[, theta_pi_sm := theta_pi]
    }
    pal <- c("Hom1" = "#1f4e79", "Het" = "#c0504d", "Hom2" = "#2c7a39")
    p8 <- decorate(
      ggplot(long, aes(mb, theta_pi_sm, color = group)) +
        geom_line(linewidth = 0.45, alpha = 0.95) +
        scale_color_manual(values = pal, name = NULL)
    ) + common_x + labs(y = "θπ invgt") + base_theme +
      theme(legend.position = "right",
            legend.key.width = unit(0.3, "cm"),
            legend.key.height = unit(0.35, "cm"),
            legend.text = element_text(size = 7)) +
      edge_labels("diverse", "low π")
  }

  fst_cols <- grep("^(fst|hudson_fst)_", names(ps), value = TRUE)
  target <- fst_cols[grepl("hom1.*hom2|hom2.*hom1", fst_cols, ignore.case = TRUE)][1]
  if (is.na(target) && length(fst_cols) > 0) target <- fst_cols[1]
  if (!is.na(target) && "mb" %in% names(ps)) {
    d <- ps[, .(mb = mb, fst = get(target))]
    setorder(d, mb)
    if (SMOOTH_WIN > 1) d[, fst_sm := rolling_median(fst, SMOOTH_WIN)] else d[, fst_sm := fst]
    p9 <- decorate(
      ggplot(d, aes(mb, fst_sm)) +
        geom_line(color = "#ffffff", linewidth = 1.2, alpha = 0.75) +
        geom_line(color = "#7b3294", linewidth = 0.45) +
        geom_hline(yintercept = 0, color = "#888888", linewidth = 0.25) +
        geom_hline(yintercept = median(d$fst_sm, na.rm = TRUE),
                   linetype = "dashed", color = "#2c3e50", linewidth = 0.25)
    ) + common_x +
      labs(y = paste0("Fst ", sub("^(fst|hudson_fst)_", "", target)),
           x = paste0(CHROM, " (Mb)")) +
      bottom_theme +
      edge_labels("differentiated", "panmictic")
  }
}

# =============================================================================
# Compose
# =============================================================================
panels  <- list(ideogram, sim_strip, p1, p2, p3, p4, p5, p6, p7)
heights <- c(0.5,          0.5,       1.3, 0.9, 0.9, 0.9, 0.6, 0.9, 1.2)
if (!is.null(p7b)) { panels <- c(panels, list(p7b)); heights <- c(heights, 0.9) }
if (!is.null(p8))  { panels <- c(panels, list(p8));  heights <- c(heights, 0.9) }
if (!is.null(p9))  { panels <- c(panels, list(p9));  heights <- c(heights, 0.9) }

# Ensure only the last panel has x-axis labels. Rebuild bottom panel theme.
last_idx <- length(panels)
# Find which of p1..p9 was used last and swap its theme to bottom_theme
swap_bottom_theme <- function(p) {
  if (is.null(p)) return(NULL)
  p + theme_classic(base_size = 9) +
      theme(plot.margin = margin(2, 46, 5, 10),
            panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.3))
}
# We already applied bottom_theme to p9 above. Only re-swap if p9 isn't the last.
if (!is.null(panels[[last_idx]]) && !identical(panels[[last_idx]], p9)) {
  panels[[last_idx]] <- swap_bottom_theme(panels[[last_idx]]) +
    labs(x = paste0(CHROM, " (Mb)"))
}

plt <- Reduce(`/`, panels) + plot_layout(heights = heights)
fig_height <- 2 + sum(heights) * 1.3
ggsave(OUT, plt, width = 11, height = fig_height, device = cairo_pdf)

# =============================================================================
# Console summary
# =============================================================================
if (is.finite(SHELF_A) && is.finite(SHELF_B)) {
  mask_shelf <- dt$mb >= SHELF_A & dt$mb <= SHELF_B
  mask_ref   <- !mask_shelf & dt$mb > 2 & dt$mb < (max(dt$mb) - 2)
  cat("\n=========  ", CHROM, " SHELF SUMMARY  =========\n", sep = "")
  cat(sprintf("Shelf: %.2f-%.2f Mb (%d windows)   Ref: %d windows\n",
              SHELF_A, SHELF_B, sum(mask_shelf), sum(mask_ref)))
  fmt <- function(lbl, s, r) {
    cat(sprintf("  %-26s shelf=%9.4g   ref=%9.4g   ratio=%6.2fx\n",
                lbl, s, r, if (!is.finite(r) || abs(r) < 1e-12) NA else s / r))
  }
  fmt("Z mean",        mean(dt$z[mask_shelf], na.rm=TRUE), mean(dt$z[mask_ref], na.rm=TRUE))
  fmt("Z sd",          sd(dt$z[mask_shelf], na.rm=TRUE),   sd(dt$z[mask_ref],   na.rm=TRUE))
  fmt("n_snps mean",   mean(dt$n_snps[mask_shelf], na.rm=TRUE),
                       mean(dt$n_snps[mask_ref],   na.rm=TRUE))
  if ("lam_1" %in% names(dt)) {
    fmt("lambda1 mean", mean(dt$lam_1[mask_shelf], na.rm=TRUE),
                        mean(dt$lam_1[mask_ref],   na.rm=TRUE))
  }
  zf <- sd(dt$z[mask_shelf], na.rm=TRUE) / abs(mean(dt$z[mask_shelf], na.rm=TRUE))
  cat(sprintf("  Z flatness (shelf SD/|mean|): %.3f   <0.10 flat / >0.15 textured\n", zf))

  if (!is.null(UNC_FILE) && file.exists(UNC_FILE)) {
    unc <- fread(UNC_FILE)
    fmt("BEAGLE uncertain",
      mean(unc[bin_mid_mb >= SHELF_A & bin_mid_mb <= SHELF_B]$mean_uncertain_frac, na.rm=TRUE),
      mean(unc[bin_mid_mb <  SHELF_A | bin_mid_mb >  SHELF_B]$mean_uncertain_frac, na.rm=TRUE))
  }
  if (!is.null(COV_FILE) && file.exists(COV_FILE)) {
    cov <- fread(COV_FILE)
    fmt("coverage mean",
      mean(cov[bin_mid_mb >= SHELF_A & bin_mid_mb <= SHELF_B]$mean_cov, na.rm=TRUE),
      mean(cov[bin_mid_mb <  SHELF_A | bin_mid_mb >  SHELF_B]$mean_cov, na.rm=TRUE))
    fmt("coverage CV",
      mean(cov[bin_mid_mb >= SHELF_A & bin_mid_mb <= SHELF_B]$cv_across_samples, na.rm=TRUE),
      mean(cov[bin_mid_mb <  SHELF_A | bin_mid_mb >  SHELF_B]$cv_across_samples, na.rm=TRUE))
  }
  if (!is.null(theta_summary) && nrow(theta_summary) > 0) {
    fmt("theta_pi mean",
      mean(theta_summary[win_mid_mb >= SHELF_A & win_mid_mb <= SHELF_B]$mean_theta, na.rm=TRUE),
      mean(theta_summary[win_mid_mb <  SHELF_A | win_mid_mb >  SHELF_B]$mean_theta, na.rm=TRUE))
    fmt("theta_pi CV",
      mean(theta_summary[win_mid_mb >= SHELF_A & win_mid_mb <= SHELF_B]$cv_theta, na.rm=TRUE),
      mean(theta_summary[win_mid_mb <  SHELF_A | win_mid_mb >  SHELF_B]$cv_theta, na.rm=TRUE))
  }
  if (!is.null(anc_dt) && nrow(anc_dt) > 0) {
    fmt("ancestry Δ12",
      mean(anc_dt[window_mid_mb >= SHELF_A & window_mid_mb <= SHELF_B]$delta12, na.rm=TRUE),
      mean(anc_dt[window_mid_mb <  SHELF_A | window_mid_mb >  SHELF_B]$delta12, na.rm=TRUE))
    if ("cv_delta12_across_samples" %in% names(anc_dt)) {
      fmt("Δ12 CV(samples)",
        mean(anc_dt[window_mid_mb >= SHELF_A & window_mid_mb <= SHELF_B]$cv_delta12_across_samples, na.rm=TRUE),
        mean(anc_dt[window_mid_mb <  SHELF_A | window_mid_mb >  SHELF_B]$cv_delta12_across_samples, na.rm=TRUE))
    }
  }
  if (!is.null(POPSTATS_INVGT_FILE) && file.exists(POPSTATS_INVGT_FILE)) {
    ps <- fread(POPSTATS_INVGT_FILE)
    setnames(ps, tolower(names(ps)))
    scol <- intersect(c("window_start", "start_bp", "start", "ws"), names(ps))[1]
    ecol <- intersect(c("window_end",   "end_bp",   "end",   "we"), names(ps))[1]
    if (!is.null(scol) && !is.null(ecol)) ps[, mb := (get(scol) + get(ecol)) / 2 / 1e6]
    fst_cols <- grep("^(fst|hudson_fst)_", names(ps), value = TRUE)
    target <- fst_cols[grepl("hom1.*hom2|hom2.*hom1", fst_cols, ignore.case = TRUE)][1]
    if (!is.na(target) && "mb" %in% names(ps)) {
      fmt(paste0("Fst ", sub("^(fst|hudson_fst)_", "", target)),
        mean(ps[mb >= SHELF_A & mb <= SHELF_B][[target]], na.rm=TRUE),
        mean(ps[mb <  SHELF_A | mb >  SHELF_B][[target]], na.rm=TRUE))
    }
  }

  cat("\n  ---- reading the profile ----\n")
  cat("  flat Z + low n_snps + low coverage + low cov CV     -> mapping / repeat artifact\n")
  cat("  flat Z + normal n_snps + stable Δ12 + low theta CV  -> uniform substructure\n")
  cat("  flat Z + high Fst Hom1-Hom2 + θπ depression in homs -> REAL INVERSION\n")
  cat("  textured Z + stable Δ12 + high theta CV             -> REAL INVERSION (classic)\n")
  cat("============================================\n")
  cat(sprintf("Figure: %s\n", OUT))
}
