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
HOBS_FILE           <- get_arg("--hobs_track", NULL)
SIMMAT_FILE <- get_arg("--sim_mat", NULL)
OUT        <- get_arg("--out")
SHELF_A    <- as.numeric(get_arg("--shelf_start_mb", NA))
SHELF_B    <- as.numeric(get_arg("--shelf_end_mb",   NA))
BP1_MB     <- as.numeric(get_arg("--breakpoint1_mb", NA))
BP2_MB     <- as.numeric(get_arg("--breakpoint2_mb", NA))
SMOOTH_WIN <- as.integer(get_arg("--smooth_win", "1"))
# --no_nodata_strips disables the gray "no data here" background bands.
# Used to render the "clean presentation" variant alongside the diagnostic.
NO_STRIPS <- !is.null(get_arg("--no_nodata_strips", NULL))
# --ref_n_bed: optional BED file of reference N / assembly-gap regions for
# the current chromosome. When supplied, the ideogram gets a dark stipple
# over those regions. Distinct from lighter stipple used for zero-SNP
# dropouts detected from the precomp.
REF_N_BED <- get_arg("--ref_n_bed", NULL)

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

# =============================================================================
# No-data regions: contiguous stretches where no SNPs exist (reference gap,
# repeat-mask, or mappability dropout). Detected on first use from dt (the
# precomp windows). Uses n_snps == 0 as the primary signal and extreme
# span_kb as a secondary signal for windows that had to stretch to collect
# the fixed SNP count.
# =============================================================================
.NODATA_CACHE <- NULL
detect_nodata_regions <- function() {
  if (!is.null(.NODATA_CACHE)) return(.NODATA_CACHE)
  if (!exists("dt", envir = .GlobalEnv)) return(data.table())
  w <- copy(dt)
  w[, is_gap := FALSE]
  if ("n_snps"  %in% names(w)) w[n_snps == 0, is_gap := TRUE]
  if ("span_kb" %in% names(w)) {
    span_hi <- stats::quantile(w$span_kb, 0.995, na.rm = TRUE)
    w[span_kb > max(span_hi, 80), is_gap := TRUE]
  }
  if (!any(w$is_gap)) { .NODATA_CACHE <<- data.table(); return(.NODATA_CACHE) }
  g <- w[is_gap == TRUE][order(start_bp)]
  # Group contiguous gap windows (runs where next gap starts within 50 kb)
  g[, run_grp := cumsum(c(TRUE, diff(start_bp) > 50e3))]
  out <- g[, .(xmin = min(start_bp) / 1e6,
               xmax = max(end_bp)   / 1e6),
           by = run_grp]
  out[, run_grp := NULL]
  .NODATA_CACHE <<- out
  out
}

# Shelf shading + breakpoint lines + no-data strips applied to every panel
decorate <- function(p) {
  nd <- if (isTRUE(NO_STRIPS)) data.table() else detect_nodata_regions()
  # 1. No-data strips behind everything (light gray, behind the data)
  if (nrow(nd) > 0) {
    p <- p + geom_rect(data = nd,
                       aes(xmin = xmin, xmax = xmax,
                           ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE,
                       fill = "#888888", alpha = 0.18)
  }
  # 2. Shelf region shading (warm amber, slightly in front of no-data)
  if (is.finite(SHELF_A) && is.finite(SHELF_B)) {
    p <- p + annotate("rect", xmin = SHELF_A, xmax = SHELF_B,
                      ymin = -Inf, ymax = Inf,
                      fill = "#f5a524", alpha = 0.10)
  }
  # 3. Breakpoint vertical lines on top
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
# Compute span_kb early so no-data detection sees it before first decorate() call
dt[, span_kb := (end_bp - start_bp) / 1e3]
z_col <- NULL
for (cand in c("max_abs_z", "robust_z", "z_robust", "z", "mds_z_robust", "mds_z1_robust", "MDS1_z")) {
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
# Ideogram layers, drawn bottom-to-top:
#   1. chromosome body (light rounded bar)
#   2. zero-SNP stipple (light gray dots, 35% density) from no-data detection
#   3. reference-N stipple (darker dots, 60% density) if --ref_n_bed supplied
#   4. shelf highlight (red rectangle)
#   5. breakpoint vertical lines
#
# Stipple is drawn as a dot grid inside the region, so it reads as "50%-ink
# printer fill" which the eye parses as "structured absence" rather than
# "data of this color".
make_stipple <- function(xmin, xmax, y_lo = 0.05, y_hi = 0.95,
                         dx = 0.05, dy = 0.15, offset = 0) {
  if (!length(xmin) || !length(xmax)) return(data.frame(x = numeric(0), y = numeric(0)))
  dots <- list()
  for (i in seq_along(xmin)) {
    xs <- seq(xmin[i] + dx / 2 + offset * dx / 2, xmax[i] - dx / 2, by = dx)
    ys <- seq(y_lo + dy / 2, y_hi - dy / 2, by = dy)
    if (!length(xs) || !length(ys)) next
    # Offset alternate rows for staggered print-like pattern
    grid_rows <- list()
    for (j in seq_along(ys)) {
      x_row <- if (j %% 2 == 1) xs else xs + dx / 2
      x_row <- x_row[x_row <= xmax[i] - dx / 4]
      if (length(x_row)) grid_rows[[j]] <- data.frame(x = x_row, y = ys[j])
    }
    if (length(grid_rows)) dots[[i]] <- do.call(rbind, grid_rows)
  }
  if (!length(dots)) return(data.frame(x = numeric(0), y = numeric(0)))
  do.call(rbind, dots)
}

load_ref_n_regions <- function(bed_path, chrom) {
  if (is.null(bed_path) || !file.exists(bed_path)) return(data.table(xmin = numeric(0), xmax = numeric(0)))
  b <- tryCatch(fread(bed_path, header = FALSE, col.names = c("chrom", "start", "end")),
                error = function(e) NULL)
  if (is.null(b) || !nrow(b)) return(data.table(xmin = numeric(0), xmax = numeric(0)))
  b <- b[chrom == CHROM][, .(xmin = start / 1e6, xmax = end / 1e6)]
  b
}

ideogram <- {
  ideo_df <- data.frame(x0 = chrom_start_mb, x1 = chrom_end_mb,
                        y0 = 0, y1 = 1)
  p <- ggplot() +
    geom_rect(data = ideo_df, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
              fill = "#e8ebef", color = "#555e69", linewidth = 0.3) +
    geom_segment(data = data.frame(x = seq(ceiling(chrom_start_mb),
                                           floor(chrom_end_mb), by = 2)),
                 aes(x = x, xend = x, y = 0.02, yend = 0.98),
                 color = "#c8cdd2", linewidth = 0.25)

  # Layer 2: zero-SNP dropouts (light stipple)
  nd_regions <- detect_nodata_regions()
  if (nrow(nd_regions) > 0) {
    dot_dx <- max(0.03, (chrom_end_mb - chrom_start_mb) * 0.003)
    light_stipple <- make_stipple(nd_regions$xmin, nd_regions$xmax,
                                  dx = dot_dx, dy = 0.2)
    if (nrow(light_stipple) > 0) {
      p <- p + geom_point(data = light_stipple, aes(x = x, y = y),
                          color = "#8a92a0", size = 0.35, alpha = 0.7,
                          shape = 19)
    }
  }

  # Layer 3: reference-N / assembly-gap regions (dark stipple)
  ref_n <- load_ref_n_regions(REF_N_BED, CHROM)
  if (nrow(ref_n) > 0) {
    dot_dx <- max(0.02, (chrom_end_mb - chrom_start_mb) * 0.002)
    dark_stipple <- make_stipple(ref_n$xmin, ref_n$xmax,
                                 dx = dot_dx, dy = 0.15, offset = 1)
    if (nrow(dark_stipple) > 0) {
      p <- p + geom_point(data = dark_stipple, aes(x = x, y = y),
                          color = "#2c3e50", size = 0.5, alpha = 0.85,
                          shape = 19)
    }
  }

  # Layer 4: shelf highlight
  if (is.finite(SHELF_A) && is.finite(SHELF_B)) {
    p <- p + geom_rect(aes(xmin = SHELF_A, xmax = SHELF_B, ymin = 0, ymax = 1),
                       fill = "#e0555c", alpha = 0.55, color = "#8a2b30",
                       linewidth = 0.3)
  }
  # Layer 5: breakpoint ticks
  if (length(BPS) > 0) {
    p <- p + geom_segment(data = data.frame(x = BPS),
                          aes(x = x, xend = x, y = -0.2, yend = 1.2),
                          color = "#e0555c", linewidth = 0.5)
  }
  p + common_x +
    scale_y_continuous(limits = c(-0.3, 1.3), expand = c(0, 0)) +
    labs(title = paste0(CHROM, " — shelf diagnostic stack"),
         subtitle = sprintf("%d local-PCA windows%s%s%s",
           nrow(dt),
           if (is.finite(SHELF_A)) sprintf("  \u00b7  shelf %g\u2013%g Mb", SHELF_A, SHELF_B) else "",
           if (SMOOTH_WIN > 1) sprintf("  \u00b7  smooth_win=%d", SMOOTH_WIN) else "",
           if (nrow(nd_regions) + nrow(ref_n) > 0) sprintf("  \u00b7  %d dropout region%s",
             nrow(nd_regions) + nrow(ref_n),
             if (nrow(nd_regions) + nrow(ref_n) == 1) "" else "s") else "")) +
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
# Track 1: Z (robust, always non-negative for max_abs_z)
# =============================================================================
# max_abs_z is always >= 0. Fix y-axis to [0, 5] for cross-chromosome
# comparability. Dashed lines at Z=2.5 (caution) and Z=3.0 (standard detection
# threshold) give the reader an anchor.
z_detect_thresh <- as.numeric(get_arg("--z_detect_thresh", "3.0"))
p1 <- decorate(
  ggplot(dt, aes(mb, z)) +
    geom_hline(yintercept = 2.5, linetype = "dashed",
               color = "#999999", linewidth = 0.25) +
    geom_hline(yintercept = z_detect_thresh, linetype = "dashed",
               color = "#c0504d", linewidth = 0.35) +
    geom_point(size = 0.22, color = "#1f4e79", alpha = 0.55) +
    scale_y_continuous(limits = c(0, 5), expand = expansion(mult = 0.02))
) + common_x + labs(y = "Robust Z") + base_theme +
  edge_labels("outlier", "typical", "#2c3e50", "#999999")

# =============================================================================
# Track 2: SNP density (per kb / 10kb / 50kb, configurable)
# =============================================================================
# n_snps is fixed by window-size config (always 100 here), so instead we show
# true local density: SNPs per SNP_DENSITY_SCALE_KB kb. Default is per-10kb.
# Override via --snp_density_scale_kb (1, 10, 50 typical).
#
# Higher values = denser SNPs = healthier variant calling in that region.
# Drops would indicate repeat/mapping dropout.
sd_scale_kb <- as.numeric(get_arg("--snp_density_scale_kb", "10"))
# span_kb already computed at load time (line 192); just derive density here
dt[, snps_per_x := n_snps * sd_scale_kb / pmax(span_kb, 1e-9)]

if (SMOOTH_WIN > 1) {
  dt[, snps_per_x_sm := rolling_median(snps_per_x, SMOOTH_WIN)]
} else {
  dt[, snps_per_x_sm := snps_per_x]
}
sd_unit <- if (sd_scale_kb == 1) "SNPs/kb"
           else sprintf("SNPs/%gkb", sd_scale_kb)

p2 <- decorate(
  ggplot(dt, aes(mb, snps_per_x_sm)) +
    geom_line(color = "#ffffff", linewidth = 1.2, alpha = 0.75) +
    geom_line(color = "#2c7a39", linewidth = 0.45) +
    geom_hline(yintercept = median(dt$snps_per_x, na.rm = TRUE),
               linetype = "dashed", color = "#2c3e50", linewidth = 0.25)
) + common_x + labs(y = sd_unit) + base_theme +
  edge_labels("dense", "sparse")

# =============================================================================

# =============================================================================
# Track 3: BEAGLE uncertainty
# =============================================================================
p3 <- NULL
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
p4 <- NULL
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
p5 <- NULL
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
p6 <- NULL
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
p7 <- NULL
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
  # Drop windows with too few informative sites — Engine F writes exact zeros
  # there, which would otherwise appear as fake drop-to-baseline spikes.
  min_sites <- as.integer(Sys.getenv("Q04_MIN_SITES", "20"))
  n_used_col <- intersect(c("n_sites_used", "n_used", "nsites"), names(ps))[1]
  if (!is.null(n_used_col)) {
    ps <- ps[get(n_used_col) >= min_sites]
  }
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
# Track 10: Engine H — HoverE_Hom1 / HoverE_Het / HoverE_Hom2 per window
# (Merot sliding-window Hobs/Hexp, patched ANGSD HWE)
# At an inversion locus: HoverE_Het -> 2, HoverE_Hom1/Hom2 -> 0.
# Outside: all three -> ~1 (HWE).
# =============================================================================
p10 <- NULL
if (!is.null(HOBS_FILE) && file.exists(HOBS_FILE)) {
  hh <- fread(HOBS_FILE)
  setnames(hh, tolower(names(hh)))
  start_h <- intersect(c("start_bp", "start"), names(hh))[1]
  end_h   <- intersect(c("end_bp",   "end"),   names(hh))[1]
  if (!is.null(start_h) && !is.null(end_h)) {
    hh[, mb := (get(start_h) + get(end_h)) / 2 / 1e6]
    hov_cols <- grep("^hovere_", names(hh), value = TRUE)
    if (length(hov_cols) > 0) {
      long <- melt(hh[, c("mb", hov_cols), with = FALSE],
                   id.vars = "mb", variable.name = "group", value.name = "hovere")
      long[, group := sub("^hovere_", "", as.character(group))]
      long[, group := toupper(group)]
      grp_order <- intersect(c("HOM1", "HET", "HOM2"), unique(long$group))
      grp_order <- c(grp_order, setdiff(unique(long$group), grp_order))
      long[, group := factor(group, levels = grp_order)]
      if (SMOOTH_WIN > 1) {
        long[, hovere_sm := rolling_median(hovere, SMOOTH_WIN), by = group]
      } else {
        long[, hovere_sm := hovere]
      }
      pal <- c("HOM1" = "#1f4e79", "HET" = "#c0504d", "HOM2" = "#2c7a39")
      p10 <- decorate(
        ggplot(long, aes(mb, hovere_sm, color = group)) +
          geom_hline(yintercept = 1, linetype = "dashed",
                     color = "#2c3e50", linewidth = 0.3) +
          geom_hline(yintercept = 2, linetype = "dotted",
                     color = "#888888", linewidth = 0.25) +
          geom_line(linewidth = 0.45, alpha = 0.9) +
          scale_color_manual(values = pal, name = NULL) +
          scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2))
      ) + common_x + labs(y = "Hobs/Hexp") + base_theme +
        theme(legend.position = "right",
              legend.key.width = unit(0.3, "cm"),
              legend.key.height = unit(0.35, "cm"),
              legend.text = element_text(size = 7)) +
        edge_labels("Het excess (~2)", "Hom deficit (~0)")
    }
  }
}

# =============================================================================
# Compose (skip NULL panels for cleaner layout)
# =============================================================================
all_panels <- list(
  list(p = ideogram, h = 0.5),
  list(p = sim_strip, h = 0.5),
  list(p = p1,  h = 1.3),
  list(p = p2,  h = 0.9),
  list(p = p3,  h = 0.9),
  list(p = p4,  h = 0.9),
  list(p = p5,  h = 0.6),
  list(p = p6,  h = 0.9),
  list(p = p7,  h = 1.2),
  list(p = p7b, h = 0.9),
  list(p = p8,  h = 0.9),
  list(p = p9,  h = 0.9),
  list(p = p10, h = 0.9)
)
keep_idx <- which(sapply(all_panels, function(x) !is.null(x$p)))
panels  <- lapply(all_panels[keep_idx], `[[`, "p")
heights <- sapply(all_panels[keep_idx], `[[`, "h")

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
