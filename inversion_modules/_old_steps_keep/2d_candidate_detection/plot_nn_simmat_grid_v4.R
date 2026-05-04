#!/usr/bin/env Rscript
# =============================================================================
# plot_nn_simmat_grid.R  (v3 -- full factorial: smoothing x axis x palette)
# =============================================================================
#
# Standalone NN sim_mat visualization for the staircase sandbox.
# Produces grid composites spanning the full factorial of:
#
#   SMOOTHING (3)
#     mds_nn          k-NN smoothing in MDS space (the sim_mat_nn<k>.rds
#                      files produced by FIX 21 precomp). Long-range
#                      structurally-similar windows pulled together.
#     spatial_band    Band crop from raw nn0 sim_mat: show sim[i,j] only
#                      where abs(i-j) <= k windows. Near-diagonal view.
#                      Cells beyond the band are NA (not plotted / white).
#     spatial_smooth  For each window i, average sim_mat rows over the
#                      spatial neighborhood [i-k/2, i+k/2], symmetrised.
#                      Full matrix output. Physical-neighborhood smoothing.
#
#   AXIS (3)
#     window     window-index geometry, uniform tiles, no gaps
#     mb_gaps    Mb position axis with white gaps where windows are sparse
#     mb_packed  window-index geometry with Mb tick labels (old style)
#
#   PALETTE (4)
#     red_blue           blue -> white -> orange -> red (diverging-ish, sim)
#                         Fixed limits 0-1.
#     orange_green_peak  dark navy -> orange -> BRIGHT GREEN at peaks,
#                         for spotting near-identical blocks above the
#                         orange field (the "glowing" peak highlight).
#                         Fixed limits 0-1.
#     adaptive           CONTRAST-STRETCHED per-matrix: maps [p5, p99]
#                         of the observed sim distribution across the
#                         full color range, with green at the top 10%
#                         of OBSERVED values. Reveals internal block
#                         structure when sim is crammed in a narrow
#                         range (typical case), which washes out in
#                         fixed-scale palettes.
#     adaptive_peak_zoom Same as adaptive, plus a LOCAL ZOOM LENS on
#                         the top 5%: dark-green -> pink gradient
#                         spans [p95, p99.5] of observed sim. Reveals
#                         internal structure WITHIN the near-identical
#                         peak blocks -- the "core of the core".
#
# For each (smoothing, axis, palette) combination a single grid composite
# is rendered showing 6 NN scales (nn0, nn20, nn40, nn80, nn160, nn320).
# That is 3 x 3 x 4 = 36 grid PNGs per chromosome (~20 min on /mnt/e/).
#
# Use --individual-plots to also save every per-scale PNG (108 files).
# By default only grids are saved.
#
# LAYOUT EXPECTED
#   <root>/
#     01_input/
#       C_gar_LG28.precomp.slim.rds           # nn0 inside $sim_mat
#       C_gar_LG28.sim_mat_nn{20,40,80,160,320}.rds
#     04_output/
#
# USAGE
#   cd ~/staircase_local/
#   Rscript 03_plot/plot_nn_simmat_grid.R
#
# CLI FLAGS
#   --input-dir       folder containing the sim_mat RDS files    [01_input]
#   --output-dir      folder to write PNGs + PDFs                 [04_output]
#   --chr             chromosome label used in filenames/titles   [C_gar_LG28]
#   --scales          comma-separated NN scales to plot          [0,20,40,80,160,320]
#   --smoothings      comma-separated smoothing modes
#                       any of: mds_nn, spatial_band, spatial_smooth
#                                                                [all three]
#   --axes            comma-separated axis reps
#                       any of: window, mb_gaps, mb_packed        [all three]
#   --palettes        comma-separated palettes
#                       any of: red_blue, orange_green_peak, adaptive,
#                       adaptive_peak_zoom                        [all four]
#   --downsample      render at approx N x N tiles                [800]
#   --individual-plots   save every per-scale PNG in addition    [off]
#   --no-grid         skip grids, only individuals (if enabled)
#
# OUTPUT
#   <chr>_grid_<smoothing>_<axis>_<palette>.png
#   <chr>_grid_<smoothing>_<axis>_<palette>.pdf
#   <chr>_<smoothing>_<axis>_<palette>_nn<k>.png   (if --individual-plots)
#   <chr>_<smoothing>_<axis>_<palette>_nn<k>.pdf
#
# DEPENDENCIES
#   data.table, ggplot2, scales
#   patchwork (for grid; falls back to facet_wrap if absent)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
if (!HAS_PATCHWORK) {
  message("[plot] patchwork not installed -- grid will use facet_wrap.\n",
          "       install.packages('patchwork') for better grids.")
}


# =============================================================================
# CLI PARSING
# =============================================================================

parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  cli <- list(
    input_dir        = "01_input",
    output_dir       = "04_output",
    chr              = "C_gar_LG28",
    scales           = c(0L, 20L, 40L, 80L, 160L, 320L),
    smoothings       = c("mds_nn", "spatial_band", "spatial_smooth"),
    axes             = c("window", "mb_gaps", "mb_packed"),
    palettes         = c("red_blue", "orange_green_peak", "adaptive",
                         "adaptive_peak_zoom"),
    downsample       = 800L,
    individual_plots = FALSE,
    grid_plot        = TRUE
  )

  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    take_next <- function() { v <- args[i + 1L]; i <<- i + 2L; v }
    if      (a == "--input-dir"  && i < length(args)) cli$input_dir  <- take_next()
    else if (a == "--output-dir" && i < length(args)) cli$output_dir <- take_next()
    else if (a == "--chr"        && i < length(args)) cli$chr        <- take_next()
    else if (a == "--scales"     && i < length(args)) {
      cli$scales <- as.integer(strsplit(take_next(), ",")[[1]])
    }
    else if (a == "--smoothings" && i < length(args)) {
      cli$smoothings <- tolower(trimws(strsplit(take_next(), ",")[[1]]))
    }
    else if (a == "--axes"       && i < length(args)) {
      cli$axes <- tolower(trimws(strsplit(take_next(), ",")[[1]]))
    }
    else if (a == "--palettes"   && i < length(args)) {
      cli$palettes <- tolower(trimws(strsplit(take_next(), ",")[[1]]))
    }
    else if (a == "--downsample" && i < length(args)) {
      cli$downsample <- as.integer(take_next())
    }
    else if (a == "--individual-plots") { cli$individual_plots <- TRUE;  i <- i + 1L }
    else if (a == "--no-grid")          { cli$grid_plot        <- FALSE; i <- i + 1L }
    else { i <- i + 1L }
  }

  valid_smoothings <- c("mds_nn", "spatial_band", "spatial_smooth")
  valid_axes       <- c("window", "mb_gaps", "mb_packed")
  valid_palettes   <- c("red_blue", "orange_green_peak", "adaptive",
                        "adaptive_peak_zoom")

  bad_s <- setdiff(cli$smoothings, valid_smoothings)
  if (length(bad_s) > 0) stop("Unknown smoothing '", paste(bad_s, collapse = ","),
                               "'. Must be subset of: ",
                               paste(valid_smoothings, collapse = ", "))
  bad_a <- setdiff(cli$axes, valid_axes)
  if (length(bad_a) > 0) stop("Unknown axis '", paste(bad_a, collapse = ","),
                               "'. Must be subset of: ",
                               paste(valid_axes, collapse = ", "))
  bad_p <- setdiff(cli$palettes, valid_palettes)
  if (length(bad_p) > 0) stop("Unknown palette '", paste(bad_p, collapse = ","),
                               "'. Must be subset of: ",
                               paste(valid_palettes, collapse = ", "))
  cli
}

cli <- parse_cli()
dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== plot_nn_simmat_grid (v3) ===\n")
cat("  input dir   :", cli$input_dir,  "\n")
cat("  output dir  :", cli$output_dir, "\n")
cat("  chromosome  :", cli$chr,        "\n")
cat("  scales      :", paste(cli$scales,     collapse = ", "), "\n")
cat("  smoothings  :", paste(cli$smoothings, collapse = ", "), "\n")
cat("  axes        :", paste(cli$axes,       collapse = ", "), "\n")
cat("  palettes    :", paste(cli$palettes,   collapse = ", "), "\n")
cat("  downsample  :", cli$downsample, "\n")
cat("  individual  :", cli$individual_plots, "\n")
cat("  grid        :", cli$grid_plot,  "\n\n")


# =============================================================================
# STYLE -- theme with thin axis lines for a more professional look
# =============================================================================

THEME_HEATMAP <- theme_minimal(base_size = 9) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "#FAFAFA", color = NA),

    # Thin axis lines -- the professional polish you asked for
    axis.line         = element_line(colour = "#333333", linewidth = 0.3),
    axis.ticks        = element_line(colour = "#333333", linewidth = 0.25),
    axis.ticks.length = unit(2, "pt"),

    # Kill the grid -- it clutters a heatmap
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),

    plot.title        = element_text(size = 11, face = "bold"),
    plot.subtitle     = element_text(size = 8,  color = "grey40"),
    plot.caption      = element_text(size = 7,  color = "grey55", hjust = 0),
    axis.text         = element_text(size = 7),
    legend.key.width  = unit(0.3, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7)
  )


# =============================================================================
# PALETTES
# =============================================================================
#
# red_blue:          dark navy -> blue -> near-white -> orange -> red
#                     diverging-ish, median anchored, the old default
#
# orange_green_peak: dark navy -> orange field -> BRIGHT GREEN at peaks
#                     near-identical (sim ~ 0.9+) pops out in green against
#                     the orange background, for Z-order "glow" effect

palette_colours <- function(name, sim_vec = NULL) {
  # sim_vec: vector of observed sim values (finite, non-NA). Used by data-driven
  # palettes (e.g. "adaptive") for percentile-based anchor computation. Fixed
  # palettes (red_blue, orange_green_peak) ignore it and derive anchors from
  # the median alone.
  #
  # Returns a list with:
  #   cols         color stops
  #   anchors      numeric vector in [0, 1] of relative stop positions
  #   limits       c(lo, hi) for scale_fill_gradientn (NULL = auto-fit to data)
  #   nice_name    label for subtitle / filename
  med_sim <- if (is.null(sim_vec)) 0.5 else {
    m <- median(sim_vec, na.rm = TRUE); if (is.finite(m)) m else 0.5
  }

  switch(name,
    red_blue = list(
      cols      = c("#0C1E3C", "#2563EB", "#F8F9FA", "#E8913A", "#C0392B"),
      anchors   = pmax(0, pmin(1, c(0, med_sim*0.7, med_sim, med_sim*1.15, 1))),
      limits    = c(0, 1),
      nice_name = "red/blue"
    ),
    orange_green_peak = list(
      # dark-navy -> orange-ish field -> green peak -> deeper green
      cols      = c("#0C1E3C", "#3B5B90", "#E8913A", "#C0392B", "#7FD64A", "#1B6B1B"),
      # Anchor the orange-to-green transition at the 88-92% sim range
      # so only the "near-identical" blocks flip into green
      anchors   = pmax(0, pmin(1, c(0, med_sim*0.6, med_sim*1.1, 0.85, 0.92, 1))),
      limits    = c(0, 1),
      nice_name = "orange->green peak"
    ),
    adaptive = {
      # ---- Contrast-stretched palette ----------------------------------
      # The actual sim distribution often occupies only a narrow band of
      # [0, 1]. Mapping the full color range across [0, 1] wastes contrast
      # on empty regions. This palette clips to the observed [p5, p99]
      # range so the visible colors span the actual data, revealing
      # internal block structure that washes out in fixed-scale palettes.
      #
      # Within the clipped range:
      #   p5 -> p50     dark navy -> blue -> near-white  (background)
      #   p50 -> p90    warm orange field                (mid-structure)
      #   p90 -> p99    bright green peak                (near-identical)
      # The "green at top 10%" is data-driven: if sim maxes out at 0.7,
      # the green band hits around 0.67. If it reaches 1.0, green kicks
      # in around 0.95. Same visual semantics, different absolute values.
      if (is.null(sim_vec) || length(sim_vec) < 10) {
        p5  <- 0;  p50 <- 0.5;  p90 <- 0.85;  p99 <- 1
      } else {
        qs  <- quantile(sim_vec, c(0.05, 0.50, 0.90, 0.99), na.rm = TRUE)
        p5  <- qs[1]; p50 <- qs[2]; p90 <- qs[3]; p99 <- qs[4]
        # Guard against pathological distributions where percentiles collide
        if (p99 - p5 < 0.01) { p5 <- p5 - 0.01; p99 <- p99 + 0.01 }
      }
      list(
        cols      = c("#0C1E3C", "#3B5B90", "#F5E6C5",  # navy -> steel -> cream
                      "#E8913A", "#C0392B",              # orange -> red
                      "#7FD64A", "#1B6B1B"),             # bright -> deep green
        anchors   = scales::rescale(c(p5, (p5 + p50) / 2, p50,
                                       (p50 + p90) / 2, p90,
                                       (p90 + p99) / 2, p99)),
        limits    = c(p5, p99),
        nice_name = sprintf("adaptive [%.2f-%.2f]", p5, p99)
      )
    },
    adaptive_peak_zoom = {
      # ---- Adaptive with a "zoom lens" on the top ---------------------
      # Extends `adaptive` by stretching the top ~5% of observed sim
      # values across their own dark-green-to-pink gradient. Within the
      # already-green peak region, you can now see a SECOND internal
      # gradient: cells at p95 read as dark green, cells at p99.5
      # render as hot pink. Equivalent to taking the range [p95, p99.5]
      # and rescaling it to 0-1 against a dedicated 2-stop palette,
      # while keeping the full p5-p95 gradient below unchanged.
      #
      # Use this to find the "core of the core" -- which cells within
      # a near-identical block are the MOST identical. Often marks the
      # tightest LD center of an inversion, or the most conserved span
      # of a near-fixed haplotype.
      if (is.null(sim_vec) || length(sim_vec) < 10) {
        p5  <- 0;   p50 <- 0.5; p90 <- 0.80; p95 <- 0.90
        p97 <- 0.95; p99 <- 0.98; p99_5 <- 1
      } else {
        qs <- quantile(sim_vec,
                       c(0.05, 0.50, 0.90, 0.95, 0.97, 0.99, 0.995),
                       na.rm = TRUE)
        p5 <- qs[1]; p50 <- qs[2]; p90 <- qs[3]
        p95 <- qs[4]; p97 <- qs[5]; p99 <- qs[6]; p99_5 <- qs[7]
        if (p99_5 - p5 < 0.01) { p5 <- p5 - 0.01; p99_5 <- p99_5 + 0.01 }
        # Ensure monotonicity in case of ties
        pts <- sort(c(p5, p50, p90, p95, p97, p99, p99_5))
        p5 <- pts[1]; p50 <- pts[2]; p90 <- pts[3]
        p95 <- pts[4]; p97 <- pts[5]; p99 <- pts[6]; p99_5 <- pts[7]
      }
      # Palette: first 5 stops span p5..p95 (the "normal" part).
      # Last 4 stops span p95..p99.5 (the "zoom"): dark green -> pink.
      # Because the top range (p95..p99.5) is tiny in sim units, the
      # zoom gradient visually blows it up to a big fraction of the
      # colorbar. Cells at p96 show dark green, p98 show reddish-pink,
      # p99+ show hot pink.
      list(
        cols    = c(
          "#0C1E3C",  # p5    dark navy
          "#3B5B90",  # p50-? steel blue
          "#F5E6C5",  # p50   cream
          "#E8913A",  # p90   warm orange
          "#1B6B1B",  # p95   dark green (zoom starts)
          "#5EB848",  # p97   bright green
          "#FF69B4",  # p99   pink
          "#FF1493"   # p99.5 hot pink
        ),
        anchors = scales::rescale(c(p5, (p5 + p50) / 2, p50, p90,
                                     p95, p97, p99, p99_5)),
        limits  = c(p5, p99_5),
        nice_name = sprintf("adaptive+zoom p95=%.2f p99.5=%.2f",
                            p95, p99_5)
      )
    },
    stop("Unknown palette: ", name)
  )
}


# =============================================================================
# DATA LOADING
# =============================================================================

load_precomp <- function(input_dir, chr) {
  f <- file.path(input_dir, paste0(chr, ".precomp.slim.rds"))
  if (!file.exists(f)) {
    f_alt <- file.path(input_dir, paste0(chr, ".precomp.rds"))
    if (file.exists(f_alt)) f <- f_alt
    else stop("precomp RDS not found in ", input_dir)
  }
  message("[plot] Loading precomp: ", f)
  pc <- readRDS(f)
  if (is.null(pc$sim_mat)) stop("precomp RDS is missing $sim_mat slot")
  pc
}

load_nn_simmat <- function(input_dir, chr, k) {
  candidates <- c(
    file.path(input_dir, sprintf("%s.sim_mat_nn%d.rds", chr, k)),
    file.path(input_dir, sprintf("%s_sim_mat_nn%d.rds", chr, k)),
    file.path(input_dir, "sim_mats", sprintf("%s.sim_mat_nn%d.rds", chr, k)),
    file.path(input_dir, "sim_mats", sprintf("%s_sim_mat_nn%d.rds", chr, k))
  )
  f <- NA_character_
  for (p in candidates) if (file.exists(p)) { f <- p; break }
  if (is.na(f)) return(NULL)
  obj <- tryCatch(readRDS(f), error = function(e) {
    message("[plot nn", k, "] read failed: ", conditionMessage(e)); NULL
  })
  if (is.null(obj)) return(NULL)
  if (is.list(obj) && !is.null(obj$sim_mat)) return(obj$sim_mat)
  if (is.matrix(obj)) return(obj)
  message("[plot nn", k, "] unrecognised RDS schema at ", f); NULL
}


# =============================================================================
# SMOOTHING TRANSFORMS applied to the raw nn0 sim_mat
# =============================================================================
#
# spatial_band:   mask sim_mat[i,j] to NA where abs(i-j) > k
# spatial_smooth: for each row i, replace with mean of rows [i-k/2, i+k/2]
#                 (symmetrised to keep matrix symmetric)

spatial_band_crop <- function(sim_mat_raw, k) {
  if (k == 0L) return(sim_mat_raw)
  n <- nrow(sim_mat_raw)
  out <- sim_mat_raw
  # Vectorised band mask: for each offset d > k, zero that diagonal
  for (d in (k + 1L):(n - 1L)) {
    # Super-diagonal d and sub-diagonal -d
    rows1 <- seq_len(n - d)
    cols1 <- rows1 + d
    out[cbind(rows1, cols1)] <- NA_real_
    out[cbind(cols1, rows1)] <- NA_real_
  }
  out
}

spatial_smooth <- function(sim_mat_raw, k) {
  if (k == 0L) return(sim_mat_raw)
  n <- nrow(sim_mat_raw)
  half <- as.integer(floor(k / 2))
  # For each row i, average rows [i-half, i+half] of the raw sim_mat.
  # Then symmetrise so the output is a valid sim_mat.
  sm <- matrix(0, n, n)
  for (i in seq_len(n)) {
    i_lo <- max(1L, i - half)
    i_hi <- min(n, i + half)
    sm[i, ] <- colMeans(sim_mat_raw[i_lo:i_hi, , drop = FALSE], na.rm = TRUE)
  }
  sm <- (sm + t(sm)) / 2
  diag(sm) <- 1
  sm
}


# =============================================================================
# AXIS VALUE COMPUTATION
# =============================================================================

compute_axis_values <- function(dt, n, axis_kind) {
  have_bp <- !is.null(dt) &&
             all(c("start_bp", "end_bp") %in% names(dt)) &&
             nrow(dt) >= n

  if (axis_kind == "window" || !have_bp) {
    return(list(x = seq_len(n), label_unit = "window index", is_mb = FALSE))
  }

  mb_center <- (dt$start_bp[seq_len(n)] + dt$end_bp[seq_len(n)]) / 2e6

  if (axis_kind == "mb_gaps") {
    return(list(x = mb_center, label_unit = "Mb", is_mb = TRUE))
  }

  if (axis_kind == "mb_packed") {
    # Window-index geometry (uniform tiles, no gaps) but with Mb tick labels.
    # We return x = window index, plus a breaks/labels override to be used
    # by scale_x_continuous.
    return(list(
      x           = seq_len(n),
      label_unit  = "Mb (packed)",
      is_mb       = FALSE,
      mb_labels   = mb_center,
      is_packed   = TRUE
    ))
  }

  stop("Unknown axis_kind: ", axis_kind)
}


# =============================================================================
# HEATMAP BUILDER
# =============================================================================

build_heatmap_tile_data <- function(sim_mat, dt, nds, axis_kind) {
  n <- nrow(sim_mat)
  step <- max(1L, n %/% nds)
  idx  <- seq(1L, n, by = step)
  ns   <- length(idx)

  ax <- compute_axis_values(dt, n, axis_kind)
  x_vals <- ax$x[idx]

  grid <- CJ(ii = seq_len(ns), jj = seq_len(ns))
  grid[, x   := x_vals[ii]]
  grid[, y   := x_vals[jj]]
  grid[, sim := sim_mat[cbind(idx[ii], idx[jj])]]
  grid <- grid[is.finite(sim), .(x, y, sim)]

  # For mb_packed, also pass through the mb labels for axis breaks
  mb_labels <- if (isTRUE(ax$is_packed)) ax$mb_labels[idx] else NULL

  list(tiles = grid, ns = ns, step = step,
       axis_label = ax$label_unit, mb_labels = mb_labels, idx = idx)
}


build_single_heatmap <- function(sim_mat, dt, chr, k,
                                   nds = 800L, axis_kind = "window",
                                   palette_name = "red_blue",
                                   title_prefix = "") {
  if (is.null(sim_mat) || !is.matrix(sim_mat) || nrow(sim_mat) < 20) return(NULL)

  td <- build_heatmap_tile_data(sim_mat, dt, nds, axis_kind)
  if (nrow(td$tiles) == 0) return(NULL)

  med_sim <- median(td$tiles$sim, na.rm = TRUE)
  if (!is.finite(med_sim)) med_sim <- 0.5

  # Pass the observed sim values so data-driven palettes (adaptive) can
  # compute percentile anchors. Fixed palettes ignore sim_vec.
  pal <- palette_colours(palette_name, sim_vec = td$tiles$sim)
  anchors   <- sort(unique(pal$anchors))
  cols_trim <- pal$cols[seq_along(anchors)]
  if (length(cols_trim) < 2L) {
    # Guard against degenerate cases (all sim values identical)
    cols_trim <- rep(pal$cols[1], 2)
    anchors   <- c(0, 1)
  }

  title_tag <- if (k == 0L) "nn0 (raw)" else paste0("nn", k)
  x_label <- paste0(chr, " (", td$axis_label, ")")

  p <- ggplot(td$tiles, aes(x = x, y = y, fill = sim)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = cols_trim,
      values  = scales::rescale(anchors),
      limits  = pal$limits,
      oob     = scales::squish,  # clamp out-of-limits values to scale range
      na.value = "white",        # white for out-of-band / missing cells
      name    = "Similarity"
    ) +
    coord_fixed(expand = FALSE) +
    labs(
      x        = x_label,
      y        = x_label,
      title    = paste0(chr, " -- ", title_prefix, title_tag,
                        "  [", toupper(axis_kind), " | ", pal$nice_name, "]"),
      subtitle = sprintf("N = %d windows | %d x %d tiles | median sim = %.3f",
                         nrow(sim_mat), td$ns, td$ns, med_sim)
    ) +
    THEME_HEATMAP

  # For mb_packed, rewrite x/y axis breaks to show Mb labels at regular intervals
  if (!is.null(td$mb_labels)) {
    n_ticks <- 6L
    tick_idx <- round(seq(1, length(td$mb_labels), length.out = n_ticks))
    breaks   <- td$mb_labels[tick_idx]  # these are Mb values
    # But our x axis is window index, not Mb. We need to map back:
    tick_positions <- td$idx[tick_idx]  # window indices matching those breaks
    labels_txt     <- sprintf("%.0f", td$mb_labels[tick_idx])
    p <- p +
      scale_x_continuous(breaks = tick_positions, labels = labels_txt,
                         expand = c(0, 0)) +
      scale_y_continuous(breaks = tick_positions, labels = labels_txt,
                         expand = c(0, 0))
  }

  p
}


# =============================================================================
# EXECUTION -- load data once, iterate over all (smoothing, axis, palette)
# =============================================================================

pc <- load_precomp(cli$input_dir, cli$chr)
sim_mat_raw <- pc$sim_mat
dt_precomp  <- pc$dt
n_windows   <- nrow(sim_mat_raw)

cat(sprintf("[plot] Raw sim_mat: %d x %d windows\n", n_windows, ncol(sim_mat_raw)))
cat(sprintf("       dt rows    : %d\n", nrow(dt_precomp %||% data.table())))

# --- Resolve matrices for each (smoothing, scale) combo. Load/compute lazily
#     and cache so we don't redo work. ---
get_matrix <- function(smoothing, k, cache) {
  key <- paste0(smoothing, "_", k)
  if (!is.null(cache[[key]])) return(list(mat = cache[[key]], cache = cache))

  if (k == 0L) {
    cache[[key]] <- sim_mat_raw
    return(list(mat = sim_mat_raw, cache = cache))
  }

  mat <- switch(smoothing,
    mds_nn         = load_nn_simmat(cli$input_dir, cli$chr, k),
    spatial_band   = spatial_band_crop(sim_mat_raw, k),
    spatial_smooth = spatial_smooth(sim_mat_raw, k),
    NULL
  )

  if (is.null(mat)) return(list(mat = NULL, cache = cache))

  cache[[key]] <- mat
  list(mat = mat, cache = cache)
}

matrix_cache <- list()

# Precompute spatial variants once (fast, avoid redoing per palette/axis)
cat("\n[plot] Pre-computing spatial variants (one-time) ...\n")
for (sm_name in intersect(cli$smoothings, c("spatial_band", "spatial_smooth"))) {
  for (k in cli$scales) {
    if (k == 0L) next
    t0 <- Sys.time()
    res <- get_matrix(sm_name, k, matrix_cache)
    matrix_cache <- res$cache
    if (!is.null(res$mat)) {
      cat(sprintf("  %-15s nn%-4d  built in %4.1f s\n",
                  sm_name, k,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
}

# Load MDS NN files up front (they're on disk, cheap to preload)
if ("mds_nn" %in% cli$smoothings) {
  cat("\n[plot] Loading mds_nn sim_mats from disk ...\n")
  for (k in cli$scales) {
    t0 <- Sys.time()
    res <- get_matrix("mds_nn", k, matrix_cache)
    matrix_cache <- res$cache
    status <- if (is.null(res$mat)) "NOT FOUND" else
      sprintf("%4.1f s (%d x %d)",
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              nrow(res$mat), ncol(res$mat))
    cat(sprintf("  mds_nn          nn%-4d  %s\n", k, status))
  }
}


# --- Main loop: render grid for every (smoothing, axis, palette) combo ----

combo_count <- length(cli$smoothings) * length(cli$axes) * length(cli$palettes)
cat(sprintf("\n[plot] Rendering %d grid combinations ...\n\n", combo_count))

combo_i <- 0L
for (sm_name in cli$smoothings) {
  for (axis_kind in cli$axes) {
    for (palette_name in cli$palettes) {
      combo_i <- combo_i + 1L
      combo_tag <- paste0(sm_name, "_", axis_kind, "_", palette_name)

      cat(sprintf("[%d/%d] %s\n", combo_i, combo_count, combo_tag))
      t_combo <- Sys.time()

      plots_for_grid <- list()

      for (k in cli$scales) {
        res <- get_matrix(sm_name, k, matrix_cache)
        matrix_cache <- res$cache
        mat <- res$mat
        if (is.null(mat)) {
          cat(sprintf("      nn%-4d  no matrix  SKIPPED\n", k))
          next
        }

        t0 <- Sys.time()
        p <- tryCatch(
          build_single_heatmap(mat, dt_precomp, cli$chr, k,
                               nds = cli$downsample,
                               axis_kind = axis_kind,
                               palette_name = palette_name),
          error = function(e) {
            message(sprintf("  [%s nn%d] build failed: %s",
                            combo_tag, k, conditionMessage(e))); NULL
          }
        )
        if (is.null(p)) next
        plots_for_grid[[as.character(k)]] <- p
        dt_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

        if (cli$individual_plots) {
          out_png <- file.path(cli$output_dir,
            sprintf("%s_%s_nn%d.png", cli$chr, combo_tag, k))
          out_pdf <- file.path(cli$output_dir,
            sprintf("%s_%s_nn%d.pdf", cli$chr, combo_tag, k))
          tryCatch({
            ggsave(out_png, p, width = 9, height = 9, dpi = 200)
            ggsave(out_pdf, p, width = 9, height = 9)
            cat(sprintf("      nn%-4d  build+save %.1f s -> %s\n",
                        k, dt_s, basename(out_png)))
          }, error = function(e) {
            message(sprintf("  [%s nn%d] save failed: %s",
                            combo_tag, k, conditionMessage(e)))
          })
        } else {
          cat(sprintf("      nn%-4d  built in %.1f s\n", k, dt_s))
        }
      }

      # ---- Grid composite ----
      if (cli$grid_plot && length(plots_for_grid) >= 2L) {
        t_grid <- Sys.time()

        # Strip per-panel titles, add lightweight sub-titles, kill legends except last
        plots_stripped <- lapply(names(plots_for_grid), function(nm) {
          k <- as.integer(nm)
          sub_title <- if (k == 0L) "nn0 (raw)" else paste0("nn", k)
          plots_for_grid[[nm]] +
            labs(title = sub_title, subtitle = NULL, caption = NULL) +
            theme(plot.title = element_text(size = 10, face = "bold",
                                             hjust = 0.5),
                  plot.subtitle = element_blank(),
                  plot.caption  = element_blank())
        })
        names(plots_stripped) <- names(plots_for_grid)

        n_plots <- length(plots_stripped)
        for (i in seq_len(n_plots - 1L)) {
          plots_stripped[[i]] <- plots_stripped[[i]] +
            theme(legend.position = "none")
        }

        if (HAS_PATCHWORK) {
          n_cols <- if (n_plots <= 3L) n_plots
                    else if (n_plots <= 6L) 3L
                    else 4L
          composite <- Reduce(`+`, plots_stripped) +
            patchwork::plot_layout(ncol = n_cols, guides = "collect") +
            patchwork::plot_annotation(
              title = paste0(cli$chr, " -- ", combo_tag),
              subtitle = sprintf(
                "smoothing = %s | axis = %s | palette = %s",
                sm_name, axis_kind, palette_name),
              theme = theme(
                plot.title    = element_text(size = 13, face = "bold"),
                plot.subtitle = element_text(size = 10, color = "grey40")
              )
            )

          out_png <- file.path(cli$output_dir,
            sprintf("%s_grid_%s.png", cli$chr, combo_tag))
          out_pdf <- file.path(cli$output_dir,
            sprintf("%s_grid_%s.pdf", cli$chr, combo_tag))
          total_w <- 9 * n_cols
          total_h <- 9 * ceiling(n_plots / n_cols)
          tryCatch({
            ggsave(out_png, composite, width = total_w, height = total_h,
                   dpi = 180, limitsize = FALSE)
            ggsave(out_pdf, composite, width = total_w, height = total_h,
                   limitsize = FALSE)
            grid_s <- as.numeric(difftime(Sys.time(), t_grid, units = "secs"))
            cat(sprintf("      [GRID] %dx%d  save=%.1fs -> %s\n",
                        total_w, total_h, grid_s, basename(out_png)))
          }, error = function(e) {
            message(sprintf("  [%s grid] save failed: %s",
                            combo_tag, conditionMessage(e)))
          })

        } else {
          # patchwork-less fallback
          fallback_list <- list()
          for (nm in names(plots_stripped)) {
            pd <- plots_stripped[[nm]]$data
            if (is.null(pd)) next
            fallback_list[[length(fallback_list) + 1L]] <- cbind(
              pd, data.table(scale = paste0("nn", as.integer(nm)))
            )
          }
          if (length(fallback_list) > 0) {
            fb_dt <- rbindlist(fallback_list)
            pal <- palette_colours(palette_name, sim_vec = fb_dt$sim)
            anchors <- sort(unique(pal$anchors))
            cols_trim <- pal$cols[seq_along(anchors)]
            if (length(cols_trim) < 2L) {
              cols_trim <- rep(pal$cols[1], 2); anchors <- c(0, 1)
            }
            p_fb <- ggplot(fb_dt, aes(x = x, y = y, fill = sim)) +
              geom_raster() +
              scale_fill_gradientn(
                colours = cols_trim, values = scales::rescale(anchors),
                limits  = pal$limits, oob = scales::squish,
                na.value = "white",
                name = "Similarity") +
              coord_fixed(expand = FALSE) +
              facet_wrap(~ scale, ncol = 3) +
              labs(
                title = paste0(cli$chr, " -- ", combo_tag),
                subtitle = "(patchwork not installed; faceted fallback)",
                x = NULL, y = NULL
              ) + THEME_HEATMAP
            out_png <- file.path(cli$output_dir,
              sprintf("%s_grid_%s.png", cli$chr, combo_tag))
            ggsave(out_png, p_fb, width = 18, height = 18, dpi = 180,
                   limitsize = FALSE)
            cat(sprintf("      [GRID facet] %s\n", basename(out_png)))
          }
        }
      }

      combo_s <- as.numeric(difftime(Sys.time(), t_combo, units = "secs"))
      cat(sprintf("      combo done in %.1f s\n\n", combo_s))
    }
  }
}


# =============================================================================
# SUMMARY
# =============================================================================

cat("============================================================\n")
cat("[DONE] Outputs in ", cli$output_dir, ":\n", sep = "")
out_files <- list.files(cli$output_dir,
                          pattern = sprintf("^%s_", cli$chr),
                          full.names = FALSE)
for (f in out_files) {
  fp <- file.path(cli$output_dir, f)
  cat(sprintf("  %8.1f KB  %s\n", file.size(fp) / 1024, f))
}
cat("============================================================\n")

cat("\n=== INTERPRETATION GUIDE ===\n\n")
cat("THE 3 SMOOTHING MODES\n")
cat("  mds_nn          each window averaged with its k most similar\n")
cat("                   windows ANYWHERE on chromosome (MDS-space).\n")
cat("                   Long-range correlated windows pulled together.\n")
cat("  spatial_band    raw sim_mat cropped to abs(i-j) <= k.\n")
cat("                   Off-diagonal cells white. Shows local structure.\n")
cat("  spatial_smooth  each row averaged over spatial neighborhood +/-k/2.\n")
cat("                   Full matrix, spatial-only smoothing.\n\n")
cat("THE 3 AXIS REPRESENTATIONS\n")
cat("  window     window-index grid, uniform tiles, no gaps. Native\n")
cat("              PCA/MDS geometry. Block that is square here is\n")
cat("              genuinely square in information space.\n")
cat("  mb_gaps    Mb position axis with WHITE gaps where windows are\n")
cat("              sparse. Shows true physical distortion.\n")
cat("  mb_packed  window-index geometry (no gaps) + Mb tick labels.\n")
cat("              Uniform tile density, Mb labels for reference.\n\n")
cat("THE 4 PALETTES\n")
cat("  red_blue            blue -> white -> orange -> red, median-anchored,\n")
cat("                      fixed scale 0-1. Standard view.\n")
cat("  orange_green_peak   dark-navy -> orange field -> BRIGHT GREEN at\n")
cat("                      sim ~0.9+. Fixed scale 0-1. Only fires when\n")
cat("                      absolute similarity reaches 0.85+.\n")
cat("  adaptive            CONTRAST-STRETCHED per-matrix: maps [p5, p99]\n")
cat("                      of the observed sim distribution across the\n")
cat("                      full color range. Green at the top 10%% of\n")
cat("                      OBSERVED values (data-driven, not absolute).\n")
cat("                      Reveals internal block structure when sim is\n")
cat("                      crammed in a narrow range (e.g. 0.5-0.75),\n")
cat("                      which washes out in fixed-scale palettes.\n")
cat("  adaptive_peak_zoom  Same as adaptive, PLUS a local zoom lens on\n")
cat("                      the top 5%%. Within the dark-green peak\n")
cat("                      region, a second dark-green -> PINK gradient\n")
cat("                      reveals internal structure. Cells at p95 read\n")
cat("                      as dark green, cells at p99.5 render as hot\n")
cat("                      pink. Use to find the 'core of the core' --\n")
cat("                      the most conserved cells WITHIN a block.\n\n")
cat("COMPARE ACROSS GRIDS\n")
cat("  Block square in WINDOW but rectangle in MB_GAPS -> density artifact.\n")
cat("  Block visible in MDS_NN but not SPATIAL_BAND -> long-range correlation\n")
cat("    (ancestry, translocation, etc.) -- NOT a local inversion.\n")
cat("  Block visible in ALL smoothings -> robust local inversion signal.\n")
cat("============================================================\n")
