#!/usr/bin/env Rscript
# =============================================================================
# STEP_D17c_overlay_plot_L1_only_v6.R
#
# v6 changes vs v5:
#   - Boundary square geometry corrected: the lower-right CORNER VERTEX of
#     each boundary square now sits exactly on the diagonal at the point
#     (i+0.5, i+0.5), and the box extends strictly into the upper triangle
#     (j > i). The diagonal touches the square at a single corner point —
#     no pixel of the box straddles the diagonal. This matches the
#     geometric meaning of "boundary at window i" for an upper-triangle
#     cross-block: rows = (i+1, ..., i+W), cols = (i-W+1, ..., i).
#   - When the new growing-W validator is in use (v6 multipass with
#     --boundary_validator_mode=grow), the validation_status column drives
#     filtering and coloring exactly as before — no change needed here.
#   - Boundary label expanded to also show grow_max_z when present.
#
# Whole-genome single-panel overlay for L1 envelope inspection. One PDF page,
# one big square showing the entire chromosome's sim_mat with envelope
# outlines drawn on top.
#
# Lower triangle: raw similarity (yellow/orange/green palette).
# Upper triangle: diagonal-distance Z (blue/white/red palette).
# Outlines:       L1 envelopes from STEP_D17_multipass_L1_only.R catalogue
#                 (status == "ENVELOPE"). Drawn as deep-blue outlined boxes
#                 with ggrepel labels showing candidate_id and Mb range.
#
# Inputs:
#   --precomp     <slim precomp .rds>
#   --sim_mat     <smoothed sim_mat .rds, default tries sim_mat_nn80.rds>
#   --catalogue   <L1 envelope TSV from STEP_D17_multipass_L1_only.R>
#   --chr         <chromosome label>
#   --outdir      <output directory>
#   --nn_label    <nn80 / nn160 / etc. — for plot title only>
#   --page_size   <inches; default 14 — single square page>
#   --outline_lw  <outline line width; default 0.6>
#   --label_size  <ggrepel label size; default 2.4>
#   --z_clip      <clip Z to +/- this; default 5>
#
# Output:
#   <outdir>/<chr>_d17L1_overlay.pdf  (single-page PDF, one panel)
#
# Author: Claude (Anthropic)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(scales)
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ---- CLI --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

precomp_f    <- get_arg("--precomp")
sim_mat_f    <- get_arg("--sim_mat")
catalogue_f  <- get_arg("--catalogue")
boundaries_f <- get_arg("--boundaries")  # optional sidecar TSV
chr_label    <- get_arg("--chr", "chr")
outdir       <- get_arg("--outdir", ".")
nn_label     <- get_arg("--nn_label", "nn80")
page_size    <- as.numeric(get_arg("--page_size", "14"))
outline_lw   <- as.numeric(get_arg("--outline_lw", "0.6"))
label_size   <- as.numeric(get_arg("--label_size", "1.2"))
z_clip       <- as.numeric(get_arg("--z_clip", "5"))
boundary_lw    <- as.numeric(get_arg("--boundary_lw", "0.4"))

# Boundary status filter & colors. Available statuses:
#   STABLE_BLUE  -- real long-range separation, persists at all val_W
#   DECAYS       -- trigger fired but val_score collapses as W grows
#   MARGINAL     -- ambiguous
#   EDGE         -- too close to chr edge to validate
# Filter values:
#   "stable"  draws only STABLE_BLUE       (recommended default)
#   "all"     draws every status (each in its own color)
#   "non_decay" draws everything except DECAYS
boundary_filter <- get_arg("--boundary_filter", "stable")
boundary_filter <- tolower(boundary_filter)
status_palette <- c(
  STABLE_BLUE = "#C8102E",   # red — real boundary
  DECAYS      = "#999999",   # gray — fake / inside an inversion
  MARGINAL    = "#E07B3F",   # orange — ambiguous
  EDGE        = "#7A4FBF"    # purple — chromosome edge
)

if (is.na(precomp_f)   || !file.exists(precomp_f))   stop("[D17L1c] --precomp required")
if (is.na(catalogue_f) || !file.exists(catalogue_f)) stop("[D17L1c] --catalogue required")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[D17L1c] precomp:   ", precomp_f, "\n")
cat("[D17L1c] sim_mat:   ", sim_mat_f %||% "(auto)", "\n")
cat("[D17L1c] catalogue: ", catalogue_f, "\n")
cat("[D17L1c] boundaries:", boundaries_f %||% "(none)", "\n")
cat("[D17L1c] boundary_filter: ", boundary_filter, "\n", sep = "")
cat("[D17L1c] chr:       ", chr_label, "\n")
cat("[D17L1c] outdir:    ", outdir, "\n")
cat("[D17L1c] page_size: ", page_size, "in\n")

# ---- Load precomp -----------------------------------------------------------

cat("[D17L1c] loading precomp\n")
pc_obj <- readRDS(precomp_f)
if (!is.null(pc_obj$dt)) {
  pc <- pc_obj
} else if (!is.null(pc_obj$pc)) {
  pc <- pc_obj$pc
} else stop("[D17L1c] cannot find pc$dt in precomp file")
dt_pc <- as.data.table(pc$dt)
n_windows_total <- nrow(dt_pc)

# Coordinate columns
if (all(c("start", "end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start
  window_end_bp   <- dt_pc$end
} else if (all(c("start_bp", "end_bp") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start_bp
  window_end_bp   <- dt_pc$end_bp
} else if (all(c("window_start", "window_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$window_start
  window_end_bp   <- dt_pc$window_end
} else if (all(c("bp_start", "bp_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$bp_start
  window_end_bp   <- dt_pc$bp_end
} else {
  stop("[D17L1c] precomp dt missing usable coordinate columns. Found: ",
       paste(names(dt_pc), collapse = ", "))
}
cat("[D17L1c] N windows: ", n_windows_total, "\n")

# ---- Load sim_mat -----------------------------------------------------------

if (is.na(sim_mat_f)) {
  candidates_path <- c(
    file.path(dirname(precomp_f), "sim_mat_nn80.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn160.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn40.rds")
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0) stop("[D17L1c] could not auto-find sim_mat_nn*.rds")
  sim_mat_f <- hit[1]
}
cat("[D17L1c] loading sim_mat: ", sim_mat_f, "\n")
sm_obj <- readRDS(sim_mat_f)
if (is.matrix(sm_obj)) {
  sim_mat <- sm_obj
} else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat)) {
  sim_mat <- sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) == 1) {
  sim_mat <- sm_obj[[1]]
} else stop("[D17L1c] sim_mat structure not recognized")
storage.mode(sim_mat) <- "double"

if (!isTRUE(nrow(sim_mat) == n_windows_total &&
            ncol(sim_mat) == n_windows_total)) {
  stop("[D17L1c] sim_mat dim ", nrow(sim_mat), "x", ncol(sim_mat),
       " does not match N windows ", n_windows_total)
}

# ---- Load catalogue ---------------------------------------------------------

cat("[D17L1c] loading catalogue\n")
cat_dt <- fread(catalogue_f)
cat("[D17L1c] catalogue rows: ", nrow(cat_dt), "\n", sep = "")

# Filter to envelope rows only
cat_dt <- cat_dt[status == "ENVELOPE"]
cat_dt <- cat_dt[!is.na(start_w) & !is.na(end_w)]
cat("[D17L1c] envelopes to draw: ", nrow(cat_dt), "\n", sep = "")

# ---- Load boundaries (optional) --------------------------------------------
bdt <- data.table()
if (!is.na(boundaries_f) && file.exists(boundaries_f)) {
  bdt <- fread(boundaries_f)
  cat("[D17L1c] boundaries to draw: ", nrow(bdt), "\n", sep = "")
} else if (!is.na(boundaries_f)) {
  cat("[D17L1c] boundaries file not found: ", boundaries_f,
      " (skipping)\n", sep = "")
}

# ---- Diagonal-distance Z ----------------------------------------------------

cat("[D17L1c] computing diagonal-distance Z\n")
N <- n_windows_total
z_mat <- matrix(0, nrow = N, ncol = N)
for (d in 0:(N - 1L)) {
  idx_i <- seq_len(N - d)
  idx_j <- idx_i + d
  idx_upper <- cbind(idx_i, idx_j)
  idx_lower <- cbind(idx_j, idx_i)
  vals_upper <- sim_mat[idx_upper]
  vals_lower <- sim_mat[idx_lower]
  vals <- if (d == 0L) vals_upper else c(vals_upper, vals_lower)
  mu_d <- mean(vals, na.rm = TRUE)
  sd_d <- sd(vals, na.rm = TRUE)
  if (!is.finite(sd_d) || sd_d < 1e-9) sd_d <- 1
  z_mat[idx_upper] <- (vals_upper - mu_d) / sd_d
  if (d > 0L) z_mat[idx_lower] <- (vals_lower - mu_d) / sd_d
  if (d > 0L && d %% 1000L == 0L) {
    cat(sprintf("[D17L1c]   diag d=%d / %d\n", d, N - 1L))
  }
}
zr <- range(z_mat, na.rm = TRUE)
cat(sprintf("[D17L1c] Z range: [%.2f, %.2f]; clipping at +/- %.1f\n",
            zr[1], zr[2], z_clip))
z_mat[z_mat >  z_clip] <-  z_clip
z_mat[z_mat < -z_clip] <- -z_clip

# ---- Build heatmap long-form (whole genome) --------------------------------

cat("[D17L1c] building long-form heatmap (",
    format(N * N, big.mark = ","), " cells)\n", sep = "")

# Note: for very large N this can be RAM-heavy. With N=4302, ~18.5M rows.
ii <- rep(seq_len(N), times = N)
jj <- rep(seq_len(N), each = N)
lower_mask <- jj <= ii
upper_mask <- jj >  ii

hm <- data.table(
  i       = ii,
  j       = jj,
  sim_val = ifelse(lower_mask, as.numeric(sim_mat), NA_real_),
  z_val   = ifelse(upper_mask, as.numeric(z_mat),   NA_real_)
)

sim_layer <- hm[!is.na(sim_val)]
z_layer   <- hm[!is.na(z_val)]

q_lo <- as.numeric(quantile(sim_layer$sim_val, 0.05, na.rm = TRUE))
q_hi <- as.numeric(quantile(sim_layer$sim_val, 0.98, na.rm = TRUE))

z_max <- max(abs(z_layer$z_val), na.rm = TRUE)
z_max <- min(z_max, z_clip)
if (!is.finite(z_max) || z_max < 1e-3) z_max <- 1

# ---- Envelope rectangles + labels -------------------------------------------

if (nrow(cat_dt) > 0L) {
  cat_dt[, draw_start := pmax(start_w, 1L)]
  cat_dt[, draw_end   := pmin(end_w, n_windows_total)]
  cat_dt[, label_text := sprintf("%s  %.2f-%.2f Mb (nW=%d)",
                                 candidate_id,
                                 start_bp / 1e6, end_bp / 1e6,
                                 n_windows)]
  # Place labels above each envelope's upper edge
  cat_dt[, lab_x := (draw_start + draw_end) / 2]
  cat_dt[, lab_y := draw_end + (n_windows_total * 0.015)]
}

# ---- Axis breaks (whole genome) --------------------------------------------

axis_breaks <- pretty(c(1L, n_windows_total), n = 8L)
axis_breaks <- unique(as.integer(round(axis_breaks)))
axis_breaks <- axis_breaks[axis_breaks >= 1L & axis_breaks <= n_windows_total]
axis_labels_mb <- sprintf("%.1f",
  window_start_bp[pmin(pmax(axis_breaks, 1L), n_windows_total)] / 1e6
)

# ---- Plot -------------------------------------------------------------------

cat("[D17L1c] building plot\n")

p <- ggplot() +
  # Lower triangle: similarity
  geom_raster(
    data = sim_layer,
    aes(x = i, y = j, fill = sim_val)
  ) +
  scale_fill_gradientn(
    colours = c("#F8F8F8", "#A8DBC2", "#F2DC78", "#E08838", "#7E1F1F"),
    values  = scales::rescale(c(0, q_lo, (q_lo + q_hi) / 2, q_hi, 1)),
    limits  = c(0, 1),
    oob     = scales::squish,
    name    = "Similarity\n(lower △)"
  ) +
  ggnewscale::new_scale_fill() +
  # Upper triangle: Z
  geom_raster(
    data = z_layer,
    aes(x = i, y = j, fill = z_val)
  ) +
  scale_fill_gradient2(
    low      = "#2C5AA0",
    mid      = "#FAFAFA",
    high     = "#B22222",
    midpoint = 0,
    limits   = c(-z_max, z_max),
    oob      = scales::squish,
    name     = "Z (diagonal)\n(upper △)"
  )

# Envelope outlines + labels
if (nrow(cat_dt) > 0L) {
  p <- p +
    geom_rect(
      data = cat_dt,
      aes(xmin = draw_start - 0.5, xmax = draw_end + 0.5,
          ymin = draw_start - 0.5, ymax = draw_end + 0.5),
      colour    = "#1F3A6E",   # deep blue
      fill      = NA,
      linewidth = outline_lw,
      linetype  = "solid",
      inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = cat_dt,
      aes(x = lab_x, y = lab_y, label = label_text),
      colour = "#1F3A6E",
      size = label_size,
      fontface = "bold",
      bg.colour = "white",
      bg.r = 0.10,
      segment.size = outline_lw,
      segment.alpha = 0.5,
      min.segment.length = 0,
      seed = 1,
      max.overlaps = Inf,
      force = 1.5,
      box.padding = 0.30,
      point.padding = 0.10,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
}

# Boundary squares (small WxW upper-triangle markers + labels)
if (nrow(bdt) > 0L) {
  # Apply status filter if validation columns are present
  if ("validation_status" %in% names(bdt)) {
    n_before <- nrow(bdt)
    keep_set <- switch(boundary_filter,
      "stable"    = c("STABLE_BLUE"),
      "non_decay" = c("STABLE_BLUE", "MARGINAL", "EDGE"),
      "all"       = c("STABLE_BLUE", "DECAYS", "MARGINAL", "EDGE"),
      c("STABLE_BLUE")
    )
    bdt <- bdt[validation_status %in% keep_set]
    cat("[D17L1c] boundary_filter='", boundary_filter,
        "': ", n_before, " -> ", nrow(bdt), " boundaries\n", sep = "")
  } else {
    cat("[D17L1c] no validation_status column — drawing all boundaries\n")
    bdt[, validation_status := "STABLE_BLUE"]   # synthetic for coloring
  }
}

if (nrow(bdt) > 0L) {
  # Trigger W (drawn square size — fixed, not tied to validation result;
  # validation is now 1D ray-based, not 2D growing-square).
  bW_trig <- if ("boundary_W" %in% names(bdt))
               as.integer(bdt$boundary_W[1]) else 5L
  bG <- if ("boundary_offset" %in% names(bdt))
          as.integer(bdt$boundary_offset[1]) else 5L

  brect <- copy(bdt)
  brect[, draw_W := bW_trig]

  # v6 square geometry: the lower-right CORNER VERTEX sits exactly on the
  # diagonal at point (i + 0.5, i + 0.5), and the box extends W cells up
  # (+y direction, columns) and W cells left (-x direction, rows). With
  # plot mapping x = row, y = column, the upper triangle is y > x.
  #
  # Pixel cells covered by the box:
  #   x-range (rows)    : i - W + 1 ... i
  #   y-range (columns) : i + 1     ... i + W
  # All cells satisfy y > x, so the box lies strictly in the upper triangle
  # and the diagonal only touches it at the single corner point (i+0.5, i+0.5).
  brect[, sq_xmin := boundary_w - draw_W + 0.5]
  brect[, sq_xmax := boundary_w + 0.5]
  brect[, sq_ymin := boundary_w + 0.5]
  brect[, sq_ymax := boundary_w + draw_W + 0.5]
  # Clip to plot area
  brect <- brect[sq_xmax <= n_windows_total + 0.5 &
                 sq_ymax <= n_windows_total + 0.5 &
                 sq_xmin >= 0.5 & sq_ymin >= 0.5]
  # Label sits just above the upper-right corner of the square
  brect[, lab_x_b := (sq_xmin + sq_xmax) / 2]
  brect[, lab_y_b := sq_ymax + (n_windows_total * 0.012)]
  # Label content: id + Mb + trigger score + grow_max_z + perp-ray stats + IN-L1
  has_perp   <- "right_frac_blue" %in% names(brect)
  has_inside <- "inside_L1"       %in% names(brect)
  has_grow   <- "grow_max_z"      %in% names(brect)
  brect[, label_text_b := {
    base <- sprintf("%s  %.2f Mb  s=%.1f",
                    boundary_idx, boundary_bp / 1e6, boundary_score)
    if (has_grow) base <- sprintf("%s  gmz=%+.2f", base, grow_max_z)
    if (has_perp) base <- sprintf("%s  rfb=%.2f lfb=%.2f rmz=%+.1f lmz=%+.1f",
                                  base, right_frac_blue, left_frac_blue,
                                  right_max_z, left_max_z)
    if (has_inside)  base <- ifelse(inside_L1,
                                    sprintf("%s  IN-L1", base),
                                    base)
    base
  }]
  # Match status_palette colors
  brect[, status_color := unname(status_palette[validation_status])]

  p <- p +
    geom_rect(
      data = brect,
      aes(xmin = sq_xmin, xmax = sq_xmax,
          ymin = sq_ymin, ymax = sq_ymax,
          colour = validation_status),
      fill      = NA,
      linewidth = boundary_lw,
      linetype  = "solid",
      inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = brect,
      aes(x = lab_x_b, y = lab_y_b, label = label_text_b,
          colour = validation_status),
      size = label_size,
      fontface = "plain",
      bg.colour = "white",
      bg.r = 0.10,
      segment.size = boundary_lw * 0.6,
      segment.alpha = 0.5,
      min.segment.length = 0,
      seed = 2,
      max.overlaps = Inf,
      force = 0.8,
      box.padding = 0.20,
      point.padding = 0.05,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = status_palette,
      breaks = c("STABLE_BLUE", "MARGINAL", "EDGE", "DECAYS"),
      name   = "Boundary status",
      drop   = TRUE
    )
}

p <- p +
  coord_equal() +
  scale_x_continuous(
    expand = c(0, 0),
    name = "window index",
    breaks = axis_breaks,
    sec.axis = sec_axis(
      ~ .,
      breaks = axis_breaks,
      labels = axis_labels_mb,
      name = "Mb"
    )
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    name = "window index",
    breaks = axis_breaks,
    sec.axis = sec_axis(
      ~ .,
      breaks = axis_breaks,
      labels = axis_labels_mb,
      name = "Mb"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.2, "cm"),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(title = sprintf(
    "%s | %s | whole genome | L1 envelopes (%d)%s",
    chr_label, nn_label, nrow(cat_dt),
    if (nrow(bdt) > 0L)
      sprintf(" | boundaries (%d)", nrow(bdt))
    else ""))

# ---- Write -----------------------------------------------------------------

outpath <- file.path(outdir, paste0(chr_label, "_d17L1_overlay.pdf"))
cat("[D17L1c] writing PDF: ", outpath, "\n")
pdf(outpath, width = page_size, height = page_size)
print(p)
dev.off()
cat("[D17L1c] done\n")
