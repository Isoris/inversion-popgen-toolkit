#!/usr/bin/env Rscript

# =============================================================================
# theme_systems_plate.R
#
# WHITE MODULAR SYSTEMS FIGURE THEME
# "Publication Systems Plate Theme"
#
# Reusable styling system for all pipeline figures.
# Source this file from any R script:
#   source("theme_systems_plate.R")
#
# Provides:
#   - theme_plate()           : base ggplot theme
#   - pal_*                   : 6 palette variants
#   - plate_scatter()         : styled scatter helper
#   - plate_track()           : chromosome profile helper
#   - plate_box()             : boxplot/violin helper
#   - plate_heatmap()         : heatmap helper
#   - plate_table()           : embedded table grob
#   - plate_callout()         : annotation/summary box
#   - plate_composite()       : patchwork assembly helper
#
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# =============================================================================
# 1. BASE THEME
# =============================================================================

theme_plate <- function(base_size = 9, base_family = "sans",
                         grid_x = FALSE, grid_y = TRUE) {

  t <- theme(
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#FAFAFA", color = NA),

    # Grid
    panel.grid.major.y = if (grid_y) element_line(color = "#E8E8E8", linewidth = 0.3)
                          else element_blank(),
    panel.grid.major.x = if (grid_x) element_line(color = "#E8E8E8", linewidth = 0.3)
                          else element_blank(),
    panel.grid.minor = element_blank(),

    # Axes
    axis.line = element_line(color = "#4A4A4A", linewidth = 0.35),
    axis.ticks = element_line(color = "#4A4A4A", linewidth = 0.25),
    axis.ticks.length = unit(0.12, "cm"),

    # Typography
    plot.title = element_text(
      size = base_size + 3, face = "bold", color = "#1A1A1A",
      margin = margin(b = 4), hjust = 0, family = base_family),
    plot.subtitle = element_text(
      size = base_size, color = "#5A5A5A",
      margin = margin(b = 8), hjust = 0, family = base_family),
    plot.caption = element_text(
      size = base_size - 2, color = "#8A8A8A", hjust = 0,
      margin = margin(t = 6), family = base_family, lineheight = 1.2),
    axis.title = element_text(
      size = base_size, color = "#3A3A3A", family = base_family),
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(
      size = base_size - 1, color = "#4A4A4A", family = base_family),

    # Facet strips
    strip.background = element_rect(fill = "#F0F0F0", color = "#D0D0D0", linewidth = 0.3),
    strip.text = element_text(
      size = base_size, face = "bold", color = "#2A2A2A",
      margin = margin(t = 3, b = 3), family = base_family),

    # Legend
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = base_size - 1, color = "#4A4A4A"),
    legend.title = element_text(size = base_size, face = "bold", color = "#3A3A3A"),
    legend.position = "bottom",
    legend.margin = margin(t = 2),

    # Spacing
    plot.margin = margin(t = 8, r = 8, b = 6, l = 6),
    panel.spacing = unit(0.6, "lines")
  )
  t
}

# =============================================================================
# 2. COLOR PALETTES (6 variants, same structure)
# =============================================================================

# Metric color families (constant across all palettes)
COL_DELTA12  <- "#2563EB"   # blue family — dominance
COL_ENTROPY  <- "#7C3AED"   # purple family — mixing
COL_ENA      <- "#059669"   # green family — effective number
COL_HET      <- "#D97706"   # amber — heterozygosity
COL_FST      <- "#DC2626"   # red — differentiation
COL_THETA    <- "#0891B2"   # cyan — diversity

# Class colors (inversion genotypes)
COL_HOM_STD  <- "#3B82F6"   # blue
COL_HET_INV  <- "#22C55E"   # green
COL_HOM_INV  <- "#EF4444"   # red
COL_RECOM    <- "#F59E0B"   # amber
COL_UNKNOWN  <- "#9CA3AF"   # gray

CLASS_COLORS <- c(
  "HOM_STD" = COL_HOM_STD, "HOMO_REF" = COL_HOM_STD,
  "HET" = COL_HET_INV,
  "HOM_INV" = COL_HOM_INV, "HOMO_INV" = COL_HOM_INV,
  "RARE_INV" = "#A855F7", "RARE_HET" = "#84CC16",
  "recombinant" = COL_RECOM, "UNKNOWN" = COL_UNKNOWN
)

BAND_COLORS <- c("1" = COL_HOM_STD, "2" = COL_HET_INV, "3" = COL_HOM_INV)

# 6 accent palette variants
pal_editorial_blue <- list(
  accent = "#2563EB", accent2 = "#1D4ED8", accent3 = "#93C5FD",
  neutral = "#64748B", bg_accent = "#EFF6FF", border = "#BFDBFE",
  name = "Editorial Blue"
)

pal_deep_indigo <- list(
  accent = "#4F46E5", accent2 = "#3730A3", accent3 = "#A5B4FC",
  neutral = "#6B7280", bg_accent = "#EEF2FF", border = "#C7D2FE",
  name = "Deep Indigo"
)

pal_slate_teal <- list(
  accent = "#0D9488", accent2 = "#115E59", accent3 = "#99F6E4",
  neutral = "#64748B", bg_accent = "#F0FDFA", border = "#99F6E4",
  name = "Slate Teal"
)

pal_steel_lavender <- list(
  accent = "#7C3AED", accent2 = "#5B21B6", accent3 = "#C4B5FD",
  neutral = "#6B7280", bg_accent = "#F5F3FF", border = "#DDD6FE",
  name = "Steel Lavender"
)

pal_warm_technical <- list(
  accent = "#D97706", accent2 = "#92400E", accent3 = "#FDE68A",
  neutral = "#78716C", bg_accent = "#FFFBEB", border = "#FDE68A",
  name = "Warm Technical"
)

pal_neutral_scientific <- list(
  accent = "#475569", accent2 = "#1E293B", accent3 = "#CBD5E1",
  neutral = "#64748B", bg_accent = "#F8FAFC", border = "#E2E8F0",
  name = "Neutral Scientific"
)

# Default palette
PAL <- pal_editorial_blue

set_palette <- function(pal) {
  PAL <<- pal
  message("[theme] Palette: ", pal$name)
}

# =============================================================================
# 3. SCALE HELPERS
# =============================================================================

# Continuous similarity scale (navy → gold → red, centered on median)
scale_fill_simmat <- function(median_val = 0.5, ...) {
  scale_fill_gradientn(
    colours = c("#0C1E3C", "#1E3A5F", "#4A7FB5", "#F0C75E", "#E8913A", "#C0392B"),
    values = scales::rescale(c(0, median_val * 0.5, median_val * 0.85,
                                median_val * 1.1, median_val * 1.3, 1)),
    name = "Similarity", ...)
}

# Diverging scale (for contrast matrices)
scale_fill_diverge <- function(name = "Value", ...) {
  scale_fill_gradient2(low = "#C0392B", mid = "white", high = "#2563EB",
                        midpoint = 0, name = name, ...)
}

# =============================================================================
# 4. HELPER: TRACK PLOT (chromosome profiles)
# =============================================================================

plate_track <- function(data, x, y, color = PAL$accent,
                         title = NULL, subtitle = NULL, caption = NULL,
                         y_lab = NULL, x_lab = NULL,
                         hline = NULL, hline_color = "#9CA3AF",
                         region_start = NULL, region_end = NULL,
                         point_size = 0.3, point_alpha = 0.45) {

  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    {if (!is.null(region_start) && !is.null(region_end))
      annotate("rect", xmin = region_start, xmax = region_end,
               ymin = -Inf, ymax = Inf, fill = PAL$bg_accent, alpha = 0.6)} +
    geom_point(size = point_size, alpha = point_alpha, color = color) +
    {if (!is.null(hline))
      geom_hline(yintercept = hline, linetype = "dashed",
                 color = hline_color, linewidth = 0.3)} +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = x_lab, y = y_lab) +
    theme_plate(grid_x = FALSE, grid_y = TRUE)
  p
}

# =============================================================================
# 5. HELPER: BOXPLOT / VIOLIN
# =============================================================================

plate_box <- function(data, x, y, fill = NULL, fill_values = CLASS_COLORS,
                       title = NULL, subtitle = NULL, caption = NULL,
                       y_lab = NULL, x_lab = NULL, legend_title = NULL) {

  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]]))

  if (!is.null(fill)) {
    p <- p + aes(fill = .data[[fill]]) +
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.75,
                    width = 0.7, color = "#4A4A4A") +
      scale_fill_manual(values = fill_values, name = legend_title)
  } else {
    p <- p + geom_boxplot(outlier.size = 0.3, linewidth = 0.3,
                           fill = PAL$accent3, color = PAL$accent, alpha = 0.7, width = 0.6)
  }

  p + labs(title = title, subtitle = subtitle, caption = caption,
           x = x_lab, y = y_lab) +
    theme_plate(grid_x = FALSE, grid_y = TRUE)
}

# =============================================================================
# 6. HELPER: SCATTER
# =============================================================================

plate_scatter <- function(data, x, y, color_by = NULL, color_values = NULL,
                           title = NULL, subtitle = NULL, caption = NULL,
                           x_lab = NULL, y_lab = NULL,
                           size = 0.8, alpha = 0.6) {

  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]]))

  if (!is.null(color_by)) {
    p <- p + aes(color = .data[[color_by]]) +
      geom_point(size = size, alpha = alpha) +
      {if (!is.null(color_values)) scale_color_manual(values = color_values)
       else scale_color_viridis_c(option = "C")}
  } else {
    p <- p + geom_point(size = size, alpha = alpha, color = PAL$accent)
  }

  p + labs(title = title, subtitle = subtitle, caption = caption,
           x = x_lab, y = y_lab) +
    theme_plate(grid_x = TRUE, grid_y = TRUE) +
    coord_fixed()
}

# =============================================================================
# 7. HELPER: HEATMAP (ggplot-based)
# =============================================================================

plate_heatmap <- function(data, x, y, fill, fill_scale = NULL,
                           title = NULL, subtitle = NULL, caption = NULL,
                           x_lab = NULL, y_lab = NULL) {

  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_tile(color = NA) +
    {if (!is.null(fill_scale)) fill_scale
     else scale_fill_simmat()} +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = x_lab, y = y_lab) +
    theme_plate(grid_x = FALSE, grid_y = FALSE) +
    theme(axis.text = element_text(size = 6),
          panel.background = element_rect(fill = "white")) +
    coord_fixed()
  p
}

# =============================================================================
# 8. HELPER: EMBEDDED TABLE (as grob for patchwork)
# =============================================================================

plate_table <- function(df, title = NULL, header_fill = "#F0F0F0",
                         header_color = "#2A2A2A", base_size = 8) {

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    message("[theme] gridExtra not installed — table as text")
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                label = paste(capture.output(print(df)), collapse = "\n"),
                                size = 2.5, family = "mono") +
             theme_void() + labs(title = title))
  }

  tt <- gridExtra::ttheme_minimal(
    base_size = base_size,
    core = list(
      bg_params = list(fill = c("white", "#F8F8F8"), col = "#E0E0E0", lwd = 0.5),
      fg_params = list(col = "#3A3A3A", fontface = "plain")
    ),
    colhead = list(
      bg_params = list(fill = header_fill, col = "#D0D0D0", lwd = 0.5),
      fg_params = list(col = header_color, fontface = "bold")
    )
  )

  tbl_grob <- gridExtra::tableGrob(df, rows = NULL, theme = tt)

  # Wrap in ggplot for patchwork compatibility
  p <- ggplot() +
    annotation_custom(tbl_grob) +
    theme_void() +
    labs(title = title) +
    theme(plot.title = element_text(size = base_size + 2, face = "bold",
                                     color = "#2A2A2A", hjust = 0, margin = margin(b = 4)))
  p
}

# =============================================================================
# 9. HELPER: CALLOUT / SUMMARY BOX
# =============================================================================

plate_callout <- function(text, title = NULL, width = 0.9,
                           bg = "#F8FAFC", border = "#D0D0D0") {

  p <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = bg, color = border, linewidth = 0.5) +
    annotate("text", x = 0.05, y = 0.5, label = text,
             hjust = 0, vjust = 0.5, size = 2.8, color = "#3A3A3A",
             family = "mono", lineheight = 1.3) +
    {if (!is.null(title))
      annotate("text", x = 0.05, y = 0.92, label = title,
               hjust = 0, size = 3.2, fontface = "bold", color = "#1A1A1A")} +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
  p
}

# =============================================================================
# 10. COMPOSITE FIGURE ASSEMBLY (patchwork)
# =============================================================================

plate_composite <- function(plots, layout = NULL, title = NULL,
                              subtitle = NULL, tag_levels = "A") {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("[theme] patchwork not installed — returning first plot")
    return(plots[[1]])
  }
  library(patchwork)

  # Combine
  combined <- Reduce(`+`, plots)
  if (!is.null(layout)) combined <- combined + plot_layout(design = layout)

  # Add labels and shared title
  combined <- combined +
    plot_annotation(
      title = title, subtitle = subtitle,
      tag_levels = tag_levels,
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", color = "#1A1A1A",
                                   margin = margin(b = 4)),
        plot.subtitle = element_text(size = 10, color = "#5A5A5A",
                                      margin = margin(b = 8))
      )
    ) &
    theme(plot.tag = element_text(size = 12, face = "bold", color = "#1A1A1A"))

  combined
}

# =============================================================================
# 11. QUICK REFERENCE
# =============================================================================
#
# Usage:
#   source("theme_systems_plate.R")
#
#   # Change palette (optional)
#   set_palette(pal_slate_teal)
#
#   # Simple track plot
#   p <- plate_track(dt, "pos_mb", "inv_likeness",
#                     title = "LG06 — Inv-likeness",
#                     caption = "Source: precomp v8.5")
#
#   # Composite figure
#   fig <- plate_composite(
#     list(p1, p2, p3, tbl),
#     layout = "AAB\nCCB",
#     title = "C_gar_LG06 — Candidate Region 36-41 Mb"
#   )
#
#   ggsave("figure.pdf", fig, width = 14, height = 10)
#
# Available palettes:
#   pal_editorial_blue    (default)
#   pal_deep_indigo
#   pal_slate_teal
#   pal_steel_lavender
#   pal_warm_technical
#   pal_neutral_scientific
#
# =============================================================================
