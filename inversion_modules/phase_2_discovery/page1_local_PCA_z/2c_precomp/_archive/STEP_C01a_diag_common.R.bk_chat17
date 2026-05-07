#!/usr/bin/env Rscript

# =============================================================================
# C01a_diag_common.R — Shared functions, theme, data loading for diag scripts
#
# Sourced by: C01a_diag_profiles.R, C01a_diag_matrices.R, C01a_diag_mds.R,
#             C01a_diag_summaries.R
# NOT run standalone.
#
# v8.5.2 — Plot overhaul:
#   G1: sys.frame bug → env vars from 00_inversion_config.sh
#   G2: theme_systems_plate.R sourced reliably
#   G3: Consistent themed palettes (no default ggplot2)
#   G4: Config integration via BASE env var
#   PC4: File naming: {chr}_{plot_type}.png
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# THEME — G1/G4: resolved from BASE env var, never sys.frame()
# =============================================================================

theme_file <- file.path(Sys.getenv("BASE", "."),
                        "inversion_codebase_v8.5/utils/theme_systems_plate.R")

if (file.exists(theme_file)) {
  source(theme_file)
  THEME_BASE <- theme_plate()
  message("[diag] Theme loaded: ", theme_file)
} else {
  # Fallback: publication-quality clean theme
  THEME_BASE <- theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "#FAFAFA", color = NA),
    panel.grid.major  = element_line(color = "#E8E8E8", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "#4A4A4A", linewidth = 0.35),
    axis.ticks        = element_line(color = "#4A4A4A", linewidth = 0.25),
    plot.title        = element_text(size = 12, face = "bold", color = "#1A1A1A",
                                     margin = margin(b = 4), hjust = 0),
    plot.subtitle     = element_text(size = 9, color = "#5A5A5A",
                                     margin = margin(b = 8), hjust = 0),
    plot.caption      = element_text(size = 7, color = "#8A8A8A", hjust = 0,
                                     margin = margin(t = 6), lineheight = 1.2),
    axis.title        = element_text(size = 9, color = "#3A3A3A"),
    axis.title.y      = element_text(margin = margin(r = 6)),
    axis.title.x      = element_text(margin = margin(t = 6)),
    axis.text         = element_text(size = 8, color = "#4A4A4A"),
    strip.background  = element_rect(fill = "#F0F0F0", color = "#D0D0D0",
                                     linewidth = 0.3),
    strip.text        = element_text(size = 9, face = "bold", color = "#2A2A2A",
                                     margin = margin(t = 3, b = 3)),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA),
    legend.key.size   = unit(0.35, "cm"),
    legend.text       = element_text(size = 8, color = "#4A4A4A"),
    legend.title      = element_text(size = 9, face = "bold", color = "#3A3A3A"),
    legend.position   = "bottom",
    legend.margin     = margin(t = 2),
    plot.margin       = margin(t = 8, r = 8, b = 6, l = 6),
    panel.spacing     = unit(0.6, "lines")
  )
  message("[diag] Fallback theme (theme_systems_plate.R not found at: ", theme_file, ")")
}

DPI <- 350

# =============================================================================
# G3: COLOR PALETTES — themed, publication-quality, consistent everywhere
# =============================================================================

# Inv-likeness status (PC14: SEED semantics, not PASS/WEAK/FAIL)
PAL_INV <- c("SEED_ELIGIBLE"   = "#2563EB",
             "SEED_MARGINAL"   = "#9CA3AF",
             "SEED_INELIGIBLE" = "#E8913A")

# Z-score classes (navy gradient)
PAL_Z <- c("z>=5"   = "#0C1E3C",
           "z>=4"   = "#1E3A5F",
           "z>=3"   = "#2563EB",
           "z>=2"   = "#4A90D9",
           "z>=1.5" = "#93C5FD",
           "z<1.5"  = "#E5E7EB")

# Adaptive NN classes
PAL_ADAPT <- c("S1S (q25)"  = "#2563EB",
               "S1M (q50)"  = "#059669",
               "S1L (q75)"  = "#D97706",
               "above q75"  = "#D1D5DB")

# k-NN block signal
PAL_BLOCK <- c("strong_block"   = "#DC2626",
               "moderate_block"  = "#D97706",
               "background"      = "#D1D5DB")

# Band colors (k=3)
PAL_BAND <- c("Band 1 (low PC1)"  = "#3B82F6",
              "Band 2 (mid/het)"   = "#22C55E",
              "Band 3 (high PC1)"  = "#EF4444")

# Family tiers
PAL_FAMILY <- c("S1S" = "#2563EB", "S1M" = "#059669", "S1L" = "#D97706")

# Fixed vs adaptive comparison
PAL_METHOD <- c("Fixed (0.70/0.80/0.90)" = "#D1D5DB",
                "Adaptive (q25/q50/q75)"  = "#3B82F6")

# Similarity matrix color ramp (navy → gold → red)
SIM_COLS <- c("#0C1E3C", "#1E3A5F", "#4A7FB5", "#F0C75E", "#E8913A", "#C0392B")

# Diverging (for contrast matrices)
DIV_LOW  <- "#2563EB"
DIV_MID  <- "#F8F9FA"
DIV_HIGH <- "#DC2626"

# Track line colors
COL_LAMBDA1 <- "#3B82F6"   # steel blue
COL_LAMBDA2 <- "#E8913A"   # warm coral
COL_PVE1    <- "#059669"   # teal green

# Regime bands
PAL_REGIME <- c("core_low" = "#2563EB", "mid" = "#86EFAC", "background" = "#E5E7EB")

# Background quantile annotation colors
COL_Q85 <- "#D97706"
COL_Q90 <- "#059669"
COL_Q95 <- "#DC2626"

# =============================================================================
# DATA LOADING
# =============================================================================

load_diag_data <- function(precomp_dir, chrom_filter = NULL) {
  message("[diag] Loading precomputed data from ", precomp_dir, " ...")
  rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
  if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

  precomp_list <- list()
  for (f in rds_files) {
    obj <- readRDS(f)
    precomp_list[[obj$chrom]] <- obj
  }

  chroms <- names(precomp_list)
  if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)
  message("[diag] ", length(chroms), " chromosomes")

  inv_like_file <- file.path(dirname(precomp_dir), "snake_inv_likeness.tsv.gz")
  inv_like_dt <- if (file.exists(inv_like_file)) fread(inv_like_file) else data.table()
  if (nrow(inv_like_dt) > 0) message("[diag] Inv-likeness: ", nrow(inv_like_dt), " windows")

  list(precomp_list = precomp_list, chroms = chroms, inv_like_dt = inv_like_dt)
}

# =============================================================================
# PLOT RENDERING ENGINE — PC4: file naming {chr}_{plot_type}.png
# =============================================================================

render_plots <- function(plot_builders, plot_dims, chroms, precomp_list, outdir) {
  use_cairo <- capabilities("cairo")

  for (pname in names(plot_builders)) {
    message("[diag] Building: ", pname, " ...")
    builder <- plot_builders[[pname]]
    dims <- plot_dims[[pname]] %||% c(w = 10, h = 7)

    plots <- list()
    for (chr in chroms) {
      p <- tryCatch(builder(chr), error = function(e) {
        message("  [WARN] ", pname, " / ", chr, ": ", e$message); NULL
      })
      if (!is.null(p)) plots[[chr]] <- p
    }

    if (length(plots) == 0) {
      message("  ", pname, ": no plots generated")
      next
    }

    # PDF multi-page
    f_pdf <- file.path(outdir, paste0(pname, ".pdf"))
    if (use_cairo) {
      cairo_pdf(f_pdf, width = dims["w"], height = dims["h"], onefile = TRUE)
    } else {
      pdf(f_pdf, width = dims["w"], height = dims["h"])
    }
    for (chr in names(plots)) {
      tryCatch(print(plots[[chr]]), error = function(e) {
        plot.new(); text(0.5, 0.5, paste0("ERROR: ", chr, "\n", e$message), cex = 0.8)
      })
    }
    dev.off()

    # PNG per chromosome — {chr}_{plot_type}.png
    png_dir <- file.path(outdir, "png", pname)
    dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
    for (chr in names(plots)) {
      f_png <- file.path(png_dir, paste0(chr, "_", pname, ".png"))
      tryCatch(
        ggsave(f_png, plots[[chr]], width = dims["w"], height = dims["h"], dpi = DPI),
        error = function(e) message("  [PNG FAIL] ", chr, ": ", e$message))
    }
    message("  ", pname, ": ", length(plots), " pages -> ", f_pdf)
  }
}

# =============================================================================
# HELPERS
# =============================================================================

rescale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0.5, length(x)))
  (x - r[1]) / diff(r)
}

parse_diag_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("Usage: Rscript <script> <precomp_dir> <outdir> [chrom]")
  list(
    precomp_dir = args[1],
    outdir = args[2],
    chrom_filter = if (length(args) >= 3 && args[3] != "all") args[3] else NULL
  )
}
