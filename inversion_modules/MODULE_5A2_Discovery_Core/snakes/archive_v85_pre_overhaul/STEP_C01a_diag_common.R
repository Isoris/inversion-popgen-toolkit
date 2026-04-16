#!/usr/bin/env Rscript

# =============================================================================
# C01a_diag_common.R — Shared functions, theme, data loading for diag scripts
#
# Sourced by: C01a_diag_profiles.R, C01a_diag_matrices.R, C01a_diag_mds.R
# NOT run standalone.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# THEME (sourced from theme_systems_plate.R)
# =============================================================================

# Try to source the systems plate theme
_theme_path <- file.path(dirname(SCRIPT_DIR %||% "."), "utils", "theme_systems_plate.R")
if (!exists("_theme_path") || !file.exists(_theme_path)) {
  _theme_path <- Sys.glob(file.path(dirname(dirname(
    sys.frame(1)$ofile %||% ".")), "utils", "theme_systems_plate.R"))[1]
}
if (!is.null(_theme_path) && file.exists(_theme_path)) {
  source(_theme_path)
  THEME_BASE <- theme_plate()
  message("[diag] Using Publication Systems Plate Theme")
} else {
  # Fallback: minimal clean theme
  DPI <- 350
  THEME_BASE <- theme_minimal(base_size = 9) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#FAFAFA", color = NA),
      plot.title = element_text(size = 10, face = "bold", color = "#1A1A1A"),
      plot.subtitle = element_text(size = 8, color = "#5A5A5A"),
      plot.caption = element_text(size = 6, color = "#8A8A8A", hjust = 0),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.3)
    )
  message("[diag] Using fallback theme (theme_systems_plate.R not found)")
}

DPI <- 350

# Color palettes
PAL_INV <- c("PASS" = "green4", "WEAK" = "orange", "FAIL" = "red3")
PAL_Z   <- c("z>=5" = "navy", "z>=4" = "blue4", "z>=3" = "royalblue",
             "z>=2" = "cornflowerblue", "z>=1.5" = "lightblue", "z<1.5" = "gray90")

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

  # Inv-likeness table
  inv_like_file <- file.path(dirname(precomp_dir), "snake_inv_likeness.tsv.gz")
  inv_like_dt <- if (file.exists(inv_like_file)) fread(inv_like_file) else data.table()
  if (nrow(inv_like_dt) > 0) message("[diag] Inv-likeness: ", nrow(inv_like_dt), " windows")

  list(precomp_list = precomp_list, chroms = chroms, inv_like_dt = inv_like_dt)
}

# =============================================================================
# PLOT RENDERING ENGINE (shared by all 3 scripts)
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

    # PDF
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

    # PNG per chromosome
    png_dir <- file.path(outdir, "png", pname)
    dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
    for (chr in names(plots)) {
      f_png <- file.path(png_dir, paste0(chr, ".png"))
      tryCatch(
        ggsave(f_png, plots[[chr]], width = dims["w"], height = dims["h"], dpi = DPI),
        error = function(e) message("  [PNG FAIL] ", chr, ": ", e$message))
    }
    message("  ", pname, ": ", length(plots), " pages -> ", f_pdf)
  }
}

# =============================================================================
# HELPER: rescale to [0,1]
# =============================================================================

rescale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0.5, length(x)))
  (x - r[1]) / diff(r)
}

# =============================================================================
# PARSE COMMON ARGS
# =============================================================================

parse_diag_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("Usage: Rscript <script> <precomp_dir> <outdir> [chrom]")
  list(
    precomp_dir = args[1],
    outdir = args[2],
    chrom_filter = if (length(args) >= 3 && args[3] != "all") args[3] else NULL
  )
}
