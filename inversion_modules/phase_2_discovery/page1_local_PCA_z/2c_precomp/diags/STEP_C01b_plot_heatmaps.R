#!/usr/bin/env Rscript

# =============================================================================
# C01b_plot_heatmaps.R  (v8.5.2)
#
# STANDALONE heatmap plotting for Snake 1 cores.
# Reads saved cores RDS + precomp RDS. No computation.
#
# v8.5.2 — Plot overhaul:
#   G1: theme loaded from env var (no sys.frame)
#   G3: consistent themed palette from theme_systems_plate.R
#   PC4: file naming {chr}_{plot_type}.{ext}
#   DPI: 350 for publication
#
# Produces per chromosome:
#   {chr}_simmat_full.{png,pdf}          — full sim_mat, no PA overlay
#   {chr}_simmat_window_idx.{png,pdf}    — sim_mat by window indices (no gaps)
#   {chr}_simmat_vs_cores.{png,pdf}      — split: upper sim_mat, lower cores PA
#
# Usage:
#   Rscript C01b_plot_heatmaps.R <precomp_dir> <cores_dir> <outdir> [chrom]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# G1: theme from env var
theme_file <- file.path(Sys.getenv("BASE", "."),
                        "inversion_codebase_v8.5/utils/theme_systems_plate.R")
if (file.exists(theme_file)) {
  source(theme_file)
  message("[C01b] Theme loaded: ", theme_file)
} else {
  theme_plate <- function(...) theme_minimal(base_size = 9) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "#FAFAFA", color = NA),
          plot.title = element_text(size = 12, face = "bold", color = "#1A1A1A"),
          plot.subtitle = element_text(size = 9, color = "#5A5A5A"),
          plot.caption = element_text(size = 7, color = "#8A8A8A", hjust = 0))
  message("[C01b] Fallback theme")
}

DPI <- 350

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript C01b_plot_heatmaps.R <precomp_dir> <cores_dir> <outdir> [chrom]")

precomp_dir <- args[1]
cores_dir   <- args[2]
outdir      <- args[3]
chrom_filter <- if (length(args) >= 4 && args[4] != "all") args[4] else NULL

dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

has_newscale <- requireNamespace("ggnewscale", quietly = TRUE)
if (has_newscale) library(ggnewscale)

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (!is.null(chrom_filter)) rds_files <- rds_files[grepl(chrom_filter, rds_files)]

# Themed PA colors
PA_COLORS <- c(
  "Background"       = "#F3F4F6",
  "Seed uncollected" = "#F59E0B",
  "1 family"         = "#3B82F6",
  "2 families"       = "#D97706",
  "All 3 families"   = "#DC2626"
)

# Sim_mat color ramp (consistent with diag_common SIM_COLS)
SIM_COLS <- c("#0C1E3C", "#1E3A5F", "#4A7FB5", "#F0C75E", "#E8913A", "#C0392B")

for (rds_file in rds_files) {
  pc <- readRDS(rds_file)
  chr <- pc$chrom
  dt <- pc$dt
  sim_mat <- pc$sim_mat
  n <- nrow(sim_mat)
  if (n < 50) next

  message("[PLOT] ", chr, " (", n, " windows)")

  state_file <- file.path(cores_dir, paste0("snake1_core_windows_", chr, ".tsv.gz"))
  if (!file.exists(state_file)) {
    message("  [SKIP] no core windows file")
    next
  }
  state_dt <- fread(state_file)

  # Build plot data (subsample max 600x600)
  step_p <- max(1L, n %/% 600)
  idx_p <- seq(1, n, by = step_p)
  ns <- length(idx_p)
  pos_mb <- (dt$start_bp[idx_p] + dt$end_bp[idx_p]) / 2e6

  # Upper triangle: sim_mat
  upper_rows <- vector("list", ns * (ns + 1) / 2)
  k <- 0L
  for (ii in seq_len(ns)) {
    for (jj in ii:ns) {
      k <- k + 1L
      v <- sim_mat[idx_p[ii], idx_p[jj]]
      upper_rows[[k]] <- list(x = pos_mb[ii], y = pos_mb[jj],
                               xi = ii, yi = jj,
                               value = if (is.finite(v)) v else NA_real_)
    }
  }
  upper_dt <- rbindlist(upper_rows)

  # Lower triangle: core PA
  lower_rows <- vector("list", ns * (ns + 1) / 2)
  k <- 0L
  for (ii in seq_len(ns)) {
    for (jj in 1:ii) {
      k <- k + 1L
      wi_i <- idx_p[ii]; wi_j <- idx_p[jj]
      st_i <- state_dt[global_window_id == dt$global_window_id[wi_i]]
      st_j <- state_dt[global_window_id == dt$global_window_id[wi_j]]
      if (nrow(st_i) > 0 && nrow(st_j) > 0) {
        nf <- min(st_i$n_families_claimed[1], st_j$n_families_claimed[1])
        pa_cat <- if (nf >= 3) "All 3 families"
                  else if (nf == 2) "2 families"
                  else if (nf == 1) "1 family"
                  else "Background"
      } else pa_cat <- "Background"
      lower_rows[[k]] <- list(x = pos_mb[ii], y = pos_mb[jj],
                               xi = ii, yi = jj, category = pa_cat)
    }
  }
  lower_dt <- rbindlist(lower_rows)
  lower_dt[, category := factor(category, levels = names(PA_COLORS))]

  sim_med <- median(upper_dt$value, na.rm = TRUE)

  sim_scale <- scale_fill_gradientn(
    colours = SIM_COLS,
    values = scales::rescale(c(0, sim_med * 0.5, sim_med * 0.85,
                                sim_med * 1.1, sim_med * 1.3, 1)),
    name = "Similarity", limits = c(0, 1))

  n_cores <- length(unique(state_dt[n_families_claimed >= 1]$core_id %||% integer(0)))

  # ── Plot A: Full sim_mat (Mb axes) ──
  mirror_dt <- copy(upper_dt[is.finite(value)])
  mirror_dt2 <- copy(mirror_dt)
  mirror_dt2[, c("x", "y") := .(y, x)]
  full_dt <- rbind(mirror_dt, mirror_dt2[x != y])

  p_sim_full <- ggplot(full_dt, aes(x = x, y = y, fill = value)) +
    geom_raster(interpolate = TRUE) +
    sim_scale + coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " \u2014 Similarity Matrix (full)"),
         subtitle = paste0(n, " windows | median sim = ", round(sim_med, 3)),
         caption = "Red blocks = shared PCA structure | Blue = dissimilar") +
    theme_plate() +
    theme(plot.caption = element_text(size = 6, hjust = 0, color = "#8A8A8A"))

  # ── Plot B: Sim_mat by window index ──
  p_idx <- ggplot(upper_dt[is.finite(value)], aes(x = xi, y = yi, fill = value)) +
    geom_raster(interpolate = TRUE) +
    sim_scale + coord_fixed() +
    labs(x = "Window index", y = "Window index",
         title = paste0(chr, " \u2014 Similarity Matrix (window index, no gaps)"),
         subtitle = paste0(ns, "\u00d7", ns, " subsampled from ", n)) +
    theme_plate()

  # ── Plot C: Split heatmap (upper sim_mat, lower cores PA) ──
  if (has_newscale) {
    p_split <- ggplot() +
      geom_raster(data = upper_dt[is.finite(value)],
                   aes(x = x, y = y, fill = value), interpolate = TRUE) +
      sim_scale +
      ggnewscale::new_scale_fill() +
      geom_tile(data = lower_dt, aes(x = x, y = y, fill = category)) +
      scale_fill_manual(values = PA_COLORS, name = "Core detection\n(lower)",
                         drop = FALSE) +
      coord_fixed()
  } else {
    p_split <- ggplot() +
      geom_raster(data = upper_dt[is.finite(value)],
                   aes(x = x, y = y, fill = value), interpolate = TRUE) +
      sim_scale +
      geom_tile(data = lower_dt[category != "Background"],
                 aes(x = x, y = y), fill = "#DC2626", alpha = 0.5) +
      coord_fixed()
  }

  p_split <- p_split +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " \u2014 Split Heatmap: Similarity (upper) vs Cores (lower)"),
         subtitle = paste0(n, " windows | ", n_cores, " cores"),
         caption = paste0("Upper: MDS similarity | Lower: Snake 1 core PA\n",
                         "Red = all 3 families | Amber = 2 | Blue = 1 | Gold = seed uncollected")) +
    theme_plate() +
    theme(legend.position = "right",
          legend.key.size = unit(0.4, "cm"),
          plot.caption = element_text(size = 6, hjust = 0, color = "#8A8A8A"))

  # ── Save (PC4: {chr}_{plot_type}.{ext}) ──
  for (ext in c("png", "pdf")) {
    tryCatch({
      ggsave(file.path(outdir, "plots", paste0(chr, "_simmat_full.", ext)),
             p_sim_full, width = 12, height = 11, dpi = DPI)
      ggsave(file.path(outdir, "plots", paste0(chr, "_simmat_window_idx.", ext)),
             p_idx, width = 12, height = 11, dpi = DPI)
      ggsave(file.path(outdir, "plots", paste0(chr, "_simmat_vs_cores.", ext)),
             p_split, width = 14, height = 12, dpi = DPI)
    }, error = function(e) message("  [", ext, " FAIL] ", e$message))
  }

  message("  [DONE] ", chr, ": 3 plot types \u00d7 2 formats")
}

message("\n[DONE] C01b_plot_heatmaps complete -> ", outdir)
