#!/usr/bin/env Rscript

# =============================================================================
# C01b_plot_heatmaps.R  (v8.5.1)
#
# STANDALONE heatmap plotting for Snake 1 cores.
# Reads saved cores RDS + precomp RDS. No computation.
#
# Produces per chromosome:
#   <chr>_A1_simmat_only.{png,pdf}      — full sim_mat, no PA overlay
#   <chr>_A1_simmat_vs_cores.{png,pdf}  — split: upper sim_mat, lower cores PA
#   <chr>_A1_simmat_window_idx.{png,pdf} — sim_mat with window indices (no gaps)
#
# Usage:
#   Rscript C01b_plot_heatmaps.R <precomp_dir> <cores_dir> <outdir> [chrom]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Try to load theme
theme_file <- file.path(Sys.getenv("BASE", "."), "inversion_codebase_v8.5/utils/theme_systems_plate.R")
if (file.exists(theme_file)) { source(theme_file) } else {
  theme_plate <- function(...) theme_minimal(base_size = 9) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "#FAFAFA", color = NA))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript C01b_plot_heatmaps.R <precomp_dir> <cores_dir> <outdir> [chrom]")

precomp_dir <- args[1]
cores_dir   <- args[2]
outdir      <- args[3]
chrom_filter <- if (length(args) >= 4 && args[4] != "all") args[4] else NULL

dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

has_newscale <- requireNamespace("ggnewscale", quietly = TRUE)
if (has_newscale) library(ggnewscale)

# ── Load precomp ──
rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (!is.null(chrom_filter)) rds_files <- rds_files[grepl(chrom_filter, rds_files)]

PA_COLORS <- c(
  "Background" = "#F3F4F6",
  "Seed uncollected" = "#FBBF24",
  "1 family" = "#3B82F6",
  "2 families" = "#F97316",
  "All 3 families" = "#DC2626"
)

for (rds_file in rds_files) {
  pc <- readRDS(rds_file)
  chr <- pc$chrom
  dt <- pc$dt
  sim_mat <- pc$sim_mat
  n <- nrow(sim_mat)
  if (n < 50) next

  message("[PLOT] ", chr, " (", n, " windows)")

  # Load cores state
  state_file <- file.path(cores_dir, paste0("snake1_core_windows_", chr, ".tsv.gz"))
  if (!file.exists(state_file)) {
    message("  [SKIP] no core windows file")
    next
  }
  state_dt <- fread(state_file)

  # ── Build plot data ──
  # Subsample: max 600×600
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
    colours = c("#0C1E3C", "#1E3A5F", "#4A7FB5", "#E8B84A", "#D4712A", "#B8282E"),
    values = scales::rescale(c(0, sim_med * 0.5, sim_med * 0.85,
                                sim_med * 1.1, sim_med * 1.3, 1)),
    name = "Similarity", limits = c(0, 1))

  n_cores <- length(unique(state_dt[n_families_claimed >= 1]$core_id %||% integer(0)))

  # ── Plot A: Sim_mat only (Mb axes) ──
  p_sim <- ggplot(upper_dt[is.finite(value)], aes(x = x, y = y, fill = value)) +
    geom_raster(interpolate = TRUE) +
    sim_scale + coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " — Similarity Matrix"),
         subtitle = paste0(n, " windows | median sim = ", round(sim_med, 3)),
         caption = "Source: MDS-space correlation of local PCA loading vectors\nRed blocks = shared sample structure | Blue = dissimilar") +
    theme_plate() +
    theme(plot.caption = element_text(size = 6, hjust = 0, color = "#8A8A8A"))

  # Mirror for full matrix display
  mirror_dt <- copy(upper_dt[is.finite(value)])
  mirror_dt2 <- copy(mirror_dt)
  mirror_dt2[, c("x", "y") := .(y, x)]
  full_dt <- rbind(mirror_dt, mirror_dt2[x != y])

  p_sim_full <- ggplot(full_dt, aes(x = x, y = y, fill = value)) +
    geom_raster(interpolate = TRUE) +
    sim_scale + coord_fixed() +
    labs(x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
         title = paste0(chr, " — Similarity Matrix (full, no gaps)"),
         subtitle = paste0(n, " windows | median = ", round(sim_med, 3))) +
    theme_plate()

  # ── Plot B: Sim_mat with window indices (no positional gaps) ──
  p_idx <- ggplot(upper_dt[is.finite(value)], aes(x = xi, y = yi, fill = value)) +
    geom_raster(interpolate = TRUE) +
    sim_scale + coord_fixed() +
    labs(x = "Window index", y = "Window index",
         title = paste0(chr, " — Similarity Matrix (window index, no gaps)"),
         subtitle = paste0(ns, "×", ns, " subsampled from ", n)) +
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
         title = paste0(chr, " — Split Heatmap: Similarity (upper) vs Cores (lower)"),
         subtitle = paste0(n, " windows | ", n_cores, " cores"),
         caption = paste0("Upper: MDS similarity | Lower: Snake 1 core PA\n",
                         "Red = all 3 families | Orange = 2 | Blue = 1 | Yellow = seed uncollected")) +
    theme_plate() +
    theme(legend.position = "right",
          legend.key.size = unit(0.4, "cm"),
          plot.caption = element_text(size = 6, hjust = 0, color = "#8A8A8A"))

  # ── Save ──
  for (ext in c("png", "pdf")) {
    tryCatch({
      ggsave(file.path(outdir, "plots", paste0(chr, "_A1_simmat_only.", ext)),
             p_sim_full, width = 12, height = 11, dpi = 350)
      ggsave(file.path(outdir, "plots", paste0(chr, "_A1_simmat_window_idx.", ext)),
             p_idx, width = 12, height = 11, dpi = 350)
      ggsave(file.path(outdir, "plots", paste0(chr, "_A1_simmat_vs_cores.", ext)),
             p_split, width = 14, height = 12, dpi = 350)
    }, error = function(e) message("  [", ext, " FAIL] ", e$message))
  }

  message("  [DONE] ", chr, ": 3 plot types × 2 formats")
}

message("\n[DONE] C01b_plot_heatmaps complete -> ", outdir)
