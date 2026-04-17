#!/usr/bin/env Rscript

# =============================================================================
# debug_core_pcas.R — Faceted PCA per core/merge region
#
# For each detected core (or merge region), extracts the PC1×PC2 scatter
# of 226 samples across that core's windows, producing a multi-page faceted
# PDF where each facet = one core. Labels show core ID, position, span,
# n_windows, family, z-score.
#
# Usage:
#   Rscript debug_core_pcas.R <precomp_dir> <core_windows_tsv> <outdir> \
#     [--chrom <chr>] [--per_page 50] [--mode cores|merge]
#
# Inputs:
#   precomp_dir:       directory with <chr>.precomp.rds files
#   core_windows_tsv:  snake1_core_windows_<chr>.tsv.gz (or merge equivalent)
#   outdir:            where to write PDF/PNG
#
# Outputs:
#   debug_pcas_<chr>.pdf — multi-page faceted PCA (50 per page default)
#   debug_pcas_<chr>.png — first page as PNG thumbnail
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript debug_core_pcas.R <precomp_dir> <core_windows_tsv> <outdir> [--chrom chr] [--per_page 50]")

precomp_dir  <- args[1]
core_file    <- args[2]
outdir       <- args[3]
chrom_filter <- NULL
PER_PAGE     <- 50L
mode         <- "cores"

i <- 4L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args))    { chrom_filter <- args[i+1]; i <- i+2 }
  else if (a == "--per_page" && i < length(args)) { PER_PAGE <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--mode" && i < length(args)) { mode <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load cores/merge windows ──
core_dt <- fread(core_file)
if (!is.null(chrom_filter)) core_dt <- core_dt[chrom == chrom_filter]
message("[DEBUG] Loaded ", nrow(core_dt), " window records from ", basename(core_file))

# Identify core/region ID column
id_col <- intersect(c("core_id", "region_id", "merge_id", "region_id"), names(core_dt))
if (length(id_col) == 0) {
  # Try to construct from available columns
  if ("scale_tier" %in% names(core_dt) && "core_idx" %in% names(core_dt)) {
    core_dt[, core_id := paste0(scale_tier, "_", core_idx)]
    id_col <- "core_id"
  } else {
    stop("Cannot find core/region ID column")
  }
} else {
  id_col <- id_col[1]
}

# Get unique cores/regions
cores <- unique(core_dt[[id_col]])
message("[DEBUG] ", length(cores), " unique ", mode, " to plot")

# ── Process per chromosome ──
chroms <- if (!is.null(chrom_filter)) chrom_filter else unique(core_dt$chrom)

for (chr in chroms) {
  message("\n[DEBUG] === ", chr, " ===")

  # Load precomp
  rds_file <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
  if (!file.exists(rds_file)) { message("  [SKIP] no precomp"); next }
  pc <- readRDS(rds_file)
  dt <- pc$dt

  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  pc2_cols <- grep("^PC_2_", names(dt), value = TRUE)
  if (length(pc1_cols) < 20 || length(pc2_cols) < 20) {
    message("  [SKIP] too few PC columns"); next
  }
  sample_names <- sub("^PC_1_", "", pc1_cols)

  chr_cores <- core_dt[chrom == chr]
  chr_core_ids <- unique(chr_cores[[id_col]])
  n_cores <- length(chr_core_ids)
  if (n_cores == 0) next
  message("  ", n_cores, " ", mode, " on ", chr)

  # ── Build facet data ──
  all_facet_data <- list()

  for (ci in seq_along(chr_core_ids)) {
    cid <- chr_core_ids[ci]
    c_wins <- chr_cores[get(id_col) == cid]
    win_ids <- c_wins$global_window_id
    win_idx <- match(win_ids, dt$global_window_id)
    win_idx <- win_idx[!is.na(win_idx)]
    if (length(win_idx) < 1) next

    # Average PC1 and PC2 across core's windows
    pc1_mat <- as.matrix(dt[win_idx, ..pc1_cols])
    pc2_mat <- as.matrix(dt[win_idx, ..pc2_cols])
    avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
    avg_pc2 <- colMeans(pc2_mat, na.rm = TRUE)

    # Core metadata for label
    start_mb <- round(min(c_wins$start_bp, na.rm = TRUE) / 1e6, 2)
    end_mb <- round(max(c_wins$end_bp, na.rm = TRUE) / 1e6, 2)
    n_win <- length(win_idx)
    fam <- if ("scale_tier" %in% names(c_wins)) c_wins$scale_tier[1] else "?"

    # k=3 for coloring
    valid <- is.finite(avg_pc1) & is.finite(avg_pc2)
    if (sum(valid) < 15) next

    km <- tryCatch(kmeans(avg_pc1[valid], centers = 3, nstart = 5),
                    error = function(e) NULL)
    band <- rep("?", length(avg_pc1))
    if (!is.null(km)) {
      co <- order(km$centers[, 1])
      band[valid][km$cluster == co[1]] <- "Band1"
      band[valid][km$cluster == co[2]] <- "Band2"
      band[valid][km$cluster == co[3]] <- "Band3"
    }

    label <- paste0(cid, " | ", start_mb, "-", end_mb, "Mb | ",
                     n_win, "w | ", fam)

    all_facet_data[[length(all_facet_data) + 1]] <- data.table(
      sample = sample_names, pc1 = avg_pc1, pc2 = avg_pc2,
      band = band, facet_label = label, core_order = ci
    )
  }

  if (length(all_facet_data) == 0) { message("  No valid cores"); next }
  pdt <- rbindlist(all_facet_data)

  # ── Paginated PDF ──
  n_pages <- ceiling(n_cores / PER_PAGE)
  pdf_file <- file.path(outdir, paste0("debug_pcas_", chr, ".pdf"))

  pdf(pdf_file, width = 24, height = 20)

  for (page in seq_len(n_pages)) {
    start_idx <- (page - 1) * PER_PAGE + 1
    end_idx <- min(page * PER_PAGE, n_cores)
    page_labels <- unique(pdt$facet_label)[start_idx:end_idx]
    page_data <- pdt[facet_label %in% page_labels]

    # Determine grid dimensions
    n_facets <- length(page_labels)
    n_cols <- ceiling(sqrt(n_facets * 1.5))
    n_rows <- ceiling(n_facets / n_cols)

    p <- ggplot(page_data, aes(x = pc1, y = pc2, color = band)) +
      geom_point(size = 0.5, alpha = 0.6) +
      scale_color_manual(values = c("Band1" = "#3b82f6", "Band2" = "#22c55e",
                                     "Band3" = "#ef4444", "?" = "grey50")) +
      facet_wrap(~ facet_label, ncol = n_cols, scales = "free") +
      labs(title = paste0(chr, " — ", mode, " debug PCA (page ", page, "/", n_pages, ")"),
           subtitle = paste0("Each facet = one ", mode, " | Color = k=3 on avg PC1 | ",
                            n_facets, " ", mode, " on this page"),
           x = "Mean PC1 across core windows",
           y = "Mean PC2 across core windows",
           caption = paste0("Source: precomp PC1/PC2 loadings averaged across each core's windows\n",
                           "3-cluster pattern = inversion (blue=hom_ref, green=het, red=hom_inv)\n",
                           "Continuous blob = family structure or noise")) +
      theme_minimal(base_size = 7) +
      theme(strip.text = element_text(size = 5, face = "bold"),
            legend.position = "bottom",
            plot.title = element_text(size = 10, face = "bold"),
            plot.caption = element_text(size = 5, color = "grey50", hjust = 0),
            panel.spacing = unit(0.3, "lines"))

    print(p)
    message("  Page ", page, ": ", n_facets, " facets")
  }

  dev.off()
  message("  → ", pdf_file)

  # First page as PNG
  png_file <- file.path(outdir, paste0("debug_pcas_", chr, "_p1.png"))
  page1_labels <- unique(pdt$facet_label)[1:min(PER_PAGE, n_cores)]
  page1_data <- pdt[facet_label %in% page1_labels]
  n_facets <- length(page1_labels)
  n_cols <- ceiling(sqrt(n_facets * 1.5))

  p1 <- ggplot(page1_data, aes(x = pc1, y = pc2, color = band)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("Band1" = "#3b82f6", "Band2" = "#22c55e",
                                   "Band3" = "#ef4444", "?" = "grey50")) +
    facet_wrap(~ facet_label, ncol = n_cols, scales = "free") +
    labs(title = paste0(chr, " — ", mode, " debug PCA (page 1)"),
         x = "Mean PC1", y = "Mean PC2") +
    theme_minimal(base_size = 7) +
    theme(strip.text = element_text(size = 5, face = "bold"),
          legend.position = "bottom",
          panel.spacing = unit(0.3, "lines"))

  ggsave(png_file, p1, width = 24, height = 20, dpi = 200)
  message("  → ", png_file)
}

message("\n[DONE]")
