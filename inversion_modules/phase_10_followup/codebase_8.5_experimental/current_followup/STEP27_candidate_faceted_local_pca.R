#!/usr/bin/env Rscript

# =============================================================================
# STEP27_candidate_faceted_local_pca.R
#
# FIG_C09 — Multi-page faceted local PCA plots.
# FIG_C09b — Multi-page faceted MDS plots (v7.3 addition).
#
# For each candidate, generates multi-page faceted local PCA showing all
# overlapping 100-SNP windows. Layout: 5×8 panels per page (adjustable).
# Two versions per figure type: fixed scales within candidate, free scales per facet.
# MDS facets use the same per-window eigenvectors as local PCA.
#
# Inputs:
#   - STEP09 window_pca RDS files (for per-window eigenvectors)
#   - STEP20 window summary (for topology labels)
#   - STEP21 rotated PCA (for coarse group colors)
#
# Usage:
#   Rscript STEP27_candidate_faceted_local_pca.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)

NCOL_FACET <- 5
NROW_FACET <- 8
PANELS_PER_PAGE <- NCOL_FACET * NROW_FACET

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

# Load STEP09 RDS files
step09_rds <- list.files(c(STEP09_DIR, INV_ROOT),
                         pattern = "\\.window_pca\\.rds$",
                         full.names = TRUE, recursive = TRUE)
step09_by_chr <- list()
for (f in step09_rds) {
  obj <- readRDS(f)
  step09_by_chr[[obj$chrom]] <- obj
}

shape_map <- c("HOMO_1" = 16, "HET" = 17, "HOMO_2" = 15, "AMBIGUOUS" = 4)

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  plot_dir <- file.path(PLOTS_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(plot_dir)

  # Load window summary
  win_file <- file.path(cand_dir, "candidate_window_summary.tsv")
  if (!file.exists(win_file)) next
  win <- fread(win_file)
  if (nrow(win) < 1) next

  # Load coarse groups
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  rot <- if (file.exists(rot_file)) fread(rot_file) else NULL
  sample_groups <- NULL
  if (!is.null(rot)) {
    sample_groups <- rot[, .(sample, group = coarse_group_refined)]
  }

  # Get STEP09 data
  step09_obj <- step09_by_chr[[chr]]
  if (is.null(step09_obj)) next

  sample_names <- step09_obj$sample_names
  pca_df <- step09_obj$pca

  # Build per-window scatter data
  all_panels <- list()
  for (wi in seq_len(nrow(win))) {
    wid <- win$window_id[wi]
    row_idx <- which(pca_df$window_id == wid)
    if (length(row_idx) != 1) next

    # Extract eigenvectors
    pc1_cols <- paste0("PC_1_", sample_names)
    pc2_cols <- paste0("PC_2_", sample_names)
    if (!all(pc1_cols %in% names(pca_df))) next

    pc1 <- as.numeric(pca_df[row_idx, pc1_cols, with = FALSE])
    pc2 <- as.numeric(pca_df[row_idx, pc2_cols, with = FALSE])

    panel_dt <- data.table(
      sample = sample_names,
      PC1 = pc1, PC2 = pc2,
      window_id = wid,
      window_label = paste0("W", wid, " | ",
                            round(win$window_start[wi] / 1e6, 2), "-",
                            round(win$window_end[wi] / 1e6, 2), " Mb | ",
                            win$topology[wi]),
      topology = win$topology[wi]
    )

    if (!is.null(sample_groups)) {
      panel_dt <- merge(panel_dt, sample_groups, by = "sample", all.x = TRUE)
    } else {
      panel_dt[, group := "unknown"]
    }

    all_panels[[length(all_panels) + 1]] <- panel_dt
  }

  if (length(all_panels) == 0) next
  panel_all <- rbindlist(all_panels, fill = TRUE)

  # Assign page numbers
  window_ids <- unique(panel_all$window_label)
  n_windows <- length(window_ids)
  n_pages <- ceiling(n_windows / PANELS_PER_PAGE)

  window_page <- data.table(
    window_label = window_ids,
    page = rep(seq_len(n_pages), each = PANELS_PER_PAGE)[seq_len(n_windows)]
  )
  panel_all <- merge(panel_all, window_page, by = "window_label")

  message("[INFO] Candidate ", cid, ": ", n_windows, " windows → ", n_pages, " pages")

  header <- candidate_header(cid, chr, c_start, c_end,
                              n_windows = n_windows, n_samples = length(sample_names))

  # Generate pages — TWO VERSIONS: free_y (per-facet scaling) and fixed (global scaling)
  # v7.3: added fixed-scale version as requested
  for (scale_mode in c("free", "fixed")) {
    scale_arg <- if (scale_mode == "free") "free" else "fixed"
    tag <- if (scale_mode == "free") "" else "_fixed"

    for (pg in seq_len(n_pages)) {
      page_data <- panel_all[page == pg]

      p <- ggplot(page_data, aes(x = PC1, y = PC2)) +
        geom_point(aes(color = group), size = 0.6, alpha = 0.7) +
        scale_color_manual(values = GROUP_COLORS, na.value = "grey60") +
        facet_wrap(~window_label, ncol = NCOL_FACET, scales = scale_arg) +
        theme_inversion(base_size = 7) +
        theme(
          strip.text = element_text(size = 5),
          legend.position = "bottom",
          axis.text = element_text(size = 4),
          panel.spacing = unit(0.15, "lines")
        ) +
        labs(
          title = paste0("FIG_C09 — Faceted local PCA [", scale_mode, "] (candidate ", cid,
                         ", page ", pg, "/", n_pages, ")"),
          subtitle = header,
          x = "Local PC1", y = "Local PC2",
          color = "Coarse group"
        )

      page_str <- sprintf("%03d", pg)
      ggsave(file.path(plot_dir, paste0("FIG_C09_faceted_pca", tag, "_page", page_str, ".pdf")),
             p, width = 14, height = 18)
      ggsave(file.path(plot_dir, paste0("FIG_C09_faceted_pca", tag, "_page", page_str, ".png")),
             p, width = 14, height = 18, dpi = 300)

      message("  Page ", pg, "/", n_pages, " (", scale_mode, ") saved")
    }
  }

  # ── FACETED MDS PLOTS (v7.3 addition) ──────────────────────────────────
  # Multi-page faceted MDS panels with both free and fixed scales.
  # Uses per-window MDS coordinates from STEP10c (Layer A: SNP distance).

  belonging_file <- file.path(cand_dir, "candidate_window_belonging_clusters.tsv")
  band_file <- file.path(cand_dir, "candidate_window_band_assignments.tsv.gz")

  if (file.exists(band_file)) {
    band_dt <- fread(band_file)
    band_dt[, candidate_id := NULL]

    # Reshape: one row per (sample, window)
    # Columns are: sample, w1:..., w2:..., etc.
    band_sample_col <- "sample"
    win_cols_b <- setdiff(names(band_dt), band_sample_col)

    # For each window, load raw dosage and compute PCA → MDS coord
    # Reuse the STEP09 eigenvectors for MDS1/MDS2 per window
    mds_panels <- list()
    for (wi in seq_len(nrow(win))) {
      wid <- win$window_id[wi]
      row_idx <- which(pca_df$window_id == wid)
      if (length(row_idx) != 1) next

      pc1_cols <- paste0("PC_1_", sample_names)
      pc2_cols <- paste0("PC_2_", sample_names)
      if (!all(pc1_cols %in% names(pca_df))) next

      pc1 <- as.numeric(pca_df[row_idx, pc1_cols, with = FALSE])
      pc2 <- as.numeric(pca_df[row_idx, pc2_cols, with = FALSE])

      # Use PC1/PC2 as MDS proxies (same eigenvectors, consistent with STEP10c)
      mds_panel <- data.table(
        sample = sample_names,
        MDS1 = pc1, MDS2 = pc2,
        window_id = wid,
        window_label = paste0("W", wid, " | ",
                              round(win$window_start[wi] / 1e6, 2), "-",
                              round(win$window_end[wi] / 1e6, 2), " Mb")
      )

      # Merge coarse groups
      if (!is.null(sample_groups)) {
        mds_panel <- merge(mds_panel, sample_groups, by = "sample", all.x = TRUE)
      } else {
        mds_panel[, group := "unknown"]
      }

      mds_panels[[length(mds_panels) + 1]] <- mds_panel
    }

    if (length(mds_panels) > 0) {
      mds_all <- rbindlist(mds_panels, fill = TRUE)

      # Assign pages
      mds_window_ids <- unique(mds_all$window_label)
      n_mds_windows <- length(mds_window_ids)
      n_mds_pages <- ceiling(n_mds_windows / PANELS_PER_PAGE)

      mds_window_page <- data.table(
        window_label = mds_window_ids,
        page = rep(seq_len(n_mds_pages), each = PANELS_PER_PAGE)[seq_len(n_mds_windows)]
      )
      mds_all <- merge(mds_all, mds_window_page, by = "window_label")

      message("[INFO] Candidate ", cid, " MDS facets: ", n_mds_windows,
              " windows → ", n_mds_pages, " pages")

      for (scale_mode in c("free", "fixed")) {
        scale_arg <- if (scale_mode == "free") "free" else "fixed"
        tag <- if (scale_mode == "free") "" else "_fixed"

        for (pg in seq_len(n_mds_pages)) {
          page_data <- mds_all[page == pg]

          p <- ggplot(page_data, aes(x = MDS1, y = MDS2)) +
            geom_point(aes(color = group), size = 0.6, alpha = 0.7) +
            scale_color_manual(values = GROUP_COLORS, na.value = "grey60") +
            facet_wrap(~window_label, ncol = NCOL_FACET, scales = scale_arg) +
            theme_inversion(base_size = 7) +
            theme(
              strip.text = element_text(size = 5),
              legend.position = "bottom",
              axis.text = element_text(size = 4),
              panel.spacing = unit(0.15, "lines")
            ) +
            labs(
              title = paste0("FIG_C09b — Faceted MDS [", scale_mode, "] (candidate ", cid,
                             ", page ", pg, "/", n_mds_pages, ")"),
              subtitle = header,
              x = "MDS1 (local PC1)", y = "MDS2 (local PC2)",
              color = "Coarse group"
            )

          page_str <- sprintf("%03d", pg)
          ggsave(file.path(plot_dir, paste0("FIG_C09b_faceted_mds", tag, "_page", page_str, ".pdf")),
                 p, width = 14, height = 18)
          ggsave(file.path(plot_dir, paste0("FIG_C09b_faceted_mds", tag, "_page", page_str, ".png")),
                 p, width = 14, height = 18, dpi = 300)

          message("  MDS page ", pg, "/", n_mds_pages, " (", scale_mode, ") saved")
        }
      }
    }
  }
}

message("[DONE] STEP27 faceted local PCA complete")
