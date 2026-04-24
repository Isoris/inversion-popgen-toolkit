#!/usr/bin/env Rscript

# =============================================================================
# STEP35_candidate_harmonized_heatmaps.R
#
# Multi-mode dosage heatmaps with polarity harmonization.
#
# Produces for each candidate:
#   1. Raw genomic-order heatmap (all samples)
#   2. Harmonized genomic-order heatmap (polarity-corrected)
#   3. Harmonized core-only heatmap
#   4. Stripe-3-only heatmap
#   5. Stripe-2-only heatmap
#   6. Stripe-2 + Stripe-3 comparison heatmap
#
# All heatmaps include:
#   - Row annotations: coarse group, coherence class, stripe quality
#   - Column annotations: genomic position, polarity, marker block
#   - Human-readable row/column numbering via lookup tables
#
# Usage:
#   Rscript STEP35_candidate_harmonized_heatmaps.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

has_CH <- suppressWarnings(require(ComplexHeatmap, quietly = TRUE))
has_circlize <- suppressWarnings(require(circlize, quietly = TRUE))

if (!has_CH) {
  message("[WARN] ComplexHeatmap not available — skipping heatmaps")
  quit("no", status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

# ── Heatmap rendering helper ────────────────────────────────────────────────
render_heatmap <- function(X, row_info, col_info, title, filepath,
                            max_markers = 300) {
  n_row <- nrow(X)
  n_col <- ncol(X)
  if (n_row < 3 || n_col < 3) return(invisible())

  # Subsample markers if too many
  if (n_col > max_markers) {
    marker_var <- apply(X, 2, var, na.rm = TRUE)
    top_idx <- order(marker_var, decreasing = TRUE)[seq_len(max_markers)]
    top_idx <- sort(top_idx)
    X <- X[, top_idx, drop = FALSE]
    col_info <- col_info[top_idx]
  }

  # Row annotations
  ha_left <- ComplexHeatmap::rowAnnotation(
    Group = row_info$coarse_group,
    Quality = row_info$stripe_quality,
    col = list(
      Group = GROUP_COLORS[intersect(names(GROUP_COLORS),
                                      unique(row_info$coarse_group))],
      Quality = c("core" = "#1B7837", "peripheral" = "#F4A582",
                   "junk" = "#D73027")[intersect(c("core", "peripheral", "junk"),
                                                  unique(row_info$stripe_quality))]
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 7)
  )

  # Column annotations
  pol_col <- c("FALSE" = "#2166AC", "TRUE" = "#D73027")
  ha_top <- NULL
  if ("polarity_reversed" %in% names(col_info)) {
    ha_top <- ComplexHeatmap::columnAnnotation(
      Polarity = as.character(col_info$polarity_reversed),
      col = list(Polarity = pol_col),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 7)
    )
  }

  # Render
  pdf(sub("\\.png$", ".pdf", filepath),
      width = max(8, ncol(X) * 0.03 + 5),
      height = max(6, nrow(X) * 0.04 + 3))

  ht <- ComplexHeatmap::Heatmap(
    X,
    name = "Dosage",
    col = circlize::colorRamp2(c(0, 1, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    left_annotation = ha_left,
    top_annotation = ha_top,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    use_raster = nrow(X) * ncol(X) > 50000
  )

  ComplexHeatmap::draw(ht)
  dev.off()

  # PNG
  png(filepath,
      width = max(8, ncol(X) * 0.03 + 5),
      height = max(6, nrow(X) * 0.04 + 3),
      units = "in", res = 300)
  ComplexHeatmap::draw(ht)
  dev.off()

  message("  Saved: ", basename(filepath))
}

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  plot_dir <- file.path(PLOTS_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(plot_dir)

  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  coh_file <- file.path(cand_dir, "candidate_sample_coherence.tsv")
  pol_file <- file.path(cand_dir, "candidate_marker_polarity.tsv")

  if (!file.exists(rot_file)) next
  rot <- fread(rot_file)
  coh <- if (file.exists(coh_file)) fread(coh_file) else NULL
  pol <- if (file.exists(pol_file)) fread(pol_file) else NULL

  # Load dosage
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) next

  dos <- fread(dos_file)
  sites <- fread(sites_file)
  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 20) next

  dos_reg <- dos[keep]
  sites_reg <- sites[keep]
  sc <- setdiff(names(dos_reg), "marker")
  if (all(grepl("^Ind", sc)) && length(sc) == nrow(rot)) {
    setnames(dos_reg, old = sc, new = rot$sample)
  }
  common <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))
  if (length(common) < 10) next

  X <- as.matrix(dos_reg[, ..common])
  storage.mode(X) <- "double"

  message("[INFO] Candidate ", cid, ": harmonized heatmaps")

  # ── Build row info ─────────────────────────────────────────────────────
  row_info_dt <- rot[match(common, sample), .(sample, coarse_group = coarse_group_refined, u)]
  if (!is.null(coh)) {
    row_info_dt <- merge(row_info_dt,
                         coh[, .(sample, coherence_class, stripe_quality)],
                         by = "sample", all.x = TRUE)
  } else {
    row_info_dt[, coherence_class := "unknown"]
    row_info_dt[, stripe_quality := "unknown"]
  }

  # Order: by group then by u
  setorder(row_info_dt, coarse_group, u)
  ordered_samples <- row_info_dt$sample

  X_ordered <- t(X)[match(ordered_samples, common), , drop = FALSE]
  rownames(X_ordered) <- ordered_samples

  # ── Build col info ─────────────────────────────────────────────────────
  col_info_dt <- data.table(
    marker_idx = seq_len(nrow(X)),
    pos = sites_reg$pos
  )
  if (!is.null(pol) && nrow(pol) == nrow(X)) {
    col_info_dt[, polarity_reversed := pol$polarity_reversed]
    col_info_dt[, delta_hom := pol$delta_hom]
  }

  # ── 1. Raw genomic-order heatmap ───────────────────────────────────────
  render_heatmap(
    X_ordered, row_info_dt, col_info_dt,
    paste0("Raw dosage heatmap | Candidate ", cid),
    file.path(plot_dir, "heatmap_raw_genomic_all.png")
  )

  # ── 2. Harmonized heatmap ──────────────────────────────────────────────
  if (!is.null(pol) && "polarity_reversed" %in% names(col_info_dt)) {
    X_harm <- X
    rev_idx <- which(col_info_dt$polarity_reversed)
    if (length(rev_idx) > 0) {
      X_harm[rev_idx, ] <- 2 - X_harm[rev_idx, ]
    }
    X_harm_ordered <- t(X_harm)[match(ordered_samples, common), , drop = FALSE]

    render_heatmap(
      X_harm_ordered, row_info_dt, col_info_dt,
      paste0("Harmonized dosage | Candidate ", cid, " | ",
             length(rev_idx), " markers flipped"),
      file.path(plot_dir, "heatmap_harmonized_genomic_all.png")
    )

    # ── 3. Core-only heatmap ─────────────────────────────────────────────
    if (!is.null(coh)) {
      core_samp <- coh[stripe_quality == "core", sample]
      core_samp <- intersect(core_samp, ordered_samples)
      if (length(core_samp) >= 10) {
        core_idx <- match(core_samp, ordered_samples)
        core_row_info <- row_info_dt[sample %in% core_samp]

        render_heatmap(
          X_harm_ordered[core_idx, , drop = FALSE],
          core_row_info, col_info_dt,
          paste0("Core-only harmonized | Candidate ", cid,
                 " | n=", length(core_samp)),
          file.path(plot_dir, "heatmap_harmonized_genomic_core.png")
        )
      }
    }

    # ── 4-6. Stripe-specific heatmaps ────────────────────────────────────
    for (stripe_name in c("HOMO_2", "HET")) {
      stripe_samp <- rot[coarse_group_refined == stripe_name, sample]
      stripe_samp <- intersect(stripe_samp, ordered_samples)
      if (length(stripe_samp) >= 5) {
        s_idx <- match(stripe_samp, ordered_samples)
        s_row <- row_info_dt[sample %in% stripe_samp]
        s_label <- if (stripe_name == "HOMO_2") "stripe3" else "stripe2"

        render_heatmap(
          X_harm_ordered[s_idx, , drop = FALSE],
          s_row, col_info_dt,
          paste0(stripe_name, " only | Candidate ", cid,
                 " | n=", length(stripe_samp)),
          file.path(plot_dir, paste0("heatmap_", s_label, "_only.png"))
        )
      }
    }

    # Stripe 2 + Stripe 3 comparison
    s23_samp <- rot[coarse_group_refined %in% c("HET", "HOMO_2"), sample]
    s23_samp <- intersect(s23_samp, ordered_samples)
    if (length(s23_samp) >= 10) {
      s23_idx <- match(s23_samp, ordered_samples)
      s23_row <- row_info_dt[sample %in% s23_samp]

      render_heatmap(
        X_harm_ordered[s23_idx, , drop = FALSE],
        s23_row, col_info_dt,
        paste0("HET + HOMO_2 comparison | Candidate ", cid),
        file.path(plot_dir, "heatmap_stripe2_stripe3_comparison.png")
      )
    }
  }
}

message("[DONE] STEP35 complete")
