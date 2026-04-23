#!/usr/bin/env Rscript
# =============================================================================
# STEP36_ordering_comparison_heatmaps.R   v1
#
# Three-version heatmap renderer for debugging sample ordering on composite
# candidates (motivated by LG28 looking scrambled despite u/v rotation).
#
# Produces, per candidate, THREE versions of the same dosage heatmap differing
# ONLY in how samples are ordered within each karyotype:
#
#   V1  RAW-PC1-NO-ROTATION     [caption: "DEBUG ONLY — not the published method"]
#       Samples ordered by raw PC1 of the shelf region, no u/v transform, no
#       Het-centering. This is the naive approach and is included for
#       comparison so we can see what the u/v rotation actually bought us.
#
#   V2  U-ROTATION-ONLY          [the established STEP21 method, as-is]
#       Samples ordered by the rotated u coordinate from candidate_pca_rotated.tsv.
#       No secondary tie-breaker. This reproduces what STEP28 and STEP35
#       currently produce today.
#
#   V3  U-THEN-HCLUST            [u primary, hclust secondary — experimental]
#       Same as V2 for the main axis, but within each karyotype, ties in u
#       are broken by hierarchical clustering on the group's dosage profile.
#       This is the proposed fix for composite candidates.
#
# Output:
#   ${PLOTS_DIR}/<chr>.candidate_<cid>/
#     heatmap_ordering_comparison.pdf      combined 3-page PDF
#     heatmap_ordering_V1_raw_pc1.pdf      + .png
#     heatmap_ordering_V2_u_only.pdf        + .png
#     heatmap_ordering_V3_u_hclust.pdf      + .png
#     candidate_ordering_audit.tsv          which samples got reordered, and why
#
# Usage:
#   Rscript STEP36_ordering_comparison_heatmaps.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

has_CH <- suppressWarnings(require(ComplexHeatmap, quietly = TRUE))
has_circlize <- suppressWarnings(require(circlize, quietly = TRUE))
if (!has_CH || !has_circlize) {
  message("[WARN] ComplexHeatmap/circlize unavailable — cannot render")
  quit("no", status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_
if (file.exists(config_file)) source(config_file)

if (!exists("FOLLOWUP_DIR"))    FOLLOWUP_DIR    <- "results/followup"
if (!exists("DOSAGE_DIR"))      DOSAGE_DIR      <- "results/dosage"
if (!exists("CANDIDATE_TABLE")) CANDIDATE_TABLE <- "results/candidates.tsv"
if (!exists("PLOTS_DIR"))       PLOTS_DIR       <- "results/plots"
if (!exists("ensure_dir"))      ensure_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); invisible(p) }
if (!exists("GROUP_COLORS"))    GROUP_COLORS <- c(HOMO_1 = "#1f4e79", HET = "#c0504d", HOMO_2 = "#2c7a39")

# Load the shared ordering helper (v2 — u-axis aware)
helper_path <- file.path(dirname(sys.frame(1)$ofile %||% "."), "within_group_ordering.R")
if (!file.exists(helper_path)) {
  # Fallback: same directory as FOLLOWUP_DIR's parent
  helper_path <- "within_group_ordering.R"
}
source(helper_path)

# ── Parameters ─────────────────────────────────────────────────────────────
MAX_MARKERS_DISPLAY <- 300L
TOP_VAR_FOR_HCLUST  <- 500L
MIN_GROUP_SIZE      <- 5L

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# Render helper — one page
# =============================================================================
render_one_page <- function(X_display, row_info, col_positions,
                             title, caption, pdf_path, png_path) {
  if (nrow(X_display) < 3 || ncol(X_display) < 3) return(invisible())

  col_fun <- circlize::colorRamp2(c(0, 1, 2), c("#2166AC", "#F7F7F7", "#B2182B"))

  # Row annotation: karyotype group
  group_v <- row_info$coarse_group
  col_map <- GROUP_COLORS[intersect(names(GROUP_COLORS), unique(group_v))]
  ha_left <- ComplexHeatmap::rowAnnotation(
    Group = group_v,
    col = list(Group = col_map),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 7),
    width = grid::unit(4, "mm")
  )

  # Column annotation: genomic position scatter
  pos_mb <- col_positions / 1e6
  ha_top <- ComplexHeatmap::columnAnnotation(
    pos_mb = ComplexHeatmap::anno_points(
      pos_mb, pch = 16, size = grid::unit(1, "mm"),
      gp = grid::gpar(col = "#555e69"),
      axis_param = list(gp = grid::gpar(fontsize = 6))
    ),
    height = grid::unit(6, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 7)
  )

  n_cell <- nrow(X_display) * ncol(X_display)
  w <- max(7, min(16, 3 + ncol(X_display) * 0.04))
  h <- max(6, min(18, 3 + nrow(X_display) * 0.035))

  ht <- ComplexHeatmap::Heatmap(
    X_display, name = "Dosage", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    left_annotation = ha_left,
    top_annotation = ha_top,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    use_raster = n_cell > 50000,
    raster_quality = 3,
    na_col = "#ffffff"
  )

  draw_with_caption <- function() {
    ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 2, 2, 2), "mm"))
    if (nzchar(caption)) {
      grid::grid.text(caption, x = 0.01, y = 0.01, just = c("left", "bottom"),
                      gp = grid::gpar(fontsize = 8, col = "#8a2b30", fontface = "italic"))
    }
  }

  tryCatch({
    pdf(pdf_path, width = w, height = h); draw_with_caption(); dev.off()
    png(png_path, width = w, height = h, units = "in", res = 150)
    draw_with_caption(); dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("  [WARN] render failed for ", basename(pdf_path), ": ",
            conditionMessage(e))
  })

  list(w = w, h = h, ht = ht, caption = caption)
}

# =============================================================================
# Main loop over candidates
# =============================================================================
cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id; chr <- row$chrom
  c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)
  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  plot_dir <- file.path(PLOTS_DIR,    paste0(chr, ".candidate_", cid))
  ensure_dir(plot_dir)

  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (!file.exists(rot_file)) {
    message("[SKIP] cid=", cid, " — no rotated PCA")
    next
  }
  rot <- fread(rot_file)
  if (nrow(rot) < 10 || !all(c("sample", "u", "coarse_group_refined") %in% names(rot))) {
    message("[SKIP] cid=", cid, " — rotated PCA missing required columns")
    next
  }

  # Load dosage
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[SKIP] cid=", cid, " — dosage files missing"); next
  }
  dos <- fread(dos_file); sites <- fread(sites_file)
  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 30) { message("[SKIP] cid=", cid, " — <30 markers"); next }
  dos_reg <- dos[keep]; sites_reg <- sites[keep]

  # Ind→sample name fixup (as in STEP28/35)
  sc <- setdiff(names(dos_reg), "marker")
  if (all(grepl("^Ind", sc)) && length(sc) == nrow(rot))
    setnames(dos_reg, old = sc, new = rot$sample)
  sample_cols <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))
  if (length(sample_cols) < 10) { message("[SKIP] cid=", cid, " — <10 samples"); next }

  X <- as.matrix(dos_reg[, ..sample_cols]); storage.mode(X) <- "double"
  message("[INFO] cid=", cid, " ", chr, " — 3-version ordering comparison")
  message("  markers=", nrow(X), " samples=", ncol(X))

  # Top-variance markers for display (same as STEP28)
  mv <- apply(X, 1L, var, na.rm = TRUE)
  mv[!is.finite(mv)] <- 0
  n_show <- min(MAX_MARKERS_DISPLAY, nrow(X))
  top_idx <- sort(order(mv, decreasing = TRUE)[seq_len(n_show)])
  X_show <- X[top_idx, , drop = FALSE]
  shown_positions <- sites_reg$pos[top_idx]

  # Row info
  ri_full <- rot[match(sample_cols, sample),
                  .(sample, coarse_group = coarse_group_refined,
                    u = u, v = if ("v" %in% names(rot)) v else NA_real_)]

  # ═══════════════════════════════════════════════════════════════════════
  # VERSION 1 — RAW PC1 (no rotation, no Het-centering)
  #
  # This is the naive baseline. We run prcomp on the full shelf dosage
  # matrix (samples × top-var markers), take PC1, and sort samples by its
  # value within each karyotype group.
  #
  # No u/v transform, no anchor geometry. Just raw PC1.
  #
  # Caption warns this is debug-only and NOT the published method.
  # ═══════════════════════════════════════════════════════════════════════
  message("  V1: raw PC1 (no rotation) ...")
  Xs <- X_show
  for (j in seq_len(ncol(Xs))) {
    bad <- !is.finite(Xs[, j])
    if (any(bad)) Xs[bad, j] <- mean(Xs[, j], na.rm = TRUE)
  }
  Xs[!is.finite(Xs)] <- 0
  raw_pc1 <- tryCatch({
    pc <- prcomp(t(Xs), center = TRUE, scale. = FALSE, rank. = 1)
    as.numeric(pc$x[, 1])
  }, error = function(e) rep(0, ncol(Xs)))
  names(raw_pc1) <- sample_cols

  v1_ord_dt <- copy(ri_full)
  v1_ord_dt[, raw_pc1 := raw_pc1[sample]]
  setorder(v1_ord_dt, coarse_group, raw_pc1)
  v1_samples <- v1_ord_dt$sample

  X_v1 <- t(X_show)[match(v1_samples, sample_cols), , drop = FALSE]
  rownames(X_v1) <- v1_samples
  v1_ri <- v1_ord_dt[match(v1_samples, sample)]

  v1_pdf <- file.path(plot_dir, "heatmap_ordering_V1_raw_pc1.pdf")
  v1_png <- file.path(plot_dir, "heatmap_ordering_V1_raw_pc1.png")
  render_one_page(
    X_v1, v1_ri, shown_positions,
    sprintf("V1 — RAW PC1 (NO ROTATION)  |  cid=%d  %s:%s-%s",
            cid, chr, format(c_start, big.mark=","), format(c_end, big.mark=",")),
    "DEBUG ONLY — naive ordering, not the published method. For visual comparison with V2/V3.",
    v1_pdf, v1_png
  )

  # ═══════════════════════════════════════════════════════════════════════
  # VERSION 2 — U ROTATION ONLY (the established STEP21 method)
  #
  # Order by coarse_group, then by u (from candidate_pca_rotated.tsv).
  # Reproduces exactly what STEP28 / STEP35 / q08 currently output.
  # ═══════════════════════════════════════════════════════════════════════
  message("  V2: u-rotation only ...")
  v2_ord_dt <- copy(ri_full)
  setorder(v2_ord_dt, coarse_group, u)
  v2_samples <- v2_ord_dt$sample

  X_v2 <- t(X_show)[match(v2_samples, sample_cols), , drop = FALSE]
  rownames(X_v2) <- v2_samples
  v2_ri <- v2_ord_dt[match(v2_samples, sample)]

  v2_pdf <- file.path(plot_dir, "heatmap_ordering_V2_u_only.pdf")
  v2_png <- file.path(plot_dir, "heatmap_ordering_V2_u_only.png")
  render_one_page(
    X_v2, v2_ri, shown_positions,
    sprintf("V2 — u ROTATION ONLY (STEP21 method)  |  cid=%d  %s:%s-%s",
            cid, chr, format(c_start, big.mark=","), format(c_end, big.mark=",")),
    "Established method: group primary, u secondary. Reproduces STEP28/STEP35/q08 output.",
    v2_pdf, v2_png
  )

  # ═══════════════════════════════════════════════════════════════════════
  # VERSION 3 — U + HCLUST WITHIN EACH GROUP
  #
  # Primary sort by u (unchanged), secondary tie-breaker by hierarchical
  # clustering within each karyotype. Reveals substructure that u alone
  # cannot capture (e.g., founder sub-haplotypes inside Hom1).
  # ═══════════════════════════════════════════════════════════════════════
  message("  V3: u then hclust within group ...")
  v3_res <- order_samples_with_u_primary(
    samples      = sample_cols,
    groups       = ri_full$coarse_group,
    u            = ri_full$u,
    v            = ri_full$v,
    dosage       = X,   # full dosage matrix, not just the displayed ones
    secondary    = "hclust",
    group_levels = c("HOMO_1", "HET", "HOMO_2"),
    min_group_n  = MIN_GROUP_SIZE,
    top_var_n    = TOP_VAR_FOR_HCLUST
  )
  v3_samples <- v3_res$ordered_samples

  X_v3 <- t(X_show)[match(v3_samples, sample_cols), , drop = FALSE]
  rownames(X_v3) <- v3_samples
  v3_ri <- ri_full[match(v3_samples, sample)]

  v3_pdf <- file.path(plot_dir, "heatmap_ordering_V3_u_hclust.pdf")
  v3_png <- file.path(plot_dir, "heatmap_ordering_V3_u_hclust.png")
  render_one_page(
    X_v3, v3_ri, shown_positions,
    sprintf("V3 — u + HCLUST WITHIN GROUP (experimental)  |  cid=%d  %s:%s-%s",
            cid, chr, format(c_start, big.mark=","), format(c_end, big.mark=",")),
    "Experimental: u primary, hierarchical clustering secondary. Intended to surface internal substructure invisible to u alone.",
    v3_pdf, v3_png
  )

  # ═══════════════════════════════════════════════════════════════════════
  # Combined 3-page PDF (one version per page, FULL SIZE)
  # Nothing is shrunk — each page is the same size as its individual file.
  # This is the primary combined output since small side-by-side panels
  # make heatmap details unreadable for heavy candidates.
  # ═══════════════════════════════════════════════════════════════════════
  combined_pdf <- file.path(plot_dir, "heatmap_ordering_comparison_3page.pdf")
  combined_ok <- FALSE

  # Try qpdf merge first (fast, no re-rendering)
  qpdf_bin <- Sys.which("qpdf")
  if (nzchar(qpdf_bin)) {
    individual <- c(v1_pdf, v2_pdf, v3_pdf)
    existing <- individual[file.exists(individual)]
    if (length(existing) >= 1) {
      cmd <- sprintf("%s --empty --pages %s -- %s",
                     qpdf_bin, paste(shQuote(existing), collapse = " "),
                     shQuote(combined_pdf))
      ret <- system(cmd, intern = FALSE)
      if (ret == 0 && file.exists(combined_pdf)) {
        message("  3-page PDF via qpdf -> ", basename(combined_pdf))
        combined_ok <- TRUE
      }
    }
  }

  # Fallback: render directly to a multi-page PDF by re-drawing
  if (!combined_ok) {
    message("  qpdf unavailable or failed — rendering combined PDF via grDevices")
    # Find the largest page size from the three individual renders
    w <- max(7, min(16, 3 + ncol(X_v1) * 0.04))
    h <- max(6, min(18, 3 + nrow(X_v1) * 0.035))
    pdf(combined_pdf, width = w, height = h, onefile = TRUE)
    tryCatch({
      col_fun <- circlize::colorRamp2(c(0, 1, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
      for (vv in list(
        list(X = X_v1, ri = v1_ri,
             title = sprintf("V1 — RAW PC1 (NO ROTATION)  |  cid=%d  %s:%s-%s",
                             cid, chr, format(c_start, big.mark=","),
                             format(c_end, big.mark=",")),
             caption = "DEBUG ONLY — naive ordering, not the published method."),
        list(X = X_v2, ri = v2_ri,
             title = sprintf("V2 — u ROTATION ONLY (STEP21 method)  |  cid=%d  %s:%s-%s",
                             cid, chr, format(c_start, big.mark=","),
                             format(c_end, big.mark=",")),
             caption = "Established method: group primary, u secondary."),
        list(X = X_v3, ri = v3_ri,
             title = sprintf("V3 — u + HCLUST WITHIN GROUP  |  cid=%d  %s:%s-%s",
                             cid, chr, format(c_start, big.mark=","),
                             format(c_end, big.mark=",")),
             caption = "Experimental: u primary, hclust secondary.")
      )) {
        col_map <- GROUP_COLORS[intersect(names(GROUP_COLORS), unique(vv$ri$coarse_group))]
        ha_left <- ComplexHeatmap::rowAnnotation(
          Group = vv$ri$coarse_group,
          col = list(Group = col_map),
          show_annotation_name = TRUE,
          annotation_name_gp = grid::gpar(fontsize = 7),
          width = grid::unit(4, "mm")
        )
        pos_mb <- shown_positions / 1e6
        ha_top <- ComplexHeatmap::columnAnnotation(
          pos_mb = ComplexHeatmap::anno_points(
            pos_mb, pch = 16, size = grid::unit(1, "mm"),
            gp = grid::gpar(col = "#555e69"),
            axis_param = list(gp = grid::gpar(fontsize = 6))
          ),
          height = grid::unit(6, "mm"),
          annotation_name_gp = grid::gpar(fontsize = 7)
        )
        ht <- ComplexHeatmap::Heatmap(
          vv$X, name = "Dosage", col = col_fun,
          cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_names = FALSE, show_column_names = FALSE,
          left_annotation = ha_left,
          top_annotation = ha_top,
          column_title = vv$title,
          column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
          use_raster = nrow(vv$X) * ncol(vv$X) > 50000,
          raster_quality = 3,
          na_col = "#ffffff"
        )
        ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 2, 2, 2), "mm"))
        grid::grid.text(vv$caption, x = 0.01, y = 0.01,
                        just = c("left", "bottom"),
                        gp = grid::gpar(fontsize = 8, col = "#8a2b30",
                                        fontface = "italic"))
      }
    }, error = function(e) message("  [WARN] combined PDF failed: ", conditionMessage(e)))
    try(dev.off(), silent = TRUE)
    message("  3-page PDF rendered directly -> ", basename(combined_pdf))
  }

  # ═══════════════════════════════════════════════════════════════════════
  # Side-by-side PDF (all 3 versions on one wide page)
  # Works best for quick visual comparison. On composite candidates with
  # lots of samples, the individual files are still the primary reference
  # — the side-by-side is for overview only.
  # ═══════════════════════════════════════════════════════════════════════
  sidebyside_pdf <- file.path(plot_dir, "heatmap_ordering_comparison_side_by_side.pdf")
  sidebyside_png <- file.path(plot_dir, "heatmap_ordering_comparison_side_by_side.png")
  tryCatch({
    col_fun_sbs <- circlize::colorRamp2(c(0, 1, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
    per_panel_w <- max(5, min(9, 2.5 + ncol(X_v1) * 0.025))
    total_w <- per_panel_w * 3 + 1.5
    total_h <- max(6, min(18, 3 + nrow(X_v1) * 0.035))

    build_ht_sbs <- function(X_dat, ri, title, cap, show_name = TRUE) {
      col_map <- GROUP_COLORS[intersect(names(GROUP_COLORS), unique(ri$coarse_group))]
      ha_left <- ComplexHeatmap::rowAnnotation(
        Group = ri$coarse_group,
        col = list(Group = col_map),
        show_annotation_name = show_name,
        annotation_name_gp = grid::gpar(fontsize = 6),
        width = grid::unit(3, "mm")
      )
      pos_mb <- shown_positions / 1e6
      ha_top <- ComplexHeatmap::columnAnnotation(
        pos_mb = ComplexHeatmap::anno_points(
          pos_mb, pch = 16, size = grid::unit(0.8, "mm"),
          gp = grid::gpar(col = "#555e69"),
          axis_param = list(gp = grid::gpar(fontsize = 5))
        ),
        height = grid::unit(5, "mm"),
        show_annotation_name = show_name,
        annotation_name_gp = grid::gpar(fontsize = 6)
      )
      ComplexHeatmap::Heatmap(
        X_dat, name = "Dosage", col = col_fun_sbs,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,
        left_annotation = ha_left,
        top_annotation = ha_top,
        column_title = title,
        column_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
        use_raster = nrow(X_dat) * ncol(X_dat) > 50000,
        raster_quality = 3,
        na_col = "#ffffff",
        show_heatmap_legend = show_name
      )
    }

    ht_v1 <- build_ht_sbs(
      X_v1, v1_ri,
      sprintf("V1 RAW PC1\ncid=%d (DEBUG)", cid),
      "not published",
      show_name = TRUE
    )
    ht_v2 <- build_ht_sbs(
      X_v2, v2_ri,
      sprintf("V2 u ROTATION\n(STEP21)"),
      "baseline",
      show_name = FALSE
    )
    ht_v3 <- build_ht_sbs(
      X_v3, v3_ri,
      sprintf("V3 u + HCLUST\n(experimental)"),
      "new",
      show_name = FALSE
    )
    ht_combined <- ht_v1 + ht_v2 + ht_v3

    sbs_title <- sprintf("Ordering comparison  |  cid=%d  %s:%s-%s",
                          cid, chr,
                          format(c_start, big.mark=","), format(c_end, big.mark=","))
    pdf(sidebyside_pdf, width = total_w, height = total_h)
    ComplexHeatmap::draw(ht_combined, column_title = sbs_title,
                          column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
                          padding = grid::unit(c(2, 2, 4, 2), "mm"))
    grid::grid.text(
      "V1 = debug baseline (naive PC1). V2 = established method. V3 = experimental hclust.",
      x = 0.01, y = 0.01, just = c("left", "bottom"),
      gp = grid::gpar(fontsize = 7, col = "#666666", fontface = "italic")
    )
    dev.off()

    png(sidebyside_png, width = total_w, height = total_h, units = "in", res = 150)
    ComplexHeatmap::draw(ht_combined, column_title = sbs_title,
                          column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
                          padding = grid::unit(c(2, 2, 4, 2), "mm"))
    grid::grid.text(
      "V1 = debug baseline (naive PC1). V2 = established method. V3 = experimental hclust.",
      x = 0.01, y = 0.01, just = c("left", "bottom"),
      gp = grid::gpar(fontsize = 7, col = "#666666", fontface = "italic")
    )
    dev.off()

    message("  Side-by-side PDF -> ", basename(sidebyside_pdf))
    message("  Side-by-side PNG -> ", basename(sidebyside_png))
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("  [WARN] side-by-side render failed: ", conditionMessage(e))
  })

  # ═══════════════════════════════════════════════════════════════════════
  # Audit TSV: which sample appears at which row in each version
  # ═══════════════════════════════════════════════════════════════════════
  audit <- data.table(
    sample = sample_cols,
    coarse_group = ri_full$coarse_group,
    u = round(ri_full$u, 6),
    v = round(ri_full$v, 6),
    raw_pc1 = round(raw_pc1, 6),
    rank_v1_raw_pc1 = match(sample_cols, v1_samples),
    rank_v2_u_only  = match(sample_cols, v2_samples),
    rank_v3_u_hclust = match(sample_cols, v3_samples)
  )
  # Add "reordering distance" — how many positions V3 moves a sample vs V2
  audit[, shift_v3_vs_v2 := rank_v3_u_hclust - rank_v2_u_only]
  audit[, shift_v1_vs_v2 := rank_v1_raw_pc1 - rank_v2_u_only]
  fwrite(audit, file.path(cand_dir, "candidate_ordering_audit.tsv"), sep = "\t")

  # Also save the V3 hclust audit from the helper
  fwrite(v3_res$within_group_audit,
         file.path(cand_dir, "candidate_within_group_ordering.tsv"),
         sep = "\t")

  # ═══════════════════════════════════════════════════════════════════════
  # Console dump
  # ═══════════════════════════════════════════════════════════════════════
  max_v3_shift <- max(abs(audit$shift_v3_vs_v2), na.rm = TRUE)
  mean_v3_shift <- mean(abs(audit$shift_v3_vs_v2), na.rm = TRUE)
  max_v1_shift <- max(abs(audit$shift_v1_vs_v2), na.rm = TRUE)
  message(sprintf("  Audit: V3 vs V2 shift — max=%d, mean=%.1f  |  V1 vs V2 shift — max=%d",
                  max_v3_shift, mean_v3_shift, max_v1_shift))
  message("  Outputs in: ", plot_dir)
  message("    Individual full-size (PDF + PNG each):")
  message("      ", basename(v1_pdf))
  message("      ", basename(v2_pdf))
  message("      ", basename(v3_pdf))
  message("    Combined:")
  message("      heatmap_ordering_comparison_3page.pdf          (3 pages, full size each)")
  message("      heatmap_ordering_comparison_side_by_side.pdf   (1 page, 3 panels)")
  message("      heatmap_ordering_comparison_side_by_side.png")
  message("    Audit:")
  message("      candidate_ordering_audit.tsv           (per-sample ranks V1/V2/V3)")
  message("      candidate_within_group_ordering.tsv    (hclust audit from V3)")
}

message("[DONE] STEP36 ordering-comparison heatmaps complete")
