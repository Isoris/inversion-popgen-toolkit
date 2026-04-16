#!/usr/bin/env Rscript
# STEP28_candidate_marker_heatmap.R  v3 — hex legend fix + quality tier track
# FIX: Removes Position columnAnnotation (root cause of hex-code legend).
#      Uses anno_simple() for Het. Adds Quality tier + polarity tracks.

suppressPackageStartupMessages({ library(data.table) })
has_CH <- suppressPackageStartupMessages(suppressWarnings(require(ComplexHeatmap, quietly = TRUE)))
has_circlize <- suppressWarnings(require(circlize, quietly = TRUE))
if (!has_CH) { message("[WARN] ComplexHeatmap unavailable"); quit("no", status = 0) }

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_
source(config_file)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]; cid <- row$candidate_id; chr <- row$chrom
  c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)
  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  plot_dir <- file.path(PLOTS_DIR, paste0(chr, ".candidate_", cid)); ensure_dir(plot_dir)

  rot <- safe_fread(file.path(cand_dir, "candidate_pca_rotated.tsv"))
  coh <- safe_fread(file.path(cand_dir, "candidate_sample_coherence.tsv"))
  pol <- safe_fread(file.path(cand_dir, "candidate_marker_polarity.tsv"))
  if (is.null(rot) || nrow(rot) < 5) next

  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) next
  dos <- fread(dos_file); sites <- fread(sites_file)
  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 20) next
  dos_reg <- dos[keep]; sites_reg <- sites[keep]

  sc <- setdiff(names(dos_reg), "marker")
  if (all(grepl("^Ind", sc)) && length(sc) == nrow(rot))
    setnames(dos_reg, old = sc, new = rot$sample)
  sample_cols <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))
  if (length(sample_cols) < 5) next

  X <- as.matrix(dos_reg[, ..sample_cols]); storage.mode(X) <- "double"

  # Top-variance markers
  mv <- apply(X, 1, var, na.rm = TRUE)
  n_show <- min(200, nrow(X))
  top_idx <- sort(order(mv, decreasing = TRUE)[seq_len(n_show)])
  X_show <- X[top_idx, , drop = FALSE]

  # Order samples: group then u
  ri <- copy(rot[sample %in% sample_cols])
  if (!is.null(coh)) ri <- merge(ri, coh[, .(sample, stripe_quality)], by = "sample", all.x = TRUE)
  if (!"stripe_quality" %in% names(ri)) ri[, stripe_quality := "unknown"]
  ri[is.na(stripe_quality), stripe_quality := "unknown"]
  setorder(ri, coarse_group_refined, u)
  so <- intersect(ri$sample, colnames(X_show))
  X_ordered <- t(X_show[, so, drop = FALSE]); rownames(X_ordered) <- so
  ri <- ri[match(so, sample)]

  # ── Row annotations (NO hex codes) ────────────────────────────────────
  het_v <- ri$regional_het; het_v[!is.finite(het_v)] <- 0
  het_cf <- circlize::colorRamp2(c(0, 0.5, 1), c("#440154", "#21918C", "#FDE725"))
  grp_c <- safe_group_colors(ri$coarse_group_refined)
  tier_c <- c("core"="#1B7837","peripheral"="#F4A582","junk"="#D73027","unknown"="#BDBDBD")

  ha_row <- ComplexHeatmap::rowAnnotation(
    Group = ri$coarse_group_refined,
    Quality = ri$stripe_quality,
    Het = ComplexHeatmap::anno_simple(het_v, col = het_cf, na_col = "grey90"),
    col = list(Group = grp_c,
               Quality = tier_c[intersect(names(tier_c), unique(ri$stripe_quality))]),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 8)
  )

  # ── Column annotation: polarity only (NO Position → no hex) ───────────
  ha_top <- NULL
  if (!is.null(pol) && nrow(pol) >= nrow(X)) {
    pol_vals <- as.character(pol$final_flip_decision[top_idx])
    if (length(pol_vals) == ncol(X_ordered)) {
      ha_top <- ComplexHeatmap::columnAnnotation(
        Polarity = pol_vals,
        col = list(Polarity = c("FALSE" = "#2166AC", "TRUE" = "#D73027")),
        show_annotation_name = TRUE, annotation_name_gp = grid::gpar(fontsize = 7))
    }
  }

  header <- candidate_header(cid, chr, c_start, c_end, n_snps = n_show, n_samples = length(so))
  w <- max(8, n_show * 0.04 + 4); h <- max(6, length(so) * 0.06 + 3)

  tryCatch({
    pdf(file.path(plot_dir, "FIG_C08_marker_heatmap.pdf"), width = w, height = h)
    ht <- ComplexHeatmap::Heatmap(X_ordered, name = "Dosage",
      col = circlize::colorRamp2(c(0,1,2), c("#2166AC","#F7F7F7","#B2182B")),
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = FALSE, show_column_names = FALSE,
      left_annotation = ha_row, top_annotation = ha_top,
      column_title = paste0("FIG_C08 | ", wrap_subtitle(header)),
      column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      use_raster = nrow(X_ordered) * ncol(X_ordered) > 50000)
    ComplexHeatmap::draw(ht); dev.off()
    png(file.path(plot_dir, "FIG_C08_marker_heatmap.png"), width = w, height = h, units = "in", res = 300)
    ComplexHeatmap::draw(ht); dev.off()
    message("  Saved: FIG_C08_marker_heatmap")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("  [WARN] FIG_C08 failed: ", conditionMessage(e))
  })
}
message("[DONE] STEP28 v3 complete")
