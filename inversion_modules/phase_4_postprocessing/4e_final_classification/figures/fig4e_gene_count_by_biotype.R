#!/usr/bin/env Rscript
# =============================================================================
# fig4e_gene_count_by_biotype.R — table-heatmap of gene counts per biotype × region
# =============================================================================
#
# Purpose
# -------
# For a user-supplied set of genomic regions (ROHs, inversions, or any
# BED-like set), count how many genes of each biotype (protein_coding,
# pseudogene, lncRNA, miRNA, etc.) fall in each region, and render it
# as a table-heatmap: rows = biotypes, columns = regions, cells shaded
# by count with the count printed inside.
#
# Inspired by the Nature paper's "Gene count in ROH on Chr 6 vs Chr 12"
# side-table (small vertical table with a blue color ramp). Same idea
# but generalized to N regions.
#
# Inputs
# ------
#   --gene_map <tsv>    Per-gene annotation, from MODULE_CONSERVATION
#                       STEP 00. Required columns:
#                         gene_id, chr, gene_start, gene_end, biotype
#                       Extra columns (strand, cds_length, is_canonical)
#                       are ignored. Multiple rows per gene (one per
#                       transcript) are de-duplicated on gene_id before
#                       counting — otherwise a gene with 5 transcripts
#                       counts 5x.
#
#   --regions <tsv>     Regions of interest. Required columns:
#                         region_id, chrom, start, end
#                       Optional columns:
#                         category {"ROH" | "inversion" | ...}
#                           Used to cluster columns visually if present.
#                         label  (short display name; defaults to
#                           region_id, e.g. "ROH Chr 6\n52.9 Mb")
#
#   --out <file>        Output PDF path. PNG is also written alongside.
#
# Optional knobs
# --------------
#   --biotype_order <csv>   Explicit biotype display order, top→bottom.
#                           Default: data-driven, sorted by total count
#                           across all regions (most abundant on top,
#                           matching the Nature figure).
#
#   --top_biotypes <int>    Cap display to top N biotypes by total count.
#                           Others lumped into "Other". Default: show all.
#
#   --min_count <int>       Cells with count below this are shown blank
#                           (white bg, no number). Default 1 (show any
#                           non-zero). Set to 0 to also show explicit
#                           zeros as faint cells.
#
#   --color_scale linear|log   Color intensity rule. Default linear.
#                              For ranges spanning 1 to hundreds, log
#                              compresses the top and makes small counts
#                              visible.
#
# Outputs
# -------
#   <out>.pdf                  The figure.
#   <out>.png                  Same, 200 dpi.
#   <out>.counts.tsv           The underlying biotype × region count
#                              matrix (wide format). Written so you can
#                              rebuild the figure later in any tool.
#
# Implementation notes
# --------------------
# Counts are done in-memory with data.table, not via bedtools, to avoid
# an external dependency. For ~40k genes × ~30 regions this is trivially
# fast. If region counts grow past a few hundred, consider precomputing
# overlaps once with foverlaps() and caching.
#
# Uses ComplexHeatmap when available (preferred — handles the cell text
# and row/col annotations cleanly), falls back to a base-grid rendering
# if not.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

has_CH <- suppressWarnings(require(ComplexHeatmap, quietly = TRUE))
has_circlize <- suppressWarnings(require(circlize, quietly = TRUE))
has_grid <- requireNamespace("grid", quietly = TRUE)

# ── CLI ──
args <- commandArgs(trailingOnly = TRUE)
gene_map_path   <- NULL
regions_path    <- NULL
out_path        <- "gene_biotype_table.pdf"
biotype_order   <- NULL   # character vector, NULL = data-driven
top_biotypes    <- 0L     # 0 = show all
min_count       <- 1L
color_scale     <- "linear"

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--gene_map"      = { gene_map_path <- args[i+1]; i <- i+2 },
    "--regions"       = { regions_path  <- args[i+1]; i <- i+2 },
    "--out"           = { out_path      <- args[i+1]; i <- i+2 },
    "--biotype_order" = { biotype_order <- strsplit(args[i+1], ",")[[1]];
                          i <- i+2 },
    "--top_biotypes"  = { top_biotypes  <- as.integer(args[i+1]); i <- i+2 },
    "--min_count"     = { min_count     <- as.integer(args[i+1]); i <- i+2 },
    "--color_scale"   = { color_scale   <- args[i+1]; i <- i+2 },
    {
      message("[gene_biotype_table] unknown arg: ", args[i]); i <- i+1
    }
  )
}

stopifnot(!is.null(gene_map_path), !is.null(regions_path))
if (!color_scale %in% c("linear", "log")) {
  stop("[gene_biotype_table] --color_scale must be 'linear' or 'log'")
}

# ── Load + validate gene map ──
gmap <- fread(gene_map_path)
required_gm <- c("gene_id", "chr", "gene_start", "gene_end", "biotype")
missing_gm <- setdiff(required_gm, names(gmap))
if (length(missing_gm) > 0) {
  stop("[gene_biotype_table] gene_map missing columns: ",
       paste(missing_gm, collapse = ", "))
}
# De-duplicate by gene_id (keep first row — the transcript identity doesn't
# matter for biotype counting, and the gene_start/end from the first
# transcript is a reasonable span).
gmap <- unique(gmap, by = "gene_id")
gmap[, biotype := ifelse(is.na(biotype) | biotype == "", "unknown", biotype)]
setkey(gmap, chr, gene_start, gene_end)

# ── Load + validate regions ──
regions <- fread(regions_path)
required_r <- c("region_id", "chrom", "start", "end")
missing_r <- setdiff(required_r, names(regions))
if (length(missing_r) > 0) {
  stop("[gene_biotype_table] regions missing columns: ",
       paste(missing_r, collapse = ", "))
}
if (!"label" %in% names(regions)) {
  # Auto-label: "ROH LG6 52.9 Mb" if category present, else just region_id
  if ("category" %in% names(regions)) {
    regions[, label := sprintf("%s\n%s\n%.1f Mb",
                               category, chrom, (end - start) / 1e6)]
  } else {
    regions[, label := region_id]
  }
}
if (!"category" %in% names(regions)) regions[, category := ""]

# ── Count genes per biotype × region ──
# One pass: for each region, subset gmap to overlapping genes and tabulate.
# "Overlap" = gene body intersects region by at least 1 bp.
count_mat_list <- list()
for (ri in seq_len(nrow(regions))) {
  r <- regions[ri]
  sub <- gmap[chr == r$chrom & gene_end >= r$start & gene_start <= r$end]
  tbl <- sub[, .N, by = biotype]
  setnames(tbl, "N", r$region_id)
  count_mat_list[[r$region_id]] <- tbl
}

# Merge all per-region counts into a wide biotype × region matrix
all_biotypes <- sort(unique(gmap$biotype))
count_mat <- matrix(0L,
                    nrow = length(all_biotypes),
                    ncol = nrow(regions),
                    dimnames = list(all_biotypes, regions$region_id))
for (rid in regions$region_id) {
  tbl <- count_mat_list[[rid]]
  if (nrow(tbl) == 0) next
  # tbl has 2 cols: biotype, <rid>
  count_mat[tbl$biotype, rid] <- tbl[[2]]
}

# ── Reorder biotypes ──
# Priority 1: user-supplied order (error on unknown entries).
# Priority 2: data-driven, by total count across regions, descending.
if (!is.null(biotype_order)) {
  unknown <- setdiff(biotype_order, rownames(count_mat))
  if (length(unknown) > 0) {
    message("[gene_biotype_table] biotype_order has unknown entries: ",
            paste(unknown, collapse = ", "), "  (dropped)")
    biotype_order <- intersect(biotype_order, rownames(count_mat))
  }
  # Append unlisted biotypes at the end
  biotype_order <- c(biotype_order,
                     setdiff(rownames(count_mat), biotype_order))
  count_mat <- count_mat[biotype_order, , drop = FALSE]
} else {
  row_totals <- rowSums(count_mat)
  count_mat <- count_mat[order(-row_totals), , drop = FALSE]
}

# ── Collapse to top-N if requested ──
if (top_biotypes > 0 && nrow(count_mat) > top_biotypes) {
  top_rows <- seq_len(top_biotypes)
  other <- colSums(count_mat[-top_rows, , drop = FALSE])
  count_mat <- rbind(count_mat[top_rows, , drop = FALSE], Other = other)
}

# ── Write the underlying counts TSV ──
counts_out <- paste0(tools::file_path_sans_ext(out_path), ".counts.tsv")
count_df <- data.frame(biotype = rownames(count_mat), count_mat,
                       check.names = FALSE)
fwrite(count_df, counts_out, sep = "\t")
message("[gene_biotype_table] Wrote counts TSV: ", counts_out)

# ── Build the display matrix ──
# Cells below min_count render as blank (NA → no fill, no text). We still
# keep count_mat for the color scale and the color extremes.
display_mat <- count_mat
display_mat[display_mat < min_count] <- NA

# Max value for color normalization (exclude NAs)
max_val <- max(count_mat, na.rm = TRUE)
if (!is.finite(max_val) || max_val < 1) max_val <- 1

# Guard: if all regions contain <=1 gene, the 3-point color ramp
# degenerates (c(1, 0.5, 1) or c(1, 1, 1) — colorRamp2 refuses). Switch
# to a 2-point ramp in that case so the plot still renders.
use_midpoint <- max_val >= 2

# Guard: empty result (no genes overlap any region). Instead of erroring
# with a blank figure, write a warning placeholder PDF and exit.
if (nrow(count_mat) == 0 || sum(count_mat) == 0) {
  pdf(out_path, width = 6, height = 3)
  plot.new()
  text(0.5, 0.5,
       "No genes overlap any of the supplied regions.\nCheck region coordinates against the gene_map.",
       cex = 0.9)
  dev.off()
  message("[gene_biotype_table] No overlaps found; wrote placeholder: ",
          out_path)
  quit(save = "no", status = 0)
}

# ── Render ──
pdf_width  <- max(5, 0.8 * ncol(count_mat) + 3)
pdf_height <- max(4, 0.35 * nrow(count_mat) + 2)

if (has_CH && has_circlize) {
  message("[gene_biotype_table] ComplexHeatmap rendering")

  # Color ramp: pale teal → deep teal, matching the Nature figure style.
  # Linear vs log affects how the midpoint is placed. If max_val < 2 the
  # 3-point ramp degenerates, so fall back to a 2-point ramp.
  if (use_midpoint) {
    col_breaks <- if (color_scale == "log") {
      c(1, sqrt(max_val), max_val)
    } else {
      c(1, max_val / 2, max_val)
    }
    col_fun <- circlize::colorRamp2(
      col_breaks,
      c("#E5F5F0", "#66C2A4", "#006D2C")
    )
  } else {
    col_fun <- circlize::colorRamp2(
      c(0, max_val),
      c("#E5F5F0", "#66C2A4")
    )
  }

  # Cell text: print the count in the cell. White text on dark fills,
  # dark text on pale fills, to keep contrast readable.
  cell_fun <- function(j, i, x, y, width, height, fill) {
    v <- display_mat[i, j]
    if (is.na(v)) return(invisible(NULL))
    # Threshold for white-on-dark. Halfway up the color ramp is a safe cut.
    text_col <- if (v > max_val * 0.5) "white" else "grey15"
    grid::grid.text(format(v, big.mark = ",", scientific = FALSE),
                    x, y,
                    gp = grid::gpar(fontsize = 9, col = text_col,
                                    fontface = "plain"))
  }

  # Column labels: from the regions$label column (supports \n for multi-line)
  col_labels <- regions$label[match(colnames(count_mat), regions$region_id)]

  # Optional top annotation by category
  top_anno <- NULL
  if ("category" %in% names(regions) && any(nzchar(regions$category))) {
    cat_vec <- regions$category[match(colnames(count_mat), regions$region_id)]
    # Robust palette: cycle a base palette when more categories than colors.
    # Previous code indexed into a 5-color vector with seq_along(unique) and
    # silently produced NAs for categories 6+. Using recycle-safe approach.
    base_cols <- c("#4575B4", "#D6604D", "#7FBC41", "#984EA3",
                   "#999999", "#E7B800", "#00A2A3", "#B2182B")
    uc <- unique(cat_vec)
    cat_pal <- setNames(base_cols[((seq_along(uc) - 1) %% length(base_cols)) + 1],
                        uc)
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      Category = cat_vec,
      col = list(Category = cat_pal),
      show_legend = TRUE,
      annotation_name_side = "left",
      simple_anno_size = grid::unit(3, "mm")
    )
  }

  ht <- ComplexHeatmap::Heatmap(
    display_mat,
    name = "Gene count",
    col = col_fun,
    na_col = "white",
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 9),
    column_names_rot = 0,
    column_labels = col_labels,
    column_names_centered = TRUE,
    cell_fun = cell_fun,
    rect_gp = grid::gpar(col = "grey80", lwd = 0.4),
    top_annotation = top_anno,
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8),
      legend_height = grid::unit(3, "cm")
    )
  )

  pdf(out_path, width = pdf_width, height = pdf_height)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "top",
                       annotation_legend_side = "top")
  dev.off()

  png_path <- paste0(tools::file_path_sans_ext(out_path), ".png")
  png(png_path, width = pdf_width, height = pdf_height,
      units = "in", res = 200)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "top",
                       annotation_legend_side = "top")
  dev.off()

  message("[gene_biotype_table] Wrote: ", out_path)
  message("[gene_biotype_table] Wrote: ", png_path)
} else {
  # ── Fallback: base-R image() + text() ──
  message("[gene_biotype_table] ComplexHeatmap not available; base fallback")
  pdf(out_path, width = pdf_width, height = pdf_height)
  par(mar = c(5, 10, 2, 4), las = 1)
  # image() plots t(matrix), so transpose for the same orientation
  breaks <- seq(0, max_val, length.out = 101)
  pal <- colorRampPalette(c("#E5F5F0", "#66C2A4", "#006D2C"))(100)
  image(x = seq_len(ncol(count_mat)), y = seq_len(nrow(count_mat)),
        z = t(count_mat)[, nrow(count_mat):1],  # flip so top row is on top
        col = pal, breaks = breaks,
        axes = FALSE, xlab = "", ylab = "")
  axis(2, at = seq_len(nrow(count_mat)),
       labels = rev(rownames(count_mat)), tick = FALSE)
  col_labels <- regions$label[match(colnames(count_mat), regions$region_id)]
  axis(1, at = seq_len(ncol(count_mat)), labels = col_labels, tick = FALSE,
       cex.axis = 0.8)
  # Cell values
  for (i in seq_len(nrow(count_mat))) {
    for (j in seq_len(ncol(count_mat))) {
      v <- count_mat[i, j]
      if (is.na(v) || v < min_count) next
      y_pos <- nrow(count_mat) - i + 1   # flipped
      text_col <- if (v > max_val * 0.5) "white" else "grey15"
      text(j, y_pos, format(v, big.mark = ","),
           col = text_col, cex = 0.8)
    }
  }
  dev.off()
  message("[gene_biotype_table] Wrote: ", out_path)
}
