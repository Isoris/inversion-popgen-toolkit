#!/usr/bin/env Rscript
# =============================================================================
# q08_shelf_heatmap.R  v2
#
# Sample × SNP genotype (dosage) heatmaps for the shelf region, rendered at 9
# different marker-selection schemes (3 schemes × 3 scales):
#
#   Scheme A (window-based): 1 top-variance SNP per local-PCA window
#     A1 = 100-SNP windows
#     A5 = 500-SNP windows
#     A10 = 1000-SNP windows
#
#   Scheme B (count-based, from STEP28/35 tradition): top-N variance markers
#     B1 = top 30
#     B5 = top 100
#     B10 = top 500
#
#   Scheme C (physical-interval): 1 top-variance SNP per fixed bp window
#     C1  = 5 kb
#     C5  = 25 kb
#     C10 = 100 kb
#
# Each panel uses ComplexHeatmap when the package is available (with row
# annotations for Hom1/Het/Hom2 assignment + column annotation for genomic
# position), and falls back to ggplot geom_raster when it is not.
#
# All 9 panels are combined into a single PDF:
#   ${QC_FIGS}/shelf_heatmap.<CHR>.pdf
#
# Individual PNG/PDF copies of each panel are also saved under:
#   ${QC_FIGS}/shelf_heatmap.<CHR>/
#
# Ported from STEP28/STEP35 conventions (top-var selection, coolwarm dosage
# palette 0/1/2 = blue/white/red, row annotations for invgt class).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
BEAGLE      <- get_arg("--beagle")
POS_FILE    <- get_arg("--pos")
PRECOMP     <- get_arg("--precomp")
CHROM       <- get_arg("--chrom")
SHELF_A     <- as.numeric(get_arg("--shelf_start_mb"))
SHELF_B     <- as.numeric(get_arg("--shelf_end_mb"))
BP1_MB      <- as.numeric(get_arg("--breakpoint1_mb", NA))
BP2_MB      <- as.numeric(get_arg("--breakpoint2_mb", NA))
SAMPLE_LIST <- get_arg("--sample_list")
INVGT_ASSIGN_FILE <- get_arg("--invgt_assign", NULL)
OUT         <- get_arg("--out")
OUT_DIR     <- get_arg("--out_dir", NULL)  # directory for individual panels

# Back-compat: old --snp_stride argument ignored in v2
invisible(get_arg("--snp_stride", NULL))

stopifnot(file.exists(BEAGLE), file.exists(POS_FILE), file.exists(PRECOMP))
if (is.null(OUT_DIR)) OUT_DIR <- sub("\\.pdf$", "", OUT)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

has_CH <- suppressWarnings(require(ComplexHeatmap, quietly = TRUE))
has_circlize <- suppressWarnings(require(circlize, quietly = TRUE))
has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))
render_mode <- if (has_CH && has_circlize) "complexheatmap" else "ggplot"
message("[q08] Render mode: ", render_mode)

# =============================================================================
# 1. Parse positions -> find shelf SNPs
# =============================================================================
message("[q08] Loading positions: ", POS_FILE)
pos <- fread(POS_FILE, header = FALSE)
pos_col <- if (ncol(pos) >= 2) 2L else 1L
pos_vec <- as.integer(pos[[pos_col]])
shelf_mask <- pos_vec >= SHELF_A * 1e6 & pos_vec <= SHELF_B * 1e6
shelf_row_idx <- which(shelf_mask)
shelf_pos <- pos_vec[shelf_row_idx]
n_shelf <- length(shelf_row_idx)
if (n_shelf < 30) stop("[q08] Only ", n_shelf, " SNPs in shelf — too few for any heatmap")
message("[q08] Shelf SNPs: ", n_shelf)

# =============================================================================
# 2. Read sample master list, stream BEAGLE for full shelf dosage matrix
# =============================================================================
sample_master <- readLines(SAMPLE_LIST)
n_samples <- length(sample_master)
message("[q08] Master samples: ", n_samples)

# Stream BEAGLE: we want shelf_row_idx rows, all columns -> dosage matrix
# Output: dosage_mat[shelf_snp x sample]
message("[q08] Streaming BEAGLE for ", n_shelf, " shelf SNPs...")
dosage_mat <- matrix(NA_real_, nrow = n_shelf, ncol = n_samples)
rownames(dosage_mat) <- as.character(shelf_pos)
colnames(dosage_mat) <- sample_master

con <- gzfile(BEAGLE, "r")
header <- readLines(con, n = 1)
header_fields <- strsplit(header, "\t", fixed = TRUE)[[1]]
n_header_samples <- (length(header_fields) - 3) / 3
if (n_header_samples != n_samples) {
  message("[q08] WARN: header says ", n_header_samples,
          " samples, master list says ", n_samples, ". Using header count.")
  n_samples <- as.integer(n_header_samples)
  dosage_mat <- matrix(NA_real_, nrow = n_shelf, ncol = n_samples)
  rownames(dosage_mat) <- as.character(shelf_pos)
  colnames(dosage_mat) <- if (length(sample_master) >= n_samples)
                           sample_master[seq_len(n_samples)]
                         else paste0("Ind", seq_len(n_samples) - 1L)
}

# Build an environment-backed hash so missing keys return NULL cleanly
# (named-integer vector subset with [[key]] throws "subscript out of bounds")
target_env <- new.env(hash = TRUE, size = length(shelf_row_idx) * 2L)
for (i in seq_along(shelf_row_idx)) {
  assign(as.character(shelf_row_idx[i]), i, envir = target_env)
}

data_row_idx <- 0L
hits <- 0L
chunk_size <- 20000L
max_target <- max(shelf_row_idx)
repeat {
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0L) break
  for (ln in lines) {
    data_row_idx <- data_row_idx + 1L
    key <- as.character(data_row_idx)
    hi <- target_env[[key]]
    if (!is.null(hi)) {
      flds <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      gl_block <- as.numeric(flds[-(1:3)])
      gl_mat <- matrix(gl_block, nrow = n_samples, byrow = TRUE)
      dosage_mat[hi, ] <- gl_mat[, 2] + 2 * gl_mat[, 3]
      hits <- hits + 1L
      if (hits %% 2000L == 0L) message("  ", hits, " / ", n_shelf)
    }
    if (data_row_idx > max_target) break
  }
  if (data_row_idx > max_target) break
}
close(con)
message("[q08] Extracted ", hits, " / ", n_shelf, " shelf rows")

# Pre-compute variance per SNP once (used by all 3 schemes)
snp_var <- apply(dosage_mat, 1L, var, na.rm = TRUE)
snp_var[!is.finite(snp_var)] <- 0

# =============================================================================
# 3. Sample ordering and class assignment
# =============================================================================
sample_order <- NULL
sample_class <- NULL   # named character Hom1/Het/Hom2

if (!is.null(INVGT_ASSIGN_FILE) && file.exists(INVGT_ASSIGN_FILE)) {
  message("[q08] Using Q07 invgt assignments: ", INVGT_ASSIGN_FILE)
  inv <- fread(INVGT_ASSIGN_FILE)
  if ("sample_master" %in% names(inv) && "pc1_mean" %in% names(inv)) {
    inv[, invgt := factor(invgt, levels = c("Hom1", "Het", "Hom2"))]
    setorder(inv, invgt, pc1_mean)
    sample_order <- inv$sample_master
    sample_class <- setNames(as.character(inv$invgt), inv$sample_master)
  }
}
if (is.null(sample_order)) {
  message("[q08] Falling back to PC1-mean ordering from precomp")
  pc <- readRDS(PRECOMP)
  dt <- as.data.table(pc$dt)
  dt[, mb := (start_bp + end_bp) / 2 / 1e6]
  shelf_mask_p <- dt$mb >= SHELF_A & dt$mb <= SHELF_B
  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  if (length(pc1_cols) > 0 && sum(shelf_mask_p) > 0) {
    pc1_means <- colMeans(as.matrix(dt[shelf_mask_p, ..pc1_cols]), na.rm = TRUE)
    pc_samples <- sub("^PC_1_", "", pc1_cols)
    if (all(grepl("^Ind[0-9]+$", pc_samples))) {
      ind_idx <- as.integer(sub("^Ind", "", pc_samples))
      if (length(sample_master) >= max(ind_idx) + 1L) {
        names(pc1_means) <- sample_master[ind_idx + 1L]
      } else {
        names(pc1_means) <- pc_samples
      }
    } else {
      names(pc1_means) <- pc_samples
    }
    sample_order <- names(sort(pc1_means))
  } else {
    sample_order <- colnames(dosage_mat)
  }
}
sample_order <- intersect(sample_order, colnames(dosage_mat))
message("[q08] Samples ordered: ", length(sample_order),
        if (!is.null(sample_class)) " (with invgt class)" else " (no class)")

# =============================================================================
# 4. Scheme definitions — 3 schemes × 3 scales = 9 configs
# =============================================================================
# Each returns the shelf-SNP indices (within dosage_mat) to display.

scheme_A <- function(window_snps) {
  # Top-var SNP per <window_snps>-SNP bin
  n_bins <- ceiling(n_shelf / window_snps)
  keep <- integer(0)
  for (b in seq_len(n_bins)) {
    lo <- (b - 1L) * window_snps + 1L
    hi <- min(n_shelf, b * window_snps)
    if (hi < lo) next
    seg_var <- snp_var[lo:hi]
    if (all(seg_var == 0)) next
    best_local <- which.max(seg_var)
    keep <- c(keep, lo + best_local - 1L)
  }
  sort(unique(keep))
}
scheme_B <- function(top_n) {
  ord <- order(snp_var, decreasing = TRUE)
  keep <- ord[seq_len(min(top_n, length(ord)))]
  keep <- keep[snp_var[keep] > 0]
  sort(unique(keep))
}
scheme_C <- function(bp_bin) {
  # 1 top-var SNP per bp-bin
  bin_id <- floor(shelf_pos / bp_bin)
  keep <- integer(0)
  for (b in unique(bin_id)) {
    idxs <- which(bin_id == b)
    if (!length(idxs)) next
    seg_var <- snp_var[idxs]
    if (all(seg_var == 0)) next
    keep <- c(keep, idxs[which.max(seg_var)])
  }
  sort(unique(keep))
}

configs <- list(
  list(id = "A1",  label = "1 top-var / 100-SNP window",  idx_fn = function() scheme_A(100)),
  list(id = "A5",  label = "1 top-var / 500-SNP window",  idx_fn = function() scheme_A(500)),
  list(id = "A10", label = "1 top-var / 1000-SNP window", idx_fn = function() scheme_A(1000)),
  list(id = "B1",  label = "top 30 var markers",          idx_fn = function() scheme_B(30)),
  list(id = "B5",  label = "top 100 var markers",         idx_fn = function() scheme_B(100)),
  list(id = "B10", label = "top 500 var markers",         idx_fn = function() scheme_B(500)),
  list(id = "C1",  label = "1 top-var / 5 kb",            idx_fn = function() scheme_C(5e3)),
  list(id = "C5",  label = "1 top-var / 25 kb",           idx_fn = function() scheme_C(25e3)),
  list(id = "C10", label = "1 top-var / 100 kb",          idx_fn = function() scheme_C(100e3))
)

# =============================================================================
# 5. Render helpers
# =============================================================================

# Palette from STEP28/35: blue / white / red dosage 0 / 1 / 2
dosage_cols <- c("#2166AC", "#F7F7F7", "#B2182B")
class_cols  <- c("Hom1" = "#1f4e79", "Het" = "#c0504d", "Hom2" = "#2c7a39")

render_ch <- function(X, positions, title, pdf_path, png_path) {
  col_fun <- circlize::colorRamp2(c(0, 1, 2), dosage_cols)

  # Row annotation (invgt class)
  ha_left <- NULL
  if (!is.null(sample_class)) {
    class_v <- sample_class[rownames(X)]
    class_v[is.na(class_v)] <- "unknown"
    col_map <- c(class_cols, "unknown" = "#CCCCCC")
    col_map <- col_map[intersect(names(col_map), unique(class_v))]
    ha_left <- ComplexHeatmap::rowAnnotation(
      invgt = class_v,
      col = list(invgt = col_map),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 7),
      width = grid::unit(4, "mm")
    )
  }

  # Column annotation: genomic position (as numeric bar) + optional breakpoints
  pos_mb <- positions / 1e6
  annot_list <- list(
    pos_mb = ComplexHeatmap::anno_points(
      pos_mb, pch = 16, size = grid::unit(1, "mm"),
      gp = grid::gpar(col = "#555e69"),
      axis_param = list(gp = grid::gpar(fontsize = 6))
    )
  )
  # If breakpoints provided, add labeled mark annotation pointing at nearest windows
  bp_keep <- c()
  bp_labels <- character()
  for (bp in c(BP1_MB, BP2_MB)) {
    if (is.finite(bp)) {
      j <- which.min(abs(pos_mb - bp))
      if (length(j) == 1L) {
        bp_keep <- c(bp_keep, j)
        bp_labels <- c(bp_labels, sprintf("BP %.3f Mb", bp))
      }
    }
  }
  if (length(bp_keep) > 0) {
    annot_list$breakpoint <- ComplexHeatmap::anno_mark(
      at = bp_keep, labels = bp_labels, which = "column", side = "top",
      labels_gp = grid::gpar(fontsize = 7, col = "#8a2b30"),
      lines_gp  = grid::gpar(col = "#e0555c", lwd = 1)
    )
  }
  ha_top <- do.call(ComplexHeatmap::columnAnnotation,
                    c(annot_list,
                      list(height = grid::unit(6, "mm"),
                           annotation_name_gp = grid::gpar(fontsize = 7))))

  n_cell <- nrow(X) * ncol(X)
  w <- max(6, min(16, 3 + ncol(X) * 0.06))
  h <- max(5, min(18, 3 + nrow(X) * 0.035))

  ht <- ComplexHeatmap::Heatmap(
    X, name = "Dosage", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    left_annotation = ha_left,
    top_annotation = ha_top,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    use_raster = n_cell > 50000,
    raster_quality = 3,
    na_col = "#ffffff"
  )

  tryCatch({
    pdf(pdf_path, width = w, height = h); ComplexHeatmap::draw(ht); dev.off()
    png(png_path, width = w, height = h, units = "in", res = 150)
    ComplexHeatmap::draw(ht); dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("  [WARN] CH render failed for ", basename(pdf_path), ": ",
            conditionMessage(e))
  })
}

render_ggplot <- function(X, positions, title, pdf_path, png_path) {
  df <- as.data.table(reshape2::melt(X, varnames = c("sample", "snp"),
                                      value.name = "dosage"))
  df[, sample := factor(sample, levels = rownames(X))]
  df[, snp := factor(snp, levels = colnames(X))]
  df[, gt := fcase(
    is.na(dosage), "NA",
    dosage < 0.5, "0/0",
    dosage < 1.5, "0/1",
    default = "1/1"
  )]
  gt_pal <- c("0/0" = dosage_cols[1], "0/1" = "#d0d0d0",
              "1/1" = dosage_cols[3], "NA" = "#ffffff")
  df[, gt := factor(gt, levels = names(gt_pal))]

  p <- ggplot(df, aes(snp, sample, fill = gt)) +
    geom_raster() +
    scale_fill_manual(values = gt_pal, name = "gt", na.value = "#ffffff") +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.key.width = unit(0.25, "cm"),
          legend.key.height = unit(0.35, "cm"))

  if (!is.null(sample_class)) {
    class_v <- sample_class[rownames(X)]
    strip_df <- data.table(sample = factor(rownames(X), levels = rownames(X)),
                            class = class_v)
    p_strip <- ggplot(strip_df, aes(x = 1, y = sample, fill = class)) +
      geom_raster() +
      scale_fill_manual(values = class_cols, na.value = "#CCCCCC") +
      theme_void() + theme(legend.position = "none")
    if (has_patchwork) {
      p <- (p_strip | p) + patchwork::plot_layout(widths = c(0.03, 1))
    }
  }

  w <- max(6, min(16, 3 + ncol(X) * 0.06))
  h <- max(5, min(18, 3 + nrow(X) * 0.035))
  ggsave(pdf_path, p, width = w, height = h, device = cairo_pdf)
  ggsave(png_path, p, width = w, height = h, dpi = 150)
}

render_empty <- function(title, pdf_path, png_path, reason) {
  p <- ggplot() + labs(title = title, subtitle = paste0("(", reason, ")")) +
    theme_void() + theme(plot.title = element_text(size = 9))
  ggsave(pdf_path, p, width = 6, height = 5, device = cairo_pdf)
  ggsave(png_path, p, width = 6, height = 5, dpi = 100)
}

# =============================================================================
# 6. Generate all 9 panels
# =============================================================================
panel_files <- character(length(configs))
summary_df <- data.table(id = character(), label = character(),
                          n_markers = integer(), ok = logical())

for (ci in seq_along(configs)) {
  cf <- configs[[ci]]
  message(sprintf("[q08] Rendering %s: %s", cf$id, cf$label))

  pdf_path <- file.path(OUT_DIR, paste0("panel_", cf$id, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0("panel_", cf$id, ".png"))
  panel_files[ci] <- pdf_path

  idx <- tryCatch(cf$idx_fn(), error = function(e) integer(0))
  n_markers <- length(idx)

  if (n_markers < 3) {
    message("  skip: only ", n_markers, " markers")
    render_empty(paste0(cf$id, ": ", cf$label), pdf_path, png_path,
                  paste0("only ", n_markers, " markers available"))
    summary_df <- rbind(summary_df, data.table(id = cf$id, label = cf$label,
                                                  n_markers = n_markers, ok = FALSE))
    next
  }

  # Build sample × SNP matrix in display orientation (rows = samples, cols = SNPs)
  X_shelf <- dosage_mat[idx, , drop = FALSE]     # [SNPs × samples]
  X <- t(X_shelf)[sample_order, , drop = FALSE]   # [samples × SNPs]
  colnames(X) <- as.character(shelf_pos[idx])

  title <- sprintf("%s: %s  (%d SNPs × %d samples)",
                   cf$id, cf$label, n_markers, nrow(X))

  if (render_mode == "complexheatmap") {
    render_ch(X, shelf_pos[idx], title, pdf_path, png_path)
  } else {
    if (!requireNamespace("reshape2", quietly = TRUE)) {
      # emergency fallback: just write a stub
      render_empty(title, pdf_path, png_path, "reshape2 missing for ggplot path")
    } else {
      render_ggplot(X, shelf_pos[idx], title, pdf_path, png_path)
    }
  }
  summary_df <- rbind(summary_df, data.table(id = cf$id, label = cf$label,
                                                n_markers = n_markers, ok = TRUE))
}

# =============================================================================
# 7. Combined summary PDF (3 × 3 grid of the 9 PDFs)
# =============================================================================
# Simplest portable approach: one PDF page per panel, ordered A1..C10
# If qpdf is available we could merge with it; otherwise we use the
# individual files plus a lightweight index PDF.

tryCatch({
  # Try qpdf merge (usually available with Cairo/Ghostscript stack)
  qpdf_bin <- Sys.which("qpdf")
  if (nzchar(qpdf_bin)) {
    existing <- panel_files[file.exists(panel_files)]
    if (length(existing) > 1) {
      cmd <- sprintf("%s --empty --pages %s -- %s",
                     qpdf_bin, paste(shQuote(existing), collapse = " "),
                     shQuote(OUT))
      ret <- system(cmd, intern = FALSE)
      if (ret == 0 && file.exists(OUT)) {
        message("[q08] Combined PDF via qpdf -> ", OUT)
      } else {
        stop("qpdf returned ", ret)
      }
    } else {
      stop("no individual panels found")
    }
  } else {
    stop("qpdf not found")
  }
}, error = function(e) {
  message("[q08] qpdf unavailable (", conditionMessage(e),
          "); writing single-page composite via cairo_pdf")
  # Fallback: open one long PDF and onePage-per-panel via grDevices
  pdf(OUT, width = 11, height = 8, onefile = TRUE)
  on.exit(try(dev.off(), silent = TRUE))
  for (pf in panel_files) {
    if (!file.exists(pf)) next
    # We can't re-open a PDF as a grob easily; instead write a stub page
    # with the panel title + pointer to the file.
    grid::grid.newpage()
    grid::grid.text(paste0("see individual file: ", basename(pf)),
                    gp = grid::gpar(fontsize = 12))
  }
  on.exit(NULL)
  try(dev.off(), silent = TRUE)
  message("[q08] Fallback composite -> ", OUT)
  message("[q08] Individual panels are under ", OUT_DIR)
})

# Summary TSV
summary_tsv <- file.path(OUT_DIR, "panel_summary.tsv")
fwrite(summary_df, summary_tsv, sep = "\t", quote = FALSE)
message("[q08] Summary TSV -> ", summary_tsv)

# Dump the console info
cat("\n========= Q08 SHELF HEATMAPS =========\n")
cat(sprintf("Chromosome : %s\n", CHROM))
cat(sprintf("Shelf      : %.2f - %.2f Mb (%d SNPs total)\n", SHELF_A, SHELF_B, n_shelf))
cat(sprintf("Samples    : %d  %s\n", length(sample_order),
            if (!is.null(sample_class)) "(invgt-annotated)" else "(unannotated)"))
cat(sprintf("Render mode: %s\n", render_mode))
cat("Panels:\n")
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %-5s %-38s n_markers=%5d  %s\n",
              summary_df$id[i], summary_df$label[i],
              summary_df$n_markers[i],
              if (summary_df$ok[i]) "OK" else "EMPTY"))
}
cat(sprintf("Summary    : %s\n", summary_tsv))
cat(sprintf("Combined   : %s\n", OUT))
cat(sprintf("Panels dir : %s\n", OUT_DIR))
cat("=======================================\n")
message("[q08] DONE")
