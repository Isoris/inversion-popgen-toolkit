#!/usr/bin/env Rscript

# =============================================================================
# STEP18_composite_candidate_figure.R
#
# Creates a publication-quality composite figure per inversion candidate,
# inspired by Todesco et al. 2020 (Nature, sunflower haploblocks, Fig. 4).
#
# Layout:
#   Row 1: [Local PCA / MDS scan]  [Regional PCA colored by group]
#   Row 2: [LD triangle heatmap]   [Het boxplot by group]
#   Row 3: [Chromosome tP curves by group, candidate region shaded]
#
# Usage:
#   Rscript STEP18_composite_candidate_figure.R \
#     <candidate_id> <candidate_table> <step12_dir> <mds_rds> \
#     <dosage_dir> <outdir> \
#     [ld_pairs_file=NONE] [sample_order_file=NONE] [bin_bp=250000]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript STEP18_... <cid> <cand_table> <step12_dir> <mds_rds> <dosage_dir> <outdir> ...")
}

cid_val        <- args[1]
candidate_file <- args[2]
step12_dir     <- args[3]
mds_rds_file   <- args[4]
dosage_dir     <- args[5]
outdir         <- args[6]
ld_pairs_file  <- if (length(args) >= 7 && args[7] != "NONE") args[7] else NA_character_
sample_order   <- if (length(args) >= 8 && args[8] != "NONE") args[8] else NA_character_
bin_bp         <- if (length(args) >= 9) as.numeric(args[9]) else 250000

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Read candidate ─────────────────────────────────────────────────────────
cand <- fread(candidate_file)
row <- cand[as.character(candidate_id) == as.character(cid_val)]
if (nrow(row) != 1) stop("Candidate ", cid_val, " not found or ambiguous")

chr <- row$chrom[1]
inv_start <- as.numeric(row$start_bp[1])
inv_end   <- as.numeric(row$end_bp[1])

# ── Panel A: MDS scan ──────────────────────────────────────────────────────
p_mds <- ggplot() + theme_void() + labs(title = "(a) MDS outlier scan")

if (file.exists(mds_rds_file)) {
  mds_obj <- readRDS(mds_rds_file)
  mds_dt <- mds_obj$dt

  # Filter to this chromosome
  mds_chr <- mds_dt[chrom == chr]
  if (nrow(mds_chr) > 0) {
    # Use MDS1 as the primary scan axis
    mds_chr[, outlier := abs(MDS1_z) >= 3]
    mds_chr[, mid_bp := (start_bp + end_bp) / 2]

    p_mds <- ggplot(mds_chr, aes(x = mid_bp / 1e6, y = MDS1)) +
      geom_point(aes(color = outlier), size = 1.2, alpha = 0.7) +
      annotate("rect",
               xmin = inv_start / 1e6, xmax = inv_end / 1e6,
               ymin = -Inf, ymax = Inf,
               alpha = 0.15, fill = "grey60") +
      scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#E69F00"),
                         guide = "none") +
      theme_bw(base_size = 10) +
      labs(title = paste0("(a) Local PCA / MDS1 — ", chr),
           x = "Position (Mb)", y = "MDS1")
  }
}

# ── Panel B: Regional PCA ─────────────────────────────────────────────────
pca_file <- file.path(step12_dir, paste0("STEP12_", chr, ".candidate_", cid_val,
                                          ".regional_pca_samples.tsv.gz"))
p_pca <- ggplot() + theme_void() + labs(title = "(b) Regional PCA")

if (file.exists(pca_file)) {
  pcs <- fread(pca_file)
  if (all(c("PC1", "PC2", "group_label", "regional_het") %in% names(pcs))) {
    group_colors <- c("Homo_1" = "#0072B2", "Het" = "#D55E00", "Homo_2" = "#009E73")
    avail_groups <- unique(pcs$group_label)
    if (!all(avail_groups %in% names(group_colors))) {
      group_colors <- setNames(c("#0072B2", "#D55E00", "#009E73")[seq_along(avail_groups)],
                               avail_groups)
    }

    p_pca <- ggplot(pcs, aes(x = PC1, y = PC2, color = group_label)) +
      geom_point(size = 1.8, alpha = 0.75) +
      scale_color_manual(values = group_colors) +
      theme_bw(base_size = 10) +
      labs(title = paste0("(b) Regional PCA — candidate ", cid_val),
           x = "PC1", y = "PC2", color = "Group") +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.4, "cm"))
  }
}

# ── Panel C: LD triangle heatmap ──────────────────────────────────────────
p_ld <- ggplot() + theme_void() +
  labs(title = "(c) LD heatmap — data not available")

if (!is.na(ld_pairs_file) && file.exists(ld_pairs_file)) {
  message("[INFO] Reading LD pairs for triangle heatmap...")

  # Compute zoom region
  inv_len <- inv_end - inv_start
  pad <- max(inv_len * 0.3, 500000)
  zoom_start <- max(1, inv_start - pad)
  zoom_end   <- inv_end + pad

  # Read and filter to zoom region
  ld <- fread(ld_pairs_file, header = FALSE, select = c(1, 2, 4),
              col.names = c("site1", "site2", "r2"))
  ld[, pos1 := as.numeric(sub(".*:", "", site1))]
  ld[, pos2 := as.numeric(sub(".*:", "", site2))]
  ld <- ld[is.finite(pos1) & is.finite(pos2) & r2 >= 0 & r2 <= 1]
  ld <- ld[pos1 >= zoom_start & pos1 <= zoom_end &
           pos2 >= zoom_start & pos2 <= zoom_end]

  if (nrow(ld) > 0) {
    # Bin
    ld[, bin1 := floor(pos1 / bin_bp) * bin_bp]
    ld[, bin2 := floor(pos2 / bin_bp) * bin_bp]
    ld[, bx := pmin(bin1, bin2)]
    ld[, by := pmax(bin1, bin2)]

    agg <- ld[, .(r2 = mean(r2, na.rm = TRUE)), by = .(bx, by)]

    apex_cols <- c("#050505", "#121212", "#241238", "#4A1D6B",
                   "#7E2E98", "#C44788", "#F1685B")

    p_ld <- ggplot(agg, aes(x = bx / 1e6, y = by / 1e6, fill = r2)) +
      geom_tile(width = bin_bp / 1e6, height = bin_bp / 1e6) +
      scale_fill_gradientn(colours = apex_cols,
                           limits = c(0, 1),
                           name = expression(r^2)) +
      annotate("rect",
               xmin = inv_start / 1e6, xmax = inv_end / 1e6,
               ymin = inv_start / 1e6, ymax = inv_end / 1e6,
               fill = NA, color = "#8B1E3F", linewidth = 0.8, linetype = "dashed") +
      coord_fixed(expand = FALSE) +
      theme_bw(base_size = 10) +
      theme(panel.grid = element_blank()) +
      labs(title = paste0("(c) LD triangle — candidate ", cid_val),
           x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"))
  }

  rm(ld)
  gc()
}

# ── Panel D: Het boxplot ─────────────────────────────────────────────────
p_het <- ggplot() + theme_void() + labs(title = "(d) Het by group")

if (exists("pcs") && "regional_het" %in% names(pcs) && "group_label" %in% names(pcs)) {
  p_het <- ggplot(pcs, aes(x = group_label, y = regional_het, fill = group_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.12, alpha = 0.4, size = 0.8) +
    scale_fill_manual(values = group_colors) +
    theme_bw(base_size = 10) +
    labs(title = "(d) Regional heterozygosity",
         x = "Inferred group", y = "Dosage het") +
    theme(legend.position = "none")
}

# ── Panel E: Chromosome tP curves ────────────────────────────────────────
plotc_file <- file.path(step12_dir, paste0("STEP12_", chr, ".candidate_", cid_val,
                                            ".plotC_summary.tsv.gz"))
p_theta <- ggplot() + theme_void() + labs(title = "(e) tP by group")

if (file.exists(plotc_file)) {
  plotc_dt <- fread(plotc_file)
  if (all(c("group_label", "WinCenter", "tP_mean") %in% names(plotc_dt))) {
    p_theta <- ggplot(plotc_dt, aes(x = WinCenter / 1e6, y = tP_mean, color = group_label)) +
      annotate("rect",
               xmin = inv_start / 1e6, xmax = inv_end / 1e6,
               ymin = -Inf, ymax = Inf,
               alpha = 0.15, fill = "grey60") +
      geom_line(linewidth = 0.7) +
      scale_color_manual(values = group_colors) +
      theme_bw(base_size = 10) +
      labs(title = paste0("(e) ", chr, " — pairwise θ by group"),
           x = "Position (Mb)", y = "Mean tP",
           color = "Group") +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.4, "cm"))
  }
}

# ── Compose ──────────────────────────────────────────────────────────────
layout <- "
AABB
CCDD
EEEE
"

combined <- p_mds + p_pca + p_ld + p_het + p_theta +
  plot_layout(design = layout, heights = c(1, 1.2, 0.7))

pdf_out <- file.path(outdir, paste0("composite_candidate_", cid_val, "_", chr, ".pdf"))
png_out <- file.path(outdir, paste0("composite_candidate_", cid_val, "_", chr, ".png"))

ggsave(pdf_out, combined, width = 12, height = 14)
ggsave(png_out, combined, width = 12, height = 14, dpi = 400)

message("[DONE] Composite figure for candidate ", cid_val, ":")
message("  ", pdf_out)
message("  ", png_out)
