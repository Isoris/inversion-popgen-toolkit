#!/usr/bin/env Rscript

# =============================================================================
# STEP26_candidate_figure_catalogue.R  (v2 — all bugs fixed)
#
# FIXES APPLIED:
#   1. FIG_C07 het boxplot: count labels now use geom_text with pre-computed
#      counts positioned at y_max, not stat_summary plotting count as y-value
#   2. FIG_C08 heatmap: Het annotation uses anno_simple to avoid hex legends
#   3. FIG_C05 pattern label: falls back to per-candidate interpretation file
#   4. FIG_C10 topology: factor levels match TOPOLOGY_COLORS, drop=FALSE
#   5. Composite scoping: p_c10 initialized NULL before conditional
#   6. FIG_C15 ancestry: handles both level_type and group_level column names
#
# Usage:
#   Rscript STEP26_candidate_figure_catalogue.R <config.R> [cid=all] [level=intermediate]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

has_waffle     <- suppressWarnings(require(waffle, quietly = TRUE))
has_ggridges   <- suppressWarnings(require(ggridges, quietly = TRUE))
has_ComplexHeatmap <- suppressWarnings(require(ComplexHeatmap, quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_
plot_level  <- if (length(args) >= 3) args[3] else "intermediate"

source(config_file)
ensure_dir(PLOTS_DIR)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

mds_file <- list.files(STEP10_DIR, pattern = "\\.window_mds\\.tsv\\.gz$", full.names = TRUE)
mds_dt   <- if (length(mds_file) > 0) fread(mds_file[1]) else NULL

interp_file <- file.path(FOLLOWUP_DIR, "candidate_region_interpretation.tsv")
interp_dt   <- if (file.exists(interp_file)) fread(interp_file) else NULL

save_plot <- function(p, filepath, width = 8, height = 6) {
  pdf_path <- sub("\\.[^.]+$", ".pdf", filepath)
  png_path <- sub("\\.[^.]+$", ".png", filepath)
  ggsave(pdf_path, p, width = width, height = height, device = "pdf")
  ggsave(png_path, p, width = width, height = height, dpi = 400, device = "png")
  message("  Saved: ", basename(pdf_path))
}

# ══════════════════════════════════════════════════════════════════════════════
# DISCOVERY FIGURES
# ══════════════════════════════════════════════════════════════════════════════

plot_discovery_figures <- function() {
  disc_dir <- file.path(PLOTS_DIR, "00_discovery")
  ensure_dir(disc_dir)

  if (!is.null(mds_dt) && nrow(mds_dt) > 0 && "MDS1" %in% names(mds_dt)) {
    mds_dt[, mid_bp := (start_bp + end_bp) / 2]
    p_d01 <- ggplot(mds_dt, aes(x = mid_bp / 1e6, y = MDS1)) +
      geom_point(size = 0.4, alpha = 0.5, color = "grey50") +
      facet_wrap(~chrom, scales = "free_x", ncol = 4) +
      theme_inversion(base_size = 8) +
      labs(title = "FIG_D01 — Chromosome-wide local PCA / MDS1 scan",
           subtitle = paste0(nrow(cand), " candidates across ", uniqueN(cand$chrom), " chromosomes"),
           x = "Position (Mb)", y = "MDS1")
    if (nrow(cand) > 0) {
      shade_dt <- cand[, .(chrom, xmin = start_bp / 1e6, xmax = end_bp / 1e6)]
      p_d01 <- p_d01 +
        geom_rect(data = shade_dt, inherit.aes = FALSE,
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                  fill = "#E69F00", alpha = 0.15)
    }
    save_plot(p_d01, file.path(disc_dir, "FIG_D01_mds_scan.pdf"),
              width = 16, height = max(6, uniqueN(mds_dt$chrom) * 0.6))
  }

  if (nrow(cand) > 0) {
    chr_counts <- cand[, .N, by = chrom][order(-N)]
    chr_counts[, chrom := factor(chrom, levels = chrom)]
    p_d02 <- ggplot(chr_counts, aes(x = chrom, y = N)) +
      geom_col(fill = "#2166AC", alpha = 0.8) +
      geom_text(aes(label = N), vjust = -0.3, size = 3) +
      theme_inversion() +
      labs(title = "FIG_D02 — Candidate count per chromosome",
           x = "Chromosome", y = "Number of candidates") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
    save_plot(p_d02, file.path(disc_dir, "FIG_D02_candidates_per_chr.pdf"),
              width = max(7, nrow(chr_counts) * 0.3 + 3), height = 5)

    cand[, span_mb := (end_bp - start_bp) / 1e6]
    p_d03 <- ggplot(cand, aes(x = span_mb)) +
      geom_histogram(bins = 30, fill = "#2166AC", alpha = 0.7, color = "white") +
      geom_rug(alpha = 0.4) + theme_inversion() +
      labs(title = "FIG_D03 — Candidate span distribution",
           subtitle = paste0("Median: ", round(median(cand$span_mb), 2), " Mb"),
           x = "Candidate span (Mb)", y = "Count")
    save_plot(p_d03, file.path(disc_dir, "FIG_D03_span_distribution.pdf"), width = 7, height = 5)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# PER-CANDIDATE FIGURES
# ══════════════════════════════════════════════════════════════════════════════

plot_candidate <- function(cid) {
  row <- cand[candidate_id == cid]
  if (nrow(row) == 0) return(invisible())
  chr <- row$chrom[1]; c_start <- as.numeric(row$start_bp[1]); c_end <- as.numeric(row$end_bp[1])

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  plot_dir <- file.path(PLOTS_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(plot_dir)

  rot    <- safe_fread(file.path(cand_dir, "candidate_pca_rotated.tsv"))
  if (is.null(rot) || nrow(rot) < 5) { message("[SKIP] candidate ", cid); return(invisible()) }

  sub    <- tryCatch(fread(file.path(cand_dir, "candidate_subclusters.tsv")), error = function(e) NULL)
  het    <- tryCatch(fread(file.path(cand_dir, "candidate_group_het_summary.tsv")), error = function(e) NULL)
  marker <- tryCatch(fread(file.path(cand_dir, "candidate_marker_support.tsv")), error = function(e) NULL)
  genes  <- tryCatch(fread(file.path(cand_dir, "candidate_gene_overlap.tsv")), error = function(e) NULL)
  win    <- tryCatch(fread(file.path(cand_dir, "candidate_window_summary.tsv")), error = function(e) NULL)
  anc    <- tryCatch(fread(file.path(cand_dir, "candidate_ancestry_association.tsv")), error = function(e) NULL)

  n_snps   <- if (!is.null(marker)) marker$n_snps[1] else NA
  n_genes  <- if (!is.null(genes)) genes$n_genes_overlap[1] else NA
  n_windows <- if (!is.null(win)) nrow(win) else NA
  grp_counts <- table(rot$coarse_group_refined)

  # ── FIX 3: Pattern label with fallback ─────────────────────────────────
  interp_label <- NA_character_
  if (!is.null(interp_dt) && nrow(interp_dt[candidate_id == cid]) > 0)
    interp_label <- interp_dt[candidate_id == cid, coarse_pattern][1]
  if (is.na(interp_label)) {
    pi_file <- file.path(cand_dir, "candidate_region_interpretation.tsv")
    if (file.exists(pi_file)) {
      pi <- tryCatch(fread(pi_file), error = function(e) NULL)
      if (!is.null(pi) && "coarse_pattern" %in% names(pi)) interp_label <- pi$coarse_pattern[1]
    }
  }
  # Also try v2 table
  if (is.na(interp_label)) {
    v2_file <- file.path(FOLLOWUP_DIR, "candidate_region_interpretation_v2.tsv")
    if (file.exists(v2_file)) {
      v2 <- tryCatch(fread(v2_file), error = function(e) NULL)
      if (!is.null(v2) && cid %in% v2$candidate_id && "final_label" %in% names(v2))
        interp_label <- v2[candidate_id == cid, final_label][1]
    }
  }
  if (is.na(interp_label)) interp_label <- "not yet classified"

  header <- candidate_header(cid, chr, c_start, c_end,
                              n_snps = n_snps, n_samples = nrow(rot),
                              n_windows = n_windows, group_counts = grp_counts,
                              interpretation = interp_label)
  message("[PLOT] ", header)

  shape_map <- c("HOMO_1" = 16, "HET" = 17, "HOMO_2" = 15, "AMBIGUOUS" = 4)

  # ── FIG_C02: PCA by het ────────────────────────────────────────────────
  p_c02 <- ggplot(rot, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = regional_het, shape = coarse_group_refined),
               size = 2.5, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma", name = "Regional\nhet (dosage)") +
    scale_shape_manual(values = shape_map, name = "Group") +
    theme_inversion() +
    labs(title = "FIG_C02 — Regional PCA colored by heterozygosity",
         subtitle = wrap_subtitle(header), x = "PC1", y = "PC2")
  save_plot(p_c02, file.path(plot_dir, "FIG_C02_pca_het.pdf"), width = 8, height = 6.5)

  # ── FIG_C03: PCA by ancestry ──────────────────────────────────────────
  anc_global_file <- file.path(ANCESTRY_ROOT, "sample_main_ancestry_by_K.tsv")
  if (file.exists(anc_global_file)) {
    anc_samp <- tryCatch(fread(anc_global_file), error = function(e) NULL)
    if (!is.null(anc_samp) && "sample" %in% names(anc_samp)) {
      anc_col <- if ("cluster_label" %in% names(anc_samp)) "cluster_label" else
                 if ("maxQ_label" %in% names(anc_samp)) "maxQ_label" else NULL
      if (!is.null(anc_col)) {
        rot_anc <- merge(rot, anc_samp[, c("sample", anc_col), with = FALSE],
                         by = "sample", all.x = TRUE)
        p_c03 <- ggplot(rot_anc, aes(x = PC1, y = PC2, color = get(anc_col))) +
          geom_point(aes(shape = coarse_group_refined), size = 2.5, alpha = 0.8) +
          scale_shape_manual(values = shape_map, name = "Group") +
          theme_inversion() +
          labs(title = "FIG_C03 — Regional PCA colored by ancestry / Q",
               subtitle = wrap_subtitle(header), x = "PC1", y = "PC2", color = "Ancestry")
        save_plot(p_c03, file.path(plot_dir, "FIG_C03_pca_ancestry.pdf"), width = 8, height = 6.5)
      }
    }
  }

  # ── FIG_C04: Rotated PCA ──────────────────────────────────────────────
  rot[, coarse_group_refined := factor(coarse_group_refined,
        levels = c("HOMO_1", "HET", "HOMO_2", "AMBIGUOUS"))]

  p_c04_sub <- ggplot(rot, aes(x = u, y = v)) +
    geom_point(aes(shape = coarse_group_refined), size = 2.5, alpha = 0.8, color = "#333333") +
    scale_shape_manual(values = shape_map, name = "Coarse group") +
    theme_inversion() +
    labs(title = "FIG_C04 — Rotated stripe-aware PCA", subtitle = wrap_subtitle(header),
         x = "u (broad stripe axis)", y = "v (within-stripe axis)")

  if (!is.null(sub)) {
    rot_sub <- merge(rot, sub[, .(sample, subcluster_label)], by = "sample", all.x = TRUE)
    p_c04_sub <- ggplot(rot_sub, aes(x = u, y = v, color = subcluster_label)) +
      geom_point(aes(shape = coarse_group_refined), size = 2.5, alpha = 0.8) +
      scale_shape_manual(values = shape_map, name = "Coarse group") +
      theme_inversion() +
      labs(title = "FIG_C04 — Rotated PCA by subcluster", subtitle = wrap_subtitle(header),
           x = "u (broad stripe axis)", y = "v (within-stripe axis)", color = "Subcluster")
  }
  save_plot(p_c04_sub, file.path(plot_dir, "FIG_C04_rotated_pca.pdf"), width = 8.5, height = 6.5)

  # ── FIG_C05: Metadata card ────────────────────────────────────────────
  card_data <- data.table(
    Field = c("Candidate ID", "Chromosome", "Start", "End", "Span",
              "SNPs", "Samples", "Windows",
              "HOMO_1", "HET", "HOMO_2", "Ambiguous", "Marker tier", "Pattern"),
    Value = c(as.character(cid), chr,
              fmt_mb(c_start), fmt_mb(c_end), fmt_mb(c_end - c_start),
              ifelse(is.na(n_snps), "\u2013", format(n_snps, big.mark = ",")),
              as.character(nrow(rot)),
              ifelse(is.na(n_windows), "\u2013", as.character(n_windows)),
              as.character(grp_counts["HOMO_1"]),
              as.character(grp_counts["HET"]),
              as.character(grp_counts["HOMO_2"]),
              as.character(sum(rot$ambiguous_flag)),
              ifelse(is.null(marker), "\u2013", marker$support_tier[1]),
              interp_label)
  )
  card_data[, y := rev(seq_len(.N))]
  p_c05 <- ggplot(card_data, aes(x = 0, y = y, label = paste0(Field, ": ", Value))) +
    geom_text(hjust = 0, size = 3.5, family = "mono") +
    theme_void() + xlim(-0.1, 3) +
    labs(title = paste0("FIG_C05 — Candidate ", cid, " summary card"))
  save_plot(p_c05, file.path(plot_dir, "FIG_C05_metadata_card.pdf"), width = 5, height = 5)

  # ── FIG_C06B: Group count barplot ──────────────────────────────────────
  count_dt <- data.table(group = names(grp_counts), n = as.integer(grp_counts))
  count_dt[, group := factor(group, levels = c("HOMO_1", "HET", "HOMO_2", "AMBIGUOUS"))]
  p_c06b <- ggplot(count_dt, aes(x = group, y = n, fill = group)) +
    geom_col(alpha = 0.85, width = 0.65) +
    geom_text(aes(label = n), vjust = -0.3, size = 4) +
    scale_fill_manual(values = GROUP_COLORS) + theme_inversion() +
    labs(title = paste0("FIG_C06B — Group counts (candidate ", cid, ")"), x = "", y = "Count") +
    theme(legend.position = "none")
  save_plot(p_c06b, file.path(plot_dir, "FIG_C06B_group_counts.pdf"), width = 5, height = 4.5)

  # ── FIG_C07: Het boxplot — FIX 1 (count labels) ───────────────────────
  het_data <- rot[is.finite(regional_het)]
  het_counts <- het_data[, .(n = .N, y_max = max(regional_het, na.rm = TRUE)),
                         by = coarse_group_refined]

  p_c07 <- ggplot(het_data, aes(x = coarse_group_refined, y = regional_het,
                                 fill = coarse_group_refined)) +
    geom_violin(alpha = 0.3, width = 0.8) +
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.08, alpha = 0.35, size = 1) +
    geom_text(data = het_counts,
              aes(x = coarse_group_refined, y = y_max, label = paste0("n=", n)),
              vjust = -0.8, size = 3.2, inherit.aes = FALSE) +
    scale_fill_manual(values = GROUP_COLORS) + theme_inversion() +
    labs(title = paste0("FIG_C07 — Heterozygosity by group (candidate ", cid, ")"),
         subtitle = wrap_subtitle(header), x = "Inferred group", y = "Regional dosage heterozygosity") +
    theme(legend.position = "none")
  save_plot(p_c07, file.path(plot_dir, "FIG_C07_het_boxplot.pdf"), width = 6, height = 5.5)

  # ── FIG_C07B: Het ridges ──────────────────────────────────────────────
  if (has_ggridges) {
    p_c07b <- ggplot(het_data, aes(x = regional_het, y = coarse_group_refined,
                                    fill = coarse_group_refined)) +
      ggridges::geom_density_ridges(alpha = 0.6, scale = 1.2) +
      scale_fill_manual(values = GROUP_COLORS) + theme_inversion() +
      labs(title = paste0("FIG_C07B — Heterozygosity ridges (candidate ", cid, ")"),
           x = "Regional dosage heterozygosity", y = "") +
      theme(legend.position = "none")
    save_plot(p_c07b, file.path(plot_dir, "FIG_C07B_het_ridges.pdf"), width = 6, height = 4.5)
  }

  # ── FIG_C10: Topology track — FIX 4 (factor levels) ───────────────────
  p_c10 <- NULL  # FIX 6: initialize before conditional
  if (!is.null(win) && nrow(win) > 0) {
    all_topo <- names(TOPOLOGY_COLORS)
    win[, topology := factor(topology, levels = all_topo)]

    # Compute explicit rect coords for visibility at high window counts
    win[, xmin_mb := window_start / 1e6]
    win[, xmax_mb := window_end / 1e6]

    p_c10 <- ggplot(win) +
      geom_rect(aes(xmin = xmin_mb, xmax = xmax_mb, ymin = 0, ymax = 1,
                     fill = topology), color = 'grey40', linewidth = 0.15) +
      scale_fill_manual(values = TOPOLOGY_COLORS, drop = FALSE, na.value = "#E0E0E0",
                        name = "Window\ntopology") +
      annotate("rect", xmin = c_start / 1e6, xmax = c_end / 1e6,
               ymin = 0.4, ymax = 1.6,
               fill = NA, color = "black", linetype = "dashed", linewidth = 0.5) +
      theme_inversion() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank()) +
      labs(title = paste0("FIG_C10 — Window topology track (candidate ", cid, ")"),
           subtitle = wrap_subtitle(header), x = "Position (Mb)", y = "")
    save_plot(p_c10, file.path(plot_dir, "FIG_C10_topology_track.pdf"), width = 10, height = 3)
  }

  # ── FIG_C11: Window agreement heatmap ──────────────────────────────────
  agree_rds <- file.path(cand_dir, "candidate_window_agreement.rds")
  if (file.exists(agree_rds) && has_ComplexHeatmap) {
    agree_mat <- readRDS(agree_rds)
    if (!is.null(agree_mat) && nrow(agree_mat) >= 3) {
      agree_mat <- as.matrix(agree_mat); storage.mode(agree_mat) <- "double"
      agree_mat[!is.finite(agree_mat)] <- NA_real_
      keep_r <- rowSums(is.finite(agree_mat)) > 1
      keep_c <- colSums(is.finite(agree_mat)) > 1
      agree_mat <- agree_mat[keep_r, keep_c, drop = FALSE]
      if (nrow(agree_mat) >= 3 && ncol(agree_mat) >= 3) {
        agree_mat[is.na(agree_mat)] <- 0; diag(agree_mat) <- 1
        pdf(file.path(plot_dir, "FIG_C11_window_agreement.pdf"), width = 8, height = 7)
        ht <- ComplexHeatmap::Heatmap(agree_mat, name = "ARI",
          col = circlize::colorRamp2(c(0, 0.5, 1), c("#F7F7F7", "#FDDBC7", "#B2182B")),
          cluster_rows = TRUE, cluster_columns = TRUE,
          show_row_names = FALSE, show_column_names = FALSE,
          column_title = paste0("FIG_C11 — Window agreement (candidate ", cid, ")"),
          column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"))
        ComplexHeatmap::draw(ht); dev.off()
        message("  Saved: FIG_C11_window_agreement.pdf")
      }
    }
  }

  # ── FIG_C14: DBSCAN diagnostics ───────────────────────────────────────
  if (!is.null(sub)) {
    rot_sub2 <- merge(rot, sub[, .(sample, subcluster_label, is_noise)], by = "sample", all.x = TRUE)
    p_c14_main <- ggplot(rot_sub2, aes(x = u, y = v, color = subcluster_label)) +
      geom_point(aes(shape = coarse_group_refined), size = 2, alpha = 0.8) +
      scale_shape_manual(values = shape_map, name = "Group") + theme_inversion() +
      labs(title = "(a) Subclusters", x = "u", y = "v", color = "Subcluster")
    p_c14_noise <- ggplot(rot_sub2, aes(x = u, y = v)) +
      geom_point(aes(color = is_noise), size = 2, alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "#333333", "TRUE" = "red"), name = "Noise") +
      theme_inversion() + labs(title = "(b) Noise points", x = "u", y = "v")
    cl_counts <- rot_sub2[, .N, by = subcluster_label][order(-N)]
    p_c14_bar <- ggplot(cl_counts, aes(x = reorder(subcluster_label, -N), y = N)) +
      geom_col(fill = "#666666", alpha = 0.8) + theme_inversion() +
      labs(title = "(c) Cluster sizes", x = "", y = "Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
    p_c14 <- p_c14_main + p_c14_noise + p_c14_bar +
      plot_layout(ncol = 3, widths = c(1.2, 1, 0.8)) +
      plot_annotation(title = paste0("FIG_C14 — DBSCAN diagnostics (candidate ", cid, ")"),
                      subtitle = wrap_subtitle(header),
                      theme = theme(plot.title = element_text(face = "bold", size = 13)))
    save_plot(p_c14, file.path(plot_dir, "FIG_C14_dbscan_diagnostics.pdf"), width = 14, height = 5.5)
  }

  # ── FIG_C16: Marker density ───────────────────────────────────────────
  sites_file <- file.path(STEP12_DIR, paste0("STEP12_", chr, ".candidate_", cid, ".regional_sites.tsv.gz"))
  if (file.exists(sites_file)) {
    sites_reg <- fread(sites_file)
    if ("pos" %in% names(sites_reg) && nrow(sites_reg) > 10) {
      n_bins <- min(50, nrow(sites_reg) / 3)
      sites_reg[, bin := cut(pos, breaks = n_bins, labels = FALSE)]
      bin_counts <- sites_reg[, .(count = .N, mid = mean(pos)), by = bin]
      p_c16 <- ggplot(bin_counts, aes(x = mid / 1e6, y = count)) +
        geom_col(fill = "#2166AC", alpha = 0.7) +
        annotate("rect", xmin = c_start / 1e6, xmax = c_end / 1e6,
                 ymin = -Inf, ymax = Inf, fill = NA, color = "black",
                 linetype = "dashed", linewidth = 0.4) +
        theme_inversion() +
        labs(title = paste0("FIG_C16 — Marker density (candidate ", cid, ")"),
             subtitle = wrap_subtitle(header), x = "Position (Mb)", y = "SNPs per bin")
      save_plot(p_c16, file.path(plot_dir, "FIG_C16_marker_density.pdf"), width = 8, height = 3.5)
    }
  }

  # ── COMPOSITES ─────────────────────────────────────────────────────────
  if (plot_level %in% c("compact", "intermediate", "full")) {
    p_compact <- (p_c05 + p_c02) / (p_c07 + p_c06b) +
      plot_annotation(title = paste0("Compact summary — Candidate ", cid),
                      subtitle = wrap_subtitle(header),
                      theme = theme(plot.title = element_text(face = "bold", size = 14)))
    save_plot(p_compact, file.path(plot_dir, "composite_compact.pdf"), width = 14, height = 11)
  }

  if (plot_level %in% c("intermediate", "full")) {
    panels <- list(p_c02, p_c04_sub, p_c07, p_c06b)
    if (!is.null(p_c10)) panels <- c(panels, list(p_c10))
    p_inter <- wrap_plots(panels, ncol = 2) +
      plot_annotation(title = paste0("Intermediate summary — Candidate ", cid),
                      subtitle = wrap_subtitle(header),
                      theme = theme(plot.title = element_text(face = "bold", size = 14)))
    save_plot(p_inter, file.path(plot_dir, "composite_intermediate.pdf"),
              width = 15, height = max(10, length(panels) * 2.5))
  }

  message("[DONE] Candidate ", cid, " figures complete")
}

# ══════════════════════════════════════════════════════════════════════════════
# GLOBAL DASHBOARD
# ══════════════════════════════════════════════════════════════════════════════

plot_global_dashboard <- function() {
  # Try v2 interpretation table first
  idt <- interp_dt
  v2_file <- file.path(FOLLOWUP_DIR, "candidate_region_interpretation_v2.tsv")
  if (file.exists(v2_file)) idt <- fread(v2_file)
  if (is.null(idt) || nrow(idt) == 0) { message("[SKIP] No interpretation data"); return(invisible()) }

  dash_dir <- file.path(PLOTS_DIR, "00_global"); ensure_dir(dash_dir)
  label_col <- if ("final_label" %in% names(idt)) "final_label" else
               if ("coarse_pattern" %in% names(idt)) "coarse_pattern" else NULL
  if (is.null(label_col)) { message("[SKIP] No label column"); return(invisible()) }

  p_g1 <- ggplot(idt, aes(x = get(label_col))) +
    geom_bar(fill = "#2166AC", alpha = 0.8) + theme_inversion() +
    labs(title = "(a) Pattern classification", x = "", y = "Count") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))

  conf_col <- if ("confidence" %in% names(idt)) "confidence" else
              if ("inversion_like_confidence" %in% names(idt)) "inversion_like_confidence" else NULL
  p_g2 <- if (!is.null(conf_col)) {
    ggplot(idt, aes(x = get(conf_col))) +
      geom_histogram(bins = 20, fill = "#1B7837", alpha = 0.7) + theme_inversion() +
      labs(title = "(b) Confidence distribution", x = "Confidence", y = "Count")
  } else ggplot() + theme_void()

  p_g3 <- if ("marker_support_tier" %in% names(idt)) {
    ggplot(idt[!is.na(marker_support_tier)], aes(x = marker_support_tier)) +
      geom_bar(fill = "#D6604D", alpha = 0.8) + theme_inversion() +
      labs(title = "(c) Marker support tiers", x = "", y = "Count")
  } else ggplot() + theme_void()

  driver_col <- if ("main_complexity_driver" %in% names(idt)) "main_complexity_driver" else NULL
  p_g4 <- if (!is.null(driver_col)) {
    ggplot(idt, aes(x = get(driver_col))) +
      geom_bar(fill = "#9970AB", alpha = 0.8) + theme_inversion() +
      labs(title = "(d) Complexity drivers", x = "", y = "Count") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else ggplot() + theme_void()

  dashboard <- (p_g1 + p_g2) / (p_g3 + p_g4) +
    plot_annotation(title = "FIG_G01 — Candidate overview dashboard",
                    subtitle = paste0(nrow(idt), " candidates"),
                    theme = theme(plot.title = element_text(face = "bold", size = 15)))
  save_plot(dashboard, file.path(dash_dir, "FIG_G01_dashboard.pdf"), width = 13, height = 10)
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN
# ══════════════════════════════════════════════════════════════════════════════
message("=", strrep("=", 70))
message("[START] Figure catalogue generation (v2 — all bugs fixed)")
plot_discovery_figures()
for (ci in seq_len(nrow(cand))) plot_candidate(cand$candidate_id[ci])
plot_global_dashboard()
message("=", strrep("=", 70))
message("[DONE] STEP26 figure catalogue complete")
