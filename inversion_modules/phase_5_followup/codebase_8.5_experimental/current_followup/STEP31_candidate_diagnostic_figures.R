#!/usr/bin/env Rscript
# STEP31_candidate_diagnostic_figures.R v2 — all 15 fixes applied
# Fail-soft: each panel guarded independently, bad panels skip not crash.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
has_ComplexHeatmap <- suppressPackageStartupMessages(
  suppressWarnings(require(ComplexHeatmap, quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_
source(config_file)

tier_colors <- c("core"="#1B7837","peripheral"="#F4A582","junk"="#D73027")
coherence_colors <- c("coherent"="#1B7837","intermediate"="#FDB863",
                       "discordant"="#D73027","insufficient"="#BDBDBD")

save_plot <- function(p, fp, width=8, height=6) {
  tryCatch({
    ggsave(sub("\\.[^.]+$",".pdf",fp), p, width=width, height=height, device="pdf")
    ggsave(sub("\\.[^.]+$",".png",fp), p, width=width, height=height, dpi=400, device="png")
    message("  Saved: ", basename(fp))
  }, error = function(e) message("  [WARN] Save failed: ", conditionMessage(e)))
}

make_state_numeric_matrix <- function(state_dt) {
  if (is.null(state_dt) || nrow(state_dt) < 1 || !"sample" %in% names(state_dt)) return(NULL)
  mat_cols <- setdiff(names(state_dt), "sample")
  if (length(mat_cols) < 1) return(NULL)
  samples <- as.character(state_dt$sample)
  char_mat <- as.matrix(state_dt[, ..mat_cols]); storage.mode(char_mat) <- "character"
  num_mat <- matrix(NA_real_, nrow(char_mat), ncol(char_mat))
  rownames(num_mat) <- samples; colnames(num_mat) <- mat_cols
  char_mat[char_mat %in% c("","NA","NaN","NULL")] <- NA_character_
  num_mat[char_mat == "AA"] <- 0; num_mat[char_mat == "AB"] <- 1; num_mat[char_mat == "BB"] <- 2
  keep_r <- rowSums(!is.na(num_mat)) > 0; keep_c <- colSums(!is.na(num_mat)) > 0
  if (!any(keep_r) || !any(keep_c)) return(NULL)
  num_mat[keep_r, keep_c, drop = FALSE]
}

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
  traj <- safe_fread(file.path(cand_dir, "candidate_trajectory_summary.tsv"))
  state_mat <- safe_fread(file.path(cand_dir, "candidate_window_state_matrix.tsv.gz"))
  core_pca <- safe_fread(file.path(cand_dir, "candidate_core_pca.tsv"))
  blocks <- safe_fread(file.path(cand_dir, "candidate_marker_blocks.tsv"))
  if (is.null(rot) || nrow(rot) < 5) next

  header <- wrap_subtitle(candidate_header(cid, chr, c_start, c_end, n_samples = nrow(rot)))
  message("[PLOT] Diagnostics for candidate ", cid)

  # ── FIG_C30: Coherence scatter ─────────────────────────────────────────
  if (!is.null(coh) && nrow(coh) > 0) {
    merged <- merge(rot, coh[, .(sample, agreement_fraction, coherence_class, stripe_quality)],
                    by = "sample", all.x = TRUE)
    merged[is.na(coherence_class), coherence_class := "insufficient"]
    merged[is.na(stripe_quality), stripe_quality := "junk"]

    p <- ggplot(merged[is.finite(agreement_fraction)],
                aes(x = u, y = agreement_fraction, color = coherence_class)) +
      geom_point(aes(shape = coarse_group_refined), size = 2.5, alpha = 0.8) +
      geom_hline(yintercept = c(0.40, 0.55, 0.70), linetype = "dashed", color = "grey50", linewidth = 0.3) +
      scale_color_manual(values = coherence_colors, name = "Coherence") +
      scale_shape_manual(values = safe_shape_values(merged$coarse_group_refined), name = "Group") +
      theme_inversion() +
      labs(title = paste0("FIG_C30 — Sample coherence (candidate ", cid, ")"),
           subtitle = header, x = "u (broad stripe axis)", y = "Agreement fraction")
    save_plot(p, file.path(plot_dir, "FIG_C30_coherence_scatter.pdf"), 8.5, 6)

    # ── FIG_C31 ──────────────────────────────────────────────────────────
    p <- ggplot(merged, aes(x = PC1, y = PC2, color = stripe_quality)) +
      geom_point(aes(shape = coarse_group_refined), size = 2.5, alpha = 0.8) +
      scale_color_manual(values = c(tier_colors, "unknown" = "#BDBDBD"), name = "Quality tier") +
      scale_shape_manual(values = safe_shape_values(merged$coarse_group_refined), name = "Group") +
      theme_inversion() +
      labs(title = paste0("FIG_C31 — Stripe quality tiers (candidate ", cid, ")"),
           subtitle = header, x = "PC1", y = "PC2")
    save_plot(p, file.path(plot_dir, "FIG_C31_stripe_quality_pca.pdf"), 8.5, 6)

    # ── FIG_C35 ──────────────────────────────────────────────────────────
    p <- ggplot(merged[is.finite(agreement_fraction)],
                aes(x = coarse_group_refined, y = agreement_fraction, fill = coarse_group_refined)) +
      geom_violin(alpha = 0.3, width = 0.8) +
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
      geom_jitter(aes(color = stripe_quality), width = 0.08, size = 1.2, alpha = 0.7) +
      scale_fill_manual(values = safe_group_colors(merged$coarse_group_refined)) +
      scale_color_manual(values = c(tier_colors, "unknown" = "#BDBDBD"), name = "Quality") +
      theme_inversion() +
      labs(title = paste0("FIG_C35 — Coherence by stripe (candidate ", cid, ")"),
           x = "Group", y = "Agreement fraction") + guides(fill = "none")
    save_plot(p, file.path(plot_dir, "FIG_C35_coherence_by_stripe.pdf"), 7, 5.5)
  }

  # ── FIG_C32: Marker polarity track ─────────────────────────────────────
  if (!is.null(pol) && nrow(pol) > 0 &&
      all(c("pos", "delta_hom_mean") %in% names(pol))) {

    flip_col <- if ("final_flip_decision" %in% names(pol)) "final_flip_decision" else
                if ("polarity_reversed" %in% names(pol)) "polarity_reversed" else NULL

    if (!is.null(flip_col)) {
      bw <- max((c_end - c_start) / 1e6 / nrow(pol) * 0.8, 1e-6)
      p <- ggplot(pol, aes(x = pos / 1e6, y = delta_hom_mean)) +
        geom_col(aes(fill = as.character(get(flip_col))), width = bw) +
        scale_fill_manual(values = c("FALSE"="#2166AC","TRUE"="#D73027"),
                          labels = c("FALSE"="Dominant","TRUE"="Reversed"), name = "Polarity") +
        geom_hline(yintercept = 0, linewidth = 0.3) + theme_inversion() +
        labs(title = paste0("FIG_C32 — Marker polarity track (candidate ", cid, ")"),
             subtitle = wrap_subtitle(paste0(header, " | ",
               sum(pol[[flip_col]] == TRUE | pol[[flip_col]] == "TRUE", na.rm = TRUE),
               "/", nrow(pol), " markers flipped")),
             x = "Position (Mb)",
             y = expression(Delta*"hom (HOMO_2 \u2212 HOMO_1 mean dosage)"))
      save_plot(p, file.path(plot_dir, "FIG_C32_polarity_track.pdf"), 10, 4)
    }
  }

  # ── FIG_C33: Trajectory heatmap (all 15 STEP31 fixes applied) ──────────
  if (!is.null(state_mat) && nrow(state_mat) > 0 && has_ComplexHeatmap) {
    num_mat <- make_state_numeric_matrix(state_mat)
    if (!is.null(num_mat) && nrow(num_mat) >= 3 && ncol(num_mat) >= 3) {
      samples_present <- rownames(num_mat)
      rod <- rot[match(samples_present, sample), .(sample, coarse_group_refined, u)]
      rod <- rod[!is.na(sample) & sample %in% samples_present]
      if (!is.null(coh)) {
        rod <- merge(rod, coh[, .(sample, stripe_quality)], by = "sample", all.x = TRUE)
      } else rod[, stripe_quality := NA_character_]
      rod[is.na(coarse_group_refined), coarse_group_refined := "AMBIGUOUS"]
      rod[is.na(u), u := Inf]
      setorder(rod, coarse_group_refined, u)
      os <- intersect(unique(rod$sample), rownames(num_mat))
      if (length(os) >= 3) {
        num_mat <- num_mat[os, , drop = FALSE]
        ri <- rod[match(os, sample)]; ri <- ri[!is.na(sample)]
        if (nrow(ri) != nrow(num_mat)) {
          common <- intersect(os, ri$sample)
          num_mat <- num_mat[common, , drop = FALSE]; ri <- ri[match(common, sample)]
        }
        if (nrow(num_mat) >= 3 && nrow(ri) >= 3) {
          gc <- safe_group_colors(ri$coarse_group_refined)
          ha <- ComplexHeatmap::rowAnnotation(
            Group = ri$coarse_group_refined, col = list(Group = gc),
            show_annotation_name = TRUE, annotation_name_gp = grid::gpar(fontsize = 7))
          w <- max(8, min(24, ncol(num_mat) * 0.03 + 4))
          h <- max(6, min(28, nrow(num_mat) * 0.045 + 3))
          pdf_f <- file.path(plot_dir, "FIG_C33_trajectory_heatmap.pdf")
          tryCatch({
            pdf(pdf_f, width = w, height = h)
            ht <- ComplexHeatmap::Heatmap(num_mat, name = "State",
              col = c("0"="#2166AC","1"="#F7F7F7","2"="#B2182B"), na_col = "#E0E0E0",
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = FALSE,
              left_annotation = ha,
              column_title = paste0("FIG_C33 — Trajectory | Candidate ", cid),
              column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
              heatmap_legend_param = list(title="State", at=c(0,1,2),
                                          labels=c("AA-like","AB-like","BB-like")),
              use_raster = nrow(num_mat) * ncol(num_mat) > 50000)
            ComplexHeatmap::draw(ht); dev.off()
            message("  Saved: FIG_C33_trajectory_heatmap.pdf")
          }, error = function(e) { try(dev.off(), silent = TRUE)
            message("  [WARN] FIG_C33: ", conditionMessage(e)) })
        }
      }
    }
  }

  # ── FIG_C34: Core vs all PCA ───────────────────────────────────────────
  if (!is.null(core_pca) && nrow(core_pca) >= 5) {
    p_all <- ggplot(rot, aes(x = PC1, y = PC2, color = coarse_group_refined)) +
      geom_point(aes(shape = coarse_group_refined), size = 2, alpha = 0.7) +
      scale_color_manual(values = safe_group_colors(rot$coarse_group_refined)) +
      scale_shape_manual(values = safe_shape_values(rot$coarse_group_refined)) +
      theme_inversion() + labs(title = paste0("(a) All samples (n=",nrow(rot),")"),
                                x = "PC1", y = "PC2") + theme(legend.position = "none")
    p_core <- ggplot(core_pca, aes(x = core_PC1, y = core_PC2, color = coarse_group_refined)) +
      geom_point(aes(shape = coarse_group_refined), size = 2, alpha = 0.7) +
      scale_color_manual(values = safe_group_colors(core_pca$coarse_group_refined)) +
      scale_shape_manual(values = safe_shape_values(core_pca$coarse_group_refined)) +
      theme_inversion() + labs(title = paste0("(b) Core-only (n=",nrow(core_pca),")"),
                                x = "Core PC1", y = "Core PC2") + theme(legend.position = "bottom")
    p <- p_all + p_core + plot_annotation(
      title = paste0("FIG_C34 — Core vs all PCA (candidate ", cid, ")"),
      subtitle = header, theme = theme(plot.title = element_text(face = "bold", size = 13)))
    save_plot(p, file.path(plot_dir, "FIG_C34_core_vs_all_pca.pdf"), 13, 6)
  }

  # ── FIG_C36: Block summary ────────────────────────────────────────────
  if (!is.null(blocks) && nrow(blocks) > 0 &&
      all(c("start_pos","end_pos") %in% names(blocks))) {
    delta_col <- if ("block_delta_hom_mean" %in% names(blocks)) "block_delta_hom_mean" else
                 if ("block_mean_delta" %in% names(blocks)) "block_mean_delta" else NULL
    dir_col <- if ("dominant_direction" %in% names(blocks)) "dominant_direction" else NULL
    if (!is.null(delta_col) && !is.null(dir_col)) {
      bw <- pmax((blocks$end_pos - blocks$start_pos) / 1e6 * 0.9, 1e-6)
      p <- ggplot(blocks, aes(x = (start_pos+end_pos)/2/1e6, y = get(delta_col), fill = get(dir_col))) +
        geom_col(width = bw) +
        scale_fill_manual(values = c("positive"="#2166AC","negative"="#D73027"), name = "Direction") +
        geom_hline(yintercept = 0, linewidth = 0.3) + theme_inversion() +
        labs(title = paste0("FIG_C36 — Marker blocks (candidate ", cid, ")"),
             subtitle = paste0(nrow(blocks), " contiguous blocks"),
             x = "Position (Mb)", y = "Block mean \u0394hom")
      save_plot(p, file.path(plot_dir, "FIG_C36_marker_blocks.pdf"), 9, 4)
    }
  }

  message("[DONE] Candidate ", cid, " diagnostic figures complete")
}
message("[DONE] STEP31 v2 complete")
