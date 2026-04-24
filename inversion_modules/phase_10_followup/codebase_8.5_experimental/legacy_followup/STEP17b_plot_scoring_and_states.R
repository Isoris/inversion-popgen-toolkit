#!/usr/bin/env Rscript

# =============================================================================
# STEP17b_plot_scoring_and_states.R
#
# Visualization suite for the integrated scoring table from STEP17.
# Per candidate:
#   1. Regional PCA colored by provisional major state
#   2. Regional PCA colored by confidence level
#   3. Evidence score distributions
# Global:
#   4. Candidate classification summary barplot
#   5. Heatmap of classification across all candidates
#
# Usage:
#   Rscript STEP17b_plot_scoring_and_states.R \
#     <scoring_table> <classification_table> <step12_dir> <outdir> \
#     [candidate_id=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript STEP17b_... <scoring.tsv.gz> <classification.tsv> <step12_dir> <outdir> [cid=all]")
}

scoring_file <- args[1]
class_file   <- args[2]
step12_dir   <- args[3]
outdir       <- args[4]
cid_filter   <- if (length(args) >= 5 && args[5] != "all") args[5] else NA_character_

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

scoring <- fread(scoring_file)
classification <- fread(class_file)

state_colors <- c(
  "FULL_A"  = "#0072B2",
  "HALF"    = "#D55E00",
  "FULL_B"  = "#009E73",
  "NOISE"   = "grey60",
  "COMPLEX" = "#CC79A7"
)

conf_colors <- c("high" = "#2ca02c", "medium" = "#ff7f0e", "low" = "#d62728")

# ── Per-candidate plots ──────────────────────────────────────────────────
candidates <- if (!is.na(cid_filter)) {
  as.integer(cid_filter)
} else {
  unique(scoring$candidate_id)
}

for (cid in candidates) {
  sub <- scoring[candidate_id == cid]
  if (nrow(sub) < 5) next
  chr <- sub$chrom[1]

  if (!all(c("PC1", "PC2") %in% names(sub))) next

  # 1. PCA by major state
  p1 <- ggplot(sub, aes(x = PC1, y = PC2, color = provisional_major_state)) +
    geom_point(size = 2, alpha = 0.75) +
    scale_color_manual(values = state_colors) +
    theme_bw(base_size = 10) +
    labs(title = paste0("Candidate ", cid, " — PCA by major state"),
         x = "PC1", y = "PC2", color = "State")

  # 2. PCA by confidence
  p2 <- ggplot(sub, aes(x = PC1, y = PC2, color = provisional_confidence)) +
    geom_point(size = 2, alpha = 0.75) +
    scale_color_manual(values = conf_colors) +
    theme_bw(base_size = 10) +
    labs(title = paste0("Candidate ", cid, " — PCA by confidence"),
         x = "PC1", y = "PC2", color = "Confidence")

  # 3. Evidence scores
  score_cols <- c("evidence_coherence", "evidence_intermediate", "evidence_noise")
  score_cols <- intersect(score_cols, names(sub))

  if (length(score_cols) > 0) {
    score_long <- melt(sub[, c("sample", "provisional_major_state", ..score_cols)],
                       id.vars = c("sample", "provisional_major_state"),
                       variable.name = "evidence_type",
                       value.name = "score")

    p3 <- ggplot(score_long[is.finite(score)],
                 aes(x = score, fill = provisional_major_state)) +
      geom_histogram(bins = 25, alpha = 0.7, position = "identity") +
      facet_wrap(~evidence_type, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = state_colors) +
      theme_bw(base_size = 10) +
      labs(title = "Evidence score distributions",
           x = "Score", y = "Count", fill = "State")

    combined <- (p1 + p2) / p3 +
      plot_layout(heights = c(1, 1.2)) +
      plot_annotation(
        title = paste0("Candidate ", cid, " (", chr, ") — state assignments"),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )
  } else {
    combined <- p1 + p2 +
      plot_annotation(
        title = paste0("Candidate ", cid, " (", chr, ") — state assignments"),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )
  }

  ggsave(file.path(outdir, paste0("states_candidate_", cid, "_", chr, ".pdf")),
         combined, width = 12, height = 10)
  ggsave(file.path(outdir, paste0("states_candidate_", cid, "_", chr, ".png")),
         combined, width = 12, height = 10, dpi = 400)
}

# ── Global classification summary ────────────────────────────────────────
if (nrow(classification) > 0 && "candidate_classification" %in% names(classification)) {
  class_colors <- c(
    "clean_inversion_like"        = "#2ca02c",
    "complex_inversion_like"      = "#ff7f0e",
    "non_inversion_or_ambiguous"  = "#d62728",
    "unclassified"                = "grey60"
  )

  p_summary <- ggplot(classification, aes(x = candidate_classification, fill = candidate_classification)) +
    geom_bar(alpha = 0.8) +
    scale_fill_manual(values = class_colors) +
    theme_bw(base_size = 12) +
    labs(title = "Candidate classification summary",
         subtitle = paste0("n = ", nrow(classification), " candidates"),
         x = "Classification", y = "Count") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(outdir, "candidate_classification_summary.pdf"),
         p_summary, width = 7, height = 5)
  ggsave(file.path(outdir, "candidate_classification_summary.png"),
         p_summary, width = 7, height = 5, dpi = 400)

  # State proportion per candidate (stacked bar)
  state_summary <- scoring[, .N, by = .(candidate_id, provisional_major_state)]
  state_summary[, candidate_id := factor(candidate_id, levels = sort(unique(candidate_id)))]

  p_stacked <- ggplot(state_summary,
                      aes(x = candidate_id, y = N, fill = provisional_major_state)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = state_colors) +
    theme_bw(base_size = 10) +
    labs(title = "State proportions across all candidates",
         x = "Candidate ID", y = "Proportion",
         fill = "State") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

  ggsave(file.path(outdir, "candidate_state_proportions.pdf"),
         p_stacked, width = max(8, nrow(classification) * 0.15 + 3), height = 5)
  ggsave(file.path(outdir, "candidate_state_proportions.png"),
         p_stacked, width = max(8, nrow(classification) * 0.15 + 3), height = 5, dpi = 400)
}

message("[DONE] Scoring and state plots written to: ", outdir)
