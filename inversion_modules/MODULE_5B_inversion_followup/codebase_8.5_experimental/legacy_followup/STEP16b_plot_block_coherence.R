#!/usr/bin/env Rscript

# =============================================================================
# STEP16b_plot_block_coherence.R
#
# Visualization suite for block coherence metrics from STEP16.
# Per candidate:
#   1. Histogram/density of block_coherence
#   2. Scatter: mean_A_support vs mean_B_support colored by PCA group
#   3. Scatter: block_coherence vs support_margin colored by provisional label
#   4. Histogram of switch_rate
#
# Usage:
#   Rscript STEP16b_plot_block_coherence.R \
#     <coherence_per_sample.tsv.gz> <group_table.tsv.gz> <outdir> \
#     [candidate_id=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript STEP16b_... <coherence.tsv.gz> <groups.tsv.gz> <outdir> [candidate_id=all]")
}

coh_file   <- args[1]
group_file <- args[2]
outdir     <- args[3]
cid_filter <- if (length(args) >= 4 && args[4] != "all") args[4] else NA_character_

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

coh <- fread(coh_file)
grp <- fread(group_file)

# Merge group info
dt <- merge(coh, grp[, .(candidate_id, sample, group_norm, group_label)],
            by = c("candidate_id", "sample"), all.x = TRUE)

if (!is.na(cid_filter)) {
  dt <- dt[as.character(candidate_id) == cid_filter]
}

candidates <- unique(dt$candidate_id)
message("[INFO] Plotting block coherence for ", length(candidates), " candidates")

label_colors <- c(
  "coherent_A"  = "#0072B2",
  "coherent_B"  = "#009E73",
  "mixed_or_het" = "#D55E00",
  "ambiguous"   = "#CC79A7",
  "low_info"    = "grey60"
)

group_colors <- c("Homo_1" = "#0072B2", "Het" = "#D55E00", "Homo_2" = "#009E73",
                  "G1" = "#0072B2", "G2" = "#D55E00", "G3" = "#009E73")

for (cid in candidates) {
  sub <- dt[candidate_id == cid]
  chr <- sub$chrom[1]
  n <- nrow(sub)

  if (n < 5) {
    message("[SKIP] Candidate ", cid, ": too few samples (", n, ")")
    next
  }

  # 1. Block coherence density
  p1 <- ggplot(sub[is.finite(block_coherence)],
               aes(x = block_coherence, fill = provisional_label)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = label_colors, na.value = "grey80") +
    theme_bw(base_size = 10) +
    labs(title = "Block coherence distribution",
         x = "Block coherence (1 - switch rate)", y = "Count",
         fill = "Label")

  # 2. A support vs B support
  p2 <- ggplot(sub[is.finite(mean_A_support) & is.finite(mean_B_support)],
               aes(x = mean_A_support, y = mean_B_support)) +
    geom_point(aes(color = group_label), size = 1.8, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    scale_color_manual(values = group_colors, na.value = "grey60") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw(base_size = 10) +
    labs(title = "Arrangement support",
         x = "Mean A support", y = "Mean B support",
         color = "PCA group")

  # 3. Coherence vs margin
  p3 <- ggplot(sub[is.finite(block_coherence) & is.finite(support_margin)],
               aes(x = support_margin, y = block_coherence, color = provisional_label)) +
    geom_point(size = 1.8, alpha = 0.7) +
    scale_color_manual(values = label_colors, na.value = "grey80") +
    theme_bw(base_size = 10) +
    labs(title = "Coherence vs support margin",
         x = "Support margin |A-B|", y = "Block coherence",
         color = "Label")

  # 4. Switch rate histogram
  p4 <- ggplot(sub[is.finite(switch_rate)],
               aes(x = switch_rate, fill = group_label)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = group_colors, na.value = "grey80") +
    theme_bw(base_size = 10) +
    labs(title = "Switch rate distribution",
         x = "Switch rate", y = "Count",
         fill = "PCA group")

  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(title = paste0("Candidate ", cid, " (", chr, ") — block coherence metrics"),
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))

  ggsave(file.path(outdir, paste0("coherence_candidate_", cid, "_", chr, ".pdf")),
         combined, width = 12, height = 10)
  ggsave(file.path(outdir, paste0("coherence_candidate_", cid, "_", chr, ".png")),
         combined, width = 12, height = 10, dpi = 400)
}

message("[DONE] Block coherence plots written to: ", outdir)
