#!/usr/bin/env Rscript
###############################################################################
# plot_relatedness_summary.R — Composite pruning summary figure
#
# 4-panel layout:
#   A. Theta distribution histogram with threshold lines
#   B. Relationship category barplot with counts
#   C. Admixture K=best with relatedness strip
#   D. NAToRA pruning summary table
#
# Usage:
#   Rscript plot_relatedness_summary.R \
#     --relatedness-file 06_relatedness/catfish_226_relatedness.res \
#     --ancestry-file structure_results/.../..._sample_ancestry.tsv \
#     --pruned-samples 06_relatedness/pruned_samples.txt \
#     --natora-summary 06_relatedness/natora_cutoff_summary.tsv \
#     --best-k 8 \
#     --theme-file utils/theme_systems_plate.R \
#     --out-dir structure_results/figures/
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(optparse)
})

option_list <- list(
  make_option("--relatedness-file", type = "character"),
  make_option("--ancestry-file", type = "character", default = NULL),
  make_option("--pruned-samples", type = "character", default = NULL),
  make_option("--natora-summary", type = "character", default = NULL),
  make_option("--best-k", type = "integer", default = 8),
  make_option("--theme-file", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--prefix", type = "character", default = "relatedness_summary"),
  make_option("--width", type = "numeric", default = 16),
  make_option("--height", type = "numeric", default = 12)
)
opt <- parse_args(OptionParser(option_list = option_list))

source(opt[["theme-file"]])

dir.create(opt[["out-dir"]], showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD DATA
# =============================================================================

res <- fread(opt[["relatedness-file"]])
ida_col <- if ("ida" %in% names(res)) "ida" else names(res)[3]
idb_col <- if ("idb" %in% names(res)) "idb" else names(res)[4]
theta_col <- if ("theta" %in% names(res)) "theta" else names(res)[5]

df <- data.table(
  ind1 = as.character(res[[ida_col]]),
  ind2 = as.character(res[[idb_col]]),
  theta = as.numeric(res[[theta_col]])
)

df[, relationship := fifelse(
  theta >= 0.354, "Duplicate/MZ",
  fifelse(theta >= 0.177, "First degree",
  fifelse(theta >= 0.0884, "Second degree",
  fifelse(theta >= 0.0442, "Third degree", "Unrelated"))))]

df[, relationship := factor(relationship,
  levels = c("Duplicate/MZ", "First degree", "Second degree",
             "Third degree", "Unrelated"))]

# Category colors matching plate theme
cat_colors <- c(
  "Duplicate/MZ"  = "#7B3294",
  "First degree"  = "#5AA5C9",
  "Second degree" = "#8ABF6A",
  "Third degree"  = "#E8A44A",
  "Unrelated"     = "#D9D9D9"
)

# =============================================================================
# PANEL A: Theta histogram
# =============================================================================

pA <- ggplot(df, aes(x = theta)) +
  geom_histogram(bins = 120, fill = .plate_accent$primary, color = "white",
                 linewidth = 0.15, alpha = 0.85) +
  geom_vline(xintercept = 0.354, linetype = "dashed", color = "#7B3294", linewidth = 0.5) +
  geom_vline(xintercept = 0.177, linetype = "dashed", color = "#5AA5C9", linewidth = 0.5) +
  geom_vline(xintercept = 0.0884, linetype = "dashed", color = "#8ABF6A", linewidth = 0.5) +
  geom_vline(xintercept = 0.0442, linetype = "dashed", color = "#E8A44A", linewidth = 0.5) +
  annotate("text", x = 0.354, y = Inf, label = "Dup/MZ", vjust = 1.5, hjust = -0.1,
           size = 2.5, color = "#7B3294", fontface = "bold") +
  annotate("text", x = 0.177, y = Inf, label = "1st", vjust = 1.5, hjust = -0.1,
           size = 2.5, color = "#5AA5C9", fontface = "bold") +
  annotate("text", x = 0.0884, y = Inf, label = "2nd", vjust = 1.5, hjust = -0.1,
           size = 2.5, color = "#8ABF6A", fontface = "bold") +
  annotate("text", x = 0.0442, y = Inf, label = "3rd", vjust = 1.5, hjust = -0.1,
           size = 2.5, color = "#E8A44A", fontface = "bold") +
  theme_plate(grid = "y") +
  labs(x = expression(paste("Kinship coefficient (", theta, ")")),
       y = "Pairs")

# =============================================================================
# PANEL B: Relationship category barplot
# =============================================================================

pair_counts <- df[, .N, by = relationship]
# Ensure all levels
all_levels <- data.table(relationship = factor(
  c("Duplicate/MZ", "First degree", "Second degree", "Third degree", "Unrelated"),
  levels = c("Duplicate/MZ", "First degree", "Second degree", "Third degree", "Unrelated")))
pair_counts <- merge(all_levels, pair_counts, by = "relationship", all.x = TRUE)
pair_counts[is.na(N), N := 0]

pB <- ggplot(pair_counts, aes(x = relationship, y = N, fill = relationship)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.4, size = 3,
            fontface = "bold", color = TEXT_BODY) +
  scale_fill_manual(values = cat_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_plate(grid = "y") +
  labs(x = NULL, y = "Pairs") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# =============================================================================
# PANEL C: Key statistics
# =============================================================================

samples <- sort(unique(c(df$ind1, df$ind2)))
n_samples <- length(samples)

n_pruned <- 0
if (!is.null(opt[["pruned-samples"]]) && file.exists(opt[["pruned-samples"]])) {
  pruned_ids <- trimws(readLines(opt[["pruned-samples"]], warn = FALSE))
  pruned_ids <- pruned_ids[nzchar(pruned_ids)]
  n_pruned <- length(pruned_ids)
}

stats_text <- paste0(
  "Total samples: ", n_samples, "\n",
  "Total pairs: ", format(nrow(df), big.mark = ","), "\n",
  "Retained after pruning: ", n_pruned, "\n",
  "Removed: ", n_samples - n_pruned, "\n\n",
  "Theta summary:\n",
  "  Mean: ", round(mean(df$theta, na.rm = TRUE), 4), "\n",
  "  Median: ", round(median(df$theta, na.rm = TRUE), 4), "\n",
  "  Max: ", round(max(df$theta, na.rm = TRUE), 4)
)

pC <- wrap_elements(
  plate_callout(stats_text, type = "finding")
)

# =============================================================================
# PANEL D: NAToRA summary (if available)
# =============================================================================

pD <- NULL
if (!is.null(opt[["natora-summary"]]) && file.exists(opt[["natora-summary"]])) {
  natora <- fread(opt[["natora-summary"]])
  if (nrow(natora) > 0) {
    pD_dt <- natora[, .(group, cutoff_label, removed_n, retained_n, removed_pct)]
    pD <- wrap_elements(
      plate_table(pD_dt, title = "NAToRA multi-cutoff summary")
    )
  }
}

# =============================================================================
# COMPOSE
# =============================================================================

if (!is.null(pD)) {
  fig <- (add_panel_label(pA, "A") | add_panel_label(pB, "B")) /
         (pC | pD) +
    plot_layout(heights = c(1.2, 1))
} else {
  fig <- (add_panel_label(pA, "A") | add_panel_label(pB, "B")) /
         pC +
    plot_layout(heights = c(1.2, 0.8))
}

fig <- fig + plot_annotation(
  title = "Relatedness and pruning summary",
  subtitle = paste0(n_samples, " samples, NAToRA + greedy first-degree pruning → ",
                    n_pruned, " retained"),
  theme = theme(
    plot.title = element_text(size = SIZE_MAIN_TITLE, face = "bold", color = TEXT_TITLE),
    plot.subtitle = element_text(size = SIZE_SUBTITLE, color = TEXT_SUBTITLE),
    plot.background = element_rect(fill = PLATE_BG, color = NA)
  )
)

# =============================================================================
# SAVE
# =============================================================================

out_pdf <- file.path(opt[["out-dir"]], paste0(opt$prefix, ".pdf"))
ggsave(out_pdf, fig, width = opt$width, height = opt$height,
       device = cairo_pdf, bg = PLATE_BG)
ggsave(sub("\\.pdf$", ".png", out_pdf), fig,
       width = opt$width, height = opt$height, dpi = 300, bg = PLATE_BG)

cat("[DONE]", out_pdf, "\n")
