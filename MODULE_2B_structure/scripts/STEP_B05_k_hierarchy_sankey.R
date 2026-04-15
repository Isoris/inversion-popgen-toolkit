#!/usr/bin/env Rscript
###############################################################################
# plot_k_hierarchy_tree.R — Sankey/alluvial founder merge tree
#
# Visualizes the K-hierarchy merge tree: K=2 at top, K=20 at bottom.
# Each node = component at K, width = n_samples. Flows show merges.
# Color = component identity at K=best, tracked through hierarchy.
#
# Usage:
#   Rscript plot_k_hierarchy_tree.R \
#     --tree-file structure_results/.../..._merge_tree.tsv \
#     --tracking-file structure_results/.../..._component_tracking.tsv \
#     --best-k 8 \
#     --theme-file utils/theme_systems_plate.R \
#     --out-dir structure_results/figures/ \
#     --file-prefix wholegenome_thin500_all226
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
  library(optparse)
})

option_list <- list(
  make_option("--tree-file", type = "character"),
  make_option("--tracking-file", type = "character"),
  make_option("--best-k", type = "integer", default = 8),
  make_option("--theme-file", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--file-prefix", type = "character"),
  make_option("--width", type = "numeric", default = 14),
  make_option("--height", type = "numeric", default = 20)
)
opt <- parse_args(OptionParser(option_list = option_list))

source(opt[["theme-file"]])

pal_k20 <- c(
  "#3B6FA0", "#CF6839", "#C44E52", "#6BA08E", "#5A8F4A",
  "#C9A83E", "#8B6DAD", "#E8919C", "#5C7A3A", "#B07850",
  "#4A8C9F", "#9E6B8A", "#7A8B3C", "#C47A5E", "#5E7FAA",
  "#A06B4F", "#6B9E7A", "#B0855A", "#7E6EA0", "#A0856B"
)

best_K <- opt[["best-k"]]
dir.create(opt[["out-dir"]], showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD COMPONENT TRACKING
# =============================================================================

tracking <- fread(opt[["tracking-file"]])
# Columns: sample_index, sample, K2_component, K3_component, ..., K20_component

k_cols <- grep("^K[0-9]+_component$", names(tracking), value = TRUE)
K_values <- as.integer(gsub("K([0-9]+)_component", "\\1", k_cols))
K_values <- sort(K_values)

cat("[INFO] K values in tracking:", paste(K_values, collapse = ","), "\n")

# =============================================================================
# BUILD ALLUVIAL DATA
# =============================================================================

# Melt to long format
id_cols <- c("sample_index", "sample")
alluv_long <- melt(tracking, id.vars = id_cols, measure.vars = k_cols,
                   variable.name = "K_col", value.name = "component")
alluv_long[, K := as.integer(gsub("K([0-9]+)_component", "\\1", K_col))]

# For coloring: use the component assignment at best_K
best_col <- paste0("K", best_K, "_component")
if (best_col %in% names(tracking)) {
  alluv_long <- merge(alluv_long,
                      tracking[, .(sample_index, best_comp = get(best_col))],
                      by = "sample_index")
} else {
  alluv_long[, best_comp := 1]
}

# Create a stratum label that's unique per K×component
alluv_long[, stratum := paste0("K", K, "_Q", component)]

# Count: for each K × component × best_comp, how many samples?
alluv_counts <- alluv_long[, .N, by = .(K, component, best_comp)]

# For ggalluvial we need wide format: one row per sample, one column per K
alluv_wide <- dcast(alluv_long, sample_index + best_comp ~ K_col,
                    value.var = "component")

# =============================================================================
# PLOT — ALLUVIAL
# =============================================================================

# ggalluvial needs the K columns as axes
# We'll use the long format with ggalluvial's to_lodes_form

lodes <- to_lodes_form(
  alluv_wide,
  axes = k_cols,
  id = "sample_index"
)
lodes[, K := as.integer(gsub("K([0-9]+)_component", "\\1", x))]
lodes[, stratum_label := paste0("Q", stratum)]

# Color by best_comp
color_map <- setNames(pal_k20[seq_len(max(alluv_wide$best_comp, na.rm = TRUE))],
                      seq_len(max(alluv_wide$best_comp, na.rm = TRUE)))

p <- ggplot(lodes,
       aes(x = x, stratum = stratum, alluvium = sample_index,
           fill = factor(best_comp), label = stratum_label)) +
  geom_flow(alpha = 0.35, color = NA, width = 0.35) +
  geom_stratum(alpha = 0.85, width = 0.35, color = "#D4D4D4", linewidth = 0.2) +
  geom_text(stat = "stratum", size = 2, color = TEXT_BODY) +
  scale_fill_manual(values = color_map,
                    name = paste0("K=", best_K, " ancestry"),
                    labels = paste0("Q", names(color_map))) +
  scale_x_discrete(labels = paste0("K=", K_values), expand = c(0.05, 0.05)) +
  theme_plate(grid = "none") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Founder lineage merge tree",
    subtitle = paste0("K=", min(K_values), " (top) to K=", max(K_values),
                      " (bottom) — color tracks K=", best_K, " identity")
  )

# Add best-K indicator line
best_K_pos <- match(paste0("K", best_K, "_component"), k_cols)
if (!is.na(best_K_pos)) {
  p <- p + annotate("segment",
    x = best_K_pos - 0.45, xend = best_K_pos + 0.45,
    y = -Inf, yend = -Inf,
    color = "#C44E52", linewidth = 2
  ) + annotate("text",
    x = best_K_pos, y = -Inf,
    label = paste0("Best K=", best_K),
    vjust = 1.5, size = 3, color = "#C44E52", fontface = "bold"
  )
}

# =============================================================================
# SAVE
# =============================================================================

out_pdf <- file.path(opt[["out-dir"]],
                     paste0(opt[["file-prefix"]], "_k_hierarchy_sankey.pdf"))
ggsave(out_pdf, p, width = opt$width, height = opt$height,
       device = cairo_pdf, bg = PLATE_BG)
ggsave(sub("\\.pdf$", ".png", out_pdf), p,
       width = opt$width, height = opt$height, dpi = 300, bg = PLATE_BG)

cat("[DONE]", out_pdf, "\n")
