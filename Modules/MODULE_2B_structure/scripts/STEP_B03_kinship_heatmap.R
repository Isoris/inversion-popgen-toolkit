#!/usr/bin/env Rscript
###############################################################################
# plot_kinship_heatmap.R — Annotated kinship heatmap (REWRITE)
#
# 226×226 matrix ordered by PCAngsd tree, with:
#   - Left/top annotation bars (dominant ancestry at K=best)
#   - Retained/pruned status bar
#   - Horizontal reference lines at theta thresholds
#   - White → warm amber → dark brown color scale
#
# Usage:
#   Rscript plot_kinship_heatmap.R \
#     --relatedness-file 06_relatedness/catfish_226_relatedness.res \
#     --ancestry-file structure_results/.../..._sample_ancestry.tsv \
#     --best-k 8 \
#     --pruned-samples 06_relatedness/pruned_samples.txt \
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
  make_option("--ancestry-file", type = "character"),
  make_option("--best-k", type = "integer", default = 8),
  make_option("--pruned-samples", type = "character", default = NULL),
  make_option("--theme-file", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--prefix", type = "character", default = "kinship_heatmap"),
  make_option("--width", type = "numeric", default = 14),
  make_option("--height", type = "numeric", default = 13)
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

# =============================================================================
# LOAD DATA
# =============================================================================

res <- fread(opt[["relatedness-file"]])
ida_col <- if ("ida" %in% names(res)) "ida" else names(res)[3]
idb_col <- if ("idb" %in% names(res)) "idb" else names(res)[4]
theta_col <- if ("theta" %in% names(res)) "theta" else names(res)[5]

pairs <- data.table(
  a = as.character(res[[ida_col]]),
  b = as.character(res[[idb_col]]),
  theta = as.numeric(res[[theta_col]])
)

samples <- sort(unique(c(pairs$a, pairs$b)))
n <- length(samples)
cat("[INFO] Building", n, "x", n, "matrix\n")

# Build matrix
theta_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))
for (i in seq_len(nrow(pairs))) {
  theta_mat[pairs$a[i], pairs$b[i]] <- pairs$theta[i]
  theta_mat[pairs$b[i], pairs$a[i]] <- pairs$theta[i]
}
diag(theta_mat) <- 0.5

# Cluster for ordering
d <- as.dist(1 - theta_mat)
hc <- hclust(d, method = "average")
sample_order <- hc$labels[hc$order]

# Load ancestry
anc_dt <- fread(opt[["ancestry-file"]])
anc_best <- anc_dt[K == best_K]

# Melt for ggplot
theta_melt <- data.table(
  Ind1 = rep(samples, n),
  Ind2 = rep(samples, each = n),
  theta = as.vector(theta_mat)
)
theta_melt[, Ind1 := factor(Ind1, levels = sample_order)]
theta_melt[, Ind2 := factor(Ind2, levels = rev(sample_order))]

# =============================================================================
# HEATMAP
# =============================================================================

p_heat <- ggplot(theta_melt, aes(x = Ind1, y = Ind2, fill = theta)) +
  geom_tile(color = NA) +
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#FFF5EB", "#FDD49E", "#FDBB84",
               "#FC8D59", "#E34A33", "#B30000"),
    name = expression(theta),
    limits = c(0, 0.5)
  ) +
  # Threshold reference lines (as geom_hline/vline on the sample axis)
  theme_plate(grid = "none") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.key.height = unit(25, "pt"),
    legend.key.width = unit(8, "pt"),
    plot.margin = margin(0, 4, 0, 0)
  ) +
  coord_fixed()

# =============================================================================
# ANCESTRY ANNOTATION BARS
# =============================================================================

anc_order <- data.table(sample = sample_order, x = seq_along(sample_order))
anc_order <- merge(anc_order, anc_best[, .(sample, cluster_index, color_hex)],
                   by = "sample", all.x = TRUE)
anc_order[is.na(cluster_index), cluster_index := 0]
anc_order[is.na(color_hex), color_hex := "#CCCCCC"]

# Top bar
p_top <- ggplot(anc_order, aes(x = x, y = 1, fill = factor(cluster_index))) +
  geom_tile(height = 1, color = NA) +
  scale_fill_manual(
    values = c("0" = "#CCCCCC", setNames(pal_k20[seq_len(best_K)], seq_len(best_K))),
    name = paste0("K=", best_K)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 4, 0, 0))

# Left bar (same data, rotated)
p_left <- ggplot(anc_order, aes(x = 1, y = n - x + 1, fill = factor(cluster_index))) +
  geom_tile(width = 1, color = NA) +
  scale_fill_manual(
    values = c("0" = "#CCCCCC", setNames(pal_k20[seq_len(best_K)], seq_len(best_K))),
    guide = "none"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# =============================================================================
# PRUNED STATUS BAR
# =============================================================================

pruned_ids <- character(0)
if (!is.null(opt[["pruned-samples"]]) && file.exists(opt[["pruned-samples"]])) {
  pruned_ids <- trimws(readLines(opt[["pruned-samples"]], warn = FALSE))
  pruned_ids <- pruned_ids[nzchar(pruned_ids)]
}

anc_order[, status := fifelse(sample %in% pruned_ids, "retained", "pruned")]

p_status <- ggplot(anc_order, aes(x = x, y = 1, fill = status)) +
  geom_tile(height = 1, color = NA) +
  scale_fill_manual(values = c(retained = "#D9D9D9", pruned = "#5AA5C9"),
                    name = "Status") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 4, 0, 0))

# =============================================================================
# COMPOSE
# =============================================================================

# Layout:
#   [spacer] [status] [anc_top]
#   [anc_left] [heatmap]

layout <- "
  ##BB
  ##CC
  AADD
"

fig <- p_left + p_status + p_top + p_heat +
  plot_layout(
    design = layout,
    widths = c(0.02, 1),
    heights = c(0.015, 0.015, 1)
  ) +
  plot_annotation(
    title = "Kinship heatmap",
    subtitle = paste0(n, " samples, ordered by hierarchical clustering, ",
                      "colored by K=", best_K, " ancestry"),
    theme = theme(
      plot.title = element_text(size = SIZE_MAIN_TITLE, face = "bold", color = TEXT_TITLE),
      plot.subtitle = element_text(size = SIZE_SUBTITLE, color = TEXT_SUBTITLE),
      plot.background = element_rect(fill = PLATE_BG, color = NA)
    )
  )

# =============================================================================
# SAVE
# =============================================================================

dir.create(opt[["out-dir"]], showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(opt[["out-dir"]], paste0(opt$prefix, ".pdf"))
ggsave(out_pdf, fig, width = opt$width, height = opt$height,
       device = cairo_pdf, bg = PLATE_BG)
ggsave(sub("\\.pdf$", ".png", out_pdf), fig,
       width = opt$width, height = opt$height, dpi = 300, bg = PLATE_BG)

cat("[DONE]", out_pdf, "\n")
