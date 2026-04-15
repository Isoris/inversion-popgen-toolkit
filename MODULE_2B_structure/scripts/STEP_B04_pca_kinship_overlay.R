#!/usr/bin/env Rscript
###############################################################################
# plot_pca_kinship_overlay.R — PCA scatter with kinship edges + ancestry hulls
#
# PC1 vs PC2 with edges for first-degree relatives, nodes colored by K=best
# ancestry, convex hull per ancestry group. Refined from existing version
# with theme_systems_plate.R styling.
#
# Usage:
#   Rscript plot_pca_kinship_overlay.R \
#     --pcangsd-cov structure_results/pcangsd_all226/pcangsd.cov \
#     --relatedness-file 06_relatedness/catfish_226_relatedness.res \
#     --ancestry-file structure_results/.../..._sample_ancestry.tsv \
#     --best-k 8 \
#     --sample-file samples.txt \
#     --theme-file utils/theme_systems_plate.R \
#     --out-dir structure_results/figures/
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--pcangsd-cov", type = "character",
              help = "PCAngsd covariance matrix file"),
  make_option("--relatedness-file", type = "character"),
  make_option("--ancestry-file", type = "character"),
  make_option("--best-k", type = "integer", default = 8),
  make_option("--sample-file", type = "character"),
  make_option("--theme-file", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--prefix", type = "character", default = "pca_kinship_overlay"),
  make_option("--theta-edge", type = "numeric", default = 0.177,
              help = "Theta threshold for drawing edges [default: 0.177]"),
  make_option("--width", type = "numeric", default = 12),
  make_option("--height", type = "numeric", default = 10)
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
# PCA FROM COVARIANCE MATRIX
# =============================================================================

cov_mat <- as.matrix(read.table(opt[["pcangsd-cov"]], header = FALSE))
samples <- trimws(readLines(opt[["sample-file"]], warn = FALSE))
samples <- samples[nzchar(samples)]
n <- length(samples)

stopifnot(nrow(cov_mat) == n)

eig <- eigen(cov_mat)
pcs <- eig$vectors
pct_var <- 100 * eig$values / sum(eig$values)

pc_dt <- data.table(
  sample = samples,
  PC1 = pcs[, 1],
  PC2 = pcs[, 2]
)

# =============================================================================
# ANCESTRY COLORS
# =============================================================================

anc_dt <- fread(opt[["ancestry-file"]])
anc_best <- anc_dt[K == best_K]
pc_dt <- merge(pc_dt, anc_best[, .(sample, cluster_index, cluster_label, qmax, color_hex)],
               by = "sample", all.x = TRUE)
pc_dt[is.na(cluster_label), cluster_label := "NA"]
pc_dt[is.na(color_hex), color_hex := "#CCCCCC"]
pc_dt[is.na(qmax), qmax := 0.5]

# =============================================================================
# KINSHIP EDGES
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
edges <- pairs[theta >= opt[["theta-edge"]]]

# Build edge segment coordinates
edge_coords <- merge(edges, pc_dt[, .(sample, x1 = PC1, y1 = PC2)],
                     by.x = "a", by.y = "sample", all.x = TRUE)
edge_coords <- merge(edge_coords, pc_dt[, .(sample, x2 = PC1, y2 = PC2)],
                     by.x = "b", by.y = "sample", all.x = TRUE)
edge_coords <- edge_coords[!is.na(x1) & !is.na(x2)]

# Within vs between ancestry
edge_coords <- merge(edge_coords,
  pc_dt[, .(sample, comp_a = cluster_label)], by.x = "a", by.y = "sample")
edge_coords <- merge(edge_coords,
  pc_dt[, .(sample, comp_b = cluster_label)], by.x = "b", by.y = "sample")
edge_coords[, edge_type := fifelse(comp_a == comp_b, "within", "between")]

n_within <- sum(edge_coords$edge_type == "within")
n_between <- sum(edge_coords$edge_type == "between")

# =============================================================================
# PLOT
# =============================================================================

color_map <- setNames(
  pal_k20[seq_len(best_K)],
  paste0("Q", seq_len(best_K))
)
color_map[["NA"]] <- "#CCCCCC"

p <- ggplot() +
  # Convex hulls
  stat_ellipse(data = pc_dt[cluster_label != "NA"],
               aes(x = PC1, y = PC2, color = cluster_label),
               type = "t", level = 0.9, linewidth = 0.3, alpha = 0.5,
               show.legend = FALSE) +
  # Edges
  geom_segment(data = edge_coords,
               aes(x = x1, y = y1, xend = x2, yend = y2, linetype = edge_type),
               alpha = 0.3, linewidth = 0.3, color = "#5C5C5C") +
  scale_linetype_manual(values = c(within = "solid", between = "dashed"),
                        name = "Edge type") +
  # Points
  geom_point(data = pc_dt,
             aes(x = PC1, y = PC2, color = cluster_label, alpha = qmax),
             size = 2.2) +
  scale_color_manual(values = color_map, name = paste0("K=", best_K, " ancestry")) +
  scale_alpha_continuous(range = c(0.4, 1.0), guide = "none") +
  theme_plate(grid = "xy", panel_border = TRUE) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pct_var[1]),
    y = sprintf("PC2 (%.1f%%)", pct_var[2]),
    title = "PCA with kinship overlay",
    subtitle = sprintf("%d within-ancestry edges, %d between-ancestry edges (theta >= %.3f)",
                       n_within, n_between, opt[["theta-edge"]])
  )

# =============================================================================
# SAVE
# =============================================================================

out_pdf <- file.path(opt[["out-dir"]], paste0(opt$prefix, ".pdf"))
ggsave(out_pdf, p, width = opt$width, height = opt$height,
       device = cairo_pdf, bg = PLATE_BG)
ggsave(sub("\\.pdf$", ".png", out_pdf), p,
       width = opt$width, height = opt$height, dpi = 300, bg = PLATE_BG)

cat("[DONE]", out_pdf, "\n")
