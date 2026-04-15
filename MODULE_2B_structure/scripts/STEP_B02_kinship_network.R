#!/usr/bin/env Rscript
###############################################################################
# plot_kinship_network.R — ggraph-based kinship network (REWRITE)
#
# Nodes colored by dominant K-best ancestry (pal_ancestry_k8/k20), sized by
# degree, darkness by max(Q). Edges: solid for Dup/MZ, dashed for 1st-degree.
# Component hulls: subtle fill. Labels only for high-degree nodes.
#
# Usage:
#   Rscript plot_kinship_network.R \
#     --relatedness-file 06_relatedness/catfish_226_relatedness.res \
#     --ancestry-file structure_results/.../wholegenome_thin500_all226_sample_ancestry.tsv \
#     --best-k 8 \
#     --theme-file utils/theme_systems_plate.R \
#     --out-dir structure_results/figures/
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggraph)
  library(tidygraph)
  library(igraph)
  library(optparse)
})

option_list <- list(
  make_option("--relatedness-file", type = "character"),
  make_option("--ancestry-file", type = "character"),
  make_option("--best-k", type = "integer", default = 8),
  make_option("--theme-file", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--prefix", type = "character", default = "kinship_network"),
  make_option("--theta-min", type = "numeric", default = 0.0884,
              help = "Min theta to draw an edge [default: 0.0884 = 2nd degree]"),
  make_option("--label-degree-min", type = "integer", default = 3,
              help = "Min degree to show node label [default: 3]"),
  make_option("--width", type = "numeric", default = 14),
  make_option("--height", type = "numeric", default = 12)
)
opt <- parse_args(OptionParser(option_list = option_list))

source(opt[["theme-file"]])

# =============================================================================
# PALETTE
# =============================================================================

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
# Determine column names
ida_col <- if ("ida" %in% names(res)) "ida" else names(res)[3]
idb_col <- if ("idb" %in% names(res)) "idb" else names(res)[4]
theta_col <- if ("theta" %in% names(res)) "theta" else names(res)[5]

pairs <- data.table(
  from = as.character(res[[ida_col]]),
  to = as.character(res[[idb_col]]),
  theta = as.numeric(res[[theta_col]])
)

# Filter to edges above threshold
edges <- pairs[theta >= opt[["theta-min"]]]

edges[, rel_class := fifelse(
  theta >= 0.354, "Dup/MZ",
  fifelse(theta >= 0.177, "First degree",
  fifelse(theta >= 0.0884, "Second degree", "Third degree")))]

cat("[INFO] Edges:", nrow(edges), "\n")

# Load ancestry for coloring
anc_dt <- fread(opt[["ancestry-file"]])
anc_best <- anc_dt[K == best_K]

# =============================================================================
# BUILD GRAPH
# =============================================================================

all_samples <- sort(unique(c(edges$from, edges$to)))

# Node attributes
node_dt <- data.table(name = all_samples)
node_dt <- merge(node_dt, anc_best[, .(sample, cluster_index, cluster_label, qmax, color_hex)],
                 by.x = "name", by.y = "sample", all.x = TRUE)
node_dt[is.na(cluster_index), cluster_index := 0]
node_dt[is.na(qmax), qmax := 0.5]
node_dt[is.na(color_hex), color_hex := "#CCCCCC"]

# Build igraph
g <- graph_from_data_frame(edges[, .(from, to)], directed = FALSE, vertices = node_dt)
E(g)$theta <- edges$theta
E(g)$rel_class <- edges$rel_class

# Node degree
V(g)$degree <- degree(g)
# Component membership
V(g)$comp_id <- components(g)$membership

# Convert to tidygraph
tg <- as_tbl_graph(g)

# =============================================================================
# PLOT
# =============================================================================

# Color nodes by ancestry component
node_colors <- setNames(node_dt$color_hex, node_dt$name)

# Edge aesthetics
edge_linetypes <- c("Dup/MZ" = "solid", "First degree" = "dashed",
                     "Second degree" = "dotted", "Third degree" = "dotted")
edge_colors <- c("Dup/MZ" = "#2D2D2D", "First degree" = "#6B6B6B",
                 "Second degree" = "#AAAAAA", "Third degree" = "#CCCCCC")

p <- ggraph(tg, layout = "stress") +
  # Component hulls (subtle)
  geom_mark_hull(
    aes(x = x, y = y, group = factor(comp_id),
        fill = factor(cluster_index)),
    alpha = 0.06, expand = unit(4, "mm"), radius = unit(3, "mm"),
    color = NA, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("0" = "#CCCCCC", setNames(pal_k20, seq_along(pal_k20))),
                    guide = "none") +
  # Edges
  geom_edge_link(
    aes(linetype = rel_class, edge_colour = rel_class, edge_alpha = theta),
    edge_width = 0.4, show.legend = TRUE
  ) +
  scale_edge_linetype_manual(values = edge_linetypes, name = "Relationship") +
  scale_edge_colour_manual(values = edge_colors, name = "Relationship") +
  scale_edge_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
  # Nodes
  geom_node_point(
    aes(size = degree, color = name, alpha = qmax),
    show.legend = FALSE
  ) +
  scale_color_manual(values = node_colors, guide = "none") +
  scale_alpha_continuous(range = c(0.4, 1.0), guide = "none") +
  scale_size_continuous(range = c(1.5, 6), guide = "none") +
  # Labels (only high-degree nodes)
  geom_node_text(
    aes(label = ifelse(degree >= opt[["label-degree-min"]], name, "")),
    size = 2, repel = TRUE, max.overlaps = 15,
    color = TEXT_BODY, family = "sans"
  ) +
  theme_plate(grid = "none") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.background = element_rect(fill = PLATE_BG, color = NA),
    legend.position = "bottom"
  ) +
  labs(
    title = "Kinship network",
    subtitle = paste0(length(all_samples), " samples, ",
                      nrow(edges), " edges (theta >= ", opt[["theta-min"]], "), ",
                      "colored by K=", best_K, " ancestry")
  )

# =============================================================================
# SAVE
# =============================================================================

dir.create(opt[["out-dir"]], showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(opt[["out-dir"]], paste0(opt$prefix, ".pdf"))
out_png <- file.path(opt[["out-dir"]], paste0(opt$prefix, ".png"))

ggsave(out_pdf, p, width = opt$width, height = opt$height,
       device = cairo_pdf, bg = PLATE_BG)
ggsave(out_png, p, width = opt$width, height = opt$height,
       dpi = 300, bg = PLATE_BG)

cat("[DONE] Saved:", out_pdf, "\n")
