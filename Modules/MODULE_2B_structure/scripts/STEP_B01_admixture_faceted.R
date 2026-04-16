#!/usr/bin/env Rscript
###############################################################################
# plot_admixture_faceted.R — Full K=2..20 admixture composite figure
#
# Produces a single tall composite per thin panel:
#   - Relatedness strip (top)
#   - For each K: evalAdmix residual strip + admixture barplot
#   - Best K indicator
#   - Side panels: loglik vs K, mean|evalAdmix| vs K
#
# Uses theme_systems_plate.R exclusively for styling.
#
# Usage:
#   Rscript plot_admixture_faceted.R \
#     --results-dir structure_results/wholegenome_thin500_all226/ \
#     --eval-dir structure_results/evaladmix/ \
#     --relatedness-file 06_relatedness/catfish_226_relatedness.res \
#     --pruned-samples 06_relatedness/pruned_samples.txt \
#     --theme-file utils/theme_systems_plate.R \
#     --file-prefix wholegenome_thin500_all226 \
#     --out-dir structure_results/figures/
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(optparse)
})

option_list <- list(
  make_option("--results-dir", type = "character"),
  make_option("--eval-dir", type = "character", default = NULL),
  make_option("--relatedness-file", type = "character", default = NULL),
  make_option("--pruned-samples", type = "character", default = NULL),
  make_option("--theme-file", type = "character"),
  make_option("--file-prefix", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--pcangsd-tree", type = "character", default = NULL,
              help = "PCAngsd .tree file for sample ordering (optional)"),
  make_option("--width", type = "numeric", default = 16),
  make_option("--height-per-k", type = "numeric", default = 0.6,
              help = "Height in inches per K row [default: 0.6]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Source theme
source(opt[["theme-file"]])

results_dir <- opt[["results-dir"]]
eval_dir    <- opt[["eval-dir"]]
file_prefix <- opt[["file-prefix"]]
out_dir     <- opt[["out-dir"]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# UNIFIED PALETTE (same as theme_systems_plate.R pal_ancestry_k8 + extension)
# =============================================================================

pal_k20 <- c(
  "#3B6FA0", "#CF6839", "#C44E52", "#6BA08E", "#5A8F4A",
  "#C9A83E", "#8B6DAD", "#E8919C", "#5C7A3A", "#B07850",
  "#4A8C9F", "#9E6B8A", "#7A8B3C", "#C47A5E", "#5E7FAA",
  "#A06B4F", "#6B9E7A", "#B0855A", "#7E6EA0", "#A0856B"
)

get_palette <- function(K) {
  if (K <= length(pal_k20)) pal_k20[seq_len(K)]
  else colorRampPalette(pal_k20)(K)
}

pal_relatedness <- c(
  retained    = "#E8E8E8",
  dup_removed = "#7B3294",
  first_removed  = "#5AA5C9",
  second_removed = "#8ABF6A",
  third_removed  = "#E8A44A"
)

# =============================================================================
# LOAD BEST SEED TABLE
# =============================================================================

best_file <- file.path(results_dir, paste0(file_prefix, "_best_seed_by_K.tsv"))
if (!file.exists(best_file)) {
  # Try auto-discovery
  candidates <- list.files(results_dir, pattern = "_best_seed_by_K\\.tsv$", full.names = TRUE)
  if (length(candidates) == 0) stop("No best_seed_by_K.tsv found in: ", results_dir)
  best_file <- candidates[1]
}
best_dt <- fread(best_file)
cat("[INFO] Best seeds loaded: K=", paste(best_dt$K, collapse = ","), "\n")

K_values <- sort(best_dt$K)

# =============================================================================
# LOAD ALL BEST QOPT FILES + BUILD LONG-FORMAT DATA
# =============================================================================

read_qopt <- function(path, K) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  do.call(rbind, lapply(strsplit(lines, "[[:space:]]+"), as.numeric))
}

# Sample order: from PCAngsd tree at best K if available, else from best-K qopt
sample_file <- file.path(results_dir, paste0(file_prefix, "_sample_order.tsv"))
if (file.exists(sample_file)) {
  sample_dt <- fread(sample_file)
  samples <- sample_dt$sample
} else {
  # Infer from first qopt
  qf <- best_dt$q_file[1]
  samples <- paste0("S", seq_len(nrow(read_qopt(qf, best_dt$K[1]))))
}
n_samples <- length(samples)

# Determine best K (lowest mean_abs_resid among K >= 4)
if ("mean_abs_resid" %in% names(best_dt) && any(!is.na(best_dt$mean_abs_resid))) {
  viable <- best_dt[K >= 4 & !is.na(mean_abs_resid)]
  if (nrow(viable) > 0) {
    best_K <- viable[which.min(mean_abs_resid), K]
  } else {
    best_K <- best_dt[which.max(loglik), K]
  }
} else {
  best_K <- best_dt[which.max(loglik), K]
}
cat("[INFO] Best K:", best_K, "\n")

# Load best-K qopt for sample ordering
best_K_qopt <- read_qopt(best_dt[K == best_K, q_file], best_K)

# Order samples by dominant component at best K, then by max Q within component
dom_comp <- max.col(best_K_qopt, ties.method = "first")
max_q <- apply(best_K_qopt, 1, max)
sample_order <- order(dom_comp, -max_q)

# Build long-format admixture data for all K
admix_long <- rbindlist(lapply(K_values, function(K_val) {
  qf <- best_dt[K == K_val, q_file]
  if (!file.exists(qf)) return(NULL)
  qmat <- read_qopt(qf, K_val)

  dt <- data.table(
    sample_idx = rep(seq_len(n_samples), K_val),
    sample = rep(samples, K_val),
    K = K_val,
    component = rep(seq_len(K_val), each = n_samples),
    Q = as.vector(qmat)
  )
  # Apply consistent sample ordering
  dt[, plot_x := match(sample_idx, sample_order)]
  dt
}), use.names = TRUE, fill = TRUE)

# =============================================================================
# LOAD EVALADMIX DATA (if available)
# =============================================================================

eval_burden <- NULL
if (!is.null(eval_dir) && dir.exists(eval_dir)) {
  eval_list <- lapply(K_values, function(K_val) {
    stem <- sprintf("%s_K%02d_best", file_prefix, K_val)
    cor_file <- file.path(eval_dir, paste0(stem, ".corres.txt"))
    if (!file.exists(cor_file)) return(NULL)

    M <- tryCatch(
      as.matrix(read.table(cor_file, header = FALSE)),
      error = function(e) NULL
    )
    if (is.null(M) || nrow(M) != n_samples) return(NULL)

    diag(M) <- NA
    # Per-sample burden = mean absolute off-diagonal residual
    burden <- rowMeans(abs(M), na.rm = TRUE)

    data.table(
      sample_idx = seq_len(n_samples),
      K = K_val,
      eval_burden = burden,
      plot_x = match(seq_len(n_samples), sample_order)
    )
  })
  eval_burden <- rbindlist(eval_list, use.names = TRUE, fill = TRUE)
  if (!is.null(eval_burden) && nrow(eval_burden) > 0) {
    cat("[INFO] evalAdmix burdens loaded for K=",
        paste(unique(eval_burden$K), collapse = ","), "\n")
  }
}

# =============================================================================
# BUILD RELATEDNESS STRIP (if available)
# =============================================================================

rel_strip <- NULL
if (!is.null(opt[["relatedness-file"]]) && file.exists(opt[["relatedness-file"]])) {
  res <- fread(opt[["relatedness-file"]])
  # Classify each sample's relationship status
  # First load pruned samples
  pruned_ids <- character(0)
  if (!is.null(opt[["pruned-samples"]]) && file.exists(opt[["pruned-samples"]])) {
    pruned_ids <- trimws(readLines(opt[["pruned-samples"]], warn = FALSE))
    pruned_ids <- pruned_ids[nzchar(pruned_ids)]
  }

  # For each sample, find its strongest relationship
  # Use column names from ngsRelate output
  ida_col <- if ("ida" %in% names(res)) "ida" else names(res)[3]
  idb_col <- if ("idb" %in% names(res)) "idb" else names(res)[4]
  theta_col <- if ("theta" %in% names(res)) "theta" else names(res)[5]

  pairs <- data.table(
    a = as.character(res[[ida_col]]),
    b = as.character(res[[idb_col]]),
    theta = as.numeric(res[[theta_col]])
  )

  # Max theta per sample
  max_theta <- rbind(
    pairs[, .(max_theta = max(theta, na.rm = TRUE)), by = .(sample = a)],
    pairs[, .(max_theta = max(theta, na.rm = TRUE)), by = .(sample = b)]
  )[, .(max_theta = max(max_theta)), by = sample]

  # Classify
  max_theta[, rel_class := fifelse(
    max_theta >= 0.354, "dup_removed",
    fifelse(max_theta >= 0.177, "first_removed",
    fifelse(max_theta >= 0.0884, "second_removed",
    fifelse(max_theta >= 0.0442, "third_removed",
    "retained"))))]

  # Override: if sample is in pruned list, it's retained
  max_theta[sample %in% pruned_ids, rel_class := "retained"]
  max_theta[!sample %in% pruned_ids & rel_class != "retained",
            rel_class := rel_class]

  # For samples not in any close pair, mark as retained
  all_samples_dt <- data.table(sample = samples)
  rel_strip <- merge(all_samples_dt, max_theta, by = "sample", all.x = TRUE)
  rel_strip[is.na(rel_class), rel_class := "retained"]
  rel_strip[, sample_idx := match(sample, samples)]
  rel_strip[, plot_x := match(sample_idx, sample_order)]
}

# =============================================================================
# BUILD THE COMPOSITE FIGURE
# =============================================================================

# Individual barplot for one K value
make_admix_bar <- function(K_val, show_x_axis = FALSE) {
  dt <- admix_long[K == K_val]
  pal <- get_palette(K_val)

  p <- ggplot(dt, aes(x = plot_x, y = Q, fill = factor(component))) +
    geom_col(width = 1, color = NA) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
    labs(y = paste0("K=", K_val)) +
    theme_plate(grid = "none", panel_border = FALSE) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 7, angle = 0, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 5),
      plot.margin = margin(0, 2, 0, 2)
    )

  # Highlight best K
  if (K_val == best_K) {
    p <- p + theme(
      axis.title.y = element_text(size = 8, face = "bold", angle = 0,
                                  vjust = 0.5, hjust = 1, color = "#C44E52")
    )
  }
  p
}

# evalAdmix strip for one K
make_eval_strip <- function(K_val) {
  if (is.null(eval_burden)) return(NULL)
  dt <- eval_burden[K == K_val]
  if (nrow(dt) == 0) return(NULL)

  ggplot(dt, aes(x = plot_x, y = 1, fill = eval_burden)) +
    geom_tile(height = 1, color = NA) +
    scale_fill_gradientn(
      colors = c("#FFFFFF", "#FFF5EB", "#FDD49E", "#FC8D59", "#D73027"),
      limits = c(0, quantile(eval_burden$eval_burden, 0.98, na.rm = TRUE)),
      guide = "none"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))
}

# Relatedness strip
make_rel_strip <- function() {
  if (is.null(rel_strip)) return(NULL)

  ggplot(rel_strip, aes(x = plot_x, y = 1, fill = rel_class)) +
    geom_tile(height = 1, color = NA) +
    scale_fill_manual(values = pal_relatedness, name = "Relatedness") +
    scale_x_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(8, "pt"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7, face = "bold"),
      plot.margin = margin(0, 2, 0, 2)
    )
}

# Side panels
make_loglik_panel <- function() {
  dt <- best_dt[, .(K, loglik)]
  dt <- dt[!is.na(loglik)]

  ggplot(dt, aes(x = K, y = loglik)) +
    geom_line(linewidth = 0.6, color = pal_k20[1]) +
    geom_point(size = 1.5, color = pal_k20[1]) +
    geom_point(data = dt[K == best_K], size = 3, color = "#C44E52", shape = 18) +
    theme_plate_compact(grid = "y", panel_border = TRUE) +
    labs(x = "K", y = "Log-likelihood", title = "Model fit") +
    theme(plot.title = element_text(size = 8))
}

make_evaladmix_panel <- function() {
  dt <- best_dt[, .(K, mean_abs_resid)]
  dt <- dt[!is.na(mean_abs_resid)]
  if (nrow(dt) == 0) return(NULL)

  ggplot(dt, aes(x = K, y = mean_abs_resid)) +
    geom_line(linewidth = 0.6, color = pal_k20[3]) +
    geom_point(size = 1.5, color = pal_k20[3]) +
    geom_point(data = dt[K == best_K], size = 3, color = "#C44E52", shape = 18) +
    theme_plate_compact(grid = "y", panel_border = TRUE) +
    labs(x = "K", y = "Mean |residual|", title = "evalAdmix") +
    theme(plot.title = element_text(size = 8))
}

# =============================================================================
# ASSEMBLE
# =============================================================================

cat("[INFO] Building composite figure...\n")

# Main column: relatedness strip + (evalAdmix strip + admixture bar) for each K
main_panels <- list()
heights_main <- numeric()

# Relatedness strip
rel_p <- make_rel_strip()
if (!is.null(rel_p)) {
  main_panels[[length(main_panels) + 1]] <- rel_p
  heights_main <- c(heights_main, 0.3)
}

# K panels
for (K_val in K_values) {
  # evalAdmix strip
  eval_p <- make_eval_strip(K_val)
  if (!is.null(eval_p)) {
    main_panels[[length(main_panels) + 1]] <- eval_p
    heights_main <- c(heights_main, 0.12)
  }
  # Admixture barplot
  main_panels[[length(main_panels) + 1]] <- make_admix_bar(K_val)
  heights_main <- c(heights_main, 1)
}

# Stack main panels
main_col <- main_panels[[1]]
for (i in 2:length(main_panels)) {
  main_col <- main_col / main_panels[[i]]
}
main_col <- main_col + plot_layout(heights = heights_main)

# Side panels
side_panels <- list(make_loglik_panel())
eval_panel <- make_evaladmix_panel()
if (!is.null(eval_panel)) {
  side_panels[[length(side_panels) + 1]] <- eval_panel
}

side_col <- side_panels[[1]]
if (length(side_panels) > 1) {
  for (i in 2:length(side_panels)) {
    side_col <- side_col / side_panels[[i]]
  }
}

# Combine main + side
fig <- main_col | side_col
fig <- fig + plot_layout(widths = c(5, 1))

fig <- fig + plot_annotation(
  title = paste0("Admixture analysis: ", file_prefix),
  subtitle = paste0(n_samples, " samples, K=", min(K_values), "-", max(K_values),
                    ", best K=", best_K, " (highlighted)"),
  theme = theme(
    plot.title = element_text(size = SIZE_MAIN_TITLE, face = "bold",
                              color = TEXT_TITLE, hjust = 0),
    plot.subtitle = element_text(size = SIZE_SUBTITLE, color = TEXT_SUBTITLE, hjust = 0),
    plot.background = element_rect(fill = PLATE_BG, color = NA)
  )
)

# =============================================================================
# SAVE
# =============================================================================

total_height <- 2 + length(K_values) * opt[["height-per-k"]] + 1.5
out_pdf <- file.path(out_dir, paste0(file_prefix, "_admixture_K", min(K_values),
                                      "to", max(K_values), ".pdf"))
out_png <- file.path(out_dir, paste0(file_prefix, "_admixture_K", min(K_values),
                                      "to", max(K_values), ".png"))

ggsave(out_pdf, fig, width = opt$width, height = total_height,
       device = cairo_pdf, bg = PLATE_BG)
ggsave(out_png, fig, width = opt$width, height = total_height,
       dpi = 300, bg = PLATE_BG)

cat("[DONE] Saved:", out_pdf, "\n")
cat("[DONE] Saved:", out_png, "\n")
