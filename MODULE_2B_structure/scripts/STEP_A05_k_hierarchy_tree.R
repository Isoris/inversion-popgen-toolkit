#!/usr/bin/env Rscript
###############################################################################
# compute_k_hierarchy_tree.R — Founder lineage merge tree from NGSadmix
#
# For consecutive K values (K, K-1), matches components by maximum Q-overlap
# to determine which components merge. This reveals the nested founder
# lineage structure in the hatchery population.
#
# Inputs:
#   --results-dir  Directory containing *_best_seed_by_K.tsv and *_best.qopt
#   --sample-file  Canonical sample order file
#   --file-prefix  e.g. "wholegenome_thin500_all226"
#   --out-dir      Output directory
#
# Outputs:
#   - {prefix}_merge_events.tsv       : K, parent, child_a, child_b, n_parent, ...
#   - {prefix}_merge_tree.tsv         : full tree as edgelist for visualization
#   - {prefix}_component_tracking.tsv : per-sample component assignment across all K
#
# Usage:
#   Rscript compute_k_hierarchy_tree.R \
#     --results-dir structure_results/wholegenome_thin500_all226/ \
#     --sample-file samples.txt \
#     --file-prefix wholegenome_thin500_all226 \
#     --out-dir structure_results/wholegenome_thin500_all226/
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--results-dir", type = "character", help = "Dir with best qopt files"),
  make_option("--sample-file", type = "character", help = "Canonical sample list"),
  make_option("--file-prefix", type = "character", help = "Prefix for file names"),
  make_option("--out-dir", type = "character", help = "Output directory"),
  make_option("--k-min", type = "integer", default = 2),
  make_option("--k-max", type = "integer", default = 20)
)
opt <- parse_args(OptionParser(option_list = option_list))

results_dir <- opt[["results-dir"]]
sample_file <- opt[["sample-file"]]
file_prefix <- opt[["file-prefix"]]
out_dir     <- opt[["out-dir"]]
K_min       <- opt[["k-min"]]
K_max       <- opt[["k-max"]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD SAMPLE LIST
# =============================================================================

samples <- trimws(readLines(sample_file, warn = FALSE))
samples <- samples[nzchar(samples)]
n_samples <- length(samples)
cat("[INFO] Samples:", n_samples, "\n")

# =============================================================================
# LOAD ALL BEST QOPT MATRICES
# =============================================================================

read_qopt <- function(path, K) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  qmat <- do.call(rbind, lapply(strsplit(lines, "[[:space:]]+"), as.numeric))
  stopifnot(nrow(qmat) == n_samples, ncol(qmat) == K)
  qmat
}

# Find all available best qopt files
qopt_list <- list()
for (K in seq(K_min, K_max)) {
  stem <- sprintf("%s_K%02d_best", file_prefix, K)
  path <- file.path(results_dir, paste0(stem, ".qopt"))
  if (file.exists(path)) {
    qopt_list[[as.character(K)]] <- read_qopt(path, K)
    cat("[INFO] Loaded K=", K, " (", path, ")\n")
  } else {
    cat("[WARN] Missing K=", K, ": ", path, "\n")
  }
}

K_available <- sort(as.integer(names(qopt_list)))
cat("[INFO] Available K values:", paste(K_available, collapse = ","), "\n")

if (length(K_available) < 2) stop("Need at least 2 K values for merge tree")

# =============================================================================
# ASSIGN DOMINANT COMPONENT PER SAMPLE AT EACH K
# =============================================================================

# For each K, assign each sample to its dominant component
assignments <- data.table(
  sample_index = rep(seq_len(n_samples), length(K_available)),
  sample = rep(samples, length(K_available)),
  K = rep(K_available, each = n_samples)
)

assignments[, component := {
  K_val <- K[1]
  qmat <- qopt_list[[as.character(K_val)]]
  max.col(qmat, ties.method = "first")
}, by = K]

assignments[, component_label := paste0("K", K, "_Q", component)]

# Write component tracking
out_tracking <- file.path(out_dir, paste0(file_prefix, "_component_tracking.tsv"))
# Wide format: one row per sample, columns for each K
tracking_wide <- dcast(assignments, sample_index + sample ~ K,
                       value.var = "component",
                       fun.aggregate = function(x) x[1])
setnames(tracking_wide,
         as.character(K_available),
         paste0("K", K_available, "_component"))
fwrite(tracking_wide, out_tracking, sep = "\t", quote = FALSE)
cat("[INFO] Wrote:", out_tracking, "\n")

# =============================================================================
# COMPUTE MERGE EVENTS
# =============================================================================

merge_events <- list()
tree_edges <- list()

for (i in seq_along(K_available)[-1]) {
  K_hi <- K_available[i]
  K_lo <- K_available[i - 1]

  if (K_hi - K_lo != 1) {
    cat("[WARN] Gap between K=", K_lo, "and K=", K_hi, "- skipping\n")
    next
  }

  qmat_hi <- qopt_list[[as.character(K_hi)]]
  qmat_lo <- qopt_list[[as.character(K_lo)]]

  # Assignment at K_hi
  assign_hi <- max.col(qmat_hi, ties.method = "first")
  # Assignment at K_lo
  assign_lo <- max.col(qmat_lo, ties.method = "first")

  # For each component at K_hi, find its best-matching component at K_lo
  # by plurality vote of its member samples
  mapping <- data.table(
    comp_hi = assign_hi,
    comp_lo = assign_lo
  )

  # For each hi component, which lo component has the most members?
  hi_to_lo <- mapping[, .(
    best_lo = {
      tab <- table(comp_lo)
      as.integer(names(tab)[which.max(tab)])
    },
    n_samples_hi = .N
  ), by = comp_hi]

  # Two hi-components mapping to the same lo-component = a merge event
  merge_groups <- hi_to_lo[, .(
    children = list(comp_hi),
    n_children = .N,
    total_samples = sum(n_samples_hi)
  ), by = best_lo]

  for (j in seq_len(nrow(merge_groups))) {
    parent_lo <- merge_groups$best_lo[j]
    children_hi <- merge_groups$children[[j]]
    n_children <- merge_groups$n_children[j]

    # Record merge event (only real merges where n_children > 1)
    if (n_children > 1) {
      merge_events[[length(merge_events) + 1]] <- data.table(
        K_transition = paste0(K_hi, "->", K_lo),
        K_from = K_hi,
        K_to = K_lo,
        parent_component = parent_lo,
        parent_label = paste0("K", K_lo, "_Q", parent_lo),
        child_components = paste(children_hi, collapse = ","),
        child_labels = paste0("K", K_hi, "_Q", paste(children_hi, collapse = paste0(",K", K_hi, "_Q"))),
        n_children = n_children,
        total_samples = merge_groups$total_samples[j]
      )
    }

    # Record all tree edges (including 1:1 mappings for the alluvial plot)
    for (ch in children_hi) {
      n_ch <- hi_to_lo[comp_hi == ch, n_samples_hi]
      tree_edges[[length(tree_edges) + 1]] <- data.table(
        K_from = K_hi,
        K_to = K_lo,
        comp_from = ch,
        comp_to = parent_lo,
        label_from = paste0("K", K_hi, "_Q", ch),
        label_to = paste0("K", K_lo, "_Q", parent_lo),
        n_samples = n_ch,
        is_merge = (n_children > 1)
      )
    }
  }
}

# Write merge events
merge_dt <- rbindlist(merge_events, use.names = TRUE, fill = TRUE)
out_merges <- file.path(out_dir, paste0(file_prefix, "_merge_events.tsv"))
fwrite(merge_dt, out_merges, sep = "\t", quote = FALSE)
cat("[INFO] Wrote:", out_merges, " (", nrow(merge_dt), " merge events)\n")

# Write tree edges (for Sankey/alluvial visualization)
tree_dt <- rbindlist(tree_edges, use.names = TRUE, fill = TRUE)
out_tree <- file.path(out_dir, paste0(file_prefix, "_merge_tree.tsv"))
fwrite(tree_dt, out_tree, sep = "\t", quote = FALSE)
cat("[INFO] Wrote:", out_tree, " (", nrow(tree_dt), " edges)\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Merge tree summary ===\n")
if (nrow(merge_dt) > 0) {
  cat("Total merge events:", nrow(merge_dt), "\n")
  cat("Deepest merge (low K):", min(merge_dt$K_to), "\n")
  cat("Shallowest merge (high K):", max(merge_dt$K_from), "\n\n")
  cat("Merge events by K transition:\n")
  print(merge_dt[, .N, by = K_transition])
} else {
  cat("No merge events found (all 1:1 mappings)\n")
}

cat("\n[DONE] K-hierarchy merge tree computation\n")
