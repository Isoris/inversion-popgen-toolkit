#!/usr/bin/env Rscript
# =============================================================================
# cli/STEP_U04_cluster_shape_profiles.R
# =============================================================================
# PCA + hierarchical (Ward.D2) + PAM (k=2..8 with silhouette) clustering on
# the candidate shape-score matrix. Writes:
#   shape_feature_matrix.tsv
#   shape_pca_coordinates.tsv
#   shape_cluster_assignments.tsv
#   shape_cluster_summary.tsv
#   shape_cluster_feature_means.tsv
#
# The atlas stores BOTH rule_based_class and unsupervised_cluster — clustering
# is descriptive, not authoritative.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(data.table)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--shape_scores",  type = "character"),
  make_option("--classes",       type = "character"),
  make_option("--out_dir",       type = "character"),
  make_option("--k_min",         type = "integer", default = 2L),
  make_option("--k_max",         type = "integer", default = 8L)
)))

scores  <- fread(opts$shape_scores)
classes <- fread(opts$classes)

feat_cols <- c("log10_length_bp",
               "dxy_inside_flank_ratio", "dxy_u_score",
               "dxy_internal_peak_score", "dxy_flatness_score",
               "abs_log2_asymmetry",
               "fst_inside_flank_ratio", "fst_u_score",
               "fst_internal_peak_score",
               "pi_imbalance_score",
               "arrangement_contrast_score",
               "oldness_score", "youngness_score",
               "neutral_u_shape_score",
               "local_adaptation_internal_score",
               "breakpoint_adaptation_score")
opt_z <- intersect(c("z_dxy_inside","z_fst_inside"), names(scores))
feat_cols <- intersect(c(feat_cols, opt_z), names(scores))

X <- as.matrix(scores[, ..feat_cols])
rownames(X) <- scores$candidate_id
# fill NAs with column medians (clustering is robust enough; we tag low-coverage)
for (c in seq_len(ncol(X))) {
  m <- median(X[, c], na.rm = TRUE)
  if (!is.finite(m)) m <- 0
  X[!is.finite(X[, c]), c] <- m
}
Xs <- scale(X)

fwrite(data.table(candidate_id = rownames(X), X),
       file.path(opts$out_dir, "shape_feature_matrix.tsv"), sep = "\t")

# PCA
pc <- prcomp(Xs, center = FALSE, scale. = FALSE)
pca_dt <- data.table(candidate_id = rownames(X),
                     PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     PC3 = if (ncol(pc$x) >= 3L) pc$x[, 3] else NA_real_)
fwrite(pca_dt, file.path(opts$out_dir, "shape_pca_coordinates.tsv"), sep = "\t")

# hierarchical (Ward.D2)
hc <- hclust(dist(Xs, method = "euclidean"), method = "ward.D2")

# PAM with silhouette
sil_avg <- function(k) {
  if (!requireNamespace("cluster", quietly = TRUE)) return(NA_real_)
  pm <- cluster::pam(Xs, k = k, diss = FALSE)
  s  <- cluster::silhouette(pm$clustering, dist(Xs))
  mean(s[, 3])
}
ks <- seq.int(opts$k_min, opts$k_max)
sils <- vapply(ks, sil_avg, numeric(1L))
best_k <- ks[which.max(sils)]
message("[U04] best k by silhouette = ", best_k,
        " (", round(max(sils, na.rm = TRUE), 3), ")")

if (requireNamespace("cluster", quietly = TRUE)) {
  pm <- cluster::pam(Xs, k = best_k, diss = FALSE)
  pam_clu <- pm$clustering
} else {
  pam_clu <- cutree(hc, k = best_k)
}
hc_clu <- cutree(hc, k = best_k)

# crude descriptive labels per cluster from feature centroids
cluster_means <- as.data.table(t(sapply(seq_len(best_k), function(k) {
  apply(X[pam_clu == k, , drop = FALSE], 2L, mean, na.rm = TRUE)
})))
setnames(cluster_means, feat_cols)
cluster_means[, cluster := seq_len(best_k)]
setcolorder(cluster_means, "cluster")

label_for <- function(row) {
  if (!is.na(row$dxy_u_score) && row$dxy_u_score > 1.5 &&
      !is.na(row$dxy_internal_peak_score) && row$dxy_internal_peak_score < 1.3)
    return("neutral_like_U_shape_cluster")
  if (!is.na(row$dxy_internal_peak_score) && row$dxy_internal_peak_score > 1.5)
    return("internal_peak_like_cluster")
  if (!is.na(row$dxy_inside_flank_ratio) && row$dxy_inside_flank_ratio > 1.5 &&
      !is.na(row$dxy_flatness_score) && row$dxy_flatness_score < 0.4)
    return("flat_deep_cluster")
  if (!is.na(row$abs_log2_asymmetry) && row$abs_log2_asymmetry > 1.0)
    return("asymmetric_edge_cluster")
  if (!is.na(row$oldness_score) && row$oldness_score < 1.0)
    return("young_weak_cluster")
  return("complex_mixed_cluster")
}
cluster_means[, cluster_label := vapply(seq_len(.N), function(i) label_for(.SD[i]),
                                        character(1L))]

assign <- data.table(candidate_id = rownames(X),
                     hc_cluster = hc_clu, pam_cluster = pam_clu)
assign <- merge(assign,
                cluster_means[, .(cluster, cluster_label)],
                by.x = "pam_cluster", by.y = "cluster", all.x = TRUE)
assign <- merge(assign, classes[, .(candidate_id, primary_class, confidence)],
                by = "candidate_id", all.x = TRUE)

fwrite(assign, file.path(opts$out_dir, "shape_cluster_assignments.tsv"), sep = "\t")
fwrite(cluster_means, file.path(opts$out_dir, "shape_cluster_feature_means.tsv"), sep = "\t")
fwrite(data.table(k = ks, silhouette_avg = sils,
                  best = (ks == best_k)),
       file.path(opts$out_dir, "shape_cluster_summary.tsv"), sep = "\t")
message("[U04] wrote cluster assignments for ", nrow(assign), " candidates")
