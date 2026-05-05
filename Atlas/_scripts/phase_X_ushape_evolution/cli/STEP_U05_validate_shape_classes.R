#!/usr/bin/env Rscript
# =============================================================================
# cli/STEP_U05_validate_shape_classes.R
# =============================================================================
# Tests whether rule-based primary_class and unsupervised PAM clusters separate
# in multivariate score space, plus a length-class control. Skips gracefully
# when vegan is missing.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(data.table)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--feature_matrix",   type = "character"),
  make_option("--class_table",      type = "character"),
  make_option("--cluster_table",    type = "character"),
  make_option("--out_dir",          type = "character"),
  make_option("--n_perm",           type = "integer", default = 999L)
)))

X <- fread(opts$feature_matrix)
classes <- fread(opts$class_table)
clusters <- fread(opts$cluster_table)

mat <- as.matrix(X[, !"candidate_id", with = FALSE])
rownames(mat) <- X$candidate_id
mat <- scale(mat)
mat[!is.finite(mat)] <- 0

meta <- merge(merge(data.table(candidate_id = rownames(mat)),
                    classes[, .(candidate_id, primary_class)], by = "candidate_id"),
              clusters[, .(candidate_id, pam_cluster, cluster_label)],
              by = "candidate_id", all.x = TRUE)

if (!requireNamespace("vegan", quietly = TRUE)) {
  message("[U05] vegan not installed; skipping PERMANOVA/ANOSIM")
  writeLines("vegan not installed",
             file.path(opts$out_dir, "validation_summary.txt"))
  quit(status = 0)
}

D <- vegan::vegdist(mat, method = "euclidean")
permanova <- list(); anosim <- list(); betadisper <- list()

run_one <- function(label, group_vec) {
  if (length(unique(group_vec)) < 2L) return(invisible(NULL))
  pad <- vegan::adonis2(D ~ group_vec, permutations = opts$n_perm)
  ano <- vegan::anosim(D, grouping = group_vec, permutations = opts$n_perm)
  bd  <- vegan::betadisper(D, group = group_vec)
  bd_p <- vegan::permutest(bd, permutations = opts$n_perm)$tab[1, "Pr(>F)"]
  permanova[[label]] <<- data.table(
    factor = label,
    R2 = pad$R2[1], F = pad$F[1], P = pad$`Pr(>F)`[1])
  anosim[[label]] <<- data.table(
    factor = label, R = ano$statistic, P = ano$signif)
  betadisper[[label]] <<- data.table(
    factor = label, dispersion_P = bd_p)
}

run_one("rule_based_primary_class", meta$primary_class)
run_one("pam_cluster",              as.factor(meta$pam_cluster))
# length class — derive from feature matrix log10_length_bp if present
if ("log10_length_bp" %in% colnames(mat)) {
  L <- 10^X$log10_length_bp
  lcl <- ifelse(L < 100e3, "small",
                ifelse(L < 1e6, "medium",
                       ifelse(L < 1e7, "large", "mega")))
  run_one("length_class", lcl)
}

if (length(permanova))
  fwrite(rbindlist(permanova), file.path(opts$out_dir, "permanova_results.tsv"), sep = "\t")
if (length(anosim))
  fwrite(rbindlist(anosim),    file.path(opts$out_dir, "anosim_results.tsv"),    sep = "\t")
if (length(betadisper))
  fwrite(rbindlist(betadisper),file.path(opts$out_dir, "betadisper_results.tsv"),sep = "\t")

writeLines(c(
  "PERMANOVA tests whether classes differ in multivariate space.",
  "ANOSIM tests whether between-class > within-class distances.",
  "betadisper checks whether significance is dispersion-driven (low dispersion P = caution).",
  paste0("n_perm=", opts$n_perm)
), file.path(opts$out_dir, "validation_summary.txt"))

message("[U05] validation tests written")
