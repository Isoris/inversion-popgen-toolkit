#!/usr/bin/env Rscript

# =============================================================================
# STEP40_architecture_classifier.R  v1
#
# Candidate Inversion Architecture Classifier
# --------------------------------------------
# Classifies each candidate dosage heatmap into a structural class:
#
#   Level 1 — architecture class
#     SIMPLE_STRONG       : clean 3-state inversion, high contrast
#     SIMPLE_WEAK         : 3-state inversion, low contrast / noisy
#     COMPOSITE_INTERNAL  : one main axis + internal marker-family substructure
#     COMPOSITE_OVERLAP   : two (or more) distinct marker-boundary systems
#     DEMOGRAPHIC_BLOCK   : local LD block driven by kinship / broodline
#     TECHNICAL           : artifact-prone (gaps, repeats, depth anomalies)
#     UNRESOLVED          : mixed evidence, cannot be separated confidently
#
#   Level 2 — confidence tier: A (strong), B, C, D (weak)
#
#   Final label e.g.: "COMPOSITE_INTERNAL-B"
#
# CORE PRINCIPLE — "separate first, correct second":
#   For composite candidates, we first partition samples using discordant
#   vertical marker patterns that encode hidden subgroup boundaries, and
#   only then apply subgroup-specific polarity harmonization. A single
#   global polarity model is wrong when the candidate contains overlapping
#   inversions, nested structure, founder haplotypes, or mixed backgrounds.
#
# Features computed per candidate:
#   F1  simple_axis_score        : variance explained by sample PC1
#   F2  boundary_consistency     : agreement of marker-level sample partitions
#   F3  marker_family_k          : # coherent marker families (silhouette-selected)
#   F4  subinterval_discordance  : left/middle/right partition agreement
#   F5  within_class_substructure: mean PC1 variance inside each karyotype
#   F6  hwe_fit                  : χ² fit of 3-band counts to HWE at inferred p
#   F7  monotonic_fraction       : fraction of markers with monotonic group means
#   F8  polarity_block_count     : # L2 polarity blocks
#   F9  technical_score          : gaps / depth / SNP-density anomaly flag
#   F10 family_ld_correlation    : corr of sample partition with kinship PCs
#
# Usage:
#   Rscript STEP40_architecture_classifier.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

if (file.exists(config_file)) source(config_file)

# ── Defaults if config missing values ──────────────────────────────────────
if (!exists("FOLLOWUP_DIR"))     FOLLOWUP_DIR     <- "results/followup"
if (!exists("DOSAGE_DIR"))       DOSAGE_DIR       <- "results/dosage"
if (!exists("CANDIDATE_TABLE"))  CANDIDATE_TABLE  <- "results/candidates.tsv"
if (!exists("PLOTS_DIR"))        PLOTS_DIR        <- "results/plots"
if (!exists("ensure_dir"))       ensure_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); invisible(p) }

# ── Parameters ─────────────────────────────────────────────────────────────
PARAMS <- list(
  min_markers_for_classification = 30L,   # skip if fewer
  min_samples_per_group           = 5L,   # skip a group if smaller
  top_marker_variance_fraction    = 0.5,  # top 50% variance markers for features
  marker_family_max_k             = 4L,   # try k=1..4 marker families
  subinterval_n_bins              = 3L,   # left / middle / right
  min_block_markers               = 5L,
  # Feature-to-class thresholds (tuned conservatively, override in config)
  thr_simple_axis_strong          = 0.35,
  thr_simple_axis_weak            = 0.15,
  thr_boundary_consistency_strong = 0.70,
  thr_boundary_consistency_weak   = 0.45,
  thr_monotonic_strong            = 0.75,
  thr_monotonic_weak              = 0.50,
  thr_within_substructure_high    = 0.30,
  thr_subinterval_discord_high    = 0.25,
  thr_technical_high              = 0.50,
  thr_hwe_p_ok                    = 0.01,
  thr_family_ld_corr_high         = 0.60
)
# Let config override
if (exists("CLASSIFIER_PARAMS") && is.list(CLASSIFIER_PARAMS)) {
  for (k in names(CLASSIFIER_PARAMS)) PARAMS[[k]] <- CLASSIFIER_PARAMS[[k]]
}

# =============================================================================
# UTILITIES
# =============================================================================

safe_var <- function(x) {
  v <- var(x, na.rm = TRUE)
  if (!is.finite(v)) 0 else v
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (safe_var(x[ok]) == 0 || safe_var(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok]))
}

# Proportion of variance explained by first PC of a matrix (samples × markers).
# Works with NA via mean imputation per marker (minimal bias, robust).
pc1_var_explained <- function(M) {
  if (is.null(M) || !is.matrix(M)) return(NA_real_)
  if (nrow(M) < 3 || ncol(M) < 2) return(NA_real_)
  # Column-wise mean imputation
  col_means <- colMeans(M, na.rm = TRUE)
  col_means[!is.finite(col_means)] <- 0
  M_imp <- M
  for (j in seq_len(ncol(M_imp))) {
    bad <- !is.finite(M_imp[, j])
    if (any(bad)) M_imp[bad, j] <- col_means[j]
  }
  # Drop zero-variance columns
  v <- apply(M_imp, 2, safe_var)
  keep <- which(v > 1e-10)
  if (length(keep) < 2) return(NA_real_)
  M_imp <- M_imp[, keep, drop = FALSE]
  pc <- tryCatch(
    prcomp(M_imp, center = TRUE, scale. = FALSE, rank. = 3),
    error = function(e) NULL
  )
  if (is.null(pc)) return(NA_real_)
  sd2 <- pc$sdev^2
  total <- sum(v)
  if (total <= 0) return(NA_real_)
  as.numeric(sd2[1] / total)
}

# Monotonicity fraction: for each marker, is the group-mean order monotonic
# across HOMO_1 → HET → HOMO_2 (either 0→1→2 or 2→1→0)?
monotonic_fraction <- function(grp_means) {
  # grp_means: markers × 3 with cols HOMO_1, HET, HOMO_2
  h1 <- grp_means[, "HOMO_1"]; he <- grp_means[, "HET"]; h2 <- grp_means[, "HOMO_2"]
  ok <- is.finite(h1) & is.finite(he) & is.finite(h2)
  if (sum(ok) < 5) return(NA_real_)
  up   <- (h1 <= he) & (he <= h2) & (h2 - h1 > 0.05)
  down <- (h1 >= he) & (he >= h2) & (h1 - h2 > 0.05)
  mean((up | down)[ok])
}

# Within-class substructure: PC1 variance explained within each karyotype group.
# High values mean the karyotype is internally heterogeneous.
within_class_substructure <- function(X_samples_markers, groups, min_n = 5) {
  # X_samples_markers: samples × markers (rows = samples)
  grp_levels <- c("HOMO_1", "HET", "HOMO_2")
  scores <- c()
  for (g in grp_levels) {
    idx <- which(groups == g)
    if (length(idx) >= min_n) {
      ve <- pc1_var_explained(X_samples_markers[idx, , drop = FALSE])
      if (is.finite(ve)) scores <- c(scores, ve)
    }
  }
  if (length(scores) == 0) return(NA_real_)
  mean(scores)
}

# Cluster markers into families by sample-dosage pattern.
# Uses correlation distance + hierarchical clustering with silhouette selection.
cluster_marker_families <- function(X_markers_samples, max_k = 4L, top_n = 200L) {
  # X_markers_samples: markers × samples
  if (nrow(X_markers_samples) < 10) return(list(k = 1L, cluster = rep(1L, nrow(X_markers_samples)), silhouette = NA_real_))

  # Take the most informative (highest-variance) markers for stability
  mv <- apply(X_markers_samples, 1, safe_var)
  keep <- order(mv, decreasing = TRUE)[seq_len(min(top_n, nrow(X_markers_samples)))]
  M <- X_markers_samples[keep, , drop = FALSE]

  # Mean-impute missing
  for (j in seq_len(ncol(M))) {
    bad <- !is.finite(M[, j])
    if (any(bad)) M[bad, j] <- mean(M[, j], na.rm = TRUE)
  }
  M[!is.finite(M)] <- 0

  # Signed correlation → we want patterns that agree OR are polarity-flipped
  # to be grouped together. So we cluster on |cor|-based distance of dosage centered rows.
  # But we also need to detect overlapping systems that have different partitions —
  # those should NOT cluster together even if correlated in magnitude.
  # So: use raw Pearson correlation (not absolute), and flip rows to a common polarity
  # only if cor sign is negative with row 1, then recompute. This way, true "opposite
  # polarity of same partition" collapses, and different partitions stay apart.
  ref_row <- M[1, ]
  for (i in seq_len(nrow(M))) {
    cr <- safe_cor(M[i, ], ref_row)
    if (is.finite(cr) && cr < 0) M[i, ] <- -M[i, ] + 2  # polarity-flip (dosage 0..2 axis)
  }

  # Now correlation distance
  cr_mat <- suppressWarnings(cor(t(M), use = "pairwise.complete.obs"))
  cr_mat[!is.finite(cr_mat)] <- 0
  d <- as.dist(1 - cr_mat)

  hc <- tryCatch(hclust(d, method = "average"), error = function(e) NULL)
  if (is.null(hc)) return(list(k = 1L, cluster = rep(1L, nrow(X_markers_samples)), silhouette = NA_real_))

  # Silhouette across k = 1..max_k
  best_k <- 1L; best_sil <- -Inf; best_cl <- rep(1L, length(keep))
  for (k in 2L:max_k) {
    cl_k <- cutree(hc, k = k)
    # Silhouette — compute only if enough members
    if (min(table(cl_k)) < 3) next
    sil <- tryCatch({
      # lightweight silhouette without extra package dependency
      # average inter-cluster - intra-cluster over sd
      cl_unique <- unique(cl_k)
      a_i <- numeric(length(cl_k))
      b_i <- numeric(length(cl_k))
      D <- as.matrix(d)
      for (i in seq_along(cl_k)) {
        same <- which(cl_k == cl_k[i] & seq_along(cl_k) != i)
        a_i[i] <- if (length(same)) mean(D[i, same]) else 0
        other_bs <- c()
        for (cc in setdiff(cl_unique, cl_k[i])) {
          oth <- which(cl_k == cc)
          if (length(oth)) other_bs <- c(other_bs, mean(D[i, oth]))
        }
        b_i[i] <- if (length(other_bs)) min(other_bs) else 0
      }
      s <- (b_i - a_i) / pmax(a_i, b_i, 1e-10)
      mean(s, na.rm = TRUE)
    }, error = function(e) -Inf)
    if (is.finite(sil) && sil > best_sil) {
      best_sil <- sil; best_k <- k; best_cl <- cl_k
    }
  }

  # Silhouette threshold for accepting k > 1
  if (best_sil < 0.15) {
    best_k <- 1L
    best_cl <- rep(1L, length(keep))
  }

  # Broadcast cluster labels back to all markers: unclustered markers get label 1
  full_cluster <- rep(1L, nrow(X_markers_samples))
  full_cluster[keep] <- best_cl

  list(k = best_k, cluster = full_cluster, silhouette = best_sil,
       keep_idx = keep)
}

# Boundary consistency: do markers agree on how samples are partitioned?
# For each informative marker, derive a sample ordering (by dosage rank).
# Compute pairwise rank correlations between markers. High mean = consistent.
boundary_consistency <- function(X_markers_samples, n_markers_use = 50L) {
  if (nrow(X_markers_samples) < 5 || ncol(X_markers_samples) < 5) return(NA_real_)
  mv <- apply(X_markers_samples, 1, safe_var)
  top <- order(mv, decreasing = TRUE)[seq_len(min(n_markers_use, nrow(X_markers_samples)))]
  M <- X_markers_samples[top, , drop = FALSE]

  # Polarity-flip negative-correlated markers to the first one
  ref <- M[1, ]
  for (i in 2:nrow(M)) {
    cr <- safe_cor(M[i, ], ref)
    if (is.finite(cr) && cr < 0) M[i, ] <- -M[i, ] + 2
  }
  # Rank-correlation matrix across markers (which samples get high vs low)
  R <- tryCatch(suppressWarnings(cor(t(M), method = "spearman", use = "pairwise.complete.obs")),
                error = function(e) NULL)
  if (is.null(R)) return(NA_real_)
  R[!is.finite(R)] <- 0
  diag(R) <- NA
  mean(abs(R), na.rm = TRUE)  # after polarity flip, high = consistent partition
}

# Subinterval discordance: split interval into n_bins, for each bin compute
# a sample PC1 score, then compute mean |1 - cor| between bin scores.
# High = overlapping systems (left bin groups samples differently than right).
subinterval_discordance <- function(X_markers_samples, positions, n_bins = 3L, min_m = 10L) {
  if (nrow(X_markers_samples) < n_bins * min_m) return(NA_real_)
  # Equal-width bins in position space
  rng <- range(positions, na.rm = TRUE)
  if (diff(rng) <= 0) return(NA_real_)
  breaks <- seq(rng[1], rng[2], length.out = n_bins + 1)
  bin_id <- cut(positions, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  scores <- matrix(NA_real_, nrow = ncol(X_markers_samples), ncol = n_bins)
  for (b in seq_len(n_bins)) {
    midx <- which(bin_id == b)
    if (length(midx) < min_m) next
    sub <- X_markers_samples[midx, , drop = FALSE]
    for (j in seq_len(ncol(sub))) {
      bad <- !is.finite(sub[, j])
      if (any(bad)) sub[bad, j] <- mean(sub[, j], na.rm = TRUE)
    }
    sub[!is.finite(sub)] <- 0
    v <- apply(sub, 1, safe_var)
    kk <- which(v > 1e-10)
    if (length(kk) < 2) next
    pc <- tryCatch(prcomp(t(sub[kk, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1),
                   error = function(e) NULL)
    if (!is.null(pc)) scores[, b] <- pc$x[, 1]
  }
  # Correlate bins pairwise
  ok_bins <- which(colSums(is.finite(scores)) >= 5)
  if (length(ok_bins) < 2) return(NA_real_)
  # Flip sign to common polarity (bin 1 as reference)
  ref_b <- ok_bins[1]
  for (b in ok_bins[-1]) {
    cr <- safe_cor(scores[, b], scores[, ref_b])
    if (is.finite(cr) && cr < 0) scores[, b] <- -scores[, b]
  }
  cors <- c()
  for (i in seq_along(ok_bins)) for (j in seq_along(ok_bins)) {
    if (i < j) {
      cr <- safe_cor(scores[, ok_bins[i]], scores[, ok_bins[j]])
      if (is.finite(cr)) cors <- c(cors, cr)
    }
  }
  if (length(cors) == 0) return(NA_real_)
  1 - mean(cors)  # 0 = perfectly consistent across bins, >0 = discordant
}

# Exact HWE χ² goodness-of-fit given 3 observed counts, allele freq = MLE
hwe_gof <- function(n_h1, n_het, n_h2) {
  n <- n_h1 + n_het + n_h2
  if (n < 10) return(list(p_value = NA_real_, chisq = NA_real_, p_allele = NA_real_))
  p <- (2 * n_h1 + n_het) / (2 * n)
  q <- 1 - p
  e_h1 <- n * p^2; e_het <- 2 * n * p * q; e_h2 <- n * q^2
  # Avoid division by zero
  e <- c(e_h1, e_het, e_h2); o <- c(n_h1, n_het, n_h2)
  keep <- e > 0
  if (sum(keep) < 2) return(list(p_value = NA_real_, chisq = NA_real_, p_allele = p))
  chi <- sum((o[keep] - e[keep])^2 / e[keep])
  pv <- 1 - pchisq(chi, df = 1)  # 1 df because p is estimated
  list(p_value = pv, chisq = chi, p_allele = p)
}

# Technical-artifact score: combines N-content, SNP-density anomaly,
# and number of low-SNP-density windows in the interval.
technical_score <- function(sites_pos, interval_start, interval_end,
                             gap_file = NULL, density_bin_bp = 50000L) {
  span <- interval_end - interval_start
  if (span <= 0) return(0)
  n_bins <- max(1L, as.integer(ceiling(span / density_bin_bp)))
  breaks <- seq(interval_start, interval_end, length.out = n_bins + 1)
  bin_ids <- cut(sites_pos, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bin_counts <- tabulate(bin_ids, nbins = n_bins)
  # Fraction of bins with 0 SNPs = dropout fraction
  dropout_frac <- mean(bin_counts == 0)
  # Coefficient of variation of nonzero bins (high CV = patchy)
  nz <- bin_counts[bin_counts > 0]
  cv <- if (length(nz) > 2) sd(nz) / mean(nz) else 0
  cv_norm <- min(cv / 2, 1)  # normalize roughly to [0,1]
  # Gap file: optional TSV with chrom,start,end of assembly gaps
  gap_frac <- 0
  if (!is.null(gap_file) && file.exists(gap_file)) {
    g <- tryCatch(fread(gap_file), error = function(e) NULL)
    if (!is.null(g) && nrow(g) > 0 && all(c("start", "end") %in% names(g))) {
      overlap_bp <- 0
      for (r in seq_len(nrow(g))) {
        o_lo <- max(interval_start, as.numeric(g$start[r]))
        o_hi <- min(interval_end, as.numeric(g$end[r]))
        if (o_hi > o_lo) overlap_bp <- overlap_bp + (o_hi - o_lo)
      }
      gap_frac <- overlap_bp / span
    }
  }
  # Aggregate
  score <- 0.5 * dropout_frac + 0.3 * cv_norm + 0.2 * gap_frac
  min(score, 1)
}

# Correlation of the candidate-level sample partition with genome-wide kinship PCs.
# High correlation means the partition looks like family/broodline structure
# rather than a real inversion signal.
family_ld_correlation <- function(sample_scores, kinship_pcs) {
  # sample_scores: named numeric vector (sample → score)
  # kinship_pcs:   data.table with columns sample, PC1, PC2, PC3, ...
  if (is.null(kinship_pcs) || !"sample" %in% names(kinship_pcs)) return(NA_real_)
  common <- intersect(names(sample_scores), kinship_pcs$sample)
  if (length(common) < 10) return(NA_real_)
  s <- sample_scores[common]
  pc_cols <- grep("^PC[0-9]+$", names(kinship_pcs), value = TRUE)
  if (length(pc_cols) == 0) return(NA_real_)
  kp <- kinship_pcs[match(common, kinship_pcs$sample), ..pc_cols]
  max_abs <- 0
  for (c in pc_cols) {
    cr <- safe_cor(s, kp[[c]])
    if (is.finite(cr) && abs(cr) > max_abs) max_abs <- abs(cr)
  }
  max_abs
}

# =============================================================================
# CORE: separate-first sample partitioning
# =============================================================================
#
# For a candidate with composite/discordant markers, partition samples into
# coherent backgrounds using the marker-family structure.
#
#   1. Cluster markers into K families (via cluster_marker_families)
#   2. For each family, compute a per-sample score (family PC1)
#   3. Stack family scores → samples × K score matrix
#   4. Cluster samples on this joint score matrix (k-means on 1..max_K+1)
#   5. Return: sample subgroup assignments + per-family sample scores
#
# Within each sample subgroup, polarity harmonization is done separately
# (delegated to the calling code via the returned subgroup labels).
#
# =============================================================================
separate_first_partition <- function(X_markers_samples, groups, max_subgroups = 4L) {
  n_markers <- nrow(X_markers_samples); n_samples <- ncol(X_markers_samples)

  # Step 1: marker families
  fam <- cluster_marker_families(X_markers_samples, max_k = PARAMS$marker_family_max_k)

  # Step 2 & 3: per-family sample scores
  fam_scores <- matrix(NA_real_, nrow = n_samples, ncol = fam$k)
  colnames(fam_scores) <- paste0("family_", seq_len(fam$k))
  for (fi in seq_len(fam$k)) {
    mi <- which(fam$cluster == fi)
    if (length(mi) < 5) next
    sub <- X_markers_samples[mi, , drop = FALSE]
    for (j in seq_len(ncol(sub))) {
      bad <- !is.finite(sub[, j])
      if (any(bad)) sub[bad, j] <- mean(sub[, j], na.rm = TRUE)
    }
    sub[!is.finite(sub)] <- 0
    v <- apply(sub, 1, safe_var)
    kk <- which(v > 1e-10)
    if (length(kk) < 2) next
    pc <- tryCatch(prcomp(t(sub[kk, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1),
                   error = function(e) NULL)
    if (!is.null(pc)) fam_scores[, fi] <- pc$x[, 1]
  }

  # Step 4: if only one family, subgroup = global PC1 sign (HOMO_1 / HET / HOMO_2 already).
  # If multiple families, cluster samples on the joint family-score space.
  subgroup <- rep(1L, n_samples)
  if (fam$k >= 2L) {
    fs <- fam_scores
    # impute any NA columns with 0, drop all-NA columns
    col_ok <- apply(fs, 2, function(z) sum(is.finite(z)) > 5)
    if (sum(col_ok) >= 2) {
      fs <- fs[, col_ok, drop = FALSE]
      for (j in seq_len(ncol(fs))) {
        m <- mean(fs[, j], na.rm = TRUE); s <- sd(fs[, j], na.rm = TRUE)
        if (!is.finite(s) || s == 0) s <- 1
        fs[, j] <- (fs[, j] - m) / s
        fs[!is.finite(fs[, j]), j] <- 0
      }
      # k-means with k=2..max_subgroups, pick by average silhouette-lite
      best_k <- 1L; best_score <- -Inf; best_labels <- rep(1L, nrow(fs))
      for (k in 2L:min(max_subgroups, nrow(fs) %/% 5)) {
        km <- tryCatch(kmeans(fs, centers = k, nstart = 10, iter.max = 50),
                       error = function(e) NULL)
        if (is.null(km)) next
        # Score = between_SS / total_SS (variance ratio)
        r <- km$betweenss / km$totss
        # Penalty for very unequal cluster sizes
        size_ok <- min(km$size) >= 5
        if (size_ok && r > best_score) { best_score <- r; best_k <- k; best_labels <- km$cluster }
      }
      if (best_k >= 2L && best_score > 0.4) subgroup <- as.integer(best_labels)
    }
  }

  list(
    marker_families = fam,
    family_scores = fam_scores,
    subgroup = subgroup,
    n_subgroups = length(unique(subgroup))
  )
}

# =============================================================================
# FEATURE EXTRACTION per candidate
# =============================================================================
compute_features <- function(X, sites_reg, rot, groups, sample_cols,
                              interval_start, interval_end,
                              gap_file = NULL, kinship_pcs = NULL) {
  # X: markers × samples, column names = sample IDs, row positions in sites_reg$pos
  # Reduce to informative markers for most features
  mv <- apply(X, 1, safe_var)
  var_thr <- quantile(mv, 1 - PARAMS$top_marker_variance_fraction, na.rm = TRUE)
  info_idx <- which(mv >= var_thr)
  if (length(info_idx) < 10) info_idx <- order(mv, decreasing = TRUE)[seq_len(min(30, nrow(X)))]
  Xi <- X[info_idx, , drop = FALSE]
  pos_i <- sites_reg$pos[info_idx]

  # ── F1: simple axis score (sample-space PC1 VE on informative markers) ──
  X_samp_mark <- t(Xi)  # samples × markers
  f1_simple_axis <- pc1_var_explained(X_samp_mark)

  # ── F2: boundary consistency (Spearman |cor| across markers post-flip) ──
  f2_boundary <- boundary_consistency(Xi, n_markers_use = min(100L, nrow(Xi)))

  # ── F3: marker family k (number of coherent families) ──
  fam <- cluster_marker_families(Xi, max_k = PARAMS$marker_family_max_k, top_n = 200L)
  f3_family_k <- fam$k
  f3_family_sil <- fam$silhouette

  # ── F4: subinterval discordance (3-bin PC1 disagreement) ──
  f4_subinterval <- subinterval_discordance(Xi, pos_i, n_bins = PARAMS$subinterval_n_bins)

  # ── F5: within-class substructure (PC1 VE inside each karyotype) ──
  # Use FULL X (all markers) restricted to samples in each group
  f5_within <- within_class_substructure(t(X), groups, min_n = PARAMS$min_samples_per_group)

  # ── F6: HWE fit ──
  n_h1 <- sum(groups == "HOMO_1"); n_het <- sum(groups == "HET"); n_h2 <- sum(groups == "HOMO_2")
  hw <- hwe_gof(n_h1, n_het, n_h2)
  f6_hwe_p <- hw$p_value

  # ── F7: monotonic fraction ──
  grp_means <- matrix(NA_real_, nrow = nrow(X), ncol = 3)
  colnames(grp_means) <- c("HOMO_1", "HET", "HOMO_2")
  for (g in colnames(grp_means)) {
    idx <- which(groups == g)
    if (length(idx) >= PARAMS$min_samples_per_group)
      grp_means[, g] <- rowMeans(X[, idx, drop = FALSE], na.rm = TRUE)
  }
  f7_monotonic <- monotonic_fraction(grp_means)

  # ── F8: polarity block count ──
  delta <- grp_means[, "HOMO_2"] - grp_means[, "HOMO_1"]
  delta[!is.finite(delta)] <- 0
  # Smoothed sign over windows of min_block_markers
  hw_win <- PARAMS$min_block_markers %/% 2L
  sm_sign <- numeric(length(delta))
  for (i in seq_along(delta)) {
    lo <- max(1, i - hw_win); hi <- min(length(delta), i + hw_win)
    local <- delta[lo:hi]
    w <- abs(local)
    pos_w <- sum(w[local > 0]); neg_w <- sum(w[local < 0])
    sm_sign[i] <- if (pos_w >= neg_w) 1 else -1
  }
  block_changes <- sum(abs(diff(sm_sign)) > 0)
  f8_block_count <- block_changes + 1L

  # ── F9: technical score ──
  f9_technical <- technical_score(sites_reg$pos, interval_start, interval_end,
                                   gap_file = gap_file)

  # ── F10: family-LD correlation ──
  # Use global PC1 score per sample as the partition summary
  pc_global <- tryCatch({
    Mi <- Xi
    for (j in seq_len(ncol(Mi))) {
      bad <- !is.finite(Mi[, j])
      if (any(bad)) Mi[bad, j] <- mean(Mi[, j], na.rm = TRUE)
    }
    Mi[!is.finite(Mi)] <- 0
    v <- apply(Mi, 1, safe_var)
    kk <- which(v > 1e-10)
    if (length(kk) < 2) return(NA_real_)
    prcomp(t(Mi[kk, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1)$x[, 1]
  }, error = function(e) rep(NA_real_, length(sample_cols)))
  names(pc_global) <- sample_cols
  f10_family_ld <- family_ld_correlation(pc_global, kinship_pcs)

  list(
    simple_axis_score          = f1_simple_axis,
    boundary_consistency       = f2_boundary,
    marker_family_k            = f3_family_k,
    marker_family_silhouette   = f3_family_sil,
    subinterval_discordance    = f4_subinterval,
    within_class_substructure  = f5_within,
    hwe_p_value                = f6_hwe_p,
    monotonic_fraction         = f7_monotonic,
    polarity_block_count       = f8_block_count,
    technical_score            = f9_technical,
    family_ld_correlation      = f10_family_ld,
    n_markers                  = nrow(X),
    n_samples                  = ncol(X),
    n_hom1                     = n_h1, n_het = n_het, n_hom2 = n_h2,
    pc_global_score            = list(pc_global)
  )
}

# =============================================================================
# CLASSIFICATION: features → (class, tier)
# =============================================================================
classify_architecture <- function(feat) {
  f <- feat
  # Technical-first gate
  if (isTRUE(f$technical_score >= PARAMS$thr_technical_high)) {
    return(list(class = "TECHNICAL", tier = "D",
                reason = sprintf("technical_score=%.2f >= %.2f",
                                 f$technical_score, PARAMS$thr_technical_high)))
  }

  # Demographic-block gate (strong family-LD, low monotonicity)
  if (isTRUE(f$family_ld_correlation >= PARAMS$thr_family_ld_corr_high) &&
      isTRUE(f$monotonic_fraction < PARAMS$thr_monotonic_weak)) {
    return(list(class = "DEMOGRAPHIC_BLOCK", tier = "C",
                reason = sprintf("family_ld=%.2f high; monotonic=%.2f low",
                                 f$family_ld_correlation, f$monotonic_fraction)))
  }

  # Core feature checks
  simple_axis_strong  <- isTRUE(f$simple_axis_score >= PARAMS$thr_simple_axis_strong)
  simple_axis_weak    <- isTRUE(f$simple_axis_score >= PARAMS$thr_simple_axis_weak)
  boundary_strong     <- isTRUE(f$boundary_consistency >= PARAMS$thr_boundary_consistency_strong)
  boundary_weak       <- isTRUE(f$boundary_consistency >= PARAMS$thr_boundary_consistency_weak)
  monotonic_strong    <- isTRUE(f$monotonic_fraction >= PARAMS$thr_monotonic_strong)
  monotonic_weak      <- isTRUE(f$monotonic_fraction >= PARAMS$thr_monotonic_weak)
  hwe_ok              <- isTRUE(f$hwe_p_value >= PARAMS$thr_hwe_p_ok)
  within_high         <- isTRUE(f$within_class_substructure >= PARAMS$thr_within_substructure_high)
  subint_high         <- isTRUE(f$subinterval_discordance >= PARAMS$thr_subinterval_discord_high)
  multifamily         <- isTRUE(f$marker_family_k >= 2L)

  # SIMPLE STRONG: everything clean
  if (simple_axis_strong && boundary_strong && monotonic_strong &&
      !multifamily && !within_high && !subint_high) {
    tier <- if (hwe_ok) "A" else "B"
    return(list(class = "SIMPLE_STRONG", tier = tier,
                reason = "axis + boundary + monotonic all strong; single family"))
  }

  # COMPOSITE OVERLAPPING: left/middle/right disagree strongly AND ≥2 families
  # with distinct sample partitions
  if (subint_high && multifamily) {
    return(list(class = "COMPOSITE_OVERLAP", tier = "B",
                reason = sprintf("subinterval_discord=%.2f + %d marker families",
                                 f$subinterval_discordance, f$marker_family_k)))
  }

  # COMPOSITE INTERNAL: main axis exists but within-class substructure OR
  # multiple marker families (without left/right disagreement). Must be
  # checked BEFORE SIMPLE_WEAK so substructure is not swallowed by a clean axis.
  if (simple_axis_weak && (within_high || multifamily)) {
    return(list(class = "COMPOSITE_INTERNAL", tier = "B",
                reason = sprintf("axis=%.2f + within_subst=%.2f + families=%d",
                                 f$simple_axis_score, f$within_class_substructure,
                                 f$marker_family_k)))
  }

  # SIMPLE WEAK: 3-state visible but low contrast / noisy (fallthrough only)
  if (simple_axis_weak && monotonic_weak && !multifamily && !subint_high) {
    return(list(class = "SIMPLE_WEAK", tier = "B",
                reason = sprintf("axis=%.2f monotonic=%.2f; 3-state weak",
                                 f$simple_axis_score, f$monotonic_fraction)))
  }

  # Fallback: unresolved
  list(class = "UNRESOLVED", tier = "C",
       reason = "mixed evidence — see feature table")
}

# =============================================================================
# MAIN: per-candidate classification
# =============================================================================
classify_candidate <- function(cid, chr, c_start, c_end,
                                cand_dir, dosage_dir,
                                gap_file = NULL, kinship_pcs = NULL) {
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (!file.exists(rot_file)) {
    return(list(status = "skip", reason = "rotated PCA missing"))
  }
  rot <- fread(rot_file)
  if (nrow(rot) < 5) return(list(status = "skip", reason = "<5 samples"))

  dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(dosage_dir, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    return(list(status = "skip", reason = "dosage files missing"))
  }
  dos <- fread(dos_file); sites <- fread(sites_file)
  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < PARAMS$min_markers_for_classification) {
    return(list(status = "skip", reason = sprintf("n_markers=%d < %d",
                                                   length(keep),
                                                   PARAMS$min_markers_for_classification)))
  }
  sites_reg <- sites[keep]; dos_reg <- dos[keep]

  sample_cols <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))
  if (length(sample_cols) < 10) {
    # Ind-prefix fallback
    ind_cols <- setdiff(names(dos_reg), "marker")
    if (all(grepl("^Ind", ind_cols)) && length(ind_cols) == nrow(rot)) {
      setnames(dos_reg, old = ind_cols, new = rot$sample)
      sample_cols <- rot$sample
    } else {
      return(list(status = "skip", reason = "<10 matched samples"))
    }
  }

  X <- as.matrix(dos_reg[, ..sample_cols]); storage.mode(X) <- "double"
  groups <- rot[match(sample_cols, sample), coarse_group_refined]

  # ── Feature extraction ────────────────────────────────────────────────
  feat <- compute_features(X, sites_reg, rot, groups, sample_cols,
                           c_start, c_end,
                           gap_file = gap_file, kinship_pcs = kinship_pcs)

  # ── Architecture classification ──────────────────────────────────────
  cls <- classify_architecture(feat)

  # ── "Separate first, correct second" — sample subgroup partitioning ──
  # Only run subgroup decomposition when the class flags composite/unresolved
  do_separate <- cls$class %in% c("COMPOSITE_INTERNAL", "COMPOSITE_OVERLAP", "UNRESOLVED")
  sep <- NULL
  if (do_separate) {
    # Restrict to informative markers for subgroup partitioning
    mv <- apply(X, 1, safe_var)
    info_idx <- which(mv >= quantile(mv, 1 - PARAMS$top_marker_variance_fraction, na.rm = TRUE))
    if (length(info_idx) < 30) info_idx <- order(mv, decreasing = TRUE)[seq_len(min(50, nrow(X)))]
    sep <- separate_first_partition(X[info_idx, , drop = FALSE], groups,
                                     max_subgroups = 4L)
  }

  list(
    status = "ok",
    candidate_id = cid, chrom = chr, start_bp = c_start, end_bp = c_end,
    class = cls$class, tier = cls$tier, label = paste0(cls$class, "-", cls$tier),
    reason = cls$reason,
    features = feat,
    separation = sep,
    sample_cols = sample_cols,
    groups = groups
  )
}

# =============================================================================
# WRITE PER-CANDIDATE OUTPUTS
# =============================================================================
write_outputs <- function(res, cand_dir) {
  if (res$status != "ok") return(invisible())
  ensure_dir(cand_dir)

  # 1. architecture label
  arch_dt <- data.table(
    candidate_id = res$candidate_id,
    chrom = res$chrom,
    start_bp = res$start_bp, end_bp = res$end_bp,
    architecture_class = res$class,
    tier = res$tier,
    label = res$label,
    reason = res$reason,
    simple_axis_score          = round(res$features$simple_axis_score, 4),
    boundary_consistency       = round(res$features$boundary_consistency, 4),
    marker_family_k            = res$features$marker_family_k,
    marker_family_silhouette   = round(res$features$marker_family_silhouette, 4),
    subinterval_discordance    = round(res$features$subinterval_discordance, 4),
    within_class_substructure  = round(res$features$within_class_substructure, 4),
    hwe_p_value                = round(res$features$hwe_p_value, 4),
    monotonic_fraction         = round(res$features$monotonic_fraction, 4),
    polarity_block_count       = res$features$polarity_block_count,
    technical_score            = round(res$features$technical_score, 4),
    family_ld_correlation      = round(res$features$family_ld_correlation, 4),
    n_markers = res$features$n_markers,
    n_samples = res$features$n_samples,
    n_hom1 = res$features$n_hom1,
    n_het  = res$features$n_het,
    n_hom2 = res$features$n_hom2
  )
  fwrite(arch_dt, file.path(cand_dir, "candidate_architecture_class.tsv"), sep = "\t")

  # 2. sample subgroup table (only when separation was performed)
  if (!is.null(res$separation)) {
    subg_dt <- data.table(
      candidate_id = res$candidate_id,
      sample = res$sample_cols,
      coarse_group = res$groups,
      subgroup_id = res$separation$subgroup
    )
    # Attach per-family scores
    fs <- res$separation$family_scores
    for (j in seq_len(ncol(fs))) {
      subg_dt[[colnames(fs)[j]]] <- round(fs[, j], 4)
    }
    fwrite(subg_dt, file.path(cand_dir, "candidate_sample_subgroups.tsv"), sep = "\t")

    # 3. marker family table
    fam <- res$separation$marker_families
    fam_dt <- data.table(
      candidate_id = res$candidate_id,
      marker_idx = seq_along(fam$cluster),
      marker_family = fam$cluster
    )
    fwrite(fam_dt, file.path(cand_dir, "candidate_marker_families.tsv"), sep = "\t")
  }

  invisible()
}

# =============================================================================
# DRIVER
# =============================================================================
main <- function() {
  if (!file.exists(CANDIDATE_TABLE)) {
    stop("CANDIDATE_TABLE not found: ", CANDIDATE_TABLE)
  }
  cand <- fread(CANDIDATE_TABLE)
  if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

  # Optional global files
  gap_file <- if (exists("ASSEMBLY_GAP_FILE") && file.exists(ASSEMBLY_GAP_FILE)) ASSEMBLY_GAP_FILE else NULL
  kinship_file <- if (exists("KINSHIP_PCS_FILE") && file.exists(KINSHIP_PCS_FILE)) KINSHIP_PCS_FILE else NULL
  kinship_pcs <- if (!is.null(kinship_file)) tryCatch(fread(kinship_file), error = function(e) NULL) else NULL

  catalog <- list()
  for (ci in seq_len(nrow(cand))) {
    row <- cand[ci]
    cid <- row$candidate_id; chr <- row$chrom
    c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)
    cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))

    message(sprintf("[INFO] Candidate %d (%s:%s-%s)", cid, chr,
                    format(c_start, big.mark = ","), format(c_end, big.mark = ",")))
    res <- tryCatch(
      classify_candidate(cid, chr, c_start, c_end, cand_dir, DOSAGE_DIR,
                         gap_file = gap_file, kinship_pcs = kinship_pcs),
      error = function(e) list(status = "error", reason = conditionMessage(e))
    )

    if (res$status != "ok") {
      message("  [SKIP] ", res$reason)
      next
    }
    write_outputs(res, cand_dir)
    catalog[[length(catalog) + 1]] <- data.table(
      candidate_id = cid, chrom = chr, start_bp = c_start, end_bp = c_end,
      label = res$label, class = res$class, tier = res$tier,
      reason = res$reason
    )
    message(sprintf("  [OK] %s — %s", res$label, res$reason))
  }

  # Catalog-level output
  if (length(catalog) > 0) {
    ensure_dir(FOLLOWUP_DIR)
    cat_dt <- rbindlist(catalog)
    fwrite(cat_dt, file.path(FOLLOWUP_DIR, "architecture_catalog.tsv"), sep = "\t")
    message("[DONE] Catalog written: ", file.path(FOLLOWUP_DIR, "architecture_catalog.tsv"))
    message("  ", nrow(cat_dt), " candidates classified")
    message("  Class distribution:")
    tbl <- cat_dt[, .N, by = .(class, tier)][order(class, tier)]
    for (i in seq_len(nrow(tbl))) {
      message(sprintf("    %-20s %s  n=%d", tbl$class[i], tbl$tier[i], tbl$N[i]))
    }
  } else {
    message("[WARN] No candidates classified")
  }
}

main()
