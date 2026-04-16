#!/usr/bin/env Rscript
# =============================================================================
# cheat22_secondary_structure.R — Secondary structure within inversion classes
#
# BIOLOGY:
#   After decomposition assigns HOM_INV, those samples may not be
#   homogeneous. They can sub-cluster by family/founder background,
#   nested sub-inversions, multiple origins, or gene conversion tracts.
#   Detecting this reveals whether "one inversion" model is sufficient.
#
# INPUT:  precomp RDS, decomposition classes, Q matrix, candidate regions
# OUTPUT: per-class subclustering, segment concordance, LD heterogeneity
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
MIN_CLASS_SIZE    <- 8L
PVE_THRESHOLD     <- 0.15   # within-class PC1 must explain > this
CRAMERS_THRESHOLD <- 0.3    # association with ancestry
K_WITHIN_RANGE    <- 1:4

# ── Within-class PCA ──────────────────────────────────────────────────

within_class_pca <- function(dt, sample_names, decomp_classes, target_class,
                              candidate_start, candidate_end) {
  empty <- list(has_secondary_structure = FALSE, n_subclusters = 1L,
                subcluster_assignments = NULL, variance_explained = NA_real_,
                ancestry_association = NA_real_)

  # Get samples in target class
  target_ids <- names(decomp_classes[decomp_classes == target_class])
  target_ids <- intersect(target_ids, sample_names)
  if (length(target_ids) < MIN_CLASS_SIZE) return(empty)

  # Extract per-window scores for target samples within region
  region_dt <- dt[start_bp >= candidate_start & end_bp <= candidate_end]
  if (nrow(region_dt) == 0) return(empty)

  # Build sample × window matrix
  # Try from PC1_loadings or inv_likeness
  if ("PC1_loadings" %in% names(region_dt) &&
      is.list(region_dt$PC1_loadings)) {
    full_mat <- do.call(rbind, region_dt$PC1_loadings)
    avail <- intersect(target_ids, colnames(full_mat))
    if (length(avail) < MIN_CLASS_SIZE) return(empty)
    mat <- t(full_mat[, avail, drop = FALSE])  # samples × windows
  } else {
    return(empty)
  }

  # Remove constant columns
  col_var <- apply(mat, 2, var, na.rm = TRUE)
  mat <- mat[, col_var > 0, drop = FALSE]
  if (ncol(mat) < 3) return(empty)

  # PCA on within-class subset
  pca <- tryCatch(prcomp(mat, center = TRUE, scale. = FALSE),
                   error = function(e) NULL)
  if (is.null(pca)) return(empty)

  pve1 <- pca$sdev[1]^2 / sum(pca$sdev^2)

  # Test k=2 on within-class PC1
  pc1_vals <- pca$x[, 1]
  if (pve1 < PVE_THRESHOLD)
    return(list(has_secondary_structure = FALSE, n_subclusters = 1L,
                subcluster_assignments = setNames(rep(1L, length(avail)), avail),
                variance_explained = round(pve1, 4),
                ancestry_association = NA_real_))

  # BIC-based k selection
  x <- matrix(pc1_vals, ncol = 1)
  n <- nrow(x)
  best_k <- 1L; best_bic <- n * log(var(pc1_vals)) + log(n)
  assignments <- setNames(rep(1L, n), avail)

  for (k in 2:min(4, floor(n / MIN_CLASS_SIZE))) {
    km <- tryCatch(kmeans(x, centers = k, nstart = 20),
                    error = function(e) NULL)
    if (is.null(km) || min(table(km$cluster)) < 3) next
    bic_k <- n * log(km$tot.withinss / n) + k * log(n)
    if (bic_k < best_bic - 6) {
      best_bic <- bic_k; best_k <- k
      assignments <- setNames(km$cluster, avail)
    }
  }

  list(has_secondary_structure = best_k > 1,
       n_subclusters = best_k,
       subcluster_assignments = assignments,
       variance_explained = round(pve1, 4),
       ancestry_association = NA_real_)
}

# ── Ancestry association (Cramér's V) ────────────────────────────────

test_ancestry_association <- function(subcluster_assignments, q_matrix) {
  if (is.null(subcluster_assignments) || is.null(q_matrix)) return(NA_real_)
  avail <- intersect(names(subcluster_assignments), rownames(q_matrix))
  if (length(avail) < 10) return(NA_real_)

  clusters <- subcluster_assignments[avail]
  if (length(unique(clusters)) < 2) return(NA_real_)

  # Assign Q-group: dominant ancestry component
  q_groups <- apply(q_matrix[avail, , drop = FALSE], 1, which.max)

  # Cramér's V
  tab <- table(clusters, q_groups)
  chi <- tryCatch(chisq.test(tab)$statistic, error = function(e) NA_real_)
  if (is.na(chi)) return(NA_real_)
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  v <- sqrt(chi / (n * (k - 1)))
  round(as.numeric(v), 3)
}

# ── Center vs flank structure ─────────────────────────────────────────

center_vs_flank_structure <- function(dt, sample_names, decomp_classes,
                                       target_class, left_bp, right_bp) {
  span <- right_bp - left_bp
  if (span <= 0) return(list(segment_concordance = NA_real_))

  # Three segments: left 25%, center 50%, right 25%
  seg_bounds <- list(
    left_flank  = c(left_bp, left_bp + span * 0.25),
    center      = c(left_bp + span * 0.25, left_bp + span * 0.75),
    right_flank = c(left_bp + span * 0.75, right_bp)
  )

  assignments_list <- list()
  for (seg_name in names(seg_bounds)) {
    bounds <- seg_bounds[[seg_name]]
    res <- within_class_pca(dt, sample_names, decomp_classes, target_class,
                             bounds[1], bounds[2])
    if (!is.null(res$subcluster_assignments) && length(res$subcluster_assignments) > 0) {
      assignments_list[[seg_name]] <- res$subcluster_assignments
    }
  }

  if (length(assignments_list) < 2)
    return(list(segment_concordance = NA_real_))

  # Compute concordance: fraction of sample pairs maintaining cluster
  all_ids <- Reduce(intersect, lapply(assignments_list, names))
  if (length(all_ids) < 5) return(list(segment_concordance = NA_real_))

  pairs_concordant <- 0L; pairs_total <- 0L
  for (i in 1:(length(all_ids)-1)) {
    for (j in (i+1):length(all_ids)) {
      sa <- all_ids[i]; sb <- all_ids[j]
      same_in_all <- TRUE
      for (seg in names(assignments_list)) {
        a_clust <- assignments_list[[seg]][sa]
        b_clust <- assignments_list[[seg]][sb]
        if (!is.na(a_clust) && !is.na(b_clust)) {
          same_in_all <- same_in_all && (a_clust == b_clust)
        }
      }
      pairs_total <- pairs_total + 1L
      if (same_in_all) pairs_concordant <- pairs_concordant + 1L
    }
  }

  concordance <- if (pairs_total > 0) pairs_concordant / pairs_total else NA_real_
  list(segment_concordance = round(concordance, 3),
       interpretation = if (is.na(concordance)) "insufficient"
         else if (concordance > 0.8) "stable_founder_haplotypes"
         else if (concordance > 0.5) "partial_recombination"
         else "mosaic_interior")
}

# ── Search mode ────────────────────────────────────────────────────────

search_secondary_structure <- function(chr, zone_start, zone_end,
                                        dt = NULL, sample_names = NULL,
                                        decomp_classes = NULL, ...) {
  empty <- data.table(method = "secondary_structure", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(dt) || is.null(decomp_classes)) return(empty)

  res <- within_class_pca(dt, sample_names, decomp_classes, "HOM_INV",
                           zone_start, zone_end)
  sc <- if (res$has_secondary_structure) 0.8 else 0.1
  data.table(method = "secondary_structure",
             best_bp = as.integer((zone_start + zone_end) / 2),
             score = round(sc, 3), is_precise = FALSE,
             detail = paste0("k=", res$n_subclusters,
                              ",pve=", res$variance_explained))
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat22 <- function(chr, candidate_start, candidate_end,
                         left_bp, right_bp, dt, sample_names,
                         decomp_classes, q_matrix = NULL) {
  message("[cheat22] ", chr, ":", round(candidate_start/1e6,1), "-",
          round(candidate_end/1e6,1), " Mb")

  # Within-class PCA for HOM_INV
  inv_pca <- within_class_pca(dt, sample_names, decomp_classes, "HOM_INV",
                               candidate_start, candidate_end)
  message("[cheat22] HOM_INV: secondary=", inv_pca$has_secondary_structure,
          " k=", inv_pca$n_subclusters,
          " PVE=", inv_pca$variance_explained)

  # Ancestry association
  cramers <- test_ancestry_association(inv_pca$subcluster_assignments, q_matrix)
  message("[cheat22] Cramér's V (ancestry): ", cramers)

  # Center vs flank
  cvf <- center_vs_flank_structure(dt, sample_names, decomp_classes,
                                    "HOM_INV", left_bp, right_bp)
  message("[cheat22] Segment concordance: ", cvf$segment_concordance,
          " → ", cvf$interpretation)

  list(within_class = inv_pca,
       ancestry_association = cramers,
       center_vs_flank = cvf,
       search_result = search_secondary_structure(chr, candidate_start,
                        candidate_end, dt, sample_names, decomp_classes))
}
