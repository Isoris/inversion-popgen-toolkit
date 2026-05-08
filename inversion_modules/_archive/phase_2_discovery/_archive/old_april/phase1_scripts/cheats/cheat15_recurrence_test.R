#!/usr/bin/env Rscript
# =============================================================================
# cheat15_recurrence_test.R — Recurrence vs single-origin test
#
# BIOLOGY:
#   Single-origin (NHEJ): ALL carriers share one ancestral haplotype at
#   breakpoint flanks, diversified only by new mutations.
#   Recurrent (NAHR): carriers can have DIFFERENT haplotype backgrounds
#   at breakpoint flanks because the inversion arose independently.
#
# INPUT:  Clair3 phased VCFs (or dosage matrix), decomposition classes,
#         boundary catalog, Cheat 14 mechanism class
# OUTPUT: per-candidate origin classification, carrier haplotype clusters
# NOTE:   At 9× coverage phase blocks are short; uses dosage IBS fallback
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
CARRIER_FLANK_BP    <- 20000L
MIN_CARRIERS        <- 5L
MIN_VARIANTS        <- 10L
BIC_DELTA_THRESHOLD <- 6

# ── Extract carrier haplotypes (dosage IBS fallback) ──────────────────

extract_carrier_haplotypes <- function(dosage_matrix, carrier_ids,
                                        window_idx = NULL) {
  if (is.null(dosage_matrix) || length(carrier_ids) < MIN_CARRIERS)
    return(NULL)
  # Subset to carriers
  avail <- intersect(carrier_ids, colnames(dosage_matrix))
  if (length(avail) < MIN_CARRIERS) return(NULL)
  mat <- dosage_matrix[, avail, drop = FALSE]
  if (!is.null(window_idx)) mat <- mat[window_idx, , drop = FALSE]
  # Filter invariant sites
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[vars > 0, , drop = FALSE]
  if (nrow(mat) < MIN_VARIANTS) return(NULL)
  t(mat)  # carriers × variants
}

# ── Pairwise IBS distance between carriers ────────────────────────────

compute_carrier_pairwise_ibs <- function(hap_mat) {
  if (is.null(hap_mat) || nrow(hap_mat) < MIN_CARRIERS) return(NULL)
  n <- nrow(hap_mat)
  d <- matrix(0, n, n, dimnames = list(rownames(hap_mat), rownames(hap_mat)))
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      valid <- !is.na(hap_mat[i,]) & !is.na(hap_mat[j,])
      if (sum(valid) < MIN_VARIANTS) { d[i,j] <- d[j,i] <- NA; next }
      ibs <- mean(abs(hap_mat[i, valid] - hap_mat[j, valid]))
      d[i,j] <- d[j,i] <- ibs
    }
  }
  as.dist(d)
}

# ── BIC-based k=1 vs k=2 clustering ──────────────────────────────────

compute_carrier_haplotype_diversity <- function(hap_mat) {
  if (is.null(hap_mat)) return(list(n_clusters = NA, bic_delta = NA,
    cluster_assignments = NULL, silhouette = NA))

  d <- compute_carrier_pairwise_ibs(hap_mat)
  if (is.null(d)) return(list(n_clusters = NA, bic_delta = NA,
    cluster_assignments = NULL, silhouette = NA))

  n <- nrow(hap_mat)
  # Convert to coordinate space for k-means
  mds <- tryCatch(cmdscale(d, k = min(3, n - 1)), error = function(e) NULL)
  if (is.null(mds)) return(list(n_clusters = NA, bic_delta = NA,
    cluster_assignments = NULL, silhouette = NA))

  # k=1: single cluster
  rss1 <- sum(scale(mds, scale = FALSE)^2)
  bic1 <- n * log(rss1 / n) + 1 * log(n)

  # k=2
  if (n < 2 * MIN_CARRIERS) {
    return(list(n_clusters = 1, bic_delta = 0,
                cluster_assignments = setNames(rep(1L, n), rownames(hap_mat)),
                silhouette = NA))
  }
  km2 <- tryCatch(kmeans(mds, centers = 2, nstart = 10),
                   error = function(e) NULL)
  if (is.null(km2))
    return(list(n_clusters = 1, bic_delta = 0,
                cluster_assignments = setNames(rep(1L, n), rownames(hap_mat)),
                silhouette = NA))

  bic2 <- n * log(km2$tot.withinss / n) + 2 * log(n)
  delta <- bic1 - bic2  # positive = k=2 better

  sil <- tryCatch({
    s <- cluster::silhouette(km2$cluster, d)
    mean(s[, "sil_width"])
  }, error = function(e) NA_real_)

  best_k <- if (delta > BIC_DELTA_THRESHOLD && sil > 0.2) 2L else 1L

  list(n_clusters = best_k, bic_delta = round(delta, 2),
       cluster_assignments = setNames(km2$cluster, rownames(hap_mat)),
       silhouette = round(sil, 3),
       within_div = if (best_k == 2) {
         c1 <- names(km2$cluster[km2$cluster == 1])
         c2 <- names(km2$cluster[km2$cluster == 2])
         list(cluster1 = length(c1), cluster2 = length(c2))
       } else NULL)
}

# ── Decision matrix: mechanism × clusters → origin class ─────────────

test_recurrence <- function(carrier_diversity, mechanism_class = "UNKNOWN") {
  nc <- carrier_diversity$n_clusters
  if (is.na(nc)) return("INSUFFICIENT_DATA")

  if (mechanism_class %in% c("NAHR_CANDIDATE", "NAHR_POSSIBLE")) {
    if (nc >= 2) return("RECURRENT")
    return("SINGLE_ORIGIN_DESPITE_REPEATS")
  }
  if (mechanism_class == "NHEJ_CANDIDATE") {
    if (nc == 1) return("SINGLE_ORIGIN")
    return("UNEXPECTED_MULTI_CLUSTER")
  }
  # Unknown mechanism
  if (nc >= 2) return("POSSIBLE_RECURRENCE")
  return("SINGLE_ORIGIN")
}

# ── Search mode ────────────────────────────────────────────────────────

search_recurrence <- function(chr, zone_start, zone_end, ...) {
  data.table(method = "recurrence_test", best_bp = NA_integer_,
             score = 0, is_precise = FALSE,
             detail = "search_not_applicable_for_recurrence")
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat15 <- function(chr, left_bp, right_bp, dosage_matrix,
                         carrier_ids, mechanism_class = "UNKNOWN") {
  message("[cheat15] ", chr, ":", round(left_bp/1e6,1), "-",
          round(right_bp/1e6,1), " Mb | ", length(carrier_ids),
          " carriers | mechanism=", mechanism_class)

  if (length(carrier_ids) < MIN_CARRIERS) {
    message("[cheat15] Too few carriers (", length(carrier_ids), " < ", MIN_CARRIERS, ")")
    return(list(origin_class = "INSUFFICIENT_DATA",
                carrier_diversity = list(n_clusters = NA),
                search_result = search_recurrence(chr, left_bp, right_bp)))
  }

  hap_mat <- extract_carrier_haplotypes(dosage_matrix, carrier_ids)
  diversity <- compute_carrier_haplotype_diversity(hap_mat)
  origin <- test_recurrence(diversity, mechanism_class)

  message("[cheat15] Clusters: ", diversity$n_clusters,
          " | BIC Δ: ", diversity$bic_delta,
          " | Silhouette: ", diversity$silhouette,
          " → ", origin)

  list(origin_class = origin, carrier_diversity = diversity,
       search_result = search_recurrence(chr, left_bp, right_bp))
}
