#!/usr/bin/env Rscript
# =============================================================================
# cheat16_multi_arrangement.R — Multi-arrangement detection (k>3)
#
# BIOLOGY:
#   Sequential inversions can accumulate over time, creating multi-
#   arrangement complexes with shared breakpoints and nested structure.
#   Signs: k>3 in PCA, breakpoint reuse, arrangement families.
#
# INPUT:  precomp RDS, boundary catalog, decomposition output
# OUTPUT: optimal k, breakpoint reuse table, arrangement graph,
#         sequential model classification
# REQUIRES: igraph (optional, for arrangement graph)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
K_RANGE_DEFAULT     <- 2:8
BIC_DELTA_MULTI     <- 10      # BIC improvement required to prefer k>3
MIN_SAMPLES_PER_K   <- 5       # each cluster must have >= this
BP_REUSE_WINDOW     <- 50000L  # ±50 kb for breakpoint matching

# ── k-selection via BIC on per-sample mean PC1 ───────────────────────

test_k_selection <- function(dt, sample_names, candidate_start, candidate_end,
                              k_range = K_RANGE_DEFAULT) {
  # Extract per-sample mean PC1 loading in region
  region_dt <- dt[start_bp >= candidate_start & end_bp <= candidate_end]
  if (nrow(region_dt) == 0 || !"PC1_loadings" %in% names(region_dt))
    return(list(best_k = 3L, bic_values = NULL, is_multi = FALSE,
                evidence_strength = 0))

  # Try to build sample × window matrix from precomp
  # Fallback: use inv_likeness per window as 1D feature
  if ("sample_PC1" %in% names(attributes(region_dt))) {
    sample_vals <- attr(region_dt, "sample_PC1")
  } else if ("PC1_loadings" %in% names(region_dt)) {
    # PC1_loadings is stored as a list column
    mat <- do.call(rbind, region_dt$PC1_loadings)
    sample_vals <- colMeans(mat, na.rm = TRUE)
  } else {
    return(list(best_k = 3L, bic_values = NULL, is_multi = FALSE,
                evidence_strength = 0))
  }

  n <- length(sample_vals)
  if (n < max(k_range) * MIN_SAMPLES_PER_K)
    k_range <- 2:floor(n / MIN_SAMPLES_PER_K)
  if (length(k_range) == 0 || max(k_range) < 2)
    return(list(best_k = 3L, bic_values = NULL, is_multi = FALSE,
                evidence_strength = 0))

  x <- matrix(sample_vals, ncol = 1)
  bic_vals <- numeric(length(k_range))
  sil_vals <- numeric(length(k_range))
  names(bic_vals) <- names(sil_vals) <- as.character(k_range)

  for (ki in seq_along(k_range)) {
    k <- k_range[ki]
    km <- tryCatch(kmeans(x, centers = k, nstart = 20, iter.max = 100),
                    error = function(e) NULL)
    if (is.null(km) || min(table(km$cluster)) < 2) {
      bic_vals[ki] <- Inf; sil_vals[ki] <- -1; next
    }
    bic_vals[ki] <- n * log(km$tot.withinss / n) + k * log(n)
    sil_vals[ki] <- tryCatch(
      mean(cluster::silhouette(km$cluster, dist(x))[, "sil_width"]),
      error = function(e) -1)
  }

  valid <- is.finite(bic_vals)
  if (sum(valid) == 0)
    return(list(best_k = 3L, bic_values = bic_vals, is_multi = FALSE,
                evidence_strength = 0))

  best_idx <- which.min(bic_vals)
  best_k   <- k_range[best_idx]

  # Evidence strength: BIC(k=3) - BIC(best_k)
  k3_idx <- which(k_range == 3)
  bic_k3 <- if (length(k3_idx) == 1 && is.finite(bic_vals[k3_idx]))
    bic_vals[k3_idx] else Inf
  strength <- bic_k3 - bic_vals[best_idx]

  list(best_k = best_k,
       bic_values = bic_vals,
       silhouette_values = sil_vals,
       is_multi = best_k > 3 && strength > BIC_DELTA_MULTI,
       evidence_strength = round(strength, 1))
}

# ── Breakpoint reuse table ────────────────────────────────────────────

build_breakpoint_reuse_table <- function(boundary_catalog, candidate_regions) {
  if (nrow(boundary_catalog) == 0 || nrow(candidate_regions) == 0)
    return(data.table(boundary_bp = integer(), n_candidates = 0L,
                       candidate_ids = character(), is_hub = FALSE))

  results <- list()
  for (i in seq_len(nrow(boundary_catalog))) {
    bp <- boundary_catalog$boundary_bp[i]
    chr_i <- boundary_catalog$chr[i]
    # Count how many candidate regions use this as a boundary (±window)
    matches <- candidate_regions[
      chr == chr_i &
      ((abs(start_bp - bp) <= BP_REUSE_WINDOW) |
       (abs(end_bp - bp) <= BP_REUSE_WINDOW))
    ]
    cand_ids <- if (nrow(matches) > 0 && "candidate_id" %in% names(matches))
      paste(matches$candidate_id, collapse = ";") else ""
    results[[i]] <- data.table(
      boundary_bp = bp, chr = chr_i,
      n_candidates = nrow(matches),
      candidate_ids = cand_ids,
      is_hub = nrow(matches) >= 2)
  }
  rbindlist(results)
}

# ── Arrangement graph (requires igraph) ───────────────────────────────

build_arrangement_graph <- function(boundary_catalog, decomp_classes,
                                     breakpoint_reuse) {
  has_igraph <- requireNamespace("igraph", quietly = TRUE)
  if (!has_igraph) {
    message("[cheat16] igraph not available — skipping graph construction")
    return(NULL)
  }

  hubs <- breakpoint_reuse[is_hub == TRUE]
  if (nrow(hubs) == 0) return(NULL)

  # Build edges: candidate pairs sharing a hub breakpoint
  edges <- list()
  for (i in seq_len(nrow(hubs))) {
    cands <- strsplit(hubs$candidate_ids[i], ";")[[1]]
    if (length(cands) < 2) next
    for (a in 1:(length(cands)-1)) {
      for (b in (a+1):length(cands)) {
        edges[[length(edges)+1]] <- data.table(
          from = cands[a], to = cands[b],
          shared_bp = hubs$boundary_bp[i])
      }
    }
  }
  if (length(edges) == 0) return(NULL)
  edge_dt <- rbindlist(edges)
  g <- igraph::graph_from_data_frame(edge_dt, directed = FALSE)
  g
}

# ── Detect sequential pattern ─────────────────────────────────────────

detect_sequential_inversions <- function(arrangement_graph) {
  if (is.null(arrangement_graph)) return("simple")
  has_igraph <- requireNamespace("igraph", quietly = TRUE)
  if (!has_igraph) return("unknown")

  n_v <- igraph::vcount(arrangement_graph)
  n_e <- igraph::ecount(arrangement_graph)
  if (n_v <= 1) return("simple")

  degrees <- igraph::degree(arrangement_graph)
  max_deg <- max(degrees)

  # Linear chain: all degrees <= 2 and graph is a path
  if (max_deg <= 2 && n_e == n_v - 1) return("linear_chain")
  # Star: one hub, rest degree 1
  if (max_deg >= 3 && sum(degrees == 1) >= n_v - 1) return("star")
  # Otherwise complex
  "complex"
}

# ── Search mode ────────────────────────────────────────────────────────

search_multi_arrangement <- function(chr, zone_start, zone_end, dt = NULL,
                                      sample_names = NULL, ...) {
  empty <- data.table(method = "multi_arrangement", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(dt) || nrow(dt) == 0) return(empty)

  ks <- test_k_selection(dt, sample_names, zone_start, zone_end)
  sc <- if (ks$is_multi) pmin(1, ks$evidence_strength / 50) else 0
  data.table(method = "multi_arrangement",
             best_bp = as.integer((zone_start + zone_end) / 2),
             score = round(sc, 3), is_precise = FALSE,
             detail = paste0("best_k=", ks$best_k, ",BIC_delta=",
                              ks$evidence_strength))
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat16 <- function(chr, candidate_start, candidate_end, dt,
                         sample_names, boundary_catalog = NULL,
                         candidate_regions = NULL) {
  message("[cheat16] ", chr, ":", round(candidate_start/1e6,1), "-",
          round(candidate_end/1e6,1), " Mb")

  # k-selection
  ks <- test_k_selection(dt, sample_names, candidate_start, candidate_end)
  message("[cheat16] Best k=", ks$best_k,
          " | Multi-arrangement: ", ks$is_multi,
          " | Evidence: ", ks$evidence_strength)

  # Breakpoint reuse
  reuse <- if (!is.null(boundary_catalog) && !is.null(candidate_regions))
    build_breakpoint_reuse_table(boundary_catalog, candidate_regions)
  else data.table()

  n_hubs <- if (nrow(reuse) > 0) sum(reuse$is_hub) else 0
  message("[cheat16] Breakpoint hubs: ", n_hubs)

  # Arrangement graph
  graph <- if (n_hubs > 0)
    build_arrangement_graph(boundary_catalog, NULL, reuse) else NULL
  model <- detect_sequential_inversions(graph)
  message("[cheat16] Sequential model: ", model)

  list(k_selection = ks, breakpoint_reuse = reuse,
       arrangement_graph = graph, sequential_model = model,
       search_result = search_multi_arrangement(chr, candidate_start,
                        candidate_end, dt, sample_names))
}
