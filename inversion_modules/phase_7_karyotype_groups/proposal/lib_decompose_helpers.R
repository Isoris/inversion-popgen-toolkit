#!/usr/bin/env Rscript
# =============================================================================
# lib_decompose_helpers.R — shared utilities for Phase 4b scripts
# =============================================================================
# Functions shared by STEP_C01i_decompose.R, STEP_C01i_b_multi_recomb.R,
# and STEP_C01i_d_seal.R. Source this at the top of each. No side effects,
# no globals — just functions.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a[1]))) b else a
}

# ─────────────────────────────────────────────────────────────────────────────
# Config / paths
# ─────────────────────────────────────────────────────────────────────────────

# Centralized path resolution. Each caller can override via env or args.
resolve_decomp_paths <- function(args = list()) {
  base <- args$base %||% Sys.getenv("BASE",
    "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04")
  registries <- args$registries %||% Sys.getenv("REGISTRIES",
    file.path(base, "inversion-popgen-toolkit/registries"))

  list(
    base         = base,
    registries   = registries,
    sample_reg   = file.path(registries, "data", "sample_registry"),
    evidence_reg = file.path(registries, "data", "evidence_registry"),
    precomp_dir  = args$precomp_dir %||% Sys.getenv("PRECOMP_DIR",
      file.path(base, "inversion_localpca_v7",
                "06_mds_candidates/snake_regions_multiscale/precomp")),
    clair3_dir   = args$clair3_dir %||% Sys.getenv("CLAIR3_DIR",
      file.path(base, "MODULE_4A_SNP_INDEL50_Clair3/postprocess_results")),
    sv_prior_dir = args$sv_prior_dir %||% Sys.getenv("FLASHLIGHT_DIR",
      file.path(base, "flashlight_v2/cache")),
    q_cache_dir = args$q_cache_dir %||% Sys.getenv("Q_CACHE_DIR",
      file.path(base, "unified_ancestry/local_Q"))
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Precomp loading — PC loadings per sample per window
# ─────────────────────────────────────────────────────────────────────────────

load_precomp <- function(chr, precomp_dir) {
  pf <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
  if (!file.exists(pf)) {
    message("[lib] no precomp for ", chr, " at ", pf)
    return(NULL)
  }
  readRDS(pf)
}

# Extract PC loadings for a candidate interval. Returns NULL on failure.
# Output: list(pc1_mat, pc2_mat, sample_names, n_windows)
extract_pc_loadings <- function(pc, start_bp, end_bp, min_windows = 5L) {
  if (is.null(pc) || is.null(pc$inv_likeness)) return(NULL)
  il <- pc$inv_likeness
  wins <- il[mid_bp >= start_bp & mid_bp <= end_bp]
  if (nrow(wins) < min_windows) return(NULL)

  pc1_cols <- grep("^PC_1_", names(wins), value = TRUE)
  pc2_cols <- grep("^PC_2_", names(wins), value = TRUE)
  if (length(pc1_cols) == 0) return(NULL)

  sample_names <- sub("^PC_1_", "", pc1_cols)
  pc1_mat <- as.matrix(wins[, ..pc1_cols])
  colnames(pc1_mat) <- sample_names

  pc2_mat <- if (length(pc2_cols) > 0) {
    m <- as.matrix(wins[, ..pc2_cols])
    colnames(m) <- sub("^PC_2_", "", pc2_cols)
    m[, sample_names, drop = FALSE]
  } else {
    matrix(0, nrow = nrow(wins), ncol = length(sample_names),
           dimnames = list(NULL, sample_names))
  }

  # FIX 38 (chat 8): previous code was
  #   win_starts = wins$start_bp %||% wins$mid_bp - 50000L,
  #   win_ends   = wins$end_bp   %||% wins$mid_bp + 50000L,
  # which parses as (wins$start_bp %||% wins$mid_bp) - 50000L because
  # %||% binds tighter than +/- in R. When wins$start_bp is present
  # (the common case) this shifted every window start back 50 kb from
  # its real value, inflating every window interval by 100 kb.
  #
  # The precomp (C01a STEP_C01a_precompute.R L849) ALWAYS writes
  # start_bp and end_bp as real window boundaries inherited from the
  # upstream per-window table `dt`. So the %||% fallback branch was
  # never supposed to fire — and the 50000 magic number was an
  # arbitrary guess (real window sizes come from C01a's multiscale
  # ladder: 20/40/80/120/160/200/240/320, not a fixed ±50 kb).
  #
  # Using the real columns directly:
  win_starts <- wins$start_bp
  win_ends   <- wins$end_bp

  list(
    pc1_mat = pc1_mat, pc2_mat = pc2_mat,
    sample_names = sample_names,
    n_windows = nrow(wins),
    win_starts = win_starts,
    win_ends   = win_ends,
    win_mids   = wins$mid_bp
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# k-means with quality metrics
# ─────────────────────────────────────────────────────────────────────────────

# Run k-means and compute silhouette score. Returns NULL if k-means fails.
kmeans_with_quality <- function(x, k = 3L, nstart = 25L, seed = 42L,
                                  compute_silhouette = TRUE) {
  set.seed(seed)
  km <- tryCatch(
    kmeans(x, centers = k, nstart = nstart, iter.max = 100),
    error = function(e) NULL
  )
  if (is.null(km)) return(NULL)

  result <- list(
    cluster = km$cluster,
    centers = km$centers,
    size = km$size,
    tot_withinss = km$tot.withinss,
    betweenss = km$betweenss
  )

  if (compute_silhouette && k >= 2 && length(unique(km$cluster)) >= 2) {
    result$silhouette_mean <- compute_silhouette_mean(x, km$cluster)
  } else {
    result$silhouette_mean <- NA_real_
  }
  result
}

# Simple silhouette mean without external cluster package dependency.
# For each point i: a = mean distance to own cluster, b = min mean distance
# to other clusters, s_i = (b-a)/max(a,b). Return mean over all points.
compute_silhouette_mean <- function(x, cluster) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n < 4) return(NA_real_)
  cl <- as.integer(cluster)
  unique_cl <- unique(cl)
  if (length(unique_cl) < 2) return(NA_real_)

  # Precompute distance matrix (may be slow for large n, but n is
  # bounded by sample count ≤ 226)
  dmat <- as.matrix(dist(x))

  s <- numeric(n)
  for (i in seq_len(n)) {
    own <- cl[i]
    own_idx <- which(cl == own & seq_len(n) != i)
    if (length(own_idx) == 0) { s[i] <- 0; next }
    a_i <- mean(dmat[i, own_idx])

    others <- setdiff(unique_cl, own)
    b_vals <- vapply(others, function(c2) {
      idx <- which(cl == c2)
      if (length(idx) == 0) return(Inf)
      mean(dmat[i, idx])
    }, numeric(1))
    b_i <- min(b_vals)

    s[i] <- if (max(a_i, b_i) > 0) (b_i - a_i) / max(a_i, b_i) else 0
  }
  mean(s)
}

# BIC-like gap between k=3 and k=2 (higher gap → k=3 is much better).
# Simple: compare tot.withinss normalized by k penalty.
compute_bic_gap <- function(x) {
  x <- as.matrix(x)
  km2 <- tryCatch(kmeans(x, centers = 2, nstart = 10), error = function(e) NULL)
  km3 <- tryCatch(kmeans(x, centers = 3, nstart = 10), error = function(e) NULL)
  if (is.null(km2) || is.null(km3)) return(NA_real_)
  n <- nrow(x)
  bic2 <- n * log(km2$tot.withinss / n) + 2 * log(n)
  bic3 <- n * log(km3$tot.withinss / n) + 3 * log(n)
  (bic2 - bic3) / abs(bic2)  # positive means k=3 fits better
}

# ─────────────────────────────────────────────────────────────────────────────
# Class label assignment
# ─────────────────────────────────────────────────────────────────────────────

# Given a k-means result on mean PC1 (or PC1+PC2), order clusters by mean
# PC1 and assign canonical labels HOM_REF (lowest), HET (middle), HOM_INV
# (highest). Returns a character vector aligned to names(x).
assign_class_labels <- function(km, x, k_use = 3L) {
  labels <- c("HOM_REF", "HET", "HOM_INV")[seq_len(k_use)]
  if (is.matrix(km$centers)) {
    # Sort by first column (PC1)
    order_by <- km$centers[, 1]
  } else {
    order_by <- km$centers
  }
  cluster_order <- order(order_by)
  cluster_to_class <- setNames(labels, as.character(cluster_order))
  class <- cluster_to_class[as.character(km$cluster)]
  names(class) <- names(x) %||% rownames(x) %||% seq_along(km$cluster)
  class
}

# ─────────────────────────────────────────────────────────────────────────────
# Loading optional data sources (fall back gracefully if absent)
# ─────────────────────────────────────────────────────────────────────────────

# Try to load sv_prior for a chromosome. Returns NULL if unavailable.
try_load_sv_prior <- function(chr, sv_prior_dir) {
  fl_path <- file.path(sv_prior_dir, paste0(chr, ".rds"))
  if (!file.exists(fl_path)) return(NULL)
  tryCatch(readRDS(fl_path), error = function(e) {
    message("[lib] sv_prior load failed for ", chr, ": ", conditionMessage(e))
    NULL
  })
}

# Load Clair3 phase info for a sample in a region. Returns data.table with
# columns POS, IS_PHASED, or empty data.table.
load_clair3_phase <- function(chr, sample_id, start_bp, end_bp, clair3_dir) {
  pf <- file.path(clair3_dir, chr, sample_id, "all_variants_with_phase.tsv")
  if (!file.exists(pf)) return(data.table())
  dt <- tryCatch(
    fread(pf, select = c(2, 41), header = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(data.table())
  setnames(dt, c("POS", "IS_PHASED"))
  dt[POS >= start_bp & POS <= end_bp]
}

# ─────────────────────────────────────────────────────────────────────────────
# Registry helpers (thin wrappers; heavy lifting in registry_loader.R)
# ─────────────────────────────────────────────────────────────────────────────

# Load reg — v10 ONLY. If the v10 registry can't be sourced, return NULL and
# let the caller fail loudly. The v9 sample_registry.R fallback was removed in
# the full-reg migration (chat B continuation, 2026-04-24). Rationale:
# * v9 → v10 migration is complete; no production scripts still depend on v9.
# * The fallback produced a list whose shape (legacy_v9 = TRUE) was never
#   actually branched on in downstream code, so the branch was dead weight.
# * Mixing v9 and v10 objects in the same session was a rare but real source
#   of silent key drops when write_block was NULL.
try_load_registry <- function() {
  for (p in c(
    "registries/api/R/registry_loader.R",
    "../registries/api/R/registry_loader.R",
    file.path(Sys.getenv("BASE", ""),
              "inversion-popgen-toolkit/registries/api/R/registry_loader.R")
  )) {
    if (file.exists(p)) {
      tryCatch({
        source(p)
        if (exists("load_registry", mode = "function")) {
          return(load_registry())
        }
      }, error = function(e) {
        message("[lib] registry source failed at ", p, ": ", conditionMessage(e))
      })
    }
  }
  NULL
}

# ─────────────────────────────────────────────────────────────────────────────
# Block write/read — REGISTRY-ONLY (full-reg migration, 2026-04-24)
# ─────────────────────────────────────────────────────────────────────────────
#
# Both write_block_safe and read_block_safe are thin wrappers around
# reg$evidence$write_block / reg$evidence$read_block. No fallback.
#
# Earlier migration stages retained an `outdir_fallback` parameter for
# backward-compat with call sites still passing it; after the call-site
# cleanup pass (2026-04-24), all 8 call sites have been updated to drop
# the argument, and the parameter has been removed from the signatures.
#
# If a new caller appears passing `outdir_fallback = ...`, R will raise
# "unused argument" — which is the right behavior. The fallback is gone;
# callers that think they need it should instead ensure registry_bridge.R
# is sourced.

write_block_safe <- function(reg, candidate_id, block_type, data,
                               source_script) {
  if (is.null(reg) || is.null(reg$evidence) ||
      is.null(reg$evidence$write_block)) {
    stop("[lib] write_block_safe: reg$evidence$write_block unavailable. ",
         "The v10 registry is now mandatory (full-reg migration). ",
         "Ensure utils/registry_bridge.R is sourced before calling this.")
  }
  Sys.setenv(CURRENT_SCRIPT = source_script)
  reg$evidence$write_block(
    candidate_id = candidate_id,
    block_type   = block_type,
    data         = data
  )
}

read_block_safe <- function(reg, candidate_id, block_type) {
  if (is.null(reg) || is.null(reg$evidence) ||
      is.null(reg$evidence$read_block)) {
    stop("[lib] read_block_safe: reg$evidence$read_block unavailable. ",
         "The v10 registry is now mandatory (full-reg migration). ",
         "Ensure utils/registry_bridge.R is sourced before calling this.")
  }
  reg$evidence$read_block(candidate_id, block_type)
}

# ─────────────────────────────────────────────────────────────────────────────
# compute_ancestry_div_hom_ref_vs_hom_inv — instant_q join for band_composition
#   fold (2026-04-24 chat C/D).
#
# For one candidate interval and decompose's per-sample class assignment,
# average the Engine B per-window Q vectors first per sample (across windows
# overlapping the interval) then per class (HOM_REF / HOM_INV), and return
# L1(mean_Q_HOM_REF, mean_Q_HOM_INV).
#
# High values → the two homozygous classes come from different ancestral
# subpopulations; candidate may be a family-driven false positive rather
# than a true inversion (see
#   _archive_superseded/bk_schemas_pre_canonical/
#   band_composition_folded_into_internal_dynamics_2026-04-24/README.md).
#
# Returns NA_real_ on ANY failure (missing file, no overlap, no samples
# in a needed class, parse trouble, etc). The key becomes NA in the
# internal_dynamics block, never breaks decompose.
#
# Arguments:
#   q_cache_dir       — directory containing Engine B per-sample local-Q
#                       caches. Typical layout: either
#                         <q_cache_dir>/<chrom>.local_Q_samples.tsv.gz
#                       or
#                         <q_cache_dir>/K<K>/<chrom>.local_Q_samples.tsv.gz
#                       or
#                         <q_cache_dir>/K<K>/<chrom>.<sample_set>.local_Q_samples.tsv.gz
#                       The function globs for all of these shapes and
#                       takes the first hit.
#   chrom             — chromosome name exactly as written in the filename.
#   start_bp, end_bp  — candidate interval in bp (integer).
#   class_assignment  — named character vector; names are sample_ids,
#                       values are one of "HOM_REF" / "HET" / "HOM_INV".
#                       (HET is ignored for this computation.)
#
# Returns: single numeric (L1 distance) or NA_real_.
#
# Implementation notes:
#   - Expects columns `start_bp`, `end_bp`, `sample_id`, and `Q1..QK`
#     (any K) in the cache file. See unified_ancestry/src/instant_q.cpp
#     write_samples() for the canonical schema.
#   - Overlap rule: window.start < end_bp AND window.end > start_bp
#     (mirrors the Python reference in
#     unified_ancestry/engines/nested_composition/internal_ancestry_composition.py
#     L204-206).
#   - Per-sample averaging is an unweighted mean over overlapping windows.
#     If a sample has no overlapping windows, it's dropped from its class.
#   - Per-class averaging is an unweighted mean over surviving samples.
#     If a class ends up empty (no samples with data), returns NA_real_.
# ─────────────────────────────────────────────────────────────────────────────
compute_ancestry_div_hom_ref_vs_hom_inv <- function(q_cache_dir, chrom,
                                                    start_bp, end_bp,
                                                    class_assignment,
                                                    verbose = FALSE) {
  if (is.null(q_cache_dir) || is.na(q_cache_dir) || q_cache_dir == "") {
    return(NA_real_)
  }
  if (!dir.exists(q_cache_dir)) {
    if (verbose) message("[ancestry_div] q_cache_dir not found: ", q_cache_dir)
    return(NA_real_)
  }

  # ── File discovery ────────────────────────────────────────────────────────
  # Try (in order): flat <chrom>.local_Q_samples.tsv.gz; K<K>/<chrom>.*; any
  # K<K>/<chrom>.<anything>.local_Q_samples.tsv.gz. Glob all K* subdirs.
  candidates <- c(
    file.path(q_cache_dir, paste0(chrom, ".local_Q_samples.tsv.gz")),
    file.path(q_cache_dir, paste0(chrom, ".local_Q_samples.tsv")),
    Sys.glob(file.path(q_cache_dir, "K*",
                        paste0(chrom, ".local_Q_samples.tsv.gz"))),
    Sys.glob(file.path(q_cache_dir, "K*",
                        paste0(chrom, ".local_Q_samples.tsv"))),
    Sys.glob(file.path(q_cache_dir, "K*",
                        paste0(chrom, ".*.local_Q_samples.tsv.gz"))),
    Sys.glob(file.path(q_cache_dir, "K*",
                        paste0(chrom, ".*.local_Q_samples.tsv")))
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) {
    if (verbose) message("[ancestry_div] no local_Q file for ", chrom,
                         " under ", q_cache_dir)
    return(NA_real_)
  }
  q_path <- candidates[1]
  if (length(candidates) > 1 && verbose) {
    message("[ancestry_div] multiple local_Q files for ", chrom,
            "; using ", q_path)
  }

  # ── Read ──────────────────────────────────────────────────────────────────
  dt <- tryCatch(
    data.table::fread(q_path, sep = "\t", header = TRUE,
                      showProgress = FALSE),
    error = function(e) {
      message("[ancestry_div] fread failed for ", q_path, ": ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NA_real_)

  # ── Schema validation ─────────────────────────────────────────────────────
  needed <- c("start_bp", "end_bp", "sample_id")
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols) > 0) {
    message("[ancestry_div] ", q_path, " missing columns: ",
            paste(missing_cols, collapse = ", "))
    return(NA_real_)
  }
  q_cols <- grep("^Q[0-9]+$", names(dt), value = TRUE)
  if (length(q_cols) < 2) {
    message("[ancestry_div] ", q_path, " has fewer than 2 Q columns — ",
            "expected Q1..QK; found: ", paste(q_cols, collapse = ", "))
    return(NA_real_)
  }
  # Sort Q columns numerically so Q10 comes after Q9, not between Q1 and Q2
  q_cols <- q_cols[order(as.integer(sub("^Q", "", q_cols)))]

  # ── Filter to candidate interval ──────────────────────────────────────────
  # Python overlap rule: window.start < parent.end && window.end > parent.start.
  # data.table uses non-standard evaluation on column names; rename the
  # function args first so the i-expression is unambiguous.
  cand_start <- start_bp
  cand_end   <- end_bp
  over <- dt[start_bp < cand_end & end_bp > cand_start]
  if (nrow(over) == 0) {
    if (verbose) message("[ancestry_div] no windows overlap ", chrom, ":",
                         cand_start, "-", cand_end, " in ", q_path)
    return(NA_real_)
  }

  # ── Filter to samples in HOM_REF ∪ HOM_INV ────────────────────────────────
  # (We don't need HET for this divergence.)
  targets <- names(class_assignment)[class_assignment %in%
                                       c("HOM_REF", "HOM_INV")]
  if (length(targets) == 0) return(NA_real_)
  over <- over[sample_id %in% targets]
  if (nrow(over) == 0) {
    if (verbose) message("[ancestry_div] no target samples have overlapping ",
                         "windows for ", chrom, ":", cand_start, "-", cand_end)
    return(NA_real_)
  }

  # ── Per-sample mean Q (across overlapping windows) ────────────────────────
  # data.table's lapply(.SD, mean) is the right tool here.
  per_sample <- over[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
                     by = sample_id, .SDcols = q_cols]

  # Tag with class. Do the named-vector lookup in base R first, then assign,
  # to avoid any NSE confusion between the `class_assignment` argument (a
  # named character vector) and a potential column of the same name.
  cls_for_rows <- unname(class_assignment[per_sample$sample_id])
  per_sample[, class_label := cls_for_rows]

  # ── Per-class mean Q ──────────────────────────────────────────────────────
  ref_q <- per_sample[class_label == "HOM_REF",
                       lapply(.SD, mean, na.rm = TRUE), .SDcols = q_cols]
  inv_q <- per_sample[class_label == "HOM_INV",
                       lapply(.SD, mean, na.rm = TRUE), .SDcols = q_cols]
  if (nrow(ref_q) == 0 || nrow(inv_q) == 0) return(NA_real_)

  ref_vec <- unlist(ref_q, use.names = FALSE)
  inv_vec <- unlist(inv_q, use.names = FALSE)
  if (length(ref_vec) != length(inv_vec) ||
      any(is.na(ref_vec)) || any(is.na(inv_vec))) {
    return(NA_real_)
  }

  # ── L1 distance ──────────────────────────────────────────────────────────
  # Sum of absolute differences across K components. Range [0, 2] for
  # probability simplices; 0 = identical ancestry composition, 2 = fully
  # disjoint.
  round(sum(abs(ref_vec - inv_vec)), 4)
}

message("[lib_decompose_helpers] loaded (full-reg mode; no JSON fallback)")
