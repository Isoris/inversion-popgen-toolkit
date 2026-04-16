#!/usr/bin/env Rscript

# =============================================================================
# STEP10c_window_sample_belonging.R
#
# THE CORE FIX: Compare windows by SAMPLE BELONGING, not window shape.
#
# Two windows can both look like clean 3-band inversions in PCA/MDS, but the
# actual samples occupying left/middle/right may be different. Window-level
# geometry alone cannot detect this. We must compare windows by:
#   - which samples are in which band
#   - whether the same samples keep grouping together across windows
#
# This script sits between STEP10 (candidate nomination) and STEP12 (regional PCA).
# It returns to marker-level dosage from STEP08 (bypassing STEP09 compression)
# and builds per-window sample-belonging vectors.
#
# PRODUCES FIVE PARALLEL MDS + HEATMAP LAYERS:
#
#   Layer A  — SNP-based window MDS (existing STEP10 logic, reproduced here
#              within each candidate for direct comparison)
#   Layer B  — Sample-belonging MDS (windows compared by which samples are in
#              which band — coarse band assignment Hamming/Jaccard)
#   Layer B2 — RAW-VECTOR INDUCED DISTANCE MATRIX COMPARISON (THE REAL TEST)
#              For each window, compute full sample×sample distance matrix from
#              the ordered per-sample marker vectors (~100 markers × 226 samples).
#              Two windows are "truly the same system" only if they produce
#              correlated sample×sample distance structures from raw profiles.
#              Three sub-layers: minor (Manhattan), major (Manhattan), 012 (Hamming).
#   Layer C  — Scalar-fingerprint MDS (lightweight: mean dosage / het fraction
#              per sample per window — fast but misses within-band substructure)
#   Layer D  — (Future) Signature-distance MDS (when Clair3 available)
#
# For each layer: MDS plot + heatmap + hclust dendrogram, all with consistent
# coloring so they can be compared side by side.
#
# DESIGN PRINCIPLES:
#   - Operates on raw dosage from STEP08, not STEP09 PCA summaries
#   - Maintains three encodings (minor, major, 012) as first-class objects
#   - Clair3 is optional: Layer D is a stub that activates when data exists
#   - All intermediate objects are saved for downstream use
#
# Usage:
#   Rscript STEP10c_window_sample_belonging.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

has_ggplot <- suppressWarnings(require(ggplot2, quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

# ── Load candidate table ─────────────────────────────────────────────────────
# Prefer snake-derived candidates (primary method) over 500kb control candidates.
# SNAKE_CANDIDATE_TABLE is set in config when STEP10e_v3 has been run.
SNAKE_CANDIDATE_TABLE <- file.path(INV_ROOT, "06_mds_candidates",
                                   "snake_regions_multiscale",
                                   "snake_candidate_regions.tsv.gz")

if (file.exists(SNAKE_CANDIDATE_TABLE)) {
  cand <- fread(SNAKE_CANDIDATE_TABLE)
  message("[INFO] Using SNAKE candidate table: ", SNAKE_CANDIDATE_TABLE)
} else if (file.exists(CANDIDATE_TABLE)) {
  cand <- fread(CANDIDATE_TABLE)
  message("[INFO] Using CONTROL candidate table (500kb): ", CANDIDATE_TABLE)
} else {
  stop("No candidate table found")
}
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

# Sort candidates by chromosome + start position for consistent ordering
cand <- cand[order(chrom, start_bp)]

message("[INFO] Processing ", nrow(cand), " candidates")

# Load STEP09 for window coordinates
step09_rds <- list.files(STEP09_DIR, pattern = "\\.window_pca\\.rds$",
                         full.names = TRUE, recursive = TRUE)
step09_by_chr <- list()
for (f in step09_rds) {
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.null(obj)) step09_by_chr[[as.character(obj$chrom)]] <- obj
}

# Clair3 check
CLAIR3_ROOT <- file.path(PROJECT_ROOT, "clair3_indel_discovery", "postprocess_results")
clair3_available <- dir.exists(CLAIR3_ROOT)

# ═══════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════

discretize_012 <- function(x) {
  out <- rep(NA_integer_, length(x))
  out[!is.na(x) & x < 0.5] <- 0L
  out[!is.na(x) & x >= 0.5 & x <= 1.5] <- 1L
  out[!is.na(x) & x > 1.5] <- 2L
  out
}

# Assign samples to bands (0/1/2) using adaptive k on PC1
assign_bands <- function(pc1_vec, max_k = 3) {
  n <- length(pc1_vec)
  if (n < 6) return(rep(NA_integer_, n))

  km2 <- tryCatch(kmeans(matrix(pc1_vec, ncol = 1), 2, nstart = 25), error = function(e) NULL)
  km3 <- tryCatch(kmeans(matrix(pc1_vec, ncol = 1), 3, nstart = 25), error = function(e) NULL)

  best_k <- 3L; km <- km3
  if (!is.null(km3) && !is.null(km2)) {
    if (km3$tot.withinss / km2$tot.withinss > 0.80) { best_k <- 2L; km <- km2 }
  } else if (!is.null(km2)) {
    best_k <- 2L; km <- km2
  }
  if (is.null(km)) return(rep(NA_integer_, n))

  # Order groups by PC1 mean → 1=left, 2=middle, 3=right (or 1=left, 2=right if k=2)
  gm <- tapply(pc1_vec, km$cluster, mean)
  go <- order(gm)
  match(km$cluster, go)
}

# Hamming distance between two integer vectors
hamming_frac <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) == 0) return(NA_real_)
  sum(a[ok] != b[ok]) / sum(ok)
}

# Jaccard distance for sample sets: 1 - |A∩B| / |A∪B|
jaccard_dist <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  union <- length(union(set_a, set_b))
  if (union == 0) return(NA_real_)
  1 - inter / union
}

# Manhattan on continuous vectors, normalized
manhattan_norm <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  n <- sum(ok)
  if (n == 0) return(NA_real_)
  sum(abs(a[ok] - b[ok])) / n
}

# Safe MDS from distance matrix
safe_mds <- function(dmat, k = 4) {
  dmat[is.na(dmat)] <- median(dmat, na.rm = TRUE)
  diag(dmat) <- 0
  dmat <- (dmat + t(dmat)) / 2  # symmetrize
  tryCatch(cmdscale(as.dist(dmat), k = min(k, nrow(dmat) - 1)), error = function(e) NULL)
}

# ── Plotting helpers ─────────────────────────────────────────────────────────
make_heatmap_png <- function(mat, filepath, title, labels = NULL) {
  if (!has_ggplot) return(invisible())
  n <- nrow(mat)
  if (n < 2) return(invisible())

  # Melt for ggplot
  dt <- data.table(
    row = rep(seq_len(n), n),
    col = rep(seq_len(n), each = n),
    value = as.vector(mat)
  )
  if (!is.null(labels)) {
    dt[, row_lab := labels[row]]
    dt[, col_lab := labels[col]]
  }

  p <- ggplot(dt, aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno", na.value = "grey80") +
    scale_y_reverse() +
    theme_minimal(base_size = 9) +
    labs(title = title, x = "Window", y = "Window", fill = "Distance") +
    coord_equal()

  ggsave(filepath, p, width = 8, height = 7, dpi = 300)
}

make_mds_plot <- function(mds_coords, color_vec, filepath, title, color_label = "Group") {
  if (!has_ggplot || is.null(mds_coords)) return(invisible())
  dt <- data.table(
    MDS1 = mds_coords[, 1],
    MDS2 = if (ncol(mds_coords) >= 2) mds_coords[, 2] else 0,
    color = color_vec
  )
  p <- ggplot(dt, aes(x = MDS1, y = MDS2, color = factor(color))) +
    geom_point(size = 3, alpha = 0.8) +
    theme_bw(base_size = 11) +
    labs(title = title, color = color_label)

  ggsave(filepath, p, width = 7, height = 5.5, dpi = 300)
}

make_dendrogram_png <- function(dmat, filepath, title, labels = NULL) {
  if (nrow(dmat) < 3) return(invisible())
  dmat[is.na(dmat)] <- median(dmat, na.rm = TRUE)
  diag(dmat) <- 0
  dmat <- (dmat + t(dmat)) / 2
  hc <- hclust(as.dist(dmat), method = "average")
  if (!is.null(labels)) hc$labels <- labels

  png(filepath, width = max(800, nrow(dmat) * 15), height = 400, res = 300)
  par(mar = c(5, 4, 3, 1), cex = 0.7)
  plot(hc, main = title, xlab = "", sub = "")
  dev.off()
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN LOOP — Per candidate
# ═══════════════════════════════════════════════════════════════════════════

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id; chr <- row$chrom
  c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)

  # Build sortable prefix: chr_candN_startMb (sorts naturally by chr then position)
  cand_prefix <- paste0(chr, "_cand", cid, "_", sprintf("%.2fMb", c_start / 1e6))

  cand_dir <- file.path(FOLLOWUP_DIR, cand_prefix)
  plot_dir <- file.path(PLOTS_DIR, cand_prefix, "window_belonging")
  ensure_dir(cand_dir); ensure_dir(plot_dir)

  message("\n[INFO] ════════════════════════════════════════════════════════")
  message("[INFO] Candidate ", cid, " — window sample-belonging analysis")

  # ── Load raw dosage (bypassing STEP09 compression) ────────────────────
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[SKIP] Missing dosage for ", chr); next
  }

  dos <- fread(dos_file); sites <- fread(sites_file)
  sample_cols <- setdiff(names(dos), "marker")
  n_samples <- length(sample_cols)

  # Sample name mapping if Ind-style
  samp_ind <- tryCatch(fread(SAMPLE_IND_FILE, header = FALSE)[[1]], error = function(e) NULL)
  if (!is.null(samp_ind) && all(grepl("^Ind", sample_cols)) && length(samp_ind) == n_samples) {
    samp_ind <- samp_ind[nchar(samp_ind) > 0]
    if (length(samp_ind) == n_samples) {
      setnames(dos, old = sample_cols, new = samp_ind)
      sample_cols <- samp_ind
    }
  }

  # Get windows overlapping candidate
  s09 <- step09_by_chr[[chr]]
  if (is.null(s09)) { message("[SKIP] No STEP09 data for ", chr); next }
  wmeta <- s09$window_meta
  woverlap <- wmeta[end_bp >= c_start & start_bp <= c_end]
  if (nrow(woverlap) < 3) { message("[SKIP] <3 windows overlap"); next }
  n_windows <- nrow(woverlap)

  message("[INFO] Windows: ", n_windows, "  Samples: ", n_samples)

  # ═════════════════════════════════════════════════════════════════════
  # PER-WINDOW: assign samples to bands using raw dosage
  # ═════════════════════════════════════════════════════════════════════

  # Store per-window: band assignments, raw vectors, PCA
  win_band_assignments <- matrix(NA_integer_, n_samples, n_windows)  # samples × windows
  win_pc1 <- matrix(NA_real_, n_samples, n_windows)
  win_labels <- character(n_windows)
  win_k <- integer(n_windows)
  win_n_markers <- integer(n_windows)

  # Per-window minor/major dosage mean vectors for scalar-fingerprint layer
  # For each window, store the sample-level mean dosage
  win_minor_means <- matrix(NA_real_, n_samples, n_windows)
  win_major_means <- matrix(NA_real_, n_samples, n_windows)

  # Per-window FULL raw marker matrices for Layer B2
  # Stored as a list of matrices (markers × samples) — one per window
  win_raw_minor <- vector("list", n_windows)   # continuous minor dosage
  win_raw_major <- vector("list", n_windows)   # continuous major dosage
  win_raw_012   <- vector("list", n_windows)   # discrete 0/1/2

  for (wi in seq_len(n_windows)) {
    wrow <- woverlap[wi]
    w_keep <- which(sites$pos >= wrow$start_bp & sites$pos <= wrow$end_bp)
    if (length(w_keep) < 10) next

    Xw <- as.matrix(dos[w_keep, ..sample_cols])
    storage.mode(Xw) <- "double"
    win_n_markers[wi] <- nrow(Xw)

    # PCA on raw dosage
    pc <- tryCatch({
      p <- prcomp(t(Xw), center = TRUE, scale. = FALSE)
      p$x[, 1]
    }, error = function(e) NULL)
    if (is.null(pc)) next

    win_pc1[, wi] <- pc

    # Assign bands
    bands <- assign_bands(pc)
    win_band_assignments[, wi] <- bands
    win_k[wi] <- max(bands, na.rm = TRUE)

    # Per-sample mean dosage in this window
    win_minor_means[, wi] <- colMeans(Xw, na.rm = TRUE)
    win_major_means[, wi] <- 2 - colMeans(Xw, na.rm = TRUE)

    # Full raw marker matrices for Layer B2
    win_raw_minor[[wi]] <- Xw                        # markers × samples
    win_raw_major[[wi]] <- 2.0 - Xw                  # markers × samples
    # v7.3 fix: use vectorized discretize_012() instead of slow element-wise apply
    X012w <- matrix(discretize_012(as.vector(Xw)), nrow = nrow(Xw), ncol = ncol(Xw))
    win_raw_012[[wi]] <- X012w                       # markers × samples

    win_labels[wi] <- paste0("w", wi, ":",
                              format(wrow$start_bp / 1e6, digits = 3), "Mb")
  }

  # Remove windows with no data
  valid_wins <- which(colSums(!is.na(win_band_assignments)) > 0)
  if (length(valid_wins) < 3) { message("[SKIP] <3 valid windows"); next }

  bands <- win_band_assignments[, valid_wins, drop = FALSE]
  labels <- win_labels[valid_wins]
  nw <- length(valid_wins)

  message("[INFO] Valid windows: ", nw)

  # ═════════════════════════════════════════════════════════════════════
  # LAYER A — SNP-BASED WINDOW DISTANCE (Frobenius on covariance)
  # This reproduces STEP10 logic within the candidate for comparison.
  # ═════════════════════════════════════════════════════════════════════

  message("[INFO] Layer A: SNP-based window distance")

  # Use per-window PC1 vectors as proxies for covariance structure
  # Frobenius approximation: ||cov_i - cov_j|| ≈ correlation of PC1 loadings
  dist_A <- matrix(NA_real_, nw, nw)
  for (i in seq_len(nw)) {
    dist_A[i, i] <- 0
    if (i < nw) for (j in (i + 1):nw) {
      # Correlation-based distance on PC1 vectors
      r <- suppressWarnings(cor(win_pc1[, valid_wins[i]], win_pc1[, valid_wins[j]],
                                 use = "pairwise.complete.obs"))
      d <- if (is.finite(r)) 1 - abs(r) else 1
      dist_A[i, j] <- d; dist_A[j, i] <- d
    }
  }

  mds_A <- safe_mds(dist_A)

  # ═════════════════════════════════════════════════════════════════════
  # LAYER B — SAMPLE-BELONGING DISTANCE (THIS IS THE CORE FIX)
  # Windows are similar iff the SAME SAMPLES occupy corresponding bands.
  # ═════════════════════════════════════════════════════════════════════

  message("[INFO] Layer B: Sample-belonging distance")

  # For each pair of windows: compare band assignment vectors
  # Use Hamming distance on integer band vectors (0-1-2-3)
  # But also: Jaccard per band (do left-band samples overlap?)
  dist_B_hamming <- matrix(NA_real_, nw, nw)
  dist_B_jaccard_mean <- matrix(NA_real_, nw, nw)

  for (i in seq_len(nw)) {
    dist_B_hamming[i, i] <- 0
    dist_B_jaccard_mean[i, i] <- 0
    if (i < nw) for (j in (i + 1):nw) {
      bi <- bands[, i]; bj <- bands[, j]

      # Hamming: fraction of samples assigned to different bands
      hd <- hamming_frac(bi, bj)

      # Also try flipped assignment (in case PC1 sign flipped)
      # If window j has bands 1,2,3 but PC1 was flipped, its "band 1" = other window's "band 3"
      bj_flip <- bj
      bj_flip[bj == 1L] <- 99L
      bj_flip[bj == max(bj, na.rm = TRUE)] <- 1L
      bj_flip[bj_flip == 99L] <- max(bj, na.rm = TRUE)
      hd_flip <- hamming_frac(bi, bj_flip)

      # Take the better (lower distance) of original and flipped
      hd_best <- min(hd, hd_flip, na.rm = TRUE)

      # Jaccard per band: average Jaccard across bands
      max_band <- max(max(bi, na.rm = TRUE), max(bj, na.rm = TRUE))
      jaccards <- numeric(max_band)
      for (b in seq_len(max_band)) {
        si <- which(bi == b); sj <- which(bj == b)
        jaccards[b] <- jaccard_dist(si, sj)
      }
      jd <- mean(jaccards, na.rm = TRUE)

      # Also Jaccard with flip
      jaccards_flip <- numeric(max_band)
      for (b in seq_len(max_band)) {
        si <- which(bi == b); sj <- which(bj_flip == b)
        jaccards_flip[b] <- jaccard_dist(si, sj)
      }
      jd_flip <- mean(jaccards_flip, na.rm = TRUE)

      jd_best <- min(jd, jd_flip, na.rm = TRUE)

      dist_B_hamming[i, j] <- hd_best; dist_B_hamming[j, i] <- hd_best
      dist_B_jaccard_mean[i, j] <- jd_best; dist_B_jaccard_mean[j, i] <- jd_best
    }
  }

  mds_B <- safe_mds(dist_B_hamming)

  # ═════════════════════════════════════════════════════════════════════
  # LAYER B2 — RAW-VECTOR INDUCED DISTANCE MATRIX COMPARISON
  #
  # THE REAL TEST: Two windows are "truly the same system" only if they
  # produce correlated sample×sample distance matrices from their full
  # ordered marker vectors.
  #
  # For each window w:
  #   1. Extract full per-sample marker vector matrix (markers × samples)
  #   2. Compute pairwise sample×sample distance matrix S_w (n_samples × n_samples)
  #   3. Flatten upper triangle → distance profile vector v_w
  # For each window pair (w_i, w_j):
  #   dist(w_i, w_j) = 1 - |Spearman_cor(v_wi, v_wj)|
  #
  # Three sub-layers: minor (Manhattan), major (Manhattan), 012 (Hamming)
  # This catches within-band subgroup differences, marker-order patterns,
  # and local profile variation that scalar fingerprints miss.
  # ═════════════════════════════════════════════════════════════════════

  message("[INFO] Layer B2: Raw-vector induced distance matrix comparison (3 encodings)")

  # Helper: compute sample×sample distance matrix from raw marker matrix
  # Returns flattened upper triangle as a numeric vector
  compute_sample_dist_profile <- function(Xmat, dist_type = "manhattan") {
    # Xmat: markers × samples
    ns <- ncol(Xmat)
    # Build sample×sample distance matrix
    sdmat <- matrix(0, ns, ns)
    for (a in seq_len(ns - 1L)) {
      for (b in (a + 1L):ns) {
        ok <- !is.na(Xmat[, a]) & !is.na(Xmat[, b])
        nok <- sum(ok)
        if (nok == 0) {
          d <- NA_real_
        } else if (dist_type == "hamming") {
          d <- sum(Xmat[ok, a] != Xmat[ok, b]) / nok
        } else {
          # manhattan, normalized
          d <- sum(abs(Xmat[ok, a] - Xmat[ok, b])) / nok
        }
        sdmat[a, b] <- d; sdmat[b, a] <- d
      }
    }
    # Flatten upper triangle (consistent ordering)
    sdmat[upper.tri(sdmat)]
  }

  # Compute distance profiles for all valid windows, for each encoding
  B2_profiles_minor <- vector("list", nw)
  B2_profiles_major <- vector("list", nw)
  B2_profiles_012   <- vector("list", nw)

  for (vi in seq_len(nw)) {
    wi <- valid_wins[vi]
    if (!is.null(win_raw_minor[[wi]])) {
      B2_profiles_minor[[vi]] <- compute_sample_dist_profile(win_raw_minor[[wi]], "manhattan")
      B2_profiles_major[[vi]] <- compute_sample_dist_profile(win_raw_major[[wi]], "manhattan")
    }
    if (!is.null(win_raw_012[[wi]])) {
      B2_profiles_012[[vi]] <- compute_sample_dist_profile(win_raw_012[[wi]], "hamming")
    }
  }

  # Compare windows by Spearman correlation of their induced distance profiles
  compute_B2_layer <- function(profiles, layer_label) {
    dmat <- matrix(NA_real_, nw, nw)
    for (i in seq_len(nw)) {
      dmat[i, i] <- 0
      if (i < nw) for (j in (i + 1):nw) {
        pi <- profiles[[i]]; pj <- profiles[[j]]
        if (is.null(pi) || is.null(pj)) {
          dmat[i, j] <- NA_real_; dmat[j, i] <- NA_real_
          next
        }
        ok <- !is.na(pi) & !is.na(pj)
        if (sum(ok) < 10) {
          dmat[i, j] <- NA_real_; dmat[j, i] <- NA_real_
          next
        }
        r <- suppressWarnings(cor(pi[ok], pj[ok], method = "spearman"))
        d <- if (is.finite(r)) 1 - abs(r) else 1
        dmat[i, j] <- d; dmat[j, i] <- d
      }
    }
    list(dist = dmat, mds = safe_mds(dmat), label = layer_label)
  }

  layer_B2_minor <- compute_B2_layer(B2_profiles_minor, "B2_rawvec_minor")
  layer_B2_major <- compute_B2_layer(B2_profiles_major, "B2_rawvec_major")
  layer_B2_012   <- compute_B2_layer(B2_profiles_012,   "B2_rawvec_012")

  message("[INFO] Layer B2 complete: ", nw, " windows × 3 encodings")

  # ═════════════════════════════════════════════════════════════════════
  # LAYER C — SCALAR-FINGERPRINT MDS (lightweight complement)
  # Windows compared by per-sample scalar summaries (mean dosage, het
  # fraction). Fast but misses within-band substructure and marker order.
  # Retained as a quick sanity check alongside the richer Layer B2.
  # ═════════════════════════════════════════════════════════════════════

  message("[INFO] Layer C: Scalar-fingerprint MDS (3 encodings)")

  # For each window, build a compact scalar fingerprint per sample
  # (mean dosage or het fraction), then compare windows by Spearman
  # correlation of these fingerprint vectors across samples

  compute_layer_C <- function(win_means_mat, label) {
    # win_means_mat: samples × windows
    vm <- win_means_mat[, valid_wins, drop = FALSE]
    dmat <- matrix(NA_real_, nw, nw)
    for (i in seq_len(nw)) {
      dmat[i, i] <- 0
      if (i < nw) for (j in (i + 1):nw) {
        # Compare: do these two windows rank samples the same way?
        r <- suppressWarnings(cor(vm[, i], vm[, j],
                                   method = "spearman", use = "pairwise.complete.obs"))
        d <- if (is.finite(r)) 1 - abs(r) else 1
        dmat[i, j] <- d; dmat[j, i] <- d
      }
    }
    list(dist = dmat, mds = safe_mds(dmat), label = label)
  }

  layer_C_minor <- compute_layer_C(win_minor_means, "minor")
  layer_C_major <- compute_layer_C(win_major_means, "major")

  # For 012: use fraction of het calls per sample per window as fingerprint
  win_het_frac <- matrix(NA_real_, n_samples, n_windows)
  for (wi in seq_len(n_windows)) {
    wrow <- woverlap[wi]
    w_keep <- which(sites$pos >= wrow$start_bp & sites$pos <= wrow$end_bp)
    if (length(w_keep) < 10) next
    Xw <- as.matrix(dos[w_keep, ..sample_cols])
    # v7.3 fix: use vectorized discretize_012() instead of slow element-wise apply
    X012 <- matrix(discretize_012(as.vector(Xw)), nrow = nrow(Xw), ncol = ncol(Xw))
    win_het_frac[, wi] <- colMeans(X012 == 1L, na.rm = TRUE)
  }
  layer_C_012 <- compute_layer_C(win_het_frac, "012_het_frac")

  # ═════════════════════════════════════════════════════════════════════
  # LAYER D — SIGNATURE DISTANCE (stub for Clair3)
  # ═════════════════════════════════════════════════════════════════════

  dist_D <- NULL
  mds_D <- NULL
  if (clair3_available) {
    message("[INFO] Layer D: Clair3 signature distance (stub — will activate with data)")
    # Future: load per-window Clair3 signature matrices, compute sample
    # belonging from indel/phase signatures, compare windows by that
  } else {
    message("[INFO] Layer D: Clair3 not available, skipping")
  }

  # ═════════════════════════════════════════════════════════════════════
  # HCLUST ON LAYER B (the core belonging layer)
  # ═════════════════════════════════════════════════════════════════════

  message("[INFO] Hierarchical clustering on sample-belonging distances")

  db <- dist_B_hamming
  db[is.na(db)] <- median(db, na.rm = TRUE)
  diag(db) <- 0
  hc_B <- hclust(as.dist(db), method = "average")

  # Cut at several heights for diagnostic coloring
  hc_k3 <- if (nw >= 3) cutree(hc_B, k = min(3, nw)) else rep(1L, nw)
  hc_k5 <- if (nw >= 5) cutree(hc_B, k = min(5, nw)) else hc_k3
  hc_h05 <- cutree(hc_B, h = 0.5)  # windows with <0.5 Hamming are "same system"

  # ═════════════════════════════════════════════════════════════════════
  # WRITE OUTPUTS
  # ═════════════════════════════════════════════════════════════════════

  message("[INFO] Writing outputs")

  # Band assignment table (samples × windows)
  band_dt <- data.table(sample = sample_cols)
  for (wi in seq_len(nw)) {
    set(band_dt, j = labels[wi], value = bands[, wi])
  }
  band_dt[, candidate_id := cid]
  fwrite(band_dt, file.path(cand_dir, paste0(cand_prefix, "_window_band_assignments.tsv.gz")), sep = "\t")

  # Distance matrices
  save_dist <- function(dmat, suffix) {
    dt <- as.data.table(dmat)
    names(dt) <- labels
    dt[, window := labels]
    setcolorder(dt, c("window", labels))
    dt[, candidate_id := cid]
    fwrite(dt, file.path(cand_dir, paste0(cand_prefix, "_", suffix)), sep = "\t")
  }

  save_dist(dist_A, "dist_A_snp.tsv")
  save_dist(dist_B_hamming, "dist_B_belonging_hamming.tsv")
  save_dist(dist_B_jaccard_mean, "dist_B_belonging_jaccard.tsv")
  save_dist(layer_B2_minor$dist, "dist_B2_rawvec_minor.tsv")
  save_dist(layer_B2_major$dist, "dist_B2_rawvec_major.tsv")
  save_dist(layer_B2_012$dist, "dist_B2_rawvec_012.tsv")
  save_dist(layer_C_minor$dist, "dist_C_scalar_minor.tsv")
  save_dist(layer_C_major$dist, "dist_C_scalar_major.tsv")
  save_dist(layer_C_012$dist, "dist_C_scalar_012.tsv")

  # Window clustering summary
  cluster_dt <- data.table(
    window_idx = seq_len(nw),
    window_label = labels,
    hclust_k3 = hc_k3,
    hclust_k5 = hc_k5,
    hclust_h05 = hc_h05,
    n_markers = win_n_markers[valid_wins],
    best_k = win_k[valid_wins],
    candidate_id = cid
  )
  fwrite(cluster_dt, file.path(cand_dir, paste0(cand_prefix, "_belonging_clusters.tsv")), sep = "\t")

  # Cross-layer agreement table
  agreement_dt <- data.table(
    layer_pair = c("A_vs_B", "A_vs_B2minor", "B_vs_B2minor",
                   "B2minor_vs_B2major", "B2minor_vs_B2_012",
                   "B2minor_vs_Cminor", "A_vs_Cminor", "B_vs_Cminor",
                   "Cminor_vs_Cmajor"),
    mantel_corr = c(
      suppressWarnings(cor(as.vector(dist_A), as.vector(dist_B_hamming), use = "complete.obs")),
      suppressWarnings(cor(as.vector(dist_A), as.vector(layer_B2_minor$dist), use = "complete.obs")),
      suppressWarnings(cor(as.vector(dist_B_hamming), as.vector(layer_B2_minor$dist), use = "complete.obs")),
      suppressWarnings(cor(as.vector(layer_B2_minor$dist), as.vector(layer_B2_major$dist), use = "complete.obs")),
      suppressWarnings(cor(as.vector(layer_B2_minor$dist), as.vector(layer_B2_012$dist), use = "complete.obs")),
      suppressWarnings(cor(as.vector(layer_B2_minor$dist), as.vector(layer_C_minor$dist), use = "complete.obs")),
      suppressWarnings(cor(as.vector(dist_A), as.vector(layer_C_minor$dist), use = "complete.obs")),
      suppressWarnings(cor(as.vector(dist_B_hamming), as.vector(layer_C_minor$dist), use = "complete.obs")),
      suppressWarnings(cor(as.vector(layer_C_minor$dist), as.vector(layer_C_major$dist), use = "complete.obs"))
    ),
    candidate_id = cid
  )
  fwrite(agreement_dt, file.path(cand_dir, paste0(cand_prefix, "_layer_agreement.tsv")), sep = "\t")

  # ═════════════════════════════════════════════════════════════════════
  # DIAGNOSTIC PLOTS
  # ═════════════════════════════════════════════════════════════════════

  if (has_ggplot) {
    message("[INFO] Generating diagnostic plots")

    # Heatmaps
    make_heatmap_png(dist_A, file.path(plot_dir, paste0(cand_prefix, "_heatmap_A_snp.png")),
                     paste0(cand_prefix, " — Layer A: SNP-based window distance"), labels)
    make_heatmap_png(dist_B_hamming, file.path(plot_dir, paste0(cand_prefix, "_heatmap_B_hamming.png")),
                     paste0(cand_prefix, " — Layer B: Sample-belonging Hamming"), labels)
    make_heatmap_png(dist_B_jaccard_mean, file.path(plot_dir, paste0(cand_prefix, "_heatmap_B_jaccard.png")),
                     paste0(cand_prefix, " — Layer B: Sample-belonging Jaccard"), labels)
    make_heatmap_png(layer_B2_minor$dist, file.path(plot_dir, paste0(cand_prefix, "_heatmap_B2_minor.png")),
                     paste0(cand_prefix, " — Layer B2: Raw-vector (minor)"), labels)
    make_heatmap_png(layer_B2_major$dist, file.path(plot_dir, paste0(cand_prefix, "_heatmap_B2_major.png")),
                     paste0(cand_prefix, " — Layer B2: Raw-vector (major)"), labels)
    make_heatmap_png(layer_B2_012$dist, file.path(plot_dir, paste0(cand_prefix, "_heatmap_B2_012.png")),
                     paste0(cand_prefix, " — Layer B2: Raw-vector (012)"), labels)
    make_heatmap_png(layer_C_minor$dist, file.path(plot_dir, paste0(cand_prefix, "_heatmap_C_minor.png")),
                     paste0(cand_prefix, " — Layer C: Scalar fingerprint (minor)"), labels)
    make_heatmap_png(layer_C_major$dist, file.path(plot_dir, paste0(cand_prefix, "_heatmap_C_major.png")),
                     paste0(cand_prefix, " — Layer C: Scalar fingerprint (major)"), labels)

    # MDS plots — colored by hclust groups from Layer B
    if (!is.null(mds_A)) {
      make_mds_plot(mds_A, hc_h05,
                    file.path(plot_dir, paste0(cand_prefix, "_mds_A_snp.png")),
                    paste0(cand_prefix, " — MDS Layer A (SNP)"),
                    "Belonging\ncluster")
    }
    if (!is.null(mds_B)) {
      make_mds_plot(mds_B, hc_h05,
                    file.path(plot_dir, paste0(cand_prefix, "_mds_B_belonging.png")),
                    paste0(cand_prefix, " — MDS Layer B (sample-belonging)"),
                    "Belonging\ncluster")
    }
    if (!is.null(layer_B2_minor$mds)) {
      make_mds_plot(layer_B2_minor$mds, hc_h05,
                    file.path(plot_dir, paste0(cand_prefix, "_mds_B2_minor.png")),
                    paste0(cand_prefix, " — MDS Layer B2 (raw-vector minor)"),
                    "Belonging\ncluster")
    }
    if (!is.null(layer_B2_major$mds)) {
      make_mds_plot(layer_B2_major$mds, hc_h05,
                    file.path(plot_dir, paste0(cand_prefix, "_mds_B2_major.png")),
                    paste0(cand_prefix, " — MDS Layer B2 (raw-vector major)"),
                    "Belonging\ncluster")
    }
    if (!is.null(layer_B2_012$mds)) {
      make_mds_plot(layer_B2_012$mds, hc_h05,
                    file.path(plot_dir, paste0(cand_prefix, "_mds_B2_012.png")),
                    paste0(cand_prefix, " — MDS Layer B2 (raw-vector 012)"),
                    "Belonging\ncluster")
    }
    if (!is.null(layer_C_minor$mds)) {
      make_mds_plot(layer_C_minor$mds, hc_h05,
                    file.path(plot_dir, paste0(cand_prefix, "_mds_C_minor.png")),
                    paste0(cand_prefix, " — MDS Layer C (scalar minor)"),
                    "Belonging\ncluster")
    }
    if (!is.null(layer_C_major$mds)) {
      make_mds_plot(layer_C_major$mds, hc_h05,
                    file.path(plot_dir, paste0(cand_prefix, "_mds_C_major.png")),
                    paste0(cand_prefix, " — MDS Layer C (scalar major)"),
                    "Belonging\ncluster")
    }

    # Dendrograms
    make_dendrogram_png(dist_A, file.path(plot_dir, paste0(cand_prefix, "_dendro_A_snp.png")),
                        paste0(cand_prefix, " — Dendrogram A (SNP)"), labels)
    make_dendrogram_png(dist_B_hamming, file.path(plot_dir, paste0(cand_prefix, "_dendro_B_belonging.png")),
                        paste0(cand_prefix, " — Dendrogram B (belonging)"), labels)
    make_dendrogram_png(layer_B2_minor$dist, file.path(plot_dir, paste0(cand_prefix, "_dendro_B2_minor.png")),
                        paste0(cand_prefix, " — Dendrogram B2 (raw-vector minor)"), labels)
    make_dendrogram_png(layer_B2_012$dist, file.path(plot_dir, paste0(cand_prefix, "_dendro_B2_012.png")),
                        paste0(cand_prefix, " — Dendrogram B2 (raw-vector 012)"), labels)
    make_dendrogram_png(layer_C_minor$dist, file.path(plot_dir, paste0(cand_prefix, "_dendro_C_minor.png")),
                        paste0(cand_prefix, " — Dendrogram C (scalar minor)"), labels)
  }

  message("[INFO] ", cand_prefix, " — window sample-belonging complete")
}

message("\n[DONE] STEP10c window sample-belonging analysis complete")