#!/usr/bin/env Rscript

# =============================================================================
# STEP10g_snake2_community_continuity.R  (v2.0)
#
# SNAKE 2: Combined hard + fuzzy neighborhood preservation with
# middle-band stability tracking.
#
# WHAT THIS SOLVES:
#   Snake 1 detects "windows with 3-stripe-like PCA structure."
#   But two different inversion systems can both look 3-stripe.
#   Current Snake 2 (v1) only asks "do neighbors persist?" — it misses
#   the case where neighbors persist but the ROLES change (different
#   samples fill the middle band in adjacent windows).
#
#   v2 adds:
#   1. Middle-band stability: are the SAME samples in the middle?
#   2. Tail-matching with flip: handle PC1 polarity flips between windows
#   3. Fuzzy neighborhood preservation: weighted overlap for bridges/peripherals
#   4. Combined score: hard + fuzzy + middle + tails in one output
#
# DESIGN DECISIONS (critical thinking applied):
#   - NO per-window u,v rotation (that's a followup tool for pooled regions,
#     not viable per 100-SNP window with noisy local PCA)
#   - NO DBSCAN per window (too fragile at small marker counts)
#   - YES k-means k=3 on PC1 for middle/tail split (simple, robust, fast)
#   - YES fuzzy weights via inverse-distance (simple, interpretable)
#   - YES one algorithm, not Snake2a/Snake2b
#   - Quality tiers are SOFT WEIGHTS, not hard discard rules
#
# OUTPUTS:
#   snake2_track.tsv.gz   — per-window: hard NN, fuzzy NN, middle stability,
#                           tail match, combined score, QC status,
#                           profile concordance + per-stripe Hamming (diagnostic)
#   snake2_summary.tsv    — per-chromosome summary
#   snake2_band_assignments.tsv.gz — per-window band membership (for Snake 4)
#   snake2_profile_layer.tsv.gz    — per-window per-stripe profile diversity (diagnostic)
#
# Usage:
#   Rscript STEP10g_snake2_community_continuity.R \
#     <step10_outprefix> <dosage_dir> <samples_ind> <outdir>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript STEP10g_snake2_community_continuity.R ",
       "<step10_outprefix> <dosage_dir> <samples_ind> <outdir>")
}

step10_prefix <- args[1]
dosage_dir    <- args[2]
samples_ind   <- args[3]
outdir        <- args[4]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PARAMETERS
# =============================================================================

NN_K                    <- 10L     # k for k-NN
NN_PASS_THRESH          <- 0.55    # combined score PASS
NN_WEAK_THRESH          <- 0.30    # combined score WEAK (below = FAIL)
FLANK_WINDOW_COUNT      <- 10L
BACKGROUND_CONTRAST_MIN <- 0.08
MIN_MARKERS_PER_WINDOW  <- 15L
MIN_SAMPLES             <- 20L

# Score weights for combined score
W_HARD_NN     <- 0.30   # hard kNN preservation
W_FUZZY_NN    <- 0.20   # fuzzy kNN preservation
W_MIDDLE      <- 0.30   # middle-band stability (most diagnostic)
W_TAIL        <- 0.20   # tail-matching with flip

# Profile concordance layer (diagnostic, does not enter combined score)
PROFILE_MIN_STRIPE_N  <- 5L     # min samples in stripe for profile analysis
PROFILE_DISCRETIZE    <- TRUE   # round(dosage) → 0/1/2

# =============================================================================
# LOAD
# =============================================================================

message("[S2v2] Loading data...")
mds_obj <- readRDS(paste0(step10_prefix, ".mds.rds"))
per_chr <- mds_obj$per_chr

sample_names <- tryCatch({
  s <- fread(samples_ind, header = FALSE)[[1]]; s[nchar(s) > 0]
}, error = function(e) NULL)
if (is.null(sample_names)) stop("Cannot read samples.ind")
n_samples <- length(sample_names)
chroms <- names(per_chr)
message("[S2v2] Chromosomes: ", length(chroms), ", Samples: ", n_samples)

# =============================================================================
# DISTANCE MATRIX (same as v1: centered correlation)
# =============================================================================

build_centered_corr_dist <- function(Xw) {
  Xc <- Xw - rowMeans(Xw, na.rm = TRUE)
  sim <- cor(Xc, use = "pairwise.complete.obs")
  sim[is.na(sim)] <- 0
  dist_mat <- 1 - abs(sim)
  diag(dist_mat) <- 0
  list(dist = dist_mat, sim = sim)
}

# =============================================================================
# HARD kNN
# =============================================================================

compute_knn <- function(dist_mat, k) {
  n <- nrow(dist_mat); k <- min(k, n - 1L)
  knn <- matrix(0L, n, k)
  for (i in seq_len(n)) {
    d <- dist_mat[i, ]; d[i] <- Inf; d[!is.finite(d)] <- Inf
    knn[i, ] <- order(d)[seq_len(k)]
  }
  knn
}

hard_nn_preservation <- function(knn_a, knn_b) {
  n <- nrow(knn_a); k <- ncol(knn_a)
  ps <- numeric(n)
  for (i in seq_len(n)) ps[i] <- length(intersect(knn_a[i,], knn_b[i,])) / k
  mean(ps)
}

# =============================================================================
# FUZZY kNN (inverse-distance weights)
# =============================================================================

compute_fuzzy_knn <- function(dist_mat, k) {
  # Returns: n × k matrix of neighbor indices + n × k matrix of weights
  n <- nrow(dist_mat); k <- min(k, n - 1L)
  idx <- matrix(0L, n, k)
  wts <- matrix(0, n, k)

  for (i in seq_len(n)) {
    d <- dist_mat[i, ]; d[i] <- Inf; d[!is.finite(d)] <- Inf
    ord <- order(d)[seq_len(k)]
    idx[i, ] <- ord
    dists <- d[ord]
    # Inverse-distance weights, add small epsilon to avoid Inf
    inv_d <- 1 / (dists + 1e-6)
    wts[i, ] <- inv_d / sum(inv_d)
  }
  list(idx = idx, wts = wts)
}

fuzzy_nn_preservation <- function(fknn_a, fknn_b) {
  # Weighted Jaccard: for each sample, compute overlap of fuzzy neighborhoods
  n <- nrow(fknn_a$idx)
  scores <- numeric(n)

  for (i in seq_len(n)) {
    # Build full weight vectors (length = n_samples)
    wa <- rep(0, n); wb <- rep(0, n)
    wa[fknn_a$idx[i, ]] <- fknn_a$wts[i, ]
    wb[fknn_b$idx[i, ]] <- fknn_b$wts[i, ]
    # Weighted Jaccard: sum(min) / sum(max)
    s_min <- sum(pmin(wa, wb))
    s_max <- sum(pmax(wa, wb))
    scores[i] <- if (s_max > 0) s_min / s_max else 0
  }
  mean(scores)
}

# =============================================================================
# MIDDLE-BAND STABILITY (the key new metric)
# =============================================================================

identify_3band <- function(dist_mat) {
  # Quick 3-way split using the first PC (cmdscale dim 1)
  # This is much simpler and more robust than DBSCAN per window
  n <- nrow(dist_mat)
  if (n < 10) return(NULL)

  # MDS dim 1
  mds <- tryCatch(cmdscale(as.dist(dist_mat), k = 1), error = function(e) NULL)
  if (is.null(mds) || length(mds) != n) return(NULL)

  pc1 <- as.numeric(mds)

  # k-means k=3 on PC1
  km <- tryCatch(kmeans(pc1, centers = 3, nstart = 5), error = function(e) NULL)
  if (is.null(km)) return(NULL)

  # Order clusters by center: leftmost, middle, rightmost
  center_order <- order(km$centers[, 1])
  role_map <- setNames(c("tailA", "middle", "tailB"), as.character(center_order))
  roles <- role_map[as.character(km$cluster)]

  # Soft middle confidence: distance to middle center / spread
  mid_center <- km$centers[center_order[2], 1]
  mid_spread <- sd(pc1[roles == "middle"])
  if (!is.finite(mid_spread) || mid_spread < 1e-8) mid_spread <- sd(pc1)
  if (!is.finite(mid_spread) || mid_spread < 1e-8) mid_spread <- 1

  # Soft membership: how middle-like is each sample (Gaussian-ish weight)
  middle_weight <- exp(-0.5 * ((pc1 - mid_center) / mid_spread)^2)
  middle_weight <- middle_weight / max(middle_weight)

  list(
    roles = roles,
    pc1 = pc1,
    middle_idx = which(roles == "middle"),
    tailA_idx = which(roles == "tailA"),
    tailB_idx = which(roles == "tailB"),
    middle_weight = middle_weight  # soft membership [0,1]
  )
}

middle_band_preservation <- function(band_a, band_b) {
  # Hard: Jaccard of middle-membership sets
  if (is.null(band_a) || is.null(band_b)) return(list(hard = NA_real_, soft = NA_real_))

  mid_a <- band_a$middle_idx
  mid_b <- band_b$middle_idx
  hard_jacc <- length(intersect(mid_a, mid_b)) / max(1, length(union(mid_a, mid_b)))

  # Soft: correlation of middle-weight vectors
  wa <- band_a$middle_weight
  wb <- band_b$middle_weight
  if (length(wa) == length(wb)) {
    soft <- cor(wa, wb, use = "complete.obs")
    if (!is.finite(soft)) soft <- 0
    soft <- max(0, soft)
  } else {
    soft <- 0
  }

  list(hard = hard_jacc, soft = soft)
}

# =============================================================================
# TAIL-MATCHING WITH FLIP
# =============================================================================

tail_match_with_flip <- function(band_a, band_b) {
  # Compare tail assignments under both orientations, keep better match
  if (is.null(band_a) || is.null(band_b)) return(NA_real_)

  tA_a <- band_a$tailA_idx; tB_a <- band_a$tailB_idx
  tA_b <- band_b$tailA_idx; tB_b <- band_b$tailB_idx

  # Direct: A→A, B→B
  jAA <- length(intersect(tA_a, tA_b)) / max(1, length(union(tA_a, tA_b)))
  jBB <- length(intersect(tB_a, tB_b)) / max(1, length(union(tB_a, tB_b)))
  direct <- (jAA + jBB) / 2

  # Swapped: A→B, B→A
  jAB <- length(intersect(tA_a, tB_b)) / max(1, length(union(tA_a, tB_b)))
  jBA <- length(intersect(tB_a, tA_b)) / max(1, length(union(tB_a, tA_b)))
  swapped <- (jAB + jBA) / 2

  max(direct, swapped)
}

# =============================================================================
# WITHIN-STRIPE PROFILE CONCORDANCE (diagnostic layer — parallel to kNN)
#
# At 9x coverage with hard filters, 100-SNP windows have clean enough dosage
# for discretized 0/1/2 Hamming distances within each stripe. This detects
# sub-haplotype structure (BB_1 vs BB_2) that kNN preservation cannot see,
# because kNN only asks "are neighbors preserved?" not "do same-stripe
# samples carry the same marker pattern?"
#
# Does NOT enter the combined score. Does NOT gate PASS/WEAK/FAIL.
# Produces diagnostic columns + a separate output file.
# =============================================================================

compute_within_stripe_profile <- function(Xw, band, min_n = PROFILE_MIN_STRIPE_N) {
  # Xw: markers × samples dosage matrix
  # band: output of identify_3band (with roles, middle_idx, tailA_idx, tailB_idx)
  # Returns: per-stripe diversity + between-window concordance metrics

  if (is.null(band)) return(NULL)
  n_mark <- nrow(Xw)

  # Discretize: round(dosage) → 0/1/2
  Gw <- round(Xw)
  Gw[is.na(Gw)] <- -1L  # explicit missing
  storage.mode(Gw) <- "integer"

  stripe_stats <- function(idx, label) {
    if (length(idx) < min_n) return(list(label = label, n = length(idx),
      mean_hamming = NA_real_, hamming_sd = NA_real_, n_informative = NA_real_))
    Gs <- Gw[, idx, drop = FALSE]  # markers × stripe_samples
    ns <- ncol(Gs)
    # Pairwise Hamming (fast: matrix ops)
    hamming_vals <- numeric(ns * (ns - 1L) / 2L)
    k <- 0L
    for (i in seq_len(ns - 1L)) {
      for (j in (i + 1L):ns) {
        valid <- Gs[, i] >= 0L & Gs[, j] >= 0L
        nv <- sum(valid)
        if (nv > 0) {
          k <- k + 1L
          hamming_vals[k] <- sum(Gs[valid, i] != Gs[valid, j]) / nv
        }
      }
    }
    hamming_vals <- hamming_vals[seq_len(k)]
    # Count markers polymorphic within this stripe
    n_poly <- sum(apply(Gs, 1, function(r) { r <- r[r >= 0L]; length(unique(r)) > 1L }))
    list(label = label, n = ns, mean_hamming = mean(hamming_vals),
         hamming_sd = sd(hamming_vals), n_informative = n_poly)
  }

  tA <- stripe_stats(band$tailA_idx, "tailA")
  tB <- stripe_stats(band$tailB_idx, "tailB")
  mid <- stripe_stats(band$middle_idx, "middle")

  # Return discretized matrix too (for between-window concordance)
  list(tailA = tA, tailB = tB, middle = mid, Gw = Gw)
}

profile_concordance_between_windows <- function(prof_a, prof_b, band_a, band_b,
                                                 min_n = PROFILE_MIN_STRIPE_N) {
  # Compare within-stripe Hamming structure between adjacent windows.
  # For each stripe: take the samples assigned to that stripe in BOTH windows,
  # compute their pairwise Hamming in window A and in window B, then correlate.
  # High correlation = same sub-structure preserved. Low = structure changed.

  if (is.null(prof_a) || is.null(prof_b) || is.null(band_a) || is.null(band_b))
    return(NA_real_)

  concordances <- numeric(0)
  for (stripe in c("tailA", "tailB", "middle")) {
    idx_a <- switch(stripe, tailA = band_a$tailA_idx, tailB = band_a$tailB_idx,
                    middle = band_a$middle_idx)
    idx_b <- switch(stripe, tailA = band_b$tailA_idx, tailB = band_b$tailB_idx,
                    middle = band_b$middle_idx)
    shared <- intersect(idx_a, idx_b)
    if (length(shared) < min_n) next

    # Hamming distances within this shared set, in window A and window B
    Ga <- prof_a$Gw[, shared, drop = FALSE]
    Gb <- prof_b$Gw[, shared, drop = FALSE]
    ns <- length(shared)
    ha <- hb <- numeric(ns * (ns - 1L) / 2L)
    k <- 0L
    for (i in seq_len(ns - 1L)) {
      for (j in (i + 1L):ns) {
        va <- Ga[, i] >= 0L & Ga[, j] >= 0L
        vb <- Gb[, i] >= 0L & Gb[, j] >= 0L
        nva <- sum(va); nvb <- sum(vb)
        k <- k + 1L
        ha[k] <- if (nva > 0) sum(Ga[va, i] != Ga[va, j]) / nva else NA_real_
        hb[k] <- if (nvb > 0) sum(Gb[vb, i] != Gb[vb, j]) / nvb else NA_real_
      }
    }
    ha <- ha[seq_len(k)]; hb <- hb[seq_len(k)]
    valid <- is.finite(ha) & is.finite(hb)
    if (sum(valid) >= 3) {
      r <- cor(ha[valid], hb[valid], method = "spearman")
      if (is.finite(r)) concordances <- c(concordances, r)
    }
  }
  if (length(concordances) > 0) mean(concordances) else NA_real_
}

# =============================================================================
# COMBINED SCORE + QC
# =============================================================================

compute_combined_score <- function(hard_nn, fuzzy_nn, middle_hard, tail_match) {
  # Replace NAs with 0 for combination
  h <- if (is.finite(hard_nn)) hard_nn else 0
  f <- if (is.finite(fuzzy_nn)) fuzzy_nn else 0
  m <- if (is.finite(middle_hard)) middle_hard else 0
  t <- if (is.finite(tail_match)) tail_match else 0
  W_HARD_NN * h + W_FUZZY_NN * f + W_MIDDLE * m + W_TAIL * t
}

classify_combined <- function(combined, contrast) {
  above_bg <- !is.na(contrast) && contrast >= BACKGROUND_CONTRAST_MIN
  if (!is.finite(combined)) return("FAIL")
  if (combined >= NN_PASS_THRESH && above_bg) return("PASS")
  if (combined >= NN_WEAK_THRESH) return("WEAK")
  "FAIL"
}

# =============================================================================
# MAIN: PER-CHROMOSOME
# =============================================================================

all_track <- list()
all_summary <- list()
all_band_rows <- list()  # NEW: per-window band assignments for Snake 4

for (chr in chroms) {
  chr_obj <- per_chr[[chr]]
  if (is.null(chr_obj)) next

  dt <- as.data.table(chr_obj$out_dt)[order(start_bp)]
  n_win <- nrow(dt)
  message("\n[S2v2] ═══════ ", chr, " (", n_win, " windows) ═══════")
  if (n_win < 3) { message("[SKIP]"); next }

  dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(dosage_dir, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) { message("[SKIP] no dosage"); next }

  dos <- fread(dos_file); sites <- fread(sites_file)
  dos_cols <- setdiff(names(dos), "marker")
  if (all(grepl("^Ind", dos_cols)) && length(sample_names) == length(dos_cols)) {
    setnames(dos, old = dos_cols, new = sample_names); dos_cols <- sample_names
  }
  if (length(dos_cols) < MIN_SAMPLES) { message("[SKIP] too few samples"); next }

  # ── Pre-compute per-window ──────────────────────────────────────────
  message("[S2v2] Pre-computing windows...")
  t0 <- proc.time()

  win_dists  <- vector("list", n_win)
  win_knns   <- vector("list", n_win)
  win_fknns  <- vector("list", n_win)
  win_bands  <- vector("list", n_win)
  win_profs  <- vector("list", n_win)  # profile layer
  win_nmark  <- integer(n_win)

  for (wi in seq_len(n_win)) {
    w_keep <- which(sites$pos >= dt$start_bp[wi] & sites$pos <= dt$end_bp[wi])
    if (length(w_keep) < MIN_MARKERS_PER_WINDOW) next

    Xw <- as.matrix(dos[w_keep, ..dos_cols]); storage.mode(Xw) <- "double"
    win_nmark[wi] <- nrow(Xw)

    res <- build_centered_corr_dist(Xw)
    win_dists[[wi]] <- res$dist
    win_knns[[wi]]  <- compute_knn(res$dist, NN_K)
    win_fknns[[wi]] <- compute_fuzzy_knn(res$dist, NN_K)
    win_bands[[wi]] <- identify_3band(res$dist)

    # Profile layer: discretize + within-stripe Hamming
    if (!is.null(win_bands[[wi]])) {
      win_profs[[wi]] <- compute_within_stripe_profile(Xw, win_bands[[wi]])
    }
  }

  n_valid <- sum(!vapply(win_dists, is.null, logical(1)))
  elapsed <- (proc.time() - t0)[3]
  message("[S2v2] ", n_valid, "/", n_win, " valid windows (", round(elapsed, 1), "s)")

  # ── Export per-window band assignments for Snake 4 ───────────────────
  for (wi in seq_len(n_win)) {
    band <- win_bands[[wi]]
    if (is.null(band)) next

    wid <- dt$global_window_id[wi]
    n_samp <- length(band$roles)

    # Store: which samples are in middle (hard) + soft weights
    all_band_rows[[length(all_band_rows) + 1]] <- data.table(
      chrom = chr,
      global_window_id = wid,
      start_bp = dt$start_bp[wi],
      end_bp = dt$end_bp[wi],
      n_middle = length(band$middle_idx),
      n_tailA = length(band$tailA_idx),
      n_tailB = length(band$tailB_idx),
      middle_idx_csv = paste(band$middle_idx, collapse = ","),
      tailA_idx_csv = paste(band$tailA_idx, collapse = ","),
      tailB_idx_csv = paste(band$tailB_idx, collapse = ","),
      middle_weight_csv = paste(round(band$middle_weight, 4), collapse = ","),
      middle_fraction = round(length(band$middle_idx) / n_samp, 4)
    )
  }

  # ── Per-window scoring ──────────────────────────────────────────────
  for (wi in seq_len(n_win)) {
    if (is.null(win_dists[[wi]])) {
      all_track[[length(all_track) + 1]] <- data.table(
        chrom = chr, global_window_id = dt$global_window_id[wi],
        start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
        n_markers = win_nmark[wi],
        hard_nn_prev = NA_real_, hard_nn_next = NA_real_, mean_hard_nn = NA_real_,
        fuzzy_nn_prev = NA_real_, fuzzy_nn_next = NA_real_, mean_fuzzy_nn = NA_real_,
        middle_stability_prev = NA_real_, middle_stability_next = NA_real_,
        mean_middle_stability = NA_real_,
        tail_match_prev = NA_real_, tail_match_next = NA_real_,
        mean_tail_match = NA_real_,
        combined_score = NA_real_,
        nn_contrast_vs_flank = NA_real_,
        profile_concordance = NA_real_,
        profile_tailA_hamming = NA_real_, profile_tailB_hamming = NA_real_,
        profile_middle_hamming = NA_real_,
        profile_tailA_poly = NA_integer_, profile_tailB_poly = NA_integer_,
        profile_middle_poly = NA_integer_,
        snake2_status = "FAIL"
      )
      next
    }

    # Hard NN
    h_prev <- if (wi > 1 && !is.null(win_knns[[wi-1]]))
      hard_nn_preservation(win_knns[[wi]], win_knns[[wi-1]]) else NA_real_
    h_next <- if (wi < n_win && !is.null(win_knns[[wi+1]]))
      hard_nn_preservation(win_knns[[wi]], win_knns[[wi+1]]) else NA_real_

    # Fuzzy NN
    f_prev <- if (wi > 1 && !is.null(win_fknns[[wi-1]]))
      fuzzy_nn_preservation(win_fknns[[wi]], win_fknns[[wi-1]]) else NA_real_
    f_next <- if (wi < n_win && !is.null(win_fknns[[wi+1]]))
      fuzzy_nn_preservation(win_fknns[[wi]], win_fknns[[wi+1]]) else NA_real_

    # Middle stability
    mb_prev <- if (wi > 1 && !is.null(win_bands[[wi-1]]) && !is.null(win_bands[[wi]]))
      middle_band_preservation(win_bands[[wi]], win_bands[[wi-1]]) else list(hard = NA_real_)
    mb_next <- if (wi < n_win && !is.null(win_bands[[wi+1]]) && !is.null(win_bands[[wi]]))
      middle_band_preservation(win_bands[[wi]], win_bands[[wi+1]]) else list(hard = NA_real_)

    # Tail matching
    tm_prev <- if (wi > 1 && !is.null(win_bands[[wi-1]]) && !is.null(win_bands[[wi]]))
      tail_match_with_flip(win_bands[[wi]], win_bands[[wi-1]]) else NA_real_
    tm_next <- if (wi < n_win && !is.null(win_bands[[wi+1]]) && !is.null(win_bands[[wi]]))
      tail_match_with_flip(win_bands[[wi]], win_bands[[wi+1]]) else NA_real_

    # Average (handling NAs)
    avg <- function(a, b) {
      vals <- c(a, b)[is.finite(c(a, b))]
      if (length(vals) > 0) mean(vals) else NA_real_
    }
    mean_h <- avg(h_prev, h_next)
    mean_f <- avg(f_prev, f_next)
    mean_m <- avg(mb_prev$hard, mb_next$hard)
    mean_t <- avg(tm_prev, tm_next)

    combined <- compute_combined_score(mean_h, mean_f, mean_m, mean_t)

    # ── Flanking contrast (on combined score, not just hard NN) ────────
    flank_scores <- c()
    left_lo <- max(2L, wi - FLANK_WINDOW_COUNT)
    left_hi <- max(1L, wi - 1L)
    if (left_lo <= left_hi) {
      for (li in left_lo:left_hi) {
        if (!is.null(win_knns[[li]]) && !is.null(win_knns[[li-1]])) {
          fh <- hard_nn_preservation(win_knns[[li]], win_knns[[li-1]])
          fm <- if (!is.null(win_bands[[li]]) && !is.null(win_bands[[li-1]]))
            middle_band_preservation(win_bands[[li]], win_bands[[li-1]])$hard else 0
          flank_scores <- c(flank_scores, 0.5 * fh + 0.5 * (fm %||% 0))
        }
      }
    }
    right_lo <- min(n_win, wi + 1L)
    right_hi <- min(n_win - 1L, wi + FLANK_WINDOW_COUNT)
    if (right_lo <= right_hi) {
      for (ri in right_lo:right_hi) {
        if (!is.null(win_knns[[ri]]) && !is.null(win_knns[[ri+1]])) {
          fh <- hard_nn_preservation(win_knns[[ri]], win_knns[[ri+1]])
          fm <- if (!is.null(win_bands[[ri]]) && !is.null(win_bands[[ri+1]]))
            middle_band_preservation(win_bands[[ri]], win_bands[[ri+1]])$hard else 0
          flank_scores <- c(flank_scores, 0.5 * fh + 0.5 * (fm %||% 0))
        }
      }
    }

    flank_base <- if (length(flank_scores) > 0) mean(flank_scores) else 0.5
    contrast <- if (is.finite(combined)) combined - flank_base else NA_real_
    status <- classify_combined(combined, contrast)

    # ── Profile concordance (between adjacent windows, within stripes) ──
    pc_prev <- if (wi > 1 && !is.null(win_profs[[wi]]) && !is.null(win_profs[[wi-1]]))
      profile_concordance_between_windows(win_profs[[wi]], win_profs[[wi-1]],
                                          win_bands[[wi]], win_bands[[wi-1]]) else NA_real_
    pc_next <- if (wi < n_win && !is.null(win_profs[[wi]]) && !is.null(win_profs[[wi+1]]))
      profile_concordance_between_windows(win_profs[[wi]], win_profs[[wi+1]],
                                          win_bands[[wi]], win_bands[[wi+1]]) else NA_real_
    prof_conc <- avg(pc_prev, pc_next)

    # Per-stripe profile stats for this window
    p <- win_profs[[wi]]
    tA_h <- if (!is.null(p) && is.finite(p$tailA$mean_hamming)) p$tailA$mean_hamming else NA_real_
    tB_h <- if (!is.null(p) && is.finite(p$tailB$mean_hamming)) p$tailB$mean_hamming else NA_real_
    mi_h <- if (!is.null(p) && is.finite(p$middle$mean_hamming)) p$middle$mean_hamming else NA_real_
    tA_p <- if (!is.null(p) && is.finite(p$tailA$n_informative)) as.integer(p$tailA$n_informative) else NA_integer_
    tB_p <- if (!is.null(p) && is.finite(p$tailB$n_informative)) as.integer(p$tailB$n_informative) else NA_integer_
    mi_p <- if (!is.null(p) && is.finite(p$middle$n_informative)) as.integer(p$middle$n_informative) else NA_integer_

    all_track[[length(all_track) + 1]] <- data.table(
      chrom = chr, global_window_id = dt$global_window_id[wi],
      start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
      n_markers = win_nmark[wi],
      hard_nn_prev = round(h_prev, 4), hard_nn_next = round(h_next, 4),
      mean_hard_nn = round(mean_h, 4),
      fuzzy_nn_prev = round(f_prev, 4), fuzzy_nn_next = round(f_next, 4),
      mean_fuzzy_nn = round(mean_f, 4),
      middle_stability_prev = round(mb_prev$hard, 4),
      middle_stability_next = round(mb_next$hard, 4),
      mean_middle_stability = round(mean_m, 4),
      tail_match_prev = round(tm_prev, 4), tail_match_next = round(tm_next, 4),
      mean_tail_match = round(mean_t, 4),
      combined_score = round(combined, 4),
      nn_contrast_vs_flank = round(contrast, 4),
      profile_concordance = round(prof_conc, 4),
      profile_tailA_hamming = round(tA_h, 4), profile_tailB_hamming = round(tB_h, 4),
      profile_middle_hamming = round(mi_h, 4),
      profile_tailA_poly = tA_p, profile_tailB_poly = tB_p,
      profile_middle_poly = mi_p,
      snake2_status = status
    )
  }

  # Summary
  chr_trk <- rbindlist(all_track[vapply(all_track, function(x) x$chrom[1] == chr, logical(1))])
  all_summary[[length(all_summary) + 1]] <- data.table(
    chrom = chr, n_windows = n_win, n_valid = n_valid,
    n_pass = sum(chr_trk$snake2_status == "PASS", na.rm = TRUE),
    n_weak = sum(chr_trk$snake2_status == "WEAK", na.rm = TRUE),
    n_fail = sum(chr_trk$snake2_status == "FAIL", na.rm = TRUE),
    mean_combined = round(mean(chr_trk$combined_score, na.rm = TRUE), 4),
    mean_hard_nn = round(mean(chr_trk$mean_hard_nn, na.rm = TRUE), 4),
    mean_middle_stability = round(mean(chr_trk$mean_middle_stability, na.rm = TRUE), 4)
  )

  message("[S2v2] ", chr, ": PASS=", sum(chr_trk$snake2_status == "PASS", na.rm = TRUE),
          " WEAK=", sum(chr_trk$snake2_status == "WEAK", na.rm = TRUE),
          " FAIL=", sum(chr_trk$snake2_status == "FAIL", na.rm = TRUE))
}

# =============================================================================
# WRITE
# =============================================================================

track_dt <- if (length(all_track) > 0) rbindlist(all_track) else {
  data.table(chrom = character(), snake2_status = character())
}
summary_dt <- if (length(all_summary) > 0) rbindlist(all_summary) else data.table()

fwrite(track_dt, file.path(outdir, "snake2_track.tsv.gz"), sep = "\t")
fwrite(summary_dt, file.path(outdir, "snake2_summary.tsv"), sep = "\t")

# Band assignments for Snake 4 consumption
band_dt <- if (length(all_band_rows) > 0) rbindlist(all_band_rows) else {
  data.table(chrom = character(), global_window_id = integer())
}
f_band <- file.path(outdir, "snake2_band_assignments.tsv.gz")
fwrite(band_dt, f_band, sep = "\t")
message("  ", f_band, " (", nrow(band_dt), " windows with band data)")

message("\n[DONE] Snake 2 v2 complete")
message("  Total: PASS=", sum(track_dt$snake2_status == "PASS", na.rm = TRUE),
        " WEAK=", sum(track_dt$snake2_status == "WEAK", na.rm = TRUE),
        " FAIL=", sum(track_dt$snake2_status == "FAIL", na.rm = TRUE))

# Profile layer: extract per-stripe detail from track for downstream
prof_cols <- c("chrom", "global_window_id", "start_bp", "end_bp", "n_markers",
               "profile_concordance",
               "profile_tailA_hamming", "profile_tailB_hamming", "profile_middle_hamming",
               "profile_tailA_poly", "profile_tailB_poly", "profile_middle_poly",
               "snake2_status")
prof_cols <- intersect(prof_cols, names(track_dt))
if (length(prof_cols) > 0 && nrow(track_dt) > 0) {
  prof_dt <- track_dt[, ..prof_cols]
  f_prof <- file.path(outdir, "snake2_profile_layer.tsv.gz")
  fwrite(prof_dt, f_prof, sep = "\t")
  message("  ", f_prof, " (", nrow(prof_dt), " windows)")
}
