#!/usr/bin/env Rscript

# =============================================================================
# fuzzy_merge_score.R
#
# Fuzzy-composition merge scoring for Snake 1 core stitching.
#
# Instead of hard binary gates (pass/fail), each merge decision is scored
# on [0,1] using fuzzy composition of two parallel relations:
#
#   Relation M (membership): How well do sample band assignments transfer
#                             from core A to core B?
#   Relation G (geometric):  How continuous is the PCA geometry between
#                             the cores (sim_mat / MDS distances)?
#
# These combine via a fuzzy composition T = M ∘ G, producing a single
# merge compatibility score per core pair.
#
# Three merge tiers are defined by the ACCEPTANCE THRESHOLD on T:
#   merge_A: T >= 0.30 (generous)
#   merge_B: T >= 0.50 (moderate)
#   merge_C: T >= 0.70 (strict)
#
# Based on: Fuzzy max-min composition (Zadeh, 1965)
#   μ_{R∘S}(x,z) = max_{y∈Y} { min( μ_R(x,y), μ_S(y,z) ) }
#
# In our case:
#   x = bands in core A (band1, band2, band3)
#   y = samples (the 226 fish)
#   z = bands in core B (band1, band2, band3)
#
#   μ_R(band_a, sample) = soft membership of sample in band_a of core A
#   μ_S(sample, band_b) = soft membership of sample in band_b of core B
#   T(band_a, band_b)   = max_sample { min(μ_R, μ_S) }
#
# If the same samples define the same bands in both cores, the diagonal
# of T will be high (~1.0) and off-diagonal low (~0). Different inversions
# will produce a diffuse T matrix.
# =============================================================================

# ── SOFT BAND MEMBERSHIP ─────────────────────────────────────────────
# Instead of hard k=3 assignment, compute fuzzy membership based on
# distance to each cluster center. Uses inverse-distance weighting.
#
# Returns: N_samples × 3 matrix of soft memberships in [0,1],
#          rows sum to 1.

soft_band_membership <- function(dt, idxs, sample_names) {
  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  n_samples <- length(available)
  if (n_samples < 10) return(NULL)

  mat <- as.matrix(dt[idxs, ..available])
  if (nrow(mat) > 1) {
    avg_loadings <- colMeans(mat, na.rm = TRUE)
  } else {
    avg_loadings <- mat[1, ]
  }

  valid <- is.finite(avg_loadings)
  if (sum(valid) < 10) return(NULL)

  vals <- avg_loadings[valid]
  snames <- sub("^PC_1_", "", names(vals))

  # k=3 clustering
  km <- tryCatch(
    kmeans(vals, centers = 3, nstart = 5),
    error = function(e) NULL
  )
  if (is.null(km)) return(NULL)

  # Order centers: low → mid → high
  center_order <- order(km$centers[, 1])
  centers <- km$centers[center_order, 1]

  # Compute fuzzy membership: inverse distance to each center
  # Add a small epsilon to avoid division by zero
  n_valid <- sum(valid)
  membership <- matrix(0, nrow = n_valid, ncol = 3)

  for (bi in 1:3) {
    d <- abs(vals - centers[bi])
    membership[, bi] <- 1 / (d + 0.01)
  }

  # Normalize rows to sum to 1
  row_sums <- rowSums(membership)
  membership <- membership / row_sums

  # Name columns and rows
  colnames(membership) <- paste0("band", 1:3)
  rownames(membership) <- snames

  membership
}


# ── GEOMETRIC CONTINUITY SCORE ───────────────────────────────────────
# Compute how geometrically continuous the transition from core A to
# core B is, using the sim_mat (precomputed similarity matrix).
#
# Uses: tail of core A → head of core B, plus gap windows if present.
# Returns: scalar in [0, 1]

geometric_continuity <- function(core_a_idxs, core_b_idxs, sim_mat, gap_idxs = integer(0)) {
  # Tail of A → first of gap/B
  tail_a <- tail(core_a_idxs, 3L)
  head_b <- head(core_b_idxs, 3L)

  scores <- numeric(0)

  # A-tail to gap
  if (length(gap_idxs) > 0) {
    for (gi in head(gap_idxs, 3L)) {
      for (ai in tail_a) {
        s <- sim_mat[ai, gi]
        if (is.finite(s)) scores <- c(scores, s)
      }
    }
    # gap to B-head
    for (gi in tail(gap_idxs, 3L)) {
      for (bi in head_b) {
        s <- sim_mat[gi, bi]
        if (is.finite(s)) scores <- c(scores, s)
      }
    }
  }

  # A-tail to B-head (direct)
  for (ai in tail_a) {
    for (bi in head_b) {
      s <- sim_mat[ai, bi]
      if (is.finite(s)) scores <- c(scores, s)
    }
  }

  if (length(scores) == 0) return(0)
  mean(scores, na.rm = TRUE)
}


# ── FUZZY MAX-MIN COMPOSITION ────────────────────────────────────────
# Given two fuzzy relations:
#   R: membership_A (N_samples × 3)
#   S: membership_B (N_samples × 3)
#
# Compute T = R^T ∘ S using max-min:
#   T[band_a, band_b] = max_sample { min(R[sample, band_a], S[sample, band_b]) }
#
# Returns: 3×3 matrix T
#
# Interpretation:
#   - Diagonal high, off-diagonal low → same inversion system
#   - Diffuse T → different or no inversion
#   - Permuted diagonal → same inversion but PC1 flipped (ok, just relabel)

fuzzy_compose <- function(membership_a, membership_b) {
  # Match samples by name
  shared <- intersect(rownames(membership_a), rownames(membership_b))
  if (length(shared) < 10) return(NULL)

  R <- membership_a[shared, , drop = FALSE]  # N_shared × 3
  S <- membership_b[shared, , drop = FALSE]  # N_shared × 3

  # T[i,j] = max_k { min(R[k,i], S[k,j]) }
  T_mat <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:3) {
    for (j in 1:3) {
      mins <- pmin(R[, i], S[, j])
      T_mat[i, j] <- max(mins, na.rm = TRUE)
    }
  }

  colnames(T_mat) <- rownames(T_mat) <- paste0("band", 1:3)
  T_mat
}


# ── COMPOSITION SCORE ────────────────────────────────────────────────
# From the 3×3 T matrix, extract a single compatibility score.
#
# Best case (same inversion): T is a permutation matrix (one high value
# per row/column). Score = mean of the row maxima.
#
# Worst case (unrelated): T is uniform (~0.33 everywhere). Score ≈ 0.33.
#
# The score measures how "crisp" the band mapping is: do bands in A
# map cleanly to bands in B, or is everything blurred?

composition_score <- function(T_mat) {
  if (is.null(T_mat)) return(0)

  # For each row (band in A), find the best match in B
  row_max <- apply(T_mat, 1, max)

  # Check that the best matches are DIFFERENT bands (permutation-like)
  best_cols <- apply(T_mat, 1, which.max)
  n_unique <- length(unique(best_cols))

  # Permutation bonus: if each A-band maps to a different B-band, add bonus
  perm_bonus <- if (n_unique == 3) 0.15 else 0

  # Base score: mean of row maxima (how well each A-band finds a match)
  base <- mean(row_max)

  # Contrast: how much better is the best match vs second-best?
  contrast <- mean(apply(T_mat, 1, function(r) {
    sorted <- sort(r, decreasing = TRUE)
    if (length(sorted) >= 2 && sorted[1] > 0) {
      (sorted[1] - sorted[2]) / sorted[1]
    } else 0
  }))

  # Combined: base match quality + contrast + permutation bonus
  score <- 0.5 * base + 0.3 * contrast + 0.2 * min(1, perm_bonus / 0.15)
  max(0, min(1, score))
}


# ── COMBINED FUZZY MERGE SCORE ───────────────────────────────────────
# Combines the two relations:
#   1. Membership composition (fuzzy max-min) → how well do bands map?
#   2. Geometric continuity (sim_mat) → how smooth is the transition?
#
# The combination itself uses a fuzzy operator:
#   - For generous merge: use MAX (either evidence is enough)
#   - For moderate merge: use weighted mean
#   - For strict merge: use MIN (both must be good)
#
# This mirrors the fuzzy philosophy: generous is optimistic (any signal
# suffices), strict is pessimistic (all signals must agree).

compute_merge_score <- function(dt, core_a, core_b, sample_names, sim_mat,
                                 combination = c("mean", "max", "min"),
                                 w_membership = 0.6, w_geometric = 0.4) {
  combination <- match.arg(combination)

  core_a_idxs <- core_a$idxs
  core_b_idxs <- core_b$idxs

  # Gap windows between cores
  gap_start <- max(core_a_idxs) + 1L
  gap_end   <- min(core_b_idxs) - 1L
  gap_idxs <- if (gap_end >= gap_start) seq(gap_start, gap_end) else integer(0)
  gap_idxs <- gap_idxs[gap_idxs >= 1 & gap_idxs <= nrow(dt)]

  # ── Relation M: membership composition ────────────────────────────
  memb_a <- soft_band_membership(dt, core_a_idxs, sample_names)
  memb_b <- soft_band_membership(dt, core_b_idxs, sample_names)

  if (!is.null(memb_a) && !is.null(memb_b)) {
    T_mat <- fuzzy_compose(memb_a, memb_b)
    m_score <- composition_score(T_mat)
  } else {
    m_score <- NA_real_
  }

  # ── Relation G: geometric continuity ──────────────────────────────
  g_score <- geometric_continuity(core_a_idxs, core_b_idxs, sim_mat, gap_idxs)

  # ── Combine ───────────────────────────────────────────────────────
  if (is.na(m_score)) {
    # No membership data — fall back to geometric only
    final <- g_score
  } else if (g_score == 0) {
    # No geometric data — fall back to membership only
    final <- m_score
  } else {
    final <- switch(combination,
      "max"  = max(m_score, g_score),
      "min"  = min(m_score, g_score),
      "mean" = w_membership * m_score + w_geometric * g_score
    )
  }

  list(
    score = max(0, min(1, final)),
    membership_score = m_score,
    geometric_score = g_score,
    T_matrix = if (exists("T_mat")) T_mat else NULL,
    gap_size = length(gap_idxs)
  )
}


# ── FUZZY MERGE ENGINE ───────────────────────────────────────────────
# Replaces the hard-gate merge with fuzzy-scored merge.
# Each core pair gets a fuzzy merge score in [0,1].
# The tier threshold determines the cutoff.

run_merge_fuzzy <- function(all_cores, dt, sample_names, sim_mat, chr, params) {
  if (length(all_cores) == 0) return(list())

  core_starts <- vapply(all_cores, function(c) min(c$idxs), integer(1))
  all_cores <- all_cores[order(core_starts)]

  # Tier-dependent combination operator
  combination <- if (params$accept_threshold < 0.40) "max"      # generous: any signal
                 else if (params$accept_threshold < 0.65) "mean" # moderate: weighted
                 else "min"                                      # strict: both must agree

  merged <- list()
  merge_scores_log <- list()  # for diagnostics
  i <- 1L

  while (i <= length(all_cores)) {
    current <- all_cores[[i]]
    cur_idxs <- current$idxs
    cur_stat <- current$statuses
    cur_scor <- current$scores

    while (i < length(all_cores)) {
      nxt <- all_cores[[i + 1]]
      nxt_start <- min(nxt$idxs)
      cur_end   <- max(cur_idxs)
      gap_size  <- nxt_start - cur_end - 1L

      # Hard proximity ceiling (still needed — can't fuzzy-merge across 10 Mb)
      if (gap_size > params$max_gap_windows || gap_size < 0) break

      # ── Compute fuzzy merge score ──────────────────────────────────
      # Build a temporary "current core" object using the tail of the
      # growing merged region (so the reference adapts as we merge)
      cur_tail_obj <- list(idxs = tail(cur_idxs, 30L))
      ms <- compute_merge_score(dt, cur_tail_obj, nxt, sample_names, sim_mat,
                                 combination = combination)

      # Log for diagnostics
      merge_scores_log[[length(merge_scores_log) + 1]] <- data.table(
        chrom = chr, merge_family = params$name,
        core_end_bp = dt$end_bp[cur_end],
        next_start_bp = dt$start_bp[nxt_start],
        gap_windows = gap_size,
        fuzzy_score = round(ms$score, 4),
        membership_score = round(ms$membership_score %||% NA_real_, 4),
        geometric_score = round(ms$geometric_score, 4),
        decision = if (ms$score >= params$accept_threshold) "merge" else "reject"
      )

      if (ms$score < params$accept_threshold) break

      # ── Merge: concatenate with gap ────────────────────────────────
      bridge_range <- if (gap_size > 0) seq(cur_end + 1L, nxt_start - 1L) else integer(0)
      bridge_range <- bridge_range[bridge_range >= 1 & bridge_range <= nrow(dt)]

      cur_idxs <- c(cur_idxs, bridge_range, nxt$idxs)
      cur_stat <- c(cur_stat,
                    rep("bridge_fuzzy", length(bridge_range)),
                    nxt$statuses)
      cur_scor <- c(cur_scor,
                    rep(NA_real_, length(bridge_range)),
                    nxt$scores)
      i <- i + 1L
    }

    if (length(cur_idxs) >= params$min_windows) {
      merged[[length(merged) + 1]] <- list(
        idxs = cur_idxs,
        statuses = cur_stat,
        scores = cur_scor,
        merge_family = params$name,
        coherence = NULL
      )
    }

    i <- i + 1L
  }

  list(regions = merged, score_log = merge_scores_log)
}
