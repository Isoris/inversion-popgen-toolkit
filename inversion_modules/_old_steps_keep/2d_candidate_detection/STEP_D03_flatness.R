# ============================================================================
# STEP_D03_flatness.R — Similarity decay profile per candidate (v2)
# ============================================================================
# For each candidate, characterize the full shape of the similarity-vs-
# distance curve within the block.
#
# Reports:
#   - Similarity at 5 quantile distances (shape descriptor)
#   - Piecewise slopes in 3 thirds (near/mid/far)
#   - Slope sign flips (artifact detector)
#   - Artifact positions (where the dip-then-recovery happens)
#   - Far/near ratio (squareness)
#   - Monotonicity score (1 = perfectly non-increasing, 0 = chaotic)
#
# Uses MEDIAN for the decay profile (robust to single-window artifacts)
# Uses Theil-Sen slopes (robust to outliers)
# ============================================================================

compute_flatness <- function(candidates, smat) {
  results <- list()

  for (i in seq_len(nrow(candidates))) {
    si <- candidates$start[i]
    ei <- candidates$end[i]
    width <- ei - si + 1
    cid <- candidates$candidate_id[i]

    if (width < 8) {
      results[[i]] <- make_empty_flatness_row(cid)
      next
    }

    # ---- Build distance profile ----
    max_dist <- width - 1
    dist_medians <- numeric(max_dist)
    dist_counts  <- integer(max_dist)

    for (d in seq_len(max_dist)) {
      vals <- numeric(0)
      for (row_i in si:(ei - d)) {
        col_j <- row_i + d
        if (col_j <= ei) {
          v <- smat[row_i, col_j]
          if (is.finite(v)) vals <- c(vals, v)
        }
      }
      dist_medians[d] <- if (length(vals) > 0) median(vals) else NA
      dist_counts[d]  <- length(vals)
    }

    valid <- !is.na(dist_medians) & dist_counts >= 2
    if (sum(valid) < 6) {
      results[[i]] <- make_empty_flatness_row(cid)
      next
    }

    dists   <- which(valid)
    profile <- dist_medians[valid]
    n_valid <- length(profile)

    # ---- Quantile decay values ----
    q_idx <- unique(pmax(1, pmin(n_valid,
      round(quantile(seq_len(n_valid), probs=c(0.1, 0.25, 0.5, 0.75, 0.9))))))

    sim_q10 <- profile[q_idx[1]]
    sim_q25 <- profile[q_idx[2]]
    sim_q50 <- profile[q_idx[3]]
    sim_q75 <- profile[q_idx[4]]
    sim_q90 <- profile[q_idx[5]]

    # ---- Far/near ratio ----
    near_n <- max(1, round(n_valid * 0.2))
    far_start <- max(1, n_valid - near_n + 1)
    near_sim <- median(profile[1:near_n])
    far_sim  <- median(profile[far_start:n_valid])
    far_near_ratio <- if (near_sim > 0) far_sim / near_sim else 0

    # ---- Piecewise slopes (thirds) ----
    third <- max(2, n_valid %/% 3)
    slope_near <- theil_sen_slope(profile[1:third])
    slope_mid  <- theil_sen_slope(profile[(third+1):min(2*third, n_valid)])
    slope_far  <- theil_sen_slope(profile[min(2*third+1, n_valid):n_valid])

    # ---- Slope sign flips (artifact detection) ----
    smoothed <- running_median(profile, span = max(3, n_valid %/% 10))
    local_slopes <- diff(smoothed)
    sign_info    <- detect_sign_flips(local_slopes)

    artifact_bins <- if (sign_info$n_neg_to_pos > 0 && length(sign_info$positions) > 0) {
      paste(dists[pmin(sign_info$positions, n_valid)], collapse=",")
    } else ""

    # ---- Monotonicity ----
    n_pairs <- length(smoothed) - 1
    monotonicity <- if (n_pairs > 0) {
      sum(diff(smoothed) <= 0.001, na.rm=TRUE) / n_pairs
    } else NA

    # ---- Flatness ----
    profile_range <- max(profile, na.rm=TRUE) - min(profile, na.rm=TRUE)
    profile_mean  <- mean(profile, na.rm=TRUE)
    flatness <- if (profile_mean > 0) max(0, 1 - profile_range / profile_mean) else NA

    results[[i]] <- data.frame(
      candidate_id    = cid,
      flatness        = round(flatness, 4),
      far_near_ratio  = round(far_near_ratio, 4),
      monotonicity    = round(monotonicity, 4),
      sim_q10 = round(sim_q10, 4),
      sim_q25 = round(sim_q25, 4),
      sim_q50 = round(sim_q50, 4),
      sim_q75 = round(sim_q75, 4),
      sim_q90 = round(sim_q90, 4),
      slope_near      = round(slope_near, 6),
      slope_mid       = round(slope_mid, 6),
      slope_far       = round(slope_far, 6),
      n_sign_flips    = sign_info$n_flips,
      n_artifact_dips = sign_info$n_neg_to_pos,
      artifact_bins   = artifact_bins,
      has_artifact    = sign_info$n_neg_to_pos > 0,
      stringsAsFactors = FALSE
    )
  }

  return(do.call(rbind, results))
}

# ---- Empty row ----
make_empty_flatness_row <- function(cid) {
  data.frame(
    candidate_id = cid, flatness = NA_real_, far_near_ratio = NA_real_,
    monotonicity = NA_real_,
    sim_q10 = NA_real_, sim_q25 = NA_real_, sim_q50 = NA_real_,
    sim_q75 = NA_real_, sim_q90 = NA_real_,
    slope_near = NA_real_, slope_mid = NA_real_, slope_far = NA_real_,
    n_sign_flips = NA_integer_, n_artifact_dips = NA_integer_,
    artifact_bins = "", has_artifact = FALSE, stringsAsFactors = FALSE
  )
}

# ---- Theil-Sen slope (median of pairwise slopes, robust) ----
theil_sen_slope <- function(seg) {
  n <- length(seg)
  if (n < 2) return(0)
  # For speed, subsample if segment is large
  if (n > 50) {
    idx <- sort(sample(n, 50))
    seg <- seg[idx]
    n <- 50
  }
  slopes <- c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      slopes <- c(slopes, (seg[j] - seg[i]) / (j - i))
    }
  }
  return(median(slopes, na.rm=TRUE))
}

# ---- Running median smoother ----
running_median <- function(x, span=5) {
  n <- length(x)
  if (n <= span) return(x)
  half <- span %/% 2
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    out[i] <- median(x[lo:hi], na.rm=TRUE)
  }
  return(out)
}

# ---- Detect slope sign flips ----
detect_sign_flips <- function(local_slopes) {
  n <- length(local_slopes)
  if (n < 2) return(list(n_flips=0L, positions=integer(0), n_neg_to_pos=0L))

  signs <- sign(local_slopes)
  # Fill zeros with previous sign
  for (i in 2:n) { if (signs[i] == 0) signs[i] <- signs[i-1] }

  if (all(signs == 0 | is.na(signs))) {
    return(list(n_flips=0L, positions=integer(0), n_neg_to_pos=0L))
  }

  flips <- 0L; neg_to_pos <- 0L; flip_pos <- c()
  for (i in 2:n) {
    if (!is.na(signs[i]) && !is.na(signs[i-1]) && signs[i] != signs[i-1]) {
      flips <- flips + 1L
      flip_pos <- c(flip_pos, i)
      if (signs[i-1] < 0 && signs[i] > 0) neg_to_pos <- neg_to_pos + 1L
    }
  }
  return(list(n_flips=flips, positions=as.integer(flip_pos), n_neg_to_pos=neg_to_pos))
}
