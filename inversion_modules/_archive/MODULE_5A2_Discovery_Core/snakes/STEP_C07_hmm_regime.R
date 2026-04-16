#!/usr/bin/env Rscript

# =============================================================================
# STEP10j_regime_hmm.R  (v1.0 — OPTIONAL LATER-STAGE MODULE)
#
# Sticky/hierarchical regime HMM built on Snake 1, Snake 2, and Snake 3
# evidence summaries. Snake 1 acts as the structural gate.
#
# STATUS: Design complete. Run ONLY after all 3 snake tracks + consensus exist.
#
# NOT the old naive HMM (raw tracks stacked into inversion-vs-not).
# This is a regime detector on snake-derived evidence summaries.
#
# HIDDEN STATES:
#   0: CORE_CONCORDANT       — all 3 snakes strong
#   1: PARTIAL_GROUP_SHIFT   — PCA strong, grouping weakening
#   2: GHSL_SUPPORTED_PARTIAL — PCA weak, GHSL still supports
#   3: TRANSITION            — between strong regime and background
#   4: WEAK_UNCERTAIN        — low support from all snakes
#   5: BACKGROUND_COMPLEX    — no snake supports inversion
#   6: QC_POOR               — data quality insufficient
#
# SNAKE 1 GATING (the key hierarchical rule):
#   If local-PCA regime collapses → strong core states MUST collapse.
#   A(STRONG) = all states
#   A(WEAK)   = all except CORE_CONCORDANT
#   A(NONE)   = only TRANSITION, WEAK_UNCERTAIN, BACKGROUND_COMPLEX, QC_POOR
#
#   Implemented as: emission probability = 0 for forbidden states.
#
# INPUTS:
#   consensus_windows.tsv.gz from STEP10i
#   snake1 outputs from STEP10e
#   snake2_track.tsv.gz from STEP10g
#   snake3_track.tsv.gz from STEP10h
#
# OUTPUTS:
#   hmm_input_features.tsv.gz     — per-window feature table
#   hmm_state_assignments.tsv.gz  — Viterbi decoded states
#   hmm_segments.tsv.gz           — contiguous segments
#   hmm_summary.tsv               — transition summary
#
# Usage:
#   Rscript STEP10j_regime_hmm.R \
#     <step10_outprefix> <snake1_dir> <snake2_dir> <snake3_dir> \
#     <consensus_dir> <outdir> [n_states=7] [sticky_alpha=0.85]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript STEP10j_regime_hmm.R ",
       "<step10_outprefix> <snake1_dir> <snake2_dir> <snake3_dir> ",
       "<consensus_dir> <outdir> [n_states=7] [sticky_alpha=0.85]")
}

step10_prefix <- args[1]
snake1_dir    <- args[2]
snake2_dir    <- args[3]
snake3_dir    <- args[4]
consensus_dir <- args[5]
outdir        <- args[6]
N_STATES      <- if (length(args) >= 7) as.integer(args[7]) else 7L
STICKY_ALPHA  <- if (length(args) >= 8) as.numeric(args[8]) else 0.85
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# STATE DEFINITIONS
# =============================================================================

STATE_LABELS <- c(
  "CORE_CONCORDANT",        # 0
  "PARTIAL_GROUP_SHIFT",    # 1
  "GHSL_SUPPORTED_PARTIAL", # 2
  "TRANSITION",             # 3
  "WEAK_UNCERTAIN",         # 4
  "BACKGROUND_COMPLEX",     # 5
  "QC_POOR"                 # 6
)

# =============================================================================
# SNAKE 1 GATING — allowed state sets
# =============================================================================

# Indices are 1-based (R convention)
ALLOWED_STRONG <- seq_len(N_STATES)       # all states
ALLOWED_WEAK   <- c(2, 3, 4, 5, 6, 7)    # no CORE_CONCORDANT
ALLOWED_NONE   <- c(4, 5, 6, 7)           # only TRANSITION + weaker

get_snake1_level <- function(pca_status, pca_score) {
  if (is.na(pca_status)) return("NONE")
  if (pca_status == "PASS" && !is.na(pca_score) && pca_score >= 0.55) return("STRONG")
  if (pca_status %in% c("PASS", "WEAK") && !is.na(pca_score) && pca_score >= 0.30) return("WEAK")
  return("NONE")
}

get_allowed_states <- function(level) {
  switch(level,
    "STRONG" = ALLOWED_STRONG,
    "WEAK"   = ALLOWED_WEAK,
    "NONE"   = ALLOWED_NONE,
    ALLOWED_NONE
  )
}

# =============================================================================
# BUILD FEATURE TABLE
# =============================================================================

message("[STEP10j] Building HMM feature table...")

# Load consensus
cons_file <- file.path(consensus_dir, "consensus_windows.tsv.gz")
if (!file.exists(cons_file)) stop("Missing consensus: ", cons_file)
cons <- fread(cons_file)

# Load snake2 track for NN preservation and matrix similarity
s2_file <- file.path(snake2_dir, "snake2_track.tsv.gz")
s2 <- if (file.exists(s2_file)) fread(s2_file) else data.table()

# Load snake3 track for GHSL metrics
s3_file <- file.path(snake3_dir, "snake3_track.tsv.gz")
s3 <- if (file.exists(s3_file)) fread(s3_file) else data.table()

# Load snake1 regions for per-window scores
s1_file <- file.path(snake1_dir, "snake_windows.tsv.gz")
s1 <- if (file.exists(s1_file)) fread(s1_file) else data.table()

# Build feature table on consensus grid
feat <- copy(cons[, .(chrom, global_window_id, start_bp, end_bp,
                       pca_status, group_status, ghsl_status,
                       n_snakes_supporting, support_signature)])

# Snake 1 features
if (nrow(s1) > 0) {
  s1_agg <- s1[, .(pca_score = mean(continuity_score, na.rm = TRUE),
                    pca_phase = snake_phase[1]),
               by = global_window_id]
  feat <- merge(feat, s1_agg, by = "global_window_id", all.x = TRUE)
} else {
  feat[, pca_score := NA_real_]
  feat[, pca_phase := "none"]
}

# Snake 2 features
if (nrow(s2) > 0) {
  s2_feat <- s2[, .(global_window_id,
                     group_nn = mean_nn_preservation,
                     group_matrix_sim = mean_matrix_sim,
                     group_nn_contrast = nn_contrast_vs_flank)]
  feat <- merge(feat, s2_feat, by = "global_window_id", all.x = TRUE)
} else {
  feat[, group_nn := NA_real_]
  feat[, group_matrix_sim := NA_real_]
  feat[, group_nn_contrast := NA_real_]
}

# Snake 3 features
if (nrow(s3) > 0) {
  s3_feat <- s3[, .(global_window_id,
                     ghsl_mean = mean_ghsl,
                     ghsl_spread = spread_ghsl,
                     ghsl_pass_frac = pass_fraction)]
  feat <- merge(feat, s3_feat, by = "global_window_id", all.x = TRUE)
} else {
  feat[, ghsl_mean := NA_real_]
  feat[, ghsl_spread := NA_real_]
  feat[, ghsl_pass_frac := NA_real_]
}

# Fill NAs with 0
num_cols <- c("pca_score", "group_nn", "group_matrix_sim", "group_nn_contrast",
              "ghsl_mean", "ghsl_spread", "ghsl_pass_frac", "n_snakes_supporting")
for (col in num_cols) {
  if (col %in% names(feat)) feat[is.na(get(col)), (col) := 0]
}

# Rolling features (±2 windows, per chromosome)
add_rolling <- function(dt, col, half_w = 2L) {
  roll_mean_col <- paste0(col, "_roll_mean")
  roll_var_col  <- paste0(col, "_roll_var")
  vals <- dt[[col]]
  n <- length(vals)
  rm <- rv <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - half_w)
    hi <- min(n, i + half_w)
    w <- vals[lo:hi]
    w <- w[is.finite(w)]
    rm[i] <- if (length(w) > 0) mean(w) else 0
    rv[i] <- if (length(w) > 1) var(w)  else 0
  }
  dt[[roll_mean_col]] <- round(rm, 4)
  dt[[roll_var_col]]  <- round(rv, 4)
  dt
}

setorder(feat, chrom, start_bp)
for (chr in unique(feat$chrom)) {
  idx <- feat$chrom == chr
  sub <- feat[idx]
  for (col in c("pca_score", "group_nn", "ghsl_mean")) {
    if (col %in% names(sub)) sub <- add_rolling(sub, col)
  }
  feat[idx] <- sub
}

# Snake 1 gating level
feat[, snake1_level := mapply(get_snake1_level, pca_status, pca_score)]

message("[STEP10j] Feature table: ", nrow(feat), " windows × ", ncol(feat), " features")

# =============================================================================
# STICKY TRANSITION MATRIX
# =============================================================================

build_transition_matrix <- function(K = N_STATES, alpha = STICKY_ALPHA) {
  TM <- matrix((1 - alpha) / (K - 1), K, K)
  diag(TM) <- alpha

  # Biologically motivated adjustments
  # CORE → PARTIAL and CORE → TRANSITION slightly favored
  TM[1, 2] <- TM[1, 2] * 1.5  # CORE → PARTIAL
  TM[1, 4] <- TM[1, 4] * 1.3  # CORE → TRANSITION
  TM[2, 4] <- TM[2, 4] * 1.5  # PARTIAL → TRANSITION
  TM[4, 6] <- TM[4, 6] * 1.3  # TRANSITION → BACKGROUND

  # Implausible jumps discouraged
  TM[1, 7] <- TM[1, 7] * 0.3  # CORE → QC_POOR
  TM[6, 1] <- TM[6, 1] * 0.3  # BACKGROUND → CORE

  # Renormalize rows
  TM <- TM / rowSums(TM)
  TM
}

# =============================================================================
# GAUSSIAN EMISSION (diagonal covariance)
# =============================================================================

# Feature columns for HMM
HMM_FEATURES <- c("pca_score", "group_nn", "group_matrix_sim",
                   "group_nn_contrast", "ghsl_mean", "ghsl_spread",
                   "ghsl_pass_frac", "n_snakes_supporting")

# Keep only features that exist
HMM_FEATURES <- intersect(HMM_FEATURES, names(feat))

log_gaussian_emission <- function(x, mu, sigma2) {
  # x, mu, sigma2: vectors of length D
  D <- length(x)
  s2 <- pmax(sigma2, 1e-6)
  -0.5 * sum((x - mu)^2 / s2 + log(2 * pi * s2))
}

# =============================================================================
# VITERBI WITH SNAKE 1 GATING
# =============================================================================

viterbi_gated <- function(obs_mat, snake1_levels, TM, mu_mat, sigma2_mat, pi0) {
  # obs_mat: T × D
  # snake1_levels: character vector length T
  # TM: K × K transition
  # mu_mat: K × D emission means
  # sigma2_mat: K × D emission variances
  # pi0: K initial probabilities

  TT <- nrow(obs_mat)
  K  <- nrow(TM)
  log_TM <- log(TM + 1e-300)

  log_delta <- matrix(-Inf, TT, K)
  psi       <- matrix(0L, TT, K)

  for (tt in seq_len(TT)) {
    allowed <- get_allowed_states(snake1_levels[tt])

    # Emission log-probabilities
    log_emit <- numeric(K)
    for (k in seq_len(K)) {
      if (k %in% allowed) {
        log_emit[k] <- log_gaussian_emission(obs_mat[tt, ], mu_mat[k, ], sigma2_mat[k, ])
      } else {
        log_emit[k] <- -Inf  # GATING: forbidden state
      }
    }

    if (tt == 1L) {
      log_delta[1, ] <- log(pi0 + 1e-300) + log_emit
    } else {
      for (k in seq_len(K)) {
        vals <- log_delta[tt - 1, ] + log_TM[, k]
        best_j <- which.max(vals)
        log_delta[tt, k] <- vals[best_j] + log_emit[k]
        psi[tt, k] <- best_j
      }
    }
  }

  # Backtrack
  states <- integer(TT)
  states[TT] <- which.max(log_delta[TT, ])
  for (tt in (TT - 1):1) {
    states[tt] <- psi[tt + 1, states[tt + 1]]
  }

  log_probs <- vapply(seq_len(TT), function(t) log_delta[t, states[t]], numeric(1))

  list(states = states, log_probs = log_probs)
}

# =============================================================================
# INITIALIZE EMISSION PARAMETERS (heuristic from data)
# =============================================================================

init_emission_params <- function(obs_mat, K = N_STATES) {
  D <- ncol(obs_mat)
  mu <- matrix(0, K, D)
  s2 <- matrix(1, K, D)

  # Global stats
  global_mean <- colMeans(obs_mat, na.rm = TRUE)
  global_sd   <- apply(obs_mat, 2, sd, na.rm = TRUE)
  global_sd[global_sd < 1e-6] <- 1e-6

  # Heuristic initialization: spread means across the feature space
  # CORE = high scores, BACKGROUND = low scores
  for (k in seq_len(K)) {
    frac <- (k - 1) / (K - 1)  # 0 for CORE, 1 for QC_POOR
    mu[k, ] <- global_mean + (0.5 - frac) * global_sd
    s2[k, ] <- global_sd^2
  }

  list(mu = mu, sigma2 = s2)
}

# =============================================================================
# RUN PER CHROMOSOME
# =============================================================================

message("[STEP10j] Running Viterbi HMM per chromosome...")

TM <- build_transition_matrix()
pi0 <- rep(1 / N_STATES, N_STATES)

all_assign <- list()
all_segs   <- list()
all_em_diag <- list()   # NEW: EM convergence diagnostics

for (chr in unique(feat$chrom)) {
  sub <- feat[chrom == chr]
  n <- nrow(sub)
  if (n < 5) {
    message("[SKIP] ", chr, ": only ", n, " windows")
    next
  }

  # Build observation matrix
  obs_mat <- as.matrix(sub[, ..HMM_FEATURES])
  obs_mat[!is.finite(obs_mat)] <- 0

  # Initialize emissions from this chromosome's data
  params <- init_emission_params(obs_mat, K = N_STATES)

  # ── EM TRAINING LOOP (v8 addition) ──────────────────────────────────
  # Hard-EM: iterate Viterbi → re-estimate emissions → repeat
  # This replaces the old heuristic-only initialization
  EM_MAX_ITER <- 15L
  EM_CONV_THRESH <- 1e-3
  prev_loglik <- -Inf
  em_trace <- list()  # track convergence per iteration

  for (em_iter in seq_len(EM_MAX_ITER)) {
    # E-step: Viterbi decoding with current parameters
    em_result <- viterbi_gated(
      obs_mat, sub$snake1_level,
      TM, params$mu, params$sigma2, pi0
    )

    # Check convergence
    current_loglik <- sum(em_result$log_probs[is.finite(em_result$log_probs)])
    improvement <- abs(current_loglik - prev_loglik) / max(1, abs(prev_loglik))

    # Track state distribution per iteration
    state_counts <- tabulate(em_result$states, nbins = N_STATES)
    em_trace[[em_iter]] <- data.table(
      chrom = chr, iteration = em_iter,
      loglik = round(current_loglik, 4),
      improvement = round(improvement, 6),
      n_states_used = sum(state_counts > 0),
      state_counts = paste(state_counts, collapse = ","),
      largest_state = STATE_LABELS[which.max(state_counts)],
      smallest_nonzero = if (any(state_counts > 0))
        STATE_LABELS[which.min(state_counts[state_counts > 0])] else NA_character_
    )

    if (em_iter > 1 && improvement < EM_CONV_THRESH) {
      message("[STEP10j] ", chr, ": EM converged at iteration ", em_iter,
              " (improvement=", round(improvement, 6), ")")
      break
    }
    prev_loglik <- current_loglik

    # M-step: re-estimate emission parameters from hard assignments
    for (k in seq_len(N_STATES)) {
      assigned <- which(em_result$states == k)
      if (length(assigned) >= 3) {
        params$mu[k, ] <- colMeans(obs_mat[assigned, , drop = FALSE], na.rm = TRUE)
        params$sigma2[k, ] <- pmax(
          apply(obs_mat[assigned, , drop = FALSE], 2, var, na.rm = TRUE),
          1e-4
        )
      }
      # If fewer than 3 assigned, keep previous parameters (regularization)
    }
  }

  if (em_iter == EM_MAX_ITER) {
    message("[STEP10j] ", chr, ": EM reached max iterations (", EM_MAX_ITER, ")")
  }

  # Save EM convergence trace
  if (length(em_trace) > 0) {
    all_em_diag[[chr]] <- rbindlist(em_trace)
  }

  # Final Viterbi with trained parameters
  result <- viterbi_gated(
    obs_mat, sub$snake1_level,
    TM, params$mu, params$sigma2, pi0
  )

  # State assignments
  sub[, hidden_state := result$states]
  sub[, state_label := STATE_LABELS[hidden_state]]
  sub[, log_prob := result$log_probs]
  sub[, posterior_confidence := exp(result$log_probs - max(result$log_probs))]

  all_assign[[chr]] <- sub

  # Build segments
  segs <- list()
  seg_start <- 1L
  for (i in 2:n) {
    if (result$states[i] != result$states[i - 1] || i == n) {
      end_i <- if (i == n && result$states[i] == result$states[i - 1]) i else i - 1
      if (i == n && result$states[i] != result$states[i - 1]) end_i <- i - 1

      segs[[length(segs) + 1]] <- data.table(
        chrom = chr,
        seg_start_bp = sub$start_bp[seg_start],
        seg_end_bp = sub$end_bp[end_i],
        n_windows = end_i - seg_start + 1L,
        dominant_state = result$states[seg_start],
        state_label = STATE_LABELS[result$states[seg_start]],
        mean_confidence = round(mean(exp(result$log_probs[seg_start:end_i] -
                                          max(result$log_probs[seg_start:end_i]))), 4)
      )

      if (i == n && result$states[i] != result$states[i - 1]) {
        segs[[length(segs) + 1]] <- data.table(
          chrom = chr, seg_start_bp = sub$start_bp[i],
          seg_end_bp = sub$end_bp[i], n_windows = 1L,
          dominant_state = result$states[i],
          state_label = STATE_LABELS[result$states[i]],
          mean_confidence = 1.0
        )
      }

      seg_start <- i
    }
  }

  all_segs[[chr]] <- rbindlist(segs)

  n_per_state <- table(factor(result$states, levels = seq_len(N_STATES)))
  msg_parts <- paste0(STATE_LABELS, "=", n_per_state, collapse = " ")
  message("[STEP10j] ", chr, ": ", msg_parts)
}

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

assign_dt <- if (length(all_assign) > 0) rbindlist(all_assign, fill = TRUE) else {
  data.table(chrom = character(), global_window_id = integer(),
             hidden_state = integer(), state_label = character())
}

seg_dt <- if (length(all_segs) > 0) rbindlist(all_segs, fill = TRUE) else {
  data.table(chrom = character(), state_label = character())
}

f1 <- file.path(outdir, "hmm_input_features.tsv.gz")
f2 <- file.path(outdir, "hmm_state_assignments.tsv.gz")
f3 <- file.path(outdir, "hmm_segments.tsv.gz")
f4 <- file.path(outdir, "hmm_summary.tsv")
f5 <- file.path(outdir, "hmm_em_convergence.tsv.gz")

fwrite(feat, f1, sep = "\t")
fwrite(assign_dt[, .(chrom, global_window_id, start_bp, end_bp,
                      hidden_state, state_label, posterior_confidence,
                      snake1_level, pca_status, group_status, ghsl_status)],
       f2, sep = "\t")
fwrite(seg_dt, f3, sep = "\t")

# Summary
summary_dt <- assign_dt[, .(n = .N), by = .(chrom, state_label)]
fwrite(summary_dt, f4, sep = "\t")

# EM convergence diagnostics
em_diag_dt <- if (length(all_em_diag) > 0) rbindlist(all_em_diag, fill = TRUE) else data.table()
fwrite(em_diag_dt, f5, sep = "\t")

message("\n[DONE] STEP10j regime HMM complete")
message("  ", f1)
message("  ", f2)
message("  ", f3)
message("  ", f4)

if (nrow(assign_dt) > 0) {
  message("\n[STEP10j] Global state distribution:")
  global_tab <- assign_dt[, .N, by = state_label][order(-N)]
  for (i in seq_len(nrow(global_tab))) {
    message("  ", global_tab$state_label[i], ": ", global_tab$N[i])
  }
}
