#!/usr/bin/env Rscript

# =============================================================================
# local_q_projector.R
#
# WINDOWED LOCAL ADMIXTURE ESTIMATION
#
# Clean-room implementation of the admixture model E-step for per-window
# local ancestry probability estimation. Uses FIXED ancestral allele
# frequencies (F matrix) from a genome-wide NGSadmix run and computes
# per-sample, per-window Q vectors (local ancestry proportions).
#
# This is NOT a copy of NGSadmix. It implements the same published EM
# algorithm (Pritchard et al. 2000; Alexander et al. 2009; Skotte et al.
# 2013) but with F fixed, reducing each window to a single E-step +
# optional Q-only EM refinement.
#
# The mathematical model (Skotte et al. 2013, Genetics 195:693-702):
#   For site j, individual i, genotype g ∈ {0,1,2}:
#     P(g_ij | Q_i, F_j) = Σ_k Q_ik × Binom(g | 2, f_jk)
#   where f_jk = allele frequency in population k at site j.
#
#   With genotype likelihoods (GL):
#     P(GL_ij | Q_i, F_j) = Σ_g GL_ijg × P(g | Q_i, F_j)
#
#   E-step: compute expected genotype given current Q,F:
#     E[g_ij] = (p1 + 2*p2) / (p0 + p1 + p2)
#     where pg = P(g) × GL_ijg
#
#   Q update (M-step for Q with F fixed):
#     Q_ik_new ∝ Q_ik × Σ_j [ E[g_ij]×f_jk/(Σ_k' Q_ik'×f_jk')
#                             + (2-E[g_ij])×(1-f_jk)/(Σ_k' Q_ik'×(1-f_jk')) ]
#
# Reference:
#   Skotte L, Korneliussen TS, Albrechtsen A (2013) Estimating individual
#   admixture proportions from next generation sequencing data. Genetics
#   195:693-702. https://doi.org/10.1534/genetics.113.154138
#
# Inputs:
#   --beagle_dir <dir>      Per-chr BEAGLE GL files (3 GL columns per sample)
#   --fmat <file>           F matrix from NGSadmix (sites × K frequencies)
#   --qinit <file>          Initial Q from genome-wide NGSadmix (samples × K)
#   --window_registry <file> Window registry with chr, start, end, window_id
#   --outdir <dir>
#   [--k <int>]             K populations (auto-detected from qinit)
#   [--em_iter <int>]       Q-only EM iterations per window (default: 5)
#   [--chrom <chr>]         Process single chromosome
#   [--ncores <int>]        Parallel cores (default: 1)
#
# Outputs:
#   local_Q_per_window.tsv.gz   — window_id × sample × K Q proportions
#   local_Q_metrics.tsv.gz      — window_id × delta12, H, ENA per sample
#   local_Q_summary.tsv.gz      — window_id × aggregated metrics
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ─── Parse arguments ──────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
beagle_dir <- NULL; fmat_file <- NULL; qinit_file <- NULL
window_file <- NULL; outdir <- "local_Q"; chrom_filter <- NULL
EM_ITER <- 5L; NCORES <- 1L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--beagle_dir" && i < length(args))      { beagle_dir <- args[i+1]; i <- i+2 }
  else if (a == "--fmat" && i < length(args))        { fmat_file <- args[i+1]; i <- i+2 }
  else if (a == "--qinit" && i < length(args))       { qinit_file <- args[i+1]; i <- i+2 }
  else if (a == "--window_registry" && i < length(args)) { window_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))      { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--chrom" && i < length(args))       { chrom_filter <- args[i+1]; i <- i+2 }
  else if (a == "--em_iter" && i < length(args))     { EM_ITER <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--ncores" && i < length(args))      { NCORES <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(beagle_dir) || is.null(qinit_file))
  stop("Required: --beagle_dir, --qinit")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
EPS <- 1e-10

message("[localQ] Local Q Projector")
message("[localQ] BEAGLE: ", beagle_dir)
message("[localQ] Q init: ", qinit_file)
message("[localQ] F mat: ", fmat_file %||% "will estimate from Q + BEAGLE")
message("[localQ] EM iter: ", EM_ITER)

# ─── Load genome-wide Q (samples × K) ────────────────────────────────
q_init <- as.matrix(fread(qinit_file, header = FALSE))
K <- ncol(q_init)
n_samples <- nrow(q_init)
message("[localQ] Q init: ", n_samples, " samples × K=", K)

# ─── Load window registry ────────────────────────────────────────────
win_dt <- NULL
if (!is.null(window_file) && file.exists(window_file)) {
  win_dt <- fread(window_file)
  if (!is.null(chrom_filter)) win_dt <- win_dt[chrom == chrom_filter]
  message("[localQ] Windows: ", nrow(win_dt))
}

# ─── Load F matrix if provided ────────────────────────────────────────
# F matrix: n_sites × K (ancestral allele frequency per population)
# If not provided, we estimate it from the BEAGLE + Q (one M-step)
f_mat <- NULL
if (!is.null(fmat_file) && file.exists(fmat_file)) {
  f_mat <- as.matrix(fread(fmat_file, header = FALSE))
  message("[localQ] F matrix: ", nrow(f_mat), " sites × K=", ncol(f_mat))
}

# =============================================================================
# CORE ALGORITHM: per-window Q estimation with fixed F
# =============================================================================

# Compute local Q for one window given:
#   gl:     matrix n_sites × (3 × n_ind), BEAGLE format (GL per genotype)
#   f_win:  matrix n_sites × K, allele frequencies per population
#   q_prior: matrix n_ind × K, prior Q (from genome-wide or previous window)
#   n_iter: number of Q-only EM iterations
#
# Returns: matrix n_ind × K of local Q proportions

compute_local_Q <- function(gl, f_win, q_prior, n_iter = EM_ITER) {
  n_sites <- nrow(gl)
  n_ind <- nrow(q_prior)
  K_loc <- ncol(q_prior)

  if (n_sites < 5 || n_ind < 10) return(q_prior)

  # Initialize Q from prior
  Q <- q_prior

  for (iter in seq_len(n_iter)) {
    # Accumulators for Q update (prodA, prodB in NGSadmix notation)
    prodA <- matrix(0, n_ind, K_loc)
    prodB <- matrix(0, n_ind, K_loc)

    for (j in seq_len(n_sites)) {
      f_j <- f_win[j, ]  # K-length vector of allele freqs
      f_j <- pmin(pmax(f_j, EPS), 1 - EPS)

      for (i in seq_len(n_ind)) {
        # Compute P(allele) for this individual at this site
        # p = Σ_k Q_ik × f_jk
        p <- sum(Q[i, ] * f_j)
        p <- pmin(pmax(p, EPS), 1 - EPS)

        # Genotype probabilities under HWE: P(g|p) × GL
        pp0 <- (1 - p)^2 * gl[j, 3*(i-1) + 1]  # P(g=0) × GL(g=0)
        pp1 <- 2 * p * (1 - p) * gl[j, 3*(i-1) + 2]  # P(g=1) × GL(g=1)
        pp2 <- p^2 * gl[j, 3*(i-1) + 3]  # P(g=2) × GL(g=2)
        denom <- pp0 + pp1 + pp2
        if (denom < EPS) next

        # Expected genotype (0-2 scale)
        exp_g <- (pp1 + 2 * pp2) / denom

        # Accumulate sufficient statistics for Q update
        for (k in seq_len(K_loc)) {
          prodA[i, k] <- prodA[i, k] + exp_g / (1 - p) * f_j[k]
          prodB[i, k] <- prodB[i, k] + (2 - exp_g) / p * (1 - f_j[k])
        }
      }
    }

    # Q update: Q_new_ik ∝ Q_ik × (prodA_ik + prodB_ik)
    for (i in seq_len(n_ind)) {
      for (k in seq_len(K_loc)) {
        Q[i, k] <- Q[i, k] * (prodA[i, k] + prodB[i, k])
      }
      # Normalize to sum to 1
      row_sum <- sum(Q[i, ])
      if (row_sum > EPS) Q[i, ] <- Q[i, ] / row_sum
      # Clamp
      Q[i, ] <- pmin(pmax(Q[i, ], EPS), 1 - EPS)
      Q[i, ] <- Q[i, ] / sum(Q[i, ])
    }
  }

  Q
}

# Estimate F from current Q and BEAGLE GL (one M-step)
# Used when F matrix is not provided
estimate_F_from_Q <- function(gl, Q, positions) {
  n_sites <- nrow(gl)
  n_ind <- nrow(Q)
  K_loc <- ncol(Q)

  F_est <- matrix(0.5, n_sites, K_loc)

  for (j in seq_len(n_sites)) {
    for (i in seq_len(n_ind)) {
      p <- sum(Q[i, ] * F_est[j, ])
      p <- pmin(pmax(p, EPS), 1 - EPS)

      pp0 <- (1 - p)^2 * gl[j, 3*(i-1) + 1]
      pp1 <- 2 * p * (1 - p) * gl[j, 3*(i-1) + 2]
      pp2 <- p^2 * gl[j, 3*(i-1) + 3]
      denom <- pp0 + pp1 + pp2
      if (denom < EPS) next

      exp_g <- (pp1 + 2 * pp2) / denom

      for (k in seq_len(K_loc)) {
        F_est[j, k] <- F_est[j, k] + Q[i, k] * exp_g
      }
    }
    # Normalize
    for (k in seq_len(K_loc)) {
      F_est[j, k] <- F_est[j, k] / (2 * sum(Q[, k]))
      F_est[j, k] <- pmin(pmax(F_est[j, k], EPS), 1 - EPS)
    }
  }
  F_est
}

# Compute metrics from Q vector
q_metrics <- function(q_vec) {
  q_vec <- q_vec[is.finite(q_vec) & q_vec > 0]
  if (length(q_vec) < 2) return(list(delta12 = NA, entropy = NA, ena = NA))
  qs <- sort(q_vec, decreasing = TRUE)
  H <- -sum(q_vec * log(q_vec + 1e-15))
  list(delta12 = qs[1] - qs[2], entropy = H, ena = exp(H))
}

# =============================================================================
# MAIN LOOP: per chromosome
# =============================================================================

# Find BEAGLE files
beagle_files <- list.files(beagle_dir, pattern = "\\.beagle\\.gz$|\\.beagle$",
                            full.names = TRUE)
if (!is.null(chrom_filter)) {
  beagle_files <- beagle_files[grepl(chrom_filter, beagle_files)]
}
message("[localQ] BEAGLE files: ", length(beagle_files))

all_q_rows <- list()
all_metric_rows <- list()
all_summary_rows <- list()

for (bf in beagle_files) {
  message("\n[localQ] Loading ", basename(bf), " ...")
  bgl <- tryCatch(fread(bf, nThread = 4L), error = function(e) NULL)
  if (is.null(bgl)) { message("  [SKIP] cannot read"); next }

  # BEAGLE format: marker, allele1, allele2, then 3 GL columns per sample
  # Columns 4:(3+3*n_ind)
  marker_col <- names(bgl)[1]
  n_gl_cols <- ncol(bgl) - 3
  if (n_gl_cols < 9 || n_gl_cols %% 3 != 0) {
    message("  [SKIP] unexpected column count: ", ncol(bgl))
    next
  }
  n_ind_bgl <- n_gl_cols %/% 3

  # Extract positions from marker names (format: chr_pos)
  markers <- bgl[[1]]
  positions <- as.integer(sub(".*_", "", markers))

  # GL matrix: n_sites × (3 × n_ind)
  gl_mat <- as.matrix(bgl[, 4:ncol(bgl)])
  n_sites_total <- nrow(gl_mat)

  chr_name <- sub("_.*", "", markers[1])
  message("  ", chr_name, ": ", n_sites_total, " sites × ", n_ind_bgl, " individuals")

  if (n_ind_bgl != n_samples) {
    message("  [WARN] sample count mismatch: BEAGLE=", n_ind_bgl, " Q=", n_samples)
    n_use <- min(n_ind_bgl, n_samples)
  } else {
    n_use <- n_samples
  }

  # Get windows for this chromosome
  if (!is.null(win_dt)) {
    chr_wins <- win_dt[chrom == chr_name]
  } else {
    # Auto-generate windows: 100 SNPs, step 20
    win_size <- 100L; win_step <- 20L
    n_wins <- max(1, (n_sites_total - win_size) %/% win_step + 1)
    chr_wins <- data.table(
      global_window_id = paste0(chr_name, "_w", seq_len(n_wins)),
      chrom = chr_name,
      start_snp_idx = seq(1, by = win_step, length.out = n_wins),
      n_snps = win_size
    )
    chr_wins[, end_snp_idx := pmin(start_snp_idx + n_snps - 1L, n_sites_total)]
    chr_wins[, start_bp := positions[start_snp_idx]]
    chr_wins[, end_bp := positions[end_snp_idx]]
    message("  Auto-generated ", nrow(chr_wins), " windows (", win_size, " SNPs, step ", win_step, ")")
  }

  message("  Processing ", nrow(chr_wins), " windows (EM_ITER=", EM_ITER, ") ...")
  t_chr <- proc.time()

  for (wi in seq_len(nrow(chr_wins))) {
    w <- chr_wins[wi]
    wid <- w$global_window_id

    # Get sites in window
    if ("start_snp_idx" %in% names(w)) {
      site_idx <- w$start_snp_idx:w$end_snp_idx
    } else {
      site_idx <- which(positions >= w$start_bp & positions <= w$end_bp)
    }
    if (length(site_idx) < 10) next

    # Extract GL sub-matrix and Q prior
    gl_win <- gl_mat[site_idx, seq_len(3 * n_use), drop = FALSE]
    q_prior <- q_init[seq_len(n_use), , drop = FALSE]

    # Estimate F for this window if not provided globally
    if (!is.null(f_mat) && max(site_idx) <= nrow(f_mat)) {
      f_win <- f_mat[site_idx, , drop = FALSE]
    } else {
      f_win <- estimate_F_from_Q(gl_win, q_prior, positions[site_idx])
    }

    # Compute local Q
    q_local <- compute_local_Q(gl_win, f_win, q_prior, n_iter = EM_ITER)

    # Per-sample metrics
    for (si in seq_len(n_use)) {
      m <- q_metrics(q_local[si, ])
      all_metric_rows[[length(all_metric_rows) + 1L]] <- data.table(
        window_id = wid, sample_idx = si,
        delta12 = round(m$delta12, 4),
        entropy = round(m$entropy, 4),
        ena = round(m$ena, 4),
        max_q = round(max(q_local[si, ]), 4),
        assigned_pop = which.max(q_local[si, ])
      )
    }

    # Window summary
    all_metrics <- rbindlist(lapply(seq_len(n_use), function(si) {
      m <- q_metrics(q_local[si, ])
      data.table(d12 = m$delta12, H = m$entropy, ena = m$ena)
    }))
    all_summary_rows[[length(all_summary_rows) + 1L]] <- data.table(
      window_id = wid, chrom = chr_name,
      start_bp = positions[site_idx[1]],
      end_bp = positions[tail(site_idx, 1)],
      n_sites = length(site_idx),
      mean_delta12 = round(mean(all_metrics$d12, na.rm = TRUE), 4),
      mean_entropy = round(mean(all_metrics$H, na.rm = TRUE), 4),
      mean_ena = round(mean(all_metrics$ena, na.rm = TRUE), 4),
      sd_delta12 = round(sd(all_metrics$d12, na.rm = TRUE), 4)
    )

    if (wi %% 500 == 0) message("    window ", wi, "/", nrow(chr_wins))
  }

  elapsed <- round((proc.time() - t_chr)[3], 1)
  message("  ", chr_name, " done in ", elapsed, "s (",
          round(nrow(chr_wins) / max(elapsed, 0.1), 0), " windows/sec)")
}

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

metric_dt <- if (length(all_metric_rows) > 0) rbindlist(all_metric_rows) else data.table()
summary_dt <- if (length(all_summary_rows) > 0) rbindlist(all_summary_rows) else data.table()

fwrite(metric_dt, file.path(outdir, "local_Q_metrics.tsv.gz"), sep = "\t")
fwrite(summary_dt, file.path(outdir, "local_Q_summary.tsv.gz"), sep = "\t")

message("\n[localQ] === OUTPUTS ===")
message("  Metrics: ", file.path(outdir, "local_Q_metrics.tsv.gz"),
        " (", nrow(metric_dt), " rows)")
message("  Summary: ", file.path(outdir, "local_Q_summary.tsv.gz"),
        " (", nrow(summary_dt), " rows)")
message("[DONE]")
