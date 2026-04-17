#!/usr/bin/env Rscript

# =============================================================================
# ancestry_bridge.R — Bridge between Module 2B/2C and inversion pipeline
#
# Provides get_local_Q(): a single function that returns per-window local
# ancestry Q vectors with metrics (delta12, entropy, ENA) for any chromosome.
#
# Design:
#   1. Check if precomputed local Q exists on disk → load and return (fast)
#   2. If not, run the local Q projector (fixed-F EM) → save → return
#   3. Downstream consumers (precomp, snakes, scoring) just call get_local_Q()
#      and don't care how Q was produced
#
# This is the ONLY interface between the ancestry modules and the inversion
# pipeline. All Q access goes through here. If we later upgrade to true
# Module 2B windowed NGSadmix, we change this file and nothing else.
#
# Precomputed files are saved per chromosome so we never recompute unless
# the user explicitly requests it with --force.
#
# Usage (standalone):
#   Rscript ancestry_bridge.R --prepare \
#     --beagle_dir <dir> --qopt <file> --fopt <file> \
#     --window_registry <file> --outdir <dir> \
#     [--chrom <chr>] [--em_iter 20] [--ncores 4]
#
# Usage (as sourced library in precomp or any R script):
#   source("ancestry_bridge.R")
#   q_dt <- get_local_Q(chr = "C_gar_LG01", cache_dir = "local_Q/")
#   # Returns data.table with: window_id, sample, Q1..QK, delta12, entropy, ena
#
#   q_summary <- get_local_Q_summary(chr = "C_gar_LG01", cache_dir = "local_Q/")
#   # Returns data.table with: window_id, mean_delta12, mean_entropy, mean_ena
#
# Cache structure:
#   <cache_dir>/
#     C_gar_LG01.local_Q_summary.tsv.gz    — per-window aggregated metrics
#     C_gar_LG01.local_Q_samples.tsv.gz    — per-window × per-sample Q + metrics
#     C_gar_LG01.local_Q_meta.tsv          — metadata (K, engine, n_windows, timestamp)
#     manifest.tsv                          — which chromosomes are cached
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
EPS_Q <- 1e-10

# =============================================================================
# CONFIGURATION (set once, used by all functions)
# =============================================================================

.ancestry_config <- new.env(parent = emptyenv())

configure_ancestry <- function(beagle_dir = NULL, qopt_file = NULL,
                                fopt_file = NULL, window_registry = NULL,
                                cache_dir = "local_Q", em_iter = 20L,
                                K = NULL) {
  .ancestry_config$beagle_dir <- beagle_dir
  .ancestry_config$qopt_file <- qopt_file
  .ancestry_config$fopt_file <- fopt_file
  .ancestry_config$window_registry <- window_registry
  .ancestry_config$cache_dir <- cache_dir
  .ancestry_config$em_iter <- em_iter

  # Load Q init
  if (!is.null(qopt_file) && file.exists(qopt_file)) {
    .ancestry_config$q_init <- as.matrix(fread(qopt_file, header = FALSE))
    .ancestry_config$K <- ncol(.ancestry_config$q_init)
    .ancestry_config$n_samples <- nrow(.ancestry_config$q_init)
    message("[ancestry] Q init loaded: ", .ancestry_config$n_samples,
            " samples × K=", .ancestry_config$K)
  }

  # Load F matrix if provided
  if (!is.null(fopt_file) && file.exists(fopt_file)) {
    .ancestry_config$f_mat <- as.matrix(fread(fopt_file, header = FALSE))
    message("[ancestry] F matrix loaded: ", nrow(.ancestry_config$f_mat),
            " sites × K=", ncol(.ancestry_config$f_mat))
  }

  if (!is.null(K)) .ancestry_config$K <- K

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

# =============================================================================
# CORE: LOCAL Q COMPUTATION (fixed-F EM, clean-room implementation)
# =============================================================================
# Published algorithm: Skotte, Korneliussen & Albrechtsen (2013) Genetics 195:693
# This implementation fixes F from genome-wide run, only updates Q per window.

compute_window_Q <- function(gl_win, f_win, q_prior, n_iter = 20L) {
  n_sites <- nrow(gl_win)
  n_ind <- nrow(q_prior)
  K_loc <- ncol(q_prior)

  if (n_sites < 5 || n_ind < 10) return(q_prior)

  Q <- q_prior
  prev_ll <- -Inf

  for (iter in seq_len(n_iter)) {
    prodA <- matrix(0, n_ind, K_loc)
    prodB <- matrix(0, n_ind, K_loc)
    ll <- 0

    for (j in seq_len(n_sites)) {
      f_j <- pmin(pmax(f_win[j, ], EPS_Q), 1 - EPS_Q)

      for (i in seq_len(n_ind)) {
        p <- sum(Q[i, ] * f_j)
        p <- pmin(pmax(p, EPS_Q), 1 - EPS_Q)

        pp0 <- (1 - p)^2     * gl_win[j, 3*(i-1) + 1]
        pp1 <- 2 * p * (1-p) * gl_win[j, 3*(i-1) + 2]
        pp2 <- p^2           * gl_win[j, 3*(i-1) + 3]
        denom <- pp0 + pp1 + pp2
        if (denom < EPS_Q) next

        ll <- ll + log(denom)
        exp_g <- (pp1 + 2 * pp2) / denom

        for (k in seq_len(K_loc)) {
          prodA[i, k] <- prodA[i, k] + exp_g / (1 - p) * f_j[k]
          prodB[i, k] <- prodB[i, k] + (2 - exp_g) / p * (1 - f_j[k])
        }
      }
    }

    # Q update
    for (i in seq_len(n_ind)) {
      for (k in seq_len(K_loc)) {
        Q[i, k] <- Q[i, k] * (prodA[i, k] + prodB[i, k])
      }
      rs <- sum(Q[i, ])
      if (rs > EPS_Q) Q[i, ] <- Q[i, ] / rs
      Q[i, ] <- pmin(pmax(Q[i, ], EPS_Q), 1 - EPS_Q)
      Q[i, ] <- Q[i, ] / sum(Q[i, ])
    }

    # Convergence check
    if (abs(ll - prev_ll) < 1e-4 * abs(prev_ll) && iter > 3) break
    prev_ll <- ll
  }
  Q
}

# Estimate F from Q + GL (one M-step, used when F matrix not provided)
estimate_F_local <- function(gl_win, Q) {
  n_sites <- nrow(gl_win)
  n_ind <- nrow(Q)
  K_loc <- ncol(Q)
  F_est <- matrix(0.25, n_sites, K_loc)

  for (j in seq_len(n_sites)) {
    numer <- rep(0, K_loc)
    denom_f <- rep(0, K_loc)
    for (i in seq_len(n_ind)) {
      p <- sum(Q[i, ] * F_est[j, ])
      p <- pmin(pmax(p, EPS_Q), 1 - EPS_Q)

      pp0 <- (1-p)^2 * gl_win[j, 3*(i-1)+1]
      pp1 <- 2*p*(1-p) * gl_win[j, 3*(i-1)+2]
      pp2 <- p^2 * gl_win[j, 3*(i-1)+3]
      d <- pp0 + pp1 + pp2
      if (d < EPS_Q) next
      exp_g <- (pp1 + 2*pp2) / d

      for (k in seq_len(K_loc)) {
        numer[k] <- numer[k] + Q[i, k] * exp_g
        denom_f[k] <- denom_f[k] + 2 * Q[i, k]
      }
    }
    for (k in seq_len(K_loc)) {
      F_est[j, k] <- if (denom_f[k] > EPS_Q) numer[k] / denom_f[k] else 0.25
      F_est[j, k] <- pmin(pmax(F_est[j, k], EPS_Q), 1 - EPS_Q)
    }
  }
  F_est
}

# Per-sample Q metrics
q_sample_metrics <- function(q_vec) {
  q_vec <- pmax(q_vec, EPS_Q)
  qs <- sort(q_vec, decreasing = TRUE)
  H <- -sum(q_vec * log(q_vec))
  list(
    max_q = qs[1],
    delta12 = qs[1] - qs[2],
    entropy = H,
    ena = exp(H),
    assigned_pop = which.max(q_vec)
  )
}

# =============================================================================
# PUBLIC API: get_local_Q() and get_local_Q_summary()
# =============================================================================

#' Get per-window × per-sample local Q for a chromosome
#' Returns cached result if available, otherwise computes and caches.
#'
#' @param chr Chromosome name (e.g. "C_gar_LG01")
#' @param cache_dir Directory for cached files (default from configure_ancestry)
#' @param force Recompute even if cache exists
#' @return data.table with: window_id, chrom, start_bp, end_bp, sample_idx,
#'         Q1..QK, max_q, delta12, entropy, ena, assigned_pop

get_local_Q <- function(chr, cache_dir = NULL, force = FALSE) {
  cache_dir <- cache_dir %||% .ancestry_config$cache_dir %||% "local_Q"
  cache_file <- file.path(cache_dir, paste0(chr, ".local_Q_samples.tsv.gz"))

  # Check cache
  if (file.exists(cache_file) && !force) {
    dt <- fread(cache_file)
    message("[ancestry] Loaded cached local Q for ", chr, ": ", nrow(dt), " rows")
    return(dt)
  }

  # Need to compute
  message("[ancestry] Computing local Q for ", chr, " ...")
  q_init <- .ancestry_config$q_init
  f_mat <- .ancestry_config$f_mat
  beagle_dir <- .ancestry_config$beagle_dir
  K <- .ancestry_config$K
  n_iter <- .ancestry_config$em_iter %||% 20L

  if (is.null(q_init) || is.null(beagle_dir)) {
    stop("[ancestry] Not configured. Call configure_ancestry() first.")
  }

  # Load BEAGLE for this chromosome
  bgl_files <- list.files(beagle_dir, pattern = paste0(".*", chr, ".*beagle"),
                           full.names = TRUE)
  if (length(bgl_files) == 0) {
    message("[ancestry] No BEAGLE found for ", chr, " in ", beagle_dir)
    return(data.table())
  }
  bgl <- fread(bgl_files[1], nThread = 4L)
  markers <- bgl[[1]]
  positions <- as.integer(sub(".*_", "", markers))
  gl_mat <- as.matrix(bgl[, 4:ncol(bgl)])
  n_sites <- nrow(gl_mat)
  n_ind <- ncol(gl_mat) %/% 3

  message("[ancestry]   ", n_sites, " sites × ", n_ind, " samples")

  # Load window registry or auto-generate
  win_reg <- .ancestry_config$window_registry
  if (!is.null(win_reg) && file.exists(win_reg)) {
    wins <- fread(win_reg)
    wins <- wins[chrom == chr]
  } else {
    # Auto: 100 SNPs, step 20
    w_size <- 100L; w_step <- 20L
    n_wins <- max(1, (n_sites - w_size) %/% w_step + 1)
    wins <- data.table(
      global_window_id = paste0(chr, "_w", seq_len(n_wins)),
      chrom = chr,
      start_idx = seq(1, by = w_step, length.out = n_wins)
    )
    wins[, end_idx := pmin(start_idx + w_size - 1L, n_sites)]
    wins[, start_bp := positions[start_idx]]
    wins[, end_bp := positions[end_idx]]
  }

  message("[ancestry]   ", nrow(wins), " windows")
  n_use <- min(n_ind, nrow(q_init))

  all_rows <- list()
  t0 <- proc.time()

  for (wi in seq_len(nrow(wins))) {
    w <- wins[wi]
    wid <- w$global_window_id

    if ("start_idx" %in% names(w)) {
      site_idx <- w$start_idx:w$end_idx
    } else {
      site_idx <- which(positions >= w$start_bp & positions <= w$end_bp)
    }
    if (length(site_idx) < 10) next

    gl_win <- gl_mat[site_idx, seq_len(3 * n_use), drop = FALSE]
    q_prior <- q_init[seq_len(n_use), , drop = FALSE]

    # Get F for this window
    if (!is.null(f_mat) && max(site_idx) <= nrow(f_mat)) {
      f_win <- f_mat[site_idx, , drop = FALSE]
    } else {
      f_win <- estimate_F_local(gl_win, q_prior)
    }

    # Compute local Q
    q_local <- compute_window_Q(gl_win, f_win, q_prior, n_iter = n_iter)

    # Per-sample output
    for (si in seq_len(n_use)) {
      m <- q_sample_metrics(q_local[si, ])
      row <- list(
        window_id = wid, chrom = chr,
        start_bp = positions[site_idx[1]],
        end_bp = positions[tail(site_idx, 1)],
        sample_idx = si
      )
      # Add Q columns
      for (k in seq_len(K)) row[[paste0("Q", k)]] <- round(q_local[si, k], 4)
      row$max_q <- round(m$max_q, 4)
      row$delta12 <- round(m$delta12, 4)
      row$entropy <- round(m$entropy, 4)
      row$ena <- round(m$ena, 4)
      row$assigned_pop <- m$assigned_pop
      all_rows[[length(all_rows) + 1L]] <- row
    }

    if (wi %% 200 == 0) {
      elapsed <- round((proc.time() - t0)[3], 1)
      rate <- round(wi / max(elapsed, 0.1), 1)
      message("[ancestry]   window ", wi, "/", nrow(wins), " (", rate, " win/s)")
    }
  }

  result <- rbindlist(all_rows, fill = TRUE)

  # Save cache
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(result, cache_file, sep = "\t")

  # Save summary
  summary_dt <- result[, .(
    mean_delta12 = round(mean(delta12, na.rm = TRUE), 4),
    mean_entropy = round(mean(entropy, na.rm = TRUE), 4),
    mean_ena = round(mean(ena, na.rm = TRUE), 4),
    sd_delta12 = round(sd(delta12, na.rm = TRUE), 4),
    n_samples = .N
  ), by = .(window_id, chrom, start_bp, end_bp)]

  summary_file <- file.path(cache_dir, paste0(chr, ".local_Q_summary.tsv.gz"))
  fwrite(summary_dt, summary_file, sep = "\t")

  # Save metadata
  meta_file <- file.path(cache_dir, paste0(chr, ".local_Q_meta.tsv"))
  fwrite(data.table(
    chrom = chr, K = K, engine = "local_projector_fixed_F",
    em_iter = n_iter, n_windows = nrow(wins), n_samples = n_use,
    n_sites_total = n_sites,
    qopt_source = .ancestry_config$qopt_file %||% "NA",
    fopt_source = .ancestry_config$fopt_file %||% "estimated",
    timestamp = Sys.time()
  ), meta_file, sep = "\t")

  # Update manifest
  manifest_file <- file.path(cache_dir, "manifest.tsv")
  manifest_row <- data.table(chrom = chr, cache_file = cache_file,
                              summary_file = summary_file,
                              n_windows = nrow(wins), timestamp = Sys.time())
  if (file.exists(manifest_file)) {
    old <- fread(manifest_file)
    old <- old[chrom != chr]  # replace existing
    fwrite(rbind(old, manifest_row), manifest_file, sep = "\t")
  } else {
    fwrite(manifest_row, manifest_file, sep = "\t")
  }

  elapsed <- round((proc.time() - t0)[3], 1)
  message("[ancestry] ", chr, " done: ", nrow(result), " rows, ",
          nrow(summary_dt), " windows, ", elapsed, "s")
  result
}

#' Get per-window summary (no per-sample detail)
get_local_Q_summary <- function(chr, cache_dir = NULL, force = FALSE) {
  cache_dir <- cache_dir %||% .ancestry_config$cache_dir %||% "local_Q"
  summary_file <- file.path(cache_dir, paste0(chr, ".local_Q_summary.tsv.gz"))

  if (file.exists(summary_file) && !force) {
    return(fread(summary_file))
  }

  # Compute full then return summary
  full <- get_local_Q(chr, cache_dir, force)
  if (nrow(full) == 0) return(data.table())

  full[, .(
    mean_delta12 = round(mean(delta12, na.rm = TRUE), 4),
    mean_entropy = round(mean(entropy, na.rm = TRUE), 4),
    mean_ena = round(mean(ena, na.rm = TRUE), 4),
    sd_delta12 = round(sd(delta12, na.rm = TRUE), 4),
    n_samples = .N
  ), by = .(window_id, chrom, start_bp, end_bp)]
}

#' List cached chromosomes
list_cached_Q <- function(cache_dir = NULL) {
  cache_dir <- cache_dir %||% .ancestry_config$cache_dir %||% "local_Q"
  manifest_file <- file.path(cache_dir, "manifest.tsv")
  if (file.exists(manifest_file)) {
    fread(manifest_file)
  } else {
    # Scan directory
    files <- list.files(cache_dir, pattern = "\\.local_Q_summary\\.tsv\\.gz$")
    data.table(chrom = sub("\\.local_Q_summary\\.tsv\\.gz$", "", files))
  }
}

# =============================================================================
# INTEGRATION WITH PRECOMP
# =============================================================================
# Call this from STEP_C01a_snake1_precompute.R to add real Q-based metrics
# to the inv_likeness table.

merge_local_Q_into_invlikeness <- function(inv_like_dt, cache_dir = "local_Q") {
  manifest <- list_cached_Q(cache_dir)
  if (nrow(manifest) == 0) {
    message("[ancestry] No cached local Q found — skipping merge")
    return(inv_like_dt)
  }

  message("[ancestry] Merging local Q from ", nrow(manifest), " chromosomes")

  for (chr in unique(inv_like_dt$chrom)) {
    summary_file <- file.path(cache_dir, paste0(chr, ".local_Q_summary.tsv.gz"))
    if (!file.exists(summary_file)) next

    q_summ <- fread(summary_file)
    if (nrow(q_summ) == 0) next

    # Match by window_id or by position overlap
    if ("window_id" %in% names(q_summ) && "global_window_id" %in% names(inv_like_dt)) {
      q_summ[, global_window_id := window_id]
      q_sub <- q_summ[, .(global_window_id,
                           localQ_delta12 = mean_delta12,
                           localQ_entropy = mean_entropy,
                           localQ_ena = mean_ena)]
      inv_like_dt <- merge(inv_like_dt, q_sub, by = "global_window_id",
                            all.x = TRUE, suffixes = c("", ".q"))
    } else {
      # Position-based matching (slower)
      for (wi in which(inv_like_dt$chrom == chr)) {
        w_start <- inv_like_dt$start_bp[wi]
        w_end <- inv_like_dt$end_bp[wi]
        match_idx <- which(q_summ$start_bp >= w_start - 500 &
                            q_summ$end_bp <= w_end + 500)
        if (length(match_idx) > 0) {
          best <- match_idx[1]
          set(inv_like_dt, wi, "localQ_delta12", q_summ$mean_delta12[best])
          set(inv_like_dt, wi, "localQ_entropy", q_summ$mean_entropy[best])
          set(inv_like_dt, wi, "localQ_ena", q_summ$mean_ena[best])
        }
      }
    }
  }

  n_filled <- sum(is.finite(inv_like_dt$localQ_delta12))
  message("[ancestry] Merged: ", n_filled, " / ", nrow(inv_like_dt), " windows have local Q")
  inv_like_dt
}

# =============================================================================
# STANDALONE MODE
# =============================================================================

if (sys.nframe() == 0) {
  # Running as standalone script
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || args[1] == "--help") {
    cat("Usage: Rscript ancestry_bridge.R --prepare \\
    --beagle_dir <dir> --qopt <file> [--fopt <file>] \\
    --outdir <dir> [--chrom <chr>] [--em_iter 20] [--window_registry <file>]\n")
    quit(status = 0)
  }

  beagle_dir <- NULL; qopt <- NULL; fopt <- NULL
  outdir <- "local_Q"; chrom <- NULL; em_iter <- 20L
  win_reg <- NULL

  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--beagle_dir") { beagle_dir <- args[i+1]; i <- i+2 }
    else if (a == "--qopt") { qopt <- args[i+1]; i <- i+2 }
    else if (a == "--fopt") { fopt <- args[i+1]; i <- i+2 }
    else if (a == "--outdir") { outdir <- args[i+1]; i <- i+2 }
    else if (a == "--chrom") { chrom <- args[i+1]; i <- i+2 }
    else if (a == "--em_iter") { em_iter <- as.integer(args[i+1]); i <- i+2 }
    else if (a == "--window_registry") { win_reg <- args[i+1]; i <- i+2 }
    else if (a == "--prepare") { i <- i+1 }
    else { i <- i+1 }
  }

  configure_ancestry(beagle_dir = beagle_dir, qopt_file = qopt,
                      fopt_file = fopt, window_registry = win_reg,
                      cache_dir = outdir, em_iter = em_iter)

  if (!is.null(chrom)) {
    get_local_Q(chrom, cache_dir = outdir)
  } else {
    # All chromosomes from FAI or BEAGLE directory
    bgl_files <- list.files(beagle_dir, pattern = "\\.beagle", full.names = TRUE)
    for (bf in bgl_files) {
      chr_name <- sub(".*?(C_gar_LG[0-9]+).*", "\\1", basename(bf))
      if (grepl("^C_gar_LG", chr_name)) {
        get_local_Q(chr_name, cache_dir = outdir)
      }
    }
  }

  message("\n[ancestry] All done. Cache: ", outdir)
}
