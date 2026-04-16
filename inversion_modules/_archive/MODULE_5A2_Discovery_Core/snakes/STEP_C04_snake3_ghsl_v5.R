#!/usr/bin/env Rscript

# =============================================================================
# STEP_C04_snake3_ghsl_v5.R
#
# SNAKE 3 v5: Within-sample haplotype divergence + rank stability
#
# FUNDAMENTAL CHANGE FROM v4:
#   v4 asked: "Do fish in the same PC1 band have concordant haplotypes?"
#             (pairwise between-sample, 25K pairs per window, founder-noisy)
#   v5 asks:  "How different are each fish's OWN two haplotypes?"
#             (within-sample, 226 measurements per window, founder-immune)
#
# BIOLOGICAL BASIS:
#   INV/INV:     both haplotypes inverted → locked together → LOW divergence
#   INV/nonINV:  one inverted, one not → can't recombine → HIGH divergence
#   nonINV/nonINV (soup): random founder haplotypes → VARIABLE divergence
#
# DETECTION LOGIC:
#   At inversion windows: sample divergence ranks are STABLE across windows
#     (each fish stays in its lane: INV/INV always low, INV/nonINV always high)
#   At founder windows: ranks SHUFFLE randomly window to window
#     (divergence depends on which founder haplotypes, varies by position)
#
# STAGES:
#   0. LOAD + PRE-INDEX (same as v4: findInterval on merged phased SNPs)
#   1. PER-SAMPLE DIVERGENCE: within-sample hap1 vs hap2 divergence per window
#   2. RANK STABILITY: Spearman correlation of ranks between adjacent windows
#   3. KARYOTYPE CALLING: stable low-rank = INV/INV, stable high = INV/nonINV
#   4. WINDOW SCORING: rank stability + divergence modality + band concordance
#   5. BLOCK DETECTION: runs of high-scoring windows = inversion domains
#
# Usage:
#   Rscript STEP_C04_snake3_ghsl_v5.R <precomp_dir> <ghsl_prep_dir> <outdir> \
#     [--chrom C_gar_LG01] [--test_windows 50] [--ncores 4]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript STEP_C04_snake3_ghsl_v5.R <precomp_dir> <ghsl_prep_dir> <outdir> [opts]")

precomp_dir   <- args[1]
ghsl_prep_dir <- args[2]
outdir        <- args[3]

chrom_filter     <- NULL
QUAL_MIN         <- 20
GQ_MIN           <- 10
MIN_PHASED_SITES <- 3L      # min phased het sites per sample per window to compute divergence
TEST_WINDOWS     <- 0L       # 0 = all, >0 = first N
NCORES           <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
RANK_WINDOW      <- 5L       # number of adjacent windows for rank stability
KARYOTYPE_QUANTILE_LO <- 0.15  # bottom 15% divergence = candidate INV/INV
KARYOTYPE_QUANTILE_HI <- 0.70  # top 30% divergence = candidate INV/nonINV
KARYOTYPE_MIN_RUN     <- 10L   # min consecutive windows to call stable karyotype

i <- 4L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args))           { chrom_filter <- args[i+1]; i <- i+2L }
  else if (a == "--test_windows" && i < length(args)){ TEST_WINDOWS <- as.integer(args[i+1]); i <- i+2L }
  else if (a == "--ncores" && i < length(args))      { NCORES <- as.integer(args[i+1]); i <- i+2L }
  else if (a == "--rank_window" && i < length(args)) { RANK_WINDOW <- as.integer(args[i+1]); i <- i+2L }
  else if (a == "--min_phased" && i < length(args))  { MIN_PHASED_SITES <- as.integer(args[i+1]); i <- i+2L }
  else { i <- i+1L }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("================================================================")
message("[S3v5] Snake 3 v5: Within-Sample Haplotype Divergence")
message("================================================================")
message("[S3v5] Precomp:    ", precomp_dir)
message("[S3v5] GHSL prep:  ", ghsl_prep_dir)
message("[S3v5] Output:     ", outdir)
message("[S3v5] MIN_PHASED_SITES=", MIN_PHASED_SITES, " RANK_WINDOW=", RANK_WINDOW)
message("[S3v5] Karyotype quantiles: INV/INV < ", KARYOTYPE_QUANTILE_LO,
        " | INV/nonINV > ", KARYOTYPE_QUANTILE_HI)
if (TEST_WINDOWS > 0) message("[S3v5] TEST MODE: first ", TEST_WINDOWS, " windows only")

# =============================================================================
# LOAD PRECOMP (for window grid and sample names)
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("[S3v5] FATAL: No .precomp.rds files in: ", precomp_dir)
message("[S3v5] Found ", length(rds_files), " precomp RDS files")

precomp_list <- list()
for (f in rds_files) { obj <- readRDS(f); precomp_list[[obj$chrom]] <- obj }
chroms <- names(precomp_list)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)

# Get sample names
sample_names <- NULL
precomp_sample_names <- NULL
for (chr_tmp in chroms) {
  pc1_cols <- grep("^PC_1_", names(precomp_list[[chr_tmp]]$dt), value = TRUE)
  if (length(pc1_cols) > 0) { precomp_sample_names <- sub("^PC_1_", "", pc1_cols); break }
}
if (is.null(precomp_sample_names)) stop("[S3v5] FATAL: No PC_1_ columns found")
n_samples <- length(precomp_sample_names)
sample_names <- precomp_sample_names

# Map Ind→CGA
if (grepl("^Ind[0-9]", sample_names[1])) {
  sf <- Sys.getenv("SAMPLES_IND", "")
  if (!nzchar(sf) || !file.exists(sf)) {
    base <- Sys.getenv("BASE", "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04")
    for (candidate in c(
      file.path(base, "het_roh/01_inputs_check/samples.ind"),
      file.path(base, "popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv")
    )) {
      if (file.exists(candidate)) { sf <- candidate; break }
    }
  }
  if (nzchar(sf) && file.exists(sf)) {
    real <- trimws(readLines(sf))
    real <- real[nzchar(real)]
    if (length(real) == n_samples) {
      sample_names <- real
      message("[S3v5] Sample mapping: ", precomp_sample_names[1], " -> ", sample_names[1],
              " (", n_samples, " samples)")
    }
  }
}

message("[S3v5] Samples: ", n_samples)
message("[S3v5] Chromosomes: ", length(chroms))

# =============================================================================
# STAGE 1: Compute per-sample within-haplotype divergence
# =============================================================================

compute_divergence_matrix <- function(ghsl_by_window, dt, n_win, available_samples) {
  # Returns TWO divergence matrices:
  #
  # 1. ghsl_div: GHSL-style haplotype divergence (primary)
  #    = n_phased_het (0|1 or 1|0) / n_total_variants_in_window
  #    Only phased hets count as divergent (we KNOW hap1 ≠ hap2).
  #    Denominator = ALL variants for this sample in this window
  #    (phased hets + unphased hets + hom_var).
  #    INV/INV: low (haplotypes locked together → few hets)
  #    INV/nonINV: high (haplotypes can't recombine → many hets)
  #
  # 2. het_div: simple het rate (secondary)
  #    = n_all_het / n_total_variants_in_window
  #    Includes unphased hets too. Less precise but more data.
  #
  # Both normalized by total sites per sample per window (sequence length proxy).

  n_samp <- length(available_samples)
  ghsl_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                     dimnames = list(available_samples, NULL))
  het_mat  <- matrix(NA_real_, nrow = n_samp, ncol = n_win,
                     dimnames = list(available_samples, NULL))
  n_sites_mat <- matrix(0L, nrow = n_samp, ncol = n_win,
                        dimnames = list(available_samples, NULL))
  n_phased_het_mat <- matrix(0L, nrow = n_samp, ncol = n_win,
                             dimnames = list(available_samples, NULL))

  for (wi in seq_len(n_win)) {
    win_dt <- ghsl_by_window[[as.character(wi)]]
    if (is.null(win_dt) || nrow(win_dt) == 0) next

    for (si in available_samples) {
      sv <- win_dt[sample_id == si]
      n_total <- nrow(sv)
      if (n_total < MIN_PHASED_SITES) next

      # Count by genotype class
      gt <- tolower(sv$gt_class)
      n_het_all <- sum(gt == "het")
      n_hom <- sum(gt %in% c("hom_ref", "hom_alt", "hom_var"))

      # Phased hets only: phase_gt contains "|" (0|1 or 1|0)
      # These are the sites where we KNOW hap1 ≠ hap2
      phased_het <- sv[gt == "het" & grepl("\\|", phase_gt)]
      n_phased_het <- nrow(phased_het)

      # GHSL divergence: phased hets / total sites
      # This is the "sequence length where haplotypes diverge" normalized by window
      ghsl_mat[si, wi] <- n_phased_het / n_total
      
      # Het rate divergence: all hets / total sites
      het_mat[si, wi] <- n_het_all / n_total

      n_sites_mat[si, wi] <- n_total
      n_phased_het_mat[si, wi] <- n_phased_het
    }

    if (wi %% 500 == 0) {
      n_scored <- sum(!is.na(ghsl_mat[, wi]))
      # Quick diagnostic: show range of ghsl values
      vals <- ghsl_mat[!is.na(ghsl_mat[, wi]), wi]
      msg_range <- if (length(vals) > 0) {
        paste0(" ghsl=[", round(min(vals), 3), "-", round(max(vals), 3),
               "] het=[", round(min(het_mat[!is.na(het_mat[, wi]), wi]), 3), "-",
               round(max(het_mat[!is.na(het_mat[, wi]), wi]), 3), "]")
      } else ""
      message("[S3v5]   Divergence: window ", wi, "/", n_win,
              " (", n_scored, " samples scored)", msg_range)
    }
  }

  list(ghsl = ghsl_mat, het = het_mat, n_sites = n_sites_mat,
       n_phased_het = n_phased_het_mat)
}

# =============================================================================
# STAGE 2: Rank stability across adjacent windows
# =============================================================================

compute_rank_stability <- function(div_mat, window_span = 5L) {
  n_samples <- nrow(div_mat)
  n_win <- ncol(div_mat)
  
  # Per-window: rank samples by divergence (1 = lowest = candidate INV/INV)
  rank_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_win,
                     dimnames = list(rownames(div_mat), NULL))
  for (wi in seq_len(n_win)) {
    vals <- div_mat[, wi]
    valid <- !is.na(vals)
    if (sum(valid) < 10) next
    rank_mat[valid, wi] <- rank(vals[valid], ties.method = "average") / sum(valid)
    # Normalized rank: 0 = lowest divergence, 1 = highest
  }

  # Rank stability: Spearman correlation between window w and w+1
  pairwise_rho <- numeric(n_win)
  pairwise_rho[] <- NA_real_

  for (wi in seq_len(n_win - 1)) {
    r1 <- rank_mat[, wi]
    r2 <- rank_mat[, wi + 1]
    valid <- !is.na(r1) & !is.na(r2)
    if (sum(valid) < 20) next
    pairwise_rho[wi] <- cor(r1[valid], r2[valid], method = "spearman")
  }

  # Smoothed rank stability: mean rho over a sliding window of width window_span
  smooth_rho <- numeric(n_win)
  smooth_rho[] <- NA_real_
  half <- floor(window_span / 2)

  for (wi in seq_len(n_win)) {
    lo <- max(1, wi - half)
    hi <- min(n_win - 1, wi + half)
    rhos <- pairwise_rho[lo:hi]
    rhos <- rhos[!is.na(rhos)]
    if (length(rhos) >= 2) smooth_rho[wi] <- mean(rhos)
  }

  list(rank_mat = rank_mat, pairwise_rho = pairwise_rho, smooth_rho = smooth_rho)
}

# =============================================================================
# STAGE 3: Per-window divergence distribution metrics
# =============================================================================

compute_window_metrics <- function(div_mat, rank_mat) {
  n_win <- ncol(div_mat)

  metrics <- data.table(
    window_idx    = seq_len(n_win),
    n_scored      = integer(n_win),
    div_mean      = numeric(n_win),
    div_median    = numeric(n_win),
    div_sd        = numeric(n_win),
    div_iqr       = numeric(n_win),
    div_skew      = numeric(n_win),   # skewness: bimodal inversions show specific shapes
    div_bimodal   = numeric(n_win),   # dip test or similar
    n_low_div     = integer(n_win),   # samples in bottom quantile (INV/INV candidates)
    n_high_div    = integer(n_win),   # samples in top quantile (INV/nonINV candidates)
    n_mid_div     = integer(n_win),   # samples in middle (soup)
    low_div_tightness  = numeric(n_win),  # IQR within low-div group
    high_div_tightness = numeric(n_win),  # IQR within high-div group
    mid_div_spread     = numeric(n_win)   # IQR within mid group
  )

  for (wi in seq_len(n_win)) {
    vals <- div_mat[, wi]
    valid <- !is.na(vals)
    n_valid <- sum(valid)
    metrics$n_scored[wi] <- n_valid
    if (n_valid < 20) next

    v <- vals[valid]
    metrics$div_mean[wi]   <- mean(v)
    metrics$div_median[wi] <- median(v)
    metrics$div_sd[wi]     <- sd(v)
    metrics$div_iqr[wi]    <- IQR(v)

    # Skewness (moment-based)
    m3 <- mean((v - mean(v))^3)
    s3 <- sd(v)^3
    metrics$div_skew[wi] <- if (s3 > 0) m3 / s3 else 0

    # Bimodality: Hartigan's dip test (if available) or simple check
    # Simple bimodality: ratio of density at median vs density at modes
    # Use kernel density
    if (n_valid >= 30 && sd(v) > 0.001) {
      d <- tryCatch({
        dens <- density(v, n = 64)
        # Count local maxima
        dv <- diff(sign(diff(dens$y)))
        n_modes <- sum(dv == -2) + 1  # peaks
        n_modes
      }, error = function(e) 1)
      metrics$div_bimodal[wi] <- d
    } else {
      metrics$div_bimodal[wi] <- 1
    }

    # Quantile-based groups
    q_lo <- quantile(v, KARYOTYPE_QUANTILE_LO)
    q_hi <- quantile(v, KARYOTYPE_QUANTILE_HI)

    low_mask  <- v <= q_lo
    high_mask <- v >= q_hi
    mid_mask  <- v > q_lo & v < q_hi

    metrics$n_low_div[wi]  <- sum(low_mask)
    metrics$n_high_div[wi] <- sum(high_mask)
    metrics$n_mid_div[wi]  <- sum(mid_mask)

    # Tightness within each group (low IQR = samples agree on divergence level)
    metrics$low_div_tightness[wi]  <- if (sum(low_mask) >= 5) IQR(v[low_mask]) else NA_real_
    metrics$high_div_tightness[wi] <- if (sum(high_mask) >= 5) IQR(v[high_mask]) else NA_real_
    metrics$mid_div_spread[wi]     <- if (sum(mid_mask) >= 5) IQR(v[mid_mask]) else NA_real_
  }

  metrics
}

# =============================================================================
# STAGE 4: Karyotype calling per sample
# =============================================================================

call_karyotypes <- function(rank_mat, min_run = 10L) {
  # For each sample, find runs of consecutive windows where rank is stable
  # and in a consistent quantile (low = INV/INV, high = INV/nonINV)
  n_samples <- nrow(rank_mat)
  n_win <- ncol(rank_mat)
  snames <- rownames(rank_mat)

  karyo_calls <- list()

  for (si in seq_len(n_samples)) {
    ranks <- rank_mat[si, ]

    # Classify each window: low / high / mid / NA
    state <- rep("NA", n_win)
    state[!is.na(ranks) & ranks <= KARYOTYPE_QUANTILE_LO] <- "LOW"   # INV/INV candidate
    state[!is.na(ranks) & ranks >= KARYOTYPE_QUANTILE_HI] <- "HIGH"  # INV/nonINV candidate
    state[!is.na(ranks) & ranks > KARYOTYPE_QUANTILE_LO & ranks < KARYOTYPE_QUANTILE_HI] <- "MID"

    # Find runs of LOW or HIGH
    rle_state <- rle(state)
    run_ends <- cumsum(rle_state$lengths)
    run_starts <- c(1L, head(run_ends, -1) + 1L)

    for (ri in seq_along(rle_state$values)) {
      if (rle_state$values[ri] %in% c("LOW", "HIGH") && rle_state$lengths[ri] >= min_run) {
        karyo_calls[[length(karyo_calls) + 1L]] <- data.table(
          sample_id = snames[si],
          window_start = run_starts[ri],
          window_end = run_ends[ri],
          n_windows = rle_state$lengths[ri],
          call = if (rle_state$values[ri] == "LOW") "INV_INV" else "INV_nonINV",
          mean_rank = mean(ranks[run_starts[ri]:run_ends[ri]], na.rm = TRUE)
        )
      }
    }
  }

  if (length(karyo_calls) > 0) rbindlist(karyo_calls) else data.table()
}

# =============================================================================
# STAGE 5: Window-level inversion scoring
# =============================================================================

score_windows <- function(metrics, smooth_rho, div_mat) {
  n_win <- nrow(metrics)

  # ── Step A: raw rank stability ──
  metrics[, rank_stability := smooth_rho[seq_len(n_win)]]

  # ── Step B: divergence contrast (low-div group mean vs high-div group mean) ──
  lo_means <- numeric(n_win)
  hi_means <- numeric(n_win)
  for (wi in seq_len(n_win)) {
    vals <- div_mat[, wi]
    valid <- !is.na(vals)
    if (sum(valid) < 20) next
    v <- vals[valid]
    q_lo <- quantile(v, KARYOTYPE_QUANTILE_LO)
    q_hi <- quantile(v, KARYOTYPE_QUANTILE_HI)
    lo_vals <- v[v <= q_lo]
    hi_vals <- v[v >= q_hi]
    lo_means[wi] <- if (length(lo_vals) > 0) mean(lo_vals) else NA_real_
    hi_means[wi] <- if (length(hi_vals) > 0) mean(hi_vals) else NA_real_
  }
  metrics[, div_contrast := hi_means - lo_means]

  # ── Step C: Z-score rank stability against chromosome-wide baseline ──
  # KEY INSIGHT: In a hatchery population, founder structure inflates rank stability
  # everywhere (median ~0.78). Raw rank_stability > 0.5 hits 100% of windows.
  # The discriminator is: how far ABOVE the chromosome baseline is this window?
  # Only windows with rank stability significantly above the founder noise floor
  # should get PASS.
  rho_valid <- metrics$rank_stability[!is.na(metrics$rank_stability)]
  rho_median <- if (length(rho_valid) > 10) median(rho_valid) else 0.5
  rho_mad    <- if (length(rho_valid) > 10) mad(rho_valid, constant = 1.4826) else 0.1
  # Floor MAD at 0.02 to avoid division by near-zero in highly uniform chromosomes
  rho_mad <- max(rho_mad, 0.02)

  metrics[, rank_stability_z := fifelse(
    !is.na(rank_stability),
    (rank_stability - rho_median) / rho_mad,
    NA_real_
  )]

  # Also z-score div_contrast
  dc_valid <- metrics$div_contrast[!is.na(metrics$div_contrast) & metrics$n_scored >= 20]
  dc_median <- if (length(dc_valid) > 10) median(dc_valid) else 0
  dc_mad    <- if (length(dc_valid) > 10) mad(dc_valid, constant = 1.4826) else 0.1
  dc_mad <- max(dc_mad, 0.02)

  metrics[, div_contrast_z := fifelse(
    !is.na(div_contrast),
    (div_contrast - dc_median) / dc_mad,
    NA_real_
  )]

  message("[S3v5]   Scoring baseline: rho_median=", round(rho_median, 4),
          " rho_MAD=", round(rho_mad, 4),
          " contrast_median=", round(dc_median, 4),
          " contrast_MAD=", round(dc_mad, 4))

  # ── Step D: Combined score (z-score based) ──
  # Components:
  #   1. rank_stability_z (weight 0.40): above-baseline rank stability
  #   2. div_contrast_z   (weight 0.25): above-baseline divergence contrast
  #   3. bimodality bonus (weight 0.20): ≥2 density modes
  #   4. tightness bonus  (weight 0.15): INV/INV group is tight (IQR < 0.05)
  #
  # Z-scores are sigmoid-squashed to [0,1] via pnorm for stable weighting:
  #   z=0 → 0.5 (at baseline), z=2 → 0.977, z=-2 → 0.023

  metrics[, ghsl_v5_score := fifelse(
    !is.na(rank_stability_z) & !is.na(div_contrast_z) & n_scored >= 20,
    pnorm(rank_stability_z) * 0.40 +
    pnorm(div_contrast_z)   * 0.25 +
    fifelse(div_bimodal >= 2, 0.20, 0) +
    fifelse(!is.na(low_div_tightness) & low_div_tightness < 0.05, 0.15, 0),
    NA_real_
  )]

  # ── Step E: Status thresholds ──
  # With z-score normalization, the expected baseline score for a founder-noise
  # window is approximately: pnorm(0)*0.40 + pnorm(0)*0.25 + 0.20 + 0.15 = 0.675
  # (assuming bimodality and tightness bonuses fire, which they often do).
  # But for windows WITHOUT bimodality: 0.5*0.40 + 0.5*0.25 = 0.325
  # PASS requires being clearly above baseline in BOTH rank_stability and contrast.
  # Threshold 0.65 means: need z > ~0.5 on both components if bimodal,
  #   or z > ~2 on rank_stability alone.
  metrics[, ghsl_v5_status := fifelse(
    is.na(ghsl_v5_score), "FAIL",
    fifelse(ghsl_v5_score > 0.65 & rank_stability_z > 1.5, "PASS",
    fifelse(ghsl_v5_score > 0.50 & rank_stability_z > 0.5, "WEAK", "FAIL"))
  )]

  metrics
}

# =============================================================================
# MAIN LOOP
# =============================================================================

all_window_metrics <- list()
all_karyo_calls   <- list()
all_summary       <- list()

for (chr in chroms) {
  pc <- precomp_list[[chr]]
  if (is.null(pc) || pc$n_windows < 20) {
    message("[S3v5] SKIP ", chr, ": only ", pc$n_windows %||% 0, " windows")
    next
  }

  dt <- pc$dt
  n_win <- nrow(dt)
  message("\n================================================================")
  message("[S3v5] ======= ", chr, " (", n_win, " windows) =======")
  message("================================================================")

  # ── Load + pre-index (same as v4) ──
  ghsl_file <- file.path(ghsl_prep_dir, paste0(chr, ".merged_phased_snps.tsv.gz"))
  if (!file.exists(ghsl_file)) {
    message("[S3v5] SKIP ", chr, ": no merged phased SNPs file")
    next
  }

  t0 <- proc.time()
  message("[S3v5] Loading ", ghsl_file, " ...")
  ghsl_dt <- fread(ghsl_file)
  message("[S3v5]   Raw: ", formatC(nrow(ghsl_dt), big.mark = ","), " variants in ",
          round((proc.time() - t0)[3], 1), "s")

  # QC filter
  if ("qual" %in% names(ghsl_dt)) ghsl_dt <- ghsl_dt[is.na(qual) | qual >= QUAL_MIN]
  if ("gq" %in% names(ghsl_dt)) ghsl_dt <- ghsl_dt[is.na(gq) | gq >= GQ_MIN]
  message("[S3v5]   After QC: ", formatC(nrow(ghsl_dt), big.mark = ","))

  available_samples <- intersect(unique(ghsl_dt$sample_id), sample_names)
  message("[S3v5]   Samples: ", length(available_samples), " / ", n_samples)

  # Pre-index by window
  message("[S3v5]   Pre-indexing...")
  t_idx <- proc.time()
  ghsl_dt <- ghsl_dt[sample_id %in% available_samples]
  window_starts <- dt$start_bp
  window_ends <- dt$end_bp
  ghsl_dt[, window_id := findInterval(pos, window_starts)]
  ghsl_dt[window_id > n_win, window_id := n_win]
  ghsl_dt[window_id < 1, window_id := 1L]
  ghsl_dt <- ghsl_dt[pos <= window_ends[window_id]]
  message("[S3v5]   Indexed: ", formatC(nrow(ghsl_dt), big.mark = ","),
          " variants in ", round((proc.time() - t_idx)[3], 1), "s")

  # Split by window
  ghsl_by_window <- split(ghsl_dt, ghsl_dt$window_id)

  n_to_process <- if (TEST_WINDOWS > 0) min(TEST_WINDOWS, n_win) else n_win

  # ── STAGE 1: Compute divergence matrix ──
  message("[S3v5] Stage 1: Computing within-sample divergence matrix...")
  t1 <- proc.time()
  div_result <- compute_divergence_matrix(ghsl_by_window, dt, n_to_process, available_samples)
  div_mat <- div_result$ghsl       # PRIMARY: phased het / total (GHSL-style)
  het_mat <- div_result$het        # SECONDARY: all het / total
  n_sites_mat <- div_result$n_sites
  elapsed1 <- round((proc.time() - t1)[3], 1)

  n_scored_total <- sum(!is.na(div_mat))
  n_possible <- length(available_samples) * n_to_process
  message("[S3v5]   GHSL divergence computed: ", formatC(n_scored_total, big.mark = ","),
          " / ", formatC(n_possible, big.mark = ","),
          " (", round(100 * n_scored_total / n_possible, 1), "%) in ", elapsed1, "s")

  # Quick distribution check on first scored windows
  ghsl_vals <- div_mat[!is.na(div_mat)]
  if (length(ghsl_vals) > 0) {
    message("[S3v5]   GHSL range: min=", round(min(ghsl_vals), 4),
            " median=", round(median(ghsl_vals), 4),
            " max=", round(max(ghsl_vals), 4))
  }
  het_vals <- het_mat[!is.na(het_mat)]
  if (length(het_vals) > 0) {
    message("[S3v5]   Het rate range: min=", round(min(het_vals), 4),
            " median=", round(median(het_vals), 4),
            " max=", round(max(het_vals), 4))
  }

  # Per-sample coverage
  per_sample_scored <- rowSums(!is.na(div_mat))
  message("[S3v5]   Per-sample windows scored: median=", median(per_sample_scored),
          " mean=", round(mean(per_sample_scored)), " min=", min(per_sample_scored))

  # Per-window coverage
  per_window_scored <- colSums(!is.na(div_mat))
  message("[S3v5]   Per-window samples scored: median=", median(per_window_scored),
          " mean=", round(mean(per_window_scored)), " min=", min(per_window_scored))

  # ── STAGE 2: Rank stability ──
  message("[S3v5] Stage 2: Computing rank stability (window span=", RANK_WINDOW, ")...")
  t2 <- proc.time()
  rank_result <- compute_rank_stability(div_mat, window_span = RANK_WINDOW)
  elapsed2 <- round((proc.time() - t2)[3], 1)

  n_stable <- sum(!is.na(rank_result$smooth_rho) & rank_result$smooth_rho > 0.5, na.rm = TRUE)
  message("[S3v5]   Rank stability computed in ", elapsed2, "s")
  message("[S3v5]   Windows with smooth_rho > 0.5: ", n_stable, " / ", n_to_process,
          " (", round(100 * n_stable / n_to_process, 1), "%)")

  # Quick quantiles of pairwise rho
  rho_vals <- rank_result$pairwise_rho[!is.na(rank_result$pairwise_rho)]
  if (length(rho_vals) > 0) {
    message("[S3v5]   Pairwise rho: median=", round(median(rho_vals), 3),
            " q90=", round(quantile(rho_vals, 0.90), 3),
            " q95=", round(quantile(rho_vals, 0.95), 3))
  }

  # ── STAGE 3: Window metrics ──
  message("[S3v5] Stage 3: Computing per-window divergence metrics...")
  t3 <- proc.time()
  metrics <- compute_window_metrics(div_mat, rank_result$rank_mat)
  
  # Add het_mat summary per window (secondary metric)
  het_div_mean <- numeric(n_to_process)
  het_div_median <- numeric(n_to_process)
  for (wi in seq_len(n_to_process)) {
    hv <- het_mat[, wi]
    hv <- hv[!is.na(hv)]
    het_div_mean[wi] <- if (length(hv) > 0) mean(hv) else NA_real_
    het_div_median[wi] <- if (length(hv) > 0) median(hv) else NA_real_
  }
  metrics[, het_div_mean := het_div_mean]
  metrics[, het_div_median := het_div_median]
  elapsed3 <- round((proc.time() - t3)[3], 1)
  message("[S3v5]   Metrics computed in ", elapsed3, "s")

  # ── STAGE 4: Karyotype calling ──
  message("[S3v5] Stage 4: Calling per-sample karyotypes...")
  t4 <- proc.time()
  karyo_dt <- call_karyotypes(rank_result$rank_mat, min_run = KARYOTYPE_MIN_RUN)
  elapsed4 <- round((proc.time() - t4)[3], 1)

  if (nrow(karyo_dt) > 0) {
    n_inv_inv <- sum(karyo_dt$call == "INV_INV")
    n_inv_non <- sum(karyo_dt$call == "INV_nonINV")
    n_unique_samples <- length(unique(karyo_dt$sample_id))
    message("[S3v5]   Karyotype calls: ", nrow(karyo_dt), " runs (",
            n_inv_inv, " INV/INV, ", n_inv_non, " INV/nonINV) across ",
            n_unique_samples, " samples in ", elapsed4, "s")
  } else {
    message("[S3v5]   No stable karyotype runs found (min_run=", KARYOTYPE_MIN_RUN, ")")
  }

  # ── STAGE 5: Scoring ──
  message("[S3v5] Stage 5: Scoring windows...")
  metrics <- score_windows(metrics, rank_result$smooth_rho, div_mat)

  n_pass <- sum(metrics$ghsl_v5_status == "PASS", na.rm = TRUE)
  n_weak <- sum(metrics$ghsl_v5_status == "WEAK", na.rm = TRUE)
  n_fail <- sum(metrics$ghsl_v5_status == "FAIL", na.rm = TRUE)
  message("[S3v5]   PASS=", n_pass, " WEAK=", n_weak, " FAIL=", n_fail)

  # Add chromosome and position info
  metrics[, `:=`(
    chrom = chr,
    global_window_id = dt$global_window_id[seq_len(n_to_process)],
    start_bp = dt$start_bp[seq_len(n_to_process)],
    end_bp = dt$end_bp[seq_len(n_to_process)],
    pos_mb = round((dt$start_bp[seq_len(n_to_process)] + dt$end_bp[seq_len(n_to_process)]) / 2e6, 4)
  )]

  if (nrow(karyo_dt) > 0) {
    karyo_dt[, chrom := chr]
    # Map window indices to genomic positions
    karyo_dt[, `:=`(
      start_bp = dt$start_bp[window_start],
      end_bp = dt$end_bp[window_end],
      start_mb = round(dt$start_bp[window_start] / 1e6, 2),
      end_mb = round(dt$end_bp[window_end] / 1e6, 2)
    )]
  }

  all_window_metrics[[length(all_window_metrics) + 1L]] <- metrics
  all_karyo_calls[[length(all_karyo_calls) + 1L]] <- karyo_dt

  all_summary[[length(all_summary) + 1L]] <- data.table(
    chrom = chr, n_windows = n_to_process,
    n_pass = n_pass, n_weak = n_weak, n_fail = n_fail,
    n_karyo_calls = nrow(karyo_dt),
    median_rho = round(median(rho_vals), 4)
  )

  message("[S3v5] ", chr, " DONE: PASS=", n_pass, " WEAK=", n_weak, " FAIL=", n_fail)
}

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

message("\n[S3v5] Writing outputs...")

metrics_dt <- if (length(all_window_metrics) > 0) rbindlist(all_window_metrics, fill = TRUE) else data.table()
karyo_dt   <- if (length(all_karyo_calls) > 0) rbindlist(all_karyo_calls, fill = TRUE) else data.table()
summ_dt    <- if (length(all_summary) > 0) rbindlist(all_summary) else data.table()

# ── TSV outputs (human-readable) ──
f_track <- file.path(outdir, "snake3v5_window_track.tsv.gz")
fwrite(metrics_dt, f_track, sep = "\t")

f_karyo <- file.path(outdir, "snake3v5_karyotype_calls.tsv.gz")
fwrite(karyo_dt, f_karyo, sep = "\t")

f_summ <- file.path(outdir, "snake3v5_summary.tsv")
fwrite(summ_dt, f_summ, sep = "\t")

# PA matrix (for cross-module integration)
if (nrow(metrics_dt) > 0) {
  pa_dt <- metrics_dt[, .(chrom, global_window_id, start_bp, end_bp, pos_mb,
                           ghsl_v5_score, ghsl_v5_status, rank_stability,
                           rank_stability_z, div_contrast, div_contrast_z,
                           div_bimodal,
                           n_low_div, n_high_div,
                           het_div_mean, het_div_median)]
  f_pa <- file.path(outdir, "snake3v5_window_pa.tsv.gz")
  fwrite(pa_dt, f_pa, sep = "\t")
}

# ── Annotation RDS (per-chromosome, for merging into precomp RDS) ──
# These small RDS files contain GHSL v5 results keyed by global_window_id.
# Any downstream script can merge them into the precomp dt:
#   annot <- readRDS("C_gar_LG01.ghsl_v5.annot.rds")
#   pc$dt <- merge(pc$dt, annot, by = "global_window_id", all.x = TRUE)
annot_dir <- file.path(outdir, "annot")
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(metrics_dt) > 0) {
  for (chr in unique(metrics_dt$chrom)) {
    chr_metrics <- metrics_dt[chrom == chr]
    
    annot_dt <- chr_metrics[, .(
      global_window_id,
      ghsl_v5_score,
      ghsl_v5_status,
      ghsl_div_mean     = div_mean,      # GHSL divergence: phased het / total
      ghsl_div_median   = div_median,
      ghsl_div_sd       = div_sd,
      ghsl_div_iqr      = div_iqr,
      ghsl_div_bimodal  = div_bimodal,
      ghsl_het_mean     = het_div_mean,   # Het rate: all het / total
      ghsl_het_median   = het_div_median,
      ghsl_rank_stability = rank_stability,
      ghsl_rank_stability_z = rank_stability_z,
      ghsl_div_contrast = div_contrast,
      ghsl_div_contrast_z = div_contrast_z,
      ghsl_n_low_div    = n_low_div,      # N samples with low divergence (INV/INV candidates)
      ghsl_n_high_div   = n_high_div,     # N samples with high divergence (INV/nonINV candidates)
      ghsl_n_scored     = n_scored
    )]
    
    f_annot <- file.path(annot_dir, paste0(chr, ".ghsl_v5.annot.rds"))
    saveRDS(annot_dt, f_annot)
    message("[S3v5] Annotation RDS: ", f_annot, " (", nrow(annot_dt), " windows)")
  }
  
  # Also save karyotype calls per chromosome
  if (nrow(karyo_dt) > 0) {
    for (chr in unique(karyo_dt$chrom)) {
      chr_karyo <- karyo_dt[chrom == chr]
      f_karyo_annot <- file.path(annot_dir, paste0(chr, ".ghsl_v5.karyotypes.rds"))
      saveRDS(chr_karyo, f_karyo_annot)
      message("[S3v5] Karyotype RDS: ", f_karyo_annot, " (", nrow(chr_karyo), " calls)")
    }
  }
}

message("\n================================================================")
message("[DONE] Snake 3 v5.3: Within-Sample Haplotype Divergence (z-score scoring)")
message("================================================================")
message("  Track:       ", f_track)
message("  Karyotype:   ", f_karyo)
message("  Summary:     ", f_summ)
message("  Annotations: ", annot_dir, "/")

if (nrow(summ_dt) > 0) {
  message("\n  Totals: PASS=", sum(summ_dt$n_pass),
          " WEAK=", sum(summ_dt$n_weak),
          " FAIL=", sum(summ_dt$n_fail))
  if (sum(summ_dt$n_karyo_calls) > 0) {
    message("  Karyotype calls: ", sum(summ_dt$n_karyo_calls), " stable runs")
  }
}

message("\n  To merge into precomp RDS:")
message("    annot <- readRDS('", annot_dir, "/<chr>.ghsl_v5.annot.rds')")
message("    pc$dt <- merge(pc$dt, annot, by='global_window_id', all.x=TRUE)")

if (TEST_WINDOWS > 0) {
  message("\n  *** TEST MODE: only processed first ", TEST_WINDOWS, " windows ***")
}
