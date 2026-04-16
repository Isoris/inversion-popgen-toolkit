#!/usr/bin/env Rscript

# =============================================================================
# STEP10h_snake3_ghsl_haplotype_contrast_v2.R
#
# SNAKE 3 v2: Adapted for Clair3 short-read phasing (~5x Illumina).
#
# Changes from v1:
#   - ADAPTIVE THRESHOLDS: computed from actual block sizes, not hardcoded
#   - TIER-AWARE: separate GHSL for TIER_1 only vs ALL tiers
#   - DATA-DRIVEN CALIBRATION: PASS/WEAK thresholds from chromosome baseline
#   - BETTER QC FLAGS: tuned for 44-819 bp blocks, not 50kb blocks
#   - phase_tier column in input (new in preparator v2)
#
# Uses the SAME common dense focal grid as Snakes 1 and 2.
# Drop-in replacement for STEP10h v1 — same interface, same output format.
#
# INPUTS:
#   <step10_outprefix>.mds.rds — per_chr data with window coordinates
#   <phased_summary_dir>/<chr>.phased_het_summary.tsv.gz — from prepare_clair3_for_ghsl.py
#   <samples_ind_file> — sample names
#
# OUTPUTS:
#   snake3_sample_window.tsv.gz   — per-sample × per-window detail
#   snake3_track.tsv.gz           — per-window summary + QC
#   snake3_summary.tsv            — per-chromosome summary
#   snake3_thresholds.tsv         — computed adaptive thresholds per chr
#
# Usage:
#   Rscript STEP10h_snake3_ghsl_haplotype_contrast_v2.R \
#     <step10_outprefix> <phased_summary_dir> <samples_ind> <outdir>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript STEP10h_snake3_ghsl_haplotype_contrast_v2.R ",
       "<step10_outprefix> <phased_summary_dir> <samples_ind> <outdir>")
}

step10_prefix   <- args[1]
phased_dir      <- args[2]
samples_ind     <- args[3]
outdir          <- args[4]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# FIXED PARAMETERS (these rarely need changing)
# =============================================================================

MIN_PHASED_VARIANT_COUNT <- 3L   # min het variants in window for a sample to count
MIN_BLOCK_COUNT          <- 1L   # min blocks overlapping window
MIN_PASS_COUNT           <- 10L  # min samples passing QC for window-level PASS
TIER1_WEIGHT             <- 1.0  # weight for TIER_1 variants in weighted GHSL
TIER2_WEIGHT             <- 0.7  # weight for TIER_2 variants (less reliable)
INDEL_MULTIPLIER         <- 2.0  # weight for indel vs SNP in GHSL

# =============================================================================
# LOAD DATA
# =============================================================================

message("[STEP10h-v2] Loading STEP10 per-chromosome outputs...")
mds_rds_file <- paste0(step10_prefix, ".mds.rds")
if (!file.exists(mds_rds_file)) stop("Missing: ", mds_rds_file)
mds_obj <- readRDS(mds_rds_file)
per_chr <- mds_obj$per_chr

sample_names <- tryCatch({
  s <- fread(samples_ind, header = FALSE)[[1]]
  s[nchar(s) > 0]
}, error = function(e) NULL)
if (is.null(sample_names)) stop("Cannot read samples.ind: ", samples_ind)

chroms <- names(per_chr)
message("[STEP10h-v2] Chromosomes: ", length(chroms), ", Samples: ", length(sample_names))

# =============================================================================
# HELPER: Compute GHSL contrast for one sample x one window
# =============================================================================

compute_sample_contrast <- function(
  sample_id, window_id, window_start, window_end,
  sp_in_window,  # data.table subset for this sample in this window
  blocks_in_window  # data.table of unique blocks overlapping window
) {
  window_bp <- window_end - window_start

  # Blocks overlapping window (clipped to window boundaries)
  b_starts <- pmax(blocks_in_window$block_start, window_start)
  b_ends   <- pmin(blocks_in_window$block_end, window_end)
  phased_bp <- sum(pmax(b_ends - b_starts, 0))
  phased_fraction <- if (window_bp > 0) phased_bp / window_bp else 0
  block_count <- nrow(blocks_in_window)

  # Variant counts
  n_het_snp <- sum(sp_in_window$type == "snp")
  n_het_indel <- sum(sp_in_window$type == "indel")
  het_indel_bp <- sum(abs(sp_in_window[type == "indel"]$length))
  total_var <- n_het_snp + n_het_indel

  # Tier breakdown
  n_tier1 <- sum(sp_in_window$phase_tier == "TIER_1_WHATSHAP")
  n_tier2 <- sum(sp_in_window$phase_tier == "TIER_2_READPAIR")
  tier2_frac <- if (total_var > 0) n_tier2 / total_var else 0

  # Contrast rates per kb of phased sequence
  if (phased_bp > 0) {
    # Tier-weighted variant counts
    snp_t1 <- sum(sp_in_window$type == "snp" & sp_in_window$phase_tier == "TIER_1_WHATSHAP")
    snp_t2 <- sum(sp_in_window$type == "snp" & sp_in_window$phase_tier == "TIER_2_READPAIR")
    indel_t1 <- sum(sp_in_window$type == "indel" & sp_in_window$phase_tier == "TIER_1_WHATSHAP")
    indel_t2 <- sum(sp_in_window$type == "indel" & sp_in_window$phase_tier == "TIER_2_READPAIR")

    # Raw rates (per kb)
    snp_rate     <- n_het_snp / phased_bp * 1000
    indel_rate   <- n_het_indel / phased_bp * 1000
    indel_bp_rate <- het_indel_bp / phased_bp * 1000
    frag_index   <- block_count / phased_bp * 1000

    # Tier-weighted rates
    snp_rate_weighted   <- (snp_t1 * TIER1_WEIGHT + snp_t2 * TIER2_WEIGHT) / phased_bp * 1000
    indel_rate_weighted <- (indel_t1 * TIER1_WEIGHT + indel_t2 * TIER2_WEIGHT) / phased_bp * 1000

    # TIER_1-only rates (for comparison)
    snp_rate_t1   <- snp_t1 / phased_bp * 1000
    indel_rate_t1 <- indel_t1 / phased_bp * 1000
  } else {
    snp_rate <- indel_rate <- indel_bp_rate <- frag_index <- 0
    snp_rate_weighted <- indel_rate_weighted <- 0
    snp_rate_t1 <- indel_rate_t1 <- 0
  }

  ghsl_total    <- snp_rate + indel_rate
  ghsl_weighted <- snp_rate_weighted + INDEL_MULTIPLIER * indel_rate_weighted
  ghsl_t1_only  <- snp_rate_t1 + INDEL_MULTIPLIER * indel_rate_t1

  data.table(
    sample_id            = sample_id,
    window_id            = window_id,
    total_window_bp      = window_bp,
    phased_bp            = phased_bp,
    phased_fraction      = round(phased_fraction, 4),
    phased_variant_count = total_var,
    phased_het_snp_count = n_het_snp,
    phased_het_indel_count = n_het_indel,
    phased_het_indel_bp  = het_indel_bp,
    block_count          = block_count,
    n_tier1              = n_tier1,
    n_tier2              = n_tier2,
    tier2_fraction       = round(tier2_frac, 4),
    ghsl_snp_contrast    = round(snp_rate, 4),
    ghsl_indel_contrast  = round(indel_rate, 4),
    ghsl_indel_bp_contrast = round(indel_bp_rate, 4),
    ghsl_total           = round(ghsl_total, 4),
    ghsl_total_weighted  = round(ghsl_weighted, 4),
    ghsl_t1_only         = round(ghsl_t1_only, 4),
    fragmentation_index  = round(frag_index, 4)
  )
}

# =============================================================================
# MAIN: PER-CHROMOSOME PROCESSING
# =============================================================================

all_sample_window <- list()
all_track_rows    <- list()
all_summary_rows  <- list()
all_thresh_rows   <- list()

for (chr in chroms) {
  chr_obj <- per_chr[[chr]]
  if (is.null(chr_obj)) next

  dt <- as.data.table(chr_obj$out_dt)
  dt <- dt[order(start_bp)]
  n_win <- nrow(dt)

  message("\n[STEP10h-v2] ======= ", chr, " (", n_win, " windows) =======")

  if (n_win < 3) {
    message("[SKIP] ", chr, ": too few windows")
    next
  }

  # ── Load phased summary ──────────────────────────────────────────
  phased_file <- file.path(phased_dir, paste0(chr, ".phased_het_summary.tsv.gz"))
  if (!file.exists(phased_file)) {
    message("[SKIP] No phased summary for ", chr)

    for (wi in seq_len(n_win)) {
      all_track_rows[[length(all_track_rows) + 1]] <- data.table(
        chrom = chr, global_window_id = dt$global_window_id[wi],
        start_bp = dt$start_bp[wi], end_bp = dt$end_bp[wi],
        mean_ghsl = 0, median_ghsl = 0, max_ghsl = 0, spread_ghsl = 0,
        mean_snp_contrast = 0, mean_indel_contrast = 0,
        mean_ghsl_t1 = 0, mean_tier2_frac = 0,
        n_pass_samples = 0L, n_weak_samples = 0L, n_fail_samples = 0L,
        pass_fraction = 0, mean_phased_fraction = 0,
        snake3_status = "FAIL"
      )
    }
    all_summary_rows[[length(all_summary_rows) + 1]] <- data.table(
      chrom = chr, n_windows = n_win, n_pass = 0L, n_weak = 0L, n_fail = n_win,
      data_available = FALSE
    )
    next
  }

  phased <- fread(phased_file)
  message("[STEP10h-v2] Loaded ", nrow(phased), " phased variant records for ", chr)

  # ── ADAPTIVE THRESHOLD CALIBRATION ───────────────────────────────
  # Compute from actual data, not hardcoded assumptions
  block_info <- unique(phased[, .(block_id, block_start, block_end, sample_id)])
  block_info[, block_span := block_end - block_start]
  block_spans <- block_info$block_span[block_info$block_span > 0]

  if (length(block_spans) > 10) {
    median_block_span <- median(block_spans)
    q25_span <- quantile(block_spans, 0.25)
  } else {
    median_block_span <- 72  # fallback from CGA097 example
    q25_span <- 30
  }

  # Adaptive per-sample QC thresholds
  # For ~5x Illumina: typical blocks are 44-819 bp
  # A window might contain 20-50 blocks totaling 1500-5000 bp
  ADAPT_MIN_PHASED_BP       <- as.integer(max(100, median_block_span * 2))
  ADAPT_LOW_PHASED_BP       <- as.integer(max(30, q25_span))
  ADAPT_MIN_PHASED_FRACTION <- 0.01  # 1% of window span is realistic for 5x

  message("[STEP10h-v2] Adaptive thresholds for ", chr, ":")
  message("  median_block_span = ", round(median_block_span), " bp")
  message("  MIN_PHASED_BP     = ", ADAPT_MIN_PHASED_BP, " bp (was 5000)")
  message("  LOW_PHASED_BP     = ", ADAPT_LOW_PHASED_BP, " bp (was 2000)")
  message("  MIN_PHASED_FRAC   = ", ADAPT_MIN_PHASED_FRACTION)

  # ── Compute baseline GHSL distribution for calibration ───────────
  # We'll collect all per-sample GHSL values first, then set thresholds
  # relative to the chromosome-wide distribution

  # Pre-index per sample
  sample_phased <- split(phased, phased$sample_id)

  # ── Process each window ──────────────────────────────────────────
  chr_sw_rows <- list()

  for (wi in seq_len(n_win)) {
    wid     <- dt$global_window_id[wi]
    w_start <- dt$start_bp[wi]
    w_end   <- dt$end_bp[wi]

    sw_rows <- list()

    for (samp in sample_names) {
      sp <- sample_phased[[samp]]

      if (is.null(sp) || nrow(sp) == 0) {
        sw_rows[[length(sw_rows) + 1]] <- data.table(
          sample_id = samp, window_id = wid,
          total_window_bp = w_end - w_start,
          phased_bp = 0L, phased_fraction = 0,
          phased_variant_count = 0L, phased_het_snp_count = 0L,
          phased_het_indel_count = 0L, phased_het_indel_bp = 0L,
          block_count = 0L, n_tier1 = 0L, n_tier2 = 0L, tier2_fraction = 0,
          ghsl_snp_contrast = 0, ghsl_indel_contrast = 0,
          ghsl_indel_bp_contrast = 0, ghsl_total = 0,
          ghsl_total_weighted = 0, ghsl_t1_only = 0, fragmentation_index = 0,
          qc_flag = "NO_DATA"
        )
        next
      }

      # Filter to window
      sp_win <- sp[pos >= w_start & pos <= w_end]

      # Blocks overlapping window
      blocks_all <- unique(sp[, .(block_start, block_end)])
      blocks_win <- blocks_all[block_end >= w_start & block_start <= w_end]

      if (nrow(sp_win) == 0 && nrow(blocks_win) == 0) {
        sw_rows[[length(sw_rows) + 1]] <- data.table(
          sample_id = samp, window_id = wid,
          total_window_bp = w_end - w_start,
          phased_bp = 0L, phased_fraction = 0,
          phased_variant_count = 0L, phased_het_snp_count = 0L,
          phased_het_indel_count = 0L, phased_het_indel_bp = 0L,
          block_count = 0L, n_tier1 = 0L, n_tier2 = 0L, tier2_fraction = 0,
          ghsl_snp_contrast = 0, ghsl_indel_contrast = 0,
          ghsl_indel_bp_contrast = 0, ghsl_total = 0,
          ghsl_total_weighted = 0, ghsl_t1_only = 0, fragmentation_index = 0,
          qc_flag = "NO_DATA"
        )
        next
      }

      row <- compute_sample_contrast(samp, wid, w_start, w_end, sp_win, blocks_win)

      # QC flags (adaptive)
      flags <- character(0)
      if (row$phased_bp < ADAPT_LOW_PHASED_BP) flags <- c(flags, "LOW_PHASED_BP")
      if (row$phased_fraction < ADAPT_MIN_PHASED_FRACTION) flags <- c(flags, "LOW_PHASED_FRACTION")
      if (row$phased_variant_count < MIN_PHASED_VARIANT_COUNT) flags <- c(flags, "LOW_VARIANT_COUNT")
      if (row$block_count < MIN_BLOCK_COUNT) flags <- c(flags, "LOW_BLOCK_SUPPORT")
      if (row$tier2_fraction > 0.80) flags <- c(flags, "TIER2_DOMINATED")
      if (row$block_count >= 1 && row$block_count < 3) flags <- c(flags, "SPARSE_BLOCKS")

      # Overall QC
      if (row$phased_bp < ADAPT_MIN_PHASED_BP || row$phased_variant_count < MIN_PHASED_VARIANT_COUNT) {
        flags <- c(flags, "LOW_CONFIDENCE")
      }

      row[, qc_flag := if (length(flags) == 0) "PASS" else paste(flags, collapse = ";")]

      sw_rows[[length(sw_rows) + 1]] <- row
    }

    sw_dt <- rbindlist(sw_rows, fill = TRUE)
    sw_dt[, chrom := chr]
    chr_sw_rows[[length(chr_sw_rows) + 1]] <- sw_dt
  }

  chr_sw_all <- rbindlist(chr_sw_rows, fill = TRUE)
  all_sample_window[[length(all_sample_window) + 1]] <- chr_sw_all

  # ── DATA-DRIVEN GHSL THRESHOLDS ────────────────────────────────
  # Compute from passing samples across all windows
  pass_ghsl <- chr_sw_all[qc_flag == "PASS"]$ghsl_total_weighted

  if (length(pass_ghsl) > 20) {
    ghsl_median <- median(pass_ghsl, na.rm = TRUE)
    ghsl_iqr    <- IQR(pass_ghsl, na.rm = TRUE)
    GHSL_PASS_THRESH <- ghsl_median + 1.5 * ghsl_iqr
    GHSL_WEAK_THRESH <- ghsl_median + 0.5 * ghsl_iqr
    # Floor: don't set thresholds below minimum useful values
    GHSL_PASS_THRESH <- max(GHSL_PASS_THRESH, 0.5)
    GHSL_WEAK_THRESH <- max(GHSL_WEAK_THRESH, 0.1)
  } else {
    # Fallback if too few passing samples
    GHSL_PASS_THRESH <- 2.0
    GHSL_WEAK_THRESH <- 0.5
  }

  MIN_PASS_FRACTION <- 0.20  # relaxed from 0.30 for short-read data

  message("[STEP10h-v2] Data-driven thresholds:")
  message("  GHSL_PASS_THRESH  = ", round(GHSL_PASS_THRESH, 3))
  message("  GHSL_WEAK_THRESH  = ", round(GHSL_WEAK_THRESH, 3))
  message("  MIN_PASS_FRACTION = ", MIN_PASS_FRACTION)

  all_thresh_rows[[length(all_thresh_rows) + 1]] <- data.table(
    chrom = chr,
    median_block_span = round(median_block_span),
    adapt_min_phased_bp = ADAPT_MIN_PHASED_BP,
    adapt_low_phased_bp = ADAPT_LOW_PHASED_BP,
    ghsl_pass_thresh = round(GHSL_PASS_THRESH, 4),
    ghsl_weak_thresh = round(GHSL_WEAK_THRESH, 4),
    n_pass_samples_baseline = length(pass_ghsl),
    ghsl_median_baseline = round(median(pass_ghsl, na.rm = TRUE), 4),
    ghsl_iqr_baseline = round(IQR(pass_ghsl, na.rm = TRUE), 4)
  )

  # ── Window-level summary ─────────────────────────────────────────
  window_ids <- unique(chr_sw_all$window_id)

  for (wid in window_ids) {
    sw_win <- chr_sw_all[window_id == wid]

    n_pass <- sum(sw_win$qc_flag == "PASS", na.rm = TRUE)
    n_fail <- sum(grepl("LOW_CONFIDENCE|NO_DATA", sw_win$qc_flag), na.rm = TRUE)
    n_weak <- nrow(sw_win) - n_pass - n_fail
    pass_frac <- n_pass / nrow(sw_win)

    ghsl_vals <- sw_win[qc_flag == "PASS"]$ghsl_total_weighted
    if (length(ghsl_vals) == 0) ghsl_vals <- c(0)

    mean_g   <- mean(ghsl_vals, na.rm = TRUE)
    median_g <- median(ghsl_vals, na.rm = TRUE)
    max_g    <- max(ghsl_vals, na.rm = TRUE)
    q75 <- quantile(ghsl_vals, 0.75, na.rm = TRUE)
    q25 <- quantile(ghsl_vals, 0.25, na.rm = TRUE)
    spread_g <- as.numeric(q75 - q25)

    # TIER1-only GHSL
    ghsl_t1_vals <- sw_win[qc_flag == "PASS"]$ghsl_t1_only
    mean_t1 <- if (length(ghsl_t1_vals) > 0) mean(ghsl_t1_vals, na.rm = TRUE) else 0

    # Mean tier2 fraction
    mean_t2f <- mean(sw_win[qc_flag == "PASS"]$tier2_fraction, na.rm = TRUE)
    if (is.nan(mean_t2f)) mean_t2f <- 0

    # Window status
    if (pass_frac >= MIN_PASS_FRACTION && n_pass >= MIN_PASS_COUNT &&
        mean_g >= GHSL_PASS_THRESH) {
      w_status <- "PASS"
    } else if (pass_frac >= MIN_PASS_FRACTION * 0.5 && mean_g >= GHSL_WEAK_THRESH) {
      w_status <- "WEAK"
    } else {
      w_status <- "FAIL"
    }

    # Look up window coordinates from the master window grid
    wi_match <- which(dt$global_window_id == wid)
    w_start_bp <- if (length(wi_match) > 0) dt$start_bp[wi_match[1]] else NA_integer_
    w_end_bp   <- if (length(wi_match) > 0) dt$end_bp[wi_match[1]]   else NA_integer_

    all_track_rows[[length(all_track_rows) + 1]] <- data.table(
      chrom = chr,
      global_window_id = wid,
      start_bp = w_start_bp,
      end_bp = w_end_bp,
      mean_ghsl = round(mean_g, 4),
      median_ghsl = round(median_g, 4),
      max_ghsl = round(max_g, 4),
      spread_ghsl = round(spread_g, 4),
      mean_snp_contrast = round(mean(sw_win[qc_flag == "PASS"]$ghsl_snp_contrast, na.rm = TRUE), 4),
      mean_indel_contrast = round(mean(sw_win[qc_flag == "PASS"]$ghsl_indel_contrast, na.rm = TRUE), 4),
      mean_ghsl_t1 = round(mean_t1, 4),
      mean_tier2_frac = round(mean_t2f, 4),
      n_pass_samples = n_pass,
      n_weak_samples = n_weak,
      n_fail_samples = n_fail,
      pass_fraction = round(pass_frac, 4),
      mean_phased_fraction = round(mean(sw_win$phased_fraction, na.rm = TRUE), 4),
      snake3_status = w_status
    )
  }

  # Per-chromosome summary
  chr_statuses <- vapply(
    all_track_rows[vapply(all_track_rows, function(x) x$chrom[1] == chr, logical(1))],
    function(x) x$snake3_status, character(1)
  )
  all_summary_rows[[length(all_summary_rows) + 1]] <- data.table(
    chrom = chr, n_windows = n_win,
    n_pass = sum(chr_statuses == "PASS"),
    n_weak = sum(chr_statuses == "WEAK"),
    n_fail = sum(chr_statuses == "FAIL"),
    data_available = TRUE
  )

  message("[STEP10h-v2] ", chr, ": PASS=", sum(chr_statuses == "PASS"),
          " WEAK=", sum(chr_statuses == "WEAK"),
          " FAIL=", sum(chr_statuses == "FAIL"))
}

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

sw_all <- if (length(all_sample_window) > 0) rbindlist(all_sample_window, fill = TRUE) else {
  data.table(sample_id = character(), window_id = integer(), chrom = character())
}
track_dt <- if (length(all_track_rows) > 0) rbindlist(all_track_rows, fill = TRUE) else {
  data.table(chrom = character(), snake3_status = character())
}
summary_dt <- if (length(all_summary_rows) > 0) rbindlist(all_summary_rows) else {
  data.table(chrom = character())
}
thresh_dt <- if (length(all_thresh_rows) > 0) rbindlist(all_thresh_rows) else {
  data.table(chrom = character())
}

f1 <- file.path(outdir, "snake3_sample_window.tsv.gz")
f2 <- file.path(outdir, "snake3_track.tsv.gz")
f3 <- file.path(outdir, "snake3_summary.tsv")
f4 <- file.path(outdir, "snake3_thresholds.tsv")

fwrite(sw_all, f1, sep = "\t")
fwrite(track_dt, f2, sep = "\t")
fwrite(summary_dt, f3, sep = "\t")
fwrite(thresh_dt, f4, sep = "\t")

message("\n[DONE] STEP10h-v2 Snake 3 GHSL haplotype contrast (adaptive)")
message("  ", f1)
message("  ", f2)
message("  ", f3)
message("  ", f4)

# Quick summary
if (nrow(summary_dt) > 0) {
  total_pass <- sum(summary_dt$n_pass, na.rm = TRUE)
  total_weak <- sum(summary_dt$n_weak, na.rm = TRUE)
  total_fail <- sum(summary_dt$n_fail, na.rm = TRUE)
  message("\n  TOTAL across chromosomes: PASS=", total_pass,
          " WEAK=", total_weak, " FAIL=", total_fail)
}
if (nrow(thresh_dt) > 0) {
  message("\n  Threshold summary:")
  message("    Median block span range: ",
          min(thresh_dt$median_block_span), "-", max(thresh_dt$median_block_span), " bp")
  message("    GHSL PASS threshold range: ",
          min(thresh_dt$ghsl_pass_thresh), "-", max(thresh_dt$ghsl_pass_thresh))
}
