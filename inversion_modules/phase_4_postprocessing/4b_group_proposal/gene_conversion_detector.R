#!/usr/bin/env Rscript
# =============================================================================
# gene_conversion_detector.R — chat-12 rewrite: per-SNP run-length detector
# =============================================================================
# WHAT THIS DETECTOR ACTUALLY MEASURES (honest framing, not the shorthand):
#
# This detector surfaces SHORT ARRANGEMENT-DISCORDANT TRACTS within a
# candidate inversion interval: stretches of ~2–10 consecutive diagnostic
# SNPs where a single sample's dosage pattern matches the OPPOSITE
# arrangement's baseline, then reverts to its own arrangement. The block
# and file are named "gene_conversion" for backward-compatibility with the
# rest of the pipeline, but that is shorthand.
#
# What we can say from the sequence data:
#   "Sample X, baseline HOM_REF, looks INV-like across SNPs at positions
#    p1..pk inside this interval."
#
# What we CANNOT say without additional experiments:
#   "Sample X underwent a gene-conversion event at this locus, with the
#    INV arrangement serving as the donor template."
#
# The classical interpretation is gene conversion (Harris 1993: tracts
# 18–774 bp, 3–26 bp boundary heteroduplex regions). But the same signal
# is also produced by:
#   - Short double crossovers
#   - Paralog mismapping (reads from a paralogous locus scoring at the
#     reference locus)
#   - Simply coincidental IBS with the other arrangement at a handful of
#     diagnostic sites — we don't have the pedigree to prove otherwise
#
# For Quentin's cohort (no trios, no mutation-accumulation panel, no
# long-read haplotypes at these tracts), the detector's output should be
# read as "arrangement-discordant IBS tracts consistent with gene
# conversion", not "confirmed gene conversion events". The downstream
# terminology in figures and manuscript text should reflect this — say
# "GC-like tracts" or "arrangement-discordant IBS tracts" rather than
# "gene conversion events" unless there's supporting independent evidence.
#
# This does NOT change the math or the gate: the detector's job is to
# surface the tract. Interpretation of cause is downstream.
#
# =============================================================================
# PRIOR STATE (deleted): A windowed-binning detector (40-SNP windows, 10-SNP
# step, CUSUM-style excursion sweep on window-mean dosage, max_tract_bins
# gate). Its lower detection limit was ~60 kb (40 SNPs × 1.6 kb/SNP at
# Quentin's density), i.e. LARGER than the biology of GC tracts (typically
# 50 bp – ~800 bp; Harris 1993 upper bound ~5 kb is atypical). The real
# tracts were invisible to that detector, and the signals it did surface
# were often SNP-column artefacts mislabelled as tracts.
#
# THIS FILE (chat-12 rewrite): per-SNP run-length detector, specified by the
# handoff under "Revised gene-conversion detector spec". High-level flow:
#
#   Step 1  SNP QC pre-filter (missingness, excess-het/paralog, depth
#            anomaly, MAF).
#   Step 2  Diagnostic-SNP mask (|AF(HOM_REF) − AF(HOM_INV)| ≥ min_delta_af).
#   Step 3  Per-sample per-SNP flag: for a HOM_REF sample, a diagnostic SNP
#            is FLAGGED when dosage ≥ 1.5 (it looks INV at a site it should
#            look REF); symmetrically for HOM_INV. HET samples are skipped
#            (ambiguous per-haplotype at 9x).
#   Step 4  Run-length aggregation in genomic order with tolerance: a tract
#            is a maximal run of flagged SNPs admitting up to max_tolerance
#            intervening non-flagged-non-consistent SNPs (HET-looking or
#            missing). Tract ends on max_tolerance+1 consecutive non-flagged
#            SNPs.
#   Step 5  Length gate: keep tracts with run_length ∈ [min_run,
#            max_flagged_snps] and span_bp ≤ max_span_bp. Longer events are
#            recombinants, handled by C01j — discard here.
#   Step 6  Confidence from run length under a geometric prior p=0.01
#            (9x per-SNP genotyping error): len=2 LOW, 3 MEDIUM, ≥4 HIGH.
#   Step 7  Direction annotation (REF_in_INV_context etc) from the sample's
#            baseline class and the tract dosage pattern.
#   Step 8  Cohort summary and snp_qc accounting block for the schema.
#
# OUTPUT (matches gene_conversion_tracts.schema.json, chat-12 update):
#   list(
#     candidate_id           = cid,
#     per_sample_summary     = data.table(sample_id, n_tracts, baseline_class,
#                                         flank_dosage = NA, detector),
#     tracts                 = data.table of per-tract fields,
#     total_tracts           = int,
#     total_samples_with_gc  = int,
#     snp_qc                 = list(…filter counts…),
#     params                 = list(…thresholds used…)
#   )
#
# CALLED FROM: STEP_C01i_b_multi_recomb.R (or a dedicated driver) with:
#   detect_cohort_gene_conversion_v2(dosage_mat, snp_pos, snp_info,
#                                    baseline_class_by_sample,
#                                    candidate_id = cid, …)
# where dosage_mat is samples × SNPs, snp_pos is a bp vector of length
# ncol(dosage_mat), snp_info is an optional data.table with per-SNP QC
# columns (call_rate, het_fraction, maf, depth_z), and
# baseline_class_by_sample is a named character vector with values in
# {HOM_REF, HET, HOM_INV, AMBIGUOUS}.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# Step 1: SNP QC pre-filter
# =============================================================================
# Inputs:
#   dosage_mat : numeric matrix, samples × SNPs (rownames=sample_ids, optional
#                 colnames=snp_ids). NA allowed for missing.
#   snp_info   : optional data.table with OPTIONAL precomputed per-SNP QC
#                 columns (call_rate, het_fraction, maf, depth_z). If absent,
#                 the detector computes them from dosage_mat on the fly.
# Returns a logical vector of length ncol(dosage_mat) — TRUE = passes QC.
# Also returns the per-reason drop counts via an attribute "drop_counts".
# =============================================================================
snp_qc_mask <- function(dosage_mat, snp_info = NULL,
                           max_het_fraction = 0.70,
                           min_call_rate    = 0.80,
                           min_maf          = 0.02,
                           depth_z_abs      = 3.0) {
  n_samples <- nrow(dosage_mat)
  n_snps    <- ncol(dosage_mat)

  call_rate <- if (!is.null(snp_info) && "call_rate" %in% names(snp_info)) {
    as.numeric(snp_info$call_rate)
  } else {
    colSums(!is.na(dosage_mat)) / n_samples
  }

  het_fraction <- if (!is.null(snp_info) && "het_fraction" %in% names(snp_info)) {
    as.numeric(snp_info$het_fraction)
  } else {
    # Count "het-looking" samples per SNP: 0.5 < dosage < 1.5
    vapply(seq_len(n_snps), function(j) {
      d <- dosage_mat[, j]
      d <- d[!is.na(d)]
      if (length(d) == 0) return(NA_real_)
      sum(d > 0.5 & d < 1.5) / length(d)
    }, numeric(1))
  }

  maf <- if (!is.null(snp_info) && "maf" %in% names(snp_info)) {
    as.numeric(snp_info$maf)
  } else {
    # MAF from mean dosage/2 under a biallelic-ish assumption
    vapply(seq_len(n_snps), function(j) {
      d <- dosage_mat[, j]
      d <- d[!is.na(d)]
      if (length(d) == 0) return(NA_real_)
      p <- mean(d) / 2
      min(p, 1 - p)
    }, numeric(1))
  }

  # Depth z only if supplied (we can't synthesize depth from dosage alone)
  depth_z <- if (!is.null(snp_info) && "depth_z" %in% names(snp_info)) {
    as.numeric(snp_info$depth_z)
  } else {
    rep(0, n_snps)
  }

  drop_missing    <- !is.na(call_rate)    & call_rate    < min_call_rate
  drop_excess_het <- !is.na(het_fraction) & het_fraction > max_het_fraction
  drop_low_maf    <- !is.na(maf)          & maf          < min_maf
  drop_depth      <- !is.na(depth_z)      & abs(depth_z) > depth_z_abs

  pass <- !(drop_missing | drop_excess_het | drop_low_maf | drop_depth)

  attr(pass, "drop_counts") <- list(
    n_dropped_missingness   = as.integer(sum(drop_missing,    na.rm = TRUE)),
    n_dropped_excess_het    = as.integer(sum(drop_excess_het, na.rm = TRUE)),
    n_dropped_low_maf       = as.integer(sum(drop_low_maf,    na.rm = TRUE)),
    n_dropped_depth_anomaly = as.integer(sum(drop_depth,      na.rm = TRUE))
  )
  pass
}

# =============================================================================
# Step 2: diagnostic-SNP mask
# =============================================================================
# For each QC-clean SNP, compute AF in HOM_REF and HOM_INV samples (using
# baseline_class_by_sample). Mark the SNP DIAGNOSTIC iff
# |AF_REF − AF_INV| ≥ min_delta_af. Returns a logical vector (same length
# as ncol(dosage_mat)); non-QC SNPs and SNPs without enough baseline-class
# samples return FALSE.
# =============================================================================
diagnostic_snp_mask <- function(dosage_mat, baseline_class_by_sample,
                                   qc_mask, min_delta_af = 0.5,
                                   min_per_class = 3L) {
  sids <- rownames(dosage_mat)
  if (is.null(sids)) stop("[gc_detector] dosage_mat needs rownames")

  cls <- baseline_class_by_sample[sids]
  ref_idx <- which(cls == "HOM_REF")
  inv_idx <- which(cls == "HOM_INV")

  if (length(ref_idx) < min_per_class || length(inv_idx) < min_per_class) {
    return(rep(FALSE, ncol(dosage_mat)))
  }

  af_ref <- colMeans(dosage_mat[ref_idx, , drop = FALSE], na.rm = TRUE) / 2
  af_inv <- colMeans(dosage_mat[inv_idx, , drop = FALSE], na.rm = TRUE) / 2

  diag_delta <- abs(af_ref - af_inv)
  diag_delta[is.nan(diag_delta) | is.na(diag_delta)] <- 0
  qc_mask & (diag_delta >= min_delta_af)
}

# =============================================================================
# Step 3 + 4: per-sample per-SNP flagging and run-length aggregation
# =============================================================================
# For one sample at the diagnostic-QC-clean subset of SNPs, build a per-SNP
# state vector with three values:
#   "flag" (flagged = inconsistent with baseline arrangement)
#   "consistent" (matches baseline arrangement)
#   "skip" (HET-looking intermediate dosage OR missing — counted toward the
#           in-tract tolerance budget but cannot start/end a tract)
# Then walk the state vector in genomic order finding maximal runs that
# start + end with "flag" and contain ≤ max_tolerance "skip" interior
# positions. Runs are broken by max_tolerance+1 consecutive "skip" (or by
# any "consistent").
#
# Returns a data.table (possibly 0 rows) with columns:
#   sample_id, tract_start_bp, tract_end_bp, tract_center_bp,
#   span_bp, tract_width_bp (=span_bp alias),
#   run_length_flagged_snps, tolerance_snps, confidence,
#   mean_dosage_tract, mean_dosage_flank, direction
# =============================================================================
scan_one_sample <- function(sid, dosage_row, snp_pos, diag_mask,
                               baseline_class,
                               min_run       = 2L,
                               max_tolerance = 1L,
                               max_span_bp   = 20000L,
                               max_flagged_snps = 10L) {

  empty <- data.table(
    sample_id = character(), tract_start_bp = integer(),
    tract_end_bp = integer(), tract_center_bp = integer(),
    span_bp = integer(), tract_width_bp = integer(),
    run_length_flagged_snps = integer(), tolerance_snps = integer(),
    confidence = character(),
    mean_dosage_tract = numeric(), mean_dosage_flank = numeric(),
    direction = character()
  )

  # HET samples: ambiguous per-haplotype at 9x — skip.
  # AMBIGUOUS / unknown: also skip.
  if (is.na(baseline_class) || !(baseline_class %in% c("HOM_REF", "HOM_INV"))) {
    return(empty)
  }

  # Restrict to diagnostic-QC-clean SNPs, in genomic order
  idx <- which(diag_mask)
  if (length(idx) < min_run) return(empty)
  ord <- order(snp_pos[idx])
  idx <- idx[ord]
  d   <- dosage_row[idx]
  p   <- snp_pos[idx]

  # Classify each position for this sample
  state <- character(length(idx))
  if (baseline_class == "HOM_REF") {
    state[!is.na(d) & d <= 0.5]                <- "consistent"   # REF-looking
    state[!is.na(d) & d >= 1.5]                <- "flag"         # INV-looking
    state[!is.na(d) & d >  0.5 & d <  1.5]     <- "skip"         # HET-looking
    state[is.na(d)]                            <- "skip"         # missing
  } else {  # HOM_INV
    state[!is.na(d) & d >= 1.5]                <- "consistent"
    state[!is.na(d) & d <= 0.5]                <- "flag"
    state[!is.na(d) & d >  0.5 & d <  1.5]     <- "skip"
    state[is.na(d)]                            <- "skip"
  }

  # Sanity: if no flags, no tracts.
  if (!any(state == "flag")) return(empty)

  # Walk: find maximal tracts that
  #   - start on a "flag"
  #   - end on a "flag"
  #   - have interior "skip" count ≤ max_tolerance
  #   - are terminated by any "consistent" OR by max_tolerance+1 consecutive
  #     "skip"
  n <- length(state)
  tracts <- list()
  i <- 1L
  while (i <= n) {
    if (state[i] != "flag") { i <- i + 1L; next }
    # Start a tract at i
    start_i <- i
    last_flag_i <- i
    tol_used <- 0L
    consec_skip <- 0L
    j <- i + 1L
    while (j <= n) {
      if (state[j] == "flag") {
        last_flag_i <- j
        consec_skip <- 0L
      } else if (state[j] == "skip") {
        consec_skip <- consec_skip + 1L
        tol_used    <- tol_used + 1L
        if (consec_skip > max_tolerance) break
      } else {  # "consistent" → hard break
        break
      }
      j <- j + 1L
    }
    # Tract spans start_i .. last_flag_i (trim any trailing skips)
    run_length <- sum(state[start_i:last_flag_i] == "flag")
    interior_skips <- sum(state[start_i:last_flag_i] == "skip")
    if (run_length >= min_run && run_length <= max_flagged_snps &&
        interior_skips <= max_tolerance) {
      span_bp <- p[last_flag_i] - p[start_i]
      if (span_bp <= max_span_bp) {
        # Confidence from run length under geometric prior p=0.01 at 9x
        conf <- if (run_length >= 4L) "HIGH"
                else if (run_length == 3L) "MEDIUM"
                else "LOW"

        tract_idx  <- start_i:last_flag_i
        flag_idx   <- tract_idx[state[tract_idx] == "flag"]
        # Flank = consistent diagnostic SNPs OUTSIDE the tract but inside
        # the provided positions
        flank_mask <- state == "consistent"
        flank_mask[tract_idx] <- FALSE
        mean_dos_tract <- mean(d[flag_idx],     na.rm = TRUE)
        mean_dos_flank <- if (any(flank_mask)) mean(d[flank_mask], na.rm = TRUE)
                          else NA_real_

        direction <- if (baseline_class == "HOM_REF") "INV_in_REF_context"
                     else                            "REF_in_INV_context"

        tracts[[length(tracts) + 1L]] <- data.table(
          sample_id               = sid,
          tract_start_bp          = as.integer(p[start_i]),
          tract_end_bp            = as.integer(p[last_flag_i]),
          tract_center_bp         = as.integer((p[start_i] + p[last_flag_i]) %/% 2),
          span_bp                 = as.integer(span_bp),
          tract_width_bp          = as.integer(span_bp),
          run_length_flagged_snps = as.integer(run_length),
          tolerance_snps          = as.integer(interior_skips),
          confidence              = conf,
          mean_dosage_tract       = round(mean_dos_tract, 3),
          mean_dosage_flank       = if (is.na(mean_dos_flank)) NA_real_
                                    else round(mean_dos_flank, 3),
          direction               = direction
        )
      }
    }
    # Advance past the last flag in the tract
    i <- last_flag_i + 1L
  }

  if (length(tracts) == 0L) empty else rbindlist(tracts)
}

# =============================================================================
# Cohort driver
# =============================================================================
# Primary public entry point. Returns the full block in the shape expected
# by gene_conversion_tracts.schema.json (chat-12 version).
# =============================================================================
detect_cohort_gene_conversion_v2 <- function(dosage_mat, snp_pos,
                                                baseline_class_by_sample,
                                                candidate_id,
                                                snp_info        = NULL,
                                                min_run         = 2L,
                                                max_tolerance   = 1L,
                                                max_span_bp     = 20000L,
                                                max_flagged_snps= 10L,
                                                min_delta_af    = 0.5,
                                                max_het_fraction= 0.70,
                                                min_maf         = 0.02,
                                                min_call_rate   = 0.80,
                                                depth_z_abs     = 3.0,
                                                min_per_class   = 3L) {

  stopifnot(is.matrix(dosage_mat),
            ncol(dosage_mat) == length(snp_pos),
            !is.null(rownames(dosage_mat)))

  sids <- rownames(dosage_mat)
  n_snps_total <- ncol(dosage_mat)

  # Step 1: SNP QC
  qc_mask <- snp_qc_mask(dosage_mat, snp_info,
                         max_het_fraction = max_het_fraction,
                         min_call_rate    = min_call_rate,
                         min_maf          = min_maf,
                         depth_z_abs      = depth_z_abs)
  qc_drops <- attr(qc_mask, "drop_counts")
  n_pass_qc <- sum(qc_mask, na.rm = TRUE)

  # Step 2: diagnostic mask
  diag_mask <- diagnostic_snp_mask(dosage_mat, baseline_class_by_sample,
                                   qc_mask, min_delta_af = min_delta_af,
                                   min_per_class = min_per_class)
  n_pass_qc_diag   <- sum(diag_mask, na.rm = TRUE)
  n_dropped_nondiag<- n_pass_qc - n_pass_qc_diag

  # Steps 3–7: per-sample scan
  summary_rows <- vector("list", length(sids))
  tract_rows   <- vector("list", length(sids))
  for (i in seq_along(sids)) {
    sid <- sids[i]
    cls <- baseline_class_by_sample[sid]
    if (is.null(cls) || is.na(cls)) cls <- NA_character_
    tr <- scan_one_sample(
      sid, dosage_mat[i, ], snp_pos, diag_mask, cls,
      min_run = min_run, max_tolerance = max_tolerance,
      max_span_bp = max_span_bp, max_flagged_snps = max_flagged_snps
    )
    tract_rows[[i]]   <- tr
    summary_rows[[i]] <- data.table(
      sample_id      = sid,
      n_tracts       = nrow(tr),
      baseline_class = cls,
      flank_dosage   = NA_real_,
      detector       = "per_snp_runlength_v1"
    )
  }

  tracts  <- rbindlist(tract_rows, fill = TRUE)
  summary <- rbindlist(summary_rows)

  snp_qc_block <- list(
    n_snps_total               = as.integer(n_snps_total),
    n_snps_pass_qc             = as.integer(n_pass_qc),
    n_snps_pass_qc_diagnostic  = as.integer(n_pass_qc_diag),
    n_dropped_missingness      = qc_drops$n_dropped_missingness,
    n_dropped_excess_het       = qc_drops$n_dropped_excess_het,
    n_dropped_depth_anomaly    = qc_drops$n_dropped_depth_anomaly,
    n_dropped_low_maf          = qc_drops$n_dropped_low_maf,
    n_dropped_non_diagnostic   = as.integer(n_dropped_nondiag)
  )

  list(
    candidate_id           = candidate_id,
    per_sample_summary     = summary,
    tracts                 = tracts,
    total_tracts           = as.integer(nrow(tracts)),
    total_samples_with_gc  = as.integer(sum(summary$n_tracts > 0L)),
    snp_qc                 = snp_qc_block,
    params                 = list(
      min_run          = as.integer(min_run),
      max_tolerance    = as.integer(max_tolerance),
      max_span_bp      = as.integer(max_span_bp),
      max_flagged_snps = as.integer(max_flagged_snps),
      min_delta_af     = as.numeric(min_delta_af),
      max_het_fraction = as.numeric(max_het_fraction),
      min_maf          = as.numeric(min_maf),
      min_call_rate    = as.numeric(min_call_rate),
      depth_z_abs      = as.numeric(depth_z_abs)
    )
  )
}

# =============================================================================
# Back-compat shim
# =============================================================================
# Old callers used detect_cohort_gene_conversion(dosage_mat, snp_pos, ...).
# Without baseline_class_by_sample we cannot run the per-SNP detector, so
# the shim emits a warning and returns an empty block with the legacy-shape
# summary columns populated. Call sites that want real output must migrate
# to detect_cohort_gene_conversion_v2().
# =============================================================================
detect_cohort_gene_conversion <- function(dosage_mat, snp_pos,
                                              window_snps = 40L,
                                              step_snps = 10L,
                                              max_tract_bins = 5L,
                                              min_magnitude = 0.5) {
  warning("[gc_detector] detect_cohort_gene_conversion() is the LEGACY shim; ",
          "the per-SNP detector requires baseline_class_by_sample. ",
          "Migrate to detect_cohort_gene_conversion_v2().")
  sids <- rownames(dosage_mat) %||% character(0)
  list(
    summary = data.table(
      sample_id    = sids,
      n_tracts     = rep(0L, length(sids)),
      flank_dosage = rep(NA_real_, length(sids)),
      detector     = rep("legacy_shim_noop", length(sids))
    ),
    tracts_by_sample = setNames(
      lapply(sids, function(x) data.table()),
      sids
    )
  )
}
