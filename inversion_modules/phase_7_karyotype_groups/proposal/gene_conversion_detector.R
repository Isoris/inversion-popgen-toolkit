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
# CALLED FROM: RUN_gene_conversion_detector.R (driver wired 2026-04-24 chat D;
# before that this detector was orphaned — built + tested but never called
# from production). The driver's per-candidate loop invokes:
#
#   detect_cohort_gene_conversion_v2(dosage_mat, snp_pos, snp_info,
#                                    baseline_class_by_sample,
#                                    candidate_id = cid, …)
#
# where dosage_mat is samples × SNPs (sliced from the per-chromosome
# DOSAGE_DIR files to the candidate interval), snp_pos is the parallel
# bp vector from <chr>.sites.tsv.gz, snp_info is an optional data.table
# with per-SNP QC columns (call_rate, het_fraction, maf, depth_z; not
# currently populated by the driver — relying on detector's fallback
# computation from dosage_mat alone), and baseline_class_by_sample is a
# named character vector with values in {HOM_REF, HET, HOM_INV, AMBIGUOUS}
# pulled from STEP_C01i_decompose's per_window_class.rds. Output is
# written both as JSON (for STEP_C01i_b_multi_recomb's load_gc_for_cid
# to consume via --gc_dir) and as a registry block (so the three flat
# keys q2_gc_total_tracts, q2_gc_total_samples_with_gc,
# q2_gc_n_snps_pass_qc_diagnostic surface via the schema's
# keys_extracted rule).
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
# =============================================================================
# For each QC-clean SNP, compute AF in HOM_REF and HOM_INV samples (using
# baseline_class_by_sample). Mark the SNP DIAGNOSTIC iff
# |AF_REF − AF_INV| ≥ min_delta_af. Returns a logical vector (same length
# as ncol(dosage_mat)); non-QC SNPs and SNPs without enough baseline-class
# samples return FALSE.
#
# Attributes on the returned mask (added 2026-04-24 chat D for the polarity
# fix in review §7):
#   attr(., "af_ref")  — per-SNP AF in HOM_REF samples, length ncol(dosage_mat).
#   attr(., "af_inv")  — per-SNP AF in HOM_INV samples, length ncol(dosage_mat).
# scan_one_sample consumes these to do the polarity-aware flagging (sample's
# dosage closer to the OTHER group's per-SNP mean than to its OWN). External
# callers that don't need the attributes just use the vector as before.
# =============================================================================
diagnostic_snp_mask <- function(dosage_mat, baseline_class_by_sample,
                                   qc_mask, min_delta_af = 0.5,
                                   min_per_class = 3L) {
  sids <- rownames(dosage_mat)
  if (is.null(sids)) stop("[gc_detector] dosage_mat needs rownames")
  nS <- ncol(dosage_mat)

  cls <- baseline_class_by_sample[sids]
  ref_idx <- which(cls == "HOM_REF")
  inv_idx <- which(cls == "HOM_INV")

  if (length(ref_idx) < min_per_class || length(inv_idx) < min_per_class) {
    out <- rep(FALSE, nS)
    attr(out, "af_ref") <- rep(NA_real_, nS)
    attr(out, "af_inv") <- rep(NA_real_, nS)
    return(out)
  }

  af_ref <- colMeans(dosage_mat[ref_idx, , drop = FALSE], na.rm = TRUE) / 2
  af_inv <- colMeans(dosage_mat[inv_idx, , drop = FALSE], na.rm = TRUE) / 2

  diag_delta <- abs(af_ref - af_inv)
  diag_delta[is.nan(diag_delta) | is.na(diag_delta)] <- 0
  mask <- qc_mask & (diag_delta >= min_delta_af)
  attr(mask, "af_ref") <- af_ref
  attr(mask, "af_inv") <- af_inv
  mask
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
                               max_flagged_snps = 10L,
                               af_ref        = NULL,
                               af_inv        = NULL) {

  # af_ref / af_inv: per-SNP allele-freq vectors (length = length(snp_pos)),
  # typically pulled off diag_mask via attr(diag_mask, "af_ref") /
  # attr(., "af_inv"). When both are provided, the scan uses POLARITY-AWARE
  # flagging: at each diagnostic SNP i, a sample is flagged iff its dosage
  # is closer to the OTHER arrangement's per-SNP group mean than to its OWN.
  # This is the fix for CODE_REVIEW_TO_FIGURE_OUT.md §7 — the previous
  # fixed-threshold rule (dosage >= 1.5 → flag for HOM_REF, dosage <= 0.5 →
  # flag for HOM_INV) assumed HOM_REF always sits at low dosage, which is
  # only true when ANGSD's major/minor call aligns with arrangement polarity.
  # In a balanced cohort a non-trivial fraction of diagnostic SNPs have
  # reverse polarity (HOM_REF at dosage ≈ 2) and the old rule produced false
  # positives at every such SNP.
  #
  # When af_ref / af_inv are NULL (no attrs, legacy caller), we fall back to
  # the old fixed-threshold rule with a one-line warning so the back-compat
  # path is auditable. Any new caller should pass them.

  empty <- data.table(
    sample_id = character(), tract_start_bp = integer(),
    tract_end_bp = integer(), tract_center_bp = integer(),
    span_bp = integer(), tract_width_bp = integer(),
    run_length_flagged_snps = integer(), tolerance_snps = integer(),
    confidence = character(),
    mean_dosage_tract = numeric(), mean_dosage_flank = numeric(),
    direction = character()
  )

  # Return shape: list(tracts = data.table, n_scanned = integer).
  empty_return <- list(tracts = empty, n_scanned = 0L)

  # HET samples: ambiguous per-haplotype at 9x — skip.
  # AMBIGUOUS / unknown: also skip.
  if (is.na(baseline_class) || !(baseline_class %in% c("HOM_REF", "HOM_INV"))) {
    return(empty_return)
  }

  # Restrict to diagnostic-QC-clean SNPs, in genomic order
  idx <- which(diag_mask)
  if (length(idx) < min_run) return(empty_return)
  ord <- order(snp_pos[idx])
  idx <- idx[ord]
  d   <- dosage_row[idx]
  p   <- snp_pos[idx]

  # Classify each position for this sample
  state <- character(length(idx))

  use_polarity_aware <- !is.null(af_ref) && !is.null(af_inv)
  if (use_polarity_aware) {
    # Per-SNP group means at the diagnostic-SNP subset (dosage scale, 0..2)
    mean_ref <- 2 * af_ref[idx]
    mean_inv <- 2 * af_inv[idx]
    # For each SNP, the sample's "own" mean is the one for its baseline class.
    own_mean   <- if (baseline_class == "HOM_REF") mean_ref else mean_inv
    other_mean <- if (baseline_class == "HOM_REF") mean_inv else mean_ref

    # Distance to each mean
    d_own   <- abs(d - own_mean)
    d_other <- abs(d - other_mean)

    # Required separation in distance terms: the two means differ by at
    # least min_delta_af * 2 (≥ 1.0 by default). A sample lies closer to
    # OTHER than to OWN when d_other < d_own, and we also require the
    # sample to be "meaningfully" on the other side — i.e. closer to OTHER
    # than halfway between the two means. This latter gate avoids flagging
    # samples with dosage exactly between the two (e.g. 1.0 at an SNP
    # where mean_ref=0.1, mean_inv=1.9 — an intermediate value is not a
    # confident "other-arrangement" signal).
    half <- abs(own_mean - other_mean) / 2
    flag_mask       <- !is.na(d) & d_other < d_own & d_own >= half
    consistent_mask <- !is.na(d) & d_own   < d_other & d_other >= half
    skip_mask       <- !is.na(d) & !flag_mask & !consistent_mask
    state[flag_mask]       <- "flag"
    state[consistent_mask] <- "consistent"
    state[skip_mask]       <- "skip"
    state[is.na(d)]        <- "skip"
  } else {
    # Legacy fixed-threshold fallback (fired when no af_ref/af_inv passed).
    # Emits a single warning per call to keep this auditable but non-noisy.
    warning("[gc_detector] scan_one_sample called without af_ref/af_inv — ",
            "falling back to fixed-polarity flagging. This is the pre-chat-D ",
            "behavior and is incorrect at reverse-polarity diagnostic SNPs ",
            "(review §7). Pass af_ref/af_inv from attr(diag_mask, ...).")
    if (baseline_class == "HOM_REF") {
      state[!is.na(d) & d <= 0.5]                <- "consistent"
      state[!is.na(d) & d >= 1.5]                <- "flag"
      state[!is.na(d) & d >  0.5 & d <  1.5]     <- "skip"
      state[is.na(d)]                            <- "skip"
    } else {
      state[!is.na(d) & d >= 1.5]                <- "consistent"
      state[!is.na(d) & d <= 0.5]                <- "flag"
      state[!is.na(d) & d >  0.5 & d <  1.5]     <- "skip"
      state[is.na(d)]                            <- "skip"
    }
  }

  # Count scanned diagnostic positions for this sample (non-skip).
  n_scanned <- as.integer(sum(state %in% c("consistent", "flag")))

  # Sanity: if no flags, no tracts.
  if (!any(state == "flag")) {
    return(list(tracts = empty, n_scanned = n_scanned))
  }

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

  tracts_dt <- if (length(tracts) == 0L) empty else rbindlist(tracts)
  list(tracts = tracts_dt, n_scanned = n_scanned)
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

  # Step 2: diagnostic mask (carries af_ref/af_inv as attributes since
  # chat D — see diagnostic_snp_mask header)
  diag_mask <- diagnostic_snp_mask(dosage_mat, baseline_class_by_sample,
                                   qc_mask, min_delta_af = min_delta_af,
                                   min_per_class = min_per_class)
  n_pass_qc_diag   <- sum(diag_mask, na.rm = TRUE)
  n_dropped_nondiag<- n_pass_qc - n_pass_qc_diag
  af_ref_v <- attr(diag_mask, "af_ref")
  af_inv_v <- attr(diag_mask, "af_inv")

  # Capture the positions behind the diagnostic mask, in genomic order.
  # Added 2026-04-24 (chat D) for auditability — schema now exposes a
  # diagnostic_snps block field so downstream consumers can see WHICH
  # positions drove the call, not just how many there were.
  diag_idx       <- which(diag_mask)
  diag_pos_ordered <- sort(as.integer(snp_pos[diag_idx]))

  # Steps 3–7: per-sample scan
  summary_rows <- vector("list", length(sids))
  tract_rows   <- vector("list", length(sids))
  per_sample_n_scanned <- integer(length(sids))
  names(per_sample_n_scanned) <- sids
  for (i in seq_along(sids)) {
    sid <- sids[i]
    cls <- baseline_class_by_sample[sid]
    if (is.null(cls) || is.na(cls)) cls <- NA_character_
    sr <- scan_one_sample(
      sid, dosage_mat[i, ], snp_pos, diag_mask, cls,
      min_run = min_run, max_tolerance = max_tolerance,
      max_span_bp = max_span_bp, max_flagged_snps = max_flagged_snps,
      af_ref = af_ref_v, af_inv = af_inv_v
    )
    # scan_one_sample returns list(tracts, n_scanned) as of 2026-04-24.
    tr <- sr$tracts
    per_sample_n_scanned[sid] <- sr$n_scanned
    tract_rows[[i]]   <- tr
    summary_rows[[i]] <- data.table(
      sample_id              = sid,
      n_tracts               = nrow(tr),
      baseline_class         = cls,
      n_diagnostic_scanned   = as.integer(sr$n_scanned),
      flank_dosage           = NA_real_,
      detector               = "per_snp_runlength_v1"
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

  # Per-sample scanned-count summary stats for the aggregates in keys.tsv.
  scanned_stats <- if (length(per_sample_n_scanned) > 0L) {
    list(
      n_samples_scannable      = as.integer(sum(per_sample_n_scanned > 0L)),
      mean_n_scanned_per_sample = round(mean(per_sample_n_scanned), 2),
      min_n_scanned_per_sample  = as.integer(min(per_sample_n_scanned)),
      max_n_scanned_per_sample  = as.integer(max(per_sample_n_scanned))
    )
  } else {
    list(n_samples_scannable = 0L, mean_n_scanned_per_sample = NA_real_,
         min_n_scanned_per_sample = NA_integer_,
         max_n_scanned_per_sample = NA_integer_)
  }

  diagnostic_snps_block <- list(
    # Scalar aggregates (also lifted into keys_extracted)
    n                         = as.integer(n_pass_qc_diag),
    n_samples_scannable       = scanned_stats$n_samples_scannable,
    mean_n_scanned_per_sample = scanned_stats$mean_n_scanned_per_sample,
    min_n_scanned_per_sample  = scanned_stats$min_n_scanned_per_sample,
    max_n_scanned_per_sample  = scanned_stats$max_n_scanned_per_sample,
    # Full position vector (integer bp, genomic order)
    positions_bp              = diag_pos_ordered
  )

  list(
    candidate_id           = candidate_id,
    per_sample_summary     = summary,
    tracts                 = tracts,
    total_tracts           = as.integer(nrow(tracts)),
    total_samples_with_gc  = as.integer(sum(summary$n_tracts > 0L)),
    snp_qc                 = snp_qc_block,
    diagnostic_snps        = diagnostic_snps_block,
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
# Legacy shim removed 2026-04-24 (chat D)
# =============================================================================
# detect_cohort_gene_conversion() (windowed-binning, no baseline_class)
# was a no-op warning stub kept for back-compat. A cohort-wide grep found
# zero live callers at the time of removal. The only off-file reference
# was a stale comment in lib_recomb_combination.R (fixed in the same
# pass). If you land here looking for the old function, the migration
# is: call detect_cohort_gene_conversion_v2() instead and pass a
# baseline_class_by_sample named character vector (values
# HOM_REF / HET / HOM_INV / AMBIGUOUS).
# =============================================================================
