#!/usr/bin/env Rscript
# =============================================================================
# _qc_shelf_reader.R — q_qc_shelf_* columns for compute_candidate_status.R
# =============================================================================
#
# Pass 13 (2026-04-24). Parallel to _axis5_final_label.R — same env-gated
# source pattern, additive columns, zero risk to existing axis logic.
#
# Reads the per-candidate summary.json files that Q10 of 4b_qc_triage writes
# (method = "phase_qc_shelf" in the evidence registry), derives a
# categorical flag from the shelf-vs-reference ratios, and appends the
# following columns to status_dt:
#
#   q_qc_shelf_ran                    logical
#   q_qc_shelf_flag                   character: clean | low_snp | high_uncertain |
#                                                coverage_artifact | messy | unknown
#   q_qc_shelf_snp_ratio              numeric
#   q_qc_shelf_uncertain_ratio        numeric
#   q_qc_shelf_coverage_ratio         numeric
#   q_qc_shelf_coverage_cv_ratio      numeric
#   q_qc_shelf_z_flatness             numeric
#   q_qc_shelf_fst_hom1_hom2_inside   numeric
#   q_qc_shelf_fst_hom1_hom2_outside  numeric
#   q_qc_shelf_fst_enrichment_fold    numeric
#   q_qc_shelf_hovere_het_inside      numeric
#   q_qc_shelf_hovere_hom1_inside     numeric
#   q_qc_shelf_hovere_hom2_inside     numeric
#
# USAGE (from compute_candidate_status.R, mirroring the axis5 pattern):
#
#   qc_shelf_dir <- Sys.getenv("QC_SHELF_EVIDENCE_DIR", "")
#   if (nzchar(qc_shelf_dir) && dir.exists(qc_shelf_dir)) {
#     source("_qc_shelf_reader.R")
#     status_dt <- append_qc_shelf_to_rows(status_dt, qc_shelf_dir)
#   }
#
# The reader finds summary JSON files by walking qc_shelf_dir for files
# matching `summary_*.json`. It expects each JSON to have:
#   $cid, $chrom
#   $summary_stats$fst_hom1_hom2_mean_inside / _outside
#   $summary_stats$fst_enrichment_fold
#   $summary_stats$hovere_HOM1_inside / hovere_HET_inside / hovere_HOM2_inside
#
# Candidates with no matching summary.json get q_qc_shelf_ran = FALSE and
# NA for all other fields. They are NOT dropped — this is a soft flag.
#
# FLAG DECISION RULES (thresholds from 4b_qc_triage/docs/interpretation.md):
#   low_snp            : snp_ratio <= 0.50
#   high_uncertain     : uncertain_ratio >= 2.00
#   coverage_artifact  : coverage_ratio <= 0.60 AND coverage_cv_ratio <= 0.70
#                        (i.e. low coverage WITH uniformity → mapping artifact)
#   messy              : any TWO of low_snp / high_uncertain / coverage_artifact
#   clean              : otherwise, with an additional requirement that
#                        fst_enrichment_fold >= 2 (the "real inversion" signal)
#                        OR hovere_HET_inside >= 1.3 (heterokaryotype signature)
#   unknown            : thresholds can't be computed (missing ratios)
#
# Thresholds are configurable via env vars so they can be tuned without
# editing this file:
#   QC_SHELF_SNP_RATIO_THRESH   (default 0.50)
#   QC_SHELF_UNC_RATIO_THRESH   (default 2.00)
#   QC_SHELF_COV_RATIO_THRESH   (default 0.60)
#   QC_SHELF_COV_CV_THRESH      (default 0.70)
#   QC_SHELF_FST_FOLD_THRESH    (default 2.00)
#   QC_SHELF_HOVERE_HET_THRESH  (default 1.30)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# -----------------------------------------------------------------------------
# Threshold config — env-overridable
# -----------------------------------------------------------------------------
.qc_shelf_threshold <- function(key, default) {
  v <- Sys.getenv(key, NA_character_)
  if (is.na(v) || !nzchar(v)) return(default)
  out <- suppressWarnings(as.numeric(v))
  if (is.finite(out)) out else default
}

# -----------------------------------------------------------------------------
# Read one summary.json file safely
# -----------------------------------------------------------------------------
.read_one_qc_shelf_summary <- function(path) {
  if (!file.exists(path)) return(NULL)
  js <- tryCatch(fromJSON(path, simplifyVector = FALSE),
                 error = function(e) NULL)
  if (is.null(js)) return(NULL)

  ss <- js$summary_stats
  if (is.null(ss)) ss <- list()

  # Access helper: pull value or NA
  .get <- function(lst, key) {
    v <- lst[[key]]
    if (is.null(v) || length(v) == 0) return(NA_real_)
    as.numeric(v[[1]])
  }

  list(
    cid                             = if (is.null(js$cid)) NA_character_ else as.character(js$cid[[1]]),
    chrom                           = if (is.null(js$chrom)) NA_character_ else as.character(js$chrom[[1]]),
    fst_hom1_hom2_inside            = .get(ss, "fst_hom1_hom2_mean_inside"),
    fst_hom1_hom2_outside           = .get(ss, "fst_hom1_hom2_mean_outside"),
    fst_enrichment_fold             = .get(ss, "fst_enrichment_fold"),
    hovere_hom1_inside              = .get(ss, "hovere_HOM1_inside"),
    hovere_het_inside               = .get(ss, "hovere_HET_inside"),
    hovere_hom2_inside              = .get(ss, "hovere_HOM2_inside"),
    # Ratios from Q04 summary: written if the q04 console-ratios companion
    # file exists (future extension). For now, default to NA — the reader
    # prioritises Fst + Hovere signals from summary.json.
    snp_ratio                       = .get(ss, "snp_ratio"),
    uncertain_ratio                 = .get(ss, "uncertain_ratio"),
    coverage_ratio                  = .get(ss, "coverage_ratio"),
    coverage_cv_ratio               = .get(ss, "coverage_cv_ratio"),
    z_flatness                      = .get(ss, "z_flatness")
  )
}

# -----------------------------------------------------------------------------
# Classify a single row into a flag
# -----------------------------------------------------------------------------
.classify_qc_shelf_flag <- function(row) {
  snp_t       <- .qc_shelf_threshold("QC_SHELF_SNP_RATIO_THRESH",  0.50)
  unc_t       <- .qc_shelf_threshold("QC_SHELF_UNC_RATIO_THRESH",  2.00)
  cov_t       <- .qc_shelf_threshold("QC_SHELF_COV_RATIO_THRESH",  0.60)
  cov_cv_t    <- .qc_shelf_threshold("QC_SHELF_COV_CV_THRESH",     0.70)
  fst_fold_t  <- .qc_shelf_threshold("QC_SHELF_FST_FOLD_THRESH",   2.00)
  hovere_t    <- .qc_shelf_threshold("QC_SHELF_HOVERE_HET_THRESH", 1.30)

  sr  <- row$snp_ratio
  ur  <- row$uncertain_ratio
  cr  <- row$coverage_ratio
  cvr <- row$coverage_cv_ratio
  ff  <- row$fst_enrichment_fold
  hh  <- row$hovere_het_inside

  # If ALL artifact-sensing ratios are NA and we don't have Fst/Hovere either,
  # we can't classify.
  if (!any(is.finite(c(sr, ur, cr, cvr, ff, hh)))) return("unknown")

  low_snp     <- is.finite(sr) && sr <= snp_t
  high_unc    <- is.finite(ur) && ur >= unc_t
  cov_artif   <- is.finite(cr) && is.finite(cvr) && cr <= cov_t && cvr <= cov_cv_t

  n_artifacts <- sum(low_snp, high_unc, cov_artif)

  if (n_artifacts >= 2)       return("messy")
  if (cov_artif)              return("coverage_artifact")
  if (high_unc)               return("high_uncertain")
  if (low_snp)                return("low_snp")

  # Otherwise: check for positive biological signal
  real_fst    <- is.finite(ff) && ff >= fst_fold_t
  real_hovere <- is.finite(hh) && hh >= hovere_t
  if (real_fst || real_hovere) return("clean")

  # No artifacts AND no positive signal → hard to call; be neutral
  "unknown"
}

# -----------------------------------------------------------------------------
# Main entry point — append q_qc_shelf_* columns to status_dt in place
# -----------------------------------------------------------------------------
append_qc_shelf_to_rows <- function(status_dt, qc_shelf_dir) {
  if (!is.data.table(status_dt)) status_dt <- as.data.table(status_dt)

  # Candidates in status_dt may use 'candidate_id' or 'cid'
  cid_col <- intersect(c("candidate_id", "cid", "interval_id"), names(status_dt))[1]
  if (is.na(cid_col)) {
    warning("[qc_shelf] no candidate_id / cid / interval_id column found on status_dt; skipping.")
    return(status_dt)
  }

  # Find all summary_*.json files under qc_shelf_dir
  summary_files <- list.files(qc_shelf_dir,
                              pattern = "^summary_.*\\.json$",
                              full.names = TRUE,
                              recursive  = TRUE)

  # Read & index by cid
  rows_by_cid <- list()
  for (f in summary_files) {
    r <- .read_one_qc_shelf_summary(f)
    if (is.null(r) || is.na(r$cid)) next
    rows_by_cid[[r$cid]] <- r
  }

  n_found <- length(rows_by_cid)
  cat(sprintf(
    "[qc_shelf] read %d summary.json files from %s\n",
    n_found, qc_shelf_dir
  ))

  # Build the append table row-by-row
  n <- nrow(status_dt)
  out <- data.table(
    q_qc_shelf_ran                   = logical(n),
    q_qc_shelf_flag                  = character(n),
    q_qc_shelf_snp_ratio             = numeric(n),
    q_qc_shelf_uncertain_ratio       = numeric(n),
    q_qc_shelf_coverage_ratio        = numeric(n),
    q_qc_shelf_coverage_cv_ratio     = numeric(n),
    q_qc_shelf_z_flatness            = numeric(n),
    q_qc_shelf_fst_hom1_hom2_inside  = numeric(n),
    q_qc_shelf_fst_hom1_hom2_outside = numeric(n),
    q_qc_shelf_fst_enrichment_fold   = numeric(n),
    q_qc_shelf_hovere_het_inside     = numeric(n),
    q_qc_shelf_hovere_hom1_inside    = numeric(n),
    q_qc_shelf_hovere_hom2_inside    = numeric(n)
  )

  # Initialise all numerics to NA (they default to 0)
  num_cols <- setdiff(names(out), c("q_qc_shelf_ran", "q_qc_shelf_flag"))
  for (nc in num_cols) out[[nc]] <- NA_real_
  out$q_qc_shelf_flag <- NA_character_
  out$q_qc_shelf_ran  <- FALSE

  n_matched <- 0L
  for (i in seq_len(n)) {
    cid_i <- status_dt[[cid_col]][i]
    r <- rows_by_cid[[as.character(cid_i)]]
    if (is.null(r)) next
    n_matched <- n_matched + 1L

    out$q_qc_shelf_ran[i]                   <- TRUE
    out$q_qc_shelf_flag[i]                  <- .classify_qc_shelf_flag(r)
    out$q_qc_shelf_snp_ratio[i]             <- r$snp_ratio
    out$q_qc_shelf_uncertain_ratio[i]       <- r$uncertain_ratio
    out$q_qc_shelf_coverage_ratio[i]        <- r$coverage_ratio
    out$q_qc_shelf_coverage_cv_ratio[i]     <- r$coverage_cv_ratio
    out$q_qc_shelf_z_flatness[i]            <- r$z_flatness
    out$q_qc_shelf_fst_hom1_hom2_inside[i]  <- r$fst_hom1_hom2_inside
    out$q_qc_shelf_fst_hom1_hom2_outside[i] <- r$fst_hom1_hom2_outside
    out$q_qc_shelf_fst_enrichment_fold[i]   <- r$fst_enrichment_fold
    out$q_qc_shelf_hovere_het_inside[i]     <- r$hovere_het_inside
    out$q_qc_shelf_hovere_hom1_inside[i]    <- r$hovere_hom1_inside
    out$q_qc_shelf_hovere_hom2_inside[i]    <- r$hovere_hom2_inside
  }

  cat(sprintf(
    "[qc_shelf] matched %d / %d status_dt rows (%d summary files had no matching candidate)\n",
    n_matched, n, max(0L, n_found - n_matched)
  ))

  # Flag breakdown summary
  if (n_matched > 0) {
    flag_counts <- table(out$q_qc_shelf_flag[out$q_qc_shelf_ran],
                         useNA = "ifany")
    cat("[qc_shelf] flag distribution among matched candidates:\n")
    for (fn in names(flag_counts)) {
      cat(sprintf("  %-20s %d\n", fn, as.integer(flag_counts[[fn]])))
    }
  }

  cbind(status_dt, out)
}

# =============================================================================
# When sourced, this file defines append_qc_shelf_to_rows() and exits.
# The caller (compute_candidate_status.R) is responsible for invoking it.
# =============================================================================
