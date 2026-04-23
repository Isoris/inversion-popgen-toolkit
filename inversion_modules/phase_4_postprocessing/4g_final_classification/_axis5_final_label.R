#!/usr/bin/env Rscript
# =============================================================================
# _axis5_final_label.R — Axis 5 extension for compute_candidate_status.R
# =============================================================================
#
# This is a STANDALONE MODULE, not a monkey-patch. Source it from
# compute_candidate_status.R (or from a wrapper) to gain Axis 5 support.
#
# Why a separate file: compute_candidate_status.R already has a tight
# 7-axis structure (build_key_spec + 4 axes). Injecting a 5th axis into
# its innards risks breaking downstream consumers that read its exact
# output columns. Instead: read the existing output, append Axis 5 as a
# new column set, and re-write.
#
# USAGE (typical, from a wrapper script):
#
#   # After compute_candidate_status.R has produced candidate_status.tsv:
#   source("inversion_modules/phase_4_postprocessing/4g_final_classification/_axis5_final_label.R")
#   add_axis5(
#     status_tsv   = "candidate_status.tsv",
#     v7_final_dir = file.path(V7_OUT, "final"),
#     out_tsv      = "candidate_status_with_axis5.tsv"
#   )
#
# OR sourcing inside compute_candidate_status.R: add at the bottom, just
# before it writes the status TSV:
#
#   axis5_dir <- Sys.getenv("V7_FINAL_DIR", "")
#   if (nzchar(axis5_dir) && dir.exists(axis5_dir)) {
#     source(file.path(dirname(sys.frame(1)$ofile), "_axis5_final_label.R"))
#     status_dt <- append_axis5_to_rows(status_dt, axis5_dir)
#   }
#
# OUTPUT COLUMNS ADDED:
#   q_overall_structural_class  — v7 final label (e.g. "supported_balanced_inversion_NAHR_like_hypothesis")
#   axis5_weakest_component     — which evidence path was weakest ("group_validation" / "mechanism" / etc.)
#   axis5_justification         — human-readable reasoning string
#   axis5_source                — "v7_final_label" (provenance)
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# -----------------------------------------------------------------------------
# read_final_label — read ONE candidate's final_label.json
# -----------------------------------------------------------------------------
read_final_label <- function(cid, v7_final_dir) {
  p <- file.path(v7_final_dir, cid, "final_label.json")
  if (!file.exists(p)) {
    return(list(
      q_overall_structural_class = NA_character_,
      axis5_weakest_component    = NA_character_,
      axis5_justification        = NA_character_,
      axis5_source               = "no_v7_output"
    ))
  }

  j <- tryCatch(
    jsonlite::fromJSON(p, simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(j)) {
    return(list(
      q_overall_structural_class = NA_character_,
      axis5_weakest_component    = NA_character_,
      axis5_justification        = sprintf("parse_error: %s", p),
      axis5_source               = "parse_error"
    ))
  }

  list(
    q_overall_structural_class = j$q_overall_structural_class %||% NA_character_,
    axis5_weakest_component    = j$weakest_component          %||% NA_character_,
    axis5_justification        = paste(j$justification %||% character(0), collapse = "; "),
    axis5_source               = "v7_final_label"
  )
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -----------------------------------------------------------------------------
# append_axis5_to_rows — given a data.table of candidate rows, append axis5
# -----------------------------------------------------------------------------
append_axis5_to_rows <- function(status_dt, v7_final_dir) {
  stopifnot("candidate_id" %in% names(status_dt))

  axis5_rows <- lapply(status_dt$candidate_id, read_final_label, v7_final_dir = v7_final_dir)
  axis5_dt <- rbindlist(axis5_rows, fill = TRUE)

  cbind(status_dt, axis5_dt)
}

# -----------------------------------------------------------------------------
# add_axis5 — top-level convenience: read existing status TSV, append, write out
# -----------------------------------------------------------------------------
add_axis5 <- function(status_tsv, v7_final_dir, out_tsv = NULL) {
  if (!file.exists(status_tsv)) {
    stop("status_tsv not found: ", status_tsv)
  }
  if (!dir.exists(v7_final_dir)) {
    stop("v7_final_dir not found: ", v7_final_dir)
  }

  status_dt <- fread(status_tsv)
  out_dt    <- append_axis5_to_rows(status_dt, v7_final_dir)

  if (is.null(out_tsv)) out_tsv <- sub("\\.tsv$", "_with_axis5.tsv", status_tsv)
  fwrite(out_dt, out_tsv, sep = "\t")

  cat(sprintf("[axis5] read %d candidates from %s\n", nrow(status_dt), status_tsv))
  cat(sprintf("[axis5] %d had v7 labels (others: NA)\n",
              sum(out_dt$axis5_source == "v7_final_label", na.rm = TRUE)))
  cat(sprintf("[axis5] wrote %s\n", out_tsv))

  invisible(out_dt)
}

# -----------------------------------------------------------------------------
# CLI entry point
# -----------------------------------------------------------------------------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) {
    status_tsv   <- args[1]
    v7_final_dir <- args[2]
    out_tsv      <- if (length(args) >= 3) args[3] else NULL
    add_axis5(status_tsv, v7_final_dir, out_tsv)
  }
}
