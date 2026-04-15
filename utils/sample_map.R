#!/usr/bin/env Rscript

# =============================================================================
# sample_map.R — Universal sample name resolver
#
# Solves the Ind0/Ind1 vs CGA009/CGA010 mismatch throughout the pipeline.
#
# BEAGLE files use Ind0, Ind1, Ind2... as column names.
# Clair3, DELLY, ngsRelate, etc. use CGA009, CGA010, CGA021...
# Precomp inherits Ind0 from BEAGLE.
# Every downstream script needs to know both.
#
# Usage (source from any R script):
#   source("utils/sample_map.R")
#   smap <- load_sample_map()    # auto-detects from config
#   smap$to_real("Ind0")         # → "CGA009"
#   smap$to_ind("CGA009")        # → "Ind0"
#   smap$to_real_vec(c("Ind0","Ind1"))  # → c("CGA009","CGA010")
#   smap$to_ind_vec(c("CGA009","CGA010"))  # → c("Ind0","Ind1")
#   smap$detect("Ind0")          # → "ind"
#   smap$detect("CGA009")        # → "real"
#   smap$ensure_real(sample_names)  # converts if needed
#   smap$ensure_ind(sample_names)   # converts if needed
#   smap$real_names               # all 226 real names in order
#   smap$ind_names                # all 226 Ind names in order
#   smap$n                        # 226
#
# The mapping is built from the sample list file (one name per line,
# in the same order as BEAGLE columns). The Ind index is 0-based
# (Ind0 = first sample).
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

load_sample_map <- function(sample_list = NULL) {

  # ── Find sample list ──
  if (is.null(sample_list) || !file.exists(sample_list)) {
    # Try environment variable from config
    sample_list <- Sys.getenv("SAMPLES_IND", "")
  }
  if (!nzchar(sample_list) || !file.exists(sample_list)) {
    # Try common locations
    candidates <- c(
      "het_roh/01_inputs_check/samples.ind",
      "popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv",
      file.path(Sys.getenv("BASE", "."), "het_roh/01_inputs_check/samples.ind"),
      file.path(Sys.getenv("BASE", "."), "popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv")
    )
    for (f in candidates) {
      if (file.exists(f)) { sample_list <- f; break }
    }
  }

  if (!file.exists(sample_list)) {
    stop("[sample_map] Cannot find sample list file. Tried:\n  ",
         paste(c(sample_list, candidates), collapse = "\n  "),
         "\n  Set SAMPLES_IND env var or pass path to load_sample_map()")
  }

  # ── Read ──
  real_names <- trimws(readLines(sample_list))
  real_names <- real_names[nzchar(real_names)]
  n <- length(real_names)
  ind_names <- paste0("Ind", seq_len(n) - 1L)

  # ── Build lookup tables ──
  real_to_ind <- setNames(ind_names, real_names)
  ind_to_real <- setNames(real_names, ind_names)

  # ── Helper functions ──

  detect_format <- function(x) {
    x <- x[1]
    if (grepl("^Ind[0-9]+$", x)) return("ind")
    if (grepl("^CGA[0-9]+$", x)) return("real")
    # Unknown format — check against both lookups
    if (x %in% ind_names) return("ind")
    if (x %in% real_names) return("real")
    return("unknown")
  }

  to_real <- function(x) {
    if (x %in% names(ind_to_real)) ind_to_real[[x]] else x
  }

  to_ind <- function(x) {
    if (x %in% names(real_to_ind)) real_to_ind[[x]] else x
  }

  to_real_vec <- function(xs) {
    ifelse(xs %in% names(ind_to_real), ind_to_real[xs], xs)
  }

  to_ind_vec <- function(xs) {
    ifelse(xs %in% names(real_to_ind), real_to_ind[xs], xs)
  }

  ensure_real <- function(xs) {
    fmt <- detect_format(xs)
    if (fmt == "ind") to_real_vec(xs)
    else xs
  }

  ensure_ind <- function(xs) {
    fmt <- detect_format(xs)
    if (fmt == "real") to_ind_vec(xs)
    else xs
  }

  # Also rename data.table columns in-place
  rename_dt_columns <- function(dt, from = "ind", prefix = "PC_1_") {
    cols <- grep(paste0("^", prefix), names(dt), value = TRUE)
    if (length(cols) == 0) return(dt)
    old_names <- sub(paste0("^", prefix), "", cols)
    fmt <- detect_format(old_names)
    if (fmt == from) {
      new_names <- if (from == "ind") to_real_vec(old_names) else to_ind_vec(old_names)
      new_cols <- paste0(prefix, new_names)
      setnames(dt, cols, new_cols)
    }
    dt
  }

  message("[sample_map] Loaded ", n, " samples from ", basename(sample_list),
          " (", real_names[1], " = ", ind_names[1], " ... ",
          real_names[n], " = ", ind_names[n], ")")

  list(
    real_names = real_names,
    ind_names = ind_names,
    n = n,
    source_file = sample_list,
    detect = detect_format,
    to_real = to_real,
    to_ind = to_ind,
    to_real_vec = to_real_vec,
    to_ind_vec = to_ind_vec,
    ensure_real = ensure_real,
    ensure_ind = ensure_ind,
    rename_dt_columns = rename_dt_columns,
    real_to_ind = real_to_ind,
    ind_to_real = ind_to_real
  )
}

# =============================================================================
# Bash helper: generate mapping TSV
# Run standalone: Rscript sample_map.R <sample_list> <output.tsv>
# =============================================================================

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    smap <- load_sample_map(args[1])
    out <- if (length(args) >= 2) args[2] else "sample_name_mapping.tsv"
    fwrite(data.table(
      ind_name = smap$ind_names,
      real_name = smap$real_names,
      index_0based = seq_len(smap$n) - 1L
    ), out, sep = "\t")
    message("[sample_map] Written: ", out, " (", smap$n, " rows)")
  } else {
    cat("Usage: Rscript sample_map.R <sample_list.txt> [output.tsv]\n")
    cat("  Generates Ind0<->CGA009 mapping table\n")
  }
}
