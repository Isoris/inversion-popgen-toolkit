#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR00b_aggregate_pestPG_on_dosage_windows.R
# =============================================================================
# Phase 2 / 2f_theta_discovery — pestPG → dosage-grid TSV bridge.
#
# Reads each sample's existing win50000.step10000 pestPG file and assigns
# its per-window θπ to the **nearest dosage scrubber window by midpoint**.
# Output is keyed on dosage window indices, ready for STEP_TR01 to pivot
# into a samples × windows matrix for local PCA.
#
# This is the v3 path — replaces v2's per-site aggregation. Both inputs
# already exist on LANTA:
#   1. ${PESTPG_DIR}/${SAMPLE}.win50000.step10000.pestPG  (from MODULE_3)
#   2. ${DOSAGE_WIN_BED_DIR}/${CHROM}/windows.bed         (from 2a_local_pca)
# No upstream rerun required.
#
# See `window_grid_alignment_v3.md` for full architectural reasoning.
#
# Output schema (long format, one row per sample × dosage window):
#   sample  chrom  window_idx  start_bp  end_bp  theta_pi  n_sites
#
# Usage (sourced 00_theta_config.sh first):
#   Rscript STEP_TR00b_aggregate_pestPG_on_dosage_windows.R <CHROM>
# Example:
#   Rscript STEP_TR00b_aggregate_pestPG_on_dosage_windows.R C_gar_LG28
#
# Walltime: ~5 min for 226 samples × 1 chromosome on LANTA scratch.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ── Args + config ─────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript STEP_TR00b_aggregate_pestPG_on_dosage_windows.R <CHROM>")
}
CHROM <- args[1]

PESTPG_DIR         <- Sys.getenv("PESTPG_DIR",         unset = NA)
DOSAGE_WIN_BED_DIR <- Sys.getenv("DOSAGE_WIN_BED_DIR", unset = NA)
SAMPLE_LIST        <- Sys.getenv("SAMPLE_LIST",        unset = NA)
THETA_TSV_DIR      <- Sys.getenv("THETA_TSV_DIR",      unset = NA)
LOG_DIR            <- Sys.getenv("LOG_DIR",            unset = ".")
PESTPG_SCALE       <- Sys.getenv("PESTPG_SCALE",       unset = "win50000.step10000")

stopifnot(!is.na(PESTPG_DIR), !is.na(DOSAGE_WIN_BED_DIR),
          !is.na(SAMPLE_LIST), !is.na(THETA_TSV_DIR))

# ── Inputs ────────────────────────────────────────────────────────────────
dosage_bed <- file.path(DOSAGE_WIN_BED_DIR, CHROM, "windows.bed")
out_tsv    <- file.path(THETA_TSV_DIR,
                        sprintf("theta_dgrid.%s.tsv.gz", CHROM))

if (!file.exists(dosage_bed)) {
  stop("[STEP_TR00b] Dosage windows BED missing: ", dosage_bed,
       "\n  Expected layout: ${DOSAGE_WIN_BED_DIR}/${CHROM}/windows.bed",
       "\n  Produced upstream by 2a_local_pca's STEP09b dense window registry.")
}

if (!file.exists(SAMPLE_LIST)) {
  stop("[STEP_TR00b] Sample list missing: ", SAMPLE_LIST)
}

samples <- readLines(SAMPLE_LIST)
samples <- samples[nchar(samples) > 0]
message("[STEP_TR00b] CHROM=", CHROM, "  n_samples=", length(samples))
message("           dosage_bed = ", dosage_bed)
message("           output     = ", out_tsv)

# ── Load dosage window grid for this chromosome ───────────────────────────
# Expected columns: chrom start_bp end_bp window_idx n_snps
# (window_idx is 0-based; matches the dosage scrubber's state.cur)
dosage <- fread(dosage_bed,
                col.names = c("chrom", "start_bp", "end_bp",
                              "window_idx", "n_snps"))
dosage <- dosage[chrom == CHROM]
if (nrow(dosage) == 0) {
  stop("[STEP_TR00b] No dosage windows found for CHROM=", CHROM,
       " in ", dosage_bed)
}
dosage[, mid_bp := as.integer((start_bp + end_bp) / 2L)]
setkey(dosage, mid_bp)
message("[STEP_TR00b] n_dosage_windows = ", nrow(dosage))

# ── Helper: load one sample's pestPG, restrict to CHROM, return mid_bp/tP ─
load_pestpg <- function(sample) {
  pestpg_file <- file.path(PESTPG_DIR,
                           sprintf("%s.%s.pestPG", sample, PESTPG_SCALE))
  if (!file.exists(pestpg_file)) {
    return(NULL)   # caller skips this sample with a warning
  }
  # pestPG header line starts with #; columns:
  #   #(indexStart,indexStop)(firstPos,lastPos)(WinStart,WinStop) Chr WinCenter
  #   tW tP tF tH tL Tajima fuf fud fayh zeng nSites
  # The first column is a parenthetical artifact we ignore.
  pp <- tryCatch(
    fread(cmd = sprintf("zcat -f %s | grep -v '^#'",
                        shQuote(pestpg_file)),
          header = FALSE,
          fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(pp) || nrow(pp) == 0) return(NULL)
  # Robust column rename — there may be 14 cols or 15 depending on ANGSD version
  if (ncol(pp) >= 14) {
    setnames(pp, 1:14,
             c("ix", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL",
               "Tajima", "fuf", "fud", "fayh", "zeng", "nSites"))
  } else {
    return(NULL)   # malformed; skip
  }
  pp <- pp[Chr == CHROM,
           .(mid_bp = as.integer(WinCenter),
             theta_pi = as.numeric(tP),
             n_sites = as.integer(nSites))]
  if (nrow(pp) == 0) return(NULL)
  setkey(pp, mid_bp)
  pp
}

# ── Loop over samples; nearest-midpoint join each one's pestPG to dosage ──
out_rows <- vector("list", length(samples))
n_ok <- 0L
n_skip <- 0L

for (i in seq_along(samples)) {
  s <- samples[i]
  pp <- load_pestpg(s)
  if (is.null(pp)) {
    warning("[STEP_TR00b] Missing/malformed pestPG for sample ", s,
            " — skipping")
    n_skip <- n_skip + 1L
    next
  }
  # Nearest-midpoint join: each dosage window picks the pestPG row
  # whose WinCenter is closest in bp space. roll = "nearest" handles
  # edge windows cleanly by clamping to the first / last pestPG row.
  joined <- pp[dosage, roll = "nearest", on = "mid_bp"]
  joined[, sample := s]
  joined[, chrom  := CHROM]
  out_rows[[i]] <- joined[, .(sample, chrom, window_idx,
                              start_bp, end_bp,
                              theta_pi, n_sites)]
  n_ok <- n_ok + 1L
  if (i %% 50 == 0) {
    message("[STEP_TR00b] processed ", i, " / ", length(samples),
            " samples")
  }
}

if (n_ok == 0) {
  stop("[STEP_TR00b] No samples produced output. Check PESTPG_DIR + sample naming.")
}

result <- rbindlist(out_rows, use.names = TRUE)
message("[STEP_TR00b] joined: n_ok=", n_ok, " n_skip=", n_skip,
        "  total_rows=", nrow(result))

# ── Write output TSV (gzip-compressed) ───────────────────────────────────
fwrite(result, out_tsv, sep = "\t", compress = "gzip", na = "")

fi <- file.info(out_tsv)
message("[STEP_TR00b] Wrote ", out_tsv,
        " (", round(fi$size / 1024 / 1024, 2), " MB)")

# Sanity: every sample × dosage window should be present
expected <- as.integer(n_ok) * nrow(dosage)
if (nrow(result) != expected) {
  warning("[STEP_TR00b] Row count mismatch: got ", nrow(result),
          " expected ", expected, " (samples × dosage_windows)")
}

message("[STEP_TR00b] DONE — chrom=", CHROM,
        " samples=", n_ok,
        " windows=", nrow(dosage))
