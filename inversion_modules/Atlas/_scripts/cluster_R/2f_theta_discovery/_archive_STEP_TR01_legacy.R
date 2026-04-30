#!/usr/bin/env Rscript
# =============================================================================
# STEP_TR01_per_window_theta_pi.R
# =============================================================================
# Phase 2 / 2f_theta_discovery — Part A input layer.
#
# Reads:
#   ${THETA_TSV_DIR}/theta.<CHROM>.win50000.step10000.tsv.gz
#     (long format: sample chrom win_start win_end win_center tP nSites)
#
# Writes:
#   ${OUT_PER_WINDOW_DIR}/theta_pi_per_window.<CHROM>.json.gz
#     (wide matrix: 226 samples × n_windows, with window grid + sample list)
#
# Schema: matches theta_pi_data_layer_spec_v0.md §2 (schema_version=1).
# Atlas reads this at pca_scrubber_v4.html line 24316 (detectSchemaAndLayers).
#
# This is the SIMPLEST of the three TR scripts — it's a pivot + JSON
# serialization. No PCA, no clustering. Just shape the matrix the atlas
# expects.
#
# Usage:
#   Rscript STEP_TR01_per_window_theta_pi.R <CHROM>
# Example:
#   Rscript STEP_TR01_per_window_theta_pi.R C_gar_LG28
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# ── Args + config ─────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript STEP_TR01_per_window_theta_pi.R <CHROM>")
CHROM <- args[1]

# Config is sourced by the launcher; here we read env-vars set by 00_theta_config.sh
THETA_TSV_DIR        <- Sys.getenv("THETA_TSV_DIR",        unset = NA)
SCALE_LABEL          <- Sys.getenv("SCALE_LABEL",          unset = NA)
WIN_BP               <- as.integer(Sys.getenv("WIN_BP",    unset = "50000"))
STEP_BP              <- as.integer(Sys.getenv("STEP_BP",   unset = "10000"))
OUT_PER_WINDOW_DIR   <- Sys.getenv("OUT_PER_WINDOW_DIR",   unset = NA)
SAMPLE_LIST          <- Sys.getenv("SAMPLE_LIST",          unset = NA)
SCHEMA_VERSION       <- as.integer(Sys.getenv("THETA_JSON_SCHEMA_VERSION",
                                              unset = "1"))
LOG_DIR              <- Sys.getenv("LOG_DIR",              unset = ".")

stopifnot(!is.na(THETA_TSV_DIR), !is.na(OUT_PER_WINDOW_DIR), !is.na(SCALE_LABEL))

input_tsv  <- file.path(THETA_TSV_DIR, sprintf("theta.%s.%s.tsv.gz",
                                               CHROM, SCALE_LABEL))
output_json <- file.path(OUT_PER_WINDOW_DIR,
                         sprintf("theta_pi_per_window.%s.json.gz", CHROM))

if (!file.exists(input_tsv)) {
  stop("[STEP_TR01] Input TSV missing: ", input_tsv,
       "\n  Run STEP_TR00_aggregate_theta.sh first.")
}

message("[STEP_TR01] CHROM=", CHROM)
message("           input  = ", input_tsv)
message("           output = ", output_json)

# ── Read sample list (canonical 226 order) ────────────────────────────────
if (file.exists(SAMPLE_LIST)) {
  samples_canonical <- readLines(SAMPLE_LIST)
  samples_canonical <- samples_canonical[nchar(samples_canonical) > 0]
} else {
  samples_canonical <- NULL
  message("[STEP_TR01] WARN: SAMPLE_LIST not found; will use samples present in TSV in their natural order.")
}

# ── Load the per-sample θπ long table ─────────────────────────────────────
dt <- fread(input_tsv)
# Expected columns: sample chrom win_start win_end win_center tP nSites
required <- c("sample", "win_start", "win_end", "tP")
miss <- setdiff(required, names(dt))
if (length(miss) > 0) {
  stop("[STEP_TR01] Input TSV missing required columns: ",
       paste(miss, collapse = ", "))
}

# nSites column may be missing in older outputs; tolerate
if (!"nSites" %in% names(dt)) {
  dt[, nSites := NA_integer_]
  message("[STEP_TR01] NOTE: nSites column missing — output will use null.")
}

# Filter to this chromosome (defensive — TSV should already be per-chrom)
if ("chrom" %in% names(dt)) {
  dt <- dt[chrom == CHROM]
}

if (nrow(dt) == 0) {
  stop("[STEP_TR01] No rows for CHROM=", CHROM, " in ", input_tsv)
}

# ── Build the canonical window grid ───────────────────────────────────────
# Each window is uniquely identified by (win_start, win_end). They should
# already be on a regular grid from thetaStat -win/-step. Sort by start.
windows_dt <- unique(dt[, .(win_start, win_end)])
setorder(windows_dt, win_start)
windows_dt[, idx := .I - 1L]  # 0-based per the spec

n_windows <- nrow(windows_dt)
message("[STEP_TR01] n_windows = ", n_windows)

# Sanity check spacing
if (n_windows >= 2) {
  step_actual <- windows_dt$win_start[2] - windows_dt$win_start[1]
  if (step_actual != STEP_BP) {
    message("[STEP_TR01] WARN: actual step (", step_actual,
            ") != configured STEP_BP (", STEP_BP, ")")
  }
  win_actual <- windows_dt$win_end[1] - windows_dt$win_start[1]
  if (win_actual != WIN_BP) {
    message("[STEP_TR01] WARN: actual window size (", win_actual,
            ") != configured WIN_BP (", WIN_BP, ")")
  }
}

# Look up nSites per window (mean across samples, since pestPG nSites is
# per-sample but should be identical across samples for the same window
# given fixed callable_sites — we record the median to be safe).
nsites_per_window <- dt[, .(n_sites = as.integer(round(median(nSites, na.rm = TRUE)))),
                       by = .(win_start, win_end)]
windows_dt <- merge(windows_dt, nsites_per_window,
                    by = c("win_start", "win_end"), all.x = TRUE)
setorder(windows_dt, idx)

# ── Build the canonical sample list ───────────────────────────────────────
samples_in_data <- sort(unique(dt$sample))
if (!is.null(samples_canonical)) {
  # Use the canonical order; warn on mismatches
  in_data_only <- setdiff(samples_in_data, samples_canonical)
  in_canon_only <- setdiff(samples_canonical, samples_in_data)
  if (length(in_data_only) > 0) {
    message("[STEP_TR01] WARN: ", length(in_data_only),
            " samples in TSV not in canonical list: ",
            paste(head(in_data_only, 3), collapse = ", "), "...")
  }
  if (length(in_canon_only) > 0) {
    message("[STEP_TR01] WARN: ", length(in_canon_only),
            " canonical samples missing from TSV (will be all-null rows).")
  }
  samples <- samples_canonical
} else {
  samples <- samples_in_data
}

n_samples <- length(samples)
message("[STEP_TR01] n_samples = ", n_samples)

# ── Pivot to wide: rows = samples, cols = windows ─────────────────────────
# Build a (sample, idx) → tP lookup, then materialize the matrix.
# Use idx (0-based window index) not win_start to keep keys compact.
dt_indexed <- merge(dt, windows_dt[, .(win_start, win_end, idx)],
                    by = c("win_start", "win_end"), all.x = TRUE)

# dcast for the matrix
mat_dt <- dcast(dt_indexed, sample ~ idx, value.var = "tP",
                fun.aggregate = mean)  # mean over duplicates if any

# Re-order rows to canonical sample order; missing samples become NA rows
mat <- matrix(NA_real_, nrow = n_samples, ncol = n_windows,
              dimnames = list(samples, NULL))
for (i in seq_len(nrow(mat_dt))) {
  s <- mat_dt$sample[i]
  ridx <- match(s, samples)
  if (!is.na(ridx)) {
    # Columns of mat_dt other than 'sample' are window indices as character
    vals <- as.numeric(mat_dt[i, -1, with = FALSE])
    # Window indices are character "0","1",...; reorder to numeric
    cidx <- as.integer(names(mat_dt)[-1]) + 1L  # +1 for R's 1-based cols
    mat[ridx, cidx] <- vals
  }
}

# Convert NA to NULL in JSON output by leaving as NA (jsonlite na='null')
# Round to 6 sig digits to keep JSON compact (θπ values are ~1e-3..1e-5 range)
mat_rounded <- signif(mat, digits = 6)

# ── Build the JSON object per the spec ────────────────────────────────────
windows_json <- lapply(seq_len(n_windows), function(i) {
  list(
    idx       = windows_dt$idx[i],          # 0-based
    start_bp  = windows_dt$win_start[i],
    end_bp    = windows_dt$win_end[i],
    n_sites   = if (is.na(windows_dt$n_sites[i])) NULL
                else windows_dt$n_sites[i]
  )
})

theta_pi_rows <- lapply(seq_len(n_samples), function(i) {
  row <- mat_rounded[i, ]
  # jsonlite handles NA → null when na='null' is set
  unname(row)
})

# NOTE on field name: atlas detectSchemaAndLayers (pca_scrubber_v4.html
# line 24319) looks for `.values` on this layer object, NOT `.theta_pi`.
# Spec v0 §2 used `theta_pi` as the values field name, but the atlas was
# written first with `values` to match its existing layer vocabulary
# (sibling layers like tracks also use `.values`). We emit `values` here
# as the canonical field name. See theta_pi_data_layer_spec_v1.md for the
# spec correction. (For a brief deprecation period, also emit
# `theta_pi` as an alias so any pre-spec-v1 code that hard-coded the old
# name keeps working — atlas reads `values` regardless.)
out_obj <- list(
  schema_version = SCHEMA_VERSION,
  layer          = "theta_pi_per_window",
  chrom          = CHROM,
  scale          = SCALE_LABEL,
  n_samples      = n_samples,
  n_windows      = n_windows,
  windows        = windows_json,
  samples        = samples,
  values         = theta_pi_rows,    # ← canonical field name (atlas-aligned)
  theta_pi       = theta_pi_rows     # ← legacy alias (spec v0 compat)
)

# ── Write gzip-compressed JSON ────────────────────────────────────────────
# auto_unbox=TRUE so single-element vectors don't render as arrays
gz <- gzfile(output_json, open = "w")
on.exit(close(gz), add = TRUE)
writeLines(toJSON(out_obj, auto_unbox = TRUE, na = "null", digits = NA), gz)

# Verify file size
fi <- file.info(output_json)
message("[STEP_TR01] Wrote ", output_json,
        " (", round(fi$size / 1024 / 1024, 2), " MB)")
message("[STEP_TR01] DONE — schema_version=", SCHEMA_VERSION,
        " layer=theta_pi_per_window chrom=", CHROM)

# Round-trip sanity check: re-read and verify shape
test <- fromJSON(output_json, simplifyVector = FALSE)
stopifnot(test$layer == "theta_pi_per_window",
          test$chrom == CHROM,
          length(test$samples) == n_samples,
          length(test$windows) == n_windows,
          length(test$values) == n_samples,           # atlas-canonical
          length(test$theta_pi) == n_samples)         # legacy alias
message("[STEP_TR01] Round-trip OK.")
