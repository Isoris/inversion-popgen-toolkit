#!/usr/bin/env Rscript

# =============================================================================
# export_ghsl_to_json_v2.R
#
# Phase 4a "chapter end" emit script. Reads the GHSL v6 classifier outputs
# and emits one JSON per chromosome at
#   <out_dir>/<chrom>/<chrom>_phase4a_ghsl.json
# for the local PCA scrubber to load as enrichment.
#
# Schema: SCHEMA_V2 (see scrubber/docs/SCHEMA_V2.md).
#
# WHAT v2 ADDS OVER v1
# --------------------
# v1 shipped a slim ghsl_heatmap (one rolling scale, optionally downsampled)
# plus karyotype runs and three tracks. v2 ships the full dense
# (sample × window × scale) panel uncompressed plus per-K stripe assignments
# at each scale for K = 2..6, all pre-computed cluster-side. Total JSON size
# is large (50-300 MB per chromosome) but the scrubber is local and can hold
# it; the wins are:
#   - K-stripe heatmaps render at any K without re-running kmeans in JS
#   - Per-sample focal±N viewer reads directly from the panel
#   - PCA recolouring by GHSL is a flat lookup
#   - L2-as-triangle-substitute aggregation happens browser-side (no
#     round trip), works for any L2 the user defines later
#
# Inputs (cluster-side, per chromosome):
#   --annot_dir   Directory with <chr>.ghsl_v6.annot.rds + .karyotypes.rds
#   --panel_dir   Directory with <chr>.ghsl_v6.per_sample.rds  (REQUIRED in v2)
#   --chrom       Chromosome name (must match the RDS files' chrom field)
#   --out_dir     Output directory. Writes
#                 <out_dir>/<chrom>/<chrom>_phase4a_ghsl.json
#                 Recommended: inversion_modules/scrubber/data/
#
# Optional:
#   --max_k 6     Maximum K for stripe partitioning (default 6, range 2..6
#                 always emitted from K=2 up to --max_k inclusive).
#   --scales s10,s20,s30,s40,s50,s100
#                 Which rolling scales to ship in the dense panel. Default
#                 is everything in the panel RDS. Use this to drop scales
#                 you don't need (cuts JSON size linearly).
#   --primary_scale s50
#                 Which scale to use for K-stripe partitioning and the
#                 catalogue summary stats. Default s50 (the classifier's
#                 primary scale).
#
# Layers this script registers in <chrom>_phase4a_ghsl.json:
#
#   tracks                       Contributes ghsl_div_median, ghsl_score,
#                                ghsl_rank_stability per window.
#   ghsl_karyotype_runs          Per-sample stable LOW/HIGH runs.
#   ghsl_panel                   Dense (sample × window × scale) divergence
#                                and rank. The big payload.
#   ghsl_kstripes                Per-K stripe assignments + summary stats
#                                at the primary scale. Small.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# PARSE ARGS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(paste(
    "Usage: Rscript export_ghsl_to_json_v2.R --annot_dir <dir> --panel_dir <dir>",
    "                                        --chrom <chr> --out_dir <dir>",
    "                                        [--max_k 6]",
    "                                        [--scales s10,s20,s30,s40,s50,s100]",
    "                                        [--primary_scale s50]",
    sep = "\n"
  ), "\n")
  quit(status = 1)
}

ANNOT_DIR      <- NULL
PANEL_DIR      <- NULL
CHROM          <- NULL
OUT_DIR        <- NULL
MAX_K          <- 6L
SCALES_FILTER  <- NULL    # NULL means "use everything available"
PRIMARY_SCALE  <- "s50"

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--annot_dir"     && i < length(args)) { ANNOT_DIR     <- args[i + 1]; i <- i + 2L }
  else if (a == "--panel_dir"     && i < length(args)) { PANEL_DIR     <- args[i + 1]; i <- i + 2L }
  else if (a == "--chrom"         && i < length(args)) { CHROM         <- args[i + 1]; i <- i + 2L }
  else if (a == "--out_dir"       && i < length(args)) { OUT_DIR       <- args[i + 1]; i <- i + 2L }
  else if (a == "--max_k"         && i < length(args)) { MAX_K         <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--scales"        && i < length(args)) { SCALES_FILTER <- strsplit(args[i + 1], ",")[[1]]; i <- i + 2L }
  else if (a == "--primary_scale" && i < length(args)) { PRIMARY_SCALE <- args[i + 1]; i <- i + 2L }
  else if (a == "-h" || a == "--help")                 { usage() }
  else { i <- i + 1L }
}

if (is.null(ANNOT_DIR) || is.null(PANEL_DIR) || is.null(CHROM) || is.null(OUT_DIR)) usage()
if (MAX_K < 2L || MAX_K > 10L) stop("--max_k must be in 2..10")

# =============================================================================
# LOAD INPUTS
# =============================================================================

annot_rds <- file.path(ANNOT_DIR, paste0(CHROM, ".ghsl_v6.annot.rds"))
karyo_rds <- file.path(ANNOT_DIR, paste0(CHROM, ".ghsl_v6.karyotypes.rds"))
panel_rds <- file.path(PANEL_DIR, paste0(CHROM, ".ghsl_v6.per_sample.rds"))

if (!file.exists(annot_rds)) stop("annot RDS not found: ", annot_rds)
if (!file.exists(karyo_rds)) stop("karyotypes RDS not found: ", karyo_rds)
if (!file.exists(panel_rds)) stop("panel RDS not found: ", panel_rds)

cat("[GHSL_EXPORT v2] Loading annot RDS: ", annot_rds, "\n", sep = "")
annot <- as.data.table(readRDS(annot_rds))
cat("[GHSL_EXPORT v2]   ", nrow(annot), " windows\n", sep = "")

cat("[GHSL_EXPORT v2] Loading karyotypes RDS: ", karyo_rds, "\n", sep = "")
karyo <- as.data.table(readRDS(karyo_rds))
cat("[GHSL_EXPORT v2]   ", nrow(karyo), " stable runs\n", sep = "")

cat("[GHSL_EXPORT v2] Loading panel RDS: ", panel_rds, "\n", sep = "")
panel <- as.data.table(readRDS(panel_rds))
panel_meta <- attr(panel, "ghsl_panel_meta")
cat("[GHSL_EXPORT v2]   ", nrow(panel), " (sample × window) rows\n", sep = "")

n_windows <- nrow(annot)

# Sample order from panel meta if present, else from data
sample_order <- if (!is.null(panel_meta$sample_order)) {
  panel_meta$sample_order
} else {
  unique(panel$sample_id)
}
n_samples <- length(sample_order)

# Window order: by window_idx, ascending
win_idx_unique <- sort(unique(panel$window_idx))
n_win_panel    <- length(win_idx_unique)

cat("[GHSL_EXPORT v2] Samples:        ", n_samples, "\n", sep = "")
cat("[GHSL_EXPORT v2] Panel windows:  ", n_win_panel, "\n", sep = "")
cat("[GHSL_EXPORT v2] Annot windows:  ", n_windows, "\n", sep = "")

# Detect which scales are present in the panel
all_div_cols   <- grep("^div_roll_s", names(panel), value = TRUE)
all_rank_cols  <- grep("^rank_in_cohort_s", names(panel), value = TRUE)
all_band_cols  <- grep("^rank_band_s", names(panel), value = TRUE)
scales_present <- sub("^div_roll_", "", all_div_cols)
cat("[GHSL_EXPORT v2] Scales in panel: ", paste(scales_present, collapse = ", "), "\n", sep = "")

# Apply scales filter if user passed one
if (!is.null(SCALES_FILTER)) {
  scales_use <- intersect(scales_present, SCALES_FILTER)
  missing    <- setdiff(SCALES_FILTER, scales_present)
  if (length(missing) > 0) {
    cat("[GHSL_EXPORT v2] WARNING: requested scales not in panel: ",
        paste(missing, collapse = ", "), "\n", sep = "")
  }
} else {
  scales_use <- scales_present
}
if (length(scales_use) == 0) stop("No usable scales after filtering")
cat("[GHSL_EXPORT v2] Scales emitted: ", paste(scales_use, collapse = ", "), "\n", sep = "")

if (!PRIMARY_SCALE %in% scales_use) {
  cat("[GHSL_EXPORT v2] WARNING: primary scale ", PRIMARY_SCALE,
      " not in emitted scales. Falling back to first: ", scales_use[1], "\n", sep = "")
  PRIMARY_SCALE <- scales_use[1]
}
cat("[GHSL_EXPORT v2] Primary scale:  ", PRIMARY_SCALE, "\n", sep = "")

# =============================================================================
# BUILD: tracks contribution
# =============================================================================
# Per-window aggregates aligned to annot's window order. The scrubber merges
# these into the existing tracks dict additively (see v3.22 scrubber update).

make_track <- function(values_raw, pos_bp = NULL) {
  values <- suppressWarnings(as.numeric(values_raw))
  finite <- values[is.finite(values)]
  rng <- if (length(finite) > 0) range(finite) else c(NA_real_, NA_real_)
  out <- list(
    min = if (is.finite(rng[1])) round(rng[1], 6) else NA_real_,
    max = if (is.finite(rng[2])) round(rng[2], 6) else NA_real_,
    values = round(values, 6)
  )
  if (!is.null(pos_bp)) out$pos_bp <- as.integer(pos_bp)
  out
}

# Per-window mid-bp grid for GHSL annot. The scrubber renders tracks at these
# positions when track.pos_bp is present (track-renderer change in v3.28+),
# decoupling the GHSL window grid from the local-PCA precomp grid.
ghsl_pos_bp <- if ("start_bp" %in% names(annot) && "end_bp" %in% names(annot)) {
  as.integer(round((annot$start_bp + annot$end_bp) / 2))
} else NULL

tracks_list <- list()
if ("div_median"     %in% names(annot)) tracks_list[["ghsl_div_median"]]     <- make_track(annot$div_median,     ghsl_pos_bp)
if ("ghsl_v6_score"  %in% names(annot)) tracks_list[["ghsl_score"]]          <- make_track(annot$ghsl_v6_score,  ghsl_pos_bp)
if ("rank_stability" %in% names(annot)) tracks_list[["ghsl_rank_stability"]] <- make_track(annot$rank_stability, ghsl_pos_bp)

cat("[GHSL_EXPORT v2] Tracks emitted: ", paste(names(tracks_list), collapse = ", "), "\n", sep = "")

# =============================================================================
# BUILD: ghsl_karyotype_runs layer
# =============================================================================

karyo_runs_list <- list()
if (nrow(karyo) > 0) {
  for (i in seq_len(nrow(karyo))) {
    r <- karyo[i]
    karyo_runs_list[[length(karyo_runs_list) + 1L]] <- list(
      sample_id     = as.character(r$sample_id),
      start_bp      = if (!is.null(r$start_bp))      as.integer(r$start_bp)      else NA_integer_,
      end_bp        = if (!is.null(r$end_bp))        as.integer(r$end_bp)        else NA_integer_,
      window_start  = if (!is.null(r$window_start))  as.integer(r$window_start)  else NA_integer_,
      window_end    = if (!is.null(r$window_end))    as.integer(r$window_end)    else NA_integer_,
      n_windows     = if (!is.null(r$n_windows))     as.integer(r$n_windows)     else NA_integer_,
      call          = as.character(r$call),
      mean_rank     = if (!is.null(r$mean_rank))     round(as.numeric(r$mean_rank), 4) else NA_real_,
      rolling_scale = if (!is.null(r$rolling_scale)) as.character(r$rolling_scale) else NA_character_
    )
  }
}
cat("[GHSL_EXPORT v2] Karyotype runs: ", length(karyo_runs_list), "\n", sep = "")

# =============================================================================
# BUILD: ghsl_panel layer (the big one)
# =============================================================================
# Layout: one nested object with samples + windows + per-scale dense matrices.
# Memory: each matrix is n_samples × n_windows × Float, 4-decimal rounded.
#
# {
#   n_samples: 226,
#   n_windows: 4302,
#   samples:        ["CGA001", ...],         length n_samples
#   window_idx:     [1, 2, ..., n_windows],  panel-window index (1-based)
#   global_window_id: [...],                 length n_windows  (alignment key)
#   start_bp:       [...],                   length n_windows
#   end_bp:         [...],                   length n_windows
#   pos_mb:         [...],                   length n_windows
#   scales: ["s10", "s20", ...],
#   div_roll: {                              one matrix per scale
#     s10:  [[s1w1, s1w2, ...], [s2w1, ...], ...],
#     s20:  ...,
#     ...
#   },
#   rank_in_cohort: {                        same shape
#     s10:  [[...]], ...
#   },
#   rank_band: {                             "LOW"|"MID"|"HIGH"|null
#     s10:  [[...]], ...
#   }
# }

cat("[GHSL_EXPORT v2] Building dense panel matrices...\n", sep = "")
t_panel <- proc.time()

# Pre-index sample / window for fast assignment
samp_to_row <- setNames(seq_along(sample_order), sample_order)
win_to_col  <- setNames(seq_along(win_idx_unique), as.character(win_idx_unique))

# Filter panel to known samples / windows once (defensive)
panel_use <- panel[sample_id %in% sample_order & window_idx %in% win_idx_unique]

# Pre-extract row / col indices once — reused across all scale columns
rows <- samp_to_row[panel_use$sample_id]
cols <- win_to_col[as.character(panel_use$window_idx)]
good <- !is.na(rows) & !is.na(cols)
rows <- rows[good]; cols <- cols[good]
fill_idx <- cbind(rows, cols)

build_matrix_list <- function(values_vec) {
  M <- matrix(NA_real_, nrow = n_samples, ncol = n_win_panel)
  M[fill_idx] <- values_vec[good]
  M_round <- round(M, 4)
  out <- vector("list", n_samples)
  for (s in seq_len(n_samples)) out[[s]] <- M_round[s, ]
  out
}

build_band_matrix_list <- function(band_vec) {
  M <- matrix(NA_character_, nrow = n_samples, ncol = n_win_panel)
  M[fill_idx] <- band_vec[good]
  out <- vector("list", n_samples)
  for (s in seq_len(n_samples)) out[[s]] <- M[s, ]
  out
}

# Window-level metadata, in panel-window order
# Pull from panel (de-duplicated to one row per window_idx) for self-consistency
win_meta <- unique(panel_use[, .(window_idx, global_window_id, start_bp, end_bp, pos_mb)])
setkey(win_meta, window_idx)
win_meta <- win_meta[J(win_idx_unique)]

div_roll_lists       <- list()
rank_in_cohort_lists <- list()
rank_band_lists      <- list()

for (sk in scales_use) {
  div_col  <- paste0("div_roll_", sk)
  rank_col <- paste0("rank_in_cohort_", sk)
  band_col <- paste0("rank_band_", sk)

  if (div_col  %in% names(panel)) div_roll_lists[[sk]]       <- build_matrix_list(panel_use[[div_col]])
  if (rank_col %in% names(panel)) rank_in_cohort_lists[[sk]] <- build_matrix_list(panel_use[[rank_col]])
  if (band_col %in% names(panel)) rank_band_lists[[sk]]      <- build_band_matrix_list(panel_use[[band_col]])
}

cat("[GHSL_EXPORT v2]   panel built in ", round((proc.time() - t_panel)[3], 1), "s\n", sep = "")

panel_layer <- list(
  n_samples         = as.integer(n_samples),
  n_windows         = as.integer(n_win_panel),
  samples           = as.character(sample_order),
  window_idx        = as.integer(win_idx_unique),
  global_window_id  = as.integer(win_meta$global_window_id),
  start_bp          = as.integer(win_meta$start_bp),
  end_bp            = as.integer(win_meta$end_bp),
  pos_mb            = round(as.numeric(win_meta$pos_mb), 4),
  scales            = as.character(scales_use),
  primary_scale     = PRIMARY_SCALE,
  div_roll          = div_roll_lists,
  rank_in_cohort    = rank_in_cohort_lists,
  rank_band         = rank_band_lists
)

# =============================================================================
# BUILD: ghsl_kstripes layer
# =============================================================================
# For K = 2..MAX_K, partition samples by their *chromosome-mean* rank at
# the primary scale, then assign each sample a stripe 1..K.
#
# We use chromosome-wide kmeans on the per-sample mean rank. This is the same
# kmeans the v6 classifier does at the candidate level (Part C), but applied
# globally so the stripes are stable across the chromosome and consistent
# with the K-stripe heatmap interpretation. Per-candidate kmeans (which the
# scrubber will eventually want for selected L2s) happens browser-side using
# the dense panel.
#
# Output:
#   ghsl_kstripes: {
#     primary_scale: "s50",
#     samples: ["CGA001", ...],
#     mean_rank_per_sample: [...],          length n_samples (input to kmeans)
#     by_k: {
#       "2": { stripe_per_sample: [1,2,1,...], stripe_means: [0.1, 0.7], stripe_medians: [...], n_per_stripe: [113, 113] },
#       "3": { stripe_per_sample: [...], ...},
#       ...
#       "<MAX_K>": ...
#     }
#   }

cat("[GHSL_EXPORT v2] Computing per-K stripe assignments at ", PRIMARY_SCALE, "...\n", sep = "")

# Mean rank per sample at the primary scale, across all panel windows
prim_rank_col <- paste0("rank_in_cohort_", PRIMARY_SCALE)
if (!prim_rank_col %in% names(panel)) {
  stop("Primary-scale rank column not in panel: ", prim_rank_col)
}

mean_rank_dt <- panel_use[, .(mean_rank = mean(get(prim_rank_col), na.rm = TRUE)),
                          by = sample_id]
# Align to sample_order
setkey(mean_rank_dt, sample_id)
mean_rank_per_sample <- mean_rank_dt[J(sample_order), mean_rank]

# Sanity: replace NaN (no valid ranks) with NA so kmeans skips them
mean_rank_per_sample[is.nan(mean_rank_per_sample)] <- NA_real_

valid_mask <- is.finite(mean_rank_per_sample)
n_valid    <- sum(valid_mask)
cat("[GHSL_EXPORT v2]   ", n_valid, "/", n_samples, " samples have valid mean rank\n", sep = "")

by_k <- list()
for (K in 2:MAX_K) {
  if (n_valid < K) {
    cat("[GHSL_EXPORT v2]   K=", K, ": skipping (not enough valid samples)\n", sep = "")
    next
  }

  vals <- mean_rank_per_sample[valid_mask]
  km <- tryCatch(kmeans(vals, centers = K, nstart = 25, iter.max = 100),
                 error = function(e) NULL)
  if (is.null(km)) {
    cat("[GHSL_EXPORT v2]   K=", K, ": kmeans failed\n", sep = "")
    next
  }

  # Reorder cluster IDs by center value: stripe 1 = lowest mean rank
  center_order <- order(km$centers[, 1])
  remap <- integer(K)
  remap[center_order] <- seq_len(K)
  stripe_valid <- remap[km$cluster]

  # Full per-sample stripe vector (NA for invalid samples)
  stripe_per_sample <- rep(NA_integer_, n_samples)
  stripe_per_sample[valid_mask] <- stripe_valid

  # Stripe summary stats — mean and median of mean_rank per stripe
  stripe_means   <- rep(NA_real_, K)
  stripe_medians <- rep(NA_real_, K)
  n_per_stripe   <- rep(0L, K)
  for (s in seq_len(K)) {
    in_stripe <- which(stripe_per_sample == s)
    if (length(in_stripe) > 0) {
      stripe_means[s]   <- round(mean(mean_rank_per_sample[in_stripe], na.rm = TRUE), 4)
      stripe_medians[s] <- round(median(mean_rank_per_sample[in_stripe], na.rm = TRUE), 4)
      n_per_stripe[s]   <- length(in_stripe)
    }
  }

  by_k[[as.character(K)]] <- list(
    stripe_per_sample = as.integer(stripe_per_sample),
    stripe_means      = stripe_means,
    stripe_medians    = stripe_medians,
    n_per_stripe      = as.integer(n_per_stripe)
  )

  cat("[GHSL_EXPORT v2]   K=", K, ": ",
      paste(sprintf("%d/%.3f", n_per_stripe, stripe_means), collapse = " | "),
      "\n", sep = "")
}

kstripes_layer <- list(
  primary_scale         = PRIMARY_SCALE,
  samples               = as.character(sample_order),
  mean_rank_per_sample  = round(mean_rank_per_sample, 4),
  by_k                  = by_k,
  max_k                 = as.integer(MAX_K)
)

# =============================================================================
# ASSEMBLE JSON
# =============================================================================

layers_present <- c()
out <- list(
  schema_version  = 2,
  chrom           = CHROM,
  n_windows       = as.integer(n_windows),
  n_samples       = as.integer(n_samples),
  `_generated_at` = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  `_generator`    = "export_ghsl_to_json_v2.R"
)

if (length(tracks_list) > 0) {
  out$tracks <- tracks_list
  layers_present <- c(layers_present, "tracks")
}

if (length(karyo_runs_list) > 0) {
  out$ghsl_karyotype_runs <- karyo_runs_list
  layers_present <- c(layers_present, "ghsl_karyotype_runs")
}

# Always emit panel + kstripes (these are the core v2 contribution)
out$ghsl_panel    <- panel_layer
out$ghsl_kstripes <- kstripes_layer
layers_present    <- c(layers_present, "ghsl_panel", "ghsl_kstripes")

out$`_layers_present` <- layers_present

# =============================================================================
# WRITE
# =============================================================================

chrom_dir <- file.path(OUT_DIR, CHROM)
dir.create(chrom_dir, recursive = TRUE, showWarnings = FALSE)
out_json  <- file.path(chrom_dir, paste0(CHROM, "_phase4a_ghsl.json"))

cat("[GHSL_EXPORT v2] Serializing JSON...\n", sep = "")
t_ser <- proc.time()
# auto_unbox=TRUE collapses length-1 vectors to scalars; na="null" emits NAs
# as JSON null which the scrubber treats correctly. digits=NA preserves
# pre-applied rounding. pretty=FALSE because the file is large enough that
# whitespace adds 20-30%.
json_str <- jsonlite::toJSON(out, auto_unbox = TRUE, na = "null",
                             digits = NA, pretty = FALSE)
writeLines(json_str, out_json)
cat("[GHSL_EXPORT v2]   serialized in ", round((proc.time() - t_ser)[3], 1), "s\n", sep = "")

f_size_mb <- round(file.info(out_json)$size / 1e6, 2)
cat("\n[GHSL_EXPORT v2] === DONE ===\n")
cat("[GHSL_EXPORT v2] Output: ", out_json, "\n", sep = "")
cat("[GHSL_EXPORT v2] Size:   ", f_size_mb, " MB\n", sep = "")
cat("[GHSL_EXPORT v2] Layers: ", paste(layers_present, collapse = ", "), "\n", sep = "")
cat("[GHSL_EXPORT v2] schema_version: 2\n")
cat("[GHSL_EXPORT v2] chrom:          ", CHROM, "\n", sep = "")
cat("[GHSL_EXPORT v2] n_windows:      ", n_windows, "\n", sep = "")
cat("[GHSL_EXPORT v2] n_samples:      ", n_samples, "\n", sep = "")
cat("[GHSL_EXPORT v2] scales emitted: ", paste(scales_use, collapse = ", "), "\n", sep = "")
cat("[GHSL_EXPORT v2] max_k:          ", MAX_K, "\n", sep = "")
cat("\n")
cat("[GHSL_EXPORT v2] Drag this file into the scrubber alongside the precomp JSON\n")
cat("[GHSL_EXPORT v2] for ", CHROM, ".\n", sep = "")
cat("[GHSL_EXPORT v2] Layers ghsl_panel and ghsl_kstripes register but do not\n")
cat("[GHSL_EXPORT v2] render yet (rendering arrives in scrubber session 2).\n")
