#!/usr/bin/env Rscript

# =============================================================================
# export_ancestry_to_json_v1.R
#
# Phase 4 "chapter end" emit script for ancestry/admixture data.
#
# Reads outputs of STEP_Q06_ancestry_tracks.sh:
#   - ancestry_window.<CHR>.tsv          (REQUIRED)
#       per-window mean summary: delta12, entropy, ena, maxQ_label,
#       cv_delta12_across_samples
#   - ancestry_sample.<CHR>.tsv.gz       (optional)
#       per-sample x per-window: maxQ, maxQ_label, delta12, delta13,
#       entropy, ena
#
# Optionally also reads upstream Engine-B caches:
#   - <LOCAL_Q_DIR>/<CHR>.local_Q_summary.tsv.gz
#       if it exposes per-K mean Q columns (Q1, Q2, ..., QK or
#       Q1_mean, Q2_mean, ...), packs them into the ancestry_q_means
#       layer for stacked-area Q-along-chromosome plots.
#
# Emits:
#   <out_dir>/<CHR>/<CHR>_phase4_ancestry.json
#
# === GRID-ALIGNMENT CONTRACT ===
# The Q06 ancestry-window grid is independent of the scrubber's local-PCA
# precomp grid (different window sizes, different boundaries). To handle
# this, every entry in the `tracks` layer carries its own `pos_bp` array
# (per-window mid-position in bp, integer). The scrubber's track renderer
# (v3.28+) prefers `track.pos_bp[i]` over `state.data.windows[i].center_mb`
# when present, and falls back to the precomp grid otherwise. This means
# Q06 tracks plot at their true chromosomal positions regardless of how
# they relate to the local-PCA window definition.
#
# JSON layers registered:
#   tracks                  (additive — auto-discovered as tracksdict_<n>
#                            entries in the popstats page chip toolbar)
#                           - ancestry_delta12
#                           - ancestry_entropy
#                           - ancestry_ena
#                           - ancestry_cv_delta12
#   ancestry_window         per-window block (numeric tracks + maxQ_label
#                            string array). Same per-window grid as tracks
#                            but carries the categorical ancestry-label too,
#                            for the regime-class strip on the eventual
#                            ancestry page.
#   ancestry_sample         (optional) dense per-sample × per-window for
#                            maxQ, delta12, delta13, entropy, ena, plus
#                            maxQ_label as a stripe-string array.
#   ancestry_q_means        (optional) per-window mean Q1, Q2, ..., QK
#                            for the stacked-area along-chromosome view.
#
# Usage:
#   Rscript export_ancestry_to_json_v1.R \
#     --window  ancestry_window.C_gar_LG28.tsv \
#     --chrom   C_gar_LG28 \
#     --out_dir inversion_modules/scrubber/data \
#    [--samples ancestry_sample.C_gar_LG28.tsv.gz] \
#    [--local_q_summary unified_ancestry/local_Q/scale_dense/C_gar_LG28.local_Q_summary.tsv.gz]
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
    "Usage: Rscript export_ancestry_to_json_v1.R --window <window.tsv> --chrom <chr> --out_dir <dir>",
    "                                            [--samples <samples.tsv.gz>]",
    "                                            [--local_q_summary <local_Q_summary.tsv.gz>]",
    "",
    "  --window           ancestry_window.<CHR>.tsv from STEP_Q06 (required)",
    "  --samples          ancestry_sample.<CHR>.tsv.gz (optional, large)",
    "  --local_q_summary  upstream local_Q_summary.tsv.gz (optional, for per-K means)",
    "  --chrom            chromosome name (must match TSV chrom column)",
    "  --out_dir          output base dir; writes <out_dir>/<chrom>/<chrom>_phase4_ancestry.json",
    sep = "\n"
  ), "\n")
  quit(status = 1)
}

WINDOW_TSV       <- NULL
SAMPLES_TSV      <- NULL
LOCAL_Q_SUMMARY  <- NULL
CHROM            <- NULL
OUT_DIR          <- NULL

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--window"          && i < length(args)) { WINDOW_TSV      <- args[i + 1]; i <- i + 2L }
  else if (a == "--samples"         && i < length(args)) { SAMPLES_TSV     <- args[i + 1]; i <- i + 2L }
  else if (a == "--local_q_summary" && i < length(args)) { LOCAL_Q_SUMMARY <- args[i + 1]; i <- i + 2L }
  else if (a == "--chrom"           && i < length(args)) { CHROM           <- args[i + 1]; i <- i + 2L }
  else if (a == "--out_dir"         && i < length(args)) { OUT_DIR         <- args[i + 1]; i <- i + 2L }
  else if (a == "-h" || a == "--help")                   { usage() }
  else { i <- i + 1L }
}

if (is.null(WINDOW_TSV) || is.null(CHROM) || is.null(OUT_DIR)) usage()
if (!file.exists(WINDOW_TSV)) stop("Window TSV not found: ", WINDOW_TSV)

cat("[ANC_EXPORT] Loading window TSV: ", WINDOW_TSV, "\n", sep = "")
win_dt <- fread(WINDOW_TSV)
cat("[ANC_EXPORT]   ", nrow(win_dt), " rows, cols: ",
    paste(names(win_dt), collapse = ", "), "\n", sep = "")

# Filter to requested chromosome
if ("chrom" %in% names(win_dt)) {
  win_dt <- win_dt[chrom == CHROM]
  cat("[ANC_EXPORT]   after chrom=", CHROM, ": ", nrow(win_dt), " rows\n", sep = "")
}
if (nrow(win_dt) == 0) stop("No rows for chrom=", CHROM, " in window TSV")

# =============================================================================
# Locate canonical columns (some downstream variants use slightly different names)
# =============================================================================

resolve_col <- function(dt, candidates) {
  hit <- intersect(candidates, names(dt))
  if (length(hit)) hit[1] else NA_character_
}

start_col   <- resolve_col(win_dt, c("window_start_bp", "window_start", "start_bp", "start"))
end_col     <- resolve_col(win_dt, c("window_end_bp",   "window_end",   "end_bp",   "end"))
midmb_col   <- resolve_col(win_dt, c("window_mid_mb",   "mid_mb",       "pos_mb"))
delta12_col <- resolve_col(win_dt, c("delta12", "delta_12"))
ent_col     <- resolve_col(win_dt, c("entropy", "shannon"))
ena_col     <- resolve_col(win_dt, c("ena", "expH", "exp_entropy", "ena_mean"))
maxlab_col  <- resolve_col(win_dt, c("maxQ_label", "max_q_label", "dom_label"))
cv12_col    <- resolve_col(win_dt, c("cv_delta12_across_samples", "cv_delta12"))

if (is.na(start_col) || is.na(end_col)) stop("Could not find window start/end columns")
# pos_mb fallback
if (is.na(midmb_col)) {
  win_dt[, `_pos_mb` := (get(start_col) + get(end_col)) / 2 / 1e6]
  midmb_col <- "_pos_mb"
}
setorderv(win_dt, midmb_col)

n_windows <- nrow(win_dt)
cat("[ANC_EXPORT] Resolved columns: start=", start_col, ", end=", end_col,
    ", mid_mb=", midmb_col, ", delta12=", delta12_col, ", entropy=", ent_col,
    ", ena=", ena_col, ", maxQ_label=", maxlab_col, ", cv12=", cv12_col, "\n", sep = "")

# =============================================================================
# BUILD: tracks contribution (auto-discovered by popstats page)
# =============================================================================
# BUILD: tracks contribution (auto-discovered by popstats page)
# =============================================================================
# Each entry: { min, max, values: [...], pos_bp: [...] }. NaN/NA -> JSON null.
#
# pos_bp is an integer per-window mid-position in bp. The scrubber's track
# renderer prefers track.pos_bp[i] for X-positioning when present and falls
# back to state.data.windows[i].center_mb otherwise. This decouples track
# data from the precomp window grid so Q06's ancestry-window grid (different
# from local-PCA windows) plots at the right chromosomal positions.

# Build the per-window mid-bp vector once (integer bp) and reuse it for every
# track in the tracks_list so they all share the same X-axis grid.
pos_bp_vec <- as.integer(round((win_dt[[start_col]] + win_dt[[end_col]]) / 2))

make_track <- function(values_raw, pos_bp = NULL) {
  v <- suppressWarnings(as.numeric(values_raw))
  finite <- v[is.finite(v)]
  rng <- if (length(finite)) range(finite) else c(NA_real_, NA_real_)
  out <- list(
    min    = if (is.finite(rng[1])) round(rng[1], 6) else NA_real_,
    max    = if (is.finite(rng[2])) round(rng[2], 6) else NA_real_,
    values = round(v, 6)
  )
  if (!is.null(pos_bp)) out$pos_bp <- as.integer(pos_bp)
  out
}

tracks_list <- list()
if (!is.na(delta12_col)) tracks_list[["ancestry_delta12"]]    <- make_track(win_dt[[delta12_col]], pos_bp_vec)
if (!is.na(ent_col))     tracks_list[["ancestry_entropy"]]    <- make_track(win_dt[[ent_col]],     pos_bp_vec)
if (!is.na(ena_col))     tracks_list[["ancestry_ena"]]        <- make_track(win_dt[[ena_col]],     pos_bp_vec)
if (!is.na(cv12_col))    tracks_list[["ancestry_cv_delta12"]] <- make_track(win_dt[[cv12_col]],    pos_bp_vec)

cat("[ANC_EXPORT] Tracks emitted: ", paste(names(tracks_list), collapse = ", "),
    " (", length(pos_bp_vec), " windows, bp range ",
    min(pos_bp_vec, na.rm=TRUE), "–", max(pos_bp_vec, na.rm=TRUE), ")\n", sep = "")

# =============================================================================
# BUILD: ancestry_window layer (rich block — includes maxQ_label categorical)
# =============================================================================
# Same per-window grid as the precomp's windows array, ordered by mid_mb.
# The numeric arrays here mirror the `tracks` entries (slightly redundant but
# convenient: a consumer can read this layer alone without needing to resolve
# the additive tracks merge). The maxQ_label string array is the new bit.

ancestry_window <- list(
  n_windows  = as.integer(n_windows),
  start_bp   = as.integer(win_dt[[start_col]]),
  end_bp     = as.integer(win_dt[[end_col]]),
  pos_mb     = round(as.numeric(win_dt[[midmb_col]]), 4),
  delta12    = if (!is.na(delta12_col)) round(as.numeric(win_dt[[delta12_col]]), 6) else NULL,
  entropy    = if (!is.na(ent_col))     round(as.numeric(win_dt[[ent_col]]),     6) else NULL,
  ena        = if (!is.na(ena_col))     round(as.numeric(win_dt[[ena_col]]),     6) else NULL,
  cv_delta12 = if (!is.na(cv12_col))    round(as.numeric(win_dt[[cv12_col]]),    6) else NULL,
  maxQ_label = if (!is.na(maxlab_col))  as.character(win_dt[[maxlab_col]]) else NULL
)
# Drop NULL entries (jsonlite would emit them as []; we want them absent)
ancestry_window <- ancestry_window[!sapply(ancestry_window, is.null)]

# =============================================================================
# BUILD: ancestry_sample layer (optional, dense)
# =============================================================================

ancestry_sample <- NULL
if (!is.null(SAMPLES_TSV)) {
  if (!file.exists(SAMPLES_TSV)) {
    cat("[ANC_EXPORT] WARNING: --samples path not found: ", SAMPLES_TSV, "\n", sep = "")
  } else {
    cat("[ANC_EXPORT] Loading samples TSV: ", SAMPLES_TSV, "\n", sep = "")
    s_dt <- fread(SAMPLES_TSV)
    cat("[ANC_EXPORT]   raw rows: ", nrow(s_dt), ", cols: ",
        paste(names(s_dt), collapse = ", "), "\n", sep = "")
    if ("chrom" %in% names(s_dt)) s_dt <- s_dt[chrom == CHROM]
    cat("[ANC_EXPORT]   after chrom filter: ", nrow(s_dt), "\n", sep = "")

    # Resolve sample-side columns
    samp_col   <- resolve_col(s_dt, c("sample", "sample_id", "iid", "ind"))
    midbp_col  <- resolve_col(s_dt, c("window_mid_bp", "mid_bp", "window_mid"))
    s_d12_col  <- resolve_col(s_dt, c("delta12", "delta_12"))
    s_d13_col  <- resolve_col(s_dt, c("delta13", "delta_13"))
    s_ent_col  <- resolve_col(s_dt, c("entropy"))
    s_ena_col  <- resolve_col(s_dt, c("ena", "expH"))
    s_maxQ_col <- resolve_col(s_dt, c("maxQ", "max_q"))
    s_lbl_col  <- resolve_col(s_dt, c("maxQ_label", "max_q_label"))

    if (is.na(samp_col) || is.na(midbp_col)) {
      cat("[ANC_EXPORT] WARNING: samples TSV missing sample/mid_bp cols, skipping\n")
    } else {
      # Build sample order and the (sample, window) panel
      sample_order <- sort(unique(as.character(s_dt[[samp_col]])))
      n_samples    <- length(sample_order)
      # Window order: align to ancestry_window's pos_mb grid via window_mid_bp
      # (each row in samples should match exactly one row in window via mid_bp).
      # We build a window-id lookup from the window TSV's start..end.
      win_idx_by_midbp <- setNames(
        seq_len(n_windows),
        sprintf("%d", round(((win_dt[[start_col]] + win_dt[[end_col]]) / 2)))
      )
      s_dt[, `_wid` := win_idx_by_midbp[sprintf("%d", get(midbp_col))]]
      # Some samples files may carry mid_bp at a different rounding; if too
      # many lookups fail, fall back to the nearest-window approach.
      n_unmapped <- sum(is.na(s_dt$`_wid`))
      if (n_unmapped > nrow(s_dt) * 0.01) {
        cat("[ANC_EXPORT]   ", n_unmapped, " samples rows have no exact window match;",
            " falling back to nearest-window\n", sep = "")
        win_mids <- (win_dt[[start_col]] + win_dt[[end_col]]) / 2
        s_dt[, `_wid` := findInterval(get(midbp_col), win_mids,
                                       all.inside = TRUE)]
      }
      s_dt <- s_dt[!is.na(`_wid`) & `_wid` >= 1L & `_wid` <= n_windows]

      build_sample_matrix_numeric <- function(col, default = NA_real_) {
        if (is.na(col)) return(NULL)
        M <- matrix(default, nrow = n_samples, ncol = n_windows)
        rownames(M) <- sample_order
        rows <- match(as.character(s_dt[[samp_col]]), sample_order)
        cols <- s_dt$`_wid`
        good <- !is.na(rows) & !is.na(cols)
        M[cbind(rows[good], cols[good])] <- as.numeric(s_dt[[col]][good])
        # Round to keep JSON compact (4 decimals)
        M <- round(M, 4)
        out <- vector("list", n_samples)
        for (i in seq_len(n_samples)) out[[i]] <- M[i, ]
        out
      }
      build_sample_matrix_chr <- function(col) {
        if (is.na(col)) return(NULL)
        M <- matrix(NA_character_, nrow = n_samples, ncol = n_windows)
        rownames(M) <- sample_order
        rows <- match(as.character(s_dt[[samp_col]]), sample_order)
        cols <- s_dt$`_wid`
        good <- !is.na(rows) & !is.na(cols)
        M[cbind(rows[good], cols[good])] <- as.character(s_dt[[col]][good])
        out <- vector("list", n_samples)
        for (i in seq_len(n_samples)) out[[i]] <- M[i, ]
        out
      }

      ancestry_sample <- list(
        n_samples  = as.integer(n_samples),
        n_windows  = as.integer(n_windows),
        samples    = sample_order,
        pos_mb     = round(as.numeric(win_dt[[midmb_col]]), 4),
        maxQ       = build_sample_matrix_numeric(s_maxQ_col),
        delta12    = build_sample_matrix_numeric(s_d12_col),
        delta13    = build_sample_matrix_numeric(s_d13_col),
        entropy    = build_sample_matrix_numeric(s_ent_col),
        ena        = build_sample_matrix_numeric(s_ena_col),
        maxQ_label = build_sample_matrix_chr(s_lbl_col)
      )
      ancestry_sample <- ancestry_sample[!sapply(ancestry_sample, is.null)]
      cat("[ANC_EXPORT] ancestry_sample: ", n_samples, " samples × ", n_windows,
          " windows\n", sep = "")
    }
  }
}

# =============================================================================
# BUILD: ancestry_q_means (optional — per-window per-K mean Q values)
# =============================================================================
# Reads the upstream local_Q_summary if provided. The expected shape is one
# row per window with columns Q1, Q2, ..., QK or Q1_mean, Q2_mean, ...
# We auto-detect either pattern. K is inferred from how many Q* columns exist.

ancestry_q_means <- NULL
if (!is.null(LOCAL_Q_SUMMARY)) {
  if (!file.exists(LOCAL_Q_SUMMARY)) {
    cat("[ANC_EXPORT] WARNING: --local_q_summary not found: ", LOCAL_Q_SUMMARY, "\n", sep = "")
  } else {
    cat("[ANC_EXPORT] Loading local_Q_summary: ", LOCAL_Q_SUMMARY, "\n", sep = "")
    q_dt <- fread(LOCAL_Q_SUMMARY)
    if ("chrom" %in% names(q_dt)) q_dt <- q_dt[chrom == CHROM]

    # Resolve the per-window mid coordinate
    qmid_col <- resolve_col(q_dt, c("window_mid_bp", "mid_bp", "window_mid"))
    qstart   <- resolve_col(q_dt, c("window_start_bp", "start_bp", "start"))
    qend     <- resolve_col(q_dt, c("window_end_bp",   "end_bp",   "end"))
    if (is.na(qmid_col) && !is.na(qstart) && !is.na(qend)) {
      q_dt[, `_mid` := (get(qstart) + get(qend)) / 2]
      qmid_col <- "_mid"
    }

    # Detect Q columns: Q1..QK or Q1_mean..QK_mean
    q_pattern_simple <- grep("^Q[0-9]+$", names(q_dt), value = TRUE)
    q_pattern_mean   <- grep("^Q[0-9]+_mean$", names(q_dt), value = TRUE)
    if (length(q_pattern_mean) > length(q_pattern_simple)) {
      q_cols <- q_pattern_mean
      k_extract <- function(s) as.integer(sub("^Q([0-9]+)_mean$", "\\1", s))
    } else {
      q_cols <- q_pattern_simple
      k_extract <- function(s) as.integer(sub("^Q([0-9]+)$", "\\1", s))
    }

    if (length(q_cols) >= 2 && !is.na(qmid_col)) {
      q_cols <- q_cols[order(k_extract(q_cols))]
      K <- length(q_cols)
      setorderv(q_dt, qmid_col)
      # Align to the same number of windows as ancestry_window if possible.
      # If the upstream cache has more rows (e.g. multiple scales merged),
      # we trust the chrom filter already applied.
      if (nrow(q_dt) != n_windows) {
        cat("[ANC_EXPORT] NOTE: q_summary rows (", nrow(q_dt),
            ") != ancestry_window rows (", n_windows, ").",
            " The two layers are emitted independently; the page should",
            " join them by pos_mb.\n", sep = "")
      }
      q_means_pos_mb <- round(as.numeric(q_dt[[qmid_col]]) / 1e6, 4)
      q_matrix <- vector("list", K)
      for (j in seq_len(K)) {
        q_matrix[[j]] <- round(as.numeric(q_dt[[q_cols[j]]]), 6)
      }
      names(q_matrix) <- paste0("Q", seq_len(K))

      ancestry_q_means <- list(
        K        = as.integer(K),
        n_windows = as.integer(nrow(q_dt)),
        pos_mb   = q_means_pos_mb,
        q_means  = q_matrix
      )
      cat("[ANC_EXPORT] ancestry_q_means: K=", K, " over ", nrow(q_dt), " windows\n", sep = "")
    } else {
      cat("[ANC_EXPORT] No Q* columns found in local_Q_summary; skipping ancestry_q_means\n")
    }
  }
}

# =============================================================================
# ASSEMBLE JSON
# =============================================================================

layers_present <- c()
out <- list(
  schema_version  = 2,
  chrom           = CHROM,
  n_windows       = as.integer(n_windows),
  `_generated_at` = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  `_generator`    = "export_ancestry_to_json_v1.R"
)

if (length(tracks_list) > 0) {
  out$tracks <- tracks_list
  layers_present <- c(layers_present, "tracks")
}
if (length(ancestry_window) > 0) {
  out$ancestry_window <- ancestry_window
  layers_present <- c(layers_present, "ancestry_window")
}
if (!is.null(ancestry_sample)) {
  out$ancestry_sample <- ancestry_sample
  layers_present <- c(layers_present, "ancestry_sample")
}
if (!is.null(ancestry_q_means)) {
  out$ancestry_q_means <- ancestry_q_means
  layers_present <- c(layers_present, "ancestry_q_means")
}
out$`_layers_present` <- layers_present

# =============================================================================
# WRITE
# =============================================================================

chrom_dir <- file.path(OUT_DIR, CHROM)
dir.create(chrom_dir, recursive = TRUE, showWarnings = FALSE)
out_json  <- file.path(chrom_dir, paste0(CHROM, "_phase4_ancestry.json"))

cat("[ANC_EXPORT] Serializing JSON...\n")
t_ser <- proc.time()
json_str <- jsonlite::toJSON(out, auto_unbox = TRUE, na = "null",
                             digits = NA, pretty = FALSE)
writeLines(json_str, out_json)
cat("[ANC_EXPORT]   serialized in ", round((proc.time() - t_ser)[3], 1), "s\n", sep = "")

f_size_kb <- round(file.info(out_json)$size / 1024, 1)
f_size_mb <- round(f_size_kb / 1024, 2)
cat("\n[ANC_EXPORT] === DONE ===\n")
cat("[ANC_EXPORT] Output: ", out_json, "\n", sep = "")
cat("[ANC_EXPORT] Size:   ",
    if (f_size_mb > 1) paste0(f_size_mb, " MB") else paste0(f_size_kb, " KB"),
    "\n", sep = "")
cat("[ANC_EXPORT] Layers: ", paste(layers_present, collapse = ", "), "\n", sep = "")
cat("[ANC_EXPORT] schema_version: 2\n")
cat("[ANC_EXPORT] chrom:          ", CHROM, "\n", sep = "")
cat("[ANC_EXPORT] n_windows:      ", n_windows, "\n", sep = "")
cat("\n")
cat("[ANC_EXPORT] Drop this file into the scrubber alongside the precomp JSON\n")
cat("[ANC_EXPORT] for ", CHROM, ". Tracks (", length(tracks_list), ") will\n", sep = "")
cat("[ANC_EXPORT] auto-appear as chips on the popstats page (page 5).\n")
cat("[ANC_EXPORT] The ancestry_window / ancestry_sample / ancestry_q_means\n")
cat("[ANC_EXPORT] layers register but do not yet have dedicated rendering;\n")
cat("[ANC_EXPORT] that arrives in the next session's ancestry page.\n")
