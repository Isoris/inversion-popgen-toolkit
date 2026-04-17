#!/usr/bin/env Rscript
# =============================================================================
# lib_ghsl_confirmation.R — Tier 3 of the chat-9 decompose upgrade design
# =============================================================================
# Independent confirmation of Tier-2 dosage changepoint calls using GHSL's
# per-sample karyotype runs.
#
# Chat 14 (2026-04-18): rewired from v5 to v6. Paths now
#   <ghsl_dir>/annot/<chr>.ghsl_v6.karyotypes.rds
#   <ghsl_dir>/annot/<chr>.ghsl_v6.annot.rds
# Column names in the karyotype RDS are unchanged (call, mean_rank,
# sample_id, chrom, window_start/_end, n_windows, start_bp/_end). The
# annot RDS carries v6-native column names (ghsl_v6_score / ghsl_v6_status
# and unprefixed rank_stability, div_contrast_z, div_bimodal).
#
# GHSL karyotype RDS: one row per run:
#   sample_id, chrom, window_start, window_end, n_windows,
#   call ∈ {INV_INV, INV_nonINV}, mean_rank
#
# A sample can have multiple rows for the same chromosome — a run is a
# consecutive stretch of windows where the sample's per-window karyotype
# indicator sits in the low-rank (INV_INV) or high-rank (INV_nonINV)
# quantile.
#
# USAGE (from STEP_C01i_b_multi_recomb.R):
#
#   source("lib_ghsl_confirmation.R")
#   ghsl_status <- ghsl_per_sample_in_interval(
#     karyo_dt        = load_ghsl_karyo(chr, paths$ghsl_dir),
#     interval        = list(chrom = chr, start_bp = s, end_bp = e),
#     annot_dt        = load_ghsl_annot(chr, paths$ghsl_dir),
#     min_run_overlap = 3L                          # require ≥ 3 windows
#   )
#   # ghsl_status is a data.table with columns:
#   #   sample_id, ghsl_call ∈ {INV_INV, INV_nonINV, SPLIT, UNINFORMATIVE},
#   #   n_runs_in_interval, n_inv_inv_windows, n_inv_non_windows
#
# INTERPRETATION (per chat-9 §Tier 3):
#   - INV_INV:       one run covering the whole interval, call=INV_INV.
#                    Confirms Tier-1 HOM_INV or INV side of recombinant.
#   - INV_nonINV:    one run covering the whole interval, call=INV_nonINV.
#                    Confirms Tier-1 HET or (possibly) recombinant-as-het.
#   - SPLIT:         ≥ 2 runs with DIFFERENT calls within the interval.
#                    DIRECT GHSL evidence of within-sample karyotype change →
#                    Tier-3 corroboration of recombinant.
#   - UNINFORMATIVE: 0 runs overlapping (sample's ranks were MID, or windows
#                    below KARYOTYPE_MIN_RUN threshold). Don't use either way.
#
# Chat 14 addition: `ghsl_per_sample_panel_in_interval()` is an
# alternative panel-based query that reads the dense per-sample panel
# RDS emitted by the v6 classifier. It is ADDITIVE to the run-based
# ghsl_per_sample_in_interval() — run-based logic captures the SPLIT
# classification, which collapses under a single interval-mean. Use
# panel-based for interval-level class summary + multi-scale context;
# use run-based for SPLIT detection. Multi_recomb currently uses
# run-based; panel-based is available for C01d Layer C / characterize_q2
# and for ad-hoc analyses.
#
# LIMITATION: for candidates with < KARYOTYPE_MIN_RUN (typically 10) windows
# total, GHSL never produces karyo runs, so every sample is UNINFORMATIVE.
# Tier 3 is then silent and Tier 2 stands alone with MEDIUM confidence.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# Load GHSL karyotype RDS for a chromosome
# =============================================================================
# Expected path (chat 14, v6): <ghsl_dir>/annot/<chrom>.ghsl_v6.karyotypes.rds
# Fallback order: annot/ subdir v6 RDS, direct v6 RDS, snake3v6 TSV,
# v5 paths (for back-compat during transition). Returns NULL if none found.
# =============================================================================
load_ghsl_karyo <- function(chrom, ghsl_dir) {
  candidates <- c(
    # v6 (chat 14 primary)
    file.path(ghsl_dir, "annot", paste0(chrom, ".ghsl_v6.karyotypes.rds")),
    file.path(ghsl_dir, paste0(chrom, ".ghsl_v6.karyotypes.rds")),
    file.path(ghsl_dir, "snake3v6_karyotype_calls.tsv.gz"),
    # v5 back-compat
    file.path(ghsl_dir, "annot", paste0(chrom, ".ghsl_v5.karyotypes.rds")),
    file.path(ghsl_dir, paste0(chrom, ".ghsl_v5.karyotypes.rds")),
    file.path(ghsl_dir, "snake3v5_karyotype_calls.tsv.gz")
  )
  for (f in candidates) {
    if (!file.exists(f)) next
    if (grepl("\\.rds$", f)) {
      return(tryCatch(as.data.table(readRDS(f)), error = function(e) NULL))
    }
    if (grepl("\\.tsv\\.gz$", f)) {
      dt <- tryCatch(fread(f), error = function(e) NULL)
      if (!is.null(dt) && "chrom" %in% names(dt)) {
        return(dt[chrom == (chrom)])
      }
      return(dt)
    }
  }
  NULL
}

# =============================================================================
# Load GHSL per-chromosome annotation RDS to map window indices → bp
# =============================================================================
load_ghsl_annot <- function(chrom, ghsl_dir) {
  candidates <- c(
    # v6 (chat 14 primary)
    file.path(ghsl_dir, "annot", paste0(chrom, ".ghsl_v6.annot.rds")),
    file.path(ghsl_dir, paste0(chrom, ".ghsl_v6.annot.rds")),
    # v5 back-compat
    file.path(ghsl_dir, "annot", paste0(chrom, ".ghsl_v5.annot.rds")),
    file.path(ghsl_dir, paste0(chrom, ".ghsl_v5.annot.rds"))
  )
  for (f in candidates) {
    if (file.exists(f)) return(tryCatch(as.data.table(readRDS(f)),
                                          error = function(e) NULL))
  }
  NULL
}

# =============================================================================
# Map window indices from karyo_dt to bp using annot_dt if needed
# =============================================================================
# karyo_dt has window_start, window_end as window INDICES (1-based).
# annot_dt has global_window_id + start_bp + end_bp per window.
# Returns karyo_dt with added start_bp, end_bp columns.
# =============================================================================
resolve_karyo_bp <- function(karyo_dt, annot_dt) {
  if (nrow(karyo_dt) == 0) return(karyo_dt)
  if (all(c("start_bp", "end_bp") %in% names(karyo_dt))) return(karyo_dt)
  if (is.null(annot_dt) || nrow(annot_dt) == 0) return(karyo_dt)
  # chat-13 Finding BB: annot_dt must have exactly one row per
  # global_window_id. If a re-run has caused two annotation passes to
  # both contribute rows, subsetting by global_window_id returns
  # multiple rows and the bp-lookup silently picks the first one —
  # leading to wrong bp assignments. Assert up front rather than
  # failing silently.
  stopifnot(uniqueN(annot_dt$global_window_id) == nrow(annot_dt))
  # Build window_index → bp mapping from annot_dt
  # Sort annot by global_window_id to get consistent ordering
  annot_sorted <- annot_dt[order(global_window_id)]
  # window_start is 1-based index into per-chromosome window list
  karyo_dt[, start_bp := annot_sorted$start_bp[pmin(window_start, nrow(annot_sorted))]]
  karyo_dt[, end_bp   := annot_sorted$end_bp[  pmin(window_end,   nrow(annot_sorted))]]
  karyo_dt
}

# =============================================================================
# Per-sample GHSL status in a given candidate interval
# =============================================================================
# Returns one row per sample in the cohort with:
#   sample_id, ghsl_call, n_runs_in_interval,
#   n_inv_inv_windows, n_inv_non_windows
# =============================================================================
ghsl_per_sample_in_interval <- function(karyo_dt,
                                           interval,
                                           annot_dt = NULL,
                                           cohort_sample_ids = NULL,
                                           min_run_overlap = 3L) {
  if (is.null(karyo_dt) || nrow(karyo_dt) == 0) {
    return(data.table(
      sample_id          = cohort_sample_ids %||% character(0),
      ghsl_call          = rep("UNINFORMATIVE",
                                length(cohort_sample_ids %||% character(0))),
      n_runs_in_interval = 0L,
      n_inv_inv_windows  = 0L,
      n_inv_non_windows  = 0L
    ))
  }

  kd <- resolve_karyo_bp(copy(karyo_dt), annot_dt)
  if (!all(c("start_bp", "end_bp") %in% names(kd))) {
    warning("[ghsl_confirmation] karyo_dt missing start_bp/end_bp after resolve")
    return(data.table())
  }

  # Restrict to the candidate's chromosome and overlapping runs
  kd <- kd[chrom == interval$chrom]
  overlap_bp <- pmax(0L,
    pmin(kd$end_bp,   interval$end_bp) -
    pmax(kd$start_bp, interval$start_bp)
  )
  # Overlap in "windows" — use median window size as heuristic
  if (!is.null(annot_dt) && nrow(annot_dt) >= 2L) {
    win_size <- as.integer(median(annot_dt$end_bp - annot_dt$start_bp,
                                     na.rm = TRUE))
    if (is.na(win_size) || win_size <= 0) win_size <- 1L
  } else win_size <- 1L
  overlap_windows <- as.integer(pmax(0L, overlap_bp %/% win_size))
  kd[, overlap_w := overlap_windows]
  kd <- kd[overlap_w >= min_run_overlap]

  # Build cohort-wide output
  cohort <- cohort_sample_ids %||% unique(kd$sample_id)
  out <- data.table(
    sample_id          = cohort,
    ghsl_call          = "UNINFORMATIVE",
    n_runs_in_interval = 0L,
    n_inv_inv_windows  = 0L,
    n_inv_non_windows  = 0L
  )

  for (sid in cohort) {
    s_runs <- kd[sample_id == sid]
    n_runs <- nrow(s_runs)
    if (n_runs == 0L) next
    n_inv_inv <- sum(s_runs[call == "INV_INV"]$overlap_w)
    n_inv_non <- sum(s_runs[call == "INV_nonINV"]$overlap_w)
    idx <- which(out$sample_id == sid)
    out[idx, `:=`(
      n_runs_in_interval = n_runs,
      n_inv_inv_windows  = as.integer(n_inv_inv),
      n_inv_non_windows  = as.integer(n_inv_non)
    )]
    # Classify per chat-9 §Tier 3:
    has_inv_inv <- n_inv_inv >= min_run_overlap
    has_inv_non <- n_inv_non >= min_run_overlap
    gcall <- if (has_inv_inv && has_inv_non) "SPLIT"
             else if (has_inv_inv) "INV_INV"
             else if (has_inv_non) "INV_nonINV"
             else "UNINFORMATIVE"
    out[idx, ghsl_call := gcall]
  }
  out
}

# =============================================================================
# Diagnose GHSL resolution for an interval (is Tier 3 applicable?)
# =============================================================================
# Returns one of:
#   "sufficient"    : annot has ≥ min_windows windows overlapping the interval
#   "insufficient"  : fewer windows → GHSL won't produce runs → Tier 3 silent
#   "no_annot"      : annot file missing — we can't tell
#
# Used by the combination rule in multi_recomb to decide whether
# "NOT G" means "Tier 2 stands alone, medium confidence" or
# "Tier 3 disagreed, flag disputed".
# =============================================================================
ghsl_interval_resolution <- function(annot_dt, interval,
                                         min_windows = 10L) {
  if (is.null(annot_dt) || nrow(annot_dt) == 0) return("no_annot")
  if (!all(c("start_bp","end_bp") %in% names(annot_dt))) return("no_annot")
  n_in <- sum(annot_dt$start_bp <= interval$end_bp &
              annot_dt$end_bp   >= interval$start_bp, na.rm = TRUE)
  if (n_in >= min_windows) "sufficient" else "insufficient"
}

# =============================================================================
# Panel-based per-sample in-interval query (chat 14, additive)
# =============================================================================
# Alternative to `ghsl_per_sample_in_interval` that reads the dense
# per-sample panel (long format) emitted by the v6 classifier. Returns
# one row per sample in `cohort_sample_ids` with:
#   sample_id, panel_call ∈ {INV_INV, INV_nonINV, MIXED, UNINFORMATIVE},
#   frac_low, frac_mid, frac_high, mean_div, rank_mean, n_windows,
#   longest_low_run_bp, longest_high_run_bp
#
# panel_call rules (applied on top of panel aggregate):
#   frac_low  ≥ low_thresh   AND frac_high < mix_thresh  → INV_INV
#   frac_high ≥ high_thresh  AND frac_low  < mix_thresh  → INV_nonINV
#   frac_low  ≥ mix_thresh   AND frac_high ≥ mix_thresh  → MIXED
#   otherwise                                            → UNINFORMATIVE
#
# Differences from run-based `ghsl_per_sample_in_interval`:
#   - Works for short candidates (< KARYOTYPE_MIN_RUN windows) because
#     it aggregates per-window rank bands rather than requiring stable
#     run emission from Part B.
#   - MIXED ≈ SPLIT of the run-based API: both flag "sample changes
#     class within the interval" but MIXED is detected by coexistence
#     of LOW and HIGH bands rather than run-boundary transitions.
#   - Uses the panel library's scale system; default scale = panel's
#     primary (s50 by default). Can be overridden via `scale` arg.
#
# Requires utils/lib_ghsl_panel.R to have been sourced; if it hasn't,
# returns NULL with a warning so the caller can fall back to the
# run-based API.
# =============================================================================
ghsl_per_sample_panel_in_interval <- function(chrom, start_bp, end_bp,
                                                 cohort_sample_ids = NULL,
                                                 ghsl_dir = NULL,
                                                 scale = NULL,
                                                 low_thresh  = 0.50,
                                                 high_thresh = 0.50,
                                                 mix_thresh  = 0.20) {
  if (!exists("ghsl_panel_aggregate", mode = "function")) {
    warning("[ghsl_confirmation] utils/lib_ghsl_panel.R not sourced — ",
            "call source('utils/lib_ghsl_panel.R') first")
    return(NULL)
  }
  agg <- ghsl_panel_aggregate(
    chrom    = chrom,
    start_bp = start_bp,
    end_bp   = end_bp,
    sample_ids = cohort_sample_ids,
    scale    = scale,
    summaries = c("mean", "frac_low", "frac_mid", "frac_high",
                  "rank_mean", "n_windows",
                  "longest_low_run_bp", "longest_high_run_bp"),
    ghsl_dir = ghsl_dir
  )
  if (is.null(agg)) return(NULL)

  # Classify into panel_call
  classify_one <- function(fl, fh) {
    if (is.na(fl) && is.na(fh)) return("UNINFORMATIVE")
    fl <- fl %||% 0
    fh <- fh %||% 0
    if (fl >= mix_thresh && fh >= mix_thresh) return("MIXED")
    if (fl >= low_thresh  && fh <  mix_thresh) return("INV_INV")
    if (fh >= high_thresh && fl <  mix_thresh) return("INV_nonINV")
    "UNINFORMATIVE"
  }
  agg[, panel_call := mapply(classify_one, frac_low, frac_high)]

  # Rename generic columns to panel_* for unambiguous downstream joins
  setnames(agg,
           c("mean"),
           c("mean_div"),
           skip_absent = TRUE)

  # If cohort_sample_ids was supplied, ensure every cohort member gets a
  # row even if they were all-NA (panel may have NA cells in this range).
  if (!is.null(cohort_sample_ids)) {
    missing <- setdiff(cohort_sample_ids, agg$sample_id)
    if (length(missing) > 0) {
      pad <- data.table(sample_id = missing, panel_call = "UNINFORMATIVE")
      agg <- rbindlist(list(agg, pad), fill = TRUE, use.names = TRUE)
    }
  }
  agg[]
}

# =============================================================================
# Chat 14: panel-based per-sample query in an interval
# =============================================================================
# Additive helper. Uses the v6 dense per-sample panel RDS via the
# utils/lib_ghsl_panel.R library. Returns per-sample interval summary
# (mean divergence, rank_mean, frac_low, frac_high, longest LOW run bp,
# longest HIGH run bp) at a chosen rolling scale. NOT a replacement for
# ghsl_per_sample_in_interval — that function captures the SPLIT
# classification via per-sample stable runs, which the panel aggregator
# collapses. Use this for C01d Layer C / characterize_q2 / ad-hoc
# interval-level reporting. Use ghsl_per_sample_in_interval for
# multi_recomb's Tier-3 SPLIT detection.
#
# Args:
#   chrom, start_bp, end_bp : interval coordinates
#   sample_ids              : optional character vector to restrict samples
#   ghsl_dir                : GHSL output directory (for the panel RDS)
#   scale                   : rolling scale, default s50 (250 kb)
#   panel_lib_path          : optional explicit path to lib_ghsl_panel.R;
#                              otherwise tries GHSL_PANEL_LIB env var and
#                              common relative locations
# Returns:
#   data.table, one row per sample, or NULL if the panel or library is
#   unavailable.
# =============================================================================
ghsl_per_sample_panel_in_interval <- function(chrom, start_bp, end_bp,
                                                sample_ids = NULL,
                                                ghsl_dir = NULL,
                                                scale = "s50",
                                                panel_lib_path = NULL) {
  if (!exists("ghsl_panel_aggregate", mode = "function")) {
    libf <- panel_lib_path
    if (is.null(libf) || !file.exists(libf)) {
      libf <- Sys.getenv("GHSL_PANEL_LIB", "")
      if (!nzchar(libf) || !file.exists(libf)) {
        for (cand in c(
          file.path(Sys.getenv("CODEBASE", ""), "utils", "lib_ghsl_panel.R"),
          file.path(Sys.getenv("BASE", ""), "inversion_codebase_v8.5",
                     "utils", "lib_ghsl_panel.R"),
          "utils/lib_ghsl_panel.R",
          "../../../utils/lib_ghsl_panel.R"
        )) {
          if (nzchar(cand) && file.exists(cand)) { libf <- cand; break }
        }
      }
    }
    if (is.null(libf) || !nzchar(libf) || !file.exists(libf)) {
      warning("[ghsl_confirmation] lib_ghsl_panel.R not found; ",
              "set GHSL_PANEL_LIB or pass panel_lib_path")
      return(NULL)
    }
    source(libf, local = FALSE)
  }
  tryCatch(
    ghsl_panel_aggregate(chrom, start_bp, end_bp,
                            sample_ids = sample_ids,
                            scale = scale,
                            summaries = c("mean", "rank_mean",
                                          "frac_low", "frac_mid", "frac_high",
                                          "longest_low_run_bp",
                                          "longest_high_run_bp",
                                          "n_windows"),
                            ghsl_dir = ghsl_dir),
    error = function(e) {
      warning("[ghsl_confirmation] panel aggregate failed: ",
              conditionMessage(e))
      NULL
    }
  )
}

