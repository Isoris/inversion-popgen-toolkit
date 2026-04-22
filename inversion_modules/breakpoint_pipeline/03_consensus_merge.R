#!/usr/bin/env Rscript

# =============================================================================
# 03_consensus_merge.R  (v1.1 — registry-wired)
#
# CONSENSUS BREAKPOINT MERGER — third stage.
#
# Loads breakpoint estimates from every available source for each candidate,
# computes weighted-median consensus, weighted-MAD-based CI, and per-method
# agreement flags. Gracefully handles missing sources.
#
# SOURCES and WEIGHTS (see METHODOLOGY.md §3.7):
#   ancestral_fragments  (02 fragment distribution mode)        weight 3.0
#   block_extension      (01 ext_left/ext_right)                weight 2.0
#   c01j_transitions     (regime_transitions.tsv)               weight 2.0
#   step40_coherence     (candidate_internal_breaks.tsv)        weight 1.0
#   step41_switches     (candidate_switching_events.tsv.gz)     weight 1.0
#   c01l_segments       (segment_summary.tsv.gz Delta_12 drops) weight 1.0
#   step37_sv           (breakpoint_support_per_candidate.tsv)  weight 0.5
#
# REGISTRY WIRING (chat-18):
#   - 01's block read via reg$evidence$read_block(cid, "dosage_blocks")
#   - 02's block read via reg$evidence$read_block(cid, "ancestral_fragments_summary")
#   - c01j / step40 / step41 / c01l / step37 still read from disk paths
#     (configured via env vars C01J_DIR / STEP40_DIR / etc.). Any of them
#     missing = that method is silently skipped (degraded consensus).
#   - OUTPUT blocks (do NOT overwrite STEP_C01g's boundary_{side} blocks):
#       boundary_refined_left  — refined LEFT consensus
#       boundary_refined_right — refined RIGHT consensus
#       breakpoints_per_method — long-format per-method audit
#   - Catalog view can be built on demand via reg$evidence$get_keys() over
#     the q3_refined_* keys that the schema's keys_extracted populates.
#
# Usage:
#   Rscript 03_consensus_merge.R [cid=all]
#   Rscript 03_consensus_merge.R LG28_1
# =============================================================================

# --- Entry: one line sources registry + ancestry stack -----------------------
Sys.setenv(CURRENT_SCRIPT = "03_consensus_merge.R")
.bridge <- Sys.getenv("REGISTRY_BRIDGE", "utils/registry_bridge.R")
if (!file.exists(.bridge)) {
  for (p in c("utils/registry_bridge.R",
              "../utils/registry_bridge.R",
              file.path(Sys.getenv("BASE", ""), "utils/registry_bridge.R"))) {
    if (file.exists(p)) { .bridge <- p; break }
  }
}
if (!file.exists(.bridge)) {
  stop("[03_consensus_merge] cannot locate utils/registry_bridge.R. ",
       "Set $REGISTRY_BRIDGE or $BASE and retry.")
}
source(.bridge)

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
config_file <- ""
cid_filter  <- NA_character_
.ai <- 1L
while (.ai <= length(args)) {
  a <- args[.ai]
  if (a == "--config" && .ai < length(args)) {
    config_file <- args[.ai + 1L]; .ai <- .ai + 2L
  } else if (a != "all") {
    cid_filter <- a; .ai <- .ai + 1L
  } else {
    .ai <- .ai + 1L
  }
}
if (nzchar(config_file) && file.exists(config_file)) source(config_file)

if (!exists("ensure_dir"))
  ensure_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); invisible(p) }

# Optional upstream dirs for the non-registry sources (env-var overridable).
# If the dir isn't set or the file isn't there, the loader returns NULL and
# that method is silently skipped.
if (!exists("C01J_DIR"))    C01J_DIR    <- Sys.getenv("C01J_DIR",    NA_character_)
if (!exists("C01L_DIR"))    C01L_DIR    <- Sys.getenv("C01L_DIR",    NA_character_)
if (!exists("STEP37_DIR"))  STEP37_DIR  <- Sys.getenv("STEP37_DIR",  NA_character_)
if (!exists("STEP40_DIR"))  STEP40_DIR  <- Sys.getenv("STEP40_DIR",  NA_character_)
if (!exists("STEP41_DIR"))  STEP41_DIR  <- Sys.getenv("STEP41_DIR",  NA_character_)
for (nm in c("C01J_DIR","C01L_DIR","STEP37_DIR","STEP40_DIR","STEP41_DIR")) {
  v <- get(nm)
  if (is.na(v) || !nzchar(v)) assign(nm, NULL)
}

BP03_PARAMS <- list(
  weights = c(
    ancestral_fragments = 3.0,
    block_extension     = 2.0,
    c01j_transitions    = 2.0,
    step40_coherence    = 1.0,
    step41_switches     = 1.0,
    c01l_segments       = 1.0,
    step37_sv           = 0.5
  ),
  agreement_mad_mult = 2.0,  # a method "agrees" if within 2*wMAD of consensus
  min_sources_required = 1L  # at least ancestral_fragments or block_extension
)
if (exists("BP03_PARAMS_OVERRIDE") && is.list(BP03_PARAMS_OVERRIDE)) {
  for (k in names(BP03_PARAMS_OVERRIDE)) {
    if (k == "weights" && is.list(BP03_PARAMS_OVERRIDE$weights)) {
      for (wk in names(BP03_PARAMS_OVERRIDE$weights))
        BP03_PARAMS$weights[[wk]] <- BP03_PARAMS_OVERRIDE$weights[[wk]]
    } else {
      BP03_PARAMS[[k]] <- BP03_PARAMS_OVERRIDE[[k]]
    }
  }
}

# =============================================================================
# Weighted median & weighted MAD
# =============================================================================
weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(NA_real_)
  if (length(x) == 1) return(x)
  ord <- order(x); x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  # smallest x with cw >= 0.5
  x[which(cw >= 0.5)[1]]
}

weighted_mad <- function(x, w, med = NULL) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(NA_real_)
  if (is.null(med)) med <- weighted_median(x, w)
  dev <- abs(x - med)
  wm <- weighted_median(dev, w)
  if (!is.finite(wm)) NA_real_ else 1.4826 * wm
}

# =============================================================================
# Source loaders — each returns a data.table with at minimum:
#   candidate_id, chrom, side ("left"|"right"), breakpoint_bp,
#   ci_low (may be NA), ci_high (may be NA)
# =============================================================================

load_ancestral_fragments <- function(cid, chr) {
  # Read 02's block from the registry
  blk <- tryCatch(reg$evidence$read_block(cid, "ancestral_fragments_summary"),
                  error = function(e) NULL)
  if (is.null(blk) || is.null(blk$data)) return(NULL)
  d <- blk$data
  if (is.null(d$status) || d$status != "ok") return(NULL)
  rows <- list()
  if (is.finite(d$frag_left_bp_mode %||% NA_real_)) {
    rows[[length(rows) + 1]] <- data.table(
      candidate_id = cid, chrom = chr, side = "left",
      method = "ancestral_fragments",
      breakpoint_bp = as.numeric(d$frag_left_bp_mode),
      ci_low  = as.numeric(d$frag_left_ci_low  %||% NA_real_),
      ci_high = as.numeric(d$frag_left_ci_high %||% NA_real_)
    )
  }
  if (is.finite(d$frag_right_bp_mode %||% NA_real_)) {
    rows[[length(rows) + 1]] <- data.table(
      candidate_id = cid, chrom = chr, side = "right",
      method = "ancestral_fragments",
      breakpoint_bp = as.numeric(d$frag_right_bp_mode),
      ci_low  = as.numeric(d$frag_right_ci_low  %||% NA_real_),
      ci_high = as.numeric(d$frag_right_ci_high %||% NA_real_)
    )
  }
  if (length(rows) == 0) NULL else rbindlist(rows)
}

load_block_extension <- function(cid, chr) {
  blk <- tryCatch(reg$evidence$read_block(cid, "dosage_blocks"),
                  error = function(e) NULL)
  if (is.null(blk) || is.null(blk$data)) return(NULL)
  d <- blk$data
  if (is.null(d$status) || d$status != "ok") return(NULL)
  if (is.null(d$ext_left_bp) || is.null(d$ext_right_bp)) return(NULL)
  rbindlist(list(
    data.table(candidate_id = cid, chrom = chr, side = "left",
                method = "block_extension",
                breakpoint_bp = as.numeric(d$ext_left_bp),
                ci_low = NA_real_, ci_high = NA_real_),
    data.table(candidate_id = cid, chrom = chr, side = "right",
                method = "block_extension",
                breakpoint_bp = as.numeric(d$ext_right_bp),
                ci_low = NA_real_, ci_high = NA_real_)
  ))
}

load_c01j_transitions <- function(cid, chr, dir_path, c_start, c_end) {
  if (is.null(dir_path)) return(NULL)
  # C01j writes regime_transitions.tsv per region / per candidate
  # Look for file patterns that likely match this candidate
  f_candidates <- c(
    file.path(dir_path, paste0("regime_transitions_", chr, "_cid", cid, ".tsv")),
    file.path(dir_path, paste0("regime_transitions_", chr, ".tsv")),
    file.path(dir_path, "regime_transitions.tsv")
  )
  f <- f_candidates[file.exists(f_candidates)][1]
  if (is.na(f)) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  # Filter to this candidate region if possible
  if ("chrom" %in% names(d)) d <- d[chrom == chr]
  # Expect columns start_mb / end_mb or start_bp / end_bp; adapt
  pos_col <- NULL
  if ("start_mb" %in% names(d))      pos_col <- function(x) x * 1e6
  else if ("start_bp" %in% names(d)) pos_col <- function(x) x
  if (is.null(pos_col)) return(NULL)

  # Enter- and exit-structured are breakpoint-like
  if ("transition_type" %in% names(d)) {
    enters <- d[grepl("enter", transition_type)]
    exits  <- d[grepl("exit",  transition_type)]
  } else {
    enters <- d; exits <- d
  }
  # Pick the enter closest to c_start (left breakpoint) and exit closest to c_end (right)
  get_col <- function(dt, col) {
    if (col %in% names(dt)) dt[[col]] else NA_real_
  }
  left_bp <- NA_real_; right_bp <- NA_real_
  if (nrow(enters) > 0) {
    pos_enter <- pos_col(get_col(enters, "start_mb") %||% get_col(enters, "start_bp"))
    left_bp <- pos_enter[which.min(abs(pos_enter - c_start))]
  }
  if (nrow(exits) > 0) {
    pos_exit <- pos_col(get_col(exits, "end_mb") %||% get_col(exits, "end_bp"))
    right_bp <- pos_exit[which.min(abs(pos_exit - c_end))]
  }
  rows <- list()
  if (is.finite(left_bp)) rows[[length(rows)+1]] <- data.table(
    candidate_id = cid, chrom = chr, side = "left",
    method = "c01j_transitions", breakpoint_bp = left_bp,
    ci_low = NA_real_, ci_high = NA_real_)
  if (is.finite(right_bp)) rows[[length(rows)+1]] <- data.table(
    candidate_id = cid, chrom = chr, side = "right",
    method = "c01j_transitions", breakpoint_bp = right_bp,
    ci_low = NA_real_, ci_high = NA_real_)
  if (length(rows) == 0) NULL else rbindlist(rows)
}

load_step40_coherence <- function(cid, chr, dir_path) {
  # STEP40 writes per-candidate internal coherence breaks. If no directory
  # configured (STEP40_DIR env unset), skip silently — it's optional.
  if (is.null(dir_path)) return(NULL)
  f_candidates <- c(
    file.path(dir_path, cid, "candidate_internal_breaks.tsv"),
    file.path(dir_path, paste0(chr, ".candidate_", cid),
              "candidate_internal_breaks.tsv"),
    file.path(dir_path, paste0("candidate_internal_breaks_", cid, ".tsv"))
  )
  f <- f_candidates[file.exists(f_candidates)][1]
  if (is.na(f)) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if (!all(c("break_start_bp", "break_end_bp") %in% names(d))) return(NULL)
  left_bp  <- min(d$break_start_bp, na.rm = TRUE)
  right_bp <- max(d$break_end_bp,   na.rm = TRUE)
  if (!is.finite(left_bp) || !is.finite(right_bp)) return(NULL)
  rbindlist(list(
    data.table(candidate_id = cid, chrom = chr, side = "left",
                method = "step40_coherence", breakpoint_bp = left_bp,
                ci_low = NA_real_, ci_high = NA_real_),
    data.table(candidate_id = cid, chrom = chr, side = "right",
                method = "step40_coherence", breakpoint_bp = right_bp,
                ci_low = NA_real_, ci_high = NA_real_)
  ))
}

load_step41_switches <- function(cid, chr, dir_path, c_start, c_end) {
  if (is.null(dir_path)) return(NULL)
  f_candidates <- c(
    file.path(dir_path, cid, "candidate_switching_events.tsv.gz"),
    file.path(dir_path, cid, "candidate_switching_events.tsv"),
    file.path(dir_path, paste0(chr, ".candidate_", cid),
              "candidate_switching_events.tsv.gz"),
    file.path(dir_path, paste0(chr, ".candidate_", cid),
              "candidate_switching_events.tsv")
  )
  f <- f_candidates[file.exists(f_candidates)][1]
  if (is.na(f)) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !"switch_bp" %in% names(d)) return(NULL)
  if (nrow(d) < 5) return(NULL)
  dens <- tryCatch(density(d$switch_bp, n = 1024), error = function(e) NULL)
  if (is.null(dens)) return(NULL)
  peak_x <- dens$x; peak_y <- dens$y
  is_peak <- peak_y[-c(1, length(peak_y))] > peak_y[-c(length(peak_y)-1, length(peak_y))] &
             peak_y[-c(1, length(peak_y))] > peak_y[-c(1, 2)]
  peak_positions <- peak_x[-c(1, length(peak_x))][is_peak]
  if (length(peak_positions) == 0) return(NULL)
  left_bp  <- peak_positions[which.min(abs(peak_positions - c_start))]
  right_bp <- peak_positions[which.min(abs(peak_positions - c_end))]
  rbindlist(list(
    data.table(candidate_id = cid, chrom = chr, side = "left",
                method = "step41_switches", breakpoint_bp = left_bp,
                ci_low = NA_real_, ci_high = NA_real_),
    data.table(candidate_id = cid, chrom = chr, side = "right",
                method = "step41_switches", breakpoint_bp = right_bp,
                ci_low = NA_real_, ci_high = NA_real_)
  ))
}

load_c01l_segments <- function(cid, chr, dir_path, c_start, c_end) {
  if (is.null(dir_path)) return(NULL)
  f <- file.path(dir_path, "segment_summary.tsv.gz")
  if (!file.exists(f)) f <- file.path(dir_path, "segment_summary.tsv")
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if ("candidate_id" %in% names(d)) d <- d[candidate_id == cid]
  if (nrow(d) == 0) return(NULL)
  # This is an indirect breakpoint signal — the breakpoint is *at* the
  # boundary between the left_flank and inv_left_half segments, and between
  # inv_right_half and right_flank. These boundaries are themselves the
  # segment boundaries (from the candidate c_start/c_end). So c01l's
  # "breakpoint estimate" is essentially the input candidate boundary —
  # but with a confidence derived from whether Delta_12 really drops there.
  # Only emit an estimate if Delta_12 is substantially lower in inv_core
  # than in the flanks (validating the candidate boundaries).
  if (!"delta_12" %in% names(d) && !"Delta_12" %in% names(d)) return(NULL)
  delta_col <- ifelse("delta_12" %in% names(d), "delta_12", "Delta_12")
  if (!"segment" %in% names(d)) return(NULL)
  core <- mean(d[grepl("core|inv_", segment), get(delta_col)], na.rm = TRUE)
  lf <- mean(d[grepl("left_flank",  segment), get(delta_col)], na.rm = TRUE)
  rf <- mean(d[grepl("right_flank", segment), get(delta_col)], na.rm = TRUE)
  # Emit only if the contrast is strong (flank substantially higher delta_12)
  if (is.finite(core) && is.finite(lf) && (lf - core) > 0.2) {
    rows_l <- data.table(candidate_id = cid, chrom = chr, side = "left",
                          method = "c01l_segments", breakpoint_bp = c_start,
                          ci_low = NA_real_, ci_high = NA_real_)
  } else rows_l <- NULL
  if (is.finite(core) && is.finite(rf) && (rf - core) > 0.2) {
    rows_r <- data.table(candidate_id = cid, chrom = chr, side = "right",
                          method = "c01l_segments", breakpoint_bp = c_end,
                          ci_low = NA_real_, ci_high = NA_real_)
  } else rows_r <- NULL
  rbindlist(list(rows_l, rows_r))
}

load_step37_sv <- function(cid, chr, dir_path) {
  if (is.null(dir_path)) return(NULL)
  f_candidates <- c(
    file.path(dir_path, "breakpoint_support_per_candidate.tsv"),
    file.path(dir_path, paste0("bp_clusters_", chr, ".tsv"))
  )
  f <- f_candidates[file.exists(f_candidates)][1]
  if (is.na(f)) return(NULL)
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if ("candidate_id" %in% names(d)) d <- d[candidate_id == cid]
  if (nrow(d) == 0) return(NULL)
  rows <- list()
  # STEP37 writes per-cluster rows with boundary_side and cluster_median
  if ("boundary_side" %in% names(d) && "cluster_median" %in% names(d)) {
    left_meds  <- d[boundary_side == "left",  cluster_median]
    right_meds <- d[boundary_side == "right", cluster_median]
    if (length(left_meds) > 0) rows[[length(rows)+1]] <- data.table(
      candidate_id = cid, chrom = chr, side = "left",
      method = "step37_sv", breakpoint_bp = median(left_meds),
      ci_low = NA_real_, ci_high = NA_real_)
    if (length(right_meds) > 0) rows[[length(rows)+1]] <- data.table(
      candidate_id = cid, chrom = chr, side = "right",
      method = "step37_sv", breakpoint_bp = median(right_meds),
      ci_low = NA_real_, ci_high = NA_real_)
  }
  if (length(rows) == 0) NULL else rbindlist(rows)
}

# =============================================================================
# Per-candidate consensus
# =============================================================================
process_candidate <- function(cid, chr, c_start, c_end) {
  # Optional: read STEP_C01g's boundary blocks for comparison (not overwritten)
  c01g_left <- tryCatch(reg$evidence$read_block(cid, "boundary_left"),
                        error = function(e) NULL)
  c01g_right <- tryCatch(reg$evidence$read_block(cid, "boundary_right"),
                         error = function(e) NULL)
  c01g_left_bp  <- if (!is.null(c01g_left)  && !is.null(c01g_left$data$boundary_bp))
                     as.numeric(c01g_left$data$boundary_bp)  else NA_real_
  c01g_right_bp <- if (!is.null(c01g_right) && !is.null(c01g_right$data$boundary_bp))
                     as.numeric(c01g_right$data$boundary_bp) else NA_real_

  # Collect all sources (registry-backed for 01/02; disk-backed for the rest)
  all_estimates <- list(
    load_ancestral_fragments(cid, chr),
    load_block_extension(cid, chr),
    load_c01j_transitions(cid, chr, C01J_DIR, c_start, c_end),
    load_step40_coherence(cid, chr, STEP40_DIR),
    load_step41_switches(cid, chr, STEP41_DIR, c_start, c_end),
    load_c01l_segments(cid, chr, C01L_DIR, c_start, c_end),
    load_step37_sv(cid, chr, STEP37_DIR)
  )
  all_estimates <- Filter(Negate(is.null), all_estimates)
  if (length(all_estimates) == 0) {
    message("[03_consensus_merge] cid=", cid, " SKIP — no sources available")
    return(invisible(NULL))
  }
  per_method <- rbindlist(all_estimates, fill = TRUE)
  per_method[, weight := BP03_PARAMS$weights[method]]
  per_method[is.na(weight), weight := 0.5]

  # Sanity: drop estimates wildly outside the scan region (10 Mb tolerance)
  per_method <- per_method[breakpoint_bp >= c_start - 10e6 &
                              breakpoint_bp <= c_end + 10e6]
  if (nrow(per_method) == 0) {
    message("[03_consensus_merge] cid=", cid, " SKIP — all sources out of range")
    return(invisible(NULL))
  }

  # Consensus per side
  compute_side <- function(side_label) {
    ps <- per_method[side == side_label]
    if (nrow(ps) < BP03_PARAMS$min_sources_required) return(NULL)
    consensus  <- weighted_median(ps$breakpoint_bp, ps$weight)
    w_mad      <- weighted_mad(ps$breakpoint_bp, ps$weight, med = consensus)
    agreement_band <- if (is.finite(w_mad))
      BP03_PARAMS$agreement_mad_mult * w_mad else NA_real_

    ps[, distance_to_consensus_kb := round((breakpoint_bp - consensus) / 1000, 2)]
    ps[, within_agreement_band := if (is.finite(agreement_band))
      abs(breakpoint_bp - consensus) <= agreement_band else NA]
    n_agree <- if (is.finite(agreement_band))
      sum(ps$within_agreement_band, na.rm = TRUE) else nrow(ps)

    list(
      per_method = ps,
      consensus_bp = consensus,
      ci_low       = if (is.finite(w_mad)) consensus - 1.96 * w_mad / sqrt(nrow(ps))
                     else NA_real_,
      ci_high      = if (is.finite(w_mad)) consensus + 1.96 * w_mad / sqrt(nrow(ps))
                     else NA_real_,
      wmad_kb      = if (is.finite(w_mad)) w_mad / 1000 else NA_real_,
      n_methods    = nrow(ps),
      n_agreeing   = n_agree,
      method_list  = unique(ps$method)
    )
  }

  left  <- compute_side("left")
  right <- compute_side("right")

  if (is.null(left) && is.null(right)) {
    message("[03_consensus_merge] cid=", cid,
            " SKIP — no sides with sufficient sources")
    return(invisible(NULL))
  }

  # Determine primary source per side (highest-weight method that was
  # present AND within the agreement band)
  best_within_band <- function(info) {
    if (is.null(info) || nrow(info$per_method) == 0) return(NA_character_)
    pm <- info$per_method
    pm_in <- pm[isTRUE(within_agreement_band) | is.na(within_agreement_band)]
    if (nrow(pm_in) == 0) pm_in <- pm
    pm_in$method[which.max(pm_in$weight)]
  }
  primary_left  <- best_within_band(left)
  primary_right <- best_within_band(right)

  # --- WRITE: boundary_refined_left block -----------------------------------
  if (!is.null(left)) {
    refined_left_shift_c01g <-
      if (is.finite(c01g_left_bp)) round((left$consensus_bp - c01g_left_bp) / 1000, 2)
      else NA_real_
    reg$evidence$write_block(cid, "boundary_refined_left", list(
      side                    = "left",
      final_bp                = as.integer(round(left$consensus_bp)),
      ci_low                  = if (is.finite(left$ci_low))
                                  as.integer(round(left$ci_low))  else NA_integer_,
      ci_high                 = if (is.finite(left$ci_high))
                                  as.integer(round(left$ci_high)) else NA_integer_,
      ci_width_kb             = if (is.finite(left$ci_low) && is.finite(left$ci_high))
                                  round((left$ci_high - left$ci_low) / 1000, 2)
                                 else NA_real_,
      n_methods_total         = left$n_methods,
      n_methods_agreeing      = left$n_agreeing,
      primary_source          = primary_left,
      sources_available       = as.list(left$method_list),
      weighted_mad_kb         = round(left$wmad_kb, 2),
      refined_vs_input_shift_kb = round((left$consensus_bp - c_start) / 1000, 2),
      refined_vs_c01g_shift_kb  = refined_left_shift_c01g,
      weights_used            = as.list(BP03_PARAMS$weights)
    ))
  }

  # --- WRITE: boundary_refined_right block ----------------------------------
  if (!is.null(right)) {
    refined_right_shift_c01g <-
      if (is.finite(c01g_right_bp)) round((right$consensus_bp - c01g_right_bp) / 1000, 2)
      else NA_real_
    reg$evidence$write_block(cid, "boundary_refined_right", list(
      side                    = "right",
      final_bp                = as.integer(round(right$consensus_bp)),
      ci_low                  = if (is.finite(right$ci_low))
                                  as.integer(round(right$ci_low))  else NA_integer_,
      ci_high                 = if (is.finite(right$ci_high))
                                  as.integer(round(right$ci_high)) else NA_integer_,
      ci_width_kb             = if (is.finite(right$ci_low) && is.finite(right$ci_high))
                                  round((right$ci_high - right$ci_low) / 1000, 2)
                                 else NA_real_,
      n_methods_total         = right$n_methods,
      n_methods_agreeing      = right$n_agreeing,
      primary_source          = primary_right,
      sources_available       = as.list(right$method_list),
      weighted_mad_kb         = round(right$wmad_kb, 2),
      refined_vs_input_shift_kb = round((right$consensus_bp - c_end) / 1000, 2),
      refined_vs_c01g_shift_kb  = refined_right_shift_c01g,
      weights_used            = as.list(BP03_PARAMS$weights)
    ))
  }

  # --- WRITE: long-format per-method block ----------------------------------
  per_method_rows <- lapply(seq_len(nrow(per_method)), function(i) {
    r <- per_method[i]
    list(
      method = r$method,
      side = r$side,
      breakpoint_bp = as.integer(round(r$breakpoint_bp)),
      weight = as.numeric(r$weight),
      ci_low  = if (is.finite(r$ci_low))  as.integer(round(r$ci_low))  else NA_integer_,
      ci_high = if (is.finite(r$ci_high)) as.integer(round(r$ci_high)) else NA_integer_,
      distance_to_consensus_kb = if (!is.null(r$distance_to_consensus_kb))
                                    r$distance_to_consensus_kb else NA_real_,
      within_agreement_band    = if (!is.null(r$within_agreement_band))
                                    r$within_agreement_band else NA
    )
  })
  sides_populated <- character(0)
  if (!is.null(left))  sides_populated <- c(sides_populated, "left")
  if (!is.null(right)) sides_populated <- c(sides_populated, "right")
  reg$evidence$write_block(cid, "breakpoints_per_method", list(
    methods         = per_method_rows,
    n_methods       = length(per_method_rows),
    sides_populated = as.list(sides_populated)
  ))

  message(sprintf(
    "[03_consensus_merge] cid=%s  left=%s (agree=%d/%d)  right=%s (agree=%d/%d)",
    cid,
    if (!is.null(left))  format(round(left$consensus_bp),  big.mark = ",") else "NA",
    if (!is.null(left))  left$n_agreeing  else 0L,
    if (!is.null(left))  left$n_methods   else 0L,
    if (!is.null(right)) format(round(right$consensus_bp), big.mark = ",") else "NA",
    if (!is.null(right)) right$n_agreeing else 0L,
    if (!is.null(right)) right$n_methods  else 0L))

  invisible(list(cid = cid, left = left, right = right))
}

# =============================================================================
# Driver
# =============================================================================
main <- function() {
  cand_path <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                          "interval_registry", "candidate_intervals.tsv")
  if (!file.exists(cand_path)) {
    stop("[03_consensus_merge] candidate_intervals.tsv not found at ", cand_path)
  }
  cand <- fread(cand_path)
  if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]
  if (nrow(cand) == 0) {
    message("[03_consensus_merge] No candidates to process",
            if (!is.na(cid_filter)) paste0(" (filter: ", cid_filter, ")") else "")
    return(invisible(NULL))
  }
  message("[03_consensus_merge] Processing ", nrow(cand),
          " candidate(s) from interval_registry")
  for (ci in seq_len(nrow(cand))) {
    row <- cand[ci]
    tryCatch(
      process_candidate(as.character(row$candidate_id),
                         as.character(row$chrom),
                         as.numeric(row$start_bp),
                         as.numeric(row$end_bp)),
      error = function(e) {
        message("[03_consensus_merge] cid=", row$candidate_id,
                " ERROR: ", conditionMessage(e))
      }
    )
  }
  message("[03_consensus_merge] DONE")
  message("[03_consensus_merge] Catalog view can be built on demand via ",
          "reg$evidence$get_keys() over q3_refined_* keys.")
}

main()
