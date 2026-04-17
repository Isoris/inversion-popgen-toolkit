#!/usr/bin/env Rscript
# =============================================================================
# lib_ghsl_panel.R — GHSL v6 per-sample panel query library
# =============================================================================
# Chat 14 (2026-04-18). Backend for on-demand queries over the dense
# per-sample panel RDS files emitted by
# phase_2_discovery/2e_ghsl/STEP_C04b_snake3_ghsl_classify.R.
#
# The panel is a long-format data.table with one row per (sample, window)
# per chromosome, carrying:
#   - window coordinates: chrom, window_idx, global_window_id, start_bp,
#     end_bp, pos_mb
#   - per-scale divergence and rank columns at every scale s in
#     {s10, s20, s30, s40, s50, s100}:
#       div_roll_<s>, rank_in_cohort_<s>, rank_band_<s>
#   - per-(sample, window) stable-run membership from Part B:
#       in_stable_run (bool), stable_run_call ∈ {INV_INV, INV_nonINV, NA}
#   - window-level (same for all samples at that window): ghsl_v6_score,
#     ghsl_v6_status, rank_stability, div_contrast_z, div_bimodal
#
# This library knows NOTHING about registries. It's a pure data-access
# layer. Registry integration lives in reg$compute$ghsl_* (see
# registries/api/R/registry_loader.R::load_compute_api). Both call into
# this library.
#
# USAGE (script-side):
#   source("utils/lib_ghsl_panel.R")
#   panel_set_default_dir("/path/to/ghsl_v6_out")
#
#   # Dense range read
#   dt <- ghsl_panel_range("C_gar_LG12", 10e6, 14e6)
#
#   # Per-sample aggregate across a range
#   agg <- ghsl_panel_aggregate("C_gar_LG12", 10e6, 14e6,
#                                 summaries = c("mean", "frac_low",
#                                               "longest_low_run_bp"))
#
#   # Multi-subblock scan (e.g. 6 sub-blocks from 4 soft boundaries)
#   subs <- data.table(subblock_id = 1:6,
#                      start_bp = c(10e6, 10.5e6, 11e6, 11.7e6, 12.4e6, 13.2e6),
#                      end_bp   = c(10.5e6, 11e6, 11.7e6, 12.4e6, 13.2e6, 14e6))
#   scan <- ghsl_panel_subblock_scan("C_gar_LG12", subs)
#
#   # Wide matrix for plotting
#   m <- ghsl_panel_wide_matrix("C_gar_LG12", 10e6, 14e6,
#                                metric = "div_roll_s50")
#
# USAGE (registry-side):
#   reg$compute$ghsl_at_candidate(cid = "LG12_17")
#   reg$compute$ghsl_at_subblocks(cid = "LG12_17", subblocks = subs)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# Default directory resolution + per-session in-memory cache
# =============================================================================

# Default GHSL directory. Overridable per call via ghsl_dir arg, or globally
# via panel_set_default_dir(). Also consults GHSL_DIR env var.
.ghsl_panel_env <- new.env(parent = emptyenv())
.ghsl_panel_env$default_dir <- NULL
.ghsl_panel_env$panel_cache <- list()   # chrom -> panel data.table
.ghsl_panel_env$meta_cache  <- list()   # chrom -> metadata attribute

panel_set_default_dir <- function(ghsl_dir) {
  .ghsl_panel_env$default_dir <- ghsl_dir
  invisible(ghsl_dir)
}

panel_get_default_dir <- function() {
  d <- .ghsl_panel_env$default_dir
  if (is.null(d) || !nzchar(d)) d <- Sys.getenv("GHSL_DIR", "")
  if (!nzchar(d)) return(NULL)
  d
}

panel_clear_cache <- function(chrom = NULL) {
  if (is.null(chrom)) {
    .ghsl_panel_env$panel_cache <- list()
    .ghsl_panel_env$meta_cache  <- list()
  } else {
    .ghsl_panel_env$panel_cache[[chrom]] <- NULL
    .ghsl_panel_env$meta_cache[[chrom]]  <- NULL
  }
  invisible(NULL)
}

# =============================================================================
# Resolve a panel RDS path for a chromosome
# =============================================================================
.panel_rds_path <- function(chrom, ghsl_dir) {
  candidates <- c(
    file.path(ghsl_dir, "per_sample", paste0(chrom, ".ghsl_v6.per_sample.rds")),
    file.path(ghsl_dir, paste0(chrom, ".ghsl_v6.per_sample.rds"))
  )
  for (f in candidates) if (file.exists(f)) return(f)
  NULL
}

.annot_rds_path <- function(chrom, ghsl_dir) {
  candidates <- c(
    file.path(ghsl_dir, "annot", paste0(chrom, ".ghsl_v6.annot.rds")),
    file.path(ghsl_dir, paste0(chrom, ".ghsl_v6.annot.rds"))
  )
  for (f in candidates) if (file.exists(f)) return(f)
  NULL
}

.karyo_rds_path <- function(chrom, ghsl_dir) {
  candidates <- c(
    file.path(ghsl_dir, "annot", paste0(chrom, ".ghsl_v6.karyotypes.rds")),
    file.path(ghsl_dir, paste0(chrom, ".ghsl_v6.karyotypes.rds"))
  )
  for (f in candidates) if (file.exists(f)) return(f)
  NULL
}

# =============================================================================
# load_ghsl_panel — memoized per-chrom loader
# =============================================================================
# Returns the per-chrom panel data.table (or NULL if the file is missing).
# Memoizes in-memory so a script looping over candidates on LG12 reads
# the RDS once, not once per candidate.
#
# Side effect: the panel's attr("ghsl_panel_meta", ...) is cached separately
# under .ghsl_panel_env$meta_cache[[chrom]] so downstream can retrieve it
# via ghsl_panel_meta(chrom) without re-reading.
# =============================================================================
load_ghsl_panel <- function(chrom, ghsl_dir = NULL, refresh = FALSE) {
  if (is.null(ghsl_dir)) ghsl_dir <- panel_get_default_dir()
  if (is.null(ghsl_dir)) {
    warning("[ghsl_panel] ghsl_dir not set (pass arg or call panel_set_default_dir)")
    return(NULL)
  }
  if (!refresh && !is.null(.ghsl_panel_env$panel_cache[[chrom]])) {
    return(.ghsl_panel_env$panel_cache[[chrom]])
  }
  f <- .panel_rds_path(chrom, ghsl_dir)
  if (is.null(f)) {
    warning("[ghsl_panel] panel RDS not found for chrom ", chrom,
            " under ", ghsl_dir)
    return(NULL)
  }
  dt <- tryCatch(readRDS(f), error = function(e) {
    warning("[ghsl_panel] readRDS failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(dt)) return(NULL)
  setDT(dt)
  # Set key for fast bp range and sample filters
  setkey(dt, chrom, start_bp)
  .ghsl_panel_env$panel_cache[[chrom]] <- dt
  .ghsl_panel_env$meta_cache[[chrom]]  <- attr(dt, "ghsl_panel_meta")
  dt
}

ghsl_panel_meta <- function(chrom, ghsl_dir = NULL) {
  if (!is.null(.ghsl_panel_env$meta_cache[[chrom]])) {
    return(.ghsl_panel_env$meta_cache[[chrom]])
  }
  # Trigger load if not cached (this populates meta_cache as a side effect)
  load_ghsl_panel(chrom, ghsl_dir)
  .ghsl_panel_env$meta_cache[[chrom]]
}

# =============================================================================
# Scale resolution: normalize user-facing --scale 50 to the column suffix
# =============================================================================
# Accepts: integer 50, "s50", "50". Returns "s50" or stops if the scale
# is not available in the panel metadata for the given chromosome.
# =============================================================================
.resolve_scale <- function(chrom, scale, ghsl_dir = NULL) {
  if (is.null(scale)) {
    meta <- ghsl_panel_meta(chrom, ghsl_dir)
    if (is.null(meta)) return("s50")  # fallback default
    return(meta$primary_scale %||% "s50")
  }
  if (is.numeric(scale)) scale <- paste0("s", as.integer(scale))
  if (!startsWith(scale, "s")) scale <- paste0("s", scale)
  scale
}

# =============================================================================
# ghsl_panel_range — dense range read
# =============================================================================
# Filter to (chrom, [start_bp, end_bp]) and optional sample subset.
# Returns the raw long-format slice with all scale columns present.
#
# Args:
#   chrom         : chromosome name, e.g. "C_gar_LG12"
#   start_bp      : inclusive lower bp (default 0)
#   end_bp        : inclusive upper bp (default Inf)
#   sample_ids    : optional character vector to restrict samples
#   ghsl_dir      : panel dir override
# Returns: data.table or NULL if panel missing.
# =============================================================================
ghsl_panel_range <- function(chrom,
                                start_bp = 0,
                                end_bp   = Inf,
                                sample_ids = NULL,
                                ghsl_dir = NULL) {
  dt <- load_ghsl_panel(chrom, ghsl_dir)
  if (is.null(dt)) return(NULL)
  # Rename args to avoid collision with the data.table's start_bp /
  # end_bp columns inside the NSE filter.
  s0 <- start_bp
  e0 <- end_bp
  out <- dt[start_bp <= e0 & end_bp >= s0]
  if (!is.null(sample_ids)) {
    out <- out[sample_id %in% sample_ids]
  }
  out
}

# =============================================================================
# ghsl_panel_aggregate — per-sample summary across a range
# =============================================================================
# Given (chrom, [start_bp, end_bp]), for each sample in the panel (or
# sample subset) return one row of summary statistics.
#
# Args:
#   chrom, start_bp, end_bp, sample_ids, ghsl_dir — as ghsl_panel_range
#   scale     : which rolling scale to summarize (default: panel primary, s50)
#   summaries : character vector of summary names to include. Available:
#     "mean"              — mean of div_roll_<scale> across the range
#     "median"            — median of div_roll_<scale>
#     "sd"                — sd of div_roll_<scale>
#     "frac_low"          — fraction of windows with rank_band_<scale> == "LOW"
#     "frac_mid"          — same for "MID"
#     "frac_high"         — same for "HIGH"
#     "longest_low_run_bp"  — longest contiguous LOW rank_band run (bp)
#     "longest_high_run_bp" — same for HIGH
#     "rank_mean"         — mean of rank_in_cohort_<scale>
#     "rank_at_peak"      — the sample's rank at the interval's
#                            maximum-contrast window (peak arrangement-
#                            diagnostic window, from div_contrast_z).
#     "stable_run_call"   — majority stable_run_call across range (NA if
#                            no cells are in a stable run)
#     "n_windows"         — total windows overlapping the range
#   ghsl_dir  : override
#
# Returns: data.table with one row per sample, columns = c("sample_id",
#          <summaries>).
# =============================================================================
ghsl_panel_aggregate <- function(chrom,
                                    start_bp = 0,
                                    end_bp   = Inf,
                                    sample_ids = NULL,
                                    scale = NULL,
                                    summaries = c("mean", "frac_low",
                                                  "frac_mid", "frac_high",
                                                  "longest_low_run_bp",
                                                  "longest_high_run_bp",
                                                  "rank_mean", "n_windows"),
                                    ghsl_dir = NULL) {
  sk <- .resolve_scale(chrom, scale, ghsl_dir)
  slice <- ghsl_panel_range(chrom, start_bp, end_bp,
                              sample_ids = sample_ids, ghsl_dir = ghsl_dir)
  if (is.null(slice) || nrow(slice) == 0) return(NULL)

  div_col  <- paste0("div_roll_",       sk)
  rank_col <- paste0("rank_in_cohort_", sk)
  band_col <- paste0("rank_band_",      sk)

  if (!div_col %in% names(slice)) {
    warning("[ghsl_panel] scale '", sk, "' not present in panel; available: ",
            paste(grep("^div_roll_", names(slice), value = TRUE), collapse = ", "))
    return(NULL)
  }

  # Longest-run helper — contiguous same-band stretch converted to bp
  longest_band_run_bp <- function(bands, starts_v, ends_v, target) {
    if (length(bands) == 0L) return(0L)
    is_t <- !is.na(bands) & bands == target
    if (!any(is_t)) return(0L)
    rl <- rle(is_t)
    idx_end <- cumsum(rl$lengths)
    idx_start <- c(1L, head(idx_end, -1) + 1L)
    best <- 0L
    for (ri in seq_along(rl$values)) {
      if (!isTRUE(rl$values[ri])) next
      run_span <- ends_v[idx_end[ri]] - starts_v[idx_start[ri]]
      if (run_span > best) best <- run_span
    }
    as.integer(best)
  }

  # Rank-at-peak: for the enclosing slice, find the window with the max
  # |div_contrast_z| — that's the most inversion-diagnostic window in the
  # range. This peak is computed ONCE across all samples since
  # div_contrast_z is window-level.
  peak_win_idx <- NA_integer_
  if ("div_contrast_z" %in% names(slice) &&
      any(!is.na(slice$div_contrast_z))) {
    # Take one row per window (they all share div_contrast_z)
    w1 <- slice[!duplicated(window_idx), .(window_idx, div_contrast_z)]
    w1 <- w1[!is.na(div_contrast_z)]
    if (nrow(w1) > 0) {
      peak_win_idx <- w1$window_idx[which.max(abs(w1$div_contrast_z))]
    }
  }

  # Aggregate per sample
  out <- slice[, {
    vals <- get(div_col)
    rks  <- if (rank_col %in% names(.SD)) get(rank_col) else rep(NA_real_, .N)
    bds  <- if (band_col %in% names(.SD)) get(band_col) else rep(NA_character_, .N)
    res <- list()
    if ("n_windows"     %in% summaries) res$n_windows     <- .N
    if ("mean"          %in% summaries) res$mean          <- mean(vals, na.rm = TRUE)
    if ("median"        %in% summaries) res$median        <- as.numeric(median(vals, na.rm = TRUE))
    if ("sd"            %in% summaries) res$sd            <- sd(vals, na.rm = TRUE)
    if ("rank_mean"     %in% summaries) res$rank_mean     <- mean(rks, na.rm = TRUE)
    if ("frac_low"      %in% summaries) res$frac_low      <- mean(bds == "LOW",  na.rm = TRUE)
    if ("frac_mid"      %in% summaries) res$frac_mid      <- mean(bds == "MID",  na.rm = TRUE)
    if ("frac_high"     %in% summaries) res$frac_high     <- mean(bds == "HIGH", na.rm = TRUE)
    if ("longest_low_run_bp"  %in% summaries)
      res$longest_low_run_bp  <- longest_band_run_bp(bds, start_bp, end_bp, "LOW")
    if ("longest_high_run_bp" %in% summaries)
      res$longest_high_run_bp <- longest_band_run_bp(bds, start_bp, end_bp, "HIGH")
    if ("rank_at_peak"  %in% summaries) {
      if (!is.na(peak_win_idx)) {
        hit <- which(window_idx == peak_win_idx)
        res$rank_at_peak <- if (length(hit) >= 1) rks[hit[1]] else NA_real_
      } else {
        res$rank_at_peak <- NA_real_
      }
    }
    if ("stable_run_call" %in% summaries) {
      # Majority non-NA stable_run_call (if none, NA)
      sc <- stable_run_call
      sc <- sc[!is.na(sc)]
      if (length(sc) == 0) res$stable_run_call <- NA_character_
      else {
        tab <- sort(table(sc), decreasing = TRUE)
        res$stable_run_call <- names(tab)[1]
      }
    }
    res
  }, by = sample_id]

  out
}

# =============================================================================
# ghsl_panel_subblock_scan — per-sample-per-subblock aggregation
# =============================================================================
# For a chromosome with a list of sub-blocks (e.g. 6 sub-blocks from 4
# soft boundaries inside a complex candidate), run ghsl_panel_aggregate
# on each sub-block and stack the results. Adds optional cross-subblock
# band-transition columns, which are the recombinant-detection signal:
# "sample X is LOW in subblocks 1..3, HIGH in subblocks 4..6."
#
# Args:
#   chrom, sample_ids, scale, summaries, ghsl_dir — as aggregate
#   subblocks : data.table with columns subblock_id, start_bp, end_bp.
#               subblock_id can be integer or character.
#   include_transitions : if TRUE (default), add a per-sample wide column
#               'band_pattern' = concatenated LOW/MID/HIGH labels across
#               sub-blocks in order, plus n_transitions = number of band
#               changes along the sub-block sequence.
# Returns: data.table with columns
#   sample_id, subblock_id, n_windows, mean, ..., and if include_transitions:
#   also an attached attr 'per_sample_patterns' data.table with one row
#   per sample giving band_pattern and n_transitions.
# =============================================================================
ghsl_panel_subblock_scan <- function(chrom,
                                        subblocks,
                                        sample_ids = NULL,
                                        scale = NULL,
                                        summaries = c("mean", "frac_low",
                                                      "frac_mid", "frac_high",
                                                      "rank_mean", "n_windows"),
                                        include_transitions = TRUE,
                                        ghsl_dir = NULL) {
  if (!is.data.table(subblocks)) subblocks <- as.data.table(subblocks)
  req <- c("subblock_id", "start_bp", "end_bp")
  if (!all(req %in% names(subblocks))) {
    stop("[ghsl_panel_subblock_scan] subblocks must have columns: ",
         paste(req, collapse = ", "))
  }
  sk <- .resolve_scale(chrom, scale, ghsl_dir)

  chunks <- vector("list", nrow(subblocks))
  for (i in seq_len(nrow(subblocks))) {
    sb <- subblocks[i]
    agg <- ghsl_panel_aggregate(chrom,
                                  start_bp = sb$start_bp,
                                  end_bp   = sb$end_bp,
                                  sample_ids = sample_ids,
                                  scale = sk,
                                  summaries = summaries,
                                  ghsl_dir = ghsl_dir)
    if (is.null(agg)) next
    agg[, subblock_id := sb$subblock_id]
    agg[, subblock_start_bp := sb$start_bp]
    agg[, subblock_end_bp   := sb$end_bp]
    chunks[[i]] <- agg
  }
  out <- rbindlist(chunks, fill = TRUE, use.names = TRUE)
  if (nrow(out) == 0) return(out)

  # Order columns nicely
  lead_cols <- c("sample_id", "subblock_id", "subblock_start_bp",
                 "subblock_end_bp")
  other_cols <- setdiff(names(out), lead_cols)
  setcolorder(out, c(lead_cols, other_cols))

  if (include_transitions && "frac_low" %in% names(out) &&
      "frac_high" %in% names(out)) {
    # Assign dominant band per (sample, subblock). Majority band wins.
    classify_band <- function(fl, fm, fh) {
      v <- c(LOW = fl %||% 0, MID = fm %||% 0, HIGH = fh %||% 0)
      v[is.na(v)] <- 0
      if (all(v == 0)) return(NA_character_)
      names(v)[which.max(v)]
    }
    out[, dominant_band := mapply(classify_band, frac_low,
                                   if ("frac_mid" %in% names(out)) frac_mid else NA,
                                   frac_high)]
    # Order subblocks by start for each sample
    setorder(out, sample_id, subblock_start_bp)
    patterns <- out[, {
      bp <- paste(dominant_band, collapse = "/")
      ntrans <- sum(
        !is.na(head(dominant_band, -1)) &
        !is.na(tail(dominant_band, -1)) &
        head(dominant_band, -1) != tail(dominant_band, -1)
      )
      list(band_pattern = bp, n_transitions = as.integer(ntrans))
    }, by = sample_id]
    setattr(out, "per_sample_patterns", patterns)
  }

  out
}

# =============================================================================
# ghsl_panel_wide_matrix — pivot to samples × windows matrix for plotting
# =============================================================================
# Args:
#   chrom, start_bp, end_bp, sample_ids, ghsl_dir — as range
#   metric : which panel column to pivot. Default "div_roll_s50".
#            For categorical columns (rank_band_s50) the returned
#            matrix is character-typed.
#   sample_order : optional character vector; if provided, rows are in
#            this order. Samples present in sample_order but not in the
#            slice get NA rows. Samples in the slice but not in
#            sample_order are appended at the end.
# Returns: matrix [samples × windows] with rownames = sample_id,
#          colnames = global_window_id (as character).
# =============================================================================
ghsl_panel_wide_matrix <- function(chrom,
                                      start_bp = 0,
                                      end_bp   = Inf,
                                      sample_ids = NULL,
                                      metric = "div_roll_s50",
                                      sample_order = NULL,
                                      ghsl_dir = NULL) {
  slice <- ghsl_panel_range(chrom, start_bp, end_bp,
                              sample_ids = sample_ids, ghsl_dir = ghsl_dir)
  if (is.null(slice) || nrow(slice) == 0) return(NULL)
  if (!metric %in% names(slice)) {
    warning("[ghsl_panel] metric '", metric, "' not in panel columns")
    return(NULL)
  }
  # Pivot wide via dcast
  form <- as.formula("sample_id ~ global_window_id")
  wide <- dcast(slice, form, value.var = metric)
  # Turn into matrix
  samps <- wide$sample_id
  mat   <- as.matrix(wide[, -"sample_id", with = FALSE])
  rownames(mat) <- samps

  if (!is.null(sample_order)) {
    # Reorder / pad
    in_order <- intersect(sample_order, samps)
    extras   <- setdiff(samps, sample_order)
    pad      <- setdiff(sample_order, samps)
    # Build desired order; NA rows for pad
    if (length(pad) > 0) {
      na_mat <- matrix(NA, nrow = length(pad), ncol = ncol(mat),
                        dimnames = list(pad, colnames(mat)))
      mat <- rbind(mat, na_mat)
    }
    mat <- mat[c(in_order, pad, extras), , drop = FALSE]
  }

  mat
}

# =============================================================================
# Plot helpers (minimal; return base-R plot or save PNG)
# =============================================================================
# These are deliberately minimal — they wrap the data access and produce
# a plot. For publication figures, grab the underlying data.table via
# ghsl_panel_range / ghsl_panel_wide_matrix and compose with figrid or
# ggplot2 separately.
# =============================================================================

# --- Per-sample multi-scale line track ---
plot_ghsl_sample_track <- function(chrom, sample_id,
                                     scales = c("s10", "s50"),
                                     start_bp = 0, end_bp = Inf,
                                     ghsl_dir = NULL,
                                     out_png = NULL,
                                     width = 1200, height = 400) {
  slice <- ghsl_panel_range(chrom, start_bp, end_bp,
                              sample_ids = sample_id, ghsl_dir = ghsl_dir)
  if (is.null(slice) || nrow(slice) == 0) {
    warning("[plot_ghsl_sample_track] no data for ", sample_id)
    return(invisible(NULL))
  }
  x <- slice$pos_mb
  cols <- setNames(grDevices::rainbow(length(scales)), scales)
  if (!is.null(out_png)) {
    grDevices::png(out_png, width = width, height = height, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  ylim <- range(unlist(lapply(scales, function(sk) {
    v <- slice[[paste0("div_roll_", sk)]]
    v[is.finite(v)]
  })), na.rm = TRUE, finite = TRUE)
  plot(x, slice[[paste0("div_roll_", scales[1])]],
       type = "l", col = cols[scales[1]],
       ylim = ylim, xlab = paste0(chrom, " (Mb)"),
       ylab = "Rolling GHSL divergence",
       main = paste0("GHSL — ", sample_id, " on ", chrom))
  for (sk in scales[-1]) {
    lines(x, slice[[paste0("div_roll_", sk)]], col = cols[sk])
  }
  legend("topright", legend = scales, col = cols, lty = 1, bty = "n")
  invisible(slice)
}

# --- Cohort heatmap at chosen scale ---
plot_ghsl_heatmap_track <- function(chrom, start_bp = 0, end_bp = Inf,
                                       metric = "rank_band_s50",
                                       sample_order = NULL,
                                       sample_ids = NULL,
                                       ghsl_dir = NULL,
                                       out_png = NULL,
                                       width = 1400, height = 900) {
  mat <- ghsl_panel_wide_matrix(chrom, start_bp, end_bp,
                                   sample_ids = sample_ids,
                                   metric = metric,
                                   sample_order = sample_order,
                                   ghsl_dir = ghsl_dir)
  if (is.null(mat)) return(invisible(NULL))
  if (!is.null(out_png)) {
    grDevices::png(out_png, width = width, height = height, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  # Categorical metrics: coerce to numeric factor codes
  if (is.character(mat)) {
    # Map LOW/MID/HIGH/NA to 1/2/3/NA for image()
    code_mat <- matrix(NA_integer_, nrow = nrow(mat), ncol = ncol(mat),
                        dimnames = dimnames(mat))
    code_mat[mat == "LOW"]  <- 1L
    code_mat[mat == "MID"]  <- 2L
    code_mat[mat == "HIGH"] <- 3L
    cols <- c("royalblue", "grey80", "firebrick")
    image(x = seq_len(ncol(code_mat)), y = seq_len(nrow(code_mat)),
          z = t(code_mat), col = cols,
          breaks = c(0.5, 1.5, 2.5, 3.5),
          xlab = "window (bp-ordered)", ylab = "sample",
          main = paste0(chrom, " : ", metric, "  (",
                         format(start_bp, big.mark = ","), "..",
                         format(end_bp,   big.mark = ","), " bp)"),
          axes = FALSE)
    legend("topright",
           legend = c("LOW", "MID", "HIGH"),
           fill = cols, bty = "n", cex = 0.8)
  } else {
    vals <- as.numeric(mat)
    rng <- range(vals, na.rm = TRUE, finite = TRUE)
    image(x = seq_len(ncol(mat)), y = seq_len(nrow(mat)),
          z = t(mat),
          col = grDevices::hcl.colors(50, "YlOrRd", rev = FALSE),
          zlim = rng,
          xlab = "window", ylab = "sample",
          main = paste0(chrom, " : ", metric),
          axes = FALSE)
  }
  invisible(mat)
}

# --- Sub-block panel ---
plot_ghsl_subblock_panel <- function(chrom, subblocks,
                                        sample_ids = NULL,
                                        scale = "s20",
                                        ghsl_dir = NULL,
                                        out_png = NULL,
                                        width = 1000, height = 900) {
  scan <- ghsl_panel_subblock_scan(chrom, subblocks,
                                     sample_ids = sample_ids, scale = scale,
                                     ghsl_dir = ghsl_dir)
  if (is.null(scan) || nrow(scan) == 0) return(invisible(NULL))
  if (!"dominant_band" %in% names(scan)) {
    warning("[plot_ghsl_subblock_panel] band classification missing; re-run with include_transitions=TRUE")
    return(invisible(scan))
  }
  # Pivot dominant_band to sample × subblock_id
  form <- as.formula("sample_id ~ subblock_id")
  wide <- dcast(scan, form, value.var = "dominant_band")
  samps <- wide$sample_id
  mat <- as.matrix(wide[, -"sample_id", with = FALSE])
  rownames(mat) <- samps
  code_mat <- matrix(NA_integer_, nrow = nrow(mat), ncol = ncol(mat),
                      dimnames = dimnames(mat))
  code_mat[mat == "LOW"]  <- 1L
  code_mat[mat == "MID"]  <- 2L
  code_mat[mat == "HIGH"] <- 3L
  cols <- c("royalblue", "grey80", "firebrick")
  if (!is.null(out_png)) {
    grDevices::png(out_png, width = width, height = height, res = 120)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  image(x = seq_len(ncol(code_mat)), y = seq_len(nrow(code_mat)),
        z = t(code_mat), col = cols,
        breaks = c(0.5, 1.5, 2.5, 3.5),
        xlab = "sub-block", ylab = "sample",
        main = paste0(chrom, " : sub-block dominant band (", scale, ")"),
        axes = FALSE)
  legend("topright", legend = c("LOW", "MID", "HIGH"),
         fill = cols, bty = "n", cex = 0.8)
  invisible(scan)
}

message("[lib_ghsl_panel] loaded — call panel_set_default_dir() to set path")
