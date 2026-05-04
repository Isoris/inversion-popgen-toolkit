#!/usr/bin/env Rscript
# =============================================================================
# STEP_D17_multipass_L1_only.R
#
# Stripped-down L1-only version of STEP_D17_multipass.R. Performs ONLY
# Level-1 envelope discovery: a chromosome-wide scan with W=800/400 (or
# whatever pass1_W is set to), with NO masking and NO composition gate.
# Output: a catalogue of envelope candidates with status="ENVELOPE".
#
# What is kept (verbatim from the original):
#   - precomp + sim_mat loading
#   - size-adaptive threshold function
#   - scan_pass (with no masking)
#   - within-scale collapse + size-biased NMS
#   - optional core trimming (trim_to_core) before envelope finalization
#   - optional post-trim NMS to merge sibling envelopes that became
#     distinct only after trimming
#
# What is removed:
#   - Level 2 / Level 3 scans (no looking inside envelopes)
#   - Ward.D2 internal clustering, FRAGMENTED salvage
#   - Adaptive boundary refiner
#   - Z-validation of splits
#   - Composition consistency gate (envelopes are allowed to be
#     heterogeneous — that is the whole point of L1)
#   - Horizontal vs nested architecture branch (only one path)
#
# Output:
#   <chr>_d17L1_envelopes.tsv       : envelope catalogue, all status=ENVELOPE
#
# Author: Claude (Anthropic)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ---- CLI --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

precomp_f      <- get_arg("--precomp")
sim_mat_f      <- get_arg("--sim_mat")
chr_label      <- get_arg("--chr", "chr")
outdir         <- get_arg("--outdir", ".")

# L1 W sizes — same as pass1_W in the full multipass
pass1_W_str    <- get_arg("--pass1_W", "800,400")
pass1_W        <- as.integer(strsplit(pass1_W_str, ",")[[1]])

# Adaptive sim_mat threshold parameters (size-dependent)
sim_thresh_min <- as.numeric(get_arg("--sim_thresh_min", "0.50"))
sim_thresh_max <- as.numeric(get_arg("--sim_thresh_max", "0.70"))
sim_thresh_tau <- as.numeric(get_arg("--sim_thresh_tau", "200"))

# L1 density floor — relaxed (0.45) by default to allow heterogeneous
# envelopes to qualify. Lower this if envelopes are missing.
l1_density_min <- as.numeric(get_arg("--l1_density_min", "0.45"))

# NMS overlap fraction (within-scale collapse + size-biased across W)
nms_overlap    <- as.numeric(get_arg("--nms_overlap", "0.50"))

# Diagonal-scan step size (windows)
scan_step      <- as.integer(get_arg("--scan_step", "5"))

# Min envelope size that survives final filter
min_n_windows  <- as.integer(get_arg("--min_n_windows", "20"))

# Optional core trimming after NMS — tightens broad envelope edges.
trim_core      <- !is.na(get_arg("--trim_core", NA_character_))
trim_drop_z    <- as.numeric(get_arg("--trim_drop_z", "1.0"))

# Optional second NMS after trim, with a tighter overlap threshold.
post_trim_nms          <- !is.na(get_arg("--post_trim_nms", NA_character_))
post_trim_nms_overlap  <- as.numeric(get_arg("--post_trim_nms_overlap", "0.30"))

# Optional diagonal boundary scan. Slides a small WxW upper-triangle cross-
# block at offset G along the diagonal. At each position i the boundary
# score is -median(Z) over the cross-block; peaks of this 1D signal mark
# positions where adjacent regions are unusually separated (i.e. the
# upper-triangle cross-block is unusually blue). Peaks are output to a
# sidecar TSV and rendered as small red squares in the overlay plot.
boundary_scan          <- !is.na(get_arg("--boundary_scan", NA_character_))
boundary_W             <- as.integer(get_arg("--boundary_W", "5"))
boundary_offset        <- as.integer(get_arg("--boundary_offset", "5"))
boundary_score_min     <- as.numeric(get_arg("--boundary_score_min", "0.7"))
boundary_min_dist      <- as.integer(get_arg("--boundary_min_dist", "10"))

# Multi-W validation of each detected peak. From each peak position i,
# scan two 1D rays of diag-Z values perpendicular to the diagonal:
#   right ray: Z[i, i+d]   for d = 1, 2, ..., d_max
#   left  ray: Z[i-d, i]   for d = 1, 2, ..., d_max
# A REAL boundary stays mostly blue (Z < 0) on both rays — the matrix
# never "rebounds" to red as we look further from the diagonal. A FAKE
# boundary (a small blue cross inside an inversion) shows blue near the
# diagonal but jumps red abruptly as we exit the cross-shape into the
# surrounding red bulk.
boundary_validate      <- get_arg("--boundary_validate", "TRUE")
boundary_validate      <- !(toupper(boundary_validate) %in%
                            c("FALSE", "F", "0", "NO"))
# Max ray length (windows)
boundary_perp_d_max    <- as.integer(get_arg("--boundary_perp_d_max", "20"))
# A boundary is STABLE_BLUE if BOTH rays have:
#   (a) frac_blue >= boundary_perp_min_blue_frac      (mostly blue along ray)
#   (b) max_z     <= boundary_perp_max_red            (no big red excursion)
# It's DECAYS if either ray crosses red abruptly:
#   first_d_red  <= boundary_perp_first_d_red_max
boundary_perp_min_blue_frac    <- as.numeric(get_arg(
  "--boundary_perp_min_blue_frac", "0.70"))
boundary_perp_max_red          <- as.numeric(get_arg(
  "--boundary_perp_max_red", "0.50"))
boundary_perp_first_d_red_max  <- as.integer(get_arg(
  "--boundary_perp_first_d_red_max", "5"))
# z threshold defining "red" per pixel
boundary_perp_red_z            <- as.numeric(get_arg(
  "--boundary_perp_red_z", "0.50"))

dry_run        <- !is.na(get_arg("--dry_run", NA_character_))

if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("[D17L1] --precomp is required and must exist")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[D17L1] precomp:    ", precomp_f, "\n")
cat("[D17L1] sim_mat:    ", sim_mat_f %||% "(auto)", "\n")
cat("[D17L1] chr:        ", chr_label, "\n")
cat("[D17L1] outdir:     ", outdir, "\n")
cat("[D17L1] pass1_W:    ", paste(pass1_W, collapse = ","), "\n")
cat("[D17L1] sim thresholds: min=", sim_thresh_min, " max=", sim_thresh_max,
    " tau=", sim_thresh_tau, "  l1_density_min=", l1_density_min, "\n", sep = "")
cat("[D17L1] nms_overlap=", nms_overlap, "  scan_step=", scan_step,
    "  min_n_windows=", min_n_windows, "\n", sep = "")
cat("[D17L1] trim_core=", trim_core, "  trim_drop_z=", trim_drop_z, "\n", sep = "")
cat("[D17L1] post_trim_nms=", post_trim_nms,
    "  post_trim_nms_overlap=", post_trim_nms_overlap, "\n", sep = "")
cat("[D17L1] boundary_scan=", boundary_scan, sep = "")
if (boundary_scan) {
  cat("  W=", boundary_W,
      "  offset=", boundary_offset,
      "  score_min=", boundary_score_min,
      "  min_dist=", boundary_min_dist, sep = "")
}
cat("\n")
if (boundary_scan && boundary_validate) {
  cat("[D17L1] boundary_validate=TRUE  (perp-ray mode)",
      "  d_max=", boundary_perp_d_max,
      "  min_blue_frac=", boundary_perp_min_blue_frac,
      "  max_red=", boundary_perp_max_red,
      "  first_d_red_max=", boundary_perp_first_d_red_max,
      "  red_z=", boundary_perp_red_z, "\n", sep = "")
}

# ---- Adaptive size-dependent threshold function ----------------------------

size_threshold <- function(W) {
  sim_thresh_min + (sim_thresh_max - sim_thresh_min) * exp(-W / sim_thresh_tau)
}

cat("[D17L1] adaptive thresholds by size:\n")
for (W in sort(unique(pass1_W), decreasing = TRUE)) {
  cat(sprintf("[D17L1]   W=%4d  threshold=%.3f\n", W, size_threshold(W)))
}

# ---- Load precomp -----------------------------------------------------------

cat("[D17L1] loading precomp\n")
pc_obj <- readRDS(precomp_f)
if (!is.null(pc_obj$dt)) {
  pc <- pc_obj
} else if (!is.null(pc_obj$pc)) {
  pc <- pc_obj$pc
} else {
  stop("[D17L1] cannot find pc$dt in precomp file")
}

dt_pc <- as.data.table(pc$dt)
n_windows_total <- nrow(dt_pc)

# ---- Normalize window coordinate columns -----------------------------------
if (all(c("start", "end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start
  window_end_bp   <- dt_pc$end
} else if (all(c("start_bp", "end_bp") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start_bp
  window_end_bp   <- dt_pc$end_bp
} else if (all(c("window_start", "window_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$window_start
  window_end_bp   <- dt_pc$window_end
} else if (all(c("bp_start", "bp_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$bp_start
  window_end_bp   <- dt_pc$bp_end
} else {
  stop("[D17L1] precomp dt missing usable coordinate columns. Found: ",
       paste(names(dt_pc), collapse = ", "))
}

cat("[D17L1] N windows: ", n_windows_total, "\n", sep = "")

# ---- Load sim_mat -----------------------------------------------------------

if (is.na(sim_mat_f)) {
  candidates_path <- c(
    file.path(dirname(precomp_f), "sim_mat_nn80.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn160.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn40.rds")
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0)
    stop("[D17L1] --sim_mat not given and could not auto-find sim_mat_nn*.rds")
  sim_mat_f <- hit[1]
  cat("[D17L1] auto-found sim_mat: ", sim_mat_f, "\n")
}

cat("[D17L1] loading sim_mat\n")
sm_obj <- readRDS(sim_mat_f)
if (is.matrix(sm_obj)) {
  sim_mat <- sm_obj
} else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat) && is.matrix(sm_obj$sim_mat)) {
  sim_mat <- sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) == 1 && is.matrix(sm_obj[[1]])) {
  sim_mat <- sm_obj[[1]]
} else {
  stop("[D17L1] sim_mat object structure not recognized; class=", class(sm_obj)[1])
}
storage.mode(sim_mat) <- "double"

if (!isTRUE(nrow(sim_mat) == n_windows_total &&
            ncol(sim_mat) == n_windows_total)) {
  stop("[D17L1] sim_mat dim ", nrow(sim_mat), "x", ncol(sim_mat),
       " does not match precomp N windows ", n_windows_total)
}

# ---- Optional candidate edge trimming --------------------------------------
trim_to_core <- function(s, e, drop_z = 1.0) {
  block <- sim_mat[s:e, s:e]
  row_means <- rowMeans(block, na.rm = TRUE)
  mu <- mean(row_means, na.rm = TRUE)
  sg <- sd(row_means, na.rm = TRUE)
  if (!is.finite(mu) || !is.finite(sg) || sg < 1e-9) return(c(s, e))
  thr <- mu - drop_z * sg
  keep <- which(row_means >= thr)
  if (length(keep) < min_n_windows) return(c(s, e))
  c(s + min(keep) - 1L, s + max(keep) - 1L)
}

apply_core_trimming <- function(dt) {
  if (nrow(dt) == 0L || !trim_core) return(dt)
  old_start <- dt$start_w
  old_end   <- dt$end_w
  for (i in seq_len(nrow(dt))) {
    tr <- trim_to_core(dt$start_w[i], dt$end_w[i], drop_z = trim_drop_z)
    dt$start_w[i] <- tr[1]
    dt$end_w[i]   <- tr[2]
  }
  dt[, trim_start_w_original := old_start]
  dt[, trim_end_w_original   := old_end]
  dt[, trim_left_dropped     := start_w - trim_start_w_original]
  dt[, trim_right_dropped    := trim_end_w_original - end_w]
  dt[, trim_total_dropped    := trim_left_dropped + trim_right_dropped]
  dt[, n_windows             := end_w - start_w + 1L]
  dt[, W                     := n_windows]
  for (i in seq_len(nrow(dt))) {
    block <- sim_mat[dt$start_w[i]:dt$end_w[i], dt$start_w[i]:dt$end_w[i]]
    dt$mean_sim[i] <- mean(block, na.rm = TRUE)
    dt$density_p70[i] <- mean(block >= 0.70, na.rm = TRUE)
    dt$density_adaptive[i] <- mean(block >= dt$threshold[i], na.rm = TRUE)
  }
  dt
}

# ---- Single-pass scanner (no masking; density floor uses l1_density_min) ----
scan_pass_l1 <- function(W_sizes, pass_index = 1L) {
  all_hits <- list()
  for (W in W_sizes) {
    thr_W <- size_threshold(W)
    half <- W %/% 2L
    centers <- seq(half + 1L, n_windows_total - half, by = scan_step)
    if (length(centers) == 0L) next

    out <- vector("list", length(centers))
    out_i <- 1L
    for (ctr in centers) {
      s <- ctr - half
      e <- s + W - 1L
      if (s < 1L || e > n_windows_total) next

      block <- sim_mat[s:e, s:e]
      msim <- mean(block, na.rm = TRUE)
      if (!is.finite(msim) || msim < thr_W) next

      dens_adaptive <- mean(block >= thr_W, na.rm = TRUE)
      if (!is.finite(dens_adaptive) || dens_adaptive < l1_density_min) next
      dens_p70 <- mean(block >= 0.70, na.rm = TRUE)

      out[[out_i]] <- data.table(
        W = W, start_w = s, end_w = e, n_windows = W,
        mean_sim = msim, density_adaptive = dens_adaptive,
        density_p70 = dens_p70, threshold = thr_W,
        discovery_pass = pass_index
      )
      out_i <- out_i + 1L
    }
    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) > 0L) all_hits[[as.character(W)]] <- rbindlist(out)
  }

  if (length(all_hits) == 0L) {
    return(data.table(
      W = integer(0), start_w = integer(0), end_w = integer(0),
      n_windows = integer(0), mean_sim = numeric(0),
      density_adaptive = numeric(0), density_p70 = numeric(0),
      threshold = numeric(0), discovery_pass = integer(0)
    ))
  }
  rbindlist(all_hits)
}

# ---- Within-scale collapse (best-first, NMS-style) -------------------------
collapse_within_scale <- function(dt) {
  if (nrow(dt) == 0L) return(dt)
  dt <- dt[order(-mean_sim)]
  kept <- list()
  for (i in seq_len(nrow(dt))) {
    row_i <- dt[i]
    is_dom <- FALSE
    if (length(kept) > 0L) {
      for (k in kept) {
        o <- max(0L, min(k$end_w, row_i$end_w) -
                 max(k$start_w, row_i$start_w) + 1L)
        smaller <- min(k$end_w - k$start_w + 1L,
                       row_i$end_w - row_i$start_w + 1L)
        if (o / smaller > nms_overlap) { is_dom <- TRUE; break }
      }
    }
    if (!is_dom) kept[[length(kept) + 1L]] <- row_i
  }
  rbindlist(kept)[order(start_w)]
}

# ---- Size-biased NMS across W (bigger wins over smaller) -------------------
size_biased_nms <- function(dt) {
  if (nrow(dt) == 0L) return(dt)
  dt <- dt[order(-W, -mean_sim)]
  kept <- list()
  for (i in seq_len(nrow(dt))) {
    row_i <- dt[i]
    is_dominated <- FALSE
    if (length(kept) > 0L) {
      for (k in kept) {
        o <- max(0L, min(k$end_w, row_i$end_w) -
                 max(k$start_w, row_i$start_w) + 1L)
        smaller <- min(k$end_w - k$start_w + 1L,
                       row_i$end_w - row_i$start_w + 1L)
        if (o / smaller > nms_overlap) { is_dominated <- TRUE; break }
      }
    }
    if (!is_dominated) kept[[length(kept) + 1L]] <- row_i
  }
  rbindlist(kept)[order(start_w)]
}

# ---- Parametric NMS used after trim (uses post-trim widths) ---------------
size_biased_nms_param <- function(dt, overlap_thr) {
  if (nrow(dt) == 0L) return(dt)
  dt <- dt[order(-mean_sim)]
  kept <- list()
  for (i in seq_len(nrow(dt))) {
    row_i <- dt[i]
    is_dominated <- FALSE
    if (length(kept) > 0L) {
      for (k in kept) {
        o <- max(0L, min(k$end_w, row_i$end_w) -
                 max(k$start_w, row_i$start_w) + 1L)
        smaller <- min(k$end_w - k$start_w + 1L,
                       row_i$end_w - row_i$start_w + 1L)
        if (smaller > 0L && o / smaller > overlap_thr) {
          is_dominated <- TRUE
          break
        }
      }
    }
    if (!is_dominated) kept[[length(kept) + 1L]] <- row_i
  }
  rbindlist(kept)[order(start_w)]
}

# ---- Run L1 -----------------------------------------------------------------

cat("\n[D17L1] === LEVEL 1 (envelopes) | W sizes: ",
    paste(pass1_W, collapse = ","), " ===\n", sep = "")

hits <- scan_pass_l1(pass1_W, 1L)
cat("[D17L1]   scan hits: ", nrow(hits), "\n", sep = "")
if (nrow(hits) == 0L) {
  cat("[D17L1] no scan hits — writing empty catalogue\n")
  empty_dt <- data.table(
    chr = character(0), candidate_id = character(0),
    start_w = integer(0), end_w = integer(0),
    start_bp = integer(0), end_bp = integer(0),
    n_windows = integer(0), mean_sim = numeric(0),
    status = character(0)
  )
  out_main <- file.path(outdir, paste0(chr_label, "_d17L1_envelopes.tsv"))
  fwrite(empty_dt, out_main, sep = "\t")
  cat("[D17L1] wrote empty catalogue: ", out_main, "\n", sep = "")
  quit(save = "no", status = 0)
}
cat("[D17L1]   by W: "); print(table(hits$W))

# Within-scale collapse (per W, then concatenate)
by_W <- split(hits, hits$W)
collapsed <- rbindlist(lapply(by_W, collapse_within_scale))
cat("[D17L1]   after within-scale collapse: ", nrow(collapsed), "\n", sep = "")

# Size-biased NMS across W
nms_dt <- size_biased_nms(collapsed)
cat("[D17L1]   after size-biased NMS:        ", nrow(nms_dt), "\n", sep = "")

if (nrow(nms_dt) == 0L) {
  cat("[D17L1] no envelopes survived NMS — writing empty catalogue\n")
  empty_dt <- data.table(
    chr = character(0), candidate_id = character(0),
    start_w = integer(0), end_w = integer(0),
    start_bp = integer(0), end_bp = integer(0),
    n_windows = integer(0), mean_sim = numeric(0),
    status = character(0)
  )
  out_main <- file.path(outdir, paste0(chr_label, "_d17L1_envelopes.tsv"))
  fwrite(empty_dt, out_main, sep = "\t")
  cat("[D17L1] wrote empty catalogue: ", out_main, "\n", sep = "")
  quit(save = "no", status = 0)
}

# Optional core trimming
if (trim_core) {
  nms_dt <- apply_core_trimming(nms_dt)
  cat("[D17L1]   after trim_to_core: median dropped windows = ",
      median(nms_dt$trim_total_dropped, na.rm = TRUE),
      " | max dropped = ",
      max(nms_dt$trim_total_dropped, na.rm = TRUE), "\n", sep = "")
}

# Optional second NMS after trim
if (post_trim_nms && nrow(nms_dt) > 1L) {
  n_before <- nrow(nms_dt)
  nms_dt <- size_biased_nms_param(nms_dt, post_trim_nms_overlap)
  cat("[D17L1]   after post-trim NMS (overlap>", post_trim_nms_overlap,
      "): ", n_before, " -> ", nrow(nms_dt), "\n", sep = "")
}

# Drop tiny envelopes
nms_dt <- nms_dt[n_windows >= min_n_windows]
cat("[D17L1]   final envelope count (n_windows >= ", min_n_windows, "): ",
    nrow(nms_dt), "\n", sep = "")

# Decorate with chr / bp / id columns and status
nms_dt <- nms_dt[order(start_w)]
nms_dt[, chr := chr_label]
nms_dt[, candidate_id := sprintf("%s_d17L1_%04d", chr_label, seq_len(.N))]
nms_dt[, start_bp := window_start_bp[start_w]]
nms_dt[, end_bp   := window_end_bp[end_w]]
nms_dt[, status := "ENVELOPE"]
setnames(nms_dt, "W", "scale_W")

# ---- Console summary --------------------------------------------------------
cat("\n[D17L1] === ENVELOPE SUMMARY for ", chr_label, " ===\n", sep = "")
cat("[D17L1] total envelopes: ", nrow(nms_dt), "\n", sep = "")
if (nrow(nms_dt) > 0L) {
  cat("[D17L1] total windows covered: ",
      sum(nms_dt$n_windows), " / ", n_windows_total,
      " (", round(100 * sum(nms_dt$n_windows) / n_windows_total, 1L), "%)\n",
      sep = "")
  cat("[D17L1] envelopes by scale_W:\n"); print(nms_dt[, .N, by = scale_W])
  cat("[D17L1] all envelopes:\n")
  for (i in seq_len(nrow(nms_dt))) {
    cat(sprintf("[D17L1]   %s  W=%d  %.2f-%.2f Mb  nW=%d  mean_sim=%.3f\n",
                nms_dt$candidate_id[i],
                nms_dt$scale_W[i],
                nms_dt$start_bp[i] / 1e6,
                nms_dt$end_bp[i] / 1e6,
                nms_dt$n_windows[i],
                nms_dt$mean_sim[i]))
  }
}

# ---- Write -----------------------------------------------------------------
out_cols <- c("chr", "candidate_id", "start_w", "end_w",
              "start_bp", "end_bp", "n_windows", "scale_W",
              "mean_sim", "density_p70", "density_adaptive",
              "threshold",
              "trim_start_w_original", "trim_end_w_original",
              "trim_left_dropped", "trim_right_dropped", "trim_total_dropped",
              "status")
for (cc in out_cols) {
  if (!cc %in% names(nms_dt)) nms_dt[, (cc) := NA]
}
out_dt <- nms_dt[, ..out_cols]

if (!dry_run) {
  out_main <- file.path(outdir, paste0(chr_label, "_d17L1_envelopes.tsv"))
  fwrite(out_dt, out_main, sep = "\t")
  cat("\n[D17L1] envelopes written: ", out_main, "\n", sep = "")
} else {
  cat("\n[D17L1] DRY RUN — no files written\n")
}

# ---- Optional: diagonal boundary scan --------------------------------------
# Slides a WxW upper-triangle cross-block at offset G along the diagonal.
# At each window position i:
#   cross = Z[(i-W+1):i, (i+G+1):(i+G+W)]
#   boundary_score[i] = -median(cross)   (bluer cross = higher score)
# Then 1D peak detection: a peak is a local maximum within +/- min_dist
# windows that exceeds boundary_score_min (in z-units).

if (boundary_scan) {
  cat("\n[D17L1] === DIAGONAL BOUNDARY SCAN ===\n", sep = "")
  cat("[D17L1]   W=", boundary_W, "  offset=", boundary_offset,
      "  score_min=", boundary_score_min,
      "  min_dist=", boundary_min_dist, "\n", sep = "")

  N <- n_windows_total
  W <- boundary_W
  G <- boundary_offset

  # Step 1: build per-diagonal mean and sd for Z normalization.
  cat("[D17L1]   precomputing diagonal mean/sd for Z\n")
  diag_mean <- numeric(N)
  diag_sd   <- numeric(N)
  for (d in 0:(N - 1L)) {
    if (d == 0L) {
      vals <- diag(sim_mat)
    } else {
      ii <- seq.int(1L, N - d)
      vals <- sim_mat[cbind(ii, ii + d)]
    }
    vals <- vals[is.finite(vals)]
    if (length(vals) >= 5L) {
      diag_mean[d + 1L] <- mean(vals)
      sg <- sd(vals)
      diag_sd[d + 1L]   <- if (is.finite(sg) && sg > 1e-9) sg else NA_real_
    } else {
      diag_mean[d + 1L] <- NA_real_
      diag_sd[d + 1L]   <- NA_real_
    }
  }

  # Step 2: at each diagonal position i, score the cross-block.
  # Valid i: cross-block must fit. Need (i - W + 1) >= 1 and (i + G + W) <= N.
  i_lo <- W
  i_hi <- N - G - W
  if (i_hi < i_lo) {
    cat("[D17L1]   chromosome too short for this W/offset — skipping scan\n")
  } else {
    centers <- seq.int(i_lo, i_hi)
    boundary_score <- rep(NA_real_, N)
    cat("[D17L1]   scoring ", length(centers),
        " positions along the diagonal\n", sep = "")
    for (ctr in centers) {
      ii <- seq.int(ctr - W + 1L, ctr)
      jj <- seq.int(ctr + G + 1L, ctr + G + W)
      block <- sim_mat[ii, jj]
      d_mat <- outer(ii, jj, FUN = function(a, b) b - a)
      idx <- d_mat + 1L
      mu <- diag_mean[idx]
      sg <- diag_sd[idx]
      z  <- (block - mu) / sg
      z[!is.finite(z)] <- NA_real_
      vals <- as.numeric(z)
      vals <- vals[is.finite(vals)]
      if (length(vals) >= 5L) {
        boundary_score[ctr] <- -median(vals)
      }
    }
    score_finite <- boundary_score[is.finite(boundary_score)]
    if (length(score_finite) > 0L) {
      cat(sprintf(
        "[D17L1]   boundary_score range: [%.2f, %.2f]   median=%.2f   p95=%.2f\n",
        min(score_finite), max(score_finite),
        median(score_finite),
        as.numeric(quantile(score_finite, 0.95))))
    }

    # Step 3: 1D peak detection.
    # A peak: (a) score >= boundary_score_min,
    #         (b) score is the max within +/- min_dist windows.
    peaks <- integer(0)
    for (ctr in centers) {
      s <- boundary_score[ctr]
      if (!is.finite(s) || s < boundary_score_min) next
      lo <- max(1L, ctr - boundary_min_dist)
      hi <- min(N, ctr + boundary_min_dist)
      win <- boundary_score[lo:hi]
      win <- win[is.finite(win)]
      if (length(win) == 0L) next
      if (s >= max(win)) peaks <- c(peaks, ctr)
    }
    cat("[D17L1]   peaks found: ", length(peaks), "\n", sep = "")

    # ---- Multi-W validation of each peak ----------------------------------
    # For each peak position i, scan two 1D rays of diag-Z values
    # perpendicular to the diagonal:
    #   right ray: Z[i, i+d]   for d = 1, 2, ..., d_max
    #   left  ray: Z[i-d, i]   for d = 1, 2, ..., d_max
    # By matrix symmetry these are the same data row/col-swapped, but we
    # keep both as separate rays in case one runs off the chromosome edge.
    #
    # A REAL boundary stays mostly blue (Z < 0) on both rays — the matrix
    # never rebounds to red as we look further from the diagonal. A FAKE
    # boundary (a small blue cross inside an inversion) shows blue near
    # the diagonal but jumps red abruptly as we exit the cross-shape.
    perp_z <- function(r, c) {
      # Single-pixel diagonal-distance Z at sim_mat[r, c]
      if (r < 1L || r > n_windows_total) return(NA_real_)
      if (c < 1L || c > n_windows_total) return(NA_real_)
      d <- abs(c - r)
      mu <- diag_mean[d + 1L]
      sg <- diag_sd[d + 1L]
      sm <- sim_mat[r, c]
      if (!is.finite(sm) || !is.finite(mu) || !is.finite(sg) || sg < 1e-9)
        return(NA_real_)
      (sm - mu) / sg
    }

    # Build output table.
    if (length(peaks) > 0L) {
      bdt <- data.table(
        chr            = chr_label,
        boundary_idx   = sprintf("%s_d17L1_b%04d", chr_label,
                                 seq_along(peaks)),
        boundary_w     = peaks,
        boundary_bp    = window_start_bp[peaks],
        boundary_score = boundary_score[peaks],
        boundary_W     = W,
        boundary_offset = G
      )
      bdt <- bdt[order(boundary_w)]

      if (boundary_validate) {
        cat("[D17L1]   validating peaks via perpendicular rays (d_max=",
            boundary_perp_d_max, ", red_z=", boundary_perp_red_z, ")\n",
            sep = "")
        # Per-peak per-ray stats
        n_peaks <- nrow(bdt)
        right_frac_blue <- numeric(n_peaks)
        left_frac_blue  <- numeric(n_peaks)
        right_max_z     <- numeric(n_peaks)
        left_max_z      <- numeric(n_peaks)
        right_first_red <- integer(n_peaks)
        left_first_red  <- integer(n_peaks)
        right_n_finite  <- integer(n_peaks)
        left_n_finite   <- integer(n_peaks)
        right_curves    <- character(n_peaks)
        left_curves     <- character(n_peaks)

        for (k in seq_len(n_peaks)) {
          i <- bdt$boundary_w[k]
          d_seq <- seq_len(boundary_perp_d_max)
          # Right ray
          r_vals <- vapply(d_seq, function(d) perp_z(i, i + d), numeric(1))
          # Left ray
          l_vals <- vapply(d_seq, function(d) perp_z(i - d, i), numeric(1))

          # Right stats
          ok_r <- is.finite(r_vals)
          right_n_finite[k] <- sum(ok_r)
          if (any(ok_r)) {
            right_frac_blue[k] <- mean(r_vals[ok_r] < 0)
            right_max_z[k]     <- max(r_vals[ok_r])
            red_idx <- which(r_vals > boundary_perp_red_z & ok_r)
            right_first_red[k] <- if (length(red_idx) > 0L)
                                    as.integer(red_idx[1]) else NA_integer_
          } else {
            right_frac_blue[k] <- NA_real_
            right_max_z[k]     <- NA_real_
            right_first_red[k] <- NA_integer_
          }
          # Left stats
          ok_l <- is.finite(l_vals)
          left_n_finite[k] <- sum(ok_l)
          if (any(ok_l)) {
            left_frac_blue[k] <- mean(l_vals[ok_l] < 0)
            left_max_z[k]     <- max(l_vals[ok_l])
            red_idx <- which(l_vals > boundary_perp_red_z & ok_l)
            left_first_red[k] <- if (length(red_idx) > 0L)
                                   as.integer(red_idx[1]) else NA_integer_
          } else {
            left_frac_blue[k] <- NA_real_
            left_max_z[k]     <- NA_real_
            left_first_red[k] <- NA_integer_
          }
          # Compact curve strings (every other pixel for readability)
          show_idx <- if (boundary_perp_d_max <= 10L) seq_len(boundary_perp_d_max)
                      else c(1L, seq.int(2L, boundary_perp_d_max - 1L, by = 2L),
                             boundary_perp_d_max)
          right_curves[k] <- paste(
            ifelse(is.finite(r_vals[show_idx]),
                   sprintf("%+.1f", r_vals[show_idx]), "NA"),
            collapse = ",")
          left_curves[k] <- paste(
            ifelse(is.finite(l_vals[show_idx]),
                   sprintf("%+.1f", l_vals[show_idx]), "NA"),
            collapse = ",")
        }

        bdt[, right_frac_blue := right_frac_blue]
        bdt[, left_frac_blue  := left_frac_blue]
        bdt[, right_max_z     := right_max_z]
        bdt[, left_max_z      := left_max_z]
        bdt[, right_first_red := right_first_red]
        bdt[, left_first_red  := left_first_red]
        bdt[, right_n_finite  := right_n_finite]
        bdt[, left_n_finite   := left_n_finite]
        bdt[, right_ray_curve := right_curves]
        bdt[, left_ray_curve  := left_curves]

        # Classify
        classify <- function(rfb, lfb, rmz, lmz, rfr, lfr) {
          # EDGE if either ray has too few finite samples
          if (!is.finite(rfb) || !is.finite(lfb)) return("EDGE")
          # DECAYS: either ray crosses red close to the diagonal
          first_red_too_close <-
            (is.finite(rfr) && rfr <= boundary_perp_first_d_red_max) ||
            (is.finite(lfr) && lfr <= boundary_perp_first_d_red_max)
          if (first_red_too_close) return("DECAYS")
          # STABLE_BLUE: both rays stay mostly blue + no big red excursion
          stable <- (rfb >= boundary_perp_min_blue_frac) &&
                    (lfb >= boundary_perp_min_blue_frac) &&
                    (is.finite(rmz) && rmz <= boundary_perp_max_red) &&
                    (is.finite(lmz) && lmz <= boundary_perp_max_red)
          if (stable) return("STABLE_BLUE")
          "MARGINAL"
        }
        bdt[, validation_status := mapply(classify,
                                          right_frac_blue, left_frac_blue,
                                          right_max_z,     left_max_z,
                                          right_first_red, left_first_red)]

        # Inside-L1 tag
        bdt[, inside_L1 := FALSE]
        bdt[, L1_id := NA_character_]
        if (nrow(nms_dt) > 0L) {
          for (k in seq_len(nrow(bdt))) {
            bw <- bdt$boundary_w[k]
            hit <- nms_dt[start_w <= bw & end_w >= bw]
            if (nrow(hit) > 0L) {
              bdt$inside_L1[k] <- TRUE
              bdt$L1_id[k]     <- hit$candidate_id[1]
            }
          }
        }

        n_stable <- sum(bdt$validation_status == "STABLE_BLUE")
        n_decay  <- sum(bdt$validation_status == "DECAYS")
        n_marg   <- sum(bdt$validation_status == "MARGINAL")
        n_edge   <- sum(bdt$validation_status == "EDGE")
        cat("[D17L1]   validation summary: STABLE_BLUE=", n_stable,
            "  DECAYS=", n_decay,
            "  MARGINAL=", n_marg,
            "  EDGE=", n_edge, "\n", sep = "")
        cat(sprintf("[D17L1]   right ray: median frac_blue=%.2f  max_z=%.2f\n",
                    median(bdt$right_frac_blue, na.rm = TRUE),
                    median(bdt$right_max_z, na.rm = TRUE)))
        cat(sprintf("[D17L1]   left  ray: median frac_blue=%.2f  max_z=%.2f\n",
                    median(bdt$left_frac_blue, na.rm = TRUE),
                    median(bdt$left_max_z, na.rm = TRUE)))

        n_stable_in_L1 <- sum(bdt$validation_status == "STABLE_BLUE" &
                              bdt$inside_L1)
        if (n_stable_in_L1 > 0L) {
          cat("[D17L1]   *** ", n_stable_in_L1,
              " STABLE_BLUE boundaries fall INSIDE L1 envelopes — ",
              "candidates for L1 split:\n", sep = "")
          show <- bdt[validation_status == "STABLE_BLUE" & inside_L1]
          for (k in seq_len(nrow(show))) {
            cat(sprintf("[D17L1]     %s  in %s  w=%d  %.2f Mb  trig=%.2f  rfb=%.2f  lfb=%.2f  rmz=%.2f  lmz=%.2f\n",
                        show$boundary_idx[k],
                        show$L1_id[k],
                        show$boundary_w[k],
                        show$boundary_bp[k] / 1e6,
                        show$boundary_score[k],
                        show$right_frac_blue[k],
                        show$left_frac_blue[k],
                        show$right_max_z[k],
                        show$left_max_z[k]))
          }
        }
      }

      cat("[D17L1]   top 10 boundaries by trigger score:\n")
      top <- bdt[order(-boundary_score)][1:min(10L, .N)]
      for (k in seq_len(nrow(top))) {
        if ("validation_status" %in% names(top)) {
          cat(sprintf("[D17L1]     %s  w=%d  %.2f Mb  trig=%.2f  status=%-11s  rfb=%.2f  lfb=%.2f  rmz=%+.2f  lmz=%+.2f\n",
                      top$boundary_idx[k], top$boundary_w[k],
                      top$boundary_bp[k] / 1e6, top$boundary_score[k],
                      top$validation_status[k],
                      top$right_frac_blue[k], top$left_frac_blue[k],
                      top$right_max_z[k],     top$left_max_z[k]))
        } else {
          cat(sprintf("[D17L1]     %s  w=%d  %.2f Mb  score=%.2f\n",
                      top$boundary_idx[k], top$boundary_w[k],
                      top$boundary_bp[k] / 1e6, top$boundary_score[k]))
        }
      }
    } else {
      bdt <- data.table(
        chr = character(0), boundary_idx = character(0),
        boundary_w = integer(0), boundary_bp = integer(0),
        boundary_score = numeric(0),
        boundary_W = integer(0), boundary_offset = integer(0)
      )
    }

    if (!dry_run) {
      out_b <- file.path(outdir, paste0(chr_label, "_d17L1_boundaries.tsv"))
      fwrite(bdt, out_b, sep = "\t")
      cat("[D17L1]   boundaries written: ", out_b, "\n", sep = "")

      # Also write the full per-window score curve, for diagnostics.
      curve_dt <- data.table(
        chr = chr_label,
        window_idx = seq_len(N),
        bp = window_start_bp,
        boundary_score = boundary_score
      )
      out_c <- file.path(outdir,
                         paste0(chr_label, "_d17L1_boundary_score_curve.tsv"))
      fwrite(curve_dt, out_c, sep = "\t")
      cat("[D17L1]   score curve written: ", out_c, "\n", sep = "")
    }
  }
}

cat("[D17L1] done\n")
