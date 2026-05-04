#!/usr/bin/env Rscript
# =============================================================================
# STEP_D17_squares_multi_scale.R
#
# Bottom-up multi-scale square detection on a smoothed sim_mat, with:
#   (1) adaptive size-dependent similarity thresholds (big blocks pass at lower
#       sim than small blocks) so we don't bias against the largest inversions
#   (2) size-biased non-maximum suppression (when blocks overlap, prefer the
#       LARGER one — the inverse of typical NMS, because we want big blocks
#       to win over their internal sub-blocks)
#   (3) per-candidate sample-composition validation (k=3 PC1 -> per-sample
#       modal band -> consistency score) to falsify candidates that are not
#       one biological signal
#
# Inputs:
#   --precomp     <slim precomp .rds with pc$dt PC_1_* columns + sample order>
#   --sim_mat     <smoothed sim_mat .rds, default tries sim_mat_nn80.rds in
#                  the precomp dir>
#   --chr         <chromosome label for output rows>
#   --outdir      <output directory>
#
# Output:
#   <outdir>/<chr>_d17_candidates.tsv with columns:
#     chr, candidate_id, start_w, end_w, start_bp, end_bp, n_windows,
#     scale_W, mean_sim, density_p70, composition_consistency,
#     n_inv_carriers, status, parent_candidate_id
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
sim_mat_f      <- get_arg("--sim_mat")          # may be NA -> auto
chr_label      <- get_arg("--chr", "chr")
outdir         <- get_arg("--outdir", ".")

# Block-size scan grid (in windows). Defaults span 50w (~0.25 Mb) to 800w (~4 Mb).
sizes_str      <- get_arg("--sizes", "50,100,200,400,800")
sizes_W        <- as.integer(strsplit(sizes_str, ",")[[1]])

# Threshold parameters
sim_thresh_min <- as.numeric(get_arg("--sim_thresh_min", "0.50"))
sim_thresh_max <- as.numeric(get_arg("--sim_thresh_max", "0.70"))
sim_thresh_tau <- as.numeric(get_arg("--sim_thresh_tau", "200"))
density_min    <- as.numeric(get_arg("--density_min", "0.55"))

# NMS overlap fraction at which two candidates are considered the "same"
nms_overlap    <- as.numeric(get_arg("--nms_overlap", "0.50"))

# Composition validation parameters
comp_min       <- as.numeric(get_arg("--comp_min", "0.65"))
comp_kmeans_k  <- as.integer(get_arg("--comp_kmeans_k", "3"))

# Step size for the diagonal scan (in windows). Smaller = denser scan = slower
scan_step      <- as.integer(get_arg("--scan_step", "5"))

# Min candidate size to keep after everything
min_n_windows  <- as.integer(get_arg("--min_n_windows", "20"))

dry_run        <- !is.na(get_arg("--dry_run", NA_character_))

if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("--precomp is required and must exist")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[D17] precomp:    ", precomp_f, "\n")
cat("[D17] sim_mat:    ", sim_mat_f %||% "(auto)", "\n")
cat("[D17] chr:        ", chr_label, "\n")
cat("[D17] outdir:     ", outdir, "\n")
cat("[D17] sizes_W:    ", paste(sizes_W, collapse = ","), "\n")
cat("[D17] thresholds: sim_min=", sim_thresh_min, " sim_max=", sim_thresh_max,
    " tau=", sim_thresh_tau, " density_min=", density_min, "\n", sep = "")
cat("[D17] NMS:        overlap=", nms_overlap, "\n", sep = "")
cat("[D17] comp:       min=", comp_min, " k=", comp_kmeans_k, "\n", sep = "")

# ---- Adaptive size-dependent similarity threshold ---------------------------
# A small block (W=50) needs higher similarity (~sim_max) to count as a real
# square. A large block (W=800) only needs sim_min, because big inversions
# at hatchery-MAF naturally have lower mean similarity than tight small ones
# (more breakpoint heterogeneity, more recombination, more dilution).
#
#   threshold(W) = sim_min + (sim_max - sim_min) * exp(-W / tau)
#
# Defaults (sim_min=0.50, sim_max=0.70, tau=200) give:
#   W=50  -> 0.66
#   W=100 -> 0.62
#   W=200 -> 0.57
#   W=400 -> 0.53
#   W=800 -> 0.51

size_threshold <- function(W) {
  sim_thresh_min + (sim_thresh_max - sim_thresh_min) * exp(-W / sim_thresh_tau)
}

cat("[D17] adaptive thresholds by size:\n")
for (W in sizes_W) {
  cat(sprintf("[D17]   W=%4d  threshold=%.3f\n", W, size_threshold(W)))
}

# ---- Load precomp -----------------------------------------------------------

cat("[D17] loading precomp\n")
pc_obj <- readRDS(precomp_f)
# Robust unpacking: precomp may have $pc or be the pc list directly
if (!is.null(pc_obj$dt)) {
  pc <- pc_obj
} else if (!is.null(pc_obj$pc)) {
  pc <- pc_obj$pc
} else {
  stop("[D17] cannot find pc$dt in precomp file")
}

dt_pc <- as.data.table(pc$dt)
n_windows_total <- nrow(dt_pc)
cat("[D17] N windows: ", n_windows_total, "\n")

# Sample IDs from PC_1_* columns
pc1_cols <- grep("^PC_1_", names(dt_pc), value = TRUE)
sample_ids <- sub("^PC_1_", "", pc1_cols)
n_samples <- length(sample_ids)
cat("[D17] N samples: ", n_samples, "\n")
if (n_samples < 50) stop("[D17] too few samples — sanity check failed")

# Window genomic coordinates
if (!all(c("start", "end") %in% names(dt_pc))) {
  stop("[D17] precomp dt is missing 'start'/'end' columns")
}
window_start_bp <- dt_pc$start
window_end_bp   <- dt_pc$end

# PC1 matrix: windows × samples
PC1 <- as.matrix(dt_pc[, ..pc1_cols])
storage.mode(PC1) <- "double"

# ---- Load smoothed sim_mat --------------------------------------------------

if (is.na(sim_mat_f)) {
  candidates_path <- c(
    file.path(dirname(precomp_f), "sim_mat_nn80.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn160.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn40.rds")
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0)
    stop("[D17] --sim_mat not given and could not auto-find sim_mat_nn*.rds in precomp dir")
  sim_mat_f <- hit[1]
  cat("[D17] auto-found sim_mat: ", sim_mat_f, "\n")
}
if (!file.exists(sim_mat_f))
  stop("[D17] --sim_mat does not exist: ", sim_mat_f)

cat("[D17] loading sim_mat\n")
sm_obj <- readRDS(sim_mat_f)
# Robust unpacking
if (is.matrix(sm_obj)) {
  sim_mat <- sm_obj
} else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat) && is.matrix(sm_obj$sim_mat)) {
  sim_mat <- sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) == 1 && is.matrix(sm_obj[[1]])) {
  sim_mat <- sm_obj[[1]]
} else {
  stop("[D17] sim_mat object structure not recognized; class=", class(sm_obj)[1])
}
storage.mode(sim_mat) <- "double"
if (!isTRUE(nrow(sim_mat) == n_windows_total &&
            ncol(sim_mat) == n_windows_total)) {
  stop("[D17] sim_mat dim ", nrow(sim_mat), "x", ncol(sim_mat),
       " does not match precomp N windows ", n_windows_total)
}

# ---- Multi-scale block scanner ----------------------------------------------
# For each block size W in sizes_W, walk the diagonal at step `scan_step`.
# At each position p, examine the W x W block on the diagonal centered on p.
# Compute mean_sim and density_p70 inside the block. Keep blocks where:
#   mean_sim   >= size_threshold(W)
#   density_p70 >= density_min          (fraction of pixels with sim >= 0.70)

cat("[D17] scanning diagonal at multiple scales\n")
scan_results <- list()

for (W in sizes_W) {
  thr_W <- size_threshold(W)
  half <- W %/% 2L
  # Start positions are window indices (1-based) of the block center.
  centers <- seq(half + 1L, n_windows_total - half, by = scan_step)
  if (length(centers) == 0L) {
    cat(sprintf("[D17]   size W=%4d: chromosome too short, skipping\n", W))
    next
  }

  out <- vector("list", length(centers))
  out_i <- 1L
  for (ctr in centers) {
    s <- ctr - half
    e <- s + W - 1L
    if (s < 1L || e > n_windows_total) next
    block <- sim_mat[s:e, s:e]
    msim <- mean(block, na.rm = TRUE)
    if (!is.finite(msim) || msim < thr_W) next
    # Density: fraction of pixels at OR above the size-adaptive threshold,
    # NOT a fixed 0.70 cutoff. This avoids biasing against big blocks whose
    # within-block similarity is genuinely diluted (mean ~0.6) but uniform.
    # A loose-envelope region will still fail because its pixel distribution
    # has a long left tail well below thr_W.
    dens <- mean(block >= thr_W, na.rm = TRUE)
    if (!is.finite(dens) || dens < density_min) next
    # Also report fraction above a fixed strong cutoff (0.70) for the output
    # table — useful for diagnostics but not for filtering.
    dens_p70 <- mean(block >= 0.70, na.rm = TRUE)
    out[[out_i]] <- data.table(
      W = W,
      start_w = s,
      end_w = e,
      n_windows = W,
      mean_sim = msim,
      density_p70 = dens_p70,
      density_adaptive = dens,
      threshold = thr_W
    )
    out_i <- out_i + 1L
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) > 0L) {
    scan_results[[as.character(W)]] <- rbindlist(out)
    cat(sprintf("[D17]   size W=%4d (thr=%.3f): %d candidate windows passed\n",
                W, thr_W, nrow(scan_results[[as.character(W)]])))
  } else {
    cat(sprintf("[D17]   size W=%4d (thr=%.3f): 0 passed\n", W, thr_W))
  }
}

if (length(scan_results) == 0L) {
  cat("[D17] no candidates at any scale; writing empty output and stopping\n")
  empty_dt <- data.table(
    chr = character(0), candidate_id = character(0),
    start_w = integer(0), end_w = integer(0),
    start_bp = integer(0), end_bp = integer(0),
    n_windows = integer(0), scale_W = integer(0),
    mean_sim = numeric(0), density_p70 = numeric(0),
    composition_consistency = numeric(0), n_inv_carriers = integer(0),
    status = character(0), parent_candidate_id = character(0)
  )
  fwrite(empty_dt, file.path(outdir, paste0(chr_label, "_d17_candidates.tsv")),
         sep = "\t")
  quit(save = "no", status = 0)
}

scan_dt <- rbindlist(scan_results)
cat("[D17] total scan hits across scales: ", nrow(scan_dt), "\n", sep = "")

# Within each scale, collapse adjacent overlapping hits into one block. Two
# scan hits at the same scale that overlap by more than 50% become one block
# spanning their union. Take the best mean_sim as representative.

collapse_within_scale <- function(dt) {
  if (nrow(dt) == 0L) return(dt)
  dt <- dt[order(start_w)]
  out <- list()
  cur <- dt[1]
  for (i in seq.int(2L, nrow(dt))) {
    row_i <- dt[i]
    overlap_w <- max(0L, min(cur$end_w, row_i$end_w) - max(cur$start_w, row_i$start_w) + 1L)
    smaller_w <- min(cur$end_w - cur$start_w + 1L,
                     row_i$end_w - row_i$start_w + 1L)
    if (overlap_w / smaller_w > nms_overlap) {
      cur$start_w <- min(cur$start_w, row_i$start_w)
      cur$end_w   <- max(cur$end_w,   row_i$end_w)
      cur$n_windows <- cur$end_w - cur$start_w + 1L
      cur$mean_sim <- max(cur$mean_sim, row_i$mean_sim)
      cur$density_p70 <- max(cur$density_p70, row_i$density_p70)
      cur$density_adaptive <- max(cur$density_adaptive, row_i$density_adaptive)
    } else {
      out[[length(out) + 1L]] <- cur
      cur <- row_i
    }
  }
  out[[length(out) + 1L]] <- cur
  rbindlist(out)
}

scan_dt_split <- split(scan_dt, scan_dt$W)
collapsed_per_scale <- lapply(scan_dt_split, collapse_within_scale)
scan_dt <- rbindlist(collapsed_per_scale)
cat("[D17] after within-scale collapse: ", nrow(scan_dt), "\n", sep = "")

# ---- Size-biased non-maximum suppression across scales ---------------------
# When candidates from different scales overlap, prefer the LARGER one. This
# is the inverse of typical NMS and is the key bias change vs. the existing
# tree builder. The reasoning: a 50w core inside a 400w real inversion is
# not a separate inversion — it's a piece of the bigger one. The bigger one
# is the right candidate.
#
# Algorithm: sort by W descending. For each candidate (largest first), if
# any already-kept candidate overlaps with it by >50%, discard. Otherwise
# keep.
# A candidate overlaps with a kept one if min(overlap_w/cand_w, overlap_w/kept_w)
# > nms_overlap. This way a 400w block "absorbs" any 50w blocks fully inside it.

cat("[D17] running size-biased NMS\n")
scan_dt <- scan_dt[order(-W, -mean_sim)]
kept <- list()
for (i in seq.int(1L, nrow(scan_dt))) {
  row_i <- scan_dt[i]
  is_dominated <- FALSE
  if (length(kept) > 0L) {
    for (k in kept) {
      o <- max(0L, min(k$end_w, row_i$end_w) - max(k$start_w, row_i$start_w) + 1L)
      smaller <- min(k$end_w - k$start_w + 1L, row_i$end_w - row_i$start_w + 1L)
      if (o / smaller > nms_overlap) {
        is_dominated <- TRUE
        break
      }
    }
  }
  if (!is_dominated) kept[[length(kept) + 1L]] <- row_i
}
nms_dt <- rbindlist(kept)
nms_dt <- nms_dt[order(start_w)]
cat("[D17] after size-biased NMS: ", nrow(nms_dt), " candidates\n", sep = "")

# ---- Sample-composition validation ------------------------------------------
# For each candidate, run k=3 k-means on each window's PC1 vector across
# samples. Each window assigns each sample to one of 3 bands. For each
# sample, compute its modal band across all windows in the candidate. The
# consistency score = mean over samples of (modal-band count) / (n_windows).
#
# Real inversion: ~0.85+ (samples consistently in same band)
# Over-merged region: ~0.5-0.65 (different windows sort samples differently)
# Pure noise: ~0.4 (random)
#
# Also report n_inv_carriers = size of the smallest of the 3 bands at the
# representative middle window of the candidate. This is a rough proxy for
# the rare-arrangement carrier set size.

per_window_kmeans_band <- function(pc1_vec, k = 3L) {
  # Returns integer vector length n_samples giving band assignment 1..k
  # ordered by cluster center (band 1 = lowest PC1 mean, band k = highest).
  if (any(!is.finite(pc1_vec))) {
    pc1_vec[!is.finite(pc1_vec)] <- 0
  }
  if (sd(pc1_vec) < 1e-9) {
    return(rep(1L, length(pc1_vec)))
  }
  km <- tryCatch(
    suppressWarnings(kmeans(pc1_vec, centers = k, nstart = 5L, iter.max = 30L)),
    error = function(e) NULL
  )
  if (is.null(km)) return(rep(1L, length(pc1_vec)))
  # Re-order so band index reflects PC1 magnitude
  ord <- order(km$centers[, 1])
  band_remap <- match(seq_len(k), ord)
  band_remap[km$cluster]
}

composition_consistency <- function(start_w, end_w) {
  ws <- seq.int(start_w, end_w)
  if (length(ws) < 5L) return(list(consistency = NA_real_, n_inv = NA_integer_))
  band_mat <- matrix(NA_integer_, nrow = n_samples, ncol = length(ws))
  for (j in seq_along(ws)) {
    band_mat[, j] <- per_window_kmeans_band(PC1[ws[j], ], k = comp_kmeans_k)
  }
  # Per-sample modal band fraction
  modal_frac <- vapply(seq_len(n_samples), function(s) {
    tab <- tabulate(band_mat[s, ], nbins = comp_kmeans_k)
    max(tab) / length(ws)
  }, numeric(1))
  consistency <- mean(modal_frac, na.rm = TRUE)
  # Carrier set size: in the middle window, smallest band cluster size
  mid <- ws[length(ws) %/% 2L + 1L]
  mid_bands <- per_window_kmeans_band(PC1[mid, ], k = comp_kmeans_k)
  band_sizes <- tabulate(mid_bands, nbins = comp_kmeans_k)
  n_inv <- min(band_sizes[band_sizes > 0])
  list(consistency = consistency, n_inv = as.integer(n_inv))
}

cat("[D17] running composition validation on ", nrow(nms_dt),
    " NMS-survived candidates\n", sep = "")
comp_consistencies <- numeric(nrow(nms_dt))
n_inv_carriers <- integer(nrow(nms_dt))
for (i in seq.int(1L, nrow(nms_dt))) {
  res <- composition_consistency(nms_dt$start_w[i], nms_dt$end_w[i])
  comp_consistencies[i] <- res$consistency
  n_inv_carriers[i] <- res$n_inv
  if (i %% 10L == 0L) {
    cat(sprintf("[D17]   ... %d/%d done\n", i, nrow(nms_dt)))
  }
}
nms_dt[, composition_consistency := comp_consistencies]
nms_dt[, n_inv_carriers := n_inv_carriers]

# ---- Two-pass replacement: failed big -> surviving smaller children --------
# A big candidate that fails composition_consistency probably represents an
# over-merge of distinct biological signals. In that case, the small
# candidates that overlap with it (and that we suppressed during NMS) might
# each be real. Re-introduce them if they pass composition.

failed_big <- nms_dt[composition_consistency < comp_min & W >= 200L]
if (nrow(failed_big) > 0L) {
  cat("[D17] ", nrow(failed_big), " big candidates failed composition; ",
      "checking smaller alternates\n", sep = "")
  # Re-scan: from the original scan_dt (pre-NMS, pre-collapse), find smaller
  # candidates that overlap any failed_big region. We need to rebuild from
  # scan_dt before collapse to give NMS-suppressed candidates a second look.
  # Use the within-scale-collapsed scan_dt as the source of "alternates".
  alternates_all <- rbindlist(collapsed_per_scale)
  alternates_all <- alternates_all[W < 200L]  # smaller scales only
  rescues <- list()
  for (i in seq.int(1L, nrow(failed_big))) {
    fb <- failed_big[i]
    inside <- alternates_all[start_w >= fb$start_w & end_w <= fb$end_w]
    if (nrow(inside) == 0L) next
    # Run composition on each
    for (j in seq.int(1L, nrow(inside))) {
      res <- composition_consistency(inside$start_w[j], inside$end_w[j])
      if (!is.na(res$consistency) && res$consistency >= comp_min) {
        rescued_row <- inside[j]
        rescued_row[, composition_consistency := res$consistency]
        rescued_row[, n_inv_carriers := res$n_inv]
        rescues[[length(rescues) + 1L]] <- rescued_row
      }
    }
  }
  if (length(rescues) > 0L) {
    rescue_dt <- rbindlist(rescues)
    # Apply size-biased NMS to the rescues themselves (in case multiple
    # alternates inside one failed_big overlap each other)
    rescue_dt <- rescue_dt[order(-W, -mean_sim)]
    kept_r <- list()
    for (i in seq.int(1L, nrow(rescue_dt))) {
      row_i <- rescue_dt[i]
      dom <- FALSE
      for (k in kept_r) {
        o <- max(0L, min(k$end_w, row_i$end_w) - max(k$start_w, row_i$start_w) + 1L)
        sm <- min(k$end_w - k$start_w + 1L, row_i$end_w - row_i$start_w + 1L)
        if (o / sm > nms_overlap) { dom <- TRUE; break }
      }
      if (!dom) kept_r[[length(kept_r) + 1L]] <- row_i
    }
    rescue_dt <- rbindlist(kept_r)
    cat("[D17] rescued ", nrow(rescue_dt), " smaller candidates\n", sep = "")
    # Mark the failed_big as REPLACED, mark rescues as RESCUED with parent
    nms_dt[, parent_candidate_id := NA_character_]
    nms_dt[, status := ifelse(composition_consistency >= comp_min, "PASS",
                              ifelse(W >= 200L, "REPLACED", "FAIL"))]
    if (nrow(rescue_dt) > 0L) {
      rescue_dt[, status := "RESCUED"]
      rescue_dt[, parent_candidate_id := NA_character_]   # filled below
      nms_dt <- rbindlist(list(nms_dt, rescue_dt), use.names = TRUE, fill = TRUE)
    }
  } else {
    nms_dt[, parent_candidate_id := NA_character_]
    nms_dt[, status := ifelse(composition_consistency >= comp_min, "PASS", "FAIL")]
  }
} else {
  nms_dt[, parent_candidate_id := NA_character_]
  nms_dt[, status := ifelse(composition_consistency >= comp_min, "PASS", "FAIL")]
}

# Drop tiny candidates
nms_dt <- nms_dt[n_windows >= min_n_windows]

# Re-sort by genomic position
nms_dt <- nms_dt[order(start_w)]

# ---- Build final output -----------------------------------------------------

nms_dt[, chr := chr_label]
nms_dt[, candidate_id := sprintf("%s_d17_%04d", chr_label, seq_len(.N))]
nms_dt[, start_bp := window_start_bp[start_w]]
nms_dt[, end_bp   := window_end_bp[end_w]]
nms_dt[, recursion_depth := 0L]
setnames(nms_dt, "W", "scale_W")

out_cols <- c(
  "chr", "candidate_id", "parent_candidate_id", "recursion_depth",
  "start_w", "end_w",
  "start_bp", "end_bp", "n_windows", "scale_W",
  "mean_sim", "density_p70", "density_adaptive",
  "composition_consistency", "n_inv_carriers",
  "status"
)
out_dt <- nms_dt[, ..out_cols]

# ---- Console summary --------------------------------------------------------

cat("\n[D17] === SUMMARY for ", chr_label, " ===\n", sep = "")
cat("[D17] total candidates: ", nrow(out_dt), "\n", sep = "")
cat("[D17] by status:\n")
print(out_dt[, .N, by = status])
cat("[D17] by scale_W (PASS only):\n")
print(out_dt[status == "PASS", .N, by = scale_W][order(scale_W)])
cat("[D17] mean_sim quartiles (PASS only):\n")
if (nrow(out_dt[status == "PASS"]) > 0L) {
  print(quantile(out_dt[status == "PASS"]$mean_sim,
                 c(0.25, 0.5, 0.75)))
}
cat("[D17] composition_consistency quartiles (PASS only):\n")
if (nrow(out_dt[status == "PASS"]) > 0L) {
  print(quantile(out_dt[status == "PASS"]$composition_consistency,
                 c(0.25, 0.5, 0.75)))
}

if (nrow(out_dt[status == "PASS"]) > 0L) {
  cat("[D17] biggest PASS candidates by n_windows:\n")
  big <- out_dt[status == "PASS"][order(-n_windows)][1:min(10L, .N)]
  for (i in seq_len(nrow(big))) {
    cat(sprintf(
      "[D17]   %s  %.2f-%.2f Mb  nW=%d  mean_sim=%.3f  comp=%.3f  carriers=%d\n",
      big$candidate_id[i],
      big$start_bp[i] / 1e6, big$end_bp[i] / 1e6,
      big$n_windows[i], big$mean_sim[i],
      big$composition_consistency[i], big$n_inv_carriers[i]
    ))
  }
}

# ---- Write -----------------------------------------------------------------

if (!dry_run) {
  outpath <- file.path(outdir, paste0(chr_label, "_d17_candidates.tsv"))
  fwrite(out_dt, outpath, sep = "\t")
  cat("\n[D17] candidates written: ", outpath, "\n", sep = "")
} else {
  cat("\n[D17] DRY RUN — no file written\n")
}

cat("[D17] done\n")
