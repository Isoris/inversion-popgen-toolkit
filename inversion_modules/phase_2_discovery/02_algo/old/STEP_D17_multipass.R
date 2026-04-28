#!/usr/bin/env Rscript
# =============================================================================
# STEP_D17_multipass.R
#
# Multi-level discovery + Ward.D2 internal clustering for catfish inversion
# candidates. Two architectures supported via --architecture:
#
# ARCHITECTURE = "nested" (default):
#
#   Level 1 (W=800/400) : envelope discovery, NO masking, NO composition gate.
#                         Output: large regions of elevated similarity that
#                         may contain one or more inversions.
#   Level 2 (W=200/100) : inversion discovery, scanned INSIDE each Level 1
#                         envelope. Composition gate applied here.
#                         Output: actual inversion candidates.
#   Level 3 (W=50, opt) : sub-feature discovery inside each Level 2 inversion.
#   Ward.D2             : runs on Level 2 inversions to characterize internal
#                         substructure (CLEAN / NESTED_ADJ / NESTED_INSIDE /
#                         FRAGMENTED / DOUBLE_CROSSOVER).
#                         FRAGMENTED is a label, not a discard. With
#                         --fragmented_emit, FRAGMENTED parents emit their
#                         large contiguous Ward runs as salvage children
#                         (status = PASS_CHILD_FRAG).
#   Catalogue           : Level 1 envelopes (status=ENVELOPE, depth=0) +
#                         Level 2 inversions (depth=1) + Ward children (depth=2).
#
# ARCHITECTURE = "horizontal" (legacy):
#   Pass 1 / 2 / 3 are mutually exclusive in space via masking. Ward runs on
#   Pass 1 + Pass 2 candidates. FRAGMENTED is a discard reason. Kept for
#   backwards compatibility and reference.
#
# OUTPUT (per chromosome):
#   <chr>_d17_catalogue.tsv         : main catalogue, all rows
#   <chr>_d17_subblocks.tsv         : per-cluster summary inside each cand
#   <chr>_d17_window_clusters.tsv   : per-window cluster ID for clustered
#                                     candidates
#   <chr>_d17_z_audit.tsv           : pre/post-merge child rows from Z
#                                     validator, when --z_validate is on
#   <chr>_phase_qc_candidates.tsv   : clean filter
#                                     (PASS + PASS_CHILD + PASS_CHILD_FRAG +
#                                      PASS_COMPLEX) only, ready for Phase QC
#
# Author: Claude (Anthropic)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)   # for silhouette()
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

# Three-pass W groups
pass1_W_str    <- get_arg("--pass1_W", "800,400")
pass2_W_str    <- get_arg("--pass2_W", "200,100")
pass3_W_str    <- get_arg("--pass3_W", "50")
pass1_W <- as.integer(strsplit(pass1_W_str, ",")[[1]])
pass2_W <- as.integer(strsplit(pass2_W_str, ",")[[1]])
pass3_W <- as.integer(strsplit(pass3_W_str, ",")[[1]])

# Adaptive sim_mat threshold parameters (size-dependent)
sim_thresh_min <- as.numeric(get_arg("--sim_thresh_min", "0.50"))
sim_thresh_max <- as.numeric(get_arg("--sim_thresh_max", "0.70"))
sim_thresh_tau <- as.numeric(get_arg("--sim_thresh_tau", "200"))
density_min    <- as.numeric(get_arg("--density_min", "0.55"))

# NMS overlap fraction
nms_overlap    <- as.numeric(get_arg("--nms_overlap", "0.50"))

# Composition validation
comp_min       <- as.numeric(get_arg("--comp_min", "0.65"))
comp_kmeans_k  <- as.integer(get_arg("--comp_kmeans_k", "3"))

# Diagonal-scan step size (windows)
scan_step      <- as.integer(get_arg("--scan_step", "5"))

# Min candidate size that survives final filter (for FAIL of any pass)
min_n_windows  <- as.integer(get_arg("--min_n_windows", "20"))

# Optional core trimming after NMS, before composition validation.
# This removes weak edge windows from a candidate using row-mean similarity.
trim_core      <- !is.na(get_arg("--trim_core", NA_character_))
trim_drop_z    <- as.numeric(get_arg("--trim_drop_z", "1.0"))

# Post-trim NMS: re-run overlap suppression after trim_to_core has tightened
# parent boundaries. Catches sibling parents that became distinct only after
# trimming and would otherwise yield duplicate Ward-clustered children.
post_trim_nms          <- !is.na(get_arg("--post_trim_nms", NA_character_))
post_trim_nms_overlap  <- as.numeric(get_arg("--post_trim_nms_overlap", "0.30"))

# Z-aware validation of NESTED_ADJ splits. After Ward emits children, we look
# at the off-diagonal cross-block between adjacent children. Cross-block Z
# uses diagonal-distance normalization (same as the overlay plotter):
#   - strongly negative cross-Z (children are LESS similar than expected at
#     their inter-window distance) ==> the split is real, keep children.
#   - strongly positive cross-Z (children are MORE similar than expected,
#     i.e. they co-segregate at long range) ==> the split is spurious; the
#     parent is one inversion that Ward over-segmented. Children are merged
#     back into a single PASS_CLEAN child spanning the parent.
# Both children always receive a z_flag column ("OK" / "OVER_SPLIT" / "REAL")
# and the original (pre-merge) child rows are also dumped to a sidecar TSV
# for audit. So the validator NEVER silently discards information.
z_validate         <- !is.na(get_arg("--z_validate", NA_character_))
z_merge_thresh     <- as.numeric(get_arg("--z_merge_thresh", "0.5"))
z_split_thresh     <- as.numeric(get_arg("--z_split_thresh", "-0.5"))

# ---- Architecture flag ------------------------------------------------------
# "nested" (default): Level 1 envelopes -> Level 2 inversions inside each
#                     envelope -> Ward.D2 inside each inversion. Levels are
#                     spatially nested. Pass 1 envelopes use NO composition
#                     gate; only Level 2 inversions do. FRAGMENTED is a
#                     pattern label, not a discard reason.
# "horizontal": legacy behaviour. Pass 1/2/3 are mutually exclusive in space
#               via masking; Ward.D2 runs on Pass 1 + Pass 2 candidates.
architecture <- get_arg("--architecture", "nested")
if (!architecture %in% c("nested", "horizontal"))
  stop("[D17mp] --architecture must be 'nested' or 'horizontal'")

# Level-1 envelope discovery: looser than Level 2.
# The envelope mean_sim threshold scales like the W=800 threshold but uses
# its own minimum to allow "elevated but heterogeneous" regions to qualify.
l1_density_min <- as.numeric(get_arg("--l1_density_min", "0.45"))

# Level-3 (W=50, optional): scan inside L2 inversions for sub-features.
l3_enabled     <- !is.na(get_arg("--l3_enabled", NA_character_))

# Allow FRAGMENTED parents to emit large contiguous Ward runs as children.
# In the nested architecture this is helpful: a FRAGMENTED L2 inversion may
# actually contain two real sub-inversions that Ward's K-best couldn't
# separate cleanly. The salvage requires runs >= fragmented_min_run_w.
fragmented_emit       <- !is.na(get_arg("--fragmented_emit", NA_character_))
fragmented_min_run_w  <- as.integer(get_arg("--fragmented_min_run_w", "40"))

# Level 2 method: "scan" (default, fixed-W scan_region) or "refiner"
# (adaptive multi-W boundary refiner). The refiner finds boundary positions
# inside an envelope by asking, at each candidate boundary b, whether the
# left half and right half are internally coherent AND well separated from
# each other (using diagonal-distance Z for the cross term). Boundaries
# stable across multiple W values become block edges; blocks are then
# composition-validated like ordinary L2 inversions.
l2_method        <- get_arg("--l2_method", "scan")
if (!l2_method %in% c("scan", "refiner"))
  stop("[D17mp] --l2_method must be 'scan' or 'refiner'")
refiner_W_str    <- get_arg("--refiner_W", "200,100,50,25")
refiner_W_vec    <- as.integer(strsplit(refiner_W_str, ",")[[1]])
refiner_step_frac <- as.numeric(get_arg("--refiner_step_frac", "0.10"))
# refiner_score_min: minimum boundary score required to call a peak. Units
# depend on cross_term:
#   cross_term="z"   : score = coherence * (-mean(diag-Z));  use ~0.05
#   cross_term="zq"  : score = -quantile(diag-Z, q);         use ~1.0 (z units)
#   cross_term="raw" : score = coherence - mean(raw cross);  use ~0.05
refiner_score_min <- as.numeric(get_arg("--refiner_score_min", "0.05"))
refiner_min_scales <- as.integer(get_arg("--refiner_min_scales", "2"))
refiner_max_depth  <- as.integer(get_arg("--refiner_max_depth", "2"))
refiner_min_block  <- as.integer(get_arg("--refiner_min_block", "30"))
# Refiner cross term:
#   "z"   : score = coherence * (-mean(diag-Z over cross-block))
#   "zq"  : score = -quantile(diag-Z over cross-block, refiner_z_quantile)
#           — robust to outliers, peaks where the BOTTOM of the cross-block
#           Z distribution dips strongly negative (the "blue valley" your
#           eye sees at boundaries).
#   "raw" : score = coherence - mean(raw cross-block sim)
refiner_cross_term <- get_arg("--refiner_cross_term", "z")
if (!refiner_cross_term %in% c("z", "zq", "raw"))
  stop("[D17mp] --refiner_cross_term must be 'z', 'zq', or 'raw'")
refiner_z_quantile <- as.numeric(get_arg("--refiner_z_quantile", "0.05"))

# Internal Ward.D2 clustering
cluster_min_W  <- as.integer(get_arg("--cluster_min_W", "100"))
k_min          <- as.integer(get_arg("--k_min", "2"))
k_max          <- as.integer(get_arg("--k_max", "5"))
silhouette_min <- as.numeric(get_arg("--silhouette_min", "0.20"))
clean_dom_frac <- as.numeric(get_arg("--clean_dom_frac", "0.85"))
crossover_runs <- as.numeric(get_arg("--crossover_runs", "0.30"))

dry_run        <- !is.na(get_arg("--dry_run", NA_character_))

if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("[D17mp] --precomp is required and must exist")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[D17mp] precomp:    ", precomp_f, "\n")
cat("[D17mp] sim_mat:    ", sim_mat_f %||% "(auto)", "\n")
cat("[D17mp] chr:        ", chr_label, "\n")
cat("[D17mp] outdir:     ", outdir, "\n")
cat("[D17mp] pass 1 (large)  W: ", paste(pass1_W, collapse = ","), "\n")
cat("[D17mp] pass 2 (medium) W: ", paste(pass2_W, collapse = ","), "\n")
cat("[D17mp] pass 3 (small)  W: ", paste(pass3_W, collapse = ","), "\n")
cat("[D17mp] sim thresholds: min=", sim_thresh_min, " max=", sim_thresh_max,
    " tau=", sim_thresh_tau, "  density_min=", density_min, "\n", sep = "")
cat("[D17mp] composition_min=", comp_min, "\n")
cat("[D17mp] trim_core=", trim_core, " trim_drop_z=", trim_drop_z, "\n", sep = "")
cat("[D17mp] post_trim_nms=", post_trim_nms,
    " post_trim_nms_overlap=", post_trim_nms_overlap, "\n", sep = "")
cat("[D17mp] z_validate=", z_validate,
    " z_merge_thresh=", z_merge_thresh,
    " z_split_thresh=", z_split_thresh, "\n", sep = "")
cat("[D17mp] architecture=", architecture,
    "  l1_density_min=", l1_density_min,
    "  l3_enabled=", l3_enabled, "\n", sep = "")
cat("[D17mp] l2_method=", l2_method, "\n", sep = "")
if (l2_method == "refiner") {
  cat("[D17mp]   refiner_W=", paste(refiner_W_vec, collapse = ","),
      "  step_frac=", refiner_step_frac,
      "  score_min=", refiner_score_min,
      "  min_scales=", refiner_min_scales,
      "  max_depth=", refiner_max_depth,
      "  min_block=", refiner_min_block,
      "  cross_term=", refiner_cross_term, sep = "")
  if (refiner_cross_term == "zq") {
    cat("  z_quantile=", refiner_z_quantile, sep = "")
  }
  cat("\n")
}
cat("[D17mp] internal clustering: cluster_min_W=", cluster_min_W,
    "  k=", k_min, "..", k_max, "\n", sep = "")

# ---- Adaptive size-dependent threshold function ----------------------------

size_threshold <- function(W) {
  sim_thresh_min + (sim_thresh_max - sim_thresh_min) * exp(-W / sim_thresh_tau)
}

cat("[D17mp] adaptive thresholds by size:\n")
for (W in sort(unique(c(pass1_W, pass2_W, pass3_W)), decreasing = TRUE)) {
  cat(sprintf("[D17mp]   W=%4d  threshold=%.3f\n", W, size_threshold(W)))
}

# ---- Load precomp -----------------------------------------------------------

cat("[D17mp] loading precomp\n")
pc_obj <- readRDS(precomp_f)
if (!is.null(pc_obj$dt)) {
  pc <- pc_obj
} else if (!is.null(pc_obj$pc)) {
  pc <- pc_obj$pc
} else {
  stop("[D17mp] cannot find pc$dt in precomp file")
}

dt_pc <- as.data.table(pc$dt)
n_windows_total <- nrow(dt_pc)

pc1_cols <- grep("^PC_1_", names(dt_pc), value = TRUE)
sample_ids <- sub("^PC_1_", "", pc1_cols)
n_samples <- length(sample_ids)

# ---- Normalize window coordinate columns -----------------------------------
# Older/newer precomp objects may store window coordinates as start/end,
# start_bp/end_bp, or window_start/window_end. D17 needs bp coordinates here.
if (all(c("start", "end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start
  window_end_bp   <- dt_pc$end
} else if (all(c("start_bp", "end_bp") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$start_bp
  window_end_bp   <- dt_pc$end_bp
  dt_pc[, `:=`(start = start_bp, end = end_bp)]
} else if (all(c("window_start", "window_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$window_start
  window_end_bp   <- dt_pc$window_end
  dt_pc[, `:=`(start = window_start, end = window_end)]
} else if (all(c("bp_start", "bp_end") %in% names(dt_pc))) {
  window_start_bp <- dt_pc$bp_start
  window_end_bp   <- dt_pc$bp_end
  dt_pc[, `:=`(start = bp_start, end = bp_end)]
} else {
  stop("[D17mp] precomp dt missing usable coordinate columns. Found: ",
       paste(names(dt_pc), collapse = ", "))
}

PC1 <- as.matrix(dt_pc[, ..pc1_cols])
storage.mode(PC1) <- "double"

cat("[D17mp] N windows: ", n_windows_total, " | N samples: ", n_samples, "\n", sep = "")

# ---- Load sim_mat -----------------------------------------------------------

if (is.na(sim_mat_f)) {
  candidates_path <- c(
    file.path(dirname(precomp_f), "sim_mat_nn80.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn160.rds"),
    file.path(dirname(precomp_f), "sim_mat_nn40.rds")
  )
  hit <- candidates_path[file.exists(candidates_path)]
  if (length(hit) == 0)
    stop("[D17mp] --sim_mat not given and could not auto-find sim_mat_nn*.rds")
  sim_mat_f <- hit[1]
  cat("[D17mp] auto-found sim_mat: ", sim_mat_f, "\n")
}

cat("[D17mp] loading sim_mat\n")
sm_obj <- readRDS(sim_mat_f)
if (is.matrix(sm_obj)) {
  sim_mat <- sm_obj
} else if (is.list(sm_obj) && !is.null(sm_obj$sim_mat) && is.matrix(sm_obj$sim_mat)) {
  sim_mat <- sm_obj$sim_mat
} else if (is.list(sm_obj) && length(sm_obj) == 1 && is.matrix(sm_obj[[1]])) {
  sim_mat <- sm_obj[[1]]
} else {
  stop("[D17mp] sim_mat object structure not recognized; class=", class(sm_obj)[1])
}
storage.mode(sim_mat) <- "double"

if (!isTRUE(nrow(sim_mat) == n_windows_total &&
            ncol(sim_mat) == n_windows_total)) {
  stop("[D17mp] sim_mat dim ", nrow(sim_mat), "x", ncol(sim_mat),
       " does not match precomp N windows ", n_windows_total)
}

# ---- Diagonal-distance Z precompute (for z_validate) -----------------------
# For each diagonal distance d (0 .. N-1), we need mean and SD of sim_mat
# values along that diagonal. Then any pixel sim_mat[i, j] has
#   Z[i, j] = (sim_mat[i, j] - diag_mean[|i - j|]) / diag_sd[|i - j|]
# Because the sim_mat is symmetric and diagonal-distance Z is well defined
# only when there are enough off-diagonal samples per d, we compute it
# once at startup and only when z_validate is requested.
diag_mean <- NULL
diag_sd   <- NULL
need_diag_z <- z_validate ||
               (l2_method == "refiner" &&
                refiner_cross_term %in% c("z", "zq"))
if (need_diag_z) {
  cat("[D17mp] precomputing diagonal-distance Z (mean + sd per d)\n")
  diag_mean <- numeric(n_windows_total)
  diag_sd   <- numeric(n_windows_total)
  for (d in 0:(n_windows_total - 1L)) {
    if (d == 0L) {
      vals <- diag(sim_mat)
    } else {
      ii <- seq.int(1L, n_windows_total - d)
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
    if (d > 0L && d %% 1000L == 0L) {
      cat("[D17mp]   diag d=", d, " / ", n_windows_total - 1L, "\n", sep = "")
    }
  }
  ok_d <- sum(!is.na(diag_sd))
  cat("[D17mp]   diag_mean/diag_sd built; usable distances: ", ok_d,
      " / ", n_windows_total, "\n", sep = "")
}

# Mean cross-block Z between two adjacent ranges [s1, e1] and [s2, e2].
# Returns NA if any of the contributing diag distances are missing.
crossblock_z <- function(s1, e1, s2, e2) {
  if (is.null(diag_mean) || is.null(diag_sd)) return(NA_real_)
  if (s2 < s1) {
    tmp <- c(s1, e1); s1 <- s2; e1 <- e2; s2 <- tmp[1]; e2 <- tmp[2]
  }
  ii <- seq.int(s1, e1)
  jj <- seq.int(s2, e2)
  # Build cross-block submatrix
  block <- sim_mat[ii, jj, drop = FALSE]
  # Distance of each cell from main diagonal: jj - ii (always positive here)
  d_mat <- outer(ii, jj, FUN = function(a, b) b - a)
  # Index into diag_mean/sd is d + 1 (because d=0 is index 1)
  idx <- d_mat + 1L
  mu  <- diag_mean[idx]
  sg  <- diag_sd[idx]
  z   <- (block - mu) / sg
  z[!is.finite(z)] <- NA_real_
  if (all(is.na(z))) return(NA_real_)
  mean(z, na.rm = TRUE)
}

# Variant that returns the full Z matrix instead of the mean. Used by the
# boundary refiner with cross_term="zq" to compute a lower-tail quantile.
crossblock_z_matrix <- function(s1, e1, s2, e2) {
  if (is.null(diag_mean) || is.null(diag_sd)) return(NULL)
  if (s2 < s1) {
    tmp <- c(s1, e1); s1 <- s2; e1 <- e2; s2 <- tmp[1]; e2 <- tmp[2]
  }
  ii <- seq.int(s1, e1)
  jj <- seq.int(s2, e2)
  block <- sim_mat[ii, jj, drop = FALSE]
  d_mat <- outer(ii, jj, FUN = function(a, b) b - a)
  idx <- d_mat + 1L
  mu  <- diag_mean[idx]
  sg  <- diag_sd[idx]
  z   <- (block - mu) / sg
  z[!is.finite(z)] <- NA_real_
  z
}


# ---- Optional candidate edge trimming ---------------------------------------
# For a candidate interval s:e, keep the dense internal core using row-mean
# similarity inside the candidate submatrix. Weak edge rows are dropped.
trim_to_core <- function(s, e, drop_z = 1.0) {
  block <- sim_mat[s:e, s:e]
  row_means <- rowMeans(block, na.rm = TRUE)

  mu <- mean(row_means, na.rm = TRUE)
  sg <- sd(row_means, na.rm = TRUE)

  if (!is.finite(mu) || !is.finite(sg) || sg < 1e-9) {
    return(c(s, e))
  }

  thr <- mu - drop_z * sg
  keep <- which(row_means >= thr)

  if (length(keep) < min_n_windows) {
    return(c(s, e))   # too aggressive, bail
  }

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

  # Recompute block metrics after trimming.
  for (i in seq_len(nrow(dt))) {
    block <- sim_mat[dt$start_w[i]:dt$end_w[i], dt$start_w[i]:dt$end_w[i]]
    dt$mean_sim[i] <- mean(block, na.rm = TRUE)
    dt$density_p70[i] <- mean(block >= 0.70, na.rm = TRUE)
    dt$density_adaptive[i] <- mean(block >= dt$threshold[i], na.rm = TRUE)
  }

  dt
}

# ---- Helper: per-window k=3 k-means band assignment -------------------------

per_window_kmeans_band <- function(pc1_vec, k = 3L) {
  if (any(!is.finite(pc1_vec))) pc1_vec[!is.finite(pc1_vec)] <- 0
  if (sd(pc1_vec) < 1e-9) return(rep(1L, length(pc1_vec)))
  km <- tryCatch(
    suppressWarnings(kmeans(pc1_vec, centers = k, nstart = 5L, iter.max = 30L)),
    error = function(e) NULL
  )
  if (is.null(km)) return(rep(1L, length(pc1_vec)))
  ord <- order(km$centers[, 1])
  band_remap <- match(seq_len(k), ord)
  band_remap[km$cluster]
}

# ---- Composition consistency for a window range ----------------------------

composition_consistency <- function(start_w, end_w) {
  ws <- seq.int(start_w, end_w)
  if (length(ws) < 5L) {
    return(list(consistency = NA_real_, n_inv = NA_integer_))
  }
  band_mat <- matrix(NA_integer_, nrow = n_samples, ncol = length(ws))
  for (j in seq_along(ws)) {
    band_mat[, j] <- per_window_kmeans_band(PC1[ws[j], ], k = comp_kmeans_k)
  }
  modal_frac <- vapply(seq_len(n_samples), function(s) {
    tab <- tabulate(band_mat[s, ], nbins = comp_kmeans_k)
    max(tab) / length(ws)
  }, numeric(1))
  consistency <- mean(modal_frac, na.rm = TRUE)
  mid <- ws[length(ws) %/% 2L + 1L]
  mid_bands <- per_window_kmeans_band(PC1[mid, ], k = comp_kmeans_k)
  band_sizes <- tabulate(mid_bands, nbins = comp_kmeans_k)
  if (sum(band_sizes > 0) == 0L) {
    n_inv <- NA_integer_
  } else {
    n_inv <- min(band_sizes[band_sizes > 0])
  }
  list(consistency = consistency, n_inv = as.integer(n_inv))
}

# ---- Single-pass scanner with masking ---------------------------------------
# Scans the diagonal at the given W sizes, but only at positions where the
# entire WxW block is in unmasked territory (no overlap with any masked range).

scan_pass <- function(W_sizes, masked_runs, pass_index, pass_color_label) {
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

      # Skip if this block overlaps any masked run
      blocked <- FALSE
      if (length(masked_runs) > 0L) {
        for (mr in masked_runs) {
          if (!(e < mr[1] || s > mr[2])) {
            blocked <- TRUE
            break
          }
        }
      }
      if (blocked) next

      block <- sim_mat[s:e, s:e]
      msim <- mean(block, na.rm = TRUE)
      if (!is.finite(msim) || msim < thr_W) next

      # Adaptive density: pixels above the size-adaptive threshold itself
      dens_adaptive <- mean(block >= thr_W, na.rm = TRUE)
      if (!is.finite(dens_adaptive) || dens_adaptive < density_min) next
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

# ---- Region-scoped scanner --------------------------------------------------
# Scans for blocks of the given W sizes whose entire WxW range lies inside
# [bound_s, bound_e]. Used by Level 2 (inside L1 envelopes) and Level 3
# (inside L2 inversions). NO masking — all positions inside the bounds are
# eligible. Pass index is recorded for downstream colouring.
scan_region <- function(W_sizes, bound_s, bound_e, pass_index) {
  all_hits <- list()
  region_w <- bound_e - bound_s + 1L
  for (W in W_sizes) {
    if (W > region_w) next
    thr_W <- size_threshold(W)
    half <- W %/% 2L
    # Centers must keep the WxW block fully inside [bound_s, bound_e]
    c_lo <- bound_s + half
    c_hi <- bound_e - half + (if (W %% 2L == 0L) 1L else 0L)
    if (c_hi < c_lo) next
    centers <- seq.int(c_lo, c_hi, by = scan_step)
    if (length(centers) == 0L) next

    out <- vector("list", length(centers))
    out_i <- 1L
    for (ctr in centers) {
      s <- ctr - half
      e <- s + W - 1L
      if (s < bound_s || e > bound_e) next

      block <- sim_mat[s:e, s:e]
      msim <- mean(block, na.rm = TRUE)
      if (!is.finite(msim) || msim < thr_W) next

      dens_adaptive <- mean(block >= thr_W, na.rm = TRUE)
      if (!is.finite(dens_adaptive) || dens_adaptive < density_min) next
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

# ---- Within-scale collapse -------------------------------------------------

collapse_within_scale <- function(dt) {
  if (nrow(dt) == 0L) return(dt)
  dt <- dt[order(-mean_sim)]      # best-first, not start_w-ordered
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

# ---- Size-biased NMS within a pass -----------------------------------------

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
        if (o / smaller > nms_overlap) {
          is_dominated <- TRUE
          break
        }
      }
    }
    if (!is_dominated) kept[[length(kept) + 1L]] <- row_i
  }
  rbindlist(kept)[order(start_w)]
}

# Parametric NMS used after trim_to_core. Sort by mean_sim (best first); a
# candidate is suppressed if any already-kept candidate overlaps it by more
# than `overlap_thr` of the smaller of the two. Uses end - start + 1 as the
# size, so it respects post-trim widths rather than the original W.
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

# ---- Run one pass: scan + collapse + NMS + composition + mask ---------------

run_pass <- function(W_sizes, masked_runs, pass_index, pass_label) {
  cat("\n[D17mp] === PASS ", pass_index, " (", pass_label,
      ") | W sizes: ", paste(W_sizes, collapse = ","), " ===\n", sep = "")

  # Scan
  hits <- scan_pass(W_sizes, masked_runs, pass_index, pass_label)
  cat("[D17mp]   scan hits: ", nrow(hits), "\n", sep = "")
  if (nrow(hits) == 0L) {
    return(list(survivors = data.table(), masked_runs = masked_runs))
  }
  cat("[D17mp]   by W: "); print(table(hits$W))

  # Within-scale collapse
  by_W <- split(hits, hits$W)
  collapsed <- rbindlist(lapply(by_W, collapse_within_scale))
  cat("[D17mp]   after within-scale collapse: ", nrow(collapsed), "\n", sep = "")

  # Size-biased NMS
  nms_dt <- size_biased_nms(collapsed)
  cat("[D17mp]   after size-biased NMS: ", nrow(nms_dt), "\n", sep = "")

  if (nrow(nms_dt) == 0L) {
    return(list(survivors = data.table(), masked_runs = masked_runs))
  }

  # Optional core trimming before composition validation.
  # This tightens broad envelopes before testing whether PC1 bands are stable.
  if (trim_core) {
    nms_dt <- apply_core_trimming(nms_dt)
    cat("[D17mp]   after trim_to_core: median dropped windows = ",
        median(nms_dt$trim_total_dropped, na.rm = TRUE),
        " | max dropped = ",
        max(nms_dt$trim_total_dropped, na.rm = TRUE), "\n", sep = "")
  }

  # Optional second NMS after trim. Trim refines edges of broad scan envelopes
  # — afterward, sibling parents that had > nms_overlap pre-trim may sit right
  # next to each other on the same biological signal. This second pass with a
  # tighter overlap threshold (default 0.30) collapses those duplicates.
  if (post_trim_nms && nrow(nms_dt) > 1L) {
    n_before <- nrow(nms_dt)
    nms_dt <- size_biased_nms_param(nms_dt, post_trim_nms_overlap)
    cat("[D17mp]   after post-trim NMS (overlap>", post_trim_nms_overlap,
        "): ", n_before, " -> ", nrow(nms_dt), "\n", sep = "")
  }

  if (nrow(nms_dt) == 0L) {
    return(list(survivors = data.table(), masked_runs = masked_runs))
  }

  # Composition validation
  cat("[D17mp]   running composition validation on ", nrow(nms_dt),
      " candidates\n", sep = "")
  comp_consistencies <- numeric(nrow(nms_dt))
  n_inv_carriers <- integer(nrow(nms_dt))
  for (i in seq_len(nrow(nms_dt))) {
    res <- composition_consistency(nms_dt$start_w[i], nms_dt$end_w[i])
    comp_consistencies[i] <- res$consistency
    n_inv_carriers[i] <- res$n_inv
  }
  nms_dt[, composition_consistency := comp_consistencies]
  nms_dt[, n_inv_carriers := n_inv_carriers]

  passed <- nms_dt[!is.na(composition_consistency) &
                   composition_consistency >= comp_min]
  failed <- nms_dt[is.na(composition_consistency) |
                   composition_consistency < comp_min]

  cat("[D17mp]   composition PASSED: ", nrow(passed),
      " | FAILED (drop, do NOT mask): ", nrow(failed), "\n", sep = "")
  if (nrow(passed) > 0L) {
    cat(sprintf("[D17mp]     PASS comp range: %.3f - %.3f (median %.3f)\n",
                min(passed$composition_consistency),
                max(passed$composition_consistency),
                median(passed$composition_consistency)))
  }

  # Update masked_runs: every PASSED candidate's window range gets masked
  new_masked <- masked_runs
  if (nrow(passed) > 0L) {
    for (i in seq_len(nrow(passed))) {
      new_masked[[length(new_masked) + 1L]] <-
        c(passed$start_w[i], passed$end_w[i])
    }
  }

  list(survivors = passed, masked_runs = new_masked)
}

# ---- Nested pipeline (Level 1 envelopes -> Level 2 inversions) -------------
# In the "nested" architecture:
#   * Level 1 (W=800/400) finds large envelopes — regions of elevated
#     similarity that may contain one or more inversions. NO composition
#     gate (envelopes are allowed to be heterogeneous, that is the whole
#     point of looking inside them) and NO masking against subsequent
#     levels (we WANT subsequent levels to look inside envelopes).
#   * Level 2 (W=200/100) scans inside each envelope for the actual
#     inversions. Composition gate applies here.
#   * Optional Level 3 (W=50) scans inside each Level 2 inversion for
#     small sub-features.
# Ward.D2 runs later on Level 2 inversions (depth 1 -> Ward children depth 2).
# Level 1 envelopes are kept in the catalogue for context but are not
# clustered themselves.

scan_collapse_nms_trim <- function(W_sizes, scan_fun, label, gate_composition) {
  # Helper that does the common: scan -> collapse -> NMS -> trim -> post-NMS
  # -> optional composition. Returns the survivor data.table (possibly empty).
  hits <- scan_fun(W_sizes)
  cat("[D17mp]   [", label, "] scan hits: ", nrow(hits), "\n", sep = "")
  if (nrow(hits) == 0L) return(data.table())
  cat("[D17mp]   [", label, "] by W: ", sep = ""); print(table(hits$W))

  by_W <- split(hits, hits$W)
  collapsed <- rbindlist(lapply(by_W, collapse_within_scale))
  cat("[D17mp]   [", label, "] after within-scale collapse: ",
      nrow(collapsed), "\n", sep = "")

  nms_dt <- size_biased_nms(collapsed)
  cat("[D17mp]   [", label, "] after size-biased NMS: ",
      nrow(nms_dt), "\n", sep = "")
  if (nrow(nms_dt) == 0L) return(data.table())

  if (trim_core) {
    nms_dt <- apply_core_trimming(nms_dt)
    cat("[D17mp]   [", label, "] after trim_to_core: median dropped = ",
        median(nms_dt$trim_total_dropped, na.rm = TRUE),
        " | max = ", max(nms_dt$trim_total_dropped, na.rm = TRUE),
        "\n", sep = "")
  }

  if (post_trim_nms && nrow(nms_dt) > 1L) {
    n_before <- nrow(nms_dt)
    nms_dt <- size_biased_nms_param(nms_dt, post_trim_nms_overlap)
    cat("[D17mp]   [", label, "] after post-trim NMS (overlap>",
        post_trim_nms_overlap, "): ", n_before, " -> ",
        nrow(nms_dt), "\n", sep = "")
  }

  if (nrow(nms_dt) == 0L) return(data.table())

  if (gate_composition) {
    comp_consistencies <- numeric(nrow(nms_dt))
    n_inv_carriers     <- integer(nrow(nms_dt))
    for (i in seq_len(nrow(nms_dt))) {
      res <- composition_consistency(nms_dt$start_w[i], nms_dt$end_w[i])
      comp_consistencies[i] <- res$consistency
      n_inv_carriers[i]     <- res$n_inv
    }
    nms_dt[, composition_consistency := comp_consistencies]
    nms_dt[, n_inv_carriers           := n_inv_carriers]
    passed <- nms_dt[!is.na(composition_consistency) &
                     composition_consistency >= comp_min]
    failed <- nms_dt[is.na(composition_consistency) |
                     composition_consistency < comp_min]
    cat("[D17mp]   [", label, "] composition PASSED: ", nrow(passed),
        " | FAILED: ", nrow(failed), "\n", sep = "")
    if (nrow(passed) > 0L) {
      cat(sprintf("[D17mp]     [%s] PASS comp range: %.3f - %.3f (median %.3f)\n",
                  label,
                  min(passed$composition_consistency),
                  max(passed$composition_consistency),
                  median(passed$composition_consistency)))
    }
    return(passed)
  } else {
    # No composition gate — fill with NA and pass everything
    nms_dt[, composition_consistency := NA_real_]
    nms_dt[, n_inv_carriers           := NA_integer_]
    return(nms_dt)
  }
}

# ---- Adaptive boundary refiner ---------------------------------------------
# Replacement for fixed-W L2 scan inside an envelope. For each candidate
# boundary b at scale W, we score:
#   left_mean  = mean(sim_mat[(b-W+1):b, (b-W+1):b])
#   right_mean = mean(sim_mat[(b+1):(b+W), (b+1):(b+W)])
#   For cross_term="z":  cross = -mean(diag-Z over sim_mat[left x right])
#   For cross_term="raw": cross = -mean(sim_mat[left x right])
#   boundary_score = (left_mean + right_mean)/2 * (-cross)   (z form)
#   boundary_score = (left_mean + right_mean)/2 - mean(cross) (raw form)
# A real boundary has high coherence on both sides AND strong negative
# diag-Z in the cross-block. Boundaries scoring above refiner_score_min at
# >= refiner_min_scales distinct W values become split points; refinement
# recurses inside each block up to refiner_max_depth.

local_boundary_score <- function(b, W, s_min, e_max, cross_term = "z") {
  l1 <- b - W + 1L; l2 <- b
  r1 <- b + 1L;     r2 <- b + W
  if (l1 < s_min || r2 > e_max) return(NA_real_)
  left_mean  <- mean(sim_mat[l1:l2, l1:l2], na.rm = TRUE)
  right_mean <- mean(sim_mat[r1:r2, r1:r2], na.rm = TRUE)
  if (!is.finite(left_mean) || !is.finite(right_mean)) return(NA_real_)
  if (cross_term == "z" && !is.null(diag_mean) && !is.null(diag_sd)) {
    cz <- crossblock_z(l1, l2, r1, r2)
    if (!is.finite(cz)) return(NA_real_)
    return(((left_mean + right_mean) / 2) * (-cz))
  } else if (cross_term == "zq" && !is.null(diag_mean) && !is.null(diag_sd)) {
    # Lower-tail quantile of cross-block Z. Robust to noisy single pixels;
    # captures "blue valleys" — boundary positions where a substantial
    # fraction of the cross-block sits well below the diag baseline.
    zmat <- crossblock_z_matrix(l1, l2, r1, r2)
    if (is.null(zmat)) return(NA_real_)
    vals <- as.numeric(zmat)
    vals <- vals[is.finite(vals)]
    if (length(vals) < 5L) return(NA_real_)
    zq <- as.numeric(quantile(vals, refiner_z_quantile, na.rm = TRUE))
    if (!is.finite(zq)) return(NA_real_)
    # Score = -zq, so blue valleys (large negative zq) -> large positive score.
    # We do NOT multiply by coherence here because the quantile itself already
    # encodes "this is below baseline" — multiplying by a near-1 coherence
    # would just rescale uniformly; the signal is in the quantile alone.
    return(-zq)
  } else {
    cross_mean <- mean(sim_mat[l1:l2, r1:r2], na.rm = TRUE)
    if (!is.finite(cross_mean)) return(NA_real_)
    return(((left_mean + right_mean) / 2) - cross_mean)
  }
}

scan_boundaries_one_W <- function(s, e, W, step, cross_term) {
  if ((e - s + 1L) < 3L * W) {
    return(data.table(boundary_w = integer(0), W = integer(0),
                      score = numeric(0)))
  }
  b_vec <- seq.int(s + W, e - W, by = step)
  if (length(b_vec) == 0L) {
    return(data.table(boundary_w = integer(0), W = integer(0),
                      score = numeric(0)))
  }
  scores <- vapply(b_vec, function(b)
    local_boundary_score(b, W, s, e, cross_term), numeric(1))
  out <- data.table(boundary_w = b_vec, W = W, score = scores)
  out[is.finite(score)]
}

scan_boundaries_multiW <- function(s, e, W_vec, step_frac, cross_term) {
  all <- list()
  for (W in W_vec) {
    if ((e - s + 1L) < 3L * W) next
    step <- max(1L, as.integer(W * step_frac))
    dt <- scan_boundaries_one_W(s, e, W, step, cross_term)
    if (nrow(dt) > 0L) all[[as.character(W)]] <- dt
  }
  if (length(all) == 0L) return(data.table())
  rbindlist(all)
}

pick_stable_boundaries <- function(bound_dt, score_min, merge_dist,
                                   min_scales) {
  if (nrow(bound_dt) == 0L) return(data.table())
  peaks <- bound_dt[score >= score_min]
  if (nrow(peaks) == 0L) return(data.table())
  setorder(peaks, boundary_w)
  groups <- list()
  current <- peaks[1L]
  if (nrow(peaks) >= 2L) {
    for (i in 2:nrow(peaks)) {
      row <- peaks[i]
      anchor <- median(current$boundary_w)
      if (abs(row$boundary_w - anchor) <= merge_dist) {
        current <- rbindlist(list(current, row), fill = TRUE)
      } else {
        groups[[length(groups) + 1L]] <- current
        current <- row
      }
    }
  }
  groups[[length(groups) + 1L]] <- current
  out <- lapply(groups, function(g) {
    data.table(
      boundary_w  = as.integer(round(weighted.mean(g$boundary_w,
                                                   pmax(g$score, 1e-6)))),
      n_scales    = uniqueN(g$W),
      W_support   = paste(sort(unique(g$W)), collapse = ","),
      max_score   = max(g$score, na.rm = TRUE),
      mean_score  = mean(g$score, na.rm = TRUE),
      n_peak_rows = nrow(g)
    )
  })
  out <- rbindlist(out)
  out[n_scales >= min_scales][order(boundary_w)]
}

boundaries_to_blocks <- function(s, e, boundaries, min_block) {
  if (nrow(boundaries) == 0L) {
    return(data.table(start_w = s, end_w = e,
                      n_windows = e - s + 1L))
  }
  b <- sort(unique(boundaries$boundary_w))
  starts <- c(s, b + 1L)
  ends   <- c(b, e)
  blocks <- data.table(start_w = starts, end_w = ends)
  blocks[, n_windows := end_w - start_w + 1L]
  blocks[n_windows >= min_block]
}

refine_region_recursive <- function(s, e, depth, max_depth,
                                    W_vec, step_frac, score_min,
                                    min_scales, min_block, cross_term,
                                    audit_records) {
  if (depth >= max_depth || (e - s + 1L) < 2L * min_block) {
    leaf <- data.table(start_w = s, end_w = e,
                       n_windows = e - s + 1L,
                       refine_depth = depth)
    return(list(leaves = list(leaf), audit = audit_records))
  }
  W_local <- W_vec[(e - s + 1L) >= 3L * W_vec]
  if (length(W_local) == 0L) {
    leaf <- data.table(start_w = s, end_w = e,
                       n_windows = e - s + 1L,
                       refine_depth = depth)
    return(list(leaves = list(leaf), audit = audit_records))
  }
  bd_raw <- scan_boundaries_multiW(s, e, W_local, step_frac, cross_term)
  merge_dist <- max(30L, as.integer(0.05 * (e - s + 1L)))
  bd <- pick_stable_boundaries(bd_raw, score_min, merge_dist, min_scales)
  if (nrow(bd_raw) > 0L) {
    audit_records[[length(audit_records) + 1L]] <-
      copy(bd_raw)[, `:=`(region_s = s, region_e = e, depth = depth,
                          picked = FALSE)]
  }
  if (nrow(bd) > 0L) {
    audit_records[[length(audit_records) + 1L]] <-
      copy(bd)[, `:=`(region_s = s, region_e = e, depth = depth,
                      picked = TRUE)]
  }
  child_blocks <- boundaries_to_blocks(s, e, bd, min_block)
  if (nrow(child_blocks) <= 1L) {
    leaf <- data.table(start_w = s, end_w = e,
                       n_windows = e - s + 1L,
                       refine_depth = depth)
    return(list(leaves = list(leaf), audit = audit_records))
  }
  all_leaves <- list()
  for (i in seq_len(nrow(child_blocks))) {
    res <- refine_region_recursive(
      s = child_blocks$start_w[i], e = child_blocks$end_w[i],
      depth = depth + 1L, max_depth = max_depth,
      W_vec = W_local, step_frac = step_frac, score_min = score_min,
      min_scales = min_scales, min_block = min_block,
      cross_term = cross_term, audit_records = audit_records)
    all_leaves <- c(all_leaves, res$leaves)
    audit_records <- res$audit
  }
  list(leaves = all_leaves, audit = audit_records)
}

refine_envelope <- function(env_s, env_e) {
  audit <- list()
  res <- refine_region_recursive(
    s = env_s, e = env_e, depth = 0L,
    max_depth = refiner_max_depth,
    W_vec = refiner_W_vec, step_frac = refiner_step_frac,
    score_min = refiner_score_min,
    min_scales = refiner_min_scales,
    min_block = refiner_min_block,
    cross_term = refiner_cross_term,
    audit_records = audit)
  leaves <- if (length(res$leaves) > 0L)
    rbindlist(res$leaves, use.names = TRUE, fill = TRUE) else data.table()
  list(blocks = leaves,
       audit  = if (length(res$audit) > 0L)
         rbindlist(res$audit, use.names = TRUE, fill = TRUE) else data.table())
}

# Decorate refiner leaf blocks to look like scan output rows so they can
# flow through trim + post-NMS + composition unchanged.
decorate_refiner_blocks <- function(leaves, pass_index) {
  if (nrow(leaves) == 0L) return(data.table())
  out <- data.table(
    W = leaves$n_windows,
    start_w = leaves$start_w,
    end_w = leaves$end_w,
    n_windows = leaves$n_windows,
    discovery_pass = pass_index
  )
  ms <- numeric(nrow(out)); dp <- numeric(nrow(out)); da <- numeric(nrow(out))
  for (i in seq_len(nrow(out))) {
    b <- sim_mat[out$start_w[i]:out$end_w[i],
                 out$start_w[i]:out$end_w[i]]
    ms[i] <- mean(b, na.rm = TRUE)
    dp[i] <- mean(b >= 0.70, na.rm = TRUE)
    da[i] <- mean(b >= 0.55, na.rm = TRUE)
  }
  out[, mean_sim := ms]
  out[, density_p70 := dp]
  out[, density_adaptive := da]
  out[, threshold := NA_real_]
  out[, refine_depth := leaves$refine_depth]
  out
}

run_nested_pipeline <- function() {
  # ---- Level 1: chromosome-wide envelope discovery ------------------------
  cat("\n[D17mp] === LEVEL 1 (envelopes) | W sizes: ",
      paste(pass1_W, collapse = ","), " ===\n", sep = "")
  # L1 uses a relaxed density floor (l1_density_min instead of density_min)
  # so envelopes can be heterogeneous regions of elevated similarity, not
  # just dense inversion blocks. Achieved by temporarily swapping the
  # density_min global for the duration of the L1 scan.
  l1_scan <- function(W_sizes) {
    saved_density <- density_min
    on.exit(density_min <<- saved_density, add = TRUE)
    density_min <<- l1_density_min
    scan_pass(W_sizes, list(), 1L, "envelope")
  }
  envelopes <- scan_collapse_nms_trim(pass1_W, l1_scan, "L1 envelope",
                                      gate_composition = FALSE)

  if (nrow(envelopes) == 0L) {
    cat("[D17mp] no Level 1 envelopes found\n")
    return(list(envelopes = data.table(), inversions = data.table()))
  }
  envelopes[, level := 1L]
  envelopes[, recursion_depth := 0L]

  cat("[D17mp] L1: ", nrow(envelopes), " envelopes (totalling ",
      sum(envelopes$end_w - envelopes$start_w + 1L), " windows = ",
      round(100 * sum(envelopes$end_w - envelopes$start_w + 1L) /
            n_windows_total, 1L), "% of chromosome)\n", sep = "")

  # ---- Level 2: per-envelope inversion discovery --------------------------
  cat("\n[D17mp] === LEVEL 2 (", l2_method,
      ") | W sizes: ", paste(pass2_W, collapse = ","), " ===\n", sep = "")
  l2_inversions <- list()
  refiner_audit_all <- list()
  for (e_idx in seq_len(nrow(envelopes))) {
    env_row <- envelopes[e_idx]
    bs <- env_row$start_w; be <- env_row$end_w
    cat("[D17mp]   envelope #", e_idx, " w=[", bs, ", ", be,
        "] (", sprintf("%.2f-%.2f Mb", window_start_bp[bs] / 1e6,
                       window_end_bp[be] / 1e6), ")\n", sep = "")
    if (l2_method == "scan") {
      l2_scan <- function(W_sizes) scan_region(W_sizes, bs, be, 2L)
      inv_dt <- scan_collapse_nms_trim(pass2_W, l2_scan,
                                       paste0("L2 env#", e_idx),
                                       gate_composition = TRUE)
    } else {
      # refiner path: produce leaf blocks, decorate, then run trim +
      # composition through the same helper. Skip NMS (refiner blocks
      # are non-overlapping by construction) by using a single-W sized
      # input that won't collapse anything.
      ref <- refine_envelope(bs, be)
      cat("[D17mp]   [L2 env#", e_idx,
          " refiner] leaf blocks: ", nrow(ref$blocks), "\n", sep = "")
      if (nrow(ref$audit) > 0L) {
        ref$audit[, envelope_idx := e_idx]
        refiner_audit_all[[length(refiner_audit_all) + 1L]] <- ref$audit
      }
      if (nrow(ref$blocks) == 0L) {
        inv_dt <- data.table()
      } else {
        deco <- decorate_refiner_blocks(ref$blocks, 2L)
        # Apply trim + composition the same way scan_collapse_nms_trim does,
        # but skip the within-scale collapse / size-biased NMS steps since
        # refiner leaves are non-overlapping.
        if (trim_core) {
          deco <- apply_core_trimming(deco)
          cat("[D17mp]   [L2 env#", e_idx,
              " refiner] after trim_to_core: median dropped = ",
              median(deco$trim_total_dropped, na.rm = TRUE),
              " | max = ", max(deco$trim_total_dropped, na.rm = TRUE),
              "\n", sep = "")
        }
        # Composition gate
        comp_consistencies <- numeric(nrow(deco))
        n_inv_carriers     <- integer(nrow(deco))
        for (i in seq_len(nrow(deco))) {
          rs <- composition_consistency(deco$start_w[i], deco$end_w[i])
          comp_consistencies[i] <- rs$consistency
          n_inv_carriers[i]     <- rs$n_inv
        }
        deco[, composition_consistency := comp_consistencies]
        deco[, n_inv_carriers           := n_inv_carriers]
        passed <- deco[!is.na(composition_consistency) &
                       composition_consistency >= comp_min]
        failed <- deco[is.na(composition_consistency) |
                       composition_consistency < comp_min]
        cat("[D17mp]   [L2 env#", e_idx,
            " refiner] composition PASSED: ", nrow(passed),
            " | FAILED: ", nrow(failed), "\n", sep = "")
        if (nrow(passed) > 0L) {
          cat(sprintf("[D17mp]     [L2 env#%d refiner] PASS comp range: %.3f - %.3f (median %.3f)\n",
                      e_idx,
                      min(passed$composition_consistency),
                      max(passed$composition_consistency),
                      median(passed$composition_consistency)))
        }
        inv_dt <- passed
      }
    }
    if (nrow(inv_dt) > 0L) {
      inv_dt[, level := 2L]
      inv_dt[, recursion_depth := 1L]
      inv_dt[, parent_envelope_idx := e_idx]
      l2_inversions[[length(l2_inversions) + 1L]] <- inv_dt
    }
  }
  inversions <- if (length(l2_inversions) > 0L)
    rbindlist(l2_inversions, use.names = TRUE, fill = TRUE) else data.table()
  cat("[D17mp] L2: ", nrow(inversions), " inversions across all envelopes\n",
      sep = "")

  # Stash refiner audit on a global so the writer can dump it later.
  if (l2_method == "refiner" && length(refiner_audit_all) > 0L) {
    refiner_audit_global <<-
      rbindlist(refiner_audit_all, use.names = TRUE, fill = TRUE)
  }

  # ---- Level 3 (optional): inside L2 inversions ---------------------------
  if (l3_enabled && nrow(inversions) > 0L) {
    cat("\n[D17mp] === LEVEL 3 (sub-features) | W sizes: ",
        paste(pass3_W, collapse = ","), " ===\n", sep = "")
    l3_features <- list()
    for (i_idx in seq_len(nrow(inversions))) {
      inv_row <- inversions[i_idx]
      bs <- inv_row$start_w; be <- inv_row$end_w
      l3_scan <- function(W_sizes) scan_region(W_sizes, bs, be, 3L)
      sub_dt <- scan_collapse_nms_trim(pass3_W, l3_scan,
                                       paste0("L3 inv#", i_idx),
                                       gate_composition = TRUE)
      if (nrow(sub_dt) > 0L) {
        sub_dt[, level := 3L]
        sub_dt[, recursion_depth := 2L]
        sub_dt[, parent_envelope_idx := inv_row$parent_envelope_idx]
        l3_features[[length(l3_features) + 1L]] <- sub_dt
      }
    }
    if (length(l3_features) > 0L) {
      l3_dt <- rbindlist(l3_features, use.names = TRUE, fill = TRUE)
      cat("[D17mp] L3: ", nrow(l3_dt), " sub-features\n", sep = "")
      # We DO want Ward to run on these too if they are >= cluster_min_W,
      # so we add them to the inversions table.
      inversions <- rbindlist(list(inversions, l3_dt),
                              use.names = TRUE, fill = TRUE)
    } else {
      cat("[D17mp] L3: 0 sub-features\n")
    }
  }

  list(envelopes = envelopes, inversions = inversions)
}

# ---- Run all three passes ---------------------------------------------------

# ---- Run pipeline -----------------------------------------------------------
# Two architectures supported:
#   "horizontal": legacy 3 mutually-exclusive masking passes; Ward.D2 runs on
#                 the union of Pass 1 and Pass 2 candidates. Pass 1 outputs
#                 are themselves treated as candidates.
#   "nested":     Level 1 envelopes -> Level 2 inversions inside each envelope
#                 -> Ward.D2 inside each inversion. Levels are spatially
#                 nested. Envelopes and inversions are kept as separate
#                 catalogue rows; Ward only runs on the inversions.

envelopes_dt <- data.table()  # populated only in nested mode
refiner_audit_global <- data.table()  # populated by L2 refiner if used

if (architecture == "horizontal") {
  masked_runs <- list()
  all_pass_survivors <- list()

  p1 <- run_pass(pass1_W, masked_runs, 1L, "large")
  all_pass_survivors[[1]] <- p1$survivors
  masked_runs <- p1$masked_runs

  p2 <- run_pass(pass2_W, masked_runs, 2L, "medium")
  all_pass_survivors[[2]] <- p2$survivors
  masked_runs <- p2$masked_runs

  p3 <- run_pass(pass3_W, masked_runs, 3L, "small")
  all_pass_survivors[[3]] <- p3$survivors
  masked_runs <- p3$masked_runs

  discovered <- rbindlist(all_pass_survivors, use.names = TRUE, fill = TRUE)
} else {
  np <- run_nested_pipeline()
  envelopes_dt <- np$envelopes
  discovered   <- np$inversions
}

if (nrow(discovered) == 0L) {
  if (architecture == "nested" && nrow(envelopes_dt) > 0L) {
    cat("[D17mp] no Level 2 inversions discovered, but ",
        nrow(envelopes_dt), " Level 1 envelopes recorded\n", sep = "")
    # Fall through — we still want to write envelopes to the catalogue
    discovered <- data.table()  # empty, Ward loop will skip
  } else {
    cat("[D17mp] no candidates discovered\n")
    empty_dt <- data.table(
      chr = character(0), candidate_id = character(0),
      start_w = integer(0), end_w = integer(0),
      start_bp = integer(0), end_bp = integer(0),
      n_windows = integer(0), discovery_pass = integer(0),
      composition_consistency = numeric(0), status = character(0)
    )
    out_pqc <- file.path(outdir, paste0(chr_label, "_phase_qc_candidates.tsv"))
    fwrite(empty_dt, out_pqc, sep = "\t")
    out_main <- file.path(outdir, paste0(chr_label, "_d17_catalogue.tsv"))
    fwrite(empty_dt, out_main, sep = "\t")
    quit(save = "no", status = 0)
  }
}

if (nrow(discovered) > 0L) {
  discovered <- discovered[order(start_w)]
  discovered[, candidate_id := sprintf("%s_d17_%04d", chr_label, seq_len(.N))]
  discovered[, parent_candidate_id := NA_character_]
  # recursion_depth: keep 1 in nested mode (set by run_nested_pipeline);
  # set to 0 in horizontal mode (legacy behaviour).
  if (architecture == "horizontal") {
    discovered[, recursion_depth := 0L]
  }
  discovered[, status := "PASS"]
}

cat("\n[D17mp] === DISCOVERY SUMMARY ===\n")
cat("[D17mp] Total inversions/candidates: ", nrow(discovered), "\n", sep = "")
if (architecture == "nested") {
  cat("[D17mp] Level 1 envelopes:           ", nrow(envelopes_dt), "\n",
      sep = "")
}
if (nrow(discovered) > 0L) {
  cat("[D17mp] By pass:\n"); print(discovered[, .N, by = discovery_pass])
}

# ---- Ward.D2 internal clustering on pass 1 + 2 candidates ------------------

cluster_subblock <- function(start_w, end_w) {
  W <- end_w - start_w + 1L
  if (W < cluster_min_W) {
    return(list(
      K_best = 1L, cluster_ids = rep(1L, W),
      silhouette = NA_real_, dom_frac = 1, max_run_w = W,
      n_runs = 1L, interleave_ratio = 0,
      pattern_class = "TOO_SMALL_TO_CLUSTER"
    ))
  }
  M_sub <- sim_mat[start_w:end_w, start_w:end_w]
  cor_mat <- suppressWarnings(cor(t(M_sub), use = "pairwise.complete.obs"))
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1
  D <- as.dist(1 - cor_mat)

  H <- tryCatch(hclust(D, method = "ward.D2"), error = function(e) NULL)
  if (is.null(H)) {
    return(list(
      K_best = 1L, cluster_ids = rep(1L, W),
      silhouette = NA_real_, dom_frac = 1, max_run_w = W,
      n_runs = 1L, interleave_ratio = 0,
      pattern_class = "CLEAN"
    ))
  }

  best_K <- 1L
  best_sil <- -Inf
  best_clusters <- rep(1L, W)
  for (K in k_min:min(k_max, W - 1L)) {
    cl <- cutree(H, k = K)
    if (length(unique(cl)) < 2L) next
    s <- tryCatch(
      mean(silhouette(cl, D)[, "sil_width"]),
      error = function(e) NA_real_
    )
    if (!is.na(s) && s > best_sil) {
      best_sil <- s
      best_K <- K
      best_clusters <- as.integer(cl)
    }
  }
  if (is.infinite(best_sil) || best_sil < silhouette_min) {
    best_K <- 1L
    best_clusters <- rep(1L, W)
    best_sil <- NA_real_
  }

  cluster_sizes <- tabulate(best_clusters, nbins = max(best_clusters))
  dom_frac <- max(cluster_sizes) / W
  rl <- rle(best_clusters)
  max_run_w <- max(rl$lengths)
  n_runs <- length(rl$lengths)
  interleave_ratio <- max(0, (n_runs - best_K) / W)

  pattern_class <- if (best_K == 1L || dom_frac >= clean_dom_frac) {
    "CLEAN"
  } else if (n_runs / W > crossover_runs) {
    "DOUBLE_CROSSOVER"
  } else if (n_runs == best_K) {
    is_nested_inside <- FALSE
    cluster_ranges <- lapply(seq_len(best_K), function(k) {
      idx <- which(best_clusters == k); c(min(idx), max(idx))
    })
    for (i in seq_len(best_K)) {
      for (j in seq_len(best_K)) {
        if (i == j) next
        ri <- cluster_ranges[[i]]; rj <- cluster_ranges[[j]]
        if (ri[1] > rj[1] && ri[2] < rj[2]) {
          is_nested_inside <- TRUE
          break
        }
      }
      if (is_nested_inside) break
    }
    if (is_nested_inside) "NESTED_INSIDE" else "NESTED_ADJ"
  } else if (n_runs <= 2 * best_K) {
    "NESTED_ADJ"
  } else {
    "FRAGMENTED"
  }

  list(
    K_best = best_K, cluster_ids = best_clusters,
    silhouette = best_sil, dom_frac = dom_frac,
    max_run_w = max_run_w, n_runs = n_runs,
    interleave_ratio = interleave_ratio,
    pattern_class = pattern_class
  )
}

# ---- Emit children based on cluster result ---------------------------------

emit_children <- function(parent_row, clust_res) {
  children <- list()
  pc_class <- clust_res$pattern_class

  if (pc_class == "NESTED_ADJ") {
    rl <- rle(clust_res$cluster_ids)
    pos <- 1L
    for (i in seq_along(rl$lengths)) {
      run_len <- rl$lengths[i]
      cluster_id <- rl$values[i]
      if (run_len >= min_n_windows) {
        child_start <- parent_row$start_w + pos - 1L
        child_end   <- parent_row$start_w + pos + run_len - 2L
        block <- sim_mat[child_start:child_end, child_start:child_end]
        comp_res <- composition_consistency(child_start, child_end)
        child_row <- copy(parent_row)
        child_row[, candidate_id := paste0(parent_row$candidate_id, "_c", i)]
        child_row[, parent_candidate_id := parent_row$candidate_id]
        child_row[, recursion_depth := (parent_row$recursion_depth %||% 0L) + 1L]
        child_row[, start_w := child_start]
        child_row[, end_w := child_end]
        child_row[, n_windows := child_end - child_start + 1L]
        child_row[, mean_sim := mean(block, na.rm = TRUE)]
        child_row[, density_p70 := mean(block >= 0.70, na.rm = TRUE)]
        child_row[, density_adaptive := mean(block >= 0.55, na.rm = TRUE)]
        child_row[, W := child_end - child_start + 1L]
        child_row[, composition_consistency := comp_res$consistency]
        child_row[, n_inv_carriers := comp_res$n_inv]
        child_row[, source_cluster_id := cluster_id]
        child_row[, status := if (!is.na(comp_res$consistency) &&
                                  comp_res$consistency >= comp_min)
                              "PASS_CHILD" else "FAIL"]
        children[[length(children) + 1L]] <- child_row
      }
      pos <- pos + run_len
    }

  } else if (pc_class == "NESTED_INSIDE") {
    cluster_ranges <- lapply(seq_len(clust_res$K_best), function(k) {
      idx <- which(clust_res$cluster_ids == k)
      c(min(idx), max(idx), length(idx))
    })
    for (i in seq_len(clust_res$K_best)) {
      ri <- cluster_ranges[[i]]
      is_inside <- FALSE
      for (j in seq_len(clust_res$K_best)) {
        if (i == j) next
        rj <- cluster_ranges[[j]]
        if (ri[1] > rj[1] && ri[2] < rj[2]) { is_inside <- TRUE; break }
      }
      if (is_inside && ri[3] >= min_n_windows) {
        child_start <- parent_row$start_w + ri[1] - 1L
        child_end   <- parent_row$start_w + ri[2] - 1L
        block <- sim_mat[child_start:child_end, child_start:child_end]
        comp_res <- composition_consistency(child_start, child_end)
        child_row <- copy(parent_row)
        child_row[, candidate_id := paste0(parent_row$candidate_id, "_n", i)]
        child_row[, parent_candidate_id := parent_row$candidate_id]
        child_row[, recursion_depth := (parent_row$recursion_depth %||% 0L) + 1L]
        child_row[, start_w := child_start]
        child_row[, end_w := child_end]
        child_row[, n_windows := child_end - child_start + 1L]
        child_row[, mean_sim := mean(block, na.rm = TRUE)]
        child_row[, density_p70 := mean(block >= 0.70, na.rm = TRUE)]
        child_row[, density_adaptive := mean(block >= 0.55, na.rm = TRUE)]
        child_row[, W := child_end - child_start + 1L]
        child_row[, composition_consistency := comp_res$consistency]
        child_row[, n_inv_carriers := comp_res$n_inv]
        child_row[, source_cluster_id := i]
        child_row[, status := if (!is.na(comp_res$consistency) &&
                                  comp_res$consistency >= comp_min)
                              "PASS_CHILD" else "FAIL"]
        children[[length(children) + 1L]] <- child_row
      }
    }

  } else if (pc_class == "FRAGMENTED" && fragmented_emit) {
    # In a FRAGMENTED parent, Ward found multiple cluster IDs but they
    # interleave too much for clean nested-adjacent splitting. We can still
    # salvage the OBVIOUS contiguous runs that are large and well-separated.
    # We use a stricter min run length (fragmented_min_run_w, default 2x the
    # normal min_n_windows) so we only emit children with enough mass to be
    # plausibly real inversions hiding inside the envelope.
    rl <- rle(clust_res$cluster_ids)
    pos <- 1L
    for (i in seq_along(rl$lengths)) {
      run_len <- rl$lengths[i]
      cluster_id <- rl$values[i]
      if (run_len >= fragmented_min_run_w) {
        child_start <- parent_row$start_w + pos - 1L
        child_end   <- parent_row$start_w + pos + run_len - 2L
        block <- sim_mat[child_start:child_end, child_start:child_end]
        comp_res <- composition_consistency(child_start, child_end)
        child_row <- copy(parent_row)
        child_row[, candidate_id := paste0(parent_row$candidate_id, "_f", i)]
        child_row[, parent_candidate_id := parent_row$candidate_id]
        child_row[, recursion_depth := (parent_row$recursion_depth %||% 0L) + 1L]
        child_row[, start_w := child_start]
        child_row[, end_w := child_end]
        child_row[, n_windows := child_end - child_start + 1L]
        child_row[, mean_sim := mean(block, na.rm = TRUE)]
        child_row[, density_p70 := mean(block >= 0.70, na.rm = TRUE)]
        child_row[, density_adaptive := mean(block >= 0.55, na.rm = TRUE)]
        child_row[, W := child_end - child_start + 1L]
        child_row[, composition_consistency := comp_res$consistency]
        child_row[, n_inv_carriers := comp_res$n_inv]
        child_row[, source_cluster_id := cluster_id]
        child_row[, status := if (!is.na(comp_res$consistency) &&
                                  comp_res$consistency >= comp_min)
                              "PASS_CHILD_FRAG" else "FAIL"]
        children[[length(children) + 1L]] <- child_row
      }
      pos <- pos + run_len
    }
  }
  children
}

# Initialize cluster-result columns on the discovered table
extra_cols <- c("pattern_K", "pattern_class", "dominant_frac",
                "max_run_w", "n_runs", "interleave_ratio", "silhouette",
                "source_cluster_id")
if (nrow(discovered) > 0L) {
  discovered[, pattern_K := NA_integer_]
  discovered[, pattern_class := NA_character_]
  discovered[, dominant_frac := NA_real_]
  discovered[, max_run_w := NA_integer_]
  discovered[, n_runs := NA_integer_]
  discovered[, interleave_ratio := NA_real_]
  discovered[, silhouette := NA_real_]
  discovered[, source_cluster_id := NA_integer_]
}

cat("\n[D17mp] === INTERNAL CLUSTERING (Ward.D2) ===\n")

window_cluster_records <- list()
subblock_records <- list()
all_children_rows <- list()
z_audit_records <- list()

for (i in seq_len(nrow(discovered))) {
  row_i <- discovered[i]
  W <- row_i$end_w - row_i$start_w + 1L
  do_cluster <- (row_i$discovery_pass <= 2L) && (W >= cluster_min_W)

  if (!do_cluster) {
    # Atomic candidate (pass 3, or too small to cluster)
    discovered[i, pattern_K := 1L]
    discovered[i, pattern_class := if (W < cluster_min_W)
      "TOO_SMALL_TO_CLUSTER" else "NOT_CLUSTERED"]
    discovered[i, dominant_frac := 1]
    discovered[i, max_run_w := W]
    discovered[i, n_runs := 1L]
    discovered[i, interleave_ratio := 0]
    next
  }

  cat(sprintf("[D17mp]   %s  pass=%d  W=%d  comp=%.3f  ",
              row_i$candidate_id, row_i$discovery_pass, W,
              row_i$composition_consistency))
  clust_res <- cluster_subblock(row_i$start_w, row_i$end_w)
  cat(sprintf("class=%s  K=%d  dom=%.2f  runs=%d\n",
              clust_res$pattern_class, clust_res$K_best,
              clust_res$dom_frac, clust_res$n_runs))

  # Update parent row with cluster info
  discovered[i, pattern_K := clust_res$K_best]
  discovered[i, pattern_class := clust_res$pattern_class]
  discovered[i, dominant_frac := clust_res$dom_frac]
  discovered[i, max_run_w := clust_res$max_run_w]
  discovered[i, n_runs := clust_res$n_runs]
  discovered[i, interleave_ratio := clust_res$interleave_ratio]
  discovered[i, silhouette := clust_res$silhouette]

  # Per-window cluster records
  for (j in seq_along(clust_res$cluster_ids)) {
    window_idx <- row_i$start_w + j - 1L
    window_cluster_records[[length(window_cluster_records) + 1L]] <- data.table(
      chr = row_i$chr %||% chr_label,
      parent_candidate_id = row_i$candidate_id,
      window_idx = window_idx,
      start_bp = window_start_bp[window_idx],
      end_bp   = window_end_bp[window_idx],
      pattern_cluster_id = clust_res$cluster_ids[j]
    )
  }

  # Per-cluster summary
  for (k in unique(clust_res$cluster_ids)) {
    idx <- which(clust_res$cluster_ids == k)
    cluster_w_idx <- row_i$start_w + idx - 1L
    cluster_block <- sim_mat[cluster_w_idx, cluster_w_idx]
    subblock_records[[length(subblock_records) + 1L]] <- data.table(
      chr = row_i$chr %||% chr_label,
      parent_candidate_id = row_i$candidate_id,
      subblock_id = paste0(row_i$candidate_id, "_sb", k),
      pattern_cluster_id = k,
      start_w = row_i$start_w + min(idx) - 1L,
      end_w   = row_i$start_w + max(idx) - 1L,
      start_bp = window_start_bp[row_i$start_w + min(idx) - 1L],
      end_bp   = window_end_bp[row_i$start_w + max(idx) - 1L],
      n_windows_in_cluster = length(idx),
      mean_intensity_in_cluster = mean(cluster_block, na.rm = TRUE),
      contiguous = length(idx) == (max(idx) - min(idx) + 1L)
    )
  }

  # Emit children for NESTED_ADJ / NESTED_INSIDE / FRAGMENTED / DOUBLE_CROSSOVER
  if (clust_res$pattern_class %in% c("NESTED_ADJ", "NESTED_INSIDE")) {
    children <- emit_children(row_i, clust_res)
    if (length(children) > 0L) {
      ch_dt <- rbindlist(children, use.names = TRUE, fill = TRUE)

      # ---- Z-aware split validation (NESTED_ADJ only) ----------------------
      # For each adjacent pair of children (sorted by start_w), compute the
      # mean Z over their cross-block. Flag pairs as REAL / OK / OVER_SPLIT
      # and, if z_validate is enabled, merge OVER_SPLIT pairs back into a
      # single child. Original child rows are dumped to z_audit regardless.
      if (clust_res$pattern_class == "NESTED_ADJ" && nrow(ch_dt) >= 2L) {
        ch_dt <- ch_dt[order(start_w)]
        # Initialise per-child fields; will be filled below.
        ch_dt[, cross_z_left  := NA_real_]
        ch_dt[, cross_z_right := NA_real_]
        ch_dt[, z_flag        := "OK"]

        for (k in seq_len(nrow(ch_dt) - 1L)) {
          z <- crossblock_z(ch_dt$start_w[k], ch_dt$end_w[k],
                            ch_dt$start_w[k + 1L], ch_dt$end_w[k + 1L])
          ch_dt[k,     cross_z_right := z]
          ch_dt[k + 1L, cross_z_left  := z]
          if (!is.na(z)) {
            if (z >= z_merge_thresh)        flag_pair <- "OVER_SPLIT"
            else if (z <= z_split_thresh)   flag_pair <- "REAL"
            else                            flag_pair <- "OK"
            # Promote stronger flag (OVER_SPLIT trumps OK trumps REAL only
            # for OK; REAL is informative). Keep the per-side direction.
            if (flag_pair == "OVER_SPLIT") {
              ch_dt[k,     z_flag := "OVER_SPLIT"]
              ch_dt[k + 1L, z_flag := "OVER_SPLIT"]
            } else if (flag_pair == "REAL") {
              if (ch_dt$z_flag[k]     == "OK") ch_dt[k,     z_flag := "REAL"]
              if (ch_dt$z_flag[k + 1L] == "OK") ch_dt[k + 1L, z_flag := "REAL"]
            }
          }
        }

        all_z <- c(ch_dt$cross_z_left, ch_dt$cross_z_right)
        all_z <- all_z[is.finite(all_z)]
        z_range_str <- if (length(all_z) == 0L) "NA, NA" else
          paste(sprintf("%.2f", range(all_z)), collapse = ", ")
        cat(sprintf(
          "[D17mp]   z_validate %s: cross_z range = [%s], over_split=%d  real=%d  ok=%d\n",
          row_i$candidate_id, z_range_str,
          sum(ch_dt$z_flag == "OVER_SPLIT"),
          sum(ch_dt$z_flag == "REAL"),
          sum(ch_dt$z_flag == "OK")))

        # Always record the pre-merge child rows in audit
        z_audit_records[[length(z_audit_records) + 1L]] <-
          copy(ch_dt)[, audit_kind := "pre_merge"]

        # If z_validate is enabled and any adjacent pair fails, do greedy
        # left-to-right merging on consecutive OVER_SPLIT pairs.
        if (z_validate && any(!is.na(ch_dt$cross_z_right) &
                              ch_dt$cross_z_right >= z_merge_thresh)) {
          merged_rows <- list()
          k <- 1L
          while (k <= nrow(ch_dt)) {
            cur <- copy(ch_dt[k])
            # Walk forward absorbing any OVER_SPLIT-linked neighbour
            while (k < nrow(ch_dt) &&
                   !is.na(cur$cross_z_right) &&
                   cur$cross_z_right >= z_merge_thresh) {
              nxt <- ch_dt[k + 1L]
              cur[, candidate_id := paste0(cur$candidate_id, "+",
                                           sub(".*_(c\\d+)$", "\\1",
                                               nxt$candidate_id))]
              cur[, end_w := nxt$end_w]
              cur[, n_windows := cur$end_w - cur$start_w + 1L]
              # Recompute simple stats over the merged span
              mb <- sim_mat[cur$start_w:cur$end_w, cur$start_w:cur$end_w]
              cur[, mean_sim := mean(mb, na.rm = TRUE)]
              cur[, density_p70 := mean(mb >= 0.70, na.rm = TRUE)]
              cur[, density_adaptive := mean(mb >= 0.55, na.rm = TRUE)]
              cmp <- composition_consistency(cur$start_w, cur$end_w)
              cur[, composition_consistency := cmp$consistency]
              cur[, n_inv_carriers := cmp$n_inv]
              # Inherit cross_z_right from the absorbed neighbour for the
              # next loop iteration's exit check
              cur[, cross_z_right := nxt$cross_z_right]
              cur[, z_flag := "MERGED"]
              k <- k + 1L
            }
            merged_rows[[length(merged_rows) + 1L]] <- cur
            k <- k + 1L
          }
          new_ch <- rbindlist(merged_rows, use.names = TRUE, fill = TRUE)
          n_pre  <- nrow(ch_dt); n_post <- nrow(new_ch)
          if (n_post < n_pre) {
            cat(sprintf("[D17mp]     -> merged %d -> %d children\n",
                        n_pre, n_post))
            ch_dt <- new_ch
            # Record the post-merge rows in audit too
            z_audit_records[[length(z_audit_records) + 1L]] <-
              copy(ch_dt)[, audit_kind := "post_merge"]
          }
        }
      }
      # ---- end of Z validation ---------------------------------------------

      all_children_rows[[length(all_children_rows) + 1L]] <- ch_dt
      # If NESTED_ADJ produced children, parent becomes PASS_PARENT
      if (clust_res$pattern_class == "NESTED_ADJ") {
        discovered[i, status := "PASS_PARENT"]
      }
      # NESTED_INSIDE keeps parent as PASS (Phase QC genotypes both)
    }
  } else if (clust_res$pattern_class == "DOUBLE_CROSSOVER") {
    discovered[i, status := "PASS_COMPLEX"]
  } else if (clust_res$pattern_class == "FRAGMENTED") {
    if (fragmented_emit) {
      children <- emit_children(row_i, clust_res)
      if (length(children) > 0L) {
        ch_dt <- rbindlist(children, use.names = TRUE, fill = TRUE)
        all_children_rows[[length(all_children_rows) + 1L]] <- ch_dt
        cat(sprintf(
          "[D17mp]   FRAGMENTED %s emitted %d salvage children (>= %d windows)\n",
          row_i$candidate_id, nrow(ch_dt), fragmented_min_run_w))
        discovered[i, status := "FRAGMENTED_PARENT"]
      } else {
        discovered[i, status := "FRAGMENTED"]
      }
    } else {
      discovered[i, status := "FRAGMENTED"]
    }
  }
  # CLEAN: parent stays PASS, no children
}

# Combine parents and children
catalogue <- discovered
if (length(all_children_rows) > 0L) {
  ch_combined <- rbindlist(all_children_rows, use.names = TRUE, fill = TRUE)
  # Children inherit pass index from parent for visualization purposes
  catalogue <- rbindlist(list(discovered, ch_combined), use.names = TRUE, fill = TRUE)
}

# In nested mode, prepend Level 1 envelopes as context rows.
if (architecture == "nested" && nrow(envelopes_dt) > 0L) {
  env_rows <- copy(envelopes_dt)
  env_rows[, candidate_id := sprintf("%s_d17_env%04d", chr_label,
                                     seq_len(.N))]
  env_rows[, parent_candidate_id := NA_character_]
  env_rows[, status := "ENVELOPE"]
  env_rows[, pattern_class := "ENVELOPE"]
  catalogue <- rbindlist(list(env_rows, catalogue),
                         use.names = TRUE, fill = TRUE)
}

catalogue <- catalogue[order(start_w)]

# Add chr + bp columns
catalogue[, chr := chr_label]
catalogue[, start_bp := window_start_bp[start_w]]
catalogue[, end_bp   := window_end_bp[end_w]]
if ("W" %in% names(catalogue)) setnames(catalogue, "W", "scale_W")

# Drop tiny candidates from the final catalogue
catalogue <- catalogue[n_windows >= min_n_windows]

# ---- Build final output -----------------------------------------------------

main_cols <- c(
  "chr", "candidate_id", "parent_candidate_id", "recursion_depth",
  "discovery_pass",
  "start_w", "end_w", "start_bp", "end_bp", "n_windows", "scale_W",
  "mean_sim", "density_p70", "density_adaptive",
  "composition_consistency", "n_inv_carriers",
  "pattern_K", "pattern_class", "dominant_frac",
  "max_run_w", "n_runs", "interleave_ratio", "silhouette",
  "source_cluster_id",
  "trim_start_w_original", "trim_end_w_original",
  "trim_left_dropped", "trim_right_dropped", "trim_total_dropped",
  "cross_z_left", "cross_z_right", "z_flag",
  "status"
)
for (cc in main_cols) {
  if (!cc %in% names(catalogue)) catalogue[, (cc) := NA]
}
catalogue_out <- catalogue[, ..main_cols]

# ---- Console summary --------------------------------------------------------

cat("\n[D17mp] === FINAL CATALOGUE SUMMARY for ", chr_label, " ===\n", sep = "")
cat("[D17mp] total rows: ", nrow(catalogue_out), "\n", sep = "")
cat("[D17mp] by status:\n"); print(catalogue_out[, .N, by = status])
cat("[D17mp] by discovery_pass:\n")
print(catalogue_out[, .N, by = discovery_pass])
cat("[D17mp] by pattern_class (PASS-like only):\n")
print(catalogue_out[status %in% c("PASS","PASS_CHILD","PASS_CHILD_FRAG",
                                    "PASS_PARENT","PASS_COMPLEX",
                                    "FRAGMENTED_PARENT"),
                    .N, by = pattern_class])

phase_qc_cands <- catalogue_out[status %in% c("PASS", "PASS_CHILD",
                                              "PASS_CHILD_FRAG",
                                              "PASS_COMPLEX")]
cat("[D17mp] => Phase QC candidates: ",
    nrow(phase_qc_cands), "\n", sep = "")

if (nrow(phase_qc_cands) > 0L) {
  cat("[D17mp] biggest Phase QC candidates:\n")
  big <- phase_qc_cands[order(-n_windows)][1:min(10L, .N)]
  for (i in seq_len(nrow(big))) {
    cat(sprintf(
      "[D17mp]   %s  pass=%d  d=%d  %.2f-%.2f Mb  nW=%d  comp=%.3f  class=%s\n",
      big$candidate_id[i],
      big$discovery_pass[i] %||% NA_integer_,
      big$recursion_depth[i],
      big$start_bp[i] / 1e6, big$end_bp[i] / 1e6,
      big$n_windows[i], big$composition_consistency[i],
      big$pattern_class[i] %||% "NA"
    ))
  }
}

# ---- Write -----------------------------------------------------------------

if (!dry_run) {
  out_main <- file.path(outdir, paste0(chr_label, "_d17_catalogue.tsv"))
  fwrite(catalogue_out, out_main, sep = "\t")
  cat("\n[D17mp] catalogue written: ", out_main, "\n", sep = "")

  if (length(subblock_records) > 0L) {
    sb_dt <- rbindlist(subblock_records, use.names = TRUE, fill = TRUE)
    out_sb <- file.path(outdir, paste0(chr_label, "_d17_subblocks.tsv"))
    fwrite(sb_dt, out_sb, sep = "\t")
    cat("[D17mp] subblocks written: ", out_sb, "  (", nrow(sb_dt), " rows)\n", sep = "")
  }

  if (length(window_cluster_records) > 0L) {
    wc_dt <- rbindlist(window_cluster_records, use.names = TRUE, fill = TRUE)
    out_wc <- file.path(outdir, paste0(chr_label, "_d17_window_clusters.tsv"))
    fwrite(wc_dt, out_wc, sep = "\t")
    cat("[D17mp] window clusters written: ", out_wc, "  (",
        nrow(wc_dt), " rows)\n", sep = "")
  }

  if (length(z_audit_records) > 0L) {
    za_dt <- rbindlist(z_audit_records, use.names = TRUE, fill = TRUE)
    out_za <- file.path(outdir, paste0(chr_label, "_d17_z_audit.tsv"))
    fwrite(za_dt, out_za, sep = "\t")
    cat("[D17mp] z audit written: ", out_za, "  (",
        nrow(za_dt), " rows)\n", sep = "")
  }

  if (nrow(refiner_audit_global) > 0L) {
    out_ra <- file.path(outdir,
                        paste0(chr_label, "_d17_refiner_audit.tsv"))
    fwrite(refiner_audit_global, out_ra, sep = "\t")
    cat("[D17mp] refiner audit written: ", out_ra, "  (",
        nrow(refiner_audit_global), " rows)\n", sep = "")
  }

  pqc_cols <- c("chr", "candidate_id", "start_bp", "end_bp", "n_windows",
                "discovery_pass", "composition_consistency", "n_inv_carriers",
                "pattern_class", "recursion_depth", "status")
  pqc_cols <- intersect(pqc_cols, names(phase_qc_cands))
  out_pqc <- file.path(outdir, paste0(chr_label, "_phase_qc_candidates.tsv"))
  fwrite(phase_qc_cands[, ..pqc_cols], out_pqc, sep = "\t")
  cat("[D17mp] Phase QC input written: ", out_pqc, "  (",
      nrow(phase_qc_cands), " rows)\n", sep = "")
} else {
  cat("\n[D17mp] DRY RUN — no files written\n")
}

cat("[D17mp] done\n")
