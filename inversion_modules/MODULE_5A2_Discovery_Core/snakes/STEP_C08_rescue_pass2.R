#!/usr/bin/env Rscript

# =============================================================================
# STEP10k_rescue_pass2.R  (v8.0)
#
# SECOND-PASS RESCUE for ambiguous/shoulder windows.
#
# Motivation: Pass 1 (STEP10e) uses conservative thresholds that miss:
#   - weaker shoulders near real candidates
#   - bridge-heavy windows
#   - partial regimes with 1.5 ≤ |z| < 3.0
#
# This module identifies rescue candidates and runs a lighter collection pass.
# ALL rescued windows are labeled distinctly — NEVER silently merged into pass 1.
#
# Labels:
#   CORE_PASS1       — from pass 1 (not touched here)
#   RESCUED_PASS2    — rescued in this pass with adequate continuity
#   SHOULDER_RESCUE  — weaker rescue adjacent to pass-1 regions
#   WEAK_EXTENSION   — minimal-confidence extension at region edges
#
# INPUTS:
#   <step10_outprefix>.mds.rds        — per_chr data from STEP10(v2)
#   <snake1_dir>/snake_windows.tsv.gz — pass-1 snake results
#   <snake1_dir>/snake_window_states.tsv.gz — pass-1 window states
#   [optional] snake2_track, snake3_track for support gating
#
# OUTPUTS:
#   rescue_windows.tsv.gz           — rescued windows with provenance labels
#   rescue_regions.tsv.gz           — rescued region summaries
#   rescue_decision_log.tsv.gz      — rescue extension decisions
#   rescue_summary.tsv              — per-chromosome counts
#
# Usage:
#   Rscript STEP10k_rescue_pass2.R <step10_outprefix> <snake1_dir> <outdir> \
#     [--rescue_z_min 1.5] [--rescue_z_max 3.0] \
#     [--rescue_proximity_bp 200000] [--rescue_min_windows 2] \
#     [--snake2_dir <path>] [--snake3_dir <path>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript STEP10k_rescue_pass2.R <step10_outprefix> <snake1_dir> <outdir> ...")
}

step10_prefix <- args[1]
snake1_dir    <- args[2]
outdir        <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Defaults
RESCUE_Z_MIN       <- 1.5
RESCUE_Z_MAX       <- 3.0
RESCUE_PROXIMITY_BP <- 200000L
RESCUE_MIN_WINDOWS <- 2L
RESCUE_ACCEPT      <- 0.40
RESCUE_TOLER       <- 0.20
RESCUE_MAX_GAP     <- 2L
RESCUE_MAX_BAD     <- 2L

snake2_dir <- NULL
snake3_dir <- NULL

i <- 4L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--rescue_z_min" && i < length(args)) {
    RESCUE_Z_MIN <- as.numeric(args[i + 1]); i <- i + 2L
  } else if (a == "--rescue_z_max" && i < length(args)) {
    RESCUE_Z_MAX <- as.numeric(args[i + 1]); i <- i + 2L
  } else if (a == "--rescue_proximity_bp" && i < length(args)) {
    RESCUE_PROXIMITY_BP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--rescue_min_windows" && i < length(args)) {
    RESCUE_MIN_WINDOWS <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--snake2_dir" && i < length(args)) {
    snake2_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--snake3_dir" && i < length(args)) {
    snake3_dir <- args[i + 1]; i <- i + 2L
  } else {
    i <- i + 1L
  }
}

message("[STEP10k] Rescue pass 2")
message("[STEP10k] z range: [", RESCUE_Z_MIN, ", ", RESCUE_Z_MAX, ")")
message("[STEP10k] Proximity: ", RESCUE_PROXIMITY_BP, " bp")

# =============================================================================
# LOAD PASS-1 RESULTS
# =============================================================================

mds_obj <- readRDS(paste0(step10_prefix, ".mds.rds"))
per_chr <- mds_obj$per_chr

# Load pass-1 window states
p1_states_file <- file.path(snake1_dir, "snake_window_states.tsv.gz")
p1_windows_file <- file.path(snake1_dir, "snake_windows.tsv.gz")

p1_states <- if (file.exists(p1_states_file)) fread(p1_states_file) else data.table()
p1_windows <- if (file.exists(p1_windows_file)) fread(p1_windows_file) else data.table()

# Load optional Snake 2/3 support
s2_track <- if (!is.null(snake2_dir) && file.exists(file.path(snake2_dir, "snake2_track.tsv.gz"))) {
  fread(file.path(snake2_dir, "snake2_track.tsv.gz"))
} else data.table()

s3_track <- if (!is.null(snake3_dir) && file.exists(file.path(snake3_dir, "snake3_track.tsv.gz"))) {
  fread(file.path(snake3_dir, "snake3_track.tsv.gz"))
} else data.table()

message("[STEP10k] Pass-1 window states: ", nrow(p1_states))
message("[STEP10k] Pass-1 collected windows: ", nrow(p1_windows))

# =============================================================================
# HELPERS (simplified from STEP10e)
# =============================================================================

make_sim_mat <- function(dmat) {
  dmax <- quantile(dmat[is.finite(dmat)], 0.95, na.rm = TRUE)
  if (!is.finite(dmax) || dmax == 0) dmax <- 1
  sim <- 1 - pmin(dmat / dmax, 1)
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  sim
}

# =============================================================================
# IDENTIFY RESCUE CANDIDATES
# =============================================================================

all_rescue_windows <- list()
all_rescue_regions <- list()
all_rescue_log     <- list()
all_rescue_summary <- list()
rescue_region_id   <- 0L

for (chr in names(per_chr)) {
  chr_obj <- per_chr[[chr]]
  if (is.null(chr_obj)) next

  dt <- as.data.table(chr_obj$out_dt)
  dt <- dt[order(start_bp)]
  n <- nrow(dt)
  if (n < 3) next

  # Get z-scores
  z_cols <- grep("^MDS[0-9]+_z$", names(dt), value = TRUE)
  if (length(z_cols) == 0) next

  dt[, max_abs_z := apply(.SD, 1, function(x) max(abs(x), na.rm = TRUE)), .SDcols = z_cols]

  # Pass-1 collected windows for this chromosome
  p1_chr <- p1_windows[chrom == chr]
  p1_wids <- if (nrow(p1_chr) > 0) unique(p1_chr$global_window_id) else integer(0)

  # Pass-1 region boundaries
  if (nrow(p1_chr) > 0) {
    p1_regions <- p1_chr[, .(region_start = min(start_bp), region_end = max(end_bp)),
                          by = snake_id]
  } else {
    p1_regions <- data.table(snake_id = integer(), region_start = numeric(), region_end = numeric())
  }

  # Identify rescue candidates:
  # 1. Not already collected by pass 1
  # 2. z in [RESCUE_Z_MIN, RESCUE_Z_MAX)
  # 3. Within RESCUE_PROXIMITY_BP of a pass-1 region
  dt[, is_p1_collected := global_window_id %in% p1_wids]
  dt[, z_in_range := max_abs_z >= RESCUE_Z_MIN & max_abs_z < RESCUE_Z_MAX]

  dt[, near_p1 := FALSE]
  for (ri in seq_len(nrow(p1_regions))) {
    rstart <- p1_regions$region_start[ri] - RESCUE_PROXIMITY_BP
    rend   <- p1_regions$region_end[ri] + RESCUE_PROXIMITY_BP
    dt[start_bp >= rstart & end_bp <= rend, near_p1 := TRUE]
  }

  # Optional: boost if Snake 2 or 3 shows WEAK+ support
  dt[, s2_support := FALSE]
  dt[, s3_support := FALSE]
  if (nrow(s2_track) > 0) {
    s2_chr <- s2_track[chrom == chr & snake2_status %in% c("PASS", "WEAK")]
    if (nrow(s2_chr) > 0) dt[global_window_id %in% s2_chr$global_window_id, s2_support := TRUE]
  }
  if (nrow(s3_track) > 0) {
    s3_chr <- s3_track[chrom == chr & snake3_status %in% c("PASS", "WEAK")]
    if (nrow(s3_chr) > 0) dt[global_window_id %in% s3_chr$global_window_id, s3_support := TRUE]
  }

  rescue_candidates <- dt[!is_p1_collected & z_in_range & (near_p1 | s2_support | s3_support)]

  if (nrow(rescue_candidates) < RESCUE_MIN_WINDOWS) {
    message("[STEP10k] ", chr, ": ", nrow(rescue_candidates), " candidates — too few, skip")
    all_rescue_summary[[chr]] <- data.table(
      chrom = chr, n_windows = n, n_rescue_candidates = nrow(rescue_candidates),
      n_rescued = 0L, n_regions = 0L
    )
    next
  }

  message("[STEP10k] ", chr, ": ", nrow(rescue_candidates), " rescue candidates")

  # ── Build sim_mat for rescue ────────────────────────────────────────
  dmat <- chr_obj$dmat
  if (nrow(dmat) != n) {
    nn <- min(nrow(dmat), n)
    dt <- dt[seq_len(nn)]
    dmat <- dmat[seq_len(nn), seq_len(nn), drop = FALSE]
  }
  sim_mat <- make_sim_mat(dmat)

  # ── Simple contiguous rescue collection ─────────────────────────────
  # Group rescue candidates into contiguous or near-contiguous runs
  rc_idx <- which(dt$global_window_id %in% rescue_candidates$global_window_id)
  if (length(rc_idx) < RESCUE_MIN_WINDOWS) next

  # Cluster by position
  runs <- list()
  cur_run <- rc_idx[1]
  for (ri in rc_idx[-1]) {
    if (ri - max(cur_run) <= RESCUE_MAX_GAP + 1L) {
      cur_run <- c(cur_run, ri)
    } else {
      if (length(cur_run) >= RESCUE_MIN_WINDOWS) runs[[length(runs) + 1]] <- cur_run
      cur_run <- ri
    }
  }
  if (length(cur_run) >= RESCUE_MIN_WINDOWS) runs[[length(runs) + 1]] <- cur_run

  # Evaluate each run
  for (run in runs) {
    # Check internal coherence
    if (length(run) >= 2) {
      coh_vals <- numeric(length(run) - 1)
      for (k in seq_along(coh_vals)) {
        s <- sim_mat[run[k], run[k + 1]]
        coh_vals[k] <- if (is.finite(s)) s else 0
      }
      mean_coh <- mean(coh_vals)
    } else {
      mean_coh <- 1.0
    }

    # Classify rescue label
    mean_z <- mean(dt$max_abs_z[run], na.rm = TRUE)
    any_s2 <- any(dt$s2_support[run])
    any_s3 <- any(dt$s3_support[run])

    if (mean_coh >= RESCUE_ACCEPT && (any_s2 || any_s3)) {
      label <- "RESCUED_PASS2"
    } else if (mean_coh >= RESCUE_TOLER && any(dt$near_p1[run])) {
      label <- "SHOULDER_RESCUE"
    } else if (mean_coh >= RESCUE_TOLER) {
      label <- "WEAK_EXTENSION"
    } else {
      next  # don't rescue
    }

    rescue_region_id <- rescue_region_id + 1L

    for (k in seq_along(run)) {
      wi <- run[k]
      all_rescue_windows[[length(all_rescue_windows) + 1]] <- data.table(
        chrom = chr,
        global_window_id = dt$global_window_id[wi],
        start_bp = dt$start_bp[wi],
        end_bp = dt$end_bp[wi],
        rescue_region_id = rescue_region_id,
        rescue_label = label,
        max_abs_z = round(dt$max_abs_z[wi], 4),
        near_p1 = dt$near_p1[wi],
        s2_support = dt$s2_support[wi],
        s3_support = dt$s3_support[wi],
        local_coherence = if (k > 1) round(sim_mat[run[k - 1], run[k]], 4) else NA_real_
      )
    }

    all_rescue_regions[[length(all_rescue_regions) + 1]] <- data.table(
      rescue_region_id = rescue_region_id,
      chrom = chr,
      start_bp = min(dt$start_bp[run]),
      end_bp = max(dt$end_bp[run]),
      n_windows = length(run),
      rescue_label = label,
      mean_coherence = round(mean_coh, 4),
      mean_z = round(mean_z, 4),
      has_s2_support = any_s2,
      has_s3_support = any_s3,
      near_p1_region = any(dt$near_p1[run])
    )
  }

  n_rescued <- sum(vapply(runs, function(r) length(r), integer(1)))
  all_rescue_summary[[chr]] <- data.table(
    chrom = chr, n_windows = n,
    n_rescue_candidates = nrow(rescue_candidates),
    n_rescued = n_rescued,
    n_regions = length(runs)
  )
  message("[STEP10k] ", chr, ": rescued ", n_rescued, " windows in ", length(runs), " regions")
}

# =============================================================================
# WRITE
# =============================================================================

rescue_windows <- if (length(all_rescue_windows) > 0) rbindlist(all_rescue_windows) else {
  data.table(chrom = character(), rescue_label = character())
}
rescue_regions <- if (length(all_rescue_regions) > 0) rbindlist(all_rescue_regions) else {
  data.table(rescue_region_id = integer())
}
rescue_summary <- if (length(all_rescue_summary) > 0) rbindlist(all_rescue_summary) else {
  data.table(chrom = character())
}

f1 <- file.path(outdir, "rescue_windows.tsv.gz")
f2 <- file.path(outdir, "rescue_regions.tsv.gz")
f3 <- file.path(outdir, "rescue_summary.tsv")

fwrite(rescue_windows, f1, sep = "\t")
fwrite(rescue_regions, f2, sep = "\t")
fwrite(rescue_summary, f3, sep = "\t")

message("\n[DONE] STEP10k rescue pass 2 complete")
message("  ", f1, " (", nrow(rescue_windows), " windows)")
message("  ", f2, " (", nrow(rescue_regions), " regions)")
message("  ", f3)

if (nrow(rescue_regions) > 0) {
  label_tab <- rescue_regions[, .N, by = rescue_label]
  for (li in seq_len(nrow(label_tab))) {
    message("  ", label_tab$rescue_label[li], ": ", label_tab$N[li], " regions")
  }
}
