#!/usr/bin/env Rscript

# =============================================================================
# STEP10l_snake4_threeband_support.R  (v2.0)
#
# SNAKE 4: Rolling 3-band support + middle-core regime identity track.
#
# READS FROM SNAKE 2 v2 (no redundant dosage recomputation):
#   Snake 2 v2 already computes per-window:
#     - 3-band k=3 assignment (middle_idx, tailA_idx, tailB_idx)
#     - fuzzy middle_weight (Gaussian soft membership)
#     - middle_band_preservation between adjacent windows
#     - combined_score, mean_middle_stability, etc.
#
#   Snake 4 reads these precomputed results and adds:
#     - Rolling block aggregation (BLOCK_SIZE windows)
#     - Middle-core regime clustering
#     - PARTIAL status for composite detection
#
# THREE CURVES PER WINDOW (evaluated on a rolling block):
#
#   Curve A — 3-band support score:
#     How many windows in the block have valid 3-band structure with
#     clean cluster separation? Uses Snake 2's combined_score + the
#     fraction of block windows where k=3 succeeded.
#
#   Curve B — Middle-stripe fraction:
#     What fraction of samples are consistently assigned to the middle
#     band across the block? Uses Snake 2's hard band assignments.
#     Safe name: "middle stripe" not "het" (may be mixed in composites).
#
#   Curve C — Middle-core membership continuity:
#     How stable is the middle set across the block? Uses both hard
#     Jaccard and fuzzy Jaccard (from Snake 2's soft middle_weight).
#     THE MOST DIAGNOSTIC curve — catches composites directly.
#
# REGIME CLUSTERING:
#   After per-window scoring, clusters all windows per chromosome by
#   their middle-membership sets into a few recurring regimes using
#   greedy single-linkage on Jaccard. Each regime = a set of windows
#   where the same samples dominate the middle stripe.
#   Output: middle_regime_id per window + regime summary table with
#   core member sample indices.
#
# COMBINED STATUS:
#   PASS    = clean 3-band + stable middle fraction + same middle core
#   PARTIAL = 3-band exists but middle core shifts (composite signal)
#   WEAK    = marginal pattern
#   FAIL    = no 3-band structure
#
# INPUTS:
#   <snake2_dir>/snake2_track.tsv.gz            — per-window scores
#   <snake2_dir>/snake2_band_assignments.tsv.gz — per-window band data
#   <step10_outprefix>.mds.rds                  — window coordinates (grid)
#
# OUTPUTS:
#   snake4_track.tsv.gz      — per-window: 3 curves + combined + status + regime
#   snake4_regimes.tsv.gz    — regime summary (which samples define each regime)
#   snake4_summary.tsv       — per-chromosome summary
#
# Usage:
#   Rscript STEP10l_snake4_threeband_support.R \
#     <step10_outprefix> <snake2_dir> <outdir> \
#     [--block_size 20] [--block_step 1]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript STEP10l_snake4_threeband_support.R ",
       "<step10_outprefix> <snake2_dir> <outdir> ",
       "[--block_size 20] [--block_step 1]")
}

step10_prefix <- args[1]
snake2_dir    <- args[2]
outdir        <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

BLOCK_SIZE            <- 20L   # windows per rolling block (20 x 100 SNPs ~ 2000 SNPs)
BLOCK_STEP            <- 1L    # step between blocks (1 = score every window)
MIN_VALID_IN_BLOCK    <- 10L   # min windows with band data in a block

# Curve A: 3-band support
THREEBAND_PASS        <- 0.35  # mean combined_score for block to count as 3-band
THREEBAND_WEAK        <- 0.20

# Curve B: middle-stripe fraction
MIDDLE_FRAC_PASS      <- 0.15  # >=15% of samples consistently in middle
MIDDLE_FRAC_WEAK      <- 0.08

# Curve C: middle membership continuity
MIDDLE_CONT_PASS      <- 0.50  # min Jaccard across block
MIDDLE_CONT_WEAK      <- 0.25

# Regime clustering
REGIME_JACCARD_THRESH <- 0.35  # windows with middle Jaccard >= this = same regime
MIN_REGIME_WINDOWS    <- 5L    # minimum windows to define a regime

# Weights for combined score
W_THREEBAND <- 0.30  # Curve A
W_MIDFRAC   <- 0.25  # Curve B
W_MIDCONT   <- 0.45  # Curve C (most diagnostic)

# Parse optional args
i <- 4L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--block_size" && i < length(args)) {
    BLOCK_SIZE <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--block_step" && i < length(args)) {
    BLOCK_STEP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--min_valid" && i < length(args)) {
    MIN_VALID_IN_BLOCK <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--regime_thresh" && i < length(args)) {
    REGIME_JACCARD_THRESH <- as.numeric(args[i + 1]); i <- i + 2L
  } else {
    i <- i + 1L
  }
}

message("[S4] Snake 4: Rolling 3-band support + regime identity")
message("[S4] Block size: ", BLOCK_SIZE, " windows, step: ", BLOCK_STEP)

# =============================================================================
# LOAD SNAKE 2 OUTPUTS
# =============================================================================

s2_track_file <- file.path(snake2_dir, "snake2_track.tsv.gz")
if (!file.exists(s2_track_file)) stop("Missing Snake 2 track: ", s2_track_file)
s2_track <- fread(s2_track_file)
message("[S4] Snake 2 track: ", nrow(s2_track), " windows")

s2_band_file <- file.path(snake2_dir, "snake2_band_assignments.tsv.gz")
if (!file.exists(s2_band_file)) stop("Missing Snake 2 band data: ", s2_band_file,
  "\n  -> Run Snake 2 v2 (STEP10g with band export) first")
s2_bands <- fread(s2_band_file)
message("[S4] Snake 2 band assignments: ", nrow(s2_bands), " windows")

mds_rds <- paste0(step10_prefix, ".mds.rds")
if (!file.exists(mds_rds)) stop("Missing: ", mds_rds)
mds_obj <- readRDS(mds_rds)

# =============================================================================
# HELPERS
# =============================================================================

parse_idx <- function(csv_str) {
  if (is.na(csv_str) || nchar(csv_str) == 0) return(integer(0))
  as.integer(strsplit(csv_str, ",")[[1]])
}

parse_wt <- function(csv_str) {
  if (is.na(csv_str) || nchar(csv_str) == 0) return(numeric(0))
  as.numeric(strsplit(csv_str, ",")[[1]])
}

jaccard <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

fuzzy_jaccard <- function(wa, wb) {
  if (length(wa) != length(wb) || length(wa) == 0) return(0)
  s_min <- sum(pmin(wa, wb))
  s_max <- sum(pmax(wa, wb))
  if (s_max == 0) return(0)
  s_min / s_max
}

# =============================================================================
# MAIN: PER-CHROMOSOME PROCESSING
# =============================================================================

chroms <- sort(unique(s2_track$chrom))
all_track   <- list()
all_summary <- list()
all_middle_sets <- list()

for (chr in chroms) {
  chr_obj <- mds_obj$per_chr[[chr]]
  if (is.null(chr_obj)) next
  grid <- as.data.table(chr_obj$out_dt)[order(start_bp)]
  n_win <- nrow(grid)

  message("\n[S4] ======= ", chr, " (", n_win, " windows) =======")
  if (n_win < BLOCK_SIZE) { message("[SKIP] fewer windows than block size"); next }

  chr_track <- s2_track[chrom == chr]
  chr_bands <- s2_bands[chrom == chr]

  # Index band data by global_window_id
  band_by_wid <- list()
  for (bi in seq_len(nrow(chr_bands))) {
    wid <- chr_bands$global_window_id[bi]
    band_by_wid[[as.character(wid)]] <- list(
      middle_idx = parse_idx(chr_bands$middle_idx_csv[bi]),
      tailA_idx  = parse_idx(chr_bands$tailA_idx_csv[bi]),
      tailB_idx  = parse_idx(chr_bands$tailB_idx_csv[bi]),
      middle_weight = parse_wt(chr_bands$middle_weight_csv[bi]),
      n_middle = chr_bands$n_middle[bi],
      middle_fraction = chr_bands$middle_fraction[bi]
    )
  }

  # Index track scores by window ID
  score_by_wid <- list()
  for (ti in seq_len(nrow(chr_track))) {
    wid <- chr_track$global_window_id[ti]
    score_by_wid[[as.character(wid)]] <- list(
      combined_score = chr_track$combined_score[ti],
      mean_middle_stability = chr_track$mean_middle_stability[ti],
      mean_hard_nn = chr_track$mean_hard_nn[ti],
      mean_fuzzy_nn = chr_track$mean_fuzzy_nn[ti],
      snake2_status = chr_track$snake2_status[ti]
    )
  }

  half_block <- BLOCK_SIZE %/% 2L

  # ── Rolling block scoring ──────────────────────────────────────────
  for (wi in seq(from = 1L, to = n_win, by = BLOCK_STEP)) {
    b_start <- max(1L, wi - half_block)
    b_end   <- min(n_win, wi + half_block - 1L)
    block_idx <- b_start:b_end
    block_wids <- grid$global_window_id[block_idx]

    has_band <- vapply(block_wids, function(w) {
      !is.null(band_by_wid[[as.character(w)]])
    }, logical(1))
    n_valid <- sum(has_band)
    wid_center <- grid$global_window_id[wi]

    if (n_valid < MIN_VALID_IN_BLOCK) {
      all_track[[length(all_track) + 1]] <- data.table(
        chrom = chr, global_window_id = wid_center,
        start_bp = grid$start_bp[wi], end_bp = grid$end_bp[wi],
        block_start_bp = grid$start_bp[b_start],
        block_end_bp = grid$end_bp[b_end],
        block_n_windows = length(block_idx),
        block_n_valid = n_valid,
        threeband_score = NA_real_,
        middle_fraction = NA_real_,
        middle_continuity = NA_real_,
        middle_continuity_fuzzy = NA_real_,
        combined_score = NA_real_,
        mean_s2_combined = NA_real_,
        mean_s2_middle_stability = NA_real_,
        snake4_status = "FAIL"
      )
      next
    }

    valid_wids <- block_wids[has_band]

    # ── Curve A: 3-band support ────────────────────────────────────
    s2_scores <- vapply(valid_wids, function(w) {
      sc <- score_by_wid[[as.character(w)]]
      if (!is.null(sc) && is.finite(sc$combined_score)) sc$combined_score else 0
    }, numeric(1))
    mean_s2_combined <- mean(s2_scores, na.rm = TRUE)
    band_frac <- n_valid / length(block_idx)
    threeband_score <- 0.70 * mean_s2_combined + 0.30 * band_frac

    # ── Curve B: middle-stripe fraction ────────────────────────────
    n_samples_max <- max(vapply(valid_wids, function(w) {
      length(band_by_wid[[as.character(w)]]$middle_weight)
    }, integer(1)))

    sample_middle_count <- integer(n_samples_max)
    for (w in valid_wids) {
      mid <- band_by_wid[[as.character(w)]]$middle_idx
      if (length(mid) > 0 && max(mid) <= n_samples_max) {
        sample_middle_count[mid] <- sample_middle_count[mid] + 1L
      }
    }
    n_consistent_middle <- sum(sample_middle_count >= n_valid * 0.50)
    middle_fraction <- n_consistent_middle / n_samples_max

    # ── Curve C: middle-core membership continuity ─────────────────
    if (length(valid_wids) >= 2) {
      consec_hard <- numeric(length(valid_wids) - 1)
      consec_fuzzy <- numeric(length(valid_wids) - 1)

      for (k in seq_len(length(valid_wids) - 1)) {
        b_a <- band_by_wid[[as.character(valid_wids[k])]]
        b_b <- band_by_wid[[as.character(valid_wids[k + 1])]]
        consec_hard[k] <- jaccard(b_a$middle_idx, b_b$middle_idx)

        if (length(b_a$middle_weight) == length(b_b$middle_weight) &&
            length(b_a$middle_weight) > 0) {
          consec_fuzzy[k] <- fuzzy_jaccard(b_a$middle_weight, b_b$middle_weight)
        } else {
          consec_fuzzy[k] <- consec_hard[k]
        }
      }
      middle_continuity       <- min(consec_hard)
      middle_continuity_fuzzy <- min(consec_fuzzy)
    } else {
      middle_continuity <- NA_real_
      middle_continuity_fuzzy <- NA_real_
    }

    # Snake 2 middle stability (adjacent-window, for reference)
    s2_mid_stab <- vapply(valid_wids, function(w) {
      sc <- score_by_wid[[as.character(w)]]
      if (!is.null(sc) && is.finite(sc$mean_middle_stability)) sc$mean_middle_stability else NA_real_
    }, numeric(1))
    mean_s2_mid_stab <- mean(s2_mid_stab, na.rm = TRUE)

    # ── Combined score ─────────────────────────────────────────────
    a_norm <- min(1, max(0, threeband_score))
    b_norm <- min(1, max(0, middle_fraction / 0.35))
    c_norm <- if (is.finite(middle_continuity)) min(1, max(0, middle_continuity)) else 0
    combined <- W_THREEBAND * a_norm + W_MIDFRAC * b_norm + W_MIDCONT * c_norm

    # ── Status ─────────────────────────────────────────────────────
    status <- if (threeband_score >= THREEBAND_PASS &&
                  middle_fraction >= MIDDLE_FRAC_PASS &&
                  is.finite(middle_continuity) && middle_continuity >= MIDDLE_CONT_PASS) {
      "PASS"
    } else if (threeband_score >= THREEBAND_PASS &&
               middle_fraction >= MIDDLE_FRAC_WEAK &&
               is.finite(middle_continuity) && middle_continuity < MIDDLE_CONT_PASS) {
      "PARTIAL"
    } else if (threeband_score >= THREEBAND_WEAK &&
               middle_fraction >= MIDDLE_FRAC_WEAK) {
      "WEAK"
    } else {
      "FAIL"
    }

    all_track[[length(all_track) + 1]] <- data.table(
      chrom = chr,
      global_window_id = wid_center,
      start_bp = grid$start_bp[wi],
      end_bp = grid$end_bp[wi],
      block_start_bp = grid$start_bp[b_start],
      block_end_bp = grid$end_bp[b_end],
      block_n_windows = length(block_idx),
      block_n_valid = n_valid,
      threeband_score = round(threeband_score, 4),
      middle_fraction = round(middle_fraction, 4),
      middle_continuity = round(middle_continuity, 4),
      middle_continuity_fuzzy = round(middle_continuity_fuzzy, 4),
      combined_score = round(combined, 4),
      mean_s2_combined = round(mean_s2_combined, 4),
      mean_s2_middle_stability = round(mean_s2_mid_stab, 4),
      snake4_status = status
    )

    # Save middle set for regime clustering
    center_band <- band_by_wid[[as.character(wid_center)]]
    if (!is.null(center_band)) {
      all_middle_sets[[length(all_middle_sets) + 1]] <- list(
        chrom = chr,
        global_window_id = wid_center,
        middle_idx = center_band$middle_idx,
        middle_weight = center_band$middle_weight
      )
    }
  }

  # ── Per-chromosome summary ─────────────────────────────────────────
  chr_rows <- all_track[vapply(all_track, function(x) x$chrom[1] == chr, logical(1))]
  chr_statuses <- vapply(chr_rows, function(x) x$snake4_status, character(1))

  all_summary[[length(all_summary) + 1]] <- data.table(
    chrom = chr,
    n_windows = n_win,
    n_scored = length(chr_statuses),
    n_pass = sum(chr_statuses == "PASS"),
    n_partial = sum(chr_statuses == "PARTIAL"),
    n_weak = sum(chr_statuses == "WEAK"),
    n_fail = sum(chr_statuses == "FAIL")
  )

  message("[S4] ", chr, ": PASS=", sum(chr_statuses == "PASS"),
          " PARTIAL=", sum(chr_statuses == "PARTIAL"),
          " WEAK=", sum(chr_statuses == "WEAK"),
          " FAIL=", sum(chr_statuses == "FAIL"))
}

# =============================================================================
# REGIME CLUSTERING (per chromosome)
# =============================================================================
# Cluster windows by shared middle-band composition into recurring regimes.
# Greedy: seed from highest-scoring windows, absorb neighbors with Jaccard >= thresh.
# This answers: "which regions share the same middle-stripe samples?"

message("\n[S4] Regime clustering...")

track_dt <- if (length(all_track) > 0) rbindlist(all_track, fill = TRUE) else {
  data.table(chrom = character(), global_window_id = integer(),
             snake4_status = character())
}

track_dt[, middle_regime_id := NA_integer_]
regime_summaries <- list()

if (length(all_middle_sets) >= MIN_REGIME_WINDOWS) {
  ms_wids <- vapply(all_middle_sets, function(x) x$global_window_id, integer(1))
  ms_chroms <- vapply(all_middle_sets, function(x) x$chrom, character(1))

  for (chr in unique(ms_chroms)) {
    chr_mask <- ms_chroms == chr
    chr_ms <- all_middle_sets[chr_mask]
    n_ms <- length(chr_ms)
    if (n_ms < MIN_REGIME_WINDOWS) next

    assigned <- rep(FALSE, n_ms)
    regime_id <- 0L

    # Sort by combined_score descending (best windows seed regimes)
    chr_wids <- vapply(chr_ms, function(x) x$global_window_id, integer(1))
    chr_scores <- track_dt[chrom == chr & global_window_id %in% chr_wids,
                           .(global_window_id, combined_score)]
    score_order <- chr_scores[order(-combined_score)]$global_window_id
    ms_order <- match(score_order, chr_wids)
    ms_order <- ms_order[!is.na(ms_order)]

    for (seed_i in ms_order) {
      if (assigned[seed_i]) next

      seed_mid <- chr_ms[[seed_i]]$middle_idx
      if (length(seed_mid) < 3) next

      members <- seed_i
      assigned[seed_i] <- TRUE

      for (j in seq_len(n_ms)) {
        if (assigned[j]) next
        j_mid <- chr_ms[[j]]$middle_idx
        if (length(j_mid) == 0) next
        if (jaccard(seed_mid, j_mid) >= REGIME_JACCARD_THRESH) {
          members <- c(members, j)
          assigned[j] <- TRUE
        }
      }

      if (length(members) >= MIN_REGIME_WINDOWS) {
        regime_id <- regime_id + 1L
        member_wids <- chr_wids[members]
        track_dt[chrom == chr & global_window_id %in% member_wids,
                 middle_regime_id := regime_id]

        # Regime core: samples in middle for >=50% of regime windows
        all_mid_idx <- unlist(lapply(chr_ms[members], function(x) x$middle_idx))
        mid_tab <- table(all_mid_idx)
        core_samples <- as.integer(names(mid_tab[mid_tab >= length(members) * 0.50]))

        regime_summaries[[length(regime_summaries) + 1]] <- data.table(
          chrom = chr,
          regime_id = regime_id,
          n_windows = length(members),
          start_bp = min(track_dt[chrom == chr & global_window_id %in% member_wids]$start_bp,
                         na.rm = TRUE),
          end_bp = max(track_dt[chrom == chr & global_window_id %in% member_wids]$end_bp,
                       na.rm = TRUE),
          n_core_middle_samples = length(core_samples),
          core_middle_samples_csv = paste(core_samples, collapse = ","),
          mean_combined = round(mean(
            track_dt[chrom == chr & global_window_id %in% member_wids]$combined_score,
            na.rm = TRUE), 4)
        )
      }
    }

    # Unassigned windows get regime 0
    unassigned_wids <- chr_wids[!assigned]
    if (length(unassigned_wids) > 0) {
      track_dt[chrom == chr & global_window_id %in% unassigned_wids &
               is.na(middle_regime_id), middle_regime_id := 0L]
    }

    message("[S4] ", chr, ": ", regime_id, " regimes identified")
  }
}

regime_dt <- if (length(regime_summaries) > 0) rbindlist(regime_summaries) else {
  data.table(chrom = character(), regime_id = integer())
}

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

summary_dt <- if (length(all_summary) > 0) rbindlist(all_summary) else data.table()

f1 <- file.path(outdir, "snake4_track.tsv.gz")
f2 <- file.path(outdir, "snake4_regimes.tsv.gz")
f3 <- file.path(outdir, "snake4_summary.tsv")

fwrite(track_dt, f1, sep = "\t")
fwrite(regime_dt, f2, sep = "\t")
fwrite(summary_dt, f3, sep = "\t")

message("\n[DONE] Snake 4 rolling 3-band support + regime identity complete")
message("  ", f1)
message("  ", f2, " (", nrow(regime_dt), " regimes)")
message("  ", f3)

if (nrow(summary_dt) > 0) {
  message("\n[S4] Total: PASS=", sum(summary_dt$n_pass),
          " PARTIAL=", sum(summary_dt$n_partial),
          " WEAK=", sum(summary_dt$n_weak),
          " FAIL=", sum(summary_dt$n_fail))
}
if (nrow(regime_dt) > 0) {
  message("[S4] Regimes: ", nrow(regime_dt), " across ",
          uniqueN(regime_dt$chrom), " chromosomes")
  message("[S4] Largest regime: ", max(regime_dt$n_windows), " windows, ",
          max(regime_dt$n_core_middle_samples), " core middle samples")
}