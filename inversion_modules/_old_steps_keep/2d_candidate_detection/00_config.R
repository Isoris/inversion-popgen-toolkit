#!/usr/bin/env Rscript
# ============================================================================
# 00_config.R — Central configuration for phase 2 / 2d candidate detection
# ============================================================================
# ALL thresholds, paths, and parameters live here.  Nothing is hidden.
# Source this file before any module.
# ============================================================================

# ---- HPC paths ----
CFG <- list()

CFG$BASE        <- "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
CFG$CODEBASE    <- file.path(CFG$BASE, "inversion_codebase_v8.5")
CFG$SNAKES      <- file.path(CFG$CODEBASE, "MODULE_5A2_Discovery_Core/snakes")
CFG$PRECOMP     <- file.path(CFG$SNAKES, "precomp")
CFG$SIM_MATS    <- file.path(CFG$PRECOMP, "sim_mats")
CFG$TRIANGLES   <- file.path(CFG$SNAKES, "triangles_v3")
CFG$OUTDIR      <- file.path(CFG$CODEBASE, "MODULE_5A2_Discovery_Core/inv_detect_out_v9.3")
CFG$SV_FILE     <- file.path(CFG$BASE, "delly_sv/merged_sv_calls.tsv")
CFG$GHSL_DIR    <- file.path(CFG$SNAKES, "ghsl_v6")   # chat 14: v5 → v6

# ---- General ----
CFG$WINDOW_SIZE_BP <- 50000L   # bin-to-bp conversion factor
CFG$N_SAMPLES      <- 226L     # total samples in the cohort

# ---- Staircase detector (01) ----
CFG$STAIR_MIN_BLOCK_WIDTH  <- 5L     # minimum block width in bins
CFG$STAIR_MIN_DROP         <- 0.03   # minimum step height (sim units) to call boundary
CFG$STAIR_PERSIST_N        <- 3L     # sustained drop must persist this many bins
CFG$STAIR_SMOOTH_SPAN      <- 3L     # running median span for depth smoothing
CFG$STAIR_EMA_ALPHA        <- 0.05   # plateau tracking rate (0=frozen, 1=instant)
CFG$STAIR_VOTE_BLUR        <- 2L     # ±bins for vote jitter tolerance

# ---- Matrix transforms (07) ----
CFG$MTX_DISTCORR_ENABLED  <- TRUE
CFG$MTX_LOCALNORM_ENABLED <- TRUE
CFG$MTX_LOCALNORM_WINDOW  <- 50L   # ±50 bins for local z-score
CFG$MTX_DENOISED_ENABLED  <- TRUE
CFG$MTX_DENOISED_LAMBDA   <- 1.0   # fused lasso lambda (use proxy if gfl unavailable)
CFG$MTX_RESIDBG_ENABLED   <- TRUE
CFG$MTX_RESIDBG_RANK      <- 3L    # SVD rank for low-rank background
CFG$MTX_RESIDBG_KERNEL    <- 200L  # large-kernel smooth window (alternative method)
CFG$MTX_EDGE_ENABLED      <- TRUE
CFG$MTX_SUPPORT_ENABLED   <- TRUE
CFG$MTX_SUPPORT_QUANTILE  <- 0.75  # adaptive threshold quantile

# ---- Block scoring (08) ----
CFG$BLOCK_NEAR_FRAC  <- 0.2   # fraction of block width for "near" zone
CFG$BLOCK_FAR_FRAC   <- 0.2   # fraction of block width for "far" zone
CFG$BLOCK_EDGE_DEPTH <- 3L    # bins inside/outside boundary for edge contrast

# ---- NN sweep / interval tree (09) ----
# 2026-04-17 FIX 21: scales expanded so D09's persistence-barcode tree
# can reach its INVERSION classifier threshold (nn_birth >= 200). Was
# c(0, 20, 40, 80) — silently broken because precomp didn't produce the
# nn-smoothed sim_mats anyway. See STEP_C01a_precompute.R FIX 21 header.
CFG$NN_SCALES_DEFAULT  <- c(0L, 20L, 40L, 80L, 120L, 160L, 200L, 240L, 320L)
CFG$NN_SCALES_EXTENDED <- seq(20L, 320L, by=20L)  # full sweep for diagnostics
CFG$NN_SWEEP_COARSE    <- seq(100L, 4000L, by=100L)  # research-only
CFG$NN_SWEEP_FINE      <- seq(20L, 100L, by=20L)     # fine range
CFG$NN_OVERLAP_THRESH  <- 0.70   # reciprocal overlap for block matching

# ---- Consensus (10) ----
CFG$CONSENSUS_MIN_VARIANTS <- 3L  # blocks must appear in ≥3 matrix variants
CFG$CONSENSUS_OVERLAP      <- 0.50  # overlap threshold for cross-variant matching

# ---- Landscape classifier (14) ----
CFG$LANDSCAPE_STRONG_SQ    <- 0.65  # squareness threshold for strong_inversion
CFG$LANDSCAPE_DIFFUSE_SQ   <- 0.40  # squareness threshold for diffuse_inversion
CFG$LANDSCAPE_FAMILY_SQ    <- 0.35  # below this + patchy = family_ld or diffuse_diag
CFG$LANDSCAPE_FIXED_HEIGHT <- 0.85  # height threshold for nested_fixed
CFG$LANDSCAPE_CONTRAST_SCALE <- 0.15  # confidence: contrast / this → 0-1
CFG$LANDSCAPE_VOTE_SCALE     <- 50    # confidence: votes / this → 0-1
CFG$LANDSCAPE_STEP_SCALE     <- 0.10  # confidence: step_sharpness / this → 0-1

# ---- Evidence modules ----
CFG$NN_SURVIVAL_RATIO  <- 0.5    # signal must keep 50% of nn0 ratio
CFG$SV_MAX_DIST_KB     <- 500    # max distance to match SV breakpoint
CFG$ICA_ENABLED        <- TRUE   # attempt ICA if fastICA available
CFG$PEEL_MAX_FRAC      <- 0.3    # never remove more than 30% of samples
CFG$PEEL_LOCAL_KIN_THRESH <- 0.7 # PC1 correlation for L1b chr-local peel
CFG$PEEL_LOCAL_K       <- 5L     # k clusters for L2 block co-segregation
CFG$PEEL_KIN_THRESHOLD <- 0.05   # ngsRelate rab threshold for L1

# ---- Known positives/negatives for calibration ----
CFG$KNOWN_POSITIVES <- list(
  list(chr="C_gar_LG01", name="I57",     start_mb=24,   end_mb=28,   notes="nested inversion system"),
  list(chr="C_gar_LG01", name="I70",     start_mb=29,   end_mb=33.5, notes="large clean inversion"),
  list(chr="C_gar_LG25", name="cand58",  start_mb=15.1, end_mb=15.7, notes="textbook 550kb inversion"),
  list(chr="C_gar_LG10", name="cand20",  start_mb=13,   end_mb=16,   notes="2.7 Mb candidate")
)
CFG$KNOWN_NEGATIVES <- list(
  list(chr="C_gar_LG01", name="I86",     start_mb=34,   end_mb=38,   notes="family LD, patchy")
)

# ---- Helper: convert Mb to bins and vice versa ----
mb_to_bin <- function(mb, window_size_bp=CFG$WINDOW_SIZE_BP) {
  as.integer(ceiling(mb * 1e6 / window_size_bp))
}
bin_to_mb <- function(bin, window_size_bp=CFG$WINDOW_SIZE_BP) {
  (bin - 1) * window_size_bp / 1e6
}
bin_to_bp <- function(bin, window_size_bp=CFG$WINDOW_SIZE_BP) {
  (bin - 1) * window_size_bp
}

cat("Config loaded: inv_detect v9.3\n")
