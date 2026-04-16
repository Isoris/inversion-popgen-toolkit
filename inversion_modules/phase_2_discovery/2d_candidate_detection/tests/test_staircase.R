#!/usr/bin/env Rscript
# ============================================================================
# test_staircase.R — Synthetic-data tests for the staircase detector
# ============================================================================
# Builds synthetic sim_mats with known blocks and verifies detection.
# Run this BEFORE testing on real data.
#
# Usage: Rscript 16_test_staircase.R
# ============================================================================

suppressPackageStartupMessages(library(data.table))

# Source staircase (with fallback config)
CFG <- list(STAIR_MIN_BLOCK_WIDTH = 5L, STAIR_MIN_DROP = 0.03,
            STAIR_PERSIST_N = 3L, STAIR_SMOOTH_SPAN = 3L,
            STAIR_EMA_ALPHA = 0.05, STAIR_VOTE_BLUR = 2L,
            WINDOW_SIZE_BP = 50000L)
bin_to_bp <- function(b, w = 50000L) (b - 1) * w
bin_to_mb <- function(b, w = 50000L) (b - 1) * w / 1e6

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(script_dir, "STEP_D01_staircase_boundaries.R"))

n_pass <- 0L; n_fail <- 0L; n_total <- 0L

check <- function(test_name, condition) {
  n_total <<- n_total + 1L
  if (isTRUE(condition)) {
    n_pass <<- n_pass + 1L
    cat("  PASS:", test_name, "\n")
  } else {
    n_fail <<- n_fail + 1L
    cat("  FAIL:", test_name, "\n")
  }
}

# ---- Helper: build synthetic sim_mat ----
make_sim <- function(N, blocks, bg = 0.2, noise_sd = 0.01) {
  smat <- matrix(bg, nrow = N, ncol = N)
  diag(smat) <- 1

  # Add distance decay to background (family LD-like)
  for (i in 1:N) for (j in 1:N) {
    d <- abs(i - j)
    smat[i, j] <- smat[i, j] + max(0, 0.1 - d * 0.001)  # mild decay
  }

  # Add blocks
  for (b in blocks) {
    for (i in b$start:b$end) {
      for (j in b$start:b$end) {
        smat[i, j] <- b$height
      }
    }
  }

  # Add noise
  noise <- matrix(rnorm(N * N, 0, noise_sd), N, N)
  noise <- (noise + t(noise)) / 2
  smat <- smat + noise
  smat <- pmax(0, pmin(1, smat))
  diag(smat) <- 1

  smat
}


# ============================================================================
# TEST 1: Single clean block
# ============================================================================
cat("\n=== TEST 1: Single clean block ===\n")

smat1 <- make_sim(200, list(
  list(start = 50, end = 100, height = 0.8)
))

r1 <- detect_blocks_staircase(smat1)

check("blocks found", nrow(r1$blocks) >= 1)
check("block near 50-100",
      any(r1$blocks$start <= 55 & r1$blocks$end >= 95))
check("no block at 150-200",
      !any(r1$blocks$start >= 140 & r1$blocks$end <= 200 & r1$blocks$height > 0.5))


# ============================================================================
# TEST 2: Two separate blocks
# ============================================================================
cat("\n=== TEST 2: Two separate blocks ===\n")

smat2 <- make_sim(300, list(
  list(start = 30, end = 80, height = 0.8),
  list(start = 180, end = 230, height = 0.7)
))

r2 <- detect_blocks_staircase(smat2)

check("at least 2 blocks", nrow(r2$blocks) >= 2)
check("block near 30-80",
      any(r2$blocks$start <= 35 & r2$blocks$end >= 75))
check("block near 180-230",
      any(r2$blocks$start <= 185 & r2$blocks$end >= 225))


# ============================================================================
# TEST 3: Nested blocks (inner brighter than outer)
# ============================================================================
cat("\n=== TEST 3: Nested blocks ===\n")

smat3 <- make_sim(200, list(
  list(start = 40, end = 120, height = 0.6),
  list(start = 60, end = 100, height = 0.85)
))

r3 <- detect_blocks_staircase(smat3)

check("at least 2 blocks", nrow(r3$blocks) >= 2)
# Should have parent-child relationship
if (nrow(r3$blocks) >= 2) {
  has_nesting <- any(!is.na(r3$blocks$parent_id))
  check("parent-child nesting detected", has_nesting)
}


# ============================================================================
# TEST 4: No blocks (pure noise)
# ============================================================================
cat("\n=== TEST 4: Pure noise (no blocks) ===\n")

smat4 <- make_sim(200, list(), noise_sd = 0.03)

r4 <- detect_blocks_staircase(smat4)

check("few or no blocks", nrow(r4$blocks) <= 3)


# ============================================================================
# TEST 5: Cascading steps (the 0.85→0.60→0.20 bug)
# ============================================================================
cat("\n=== TEST 5: Cascading steps (regression test) ===\n")

smat5 <- make_sim(300, list(
  list(start = 50, end = 152, height = 0.6),   # outer rare
  list(start = 52, end = 150, height = 0.85)    # inner common
))

r5 <- detect_blocks_staircase(smat5)

# The inner block should be detected
check("inner block found",
      any(r5$blocks$start <= 55 & r5$blocks$end >= 145 & r5$blocks$height > 0.7))

# Check step table has steps at both boundaries
if (nrow(r5$step_table) > 0) {
  steps_near_150 <- r5$step_table[step_position >= 148 & step_position <= 155]
  check("steps found near inner boundary (~150)",
        nrow(steps_near_150) > 0)
}


# ============================================================================
# TEST 6: Tiny block inside large (the 3-window bug)
# ============================================================================
cat("\n=== TEST 6: Tiny block inside large ===\n")

smat6 <- make_sim(200, list(
  list(start = 40, end = 120, height = 0.6),
  list(start = 75, end = 78, height = 0.9)   # 4 windows, very bright
))

r6 <- detect_blocks_staircase(smat6)

# The outer block should be found
check("outer block found",
      any(r6$blocks$start <= 45 & r6$blocks$end >= 115))

# The tiny inner block: check if pass 2 catches it
tiny_found <- any(r6$blocks$start >= 73 & r6$blocks$end <= 80 & r6$blocks$width <= 6)
check("tiny inner block found (pass 2)", tiny_found)


# ============================================================================
# TEST 7: Family LD flag
# ============================================================================
cat("\n=== TEST 7: Family LD detection ===\n")

smat7 <- make_sim(200, list(), bg = 0.2, noise_sd = 0.01)
# Make some windows have elevated background (family LD)
for (i in 80:120) {
  smat7[i, ] <- smat7[i, ] + 0.15
  smat7[, i] <- smat7[, i] + 0.15
}

r7 <- detect_blocks_staircase(smat7)

if ("family_ld" %in% names(r7$vote_profile)) {
  n_ld <- sum(r7$vote_profile$family_ld, na.rm = TRUE)
  check("family LD windows detected", n_ld >= 20)
  check("family LD windows in right range",
        any(r7$vote_profile$family_ld[80:120]))
}


# ============================================================================
# SUMMARY
# ============================================================================
cat("\n==========================================\n")
cat("RESULTS:", n_pass, "passed,", n_fail, "failed, of", n_total, "tests\n")
cat("==========================================\n")

if (n_fail > 0) {
  cat("SOME TESTS FAILED — check output above\n")
  quit(status = 1)
} else {
  cat("ALL TESTS PASSED\n")
}
