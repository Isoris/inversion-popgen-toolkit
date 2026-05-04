#!/usr/bin/env Rscript
# =============================================================================
# engine_b_smoke_test.R — fail-fast self-test before long compute
# =============================================================================
# Problem being solved:
#   A silent Engine B failure (binary not compiled, dispatcher misconfigured,
#   Q matrix missing, BAM list wrong) turns a 6-hour precomp into 6 hours of
#   NA values. You don't notice until the catalog is empty. This takes 5
#   seconds to catch.
#
# Usage:
#   source("utils/engine_b_smoke_test.R")
#   engine_b_ok <- run_engine_b_smoke_test(
#     test_chr = "LG01",     # any chromosome that has a BAM
#     test_start = 1000000,
#     test_end = 1100000,    # 100kb region, cheap
#     min_samples_per_group = 5
#   )
#   if (!engine_b_ok) {
#     stop("Engine B smoke test failed. Fix before running precomp.")
#   }
#
# Exits fast (≤5s typical) with a clear error if anything is wrong.
# Exits silently if everything works.
#
# What it actually tests:
#   1. load_bridge.R is sourceable
#   2. get_region_stats / get_Q functions exist and are callable
#   3. Sample registry has ≥2 populated groups (or falls back to first 10
#      samples split in half)
#   4. One Fst call returns a finite number (not NA)
#   5. If Q matrix is loaded: one get_Q call returns a matrix with the
#      expected K ancestry components
#
# If the smoke test passes, Engine B is REACHABLE. Whether every downstream
# call will produce good numbers is a different question — that requires
# real data. But the most common failure modes (unreachable, miswired, no
# data) are caught here.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

run_engine_b_smoke_test <- function(
  test_chr = NULL,
  test_start = NULL,
  test_end = NULL,
  min_samples_per_group = 5L,
  verbose = TRUE
) {
  msg <- function(...) if (verbose) message("[smoke_test] ", ...)
  fail <- function(reason) {
    message("[smoke_test] ✗ FAILED: ", reason)
    return(FALSE)
  }

  t_start <- Sys.time()
  msg("starting Engine B smoke test")

  # ── Check 1: load_bridge.R sourceable ───────────────────────────────
  bridge_path <- NULL
  for (p in c("utils/load_bridge.R", "../utils/load_bridge.R",
              Sys.getenv("LOAD_BRIDGE", ""),
              file.path(Sys.getenv("BASE", ""),
                        "inversion-popgen-toolkit/utils/load_bridge.R"))) {
    if (nzchar(p) && file.exists(p)) { bridge_path <- p; break }
  }
  if (is.null(bridge_path)) {
    return(fail("load_bridge.R not found in any standard location"))
  }
  msg("check 1: load_bridge.R at ", bridge_path)

  ok <- tryCatch({
    source(bridge_path); TRUE
  }, error = function(e) {
    msg("error sourcing load_bridge.R: ", conditionMessage(e))
    FALSE
  })
  if (!ok) return(fail("load_bridge.R source error"))

  # ── Check 2: key functions exist ────────────────────────────────────
  if (!exists("get_region_stats", mode = "function")) {
    return(fail("get_region_stats() not available after sourcing bridge"))
  }
  if (!exists(".bridge_available")) {
    return(fail(".bridge_available flag not set; bridge loaded incorrectly"))
  }
  if (isFALSE(.bridge_available)) {
    return(fail(".bridge_available = FALSE; Engine B binary likely not compiled. ",
                "Run: cd unified_ancestry/src && make && cd ../engines/fst_dxy && make"))
  }
  msg("check 2: get_region_stats available, .bridge_available = TRUE")

  # ── Check 3: sample groups available ────────────────────────────────
  test_groups <- NULL
  if (exists("reg") && !is.null(reg) &&
      !is.null(reg$list_groups) && nrow(reg$list_groups()) >= 2) {
    groups_dt <- reg$list_groups()
    # Pick two groups with ≥ min_samples_per_group members
    valid <- groups_dt[n >= min_samples_per_group]
    if (nrow(valid) >= 2) {
      g1_id <- valid$group_id[1]; g2_id <- valid$group_id[2]
      g1 <- reg$get_group(g1_id); g2 <- reg$get_group(g2_id)
      test_groups <- list(g1 = g1, g2 = g2)
      msg("check 3: using registry groups ", g1_id, " (n=", length(g1),
          ") and ", g2_id, " (n=", length(g2), ")")
    }
  }

  if (is.null(test_groups)) {
    # Fallback: split first 10 samples in half
    sample_list_path <- NULL
    for (p in c(
      Sys.getenv("SAMPLES_IND", ""),
      file.path(Sys.getenv("BASE", "."),
                "het_roh/01_inputs_check/samples.ind"),
      file.path(Sys.getenv("BASE", "."),
                "popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv")
    )) {
      if (nzchar(p) && file.exists(p)) { sample_list_path <- p; break }
    }
    if (is.null(sample_list_path)) {
      return(fail("no registry groups AND no sample list found for fallback"))
    }
    all_samps <- trimws(readLines(sample_list_path))
    all_samps <- all_samps[nzchar(all_samps)]
    if (length(all_samps) < 2 * min_samples_per_group) {
      return(fail(paste0("only ", length(all_samps), " samples in ", sample_list_path,
                          "; need at least ", 2 * min_samples_per_group)))
    }
    n_per <- min_samples_per_group
    test_groups <- list(
      g1 = all_samps[1:n_per],
      g2 = all_samps[(n_per + 1):(2 * n_per)]
    )
    msg("check 3: fallback — split first ", 2 * n_per,
        " samples from ", sample_list_path)
  }

  # ── Check 4: test region defaults ───────────────────────────────────
  if (is.null(test_chr)) {
    # Pick any chromosome that's been defined; otherwise "LG01"
    test_chr <- if (!is.null(Sys.getenv("TEST_CHR_SMOKE", ""))) {
      Sys.getenv("TEST_CHR_SMOKE", "LG01")
    } else "LG01"
  }
  if (is.null(test_start)) test_start <- 1000000L
  if (is.null(test_end))   test_end   <- 1100000L  # 100 kb, cheap

  msg("check 4: test region = ", test_chr, ":", test_start, "-", test_end)

  # ── Check 5: one Fst call ───────────────────────────────────────────
  t0 <- Sys.time()
  fst_result <- tryCatch(
    get_region_stats(
      test_chr,
      start_bp = as.integer(test_start),
      end_bp = as.integer(test_end),
      what = "Fst",
      groups = test_groups
    ),
    error = function(e) {
      msg("exception in get_region_stats(Fst): ", conditionMessage(e))
      NULL
    }
  )
  t_fst <- round(as.numeric(Sys.time() - t0), 2)
  if (is.null(fst_result)) {
    return(fail("get_region_stats(Fst) threw an exception; Engine B is broken"))
  }

  # Extract the numeric Fst — the dispatcher returns a list
  fst_val <- NULL
  if (is.list(fst_result) && !is.null(fst_result$Fst)) {
    fst_val <- as.numeric(fst_result$Fst[[1]])
  } else if (is.numeric(fst_result)) {
    fst_val <- as.numeric(fst_result[1])
  }
  if (is.null(fst_val) || !is.finite(fst_val)) {
    return(fail(paste0("get_region_stats(Fst) returned non-finite: ",
                        capture.output(str(fst_result))[1],
                        ". Check BAM list, Q matrix, region bounds.")))
  }
  msg("check 5: Fst = ", round(fst_val, 4), " (", t_fst, "s) — finite, reachable")

  # ── Check 6 (optional): one get_Q call ──────────────────────────────
  if (exists("get_Q", mode = "function")) {
    t0 <- Sys.time()
    q_result <- tryCatch(
      get_Q(test_chr, start_bp = as.integer(test_start),
                       end_bp = as.integer(test_end)),
      error = function(e) {
        msg("exception in get_Q(): ", conditionMessage(e))
        NULL
      }
    )
    t_q <- round(as.numeric(Sys.time() - t0), 2)
    if (is.null(q_result)) {
      msg("check 6: get_Q returned NULL — Q computation unavailable, but Fst works")
    } else if (is.matrix(q_result) || is.data.frame(q_result)) {
      K <- ncol(q_result)
      msg("check 6: get_Q returned ", nrow(q_result), "x", K, " matrix (", t_q, "s)")
      # Sanity: proportions should sum to ≈ 1 per row
      row_sums <- rowSums(q_result, na.rm = TRUE)
      if (any(abs(row_sums - 1) > 0.1, na.rm = TRUE)) {
        msg("  WARNING: row sums deviate from 1 by > 0.1 — Q matrix may be malformed")
      }
    } else {
      msg("check 6: get_Q returned unexpected type: ", class(q_result)[1])
    }
  } else {
    msg("check 6: get_Q not present — skipped (Q scan unavailable, Fst OK)")
  }

  t_total <- round(as.numeric(Sys.time() - t_start), 2)
  message("[smoke_test] ✓ Engine B OK (", t_total, "s total)")
  return(TRUE)
}

# ── Convenience: if sourced with RUN_SMOKE=1, run immediately ──
if (nzchar(Sys.getenv("RUN_SMOKE_TEST", ""))) {
  ok <- run_engine_b_smoke_test()
  if (!ok) {
    message("")
    message("EXAMPLES of what to check:")
    message("  1. Is Engine B compiled?   cd unified_ancestry/src && make")
    message("  2. Is Q matrix present?    ls ${PRECOMP_OUT}/*.qopt")
    message("  3. Is BAM list correct?    head ${BAM_LIST}")
    message("  4. Is dispatcher wired?    grep get_region_stats utils/region_stats_dispatcher*.R")
    message("")
    quit(status = 1)
  }
}
