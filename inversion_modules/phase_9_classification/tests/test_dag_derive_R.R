#!/usr/bin/env Rscript
# =============================================================================
# test_dag_derive_R.R — fixture for chat-12 DAG rewrite of derive_R_from_regime
# =============================================================================
# Handoff DoD item 8: "Add a new fixture for the DAG shape check (verify on
# A-A-A-B-B-A and A-A-A-B-A-B-B-A patterns)."
#
# Also covers:
#   - single-run case (no deviation at all)
#   - all-different case (maximum alternation)
#   - empty/degenerate input
#   - interval filtering (windows outside candidate are dropped)
#   - window-size inference from pos_mid_mb
#   - R_fired gate semantics: fraction AND bp span must both clear threshold
#   - back-compat: as_regime_dt() unwraps the list to the per_sample table
#   - back-compat: combine_cohort_recomb() accepts the new list shape
#
# Run:  Rscript tests/test_dag_derive_R.R
# Pass: all CHECK lines print PASS; final line "ALL TESTS PASSED".
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ── Locate the library relative to this test file ───────────────────────────
this_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
  error = function(e) NULL
)
if (is.null(this_file) || !nzchar(this_file) || !file.exists(this_file)) {
  # Running via Rscript: use commandArgs
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  this_file <- if (length(file_arg)) normalizePath(file_arg, mustWork = FALSE) else NA_character_
}
repo_root <- if (!is.na(this_file)) {
  # pass 15: test is now at inversion_modules/phase_9_classification/tests/
  # so repo root is 3 levels up (was 2 in the old phase_4_postprocessing/tests/ layout)
  dirname(dirname(dirname(this_file)))
} else {
  getwd()
}
lib_path <- file.path(repo_root, "inversion_modules",
                      "phase_7_karyotype_groups", "proposal",
                      "lib_recomb_combination.R")
if (!file.exists(lib_path)) {
  # Fallback: running from the proposal/ directory or repo root
  alts <- c("lib_recomb_combination.R",
            "../inversion_modules/phase_7_karyotype_groups/proposal/lib_recomb_combination.R",
            file.path(getwd(), "inversion_modules", "phase_7_karyotype_groups",
                      "proposal", "lib_recomb_combination.R"))
  for (a in alts) if (file.exists(a)) { lib_path <- a; break }
}
if (!file.exists(lib_path)) stop("cannot locate lib_recomb_combination.R — cwd=", getwd())
source(lib_path)

# ── Micro assert helper ─────────────────────────────────────────────────────
n_pass <- 0L; n_fail <- 0L
CHECK <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    cat("  FAIL  ", label, "  [error: ", conditionMessage(e), "]\n", sep = "")
    n_fail <<- n_fail + 1L
    return(FALSE)
  })
  if (isTRUE(ok)) {
    cat("  PASS  ", label, "\n", sep = "")
    n_pass <<- n_pass + 1L
  } else {
    cat("  FAIL  ", label, "\n", sep = "")
    n_fail <<- n_fail + 1L
  }
}
section <- function(title) cat("\n── ", title, " ──\n", sep = "")

# ── Helper: build a regime_memberships table for one sample from a
#    pattern string like "AAABBA". Windows are placed at 1, 2, …, n Mb.
# ───────────────────────────────────────────────────────────────────────────
mk_memb <- function(pattern_by_sample, window_start_mb = 1.0,
                    window_step_mb = 0.1) {
  rows <- rbindlist(lapply(names(pattern_by_sample), function(sid) {
    pat <- strsplit(pattern_by_sample[[sid]], "")[[1]]
    n   <- length(pat)
    data.table(
      sample     = sid,
      pos_mid_mb = window_start_mb + (seq_len(n) - 1L) * window_step_mb,
      group_id   = pat
    )
  }))
  rows
}

# ─────────────────────────────────────────────────────────────────────────────
section("Pattern 1: A-A-A-B-B-A (the primary handoff case)")
# Old counter was going to call long_A -> long_B -> long_A if min_run=3, or drop
# the B-burst as noise if min_run=2 — both failure modes in the chat-12 spec.
# DAG expects:
#   rle(AAABBA)   = values(A,B,A)  lengths(3,2,1)
#   dominant      = A  (3+1 = 4 windows; the lone B-run totals 2)
#   n_distinct    = 2  (A, B)
#   n_edges       = 2  (A->B, B->A)
#   deviation_fraction      = 2/6 = 0.3333...
#   longest_deviation_windows = 2 (the contiguous B-run)
#   longest_deviation_bp    = 2 * window_size_bp
#     — with 0.1 Mb step, window_size_bp should infer to 100_000 → 200_000 bp
#   terminal_matches_start  = TRUE (starts A, ends A)
# ─────────────────────────────────────────────────────────────────────────────
memb1 <- mk_memb(list(S1 = "AAABBA"))
res1  <- derive_R_from_regime(memb1, interval = list(start_mb = 0, end_mb = 10),
                              min_deviation_fraction = 0.15,
                              min_deviation_bp       = 150000L)
ps1 <- res1$per_sample
CHECK("AAABBA: one row",            nrow(ps1) == 1L)
CHECK("AAABBA: dominant = A",       ps1$dominant_regime == "A")
CHECK("AAABBA: n_distinct_nodes=2", ps1$n_distinct_nodes == 2L)
CHECK("AAABBA: n_edges=2",          ps1$n_edges == 2L)
CHECK("AAABBA: n_regime_changes=2 (back-compat alias)",
      ps1$n_regime_changes == 2L)
CHECK("AAABBA: deviation_fraction ~ 1/3",
      abs(ps1$deviation_fraction - (2/6)) < 1e-9)
CHECK("AAABBA: longest_deviation_windows = 2",
      ps1$longest_deviation_windows == 2L)
CHECK("AAABBA: window_size inferred to 100_000 bp",
      res1$params$window_size_bp == 100000L)
CHECK("AAABBA: longest_deviation_bp = 200_000",
      ps1$longest_deviation_bp == 200000L)
CHECK("AAABBA: terminal_matches_start = TRUE",
      isTRUE(ps1$terminal_matches_start))
CHECK("AAABBA: R_fired (frac=0.33>=0.15 AND 200kb>=150kb)",
      isTRUE(ps1$R_fired))
CHECK("AAABBA: dag_nodes_compact = A->B->A",
      ps1$dag_nodes_compact == "A->B->A")
CHECK("AAABBA: dag_node_weights = c(3,2,1)",
      identical(ps1$dag_node_weights[[1]], c(3L, 2L, 1L)))

# ─────────────────────────────────────────────────────────────────────────────
section("Pattern 2: A-A-A-B-A-B-B-A (the alternating-noise case)")
# rle(AAABABBA) = values(A,B,A,B,A)  lengths(3,1,1,2,1)
#   A-windows total = 3+1+1 = 5
#   B-windows total = 1+2   = 3
#   dominant = A, dom_w = 5, n = 8, deviation_fraction = 3/8 = 0.375
#   longest non-dominant run = 2 (the middle BB)
#   n_distinct_nodes = 2, n_edges = 4
#   terminal_matches_start = TRUE (A ... A)
# With window_size_bp = 100_000 and min_deviation_bp = 150_000,
#   longest_dev_bp = 200_000 -> R_fired = TRUE (fraction 0.375 >= 0.15).
# With min_deviation_bp = 250_000, R_fired = FALSE despite high fraction —
# the bp gate protects against alternating-noise false positives of the kind
# the handoff specifically called out.
# ─────────────────────────────────────────────────────────────────────────────
memb2 <- mk_memb(list(S2 = "AAABABBA"))
res2a <- derive_R_from_regime(memb2, list(start_mb = 0, end_mb = 10),
                              min_deviation_fraction = 0.15,
                              min_deviation_bp       = 150000L)
ps2a <- res2a$per_sample
CHECK("AAABABBA: dominant = A",             ps2a$dominant_regime == "A")
CHECK("AAABABBA: deviation_fraction = 3/8", abs(ps2a$deviation_fraction - 3/8) < 1e-9)
CHECK("AAABABBA: longest_dev_windows = 2",  ps2a$longest_deviation_windows == 2L)
CHECK("AAABABBA: n_edges = 4",              ps2a$n_edges == 4L)
CHECK("AAABABBA: n_distinct_nodes = 2",     ps2a$n_distinct_nodes == 2L)
CHECK("AAABABBA: terminal_matches_start = TRUE",
      isTRUE(ps2a$terminal_matches_start))
CHECK("AAABABBA: R_fired at min_bp=150kb",  isTRUE(ps2a$R_fired))

res2b <- derive_R_from_regime(memb2, list(start_mb = 0, end_mb = 10),
                              min_deviation_fraction = 0.15,
                              min_deviation_bp       = 250000L)
CHECK("AAABABBA: R_fired = FALSE at min_bp=250kb (bp gate protects)",
      isFALSE(res2b$per_sample$R_fired))

# ─────────────────────────────────────────────────────────────────────────────
section("Single-run case AAAAAA")
# Only one regime visited, zero deviation.
# ─────────────────────────────────────────────────────────────────────────────
memb3 <- mk_memb(list(S3 = "AAAAAA"))
res3  <- derive_R_from_regime(memb3, list(start_mb = 0, end_mb = 10))
ps3   <- res3$per_sample
CHECK("AAAAAA: deviation_fraction = 0",     ps3$deviation_fraction == 0)
CHECK("AAAAAA: longest_dev_windows = 0",    ps3$longest_deviation_windows == 0L)
CHECK("AAAAAA: n_distinct_nodes = 1",       ps3$n_distinct_nodes == 1L)
CHECK("AAAAAA: n_edges = 0",                ps3$n_edges == 0L)
CHECK("AAAAAA: terminal_matches_start = TRUE",
      isTRUE(ps3$terminal_matches_start))
CHECK("AAAAAA: R_fired = FALSE",            isFALSE(ps3$R_fired))

# ─────────────────────────────────────────────────────────────────────────────
section("Boundary recombinant AAABBB (dominant-tie broken, not-matching terminals)")
# rle = (A,B) lengths (3,3). Dominant: both tied at 3 windows — our
# tie-break rule is "first-occurrence wins", so dominant = A.
# deviation_fraction = 3/6 = 0.5, longest_dev = 3, terminals do NOT match.
# ─────────────────────────────────────────────────────────────────────────────
memb4 <- mk_memb(list(S4 = "AAABBB"))
res4  <- derive_R_from_regime(memb4, list(start_mb = 0, end_mb = 10))
ps4   <- res4$per_sample
CHECK("AAABBB: dominant = A (first-occurrence tie-break)",
      ps4$dominant_regime == "A")
CHECK("AAABBB: deviation_fraction = 0.5",   ps4$deviation_fraction == 0.5)
CHECK("AAABBB: terminal_matches_start = FALSE",
      isFALSE(ps4$terminal_matches_start))

# ─────────────────────────────────────────────────────────────────────────────
section("Cohort aggregation with three samples of mixed shapes")
# ─────────────────────────────────────────────────────────────────────────────
memb5 <- mk_memb(list(
  SA = "AAAAAA",     # no deviation
  SB = "AAABBA",     # deviation_fraction 1/3
  SC = "AAAAAB"      # tiny tail; dev_frac = 1/6 = 0.1667
))
res5 <- derive_R_from_regime(memb5, list(start_mb = 0, end_mb = 10),
                             min_deviation_fraction = 0.15,
                             min_deviation_bp       = 150000L)
ps5 <- res5$per_sample[order(sample_id)]
CHECK("cohort: 3 samples present",
      nrow(ps5) == 3L && setequal(ps5$sample_id, c("SA","SB","SC")))
CHECK("cohort: SA R_fired = FALSE",
      isFALSE(ps5[sample_id == "SA"]$R_fired))
CHECK("cohort: SB R_fired = TRUE (0.33 frac, 200kb)",
      isTRUE(ps5[sample_id == "SB"]$R_fired))
CHECK("cohort: SC R_fired = FALSE (frac 0.17 passes, longest_bp=100kb fails)",
      isFALSE(ps5[sample_id == "SC"]$R_fired))
CHECK("cohort: n_samples_total = 3",        res5$cohort$n_samples_total == 3L)
CHECK("cohort: n_samples_with_deviation = 2",
      res5$cohort$n_samples_with_deviation == 2L)
CHECK("cohort: fraction_samples_R_fired = 1/3",
      abs(res5$cohort$fraction_samples_R_fired - 1/3) < 1e-9)
CHECK("cohort: dominant_regime_cohort = A",
      res5$cohort$dominant_regime_cohort == "A")
CHECK("cohort: dominant_regime_support = 1",
      res5$cohort$dominant_regime_support == 1)

# ─────────────────────────────────────────────────────────────────────────────
section("Interval filter drops out-of-range windows")
# ─────────────────────────────────────────────────────────────────────────────
memb6 <- mk_memb(list(S = "AABBAA"), window_start_mb = 0.5, window_step_mb = 0.5)
# Windows at 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 Mb. Labels A A B B A A.
# Interval [1.2, 2.7] keeps windows at 1.5, 2.0, 2.5 → labels B,B,A → 3 windows.
res6 <- derive_R_from_regime(memb6, list(start_mb = 1.2, end_mb = 2.7))
ps6  <- res6$per_sample
CHECK("interval filter: 3 windows retained",  ps6$n_windows == 3L)
CHECK("interval filter: dominant = B (2 windows beats 1)",
      ps6$dominant_regime == "B")
CHECK("interval filter: deviation_fraction = 1/3",
      abs(ps6$deviation_fraction - 1/3) < 1e-9)

# ─────────────────────────────────────────────────────────────────────────────
section("Empty/degenerate inputs return the empty result shape")
# ─────────────────────────────────────────────────────────────────────────────
res7 <- derive_R_from_regime(
  data.table(sample = character(), pos_mid_mb = numeric(), group_id = character()),
  list(start_mb = 0, end_mb = 10)
)
CHECK("empty input: per_sample is 0-row data.table",
      is.data.table(res7$per_sample) && nrow(res7$per_sample) == 0L)
CHECK("empty input: cohort n_samples_total = 0",
      res7$cohort$n_samples_total == 0L)

res8 <- derive_R_from_regime(NULL, list(start_mb = 0, end_mb = 10))
CHECK("NULL input: returns empty shape",
      is.data.table(res8$per_sample) && nrow(res8$per_sample) == 0L)

# ─────────────────────────────────────────────────────────────────────────────
section("Back-compat: as_regime_dt() unwrap")
# ─────────────────────────────────────────────────────────────────────────────
dt5 <- as_regime_dt(res5)
CHECK("as_regime_dt(): returns data.table with R_fired column",
      is.data.table(dt5) && "R_fired" %in% names(dt5) && nrow(dt5) == 3L)
CHECK("as_regime_dt(): back-compat n_regime_changes column present",
      "n_regime_changes" %in% names(dt5))

# ─────────────────────────────────────────────────────────────────────────────
section("Back-compat: combine_cohort_recomb() accepts list OR data.table")
# ─────────────────────────────────────────────────────────────────────────────
g_dt <- data.table(
  sample_id = c("SA","SB","SC"),
  ghsl_call = c("UNINFORMATIVE","SPLIT","UNINFORMATIVE")
)
# Pass the DAG list directly (old call sites that haven't migrated):
out_list <- combine_cohort_recomb(r_dt = res5, g_dt = g_dt,
                                  g_resolution = "sufficient")
CHECK("combine_cohort_recomb(list): 3 rows produced",
      is.data.table(out_list) && nrow(out_list) == 3L)
CHECK("combine_cohort_recomb(list): SB = RECOMBINANT HIGH (R & G)",
      out_list[sample_id == "SB"]$recomb_status == "RECOMBINANT" &&
      out_list[sample_id == "SB"]$confidence    == "HIGH")

# Pass the unwrapped data.table (the explicit migration path):
out_dt <- combine_cohort_recomb(r_dt = res5$per_sample, g_dt = g_dt,
                                g_resolution = "sufficient")
CHECK("combine_cohort_recomb(data.table): same result",
      identical(out_list[order(sample_id)]$recomb_status,
                out_dt  [order(sample_id)]$recomb_status))

# ─────────────────────────────────────────────────────────────────────────────
# Finalise
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== RESULTS: ", n_pass, " passed, ", n_fail, " failed ===\n", sep = "")
if (n_fail > 0L) {
  quit(status = 1L)
} else {
  cat("ALL TESTS PASSED\n")
}
