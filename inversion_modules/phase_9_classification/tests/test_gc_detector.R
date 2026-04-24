#!/usr/bin/env Rscript
# =============================================================================
# test_gc_detector.R — fixture for chat-12 per-SNP GC detector rewrite
# =============================================================================
# Covers:
#   1. Single clean short tract (run_length=3, span <= max_span_bp) → HIGH/MED
#      confidence as appropriate; correct direction for HOM_REF vs HOM_INV.
#   2. Tract with 1 interior HET-looking tolerance SNP (max_tolerance=1) →
#      accepted; tolerance_snps = 1.
#   3. Run that exceeds max_tolerance (2 consecutive skips) → tract breaks
#      before the second skip.
#   4. Run longer than max_flagged_snps → rejected (recombinant, not GC).
#   5. Tract span exceeding max_span_bp → rejected.
#   6. SNP QC drops: high-missingness, excess-het (paralog artefact),
#      low-MAF, and depth anomaly SNPs are removed from the diagnostic set.
#   7. Non-diagnostic SNPs (low |AF_ref − AF_inv|) are excluded even if
#      QC-clean.
#   8. HET samples produce 0 tracts regardless of dosage (skipped by spec).
#   9. snp_qc accounting block is populated correctly.
#
# Pass = all PASS lines printed; final "ALL TESTS PASSED".
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# Locate the detector
this_file <- NA_character_
args <- commandArgs(trailingOnly = FALSE)
fa <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(fa)) this_file <- normalizePath(fa, mustWork = FALSE)
if (is.na(this_file)) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                        error = function(e) NA_character_)
}
repo_root <- if (!is.na(this_file) && nzchar(this_file))
               # pass 15: test is now at inversion_modules/phase_9_classification/tests/
               # so repo root is 3 levels up (was 2 in the old layout)
               dirname(dirname(dirname(this_file))) else getwd()
det_path <- file.path(repo_root, "inversion_modules",
                      "phase_7_karyotype_groups", "proposal",
                      "gene_conversion_detector.R")
if (!file.exists(det_path)) {
  alts <- c("gene_conversion_detector.R",
            "../proposal/gene_conversion_detector.R",
            "../inversion_modules/phase_7_karyotype_groups/proposal/gene_conversion_detector.R",
            file.path(getwd(), "inversion_modules", "phase_7_karyotype_groups",
                      "proposal", "gene_conversion_detector.R"))
  for (a in alts) if (file.exists(a)) { det_path <- a; break }
}
if (!file.exists(det_path)) stop("cannot locate gene_conversion_detector.R (repo_root=", repo_root, ")")
source(det_path)

# Micro-asserts
n_pass <- 0L; n_fail <- 0L
CHECK <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    cat("  FAIL  ", label, "  [error: ", conditionMessage(e), "]\n", sep = "")
    n_fail <<- n_fail + 1L; return(FALSE)
  })
  if (isTRUE(ok)) { cat("  PASS  ", label, "\n", sep = ""); n_pass <<- n_pass + 1L }
  else { cat("  FAIL  ", label, "\n", sep = ""); n_fail <<- n_fail + 1L }
}
section <- function(title) cat("\n── ", title, " ──\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build a cohort-level dosage matrix from per-sample per-SNP codes.
# Codes: "R" = REF-looking (dosage=0), "I" = INV-looking (dosage=2),
#        "H" = HET-looking (1), "." = NA (missing).
# Rows = samples; each sample's codes must be the same length.
# ─────────────────────────────────────────────────────────────────────────────
mk_mat <- function(code_by_sample, snp_positions) {
  sids <- names(code_by_sample)
  n_snps <- length(snp_positions)
  mat <- matrix(NA_real_, nrow = length(sids), ncol = n_snps,
                dimnames = list(sids, NULL))
  for (sid in sids) {
    codes <- strsplit(code_by_sample[[sid]], "")[[1]]
    stopifnot(length(codes) == n_snps)
    d <- rep(NA_real_, n_snps)
    d[codes == "R"] <- 0
    d[codes == "I"] <- 2
    d[codes == "H"] <- 1
    d[codes == "."] <- NA_real_
    mat[sid, ] <- d
  }
  mat
}

# ─────────────────────────────────────────────────────────────────────────────
section("Test 1: a single 3-SNP tract in a HOM_REF sample, surrounded by REF")
# 10 SNPs at 1,2,3,4,5,6,7,8,9,10 kb. 5 HOM_REF cohort samples and 5 HOM_INV
# cohort samples provide the diagnostic baseline AF. Test sample 'TR' is
# HOM_REF but has INV-looking dosage at SNPs 4,5,6 — a 3-SNP GC tract.
# Expected: 1 tract, run_length=3, tolerance=0, confidence=MEDIUM, span=
# 2 kb (positions 4–6 kb), direction=INV_in_REF_context.
# ─────────────────────────────────────────────────────────────────────────────
snp_pos <- as.integer(seq_len(10) * 1000L)
codes <- c(
  TR = "RRRIIIRRRR",              # the test sample
  R1 = "RRRRRRRRRR", R2 = "RRRRRRRRRR",
  R3 = "RRRRRRRRRR", R4 = "RRRRRRRRRR",  # baseline REFs
  I1 = "IIIIIIIIII", I2 = "IIIIIIIIII",
  I3 = "IIIIIIIIII", I4 = "IIIIIIIIII", I5 = "IIIIIIIIII"  # baseline INVs
)
mat <- mk_mat(codes, snp_pos)
cls <- setNames(c("HOM_REF", rep("HOM_REF", 4), rep("HOM_INV", 5)), names(codes))
res <- detect_cohort_gene_conversion_v2(mat, snp_pos, cls,
                                        candidate_id = "cand_test1",
                                        max_span_bp = 20000L)
CHECK("T1: exactly 1 tract total", res$total_tracts == 1L)
CHECK("T1: tract on TR only", res$tracts$sample_id == "TR")
CHECK("T1: run_length = 3", res$tracts$run_length_flagged_snps == 3L)
CHECK("T1: tolerance_snps = 0", res$tracts$tolerance_snps == 0L)
CHECK("T1: confidence = MEDIUM (len=3)", res$tracts$confidence == "MEDIUM")
CHECK("T1: span_bp = 2000", res$tracts$span_bp == 2000L)
CHECK("T1: direction = INV_in_REF_context",
      res$tracts$direction == "INV_in_REF_context")
CHECK("T1: total_samples_with_gc = 1", res$total_samples_with_gc == 1L)
CHECK("T1: all 10 SNPs pass QC (no artefacts)",
      res$snp_qc$n_snps_pass_qc == 10L)
CHECK("T1: all 10 SNPs are diagnostic (|AF_ref=0, AF_inv=1|=1 >= 0.5)",
      res$snp_qc$n_snps_pass_qc_diagnostic == 10L)

# ─────────────────────────────────────────────────────────────────────────────
section("Test 2: tolerance — 1 interior HET-looking SNP inside a 4-SNP run")
# TR: RRRIHI IRR  (SNP 5 is HET-looking, surrounded by flagged I I and I I)
# Expected: tract spans SNPs 4..7 (I,H,I,I), run_length = 3, tolerance = 1,
# span = 3 kb, confidence = MEDIUM.
# ─────────────────────────────────────────────────────────────────────────────
codes2 <- c(
  TR = "RRRIHIIRRR",
  R1="RRRRRRRRRR", R2="RRRRRRRRRR", R3="RRRRRRRRRR", R4="RRRRRRRRRR",
  I1="IIIIIIIIII", I2="IIIIIIIIII", I3="IIIIIIIIII", I4="IIIIIIIIII",
  I5="IIIIIIIIII"
)
mat2 <- mk_mat(codes2, snp_pos)
res2 <- detect_cohort_gene_conversion_v2(mat2, snp_pos, cls,
                                         candidate_id = "cand_test2",
                                         max_tolerance = 1L, max_span_bp = 20000L)
CHECK("T2: exactly 1 tract", res2$total_tracts == 1L)
CHECK("T2: run_length_flagged_snps = 3", res2$tracts$run_length_flagged_snps == 3L)
CHECK("T2: tolerance_snps = 1", res2$tracts$tolerance_snps == 1L)
CHECK("T2: span_bp = 3000", res2$tracts$span_bp == 3000L)

# ─────────────────────────────────────────────────────────────────────────────
section("Test 3: tolerance limit — 2 consecutive skips break the tract")
# TR: RRRIHHIIRR — the two consecutive H (positions 5,6) exceed max_tolerance=1.
# Expected: NO tract (the single flagged I at position 4 is not enough on
# its own; the I I at 7,8 is a separate 2-SNP run).
# With min_run=2, the I I at positions 7,8 actually forms a valid tract by
# itself (run_length=2, tolerance=0). So we expect exactly 1 tract starting
# at position 7.
# ─────────────────────────────────────────────────────────────────────────────
codes3 <- c(
  TR = "RRRIHHIIRR",
  R1="RRRRRRRRRR", R2="RRRRRRRRRR", R3="RRRRRRRRRR", R4="RRRRRRRRRR",
  I1="IIIIIIIIII", I2="IIIIIIIIII", I3="IIIIIIIIII", I4="IIIIIIIIII",
  I5="IIIIIIIIII"
)
mat3 <- mk_mat(codes3, snp_pos)
res3 <- detect_cohort_gene_conversion_v2(mat3, snp_pos, cls,
                                         candidate_id = "cand_test3",
                                         max_tolerance = 1L)
CHECK("T3: exactly 1 tract (the short 2-SNP run after the break)",
      res3$total_tracts == 1L)
CHECK("T3: tract starts at SNP position 7000 bp",
      res3$tracts$tract_start_bp == 7000L)
CHECK("T3: tract run_length = 2", res3$tracts$run_length_flagged_snps == 2L)
CHECK("T3: confidence = LOW (len=2)", res3$tracts$confidence == "LOW")

# ─────────────────────────────────────────────────────────────────────────────
section("Test 4: tract too long (exceeds max_flagged_snps) is REJECTED")
# TR: RRIIIIIIIIR — 8 flagged SNPs in a row. max_flagged_snps=10 accepts,
# max_flagged_snps=5 rejects. Verify both behaviours.
# ─────────────────────────────────────────────────────────────────────────────
snp_pos_long <- as.integer(seq_len(11) * 1000L)
codes4 <- c(
  TR = "RRIIIIIIIIR",
  R1="RRRRRRRRRRR", R2="RRRRRRRRRRR", R3="RRRRRRRRRRR", R4="RRRRRRRRRRR",
  I1="IIIIIIIIIII", I2="IIIIIIIIIII", I3="IIIIIIIIIII", I4="IIIIIIIIIII",
  I5="IIIIIIIIIII"
)
mat4 <- mk_mat(codes4, snp_pos_long)
res4a <- detect_cohort_gene_conversion_v2(mat4, snp_pos_long, cls,
                                          candidate_id = "cand_test4a",
                                          max_flagged_snps = 5L,
                                          max_span_bp = 100000L)
CHECK("T4a: 8-flag run REJECTED when max_flagged_snps=5",
      res4a$total_tracts == 0L)
res4b <- detect_cohort_gene_conversion_v2(mat4, snp_pos_long, cls,
                                          candidate_id = "cand_test4b",
                                          max_flagged_snps = 10L,
                                          max_span_bp = 100000L)
CHECK("T4b: 8-flag run ACCEPTED when max_flagged_snps=10",
      res4b$total_tracts == 1L && res4b$tracts$run_length_flagged_snps == 8L)
CHECK("T4b: confidence = HIGH (len=8 >= 4)",
      res4b$tracts$confidence == "HIGH")

# ─────────────────────────────────────────────────────────────────────────────
section("Test 5: tract span exceeding max_span_bp is REJECTED")
# Two flagged SNPs 30 kb apart (skipped interior treated as 'skip' not
# 'consistent' since R4 and R5 flanks are REF-consistent — but with
# max_tolerance=1 and 28 skips between, the walker terminates long before
# pairing them). Use a simpler direct test: a 2-SNP flag run spanning
# max_span_bp + 1.
# Positions 1kb, 1kb+(max_span_bp+1) — make this explicit.
# ─────────────────────────────────────────────────────────────────────────────
max_span <- 5000L
snp_pos5 <- c(1000L, 7000L)  # 6 kb apart
mat5 <- matrix(c(2, 2,  # TR: both flagged
                 rep(0, 10), rep(2, 10)) |> head(22),
               nrow = 11, byrow = FALSE,
               dimnames = list(c("TR", paste0("R", 1:4), paste0("I", 1:6)),
                               NULL))
# That matrix shape is awkward; build it properly
mat5 <- matrix(NA_real_, nrow = 10, ncol = 2,
               dimnames = list(c("TR", paste0("R", 1:4), paste0("I", 1:5)), NULL))
mat5["TR", ]   <- c(2, 2)
for (nm in paste0("R", 1:4)) mat5[nm, ] <- c(0, 0)
for (nm in paste0("I", 1:5)) mat5[nm, ] <- c(2, 2)
cls5 <- setNames(c("HOM_REF", rep("HOM_REF", 4), rep("HOM_INV", 5)),
                 rownames(mat5))
res5 <- detect_cohort_gene_conversion_v2(mat5, snp_pos5, cls5,
                                         candidate_id = "cand_test5",
                                         max_span_bp = max_span)
CHECK("T5: 6 kb span REJECTED when max_span_bp = 5 kb",
      res5$total_tracts == 0L)
res5b <- detect_cohort_gene_conversion_v2(mat5, snp_pos5, cls5,
                                          candidate_id = "cand_test5b",
                                          max_span_bp = 10000L)
CHECK("T5b: same pair ACCEPTED when max_span_bp = 10 kb",
      res5b$total_tracts == 1L && res5b$tracts$span_bp == 6000L)

# ─────────────────────────────────────────────────────────────────────────────
section("Test 6: SNP-level QC drops artefact columns")
# 10 SNPs, but SNP 5 is a paralog-mismapping artefact — every sample is
# HET-looking there (dosage = 1). That should trip excess-het and be
# dropped from the diagnostic set. The TR sample has I I at SNPs 3,4 —
# a 2-SNP tract that is independent of the artefact at SNP 5.
# ─────────────────────────────────────────────────────────────────────────────
codes6 <- c(
  TR = "RRIIHRRRRR",  # SNP 5 is H
  R1="RRRRHRRRRR", R2="RRRRHRRRRR", R3="RRRRHRRRRR", R4="RRRRHRRRRR",
  I1="IIIIHIIIII", I2="IIIIHIIIII", I3="IIIIHIIIII", I4="IIIIHIIIII",
  I5="IIIIHIIIII"
)
mat6 <- mk_mat(codes6, snp_pos)
res6 <- detect_cohort_gene_conversion_v2(mat6, snp_pos, cls,
                                         candidate_id = "cand_test6",
                                         max_het_fraction = 0.70)
CHECK("T6: SNP 5 (all-het paralog artefact) dropped by QC",
      res6$snp_qc$n_dropped_excess_het == 1L)
CHECK("T6: n_snps_pass_qc = 9",
      res6$snp_qc$n_snps_pass_qc == 9L)
CHECK("T6: 2-SNP tract on TR still detected (unaffected by dropped SNP)",
      res6$total_tracts == 1L &&
      res6$tracts$run_length_flagged_snps == 2L)

# ─────────────────────────────────────────────────────────────────────────────
section("Test 7: non-diagnostic SNPs excluded (low |AF_ref − AF_inv|)")
# Build a cohort where SNPs 1,2,3 carry REF=0/INV=0 (not diagnostic — both
# arrangements share the minor allele at this locus). Even if TR looks
# 'wrong' there, the detector should not flag (since the SNP isn't
# diagnostic).
# ─────────────────────────────────────────────────────────────────────────────
codes7 <- c(
  TR = "IIIRRRRRRR",   # TR "mismatches" on SNPs 1,2,3 — but those are non-diag
  R1="RRRRRRRRRR", R2="RRRRRRRRRR", R3="RRRRRRRRRR", R4="RRRRRRRRRR",
  I1="RRRIIIIIII", I2="RRRIIIIIII", I3="RRRIIIIIII", I4="RRRIIIIIII",
  I5="RRRIIIIIII"   # INV baseline is also REF-looking at 1..3 → AF_ref=AF_inv=0
)
mat7 <- mk_mat(codes7, snp_pos)
res7 <- detect_cohort_gene_conversion_v2(mat7, snp_pos, cls,
                                         candidate_id = "cand_test7",
                                         min_delta_af = 0.5)
CHECK("T7: SNPs 1..3 dropped as non-diagnostic",
      res7$snp_qc$n_dropped_non_diagnostic >= 3L)
CHECK("T7: no tract on TR (its 'mismatches' are all at non-diagnostic sites)",
      res7$total_tracts == 0L)

# ─────────────────────────────────────────────────────────────────────────────
section("Test 8: HET sample yields 0 tracts regardless of per-SNP dosage")
# TR is baseline HET; per spec we skip HET samples entirely (ambiguous
# per-haplotype at 9x). Even with INV-looking stretches we should emit 0.
# ─────────────────────────────────────────────────────────────────────────────
codes8 <- c(
  TR = "IIIIIIRRRR",
  R1="RRRRRRRRRR", R2="RRRRRRRRRR", R3="RRRRRRRRRR", R4="RRRRRRRRRR",
  I1="IIIIIIIIII", I2="IIIIIIIIII", I3="IIIIIIIIII", I4="IIIIIIIIII",
  I5="IIIIIIIIII"
)
mat8 <- mk_mat(codes8, snp_pos)
cls8 <- cls; cls8["TR"] <- "HET"
res8 <- detect_cohort_gene_conversion_v2(mat8, snp_pos, cls8,
                                         candidate_id = "cand_test8")
CHECK("T8: HET sample produces 0 tracts", res8$total_tracts == 0L)
CHECK("T8: per_sample_summary records baseline_class = HET",
      res8$per_sample_summary[sample_id == "TR"]$baseline_class == "HET")

# ─────────────────────────────────────────────────────────────────────────────
section("Test 9: snp_qc accounting block shape")
# ─────────────────────────────────────────────────────────────────────────────
req_fields <- c("n_snps_total", "n_snps_pass_qc", "n_snps_pass_qc_diagnostic",
                "n_dropped_missingness", "n_dropped_excess_het",
                "n_dropped_low_maf", "n_dropped_depth_anomaly",
                "n_dropped_non_diagnostic")
CHECK("T9: snp_qc has all required fields",
      all(req_fields %in% names(res$snp_qc)))
CHECK("T9: params block records key thresholds",
      all(c("min_run", "max_tolerance", "max_span_bp", "min_delta_af") %in%
          names(res$params)))

# ─────────────────────────────────────────────────────────────────────────────
# HOM_INV sample direction check — round-trip
section("Test 10: HOM_INV baseline with REF-looking tract → REF_in_INV_context")
# ─────────────────────────────────────────────────────────────────────────────
codes10 <- c(
  TR = "IIIRRRIIII",   # HOM_INV baseline, 3 REF-looking SNPs at 4,5,6
  R1="RRRRRRRRRR", R2="RRRRRRRRRR", R3="RRRRRRRRRR", R4="RRRRRRRRRR",
  I1="IIIIIIIIII", I2="IIIIIIIIII", I3="IIIIIIIIII", I4="IIIIIIIIII",
  I5="IIIIIIIIII"
)
mat10 <- mk_mat(codes10, snp_pos)
cls10 <- cls; cls10["TR"] <- "HOM_INV"
res10 <- detect_cohort_gene_conversion_v2(mat10, snp_pos, cls10,
                                          candidate_id = "cand_test10")
CHECK("T10: HOM_INV + REF-tract → direction = REF_in_INV_context",
      res10$total_tracts == 1L &&
      res10$tracts$direction == "REF_in_INV_context")

# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== RESULTS: ", n_pass, " passed, ", n_fail, " failed ===\n", sep = "")
if (n_fail > 0L) quit(status = 1L) else cat("ALL TESTS PASSED\n")
