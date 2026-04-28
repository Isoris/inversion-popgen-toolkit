#!/usr/bin/env Rscript
# =============================================================================
# smoke_test_M06.R
#
# Builds synthetic inputs for STEP_M06_emit_boundary_evidence.R, runs the
# script, and validates the output JSON against schema 2.8 §12.
#
# Synthetic scenario: one candidate at C_gar_LG28:12,000,000..14,000,000.
# Three regimes (g0=homo1, g1=het, g2=homo2). Markers in scan range
# (10.5..15.5 Mb at 5 kb resolution → ~1000 markers). Engineered to
# produce visible FST + theta-pi signal at both edges (markers near
# 11.95 Mb and 14.05 Mb show strong allele-freq differences between
# regimes; markers in the candidate body show heterozygosity-elevated
# pattern in g1).
#
# Exits non-zero on any failure. Suitable for CI.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# Test helpers
.pass <- 0L; .fail <- 0L
ok <- function(cond, msg) {
  if (isTRUE(cond)) {
    cat("ok:   ", msg, "\n", sep = "")
    .pass <<- .pass + 1L
  } else {
    cat("FAIL: ", msg, "\n", sep = "")
    .fail <<- .fail + 1L
  }
}

# =============================================================================
# Working dir — stable location so cross_check_M06.js can find it
# =============================================================================
work_root <- "/tmp/m06_smoke_latest"
unlink(work_root, recursive = TRUE, force = TRUE)
dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
cat("[SMOKE] work dir: ", work_root, "\n", sep = "")

CHROM   <- "C_gar_LG28"
N_FISH  <- 60L  # 20 × g0, 20 × g1, 20 × g2
CAND_ID <- 1L
START_BP <- 12000000L
END_BP   <- 14000000L

# =============================================================================
# Synthetic candidates_registry.tsv
# =============================================================================
cands_path <- file.path(work_root, "candidates_registry.tsv")
cands <- data.table(
  candidate_id = CAND_ID,
  chrom = CHROM,
  start_bp = START_BP,
  end_bp = END_BP,
  start_w = 600L,
  end_w   = 700L,
  source = "manual",
  K = 3L
)
fwrite(cands, cands_path, sep = "\t")
cat("[SMOKE] wrote candidates_registry.tsv\n")

# =============================================================================
# Synthetic samples list
# =============================================================================
samples_path <- file.path(work_root, "samples.txt")
sample_ids <- sprintf("F%04d", seq_len(N_FISH))
writeLines(sample_ids, samples_path)
cat("[SMOKE] wrote ", N_FISH, " samples\n", sep = "")

# Regime assignments — first 20 → g0, next 20 → g1, last 20 → g2
regime <- rep(c("g0", "g1", "g2"), each = 20L)

# =============================================================================
# Synthetic fish_regime_calls.tsv
# =============================================================================
fr_path <- file.path(work_root, "fish_regime_calls.tsv")
fr <- data.table(
  candidate_id = CAND_ID,
  sample = sample_ids,
  regime = regime
)
fwrite(fr, fr_path, sep = "\t")
cat("[SMOKE] wrote fish_regime_calls.tsv\n")

# =============================================================================
# Synthetic dosage matrix
# =============================================================================
# Markers at 5 kb spacing across 10 Mb..16 Mb → 1200 markers
# - "edge" markers (near 11.95 Mb and 14.05 Mb): strong FST signal
#   (g0 dominantly homo-ref, g2 dominantly homo-alt)
# - "body" markers (12..14 Mb): elevated het signal in g1 (hetero advantage)
# - "outside" markers: neutral, low signal
set.seed(42)
n_markers <- 1200L
marker_pos <- as.integer(seq(10000000L, 16000000L, length.out = n_markers))

# Build dosage matrix
dosage <- matrix(NA_integer_, nrow = n_markers, ncol = N_FISH)
for (mi in seq_len(n_markers)) {
  pos <- marker_pos[mi]
  # Determine "zone" — outside / left-edge / body / right-edge / outside
  is_left_edge  <- pos >= 11900000L && pos <= 12000000L
  is_right_edge <- pos >= 14000000L && pos <= 14100000L
  is_body       <- pos > 12000000L && pos < 14000000L

  # Per-regime allele freq for ALT
  if (is_left_edge || is_right_edge) {
    # Strong FST signal — g0 mostly REF (low p), g2 mostly ALT (high p)
    p_g0 <- 0.10
    p_g1 <- 0.50
    p_g2 <- 0.90
  } else if (is_body) {
    # Body: g1 elevated het (not so different freq; just hetero biased)
    p_g0 <- 0.10
    p_g1 <- 0.50
    p_g2 <- 0.90
  } else {
    # Neutral region — all regimes share similar freq
    p_g0 <- 0.30
    p_g1 <- 0.30
    p_g2 <- 0.30
  }

  for (si in seq_len(N_FISH)) {
    p <- if (regime[si] == "g0") p_g0
         else if (regime[si] == "g1") p_g1
         else p_g2
    if (regime[si] == "g1" && is_body) {
      # Body markers: g1 forced heterozygous (matches IH-W expectation)
      dosage[mi, si] <- 1L
    } else {
      # Sample dosage as binomial(2, p) — 0/1/2
      dosage[mi, si] <- rbinom(1L, 2L, p)
    }
  }
}

# Add a few NA cells (~1% of cells) to test NA handling
n_na <- as.integer(0.01 * n_markers * N_FISH)
na_idx <- sample(seq_len(n_markers * N_FISH), n_na, replace = FALSE)
dosage[na_idx] <- -1L   # -1 = NA encoding, M06 should map to NA

# Build dosage TSV: marker_id + pos_bp + N_FISH columns
dos_path <- file.path(work_root, paste0(CHROM, ".dos.tsv"))
dos_dt <- data.table(
  marker_id = paste0(CHROM, "_", marker_pos),
  pos_bp = marker_pos
)
for (si in seq_len(N_FISH)) {
  dos_dt[[sample_ids[si]]] <- dosage[, si]
}
fwrite(dos_dt, dos_path, sep = "\t")
cat("[SMOKE] wrote dosage matrix: ", n_markers, " markers x ", N_FISH, " samples\n", sep = "")

# =============================================================================
# Synthetic SV VCFs — 1 DELLY (INV at left edge), 1 Manta (INV5 at right edge)
# =============================================================================
delly_vcf <- file.path(work_root, "delly_inv.vcf")
delly_lines <- c(
  "##fileformat=VCFv4.2",
  "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
  "##INFO=<ID=CT,Number=1,Type=String,Description=\"Connection type\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  paste0(CHROM, "\t11950000\tINV0001\tN\t<INV>\t80\tPASS\tSVTYPE=INV;CT=3to3;DELLY=v1"),
  paste0(CHROM, "\t14025000\tINV0002\tN\t<INV>\t60\tPASS\tSVTYPE=INV;CT=5to5;DELLY=v1")
)
writeLines(delly_lines, delly_vcf)

manta_vcf <- file.path(work_root, "manta_inv.vcf")
manta_lines <- c(
  "##fileformat=VCFv4.2",
  "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  paste0(CHROM, "\t14050000\tMantaINV1\tN\t<INV>\t95\tPASS\tSVTYPE=INV;INV5"),
  # OFF-TARGET — out of scan region; should NOT appear in output
  paste0(CHROM, "\t20000000\tMantaINV2\tN\t<INV>\t75\tPASS\tSVTYPE=INV;INV3"),
  # OFF-CHROM — should NOT appear
  paste0("C_gar_LG27\t12000000\tMantaINV3\tN\t<INV>\t60\tPASS\tSVTYPE=INV;INV3")
)
writeLines(manta_lines, manta_vcf)
cat("[SMOKE] wrote DELLY VCF (2 INVs) and Manta VCF (1 in-range, 2 out-of-range)\n")

# =============================================================================
# Synthetic discordant-pair BED
# =============================================================================
# Spike near both edges, low background everywhere else
disc_bed <- file.path(work_root, "discordant.bed")
bed_rows <- list()
# Background — small piles every 50 kb across entire scan range
for (pos in seq(10500000L, 15500000L, by = 50000L)) {
  bed_rows[[length(bed_rows) + 1L]] <- data.frame(
    chrom = CHROM, start = pos, end = pos + 1000L, count = 2L
  )
}
# Edge spikes — left edge ~11.95 Mb, right edge ~14.05 Mb, big counts
for (pos in seq(11900000L, 12000000L, by = 5000L)) {
  bed_rows[[length(bed_rows) + 1L]] <- data.frame(
    chrom = CHROM, start = pos, end = pos + 1000L, count = 50L
  )
}
for (pos in seq(14000000L, 14100000L, by = 5000L)) {
  bed_rows[[length(bed_rows) + 1L]] <- data.frame(
    chrom = CHROM, start = pos, end = pos + 1000L, count = 80L
  )
}
# Off-chrom row (should be filtered out)
bed_rows[[length(bed_rows) + 1L]] <- data.frame(
  chrom = "C_gar_LG27", start = 5000000L, end = 5001000L, count = 999L
)
bed_dt <- rbindlist(bed_rows)
fwrite(bed_dt, disc_bed, sep = "\t", col.names = FALSE)
cat("[SMOKE] wrote discordant BED: ", nrow(bed_dt), " rows\n", sep = "")

# =============================================================================
# RUN STEP_M06
# =============================================================================
out_dir <- file.path(work_root, "scrubber_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

m06_path <- "/home/claude/work/r_module6/STEP_M06_emit_boundary_evidence.R"
cat("\n[SMOKE] running STEP_M06 ...\n")
cmd <- paste(
  "Rscript", shQuote(m06_path),
  "--candidates_registry", shQuote(cands_path),
  "--chrom", shQuote(CHROM),
  "--out_dir", shQuote(out_dir),
  "--scan_radius_bp 1500000",
  "--scan_window_bp 5000",
  "--dosage", shQuote(dos_path),
  "--fish_regimes", shQuote(fr_path),
  "--samples", shQuote(samples_path),
  "--sv_vcf", shQuote(delly_vcf), shQuote(manta_vcf),
  "--discordant_bed", shQuote(disc_bed)
)
status <- system(cmd)
ok(status == 0L, "STEP_M06 exited with status 0")
if (status != 0L) {
  cat("\n[SMOKE] M06 failed — aborting\n")
  quit(status = 1)
}

# =============================================================================
# Validate output JSON
# =============================================================================
out_json <- file.path(out_dir, CHROM, paste0(CHROM, "_boundary_evidence.json"))
ok(file.exists(out_json), paste0("output JSON exists: ", out_json))

j <- jsonlite::fromJSON(out_json, simplifyVector = FALSE)

# Schema-level checks
ok(j$schema_version == 2L, "schema_version = 2")
ok("boundary_evidence" %in% j$`_layers_present`,
   "_layers_present includes boundary_evidence")
ok(grepl("STEP_M06", j$`_generator`),
   "_generator is STEP_M06_emit_boundary_evidence.R")

# boundary_evidence array
be <- j$boundary_evidence
ok(length(be) == 1L, paste0("one candidate in boundary_evidence (got ", length(be), ")"))

cand <- be[[1]]
ok(cand$candidate_id == CAND_ID, paste0("candidate_id matches: ", cand$candidate_id))
ok(cand$chrom == CHROM, "chrom matches")
ok(cand$scan_window_bp == 5000L, "scan_window_bp = 5000")
ok(cand$scan_start_bp <= START_BP, paste0("scan_start_bp <= start_bp: ", cand$scan_start_bp, " <= ", START_BP))
ok(cand$scan_end_bp >= END_BP, "scan_end_bp >= end_bp")

# Tracks present
tracks <- cand$tracks
ok(!is.null(tracks$fst), "tracks.fst present")
ok(!is.null(tracks$theta_pi_homo1), "tracks.theta_pi_homo1 present")
ok(!is.null(tracks$theta_pi_het), "tracks.theta_pi_het present")
ok(!is.null(tracks$theta_pi_homo2), "tracks.theta_pi_homo2 present")
ok(!is.null(tracks$discordant_pair_pileup), "tracks.discordant_pair_pileup present")
ok(!is.null(tracks$sv_anchors), "tracks.sv_anchors present")

# Track lengths consistent
expected_n_win <- as.integer((cand$scan_end_bp - cand$scan_start_bp) / 5000L)
ok(length(tracks$fst) == expected_n_win,
   paste0("fst length matches scan window count (", length(tracks$fst), " == ", expected_n_win, ")"))
ok(length(tracks$theta_pi_homo1) == expected_n_win, "theta_pi_homo1 length matches")
ok(length(tracks$theta_pi_het) == expected_n_win,   "theta_pi_het length matches")
ok(length(tracks$theta_pi_homo2) == expected_n_win, "theta_pi_homo2 length matches")
ok(length(tracks$discordant_pair_pileup) == expected_n_win,
   paste0("discordant_pair_pileup length matches (", length(tracks$discordant_pair_pileup), ")"))

# FST signal — should be high near edges, low elsewhere
# Convert "null" entries (which jsonlite gives as NULL list elements) to NA
unlist_safe <- function(lst) {
  vapply(lst, function(x) if (is.null(x)) NA_real_ else as.numeric(x), numeric(1))
}
fst_v <- unlist_safe(tracks$fst)
left_edge_idx  <- which(seq_along(fst_v) >= floor((11900000 - cand$scan_start_bp) / 5000) + 1L
                        & seq_along(fst_v) <= floor((12000000 - cand$scan_start_bp) / 5000) + 1L)
right_edge_idx <- which(seq_along(fst_v) >= floor((14000000 - cand$scan_start_bp) / 5000) + 1L
                        & seq_along(fst_v) <= floor((14100000 - cand$scan_start_bp) / 5000) + 1L)

cat("[SMOKE] FST stats:\n")
cat("[SMOKE]   left edge (11.9..12.0 Mb): mean = ", round(mean(fst_v[left_edge_idx], na.rm = TRUE), 4), "\n", sep = "")
cat("[SMOKE]   right edge (14.0..14.1 Mb): mean = ", round(mean(fst_v[right_edge_idx], na.rm = TRUE), 4), "\n", sep = "")
cat("[SMOKE]   global mean: ", round(mean(fst_v, na.rm = TRUE), 4), "\n", sep = "")
ok(mean(fst_v[left_edge_idx], na.rm = TRUE) > 0.5,
   paste0("FST elevated at left edge (>0.5; got ",
          round(mean(fst_v[left_edge_idx], na.rm = TRUE), 3), ")"))
ok(mean(fst_v[right_edge_idx], na.rm = TRUE) > 0.5,
   "FST elevated at right edge (>0.5)")

# Discordant pile: should spike at both edges
disc_v <- unlist_safe(tracks$discordant_pair_pileup)
left_disc <- max(disc_v[left_edge_idx], na.rm = TRUE)
right_disc <- max(disc_v[right_edge_idx], na.rm = TRUE)
cat("[SMOKE]   discordant left edge max: ", left_disc, "\n", sep = "")
cat("[SMOKE]   discordant right edge max: ", right_disc, "\n", sep = "")
ok(left_disc >= 50L, paste0("discordant pile spike at left edge (>=50; got ", left_disc, ")"))
ok(right_disc >= 80L, "discordant pile spike at right edge (>=80)")

# SV anchors: 3 expected (DELLY left, DELLY right, Manta_INV5 right). The
# off-target Manta calls (out-of-scan and off-chrom) must be filtered.
sv <- tracks$sv_anchors
ok(length(sv) == 3L, paste0("3 SV anchors in scan range (got ", length(sv), ")"))
sv_kinds <- vapply(sv, function(x) x$kind, character(1))
sv_pos   <- vapply(sv, function(x) x$pos_bp, integer(1))
ok(any(grepl("DELLY", sv_kinds)), "DELLY anchor present")
ok(any(grepl("Manta_INV5", sv_kinds)), "Manta_INV5 anchor present")
ok(all(sv_pos >= cand$scan_start_bp & sv_pos <= cand$scan_end_bp),
   "all SV anchors are inside scan range")
# Off-target Manta_INV3 at 20 Mb should NOT be present
ok(!any(sv_kinds == "Manta_INV3"),
   "off-target Manta_INV3 (out of scan range) was filtered out")

# CT field — DELLY anchors should have CT = 3to3 or 5to5
delly_idx <- which(grepl("DELLY", sv_kinds))
ok(all(vapply(sv[delly_idx], function(x) !is.null(x$ct) && x$ct %in% c("3to3", "5to5"), logical(1))),
   "DELLY anchors carry CT field with valid values")

# =============================================================================
# Summary
# =============================================================================
cat("\n[SMOKE] === SUMMARY ===\n")
cat("[SMOKE] Pass: ", .pass, "\n", sep = "")
cat("[SMOKE] Fail: ", .fail, "\n", sep = "")
cat("[SMOKE] Output JSON: ", out_json, "\n", sep = "")
cat("[SMOKE]   size: ", round(file.info(out_json)$size / 1024, 1), " KB\n", sep = "")

if (.fail > 0L) quit(status = 1)
cat("[SMOKE] OK — STEP_M06 produces a valid boundary_evidence JSON.\n")
