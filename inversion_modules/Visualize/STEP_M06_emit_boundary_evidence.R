#!/usr/bin/env Rscript

# =============================================================================
# STEP_M06_emit_boundary_evidence.R
#
# Per-chromosome emit of the optional `boundary_evidence` layer for the
# scrubber's v3.97 boundary zone refinement module (schema 2.8 §12).
#
# This script feeds the four optional tracks that the scrubber's
# `_buildBoundaryTrackScores` already detects and consumes when present:
#
#   - fst                        : per-window Hudson FST between regimes
#                                  (homo1 vs homo2 of the K=3 candidate)
#   - theta_pi_homo1/het/homo2   : per-window regime-conditional pairwise
#                                  nucleotide diversity (sample-size corrected)
#   - discordant_pair_pileup     : per-window count of discordant read pairs
#                                  (from a pre-aggregated BED — typically
#                                   produced by samtools view -F 0x2 + bedtools)
#   - sv_anchors                 : DELLY/Manta INV/BND/DEL calls inside the
#                                  candidate scan region
#
# Each track is independently optional — when its source input is missing,
# the track is simply absent from `tracks: { ... }`. The scrubber skips
# absent tracks gracefully.
#
# Output: ONE JSON per chromosome:
#   <out_dir>/<chrom>/<chrom>_boundary_evidence.json
#
# JSON structure (schema 2.8 §12):
# {
#   "schema_version": 2,
#   "_layers_present": ["boundary_evidence"],
#   "boundary_evidence": [
#     {
#       "candidate_id": <int>,
#       "chrom": "<chrom>",
#       "scan_start_bp": <int>,           # candidate ± scan_radius_bp, clamped
#       "scan_end_bp":   <int>,
#       "scan_window_bp": <int>,          # bp resolution of arrays below
#       "tracks": {
#         "fst": [ ... ],                 # length = floor((end-start)/win)
#         "theta_pi_homo1": [ ... ],
#         "theta_pi_het":   [ ... ],
#         "theta_pi_homo2": [ ... ],
#         "discordant_pair_pileup": [ ... ],
#         "sv_anchors": [ {kind, ct?, pos_bp, qual?}, ... ]
#       }
#     },
#     ...
#   ]
# }
#
# Usage:
#   Rscript STEP_M06_emit_boundary_evidence.R \
#     --candidates_registry candidates_registry.tsv \
#     --chrom C_gar_LG28 \
#     --out_dir scrubber_data/ \
#     [--scan_radius_bp 1500000] \
#     [--scan_window_bp 5000] \
#     [--dosage  C_gar_LG28.dos.tsv.gz] \
#     [--fish_regimes fish_regime_calls.tsv] \
#     [--sv_vcf delly.inv.vcf manta.inv.vcf] \
#     [--discordant_bed discordant_pairs.bed] \
#     [--samples sample_list.txt]
#
# Conventions:
#   - All integers are emitted as JSON integers.
#   - All floats are rounded to 6dp (FST/θπ are typically 0..0.5 range).
#   - NA in any track maps to JSON null. Scrubber treats null as
#     non-finite and excludes from per-window aggregations.
#   - The script aborts ONLY on missing required inputs or candidates_registry
#     parse failure. Optional-input failures (missing VCF, malformed BED,
#     etc.) emit warnings and skip that track.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# =============================================================================
# CLI parsing
# =============================================================================
usage <- function() {
  cat(paste(
    "Usage: Rscript STEP_M06_emit_boundary_evidence.R \\",
    "  --candidates_registry <file.tsv> --chrom <chrom> --out_dir <dir>",
    "  [--scan_radius_bp 1500000]",
    "  [--scan_window_bp 5000]",
    "  [--dosage <chrom>.dos.tsv.gz]",
    "  [--fish_regimes fish_regime_calls.tsv]",
    "  [--sv_vcf vcf1 vcf2 ...]",
    "  [--discordant_bed discordant_pairs.bed]",
    "  [--samples sample_list.txt]",
    "",
    "Required:",
    "  --candidates_registry   v3.97 candidate registry TSV (must include candidate_id,",
    "                          chrom, start_bp, end_bp).",
    "  --chrom                 Restrict to this chromosome (must match registry's chrom column).",
    "  --out_dir               Output dir; writes <out_dir>/<chrom>/<chrom>_boundary_evidence.json",
    "",
    "Optional inputs (each enables a different track):",
    "  --dosage + --fish_regimes  → enables fst, theta_pi_homo1, theta_pi_het, theta_pi_homo2",
    "  --sv_vcf                   → enables sv_anchors (1+ VCF files)",
    "  --discordant_bed           → enables discordant_pair_pileup (4-col BED: chrom start end count)",
    "  --samples                  → cohort sample order (defaults to dosage TSV header order)",
    "",
    "Optional tuning:",
    "  --scan_radius_bp N      Default 1500000 (1.5 Mb); matches scrubber default",
    "  --scan_window_bp N      Default 5000 (5 kb); matches schema 2.8 §12 default",
    sep = "\n"), "\n")
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)

CANDIDATES_TSV <- NULL
CHROM          <- NULL
OUT_DIR        <- NULL
SCAN_RADIUS_BP <- 1500000L
SCAN_WINDOW_BP <- 5000L
DOSAGE_FILE    <- NULL
FISH_REGIMES   <- NULL
SV_VCFS        <- character(0)
DISCORDANT_BED <- NULL
SAMPLES_FILE   <- NULL

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--candidates_registry" && i < length(args)) { CANDIDATES_TSV <- args[i + 1]; i <- i + 2L }
  else if (a == "--chrom"               && i < length(args)) { CHROM          <- args[i + 1]; i <- i + 2L }
  else if (a == "--out_dir"             && i < length(args)) { OUT_DIR        <- args[i + 1]; i <- i + 2L }
  else if (a == "--scan_radius_bp"      && i < length(args)) { SCAN_RADIUS_BP <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--scan_window_bp"      && i < length(args)) { SCAN_WINDOW_BP <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--dosage"              && i < length(args)) { DOSAGE_FILE    <- args[i + 1]; i <- i + 2L }
  else if (a == "--fish_regimes"        && i < length(args)) { FISH_REGIMES   <- args[i + 1]; i <- i + 2L }
  else if (a == "--discordant_bed"      && i < length(args)) { DISCORDANT_BED <- args[i + 1]; i <- i + 2L }
  else if (a == "--samples"             && i < length(args)) { SAMPLES_FILE   <- args[i + 1]; i <- i + 2L }
  else if (a == "--sv_vcf"              && i < length(args)) {
    # Consume all remaining args until the next flag (--*) — supports multiple VCFs
    j <- i + 1L
    while (j <= length(args) && !startsWith(args[j], "--")) {
      SV_VCFS <- c(SV_VCFS, args[j])
      j <- j + 1L
    }
    i <- j
  }
  else if (a == "-h" || a == "--help") { usage() }
  else { i <- i + 1L }
}

if (is.null(CANDIDATES_TSV) || is.null(CHROM) || is.null(OUT_DIR)) usage()
if (!file.exists(CANDIDATES_TSV))
  stop("candidates registry not found: ", CANDIDATES_TSV)
if (SCAN_WINDOW_BP < 1L)  stop("--scan_window_bp must be >= 1")
if (SCAN_RADIUS_BP < 0L)  stop("--scan_radius_bp must be >= 0")

cat("[M06] STEP_M06_emit_boundary_evidence.R\n")
cat("[M06]   chrom:           ", CHROM, "\n", sep = "")
cat("[M06]   scan_radius_bp:  ", SCAN_RADIUS_BP, "\n", sep = "")
cat("[M06]   scan_window_bp:  ", SCAN_WINDOW_BP, "\n", sep = "")

# =============================================================================
# Load candidates registry — filter to this chrom
# =============================================================================
cat("[M06] Loading candidates registry: ", CANDIDATES_TSV, "\n", sep = "")
reg <- fread(CANDIDATES_TSV)
need_cols <- c("candidate_id", "chrom", "start_bp", "end_bp")
miss_cols <- setdiff(need_cols, names(reg))
if (length(miss_cols) > 0)
  stop("candidates registry missing required columns: ",
       paste(miss_cols, collapse = ", "))

reg_chrom <- reg[chrom == CHROM]
n_cand <- nrow(reg_chrom)
cat("[M06]   ", n_cand, " candidates on this chrom (of ", nrow(reg), " total)\n", sep = "")
if (n_cand == 0) {
  cat("[M06] No candidates for ", CHROM, " — emitting empty layer.\n", sep = "")
}

# Coerce integer
reg_chrom[, candidate_id := as.integer(candidate_id)]
reg_chrom[, start_bp     := as.integer(start_bp)]
reg_chrom[, end_bp       := as.integer(end_bp)]

# =============================================================================
# Optional: load dosage matrix + sample list
# =============================================================================
dosage_mat   <- NULL    # numeric matrix [n_markers × n_samples]; NA permitted
marker_pos   <- NULL    # integer vector aligned to dosage_mat rows
sample_ids   <- NULL    # character vector; column order in dosage_mat

if (!is.null(DOSAGE_FILE)) {
  if (!file.exists(DOSAGE_FILE)) stop("dosage file not found: ", DOSAGE_FILE)
  cat("[M06] Loading dosage matrix: ", DOSAGE_FILE, "\n", sep = "")
  dos_dt <- fread(DOSAGE_FILE)
  cat("[M06]   raw dim: ", nrow(dos_dt), " markers x ", ncol(dos_dt) - 1, " sample cols\n", sep = "")

  # Format detection — same as STEP_M04: marker / pos / pos_bp / marker_id
  has_marker_col <- any(c("marker", "marker_id") %in% names(dos_dt))
  has_pos_col    <- any(c("pos", "pos_bp")      %in% names(dos_dt))

  if (has_pos_col) {
    pos_col <- if ("pos_bp" %in% names(dos_dt)) "pos_bp" else "pos"
    marker_pos <- as.integer(dos_dt[[pos_col]])
  } else if (has_marker_col) {
    # Parse pos from marker IDs like "<chrom>_<pos>" or "<chrom>:<pos>"
    mcol <- if ("marker_id" %in% names(dos_dt)) "marker_id" else "marker"
    mid  <- as.character(dos_dt[[mcol]])
    # Try _pos pattern first, then :pos
    parts <- strsplit(mid, "_", fixed = TRUE)
    last  <- vapply(parts, function(p) p[length(p)], character(1))
    pos_n <- suppressWarnings(as.integer(last))
    if (any(is.na(pos_n))) {
      parts2 <- strsplit(mid, ":", fixed = TRUE)
      last2  <- vapply(parts2, function(p) p[length(p)], character(1))
      pos_n  <- suppressWarnings(as.integer(last2))
    }
    if (any(is.na(pos_n)))
      stop("Could not parse positions from marker column. Use --dosage with explicit pos_bp column.")
    marker_pos <- pos_n
  } else {
    stop("Dosage TSV needs either pos/pos_bp column or marker/marker_id column.")
  }

  # Sample columns = everything except metadata
  meta_cols <- intersect(c("marker", "marker_id", "pos", "pos_bp", "chrom",
                             "missingness", "diagnostic_score"), names(dos_dt))
  sample_cols <- setdiff(names(dos_dt), meta_cols)

  if (!is.null(SAMPLES_FILE)) {
    if (!file.exists(SAMPLES_FILE)) stop("samples file not found: ", SAMPLES_FILE)
    sample_ids <- as.character(fread(SAMPLES_FILE, header = FALSE)[[1]])
    sample_ids <- sample_ids[nchar(sample_ids) > 0]
    miss_in_dosage <- setdiff(sample_ids, sample_cols)
    if (length(miss_in_dosage) > 0)
      stop("samples not found in dosage TSV header: ",
           paste(head(miss_in_dosage, 5), collapse = ", "),
           if (length(miss_in_dosage) > 5) ", ..." else "")
    sample_cols <- sample_ids
  } else {
    sample_ids <- sample_cols
  }

  # Build matrix; NA-permissive integer (FST formula uses doubles so coerce)
  dosage_mat <- as.matrix(dos_dt[, sample_cols, with = FALSE])
  storage.mode(dosage_mat) <- "double"
  # NA encoding: many TSVs use -1 for NA; respect both
  dosage_mat[dosage_mat == -1] <- NA_real_
  cat("[M06]   dosage matrix: ", nrow(dosage_mat), " x ", ncol(dosage_mat), "  (",
      sum(is.na(dosage_mat)), " NA cells)\n", sep = "")
}

# =============================================================================
# Optional: fish_regimes — per-candidate × per-sample regime assignments
# =============================================================================
# Schema (v3.80 fish_regime_calls.tsv): candidate_id, sample, regime, ...
# Regime values for K=3 are: g0 / g1 / g2  (homo1 / het / homo2 in order of
# PC1 rank). For K=2 candidates we still emit a homo1+homo2 split with no
# theta_pi_het.
regimes_dt <- NULL
if (!is.null(FISH_REGIMES)) {
  if (!file.exists(FISH_REGIMES)) stop("fish_regimes file not found: ", FISH_REGIMES)
  cat("[M06] Loading fish_regimes: ", FISH_REGIMES, "\n", sep = "")
  regimes_dt <- fread(FISH_REGIMES)
  need <- c("candidate_id", "sample", "regime")
  miss <- setdiff(need, names(regimes_dt))
  if (length(miss) > 0)
    stop("fish_regimes missing columns: ", paste(miss, collapse = ", "))
  regimes_dt[, candidate_id := as.integer(candidate_id)]
  regimes_dt[, sample := as.character(sample)]
  regimes_dt[, regime := as.character(regime)]
  cat("[M06]   ", nrow(regimes_dt), " fish-regime calls\n", sep = "")
}

# =============================================================================
# Optional: load SV VCFs
# =============================================================================
# Helper: parse a VCF for SV anchors. Returns data.table with columns
# kind, ct, pos_bp, qual. Handles DELLY (CT=...) and Manta (INV3/INV5
# in INFO or in ALT field) — all SVTYPE in {INV, DEL, BND}.
# Defined BEFORE its use site below — R hoists top-level function defs
# only when written in the same evaluation order, so we put the def
# first to be safe.
parse_sv_vcf <- function(path, chrom_filter) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  recs <- list()
  while (length(line <- readLines(con, n = 1L, warn = FALSE)) == 1L) {
    if (nchar(line) == 0L) next
    if (substr(line, 1L, 1L) == "#") next   # header
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 8L) next
    chr <- f[1]
    if (chr != chrom_filter) next
    pos <- suppressWarnings(as.integer(f[2]))
    if (is.na(pos)) next
    qual_str <- f[6]
    qual <- suppressWarnings(as.numeric(qual_str))
    if (is.na(qual)) qual <- NA_real_
    info <- f[8]
    alt  <- f[5]
    # SVTYPE
    svtype <- NA_character_
    m <- regmatches(info, regexpr("SVTYPE=[^;]+", info))
    if (length(m) > 0) svtype <- sub("SVTYPE=", "", m)
    if (is.na(svtype)) {
      # Manta INV3/INV5 — parse from ALT or INFO
      if (grepl("INV", alt) || grepl("INV3|INV5", info)) svtype <- "INV"
    }
    if (is.na(svtype) || !(svtype %in% c("INV", "DEL", "BND"))) next
    # CT (DELLY) — captures connection type for inversions
    ct <- NA_character_
    mc <- regmatches(info, regexpr("CT=[^;]+", info))
    if (length(mc) > 0) ct <- sub("CT=", "", mc)
    # Manta INV3/INV5 detection
    kind <- if (grepl("INV3", info, fixed = TRUE)) "Manta_INV3"
            else if (grepl("INV5", info, fixed = TRUE)) "Manta_INV5"
            else if (svtype == "INV" && grepl("DELLY", info, fixed = TRUE)) "DELLY_INV"
            else if (svtype == "INV") paste0(svtype)
            else if (svtype == "DEL") "DEL"
            else if (svtype == "BND") "BND"
            else svtype
    recs[[length(recs) + 1L]] <- list(
      kind   = kind,
      ct     = ct,
      pos_bp = pos,
      qual   = qual
    )
  }
  if (length(recs) == 0L) return(data.table())
  rbindlist(recs, use.names = TRUE, fill = TRUE)
}

sv_anchors_all <- list()
if (length(SV_VCFS) > 0L) {
  cat("[M06] Loading SV VCFs (", length(SV_VCFS), " files)\n", sep = "")
  for (vcf in SV_VCFS) {
    if (!file.exists(vcf)) {
      cat("[M06]   WARNING: SV VCF not found, skipping: ", vcf, "\n", sep = "")
      next
    }
    parsed <- tryCatch(parse_sv_vcf(vcf, CHROM),
                       error = function(e) {
                         cat("[M06]   WARNING: failed to parse ", vcf, ": ", conditionMessage(e), "\n", sep = "")
                         data.table()
                       })
    if (nrow(parsed) > 0)
      sv_anchors_all[[vcf]] <- parsed
  }
  total_anchors <- sum(vapply(sv_anchors_all, nrow, integer(1)))
  cat("[M06]   ", total_anchors, " SV anchors across ", length(sv_anchors_all),
      " VCF(s) for ", CHROM, "\n", sep = "")
}

# =============================================================================
# Optional: discordant-pair BED
# =============================================================================
disc_bed_dt <- NULL
if (!is.null(DISCORDANT_BED)) {
  if (!file.exists(DISCORDANT_BED)) {
    cat("[M06]   WARNING: discordant BED not found, skipping: ", DISCORDANT_BED, "\n", sep = "")
  } else {
    cat("[M06] Loading discordant-pair BED: ", DISCORDANT_BED, "\n", sep = "")
    bed_raw <- fread(DISCORDANT_BED, header = FALSE)
    if (ncol(bed_raw) < 4L)
      stop("discordant BED must be 4-col: chrom start end count")
    setnames(bed_raw, 1:4, c("chr", "start", "end", "count"))
    disc_bed_dt <- bed_raw[chr == CHROM]
    cat("[M06]   ", nrow(disc_bed_dt), " discordant rows on this chrom\n", sep = "")
  }
}

# =============================================================================
# Per-window FST (Hudson 1992) and per-window θπ helpers
# =============================================================================
# Hudson FST per marker comparing two regimes A and B with allele freqs pA, pB:
#   numer = (pA - pB)^2 - pA*(1-pA)/(nA-1) - pB*(1-pB)/(nB-1)
#   denom = pA*(1-pB) + pB*(1-pA)
#   fst   = numer / denom
# Window FST = mean(numer) / mean(denom)  (ratio of averages, not avg of ratios).
hudson_fst_window <- function(dos_window, idx_A, idx_B) {
  if (length(idx_A) < 2 || length(idx_B) < 2) return(NA_real_)
  numer_sum <- 0.0
  denom_sum <- 0.0
  n_used <- 0L
  for (mi in seq_len(nrow(dos_window))) {
    rowA <- dos_window[mi, idx_A]
    rowB <- dos_window[mi, idx_B]
    rowA <- rowA[!is.na(rowA)]
    rowB <- rowB[!is.na(rowB)]
    nA <- length(rowA); nB <- length(rowB)
    if (nA < 2 || nB < 2) next
    pA <- mean(rowA) / 2
    pB <- mean(rowB) / 2
    if (!is.finite(pA) || !is.finite(pB)) next
    numer_marker <- (pA - pB)^2 - pA * (1 - pA) / (nA - 1) - pB * (1 - pB) / (nB - 1)
    denom_marker <- pA * (1 - pB) + pB * (1 - pA)
    if (!is.finite(denom_marker) || denom_marker <= 0) next
    numer_sum <- numer_sum + numer_marker
    denom_sum <- denom_sum + denom_marker
    n_used <- n_used + 1L
  }
  if (n_used < 1 || denom_sum <= 0) return(NA_real_)
  numer_sum / denom_sum
}

# θπ within a regime for a window — mean per-marker pairwise diversity with
# sample-size correction. Per marker: 2 * p * (1-p) * n / (n-1) where p is
# allele freq, n is sample count after NA-drop. Window θπ is mean across
# markers in the window.
theta_pi_window <- function(dos_window, idx_set) {
  if (length(idx_set) < 2) return(NA_real_)
  vals <- numeric(0)
  for (mi in seq_len(nrow(dos_window))) {
    row <- dos_window[mi, idx_set]
    row <- row[!is.na(row)]
    n <- length(row)
    if (n < 2) next
    p <- mean(row) / 2
    if (!is.finite(p)) next
    vals <- c(vals, 2 * p * (1 - p) * n / (n - 1))
  }
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

# =============================================================================
# Build boundary_evidence list — ONE entry per candidate
# =============================================================================
chrom_max_bp <- if (!is.null(marker_pos)) max(marker_pos, na.rm = TRUE) else NA_integer_
boundary_evidence <- vector("list", n_cand)

for (ci in seq_len(n_cand)) {
  cand <- reg_chrom[ci]
  cand_id <- cand$candidate_id
  span <- cand$end_bp - cand$start_bp
  # Schema "huge candidate" expansion rule (matches scrubber)
  radius <- SCAN_RADIUS_BP
  if (span > 3000000L) radius <- max(radius, as.integer(span * 0.5))
  scan_start <- max(0L, cand$start_bp - radius)
  scan_end   <- if (is.na(chrom_max_bp)) cand$end_bp + radius else min(chrom_max_bp, cand$end_bp + radius)
  if (scan_end <= scan_start) {
    cat("[M06]   cand ", cand_id, ": skipping (degenerate scan range)\n", sep = "")
    next
  }
  # Number of windows
  n_win <- as.integer(floor((scan_end - scan_start) / SCAN_WINDOW_BP))
  if (n_win < 1L) n_win <- 1L
  # Window starts (inclusive): scan_start, scan_start+win, ..., scan_start+(n-1)*win
  win_starts <- scan_start + (seq_len(n_win) - 1L) * SCAN_WINDOW_BP
  win_ends   <- win_starts + SCAN_WINDOW_BP - 1L

  cat("[M06]   cand ", cand_id, ": scan ", scan_start, "..", scan_end,
      " (", n_win, " windows of ", SCAN_WINDOW_BP, " bp)\n", sep = "")

  tracks <- list()

  # -------- FST + θπ tracks (require dosage + regime assignments) --------
  if (!is.null(dosage_mat) && !is.null(regimes_dt)) {
    cand_regimes <- regimes_dt[candidate_id == cand_id]
    if (nrow(cand_regimes) > 0L) {
      # Map sample names to dosage column indices
      sm <- match(cand_regimes$sample, sample_ids)
      cand_regimes[, idx := sm]
      # Drop fish whose sample isn't in the dosage matrix
      cand_regimes <- cand_regimes[!is.na(idx)]
      idx_homo1 <- cand_regimes[regime == "g0", idx]
      idx_het   <- cand_regimes[regime == "g1", idx]
      idx_homo2 <- cand_regimes[regime == "g2", idx]
      n0 <- length(idx_homo1); n1 <- length(idx_het); n2 <- length(idx_homo2)
      cat("[M06]     regimes: g0 n=", n0, "  g1 n=", n1, "  g2 n=", n2, "\n", sep = "")

      # Marker filter: in scan range
      m_in <- which(marker_pos >= scan_start & marker_pos <= scan_end)
      cat("[M06]     ", length(m_in), " markers in scan range\n", sep = "")

      if (length(m_in) > 0L && n0 >= 2 && n2 >= 2) {
        # Allocate output vectors
        fst_v        <- rep(NA_real_, n_win)
        tp_homo1_v   <- if (n0 >= 2) rep(NA_real_, n_win) else NULL
        tp_het_v     <- if (n1 >= 2) rep(NA_real_, n_win) else NULL
        tp_homo2_v   <- if (n2 >= 2) rep(NA_real_, n_win) else NULL

        # Bin markers by window. Sort ascending first (typically already so).
        mp <- marker_pos[m_in]
        ord <- order(mp)
        m_in_ord <- m_in[ord]
        mp_ord   <- mp[ord]

        # Compute window indices for each marker (1-based)
        win_idx_for_marker <- as.integer(floor((mp_ord - scan_start) / SCAN_WINDOW_BP)) + 1L
        win_idx_for_marker[win_idx_for_marker > n_win] <- n_win

        # Walk window by window
        for (wi in seq_len(n_win)) {
          mi_in_win <- m_in_ord[win_idx_for_marker == wi]
          if (length(mi_in_win) == 0L) next
          dos_w <- dosage_mat[mi_in_win, , drop = FALSE]
          fst_v[wi] <- hudson_fst_window(dos_w, idx_homo1, idx_homo2)
          if (!is.null(tp_homo1_v)) tp_homo1_v[wi] <- theta_pi_window(dos_w, idx_homo1)
          if (!is.null(tp_het_v))   tp_het_v[wi]   <- theta_pi_window(dos_w, idx_het)
          if (!is.null(tp_homo2_v)) tp_homo2_v[wi] <- theta_pi_window(dos_w, idx_homo2)
        }

        tracks$fst <- round(fst_v, 6)
        if (!is.null(tp_homo1_v)) tracks$theta_pi_homo1 <- round(tp_homo1_v, 6)
        if (!is.null(tp_het_v))   tracks$theta_pi_het   <- round(tp_het_v, 6)
        if (!is.null(tp_homo2_v)) tracks$theta_pi_homo2 <- round(tp_homo2_v, 6)
      } else {
        cat("[M06]     skipping FST/θπ: insufficient markers or regime sizes\n")
      }
    } else {
      cat("[M06]     no fish_regime calls for cand ", cand_id, "; skipping FST/θπ\n", sep = "")
    }
  }

  # -------- discordant_pair_pileup --------
  if (!is.null(disc_bed_dt) && nrow(disc_bed_dt) > 0L) {
    # Bin BED rows into scan windows. Each BED row contributes its `count`
    # to the window covering its midpoint (or split proportionally? — the
    # simpler rule is: midpoint binning, which matches typical RD/coverage
    # tooling).
    in_scan <- disc_bed_dt[start <= scan_end & end >= scan_start]
    if (nrow(in_scan) > 0L) {
      mid <- as.integer((in_scan$start + in_scan$end) / 2L)
      win_idx <- as.integer(floor((mid - scan_start) / SCAN_WINDOW_BP)) + 1L
      win_idx[win_idx < 1L] <- 1L
      win_idx[win_idx > n_win] <- n_win
      pile_v <- integer(n_win)
      for (k in seq_along(win_idx)) {
        pile_v[win_idx[k]] <- pile_v[win_idx[k]] + as.integer(in_scan$count[k])
      }
      tracks$discordant_pair_pileup <- pile_v
    }
  }

  # -------- sv_anchors --------
  if (length(sv_anchors_all) > 0L) {
    all_anchors <- rbindlist(sv_anchors_all, use.names = TRUE, fill = TRUE)
    if (nrow(all_anchors) > 0L) {
      in_scan <- all_anchors[pos_bp >= scan_start & pos_bp <= scan_end]
      if (nrow(in_scan) > 0L) {
        anchor_list <- vector("list", nrow(in_scan))
        for (k in seq_len(nrow(in_scan))) {
          ak <- list(kind = in_scan$kind[k], pos_bp = as.integer(in_scan$pos_bp[k]))
          if (!is.na(in_scan$ct[k])) ak$ct <- in_scan$ct[k]
          if (!is.na(in_scan$qual[k])) ak$qual <- as.numeric(in_scan$qual[k])
          anchor_list[[k]] <- ak
        }
        tracks$sv_anchors <- anchor_list
      }
    }
  }

  # Assemble candidate entry. If tracks{} is empty, skip the candidate
  # — emitting an empty boundary_evidence row would still register the
  # layer with the scrubber but provide no signal.
  if (length(tracks) == 0L) {
    cat("[M06]     no tracks computable for cand ", cand_id, "; skipping\n", sep = "")
    next
  }

  boundary_evidence[[ci]] <- list(
    candidate_id   = as.integer(cand_id),
    chrom          = CHROM,
    scan_start_bp  = as.integer(scan_start),
    scan_end_bp    = as.integer(scan_end),
    scan_window_bp = as.integer(SCAN_WINDOW_BP),
    tracks         = tracks
  )
}

# Drop NULL entries (candidates we skipped)
boundary_evidence <- boundary_evidence[!vapply(boundary_evidence, is.null, logical(1))]

# =============================================================================
# Assemble + write JSON
# =============================================================================
chrom_dir <- file.path(OUT_DIR, CHROM)
dir.create(chrom_dir, recursive = TRUE, showWarnings = FALSE)
out_json <- file.path(chrom_dir, paste0(CHROM, "_boundary_evidence.json"))

out <- list(
  schema_version    = 2L,
  `_layers_present` = list("boundary_evidence"),
  boundary_evidence = boundary_evidence,
  `_generated_at`   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  `_generator`      = "STEP_M06_emit_boundary_evidence.R"
)

cat("[M06] Serializing JSON ...\n")
t_ser <- proc.time()
# auto_unbox=TRUE → length-1 vectors emit as scalars; na="null" → NA → JSON null
# digits=NA preserves the round() we applied above
json_str <- jsonlite::toJSON(out, auto_unbox = TRUE, na = "null",
                             digits = NA, pretty = FALSE)
writeLines(json_str, out_json)
cat("[M06]   serialized in ", round((proc.time() - t_ser)[3], 1), "s\n", sep = "")

f_size_kb <- round(file.info(out_json)$size / 1024, 1)
cat("\n[M06] === DONE ===\n")
cat("[M06] Output: ", out_json, "\n", sep = "")
cat("[M06] Size:   ", f_size_kb, " KB\n", sep = "")
cat("[M06] Candidates with at least one track: ", length(boundary_evidence), " of ", n_cand, "\n", sep = "")
cat("[M06] schema_version: 2\n")
cat("\n")
cat("[M06] Drop this JSON into the scrubber alongside other enrichment files.\n")
cat("[M06] The scrubber registers it as the 'boundary_evidence' layer; v3.97's\n")
cat("[M06] auto-propose algorithm picks up the additional fst / theta_pi /\n")
cat("[M06] discordant_pile / sv_anchor tracks automatically.\n")
