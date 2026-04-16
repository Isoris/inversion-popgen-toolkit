#!/usr/bin/env Rscript

# =============================================================================
# STEP17_build_candidate_scoring_table.R
#
# Integrated per-sample per-candidate scoring table combining:
#   - Regional PCA groups + coordinates (STEP15)
#   - Block coherence metrics (STEP16)
#   - Theta/tP summary within candidate (from STEP12 plotC or STEP14)
#   - Optional ancestry/Q data
#
# Produces evidence scores and provisional major-state labels:
#   FULL_A / HALF / FULL_B / NOISE / COMPLEX
#
# Usage:
#   Rscript STEP17_build_candidate_scoring_table.R \
#     <group_table> <coherence_table> <candidate_table> \
#     <outdir> [ancestry_table=NONE] [theta_overlap_dir=NONE]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript STEP17_... <group_table> <coherence_table> <candidate_table> <outdir> [ancestry=NONE] [theta_dir=NONE]")
}

group_file     <- args[1]
coherence_file <- args[2]
candidate_file <- args[3]
outdir         <- args[4]
ancestry_file  <- if (length(args) >= 5 && args[5] != "NONE") args[5] else NA_character_
theta_dir      <- if (length(args) >= 6 && args[6] != "NONE") args[6] else NA_character_

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Read inputs ────────────────────────────────────────────────────────────
groups <- fread(group_file)
coherence <- fread(coherence_file)
cand <- fread(candidate_file)

stopifnot("candidate_id" %in% names(groups))
stopifnot("candidate_id" %in% names(coherence))

message("[INFO] Group records: ", nrow(groups))
message("[INFO] Coherence records: ", nrow(coherence))

# ── Merge groups + coherence ──────────────────────────────────────────────
scoring <- merge(groups, coherence,
                 by = c("candidate_id", "chrom", "sample"),
                 all = TRUE)

# ── Add candidate coordinates ────────────────────────────────────────────
if (all(c("start_bp", "end_bp") %in% names(cand))) {
  scoring <- merge(scoring,
                   cand[, .(candidate_id, start_bp, end_bp)],
                   by = "candidate_id", all.x = TRUE)
}

# ── Optional ancestry merge ──────────────────────────────────────────────
if (!is.na(ancestry_file) && file.exists(ancestry_file)) {
  ancestry <- fread(ancestry_file)
  message("[INFO] Ancestry records: ", nrow(ancestry))

  if ("sample" %in% names(ancestry)) {
    # Keep Q columns and dominant ancestry label
    q_cols <- grep("^Q[0-9]", names(ancestry), value = TRUE)
    keep <- unique(c("sample", "cluster_label", "qmax", q_cols))
    keep <- intersect(keep, names(ancestry))
    if (length(keep) > 1) {
      scoring <- merge(scoring, ancestry[, ..keep], by = "sample", all.x = TRUE)
    }
  }
}

# ── Compute evidence scores ──────────────────────────────────────────────
# Evidence_coherence: high if block_coherence high, switch_rate low, support_margin high
scoring[, evidence_coherence := fifelse(
  is.finite(block_coherence) & is.finite(support_margin),
  (block_coherence * 0.4 + support_margin * 0.4 +
   fifelse(is.finite(long_range_agreement), long_range_agreement * 0.2, 0)),
  NA_real_
)]

# Evidence_intermediate: high if support is balanced (low margin, intermediate het)
scoring[, evidence_intermediate := fifelse(
  is.finite(support_margin) & is.finite(regional_het),
  (1 - support_margin) * 0.5 + pmin(regional_het / max(regional_het, na.rm = TRUE), 1) * 0.5,
  NA_real_
)]

# Evidence_noise: high if low coherence, high switch rate
scoring[, evidence_noise := fifelse(
  is.finite(block_coherence) & is.finite(switch_rate),
  (1 - block_coherence) * 0.5 + switch_rate * 0.5,
  NA_real_
)]

# ── Provisional major state ─────────────────────────────────────────────
# Logic:
#   Group 1 or 3 + high coherence → FULL_A or FULL_B
#   Group 2 + high intermediate score → HALF
#   High noise → NOISE
#   Otherwise → COMPLEX
scoring[, provisional_major_state := {
  state <- "COMPLEX"

  if (!is.na(evidence_coherence) && !is.na(group_norm)) {
    if (group_norm == 1 && evidence_coherence >= 0.5 && !is.na(mean_A_support) && mean_A_support > 0.6) {
      state <- "FULL_A"
    } else if (group_norm == 3 && evidence_coherence >= 0.5 && !is.na(mean_B_support) && mean_B_support > 0.6) {
      state <- "FULL_B"
    } else if (group_norm == 2 && !is.na(evidence_intermediate) && evidence_intermediate >= 0.4) {
      state <- "HALF"
    }
  }

  if (!is.na(evidence_noise) && evidence_noise >= 0.7) {
    state <- "NOISE"
  }

  state
}, by = seq_len(nrow(scoring))]

# Confidence
scoring[, provisional_confidence := fifelse(
  provisional_major_state %in% c("FULL_A", "FULL_B") & evidence_coherence >= 0.7, "high",
  fifelse(
    provisional_major_state == "HALF" & evidence_intermediate >= 0.5, "medium",
    fifelse(
      provisional_major_state == "NOISE", "low",
      "low"
    )
  )
)]

# Complexity flag
scoring[, complexity_flag := fifelse(
  provisional_major_state == "COMPLEX" |
  (!is.na(evidence_noise) & evidence_noise >= 0.5 & provisional_major_state != "NOISE"),
  TRUE, FALSE
)]

# ── Write outputs ────────────────────────────────────────────────────────
fwrite(scoring,
       file.path(outdir, "inversion_candidate_scoring_table.tsv.gz"),
       sep = "\t")

# Per-candidate summary
cand_summary <- scoring[, .(
  n_samples    = .N,
  n_FULL_A     = sum(provisional_major_state == "FULL_A"),
  n_HALF       = sum(provisional_major_state == "HALF"),
  n_FULL_B     = sum(provisional_major_state == "FULL_B"),
  n_NOISE      = sum(provisional_major_state == "NOISE"),
  n_COMPLEX    = sum(provisional_major_state == "COMPLEX"),
  mean_coherence = round(mean(evidence_coherence, na.rm = TRUE), 4),
  frac_high_confidence = round(mean(provisional_confidence == "high"), 4),
  complexity_rate = round(mean(complexity_flag), 4)
), by = .(candidate_id, chrom)]

# Candidate-level classification
cand_summary[, candidate_classification := fifelse(
  n_FULL_A >= 3 & n_HALF >= 3 & n_FULL_B >= 3 & mean_coherence >= 0.5 & complexity_rate < 0.3,
  "clean_inversion_like",
  fifelse(
    (n_FULL_A + n_FULL_B) >= 5 & mean_coherence >= 0.3,
    "complex_inversion_like",
    fifelse(
      complexity_rate >= 0.5 | mean_coherence < 0.3,
      "non_inversion_or_ambiguous",
      "unclassified"
    )
  )
)]

fwrite(cand_summary,
       file.path(outdir, "inversion_candidate_classification.tsv"),
       sep = "\t")

message("[DONE] Wrote:")
message("  ", file.path(outdir, "inversion_candidate_scoring_table.tsv.gz"))
message("  ", file.path(outdir, "inversion_candidate_classification.tsv"))
