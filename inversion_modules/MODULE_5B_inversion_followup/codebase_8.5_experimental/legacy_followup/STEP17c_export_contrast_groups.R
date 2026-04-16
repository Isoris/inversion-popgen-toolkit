#!/usr/bin/env Rscript

# =============================================================================
# STEP17c_export_contrast_groups.R
#
# Reads the STEP17 scoring table and exports per-candidate group files
# that serve as the SHARED grouping layer for both LD and FST modules.
#
# Outputs per candidate:
#   candidate_<id>.contrast_groups.tsv   — full table
#   candidate_<id>.FULL_A.samples.txt    — one sample per line
#   candidate_<id>.FULL_B.samples.txt
#   candidate_<id>.HALF.samples.txt
#   candidate_<id>.NOISE.samples.txt
#
# Also produces:
#   all_candidates_contrast_summary.tsv  — counts per state per candidate
#
# Usage:
#   Rscript STEP17c_export_contrast_groups.R \
#     <scoring_table.tsv.gz> <candidate_table> <outdir> \
#     [min_per_group=3]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript STEP17c_export_contrast_groups.R <scoring.tsv.gz> <candidate.tsv> <outdir> [min_per_group=3]")
}

scoring_file   <- args[1]
candidate_file <- args[2]
outdir         <- args[3]
min_per_group  <- if (length(args) >= 4) as.integer(args[4]) else 3L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

scoring <- fread(scoring_file)
cand    <- fread(candidate_file)

stopifnot(all(c("candidate_id", "sample", "provisional_major_state") %in% names(scoring)))
stopifnot(all(c("candidate_id", "chrom", "start_bp", "end_bp") %in% names(cand)))

message("[INFO] Scoring records: ", nrow(scoring))
message("[INFO] Candidates: ", nrow(cand))

# Merge candidate coordinates into scoring
if (!("start_bp" %in% names(scoring))) {
  scoring <- merge(scoring, cand[, .(candidate_id, start_bp, end_bp)],
                   by = "candidate_id", all.x = TRUE)
}

# Ensure chrom is present
if (!("chrom" %in% names(scoring))) {
  scoring <- merge(scoring, cand[, .(candidate_id, chrom)],
                   by = "candidate_id", all.x = TRUE)
}

summary_list <- list()

for (cid_val in unique(scoring$candidate_id)) {
  sub <- scoring[candidate_id == cid_val]
  chr <- sub$chrom[1]

  # Build contrast groups table
  keep_cols <- intersect(
    c("candidate_id", "chrom", "sample", "provisional_major_state",
      "provisional_confidence", "group_norm", "group_label",
      "evidence_coherence", "evidence_intermediate", "evidence_noise",
      "complexity_flag", "PC1", "PC2", "regional_het",
      "block_coherence", "support_margin", "mean_A_support", "mean_B_support"),
    names(sub)
  )

  contrast_dt <- sub[, ..keep_cols]

  # Add keep_for_contrast flag: TRUE if state is FULL_A, FULL_B, or HALF
  contrast_dt[, keep_for_contrast := provisional_major_state %in% c("FULL_A", "FULL_B", "HALF")]

  # Write contrast groups table
  fwrite(contrast_dt,
         file.path(outdir, paste0("candidate_", cid_val, ".contrast_groups.tsv")),
         sep = "\t")

  # Write per-state sample lists
  states <- c("FULL_A", "FULL_B", "HALF", "NOISE", "COMPLEX")
  state_counts <- list()

  for (st in states) {
    samples <- sub[provisional_major_state == st, sample]
    state_counts[[st]] <- length(samples)

    outfile <- file.path(outdir, paste0("candidate_", cid_val, ".", st, ".samples.txt"))
    if (length(samples) > 0) {
      writeLines(samples, outfile)
    } else {
      writeLines(character(0), outfile)
    }
  }

  # Determine if contrast is viable
  n_A <- state_counts[["FULL_A"]]
  n_B <- state_counts[["FULL_B"]]
  n_H <- state_counts[["HALF"]]

  contrast_viable <- n_A >= min_per_group & n_B >= min_per_group
  half_viable     <- n_H >= min_per_group

  summary_list[[length(summary_list) + 1]] <- data.table(
    candidate_id     = cid_val,
    chrom            = chr,
    n_FULL_A         = n_A,
    n_FULL_B         = n_B,
    n_HALF           = n_H,
    n_NOISE          = state_counts[["NOISE"]],
    n_COMPLEX        = state_counts[["COMPLEX"]],
    n_total          = nrow(sub),
    contrast_viable  = contrast_viable,
    half_viable      = half_viable
  )
}

summary_dt <- rbindlist(summary_list)

# Merge candidate coordinates
summary_dt <- merge(summary_dt,
                    cand[, .(candidate_id, start_bp, end_bp)],
                    by = "candidate_id", all.x = TRUE)

fwrite(summary_dt,
       file.path(outdir, "all_candidates_contrast_summary.tsv"),
       sep = "\t")

message("[DONE] Exported groups for ", length(unique(scoring$candidate_id)), " candidates to: ", outdir)
message("  Viable contrasts (FULL_A vs FULL_B, n>=", min_per_group, "): ",
        sum(summary_dt$contrast_viable), " / ", nrow(summary_dt))
message("  Viable HALF groups: ", sum(summary_dt$half_viable), " / ", nrow(summary_dt))
