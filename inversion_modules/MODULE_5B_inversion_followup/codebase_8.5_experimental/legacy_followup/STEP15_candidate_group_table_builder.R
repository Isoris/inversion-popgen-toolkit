#!/usr/bin/env Rscript

# =============================================================================
# STEP15_candidate_group_table_builder.R
#
# Parse all STEP12 regional_pca_samples.tsv.gz outputs into unified tables:
#   1. Long-format: sample × candidate → group assignment + PCA coords
#   2. Group summary: per-candidate group counts
#   3. Wide matrix: rows=samples, cols=candidates, values=group code
#
# Usage:
#   Rscript STEP15_candidate_group_table_builder.R \
#     <step12_dir> <candidate_table> <outdir>
#
# Inputs:
#   step12_dir       — directory containing STEP12_*.regional_pca_samples.tsv.gz
#   candidate_table  — candidate_regions.tsv.gz or supp_candidate_master_table.tsv
#   outdir           — output directory
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript STEP15_candidate_group_table_builder.R <step12_dir> <candidate_table> <outdir>")
}

step12_dir     <- args[1]
candidate_file <- args[2]
outdir         <- args[3]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Read candidate table ───────────────────────────────────────────────────
cand <- fread(candidate_file)
stopifnot("candidate_id" %in% names(cand), "chrom" %in% names(cand))
message("[INFO] Candidates: ", nrow(cand))

# ── Find all STEP12 regional_pca_samples files ────────────────────────────
pca_files <- list.files(step12_dir, pattern = "\\.regional_pca_samples\\.tsv\\.gz$",
                        full.names = TRUE, recursive = FALSE)
if (length(pca_files) == 0) stop("No STEP12 .regional_pca_samples.tsv.gz found in: ", step12_dir)
message("[INFO] Found ", length(pca_files), " STEP12 PCA files")

# ── Parse each file ──────────────────────────────────────────────────────
all_long <- list()

for (f in pca_files) {
  bn <- basename(f)

  # Extract candidate_id from filename: STEP12_<chrom>.candidate_<id>.regional_pca_samples.tsv.gz
  m <- regmatches(bn, regexpr("candidate_([0-9]+)", bn))
  if (length(m) == 0) {
    message("[WARN] Cannot parse candidate_id from: ", bn)
    next
  }
  cid <- as.integer(sub("candidate_", "", m))

  # Extract chromosome
  m_chr <- regmatches(bn, regexpr("STEP12_([^.]+)", bn))
  chr_val <- if (length(m_chr) > 0) sub("STEP12_", "", m_chr) else NA_character_

  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) {
    message("[WARN] Empty or unreadable: ", bn)
    next
  }

  # Required columns
  req <- c("sample", "PC1", "PC2")
  if (!all(req %in% names(dt))) {
    message("[WARN] Missing required columns in: ", bn)
    next
  }

  dt[, candidate_id := cid]
  dt[, chrom := chr_val]
  dt[, source_file := bn]

  # Normalize group columns
  if ("ordered_group" %in% names(dt)) {
    dt[, group_norm := ordered_group]
  } else if ("raw_group" %in% names(dt)) {
    # Order by PC1 centroid
    grp_ord <- dt[, .(pc1_mean = mean(PC1)), by = raw_group][order(pc1_mean)]
    grp_ord[, group_norm := seq_len(.N)]
    dt <- merge(dt, grp_ord[, .(raw_group, group_norm)], by = "raw_group", all.x = TRUE)
  } else {
    dt[, group_norm := NA_integer_]
  }

  # Keep group_label if present
  if (!("group_label" %in% names(dt))) {
    dt[, group_label := paste0("G", group_norm)]
  }

  # Standardize regional_het if present
  if (!("regional_het" %in% names(dt))) dt[, regional_het := NA_real_]
  if (!("regional_hap_score" %in% names(dt))) dt[, regional_hap_score := NA_real_]

  keep_cols <- c("candidate_id", "chrom", "sample", "group_norm", "group_label",
                 "PC1", "PC2", "regional_het", "regional_hap_score", "source_file")
  keep_cols <- intersect(keep_cols, names(dt))

  all_long[[length(all_long) + 1]] <- dt[, ..keep_cols]
}

if (length(all_long) == 0) stop("No STEP12 files parsed successfully")

long_dt <- rbindlist(all_long, fill = TRUE)
message("[INFO] Parsed ", nrow(long_dt), " sample-candidate records across ",
        uniqueN(long_dt$candidate_id), " candidates")

# ── 1. Long-format table ─────────────────────────────────────────────────
fwrite(long_dt,
       file.path(outdir, "inversion_candidate_sample_groups.tsv.gz"),
       sep = "\t")

# ── 2. Group summary ────────────────────────────────────────────────────
group_summary <- long_dt[, {
  tab <- table(group_norm)
  n_groups <- length(tab)
  list(
    group1_n = as.integer(tab["1"]),
    group2_n = as.integer(tab["2"]),
    group3_n = as.integer(tab["3"]),
    n_samples = .N,
    n_groups_detected = n_groups
  )
}, by = .(candidate_id, chrom)]

# Replace NA counts with 0
for (col in c("group1_n", "group2_n", "group3_n")) {
  group_summary[is.na(get(col)), (col) := 0L]
}

# Merge with candidate coordinates
if (all(c("start_bp", "end_bp") %in% names(cand))) {
  group_summary <- merge(group_summary,
                         cand[, .(candidate_id, start_bp, end_bp)],
                         by = "candidate_id", all.x = TRUE)
}

fwrite(group_summary,
       file.path(outdir, "inversion_candidate_group_summary.tsv"),
       sep = "\t")

# ── 3. Wide matrix (rows=samples, cols=candidate_id, values=group_norm) ──
wide_dt <- dcast(long_dt, sample ~ candidate_id, value.var = "group_norm")
fwrite(wide_dt,
       file.path(outdir, "inversion_candidate_group_matrix_wide.tsv.gz"),
       sep = "\t")

# ── 4. Numeric matrix for R (.Rtab) ────────────────────────────────────
mat <- as.matrix(wide_dt[, -1])
rownames(mat) <- wide_dt$sample
write.table(mat,
            file.path(outdir, "inversion_candidate_group_matrix.Rtab"),
            sep = "\t", quote = FALSE)

message("[DONE] Wrote to: ", outdir)
message("  inversion_candidate_sample_groups.tsv.gz")
message("  inversion_candidate_group_summary.tsv")
message("  inversion_candidate_group_matrix_wide.tsv.gz")
message("  inversion_candidate_group_matrix.Rtab")
