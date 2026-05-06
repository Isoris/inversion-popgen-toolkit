#!/usr/bin/env Rscript
# =============================================================================
# cli/STEP_U01_compute_window_stats.R
# =============================================================================
# Offline batch entry point for U-shape window stats.
#
# For each candidate, reads <chrom>.dosage.tsv.gz, intersects with the
# karyotype-group table, and emits per-window pi/dXY/FST/AFD/SNP-class counts.
#
# When the popstats server is up, the server endpoint reuses the same
# scoring/classification on top of region_popstats output — this CLI
# is the offline fallback (dosage-only mode).
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# locate the lib relative to this script
.this_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(),
                                             value = TRUE)[1L]))
source(file.path(.this_dir, "..", "R", "ushape_lib.R"))


opts <- parse_args(OptionParser(option_list = list(
  make_option("--candidate_table", type = "character"),
  make_option("--group_table",     type = "character"),
  make_option("--dosage_dir",      type = "character"),
  make_option("--out_dir",         type = "character"),
  make_option("--n_windows_inside", type = "integer", default = 100L),
  make_option("--max_flank_bp",     type = "integer", default = 2000000L),
  make_option("--min_flank_bp",     type = "integer", default = 100000L),
  make_option("--min_window_bp",    type = "integer", default = 1000L),
  make_option("--edge_fraction",    type = "double",  default = 0.20),
  make_option("--center_fraction",  type = "double",  default = 0.60),
  make_option("--min_group_n",      type = "integer", default = 3L),
  make_option("--use_extra_groups", type = "logical", default = FALSE)
)))

stopifnot(file.exists(opts$candidate_table),
          file.exists(opts$group_table),
          dir.exists(opts$dosage_dir))
dir.create(file.path(opts$out_dir, "window_stats"),
           recursive = TRUE, showWarnings = FALSE)

cands  <- fread(opts$candidate_table)
groups <- fread(opts$group_table)
stopifnot(all(c("candidate_id","chrom","start_bp","end_bp") %in% names(cands)))
stopifnot(all(c("candidate_id","sample_id","karyotype_group") %in% names(groups)))

params <- list(
  n_windows_inside = opts$n_windows_inside,
  max_flank_bp = opts$max_flank_bp,  min_flank_bp = opts$min_flank_bp,
  min_window_bp = opts$min_window_bp,
  edge_fraction = opts$edge_fraction, center_fraction = opts$center_fraction,
  min_group_n = opts$min_group_n
)

raw_summaries <- list()
for (i in seq_len(nrow(cands))) {
  cand <- as.list(cands[i])
  cand$length_bp <- cand$end_bp - cand$start_bp + 1L

  g_rows <- groups[candidate_id == cand$candidate_id]
  if (nrow(g_rows) == 0L) {
    message("[U01] no groups for ", cand$candidate_id, "; skip")
    next
  }
  this_groups <- list(
    HOMO_1 = g_rows[karyotype_group == "HOMO_1", sample_id],
    HET    = g_rows[karyotype_group == "HET",    sample_id],
    HOMO_2 = g_rows[karyotype_group == "HOMO_2", sample_id]
  )
  if (length(this_groups$HOMO_1) < opts$min_group_n ||
      length(this_groups$HOMO_2) < opts$min_group_n) {
    message("[U01] ", cand$candidate_id, ": HOMO group < min_group_n; skip")
    next
  }

  dosage_path <- file.path(opts$dosage_dir, paste0(cand$chrom, ".dosage.tsv.gz"))
  if (!file.exists(dosage_path)) {
    message("[U01] missing dosage: ", dosage_path, "; skip")
    next
  }
  dosage_dt <- fread(dosage_path)
  win <- compute_window_stats_from_dosage(dosage_dt, cand, this_groups, params)
  if (nrow(win) == 0L) next

  fwrite(win,
         file.path(opts$out_dir, "window_stats",
                   paste0(cand$candidate_id, ".window_stats.tsv")),
         sep = "\t")
  raw_summaries[[length(raw_summaries) + 1L]] <-
    candidate_raw_summary(win, cand, this_groups)
  message("[U01] done ", cand$candidate_id,
          " (", nrow(win), " windows)")
}

if (length(raw_summaries) > 0L) {
  fwrite(rbindlist(raw_summaries, fill = TRUE),
         file.path(opts$out_dir, "candidate_raw_summary.tsv"),
         sep = "\t")
}
message("[U01] wrote ", length(raw_summaries), " candidate summaries")
