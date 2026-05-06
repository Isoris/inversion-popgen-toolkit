#!/usr/bin/env Rscript
# =============================================================================
# cli/STEP_U02_score_and_classify_shapes.R
# =============================================================================
# Reads candidate_raw_summary.tsv + per-candidate window_stats files,
# computes shape scores + age-gated classification, writes:
#   output/candidate_shape_scores.tsv
#   output/candidate_shape_classes.tsv
#   output/classification_summary.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(data.table)
})

.this_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1L]))
source(file.path(.this_dir, "..", "R", "ushape_lib.R"))
source(file.path(.this_dir, "..", "R", "ushape_classify.R"))


opts <- parse_args(OptionParser(option_list = list(
  make_option("--raw_summary",        type = "character"),
  make_option("--window_stats_dir",   type = "character"),
  make_option("--out_dir",            type = "character"),
  make_option("--u_score_high",        type = "double",  default = 1.5),
  make_option("--inside_flank_high",   type = "double",  default = 1.5),
  make_option("--internal_peak_high",  type = "double",  default = 1.5),
  make_option("--asymmetry_log2_high", type = "double",  default = 1.0),
  make_option("--fst_enrichment_high", type = "double",  default = 1.5),
  make_option("--oldness_min",         type = "double",  default = 1.2),
  make_option("--dxy_min_inside",      type = "double",  default = 1e-4),
  make_option("--min_inside_windows_for_shape", type = "integer", default = 8L),
  make_option("--flatness_max_for_flat", type = "double", default = 0.35)
)))

raw <- fread(opts$raw_summary)
stopifnot("candidate_id" %in% names(raw))

scores_list <- list(); class_list <- list()
for (i in seq_len(nrow(raw))) {
  rr <- raw[i]
  win_path <- file.path(opts$window_stats_dir,
                        paste0(rr$candidate_id, ".window_stats.tsv"))
  if (!file.exists(win_path)) {
    message("[U02] missing windows for ", rr$candidate_id); next
  }
  win <- fread(win_path)
  sc <- score_candidate(rr, win, opts)
  scores_list[[length(scores_list) + 1L]] <- sc

  inside_n <- nrow(win[zone %in% c("left_edge","center","right_edge")])
  cl <- classify_candidate(sc, rr, opts, inside_n)
  class_list[[length(class_list) + 1L]] <- data.table(
    candidate_id = rr$candidate_id,
    primary_class = cl$primary_class,
    secondary_class = cl$secondary_class,
    confidence = cl$confidence,
    reason = cl$reason,
    flags = paste(cl$flags, collapse = ";")
  )
}

if (length(scores_list)) {
  fwrite(rbindlist(scores_list, fill = TRUE),
         file.path(opts$out_dir, "candidate_shape_scores.tsv"), sep = "\t")
  fwrite(rbindlist(class_list, fill = TRUE),
         file.path(opts$out_dir, "candidate_shape_classes.tsv"), sep = "\t")

  # tiny summary
  cls <- rbindlist(class_list)
  smry <- cls[, .N, by = primary_class][order(-N)]
  fwrite(smry, file.path(opts$out_dir, "classification_summary.tsv"), sep = "\t")
}
message("[U02] scored ", length(scores_list), " candidates")
