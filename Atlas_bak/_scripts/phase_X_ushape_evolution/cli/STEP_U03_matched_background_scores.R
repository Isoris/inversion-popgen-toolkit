#!/usr/bin/env Rscript
# =============================================================================
# cli/STEP_U03_matched_background_scores.R
# =============================================================================
# Local-flank-resampling matched background. NO genomewide windows needed.
# This script is a thin wrapper around flank_resample_zscores() in
# R/ushape_lib.R — the same helper is called inline by the popstats
# server endpoint, so offline batch and dynamic atlas fetches give
# identical numbers.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(data.table)
})

.this_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(),
                                             value = TRUE)[1L]))
source(file.path(.this_dir, "..", "R", "ushape_lib.R"))


opts <- parse_args(OptionParser(option_list = list(
  make_option("--candidate_table",  type = "character"),
  make_option("--window_stats_dir", type = "character"),
  make_option("--out_dir",          type = "character"),
  make_option("--n_random",         type = "integer", default = 1000L),
  make_option("--seed",             type = "integer", default = 42L)
)))

cands <- fread(opts$candidate_table)
scores_path <- file.path(opts$out_dir, "candidate_shape_scores.tsv")
if (!file.exists(scores_path)) stop("run U02 first: ", scores_path, " missing")
scores <- fread(scores_path)

# length-class assignment for downstream filters
.lclass <- function(L) {
  if (L < 100e3)   "small"
  else if (L < 1e6)   "medium"
  else if (L < 1e7)   "large"
  else "mega"
}
cands[, length_bp := end_bp - start_bp + 1L]
cands[, length_class := vapply(length_bp, .lclass, character(1L))]

z_rows <- list()
for (i in seq_len(nrow(cands))) {
  cid <- cands$candidate_id[i]
  win_path <- file.path(opts$window_stats_dir, paste0(cid, ".window_stats.tsv"))
  if (!file.exists(win_path)) next
  cw <- fread(win_path)
  z <- flank_resample_zscores(cw, n_random = opts$n_random, seed = opts$seed)
  z[, candidate_id := cid]
  z_rows[[length(z_rows) + 1L]] <- z
}

zdt <- if (length(z_rows)) rbindlist(z_rows, fill = TRUE) else
        data.table(candidate_id = character(0))

if (nrow(zdt) > 0L) {
  setcolorder(zdt, c("candidate_id",
                     setdiff(names(zdt), "candidate_id")))
  scores <- merge(scores, zdt, by = "candidate_id", all.x = TRUE)
  scores[, matched_background_mode := "flank_resample"]
  fwrite(scores, scores_path, sep = "\t")
  fwrite(zdt, file.path(opts$out_dir, "candidate_matched_background.tsv"),
         sep = "\t")
}

fwrite(cands[, .(candidate_id, length_class, length_bp)],
       file.path(opts$out_dir, "candidate_length_class.tsv"), sep = "\t")

message("[U03] flank-resample z-scored ", nrow(zdt), " candidates")
