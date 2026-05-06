#!/usr/bin/env Rscript
# =============================================================================
# cli/STEP_U06_export_ushape_json.R
# =============================================================================
# Combines window stats + raw summaries + shape scores + classes + cluster
# assignments + thresholds into output/ushape_evolution_v1.json.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(jsonlite)
})

.this_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1L]))
source(file.path(.this_dir, "..", "R", "ushape_io.R"))


opts <- parse_args(OptionParser(option_list = list(
  make_option("--candidate_table",    type = "character"),
  make_option("--window_stats_dir",   type = "character"),
  make_option("--raw_summary",        type = "character"),
  make_option("--shape_scores",       type = "character"),
  make_option("--shape_classes",      type = "character"),
  make_option("--cluster_assignments",type = "character", default = NA),
  make_option("--out_json",           type = "character"),
  make_option("--n_windows_inside",   type = "integer", default = 100L),
  make_option("--max_flank_bp",       type = "integer", default = 2000000L),
  make_option("--min_flank_bp",       type = "integer", default = 100000L),
  make_option("--min_window_bp",      type = "integer", default = 1000L),
  make_option("--edge_fraction",      type = "double",  default = 0.20),
  make_option("--center_fraction",    type = "double",  default = 0.60),
  make_option("--u_score_high",        type = "double", default = 1.5),
  make_option("--inside_flank_high",   type = "double", default = 1.5),
  make_option("--internal_peak_high",  type = "double", default = 1.5),
  make_option("--asymmetry_log2_high", type = "double", default = 1.0),
  make_option("--fst_enrichment_high", type = "double", default = 1.5),
  make_option("--oldness_min",         type = "double", default = 1.2),
  make_option("--dxy_min_inside",      type = "double", default = 1e-4),
  make_option("--flatness_max_for_flat",type = "double", default = 0.35)
)))

raw <- fread(opts$raw_summary)
scores <- fread(opts$shape_scores)
classes <- fread(opts$shape_classes)
clu <- if (!is.na(opts$cluster_assignments) && file.exists(opts$cluster_assignments))
         fread(opts$cluster_assignments) else
         data.table(candidate_id = character(0))

blocks <- list()
for (i in seq_len(nrow(raw))) {
  rr <- raw[i]
  win_path <- file.path(opts$window_stats_dir,
                        paste0(rr$candidate_id, ".window_stats.tsv"))
  if (!file.exists(win_path)) next
  win <- fread(win_path)

  sc <- scores[candidate_id == rr$candidate_id]
  if (nrow(sc) == 0L) next

  cl_row <- classes[candidate_id == rr$candidate_id]
  if (nrow(cl_row) == 0L) next
  classification <- list(
    primary_class   = cl_row$primary_class,
    secondary_class = if (is.na(cl_row$secondary_class)) NA_character_ else cl_row$secondary_class,
    confidence      = cl_row$confidence,
    reason          = cl_row$reason,
    flags           = strsplit(cl_row$flags %||% "", ";")[[1]]
  )
  if (nrow(clu) > 0L) {
    cr <- clu[candidate_id == rr$candidate_id]
    if (nrow(cr)) {
      classification$unsupervised_cluster <- as.character(cr$pam_cluster[1])
      classification$cluster_label <- cr$cluster_label[1]
    }
  }

  groups <- list(HOMO_1 = rep("", rr$n_homo1),
                 HET    = rep("", rr$n_het %||% 0),
                 HOMO_2 = rep("", rr$n_homo2))

  # matched-background z-scores (if U03 was run and merged into shape_scores)
  mbg_cols <- intersect(c("z_dxy_inside","z_fst_inside","z_dxy_u_score",
                          "z_dxy_internal_peak","bg_n_pool","bg_mode"),
                        names(sc))
  mbg <- if (length(mbg_cols)) sc[, ..mbg_cols] else NULL

  blocks[[length(blocks) + 1L]] <-
    ushape_candidate_block(win, rr, sc, classification, groups, matched_bg = mbg)
}

write_ushape_json(opts$out_json, blocks, params = opts,
                  created_by = "STEP_U06_export_ushape_json.R")
message("[U06] wrote ", length(blocks), " candidates to ", opts$out_json)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
