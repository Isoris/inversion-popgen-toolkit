# =============================================================================
# R/ushape_io.R
# =============================================================================
# JSON writer for the atlas. Output schema: ushape_evolution_v1.
#
# Two surfaces:
#   write_ushape_json(out_path, candidates_list, params, score_defs, ...)
#       Used by STEP_U06_export_ushape_json.R for offline batch.
#   ushape_candidate_block(window_dt, raw_row, scores, classification)
#       Used by the popstats-server endpoint to return ONE candidate's
#       payload inline (no file write). The atlas merges it with any
#       cached layer it already has client-side.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})


SCHEMA_VERSION <- "ushape_evolution_v1"


# Score equations as plain text — embedded in every JSON so the atlas can show
# them in tooltips without an out-of-band lookup.
SCORE_DEFS <- list(
  dxy_inside_flank_ratio   = "mean(dXY in inside windows) / mean(dXY in flank windows)",
  dxy_u_score              = "mean(dXY at left+right edge) / mean(dXY at center)",
  dxy_asymmetry_score      = "mean(dXY left edge) / mean(dXY right edge)",
  abs_log2_asymmetry       = "abs(log2(dxy_asymmetry_score))",
  dxy_internal_peak_score  = "max(dXY in center windows) / mean(dXY at edge windows)",
  dxy_flatness_score       = "sd(dXY inside) / mean(dXY inside)",
  fst_inside_flank_ratio   = "mean(FST in inside) / mean(FST in flanks)",
  fst_u_score              = "mean(FST at edges) / mean(FST at center)",
  fst_internal_peak_score  = "max(FST center) / mean(FST edges)",
  pi_ratio_homo1_homo2     = "pi_HOMO1_inside / pi_HOMO2_inside",
  pi_imbalance_score       = "abs(log2(pi_ratio_homo1_homo2))",
  arrangement_contrast_score = "mean of capped(dxy_inside_flank, fst_inside_flank, 5*afd_inside)",
  oldness_score            = paste(
    "mean of capped(dxy_inside_flank, dxy_inside/pi_avg, fst_inside_flank,",
    "private_homo1+private_homo2 / shared_poly)"),
  youngness_score          = "max(0, 2 - oldness_score)",
  neutral_u_shape_score    = "mean(dxy_u_score, fst_u_score, -dxy_internal_peak_score)",
  local_adaptation_internal_score = "mean(dxy_internal_peak_score, fst_internal_peak_score)",
  breakpoint_adaptation_score = paste(
    "mean of (mean(edge dXY)/flank dXY, mean(edge FST)/flank FST, dxy_u_score)")
)


# ---- a tiny helper: drop NaN/Inf so jsonlite doesn't emit them ---------------
.jsafe <- function(x) {
  if (is.list(x)) return(lapply(x, .jsafe))
  if (is.atomic(x)) {
    x[!is.finite(x)] <- NA
    return(x)
  }
  x
}


ushape_window_block <- function(window_dt) {
  if (nrow(window_dt) == 0L) return(list())
  apply(window_dt, 1L, function(r) {
    list(
      window_id = as.integer(r["window_id"]),
      start_bp  = as.integer(r["start_bp"]),
      end_bp    = as.integer(r["end_bp"]),
      relative_position = as.numeric(r["relative_position"]),
      zone      = as.character(r["zone"]),
      n_snps    = as.integer(r["n_snps"]),
      pi_homo1  = as.numeric(r["pi_homo1"]),
      pi_homo2  = as.numeric(r["pi_homo2"]),
      dxy_homo1_homo2 = as.numeric(r["dxy_homo1_homo2"]),
      fst_homo1_homo2 = as.numeric(r["fst_homo1_homo2"]),
      allele_freq_delta = as.numeric(r["allele_freq_delta"])
    )
  })
}


ushape_candidate_block <- function(window_dt, raw_row, scores,
                                   classification, groups,
                                   matched_bg = NULL) {
  list(
    candidate_id = raw_row$candidate_id,
    chrom        = raw_row$chrom,
    start_bp     = as.integer(raw_row$start_bp),
    end_bp       = as.integer(raw_row$end_bp),
    length_bp    = as.integer(raw_row$length_bp),

    groups_used = list(
      homo1_label = "HOMO_1", homo2_label = "HOMO_2", het_label = "HET",
      n_homo1 = as.integer(raw_row$n_homo1),
      n_homo2 = as.integer(raw_row$n_homo2),
      n_het   = as.integer(raw_row$n_het)
    ),

    raw_summary = .jsafe(list(
      dxy_inside_mean      = raw_row$dxy_inside_mean,
      dxy_flank_mean       = raw_row$dxy_flank_mean,
      dxy_left_edge_mean   = raw_row$dxy_left_edge_mean,
      dxy_right_edge_mean  = raw_row$dxy_right_edge_mean,
      dxy_center_mean      = raw_row$dxy_center_mean,
      fst_inside_mean      = raw_row$fst_inside_mean,
      fst_flank_mean       = raw_row$fst_flank_mean,
      fst_left_edge_mean   = raw_row$fst_left_edge_mean,
      fst_right_edge_mean  = raw_row$fst_right_edge_mean,
      fst_center_mean      = raw_row$fst_center_mean,
      pi_homo1_inside_mean = raw_row$pi_homo1_inside_mean,
      pi_homo2_inside_mean = raw_row$pi_homo2_inside_mean,
      pi_homo1_flank_mean  = raw_row$pi_homo1_flank_mean,
      pi_homo2_flank_mean  = raw_row$pi_homo2_flank_mean,
      allele_freq_delta_inside_mean = raw_row$allele_freq_delta_inside_mean,
      private_snp_homo1_count = as.integer(raw_row$private_snp_homo1_count %||% NA),
      private_snp_homo2_count = as.integer(raw_row$private_snp_homo2_count %||% NA),
      fixed_diff_count        = as.integer(raw_row$fixed_diff_count %||% NA),
      shared_snp_count        = as.integer(raw_row$shared_snp_count %||% NA)
    )),

    shape_scores = .jsafe(as.list(scores[, !"candidate_id", with = FALSE])),

    matched_background = .jsafe(if (is.null(matched_bg) || nrow(matched_bg) == 0L) {
      list(bg_mode = "absent")
    } else {
      as.list(matched_bg[1L])
    }),

    classification = list(
      primary_class      = classification$primary_class,
      secondary_class    = classification$secondary_class,
      confidence         = classification$confidence,
      reason             = classification$reason,
      flags              = classification$flags,
      rule_based_class   = classification$primary_class,
      unsupervised_cluster = classification$unsupervised_cluster %||% NA_character_,
      cluster_label      = classification$cluster_label %||% NA_character_
    ),

    windows = ushape_window_block(window_dt)
  )
}


write_ushape_json <- function(out_path, candidate_blocks, params,
                              created_by = "STEP_U06_export_ushape_json.R") {
  payload <- list(
    format_version  = SCHEMA_VERSION,
    analysis_name   = "inversion_ushape_evolution",
    created_by      = created_by,
    candidate_count = length(candidate_blocks),
    window_definition = list(
      n_windows_inside = params$n_windows_inside,
      edge_fraction    = params$edge_fraction,
      center_fraction  = params$center_fraction,
      max_flank_bp     = params$max_flank_bp,
      min_flank_bp     = params$min_flank_bp,
      min_window_bp    = params$min_window_bp
    ),
    score_definitions = SCORE_DEFS,
    thresholds = list(
      u_score_high            = params$u_score_high,
      inside_flank_high       = params$inside_flank_high,
      internal_peak_high      = params$internal_peak_high,
      asymmetry_log2_high     = params$asymmetry_log2_high,
      fst_enrichment_high     = params$fst_enrichment_high,
      oldness_min             = params$oldness_min,
      dxy_min_inside          = params$dxy_min_inside,
      flatness_max_for_flat   = params$flatness_max_for_flat
    ),
    candidates = candidate_blocks
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write(jsonlite::toJSON(payload, auto_unbox = TRUE, na = "null",
                         pretty = TRUE, digits = 6), out_path)
  invisible(out_path)
}


`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
