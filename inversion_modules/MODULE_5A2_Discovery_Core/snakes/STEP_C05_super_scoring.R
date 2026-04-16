#!/usr/bin/env Rscript

# =============================================================================
# STEP_C05_super_scoring.R  (v1.0)
#
# THREE-LEVEL INVERSION SCORING SYSTEM
#
# Consumes outputs from all upstream modules and produces:
#   TABLE A: inversion_master.tsv.gz    — one row per candidate inversion
#   TABLE B: breakpoint_master.tsv.gz   — one row per breakpoint
#   TABLE C: inversion_segments.tsv.gz  — one row per internal segment
#
# INPUTS CONSUMED:
#   - snake_candidate_regions.tsv.gz        (from C01b_2_merge)
#   - snake_merge_scores.tsv.gz             (from C01b_2_merge)
#   - snake_merge_regions.tsv.gz            (from C01b_2_merge)
#   - triangle_intervals.tsv.gz             (from C01c_triangle_regimes)
#   - triangle_offdiag_linkage.tsv.gz       (from C01c)
#   - triangle_subregimes.tsv.gz            (from C01c)
#   - snake3v5_window_track.tsv.gz          (from C04_snake3_ghsl)
#   - snake3v5_karyotype_calls.tsv.gz       (from C04_snake3_ghsl)
#   - snake_inv_likeness.tsv.gz             (from C01a_precompute)
#   - precomp/*.precomp.rds                 (from C01a_precompute)
#   - boundary_catalog_*.tsv.gz             (from landscape, if available)
#
# SCORING PHILOSOPHY:
#   1. Raw variables first (measurements)
#   2. Interpreted scores [0-1] (per evidence channel)
#   3. Category scores (structural / breakpoint / population / complexity)
#   4. Final combined score + decision label
#   Never combine too early — preserve interpretability.
#
# Usage:
#   Rscript STEP_C05_super_scoring.R <pipeline_outdir> <scoring_outdir> \
#     [--chrom C_gar_LG01] [--repeat_bed <bed>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# ARGS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript STEP_C05_super_scoring.R <pipeline_outdir> <scoring_outdir> [opts]")

pipeline_dir <- args[1]
scoring_dir  <- args[2]
chrom_filter <- NULL
repeat_bed   <- NULL

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args))       { chrom_filter <- args[i+1]; i <- i+2L }
  else if (a == "--repeat_bed" && i < length(args)) { repeat_bed <- args[i+1]; i <- i+2L }
  else { i <- i+1L }
}

dir.create(scoring_dir, recursive = TRUE, showWarnings = FALSE)

message("================================================================")
message("[SCORE] Super Scoring System v1.0")
message("================================================================")
message("[SCORE] Pipeline dir: ", pipeline_dir)
message("[SCORE] Output dir:   ", scoring_dir)

# =============================================================================
# HELPER: safe file loader
# =============================================================================

safe_load <- function(path, label) {
  if (file.exists(path)) {
    dt <- tryCatch(fread(path), error = function(e) { message("[SCORE] WARN: ", label, " unreadable"); data.table() })
    message("[SCORE] Loaded ", label, ": ", nrow(dt), " rows")
    dt
  } else {
    message("[SCORE] SKIP ", label, ": not found (", path, ")")
    data.table()
  }
}

# =============================================================================
# LOAD ALL UPSTREAM OUTPUTS
# =============================================================================

# Snake 1 merge candidates — the primary candidate list
cand_dt  <- safe_load(file.path(pipeline_dir, "snake_candidate_regions.tsv.gz"), "candidates")
merge_dt <- safe_load(file.path(pipeline_dir, "snake_merge_regions.tsv.gz"), "merge_regions")
merge_scores_dt <- safe_load(file.path(pipeline_dir, "snake_merge_scores.tsv.gz"), "merge_scores")

# Triangle regimes
tri_dt   <- safe_load(file.path(pipeline_dir, "triangle_intervals.tsv.gz"), "triangles")
tri_od   <- safe_load(file.path(pipeline_dir, "triangle_offdiag_linkage.tsv.gz"), "offdiag_linkage")
tri_sub  <- safe_load(file.path(pipeline_dir, "triangle_subregimes.tsv.gz"), "subregimes")

# Snake 3 GHSL
ghsl_dt  <- safe_load(file.path(pipeline_dir, "snake3v5_window_track.tsv.gz"), "ghsl_track")
karyo_dt <- safe_load(file.path(pipeline_dir, "snake3v5_karyotype_calls.tsv.gz"), "karyotype_calls")

# Inv-likeness (from precompute)
invlike_dt <- safe_load(file.path(pipeline_dir, "snake_inv_likeness.tsv.gz"), "inv_likeness")

# Precomp RDS (for sim_mat, boundaries)
precomp_dir <- file.path(pipeline_dir, "precomp")
precomp_list <- list()
if (dir.exists(precomp_dir)) {
  rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
  for (f in rds_files) {
    obj <- readRDS(f)
    precomp_list[[obj$chrom]] <- obj
  }
  message("[SCORE] Precomp RDS: ", length(precomp_list), " chromosomes")
}

# Repeat annotations (optional)
repeat_dt <- data.table()
if (!is.null(repeat_bed) && file.exists(repeat_bed)) {
  repeat_dt <- tryCatch(fread(repeat_bed, header = FALSE,
    col.names = c("chrom", "start", "end", "repeat_class", "repeat_family", "strand")[
      seq_len(min(6, ncol(fread(repeat_bed, nrows = 1, header = FALSE))))]),
    error = function(e) data.table())
  if (nrow(repeat_dt) > 0) message("[SCORE] Repeats: ", nrow(repeat_dt), " annotations")
}

# Boundary catalog (from landscape, if available)
bound_files <- list.files(pipeline_dir, pattern = "boundary_catalog_.*\\.tsv\\.gz$", full.names = TRUE)
bound_dt <- if (length(bound_files) > 0) {
  rbindlist(lapply(bound_files, function(f) tryCatch(fread(f), error = function(e) data.table())), fill = TRUE)
} else data.table()
if (nrow(bound_dt) > 0) message("[SCORE] Boundary catalog: ", nrow(bound_dt), " entries")

# =============================================================================
# CANDIDATE FILTERING
# =============================================================================

if (nrow(cand_dt) == 0) {
  message("[SCORE] No candidates found — writing empty tables")
  fwrite(data.table(), file.path(scoring_dir, "inversion_master.tsv.gz"), sep = "\t")
  fwrite(data.table(), file.path(scoring_dir, "breakpoint_master.tsv.gz"), sep = "\t")
  fwrite(data.table(), file.path(scoring_dir, "inversion_segments.tsv.gz"), sep = "\t")
  message("[DONE]")
  quit(save = "no")
}

if (!is.null(chrom_filter)) cand_dt <- cand_dt[chrom == chrom_filter]
message("[SCORE] Scoring ", nrow(cand_dt), " candidate inversions")

# =============================================================================
# TABLE A: INVERSION MASTER
# =============================================================================

score_one_inversion <- function(ri, cand) {
  chr     <- cand$chrom
  inv_s   <- cand$start_bp
  inv_e   <- cand$end_bp
  span_bp <- inv_e - inv_s
  span_mb <- round(span_bp / 1e6, 3)

  # ── RAW SIGNALS ──

  # 1. Snake 1 merge support
  merge_nwin   <- cand$n_windows %||% NA_integer_
  merge_density <- cand$density %||% NA_real_
  merge_inv_mean <- cand$mean_inv %||% NA_real_
  merge_coherence <- cand$coherence %||% NA_real_

  # 2. Inv-likeness within region
  il_region <- invlike_dt[chrom == chr & start_bp >= inv_s & end_bp <= inv_e]
  il_median <- if (nrow(il_region) > 0) median(il_region$inv_likeness, na.rm = TRUE) else NA_real_
  il_q90    <- if (nrow(il_region) > 0) quantile(il_region$inv_likeness, 0.90, na.rm = TRUE) else NA_real_
  struct_median <- if (nrow(il_region) > 0 && "structure_likeness" %in% names(il_region))
    median(il_region$structure_likeness, na.rm = TRUE) else NA_real_
  family_median <- if (nrow(il_region) > 0 && "family_likeness" %in% names(il_region))
    median(il_region$family_likeness, na.rm = TRUE) else NA_real_

  # 3. Triangle support
  tri_overlap <- tri_dt[chrom == chr &
    !is.na(start_mb) & !is.na(end_mb) &
    start_mb * 1e6 < inv_e & end_mb * 1e6 > inv_s]
  tri_best_type <- if (nrow(tri_overlap) > 0) tri_overlap$type[1] else "none"
  tri_contrast <- if (nrow(tri_overlap) > 0) max(tri_overlap$contrast, na.rm = TRUE) else 0
  tri_sharpness <- if (nrow(tri_overlap) > 0 && "sharpness" %in% names(tri_overlap))
    max(tri_overlap$sharpness, na.rm = TRUE) else 0
  tri_squareness <- if (nrow(tri_overlap) > 0 && "squareness" %in% names(tri_overlap))
    max(tri_overlap$squareness, na.rm = TRUE) else 0

  # 4. Off-diagonal linkage (distant windows sharing structure)
  od_count <- 0L
  if (nrow(tri_od) > 0 && all(c("chrom", "region_a_start_mb") %in% names(tri_od))) {
    od_hits <- tri_od[chrom == chr &
      ((region_a_start_mb * 1e6 >= inv_s & region_a_start_mb * 1e6 <= inv_e) |
       (region_b_start_mb * 1e6 >= inv_s & region_b_start_mb * 1e6 <= inv_e))]
    od_count <- nrow(od_hits)
  }

  # 5. GHSL (Snake 3) support
  ghsl_region <- ghsl_dt[chrom == chr & start_bp >= inv_s & end_bp <= inv_e]
  ghsl_pass_frac <- if (nrow(ghsl_region) > 0)
    sum(ghsl_region$ghsl_v5_status == "PASS", na.rm = TRUE) / nrow(ghsl_region) else 0
  ghsl_rho_med <- if (nrow(ghsl_region) > 0 && "rank_stability" %in% names(ghsl_region))
    median(ghsl_region$rank_stability, na.rm = TRUE) else NA_real_
  ghsl_contrast_med <- if (nrow(ghsl_region) > 0 && "div_contrast" %in% names(ghsl_region))
    median(ghsl_region$div_contrast, na.rm = TRUE) else NA_real_

  # 6. Karyotype calls within region
  karyo_region <- karyo_dt[chrom == chr & start_bp >= inv_s & end_bp <= inv_e]
  n_karyo_inv_inv <- sum(karyo_region$call == "INV_INV", na.rm = TRUE)
  n_karyo_inv_non <- sum(karyo_region$call == "INV_nonINV", na.rm = TRUE)
  n_karyo_samples <- length(unique(karyo_region$sample_id))

  # 7. Sim_mat boundary hardness (from precomp)
  boundary_hardness_left  <- NA_real_
  boundary_hardness_right <- NA_real_
  center_vs_flank <- NA_real_
  pc <- precomp_list[[chr]]
  if (!is.null(pc)) {
    dt_pc <- pc$dt
    sim_mat <- pc$sim_mat
    n_w <- nrow(sim_mat)

    # Map bp to window indices
    wi_start <- which.min(abs(dt_pc$start_bp - inv_s))
    wi_end   <- which.min(abs(dt_pc$end_bp - inv_e))
    if (wi_start > 0 && wi_end > 0 && wi_end > wi_start) {
      # Center similarity
      center_idx <- wi_start:wi_end
      if (length(center_idx) >= 4) {
        center_sim <- mean(sim_mat[center_idx, center_idx], na.rm = TRUE)

        # Left flank similarity
        flank_w <- min(20L, length(center_idx) %/% 2)
        left_flank <- max(1, wi_start - flank_w):(wi_start - 1)
        right_flank <- (wi_end + 1):min(n_w, wi_end + flank_w)

        flank_sim <- NA_real_
        if (length(left_flank) >= 2 && length(right_flank) >= 2) {
          flank_sim <- mean(c(
            mean(sim_mat[left_flank, left_flank], na.rm = TRUE),
            mean(sim_mat[right_flank, right_flank], na.rm = TRUE)
          ))
        }
        center_vs_flank <- if (is.finite(flank_sim) && flank_sim > 0) center_sim / flank_sim else NA_real_

        # Boundary hardness: cross-boundary sim drop
        if (length(left_flank) >= 2) {
          cross_left <- mean(sim_mat[left_flank, center_idx[1:min(5, length(center_idx))]], na.rm = TRUE)
          boundary_hardness_left <- center_sim - cross_left
        }
        if (length(right_flank) >= 2) {
          cross_right <- mean(sim_mat[center_idx[max(1, length(center_idx)-4):length(center_idx)], right_flank], na.rm = TRUE)
          boundary_hardness_right <- center_sim - cross_right
        }
      }
    }
  }

  # 8. Repeat hotspot (if available)
  repeat_overlap_bp <- 0L
  repeat_classes <- ""
  if (nrow(repeat_dt) > 0) {
    rep_hits <- repeat_dt[chrom == chr & start < inv_e & end > inv_s]
    if (nrow(rep_hits) > 0) {
      repeat_overlap_bp <- sum(pmin(rep_hits$end, inv_e) - pmax(rep_hits$start, inv_s))
      if ("repeat_class" %in% names(rep_hits))
        repeat_classes <- paste(unique(rep_hits$repeat_class), collapse = ";")
    }
  }

  # ── INTERPRETED SCORES [0-1] ──

  # Structural support: inv_likeness + triangle + merge coherence
  s_invlike <- if (is.finite(il_median)) pmin(1, pmax(0, il_median / 0.7)) else 0
  s_triangle <- pmin(1, pmax(0, tri_contrast / 0.08))
  s_merge <- if (is.finite(merge_coherence)) pmin(1, merge_coherence) else 0
  structural_support <- round(0.40 * s_invlike + 0.35 * s_triangle + 0.25 * s_merge, 4)

  # Breakpoint support: boundary hardness + sharpness + squareness
  s_hard_l <- if (is.finite(boundary_hardness_left)) pmin(1, pmax(0, boundary_hardness_left / 0.15)) else 0
  s_hard_r <- if (is.finite(boundary_hardness_right)) pmin(1, pmax(0, boundary_hardness_right / 0.15)) else 0
  s_sharp <- pmin(1, pmax(0, tri_sharpness / 0.05))
  s_square <- pmin(1, tri_squareness)
  breakpoint_support <- round(0.30 * (s_hard_l + s_hard_r) / 2 + 0.35 * s_sharp + 0.35 * s_square, 4)

  # Population support: GHSL rank stability + karyotype calls
  s_ghsl <- if (is.finite(ghsl_rho_med)) pmin(1, pmax(0, ghsl_rho_med)) else 0
  s_karyo <- if (n_karyo_samples > 0) pmin(1, n_karyo_samples / 50) else 0
  s_ghsl_pass <- pmin(1, ghsl_pass_frac / 0.3)
  population_support <- round(0.50 * s_ghsl + 0.25 * s_karyo + 0.25 * s_ghsl_pass, 4)

  # Repeat hotspot
  repeat_hotspot <- round(pmin(1, repeat_overlap_bp / span_bp), 4)

  # Complexity: family contamination + nestedness hints
  s_family <- if (is.finite(family_median)) pmin(1, family_median) else 0
  complexity_score <- round(s_family, 4)

  # ── OVERALL CONFIDENCE ──
  overall_confidence <- round(
    0.35 * structural_support +
    0.25 * breakpoint_support +
    0.25 * population_support +
    0.15 * (1 - complexity_score),  # low family = good
    4)

  # ── MECHANISM CLASS ──
  mechanism_class <- "unknown"
  if (structural_support > 0.6 && breakpoint_support > 0.5 && s_square > 0.4) {
    mechanism_class <- "strong_simple_inversion"
  } else if (structural_support > 0.4 && breakpoint_support > 0.3) {
    mechanism_class <- "likely_real_moderate"
  } else if (s_family > 0.5 && structural_support > 0.3) {
    mechanism_class <- "family_mimic_suspect"
  } else if (repeat_hotspot > 0.3 && breakpoint_support > 0.3) {
    mechanism_class <- "repeat_mediated_suspect"
  } else if (structural_support > 0.3 && breakpoint_support < 0.2) {
    mechanism_class <- "weak_or_eroded"
  } else if (overall_confidence < 0.2) {
    mechanism_class <- "probable_artifact"
  }

  # ── FINAL LABEL ──
  final_label <- if (overall_confidence >= 0.6) "HIGH_CONFIDENCE"
    else if (overall_confidence >= 0.35) "MODERATE"
    else if (overall_confidence >= 0.15) "WEAK"
    else "ARTIFACT"

  data.table(
    inv_id = ri,
    chrom = chr,
    start_bp = inv_s, end_bp = inv_e, span_bp = span_bp, span_mb = span_mb,
    # Candidate source
    source = cand$source %||% "snake",
    n_windows = merge_nwin,
    merge_density = round(merge_density %||% NA_real_, 4),
    merge_coherence = round(merge_coherence %||% NA_real_, 4),
    # Raw signals
    invlike_median = round(il_median %||% NA_real_, 4),
    invlike_q90 = round(il_q90 %||% NA_real_, 4),
    structure_median = round(struct_median %||% NA_real_, 4),
    family_median = round(family_median %||% NA_real_, 4),
    tri_type = tri_best_type,
    tri_contrast = round(tri_contrast, 4),
    tri_sharpness = round(tri_sharpness, 4),
    tri_squareness = round(tri_squareness, 4),
    offdiag_links = od_count,
    ghsl_pass_frac = round(ghsl_pass_frac, 4),
    ghsl_rho_median = round(ghsl_rho_med %||% NA_real_, 4),
    ghsl_contrast_median = round(ghsl_contrast_med %||% NA_real_, 4),
    n_karyo_inv_inv = n_karyo_inv_inv,
    n_karyo_inv_non = n_karyo_inv_non,
    n_karyo_samples = n_karyo_samples,
    boundary_hardness_left = round(boundary_hardness_left %||% NA_real_, 4),
    boundary_hardness_right = round(boundary_hardness_right %||% NA_real_, 4),
    center_vs_flank = round(center_vs_flank %||% NA_real_, 4),
    repeat_overlap_bp = repeat_overlap_bp,
    repeat_classes = repeat_classes,
    # Interpreted scores [0-1]
    structural_support = structural_support,
    breakpoint_support = breakpoint_support,
    population_support = population_support,
    repeat_hotspot_score = repeat_hotspot,
    complexity_score = complexity_score,
    # Final
    overall_confidence = overall_confidence,
    mechanism_class = mechanism_class,
    final_label = final_label
  )
}

message("[SCORE] Building TABLE A: inversion_master...")
inv_rows <- lapply(seq_len(nrow(cand_dt)), function(ri) {
  tryCatch(score_one_inversion(ri, cand_dt[ri]),
           error = function(e) { message("[SCORE] WARN inv ", ri, ": ", e$message); NULL })
})
inv_master <- rbindlist(Filter(Negate(is.null), inv_rows), fill = TRUE)
message("[SCORE] TABLE A: ", nrow(inv_master), " inversions scored")

# =============================================================================
# TABLE B: BREAKPOINT MASTER
# =============================================================================

message("[SCORE] Building TABLE B: breakpoint_master...")
bp_rows <- list()
bp_id <- 0L

for (ri in seq_len(nrow(inv_master))) {
  inv <- inv_master[ri]
  chr <- inv$chrom
  pc <- precomp_list[[chr]]

  for (side in c("left", "right")) {
    bp_id <- bp_id + 1L
    bp_pos <- if (side == "left") inv$start_bp else inv$end_bp
    hardness <- if (side == "left") inv$boundary_hardness_left else inv$boundary_hardness_right

    # Boundary type from landscape catalog
    bp_class <- "unknown"
    if (nrow(bound_dt) > 0 && "chrom" %in% names(bound_dt)) {
      near <- bound_dt[chrom == chr & abs(start - bp_pos) < 50000]
      if (nrow(near) > 0 && "boundary_type" %in% names(near))
        bp_class <- near$boundary_type[1]
    }

    # Repeat overlap at breakpoint (±5kb)
    bp_repeat_class <- ""
    bp_repeat_bp <- 0L
    if (nrow(repeat_dt) > 0) {
      rep_near <- repeat_dt[chrom == chr & start < bp_pos + 5000 & end > bp_pos - 5000]
      if (nrow(rep_near) > 0) {
        bp_repeat_bp <- sum(pmin(rep_near$end, bp_pos + 5000) - pmax(rep_near$start, bp_pos - 5000))
        if ("repeat_class" %in% names(rep_near))
          bp_repeat_class <- paste(unique(rep_near$repeat_class), collapse = ";")
      }
    }

    # Sim_mat drop at boundary
    bp_sim_drop <- NA_real_
    if (!is.null(pc)) {
      dt_pc <- pc$dt
      sim_mat <- pc$sim_mat
      wi_bp <- which.min(abs(dt_pc$start_bp - bp_pos))
      if (wi_bp > 5 && wi_bp < nrow(sim_mat) - 5) {
        inside <- if (side == "left") sim_mat[wi_bp, (wi_bp+1):(wi_bp+5)]
                  else sim_mat[wi_bp, (wi_bp-5):(wi_bp-1)]
        outside <- if (side == "left") sim_mat[wi_bp, (wi_bp-5):(wi_bp-1)]
                   else sim_mat[wi_bp, (wi_bp+1):(wi_bp+5)]
        bp_sim_drop <- mean(inside, na.rm = TRUE) - mean(outside, na.rm = TRUE)
      }
    }

    # Scores
    bp_cleanliness <- pmin(1, pmax(0, (hardness %||% 0) / 0.15))
    bp_mechanism_signal <- if (bp_repeat_bp > 500) 0.8 else if (bp_repeat_bp > 100) 0.4 else 0
    bp_confidence <- round(0.6 * bp_cleanliness + 0.2 * pmin(1, pmax(0, (bp_sim_drop %||% 0) / 0.1)) + 0.2 * (1 - bp_mechanism_signal * 0.3), 4)

    bp_rows[[bp_id]] <- data.table(
      bp_id = bp_id,
      inv_id = inv$inv_id,
      chrom = chr,
      pos_bp = bp_pos,
      pos_mb = round(bp_pos / 1e6, 3),
      bp_side = side,
      bp_class = bp_class,
      hardness_score = round(hardness %||% NA_real_, 4),
      sim_drop = round(bp_sim_drop %||% NA_real_, 4),
      repeat_class = bp_repeat_class,
      repeat_overlap_bp = bp_repeat_bp,
      bp_cleanliness_score = round(bp_cleanliness, 4),
      bp_mechanism_signal_score = round(bp_mechanism_signal, 4),
      bp_confidence_score = bp_confidence
    )
  }
}

bp_master <- if (length(bp_rows) > 0) rbindlist(bp_rows, fill = TRUE) else data.table()
message("[SCORE] TABLE B: ", nrow(bp_master), " breakpoints scored")

# Link breakpoint IDs back to inversion master
if (nrow(bp_master) > 0 && nrow(inv_master) > 0) {
  for (ri in seq_len(nrow(inv_master))) {
    inv_id_val <- inv_master$inv_id[ri]
    left_bp <- bp_master[inv_id == inv_id_val & bp_side == "left"]
    right_bp <- bp_master[inv_id == inv_id_val & bp_side == "right"]
    if (nrow(left_bp) > 0) inv_master[ri, left_bp_id := left_bp$bp_id[1]]
    if (nrow(right_bp) > 0) inv_master[ri, right_bp_id := right_bp$bp_id[1]]
  }
}

# =============================================================================
# TABLE C: INVERSION SEGMENTS
# =============================================================================

message("[SCORE] Building TABLE C: inversion_segments...")
seg_rows <- list()
seg_id <- 0L

for (ri in seq_len(nrow(inv_master))) {
  inv <- inv_master[ri]
  chr <- inv$chrom
  inv_s <- inv$start_bp; inv_e <- inv$end_bp
  span <- inv_e - inv_s

  # Define 5 segments: outside_left, inside_left, center, inside_right, outside_right
  flank_bp <- min(span * 0.15, 500000)  # 15% of span or 500kb max
  edge_bp  <- min(span * 0.10, 200000)  # 10% of span or 200kb max

  segments <- list(
    list(label = "outside_left",  s = inv_s - flank_bp, e = inv_s),
    list(label = "inside_left",   s = inv_s, e = inv_s + edge_bp),
    list(label = "center",        s = inv_s + edge_bp, e = inv_e - edge_bp),
    list(label = "inside_right",  s = inv_e - edge_bp, e = inv_e),
    list(label = "outside_right", s = inv_e, e = inv_e + flank_bp)
  )

  for (seg in segments) {
    seg_id <- seg_id + 1L

    # Inv-likeness in segment
    il_seg <- invlike_dt[chrom == chr & start_bp >= seg$s & end_bp <= seg$e]
    il_seg_med <- if (nrow(il_seg) > 0) median(il_seg$inv_likeness, na.rm = TRUE) else NA_real_

    # Het signal in segment
    het_seg <- NA_real_
    if (nrow(il_seg) > 0 && "dosage_het_rate_cv" %in% names(il_seg))
      het_seg <- median(il_seg$dosage_het_rate_cv, na.rm = TRUE)

    # GHSL in segment
    ghsl_seg <- ghsl_dt[chrom == chr & start_bp >= seg$s & end_bp <= seg$e]
    ghsl_rho_seg <- if (nrow(ghsl_seg) > 0 && "rank_stability" %in% names(ghsl_seg))
      median(ghsl_seg$rank_stability, na.rm = TRUE) else NA_real_

    # Repeat density
    rep_seg <- 0
    if (nrow(repeat_dt) > 0) {
      rh <- repeat_dt[chrom == chr & start < seg$e & end > seg$s]
      if (nrow(rh) > 0) rep_seg <- sum(pmin(rh$end, seg$e) - pmax(rh$start, seg$s)) / max(1, seg$e - seg$s)
    }

    # Subregime patchiness (from triangle)
    patch_score <- 0
    if (nrow(tri_sub) > 0 && all(c("chrom", "start_mb", "end_mb") %in% names(tri_sub))) {
      sub_hits <- tri_sub[chrom == chr & start_mb * 1e6 < seg$e & end_mb * 1e6 > seg$s]
      if (nrow(sub_hits) > 0) {
        n_hot <- sum(sub_hits$sub_type == "hot", na.rm = TRUE)
        patch_score <- n_hot / nrow(sub_hits)
      }
    }

    seg_rows[[seg_id]] <- data.table(
      inv_id = inv$inv_id,
      segment_id = seg_id,
      segment_class = seg$label,
      seg_start = as.integer(seg$s),
      seg_end = as.integer(seg$e),
      seg_span_bp = as.integer(seg$e - seg$s),
      invlike_median = round(il_seg_med %||% NA_real_, 4),
      het_cv_median = round(het_seg %||% NA_real_, 4),
      ghsl_rho_median = round(ghsl_rho_seg %||% NA_real_, 4),
      repeat_density = round(rep_seg, 4),
      patchiness_score = round(patch_score, 4),
      n_windows = nrow(il_seg)
    )
  }
}

seg_master <- if (length(seg_rows) > 0) rbindlist(seg_rows, fill = TRUE) else data.table()
message("[SCORE] TABLE C: ", nrow(seg_master), " segments scored")

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

fwrite(inv_master, file.path(scoring_dir, "inversion_master.tsv.gz"), sep = "\t")
fwrite(bp_master, file.path(scoring_dir, "breakpoint_master.tsv.gz"), sep = "\t")
fwrite(seg_master, file.path(scoring_dir, "inversion_segments.tsv.gz"), sep = "\t")

# Summary
message("\n================================================================")
message("[DONE] Super Scoring System v1.0")
message("================================================================")
message("  TABLE A: ", file.path(scoring_dir, "inversion_master.tsv.gz"),
        " (", nrow(inv_master), " inversions)")
message("  TABLE B: ", file.path(scoring_dir, "breakpoint_master.tsv.gz"),
        " (", nrow(bp_master), " breakpoints)")
message("  TABLE C: ", file.path(scoring_dir, "inversion_segments.tsv.gz"),
        " (", nrow(seg_master), " segments)")

if (nrow(inv_master) > 0) {
  message("\n  Confidence distribution:")
  for (lab in c("HIGH_CONFIDENCE", "MODERATE", "WEAK", "ARTIFACT")) {
    n_lab <- sum(inv_master$final_label == lab)
    message("    ", lab, ": ", n_lab)
  }
  message("\n  Mechanism classes:")
  for (mc in unique(inv_master$mechanism_class)) {
    message("    ", mc, ": ", sum(inv_master$mechanism_class == mc))
  }
}
