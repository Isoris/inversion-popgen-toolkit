#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01d_candidate_scoring.R  (v9.3.4 — Engine B popgen stream + Cheat 25)
#
# TIERED EVIDENCE FRAMEWORK for inversion candidate scoring.
#
# v9.3.4 changes:
#   - Sources load_bridge.R for Engine B access.
#   - NEW: Popgen annotation stream via get_region_stats():
#     Fst(band1,band3), theta_pi, theta_W, Tajima_D per candidate.
#     These are NOT scoring dimensions — they are annotation columns
#     for downstream hypothesis tests and manuscript tables.
#     The C binary computes all from BEAGLE dosage in ~2s per call.
#   - Cheat 5 family_fst_ratio from precomp (if available from v9.3.4 precomp).
#
# v9.3.2 changes:
#   - D11 boundary concordance: reads C01g boundary catalog, scores each
#     candidate by the concordance of its left + right boundaries.
#   - D12 snake concordance: reads C01b core_regions, checks overlap with
#     staircase blocks. Blocks without snake core → capped Tier 4.
#   - Cheat 25 block viability: 4-test battery using SV + boundary + peel.
#     DEAD candidates → forced Tier 4.
#   - New args: --cores_dir, --boundary_dir, --flashlight_dir
#
# v9.3: Rebuilt to consume inv_detect_v9.3 scoring tables DIRECTLY.
# The old C01c triangle/tube/GHSL pipeline is superseded.
#
# OLD (v8.4):  C01c triangles → bridge → C01d (10 dimensions from triangles)
# NEW (v9.3):  inv_detect scoring_table_*.tsv → C01d (12 dimensions from detector)
#                                                 [chat-13 Finding AO: was
#                                                  documented as 10; actual
#                                                  count is 12 per the mapping
#                                                  below, D1–D12.]
#
# Dimension mapping (12 dimensions after v9.3.2):
#   D1  Block strength    — contrast × squareness × sharpness (from 08_bloc_scoring)
#   D2  Block shape       — shape_class + occupancy + homogeneity (from 08 + 04)
#   D3  NN persistence    — survives_nn40/nn80 + nn_birth (from 02 + 09)
#   D4  Decay flatness    — far_near_ratio + flatness + monotonicity (from 03)
#   D5  Interior quality  — homogeneity + low patchiness + low stripe_count (from 04)
#   D6  Consensus support — n_variants found in + confidence label (from 10)
#   D7  SV breakpoint     — sv_overlap_pct + sv_left/right_dist (from 06)
#   D8  Peel diagnostic   — L1b + L2 effect classes (from C01n peeling)
#   D9  PCA/ICA clusters  — pc1_silhouette + ica_kurtosis (from 08_pca)
#   D10 Partition (GHSL)  — partition_stability + entropy (from 05)
#   D11 Boundary concord  — C01g boundary catalog (v9.3.2, filled after loop)
#   D12 Snake concordance — C01b core_regions overlap (v9.3.2, filled after loop)
#
# D3 (ancestry diversity / eff_K) is GONE — replaced by D3 (NN persistence)
# which is a genuinely independent structural test.
#
# Tier assignment: count-based matching manuscript.
#   Tier 1: ≥8/12 dimensions positive
#   Tier 2: 6–7 dimensions positive
#   Tier 3: 4–5 dimensions positive
#   Tier 4: <4 dimensions positive
#
# Usage:
#   Rscript STEP_C01d_candidate_scoring.R <detector_dir> <outdir> \
#     [--precomp_dir <dir>] [--hyp_dir <dir>]
#
# Where <detector_dir> contains scoring_table_*.tsv files from inv_detect_v9.3.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ── Source load_bridge.R (provides smap, reg, get_Q, get_region_stats) ───────
.bridge_file <- Sys.getenv("LOAD_BRIDGE", "")
if (!nzchar(.bridge_file)) {
  for (.bp in c("utils/load_bridge.R", "../utils/load_bridge.R",
                 file.path(Sys.getenv("BASE", ""), "inversion_codebase_v8.5/utils/load_bridge.R"))) {
    if (file.exists(.bp)) { .bridge_file <- .bp; break }
  }
}
.bridge_available <- FALSE
if (nzchar(.bridge_file) && file.exists(.bridge_file)) {
  tryCatch({
    source(.bridge_file)
    .bridge_available <- TRUE
    message("[C01d] load_bridge.R sourced — Engine B available for popgen annotations")
  }, error = function(e) {
    message("[C01d] WARNING: load_bridge.R failed: ", conditionMessage(e))
  })
} else {
  message("[C01d] load_bridge.R not found — popgen annotations will be skipped")
}
safe_num <- function(x, default = NA_real_) {
  x <- as.numeric(x)
  fifelse(is.finite(x), x, default)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript STEP_C01d_candidate_scoring.R <detector_dir> <outdir> [opts]")

detector_dir <- args[1]; outdir <- args[2]
precomp_dir <- NULL; hyp_dir <- NULL
cores_dir <- NULL; boundary_dir <- NULL

# BUGFIX 2026-04-17: --flashlight_dir was parsed but never read in this
# script. D7 (sv_breakpoint) already gets its SV info from the staircase
# scoring table's sv_overlap_pct / n_sv_hits columns (populated by
# phase_2/2d/STEP_D06_sv_overlap.R). The flag is accepted for backward
# compatibility with existing launchers but is silently ignored.
#
# 2026-04-17 (chat 5): an earlier attempt (FIX 29 first draft) added a
# --phase3_dir flag + D13 scoring dimension to pull phase_3 OR-test
# results into Layer A's composite. That was the wrong abstraction —
# phase_3's OR test is Layer D of the 4-layer evidence model, not a
# sub-dimension of Layer A. Layer D is written into the evidence
# registry directly by phase_3/STEP03 via the `existence_layer_d`
# block, and consumed by C01f's compute_group_validation() through
# the `q7_layer_d_fisher_p` / `q7_layer_d_fisher_or` flat keys.
# C01d stays Layer-A-only. The --phase3_dir flag was reverted.
i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--precomp_dir" && i < length(args))    { precomp_dir <- args[i+1]; i <- i+2L }
  else if (a == "--hyp_dir" && i < length(args))    { hyp_dir <- args[i+1]; i <- i+2L }
  else if (a == "--cores_dir" && i < length(args))  { cores_dir <- args[i+1]; i <- i+2L }
  else if (a == "--boundary_dir" && i < length(args)) { boundary_dir <- args[i+1]; i <- i+2L }
  else if (a == "--flashlight_dir" && i < length(args)) {
    # accepted but ignored — D7 uses scoring_table SV columns instead
    i <- i+2L
  }
  else { i <- i+1L }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

DPI <- 350
THEME_BASE <- theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        plot.caption = element_text(size = 6, color = "grey60", hjust = 0))

# =============================================================================
# LOAD inv_detect_v9.3 SCORING TABLES
# =============================================================================

score_files <- sort(list.files(detector_dir, pattern = "^scoring_table_.*\\.tsv$",
                                full.names = TRUE))
if (length(score_files) == 0) {
  # Try .tsv.gz
  score_files <- sort(list.files(detector_dir, pattern = "^scoring_table_.*\\.tsv\\.gz$",
                                  full.names = TRUE))
}
if (length(score_files) == 0) stop("No scoring_table_*.tsv files in ", detector_dir)

message("[C01d] Loading ", length(score_files), " scoring tables from: ", detector_dir)
all_tabs <- lapply(score_files, fread, fill = TRUE)
iv_dt <- rbindlist(all_tabs, fill = TRUE)

# Standardize chromosome column
if (!"chrom" %in% names(iv_dt) && "chr" %in% names(iv_dt)) {
  setnames(iv_dt, "chr", "chrom")
}

# Ensure start_mb / end_mb exist
if (!"start_mb" %in% names(iv_dt) && "start" %in% names(iv_dt)) {
  window_bp <- 50000  # default
  iv_dt[, start_mb := round((start - 1) * window_bp / 1e6, 3)]
  iv_dt[, end_mb := round(end * window_bp / 1e6, 3)]
}

# Ensure interval_id exists
if (!"interval_id" %in% names(iv_dt)) {
  if ("block_id" %in% names(iv_dt)) {
    iv_dt[, interval_id := block_id]
  } else if ("candidate_id" %in% names(iv_dt)) {
    iv_dt[, interval_id := candidate_id]
  } else {
    iv_dt[, interval_id := seq_len(.N)]
  }
}

message("[C01d] Total candidates: ", nrow(iv_dt), " across ",
        length(unique(iv_dt$chrom)), " chromosomes")

# Load hypothesis verdicts if available
hyp_verd <- data.table()
if (!is.null(hyp_dir)) {
  f <- file.path(hyp_dir, "hypothesis_verdicts.tsv")
  if (file.exists(f)) { hyp_verd <- fread(f); message("[C01d] Hypothesis verdicts: ", nrow(hyp_verd)) }
}

# =============================================================================
# SCORING DIMENSIONS (10 — all from inv_detect_v9.3 evidence stack)
# =============================================================================

score_candidates <- function(iv_dt, hyp_verd) {
  scores <- list()

  for (ri in seq_len(nrow(iv_dt))) {
    iv <- iv_dt[ri]
    chr <- iv$chrom; iid <- iv$interval_id

    # --- D1: Block strength (contrast + squareness + sharpness) ---
    d1_contrast  <- pmin(1, safe_num(iv$contrast, 0) / 0.08)
    d1_sq        <- pmin(1, safe_num(iv$squareness, 0))
    d1_sharp     <- pmin(1, safe_num(iv$sharpness, 0) / 0.04)
    d1 <- 0.40 * d1_contrast + 0.35 * d1_sq + 0.25 * d1_sharp

    # --- D2: Block shape (shape_class + occupancy + homogeneity) ---
    shape <- iv$shape_class %||% "unknown"
    shape_score <- if (shape == "strong_square") 1.0
                   else if (shape == "diffuse_square") 0.6
                   else if (shape == "diagonal_band") 0.3
                   else if (shape == "ambiguous") 0.15
                   else 0  # noise, unknown
    occ <- pmin(1, safe_num(iv$occupancy, 0))
    homog <- pmin(1, safe_num(iv$homogeneity, 0))
    d2 <- 0.50 * shape_score + 0.25 * occ + 0.25 * homog

    # --- D3: NN persistence (survives nn40/80 + nn_birth) ---
    nn40 <- if ("survives_nn40" %in% names(iv) && !is.na(iv$survives_nn40) &&
                iv$survives_nn40 == TRUE) 1.0 else 0
    nn80 <- if ("survives_nn80" %in% names(iv) && !is.na(iv$survives_nn80) &&
                iv$survives_nn80 == TRUE) 1.0 else 0
    nn_birth <- safe_num(iv$nn_birth, 0)
    nn_birth_score <- pmin(1, nn_birth / 200)  # saturates at nn_birth=200
    d3 <- 0.30 * nn40 + 0.30 * nn80 + 0.40 * nn_birth_score

    # --- D4: Decay flatness (far_near_ratio + flatness) ---
    fnr <- pmin(1, safe_num(iv$far_near_ratio, 0))
    flat <- pmin(1, safe_num(iv$flatness, 0))
    mono <- pmin(1, safe_num(iv$monotonicity, 0))
    d4 <- 0.40 * fnr + 0.35 * flat + 0.25 * mono

    # --- D5: Interior quality (low cv, low patchiness, low stripes) ---
    hom_score <- pmin(1, safe_num(iv$homogeneity, 0))
    patch_score <- pmin(1, pmax(0, 1 - safe_num(iv$patchiness, 0) * 3))
    stripe_n <- safe_num(iv$stripe_count, 0)
    stripe_score <- pmin(1, pmax(0, 1 - stripe_n / 5))
    d5 <- 0.40 * hom_score + 0.35 * patch_score + 0.25 * stripe_score

    # --- D6: Consensus support (how many matrix variants found this block) ---
    n_var <- safe_num(iv$n_variants, 1)
    consensus_label <- iv$consensus_confidence %||% iv$confidence %||% "LOW"
    d6 <- pmin(1, (n_var - 1) / 4)  # saturates at 5 variants

    # --- D7: SV breakpoint evidence ---
    sv_pct <- safe_num(iv$sv_overlap_pct, 0)
    n_sv <- safe_num(iv$n_sv_hits, 0)
    d7 <- pmin(1, sv_pct * 0.6 + pmin(1, n_sv / 3) * 0.4)

    # --- D8: Peel diagnostic (L1b + L2 effect classes) ---
    d8 <- 0.5  # neutral default
    l1b_eff <- iv$l1b_effect %||% iv$l1b_peel_effect %||% "unknown"
    l2_eff  <- iv$l2_effect %||% "unknown"

    # L1b score
    l1b_score <- if (l1b_eff == "stable") 1.0
                 else if (l1b_eff == "revealed_child") 0.9
                 else if (l1b_eff == "weakened") 0.5
                 else if (l1b_eff == "disappeared") 0.15
                 else 0.5  # unknown/ambiguous
    # L2 score
    l2_score <- if (l2_eff == "stable") 1.0
                else if (l2_eff == "revealed_child") 0.9
                else if (l2_eff == "weakened") 0.6
                else if (l2_eff == "disappeared") 0.1
                else 0.5
    d8 <- 0.55 * l1b_score + 0.45 * l2_score

    # Override with hypothesis verdict if available
    hyp_verdict <- "not_tested"
    if (nrow(hyp_verd) > 0) {
      hv <- hyp_verd[chrom == chr & interval_id == iid]
      if (nrow(hv) > 0) {
        hyp_verdict <- hv$verdict[1]
        hyp_d8 <- if (hyp_verdict %in% c("H2_or_H3_real_inversion",
                       "H2_broad_inversion_hot_core", "H3_nested_composite")) 1.0
                  else if (hyp_verdict %in% c("likely_real_inversion",
                            "H4_sub_haplotype_founder")) 0.7
                  else if (hyp_verdict %in% c("mixed_real_plus_family",
                            "mixed_family_and_signal")) 0.4
                  else if (hyp_verdict == "ambiguous_family_linked") 0.35
                  else if (hyp_verdict == "H1_family_structure") 0.1
                  else if (hyp_verdict == "H5_technical_possible") 0.15
                  else 0.5
        # Take the more informative of peel vs hypothesis
        d8 <- max(d8, hyp_d8)
      }
    }

    # --- D9: PCA/ICA cluster quality ---
    sil <- safe_num(iv$pc1_silhouette, 0)
    ica_k <- safe_num(iv$ica_kurtosis, 0)
    d9 <- pmin(1, pmax(0, 0.60 * sil + 0.40 * pmin(1, ica_k / 5)))

    # --- D10: Partition stability (GHSL-like) ---
    # BUGFIX 2026-04-17 (FIX 22, DESIGN): D10 now consumes 2e/C04's
    # GHSL v5 phased-Clair3 haplotype divergence when available, not
    # just D05's sim_mat-based partition stability. The two measure
    # different things and are complementary:
    #   - D05 sim_mat partition_stability: "do sample-PC1 cluster
    #     assignments stay consistent across the block's windows?"
    #     Population-structure signal; always available.
    #   - C04 ghsl_v5_score_max: "do within-sample haplotype
    #     divergences form a stable bimodal pattern?" Direct
    #     biological karyotype signal; available only where Clair3
    #     phasing has finished.
    # When both exist, blend with preference for the more direct
    # biological signal (ghsl). When only sim_mat, use it alone
    # (preserves pre-FIX-22 behaviour).
    part_stab <- safe_num(iv$partition_stability, 0)
    part_ent  <- safe_num(iv$partition_entropy, 0)
    ent_score <- pmin(1, pmax(0, 1 - part_ent / 2))
    d10_simmat <- 0.60 * part_stab + 0.40 * ent_score

    # GHSL v5 contribution (from C04 via FIX 22 merge in run_all.R)
    ghsl_score_max    <- safe_num(iv$ghsl_v5_score_max, NA_real_)
    ghsl_bimodal_frac <- safe_num(iv$ghsl_div_bimodal_frac, NA_real_)
    ghsl_pass_frac    <- safe_num(iv$ghsl_pass_frac, NA_real_)
    ghsl_n_scored     <- safe_num(iv$ghsl_n_scored_windows, 0)

    if (ghsl_n_scored >= 3 && is.finite(ghsl_score_max)) {
      # C04 has data for this block — build phased-Clair3 D10 component
      # ghsl_v5_score_max is already sigmoid-squashed to [0,1] in C04
      # (see STEP_C04_snake3_ghsl_v5.R L461-468). We take it as a
      # direct partition-quality proxy.
      d10_ghsl <- 0.50 * ghsl_score_max +
                  0.25 * safe_num(ghsl_bimodal_frac, 0) +
                  0.25 * safe_num(ghsl_pass_frac, 0)
      d10_ghsl <- pmin(1, pmax(0, d10_ghsl))
      # Blend: ghsl gets 0.6 weight when both are present (direct
      # biological evidence preferred over structural proxy)
      d10 <- 0.60 * d10_ghsl + 0.40 * d10_simmat
      d10_source <- "ghsl_and_simmat"
    } else {
      # Fall back to sim_mat partition stability alone
      d10 <- d10_simmat
      d10_source <- "simmat_only"
    }

    # --- D11: Boundary concordance (from C01g boundary catalog) ---
    # How well are this candidate's edges supported by independent evidence?
    d11 <- 0  # default: no boundary data
    boundary_verdict_left <- "unknown"
    boundary_verdict_right <- "unknown"
    n_boundary_cheats <- 0L

    # --- D12: Snake-staircase concordance ---
    # Does a snake core overlap this staircase block?
    d12 <- 0  # default: no snake data
    snake_overlap_status <- "no_data"

    # --- COMPOSITE SCORE (12 dimensions) ---
    # D1 (block) and D3 (NN) are strongest structural evidence.
    # D8 (peel) is the confound test. D11 (boundary) and D12 (snake) are new.
    final_score <- 0.14 * d1 + 0.08 * d2 + 0.13 * d3 +
                   0.06 * d4 + 0.05 * d5 + 0.06 * d6 +
                   0.06 * d7 + 0.12 * d8 + 0.06 * d9 + 0.06 * d10 +
                   0.09 * d11 + 0.09 * d12

    # --- TIER ASSIGNMENT (count-based, matches manuscript) ---
    dim_positive <- sum(c(
      d1  >= 0.30,   # D1: block has meaningful contrast + shape
      d2  >= 0.30,   # D2: good shape class
      d3  >= 0.30,   # D3: NN persistence present
      d4  >= 0.40,   # D4: flat decay (inversion-like)
      d5  >= 0.40,   # D5: homogeneous interior
      d6  >= 0.25,   # D6: found in ≥2 matrix variants
      d7  >= 0.20,   # D7: SV breakpoint match
      d8  >= 0.60,   # D8: peel says stable / hypothesis says real
      d9  >= 0.30,   # D9: PCA clustering quality
      d10 >= 0.30,   # D10: partition stability
      d11 >= 0.30,   # D11: boundary concordance
      d12 >= 0.30    # D12: snake concordance
    ))

    tier <- if (dim_positive >= 8) 1L
            else if (dim_positive >= 6) 2L
            else if (dim_positive >= 4) 3L
            else 4L
    # Downgrade if peel says disappeared on both levels
    if (l1b_eff == "disappeared" && l2_eff == "disappeared" && tier <= 2) tier <- 3L

    # --- PATTERN CLASSIFICATION ---
    n_children <- safe_num(iv$n_children %||% iv$n_nested, 0)
    shape_cl <- iv$shape_class %||% "unknown"
    landscape_cat <- iv$category %||% iv$landscape_category %||% "unclassified"

    pattern <- if (landscape_cat == "complex_system" || n_children >= 2) {
      "complex_system"
    } else if (landscape_cat == "nested_fixed") {
      "nested_fixed"
    } else if (landscape_cat == "nested_rare") {
      "nested_rare"
    } else if (landscape_cat == "family_ld_band") {
      "family_ld"
    } else if (shape_cl == "strong_square" && d1 >= 0.5) {
      "strong_inversion"
    } else if (shape_cl == "diffuse_square" && d1 >= 0.3) {
      "diffuse_inversion"
    } else if (shape_cl == "diagonal_band") {
      "diagonal_band"
    } else if (d1 < 0.2 && d2 < 0.2) {
      "noise"
    } else {
      "unclassified"
    }

    span_mb <- round(safe_num(iv$end_mb, 0) - safe_num(iv$start_mb, 0), 3)

    scores[[ri]] <- data.table(
      chrom = chr, interval_id = iid,
      start_mb = round(safe_num(iv$start_mb, 0), 3),
      end_mb = round(safe_num(iv$end_mb, 0), 3),
      span_mb = span_mb,
      # 12 scoring dimensions (D11 + D12 filled after loop at L819+)
      d1_block_strength = round(d1, 3),
      d2_block_shape = round(d2, 3),
      d3_nn_persistence = round(d3, 3),
      d4_decay_flatness = round(d4, 3),
      d5_interior_quality = round(d5, 3),
      d6_consensus = round(d6, 3),
      d7_sv_breakpoint = round(d7, 3),
      d8_peel_or_hyp = round(d8, 3),
      d9_pca_clusters = round(d9, 3),
      d10_partition = round(d10, 3),
      d10_source = d10_source,  # BUGFIX 2026-04-17 (FIX 22): simmat_only vs ghsl_and_simmat
      d11_boundary_concordance = round(d11, 3),
      d12_snake_concordance = round(d12, 3),
      # Summary
      final_score = round(final_score, 3),
      dim_positive = dim_positive,
      tier = tier,
      pattern = pattern,
      # Source evidence
      shape_class = shape_cl,
      landscape_category = landscape_cat,
      hyp_verdict = hyp_verdict,
      l1b_peel = l1b_eff,
      l2_peel = l2_eff,
      survives_nn40 = safe_num(iv$survives_nn40, NA),
      nn_birth = safe_num(iv$nn_birth, NA),
      n_variants = safe_num(iv$n_variants %||% n_var, NA),
      n_sv_hits = safe_num(iv$n_sv_hits, NA),
      snake_overlap = snake_overlap_status,
      boundary_verdict_left = boundary_verdict_left,
      boundary_verdict_right = boundary_verdict_right,
      n_boundary_cheats = n_boundary_cheats
    )
  }
  rbindlist(scores, fill = TRUE)
}

# =============================================================================
# SCORE ALL CANDIDATES
# =============================================================================

message("[C01d] Scoring candidates...")
cand_dt <- score_candidates(iv_dt, hyp_verd)
message("[C01d] Initial scoring: ", nrow(cand_dt), " candidates")

# =============================================================================
# D12: SNAKE-STAIRCASE CONCORDANCE (post-scoring adjustment)
# =============================================================================
# For each staircase-derived candidate, check if a snake core overlaps.
# Overlap → D12 high, eligible for Tier 1-2.
# No overlap → D12 = 0, cap at Tier 4.

if (!is.null(cores_dir) && dir.exists(cores_dir)) {
  message("[C01d] ── D12: Seeded-region / staircase concordance ──")
  # C01b_1 (phase_2/2c) writes seeded_regions_summary_<chr>.tsv.gz.
  # Legacy snake1_core_regions_<chr>.tsv.gz pattern kept as fallback for
  # pre-rename output directories.
  core_files <- list.files(cores_dir,
                            pattern = "^seeded_regions_summary_.*\\.tsv\\.gz$",
                            full.names = TRUE)
  if (length(core_files) == 0) {
    core_files <- list.files(cores_dir,
                              pattern = "^snake1_core_regions_.*\\.tsv\\.gz$",
                              full.names = TRUE)
    if (length(core_files) > 0) {
      message("[C01d] Reading legacy snake1_core_regions_*.tsv.gz (pre-rename)")
    }
  }
  if (length(core_files) > 0) {
    all_cores <- rbindlist(lapply(core_files, fread), fill = TRUE)
    # Back-compat: legacy input had core_family/snake_id/cheat26_status; new
    # input has scale_tier/region_id/test26_status. Normalise on read so the
    # overlap logic below doesn't need to know which it got.
    if ("scale_tier"   %in% names(all_cores) && !"core_family"    %in% names(all_cores)) all_cores[, core_family   := scale_tier]
    if ("region_id"    %in% names(all_cores) && !"snake_id"       %in% names(all_cores)) all_cores[, snake_id      := region_id]
    if ("test26_status"%in% names(all_cores) && !"cheat26_status" %in% names(all_cores)) all_cores[, cheat26_status := test26_status]
    message("[C01d] Loaded ", nrow(all_cores), " seeded regions from ", length(core_files), " chromosomes")

    for (ci in seq_len(nrow(cand_dt))) {
      cd <- cand_dt[ci]
      s_bp <- cd$start_mb * 1e6; e_bp <- cd$end_mb * 1e6

      # Find seeded regions that overlap this candidate (≥30% reciprocal overlap)
      chr_cores <- all_cores[chrom == cd$chrom]
      if (nrow(chr_cores) == 0) next

      for (ki in seq_len(nrow(chr_cores))) {
        cr <- chr_cores[ki]
        overlap_start <- max(s_bp, cr$start_bp)
        overlap_end <- min(e_bp, cr$end_bp)
        overlap_len <- max(0, overlap_end - overlap_start)
        cand_len <- e_bp - s_bp
        core_len <- cr$end_bp - cr$start_bp

        if (cand_len > 0 && core_len > 0) {
          recip_cand <- overlap_len / cand_len
          recip_core <- overlap_len / core_len

          if (recip_cand > 0.3 || recip_core > 0.3) {
            # Seeded region overlaps this staircase block
            overlap_frac <- max(recip_cand, recip_core)
            d12_val <- pmin(1, overlap_frac * 1.2)  # slight boost for strong overlap

            # Check test_26 kin-pruned retention status if available
            c26 <- cr$cheat26_status %||% "untested"
            if (c26 == "persists") d12_val <- pmin(1, d12_val + 0.2)
            else if (c26 == "collapsed") d12_val <- d12_val * 0.5

            if (d12_val > cand_dt$d12_snake_concordance[ci]) {
              set(cand_dt, ci, "d12_snake_concordance", round(d12_val, 3))
              set(cand_dt, ci, "snake_overlap", paste0(
                cr$core_family, "_region", cr$snake_id,
                "_", round(overlap_frac * 100), "pct",
                if (c26 != "untested") paste0("_t26=", c26) else ""))
            }
          }
        }
      }
    }

    # Tier cap: staircase blocks without any seeded-region overlap → Tier 4
    no_snake <- cand_dt$d12_snake_concordance == 0 & cand_dt$snake_overlap == "no_data"
    n_capped <- sum(no_snake & cand_dt$tier <= 3)
    cand_dt[no_snake & tier <= 3, tier := 4L]
    message("[C01d] Seeded-region concordance: ",
            sum(cand_dt$d12_snake_concordance > 0), " with overlap, ",
            n_capped, " capped to Tier 4 (no seeded region)")
  } else {
    message("[C01d] No seeded_regions_summary_*.tsv.gz files found in ", cores_dir)
  }
} else {
  message("[C01d] No --cores_dir — D12 seeded-region concordance skipped")
}

# =============================================================================
# D11: BOUNDARY CONCORDANCE (from C01g boundary catalog)
# =============================================================================
# For each candidate, find its left and right boundaries in the catalog,
# read their concordance scores, and compute D11.

if (!is.null(boundary_dir) && dir.exists(boundary_dir)) {
  message("[C01d] ── D11: Boundary concordance ──")
  bcat_file <- file.path(boundary_dir, "boundary_catalog_unified.tsv.gz")
  if (file.exists(bcat_file)) {
    bcat <- fread(bcat_file)
    message("[C01d] Loaded boundary catalog: ", nrow(bcat), " boundaries")

    for (ci in seq_len(nrow(cand_dt))) {
      cd <- cand_dt[ci]
      s_bp <- cd$start_mb * 1e6; e_bp <- cd$end_mb * 1e6

      chr_bounds <- bcat[chrom == cd$chrom]
      if (nrow(chr_bounds) == 0) next

      # Find nearest boundary to left edge
      left_dist <- abs(chr_bounds$boundary_bp - s_bp)
      left_idx <- which.min(left_dist)
      left_match <- if (left_dist[left_idx] < 100000) chr_bounds[left_idx] else NULL

      # Find nearest boundary to right edge
      right_dist <- abs(chr_bounds$boundary_bp - e_bp)
      right_idx <- which.min(right_dist)
      right_match <- if (right_dist[right_idx] < 100000) chr_bounds[right_idx] else NULL

      # Score D11 from boundary verdicts
      verdict_to_score <- function(v) {
        if (is.null(v)) return(0)
        vv <- v$boundary_verdict %||% "unresolved"
        if (vv == "confirmed_structural") 1.0
        else if (vv == "likely_structural") 0.8
        else if (vv == "probable_structural") 0.6
        else if (vv == "single_evidence_structural") 0.35
        else if (vv == "internal_feature") 0.25
        else if (vv == "system_junction") 0.5
        else 0.1  # unresolved / background
      }

      left_score <- verdict_to_score(left_match)
      right_score <- verdict_to_score(right_match)
      d11_val <- (left_score + right_score) / 2

      n_cheats_total <- 0L
      if (!is.null(left_match)) {
        n_cheats_total <- n_cheats_total + (left_match$n_cheats_supporting %||% 0L)
        set(cand_dt, ci, "boundary_verdict_left", left_match$boundary_verdict %||% "unknown")
      }
      if (!is.null(right_match)) {
        n_cheats_total <- n_cheats_total + (right_match$n_cheats_supporting %||% 0L)
        set(cand_dt, ci, "boundary_verdict_right", right_match$boundary_verdict %||% "unknown")
      }

      set(cand_dt, ci, "d11_boundary_concordance", round(d11_val, 3))
      set(cand_dt, ci, "n_boundary_cheats", n_cheats_total)
    }

    message("[C01d] Boundary concordance: ",
            sum(cand_dt$d11_boundary_concordance > 0.3, na.rm = TRUE),
            " candidates with supported boundaries")
  } else {
    message("[C01d] boundary_catalog_unified.tsv.gz not found in ", boundary_dir)
  }
} else {
  message("[C01d] No --boundary_dir — D11 boundary concordance skipped")
}

# =============================================================================
# CHEAT 25: BLOCK VIABILITY TEST (4-test battery → DEAD = Tier 4)
# =============================================================================
# Tests whether a candidate has ANY independent evidence beyond sim_mat.
# If all 4 tests are negative → the block is DEAD → forced Tier 4.
#
# Test 1: SV overlap (from D7) — any SV caller sees this region?
# Test 2: Boundary support (from D11) — are edges confirmed?
# Test 3: Peel survival (from D8) — does the block survive kin-pruning?
# Test 4: Snake core overlap (from D12) — did the snake find this region?

message("[C01d] ── Cheat 25: Block viability test ──")

cand_dt[, cheat25_status := "ALIVE"]  # default

for (ci in seq_len(nrow(cand_dt))) {
  cd <- cand_dt[ci]

  test1_sv    <- cd$d7_sv_breakpoint >= 0.20
  test2_bound <- cd$d11_boundary_concordance >= 0.30
  test3_peel  <- cd$d8_peel_or_hyp >= 0.50
  test4_snake <- cd$d12_snake_concordance >= 0.20

  n_alive <- sum(c(test1_sv, test2_bound, test3_peel, test4_snake))

  status <- if (n_alive >= 3) "ALIVE"
            else if (n_alive == 2) "UNCERTAIN"
            else if (n_alive == 1) "WEAK"
            else "DEAD"

  set(cand_dt, ci, "cheat25_status", status)

  # DEAD → force Tier 4
  if (status == "DEAD" && cand_dt$tier[ci] <= 3) {
    set(cand_dt, ci, "tier", 4L)
  }
  # WEAK → cap at Tier 3
  if (status == "WEAK" && cand_dt$tier[ci] <= 2) {
    set(cand_dt, ci, "tier", 3L)
  }
}

viab_tab <- table(cand_dt$cheat25_status)
message("[C01d] Viability: ", paste(names(viab_tab), viab_tab, sep = "=", collapse = " "))

# =============================================================================
# POPGEN ANNOTATION STREAM (Engine B dispatcher)
# =============================================================================
# Per-candidate Fst, theta_pi, theta_W, Tajima_D from region_popstats C binary.
# These are ANNOTATIONS for downstream use (C01f hypothesis, manuscript tables).
# NOT used as scoring dimensions — the 12 dimensions are structural evidence.
# Cost: ~2 seconds per candidate (100 SNPs × 226 samples).

if (.bridge_available && exists("get_region_stats", mode = "function")) {
  message("[C01d] ── Popgen annotation stream (Engine B) ──")

  # Initialize columns
  cand_dt[, `:=`(
    popgen_fst_b1b3 = NA_real_,
    popgen_theta_pi = NA_real_,
    popgen_theta_W  = NA_real_,
    popgen_Tajima_D = NA_real_,
    popgen_method   = NA_character_
  )]

  # Helper
  extract_fst_d <- function(s) {
    if (is.null(s$Fst) || length(s$Fst) == 0) return(NA_real_)
    as.numeric(s$Fst[[1]])
  }

  n_annotated <- 0L
  for (ci in seq_len(nrow(cand_dt))) {
    cd <- cand_dt[ci]
    chr <- cd$chrom
    s_bp <- cd$start_mb * 1e6
    e_bp <- cd$end_mb * 1e6
    if (!is.finite(s_bp) || !is.finite(e_bp) || e_bp <= s_bp) next

    # Get PC1 band assignments for this candidate from precomp
    pc_obj <- NULL
    if (!is.null(precomp_dir)) {
      pf <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
      if (file.exists(pf)) pc_obj <- readRDS(pf)
    }

    b1_cga <- character(0); b3_cga <- character(0)
    if (!is.null(pc_obj)) {
      dt_chr <- pc_obj$dt
      pc1_cols <- grep("^PC_1_", names(dt_chr), value = TRUE)
      inner_w <- which(dt_chr$start_bp >= s_bp & dt_chr$end_bp <= e_bp)
      if (length(inner_w) >= 3 && length(pc1_cols) >= 20) {
        avg_pc1 <- colMeans(as.matrix(dt_chr[inner_w, ..pc1_cols]), na.rm = TRUE)
        valid <- is.finite(avg_pc1)
        if (sum(valid) >= 20) {
          km <- tryCatch(kmeans(avg_pc1[valid], centers = 3, nstart = 5), error = function(e) NULL)
          if (!is.null(km)) {
            co <- order(km$centers[, 1])
            bands <- integer(sum(valid))
            for (b in 1:3) bands[km$cluster == co[b]] <- b
            snames <- sub("^PC_1_", "", names(avg_pc1)[valid])
            names(bands) <- snames

            b1_names <- snames[bands == 1]
            b3_names <- snames[bands == 3]

            # Map to CGA
            if (!is.null(smap) && grepl("^Ind[0-9]", b1_names[1])) {
              b1_cga <- smap$to_real_vec(b1_names)
              b3_cga <- smap$to_real_vec(b3_names)
            } else if (!is.null(real_names) && grepl("^Ind[0-9]", b1_names[1]) &&
                       length(real_names) > 0) {
              ind_to_real <- setNames(real_names, paste0("Ind", seq_along(real_names) - 1))
              b1_cga <- ind_to_real[b1_names]; b1_cga <- b1_cga[!is.na(b1_cga)]
              b3_cga <- ind_to_real[b3_names]; b3_cga <- b3_cga[!is.na(b3_cga)]
            } else {
              b1_cga <- b1_names; b3_cga <- b3_names
            }
          }
        }
      }
    }

    tryCatch({
      # Request Fst + theta in one call (C binary computes all together)
      what_req <- c("theta_pi", "theta_W", "Tajima_D")
      groups_arg <- NULL
      if (length(b1_cga) >= 5 && length(b3_cga) >= 5) {
        what_req <- c("Fst", what_req)
        groups_arg <- list(b1 = b1_cga, b3 = b3_cga)
      }

      stats <- get_region_stats(chr, s_bp, e_bp,
                                 what = what_req, groups = groups_arg)

      if (!is.null(stats$Fst))
        set(cand_dt, ci, "popgen_fst_b1b3", round(extract_fst_d(stats), 6))
      if (!is.null(stats$theta_pi))
        set(cand_dt, ci, "popgen_theta_pi", stats$theta_pi)
      if (!is.null(stats$theta_W))
        set(cand_dt, ci, "popgen_theta_W", stats$theta_W)
      if (!is.null(stats$Tajima_D))
        set(cand_dt, ci, "popgen_Tajima_D", stats$Tajima_D)
      set(cand_dt, ci, "popgen_method", "engine_b")
      n_annotated <- n_annotated + 1L
    }, error = function(e) NULL)

    if (ci %% 20 == 0) message("  ", ci, "/", nrow(cand_dt), " candidates")
  }
  message("[C01d] Popgen annotation: ", n_annotated, "/", nrow(cand_dt), " candidates annotated")

  # Also pull test_05 family_fst_ratio from precomp if available.
  # BUGFIX 2026-04-17: three issues fixed —
  #   (1) column renamed cheat5_family_fst_ratio → test05_family_fst_ratio
  #       in C01a during the cheat→test rename; this reader still looked
  #       for the old name. Now reads the new name with legacy fallback.
  #   (2) pc_obj does NOT have a top-level `inv_likeness` element — the
  #       column lives at pc_obj$dt$<col>. Was always NULL, so the loop
  #       silently did nothing. Now reads pc_obj$dt directly.
  #   (3) pc_obj$dt has start_bp/end_bp, NOT mid_bp. Use window-overlap
  #       against the candidate's bp range.
  # Output column name in cand_dt kept as cheat5_family_fst_ratio so
  # downstream readers (if any) don't need changing in this fix pass.
  if (!is.null(precomp_dir)) {
    cand_dt[, cheat5_family_fst_ratio := NA_real_]
    for (ci in seq_len(nrow(cand_dt))) {
      cd <- cand_dt[ci]; chr <- cd$chrom
      pf <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
      if (!file.exists(pf)) next
      pc_obj <- tryCatch(readRDS(pf), error = function(e) NULL)
      if (is.null(pc_obj) || is.null(pc_obj$dt)) next
      il <- pc_obj$dt
      # Accept both new name (test05) and legacy (cheat5) for backward
      # compatibility with pre-rename precomp caches.
      fst_col <- if ("test05_family_fst_ratio" %in% names(il)) {
        "test05_family_fst_ratio"
      } else if ("cheat5_family_fst_ratio" %in% names(il)) {
        "cheat5_family_fst_ratio"
      } else NULL
      if (is.null(fst_col)) next
      s_bp <- cd$start_mb * 1e6; e_bp <- cd$end_mb * 1e6
      # Windows whose center lies inside the candidate region
      mid_bp_vec <- (il$start_bp + il$end_bp) / 2
      region <- il[mid_bp_vec >= s_bp & mid_bp_vec <= e_bp]
      if (nrow(region) > 0) {
        set(cand_dt, ci, "cheat5_family_fst_ratio",
            round(mean(region[[fst_col]], na.rm = TRUE), 4))
      }
    }
    n_ratio <- sum(is.finite(cand_dt$cheat5_family_fst_ratio))
    if (n_ratio > 0) message("[C01d] test_05 family_fst_ratio: ", n_ratio, " candidates annotated")
  }
} else {
  message("[C01d] Popgen annotation skipped (no Engine B)")
}

# =============================================================================
# RECOMPUTE FINAL SCORE WITH D11 + D12
# =============================================================================
# Now that D11 and D12 are filled, recompute the composite score

cand_dt[, final_score := round(
  0.14 * d1_block_strength + 0.08 * d2_block_shape + 0.13 * d3_nn_persistence +
  0.06 * d4_decay_flatness + 0.05 * d5_interior_quality + 0.06 * d6_consensus +
  0.06 * d7_sv_breakpoint + 0.12 * d8_peel_or_hyp + 0.06 * d9_pca_clusters +
  0.06 * d10_partition + 0.09 * d11_boundary_concordance + 0.09 * d12_snake_concordance,
3)]

# Recompute dim_positive with D11+D12
cand_dt[, dim_positive := rowSums(cbind(
  d1_block_strength >= 0.30, d2_block_shape >= 0.30, d3_nn_persistence >= 0.30,
  d4_decay_flatness >= 0.40, d5_interior_quality >= 0.40, d6_consensus >= 0.25,
  d7_sv_breakpoint >= 0.20, d8_peel_or_hyp >= 0.60, d9_pca_clusters >= 0.30,
  d10_partition >= 0.30, d11_boundary_concordance >= 0.30, d12_snake_concordance >= 0.30
), na.rm = TRUE)]

# BUGFIX 2026-04-17 (chat 7, FIX 33): recompute tier after dim_positive
# recomputation. The tier assignment at L373–376 runs inside the
# per-candidate loop BEFORE D11 and D12 are filled (both default to 0
# in the loop). Since tier is count-based on dim_positive, every
# candidate with non-zero D11/D12 was getting a tier computed on a
# dim_positive that was up to 2 lower than the final value. Downstream
# readers saw dim_positive = 8, tier = 2 — internally inconsistent.
# Apply the same ≥8/≥6/≥4 ladder as L373–376 and re-apply the
# peel-disappeared downgrade so tier matches the final dim_positive.
cand_dt[, tier := fifelse(dim_positive >= 8L, 1L,
                    fifelse(dim_positive >= 6L, 2L,
                      fifelse(dim_positive >= 4L, 3L, 4L)))]
# Peel-disappeared downgrade (matches L378): if both peel effects say
# "disappeared", a candidate can't sit at Tier 1–2.
if ("l1b_peel" %in% names(cand_dt) && "l2_peel" %in% names(cand_dt)) {
  cand_dt[l1b_peel == "disappeared" & l2_peel == "disappeared" & tier <= 2L,
          tier := 3L]
}

# Re-sort by final_score
cand_dt <- cand_dt[order(-final_score)]

message("[C01d] Tiers: T1=", sum(cand_dt$tier == 1),
        " T2=", sum(cand_dt$tier == 2),
        " T3=", sum(cand_dt$tier == 3),
        " T4=", sum(cand_dt$tier == 4))
message("[C01d] Patterns: ", paste(names(table(cand_dt$pattern)),
        table(cand_dt$pattern), sep = "=", collapse = " "))

# =============================================================================
# MORPHOLOGY PA SUMMARY (from precomp, if available)
# =============================================================================

if (!is.null(precomp_dir) && dir.exists(precomp_dir)) {
  message("[C01d] Loading precomp for morphology PA from: ", precomp_dir)
  morph_rows <- list()

  for (ci in seq_len(nrow(cand_dt))) {
    cd <- cand_dt[ci]; chr <- cd$chrom
    start_bp <- cd$start_mb * 1e6; end_bp <- cd$end_mb * 1e6

    rds_f <- file.path(precomp_dir, paste0(chr, ".precomp.rds"))
    if (!file.exists(rds_f)) {
      morph_rows[[ci]] <- data.table(chrom = chr, interval_id = cd$interval_id)
      next
    }

    pc <- tryCatch(readRDS(rds_f), error = function(e) NULL)
    if (is.null(pc) || is.null(pc$dt)) {
      morph_rows[[ci]] <- data.table(chrom = chr, interval_id = cd$interval_id)
      next
    }

    dt_pc <- pc$dt
    win_idx <- which(dt_pc$start_bp >= start_bp & dt_pc$end_bp <= end_bp)
    if (length(win_idx) < 3) {
      morph_rows[[ci]] <- data.table(chrom = chr, interval_id = cd$interval_id)
      next
    }

    sub <- dt_pc[win_idx]
    sm <- function(x) { x <- x[is.finite(x)]; if (length(x) > 0) round(mean(x), 4) else NA_real_ }
    sx <- function(x) { x <- x[is.finite(x)]; if (length(x) > 0) round(max(x), 4) else NA_real_ }

    morph_rows[[ci]] <- data.table(
      chrom = chr, interval_id = cd$interval_id,
      pa_flat_inv_mean = sm(sub$flat_inv_score),
      pa_spiky_inv_mean = sm(sub$spiky_inv_score),
      pa_frag_mean = sm(sub$fragmentation_score),
      pa_family_like_mean = sm(sub$family_likeness),
      pa_jaggedness_mean = sm(sub$local_jaggedness),
      pa_block_compactness_mean = sm(sub$local_block_compactness)
    )
  }

  if (length(morph_rows) > 0) {
    morph_summ <- rbindlist(morph_rows, fill = TRUE)
    cand_dt <- merge(cand_dt, morph_summ, by = c("chrom", "interval_id"), all.x = TRUE)
    message("[C01d] Morphology PA added: ", ncol(morph_summ) - 2, " columns")
  }
} else {
  message("[C01d] No precomp_dir — morphology PA skipped")
}

# =============================================================================
# PLOTS
# =============================================================================

# Evidence heatmap
score_cols <- c("d1_block_strength", "d2_block_shape", "d3_nn_persistence",
                "d4_decay_flatness", "d5_interior_quality", "d6_consensus",
                "d7_sv_breakpoint", "d8_peel_or_hyp", "d9_pca_clusters",
                "d10_partition", "d11_boundary_concordance", "d12_snake_concordance")

top_n <- min(40, nrow(cand_dt))
top_cand <- cand_dt[order(-final_score)][seq_len(top_n)]
top_cand[, cand_label := paste0(chrom, ":", start_mb, "-", end_mb)]
top_cand[, cand_label := factor(cand_label, levels = rev(cand_label))]

melt_dt <- melt(top_cand, id.vars = c("cand_label", "tier", "pattern", "final_score"),
                measure.vars = score_cols, variable.name = "dimension", value.name = "score")

tryCatch({
  pC <- ggplot(melt_dt, aes(x = dimension, y = cand_label, fill = score)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colours = c("grey95", "lightyellow", "gold", "darkorange", "red3", "darkred"),
      name = "Score", limits = c(0, 1)
    ) +
    geom_text(aes(label = round(score, 2)), size = 2, color = "grey20") +
    labs(title = "Evidence Scoring Heatmap (inv_detect v9.3)",
         subtitle = paste0("Top ", top_n, " candidates by composite score"),
         x = NULL, y = NULL,
         caption = paste(
           "D1=block strength | D2=shape | D3=NN persistence | D4=decay flatness",
           "D5=interior quality | D6=consensus | D7=SV breakpoints",
           "D8=peel/hypothesis | D9=PCA clusters | D10=partition stability",
           sep = "\n")) +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 5))
  ggsave(file.path(outdir, "plots", "evidence_heatmap.png"),
         pC, width = 10, height = max(5, top_n * 0.3), dpi = DPI)
}, error = function(e) message("[PLOT] heatmap: ", e$message))

# Genome-wide ideogram
tier_pal <- c("1" = "red3", "2" = "darkorange", "3" = "steelblue", "4" = "grey60")
chr_sizes <- cand_dt[, .(max_mb = max(end_mb, na.rm = TRUE)), by = chrom]
chr_sizes[, chrom := factor(chrom, levels = sort(unique(as.character(chrom))))]

tryCatch({
  pA <- ggplot() +
    geom_segment(data = chr_sizes,
                 aes(x = 0, xend = max_mb, y = chrom, yend = chrom),
                 color = "grey80", linewidth = 2) +
    geom_segment(data = cand_dt,
                 aes(x = start_mb, xend = end_mb,
                     y = factor(chrom, levels = levels(chr_sizes$chrom)),
                     yend = factor(chrom, levels = levels(chr_sizes$chrom)),
                     color = factor(tier)),
                 linewidth = 4, alpha = 0.8) +
    scale_color_manual(values = tier_pal, name = "Tier") +
    labs(title = "Genome-wide Candidate Overview",
         subtitle = paste0(nrow(cand_dt), " candidates across ",
                          length(unique(cand_dt$chrom)), " chromosomes"),
         x = "Position (Mb)", y = NULL) +
    THEME_BASE +
    theme(axis.text.y = element_text(size = 6))
  ggsave(file.path(outdir, "plots", "genome_ideogram.png"),
         pA, width = 14, height = 8, dpi = DPI)
}, error = function(e) message("[PLOT] ideogram: ", e$message))

# Tier counts per chromosome
tier_counts <- cand_dt[, .(T1 = sum(tier == 1), T2 = sum(tier == 2),
                            T3 = sum(tier == 3), T4 = sum(tier == 4),
                            total = .N), by = chrom]

# =============================================================================
# WRITE
# =============================================================================

message("[C01d] Writing outputs...")
fwrite(cand_dt, file.path(outdir, "candidate_scores.tsv.gz"), sep = "\t")
fwrite(cand_dt[, .(chrom, interval_id, start_mb, end_mb, span_mb,
                    final_score, dim_positive, tier, pattern, shape_class)],
       file.path(outdir, "candidate_summary.tsv"), sep = "\t")
fwrite(tier_counts, file.path(outdir, "candidate_tier_counts.tsv"), sep = "\t")

message("\n[C01d] TIER SUMMARY:")
message("  Tier 1 (strong):  ", sum(cand_dt$tier == 1), " candidates")
message("  Tier 2 (moderate):", sum(cand_dt$tier == 2), " candidates")
message("  Tier 3 (complex): ", sum(cand_dt$tier == 3), " candidates")
message("  Tier 4 (weak):    ", sum(cand_dt$tier == 4), " candidates")


# ═══════════════════════════════════════════════════════════════════════
# REGISTER IN EVIDENCE REGISTRY (catalog birth)
# ═══════════════════════════════════════════════════════════════════════
tryCatch({
  # Source helpers if not already loaded
  for (.hf in c("utils/registry_key_helpers.R", "../utils/registry_key_helpers.R")) {
    if (file.exists(.hf)) { source(.hf); break }
  }

  if (.bridge_available && exists("reg") && !is.null(reg) && !is.null(reg$register_candidate)) {
    message("[C01d] Registering ", nrow(cand_dt), " candidates in evidence registry...")
    n_keys_total <- 0L
    for (ci in seq_len(nrow(cand_dt))) {
      cd <- cand_dt[ci]
      cid <- paste0(cd$chrom, "_", cd$interval_id)

      # Birth event: register the candidate
      reg$register_candidate(
        candidate_id = cid, chrom = cd$chrom,
        start_bp = as.integer(cd$start_mb * 1e6),
        end_bp = as.integer(cd$end_mb * 1e6),
        tier = cd$tier, score = cd$final_score
      )

      # Bulk register all Q1/Q7/popgen keys from this candidate row
      if (exists("register_C01d_keys", mode = "function")) {
        n_keys_total <- n_keys_total + register_C01d_keys(cd, cid, outdir)
      }
      # Write full scoring data to candidate folder
      if (exists("store_C01d_results", mode = "function")) {
        store_C01d_results(cd, cid, outdir)
      }
    }
    message("[C01d] Registered ", nrow(cand_dt), " candidates, ", n_keys_total, " evidence keys")
    reg$status(max_priority = 1)
  }
}, error = function(e) message("[C01d] Registry wiring: ", conditionMessage(e)))

message("\n[DONE] -> ", outdir)
