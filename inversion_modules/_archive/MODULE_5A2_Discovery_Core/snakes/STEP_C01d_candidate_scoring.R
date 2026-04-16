#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01d_candidate_scoring.R  (v8.4)
#
# TIERED EVIDENCE FRAMEWORK for inversion candidate scoring.
#
# Integrates outputs from:
#   C01a (precomp diagnostics)
#   C01c (triangle regimes, bridges, sub-regimes, off-diagonal)
#   C01a2 (network diagnostics, ancestry diversity, relatedness)
#   C01b (snake cores, merge regions)
#
# For each candidate region:
#   1. Triangle evidence (type, contrast, sharpness, sub-regime structure)
#   2. PCA band structure (trimodality, band sizes, symmetry)
#   3. Ancestry diversity (effective K per band)
#   4. Relatedness check (within vs between band theta)
#   5. Off-diagonal linkage (connected to other regions?)
#   6. Bridge status (part of a larger system?)
#   7. Snake core support (how many cores overlap this region?)
#
# Assigns each candidate:
#   - Per-dimension scores (0-1)
#   - Final composite score
#   - Tier (1=strong, 2=moderate, 3=complex, 4=weak/artifact)
#   - Pattern class (canonical_3band, asymmetric, split_het, diffuse, ancestry_like)
#
# Outputs:
#   candidate_scores.tsv.gz       -- per-candidate scoring matrix
#   candidate_summary.tsv         -- summary with tier + pattern
#   candidate_tier_counts.tsv     -- counts per tier per chr
#   plots/evidence_heatmap.png    -- Panel C of manuscript figure
#   plots/genome_ideogram.png     -- Panel A
#   plots/tier_summary.png        -- Panel B
#
# Usage:
#   Rscript STEP_C01d_candidate_scoring.R <triangle_dir> <outdir> \
#     [--cores_dir <dir>] [--network_dir <dir>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: see header")

triangle_dir <- args[1]; outdir <- args[2]
cores_dir <- NULL; network_dir <- NULL
i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--cores_dir" && i < length(args))   { cores_dir <- args[i+1]; i <- i+2L }
  else if (a == "--network_dir" && i < length(args)) { network_dir <- args[i+1]; i <- i+2L }
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
# LOAD UPSTREAM DATA
# =============================================================================

message("[C01d] Loading triangle data from: ", triangle_dir)

iv_dt <- fread(file.path(triangle_dir, "triangle_intervals.tsv.gz"))
comp_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_sample_composition.tsv.gz")),
                     error = function(e) data.table())
bridge_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_bridges.tsv.gz")),
                       error = function(e) data.table())
offdiag_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_offdiag_linkage.tsv.gz")),
                        error = function(e) data.table())
subreg_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_subregimes.tsv.gz")),
                       error = function(e) data.table())
subtri_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_subtriangles.tsv.gz")),
                       error = function(e) data.table())
cmp_dt <- tryCatch(fread(file.path(triangle_dir, "triangle_interval_comparison.tsv.gz")),
                    error = function(e) data.table())

message("[C01d] Intervals: ", nrow(iv_dt), " | Composition: ", nrow(comp_dt),
        " | Bridges: ", nrow(bridge_dt))

# Load hypothesis test verdicts (from C01f) if available
hyp_dir <- NULL; hyp_verd <- data.table()
hi <- match("--hyp_dir", args)
if (!is.na(hi) && hi < length(args)) hyp_dir <- args[hi + 1]
if (!is.null(hyp_dir)) {
  f <- file.path(hyp_dir, "hypothesis_verdicts.tsv")
  if (file.exists(f)) { hyp_verd <- fread(f); message("[C01d] Hypothesis verdicts: ", nrow(hyp_verd)) }
}

# Load tube graph verdicts (from C01g) if available
tube_dir <- NULL; tube_verd <- data.table()
ti <- match("--tube_dir", args)
if (!is.na(ti) && ti < length(args)) tube_dir <- args[ti + 1]
if (!is.null(tube_dir)) {
  f <- file.path(tube_dir, "tube_verdicts.tsv")
  if (file.exists(f)) { tube_verd <- fread(f); message("[C01d] Tube verdicts: ", nrow(tube_verd)) }
}

# Load GHSL Snake 3 track summary if available
ghsl_dir2 <- NULL; ghsl_track <- data.table()
gi <- match("--ghsl_dir", args)
if (!is.na(gi) && gi < length(args)) ghsl_dir2 <- args[gi + 1]
if (!is.null(ghsl_dir2)) {
  f <- file.path(ghsl_dir2, "snake3_track.tsv.gz")
  if (file.exists(f)) { ghsl_track <- fread(f); message("[C01d] GHSL track: ", nrow(ghsl_track)) }
}

# =============================================================================
# SCORING DIMENSIONS (10 total)
# =============================================================================

score_candidates <- function(iv_dt, comp_dt, bridge_dt, offdiag_dt, subreg_dt, subtri_dt,
                              hyp_verd, tube_verd, ghsl_track) {
  scores <- list()

  for (ri in seq_len(nrow(iv_dt))) {
    iv <- iv_dt[ri]
    chr <- iv$chrom; iid <- iv$interval_id

    # --- D1: Triangle strength (contrast + sharpness) ---
    d1_contrast <- pmin(1, iv$contrast / 0.08)  # saturates at 0.08
    d1_sharpness <- pmin(1, iv$sharpness / 0.04)
    d1 <- 0.6 * d1_contrast + 0.4 * d1_sharpness

    # --- D2: Band structure (trimodality from composition) ---
    d2 <- 0
    band_sizes <- c(0, 0, 0)
    band_symmetry <- 0
    if (nrow(comp_dt) > 0) {
      iv_comp <- comp_dt[chrom == chr & interval_id == iid]
      if (nrow(iv_comp) > 0) {
        for (b in 1:3) band_sizes[b] <- sum(iv_comp$band == paste0("band", b))
        total <- sum(band_sizes)
        if (total > 0) {
          # Trimodality: all 3 bands should have samples
          n_nonempty <- sum(band_sizes > 0)
          # Minimum band fraction (smallest band)
          min_frac <- min(band_sizes) / max(total, 1)
          # Symmetry: how balanced are band1 vs band3?
          if (band_sizes[1] > 0 && band_sizes[3] > 0) {
            band_symmetry <- 1 - abs(band_sizes[1] - band_sizes[3]) / max(band_sizes[1] + band_sizes[3], 1)
          }
          # Het fraction (band2 / total)
          het_frac <- band_sizes[2] / max(total, 1)

          d2 <- if (n_nonempty == 3 && min_frac > 0.05) {
            0.3 + 0.3 * min_frac / 0.2 + 0.2 * band_symmetry + 0.2 * pmin(1, het_frac / 0.5)
          } else if (n_nonempty == 3) {
            0.3 + 0.2 * band_symmetry
          } else if (n_nonempty == 2) {
            0.2  # two bands only
          } else {
            0.0  # unimodal
          }
          d2 <- pmin(1, d2)
        }
      }
    }

    # --- D3: Ancestry diversity (effective K per band) ---
    d3 <- 0
    if (nrow(comp_dt) > 0 && "eff_K" %in% names(comp_dt)) {
      iv_comp <- comp_dt[chrom == chr & interval_id == iid & is.finite(eff_K)]
      if (nrow(iv_comp) > 0) {
        mean_eff_k <- mean(iv_comp$eff_K, na.rm = TRUE)
        # High eff_K = multiple families in same band = real inversion
        d3 <- pmin(1, (mean_eff_k - 1) / 4)  # saturates at eff_K = 5
      }
    }

    # --- D4: Off-diagonal linkage ---
    d4 <- 0
    if (nrow(offdiag_dt) > 0) {
      od_links <- offdiag_dt[(interval_id_1 == iid | interval_id_2 == iid) &
                              chrom == chr]
      if (nrow(od_links) > 0) {
        max_contrast <- max(od_links$cross_contrast, na.rm = TRUE)
        d4 <- pmin(1, max_contrast / 0.05)
      }
    }

    # --- D5: Sub-regime structure (internal heterogeneity) ---
    d5 <- 0
    if ("n_hot_subregimes" %in% names(iv) && "frac_hot" %in% names(iv)) {
      d5 <- pmin(1, iv$frac_hot * 1.5)  # reward hot sub-regimes
    }

    # --- D6: Bridge membership (part of larger system?) ---
    d6_bridge <- 0
    if (nrow(bridge_dt) > 0) {
      br <- bridge_dt[(id_a == iid | id_b == iid) & chrom == chr]
      if (nrow(br) > 0) {
        max_br <- max(br$bridge_score, na.rm = TRUE)
        d6_bridge <- pmin(1, max_br)
      }
    }

    # --- D7: Recursive depth (nested sub-triangles) ---
    d7 <- 0
    if (nrow(subtri_dt) > 0 && "parent_interval" %in% names(subtri_dt)) {
      n_nested <- nrow(subtri_dt[parent_interval == iid & chrom == chr])
      d7 <- pmin(1, n_nested / 5)
    }

    # --- D8: Hypothesis test verdict (from C01f) ---
    d8 <- 0.5  # neutral default if no hypothesis tests run
    hyp_verdict <- "not_tested"
    if (nrow(hyp_verd) > 0) {
      hv <- hyp_verd[chrom == chr & interval_id == iid]
      if (nrow(hv) > 0) {
        hyp_verdict <- hv$verdict[1]
        d8 <- if (hyp_verdict %in% c("H2_or_H3_real_inversion", "H2_broad_inversion_hot_core",
                                       "H3_nested_composite")) 1.0
              else if (hyp_verdict %in% c("likely_real_inversion", "H4_sub_haplotype_founder")) 0.7
              else if (hyp_verdict %in% c("mixed_real_plus_family", "mixed_family_and_signal")) 0.4
              else if (hyp_verdict == "H1_family_structure") 0.1
              else if (hyp_verdict == "H5_technical_possible") 0.15
              else 0.5  # unresolved
      }
    }

    # --- D9: Tube graph Stage E verdict (from C01g) ---
    d9 <- 0.5  # neutral default
    tube_stage_e <- "not_tested"
    tube_stage_d <- "not_tested"
    if (nrow(tube_verd) > 0) {
      tv <- tube_verd[chrom == chr & interval_id == iid]
      if (nrow(tv) > 0) {
        tube_stage_e <- tv$stage_e[1]
        tube_stage_d <- tv$stage_d[1]
        d9 <- if (tube_stage_e %in% c("broad_clean_inversion", "hot_core_inversion")) 1.0
              else if (tube_stage_e %in% c("recombinant_confirmed", "recombinant_ghsl_confirmed",
                                            "double_crossover_confirmed", "double_crossover_ghsl_confirmed")) 0.9
              else if (tube_stage_e %in% c("simple_recombinant", "double_crossover_like")) 0.7
              else if (tube_stage_e == "nested_composite") 0.6
              else if (tube_stage_e == "complex_mosaic") 0.4
              else if (tube_stage_e == "adjacent_systems") 0.3
              else 0.5
        # Penalize if Stage D says confounded
        if (tube_stage_d == "family_dominated") d9 <- d9 * 0.5
        if (tube_stage_d == "unstable_middle") d9 <- d9 * 0.7
      }
    }

    # --- D10: GHSL Snake 3 support (independent haplotype contrast) ---
    d10 <- 0
    ghsl_mean <- NA_real_
    if (nrow(ghsl_track) > 0) {
      # Match GHSL windows overlapping this interval
      start_bp <- iv$start_mb * 1e6; end_bp <- iv$end_mb * 1e6
      gw <- ghsl_track[chrom == chr & start_bp >= start_bp & end_bp <= end_bp]
      if (nrow(gw) > 0) {
        n_pass <- sum(gw$snake3_status == "PASS")
        ghsl_mean <- mean(gw$mean_ghsl, na.rm = TRUE)
        pass_frac <- n_pass / nrow(gw)
        d10 <- pmin(1, pass_frac * 0.6 + pmin(1, ghsl_mean / 5) * 0.4)
      }
    }

    # --- COMPOSITE SCORE (10 dimensions) ---
    final_score <- 0.18 * d1 + 0.15 * d2 + 0.10 * d3 +
                   0.07 * d4 + 0.05 * d5 + 0.05 * d6_bridge + 0.05 * d7 +
                   0.15 * d8 + 0.12 * d9 + 0.08 * d10

    # --- TIER ASSIGNMENT (uses hypothesis + tube when available) ---
    tier <- if (final_score >= 0.6 && d1 >= 0.4 && d8 >= 0.7) {
      1L  # Strong: good triangle + hypothesis confirmed
    } else if (final_score >= 0.5 && d1 >= 0.3 && d2 >= 0.3) {
      1L  # Strong: good triangle + good bands even without hypothesis test
    } else if (final_score >= 0.4 && d1 >= 0.3) {
      2L  # Moderate
    } else if (final_score >= 0.25 || d4 >= 0.5 || d9 >= 0.7) {
      3L  # Complex: supported by tube or off-diagonal
    } else {
      4L  # Weak
    }
    # Downgrade if hypothesis says family structure
    if (d8 <= 0.15 && tier <= 2) tier <- 3L

    # --- PATTERN CLASSIFICATION (enriched with tube verdicts) ---
    pattern <- if (tube_stage_e %in% c("recombinant_confirmed", "recombinant_ghsl_confirmed",
                                        "simple_recombinant")) {
      "recombinant"
    } else if (tube_stage_e %in% c("double_crossover_confirmed", "double_crossover_ghsl_confirmed",
                                    "double_crossover_like")) {
      "double_crossover"
    } else if (tube_stage_e == "complex_mosaic") {
      "complex_mosaic"
    } else if (d2 >= 0.7 && band_symmetry >= 0.6) {
      "canonical_3band"
    } else if (d2 >= 0.5 && band_symmetry < 0.4) {
      "asymmetric"
    } else if (d2 >= 0.3 && band_sizes[2] > band_sizes[1] + band_sizes[3]) {
      "split_het"
    } else if (d3 < 0.2 && d1 >= 0.3) {
      "ancestry_like"
    } else if (d1 < 0.3 && d2 < 0.3) {
      "diffuse"
    } else {
      "unclassified"
    }

    scores[[ri]] <- data.table(
      chrom = chr, interval_id = iid,
      start_mb = iv$start_mb, end_mb = iv$end_mb,
      span_mb = round(iv$end_mb - iv$start_mb, 3),
      interval_type = iv$interval_type,
      d1_triangle = round(d1, 3), d2_bands = round(d2, 3),
      d3_ancestry = round(d3, 3), d4_offdiag = round(d4, 3),
      d5_subregime = round(d5, 3), d6_bridge = round(d6_bridge, 3),
      d7_nesting = round(d7, 3), d8_hypothesis = round(d8, 3),
      d9_tube = round(d9, 3), d10_ghsl = round(d10, 3),
      final_score = round(final_score, 3),
      tier = tier, pattern = pattern,
      hyp_verdict = hyp_verdict, tube_stage_d = tube_stage_d, tube_stage_e = tube_stage_e,
      ghsl_mean = round(ghsl_mean, 4),
      band1_n = band_sizes[1], band2_n = band_sizes[2], band3_n = band_sizes[3],
      band_symmetry = round(band_symmetry, 3)
    )
  }
  rbindlist(scores)
}

# =============================================================================
# SCORE ALL CANDIDATES
# =============================================================================

message("[C01d] Scoring candidates...")
cand_dt <- score_candidates(iv_dt, comp_dt, bridge_dt, offdiag_dt, subreg_dt, subtri_dt,
                             hyp_verd, tube_verd, ghsl_track)
message("[C01d] Candidates scored: ", nrow(cand_dt))
message("[C01d] Tiers: T1=", sum(cand_dt$tier == 1),
        " T2=", sum(cand_dt$tier == 2),
        " T3=", sum(cand_dt$tier == 3),
        " T4=", sum(cand_dt$tier == 4))
message("[C01d] Patterns: ", paste(names(table(cand_dt$pattern)),
        table(cand_dt$pattern), sep = "=", collapse = " "))

# =============================================================================
# PA AGGREGATION LAYER (addon — does not replace tier gating)
# =============================================================================
# Load all PA matrices if available, join by window, produce per-candidate
# cross-layer agreement summary. Adds columns to cand_dt.

pa_dir <- NULL
pi <- match("--pa_dir", args)
if (!is.na(pi) && pi < length(args)) pa_dir <- args[pi + 1]

# Also check for individual PA files via other dirs
core_pa_file <- NULL; merge_pa_file <- NULL; tri_pa_file <- NULL; regime_pa_file <- NULL
if (!is.null(pa_dir)) {
  core_pa_file <- file.path(pa_dir, "snake1_window_states.tsv.gz")
  merge_pa_file <- file.path(pa_dir, "snake_merge_window_pa.tsv.gz")
  tri_pa_file <- file.path(pa_dir, "triangle_window_pa.tsv.gz")
  regime_pa_file <- file.path(pa_dir, "regime_window_pa.tsv.gz")
}
# Also try from individual dirs
ci <- match("--core_pa", args); if (!is.na(ci) && ci < length(args)) core_pa_file <- args[ci + 1]
mi <- match("--merge_pa", args); if (!is.na(mi) && mi < length(args)) merge_pa_file <- args[mi + 1]
tp <- match("--tri_pa", args); if (!is.na(tp) && tp < length(args)) tri_pa_file <- args[tp + 1]
rp <- match("--regime_pa", args); if (!is.na(rp) && rp < length(args)) regime_pa_file <- args[rp + 1]

core_pa <- if (!is.null(core_pa_file) && file.exists(core_pa_file)) fread(core_pa_file) else data.table()
merge_pa <- if (!is.null(merge_pa_file) && file.exists(merge_pa_file)) fread(merge_pa_file) else data.table()
tri_pa <- if (!is.null(tri_pa_file) && file.exists(tri_pa_file)) fread(tri_pa_file) else data.table()
regime_pa <- if (!is.null(regime_pa_file) && file.exists(regime_pa_file)) fread(regime_pa_file) else data.table()

has_pa <- nrow(core_pa) > 0 || nrow(merge_pa) > 0 || nrow(tri_pa) > 0 || nrow(regime_pa) > 0

if (has_pa) {
  message("[C01d] PA aggregation: core=", nrow(core_pa), " merge=", nrow(merge_pa),
          " tri=", nrow(tri_pa), " regime=", nrow(regime_pa))

  pa_summary_rows <- list()
  for (ci in seq_len(nrow(cand_dt))) {
    cd <- cand_dt[ci]; chr <- cd$chrom
    start_mb <- cd$start_mb; end_mb <- cd$end_mb

    # Core PA: windows in this region
    core_in <- if (nrow(core_pa) > 0) {
      core_pa[chrom == chr & pos_mb >= start_mb & pos_mb <= end_mb]
    } else data.table()

    # Merge PA
    merge_in <- if (nrow(merge_pa) > 0) {
      merge_pa[chrom == chr & pos_mb >= start_mb & pos_mb <= end_mb]
    } else data.table()

    # Triangle PA
    tri_in <- if (nrow(tri_pa) > 0) {
      tri_pa[chrom == chr & pos_mb >= start_mb & pos_mb <= end_mb]
    } else data.table()

    # Regime PA
    regime_in <- if (nrow(regime_pa) > 0) {
      regime_pa[chrom == chr & pos_mid_mb >= start_mb & pos_mid_mb <= end_mb]
    } else data.table()

    n_windows_region <- max(1, nrow(core_in))

    # Core summary
    pa_core_SML <- if (nrow(core_in) > 0 && "pa_pattern" %in% names(core_in))
      sum(core_in$pa_pattern == "SML") / n_windows_region else NA
    pa_core_any <- if (nrow(core_in) > 0 && "n_families_claimed" %in% names(core_in))
      sum(core_in$n_families_claimed > 0) / n_windows_region else NA
    pa_seed_unc <- if (nrow(core_in) > 0 && "collector_state" %in% names(core_in))
      sum(core_in$collector_state == "seed_uncollected") / n_windows_region else NA

    # Merge summary
    pa_merge_ABC <- if (nrow(merge_in) > 0 && "merge_pa" %in% names(merge_in))
      sum(merge_in$merge_pa == "ABC") / max(1, nrow(merge_in)) else NA
    pa_merge_bridge <- if (nrow(merge_in) > 0 && "merge_role" %in% names(merge_in))
      sum(merge_in$merge_role == "bridge") / max(1, nrow(merge_in)) else NA
    pa_merge_gap_sig <- if (nrow(merge_in) > 0 && "merge_role" %in% names(merge_in))
      sum(merge_in$merge_role == "gap_with_signal") / max(1, nrow(merge_in)) else NA

    # Triangle summary
    pa_tri_strong <- if (nrow(tri_in) > 0 && "tri_state" %in% names(tri_in))
      sum(tri_in$tri_state == "strong_triangle") / max(1, nrow(tri_in)) else NA
    pa_tri_outside <- if (nrow(tri_in) > 0 && "tri_state" %in% names(tri_in))
      sum(tri_in$tri_state == "outside") / max(1, nrow(tri_in)) else NA

    # Regime summary
    pa_regime_structured <- if (nrow(regime_in) > 0 && "regime_state" %in% names(regime_in))
      sum(regime_in$regime_state %in% c("clean_inversion", "structured_moderate",
          "structured_complex")) / max(1, nrow(regime_in)) else NA
    pa_regime_soup <- if (nrow(regime_in) > 0 && "regime_state" %in% names(regime_in))
      sum(regime_in$regime_state == "background_soup") / max(1, nrow(regime_in)) else NA
    pa_regime_mean_struct <- if (nrow(regime_in) > 0 && "regime_structure_score" %in% names(regime_in))
      mean(regime_in$regime_structure_score, na.rm = TRUE) else NA

    # Cross-layer agreement
    layers_agree <- sum(c(
      !is.na(pa_core_any) && pa_core_any > 0.5,
      !is.na(pa_merge_ABC) && pa_merge_ABC > 0.3,
      !is.na(pa_tri_strong) && pa_tri_strong > 0.5,
      !is.na(pa_regime_structured) && pa_regime_structured > 0.5
    ))

    pa_summary_rows[[ci]] <- data.table(
      chrom = chr, interval_id = cd$interval_id,
      pa_core_frac_SML = round(pa_core_SML, 3),
      pa_core_frac_any = round(pa_core_any, 3),
      pa_seed_uncollected = round(pa_seed_unc, 3),
      pa_merge_frac_ABC = round(pa_merge_ABC, 3),
      pa_merge_frac_bridge = round(pa_merge_bridge, 3),
      pa_merge_gap_signal = round(pa_merge_gap_sig, 3),
      pa_tri_frac_strong = round(pa_tri_strong, 3),
      pa_tri_frac_outside = round(pa_tri_outside, 3),
      pa_regime_frac_structured = round(pa_regime_structured, 3),
      pa_regime_frac_soup = round(pa_regime_soup, 3),
      pa_regime_mean_structure = round(pa_regime_mean_struct, 3),
      pa_layers_agree = layers_agree
    )
  }

  if (length(pa_summary_rows) > 0) {
    pa_summ <- rbindlist(pa_summary_rows, fill = TRUE)
    cand_dt <- merge(cand_dt, pa_summ, by = c("chrom", "interval_id"), all.x = TRUE)
    message("[C01d] PA columns added: ", ncol(pa_summ) - 2, " new columns")
    message("[C01d] Layer agreement: ",
            sum(cand_dt$pa_layers_agree == 4, na.rm = TRUE), " all-4, ",
            sum(cand_dt$pa_layers_agree == 3, na.rm = TRUE), " 3-of-4, ",
            sum(cand_dt$pa_layers_agree <= 2, na.rm = TRUE), " <=2")
  }
}

# =============================================================================
# PLOTS
# =============================================================================

# --- Panel C: Evidence scoring heatmap ---
# Top candidates sorted by final_score
top_n <- min(40, nrow(cand_dt))
top_cand <- cand_dt[order(-final_score)][seq_len(top_n)]
top_cand[, cand_label := paste0(chrom, ":", start_mb, "-", end_mb)]
top_cand[, cand_label := factor(cand_label, levels = rev(cand_label))]

# Melt score dimensions for heatmap
score_cols <- c("d1_triangle", "d2_bands", "d3_ancestry", "d4_offdiag",
                "d5_subregime", "d6_bridge", "d7_nesting",
                "d8_hypothesis", "d9_tube", "d10_ghsl")
melt_dt <- melt(top_cand, id.vars = c("cand_label", "tier", "pattern", "final_score"),
                measure.vars = score_cols, variable.name = "dimension", value.name = "score")

pC <- ggplot(melt_dt, aes(x = dimension, y = cand_label, fill = score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("grey95", "lightyellow", "gold", "darkorange", "red3", "darkred"),
    name = "Score", limits = c(0, 1)
  ) +
  geom_text(aes(label = round(score, 2)), size = 2, color = "grey20") +
  # Add tier + pattern annotation columns
  geom_tile(data = top_cand, aes(x = "Tier", y = cand_label,
            fill = fifelse(tier == 1, 0.9, fifelse(tier == 2, 0.6, fifelse(tier == 3, 0.3, 0.1)))),
            color = "white", linewidth = 0.3) +
  geom_text(data = top_cand, aes(x = "Tier", y = cand_label, label = tier),
            size = 2.5, fontface = "bold") +
  scale_x_discrete(labels = c(score_cols, "Tier")) +
  labs(title = "Evidence Scoring Heatmap",
       subtitle = paste0("Top ", top_n, " candidates by composite score"),
       x = NULL, y = NULL,
       caption = "D1=triangle | D2=bands | D3=ancestry | D4=off-diagonal\nD5=sub-regime | D6=bridge | D7=nesting | D8=hypothesis test\nD9=tube graph | D10=GHSL haplotype contrast") +
  THEME_BASE +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 5))

tryCatch(ggsave(file.path(outdir, "plots", "evidence_heatmap.png"),
                 pC, width = 10, height = max(5, top_n * 0.3), dpi = DPI),
         error = function(e) message("[PLOT] heatmap: ", e$message))

# --- Panel A: Genome-wide ideogram ---
tier_pal <- c("1" = "red3", "2" = "darkorange", "3" = "steelblue", "4" = "grey60")

# Chromosome sizes from intervals
chr_sizes <- cand_dt[, .(max_mb = max(end_mb, na.rm = TRUE)), by = chrom]
chr_sizes[, chrom := factor(chrom, levels = sort(unique(as.character(chrom))))]

pA <- ggplot() +
  # Chromosome backbones
  geom_segment(data = chr_sizes,
               aes(x = 0, xend = max_mb, y = chrom, yend = chrom),
               color = "grey80", linewidth = 2) +
  # Candidate regions
  geom_segment(data = cand_dt,
               aes(x = start_mb, xend = end_mb,
                   y = factor(chrom, levels = levels(chr_sizes$chrom)),
                   yend = factor(chrom, levels = levels(chr_sizes$chrom)),
                   color = factor(tier)),
               linewidth = 4, alpha = 0.8) +
  scale_color_manual(values = tier_pal, name = "Tier") +
  labs(title = "Genome-wide Candidate Overview",
       subtitle = paste0(nrow(cand_dt), " candidates across ", length(unique(cand_dt$chrom)),
                        " chromosomes"),
       x = "Position (Mb)", y = NULL) +
  THEME_BASE +
  theme(axis.text.y = element_text(size = 6))

tryCatch(ggsave(file.path(outdir, "plots", "genome_ideogram.png"),
                 pA, width = 14, height = 8, dpi = DPI),
         error = function(e) message("[PLOT] ideogram: ", e$message))

# --- Panel B: Tier summary ---
if (nrow(cand_dt) > 0) {
  pattern_pal <- c("canonical_3band" = "green4", "asymmetric" = "purple3",
                    "split_het" = "darkorange", "ancestry_like" = "grey60",
                    "diffuse" = "plum3", "recombinant" = "red3",
                    "double_crossover" = "darkred", "complex_mosaic" = "brown",
                    "unclassified" = "grey80")

  pB <- ggplot(cand_dt, aes(x = d1_triangle, y = d2_bands,
                              color = pattern, shape = factor(tier))) +
    geom_point(aes(size = span_mb), alpha = 0.7) +
    scale_color_manual(values = pattern_pal, name = "Pattern") +
    scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15, "4" = 1),
                       name = "Tier") +
    scale_size_continuous(range = c(1, 6), name = "Span (Mb)") +
    facet_wrap(~ paste0("Tier ", tier), ncol = 4) +
    labs(title = "Candidate Classification Summary",
         subtitle = "D1 (triangle) vs D2 (bands) colored by pattern class",
         x = "D1: Triangle strength", y = "D2: Band structure") +
    THEME_BASE

  tryCatch(ggsave(file.path(outdir, "plots", "tier_summary.png"),
                   pB, width = 14, height = 5, dpi = DPI),
           error = function(e) message("[PLOT] tier: ", e$message))
}

# --- Tier counts per chromosome ---
tier_counts <- cand_dt[, .(T1 = sum(tier == 1), T2 = sum(tier == 2),
                            T3 = sum(tier == 3), T4 = sum(tier == 4),
                            total = .N), by = chrom]

# =============================================================================
# WRITE
# =============================================================================

message("[C01d] Writing outputs...")
fwrite(cand_dt, file.path(outdir, "candidate_scores.tsv.gz"), sep = "\t")
fwrite(cand_dt[, .(chrom, interval_id, start_mb, end_mb, span_mb,
                    final_score, tier, pattern, interval_type)],
       file.path(outdir, "candidate_summary.tsv"), sep = "\t")
fwrite(tier_counts, file.path(outdir, "candidate_tier_counts.tsv"), sep = "\t")

message("\n[C01d] TIER SUMMARY:")
message("  Tier 1 (strong):  ", sum(cand_dt$tier == 1), " candidates")
message("  Tier 2 (moderate):", sum(cand_dt$tier == 2), " candidates")
message("  Tier 3 (complex): ", sum(cand_dt$tier == 3), " candidates")
message("  Tier 4 (weak):    ", sum(cand_dt$tier == 4), " candidates")
message("\n[DONE] -> ", outdir)
