#!/usr/bin/env Rscript
# =============================================================================
# cheat25_block_viability.R — Block viability test (4-test battery)
#
# BIOLOGY:
#   Not every candidate block with elevated inv_likeness is a real
#   inversion. This test applies four orthogonal checks to determine
#   whether a block is ALIVE (real signal), UNCERTAIN, or DEAD (noise).
#
# TEST A: SV concordance — flashlight genotypes match PCA bands?
# TEST B: Theta heterozygosity — HET-band elevated θ_P?
# TEST C: Indel slope — het/hom ratio differs between bands?
# TEST D: Band stability — same samples in same bands across >80% windows?
#
# ≥3 pass → ALIVE, 1-2 → UNCERTAIN, 0 → DEAD
#
# INPUT:  precomp dt, sample names, flashlight data, theta data
# OUTPUT: viability status, per-test results
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
SV_CONC_THRESHOLD     <- 0.60   # Test A: concordance rate
THETA_SEP_THRESHOLD   <- 0.50   # Test B: band separation z-score
INDEL_RATIO_THRESHOLD <- 1.3    # Test C: het/hom indel ratio
BAND_STABILITY_FRAC   <- 0.80   # Test D: fraction of windows stable
MIN_WINDOWS_VIABILITY <- 5L

# ── Test A: SV concordance ────────────────────────────────────────────

test_sv_concordance <- function(chr, candidate_start, candidate_end,
                                 band_assignments, fl = NULL) {
  if (is.null(fl)) {
    # Try loading via flashlight functions if available
    fl_data <- tryCatch(
      if (exists("load_flashlight", mode = "function"))
        load_flashlight(chr) else NULL,
      error = function(e) NULL)
  } else {
    fl_data <- fl
  }

  if (is.null(fl_data)) return(list(pass = NA, concordance = NA_real_,
                                      detail = "no_flashlight"))

  anch <- tryCatch(
    get_sv_anchors(chr, candidate_start, candidate_end),
    error = function(e) data.table())

  if (nrow(anch) == 0) return(list(pass = NA, concordance = NA_real_,
                                     detail = "no_anchors"))

  conc <- tryCatch(
    get_anchor_concordance(chr, candidate_start, candidate_end, band_assignments),
    error = function(e) list(concordance = NA_real_))

  pass <- !is.na(conc$concordance) && conc$concordance >= SV_CONC_THRESHOLD

  list(pass = pass, concordance = round(conc$concordance, 3),
       detail = paste0("n_anchors=", nrow(anch),
                        ",concordance=", round(conc$concordance, 3)))
}

# ── Test B: Theta heterozygosity ──────────────────────────────────────

test_theta_het <- function(chr, candidate_start, candidate_end,
                            band_assignments, theta_dir = NULL) {
  if (is.null(theta_dir) || !nzchar(theta_dir))
    theta_dir <- Sys.getenv("THETA_DIR", "")
  if (!nzchar(theta_dir))
    return(list(pass = NA, separation = NA_real_, detail = "no_theta_dir"))

  # Attempt to use cheat12 infrastructure
  if (!exists("compute_theta_zscores", mode = "function"))
    return(list(pass = NA, separation = NA_real_,
                detail = "cheat12_not_loaded"))

  sample_ids <- names(band_assignments)
  theta_z <- tryCatch(
    compute_theta_zscores(theta_dir, sample_ids, chr,
                           candidate_start, candidate_end),
    error = function(e) data.table())

  if (nrow(theta_z) == 0)
    return(list(pass = NA, separation = NA_real_, detail = "no_theta_data"))

  band_result <- tryCatch(
    compute_band_theta(theta_z, band_assignments),
    error = function(e) list(per_band = data.table()))

  pb <- band_result$per_band
  if (nrow(pb) < 2)
    return(list(pass = NA, separation = NA_real_, detail = "too_few_bands"))

  het_band <- pb[band == 2]$band_mean_tP_z
  hom_bands <- pb[band %in% c(1, 3)]$band_mean_tP_z
  sep <- if (length(het_band) > 0 && length(hom_bands) > 0)
    het_band - mean(hom_bands) else NA_real_

  pass <- !is.na(sep) && sep > THETA_SEP_THRESHOLD

  list(pass = pass, separation = round(sep, 3),
       detail = paste0("het_z=", round(het_band, 3),
                        ",hom_z=", round(mean(hom_bands), 3)))
}

# ── Test C: Indel slope ───────────────────────────────────────────────

test_indel_slope <- function(chr, candidate_start, candidate_end,
                              band_assignments) {
  # Attempt to use cheat9 infrastructure
  if (!exists("classify_by_indel_slope", mode = "function"))
    return(list(pass = NA, ratio = NA_real_,
                detail = "cheat9_not_loaded"))

  indel_result <- tryCatch(
    classify_by_indel_slope(chr, candidate_start, candidate_end,
                             names(band_assignments)),
    error = function(e) NULL)

  if (is.null(indel_result))
    return(list(pass = NA, ratio = NA_real_, detail = "indel_failed"))

  # Check if het band has elevated indel het/hom ratio
  if (is.data.table(indel_result) && "het_hom_ratio" %in% names(indel_result)) {
    band2_ids <- names(band_assignments[band_assignments == 2])
    band13_ids <- names(band_assignments[band_assignments %in% c(1, 3)])

    het_ratio <- mean(indel_result[sample_id %in% band2_ids]$het_hom_ratio,
                       na.rm = TRUE)
    hom_ratio <- mean(indel_result[sample_id %in% band13_ids]$het_hom_ratio,
                       na.rm = TRUE)

    ratio <- if (hom_ratio > 0) het_ratio / hom_ratio else NA_real_
    pass <- !is.na(ratio) && ratio > INDEL_RATIO_THRESHOLD

    return(list(pass = pass, ratio = round(ratio, 3),
                detail = paste0("het_ratio=", round(het_ratio, 3),
                                 ",hom_ratio=", round(hom_ratio, 3))))
  }

  list(pass = NA, ratio = NA_real_, detail = "incompatible_output")
}

# ── Test D: Band stability ────────────────────────────────────────────

test_band_stability <- function(dt, candidate_start, candidate_end,
                                 sample_names) {
  region_dt <- dt[start_bp >= candidate_start & end_bp <= candidate_end]
  if (nrow(region_dt) < MIN_WINDOWS_VIABILITY)
    return(list(pass = NA, stability = NA_real_, detail = "too_few_windows"))

  # Check if per-window band assignments are stored
  if (!"band_assignments" %in% names(attributes(region_dt)) &&
      !"fixed_bands" %in% names(region_dt))
    return(list(pass = NA, stability = NA_real_,
                detail = "no_per_window_bands"))

  # Fallback: use consistency of PC1 ranking across windows
  # Build sample × window PC1 matrix
  if ("PC1_loadings" %in% names(region_dt) &&
      is.list(region_dt$PC1_loadings)) {
    mat <- do.call(rbind, region_dt$PC1_loadings)
    # For each pair of windows, count fraction where sample order is preserved
    n_win <- nrow(mat)
    if (n_win < MIN_WINDOWS_VIABILITY)
      return(list(pass = NA, stability = NA_real_, detail = "too_few_windows"))

    n_agree <- 0L; n_total <- 0L
    # Sample pairs: do they maintain relative position?
    ref_ranks <- rank(mat[1, ])
    for (w in 2:n_win) {
      w_ranks <- rank(mat[w, ])
      cor_val <- cor(ref_ranks, w_ranks, method = "spearman",
                      use = "pairwise.complete.obs")
      if (is.finite(cor_val)) {
        n_total <- n_total + 1L
        if (cor_val > 0.7) n_agree <- n_agree + 1L
      }
    }

    stability <- if (n_total > 0) n_agree / n_total else NA_real_
    pass <- !is.na(stability) && stability >= BAND_STABILITY_FRAC

    return(list(pass = pass, stability = round(stability, 3),
                detail = paste0("agree_windows=", n_agree, "/", n_total)))
  }

  list(pass = NA, stability = NA_real_, detail = "no_loadings")
}

# ── Search mode ────────────────────────────────────────────────────────

search_block_viability <- function(chr, zone_start, zone_end, ...) {
  data.table(method = "block_viability", best_bp = NA_integer_,
             score = 0, is_precise = FALSE,
             detail = "search_not_applicable")
}

# ── Main viability assessment ─────────────────────────────────────────

run_block_viability <- function(chr, candidate_start, candidate_end,
                                 dt, sample_names, band_assignments,
                                 fl = NULL, theta_dir = NULL) {
  message("[cheat25] ", chr, ":", round(candidate_start/1e6,1), "-",
          round(candidate_end/1e6,1), " Mb")

  # Run all four tests
  ta <- test_sv_concordance(chr, candidate_start, candidate_end,
                             band_assignments, fl)
  tb <- test_theta_het(chr, candidate_start, candidate_end,
                        band_assignments, theta_dir)
  tc <- test_indel_slope(chr, candidate_start, candidate_end,
                          band_assignments)
  td <- test_band_stability(dt, candidate_start, candidate_end, sample_names)

  # Count passes (NA = not tested, doesn't count)
  tests <- list(A_sv = ta, B_theta = tb, C_indel = tc, D_stability = td)
  n_pass <- sum(sapply(tests, function(x) isTRUE(x$pass)))
  n_tested <- sum(sapply(tests, function(x) !is.na(x$pass)))

  status <- if (n_pass >= 3) "ALIVE"
    else if (n_pass >= 1) "UNCERTAIN"
    else if (n_tested == 0) "UNTESTED"
    else "DEAD"

  message("[cheat25] Tests: A=", ta$pass, " B=", tb$pass,
          " C=", tc$pass, " D=", td$pass,
          " | Pass: ", n_pass, "/", n_tested,
          " → ", status)

  list(status = status, n_pass = n_pass, n_tested = n_tested,
       tests = tests,
       search_result = search_block_viability(chr, candidate_start,
                        candidate_end))
}

# Alias for wiring guide compatibility
run_cheat25 <- run_block_viability
