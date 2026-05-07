#!/usr/bin/env Rscript
# =============================================================================
# PATCH: Flashlight integration for STEP_C04_snake3_ghsl_v3.R
#
# After GHSL computes within-band vs between-band concordance, check:
#   Do DELLY-confirmed HOM_INV samples show higher within-band concordance
#   than DELLY-confirmed HOM_REF samples?
#
# This is a DIRECT test that phased haplotypes agree with SV genotypes.
# If GHSL says "band 3 samples share haplotypes" AND DELLY says "band 3
# samples are HOM_INV" → strong independent confirmation.
#
# Also: for each SV anchor, compare its GHSL assignment_confidence with
# its SV genotype. Concordant = both agree. Discordant = PCA/GHSL and
# SV callers disagree, which could mean:
#   - SV genotype error (~10-15% at 9×)
#   - GHSL band assignment error (phase noise)
#   - Recombinant sample (true biological discordance)
#
# INSERT: After process_window() returns sample_scores, add sv_prior
# cross-check as a post-processing step per chromosome.
# =============================================================================

# ─── Source sv_prior ───────────────────────────────────────────────

fl_loader <- Sys.getenv("SV_PRIOR_LOADER", "")
if (!nzchar(fl_loader)) {
  for (p in c(
    file.path(dirname(outdir), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path),
    file.path(dirname(dirname(outdir)), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path)
  )) if (file.exists(p)) { fl_loader <- p; break }
}
.ghsl_has_sv_prior <- FALSE
if (nzchar(fl_loader) && file.exists(fl_loader)) {
  tryCatch({
    source(fl_loader)
    .ghsl_has_sv_prior <- TRUE
    message("[S3v3] Flashlight sourced — SV anchor cross-check available")
  }, error = function(e) message("[S3v3] Flashlight: ", e$message))
}

# ─── GHSL × FLASHLIGHT CROSS-CHECK ──────────────────────────────────
# Called once per chromosome, after all windows are processed.
# Takes the aggregated sample_scores and cross-references with SV anchors.

ghsl_sv_prior_crosscheck <- function(chr_scores_dt, chr, bands) {
  # chr_scores_dt: data.table from all_scores for this chromosome
  #   columns: sample_id, pc1_band, within_mean, between_mean, sample_ghsl, ...
  # bands: named integer vector (sample → band) from chromosome-wide k=3

  empty_result <- list(
    crosscheck = data.table(),
    summary = data.table(
      chrom = chr, n_anchors_tested = 0L, n_concordant = 0L,
      n_discordant = 0L, mean_ghsl_hom_inv = NA_real_,
      mean_ghsl_hom_ref = NA_real_, ghsl_sv_agreement = NA_real_
    )
  )

  if (!.ghsl_has_sv_prior || nrow(chr_scores_dt) == 0) return(empty_result)

  fl <- load_sv_prior(chr)
  if (is.null(fl) || nrow(fl$inv_calls) == 0) return(empty_result)

  # Get all SV anchors for this chromosome
  chr_start <- 0; chr_end <- max(fl$inv_calls$bp2, na.rm = TRUE) + 1e6
  anchors <- get_sv_anchors(chr, chr_start, chr_end, min_confidence = "LOW")
  if (nrow(anchors) < 5) return(empty_result)

  # Aggregate GHSL scores per sample (mean across all windows)
  if (!"sample_id" %in% names(chr_scores_dt)) return(empty_result)

  ghsl_per_sample <- chr_scores_dt[
    !is.na(sample_ghsl),
    .(mean_ghsl = mean(sample_ghsl, na.rm = TRUE),
      mean_within = mean(within_mean, na.rm = TRUE),
      mean_between = mean(between_mean, na.rm = TRUE),
      n_windows = .N,
      n_confirms_pc1 = sum(ghsl_confirms_pc1 == TRUE, na.rm = TRUE)),
    by = sample_id
  ]

  # Cross-reference with SV anchors
  anchors_unique <- anchors[!duplicated(sample_id)]
  merged <- merge(ghsl_per_sample, anchors_unique,
                  by = "sample_id", all = FALSE)

  if (nrow(merged) < 5) return(empty_result)

  # Add PCA band
  if (!is.null(bands) && length(bands) > 0) {
    merged[, pca_band := bands[sample_id]]
  }

  # Check concordance: SV genotype vs PCA band
  # Expected: HOM_REF → band 1, HET → band 2, HOM_INV → band 3
  expected_fwd <- c(HOM_REF = 1L, HET = 2L, HOM_INV = 3L)
  expected_rev <- c(HOM_REF = 3L, HET = 2L, HOM_INV = 1L)

  merged[, exp_fwd := expected_fwd[sv_genotype]]
  merged[, exp_rev := expected_rev[sv_genotype]]
  n_fwd <- sum(merged$pca_band == merged$exp_fwd, na.rm = TRUE)
  n_rev <- sum(merged$pca_band == merged$exp_rev, na.rm = TRUE)
  orientation <- if (n_fwd >= n_rev) "forward" else "reversed"

  if (orientation == "forward") {
    merged[, sv_pca_concordant := (pca_band == exp_fwd)]
  } else {
    merged[, sv_pca_concordant := (pca_band == exp_rev)]
  }

  # Key comparison: do HOM_INV anchors have higher within-band GHSL
  # than HOM_REF anchors?
  ghsl_hom_inv <- merged[sv_genotype == "HOM_INV"]$mean_ghsl
  ghsl_hom_ref <- merged[sv_genotype == "HOM_REF"]$mean_ghsl
  ghsl_het     <- merged[sv_genotype == "HET"]$mean_ghsl

  mean_ghsl_inv <- if (length(ghsl_hom_inv) > 0) mean(ghsl_hom_inv, na.rm = TRUE) else NA_real_
  mean_ghsl_ref <- if (length(ghsl_hom_ref) > 0) mean(ghsl_hom_ref, na.rm = TRUE) else NA_real_
  mean_ghsl_het <- if (length(ghsl_het) > 0) mean(ghsl_het, na.rm = TRUE) else NA_real_

  # Wilcoxon test: HOM_INV within-band GHSL > HOM_REF within-band GHSL?
  wt_p <- NA_real_
  if (length(ghsl_hom_inv) >= 3 && length(ghsl_hom_ref) >= 3) {
    wt <- tryCatch(
      wilcox.test(ghsl_hom_inv, ghsl_hom_ref, alternative = "greater"),
      error = function(e) NULL
    )
    if (!is.null(wt)) wt_p <- wt$p.value
  }

  # GHSL-SV agreement rate
  n_concordant <- sum(merged$sv_pca_concordant, na.rm = TRUE)
  n_total <- sum(!is.na(merged$sv_pca_concordant))
  agreement <- if (n_total > 0) n_concordant / n_total else NA_real_

  # Flag discordant samples (candidates for recombination detection)
  discordant_samples <- merged[sv_pca_concordant == FALSE]$sample_id

  # Per-anchor detail table
  crosscheck <- merged[, .(
    sample_id, sv_genotype, sv_confidence, pca_band,
    mean_ghsl, mean_within, mean_between,
    n_windows, n_confirms_pc1,
    sv_pca_concordant
  )]

  summary_dt <- data.table(
    chrom = chr,
    n_anchors_tested = nrow(merged),
    n_concordant = n_concordant,
    n_discordant = sum(!merged$sv_pca_concordant, na.rm = TRUE),
    mean_ghsl_hom_inv = round(mean_ghsl_inv, 4),
    mean_ghsl_hom_ref = round(mean_ghsl_ref, 4),
    mean_ghsl_het = round(mean_ghsl_het, 4),
    ghsl_sv_agreement = round(agreement, 3),
    wilcox_p = round(wt_p, 6),
    orientation = orientation,
    n_discordant_samples = length(discordant_samples)
  )

  if (agreement > 0.7 && !is.na(wt_p) && wt_p < 0.05) {
    message("[S3v3] [sv_prior] ", chr, ": GHSL confirms SV genotypes!",
            " agreement=", round(agreement, 2),
            " Wilcoxon p=", format(wt_p, digits = 3),
            " (HOM_INV GHSL=", round(mean_ghsl_inv, 3),
            " > HOM_REF GHSL=", round(mean_ghsl_ref, 3), ")")
  } else if (!is.na(agreement) && agreement < 0.4) {
    message("[S3v3] [sv_prior] ", chr, ": GHSL discordant with SV",
            " agreement=", round(agreement, 2))
  }

  list(crosscheck = crosscheck, summary = summary_dt)
}

# ─── USAGE: After the main per-chromosome loop ──────────────────────
# At the end of the `for (chr in chroms)` loop, after collecting
# chr_scores, add:
#
#   # Flashlight cross-check
#   if (.ghsl_has_sv_prior && length(chr_scores) > 0) {
#     chr_scores_dt <- rbindlist(chr_scores, fill = TRUE)
#     fl_check <- ghsl_sv_prior_crosscheck(chr_scores_dt, chr, bands)
#     if (nrow(fl_check$crosscheck) > 0) {
#       fl_check$crosscheck[, chrom := chr]
#       all_fl_crosscheck <- c(all_fl_crosscheck, list(fl_check$crosscheck))
#     }
#     all_fl_summary <- c(all_fl_summary, list(fl_check$summary))
#   }
#
# At script end, write sv_prior cross-check outputs:
#
#   fl_cross_dt <- if (length(all_fl_crosscheck) > 0) rbindlist(all_fl_crosscheck, fill = TRUE) else data.table()
#   fl_summ_dt  <- if (length(all_fl_summary) > 0) rbindlist(all_fl_summary, fill = TRUE) else data.table()
#   fwrite(fl_cross_dt, file.path(outdir, "snake3v3_sv_prior_crosscheck.tsv.gz"), sep = "\t")
#   fwrite(fl_summ_dt, file.path(outdir, "snake3v3_sv_prior_summary.tsv"), sep = "\t")
