#!/usr/bin/env Rscript
# =============================================================================
# cheat30_gds_by_genotype.R — Pairwise GDS (Genotype Dosage Similarity) stratified by inversion genotype
#
# Based on Porubsky et al. 2022 Cell, Figure 3C-D.
#
# THE KEY INSIGHT:
#   For a REAL inversion, pairwise GDS (identity-by-state) distributions
#   DIFFER between genotype pairs:
#     I/I (inv/inv):   high GDS (shared inverted haplotype)
#     D/D (dir/dir):   high GDS (shared direct haplotype)
#     D/I (het pairs): LOWER GDS (different arrangements diverged)
#
#   For a single-origin inversion: I/I GDS is very high (one founder)
#   For a recurrent inversion: I/I GDS shows a BIMODAL distribution
#     (multiple independent origins → not all carriers share a haplotype)
#
#   This is Porubsky's Fig 3C (8p23.1) vs 3D (17q21.31):
#     8p23.1: D/D bimodal → multiple direct haplotype backgrounds
#     17q21.31: clean separation → ancient single-origin inversion
#
# FOR CATFISH:
#   We have decomposition genotypes (HOM_REF, HET, HOM_INV) from C01i.
#   We have dosage matrices from BEAGLE.
#   We can compute pairwise GDS within the inversion region and stratify
#   by genotype pair. The distribution shape tells us:
#     - Is the inversion real? (D/I lower than D/D and I/I)
#     - Is it single-origin? (I/I unimodal, tight)
#     - Is it recurrent? (I/I bimodal or wide)
#     - How old is it? (wider D/I gap = older divergence)
#
# INPUT:
#   dosage_matrix:    from BEAGLE (sites × samples) for the candidate region
#   decomp_classes:   named vector (sample_id → "HOM_REF"/"HET"/"HOM_INV")
#   candidate region: chr, start_bp, end_bp
#
# OUTPUT:
#   gds_summary:      per-genotype-pair mean, sd, median GDS
#   separation_test:  Wilcoxon D/I vs (D/D + I/I), p-value
#   bimodality_test:  Hartigan's dip test on I/I distribution
#   origin_class:     "single_origin" / "recurrent" / "inconclusive"
#   age_proxy:        D/I vs D/D GDS gap (larger = older)
#
# INTEGRATION:
#   C01f: as T14_ibs_genotype test (independent confirmation of inversion)
#   Cheat 15: bimodality of I/I feeds recurrence classification
#   Cheat 29: mechanism + GDS = full origin model
#
# REQUIRES: BEAGLE dosage (from dispatcher cache or direct), decomposition
# DISPATCHER: No (pure R computation on dosage matrix)
# DIFFICULTY: Medium
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
MIN_SAMPLES_PER_CLASS <- 5L
MIN_SITES             <- 50L
IBS_SUBSAMPLE_SITES   <- 5000L   # subsample sites for speed
DIP_PVAL_THRESHOLD    <- 0.05    # Hartigan's dip test for bimodality

# ── Compute pairwise GDS from dosage matrix ───────────────────────────

compute_pairwise_gds <- function(dosage_mat) {
# Backward-compatible alias
compute_pairwise_ibs <- compute_pairwise_gds
  # dosage_mat: sites × samples (values 0-2, expected dosage)
  n <- ncol(dosage_mat)
  ns <- nrow(dosage_mat)
  if (n < 2 || ns < MIN_SITES) return(NULL)

  # Subsample sites if too many
  if (ns > IBS_SUBSAMPLE_SITES) {
    idx <- sort(sample(ns, IBS_SUBSAMPLE_SITES))
    dosage_mat <- dosage_mat[idx, , drop = FALSE]
    ns <- nrow(dosage_mat)
  }

  # GDS = 1 - |d_i - d_j| / 2, averaged across sites
  ibs <- matrix(NA_real_, n, n, dimnames = list(colnames(dosage_mat),
                                                  colnames(dosage_mat)))
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      diff <- abs(dosage_mat[, i] - dosage_mat[, j])
      ok <- is.finite(diff)
      if (sum(ok) >= MIN_SITES) {
        ibs[i, j] <- 1 - mean(diff[ok]) / 2
        ibs[j, i] <- ibs[i, j]
      }
    }
  }
  diag(ibs) <- 1
  ibs
}

# ── Stratify GDS by genotype pairs ────────────────────────────────────

stratify_ibs_by_genotype <- function(ibs_mat, decomp_classes) {
  if (is.null(ibs_mat)) return(NULL)

  samples <- colnames(ibs_mat)
  classes <- decomp_classes[samples]
  classes <- classes[!is.na(classes)]
  samples <- names(classes)

  if (length(samples) < 10) return(NULL)

  # Build pair table
  pairs <- list()
  for (i in 1:(length(samples)-1)) {
    for (j in (i+1):length(samples)) {
      si <- samples[i]; sj <- samples[j]
      ci <- classes[si]; cj <- classes[sj]
      ibs_val <- ibs_mat[si, sj]
      if (!is.finite(ibs_val)) next

      # Pair type (order-independent)
      pair_type <- paste(sort(c(ci, cj)), collapse = "_")
      pairs[[length(pairs) + 1]] <- list(
        sample_a = si, sample_b = sj,
        class_a = ci, class_b = cj,
        pair_type = pair_type, ibs = ibs_val
      )
    }
  }

  if (length(pairs) == 0) return(NULL)
  rbindlist(pairs)
}

# ── Run separation test ───────────────────────────────────────────────

test_genotype_separation <- function(pair_dt) {
  if (is.null(pair_dt) || nrow(pair_dt) == 0)
    return(list(separation_p = NA, separation_effect = NA,
                age_proxy = NA))

  # D/I (heterozygous pairs) vs same-genotype pairs
  same_geno <- pair_dt[pair_type %in% c("HOM_REF_HOM_REF", "HOM_INV_HOM_INV")]
  diff_geno <- pair_dt[pair_type == "HOM_INV_HOM_REF" |
                         pair_type == "HOM_REF_HOM_INV"]

  if (nrow(same_geno) < 5 || nrow(diff_geno) < 5)
    return(list(separation_p = NA, separation_effect = NA, age_proxy = NA))

  # Wilcoxon: same-geno GDS should be HIGHER than diff-geno GDS
  wt <- tryCatch(
    wilcox.test(same_geno$ibs, diff_geno$ibs, alternative = "greater"),
    error = function(e) list(p.value = NA)
  )

  effect <- mean(same_geno$ibs, na.rm = TRUE) - mean(diff_geno$ibs, na.rm = TRUE)

  # Age proxy: bigger gap = older divergence between arrangements
  ref_ref <- pair_dt[pair_type == "HOM_REF_HOM_REF"]
  inv_inv <- pair_dt[pair_type == "HOM_INV_HOM_INV"]
  age_proxy <- NA_real_
  if (nrow(ref_ref) > 0 && nrow(diff_geno) > 0) {
    age_proxy <- mean(ref_ref$ibs, na.rm = TRUE) - mean(diff_geno$ibs, na.rm = TRUE)
  }

  list(separation_p = wt$p.value,
       separation_effect = round(effect, 4),
       age_proxy = round(age_proxy, 4))
}

# ── Test bimodality of I/I distribution (recurrence signal) ───────────

test_ibs_bimodality <- function(pair_dt) {
  inv_inv <- pair_dt[pair_type == "HOM_INV_HOM_INV"]
  if (nrow(inv_inv) < 10)
    return(list(dip_stat = NA, dip_p = NA, is_bimodal = NA))

  # Hartigan's dip test
  dip_result <- tryCatch({
    if (requireNamespace("diptest", quietly = TRUE)) {
      d <- diptest::dip.test(inv_inv$ibs)
      list(dip_stat = round(d$statistic, 4),
           dip_p = d$p.value,
           is_bimodal = d$p.value < DIP_PVAL_THRESHOLD)
    } else {
      # Fallback: use coefficient of variation as bimodality proxy
      cv <- sd(inv_inv$ibs) / mean(inv_inv$ibs)
      list(dip_stat = round(cv, 4), dip_p = NA,
           is_bimodal = cv > 0.15)  # high CV suggests multimodality
    }
  }, error = function(e) list(dip_stat = NA, dip_p = NA, is_bimodal = NA))

  dip_result
}

# ── Runner ────────────────────────────────────────────────────────────

run_cheat30 <- function(chr, start_bp, end_bp,
                         dosage_matrix, decomp_classes) {
  message("[cheat30] ", chr, ":", round(start_bp/1e6, 1), "-",
          round(end_bp/1e6, 1), " Mb")

  default <- list(
    separation_p = NA, separation_effect = NA, age_proxy = NA,
    dip_stat = NA, dip_p = NA, is_bimodal = NA,
    origin_class = "inconclusive",
    n_ref = 0L, n_het = 0L, n_inv = 0L,
    mean_ibs_same = NA, mean_ibs_diff = NA
  )

  if (is.null(dosage_matrix) || is.null(decomp_classes))
    return(default)

  # Check sample counts per class
  n_ref <- sum(decomp_classes == "HOM_REF")
  n_het <- sum(decomp_classes == "HET")
  n_inv <- sum(decomp_classes == "HOM_INV")

  if (n_ref < MIN_SAMPLES_PER_CLASS || n_inv < MIN_SAMPLES_PER_CLASS) {
    message("[cheat30] Too few samples per class (REF=", n_ref,
            " HET=", n_het, " INV=", n_inv, ")")
    return(modifyList(default, list(n_ref = n_ref, n_het = n_het, n_inv = n_inv)))
  }

  # Compute GDS
  ibs <- compute_pairwise_ibs(dosage_matrix)
  if (is.null(ibs)) {
    message("[cheat30] GDS computation failed")
    return(modifyList(default, list(n_ref = n_ref, n_het = n_het, n_inv = n_inv)))
  }

  # Stratify
  pair_dt <- stratify_ibs_by_genotype(ibs, decomp_classes)
  if (is.null(pair_dt) || nrow(pair_dt) < 20) {
    message("[cheat30] Too few pairs")
    return(modifyList(default, list(n_ref = n_ref, n_het = n_het, n_inv = n_inv)))
  }

  # Tests
  sep <- test_genotype_separation(pair_dt)
  bim <- test_ibs_bimodality(pair_dt)

  # Classification
  origin_class <- "inconclusive"
  if (!is.na(sep$separation_p) && sep$separation_p < 0.001 &&
      sep$separation_effect > 0.02) {
    # Real inversion (significant genotype-GDS association)
    if (!is.na(bim$is_bimodal) && bim$is_bimodal) {
      origin_class <- "recurrent"  # I/I bimodal → multiple origins
    } else {
      origin_class <- "single_origin"  # I/I unimodal → one founder
    }
  } else if (!is.na(sep$separation_p) && sep$separation_p < 0.05) {
    origin_class <- "weak_signal"
  }

  same_geno <- pair_dt[pair_type %in% c("HOM_REF_HOM_REF", "HOM_INV_HOM_INV")]
  diff_geno <- pair_dt[grepl("HOM_REF.*HOM_INV|HOM_INV.*HOM_REF", pair_type)]

  message("[cheat30] Separation p=", signif(sep$separation_p, 3),
          " effect=", sep$separation_effect,
          " | Bimodal=", bim$is_bimodal,
          " | Origin: ", origin_class)

  list(
    separation_p = sep$separation_p,
    separation_effect = sep$separation_effect,
    age_proxy = sep$age_proxy,
    dip_stat = bim$dip_stat,
    dip_p = bim$dip_p,
    is_bimodal = bim$is_bimodal,
    origin_class = origin_class,
    n_ref = n_ref, n_het = n_het, n_inv = n_inv,
    mean_ibs_same = round(mean(same_geno$ibs, na.rm = TRUE), 4),
    mean_ibs_diff = round(mean(diff_geno$ibs, na.rm = TRUE), 4),
    pair_table = pair_dt  # for diagnostic plots
  )
}
