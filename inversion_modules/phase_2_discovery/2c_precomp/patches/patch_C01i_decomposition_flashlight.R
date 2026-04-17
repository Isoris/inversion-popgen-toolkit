#!/usr/bin/env Rscript
# =============================================================================
# PATCH: Flashlight integration for STEP_C01i_multi_inversion_decomposition.R
#
# THE KEY CONNECTION — this is where the sv_prior has the biggest impact.
#
# THREE modifications:
#
# 1. SEEDED INITIALIZATION (Cheat 1)
#    If ≥10 anchor samples exist for a candidate, use their SV genotypes
#    as cluster initialization instead of random k-means. Known HOM_REF →
#    cluster 1 centroid, known HET → cluster 2, known HOM_INV → cluster 3.
#    Remaining samples assigned by distance to seed cluster means.
#    Discordant samples (PCA says X, SV says Y) are flagged.
#
# 2. HET-DEL CONSTRAINT (Cheat 2)
#    If a sample has a het-DEL at the inversion breakpoint, it carries
#    the reference arrangement on ≥1 haplotype. Therefore it CANNOT be
#    HOM_INV. This constrains the clustering: any sample with het-DEL at
#    breakpoint is forced into HOM_REF or HET (never HOM_INV).
#
# 3. HEMIZYGOUS SEGMENT INTERPRETATION (Cheat 3)
#    For samples with het-DELs INSIDE the inversion, the deleted region
#    shows only one haplotype. If that haplotype's allele state is pure
#    (no mixing), we get UNAMBIGUOUS haplotype assignment for that segment.
#    This resolves heterozygotes that PCA struggles with at 9× coverage.
#
# INSERT POINTS:
#   - Source sv_prior_loader.R at script top
#   - Replace classify_samples_per_block() with sv_prior-aware version
#   - Add Cheat 3 interpretation after initial classification
# =============================================================================

# ─── Source sv_prior (add at script top) ───────────────────────────

fl_loader <- Sys.getenv("SV_PRIOR_LOADER", "")
if (!nzchar(fl_loader)) {
  for (p in c(
    file.path(dirname(outdir), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path),
    file.path(dirname(dirname(outdir)), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path)
  )) if (file.exists(p)) { fl_loader <- p; break }
}
.decomp_has_sv_prior <- FALSE
if (nzchar(fl_loader) && file.exists(fl_loader)) {
  tryCatch({
    source(fl_loader)
    .decomp_has_sv_prior <- TRUE
    message("[C01i] Flashlight loader sourced — seeded decomposition available")
  }, error = function(e) message("[C01i] Flashlight: ", e$message))
}

# ─── SEEDED k-MEANS INITIALIZATION ──────────────────────────────────
# Uses SV genotypes as cluster seeds instead of random initialization.
# Falls back to standard k-means if <10 anchor samples available.

seeded_kmeans <- function(x, seeds, k = 3L, nstart_fallback = 10L) {
  # x: numeric vector (dosage_mean per sample), named
  # seeds: named list with $HOM_REF, $HET, $HOM_INV sample names
  # Returns: kmeans-compatible object with $cluster, $centers

  seed_ref <- intersect(seeds$HOM_REF, names(x))
  seed_het <- intersect(seeds$HET, names(x))
  seed_inv <- intersect(seeds$HOM_INV, names(x))

  n_seeds <- length(seed_ref) + length(seed_het) + length(seed_inv)

  if (n_seeds < 10 || length(seed_ref) < 3 || length(seed_inv) < 2) {
    # Fallback: standard k-means
    km <- tryCatch(kmeans(x, centers = k, nstart = nstart_fallback),
                   error = function(e) NULL)
    if (!is.null(km)) km$seeded <- FALSE
    return(km)
  }

  # Compute seed centroids
  center_ref <- mean(x[seed_ref], na.rm = TRUE)
  center_het <- if (length(seed_het) >= 3) mean(x[seed_het], na.rm = TRUE) else {
    (center_ref + mean(x[seed_inv], na.rm = TRUE)) / 2
  }
  center_inv <- mean(x[seed_inv], na.rm = TRUE)

  # Order centroids low → high
  centers <- sort(c(center_ref, center_het, center_inv))

  # Assign ALL samples by distance to seed centroids
  dists <- matrix(0, length(x), k)
  for (ci in seq_len(k)) dists[, ci] <- abs(x - centers[ci])
  cluster <- apply(dists, 1, which.min)
  names(cluster) <- names(x)

  # Refine with one round of k-means (seeded start, not random)
  km <- tryCatch(
    kmeans(x, centers = matrix(centers, ncol = 1), iter.max = 20),
    error = function(e) NULL
  )

  if (!is.null(km)) {
    km$seeded <- TRUE
    km$n_seed_ref <- length(seed_ref)
    km$n_seed_het <- length(seed_het)
    km$n_seed_inv <- length(seed_inv)
  } else {
    # Manual fallback
    km <- list(
      cluster = cluster,
      centers = matrix(centers, ncol = 1),
      size = as.integer(table(factor(cluster, levels = 1:k))),
      seeded = TRUE,
      n_seed_ref = length(seed_ref),
      n_seed_het = length(seed_het),
      n_seed_inv = length(seed_inv)
    )
  }
  km
}

# ─── FLASHLIGHT-AWARE SAMPLE CLASSIFICATION ─────────────────────────
# Replaces the unsupervised classify_samples_per_block() with one that:
# 1. Checks sv_prior for anchor genotypes (Cheat 1)
# 2. Uses seeded k-means if enough anchors exist
# 3. Applies Cheat 2 constraint (het-DEL → not HOM_INV)
# 4. Flags discordant samples

classify_samples_sv_prior <- function(gt, block, samples, chr, candidate_start, candidate_end) {
  # First: standard dosage computation (same as original)
  b_gt <- gt[block$marker_idx, , drop = FALSE]
  ns <- ncol(b_gt)
  n0 <- colSums(b_gt == 0, na.rm = TRUE)
  n1 <- colSums(b_gt == 1, na.rm = TRUE)
  n2 <- colSums(b_gt == 2, na.rm = TRUE)
  n_valid <- n0 + n1 + n2
  het_rate <- n1 / pmax(n_valid, 1)
  dosage_mean <- (n1 + 2 * n2) / pmax(n_valid, 1)
  names(dosage_mean) <- samples

  # ── Cheat 1: Query sv_prior for anchor samples ──
  seeds <- list(HOM_REF = character(0), HET = character(0), HOM_INV = character(0))
  cheat2_constrained <- character(0)
  sv_prior_mode <- "unsupervised"

  if (.decomp_has_sv_prior) {
    fl <- load_sv_prior(chr)
    if (!is.null(fl)) {
      # Get SV anchors overlapping this candidate region
      anchors <- get_sv_anchors(chr, candidate_start, candidate_end,
                                min_confidence = "MEDIUM")

      if (nrow(anchors) >= 10) {
        seeds$HOM_REF <- unique(anchors[sv_genotype == "HOM_REF"]$sample_id)
        seeds$HET     <- unique(anchors[sv_genotype == "HET"]$sample_id)
        seeds$HOM_INV <- unique(anchors[sv_genotype == "HOM_INV"]$sample_id)
        sv_prior_mode <- "seeded"
        message("    [sv_prior] Seeded mode: REF=", length(seeds$HOM_REF),
                " HET=", length(seeds$HET), " INV=", length(seeds$HOM_INV))
      }

      # ── Cheat 2: Het-DEL constraint ──
      # Samples with het-DEL at breakpoint → cannot be HOM_INV
      bp_dels_left  <- get_breakpoint_dels(chr, candidate_start, window = 50000L)
      bp_dels_right <- get_breakpoint_dels(chr, candidate_end, window = 50000L)
      all_bp_del_carriers <- unique(c(
        unlist(bp_dels_left$het_carriers),
        unlist(bp_dels_right$het_carriers)
      ))
      cheat2_constrained <- intersect(all_bp_del_carriers, samples)
      if (length(cheat2_constrained) > 0) {
        message("    [sv_prior] Cheat 2: ", length(cheat2_constrained),
                " samples constrained (het-DEL at breakpoint → not HOM_INV)")
      }
    }
  }

  # ── Clustering ──
  if (sv_prior_mode == "seeded") {
    km <- seeded_kmeans(dosage_mean, seeds, k = 3L)
  } else {
    km <- tryCatch(kmeans(dosage_mean, centers = 3, nstart = 10),
                   error = function(e) NULL)
  }

  if (!is.null(km)) {
    co <- order(km$centers[, 1])
    group <- character(ns)
    group[km$cluster == co[1]] <- "HOMO_REF"
    group[km$cluster == co[2]] <- "HET"
    group[km$cluster == co[3]] <- "HOMO_INV"
  } else {
    group <- fifelse(dosage_mean < 0.5, "HOMO_REF",
             fifelse(dosage_mean > 1.5, "HOMO_INV", "HET"))
  }

  # ── Apply Cheat 2 constraint ──
  # Force het-DEL-at-breakpoint samples out of HOM_INV
  n_cheat2_reclassified <- 0L
  if (length(cheat2_constrained) > 0) {
    for (si in seq_len(ns)) {
      if (samples[si] %in% cheat2_constrained && group[si] == "HOMO_INV") {
        # Reclassify as HET (het-DEL means at least one REF haplotype)
        group[si] <- "HET"
        n_cheat2_reclassified <- n_cheat2_reclassified + 1L
      }
    }
    if (n_cheat2_reclassified > 0) {
      message("    [sv_prior] Cheat 2 reclassified ", n_cheat2_reclassified,
              " samples: HOM_INV → HET (het-DEL at breakpoint)")
    }
  }

  # ── Flag discordant samples (SV says X, clustering says Y) ──
  discordant <- character(0)
  if (sv_prior_mode == "seeded") {
    for (si in seq_len(ns)) {
      sid <- samples[si]
      sv_gt <- NULL
      if (sid %in% seeds$HOM_REF) sv_gt <- "HOMO_REF"
      else if (sid %in% seeds$HET) sv_gt <- "HET"
      else if (sid %in% seeds$HOM_INV) sv_gt <- "HOMO_INV"

      if (!is.null(sv_gt) && sv_gt != group[si]) {
        discordant <- c(discordant, sid)
      }
    }
    if (length(discordant) > 0) {
      message("    [sv_prior] Discordant: ", length(discordant),
              " samples (SV ≠ clustering)")
    }
  }

  # Frequency
  freq <- (sum(group %in% c("HOMO_INV", "RARE_INV")) * 2 +
            sum(group %in% c("HET", "RARE_HET"))) / (2 * ns)

  # Build output
  profiles <- character(ns)
  for (si in seq_len(ns)) profiles[si] <- encode_012(b_gt[, si])

  out_dt <- data.table(
    sample = samples, group = group,
    dosage_mean = round(dosage_mean, 3), het_rate = round(het_rate, 3),
    n_valid = n_valid, profile = profiles,
    frequency = round(freq, 4),
    sv_prior_mode = sv_prior_mode,
    is_discordant = samples %in% discordant,
    is_cheat2_constrained = samples %in% cheat2_constrained
  )

  out_dt
}

# ─── CHEAT 3: HEMIZYGOUS SEGMENT INTERPRETATION ─────────────────────
# After initial classification, for samples with het-DELs INSIDE the
# inversion, check the allele state in the deleted region.
# The deleted region shows only one haplotype = UNAMBIGUOUS.

interpret_hemizygous_segments <- function(gt, pos, samples, sample_class,
                                           chr, candidate_start, candidate_end) {
  if (!.decomp_has_sv_prior) return(data.table())

  int_dels <- get_internal_dels(chr, candidate_start, candidate_end)
  if (nrow(int_dels) == 0) return(data.table())

  hemi_rows <- list()

  for (di in seq_len(nrow(int_dels))) {
    del_start <- int_dels$del_start[di]
    del_end   <- int_dels$del_end[di]
    carriers  <- int_dels$het_carriers[[di]]
    carriers  <- intersect(carriers, samples)
    if (length(carriers) == 0) next

    # Find markers in the deleted segment
    del_marker_idx <- which(pos >= del_start & pos <= del_end)
    if (length(del_marker_idx) < 5) next

    del_gt <- gt[del_marker_idx, , drop = FALSE]

    for (sid in carriers) {
      si <- match(sid, samples)
      if (is.na(si)) next

      # In the deleted segment, this sample is hemizygous
      # Reads come from only one haplotype (~4.5× of one, ~0× of other)
      # The genotype should be nearly pure: all 0 or all 2
      s_gt <- del_gt[, si]
      s_gt <- s_gt[!is.na(s_gt)]
      if (length(s_gt) < 5) next

      frac_0 <- mean(s_gt == 0)
      frac_2 <- mean(s_gt == 2)
      frac_1 <- mean(s_gt == 1)  # should be very low if truly hemizygous

      # Determine surviving haplotype
      if (frac_0 > 0.7) {
        surviving <- "REF"  # reference haplotype survived
      } else if (frac_2 > 0.7) {
        surviving <- "INV"  # inverted haplotype survived
      } else {
        surviving <- "AMBIGUOUS"
      }

      # Check purity: low het rate = truly hemizygous = high confidence
      purity <- 1 - frac_1
      confidence <- if (purity > 0.85 && frac_1 < 0.10) "HIGH"
                    else if (purity > 0.70) "MODERATE"
                    else "LOW"

      hemi_rows[[length(hemi_rows) + 1]] <- data.table(
        sample_id = sid,
        del_id = int_dels$del_id[di],
        del_start = del_start, del_end = del_end,
        del_svlen = int_dels$del_svlen[di],
        n_markers_in_del = length(s_gt),
        frac_hom_ref = round(frac_0, 3),
        frac_het = round(frac_1, 3),
        frac_hom_inv = round(frac_2, 3),
        purity = round(purity, 3),
        surviving_haplotype = surviving,
        confidence = confidence,
        overall_group = sample_class[sample == sid]$group[1]
      )
    }
  }

  hemi_dt <- if (length(hemi_rows) > 0) rbindlist(hemi_rows) else data.table()

  if (nrow(hemi_dt) > 0) {
    n_high <- sum(hemi_dt$confidence == "HIGH")
    n_ref  <- sum(hemi_dt$surviving_haplotype == "REF")
    n_inv  <- sum(hemi_dt$surviving_haplotype == "INV")
    message("    [sv_prior] Cheat 3: ", nrow(hemi_dt), " hemizygous segments",
            " (", n_high, " HIGH conf, REF=", n_ref, " INV=", n_inv, ")")
  }

  hemi_dt
}

# ─── USAGE IN MAIN LOOP ─────────────────────────────────────────────
# Replace the call to classify_samples_per_block() with:
#
#   samp_class <- classify_samples_sv_prior(
#     gt_data$gt, blk, gt_data$samples,
#     chr, cand$start_mb * 1e6, cand$end_mb * 1e6
#   )
#   samp_class[, `:=`(chrom = chr, interval_id = iid, system_id = bi)]
#
# After classification, add Cheat 3 interpretation:
#
#   hemi_dt <- interpret_hemizygous_segments(
#     gt_data$gt, gt_data$pos, gt_data$samples, samp_class,
#     chr, cand$start_mb * 1e6, cand$end_mb * 1e6
#   )
#   if (nrow(hemi_dt) > 0) {
#     hemi_dt[, `:=`(chrom = chr, interval_id = iid, system_id = bi)]
#     all_hemi[[length(all_hemi) + 1]] <- hemi_dt
#   }
#
# At script end, write hemizygous segments:
#
#   hemi_all <- if (length(all_hemi) > 0) rbindlist(all_hemi, fill = TRUE) else data.table()
#   fwrite(hemi_all, file.path(outdir, "hemizygous_segments.tsv.gz"), sep = "\t")
