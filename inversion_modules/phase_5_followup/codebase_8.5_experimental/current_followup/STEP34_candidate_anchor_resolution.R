#!/usr/bin/env Rscript

# =============================================================================
# STEP34_candidate_anchor_resolution.R
#
# Anchor-based inversion-state resolution.
#
# If one homokaryotype stripe is sufficiently clean (compact, coherent), use
# it as an anchor to decompose the HET band into directional subgroups.
#
# For each candidate:
#   1. Evaluate anchor eligibility per homo stripe
#   2. Define anchor markers (stable in the clean homo stripe)
#   3. For each HET sample: compute anchor-side composition
#   4. Classify HET samples by AA-side vs BB-side similarity
#   5. Detect AB subgroups and directional relationships
#   6. Export LD/FST-ready contrast group manifests
#
# Inputs:
#   - STEP21 rotated PCA
#   - STEP29 coherence scores + marker polarity
#   - Dosage matrix
#
# Outputs:
#   - candidate_anchor_eligibility.tsv
#   - candidate_het_anchor_composition.tsv
#   - candidate_contrast_groups.tsv           (for downstream LD/FST/Hobs)
#   - candidate_contrast_manifests/           (per-contrast sample lists)
#
# Usage:
#   Rscript STEP34_candidate_anchor_resolution.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  coh_file <- file.path(cand_dir, "candidate_sample_coherence.tsv")
  if (!file.exists(rot_file)) next

  rot <- fread(rot_file)
  coh <- if (file.exists(coh_file)) fread(coh_file) else NULL

  # Load dosage
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) next

  dos <- fread(dos_file)
  sites <- fread(sites_file)
  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 20) next

  dos_reg <- dos[keep]
  sc <- setdiff(names(dos_reg), "marker")
  if (all(grepl("^Ind", sc)) && length(sc) == nrow(rot)) {
    setnames(dos_reg, old = sc, new = rot$sample)
  }
  common <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))
  if (length(common) < 10) next

  X <- as.matrix(dos_reg[, ..common])
  storage.mode(X) <- "double"

  message("[INFO] Candidate ", cid, ": anchor resolution")

  # ── Evaluate anchor eligibility ────────────────────────────────────────
  anchor_results <- list()
  for (g in c("HOMO_1", "HOMO_2")) {
    g_samples <- rot[coarse_group_refined == g, sample]
    g_samples <- intersect(g_samples, common)
    ng <- length(g_samples)

    if (ng < 5) {
      anchor_results[[g]] <- data.table(
        candidate_id = cid, stripe = g, n_samples = ng,
        eligible = FALSE, reason = "too_few_samples"
      )
      next
    }

    # Core samples only if coherence available
    if (!is.null(coh)) {
      core_g <- coh[coarse_group == g & stripe_quality == "core", sample]
      core_g <- intersect(core_g, common)
      n_core <- length(core_g)
    } else {
      core_g <- g_samples
      n_core <- ng
    }

    # Compactness: coefficient of variation of PCA distances
    g_rot <- rot[sample %in% core_g]
    if (nrow(g_rot) < 3) {
      anchor_results[[g]] <- data.table(
        candidate_id = cid, stripe = g, n_samples = ng,
        n_core = n_core, eligible = FALSE, reason = "too_few_core"
      )
      next
    }

    dists <- sqrt((g_rot$u - mean(g_rot$u))^2 + (g_rot$v - mean(g_rot$v))^2)
    cv <- sd(dists) / mean(dists)
    compact <- cv < 1.5

    # Coherence: fraction coherent
    if (!is.null(coh)) {
      frac_coherent <- mean(coh[coarse_group == g, coherence_class] == "coherent")
    } else {
      frac_coherent <- NA_real_
    }

    eligible <- compact && n_core >= 5 && (is.na(frac_coherent) || frac_coherent >= 0.5)
    reason <- if (eligible) "eligible" else
              if (!compact) "not_compact" else
              if (n_core < 5) "core_too_small" else "low_coherence"

    anchor_results[[g]] <- data.table(
      candidate_id = cid, stripe = g, n_samples = ng,
      n_core = n_core, cv_distance = round(cv, 4),
      frac_coherent = round(frac_coherent, 4),
      eligible = eligible, reason = reason
    )
  }

  anch_dt <- rbindlist(anchor_results, fill = TRUE)
  fwrite(anch_dt, file.path(cand_dir, "candidate_anchor_eligibility.tsv"), sep = "\t")

  # ── Anchor marker identification ───────────────────────────────────────
  h1_core <- if (anch_dt[stripe == "HOMO_1", eligible]) {
    if (!is.null(coh)) coh[coarse_group == "HOMO_1" & stripe_quality == "core", sample]
    else rot[coarse_group_refined == "HOMO_1", sample]
  } else NULL
  h2_core <- if (anch_dt[stripe == "HOMO_2", eligible]) {
    if (!is.null(coh)) coh[coarse_group == "HOMO_2" & stripe_quality == "core", sample]
    else rot[coarse_group_refined == "HOMO_2", sample]
  } else NULL

  h1_core <- intersect(h1_core, common)
  h2_core <- intersect(h2_core, common)

  has_both_anchors <- length(h1_core) >= 3 && length(h2_core) >= 3

  if (has_both_anchors) {
    h1_mean <- rowMeans(X[, h1_core, drop = FALSE], na.rm = TRUE)
    h2_mean <- rowMeans(X[, h2_core, drop = FALSE], na.rm = TRUE)
    delta <- h2_mean - h1_mean

    # Anchor markers: strong contrast between core homos
    anchor_thresh <- quantile(abs(delta), 0.7, na.rm = TRUE)
    anchor_idx <- which(abs(delta) >= anchor_thresh)

    message("[INFO]   Anchors: HOMO_1 core=", length(h1_core),
            " HOMO_2 core=", length(h2_core),
            " anchor markers=", length(anchor_idx))

    # ── HET anchor composition ───────────────────────────────────────────
    het_samples <- intersect(rot[coarse_group_refined == "HET", sample], common)
    if (length(het_samples) >= 3 && length(anchor_idx) >= 10) {
      het_comp <- list()
      for (sid in het_samples) {
        x_s <- X[anchor_idx, sid]
        d_h1 <- abs(x_s - h1_mean[anchor_idx])
        d_h2 <- abs(x_s - h2_mean[anchor_idx])

        valid <- is.finite(d_h1) & is.finite(d_h2)
        frac_aa <- mean(d_h1[valid] < d_h2[valid])
        frac_bb <- mean(d_h2[valid] < d_h1[valid])
        frac_mid <- 1 - frac_aa - frac_bb

        het_comp[[length(het_comp) + 1]] <- data.table(
          candidate_id = cid, sample = sid,
          frac_aa_like = round(frac_aa, 4),
          frac_bb_like = round(frac_bb, 4),
          frac_intermediate = round(frac_mid, 4),
          ab_balance = round(frac_aa - frac_bb, 4),
          n_anchor_markers = length(anchor_idx)
        )
      }
      het_comp_dt <- rbindlist(het_comp)
      fwrite(het_comp_dt, file.path(cand_dir, "candidate_het_anchor_composition.tsv"),
             sep = "\t")
    }
  }

  # ── Export contrast group manifests ─────────────────────────────────────
  manifest_dir <- file.path(cand_dir, "candidate_contrast_manifests")
  ensure_dir(manifest_dir)

  # Build group definitions
  contrasts <- list()

  for (g in c("HOMO_1", "HET", "HOMO_2")) {
    g_all <- rot[coarse_group_refined == g, sample]
    g_core <- if (!is.null(coh)) {
      coh[coarse_group == g & stripe_quality == "core", sample]
    } else g_all

    # Write sample lists
    writeLines(g_all, file.path(manifest_dir, paste0(g, "_all.txt")))
    writeLines(g_core, file.path(manifest_dir, paste0(g, "_core.txt")))

    contrasts[[length(contrasts) + 1]] <- data.table(
      candidate_id = cid,
      group_name = paste0(g, "_all"),
      n_samples = length(g_all),
      type = "all"
    )
    contrasts[[length(contrasts) + 1]] <- data.table(
      candidate_id = cid,
      group_name = paste0(g, "_core"),
      n_samples = length(g_core),
      type = "core"
    )
  }

  # Standard contrasts
  contrast_pairs <- data.table(
    candidate_id = cid,
    contrast = c("HOMO_1_core_vs_HOMO_2_core",
                 "HOMO_1_all_vs_HOMO_2_all",
                 "HOMO_1_core_vs_HET_core",
                 "HET_core_vs_HOMO_2_core"),
    group_A = c("HOMO_1_core", "HOMO_1_all", "HOMO_1_core", "HET_core"),
    group_B = c("HOMO_2_core", "HOMO_2_all", "HET_core", "HOMO_2_core"),
    chrom = chr,
    start_bp = c_start,
    end_bp = c_end
  )

  fwrite(rbindlist(contrasts), file.path(cand_dir, "candidate_contrast_groups.tsv"), sep = "\t")
  fwrite(contrast_pairs, file.path(cand_dir, "candidate_contrast_pairs.tsv"), sep = "\t")

  message("[INFO]   Contrast manifests written to ", manifest_dir)
}

message("[DONE] STEP34 complete")
