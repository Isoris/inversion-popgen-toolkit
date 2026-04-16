#!/usr/bin/env Rscript

# =============================================================================
# STEP36_candidate_clair3_local_signatures.R
#
# Clair3 local-signature support layer for inversion subgroup resolution.
#
# DESIGN PRINCIPLE: Clair3 outputs are NOT the primary source of broad stripe
# definitions. They are a local support layer for:
#   - subgroup resolution within stripes
#   - regime consistency across windows
#   - phase/read support from Clair3-derived local evidence
#   - stabilizing subgroup labels across many markers
#
# For each candidate, integrates:
#   A. Phase-aware mini-haplotype support (from STEP02B)
#   B. Weak-indel local signatures (from STEP04/STEP05)
#   C. BAM-backed regenotyped support (from STEP07)
#   D. Final classified variant markers (from STEP08)
#
# Produces:
#   - candidate_clair3_phase_support.tsv      local phase block summary
#   - candidate_clair3_weak_indel_sharing.tsv  cross-sample weak-indel sharing
#   - candidate_clair3_regen_support.tsv       regenotyped support by group
#   - candidate_clair3_marker_by_group.tsv     variant markers by inversion group
#   - candidate_clair3_local_signature_matrix.tsv.gz  sample × marker local state
#
# Usage:
#   Rscript STEP36_candidate_clair3_local_signatures.R \
#     <config.R> [cid=all] \
#     [--clair3_root <path>]   # root of postprocess_results/<chrom>/<sample>/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

# Parse --clair3_root if present
clair3_root <- NULL
for (i in seq_along(args)) {
  if (args[i] == "--clair3_root" && i < length(args)) {
    clair3_root <- args[i + 1]
  }
}

source(config_file)
ensure_dir(FOLLOWUP_DIR)

# Default Clair3 root
if (is.null(clair3_root)) {
  clair3_root <- file.path(PROJECT_ROOT, "clair3_indel_discovery", "postprocess_results")
}

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

# ── Helper: find Clair3 file for a sample/chrom ──────────────────────────────
find_clair3_file <- function(sample_id, chrom, filename) {
  # Try multiple path patterns
  paths <- c(
    file.path(clair3_root, chrom, sample_id, filename),
    file.path(clair3_root, sample_id, chrom, filename),
    file.path(clair3_root, sample_id, filename)
  )
  for (p in paths) {
    if (file.exists(p)) return(p)
  }
  return(NULL)
}

# ── Helper: find population-level Clair3 file ────────────────────────────────
find_pop_file <- function(chrom, filename) {
  paths <- c(
    file.path(clair3_root, chrom, "_population", filename),
    file.path(clair3_root, "_population", chrom, filename),
    file.path(clair3_root, "_population", filename)
  )
  for (p in paths) {
    if (file.exists(p)) return(p)
  }
  return(NULL)
}

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  if (!dir.exists(cand_dir)) next

  # Load group assignments
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (!file.exists(rot_file)) next
  rot <- fread(rot_file)
  if (nrow(rot) < 5) next

  # Load coherence if available
  coh_file <- file.path(cand_dir, "candidate_sample_coherence.tsv")
  coh <- if (file.exists(coh_file)) fread(coh_file) else NULL

  sample_list <- rot$sample
  n_samples <- length(sample_list)

  message("[INFO] Candidate ", cid, " (", chr, "): Clair3 local signatures")

  # ══════════════════════════════════════════════════════════════════════
  # A. PHASE-AWARE MINI-HAPLOTYPE SUPPORT
  # ══════════════════════════════════════════════════════════════════════

  phase_results <- list()
  for (sid in sample_list) {
    pf <- find_clair3_file(sid, chr, "all_variants_with_phase.tsv")
    if (is.null(pf)) next

    ph <- tryCatch(fread(pf), error = function(e) NULL)
    if (is.null(ph) || nrow(ph) == 0) next

    # Filter to candidate region
    ph_reg <- ph[CHROM == chr & POS >= c_start & POS <= c_end]
    if (nrow(ph_reg) == 0) next

    # Count phase blocks and phased variants
    n_phased <- sum(ph_reg$IS_PHASED == TRUE | ph_reg$IS_PHASED == 1, na.rm = TRUE)
    n_total <- nrow(ph_reg)

    # Phase block sizes
    if ("PHASE_BLOCK_ID" %in% names(ph_reg)) {
      block_sizes <- ph_reg[!is.na(PHASE_BLOCK_ID), .N, by = PHASE_BLOCK_ID]
      n_blocks <- nrow(block_sizes)
      mean_block_size <- if (nrow(block_sizes) > 0) round(mean(block_sizes$N), 2) else 0
      max_block_size <- if (nrow(block_sizes) > 0) max(block_sizes$N) else 0
    } else {
      n_blocks <- 0; mean_block_size <- 0; max_block_size <- 0
    }

    # Phase tier distribution
    tier_tab <- if ("PHASE_TIER" %in% names(ph_reg)) table(ph_reg$PHASE_TIER) else NULL
    n_tier1 <- as.integer(tier_tab["TIER_1_WHATSHAP"]); if (is.na(n_tier1)) n_tier1 <- 0L
    n_tier2 <- as.integer(tier_tab["TIER_2_READPAIR"]); if (is.na(n_tier2)) n_tier2 <- 0L

    phase_results[[length(phase_results) + 1]] <- data.table(
      candidate_id = cid, sample = sid,
      n_variants_in_region = n_total,
      n_phased = n_phased,
      frac_phased = round(n_phased / n_total, 4),
      n_phase_blocks = n_blocks,
      mean_block_size = mean_block_size,
      max_block_size = max_block_size,
      n_tier1_whatshap = n_tier1,
      n_tier2_readpair = n_tier2
    )
  }

  if (length(phase_results) > 0) {
    phase_dt <- rbindlist(phase_results)
    # Merge group labels
    phase_dt <- merge(phase_dt, rot[, .(sample, coarse_group_refined)],
                      by = "sample", all.x = TRUE)
    fwrite(phase_dt, file.path(cand_dir, "candidate_clair3_phase_support.tsv"), sep = "\t")
    message("[INFO]   Phase support: ", nrow(phase_dt), " samples with phase data")
  }

  # ══════════════════════════════════════════════════════════════════════
  # B. WEAK-INDEL SHARING ACROSS SAMPLES
  # ══════════════════════════════════════════════════════════════════════

  pop_clusters_file <- find_pop_file(chr, "weak_indel_population_clusters.tsv")
  if (!is.null(pop_clusters_file)) {
    pop_cl <- tryCatch(fread(pop_clusters_file), error = function(e) NULL)
    if (!is.null(pop_cl) && nrow(pop_cl) > 0) {
      # Filter to candidate region
      pop_reg <- pop_cl[CHROM == chr & POS >= c_start & POS <= c_end]
      pop_reg <- pop_reg[SAMPLE_ID %in% sample_list]

      if (nrow(pop_reg) > 0) {
        # For each cluster: which inversion groups share it?
        cluster_col <- if ("POP_CLUSTER_ID" %in% names(pop_reg)) "POP_CLUSTER_ID" else
                       if ("LOCAL_CLUSTER_ID" %in% names(pop_reg)) "LOCAL_CLUSTER_ID" else NULL

        if (!is.null(cluster_col)) {
          # Merge group labels
          pop_reg <- merge(pop_reg,
                           rot[, .(sample = sample, inv_group = coarse_group_refined)],
                           by.x = "SAMPLE_ID", by.y = "sample", all.x = TRUE)

          sharing <- pop_reg[, .(
            n_samples = uniqueN(SAMPLE_ID),
            n_HOMO_1 = sum(inv_group == "HOMO_1", na.rm = TRUE),
            n_HET = sum(inv_group == "HET", na.rm = TRUE),
            n_HOMO_2 = sum(inv_group == "HOMO_2", na.rm = TRUE),
            groups_present = paste(sort(unique(na.omit(inv_group))), collapse = ","),
            representative_pos = POS[1],
            var_type = VAR_TYPE[1]
          ), by = cluster_col]

          sharing[, candidate_id := cid]
          # Is this signature group-specific?
          sharing[, group_specific := n_samples >= 2 &
                    (n_HOMO_1 == 0 | n_HET == 0 | n_HOMO_2 == 0)]

          fwrite(sharing, file.path(cand_dir, "candidate_clair3_weak_indel_sharing.tsv"),
                 sep = "\t")
          n_specific <- sum(sharing$group_specific)
          message("[INFO]   Weak-indel clusters in region: ", nrow(sharing),
                  " (", n_specific, " group-specific)")
        }
      }
    }
  }

  # ══════════════════════════════════════════════════════════════════════
  # C. BAM-BACKED REGENOTYPED SUPPORT BY GROUP
  # ══════════════════════════════════════════════════════════════════════

  regen_file <- find_pop_file(chr, "regenotyped_rescued_indels.tsv")
  if (!is.null(regen_file)) {
    regen <- tryCatch(fread(regen_file), error = function(e) NULL)
    if (!is.null(regen) && nrow(regen) > 0) {
      regen_reg <- regen[CHROM == chr & POS >= c_start & POS <= c_end]
      regen_reg <- regen_reg[SAMPLE_ID %in% sample_list]

      if (nrow(regen_reg) > 0) {
        # Merge group labels
        regen_reg <- merge(regen_reg,
                           rot[, .(sample = sample, inv_group = coarse_group_refined)],
                           by.x = "SAMPLE_ID", by.y = "sample", all.x = TRUE)

        # Per-marker per-group summary
        regen_summary <- regen_reg[, .(
          n_samples = uniqueN(SAMPLE_ID),
          n_called_alt = sum(REGEN_GT %in% c("0/1", "1/1")),
          n_hom_alt = sum(REGEN_GT == "1/1"),
          n_het_call = sum(REGEN_GT == "0/1"),
          n_hom_ref = sum(REGEN_GT == "0/0"),
          mean_support_frac = round(mean(SUPPORT_FRACTION, na.rm = TRUE), 4),
          mean_mapq = round(mean(MEAN_MAPQ_SUPPORT, na.rm = TRUE), 2)
        ), by = .(REGEN_ID, POS, inv_group)]

        regen_summary[, candidate_id := cid]
        fwrite(regen_summary, file.path(cand_dir, "candidate_clair3_regen_support.tsv"),
               sep = "\t")
        message("[INFO]   Regenotyped markers in region: ",
                uniqueN(regen_summary$REGEN_ID), " × ",
                uniqueN(regen_summary$inv_group), " groups")
      }
    }
  }

  # ══════════════════════════════════════════════════════════════════════
  # D. FINAL CLASSIFIED VARIANTS AS LOCAL MARKERS
  # ══════════════════════════════════════════════════════════════════════

  # Build sample × marker local signature matrix from final classification
  sig_rows <- list()
  marker_positions <- c()

  for (sid in sample_list) {
    fc_file <- find_clair3_file(sid, chr, "final_variant_classification.tsv")
    if (is.null(fc_file)) next

    fc <- tryCatch(fread(fc_file), error = function(e) NULL)
    if (is.null(fc) || nrow(fc) == 0) next

    # Filter to candidate region, keep PASS + rescued
    fc_reg <- fc[CHROM == chr & POS >= c_start & POS <= c_end]
    fc_reg <- fc_reg[FINAL_CLASS %in% c("STRICT_PASS", "RESCUED_STRONG_SINGLE_SAMPLE",
                                         "RESCUED_POPULATION_SIGNATURE")]

    if (nrow(fc_reg) == 0) next

    # Build local state: use GT_CLASS if available, else derive from GT
    gt_col <- if ("GT_CLASS" %in% names(fc_reg)) "GT_CLASS" else
              if ("GT" %in% names(fc_reg)) "GT" else NULL

    if (!is.null(gt_col)) {
      for (ri in seq_len(nrow(fc_reg))) {
        pos_key <- fc_reg$POS[ri]
        gt_val <- fc_reg[[gt_col]][ri]
        # Encode: HOM_REF=0, HET=1, HOM_VAR=2
        state <- if (gt_val %in% c("HOM_REF", "0/0", "0|0")) 0L else
                 if (gt_val %in% c("HET", "0/1", "0|1", "1|0")) 1L else
                 if (gt_val %in% c("HOM_VAR", "1/1", "1|1")) 2L else NA_integer_

        sig_rows[[length(sig_rows) + 1]] <- data.table(
          sample = sid, pos = pos_key, state = state,
          var_type = fc_reg$VAR_TYPE[ri],
          qual = fc_reg$QUAL[ri],
          final_class = fc_reg$FINAL_CLASS[ri]
        )
        marker_positions <- unique(c(marker_positions, pos_key))
      }
    }
  }

  if (length(sig_rows) > 0) {
    sig_dt <- rbindlist(sig_rows)

    # Build wide matrix: samples × positions
    sig_wide <- dcast(sig_dt, sample ~ pos, value.var = "state",
                      fun.aggregate = function(x) {
                        x <- x[!is.na(x)]
                        if (length(x) == 0) return(NA_integer_)
                        as.integer(names(sort(table(x), decreasing = TRUE))[1])
                      })

    # Merge group labels
    sig_wide <- merge(
      rot[, .(sample, coarse_group_refined)],
      sig_wide, by = "sample", all.x = TRUE
    )

    fwrite(sig_wide, file.path(cand_dir, "candidate_clair3_local_signature_matrix.tsv.gz"),
           sep = "\t")

    # Per-marker per-group summary
    marker_cols <- setdiff(names(sig_wide), c("sample", "coarse_group_refined"))
    if (length(marker_cols) >= 5) {
      grp_marker_list <- list()
      for (g in c("HOMO_1", "HET", "HOMO_2")) {
        g_rows <- sig_wide[coarse_group_refined == g]
        if (nrow(g_rows) < 2) next
        for (mc in marker_cols) {
          vals <- g_rows[[mc]]
          vals <- vals[!is.na(vals)]
          if (length(vals) < 2) next
          dom_state <- as.integer(names(sort(table(vals), decreasing = TRUE))[1])
          dom_freq <- max(table(vals)) / length(vals)
          grp_marker_list[[length(grp_marker_list) + 1]] <- data.table(
            candidate_id = cid, pos = as.integer(mc), group = g,
            n_samples = length(vals),
            dominant_state = dom_state,
            dominant_freq = round(dom_freq, 4),
            mean_state = round(mean(vals), 4),
            median_state = round(median(vals), 4)
          )
        }
      }

      if (length(grp_marker_list) > 0) {
        grp_mk <- rbindlist(grp_marker_list)
        fwrite(grp_mk, file.path(cand_dir, "candidate_clair3_marker_by_group.tsv"), sep = "\t")
        message("[INFO]   Clair3 signature matrix: ",
                ncol(sig_wide) - 2, " markers × ", nrow(sig_wide), " samples")
      }
    }
  }

  message("[INFO]   Clair3 local signatures done for candidate ", cid)
}

message("[DONE] STEP36 complete")
