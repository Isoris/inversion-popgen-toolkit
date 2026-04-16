#!/usr/bin/env Rscript

# =============================================================================
# STEP25B_candidate_interpretation_table_v2.R
#
# Upgraded master interpretation table integrating all new analytical layers.
#
# Reads outputs from STEP20-STEP34 and produces an enriched one-row-per-
# candidate interpretation table with fields for:
#   - stripe resolvability
#   - coherent core yes/no
#   - trajectory stability
#   - marker support + polarity summary
#   - multiscale persistence
#   - anchor eligibility
#   - AB asymmetry
#   - expanded final labels
#   - expanded confidence scoring
#
# Usage:
#   Rscript STEP25B_candidate_interpretation_table_v2.R <config.R> [cid=all]
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

interp_list <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  if (!dir.exists(cand_dir)) next

  # ── Load all available component tables ────────────────────────────────
  load_tsv <- function(fn) tryCatch(fread(file.path(cand_dir, fn)), error = function(e) NULL)

  rot    <- load_tsv("candidate_pca_rotated.tsv")
  sub    <- load_tsv("candidate_subclusters.tsv")
  geom   <- load_tsv("candidate_stripe_geometry.tsv")
  het    <- load_tsv("candidate_group_het_summary.tsv")
  marker <- load_tsv("candidate_marker_support.tsv")
  genes  <- load_tsv("candidate_gene_overlap.tsv")
  win    <- load_tsv("candidate_window_summary.tsv")
  coh    <- load_tsv("candidate_sample_coherence.tsv")
  traj   <- load_tsv("candidate_trajectory_summary.tsv")
  mscale <- load_tsv("candidate_multiscale_summary.tsv")
  sg_det <- load_tsv("candidate_stripe_geometry_detailed.tsv")
  anchor <- load_tsv("candidate_anchor_eligibility.tsv")
  ab_asym <- load_tsv("candidate_ab_asymmetry.tsv")
  blocks <- load_tsv("candidate_marker_blocks.tsv")

  if (is.null(rot)) next
  n_samples <- nrow(rot)

  # ── Group counts ───────────────────────────────────────────────────────
  grp_tab <- table(rot$coarse_group_refined)
  n_homo1 <- as.integer(grp_tab["HOMO_1"]); if (is.na(n_homo1)) n_homo1 <- 0L
  n_het <- as.integer(grp_tab["HET"]); if (is.na(n_het)) n_het <- 0L
  n_homo2 <- as.integer(grp_tab["HOMO_2"]); if (is.na(n_homo2)) n_homo2 <- 0L
  n_ambig <- sum(rot$ambiguous_flag, na.rm = TRUE)

  # ── Coherence summary ──────────────────────────────────────────────────
  n_coherent <- n_intermediate <- n_discordant <- NA_integer_
  n_core <- n_peripheral <- n_junk <- NA_integer_
  frac_core <- NA_real_
  if (!is.null(coh)) {
    ct <- table(coh$coherence_class)
    n_coherent <- as.integer(ct["coherent"]); if (is.na(n_coherent)) n_coherent <- 0L
    n_intermediate <- as.integer(ct["intermediate"]); if (is.na(n_intermediate)) n_intermediate <- 0L
    n_discordant <- as.integer(ct["discordant"]); if (is.na(n_discordant)) n_discordant <- 0L

    qt <- table(coh$stripe_quality)
    n_core <- as.integer(qt["core"]); if (is.na(n_core)) n_core <- 0L
    n_peripheral <- as.integer(qt["peripheral"]); if (is.na(n_peripheral)) n_peripheral <- 0L
    n_junk <- as.integer(qt["junk"]); if (is.na(n_junk)) n_junk <- 0L
    frac_core <- round(n_core / n_samples, 4)
  }

  # ── Trajectory stability ───────────────────────────────────────────────
  traj_mean_stability <- traj_frac_stable <- NA_real_
  if (!is.null(traj)) {
    traj_mean_stability <- round(mean(traj$stability, na.rm = TRUE), 4)
    traj_frac_stable <- round(mean(traj$trajectory_class == "stable", na.rm = TRUE), 4)
  }

  # ── Marker polarity ───────────────────────────────────────────────────
  n_markers_reversed <- n_blocks <- NA_integer_
  frac_reversed <- NA_real_
  if (!is.null(blocks)) {
    n_blocks <- nrow(blocks)
  }
  pol <- load_tsv("candidate_marker_polarity.tsv")
  if (!is.null(pol)) {
    n_markers_reversed <- sum(pol$polarity_reversed, na.rm = TRUE)
    frac_reversed <- round(n_markers_reversed / nrow(pol), 4)
  }

  # ── Multiscale persistence ────────────────────────────────────────────
  multiscale_persistent <- NA
  if (!is.null(mscale) && nrow(mscale) >= 2) {
    w500 <- mscale[window_size == 500]
    if (nrow(w500) > 0) {
      multiscale_persistent <- w500$frac_clean[1] >= 0.3
    }
  }

  # ── Anchor eligibility ────────────────────────────────────────────────
  h1_anchor <- h2_anchor <- FALSE
  if (!is.null(anchor)) {
    h1_anchor <- any(anchor[stripe == "HOMO_1", eligible])
    h2_anchor <- any(anchor[stripe == "HOMO_2", eligible])
  }

  # ── AB asymmetry ──────────────────────────────────────────────────────
  ab_mean_balance <- NA_real_
  if (!is.null(ab_asym)) {
    ab_mean_balance <- round(mean(abs(ab_asym$ab_balance), na.rm = TRUE), 4)
  }

  # ── Previous fields (from STEP25 v1) ──────────────────────────────────
  n_snps <- if (!is.null(marker)) marker$n_snps[1] else NA_integer_
  support_tier <- if (!is.null(marker)) marker$support_tier[1] else NA_character_
  n_genes <- if (!is.null(genes)) genes$n_genes_overlap[1] else NA_integer_
  n_windows <- if (!is.null(win)) nrow(win) else NA_integer_
  window_coherence <- if (!is.null(win)) win$candidate_window_coherence[1] else NA_character_

  # ── Stripe geometry ───────────────────────────────────────────────────
  stripe_geom_main <- "unknown"
  if (!is.null(geom)) {
    stripe_geom_main <- paste(geom$geometry_label, collapse = "/")
  }

  # Stripe resolvability from detailed geometry
  stripe_resolvability <- "unknown"
  if (!is.null(sg_det) && nrow(sg_det) > 0 && "elongation_ratio" %in% names(sg_det)) {
    max_elong <- max(sg_det$elongation_ratio, na.rm = TRUE)
    if (is.finite(max_elong)) {
      stripe_resolvability <- if (max_elong < 3) "resolvable" else
                               if (max_elong < 8) "partially_resolvable" else "gradient_like"
    }
  }

  # ── Derive final label ────────────────────────────────────────────────
  has_clean_het_support <- !is.null(het) && nrow(het) >= 3
  has_coherent_core <- !is.na(frac_core) && frac_core >= 0.4
  is_trajectory_stable <- !is.na(traj_frac_stable) && traj_frac_stable >= 0.5
  is_polarity_complex <- !is.na(frac_reversed) && frac_reversed > 0.2
  is_multiscale_ok <- !is.na(multiscale_persistent) && multiscale_persistent

  final_label <- "unclassified"
  if (has_coherent_core && is_trajectory_stable && !is_polarity_complex) {
    if (stripe_resolvability == "resolvable") {
      final_label <- "clean_inversion_like"
    } else {
      final_label <- "inversion_like_with_background_diffusion"
    }
  } else if (is_polarity_complex && !is.na(n_blocks) && n_blocks >= 4) {
    final_label <- "merged_multiblock"
  } else if (!is.na(frac_core) && frac_core < 0.3) {
    final_label <- "diffuse_ambiguous"
  } else if (stripe_resolvability == "gradient_like") {
    final_label <- "gradient_like_nonresolvable"
  } else if (has_coherent_core) {
    final_label <- "inversion_like_with_background_diffusion"
  }

  # ── Confidence ─────────────────────────────────────────────────────────
  conf <- 0
  if (has_coherent_core) conf <- conf + 0.2
  if (is_trajectory_stable) conf <- conf + 0.15
  if (!is_polarity_complex) conf <- conf + 0.1
  if (is_multiscale_ok) conf <- conf + 0.15
  if (n_homo1 >= 5 && n_homo2 >= 5) conf <- conf + 0.1
  if (n_het >= n_homo1 && n_het >= n_homo2) conf <- conf + 0.05
  if (!is.na(support_tier) && support_tier %in% c("high", "medium")) conf <- conf + 0.1
  if (h1_anchor || h2_anchor) conf <- conf + 0.1
  if (!is.na(window_coherence) && window_coherence == "coherent_block") conf <- conf + 0.05
  conf <- min(conf, 1)

  # ── Summary sentence ──────────────────────────────────────────────────
  parts <- c(final_label)
  parts <- c(parts, paste0("H1:", n_homo1, " Het:", n_het, " H2:", n_homo2))
  if (!is.na(frac_core)) parts <- c(parts, paste0("core:", round(frac_core * 100), "%"))
  if (!is.na(traj_frac_stable)) parts <- c(parts, paste0("stable_traj:", round(traj_frac_stable * 100), "%"))
  if (is_polarity_complex) parts <- c(parts, paste0("rev_markers:", n_markers_reversed))
  if (h1_anchor) parts <- c(parts, "H1_anchor")
  if (h2_anchor) parts <- c(parts, "H2_anchor")
  summary_sentence <- paste(parts, collapse = "; ")

  # ── Build row ──────────────────────────────────────────────────────────
  interp_list[[length(interp_list) + 1]] <- data.table(
    candidate_id = cid, chrom = chr,
    start_bp = c_start, end_bp = c_end,
    candidate_length_bp = c_end - c_start,
    n_snps = n_snps, n_samples = n_samples, n_windows = n_windows,
    n_genes_overlap = n_genes,
    n_HOMO_1 = n_homo1, n_HET = n_het, n_HOMO_2 = n_homo2, n_AMBIGUOUS = n_ambig,
    n_coherent = n_coherent, n_intermediate = n_intermediate, n_discordant = n_discordant,
    n_core = n_core, n_peripheral = n_peripheral, n_junk = n_junk,
    frac_core = frac_core,
    trajectory_mean_stability = traj_mean_stability,
    trajectory_frac_stable = traj_frac_stable,
    n_markers_reversed = n_markers_reversed, frac_markers_reversed = frac_reversed,
    n_polarity_blocks = n_blocks,
    multiscale_persistent = multiscale_persistent,
    stripe_resolvability = stripe_resolvability,
    stripe_geometry = stripe_geom_main,
    anchor_HOMO_1 = h1_anchor, anchor_HOMO_2 = h2_anchor,
    ab_mean_balance = ab_mean_balance,
    marker_support_tier = support_tier,
    window_coherence = window_coherence,
    final_label = final_label,
    confidence = round(conf, 3),
    summary_sentence = summary_sentence
  )

  message("[INFO] Candidate ", cid, ": ", final_label,
          " | conf=", round(conf, 2),
          " | core=", n_core, "/", n_samples)
}

if (length(interp_list) > 0) {
  interp_dt <- rbindlist(interp_list, fill = TRUE)
  fwrite(interp_dt,
         file.path(FOLLOWUP_DIR, "candidate_region_interpretation_v2.tsv"),
         sep = "\t")
  message("[DONE] Master interpretation table v2: ", nrow(interp_dt), " candidates")
}

message("[DONE] STEP25B complete")
