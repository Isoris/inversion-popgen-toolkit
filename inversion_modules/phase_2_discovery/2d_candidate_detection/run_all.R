#!/usr/bin/env Rscript
# ============================================================================
# run_all.R — Orchestrator for phase 2 / 2d candidate detection
# ============================================================================
#
# Codebase:    inversion_modules v8.5 / detector v9.3
# Upstream:    phase_2/2c precompute output — $SIM_MATS_DIR/<chr>.sim_mat_nn*.rds
#              (and $PRECOMP_DIR/<chr>.precomp.rds for Phase 9 peeling)
# Downstream:  STEP_D12_bridge_to_C01d.R → triangle_intervals.tsv.gz → C01d
# Formerly:    00_run_all.R (in the inv_detect_v9.3 drop)
#
# Purpose
# -------
# Matrix-based candidate detection track. Runs after 2c precompute, parallel
# to the seed-based track (2c STEP_C01b_1 seeded region-growing). Produces
# the same kind of candidate blocks but from the similarity matrix directly
# rather than from MDS z-outlier seeds. Both tracks feed the same downstream
# scoring (C01d) via the bridge script.
#
# Phases
# ------
#   1. Staircase boundary detection on raw sim_mat                (D01)
#   2. Matrix transforms → 6 treated variants                     (D07)
#   3. Block scoring on each variant                              (D08)
#   4. NN persistence across scales                               (D02)
#   5. NN sweep interval tree                                     (D09)  [optional]
#   6. Consensus across matrix variants                           (D10)
#   7. Evidence stack — flatness / CV / GHSL / SV / local PCA    (D03..D06, D08b)
#   8. Final scoring table assembly
#   9. Blockwise peeling diagnostic                               (D09n) [optional]
#
# Usage
# -----
#   Rscript run_all.R --chr C_gar_LG01 --sim-mat-dir $SIM_MATS_DIR
#   Rscript run_all.R --chr C_gar_LG01 --sim-mat-dir <dir> --phases 1:3
#   Rscript run_all.R --chr C_gar_LG01 --sim-mat-dir <dir> --skip-transforms
#
# Interactive
#   CHR <- "C_gar_LG01"; SIM_MAT_DIR <- "../precomp"
#   source("run_all.R")
#
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# ---- Parse arguments ----
option_list <- list(
  make_option("--chr",         type="character", default=NULL),
  make_option("--sim-mat-dir", type="character", default=NULL),
  make_option("--precomp-dir", type="character", default="",
              help="Precomp RDS dir (for peeling diagnostic, Phase 9)"),
  make_option("--sv-file",     type="character", default=""),
  make_option("--ghsl-dir",    type="character", default=""),
  make_option("--pruned-list", type="character", default="",
              help="pruned_samples.txt or per-chr pruning table (Phase 9)"),
  make_option("--coseg-dir",   type="character", default="",
              help="C01i multi_inversion output dir (Phase 9)"),
  make_option("--outdir",      type="character", default="inv_detect_out_v9.3"),
  make_option("--nn-scales",   type="character", default="0,20,40,80,120,160,200,240,320"),
  make_option("--phases",      type="character", default="all",
              help="Which phases: 'all', '1', '1:3', '2:5', etc."),
  make_option("--skip-transforms", type="logical",   default=FALSE, action="store_true",
              help="Skip matrix transforms (Phase 2) for faster initial test"),
  make_option("--skip-tree",   type="logical",   default=FALSE, action="store_true",
              help="Skip full NN sweep tree (Phase 5)")
)

if (interactive()) {
  ARGS <- list(
    chr         = if (exists("CHR")) CHR else Sys.getenv("CHR", "C_gar_LG01"),
    sim_mat_dir = if (exists("SIM_MAT_DIR")) SIM_MAT_DIR else Sys.getenv("SIM_MAT_DIR", ""),
    precomp_dir = if (exists("PRECOMP_DIR")) PRECOMP_DIR else Sys.getenv("PRECOMP_DIR", ""),
    sv_file     = if (exists("SV_FILE")) SV_FILE else "",
    ghsl_dir    = if (exists("GHSL_DIR")) GHSL_DIR else "",
    pruned_list = if (exists("PRUNED_LIST")) PRUNED_LIST else "",
    coseg_dir   = if (exists("COSEG_DIR")) COSEG_DIR else "",
    outdir      = if (exists("OUTDIR")) OUTDIR else "inv_detect_out_v9.3",
    nn_scales   = if (exists("NN_SCALES")) NN_SCALES else "0,20,40,80,120,160,200,240,320",
    phases      = if (exists("PHASES")) PHASES else "all",
    skip_transforms = if (exists("SKIP_TRANSFORMS")) SKIP_TRANSFORMS else FALSE,
    skip_tree   = if (exists("SKIP_TREE")) SKIP_TREE else FALSE
  )
} else {
  ARGS <- parse_args(OptionParser(option_list=option_list))
}

nn_scales <- as.integer(strsplit(ARGS$nn_scales, ",")[[1]])
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error=function(e) ".")
if (script_dir == "." || script_dir == "") script_dir <- getwd()

# Parse phases
if (ARGS$phases == "all") {
  run_phases <- 1:9
} else if (grepl(":", ARGS$phases)) {
  parts <- as.integer(strsplit(ARGS$phases, ":")[[1]])
  run_phases <- parts[1]:parts[2]
} else {
  run_phases <- as.integer(strsplit(ARGS$phases, ",")[[1]])
}

dir.create(ARGS$outdir, recursive=TRUE, showWarnings=FALSE)

# ---- Source all modules ----
source(file.path(script_dir, "00_config.R"))
source(file.path(script_dir, "STEP_D01_staircase_boundaries.R"))
source(file.path(script_dir, "STEP_D02_nn_persistence.R"))
source(file.path(script_dir, "STEP_D03_flatness.R"))
source(file.path(script_dir, "STEP_D04_interior_cv.R"))
source(file.path(script_dir, "STEP_D05_ghsl_stability.R"))
source(file.path(script_dir, "STEP_D06_sv_overlap.R"))
source(file.path(script_dir, "STEP_D07_matrix_transforms.R"))
source(file.path(script_dir, "STEP_D08_block_scoring.R"))
source(file.path(script_dir, "STEP_D08b_local_pca.R"))
source(file.path(script_dir, "STEP_D09_nn_sweep_tree.R"))
source(file.path(script_dir, "STEP_D10_variant_consensus.R"))

cat("============================================================\n")
cat("INVERSION DETECTOR v9.3 — Full Pipeline\n")
cat("Chromosome:", ARGS$chr, "\n")
cat("Phases:", paste(run_phases, collapse=", "), "\n")
cat("NN scales:", paste(nn_scales, collapse=", "), "\n")
cat("Output:", ARGS$outdir, "\n")
cat("============================================================\n\n")


# ---- Helper: load sim_mat ----
load_sim_mat <- function(chr, nn, sim_mat_dir) {
  paths <- c(
    file.path(sim_mat_dir, sprintf("%s.sim_mat_nn%d.rds", chr, nn)),
    file.path(sim_mat_dir, sprintf("%s_sim_mat_nn%d.rds", chr, nn)),
    file.path(sim_mat_dir, "..", "triangles_v3", sprintf("nn%d", nn),
              sprintf("%s_sim_mat_nn%d.rds", chr, nn))
  )
  if (nn == 0) {
    paths <- c(paths,
      file.path(sim_mat_dir, sprintf("%s.precomp.rds", chr)),
      file.path(sim_mat_dir, sprintf("%s_precomp.rds", chr)),
      file.path(sim_mat_dir, "..", sprintf("%s.precomp.rds", chr))
    )
  }
  for (p in paths) {
    if (file.exists(p)) {
      cat("  Loading:", p, "\n")
      obj <- readRDS(p)
      if (is.matrix(obj)) return(obj)
      if (is.list(obj) && "sim_mat" %in% names(obj)) return(obj$sim_mat)
      return(obj)
    }
  }
  return(NULL)
}

# ---- Load similarity matrices ----
cat("Loading similarity matrices...\n")
sim_mats <- list()
for (nn in nn_scales) {
  sm <- load_sim_mat(ARGS$chr, nn, ARGS$sim_mat_dir)
  if (!is.null(sm)) sim_mats[[as.character(nn)]] <- sm
}
if (length(sim_mats) == 0) stop("No sim_mats loaded. Check --sim-mat-dir.")
smat <- sim_mats[[as.character(nn_scales[1])]]
if (is.null(smat)) smat <- sim_mats[[1]]
cat("  Loaded", length(sim_mats), "NN scales. Primary:", nrow(smat), "x", ncol(smat), "\n\n")


# =====================================================
# PHASE 1: Staircase detector on raw sim_mat
# =====================================================
blocks_file    <- file.path(ARGS$outdir, sprintf("blocks_%s.tsv", ARGS$chr))
staircase_file <- file.path(ARGS$outdir, sprintf("staircase_%s.rds", ARGS$chr))

if (1 %in% run_phases) {
  cat("=== PHASE 1: Staircase detector ===\n")
  staircase_result <- detect_blocks_staircase(smat)

  fwrite(staircase_result$blocks, blocks_file, sep="\t")
  saveRDS(staircase_result, staircase_file)
  cat("Blocks found:", nrow(staircase_result$blocks), "\n")

  # Save artifact registry
  if (nrow(staircase_result$artifacts) > 0) {
    art_file <- file.path(ARGS$outdir, sprintf("artifacts_%s.tsv", ARGS$chr))
    fwrite(staircase_result$artifacts, art_file, sep="\t")
    cat("Artifacts logged:", nrow(staircase_result$artifacts), "\n")
  }
  cat("\n")
} else {
  cat("Phase 1: Loading from disk...\n")
  staircase_result <- readRDS(staircase_file)
}

blocks <- staircase_result$blocks
if (nrow(blocks) == 0) {
  cat("No blocks found. Exiting.\n")
  quit(save="no")
}

# Compatibility: ensure candidate_id column exists
if (!"candidate_id" %in% names(blocks)) {
  blocks[, candidate_id := block_id]
}
candidates <- copy(blocks)


# =====================================================
# PHASE 2: Matrix transforms → generate treated variants
# =====================================================
variants_file <- file.path(ARGS$outdir, sprintf("variants_%s.rds", ARGS$chr))

if (2 %in% run_phases && !ARGS$skip_transforms) {
  cat("=== PHASE 2: Matrix transforms ===\n")
  variants <- generate_all_variants(smat)
  saveRDS(variants, variants_file)
  cat("Saved", length(variants), "variants.\n\n")
} else if (file.exists(variants_file)) {
  cat("Phase 2: Loading variants from disk...\n")
  variants <- readRDS(variants_file)
} else {
  cat("Phase 2: Skipped. Using raw only.\n")
  variants <- list(raw = smat)
}


# =====================================================
# PHASE 3: Block scoring on each variant
# =====================================================
scores_file <- file.path(ARGS$outdir, sprintf("block_scores_%s.tsv", ARGS$chr))

if (3 %in% run_phases) {
  cat("=== PHASE 3: Block scoring ===\n")
  scores_all <- score_blocks_all_variants(blocks, variants)
  fwrite(scores_all, scores_file, sep="\t")
  cat("Scores written:", nrow(scores_all), "rows\n\n")
} else if (file.exists(scores_file)) {
  scores_all <- fread(scores_file)
} else {
  scores_all <- NULL
}


# =====================================================
# PHASE 4: NN persistence (lightweight 4-scale)
# =====================================================
nn_file <- file.path(ARGS$outdir, sprintf("nn_persistence_%s.tsv", ARGS$chr))

if (4 %in% run_phases && length(sim_mats) >= 2) {
  cat("=== PHASE 4: NN persistence ===\n")
  nn_evidence <- compute_nn_persistence_v2(sim_mats)
  fwrite(nn_evidence, nn_file, sep="\t")
  cat("NN evidence written\n\n")
} else if (file.exists(nn_file)) {
  nn_evidence <- fread(nn_file)
} else {
  nn_evidence <- NULL
}


# =====================================================
# PHASE 5: Full NN sweep interval tree (optional)
# =====================================================
tree_file <- file.path(ARGS$outdir, sprintf("nn_tree_%s.tsv", ARGS$chr))

if (5 %in% run_phases && !ARGS$skip_tree && length(sim_mats) >= 3) {
  cat("=== PHASE 5: NN sweep interval tree ===\n")
  tree <- build_nn_sweep_tree(sim_mats)
  fwrite(tree, tree_file, sep="\t")
  print_tree_summary(tree)
  cat("\n")
} else if (file.exists(tree_file)) {
  tree <- fread(tree_file)
} else {
  tree <- NULL
}


# =====================================================
# PHASE 6: Consensus across matrix variants
# =====================================================
consensus_file <- file.path(ARGS$outdir, sprintf("consensus_%s.tsv", ARGS$chr))

if (6 %in% run_phases && length(variants) > 1) {
  cat("=== PHASE 6: Consensus across variants ===\n")

  # Run staircase on each variant
  blocks_by_variant <- list()
  blocks_by_variant[["raw"]] <- blocks

  for (vname in setdiff(names(variants), "raw")) {
    cat("  Staircase on", vname, "...\n")
    vresult <- detect_blocks_staircase(variants[[vname]])
    blocks_by_variant[[vname]] <- vresult$blocks
    cat("    Blocks:", nrow(vresult$blocks), "\n")
  }

  consensus <- build_consensus(blocks_by_variant, scores_all)
  fwrite(consensus, consensus_file, sep="\t")
  print_consensus(consensus)
  cat("\n")
} else if (file.exists(consensus_file)) {
  consensus <- fread(consensus_file)
} else {
  consensus <- NULL
}


# =====================================================
# PHASE 7: Evidence stack for survivors
# =====================================================

if (7 %in% run_phases) {
  cat("=== PHASE 7: Evidence stack ===\n")

  # 7a: Flatness
  cat("  7a: Flatness / decay profile...\n")
  flatness_ev <- compute_flatness(candidates, smat)
  fwrite(as.data.table(flatness_ev),
         file.path(ARGS$outdir, sprintf("ev_flatness_%s.tsv", ARGS$chr)),
         sep="\t")

  # 7b: Interior CV
  cat("  7b: Interior CV...\n")
  cv_ev <- compute_interior_cv(candidates, smat)
  fwrite(as.data.table(cv_ev),
         file.path(ARGS$outdir, sprintf("ev_cv_%s.tsv", ARGS$chr)),
         sep="\t")

  # 7c: GHSL partition stability
  cat("  7c: GHSL partition stability...\n")
  ghsl_ev <- compute_ghsl_from_simmat(candidates, smat)
  fwrite(as.data.table(ghsl_ev),
         file.path(ARGS$outdir, sprintf("ev_ghsl_%s.tsv", ARGS$chr)),
         sep="\t")

  # 7d: SV overlap
  if (ARGS$sv_file != "" && file.exists(ARGS$sv_file)) {
    cat("  7d: SV breakpoint overlap...\n")
    sv_ev <- compute_sv_overlap(candidates, ARGS$sv_file, ARGS$chr,
                                 window_size_bp=CFG$WINDOW_SIZE_BP)
    fwrite(as.data.table(sv_ev),
           file.path(ARGS$outdir, sprintf("ev_sv_%s.tsv", ARGS$chr)),
           sep="\t")
  } else {
    cat("  7d: SV overlap skipped (no SV file)\n")
    sv_ev <- NULL
  }

  # 7e: Local PCA/ICA
  cat("  7e: Local PCA/ICA...\n")
  pca_ev <- compute_local_pca(candidates, smat, do_ica=CFG$ICA_ENABLED)
  fwrite(as.data.table(pca_ev),
         file.path(ARGS$outdir, sprintf("ev_pca_%s.tsv", ARGS$chr)),
         sep="\t")

  cat("  Evidence stack complete.\n\n")
} else {
  # Try loading
  flat_f <- file.path(ARGS$outdir, sprintf("ev_flatness_%s.tsv", ARGS$chr))
  flatness_ev <- if (file.exists(flat_f)) fread(flat_f) else NULL

  cv_f <- file.path(ARGS$outdir, sprintf("ev_cv_%s.tsv", ARGS$chr))
  cv_ev <- if (file.exists(cv_f)) fread(cv_f) else NULL

  ghsl_f <- file.path(ARGS$outdir, sprintf("ev_ghsl_%s.tsv", ARGS$chr))
  ghsl_ev <- if (file.exists(ghsl_f)) fread(ghsl_f) else NULL

  sv_f <- file.path(ARGS$outdir, sprintf("ev_sv_%s.tsv", ARGS$chr))
  sv_ev <- if (file.exists(sv_f)) fread(sv_f) else NULL

  pca_f <- file.path(ARGS$outdir, sprintf("ev_pca_%s.tsv", ARGS$chr))
  pca_ev <- if (file.exists(pca_f)) fread(pca_f) else NULL
}


# =====================================================
# PHASE 8: Final scoring table assembly
# =====================================================
if (8 %in% run_phases) {
  cat("=== PHASE 8: Final scoring table ===\n")

  tab <- copy(candidates)
  tab[, chr := ARGS$chr]

  # Join evidence
  join_ev <- function(tab, ev) {
    if (!is.null(ev) && nrow(ev) > 0) {
      ev_dt <- as.data.table(ev)
      if ("candidate_id" %in% names(ev_dt)) {
        tab <- merge(tab, ev_dt, by="candidate_id", all.x=TRUE)
      }
    }
    return(tab)
  }

  tab <- join_ev(tab, flatness_ev)
  tab <- join_ev(tab, cv_ev)
  tab <- join_ev(tab, ghsl_ev)
  tab <- join_ev(tab, sv_ev)
  tab <- join_ev(tab, pca_ev)

  # Join NN persistence
  if (!is.null(nn_evidence) && nrow(nn_evidence) > 0) {
    tab <- merge(tab, nn_evidence, by="candidate_id", all.x=TRUE)
  }

  # BUGFIX 2026-04-17 (FIX 16, DESIGN): merge nn_birth from the NN sweep
  # tree (Phase 5, D09). C01d's D3 scoring reads iv$nn_birth (L210) and
  # expects it to saturate at 200. Before this fix nn_birth was never
  # merged into the scoring table — D09 wrote nn_tree_<chr>.tsv to disk
  # but Phase 8 never joined it in. Result: iv$nn_birth was always NA,
  # safe_num(NA, 0) = 0, nn_birth_score = 0, and D3 capped at 0.6 instead
  # of 1.0 on every candidate, every chromosome.
  #
  # Mapping strategy: tree nodes are intervals in bin coordinates
  # (matching blocks$start/end). For each block, pick the tree node whose
  # interval has the highest reciprocal overlap with the block. Take that
  # node's nn_birth. If no node overlaps, leave as NA (not 0 — the
  # downstream safe_num(NA, 0) already handles that correctly).
  if (!is.null(tree) && nrow(tree) > 0 &&
      all(c("start", "end", "nn_birth") %in% names(tree))) {
    nn_birth_per_block <- numeric(nrow(tab))
    nn_birth_per_block[] <- NA_real_
    for (ti in seq_len(nrow(tab))) {
      bs <- tab$start[ti]; be <- tab$end[ti]
      if (!is.finite(bs) || !is.finite(be)) next
      bw <- be - bs + 1
      # Reciprocal overlap with every tree node
      ov_num <- pmax(0, pmin(tree$end, be) - pmax(tree$start, bs) + 1)
      tw <- tree$end - tree$start + 1
      recip <- pmin(ov_num / bw, ov_num / pmax(tw, 1))
      best <- which.max(recip)
      if (length(best) > 0 && recip[best] >= 0.5) {
        nn_birth_per_block[ti] <- tree$nn_birth[best]
      }
    }
    tab[, nn_birth := nn_birth_per_block]
    cat("  Merged nn_birth from tree for",
        sum(!is.na(nn_birth_per_block)), "/", nrow(tab), "blocks\n")
  }

  # BUGFIX 2026-04-17 (FIX 22, DESIGN): wire 2e GHSL (C04 snake3
  # phased-Clair3 haplotype divergence) into the scoring table.
  # Chat 14 (2026-04-18): rewired from v5 to v6. v6 annot RDS uses
  # v6-native column names — ghsl_v6_score / ghsl_v6_status plus
  # unprefixed rank_stability, div_contrast_z, div_bimodal. Paths now
  # <chr>.ghsl_v6.annot.rds.
  #
  # Before this fix, C04's per-window annot.rds shards existed but were
  # read by nobody in the main C01d pipeline — they were consumed only
  # by two phase_5_followup supplementary-figure scripts. D10 in C01d
  # read partition_stability from D05's sim_mat clustering, never from
  # C04. The 2c_precomp README's "Layer C = GHSL" 4-layer independence
  # framework was aspirational on the scoring path.
  #
  # Wiring strategy: C04 writes <chr>.ghsl_v6.annot.rds keyed by
  # global_window_id, where global_window_id = bin index matching
  # D01/D09's start/end columns. For each block, aggregate across the
  # block's windows (cols the interval [start, end]):
  #   ghsl_v6_score_max       : peak composite score within block
  #   ghsl_rank_stability_max : peak inter-window rank stability
  #   ghsl_div_contrast_z_max : peak bimodality-contrast z-score
  #   ghsl_div_bimodal_frac   : fraction of windows with >=2 density modes
  #   ghsl_pass_frac          : fraction of windows with status=="PASS"
  #   ghsl_n_scored_windows   : windows in block that C04 scored
  #
  # If --ghsl-dir is empty or no <chr>.ghsl_v6.annot.rds exists (Clair3
  # still running on this chromosome), the merge is a no-op. C01d handles
  # missing columns via safe_num()/%||% fallbacks.
  if (ARGS$ghsl_dir != "") {
    ghsl_annot_file <- file.path(ARGS$ghsl_dir, "annot",
                                  paste0(ARGS$chr, ".ghsl_v6.annot.rds"))
    # Fallback: accept annot.rds directly in ghsl_dir (no annot/ subdir)
    if (!file.exists(ghsl_annot_file)) {
      ghsl_annot_file <- file.path(ARGS$ghsl_dir,
                                    paste0(ARGS$chr, ".ghsl_v6.annot.rds"))
    }
    # Chat 14 back-compat: also accept v5 annot if v6 not yet emitted on
    # this chromosome (allows a partial transition when not every chrom
    # has been re-classified with v6 yet). v5 column names mapped below.
    using_v5_compat <- FALSE
    if (!file.exists(ghsl_annot_file)) {
      v5_candidate <- file.path(ARGS$ghsl_dir, "annot",
                                  paste0(ARGS$chr, ".ghsl_v5.annot.rds"))
      if (!file.exists(v5_candidate)) {
        v5_candidate <- file.path(ARGS$ghsl_dir,
                                    paste0(ARGS$chr, ".ghsl_v5.annot.rds"))
      }
      if (file.exists(v5_candidate)) {
        ghsl_annot_file <- v5_candidate
        using_v5_compat <- TRUE
      }
    }
    if (file.exists(ghsl_annot_file)) {
      cat(sprintf("  Loading GHSL %s annot: %s\n",
                   if (using_v5_compat) "v5 (compat)" else "v6",
                   ghsl_annot_file))
      ghsl_annot <- tryCatch(readRDS(ghsl_annot_file),
                              error = function(e) NULL)
      # Column-name resolution: v6 uses ghsl_v6_score + unprefixed
      # rank_stability, div_contrast_z, div_bimodal. v5 used ghsl_v5_score
      # + ghsl_rank_stability + ghsl_div_contrast_z + ghsl_div_bimodal.
      if (!is.null(ghsl_annot) && nrow(ghsl_annot) > 0 &&
          "global_window_id" %in% names(ghsl_annot)) {
        if (using_v5_compat) {
          score_col    <- "ghsl_v5_score"
          rstab_col    <- "ghsl_rank_stability"
          divcz_col    <- "ghsl_div_contrast_z"
          divbim_col   <- "ghsl_div_bimodal"
        } else {
          score_col    <- "ghsl_v6_score"
          rstab_col    <- "rank_stability"
          divcz_col    <- "div_contrast_z"
          divbim_col   <- "div_bimodal"
        }
        # Aggregate per block
        ghsl_per_block <- data.table(
          block_id                 = integer(nrow(tab)),
          ghsl_v6_score_max        = NA_real_,
          ghsl_rank_stability_max  = NA_real_,
          ghsl_div_contrast_z_max  = NA_real_,
          ghsl_div_bimodal_frac    = NA_real_,
          ghsl_pass_frac           = NA_real_,
          ghsl_n_scored_windows    = 0L
        )
        for (ti in seq_len(nrow(tab))) {
          bs <- tab$start[ti]; be <- tab$end[ti]
          bid <- tab$block_id[ti]
          ghsl_per_block$block_id[ti] <- bid
          if (!is.finite(bs) || !is.finite(be)) next
          win_rows <- ghsl_annot[global_window_id >= bs & global_window_id <= be]
          score_vals <- if (score_col %in% names(win_rows)) win_rows[[score_col]] else NA_real_
          n_scored <- sum(!is.na(score_vals))
          ghsl_per_block$ghsl_n_scored_windows[ti] <- n_scored
          if (n_scored >= 3) {
            ghsl_per_block$ghsl_v6_score_max[ti] <-
              max(score_vals, na.rm = TRUE)
            if (rstab_col %in% names(win_rows)) {
              ghsl_per_block$ghsl_rank_stability_max[ti] <-
                max(win_rows[[rstab_col]], na.rm = TRUE)
            }
            if (divcz_col %in% names(win_rows)) {
              ghsl_per_block$ghsl_div_contrast_z_max[ti] <-
                max(win_rows[[divcz_col]], na.rm = TRUE)
            }
            if (divbim_col %in% names(win_rows)) {
              ghsl_per_block$ghsl_div_bimodal_frac[ti] <-
                mean(win_rows[[divbim_col]] >= 2, na.rm = TRUE)
            }
            # PASS fraction — v6 emits ghsl_v6_status directly; v5 did
            # not carry status on the annot RDS and we approximated with
            # score > 0.65. Honor that distinction.
            if ("ghsl_v6_status" %in% names(win_rows)) {
              ghsl_per_block$ghsl_pass_frac[ti] <-
                mean(win_rows$ghsl_v6_status == "PASS", na.rm = TRUE)
            } else {
              ghsl_per_block$ghsl_pass_frac[ti] <-
                mean(score_vals > 0.65, na.rm = TRUE)
            }
          }
        }
        # Replace infinite max()/mean() on empty-after-na.rm with NA
        for (col in c("ghsl_v6_score_max", "ghsl_rank_stability_max",
                       "ghsl_div_contrast_z_max", "ghsl_div_bimodal_frac",
                       "ghsl_pass_frac")) {
          vals <- ghsl_per_block[[col]]
          vals[!is.finite(vals)] <- NA_real_
          ghsl_per_block[[col]] <- vals
        }
        tab <- merge(tab, ghsl_per_block, by = "block_id", all.x = TRUE)
        cat("  Merged GHSL v6 annotations for",
            sum(ghsl_per_block$ghsl_n_scored_windows > 0), "/",
            nrow(tab), "blocks\n")
      } else {
        cat("  GHSL annot file present but empty or malformed — skipping\n")
      }
    } else {
      cat("  Phase 8 GHSL merge: no annot.rds at ", ghsl_annot_file,
          " (Clair3 still running?) — skipping\n")
    }
  }

  # Join block scores (raw variant)
  # BUGFIX 2026-04-17 (FIX 15, SILENT): D01 writes a `contrast` column
  # (block_height - bg_sim), and D08 also writes `contrast`
  # (inside_mean - flank_outside_mean). A plain merge on block_id produced
  # `contrast.x` (D01) and `contrast.y` (D08), so `iv$contrast` in C01d
  # was NULL → safe_num(NULL,0)=0 → D1_block_strength silently collapsed
  # to 0.35*squareness + 0.25*sharpness. Drop D01's contrast before the
  # merge so D08's (more informative, flank-based, variant-aware) version
  # is what C01d sees. D01's original `height` and `bg_sim` columns still
  # survive for anyone who wants the global-background contrast.
  if (!is.null(scores_all) && nrow(scores_all) > 0) {
    raw_scores <- scores_all[variant == "raw"]
    if (nrow(raw_scores) > 0) {
      if ("contrast" %in% names(tab)) tab[, contrast := NULL]
      tab <- merge(tab, raw_scores[, .(block_id, contrast, squareness,
                                        occupancy, patchiness, sharpness,
                                        shape_class)],
                   by="block_id", all.x=TRUE)
    }
  }

  setorder(tab, start)

  final_file <- file.path(ARGS$outdir, sprintf("scoring_table_%s.tsv", ARGS$chr))
  fwrite(tab, final_file, sep="\t")

  cat("Final scoring table:", nrow(tab), "candidates,", ncol(tab), "columns\n")
  cat("Wrote:", final_file, "\n")

  # Quick calibration check against known positives/negatives
  cat("\n--- Calibration Check ---\n")
  for (kp in CFG$KNOWN_POSITIVES) {
    if (kp$chr != ARGS$chr) next
    kp_start <- mb_to_bin(kp$start_mb)
    kp_end   <- mb_to_bin(kp$end_mb)
    hits <- tab[start <= kp_end & end >= kp_start]
    cat(sprintf("  %s (%s, %.1f-%.1f Mb): %d blocks overlap",
                kp$name, kp$chr, kp$start_mb, kp$end_mb, nrow(hits)))
    if (nrow(hits) > 0 && "shape_class" %in% names(hits)) {
      cat(" | shapes:", paste(unique(hits$shape_class), collapse=","))
    }
    if (nrow(hits) > 0 && "survives_nn40" %in% names(hits)) {
      cat(" | nn40:", sum(hits$survives_nn40, na.rm=TRUE), "/", nrow(hits))
    }
    cat("\n")
  }

  for (kn in CFG$KNOWN_NEGATIVES) {
    if (kn$chr != ARGS$chr) next
    kn_start <- mb_to_bin(kn$start_mb)
    kn_end   <- mb_to_bin(kn$end_mb)
    hits <- tab[start <= kn_end & end >= kn_start]
    cat(sprintf("  %s (%s, %.1f-%.1f Mb): %d blocks overlap (should be 0 or weak)",
                kn$name, kn$chr, kn$start_mb, kn$end_mb, nrow(hits)))
    if (nrow(hits) > 0 && "shape_class" %in% names(hits)) {
      cat(" | shapes:", paste(unique(hits$shape_class), collapse=","))
    }
    cat("\n")
  }
}

# =====================================================
# PHASE 9: Blockwise peeling diagnostic (STEP_C01n)
# =====================================================
peel_file <- file.path(ARGS$outdir, sprintf("peel_diagnostic_%s.tsv", ARGS$chr))

if (9 %in% run_phases) {
  cat("=== PHASE 9: Blockwise peeling diagnostic ===\n")

  # Need precomp RDS (which has dt with PC_1_<sample> columns)
  precomp_dir_use <- ARGS$precomp_dir
  if (precomp_dir_use == "") precomp_dir_use <- ARGS$sim_mat_dir

  precomp_file <- file.path(precomp_dir_use,
                             sprintf("%s.precomp.rds", ARGS$chr))
  # Try alternative paths
  if (!file.exists(precomp_file)) {
    precomp_file <- file.path(precomp_dir_use, "..",
                               sprintf("%s.precomp.rds", ARGS$chr))
  }

  if (file.exists(precomp_file)) {
    source(file.path(script_dir, "STEP_D09n_peeling_diagnostic.R"))

    precomp <- readRDS(precomp_file)
    peel_result <- run_peel_diagnostic(
      precomp, blocks,
      pruned_list = if (ARGS$pruned_list != "") ARGS$pruned_list else NULL,
      coseg_dir   = if (ARGS$coseg_dir != "") ARGS$coseg_dir else NULL,
      chr_name    = ARGS$chr,
      max_peel_frac = 0.3,
      local_k = 5L,
      n_pcs_use = 2L
    )

    fwrite(peel_result$diagnostics, peel_file, sep = "\t")
    fwrite(peel_result$sample_groups,
           file.path(ARGS$outdir, sprintf("peel_groups_%s.tsv", ARGS$chr)),
           sep = "\t")
    print_interpretation(peel_result$diagnostics)

    # Merge key peel results into final scoring table if Phase 8 ran
    if (8 %in% run_phases && nrow(peel_result$diagnostics) > 0) {
      # Add L1b effect as a column (the most informative peel)
      l1b <- peel_result$diagnostics[peel_mode == "L1b_chrlocal_kin",
                                      .(block_id, l1b_effect = effect_class,
                                        l1b_contrast_ratio = round(after_contrast / pmax(before_contrast, 0.001), 3))]
      l2 <- peel_result$diagnostics[peel_mode == "L2_local_coseg",
                                     .(block_id, l2_effect = effect_class,
                                       l2_contrast_ratio = round(after_contrast / pmax(before_contrast, 0.001), 3))]

      if (nrow(l1b) > 0) tab <- merge(tab, l1b, by = "block_id", all.x = TRUE)
      if (nrow(l2) > 0)  tab <- merge(tab, l2, by = "block_id", all.x = TRUE)

      # Re-write final scoring table with peel columns
      fwrite(tab, final_file, sep = "\t")
      cat("  Peel columns added to scoring table\n")
    }
  } else {
    cat("  Phase 9 skipped: no precomp RDS found at ", precomp_file, "\n")
    cat("  (Need --precomp-dir pointing to precomp/<chr>.precomp.rds)\n")
  }
}

cat("\n============================================================\n")
cat("DONE. v9.3 pipeline complete for", ARGS$chr, "\n")
cat("============================================================\n")
