#!/usr/bin/env Rscript

# =============================================================================
# STEP_C64_cross_inversion_synthesis.R — 6E synthesis layer
#
# Cross-inversion convergence analysis. Two main outputs:
#
# 1. Family / GO convergence table — for each family or GO term, how many
#    inversions independently captured it? With null comparison: how many
#    inversions would be expected to capture it by chance given the family's
#    genome-wide frequency?
#
# 2. Top-candidate gene table — across all inversions, which genes show the
#    strongest combined evidence: large |delta_burden| + significant config_mi
#    + arrangement-divergent FST + (optionally) high cargo_lof_count?
#
# Inputs:
#   ${CARGO_INVENTORY_DIR}/<cid>/genes.tsv
#   ${CARGO_INVENTORY_DIR}/diagnostic_table.tsv
#   ${CARGO_ENRICH_DIR}/<cid>/{go,kegg,family}_enrichment.tsv
#   results_registry pairwise + interval_summary rows for cargo_* stats
#
# Outputs:
#   ${CARGO_SYNTH_DIR}/family_convergence.tsv
#   ${CARGO_SYNTH_DIR}/go_convergence.tsv
#   ${CARGO_SYNTH_DIR}/top_candidate_genes.tsv
#   ${CARGO_SYNTH_DIR}/cross_inversion_summary.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

CARGO_INVENTORY_DIR <- Sys.getenv("CARGO_INVENTORY_DIR")
CARGO_ENRICH_DIR    <- Sys.getenv("CARGO_ENRICH_DIR")
CARGO_BURDEN_DIR    <- Sys.getenv("CARGO_BURDEN_DIR")
CARGO_CONFIG_DIR    <- Sys.getenv("CARGO_CONFIG_DIR")
CARGO_SYNTH_DIR     <- Sys.getenv("CARGO_SYNTH_DIR")
RESULTS_REGISTRY_DIR <- Sys.getenv("RESULTS_REGISTRY_DIR")
REGISTRY_LOADER_R    <- Sys.getenv("REGISTRY_LOADER_R")

stopifnot(nzchar(CARGO_INVENTORY_DIR), nzchar(CARGO_SYNTH_DIR))
dir.create(CARGO_SYNTH_DIR, recursive = TRUE, showWarnings = FALSE)

# Try to source the R registry; fall back to direct file reads if not available.
reg <- NULL
if (nzchar(REGISTRY_LOADER_R) && file.exists(REGISTRY_LOADER_R)) {
  source(REGISTRY_LOADER_R)
  reg <- tryCatch(load_registry(), error = function(e) {
    message("[C64] R registry unavailable: ", conditionMessage(e))
    NULL
  })
}

# ── Load diagnostic table ──
diag <- fread(file.path(CARGO_INVENTORY_DIR, "diagnostic_table.tsv"))
candidates <- diag[n_genes_inside > 0, candidate_id]
message("[C64] ", length(candidates), " candidates with cargo")

# ── 1. Family convergence ──
all_inv_genes <- rbindlist(lapply(candidates, function(cid) {
  gf <- file.path(CARGO_INVENTORY_DIR, cid, "genes.tsv")
  if (!file.exists(gf)) return(NULL)
  g <- fread(gf)
  g[, candidate_id := cid]
  g
}), fill = TRUE)

if (nrow(all_inv_genes) == 0) {
  message("[C64] No inventory data found")
  quit(status = 0)
}

# Family per (cid, family) → 1 (presence)
fam_pres <- all_inv_genes[nzchar(family),
                          .(n_genes = .N), by = .(candidate_id, family)]
fam_summary <- fam_pres[, .(n_inversions_captured = uniqueN(candidate_id),
                              total_genes_across_inversions = sum(n_genes),
                              max_in_single_inversion = max(n_genes)),
                         by = family]
setorder(fam_summary, -n_inversions_captured, -total_genes_across_inversions)

# Quick null: shuffle gene-to-family labels among all inversion genes,
# recompute n_inversions_captured per family. Repeat 200 times.
shuffle_null <- function(dt, n_perm = 200) {
  set.seed(42)
  observed <- dt[nzchar(family),
                  .(n_inv = uniqueN(candidate_id)),
                  by = family]
  null_max <- replicate(n_perm, {
    shuffled <- copy(dt)
    shuffled[, family := sample(family)]
    shuffled[nzchar(family), uniqueN(candidate_id), by = family][, max(V1)]
  })
  list(observed = observed, null_max_distribution = null_max)
}
nul <- shuffle_null(all_inv_genes, n_perm = 200)
fam_summary[, exceeds_null_95pct :=
              n_inversions_captured > quantile(nul$null_max_distribution, 0.95)]

fwrite(fam_summary, file.path(CARGO_SYNTH_DIR, "family_convergence.tsv"), sep = "\t")
message("[C64] Wrote family_convergence.tsv (", nrow(fam_summary), " families)")

# ── 2. GO convergence ──
# For each GO term, how many inversions had it significantly enriched?
go_long <- rbindlist(lapply(candidates, function(cid) {
  ef <- file.path(CARGO_ENRICH_DIR, cid, "go_enrichment.tsv")
  if (!file.exists(ef)) return(NULL)
  e <- fread(ef)
  if (nrow(e) == 0) return(NULL)
  e[, candidate_id := cid]
  e
}), fill = TRUE)

if (nrow(go_long) > 0 && "fdr_bh" %in% names(go_long)) {
  sig_go <- go_long[!is.na(fdr_bh) & fdr_bh < 0.05 & background == "matched_collinear"]
  go_summary <- sig_go[, .(n_inversions_enriched = uniqueN(candidate_id),
                            min_fdr = min(fdr_bh, na.rm = TRUE),
                            mean_odds = mean(odds, na.rm = TRUE)),
                       by = term]
  setorder(go_summary, -n_inversions_enriched, min_fdr)
  fwrite(go_summary, file.path(CARGO_SYNTH_DIR, "go_convergence.tsv"), sep = "\t")
  message("[C64] Wrote go_convergence.tsv (", nrow(go_summary), " GO terms)")
} else {
  message("[C64] No GO enrichment results available for convergence analysis")
}

# ── 3. Top candidate genes (combined-evidence ranking) ──
# Pull cargo_burden_diff (pairwise), cargo_config_mi (interval_summary).
# If registry unavailable, read fallback files written by C61/C62.

read_burden_diff_one <- function(cid) {
  fp_fallback <- file.path(CARGO_BURDEN_DIR, cid, "burden_diff.tsv")
  # Try registry first
  if (!is.null(reg)) {
    g_REF <- paste0("inv_", cid, "_HOM_REF")
    g_INV <- paste0("inv_", cid, "_HOM_INV")
    rows <- tryCatch(reg$ask(who = c(g_REF, g_INV), what = "cargo_burden_diff"),
                     error = function(e) NULL)
    if (!is.null(rows) && length(rows) > 0 && nrow(rows) > 0) {
      f <- file.path(RESULTS_REGISTRY_DIR, rows$file[1])
      if (file.exists(f)) return(fread(f)[, candidate_id := cid])
    }
  }
  if (file.exists(fp_fallback)) return(fread(fp_fallback)[, candidate_id := cid])
  NULL
}

read_config_mi_one <- function(cid, kary) {
  fp_fallback <- file.path(CARGO_CONFIG_DIR, cid, paste0("config_mi__inv_", cid, "_", kary, ".tsv"))
  if (!is.null(reg)) {
    g <- paste0("inv_", cid, "_", kary)
    rows <- tryCatch(reg$ask(who = g, what = "cargo_config_mi"),
                     error = function(e) NULL)
    if (!is.null(rows) && length(rows) > 0 && nrow(rows) > 0) {
      f <- file.path(RESULTS_REGISTRY_DIR, rows$file[1])
      if (file.exists(f)) return(fread(f)[, `:=`(candidate_id = cid, karyotype = kary)])
    }
  }
  if (file.exists(fp_fallback)) return(fread(fp_fallback)[, `:=`(candidate_id = cid, karyotype = kary)])
  NULL
}

bd_all <- rbindlist(lapply(candidates, read_burden_diff_one), fill = TRUE)
mi_REF <- rbindlist(lapply(candidates, function(c) read_config_mi_one(c, "HOM_REF")), fill = TRUE)
mi_INV <- rbindlist(lapply(candidates, function(c) read_config_mi_one(c, "HOM_INV")), fill = TRUE)

# Combined ranking
top_rows <- list()
if (nrow(bd_all) > 0) {
  bd_all[, abs_delta := abs(as.numeric(delta_burden_INV_minus_REF))]
  bd_all[, perm_p := suppressWarnings(as.numeric(permutation_p))]
  # Join MI Z-score from the more interesting arrangement (HOM_INV)
  if (nrow(mi_INV) > 0) {
    bd_all <- merge(bd_all,
                    mi_INV[, .(candidate_id, gene_id,
                               mi_z_INV = suppressWarnings(as.numeric(mi_z_independent)),
                               mi_p_INV = suppressWarnings(as.numeric(mi_p_independent)))],
                    by = c("candidate_id", "gene_id"), all.x = TRUE)
  }
  # Combined rank: standardize each column → average rank
  rank_safe <- function(x) {
    x[!is.finite(x)] <- 0
    rank(-x, ties.method = "average")
  }
  bd_all[, score_burden_rank := rank_safe(abs_delta)]
  if ("mi_z_INV" %in% names(bd_all)) {
    bd_all[, score_mi_rank := rank_safe(mi_z_INV)]
  } else {
    bd_all[, score_mi_rank := NA_real_]
  }
  bd_all[, n_lof_total := as.integer(n_LoF_HOM_REF) + as.integer(n_LoF_HOM_INV)]
  bd_all[, score_lof_rank := rank_safe(n_lof_total)]
  bd_all[, combined_rank :=
           rowMeans(.SD, na.rm = TRUE),
         .SDcols = c("score_burden_rank", "score_mi_rank", "score_lof_rank")]
  top <- bd_all[order(combined_rank)]
  fwrite(top, file.path(CARGO_SYNTH_DIR, "top_candidate_genes.tsv"), sep = "\t")
  message("[C64] Wrote top_candidate_genes.tsv (", nrow(top), " genes)")
}

# ── 4. Cross-inversion summary ──
summary_dt <- diag[, .(candidate_id, chrom, start_bp, end_bp, length_bp,
                       n_HOM_REF, n_HOM_INV,
                       n_genes_inside, n_missense_sites_inside,
                       n_lof_sites_inside,
                       level_1_ok, level_2_ok)]
fwrite(summary_dt, file.path(CARGO_SYNTH_DIR, "cross_inversion_summary.tsv"), sep = "\t")
message("[C64] Wrote cross_inversion_summary.tsv")

message("[C64] Done")
