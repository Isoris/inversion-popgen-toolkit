#!/usr/bin/env Rscript

# =============================================================================
# STEP_C65_per_inversion_profile.R — 6E figures
#
# Per-candidate "cargo profile" figure. One PDF per candidate:
#   Panel A: gene density along inversion (lollipop, colored by family)
#   Panel B: per-site FST between arrangements (manhattan, colored by VESM class)
#   Panel C: per-arrangement burden — top 20 genes ranked by |delta_burden|
#   Panel D: configuration-MI Z-score per gene (HOM_REF vs HOM_INV)
#
# Inputs (per candidate):
#   ${CARGO_INVENTORY_DIR}/<cid>/genes.tsv
#   FST table (registry: cargo_fst_per_site, or fallback file)
#   burden_diff table (registry: cargo_burden_diff, or fallback)
#   config_mi tables (registry: cargo_config_mi, or fallback)
#
# Output: ${CARGO_FIG_DIR}/<cid>__profile.pdf
#
# Usage:
#   Rscript STEP_C65_per_inversion_profile.R [candidate_id|all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

CARGO_INVENTORY_DIR <- Sys.getenv("CARGO_INVENTORY_DIR")
CARGO_BURDEN_DIR    <- Sys.getenv("CARGO_BURDEN_DIR")
CARGO_CONFIG_DIR    <- Sys.getenv("CARGO_CONFIG_DIR")
CARGO_FIG_DIR       <- Sys.getenv("CARGO_FIG_DIR")
RESULTS_REGISTRY_DIR <- Sys.getenv("RESULTS_REGISTRY_DIR")
REGISTRY_LOADER_R    <- Sys.getenv("REGISTRY_LOADER_R")

dir.create(CARGO_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

reg <- NULL
if (nzchar(REGISTRY_LOADER_R) && file.exists(REGISTRY_LOADER_R)) {
  source(REGISTRY_LOADER_R)
  reg <- tryCatch(load_registry(), error = function(e) NULL)
}

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[1] else "all"

# ── Loaders with registry-then-fallback pattern ──
load_via_registry_or_file <- function(cid, stat_name, fallback_file) {
  if (!is.null(reg)) {
    rows <- tryCatch(reg$ask(what = stat_name,
                              who = paste0("inv_", cid, c("_HOM_REF", "_HOM_INV"))),
                     error = function(e) NULL)
    if (!is.null(rows) && nrow(rows) > 0) {
      f <- file.path(RESULTS_REGISTRY_DIR, rows$file[1])
      if (file.exists(f)) return(fread(f))
    }
  }
  if (file.exists(fallback_file)) return(fread(fallback_file))
  NULL
}

# ── Panel builders ──

panel_A_gene_track <- function(genes_dt, cand_chrom, start_bp, end_bp) {
  if (nrow(genes_dt) == 0) return(NULL)
  dt <- copy(genes_dt)
  dt[, mid := (start + end) / 2]
  dt[, fam_simple := ifelse(nzchar(family),
                              substr(family, 1, 20),
                              "unknown")]
  ggplot(dt, aes(x = mid, y = 1, color = fam_simple)) +
    geom_segment(aes(xend = mid, yend = 0), alpha = 0.6) +
    geom_point(size = 1.2) +
    scale_y_continuous(limits = c(0, 1.1), breaks = NULL) +
    scale_x_continuous(name = paste0(cand_chrom, " position (bp)"),
                       limits = c(start_bp, end_bp),
                       labels = scales::comma) +
    labs(title = "A. Gene track inside inversion",
         color = "family") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.text = element_text(size = 7),
          panel.grid.major.y = element_blank())
}

panel_B_fst <- function(fst_dt, cand_chrom, start_bp, end_bp) {
  if (is.null(fst_dt) || nrow(fst_dt) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                label = "no FST data") +
                       labs(title = "B. per-site FST") +
                       theme_void())
  }
  dt <- copy(fst_dt)
  dt[, pos := as.integer(sub("^[^:]+:([0-9]+):.*", "\\1", var_key))]
  dt[, fst := suppressWarnings(as.numeric(fst_hudson))]
  dt[, vesm := suppressWarnings(as.numeric(vesm_llr))]
  dt[, vesm_class := fcase(
    is.na(vesm), "unknown",
    vesm < -7,    "likely_damaging",
    vesm < -3,    "possibly_damaging",
    default       = "benign")]
  ggplot(dt[!is.na(fst)], aes(x = pos, y = fst, color = vesm_class)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c(likely_damaging = "#d62728",
                                    possibly_damaging = "#ff7f0e",
                                    benign = "#1f77b4",
                                    unknown = "grey60"),
                       drop = FALSE) +
    scale_x_continuous(name = paste0(cand_chrom, " position (bp)"),
                       limits = c(start_bp, end_bp),
                       labels = scales::comma) +
    scale_y_continuous(name = "Hudson FST (HOM_REF vs HOM_INV)",
                       limits = c(0, 1)) +
    labs(title = "B. Per-site missense FST between arrangements",
         color = "VESM class") +
    theme_minimal()
}

panel_C_burden_top20 <- function(burden_dt) {
  if (is.null(burden_dt) || nrow(burden_dt) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                label = "no burden data") +
                       theme_void())
  }
  dt <- copy(burden_dt)
  dt[, dB := suppressWarnings(as.numeric(delta_burden_INV_minus_REF))]
  dt[, abs_dB := abs(dB)]
  top <- dt[order(-abs_dB)][1:min(20, nrow(dt))]
  top[, gene_id := factor(gene_id, levels = rev(top$gene_id))]
  ggplot(top, aes(x = dB, y = gene_id, fill = dB > 0)) +
    geom_col() +
    scale_fill_manual(values = c(`TRUE` = "#d62728", `FALSE` = "#1f77b4"),
                       labels = c(`TRUE` = "higher on HOM_INV",
                                   `FALSE` = "higher on HOM_REF"),
                       name = "") +
    labs(title = "C. Top 20 genes by burden asymmetry",
         x = "Δ VESM burden  (HOM_INV − HOM_REF)",
         y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7),
          legend.position = "bottom")
}

panel_D_config_mi <- function(mi_REF, mi_INV) {
  has_R <- !is.null(mi_REF) && nrow(mi_REF) > 0
  has_I <- !is.null(mi_INV) && nrow(mi_INV) > 0
  if (!has_R && !has_I) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                label = "no config-MI data") +
                       theme_void())
  }
  if (!has_R) mi_REF <- mi_INV[0]
  if (!has_I) mi_INV <- mi_REF[0]
  d <- merge(
    mi_REF[, .(gene_id, z_REF = suppressWarnings(as.numeric(mi_z_independent)))],
    mi_INV[, .(gene_id, z_INV = suppressWarnings(as.numeric(mi_z_independent)))],
    by = "gene_id", all = TRUE)
  d[is.na(z_REF), z_REF := 0]
  d[is.na(z_INV), z_INV := 0]
  d[, asym := abs(z_INV - z_REF)]
  ggplot(d, aes(x = z_REF, y = z_INV, color = asym)) +
    geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
    geom_hline(yintercept = 2, color = "grey80") +
    geom_vline(xintercept = 2, color = "grey80") +
    geom_point(alpha = 0.8) +
    scale_color_viridis_c(option = "magma") +
    labs(title = "D. Configuration-MI Z-score per gene",
         x = "Z(MI) on HOM_REF",
         y = "Z(MI) on HOM_INV",
         color = "|asym|") +
    theme_minimal()
}

# ── Assemble per candidate ──
process_one <- function(cid) {
  ginfo <- file.path(CARGO_INVENTORY_DIR, cid, "genes.tsv")
  if (!file.exists(ginfo)) {
    message("[C65] No inventory for ", cid); return(invisible())
  }
  gdt <- fread(ginfo)
  if (nrow(gdt) == 0) {
    message("[C65] Empty gene table for ", cid); return(invisible())
  }
  cand_chrom <- gdt$chrom[1]
  start_bp <- min(gdt$start)
  end_bp   <- max(gdt$end)

  fst_dt    <- load_via_registry_or_file(
    cid, "cargo_fst_per_site",
    file.path(CARGO_BURDEN_DIR, cid, "fst_per_site.tsv"))
  burden_dt <- load_via_registry_or_file(
    cid, "cargo_burden_diff",
    file.path(CARGO_BURDEN_DIR, cid, "burden_diff.tsv"))
  mi_REF <- load_via_registry_or_file(
    cid, "cargo_config_mi",
    file.path(CARGO_CONFIG_DIR, cid, paste0("config_mi__inv_", cid, "_HOM_REF.tsv")))
  mi_INV <- load_via_registry_or_file(
    cid, "cargo_config_mi",
    file.path(CARGO_CONFIG_DIR, cid, paste0("config_mi__inv_", cid, "_HOM_INV.tsv")))

  pA <- panel_A_gene_track(gdt, cand_chrom, start_bp, end_bp)
  pB <- panel_B_fst(fst_dt, cand_chrom, start_bp, end_bp)
  pC <- panel_C_burden_top20(burden_dt)
  pD <- panel_D_config_mi(mi_REF, mi_INV)

  layout <- (pA / pB) | (pC / pD)
  combined <- layout +
    plot_annotation(title = paste0("Cargo profile — ", cid,
                                     "  (", cand_chrom, ":",
                                     format(start_bp, big.mark = ","), "–",
                                     format(end_bp, big.mark = ","), ")"))

  out <- file.path(CARGO_FIG_DIR, paste0(cid, "__profile.pdf"))
  ggsave(out, combined, width = 14, height = 10, device = cairo_pdf)
  message("[C65] Wrote ", out)
}

if (target == "all") {
  diag <- fread(file.path(CARGO_INVENTORY_DIR, "diagnostic_table.tsv"))
  cids <- diag[n_genes_inside > 0, candidate_id]
} else {
  cids <- target
}
message("[C65] Generating profiles for ", length(cids), " candidate(s)")
for (cid in cids) process_one(cid)
message("[C65] Done")
