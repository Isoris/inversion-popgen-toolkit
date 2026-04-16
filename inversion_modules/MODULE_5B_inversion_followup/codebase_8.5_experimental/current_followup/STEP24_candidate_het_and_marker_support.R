#!/usr/bin/env Rscript

# =============================================================================
# STEP24_candidate_het_and_marker_support.R
#
# Heterozygosity support layer + marker density / support metrics.
#
# For each candidate:
#   1. Summarize regional heterozygosity by coarse group and subcluster
#   2. Compute marker density, inter-SNP gap distribution, support tier
#   3. Count overlapping genes from GFF
#
# Inputs:
#   - STEP21 rotated PCA tables (have regional_het)
#   - STEP22 subcluster tables
#   - STEP12 regional sites tables
#   - GFF annotation
#
# Outputs per candidate:
#   - candidate_group_het_summary.tsv
#   - candidate_marker_support.tsv
#   - candidate_gene_overlap.tsv
#
# Usage:
#   Rscript STEP24_candidate_het_and_marker_support.R <config.R> [cid=all]
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

# ── Load GFF for gene overlap ───────────────────────────────────────────────
gff_dt <- NULL
if (file.exists(GFF_FILE)) {
  gff_raw <- fread(
    GFF_FILE,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    comment.char = "#"
  )

  if (ncol(gff_raw) < 9) {
    stop("GFF file does not have at least 9 columns after removing comments: ", GFF_FILE)
  }

  gff_raw <- gff_raw[, .(chrom = V1, feature = V3, start = V4, end = V5, attr = V9)]


  gff_dt <- gff_raw[feature %in% c("gene", "mRNA")]
  # Extract gene name/ID
  gff_dt[, gene_id := sub(".*ID=([^;]+).*", "\\1", attr)]
  gff_dt[, gene_name := fifelse(
    grepl("Name=", attr),
    sub(".*Name=([^;]+).*", "\\1", attr),
    gene_id
  )]
  message("[INFO] Loaded GFF: ", nrow(gff_dt), " gene/mRNA features")
}

# ── Heterozygosity summary function ──────────────────────────────────────────
het_summary <- function(het_vals, group_label) {
  data.table(
    group_label = group_label,
    n = length(het_vals),
    mean_het = round(mean(het_vals, na.rm = TRUE), 6),
    median_het = round(median(het_vals, na.rm = TRUE), 6),
    sd_het = round(sd(het_vals, na.rm = TRUE), 6),
    iqr_het = round(IQR(het_vals, na.rm = TRUE), 6),
    min_het = round(min(het_vals, na.rm = TRUE), 6),
    max_het = round(max(het_vals, na.rm = TRUE), 6)
  )
}

# ── Main loop ────────────────────────────────────────────────────────────────
all_het <- list()
all_marker <- list()
all_genes <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end <- as.numeric(row$end_bp)
  c_len <- c_end - c_start

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  sub_file <- file.path(cand_dir, "candidate_subclusters.tsv")

  if (!file.exists(rot_file)) next
  rot <- fread(rot_file)
  sub <- if (file.exists(sub_file)) fread(sub_file) else NULL

  # ── Heterozygosity by coarse group ─────────────────────────────────────
  het_results <- list()
  if ("regional_het" %in% names(rot) && sum(is.finite(rot$regional_het)) > 3) {
    for (g in sort(unique(rot$coarse_group_refined))) {
      idx <- rot$coarse_group_refined == g & is.finite(rot$regional_het)
      if (sum(idx) >= 2) {
        h <- het_summary(rot$regional_het[idx], g)
        h[, candidate_id := cid]
        h[, group_level := "coarse_group"]
        het_results[[length(het_results) + 1]] <- h
      }
    }

    # By subcluster
    if (!is.null(sub) && "subcluster_label" %in% names(sub)) {
      merged_sub <- merge(rot[, .(sample, regional_het)],
                          sub[, .(sample, subcluster_label)],
                          by = "sample")
      for (sl in unique(merged_sub$subcluster_label)) {
        idx <- merged_sub$subcluster_label == sl & is.finite(merged_sub$regional_het)
        if (sum(idx) >= 2) {
          h <- het_summary(merged_sub$regional_het[idx], sl)
          h[, candidate_id := cid]
          h[, group_level := "subcluster"]
          het_results[[length(het_results) + 1]] <- h
        }
      }
    }
  }

  if (length(het_results) > 0) {
    het_dt <- rbindlist(het_results, fill = TRUE)
    fwrite(het_dt, file.path(cand_dir, "candidate_group_het_summary.tsv"), sep = "\t")
    all_het[[length(all_het) + 1]] <- het_dt
  }

  # ── Marker density / support ───────────────────────────────────────────
  sites_file <- file.path(STEP12_DIR,
                          paste0("STEP12_", chr, ".candidate_", cid, ".regional_sites.tsv.gz"))
  if (!file.exists(sites_file)) {
    # Try dosage dir
    sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  }

  if (file.exists(sites_file)) {
    sites <- fread(sites_file)
    if ("pos" %in% names(sites)) {
      # Filter to candidate region
      if ("chrom" %in% names(sites)) {
        sites_reg <- sites[chrom == chr & pos >= c_start & pos <= c_end]
      } else {
        sites_reg <- sites[pos >= c_start & pos <= c_end]
      }

      n_snps <- nrow(sites_reg)
      if (n_snps >= 2) {
        gaps <- diff(sort(sites_reg$pos))
        snps_per_kb <- n_snps / (c_len / 1000)
        bp_per_snp <- c_len / n_snps

        marker_dt <- data.table(
          candidate_id = cid,
          chrom = chr,
          start_bp = c_start,
          end_bp = c_end,
          length_bp = c_len,
          n_snps = n_snps,
          snps_per_kb = round(snps_per_kb, 2),
          bp_per_snp = round(bp_per_snp, 0),
          median_gap_bp = round(median(gaps), 0),
          p90_gap_bp = round(quantile(gaps, 0.9), 0),
          max_gap_bp = max(gaps),
          support_tier = fifelse(
            snps_per_kb >= 5 & max(gaps) < c_len * 0.1, "high",
            fifelse(snps_per_kb >= 2, "medium", "low")
          )
        )

        fwrite(marker_dt, file.path(cand_dir, "candidate_marker_support.tsv"), sep = "\t")
        all_marker[[length(all_marker) + 1]] <- marker_dt

        # Save gap distribution for plotting
        gap_dt <- data.table(candidate_id = cid, gap_bp = gaps,
                             gap_start = sort(sites_reg$pos)[-length(sort(sites_reg$pos))])
        fwrite(gap_dt, file.path(cand_dir, "candidate_marker_gaps.tsv.gz"), sep = "\t")
      }
    }
  }

  # ── Gene overlap ──────────────────────────────────────────────────────
  if (!is.null(gff_dt)) {
    genes_overlap <- gff_dt[chrom == chr & feature == "gene" &
                            start <= c_end & end >= c_start]
    n_genes <- nrow(genes_overlap)

    # Also check flanks (20% each side)
    flank <- max(c_len * 0.2, 50000)
    genes_left <- gff_dt[chrom == chr & feature == "gene" &
                         start <= c_start & end >= (c_start - flank)]
    genes_right <- gff_dt[chrom == chr & feature == "gene" &
                          start <= (c_end + flank) & end >= c_end]

    gene_dt <- data.table(
      candidate_id = cid,
      chrom = chr,
      n_genes_overlap = n_genes,
      n_genes_left_flank = nrow(genes_left),
      n_genes_right_flank = nrow(genes_right),
      top_genes = paste(head(genes_overlap$gene_name, 5), collapse = ", ")
    )

    fwrite(gene_dt, file.path(cand_dir, "candidate_gene_overlap.tsv"), sep = "\t")
    all_genes[[length(all_genes) + 1]] <- gene_dt
  }

  message("[INFO] Candidate ", cid, " het/marker/gene analysis done")
}

# ── Global outputs ───────────────────────────────────────────────────────────
if (length(all_het) > 0) {
  fwrite(rbindlist(all_het, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_het_summary.tsv.gz"), sep = "\t")
}
if (length(all_marker) > 0) {
  fwrite(rbindlist(all_marker, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_marker_support.tsv.gz"), sep = "\t")
}
if (length(all_genes) > 0) {
  fwrite(rbindlist(all_genes, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_gene_overlap.tsv.gz"), sep = "\t")
}

message("[DONE] STEP24 complete")
