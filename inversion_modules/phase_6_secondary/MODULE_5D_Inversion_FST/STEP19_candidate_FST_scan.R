#!/usr/bin/env Rscript

# =============================================================================
# STEP19_candidate_FST_scan.R  (v2)
#
# Candidate-region FST scans using STEP17c-exported contrast groups.
# Uses the SAME grouping layer as the LD module.
#
# Reads: <group_dir>/candidate_<id>.contrast_groups.tsv
# Uses: provisional_major_state column for sample assignment
#
# Computes Hudson FST from dosage for:
#   FULL_A vs FULL_B  (main contrast)
#   FULL_A vs HALF
#   FULL_B vs HALF
#
# Usage:
#   Rscript STEP19_candidate_FST_scan.R \
#     <group_dir> <candidate_table> <dosage_dir> <outdir> \
#     [candidate_id=all] [flank_frac=0.10] [win_snps=50] [min_per_group=5]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript STEP19_... <group_dir> <candidate_table> <dosage_dir> <outdir> ...")
}

group_dir      <- args[1]
candidate_file <- args[2]
dosage_dir     <- args[3]
outdir         <- args[4]
cid_filter     <- if (length(args) >= 5 && args[5] != "all") args[5] else NA_character_
flank_frac     <- if (length(args) >= 6) as.numeric(args[6]) else 0.10
win_snps       <- if (length(args) >= 7) as.integer(args[7]) else 50L
min_per_group  <- if (length(args) >= 8) as.integer(args[8]) else 5L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cand <- fread(candidate_file)
stopifnot(all(c("candidate_id", "chrom", "start_bp", "end_bp") %in% names(cand)))

if (!is.na(cid_filter)) cand <- cand[as.character(candidate_id) == cid_filter]
if (nrow(cand) == 0) stop("No candidates")
message("[INFO] Candidates: ", nrow(cand))

# Hudson FST from dosage
hudson_fst_vec <- function(dos_mat, idx1, idx2) {
  p1 <- rowMeans(dos_mat[, idx1, drop = FALSE], na.rm = TRUE) / 2
  p2 <- rowMeans(dos_mat[, idx2, drop = FALSE], na.rm = TRUE) / 2
  num <- (p1 - p2)^2
  den <- p1 * (1 - p2) + p2 * (1 - p1)
  den[den == 0] <- NA_real_
  fst <- num / den
  fst[!is.finite(fst)] <- NA_real_
  fst
}

all_fst <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- as.character(row$candidate_id)
  chr <- row$chrom
  inv_start <- as.numeric(row$start_bp); inv_end <- as.numeric(row$end_bp)
  inv_len <- inv_end - inv_start

  # ── Read STEP17c contrast groups ───────────────────────────────────
  grp_file <- file.path(group_dir, paste0("candidate_", cid, ".contrast_groups.tsv"))
  if (!file.exists(grp_file)) {
    message("[WARN] No contrast groups for candidate ", cid); next
  }

  grp <- fread(grp_file)
  if (!("provisional_major_state" %in% names(grp))) {
    message("[WARN] Missing provisional_major_state in ", grp_file); next
  }

  samp_A <- grp[provisional_major_state == "FULL_A", sample]
  samp_B <- grp[provisional_major_state == "FULL_B", sample]
  samp_H <- grp[provisional_major_state == "HALF", sample]

  if (length(samp_A) < min_per_group || length(samp_B) < min_per_group) {
    message("[INFO] Candidate ", cid, ": too few FULL_A (", length(samp_A),
            ") or FULL_B (", length(samp_B), ") for FST")
    next
  }

  # ── Read dosage ──────────────────────────────────────────────────
  dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(dosage_dir, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[WARN] Missing dosage for ", chr); next
  }

  dos <- fread(dos_file); sites <- fread(sites_file)

  flank_bp <- max(50000, inv_len * flank_frac)
  reg_start <- max(1, inv_start - flank_bp)
  reg_end <- inv_end + flank_bp

  keep <- which(sites$chrom == chr & sites$pos >= reg_start & sites$pos <= reg_end)
  if (length(keep) < win_snps) next

  dos_reg <- dos[keep]; sites_reg <- sites[keep]
  sample_cols <- setdiff(names(dos_reg), "marker")

  idx_A <- which(sample_cols %in% samp_A)
  idx_B <- which(sample_cols %in% samp_B)
  idx_H <- which(sample_cols %in% samp_H)

  X <- as.matrix(dos_reg[, ..sample_cols]); storage.mode(X) <- "double"

  comparisons <- list()
  if (length(idx_A) >= min_per_group && length(idx_B) >= min_per_group)
    comparisons[["FULL_A_vs_FULL_B"]] <- list(idx_A, idx_B)
  if (length(idx_A) >= min_per_group && length(idx_H) >= min_per_group)
    comparisons[["FULL_A_vs_HALF"]] <- list(idx_A, idx_H)
  if (length(idx_B) >= min_per_group && length(idx_H) >= min_per_group)
    comparisons[["FULL_B_vs_HALF"]] <- list(idx_B, idx_H)

  if (length(comparisons) == 0) next

  n_win <- floor(nrow(X) / win_snps)
  if (n_win < 1) next

  for (comp_name in names(comparisons)) {
    i1 <- comparisons[[comp_name]][[1]]; i2 <- comparisons[[comp_name]][[2]]
    fst_vec <- hudson_fst_vec(X, i1, i2)

    for (w in seq_len(n_win)) {
      s1 <- (w-1L)*win_snps + 1L; s2 <- w*win_snps
      all_fst[[length(all_fst)+1]] <- data.table(
        candidate_id = as.integer(cid), chrom = chr,
        comparison = comp_name, window = w,
        win_center = mean(sites_reg$pos[s1:s2]),
        fst_mean = round(mean(fst_vec[s1:s2], na.rm=TRUE), 6),
        in_candidate = mean(sites_reg$pos[s1:s2]) >= inv_start & mean(sites_reg$pos[s1:s2]) <= inv_end,
        n1 = length(i1), n2 = length(i2)
      )
    }
  }

  message("[INFO] Candidate ", cid, ": FST for ", length(comparisons), " comparisons (",
          paste(names(comparisons), collapse=", "), ")")
}

if (length(all_fst) == 0) { message("[WARN] No FST results"); quit("no", status=0) }

fst_dt <- rbindlist(all_fst, fill = TRUE)
fwrite(fst_dt, file.path(outdir, "inversion_candidate_fst_scan.tsv.gz"), sep = "\t")

# Plots
for (cid in unique(fst_dt$candidate_id)) {
  sub <- fst_dt[candidate_id == cid]; chr <- sub$chrom[1]
  cr <- cand[candidate_id == cid]
  inv_s <- as.numeric(cr$start_bp[1]); inv_e <- as.numeric(cr$end_bp[1])

  p <- ggplot(sub, aes(x = win_center / 1e6, y = fst_mean, color = comparison)) +
    annotate("rect", xmin=inv_s/1e6, xmax=inv_e/1e6, ymin=-Inf, ymax=Inf, alpha=0.15, fill="grey60") +
    geom_line(linewidth=0.7) + geom_point(size=0.8, alpha=0.6) +
    facet_wrap(~comparison, ncol=1, scales="free_y") +
    theme_bw(base_size=11) +
    labs(title=paste0("Candidate ", cid, " — FST (", chr, ")"),
         subtitle="Hudson FST from dosage | groups from STEP17c",
         x="Position (Mb)", y=expression(F[ST])) +
    theme(legend.position="none", strip.text=element_text(size=9))

  ht <- 3 * length(unique(sub$comparison)) + 1
  ggsave(file.path(outdir, paste0("fst_scan_candidate_", cid, "_", chr, ".pdf")), p, width=8, height=ht)
  ggsave(file.path(outdir, paste0("fst_scan_candidate_", cid, "_", chr, ".png")), p, width=8, height=ht, dpi=400)
}

message("[DONE] FST scan complete")
