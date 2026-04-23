#!/usr/bin/env Rscript
# =============================================================================
# q09b_shelf_ld_check.R
# =============================================================================
# For a suspected inversion, test whether the alleles distinguishing Hom1 vs
# Hom2 at different positions in the shelf are in LD. For a single inversion,
# the inverted haplotype is one linked block, so "diagnostic SNPs" at the left
# edge should be highly correlated with diagnostic SNPs at the right edge.
#
# Method:
# 1. Load dosage matrix for shelf SNPs x samples (from BEAGLE stream or from
#    Q08's intermediate output if cached).
# 2. Use existing invgt assignments from Q07: Hom1, Het, Hom2 groups.
# 3. For each SNP, compute the "diagnostic score" = |p_Hom1 - p_Hom2|.
#    SNPs with score > 0.4 are "diagnostic" for the invgt axis.
# 4. Correlate dosage matrices across subsets of diagnostic SNPs from
#    different parts of the shelf. If sub-regions A and B are in the same
#    inversion, the sample-level dosages at A-SNPs and B-SNPs will be highly
#    correlated. If they're separate inversions, correlation will be low.
#
# Output: heatmap of between-window dosage correlations + summary table of
# within-vs-between-region LD.
#
# Usage:
#   Rscript q09b_shelf_ld_check.R \
#     --beagle path/to/beagle.gz --pos path/to/pos.fixed \
#     --invgt path/to/invgt_assignments.tsv \
#     --chrom C_gar_LG28 --shelf_start_mb 15 --shelf_end_mb 18 \
#     --out_pdf path/to/shelf_ld.pdf \
#     --out_tsv path/to/shelf_ld_summary.tsv
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(f, d = NULL) {
  i <- which(args == f)[1]
  if (is.na(i) || i >= length(args)) d else args[i + 1]
}

BEAGLE      <- get_arg("--beagle")
POS         <- get_arg("--pos")
INVGT_FILE  <- get_arg("--invgt")
CHROM       <- get_arg("--chrom")
SHELF_A     <- as.numeric(get_arg("--shelf_start_mb", "0"))
SHELF_B     <- as.numeric(get_arg("--shelf_end_mb",   "Inf"))
OUT_PDF     <- get_arg("--out_pdf")
OUT_TSV     <- get_arg("--out_tsv")
SAMPLE_LIST <- get_arg("--sample_list", NULL)
N_BINS      <- as.integer(get_arg("--n_bins", "20"))

message("[q09b] loading pos file")
pos <- fread(POS, header = FALSE, col.names = c("chrom", "bp", "major", "minor"),
             fill = TRUE)
pos[, idx := .I]
shelf_pos <- pos[chrom == CHROM & bp >= SHELF_A * 1e6 & bp <= SHELF_B * 1e6]
if (nrow(shelf_pos) == 0) stop("No shelf SNPs found")
message("[q09b] shelf SNPs: ", nrow(shelf_pos))

# Load sample master list
sample_master <- if (!is.null(SAMPLE_LIST) && file.exists(SAMPLE_LIST)) {
  readLines(SAMPLE_LIST)
} else character(0)

# Load invgt assignments
message("[q09b] loading invgt assignments")
inv <- fread(INVGT_FILE)
setnames(inv, tolower(names(inv)))
id_col <- intersect(c("sample", "sample_id", "id", "iid"), names(inv))[1]
grp_col <- intersect(c("invgt", "group", "assigned_group", "pca_group"), names(inv))[1]
if (is.null(id_col) || is.null(grp_col)) stop("invgt file needs sample and invgt columns")
setnames(inv, c(id_col, grp_col), c("sample", "invgt"))

# Stream BEAGLE for shelf SNPs
message("[q09b] streaming BEAGLE for ", nrow(shelf_pos), " shelf SNPs")
con <- gzfile(BEAGLE, "r")
hdr <- readLines(con, n = 1L)
hdr_flds <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
n_samples <- (length(hdr_flds) - 3L) / 3L

dosage_mat <- matrix(NA_real_, nrow = nrow(shelf_pos), ncol = n_samples)
if (length(sample_master) >= n_samples) {
  colnames(dosage_mat) <- sample_master[seq_len(n_samples)]
} else {
  colnames(dosage_mat) <- paste0("Ind", seq_len(n_samples) - 1L)
}

target_env <- new.env(hash = TRUE, size = nrow(shelf_pos) * 2L)
for (i in seq_len(nrow(shelf_pos))) assign(as.character(shelf_pos$idx[i]), i, envir = target_env)

row_idx <- 0L; hits <- 0L; chunk_size <- 20000L
max_target <- max(shelf_pos$idx)
repeat {
  lines <- readLines(con, n = chunk_size)
  if (!length(lines)) break
  for (ln in lines) {
    row_idx <- row_idx + 1L
    hi <- target_env[[as.character(row_idx)]]
    if (!is.null(hi)) {
      flds <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      gl <- as.numeric(flds[-(1:3)])
      m <- matrix(gl, nrow = n_samples, byrow = TRUE)
      dosage_mat[hi, ] <- m[, 2] + 2 * m[, 3]
      hits <- hits + 1L
      if (hits %% 2000L == 0L) message("  ", hits, " / ", nrow(shelf_pos))
    }
    if (row_idx > max_target) break
  }
  if (row_idx > max_target) break
}
close(con)
message("[q09b] extracted ", hits, " shelf SNPs")

# Match invgt assignments to columns
samp_ok <- colnames(dosage_mat) %in% inv$sample
if (!any(samp_ok)) stop("No samples match between BEAGLE header and invgt file")
dosage_mat <- dosage_mat[, samp_ok, drop = FALSE]
inv_lookup <- setNames(inv$invgt, inv$sample)
sample_invgt <- inv_lookup[colnames(dosage_mat)]
message("[q09b] per-group sizes: ",
        paste(names(table(sample_invgt)), table(sample_invgt), sep = "=", collapse = ", "))

hom1_idx <- which(sample_invgt == "Hom1")
hom2_idx <- which(sample_invgt == "Hom2")
het_idx  <- which(sample_invgt == "Het")

# =============================================================================
# Per-SNP allele frequency and diagnostic score
# =============================================================================
freq_in_group <- function(mat, cols) {
  if (!length(cols)) return(rep(NA_real_, nrow(mat)))
  rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE) / 2
}
p_hom1 <- freq_in_group(dosage_mat, hom1_idx)
p_hom2 <- freq_in_group(dosage_mat, hom2_idx)
p_het  <- freq_in_group(dosage_mat, het_idx)

snp_dt <- data.table(
  idx = shelf_pos$idx,
  bp  = shelf_pos$bp,
  mb  = shelf_pos$bp / 1e6,
  p_hom1 = p_hom1,
  p_hom2 = p_hom2,
  p_het  = p_het,
  diagnostic = abs(p_hom1 - p_hom2)
)
message("[q09b] diagnostic score > 0.4: ",
        sum(snp_dt$diagnostic > 0.4, na.rm = TRUE), " / ", nrow(snp_dt))

# =============================================================================
# Per-sample inversion-haplotype score PER BIN
# Partition the shelf into N_BINS bins. Within each bin, score each sample
# as mean dosage of diagnostic SNPs weighted by p_hom1 - p_hom2 sign. This
# gives a per-sample "which arrangement in this bin" profile.
# Then: correlate the per-sample bin profiles across all pairs of bins.
# =============================================================================
snp_dt[, bin := cut(mb, breaks = seq(SHELF_A, SHELF_B, length.out = N_BINS + 1),
                    include.lowest = TRUE, labels = FALSE)]

score_per_bin <- matrix(NA_real_, nrow = ncol(dosage_mat), ncol = N_BINS)
rownames(score_per_bin) <- colnames(dosage_mat)

for (b in seq_len(N_BINS)) {
  snps_in_bin <- snp_dt[bin == b & diagnostic > 0.3, .(idx, p_hom1, p_hom2)]
  if (nrow(snps_in_bin) < 5) next
  # Polarize: for each SNP, the arrangement signal is the difference
  #   sign(p_hom2 - p_hom1) * (dosage - mean)
  # Sum it across all diagnostic SNPs in the bin and normalize.
  row_slice <- match(snps_in_bin$idx, shelf_pos$idx)
  sub <- dosage_mat[row_slice, , drop = FALSE]
  pol <- sign(snps_in_bin$p_hom2 - snps_in_bin$p_hom1)
  # Rescale dosage to -1..1 with sign of polarity
  scored <- sweep(sub, 1, pol, `*`)
  score_per_bin[, b] <- colMeans(scored, na.rm = TRUE)
}

# =============================================================================
# Between-bin correlation matrix — the headline diagnostic
# If all bins correlate > 0.8 with each other -> one inversion
# If bins fall into two blocks with high within-block correlation and low
# between-block -> two inversions with different sample assignments
# =============================================================================
cor_mat <- suppressWarnings(cor(score_per_bin, use = "pairwise.complete.obs"))

# Render heatmap
message("[q09b] rendering correlation heatmap")
cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("bin_a", "bin_b", "r"))
cor_dt[, bin_a := as.integer(bin_a)]
cor_dt[, bin_b := as.integer(bin_b)]
# Bin centers for axis labels
bin_edges <- seq(SHELF_A, SHELF_B, length.out = N_BINS + 1)
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2
cor_dt[, mb_a := bin_centers[bin_a]]
cor_dt[, mb_b := bin_centers[bin_b]]

p_cor <- ggplot(cor_dt, aes(mb_a, mb_b, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "#c0504d", mid = "white", high = "#1f4e79",
                       midpoint = 0, limits = c(-1, 1),
                       name = "sample-\narrangement\ncorr") +
  coord_fixed() +
  labs(x = paste0(CHROM, " (Mb)"), y = paste0(CHROM, " (Mb)"),
       title = sprintf("Between-bin sample-arrangement correlation  (%d bins)", N_BINS),
       subtitle = "single inversion = all blue / light red\ntwo inversions = two blue blocks, off-diagonal white/red") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#555e69"))

# Plot of per-SNP diagnostic score across the shelf too
p_diag <- ggplot(snp_dt, aes(mb, diagnostic)) +
  geom_point(alpha = 0.3, size = 0.4, color = "#1f4e79") +
  geom_hline(yintercept = c(0.3, 0.6), linetype = "dashed", color = "#c0504d") +
  labs(x = paste0(CHROM, " (Mb)"), y = "|p_Hom1 − p_Hom2|",
       title = "Per-SNP diagnostic score",
       subtitle = "dashed lines: 0.3 (loose) / 0.6 (strong) diagnostic threshold") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# Summary table: within-region mean correlation vs between-region
# Heuristic: split shelf in half, compare within-halves to across-halves
mid <- N_BINS %/% 2L
a_bins <- 1:mid
b_bins <- (mid + 1):N_BINS
within_a  <- mean(cor_mat[a_bins, a_bins][upper.tri(cor_mat[a_bins, a_bins])], na.rm = TRUE)
within_b  <- mean(cor_mat[b_bins, b_bins][upper.tri(cor_mat[b_bins, b_bins])], na.rm = TRUE)
between   <- mean(cor_mat[a_bins, b_bins], na.rm = TRUE)
overall   <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)

summary_dt <- data.table(
  chrom = CHROM,
  shelf_start_mb = SHELF_A,
  shelf_end_mb = SHELF_B,
  n_bins = N_BINS,
  n_snps_shelf = nrow(snp_dt),
  n_diagnostic_03 = sum(snp_dt$diagnostic > 0.3, na.rm = TRUE),
  n_diagnostic_06 = sum(snp_dt$diagnostic > 0.6, na.rm = TRUE),
  within_first_half_corr = within_a,
  within_second_half_corr = within_b,
  between_halves_corr = between,
  overall_mean_corr = overall,
  verdict = ifelse(between > 0.7, "SINGLE_INVERSION",
            ifelse(between < 0.3, "LIKELY_MULTIPLE_ARRANGEMENTS",
                                  "AMBIGUOUS_OR_PARTIAL_LD"))
)

fwrite(summary_dt, OUT_TSV, sep = "\t")
cat("\n=========  SHELF LD SUMMARY  =========\n")
for (col in names(summary_dt)) {
  cat(sprintf("  %-30s: %s\n", col, summary_dt[[col]]))
}
cat("=====================================\n\n")

# Combined PDF
pdf(OUT_PDF, width = 10, height = 12)
print(p_diag)
print(p_cor)
dev.off()
message("[q09b] wrote ", OUT_PDF)
message("[q09b] wrote ", OUT_TSV)
