#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop(paste(
    "Usage: Rscript STEP14_candidate_followup.R",
    "<candidate_regions.tsv.gz>",
    "<candidate_id>",
    "<step12_regional_pca_samples.tsv.gz>",
    "<sample_tP_windows.tsv.gz>",
    "<outprefix>",
    "<chrom_sizes.tsv>",
    "<family_table_or_NONE>",
    "[flank_frac=0.10]",
    "[min_flank_bp=20000]"
  ))
}

candidate_file   <- args[1]
candidate_id_val <- args[2]
pcs_file         <- args[3]
tp_file          <- args[4]
outprefix        <- args[5]
chrom_sizes_file <- args[6]
family_file      <- args[7]
flank_frac       <- if (length(args) >= 8) as.numeric(args[8]) else 0.10
min_flank_bp     <- if (length(args) >= 9) as.numeric(args[9]) else 20000

# -----------------------------
# Read candidate
# -----------------------------
cand <- fread(candidate_file)
stopifnot("candidate_id" %in% names(cand))
row <- cand[as.character(candidate_id) == as.character(candidate_id_val)]

if (nrow(row) != 1) {
  stop("Expected exactly 1 candidate row for candidate_id=", candidate_id_val,
       " but found ", nrow(row))
}

cand_chrom <- row$chrom[1]
inv_start  <- as.numeric(row$start_bp[1])
inv_end    <- as.numeric(row$end_bp[1])
inv_len    <- inv_end - inv_start + 1

# -----------------------------
# Chrom sizes
# -----------------------------
chr_sizes <- fread(chrom_sizes_file)
chr_col  <- names(chr_sizes)[1]
size_col <- names(chr_sizes)[2]
setnames(chr_sizes, c(chr_col, size_col), c("chrom", "chrom_size"))

chr_size <- chr_sizes[chrom == cand_chrom, chrom_size]
if (length(chr_size) != 1) {
  stop("Could not find chromosome size for ", cand_chrom)
}

flank_bp   <- max(min_flank_bp, round(inv_len * flank_frac))
zoom_start <- max(1, inv_start - flank_bp)
zoom_end   <- min(chr_size, inv_end + flank_bp)

# -----------------------------
# Read STEP12 PCA samples
# -----------------------------
pcs <- fread(pcs_file)
req_pcs <- c("sample", "PC1", "PC2", "group_label", "regional_het")
miss_pcs <- setdiff(req_pcs, names(pcs))
if (length(miss_pcs) > 0) {
  stop("STEP12 pcs file missing columns: ", paste(miss_pcs, collapse = ", "))
}

# -----------------------------
# Optional family / ancestry table
# -----------------------------
if (!identical(family_file, "NONE")) {
  family_dt <- fread(family_file)
  if (!("sample" %in% names(family_dt))) {
    stop("Family table must contain a 'sample' column")
  }

  possible <- setdiff(names(family_dt), "sample")
  if (length(possible) == 0) {
    stop("Family table has no label column besides sample")
  }

  family_col <- possible[1]
  fam_keep <- family_dt[, c("sample", family_col), with = FALSE]
  setnames(fam_keep, family_col, "family_label")

  pcs <- merge(pcs, fam_keep, by = "sample", all.x = TRUE)
  pcs[is.na(family_label), family_label := "unknown"]
} else {
  pcs[, family_label := "all"]
}

# -----------------------------
# Read tP windows
# -----------------------------
tp <- fread(tp_file)
req_tp <- c("sample", "chrom", "WinStart", "WinStop", "WinCenter", "tP")
miss_tp <- setdiff(req_tp, names(tp))
if (length(miss_tp) > 0) {
  stop("tP file missing columns: ", paste(miss_tp, collapse = ", "))
}

tp <- tp[chrom == cand_chrom]
tp <- tp[sample %in% pcs$sample]

if (nrow(tp) == 0) {
  stop("No tP rows left after filtering to chromosome and STEP12 samples")
}

tp <- merge(
  tp,
  pcs[, .(sample, group_label)],
  by = "sample",
  all.x = FALSE
)

# -----------------------------
# Normalization
# -----------------------------
chr_mean <- mean(tp$tP, na.rm = TRUE)
chr_sd   <- sd(tp$tP, na.rm = TRUE)
tp[, tP_chr_z := if (is.finite(chr_sd) && chr_sd > 0) (tP - chr_mean) / chr_sd else NA_real_]

tp_zoom <- tp[WinCenter >= zoom_start & WinCenter <= zoom_end]
if (nrow(tp_zoom) == 0) {
  stop("No tP windows in zoom interval")
}

# Flanks = zoom interval excluding candidate overlap
flank_rows <- tp_zoom[!(WinStart <= inv_end & WinStop >= inv_start)]
flank_mean <- mean(flank_rows$tP, na.rm = TRUE)

if (!is.finite(flank_mean) || flank_mean <= 0) {
  flank_mean <- mean(tp_zoom$tP, na.rm = TRUE)
}

tp_zoom[, tP_flank_norm := tP / flank_mean]

# -----------------------------
# Summaries for plotting
# -----------------------------
sum_raw <- tp_zoom[
  is.finite(tP),
  .(
    tP_mean   = mean(tP, na.rm = TRUE),
    tP_median = median(tP, na.rm = TRUE),
    n         = .N
  ),
  by = .(group_label, WinCenter)
]

sum_z <- tp_zoom[
  is.finite(tP_chr_z),
  .(
    tP_chr_z_mean   = mean(tP_chr_z, na.rm = TRUE),
    tP_chr_z_median = median(tP_chr_z, na.rm = TRUE),
    n               = .N
  ),
  by = .(group_label, WinCenter)
]

sum_norm <- tp_zoom[
  is.finite(tP_flank_norm),
  .(
    tP_flank_norm_mean   = mean(tP_flank_norm, na.rm = TRUE),
    tP_flank_norm_median = median(tP_flank_norm, na.rm = TRUE),
    n                    = .N
  ),
  by = .(group_label, WinCenter)
]

zoom_table <- merge(sum_raw, sum_z, by = c("group_label", "WinCenter", "n"), all = TRUE)
zoom_table <- merge(zoom_table, sum_norm, by = c("group_label", "WinCenter", "n"), all = TRUE)

zoom_table[, candidate_id := candidate_id_val]
zoom_table[, chrom := cand_chrom]
zoom_table[, inv_start := inv_start]
zoom_table[, inv_end := inv_end]
zoom_table[, zoom_start := zoom_start]
zoom_table[, zoom_end := zoom_end]
zoom_table[, flank_mean := flank_mean]

fwrite(
  zoom_table,
  paste0(outprefix, ".candidate_", candidate_id_val, ".theta_zoom.tsv.gz"),
  sep = "\t"
)

# -----------------------------
# Plot 1: theta zoom raw
# -----------------------------
p1 <- ggplot(zoom_table, aes(x = WinCenter / 1e6, y = tP_mean, color = group_label)) +
  annotate(
    "rect",
    xmin = inv_start / 1e6, xmax = inv_end / 1e6,
    ymin = -Inf, ymax = Inf,
    alpha = 0.12, fill = "grey60"
  ) +
  geom_line(linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Candidate ", candidate_id_val, " - theta zoom raw"),
    subtitle = paste0(cand_chrom, ": ", format(zoom_start, big.mark = ","), "-", format(zoom_end, big.mark = ",")),
    x = "Chromosome position (Mb)",
    y = "Mean tP",
    color = "Group"
  )

# -----------------------------
# Plot 2: theta zoom flank-normalized
# -----------------------------
p2 <- ggplot(zoom_table, aes(x = WinCenter / 1e6, y = tP_flank_norm_mean, color = group_label)) +
  annotate(
    "rect",
    xmin = inv_start / 1e6, xmax = inv_end / 1e6,
    ymin = -Inf, ymax = Inf,
    alpha = 0.12, fill = "grey60"
  ) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Candidate ", candidate_id_val, " - theta zoom flank normalized"),
    subtitle = "tP divided by mean flank tP",
    x = "Chromosome position (Mb)",
    y = "Mean tP / flank mean",
    color = "Group"
  )

# -----------------------------
# Plot 3: theta zoom chromosome z-score
# -----------------------------
p3 <- ggplot(zoom_table, aes(x = WinCenter / 1e6, y = tP_chr_z_mean, color = group_label)) +
  annotate(
    "rect",
    xmin = inv_start / 1e6, xmax = inv_end / 1e6,
    ymin = -Inf, ymax = Inf,
    alpha = 0.12, fill = "grey60"
  ) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Candidate ", candidate_id_val, " - theta zoom chromosome z-score"),
    subtitle = "tP standardized within chromosome",
    x = "Chromosome position (Mb)",
    y = "Mean chromosome z(tP)",
    color = "Group"
  )

# -----------------------------
# Plot 4: family-colored PCA
# -----------------------------
p4 <- ggplot(pcs, aes(x = PC1, y = PC2, color = family_label, shape = group_label)) +
  geom_point(size = 2.2, alpha = 0.85) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Candidate ", candidate_id_val, " - regional PCA colored by family/ancestry"),
    subtitle = cand_chrom,
    x = "PC1",
    y = "PC2",
    color = "Family",
    shape = "Group"
  )

# -----------------------------
# Plot 5: heterozygosity by group
# -----------------------------
p5 <- ggplot(pcs, aes(x = group_label, y = regional_het, fill = group_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Candidate ", candidate_id_val, " - regional heterozygosity by inferred group"),
    x = "Group",
    y = "Regional dosage heterozygosity"
  ) +
  theme(legend.position = "none")

# -----------------------------
# Save plots: PDF + PNG
# -----------------------------
ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".theta_zoom_raw.pdf"), p1, width = 8, height = 4.8)
ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".theta_zoom_raw.png"), p1, width = 8, height = 4.8, dpi = 400)

ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".theta_zoom_flanknorm.pdf"), p2, width = 8, height = 4.8)
ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".theta_zoom_flanknorm.png"), p2, width = 8, height = 4.8, dpi = 400)

ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".theta_zoom_chrZ.pdf"), p3, width = 8, height = 4.8)
ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".theta_zoom_chrZ.png"), p3, width = 8, height = 4.8, dpi = 400)

ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".pca_family.pdf"), p4, width = 7, height = 5.5)
ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".pca_family.png"), p4, width = 7, height = 5.5, dpi = 400)

ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".het_by_group.pdf"), p5, width = 6.2, height = 4.8)
ggsave(paste0(outprefix, ".candidate_", candidate_id_val, ".het_by_group.png"), p5, width = 6.2, height = 4.8, dpi = 400)

message("[DONE] Wrote follow-up files for candidate ", candidate_id_val)
