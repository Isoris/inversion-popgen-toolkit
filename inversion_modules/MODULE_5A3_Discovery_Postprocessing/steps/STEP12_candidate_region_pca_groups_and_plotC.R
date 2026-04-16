#!/usr/bin/env Rscript

# =============================================================================
# STEP12_candidate_region_pca_groups_and_plotC.R  (v2 — adaptive-k)
#
# Regional PCA per candidate with BIC-optimal k selection.
#
# Changes from v1:
#   - k is selected from range [2,5] by BIC (using mclust if available, else
#     within-cluster SS elbow heuristic)
#   - Gradient detection: if variance explained by k=2 grouping on PC1 is low,
#     the candidate is flagged as gradient rather than forced into discrete groups
#   - Output includes best_k, BIC values for all k, and model_type
#   - Backward compatible: default still outputs the same columns
#   - The n_groups argument is now interpreted as max_k (default 5)
#
# Usage:
#   Rscript STEP12_candidate_region_pca_groups_and_plotC.R \
#     <candidate_regions.tsv.gz> <candidate_id> <sites.tsv.gz> \
#     <dosage.tsv.gz> <sample_tP_windows.tsv.gz> <sample_order_file> <outprefix> \
#     [max_k=5]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

use_mclust <- suppressWarnings(require(mclust, quietly = TRUE))
if (use_mclust) {
  message("[INFO] Using mclust for BIC-based k selection")
} else {
  message("[INFO] mclust unavailable; using WSS elbow heuristic for k selection")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop(paste(
    "Usage: Rscript STEP12_candidate_region_pca_groups_and_plotC.R",
    "<candidate_regions.tsv.gz> <candidate_id> <sites.tsv.gz>",
    "<dosage.tsv.gz> <sample_tP_windows.tsv.gz> <sample_order_file> <outprefix>",
    "[max_k=5]"
  ))
}

candidate_file     <- args[1]
candidate_id_value <- args[2]
sites_file         <- args[3]
dosage_file        <- args[4]
tP_windows_file    <- args[5]
sample_order_file  <- args[6]
outprefix          <- args[7]
max_k              <- if (length(args) >= 8) as.integer(args[8]) else 5L
min_k              <- 2L

# ── Read candidate region ──────────────────────────────────────────────────
message("[INFO] Reading candidate regions")
cand <- fread(candidate_file)

if (!("candidate_id" %in% names(cand))) {
  stop("candidate_regions file must contain 'candidate_id' column")
}

cand_row <- cand[as.character(candidate_id) == as.character(candidate_id_value)]
if (nrow(cand_row) != 1) {
  stop("Expected exactly one candidate row for candidate_id=", candidate_id_value,
       "; found ", nrow(cand_row))
}

region_chrom <- cand_row$chrom[1]
region_start <- as.numeric(cand_row$start_bp[1])
region_end   <- as.numeric(cand_row$end_bp[1])

message("[INFO] Candidate region: ", region_chrom, ":",
        format(region_start, big.mark = ","), "-",
        format(region_end, big.mark = ","))

# ── Read sites + dosage ───────────────────────────────────────────────────
message("[INFO] Reading sites")
sites <- fread(sites_file)
req_sites <- c("marker", "chrom", "pos", "allele1", "allele2")
miss_sites <- setdiff(req_sites, names(sites))
if (length(miss_sites) > 0) {
  stop("Sites file missing required columns: ", paste(miss_sites, collapse = ", "))
}

message("[INFO] Reading dosage")
dos <- fread(dosage_file)
if (!("marker" %in% names(dos))) stop("Dosage file must contain 'marker' column")

if (!identical(sites$marker, dos$marker)) {
  setkeyv(sites, "marker")
  setkeyv(dos, "marker")
  dos <- dos[sites$marker]
  if (!identical(sites$marker, dos$marker)) stop("Sites and dosage markers do not match")
}

# ── Map Ind0..IndN to real sample IDs if needed ───────────────────────────
sample_names <- setdiff(names(dos), "marker")

if (all(grepl("^Ind[0-9]+$", sample_names))) {
  message("[INFO] Dosage columns are Ind-style; mapping to real sample IDs")
  samp_raw <- fread(sample_order_file, header = FALSE, sep = "\t", fill = TRUE)
  sample_order <- samp_raw[[1]]
  sample_order <- sample_order[nchar(sample_order) > 0]
  if (length(sample_order) != length(sample_names)) {
    stop("sample_order_file length mismatch")
  }
  setnames(dos, old = sample_names, new = sample_order)
  sample_names <- sample_order
}

message("[INFO] Samples: ", length(sample_names))

# ── Filter to candidate region ─────────────────────────────────────────────
keep <- which(sites$chrom == region_chrom &
              sites$pos >= region_start &
              sites$pos <= region_end)
if (length(keep) < 10) {
  stop("Too few SNPs in candidate region: ", length(keep))
}
message("[INFO] SNPs in candidate region: ", length(keep))

sites_reg <- sites[keep]
dos_reg   <- dos[keep]
X <- as.matrix(dos_reg[, ..sample_names])
storage.mode(X) <- "double"

# ── Regional PCA on individuals ────────────────────────────────────────────
pca <- prcomp(t(X), center = TRUE, scale. = FALSE)

pcs <- data.table(
  sample = sample_names,
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2]
)

var_expl <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_pct <- round(var_expl[1], 1)
pc2_pct <- round(var_expl[2], 1)

# ── Adaptive k selection ─────────────────────────────────────────────────
pc1_data <- matrix(pcs$PC1, ncol = 1)

if (use_mclust) {
  # BIC-based selection using Gaussian mixture models
  mc <- tryCatch(
    Mclust(pc1_data, G = min_k:max_k, verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(mc)) {
    best_k <- mc$G
    bic_values <- mc$BIC
    # Extract BIC for each k (mclust returns matrix: models × k)
    # Take best model per k
    bic_per_k <- apply(bic_values, 2, max, na.rm = TRUE)
    bic_per_k[!is.finite(bic_per_k)] <- NA_real_
    model_type <- if (best_k == 1) "gradient" else paste0("discrete_k", best_k)
    cluster_assignments <- mc$classification
  } else {
    message("[WARN] mclust failed; falling back to k-means WSS")
    use_mclust <- FALSE
  }
}

if (!use_mclust) {
  # WSS elbow heuristic
  wss <- numeric(max_k)
  km_results <- list()
  for (k in seq_len(max_k)) {
    set.seed(1)
    km_k <- kmeans(pc1_data, centers = k, nstart = 50)
    wss[k] <- km_k$tot.withinss
    km_results[[k]] <- km_k
  }

  # Elbow: largest drop in WSS (second derivative)
  if (max_k >= 3) {
    dwss <- diff(wss)  # first differences (negative)
    ddwss <- diff(dwss) # second differences
    # Best k = where the drop in WSS reduction is largest (biggest elbow)
    # Candidates: k = 2, 3, ..., max_k
    # ddwss[i] corresponds to k = i+2
    elbow_idx <- which.max(ddwss) + 1  # +1 because ddwss[1] = k=3
    best_k <- max(min_k, min(elbow_idx + 1, max_k))
  } else {
    best_k <- min_k
  }

  bic_per_k <- -wss[min_k:max_k]  # negative WSS as pseudo-BIC (higher = better)
  names(bic_per_k) <- paste0("k", min_k:max_k)
  cluster_assignments <- km_results[[best_k]]$cluster
  model_type <- paste0("discrete_k", best_k)
}

# Gradient detection: if PC1 variance explained is low or best_k grouping
# has poor separation, flag as gradient
if (best_k >= 2) {
  grp_means <- tapply(pcs$PC1, cluster_assignments, mean)
  grp_sds <- tapply(pcs$PC1, cluster_assignments, sd)
  overall_sd <- sd(pcs$PC1)
  # Between-group variance fraction
  between_var <- var(grp_means[cluster_assignments]) / var(pcs$PC1)
  if (!is.finite(between_var)) between_var <- 0
  if (between_var < 0.3) {
    model_type <- "gradient"
    message("[INFO] Low between-group variance (", round(between_var, 3),
            "); flagging as gradient")
  }
}

message("[INFO] Best k: ", best_k, " (model_type: ", model_type, ")")

# ── Apply grouping ─────────────────────────────────────────────────────────
set.seed(1)
km <- kmeans(pc1_data, centers = best_k, nstart = 50)
pcs[, raw_group := km$cluster]

grp_order <- pcs[, .(pc1_mean = mean(PC1)), by = raw_group][order(pc1_mean)]
grp_order[, ordered_group := seq_len(.N)]
pcs <- merge(pcs, grp_order[, .(raw_group, ordered_group)], by = "raw_group", all.x = TRUE)
pcs <- pcs[match(sample_names, sample)]

# Dynamic group labels
if (best_k == 3) {
  group_labels <- c("Homo_1", "Het", "Homo_2")
} else if (best_k == 2) {
  group_labels <- c("Group_1", "Group_2")
} else {
  group_labels <- paste0("G", seq_len(best_k))
}
pcs[, group_label := group_labels[ordered_group]]
pcs[, group_label := factor(group_label, levels = group_labels)]

# Add model selection metadata
pcs[, best_k := best_k]
pcs[, model_type := model_type]

message("[INFO] Group sizes: ",
        paste(pcs[, .N, by = group_label][order(group_label), paste0(group_label, "=", N)], collapse = ", "))

# ── Regional heterozygosity from dosage ────────────────────────────────────
H <- X * (2 - X)
regional_het <- colMeans(H, na.rm = TRUE)
pcs[, regional_het := regional_het[match(sample, sample_names)]]

pc1_vec <- pcs$PC1[match(sample_names, pcs$sample)]

flip_flag <- rep(FALSE, nrow(X))
for (i in seq_len(nrow(X))) {
  xi <- X[i, ]
  ok <- is.finite(xi) & is.finite(pc1_vec)
  if (sum(ok) >= 3) {
    rr <- suppressWarnings(cor(xi[ok], pc1_vec[ok]))
    if (is.finite(rr) && rr < 0) flip_flag[i] <- TRUE
  }
}

X_oriented <- X
X_oriented[flip_flag, ] <- 2 - X_oriented[flip_flag, ]
regional_hap_score <- colMeans(X_oriented, na.rm = TRUE)
pcs[, regional_hap_score := regional_hap_score[match(sample, sample_names)]]

# ── Representative individuals ─────────────────────────────────────────────
centers <- pcs[, .(PC1c = mean(PC1), PC2c = mean(PC2)), by = group_label]
pcs2 <- merge(pcs, centers, by = "group_label", all.x = TRUE)
pcs2[, dist_center := sqrt((PC1 - PC1c)^2 + (PC2 - PC2c)^2)]
reps <- pcs2[order(dist_center), .SD[1], by = group_label]

reps[, candidate_id       := as.integer(candidate_id_value)]
reps[, candidate_chrom    := region_chrom]
reps[, candidate_start_bp := region_start]
reps[, candidate_end_bp   := region_end]

# ── Read per-sample tP windows ────────────────────────────────────────────
message("[INFO] Reading per-sample tP windows")
tP_all <- fread(tP_windows_file)

req_tP <- c("sample", "chrom", "WinStart", "WinStop", "WinCenter", "tP")
miss_tP <- setdiff(req_tP, names(tP_all))
if (length(miss_tP) > 0) {
  stop("tP windows file missing required columns: ", paste(miss_tP, collapse = ", "))
}

tP_chr <- tP_all[sample %in% sample_names & chrom == region_chrom]

if (nrow(tP_chr) == 0) {
  stop("No tP data for chromosome ", region_chrom, " and the analysis samples")
}

message("[INFO] tP windows for this chromosome: ", nrow(tP_chr),
        " (", uniqueN(tP_chr$sample), " samples)")

tP_chr <- merge(
  tP_chr,
  pcs[, .(sample, group_label, ordered_group, regional_het, regional_hap_score)],
  by = "sample",
  all.x = FALSE
)

# ── Summarize tP by group and window ──────────────────────────────────────
plotc_dt <- tP_chr[
  is.finite(tP),
  .(
    tP_mean = mean(tP, na.rm = TRUE),
    tP_sd   = sd(tP, na.rm = TRUE),
    n       = .N
  ),
  by = .(group_label, ordered_group, WinCenter)
][order(ordered_group, WinCenter)]

# ── Write outputs ─────────────────────────────────────────────────────────
cid <- candidate_id_value
pcs_out   <- paste0(outprefix, ".candidate_", cid, ".regional_pca_samples.tsv.gz")
sites_out <- paste0(outprefix, ".candidate_", cid, ".regional_sites.tsv.gz")
plotc_out <- paste0(outprefix, ".candidate_", cid, ".plotC_summary.tsv.gz")
rep_out   <- paste0(outprefix, ".candidate_", cid, ".representatives.tsv.gz")

# Model selection summary
model_out <- paste0(outprefix, ".candidate_", cid, ".model_selection.tsv")
model_dt <- data.table(
  candidate_id = as.integer(cid),
  best_k = best_k,
  model_type = model_type,
  pc1_var_pct = pc1_pct,
  pc2_var_pct = pc2_pct,
  n_samples = length(sample_names),
  n_snps = length(keep)
)
# Add BIC columns
for (k in min_k:max_k) {
  kname <- paste0("bic_k", k)
  kidx <- k - min_k + 1
  model_dt[, (kname) := if (kidx <= length(bic_per_k)) bic_per_k[kidx] else NA_real_]
}

fwrite(pcs, pcs_out, sep = "\t")
fwrite(sites_reg, sites_out, sep = "\t")
fwrite(plotc_dt, plotc_out, sep = "\t")
fwrite(reps, rep_out, sep = "\t")
fwrite(model_dt, model_out, sep = "\t")

# ── Plots ─────────────────────────────────────────────────────────────────
p_a <- ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(shape = 1, alpha = 0.6, size = 2) +
  theme_bw(base_size = 12) +
  labs(
    title    = paste0("Candidate ", cid, " — regional PCA (best k=", best_k, ", ", model_type, ")"),
    subtitle = paste0(region_chrom, ": ",
                      format(region_start, big.mark = ","), " – ",
                      format(region_end, big.mark = ","),
                      "  (", length(keep), " SNPs)"),
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  )

p_b <- ggplot(pcs, aes(x = PC1, y = PC2, color = regional_het)) +
  geom_point(aes(shape = group_label), size = 2.5, alpha = 0.8) +
  geom_point(
    data = reps,
    aes(x = PC1, y = PC2),
    inherit.aes = FALSE,
    shape = 1, size = 6, stroke = 1.2, color = "black"
  ) +
  geom_text(
    data = reps,
    aes(x = PC1, y = PC2, label = group_label),
    inherit.aes = FALSE,
    nudge_y = diff(range(pcs$PC2)) * 0.06,
    color = "black", size = 3.5, fontface = "bold"
  ) +
  theme_bw(base_size = 12) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title    = paste0("Candidate ", cid, " — PCA colored by het (k=", best_k, ")"),
    subtitle = "Circled = representative individual nearest each group center",
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)"),
    color = "Regional\nhet (dosage)",
    shape = "Group"
  )

p_c <- ggplot(plotc_dt, aes(x = WinCenter / 1e6, y = tP_mean, color = group_label)) +
  annotate(
    "rect",
    xmin = region_start / 1e6, xmax = region_end / 1e6,
    ymin = -Inf, ymax = Inf,
    alpha = 0.15, fill = "grey60"
  ) +
  geom_line(linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title    = paste0(region_chrom, " — pairwise θ (tP) by inferred PCA group"),
    subtitle = paste0("Shaded region = candidate ", cid,
                      " (k=", best_k, ", ", model_type, ")"),
    x = "Chromosome position (Mb)",
    y = "Mean tP (pairwise θ)",
    color = "Group"
  )

pdf_a <- paste0(outprefix, ".candidate_", cid, ".panelA_regional_PCA.pdf")
pdf_b <- paste0(outprefix, ".candidate_", cid, ".panelB_regional_PCA_het.pdf")
pdf_c <- paste0(outprefix, ".candidate_", cid, ".panelC_chrom_tP_by_group.pdf")

ggsave(pdf_a, p_a, width = 6, height = 5)
ggsave(pdf_b, p_b, width = 7.5, height = 5.5)
ggsave(pdf_c, p_c, width = 8, height = 4.8)

message("[DONE] Wrote:")
message("  ", pcs_out)
message("  ", sites_out)
message("  ", plotc_out)
message("  ", rep_out)
message("  ", model_out)
message("  ", pdf_a)
message("  ", pdf_b)
message("  ", pdf_c)
