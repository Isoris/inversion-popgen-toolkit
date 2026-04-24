#!/usr/bin/env Rscript
# =============================================================================
# q07_build_groups.R
#
# Builds two sets of sample group files for Engine F (region_popstats):
#
#   Strategy A — ancestry groups (genome-wide / chromosome-wide):
#     For each sample, take its dominant maxQ_label from Engine B's
#     local_Q_samples (mode across all windows on this chromosome, or
#     majority-vote). Bucket samples by dominant label.
#
#   Strategy B — inversion-genotype groups (region-local):
#     For samples' PC1 loadings averaged across windows inside the
#     shelf region, k-means(k=3). Label by PC1 centroid sign: lowest
#     = Hom1, middle = Het, highest = Hom2.
#     Also emit assignment TSV for downstream use + Q04 panel.
#
# Writes:
#   <groups_dir>/ancestry/<label>.txt      one sample ID per line
#   <groups_dir>/invgt/<Hom1|Het|Hom2>.txt
#   <invgt_out>                            assignment table
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
PRECOMP        <- get_arg("--precomp")
CHROM          <- get_arg("--chrom")
SAMPLE_LIST    <- get_arg("--sample_list")
LOCALQ_SAMPLES <- get_arg("--localq_samples", NULL)
GROUPS_DIR     <- get_arg("--groups_dir")
INVGT_OUT      <- get_arg("--invgt_out")
MIN_GROUP_N    <- as.integer(get_arg("--min_group_n", "10"))
INV_K          <- as.integer(get_arg("--inv_k", "3"))
SHELF_A        <- as.numeric(get_arg("--shelf_start_mb", NA))
SHELF_B        <- as.numeric(get_arg("--shelf_end_mb",   NA))

stopifnot(!is.null(PRECOMP), file.exists(PRECOMP))
stopifnot(!is.null(SAMPLE_LIST), file.exists(SAMPLE_LIST))
stopifnot(!is.null(GROUPS_DIR))

dir.create(file.path(GROUPS_DIR, "ancestry"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(GROUPS_DIR, "invgt"),    recursive = TRUE, showWarnings = FALSE)

# ---- Load master sample list in BAM order -----------------------------------
# Engine F needs sample lists that are SUBSETS of this master list (any order)
sample_master <- fread(SAMPLE_LIST, header = FALSE)[[1]]
message("[q07] Master samples: ", length(sample_master))

# ---- Load precomp for PC loadings -------------------------------------------
pc <- readRDS(PRECOMP)
dt <- as.data.table(pc$dt)
dt[, mb := (start_bp + end_bp) / 2 / 1e6]
pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
# Normalize sample IDs used by precomp (Ind0..Ind225 OR CGA...)
precomp_samples <- sub("^PC_1_", "", pc1_cols)

# Mapping precomp samples <-> master samples
# Case A: precomp uses Ind0..IndN → map by positional index to sample_master
# Case B: precomp uses CGAxxx → identity
if (all(grepl("^Ind[0-9]+$", precomp_samples))) {
  ind_idx <- as.integer(sub("^Ind", "", precomp_samples))
  # Ind0 = sample_master[1], Ind1 = sample_master[2], etc.
  precomp_to_master <- sample_master[ind_idx + 1L]
} else {
  precomp_to_master <- precomp_samples
}
# Names of samples in precomp, mapped to master IDs
names(precomp_to_master) <- precomp_samples
message("[q07] Precomp sample format: ",
        if (all(grepl("^Ind[0-9]+$", precomp_samples))) "Ind-indexed" else "CGA-direct")

# =============================================================================
# STRATEGY A: Ancestry groups
# =============================================================================
anc_files_written <- 0
if (!is.null(LOCALQ_SAMPLES) && file.exists(LOCALQ_SAMPLES)) {
  message("[q07] Loading local_Q_samples: ", LOCALQ_SAMPLES)
  ls <- fread(LOCALQ_SAMPLES)
  setnames(ls, tolower(names(ls)))
  if ("chrom" %in% names(ls)) ls <- ls[chrom == CHROM]

  # Expect column 'maxq_label' (K1..K8) and 'sample'
  if (all(c("sample", "maxq_label") %in% names(ls))) {
    # Take the most common maxQ across windows per sample (majority vote)
    mode_label <- function(x) {
      tt <- table(x, useNA = "no")
      if (length(tt) == 0) return(NA_character_)
      names(tt)[which.max(tt)]
    }
    dom <- ls[, .(dom_anc = mode_label(maxq_label)), by = sample]
    setkey(dom, sample)
    # Bucket: unique levels
    for (lbl in unique(na.omit(dom$dom_anc))) {
      members <- dom[dom_anc == lbl, sample]
      # Members are already master IDs (Engine B uses master sample names)
      # Filter to master list for sanity
      members <- intersect(members, sample_master)
      if (length(members) >= MIN_GROUP_N) {
        out_f <- file.path(GROUPS_DIR, "ancestry", paste0(lbl, ".txt"))
        writeLines(members, out_f)
        message("[q07] ancestry ", lbl, ": n=", length(members), " -> ", basename(out_f))
        anc_files_written <- anc_files_written + 1
      } else {
        message("[q07] ancestry ", lbl, ": n=", length(members),
                " < MIN_GROUP_N=", MIN_GROUP_N, ", skipping")
      }
    }
  } else {
    message("[q07] local_Q_samples lacks required columns (sample, maxq_label); ",
            "no ancestry groups built")
  }
} else {
  message("[q07] No local_Q_samples provided; skipping Strategy A")
}
message("[q07] Strategy A: ", anc_files_written, " group file(s) written")

# =============================================================================
# STRATEGY B: Inversion-genotype groups via PCA k-means(k=3) at shelf region
# =============================================================================
invgt_files_written <- 0
assign_dt <- NULL

if (is.finite(SHELF_A) && is.finite(SHELF_B)) {
  message(sprintf("[q07] Strategy B: k=%d clustering on PC1 at %.2f-%.2f Mb",
                  INV_K, SHELF_A, SHELF_B))
  shelf_mask <- dt$mb >= SHELF_A & dt$mb <= SHELF_B
  n_shelf <- sum(shelf_mask)
  if (n_shelf < 5) {
    message("[q07] Shelf region has only ", n_shelf, " windows, need >=5. Skipping Strategy B")
  } else {
    # Mean PC1 loading per sample across shelf windows
    pc1_mat <- as.matrix(dt[shelf_mask, ..pc1_cols])
    pc1_mean <- colMeans(pc1_mat, na.rm = TRUE)
    # k-means
    set.seed(1L)
    km <- kmeans(pc1_mean, centers = INV_K, nstart = 25L, iter.max = 50L)
    centroids <- km$centers[, 1]
    # Rank centroids; label by rank: 1 = Hom1, middle = Het, last = Hom2
    ord <- order(centroids)
    centroid_rank <- match(seq_along(centroids), ord)
    if (INV_K == 3) {
      lbl_by_rank <- c("Hom1", "Het", "Hom2")
    } else {
      lbl_by_rank <- paste0("G", seq_len(INV_K))
    }
    labels_per_cluster <- lbl_by_rank[centroid_rank]
    sample_labels <- labels_per_cluster[km$cluster]

    # Map precomp IDs -> master IDs
    assign_dt <- data.table(
      sample_precomp = precomp_samples,
      sample_master  = precomp_to_master,
      pc1_mean       = round(pc1_mean, 5),
      invgt          = sample_labels
    )
    # Write assignment TSV
    fwrite(assign_dt, INVGT_OUT, sep = "\t", quote = FALSE)
    message("[q07] invgt assignments -> ", INVGT_OUT)

    # Cluster diagnostics
    msg_parts <- paste(lbl_by_rank, "n=", table(factor(sample_labels, levels = lbl_by_rank))[lbl_by_rank],
                       "centroid=", sprintf("%+.4f", centroids[ord]))
    message("[q07] clusters: ", paste(msg_parts, collapse = "; "))

    # Write group files (master IDs only, matched to sample_master)
    for (lbl in unique(sample_labels)) {
      members <- assign_dt[invgt == lbl, sample_master]
      members <- intersect(members, sample_master)
      if (length(members) >= 3) {
        out_f <- file.path(GROUPS_DIR, "invgt", paste0(lbl, ".txt"))
        writeLines(members, out_f)
        message("[q07] invgt ", lbl, ": n=", length(members), " -> ", basename(out_f))
        invgt_files_written <- invgt_files_written + 1
      }
    }
  }
} else {
  message("[q07] No shelf region specified; skipping Strategy B")
  message("      Pass --shelf_start_mb X --shelf_end_mb Y to enable inv-genotype grouping")
}

message("[q07] Strategy B: ", invgt_files_written, " group file(s) written")
message("[q07] DONE.")
