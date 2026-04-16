#!/usr/bin/env Rscript

# =============================================================================
# STEP09b_stage1_perchr_pca.R  (v8.3-parallel)
#
# PARALLEL STAGE 1: Compute dense-window PCA for ONE chromosome.
#
# This replaces the inner loop of monolithic STEP09b. Each SLURM array task
# processes one chromosome, producing per-chr outputs with LOCAL indices only.
# Global window_id assignment is deferred to Stage 2 (merge).
#
# FEATURES:
#   - Progress bar with ETA for the PCA window loop
#   - Skip-if-done logic (safe re-runs)
#   - Writes temp outputs to <outdir>/tmp/ (Stage 2 reads from there)
#
# INPUTS:
#   --dosage_dir  — directory with per-chr .dosage.tsv.gz + .sites.tsv.gz
#   --outdir      — output directory (will write to <outdir>/tmp/)
#   --chrom       — chromosome name (e.g. C_gar_LG01)
#   [--winsize 100] [--step 20] [--npc 2]
#
# OUTPUTS (in <outdir>/tmp/):
#   <chrom>.chr_meta.tsv.gz       — per-window metadata (local indices)
#   <chrom>.window_pca_tmp.rds    — PCA data (window_id = NA, patched in Stage 2)
#   <chrom>.pca_table.tsv.gz      — PCA loadings table (window_id = NA)
#   <chrom>.summary.tsv           — single-row summary for this chromosome
#
# Usage:
#   Rscript STEP09b_stage1_perchr_pca.R \
#     --dosage_dir /path/to/04_dosage_by_chr \
#     --outdir /path/to/05b_dense_registry \
#     --chrom C_gar_LG01 \
#     [--winsize 100] [--step 20] [--npc 2]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(RSpectra)
})

# =============================================================================
# PARSE ARGS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

dosage_dir <- NULL
outdir     <- NULL
CHROM      <- NULL
WINSIZE    <- 100L
STEP       <- NA_integer_
NPC        <- 2L

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--dosage_dir" && i < length(args)) {
    dosage_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]; i <- i + 2L
  } else if (a == "--chrom" && i < length(args)) {
    CHROM <- args[i + 1]; i <- i + 2L
  } else if (a == "--winsize" && i < length(args)) {
    WINSIZE <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--step" && i < length(args)) {
    STEP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--npc" && i < length(args)) {
    NPC <- as.integer(args[i + 1]); i <- i + 2L
  } else {
    i <- i + 1L
  }
}

if (is.null(dosage_dir) || is.null(outdir) || is.null(CHROM)) {
  stop("Usage: Rscript STEP09b_stage1_perchr_pca.R --dosage_dir <dir> --outdir <dir> --chrom <chr> ",
       "[--winsize 100] [--step 20] [--npc 2]")
}

if (is.na(STEP)) STEP <- WINSIZE
WINDOW_SET_ID <- paste0("snp", WINSIZE, "_step", STEP)

tmpdir <- file.path(outdir, "tmp")
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SKIP-IF-DONE
# =============================================================================

meta_out <- file.path(tmpdir, paste0(CHROM, ".chr_meta.tsv.gz"))
rds_out  <- file.path(tmpdir, paste0(CHROM, ".window_pca_tmp.rds"))

if (file.exists(meta_out) && file.exists(rds_out)) {
  message("[STEP09b-S1] ", CHROM, ": outputs already exist, skipping")
  message("  ", meta_out)
  message("  ", rds_out)
  quit(status = 0)
}

# =============================================================================
# PCA HELPER (identical to monolithic STEP09b)
# =============================================================================

cov_pca_one_window <- function(x, k) {
  rm <- rowMeans(x, na.rm = TRUE)
  xc <- x - rm
  covmat <- suppressWarnings(cov(xc, use = "pairwise.complete.obs"))
  if (anyNA(covmat) || nrow(covmat) < k) {
    return(list(total = NA_real_, eigval = rep(NA_real_, k),
                eigvec = matrix(NA_real_, nrow = ncol(x), ncol = k)))
  }
  total <- sum(covmat^2)
  if (k >= ncol(covmat)) {
    ee <- eigen(covmat, symmetric = TRUE)
    vals <- ee$values[seq_len(k)]
    vecs <- ee$vectors[, seq_len(k), drop = FALSE]
  } else {
    ee <- tryCatch(
      RSpectra::eigs_sym(covmat, k = k, which = "LM"),
      error = function(e) {
        ee2 <- eigen(covmat, symmetric = TRUE)
        list(values = ee2$values[seq_len(k)], vectors = ee2$vectors[, seq_len(k), drop = FALSE])
      }
    )
    vals <- ee$values; vecs <- ee$vectors
  }
  list(total = total, eigval = vals, eigvec = vecs)
}

# =============================================================================
# PROGRESS BAR HELPER
# =============================================================================

make_progress <- function(total, label = "STEP09b-S1") {
  t0 <- proc.time()[3]
  last_print <- 0
  PRINT_INTERVAL <- 2  # seconds between progress updates

  function(current) {
    now <- proc.time()[3]
    if (now - last_print < PRINT_INTERVAL && current < total) return(invisible())
    last_print <<- now

    elapsed <- now - t0
    pct <- 100 * current / total
    if (current > 0) {
      eta_sec <- elapsed * (total - current) / current
      eta_str <- if (eta_sec >= 3600) {
        sprintf("%dh%02dm", as.integer(eta_sec / 3600), as.integer((eta_sec %% 3600) / 60))
      } else if (eta_sec >= 60) {
        sprintf("%dm%02ds", as.integer(eta_sec / 60), as.integer(eta_sec %% 60))
      } else {
        sprintf("%ds", as.integer(eta_sec))
      }
    } else {
      eta_str <- "?"
    }

    # Build bar
    bar_width <- 30
    filled <- as.integer(bar_width * current / total)
    bar <- paste0("[", strrep("=", filled), strrep(" ", bar_width - filled), "]")

    message(sprintf("\r[%s] %s  %s  %d/%d  (%.1f%%)  ETA %s",
                    label, CHROM, bar, current, total, pct, eta_str))
  }
}

# =============================================================================
# LOAD DATA
# =============================================================================

message("[STEP09b-S1] Dense window PCA — ", CHROM)
message("[STEP09b-S1] Window: ", WINSIZE, " SNPs, Step: ", STEP, " SNPs")
message("[STEP09b-S1] Set ID: ", WINDOW_SET_ID)

sf <- file.path(dosage_dir, paste0(CHROM, ".sites.tsv.gz"))
dos_file <- file.path(dosage_dir, paste0(CHROM, ".dosage.tsv.gz"))

if (!file.exists(sf))       stop("Missing sites file: ", sf)
if (!file.exists(dos_file)) stop("Missing dosage file: ", dos_file)

message("[STEP09b-S1] Reading sites + dosage...")
t_load <- proc.time()[3]
sites <- fread(sf)
dos   <- fread(dos_file)
message("[STEP09b-S1] Loaded in ", round(proc.time()[3] - t_load, 1), "s")

# Align
if (!identical(sites$marker, dos$marker)) {
  setkeyv(sites, "marker"); setkeyv(dos, "marker")
  dos <- dos[sites$marker]
}

sample_names <- setdiff(names(dos), "marker")
X <- as.matrix(dos[, ..sample_names])
storage.mode(X) <- "double"
n_snps <- nrow(X)

starts <- seq(1L, n_snps - WINSIZE + 1L, by = STEP)
n_windows <- length(starts)

if (n_windows < 1) {
  message("[STEP09b-S1] ", CHROM, " — ", n_snps, " SNPs, not enough for one window")
  # Write empty sentinel so Stage 2 knows this chr was processed
  fwrite(data.table(chrom = CHROM, n_snps = n_snps, n_windows = 0L),
         file.path(tmpdir, paste0(CHROM, ".summary.tsv")), sep = "\t")
  quit(status = 0)
}

message("[STEP09b-S1] ", CHROM, ": ", n_snps, " SNPs → ", n_windows, " windows")

# =============================================================================
# PCA LOOP WITH PROGRESS
# =============================================================================

meta_list <- vector("list", n_windows)
pca_rows  <- vector("list", n_windows)

pb <- make_progress(n_windows)

for (wi in seq_len(n_windows)) {
  i1 <- starts[wi]
  i2 <- min(i1 + WINSIZE - 1L, n_snps)
  xw <- X[i1:i2, , drop = FALSE]
  sw <- sites[i1:i2]
  pca <- cov_pca_one_window(xw, NPC)

  # Local metadata — window_id = NA (assigned by Stage 2)
  meta_list[[wi]] <- data.table(
    window_id        = NA_integer_,
    window_set_id    = WINDOW_SET_ID,
    chrom            = CHROM,
    window_index_chr = wi - 1L,
    start_bp         = sw$pos[1],
    end_bp           = sw$pos[nrow(sw)],
    center_bp        = as.integer((sw$pos[1] + sw$pos[nrow(sw)]) / 2),
    start_snp_idx    = i1,
    end_snp_idx      = i2,
    n_snps           = nrow(xw),
    step_snps        = STEP,
    window_snps      = WINSIZE,
    focal_mode       = "focal_dense"
  )

  rowvec <- c(total = pca$total, setNames(pca$eigval, paste0("lam_", seq_len(NPC))))
  for (pc in seq_len(NPC)) {
    rowvec <- c(rowvec, setNames(pca$eigvec[, pc], paste0("PC_", pc, "_", sample_names)))
  }
  pca_rows[[wi]] <- rowvec

  pb(wi)
}
message("")  # newline after progress bar

# =============================================================================
# ASSEMBLE + WRITE
# =============================================================================

chr_meta <- rbindlist(meta_list)
pca_mat  <- as.data.frame(do.call(rbind, pca_rows))
pca_mat$window_id <- NA_integer_
pca_mat <- pca_mat[, c("window_id", setdiff(names(pca_mat), "window_id"))]

# Per-chr meta
fwrite(chr_meta, meta_out, sep = "\t")

# Per-chr PCA (temp — Stage 2 will patch window_id and re-save to final location)
saveRDS(list(
  chrom = CHROM, sample_names = sample_names, npc = NPC,
  winsize = WINSIZE, step = STEP,
  window_meta = chr_meta, pca = pca_mat
), rds_out)

# PCA table
pca_tsv <- file.path(tmpdir, paste0(CHROM, ".pca_table.tsv.gz"))
fwrite(as.data.table(pca_mat), pca_tsv, sep = "\t")

# Summary
fwrite(data.table(
  chrom = CHROM, n_snps = n_snps, n_windows = n_windows,
  step_snps = STEP, window_snps = WINSIZE
), file.path(tmpdir, paste0(CHROM, ".summary.tsv")), sep = "\t")

message("[DONE] STEP09b Stage 1 — ", CHROM, ": ", n_windows, " windows")
message("  ", meta_out)
message("  ", rds_out)
