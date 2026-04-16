#!/usr/bin/env Rscript
# =============================================================================
# STEP_D11b_perchr_pruning_table.R — Build per-chromosome kin-pruning table
# =============================================================================
#
# Reads per-chromosome ngsRelate .res files (from 11a), applies greedy
# pruning independently per chromosome, writes one unified table.
#
# NaToRA logic: iteratively remove the sample with the most first-degree
# connections until no first-degree pairs remain. "Most informative" =
# the sample whose removal breaks the most high-kinship pairs.
#
# Output:
#   per_chr_pruned.tsv — columns: chr | sample | status | n_connections | rab_max
#
# INVERSION CONFOUND NOTE:
#   On chromosomes with large inversions, ngsRelate inflates kinship between
#   carriers. This means inversion carriers get pruned preferentially.
#   This is FINE for the peeling diagnostic (it's what we want to test).
#   It's NOT fine for population genetics estimates (use genome-wide pruning
#   or mask inversions first for that).
#
# Usage:
#   Rscript 11b_build_perchr_pruning_table.R \
#     --ngsrelate_dir ngsrelate_perchr/ \
#     --sample_list samples.txt \
#     --outfile per_chr_pruned.tsv \
#     [--kin_threshold 0.05] \
#     [--max_remove_frac 0.4]
#
# Or: also compute the FAST PC1-based pruning as L1b backup:
#   Rscript 11b_build_perchr_pruning_table.R \
#     --ngsrelate_dir ngsrelate_perchr/ \
#     --precomp_dir precomp/ \
#     --sample_list samples.txt \
#     --outfile per_chr_pruned.tsv \
#     --also_pc1
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
ngsrel_dir <- NULL; sample_list_file <- NULL; precomp_dir <- NULL
outfile <- "per_chr_pruned.tsv"
kin_threshold <- 0.05; max_remove_frac <- 0.4
also_pc1 <- FALSE; pc1_cor_threshold <- 0.7

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--ngsrelate_dir" && i < length(args))   { ngsrel_dir <- args[i+1]; i <- i+2 }
  else if (a == "--sample_list" && i < length(args)) { sample_list_file <- args[i+1]; i <- i+2 }
  else if (a == "--precomp_dir" && i < length(args)) { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--outfile" && i < length(args))     { outfile <- args[i+1]; i <- i+2 }
  else if (a == "--kin_threshold" && i < length(args)) { kin_threshold <- as.numeric(args[i+1]); i <- i+2 }
  else if (a == "--max_remove_frac" && i < length(args)) { max_remove_frac <- as.numeric(args[i+1]); i <- i+2 }
  else if (a == "--also_pc1")                        { also_pc1 <- TRUE; i <- i+1 }
  else if (a == "--pc1_cor_threshold" && i < length(args)) { pc1_cor_threshold <- as.numeric(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(ngsrel_dir)) stop("--ngsrelate_dir required")

# ---- Load sample list ----
sample_names <- NULL
if (!is.null(sample_list_file) && file.exists(sample_list_file)) {
  raw_names <- readLines(sample_list_file)
  sample_names <- sapply(raw_names, function(x) {
    x <- trimws(basename(x))
    x <- sub("\\.sorted\\.markdup\\.bam$", "", x)
    sub("\\.bam$", "", x)
  }, USE.NAMES = FALSE)
  sample_names <- sample_names[nchar(sample_names) > 0]
  message("[PRUNE] Sample list: ", length(sample_names), " samples")
}

# ---- Find ngsRelate result files ----
res_files <- list.files(ngsrel_dir, pattern = "\\.ngsrelate\\.res$", full.names = TRUE)
if (length(res_files) == 0) stop("No .ngsrelate.res files in ", ngsrel_dir)
message("[PRUNE] Found ", length(res_files), " chromosome results")

# ---- NaToRA-style greedy pruning per chromosome ----

natora_prune <- function(pairs_dt, all_samples, kin_thresh, max_remove_frac) {
  # pairs_dt: data.table with columns sample1, sample2, rab
  # Returns: character vector of samples to REMOVE

  high_kin <- pairs_dt[rab >= kin_thresh]
  if (nrow(high_kin) == 0) return(character(0))

  to_remove <- character(0)
  max_remove <- floor(length(all_samples) * max_remove_frac)
  remaining <- copy(high_kin)

  while (nrow(remaining) > 0 && length(to_remove) < max_remove) {
    # Count connections per sample
    counts <- table(c(remaining$sample1, remaining$sample2))
    # Remove the most-connected sample
    worst <- names(which.max(counts))
    to_remove <- c(to_remove, worst)
    remaining <- remaining[sample1 != worst & sample2 != worst]
  }

  intersect(to_remove, all_samples)
}

all_rows <- list()

for (res_file in sort(res_files)) {
  chr <- sub("\\.ngsrelate\\.res$", "", basename(res_file))
  message("\n[PRUNE] ", chr)

  res <- tryCatch(fread(res_file), error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) {
    message("  Empty result — keeping all samples")
    if (!is.null(sample_names)) {
      all_rows[[length(all_rows) + 1]] <- data.table(
        chr = chr, sample = sample_names, status = "keep",
        n_connections = 0L, rab_max = 0, prune_method = "ngsRelate"
      )
    }
    next
  }

  # Parse ngsRelate output
  # Standard columns: a, b, rab (0-based indices), plus KING, fratij, R0, R1, etc.
  if (!"rab" %in% names(res)) {
    # Try column 3 as rab
    if (ncol(res) >= 3) {
      names(res)[1:3] <- c("a", "b", "rab")
    } else {
      message("  Cannot parse — skipping"); next
    }
  }

  # Convert 0-based indices to sample names
  if (!is.null(sample_names) && is.numeric(res$a)) {
    res[, sample1 := sample_names[a + 1L]]
    res[, sample2 := sample_names[b + 1L]]
  } else if ("ida" %in% names(res)) {
    res[, sample1 := ida]
    res[, sample2 := idb]
  } else {
    res[, sample1 := as.character(a)]
    res[, sample2 := as.character(b)]
  }

  # Keep only valid pairs
  res <- res[is.finite(rab) & !is.na(sample1) & !is.na(sample2)]

  n_high <- sum(res$rab >= kin_threshold)
  message("  Pairs: ", nrow(res), " total, ", n_high, " above threshold ", kin_threshold)

  # Get all samples for this chromosome
  chr_samples <- if (!is.null(sample_names)) sample_names else {
    unique(c(res$sample1, res$sample2))
  }

  # NaToRA-style pruning
  to_remove <- natora_prune(res[, .(sample1, sample2, rab)],
                             chr_samples, kin_threshold, max_remove_frac)

  # Per-sample stats
  conn_counts <- integer(length(chr_samples))
  rab_maxes <- numeric(length(chr_samples))
  names(conn_counts) <- chr_samples
  names(rab_maxes) <- chr_samples

  high_pairs <- res[rab >= kin_threshold]
  if (nrow(high_pairs) > 0) {
    for (s in chr_samples) {
      s_pairs <- high_pairs[sample1 == s | sample2 == s]
      conn_counts[s] <- nrow(s_pairs)
      rab_maxes[s] <- if (nrow(s_pairs) > 0) max(s_pairs$rab) else 0
    }
  }

  message("  Remove: ", length(to_remove), " / Keep: ",
          length(chr_samples) - length(to_remove))

  all_rows[[length(all_rows) + 1]] <- data.table(
    chr = chr,
    sample = chr_samples,
    status = ifelse(chr_samples %in% to_remove, "remove", "keep"),
    n_connections = conn_counts[chr_samples],
    rab_max = round(rab_maxes[chr_samples], 4),
    prune_method = "ngsRelate"
  )
}

# ---- Optional: also compute PC1-based pruning ----
if (also_pc1 && !is.null(precomp_dir)) {
  message("\n[PRUNE] Also computing PC1-based chromosome-local pruning...")

  rds_files <- list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE)
  for (rds_file in sort(rds_files)) {
    chr <- sub("\\.precomp\\.rds$", "", basename(rds_file))

    pc <- readRDS(rds_file)
    dt <- pc$dt
    pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
    snames <- sub("^PC_1_", "", pc1_cols)

    if (length(snames) < 20) next

    # Subsample windows
    n_win <- nrow(dt)
    win_idx <- if (n_win > 2000) sort(sample(n_win, 2000)) else seq_len(n_win)

    sub <- as.matrix(dt[win_idx, ..pc1_cols])
    sub[is.na(sub)] <- 0
    samp_cor <- cor(sub, use = "pairwise.complete.obs")
    samp_cor[!is.finite(samp_cor)] <- 0
    diag(samp_cor) <- 0

    # Greedy removal
    adj <- abs(samp_cor) >= pc1_cor_threshold
    pc1_remove <- character(0)
    max_rm <- floor(length(snames) * max_remove_frac)

    while (any(adj) && length(pc1_remove) < max_rm) {
      conn <- colSums(adj)
      if (max(conn) == 0) break
      worst_idx <- which.max(conn)
      pc1_remove <- c(pc1_remove, snames[worst_idx])
      adj[worst_idx, ] <- FALSE
      adj[, worst_idx] <- FALSE
    }

    message("  [PC1] ", chr, ": remove ", length(pc1_remove))

    all_rows[[length(all_rows) + 1]] <- data.table(
      chr = chr,
      sample = snames,
      status = ifelse(snames %in% pc1_remove, "remove", "keep"),
      n_connections = 0L,
      rab_max = 0,
      prune_method = "PC1_correlation"
    )
  }
}

# ---- Write ----
result <- rbindlist(all_rows, fill = TRUE)
fwrite(result, outfile, sep = "\t")

# Summary
summary_dt <- result[prune_method == "ngsRelate",
                      .(n_keep = sum(status == "keep"),
                        n_remove = sum(status == "remove"),
                        pct_keep = round(100 * sum(status == "keep") / .N, 1)),
                      by = chr]

message("\n=== Per-chromosome ngsRelate pruning summary ===")
for (r in seq_len(nrow(summary_dt))) {
  message("  ", summary_dt$chr[r], ": keep ", summary_dt$n_keep[r],
          " / remove ", summary_dt$n_remove[r],
          " (", summary_dt$pct_keep[r], "% kept)")
}

if (also_pc1) {
  pc1_summary <- result[prune_method == "PC1_correlation",
                         .(n_remove = sum(status == "remove")), by = chr]
  message("\n=== PC1 pruning comparison ===")
  for (r in seq_len(nrow(pc1_summary))) {
    message("  ", pc1_summary$chr[r], ": PC1 removes ", pc1_summary$n_remove[r])
  }
}

message("\n[DONE] Wrote: ", outfile, " (", nrow(result), " rows, ",
        length(unique(result$chr)), " chromosomes)")
