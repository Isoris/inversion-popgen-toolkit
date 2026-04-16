#!/usr/bin/env Rscript

# =============================================================================
# STEP10i_three_snake_consensus.R  (v1.0)
#
# OVERLAP-BASED CONSENSUS across Snake 1, Snake 2, and Snake 3.
#
# Do NOT combine the 3 snakes by one weighted average.
# Instead, record which snakes support each window, then build meta-regions
# from contiguous supported spans.
#
# Support signatures:
#   PCA+GROUP+GHSL   — 3-way core (strongest)
#   PCA+GROUP         — 2-way
#   PCA+GHSL          — 2-way
#   GROUP+GHSL        — 2-way
#   PCA_ONLY          — 1-way
#   GROUP_ONLY        — 1-way
#   GHSL_ONLY         — 1-way
#   NONE              — no support
#
# Consensus classes for meta-regions:
#   CORE_3WAY, SHOULDER_2WAY, WEAK_1WAY, COMPLEX_SPLIT,
#   PCA_DOMINANT, GROUP_DOMINANT, GHSL_REFINED_CORE, DISCORDANT_COMPLEX
#
# INPUTS:
#   snake1: snake_windows.tsv.gz from STEP10e
#   snake2: snake2_track.tsv.gz from STEP10g
#   snake3: snake3_track.tsv.gz from STEP10h
#   step10: .mds.rds for window coordinates
#
# OUTPUTS:
#   consensus_windows.tsv.gz   — per-window consensus record
#   consensus_regions.tsv.gz   — meta-region summaries
#   consensus_summary.tsv      — per-chromosome summary
#
# Usage:
#   Rscript STEP10i_three_snake_consensus.R \
#     <step10_outprefix> <snake1_dir> <snake2_dir> <snake3_dir> <outdir>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript STEP10i_three_snake_consensus.R ",
       "<step10_outprefix> <snake1_dir> <snake2_dir> <snake3_dir> <outdir>")
}

step10_prefix <- args[1]
snake1_dir    <- args[2]
snake2_dir    <- args[3]
snake3_dir    <- args[4]
outdir        <- args[5]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PARAMETERS
# =============================================================================

GAP_TOLERANCE    <- 2L     # max gap windows within a region
MIN_REGION_WINS  <- 3L     # minimum windows for a consensus region

# =============================================================================
# LOAD ALL THREE SNAKE TRACKS
# =============================================================================

message("[STEP10i] Loading 3 snake tracks...")

# ── STEP10 window coordinates (the common grid) ─────────────────────
mds_rds <- paste0(step10_prefix, ".mds.rds")
if (!file.exists(mds_rds)) stop("Missing: ", mds_rds)
mds_obj <- readRDS(mds_rds)

# Build master window table
all_windows <- rbindlist(lapply(mds_obj$per_chr, function(x) {
  dt <- as.data.table(x$out_dt)
  dt[, .(chrom, global_window_id, start_bp, end_bp)]
}), fill = TRUE)
setkey(all_windows, global_window_id)
message("[STEP10i] Master grid: ", nrow(all_windows), " windows")

# ── Snake 1 (STEP10e) ───────────────────────────────────────────────
snake1_file <- file.path(snake1_dir, "snake_windows.tsv.gz")
if (file.exists(snake1_file)) {
  s1 <- fread(snake1_file)
  # A window is PASS if it was included in ANY snake (core or merged)
  s1_pass_wids <- unique(s1$global_window_id)
  message("[STEP10i] Snake 1: ", length(s1_pass_wids), " supported windows")
} else {
  message("[WARN] Snake 1 file not found: ", snake1_file)
  s1_pass_wids <- integer(0)
}

# ── Snake 2 (STEP10g) ───────────────────────────────────────────────
snake2_file <- file.path(snake2_dir, "snake2_track.tsv.gz")
if (file.exists(snake2_file)) {
  s2 <- fread(snake2_file)
  s2_status <- s2[, .(global_window_id, snake2_status)]
  message("[STEP10i] Snake 2: ", sum(s2$snake2_status == "PASS", na.rm = TRUE), " PASS windows")
} else {
  message("[WARN] Snake 2 file not found: ", snake2_file)
  s2_status <- data.table(global_window_id = integer(), snake2_status = character())
}

# ── Snake 3 (STEP10h) ───────────────────────────────────────────────
snake3_file <- file.path(snake3_dir, "snake3_track.tsv.gz")
if (file.exists(snake3_file)) {
  s3 <- fread(snake3_file)
  s3_status <- s3[, .(global_window_id, snake3_status)]
  message("[STEP10i] Snake 3: ", sum(s3$snake3_status == "PASS", na.rm = TRUE), " PASS windows")
} else {
  message("[WARN] Snake 3 file not found: ", snake3_file)
  s3_status <- data.table(global_window_id = integer(), snake3_status = character())
}

# =============================================================================
# BUILD PER-WINDOW CONSENSUS
# =============================================================================

message("[STEP10i] Building per-window consensus...")

# Merge all tracks onto master grid
cons <- copy(all_windows)

# Snake 1: PASS if window is in snake_windows
cons[, pca_status := fifelse(global_window_id %in% s1_pass_wids, "PASS", "FAIL")]

# Snake 2: merge
cons <- merge(cons, s2_status, by = "global_window_id", all.x = TRUE)
cons[is.na(snake2_status), snake2_status := "FAIL"]
setnames(cons, "snake2_status", "group_status")

# Snake 3: merge
cons <- merge(cons, s3_status, by = "global_window_id", all.x = TRUE)
cons[is.na(snake3_status), snake3_status := "FAIL"]
setnames(cons, "snake3_status", "ghsl_status")

# Support flags
cons[, pca_support := pca_status == "PASS"]
cons[, group_support := group_status == "PASS"]
cons[, ghsl_support := ghsl_status == "PASS"]
cons[, n_snakes_supporting := as.integer(pca_support) + as.integer(group_support) + as.integer(ghsl_support)]

# Support signature
cons[, support_signature := {
  sig <- character(.N)
  sig[pca_support & group_support & ghsl_support] <- "PCA+GROUP+GHSL"
  sig[pca_support & group_support & !ghsl_support] <- "PCA+GROUP"
  sig[pca_support & !group_support & ghsl_support] <- "PCA+GHSL"
  sig[!pca_support & group_support & ghsl_support] <- "GROUP+GHSL"
  sig[pca_support & !group_support & !ghsl_support] <- "PCA_ONLY"
  sig[!pca_support & group_support & !ghsl_support] <- "GROUP_ONLY"
  sig[!pca_support & !group_support & ghsl_support] <- "GHSL_ONLY"
  sig[!pca_support & !group_support & !ghsl_support] <- "NONE"
  sig
}]

setorder(cons, chrom, start_bp)

message("[STEP10i] Signature distribution:")
sig_tab <- cons[, .N, by = support_signature][order(-N)]
for (i in seq_len(nrow(sig_tab))) {
  message("  ", sig_tab$support_signature[i], ": ", sig_tab$N[i])
}

# =============================================================================
# BUILD CONSENSUS REGIONS
# =============================================================================

message("[STEP10i] Building consensus regions...")

build_regions_chr <- function(cdt) {
  n <- nrow(cdt)
  if (n == 0) return(NULL)

  regions <- list()
  current <- NULL
  gap <- 0L

  for (i in seq_len(n)) {
    if (cdt$n_snakes_supporting[i] > 0) {
      if (is.null(current)) {
        current <- list(start_idx = i, end_idx = i)
      } else {
        current$end_idx <- i
      }
      gap <- 0L
    } else {
      if (!is.null(current)) {
        gap <- gap + 1L
        if (gap > GAP_TOLERANCE) {
          # Close region
          if ((current$end_idx - current$start_idx + 1) >= MIN_REGION_WINS) {
            regions[[length(regions) + 1]] <- current
          }
          current <- NULL
          gap <- 0L
        }
      }
    }
  }

  # Close final region
  if (!is.null(current) &&
      (current$end_idx - current$start_idx + 1) >= MIN_REGION_WINS) {
    regions[[length(regions) + 1]] <- current
  }

  regions
}

region_rows <- list()
region_id <- 0L

for (chr in unique(cons$chrom)) {
  cdt <- cons[chrom == chr]
  regs <- build_regions_chr(cdt)

  if (is.null(regs) || length(regs) == 0) next

  for (reg in regs) {
    region_id <- region_id + 1L
    sub <- cdt[reg$start_idx:reg$end_idx]

    # Span breakdown
    three_way <- sub[n_snakes_supporting == 3]
    two_way   <- sub[n_snakes_supporting >= 2]
    one_way   <- sub[n_snakes_supporting >= 1]

    # Dominant signature
    sig_counts <- sub[, .N, by = support_signature]
    dominant <- sig_counts[which.max(N)]$support_signature

    # Internal discordance
    supported <- sub[n_snakes_supporting > 0]
    discordance <- if (nrow(supported) > 0) {
      1 - sum(supported$n_snakes_supporting) / (3 * nrow(supported))
    } else 1

    # Classification
    n_3way <- nrow(three_way)
    n_total <- nrow(sub)

    if (n_3way >= n_total * 0.3) {
      support_class <- "CORE_3WAY"
    } else if (nrow(two_way) >= n_total * 0.5) {
      support_class <- "SHOULDER_2WAY"
    } else if (dominant == "PCA_ONLY") {
      support_class <- "PCA_DOMINANT"
    } else if (dominant == "GROUP_ONLY") {
      support_class <- "GROUP_DOMINANT"
    } else if (discordance > 0.5) {
      support_class <- "DISCORDANT_COMPLEX"
    } else {
      support_class <- "WEAK_1WAY"
    }

    # Per-track spans
    pca_wins  <- sub[pca_support == TRUE]
    grp_wins  <- sub[group_support == TRUE]
    ghsl_wins <- sub[ghsl_support == TRUE]

    region_rows[[length(region_rows) + 1]] <- data.table(
      region_id = region_id,
      chrom = chr,
      start_bp = min(sub$start_bp),
      end_bp = max(sub$end_bp),
      n_windows = n_total,
      n_3way = n_3way,
      n_2way = nrow(two_way),
      n_1way = nrow(one_way),

      core_3way_start = if (n_3way > 0) min(three_way$start_bp) else NA_real_,
      core_3way_end   = if (n_3way > 0) max(three_way$end_bp)   else NA_real_,

      support_class = support_class,
      dominant_signature = dominant,
      internal_discordance = round(discordance, 4),

      pca_span_start  = if (nrow(pca_wins) > 0) min(pca_wins$start_bp)   else NA_real_,
      pca_span_end    = if (nrow(pca_wins) > 0) max(pca_wins$end_bp)     else NA_real_,
      group_span_start = if (nrow(grp_wins) > 0) min(grp_wins$start_bp)  else NA_real_,
      group_span_end   = if (nrow(grp_wins) > 0) max(grp_wins$end_bp)    else NA_real_,
      ghsl_span_start  = if (nrow(ghsl_wins) > 0) min(ghsl_wins$start_bp) else NA_real_,
      ghsl_span_end    = if (nrow(ghsl_wins) > 0) max(ghsl_wins$end_bp)   else NA_real_
    )
  }
}

region_dt <- if (length(region_rows) > 0) rbindlist(region_rows) else {
  data.table(region_id = integer(), chrom = character(), support_class = character())
}

message("[STEP10i] Consensus regions: ", nrow(region_dt))
if (nrow(region_dt) > 0) {
  class_tab <- region_dt[, .N, by = support_class][order(-N)]
  for (i in seq_len(nrow(class_tab))) {
    message("  ", class_tab$support_class[i], ": ", class_tab$N[i])
  }
}

# =============================================================================
# PER-CHROMOSOME SUMMARY
# =============================================================================

chr_summary <- cons[, .(
  n_windows = .N,
  n_pca_pass = sum(pca_support),
  n_group_pass = sum(group_support),
  n_ghsl_pass = sum(ghsl_support),
  n_3way = sum(n_snakes_supporting == 3),
  n_2way = sum(n_snakes_supporting == 2),
  n_1way = sum(n_snakes_supporting == 1),
  n_none = sum(n_snakes_supporting == 0)
), by = chrom]

# =============================================================================
# WRITE OUTPUTS
# =============================================================================

f1 <- file.path(outdir, "consensus_windows.tsv.gz")
f2 <- file.path(outdir, "consensus_regions.tsv.gz")
f3 <- file.path(outdir, "consensus_summary.tsv")

fwrite(cons, f1, sep = "\t")
fwrite(region_dt, f2, sep = "\t")
fwrite(chr_summary, f3, sep = "\t")

message("\n[DONE] STEP10i 3-snake consensus complete")
message("  ", f1)
message("  ", f2)
message("  ", f3)
