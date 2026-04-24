#!/usr/bin/env Rscript
# =============================================================================
# sd_substrate_biser2.R  —  Angle B of SD substrate detection
# =============================================================================
#
# ROLE:
#   Catalog lookup into the pre-computed BISER2 segmental duplication
#   catalog. For the two breakpoints of a candidate, finds SD pairs where
#   one copy is within ±flank_bp of the left breakpoint AND the other
#   copy is within ±flank_bp of the right breakpoint. Classifies based on
#   orientation (inverted vs direct) and size/identity thresholds.
#
#   Lives at:
#     q4_mechanism/sd_substrate/sd_substrate_biser2.R
#
#   Companion scripts (same directory):
#     sd_substrate_minimap2.R    — angle A (de novo alignment)
#     sd_substrate_concordance.R — reads A + B, writes joint verdict
#
# =============================================================================
# INPUT — BISER2 catalog TSV
# =============================================================================
# BISER2 output format (per the wiki): 9 columns, tab-separated
#   chr1 start1 end1 chr2 start2 end2 strand_orientation identity SD_type
# or 8 columns without trailing type. This script handles both by
# checking column count.
#
# For NAHR-substrate detection we need:
#   - intrachromosomal SD pairs (chr1 == chr2 == candidate's chrom)
#   - one copy in the left breakpoint's flank, other in the right breakpoint's flank
#   - "inverted" orientation tag is what matters; direct SDs are logged
#     but don't support NAHR
#
# =============================================================================
# OUTPUTS (all in <outdir>/<candidate_id>/sd_substrate/biser2/)
# =============================================================================
#   flanking_sds.tsv   BISER2 rows matching the straddle criterion
#   flanking_sds.bed   BED9 with itemRgb:
#                        red (255,0,0) = inverted SD (NAHR substrate)
#                        green (0,128,0) = direct SD
#
# And:
#   <outdir>/<candidate_id>/structured/sd_substrate_biser2.json
#     Tier-2 block
#
# =============================================================================
# SCHEMA BLOCK WRITTEN: sd_substrate_biser2
# =============================================================================
# Flat keys:
#   q4b_biser2_ran                bool   (catalog was parseable)
#   q4b_biser2_n_flanking_sds     int    straddling SDs (both sides)
#   q4b_biser2_n_inverted         int    inverted subset
#   q4b_biser2_n_direct           int    direct subset
#   q4b_biser2_best_inv_length    int    bp — largest inverted SD
#   q4b_biser2_best_inv_identity  float  % identity of that SD
#   q4b_biser2_nahr_class         strong | weak | subthreshold |
#                                 direct_only | none
#   q4b_biser2_catalog_path       string (for audit)
#
# =============================================================================
# CLI
# =============================================================================
#   Rscript sd_substrate_biser2.R \
#     --candidate LG28_cand_1 \
#     --chrom C_gar_LG28 \
#     --left_bp 15115243 \
#     --right_bp 18005891 \
#     --biser2_tsv /path/to/biser2_results.tsv \
#     --outdir /per-candidate/output/root \
#     [--flank_bp 50000]        # default
#     [--overwrite]
#     [--registries_root ...]
#
# REQUIRES: nothing external (pure TSV + interval overlap)
# DIFFICULTY: easy
# DISPATCHER: no
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# --------------------------------------------------------------------------
# Tunables (override via CLI)
# --------------------------------------------------------------------------
DEFAULT_FLANK_BP  <- 50000L
SD_MIN_LEN_NAHR   <- 1000L
SD_MIN_LEN_WEAK   <- 500L
SD_MIN_IDENT      <- 85
SD_STRONG_IDENT   <- 90

log_msg <- function(...) message("[sd_substrate_b2] ", ...)

# --------------------------------------------------------------------------
# BISER2 loader — tolerates minor format variants
# --------------------------------------------------------------------------

load_biser2 <- function(biser2_file) {
  if (!file.exists(biser2_file)) {
    stop("[sd_substrate_biser2] BISER2 file not found: ", biser2_file)
  }
  dt <- fread(biser2_file, header = FALSE)
  nc <- ncol(dt)

  if (nc == 9) {
    setnames(dt, c("chr1","start1","end1","chr2","start2","end2",
                    "orientation","identity","sd_type"))
  } else if (nc == 8) {
    setnames(dt, c("chr1","start1","end1","chr2","start2","end2",
                    "orientation","identity"))
    dt[, sd_type := NA_character_]
  } else {
    # Try reading with header
    dt <- fread(biser2_file)
    if (!all(c("chr1","start1","end1","chr2","start2","end2") %in% names(dt))) {
      stop("[sd_substrate_biser2] cannot parse BISER2 format: ", nc, " columns ",
           "and no standard header — check the file")
    }
  }

  # Coerce types defensively
  for (col in c("start1","end1","start2","end2")) {
    dt[, (col) := as.integer(get(col))]
  }
  dt[, identity := as.numeric(identity)]

  dt[, sd_length := pmax(abs(end1 - start1), abs(end2 - start2))]
  # Normalize orientation: BISER2 sometimes writes "+/-" or "inv" / "dir"
  dt[, is_inverted := grepl("inv|^-$|^-/|-/-", tolower(orientation))]

  dt
}

# --------------------------------------------------------------------------
# Find SDs straddling the candidate's two breakpoints
# --------------------------------------------------------------------------

find_straddling_sds <- function(sd_dt, chr, left_bp, right_bp, flank) {
  if (nrow(sd_dt) == 0) return(data.table())

  # For an inversion to have an SD-NAHR substrate, we need:
  #   one SD copy in the left-breakpoint flank zone   [left-flank, left+flank]
  #   other copy in the right-breakpoint flank zone   [right-flank, right+flank]
  # The BISER2 row has two copies (1, 2); either order is possible.

  left_lo  <- max(1L, left_bp - flank);  left_hi  <- left_bp + flank
  right_lo <- max(1L, right_bp - flank); right_hi <- right_bp + flank

  intrachrom <- sd_dt[chr1 == chr & chr2 == chr]
  if (nrow(intrachrom) == 0) return(data.table())

  # One copy in left zone, other in right zone (either assignment)
  hits <- intrachrom[
    ((start1 <= left_hi  & end1 >= left_lo  & start2 <= right_hi & end2 >= right_lo) |
     (start2 <= left_hi  & end2 >= left_lo  & start1 <= right_hi & end1 >= right_lo))
  ]
  hits
}

classify_straddling <- function(flanking) {
  if (nrow(flanking) == 0) {
    return(list(nahr_class = "none",
                n_flanking = 0L, n_inverted = 0L, n_direct = 0L,
                best_inv_length = NA_integer_, best_inv_identity = NA_real_))
  }
  inv <- flanking[is_inverted == TRUE]
  dir <- flanking[is_inverted == FALSE]

  if (nrow(inv) == 0) {
    return(list(nahr_class = "direct_only",
                n_flanking = nrow(flanking),
                n_inverted = 0L, n_direct = nrow(dir),
                best_inv_length = NA_integer_, best_inv_identity = NA_real_))
  }

  best <- inv[which.max(sd_length)]
  len <- best$sd_length[1]
  ident <- best$identity[1]
  nahr_class <- if (len >= SD_MIN_LEN_NAHR && ident >= SD_STRONG_IDENT) {
    "strong"
  } else if (len >= SD_MIN_LEN_WEAK && ident >= SD_MIN_IDENT) {
    "weak"
  } else {
    "subthreshold"
  }

  list(nahr_class = nahr_class,
       n_flanking = nrow(flanking),
       n_inverted = nrow(inv),
       n_direct = nrow(dir),
       best_inv_length = as.integer(len),
       best_inv_identity = round(as.numeric(ident), 2))
}

# --------------------------------------------------------------------------
# BED writer
# --------------------------------------------------------------------------

write_flanking_bed <- function(flanking, chr, bed_out) {
  lines <- c(
    "track name=\"sd_substrate_biser2\" description=\"BISER2 SDs straddling the candidate; red=inverted green=direct\" itemRgb=On"
  )
  if (nrow(flanking) > 0) {
    for (i in seq_len(nrow(flanking))) {
      r <- flanking[i]
      color <- if (r$is_inverted) "255,0,0" else "0,128,0"
      strand <- if (r$is_inverted) "-" else "+"
      score <- as.integer(pmin(1000, round(r$identity * 10)))
      # Emit two BED records per SD pair — one per copy
      name1 <- sprintf("SD%d_copy1_%.1f%%", i, r$identity)
      name2 <- sprintf("SD%d_copy2_%.1f%%", i, r$identity)
      lines <- c(lines, sprintf("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s",
                                chr, r$start1, r$end1, name1, score, strand,
                                r$start1, r$end1, color))
      lines <- c(lines, sprintf("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s",
                                chr, r$start2, r$end2, name2, score, strand,
                                r$start2, r$end2, color))
    }
  }
  writeLines(lines, bed_out)
}

# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

parse_args <- function(argv) {
  opts <- list(
    candidate = NA_character_,
    chrom = NA_character_,
    left_bp = NA_integer_,
    right_bp = NA_integer_,
    biser2_tsv = NA_character_,
    outdir = NA_character_,
    flank_bp = DEFAULT_FLANK_BP,
    overwrite = FALSE,
    registries_root = NA_character_
  )
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]; v <- if (i < length(argv)) argv[i + 1] else NA_character_
    switch(a,
      "--candidate"       = { opts$candidate <- v; i <- i + 2 },
      "--chrom"           = { opts$chrom <- v; i <- i + 2 },
      "--left_bp"         = { opts$left_bp <- as.integer(v); i <- i + 2 },
      "--right_bp"        = { opts$right_bp <- as.integer(v); i <- i + 2 },
      "--biser2_tsv"      = { opts$biser2_tsv <- v; i <- i + 2 },
      "--outdir"          = { opts$outdir <- v; i <- i + 2 },
      "--flank_bp"        = { opts$flank_bp <- as.integer(v); i <- i + 2 },
      "--overwrite"       = { opts$overwrite <- TRUE; i <- i + 1 },
      "--registries_root" = { opts$registries_root <- v; i <- i + 2 },
      { i <- i + 1 }
    )
  }
  opts
}

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  opts <- parse_args(argv)
  for (k in c("candidate", "chrom", "left_bp", "right_bp",
              "biser2_tsv", "outdir")) {
    if (is.na(opts[[k]])) stop("[sd_substrate_biser2] missing required --", k)
  }

  cand_dir <- file.path(opts$outdir, opts$candidate)
  b2_dir   <- file.path(cand_dir, "sd_substrate", "biser2")
  str_dir  <- file.path(cand_dir, "structured")
  dir.create(b2_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(str_dir, recursive = TRUE, showWarnings = FALSE)

  flanking_tsv <- file.path(b2_dir, "flanking_sds.tsv")
  flanking_bed <- file.path(b2_dir, "flanking_sds.bed")

  log_msg(opts$candidate, " — loading BISER2 catalog: ", opts$biser2_tsv)
  sd_dt <- load_biser2(opts$biser2_tsv)
  log_msg("loaded ", nrow(sd_dt), " SD pairs total (",
          sum(sd_dt$is_inverted), " inverted, ",
          sum(!sd_dt$is_inverted), " direct)")

  flanking <- find_straddling_sds(sd_dt, opts$chrom, opts$left_bp,
                                   opts$right_bp, opts$flank_bp)
  log_msg(opts$candidate, " — ", nrow(flanking),
          " SD pairs straddle the ±", opts$flank_bp, " bp flanks (",
          sum(flanking$is_inverted %||% FALSE), " inverted)")

  fwrite(flanking, flanking_tsv, sep = "\t")
  write_flanking_bed(flanking, opts$chrom, flanking_bed)

  cls <- classify_straddling(flanking)

  block <- list(
    block_type    = "sd_substrate_biser2",
    candidate_id  = opts$candidate,
    source_script = "sd_substrate_biser2.R",
    data = list(
      q4b_biser2_ran                = TRUE,
      q4b_biser2_flank_bp           = opts$flank_bp,
      q4b_biser2_n_flanking_sds     = as.integer(cls$n_flanking),
      q4b_biser2_n_inverted         = as.integer(cls$n_inverted),
      q4b_biser2_n_direct           = as.integer(cls$n_direct),
      q4b_biser2_best_inv_length    = cls$best_inv_length,
      q4b_biser2_best_inv_identity  = cls$best_inv_identity,
      q4b_biser2_nahr_class         = cls$nahr_class,
      q4b_biser2_catalog_path       = normalizePath(opts$biser2_tsv, mustWork = FALSE),
      q4b_biser2_flanking_tsv_path  = normalizePath(flanking_tsv, mustWork = FALSE),
      q4b_biser2_flanking_bed_path  = normalizePath(flanking_bed, mustWork = FALSE)
    )
  )
  out_json <- file.path(str_dir, "sd_substrate_biser2.json")
  writeLines(toJSON(block, auto_unbox = TRUE, null = "null", na = "null",
                    pretty = TRUE), out_json)
  log_msg("wrote ", out_json)

  if (!is.na(opts$registries_root)) {
    reg_r <- file.path(opts$registries_root, "api", "R", "registry_loader.R")
    if (file.exists(reg_r)) {
      tryCatch({
        source(reg_r)
        reg <- load_registry()
        reg$evidence$write_block(
          candidate_id  = opts$candidate,
          block_type    = "sd_substrate_biser2",
          data          = block$data,
          source_script = "sd_substrate_biser2.R"
        )
        log_msg("registered block for ", opts$candidate)
      }, error = function(e) {
        log_msg("registry write failed: ", e$message)
      })
    }
  }

  log_msg(opts$candidate,
          " — nahr_class=", cls$nahr_class,
          " n_inv=", cls$n_inverted,
          " n_dir=", cls$n_direct)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

invoked_as_script <- any(grepl("--file=", commandArgs(trailingOnly = FALSE)))
if (invoked_as_script) main()
