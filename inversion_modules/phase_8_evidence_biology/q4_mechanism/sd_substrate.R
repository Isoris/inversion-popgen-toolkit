#!/usr/bin/env Rscript
# =============================================================================
# sd_substrate.R — NAHR substrate detection from minimap2 + BISER2 concordance
#
# Merged from (per docs/toolkit_audit.md §D.1, acted on pass 18, 2026-04-24):
#   - cheat14_repeat_architecture.R  (archived, minimap2 self-alignment)
#   - cheat14_self_align.sh          (archived, bash wrapper)
#   - cheat27_sd_nahr_substrate.R    (was function-library-only)
#
# The two cheats were always designed to produce a joint "concordance"
# verdict. Keeping them in separate files meant you had to either run
# them with external plumbing or manually re-implement the comparison.
# The merge keeps both angles and produces the agree/disagree output in
# one pass.
#
# =============================================================================
# WHAT THIS SCRIPT ANSWERS (Q4 — formation mechanism)
# =============================================================================
# "Is there an inverted segmental duplication pair flanking this candidate's
#  breakpoints?" If yes → NAHR mechanism (recurrent, can form on multiple
#  haplotype backgrounds). If no → NHEJ mechanism (single-origin, two
#  independent double-strand breaks joined).
#
# Two independent angles:
#
#   ANGLE A (cheat14 legacy): minimap2 self-alignment on ±50 kb flanks.
#     Runs de novo on the reference FASTA. Finds inverted repeats
#     directly. Output: per-breakpoint mechanism class
#     (NAHR_CANDIDATE / NAHR_POSSIBLE / NHEJ_CANDIDATE / COMPLEX).
#     Cost: ~30 s per candidate on a laptop, ~5 min for 200 candidates.
#
#   ANGLE B (cheat27 legacy): lookup into pre-computed BISER2 catalog.
#     BISER2 already catalogued all SD pairs genome-wide. This angle
#     is just interval overlap. Output: per-breakpoint BISER2 SD
#     class (strong / weak / direct_only / subthreshold / none).
#     Cost: negligible.
#
# The two angles are combined into a `concordance` column:
#   agree_nahr            both say NAHR            → high confidence
#   agree_nhej            both say no NAHR         → high confidence
#   disagree_minimap2_nahr minimap2 yes, BISER2 no → check (BISER2 miss?)
#   disagree_biser2_nahr   BISER2 yes, minimap2 no → check (minimap2 miss?)
#   minimap2_only          BISER2 catalog missing  → moderate
#   biser2_only_nahr       --no_minimap2 mode      → moderate
#   no_data                both missing            → low
#
# =============================================================================
# KEYS WRITTEN TO mechanism.schema.json (Tier-2 block)
# =============================================================================
# Direct fields:
#   sd_present_left, sd_present_right        boolean — per-side minimap2 hit
#   sd_identity_left, sd_identity_right      number  — best inverted hit identity
#   sd_length_left_bp, sd_length_right_bp    integer — best inverted hit length
#   sd_concordance                           enum    — agree_nahr / disagree / unresolved
#   mechanism_confidence                     enum    — high / medium / low
#
# Flat keys (via schema keys_extracted):
#   q4_sd_identity_left, q4_sd_identity_right
#   q4_sd_length_left_bp, q4_sd_length_right_bp
#   q4_sd_concordance
#   q4_mechanism_confidence
#
# Additional structured fields (not in flat-key spec yet, saved in block):
#   minimap2_mechanism_per_bp     per-breakpoint NAHR_CANDIDATE / NHEJ_CANDIDATE / ...
#   biser2_nahr_class_per_bp      per-breakpoint strong / weak / direct_only / none
#   concordance_per_bp            per-breakpoint concordance verdict
#   paf_path                      where the raw minimap2 PAF was saved (for audit)
#
# =============================================================================
# REGISTRY_CONTRACT
#   BLOCKS_WRITTEN:
#     - mechanism: registries/schemas/structured_block_schemas/mechanism.schema.json
#       keys: q4_sd_identity_left, q4_sd_identity_right, q4_sd_length_left_bp,
#             q4_sd_length_right_bp, q4_sd_concordance, q4_mechanism_confidence
#       status: WIRED
#       note: Fills only the SD portion of mechanism.schema.json. The junction
#             portion (q4_mechanism_class, q4_microhomology_length_bp,
#             q4_junction_orientation) is filled by junction_ref / junction_asm
#             scripts. mechanism_class is NOT written here — junction analysis
#             makes the final call.
#   KEYS_IN:
#     - q3_final_left_bp, q3_final_right_bp (from phase_6 via bp_bridge)
# =============================================================================
#
# USAGE (single candidate):
#   Rscript sd_substrate.R \
#     --candidate LG28_cand_1 \
#     --chrom C_gar_LG28 \
#     --left_bp 15115243 \
#     --right_bp 18005891 \
#     --ref_fasta /path/to/ref.fasta \
#     --biser2_tsv /path/to/biser2.tsv \
#     --outdir /per-candidate/output/dir
#
# USAGE (--no_minimap2 for fast BISER2-only mode when LANTA is busy):
#   Rscript sd_substrate.R --candidate ... --biser2_tsv ... --no_minimap2 \
#     --outdir ...
#
# REQUIRES: samtools, minimap2 on PATH (unless --no_minimap2)
# DIFFICULTY: medium (CLI + two data sources + concordance logic)
# DISPATCHER: no (pure sequence analysis, no popgen stats)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# =============================================================================
# PARAMETERS (tuned for catfish subgenomes; override via CLI)
# =============================================================================
FLANK_BP         <- 50000L   # ±50 kb flanks for both angles
NAHR_MIN_LEN     <- 500L     # minimap2: min inverted-repeat length for NAHR_CANDIDATE
NAHR_MIN_IDENT   <- 85       # minimap2: min % identity for NAHR_CANDIDATE
NAHR_POSS_MIN    <- 100L     # minimap2: min length for NAHR_POSSIBLE
SD_MIN_LEN_NAHR  <- 1000L    # BISER2:   min SD length for strong NAHR
SD_MIN_LEN_WEAK  <- 500L     # BISER2:   min SD length for weak NAHR
SD_MIN_IDENT     <- 85       # BISER2:   min % identity
SD_STRONG_IDENT  <- 90       # BISER2:   threshold for "strong" substrate

# =============================================================================
# ANGLE A — minimap2 self-alignment (was cheat14)
# =============================================================================

#' Extract ±flank regions around breakpoints from the reference FASTA
#' @return path to fasta, or NULL on failure
extract_breakpoint_flanks <- function(ref_fasta, boundary_bp, chr,
                                       flank_size = FLANK_BP,
                                       out_fasta = tempfile(fileext = ".fa")) {
  if (!file.exists(ref_fasta)) {
    message("[sd_substrate] Reference FASTA not found: ", ref_fasta)
    return(NULL)
  }
  regions <- character()
  for (bp in boundary_bp) {
    ls <- max(1, bp - flank_size); le <- bp
    rs <- bp;                       re <- bp + flank_size
    regions <- c(regions,
                 paste0(chr, ":", ls, "-", le),
                 paste0(chr, ":", rs, "-", re))
  }
  cmd <- paste0("samtools faidx ", shQuote(ref_fasta), " ",
                paste(regions, collapse = " "),
                " > ", shQuote(out_fasta), " 2>/dev/null")
  if (system(cmd) != 0 || file.size(out_fasta) == 0) {
    message("[sd_substrate] samtools faidx failed"); return(NULL)
  }
  out_fasta
}

#' Run minimap2 self-alignment and return the PAF as data.table
self_align_flanks <- function(flank_fasta, save_paf_to = NULL) {
  if (is.null(flank_fasta) || !file.exists(flank_fasta)) return(data.table())
  paf <- if (is.null(save_paf_to)) tempfile(fileext = ".paf") else save_paf_to
  cmd <- paste0("minimap2 -X -c --eqx ", shQuote(flank_fasta), " ",
                shQuote(flank_fasta), " > ", shQuote(paf), " 2>/dev/null")
  if (system(cmd) != 0 || !file.exists(paf) || file.size(paf) == 0) {
    message("[sd_substrate] minimap2 self-alignment failed"); return(data.table())
  }
  dt <- tryCatch(
    fread(paf, header = FALSE, fill = TRUE, select = 1:12,
          col.names = c("qname","qlen","qstart","qend","strand",
                        "tname","tlen","tstart","tend","nmatch","alen","mapq")),
    error = function(e) data.table())
  if (is.null(save_paf_to)) unlink(paf)
  if (nrow(dt) == 0) return(dt)
  dt <- dt[qname != tname | abs(qstart - tstart) > 100]  # remove self-hits
  dt <- dt[alen >= 50]                                   # remove tiny hits
  dt[, identity := round(nmatch / alen * 100, 1)]
  dt[, orientation := fifelse(strand == "+", "direct", "inverted")]
  dt
}

#' Classify mechanism per breakpoint from minimap2 PAF
#' @return data.table with one row per breakpoint
classify_mechanism_minimap2 <- function(paf_results, boundary_bp) {
  if (length(boundary_bp) == 0) {
    return(data.table(bp_pair = character(), mechanism = character(),
                      repeat_length = integer(), repeat_identity = numeric(),
                      orientation = character()))
  }
  if (nrow(paf_results) == 0) {
    return(data.table(
      bp_pair = paste0("bp_", seq_along(boundary_bp)),
      mechanism = rep("NHEJ_CANDIDATE", length(boundary_bp)),
      repeat_length = 0L, repeat_identity = 0, orientation = "none"))
  }
  inv_hits <- paf_results[orientation == "inverted"]
  results <- list()
  for (i in seq_along(boundary_bp)) {
    label <- paste0("bp_", i)
    hits_i <- inv_hits[grepl(paste0("bp", i), qname) |
                         grepl(paste0("bp", i), tname)]
    if (nrow(hits_i) == 0) {
      results[[i]] <- data.table(bp_pair = label, mechanism = "NHEJ_CANDIDATE",
                                 repeat_length = 0L, repeat_identity = 0,
                                 orientation = "none")
    } else {
      best <- hits_i[which.max(alen)]
      mech <- if (best$alen >= NAHR_MIN_LEN && best$identity >= NAHR_MIN_IDENT) {
        "NAHR_CANDIDATE"
      } else if (best$alen >= NAHR_POSS_MIN) {
        "NAHR_POSSIBLE"
      } else {
        "COMPLEX_ARCHITECTURE"
      }
      results[[i]] <- data.table(bp_pair = label, mechanism = mech,
                                 repeat_length = as.integer(best$alen),
                                 repeat_identity = best$identity,
                                 orientation = best$orientation)
    }
  }
  rbindlist(results)
}

# =============================================================================
# ANGLE B — BISER2 catalog lookup (was cheat27)
# =============================================================================

#' Load BISER2 output TSV and normalize column names
load_biser2 <- function(biser2_file) {
  if (!file.exists(biser2_file)) {
    message("[sd_substrate] BISER2 file not found: ", biser2_file)
    return(data.table())
  }
  dt <- fread(biser2_file)
  expected <- c("chr1", "start1", "end1", "chr2", "start2", "end2",
                "identity", "orientation")
  if (!all(expected[1:6] %in% names(dt))) {
    if (ncol(dt) >= 8) {
      setnames(dt, seq_len(min(8, ncol(dt))),
               expected[seq_len(min(8, ncol(dt)))])
    } else {
      message("[sd_substrate] Cannot parse BISER2 format: ", ncol(dt), " columns")
      return(data.table())
    }
  }
  dt[, sd_length := pmax(abs(end1 - start1), abs(end2 - start2))]
  dt
}

#' Find BISER2 SD pairs where one copy sits in the ±flank of each breakpoint
find_flanking_sds <- function(sd_dt, chr, boundary_bp, flank = FLANK_BP) {
  if (nrow(sd_dt) == 0) return(data.table())
  left_start  <- boundary_bp - flank
  left_end    <- boundary_bp
  right_start <- boundary_bp
  right_end   <- boundary_bp + flank
  hits <- sd_dt[
    chr1 == chr & chr2 == chr &
      ((start1 >= left_start & end1 <= left_end & start2 >= right_start & end2 <= right_end) |
       (start2 >= left_start & end2 <= left_end & start1 >= right_start & end1 <= right_end))
  ]
  if (nrow(hits) == 0) return(data.table())
  hits[, `:=`(
    boundary_bp = boundary_bp,
    is_inverted = grepl("inv", tolower(orientation)),
    straddles_breakpoint = TRUE
  )]
  hits
}

#' Classify the BISER2 SD context per breakpoint
classify_sd_mechanism_biser2 <- function(flanking_sds) {
  if (nrow(flanking_sds) == 0) {
    return(list(biser2_has_sd = FALSE, biser2_orientation = "none",
                biser2_identity = NA_real_, biser2_length = NA_integer_,
                biser2_nahr_class = "none"))
  }
  inv_sds <- flanking_sds[is_inverted == TRUE]
  if (nrow(inv_sds) > 0) {
    best <- inv_sds[which.max(sd_length)]
    ident <- best$identity[1]
    len <- best$sd_length[1]
    substrate <- if (len >= SD_MIN_LEN_NAHR && ident >= SD_STRONG_IDENT) "strong"
                 else if (len >= SD_MIN_LEN_WEAK && ident >= SD_MIN_IDENT) "weak"
                 else "subthreshold"
    return(list(biser2_has_sd = TRUE,
                biser2_orientation = "inverted",
                biser2_identity = round(ident, 1),
                biser2_length = as.integer(len),
                biser2_nahr_class = substrate))
  }
  best_dir <- flanking_sds[which.max(sd_length)]
  list(biser2_has_sd = TRUE,
       biser2_orientation = "direct",
       biser2_identity = round(best_dir$identity[1], 1),
       biser2_length = as.integer(best_dir$sd_length[1]),
       biser2_nahr_class = "direct_only")
}

# =============================================================================
# CONCORDANCE — combine ANGLE A and ANGLE B
# =============================================================================

#' Decide concordance between minimap2 verdict and BISER2 verdict per breakpoint
#' @return concordance string + confidence string
combine_angles <- function(minimap2_mech, biser2_nahr_class) {
  # Normalize to "nahr" / "nhej" / "unknown"
  mm_says <- if (is.na(minimap2_mech) || is.null(minimap2_mech)) "unknown"
             else if (grepl("NAHR", minimap2_mech)) "nahr"
             else if (minimap2_mech == "NHEJ_CANDIDATE") "nhej"
             else "unknown"
  b2_says <- if (is.na(biser2_nahr_class) || is.null(biser2_nahr_class)) "unknown"
             else if (biser2_nahr_class %in% c("strong", "weak")) "nahr"
             else if (biser2_nahr_class %in% c("none", "direct_only", "subthreshold")) "nhej"
             else "unknown"

  concordance <- if (mm_says == "unknown" && b2_says == "unknown") "no_data"
                 else if (mm_says == "unknown") {
                   if (b2_says == "nahr") "biser2_only_nahr" else "biser2_only_nhej"
                 } else if (b2_says == "unknown") {
                   if (mm_says == "nahr") "minimap2_only_nahr" else "minimap2_only_nhej"
                 } else if (mm_says == "nahr" && b2_says == "nahr") "agree_nahr"
                 else if (mm_says == "nhej" && b2_says == "nhej") "agree_nhej"
                 else if (mm_says == "nahr" && b2_says == "nhej") "disagree_minimap2_nahr"
                 else if (mm_says == "nhej" && b2_says == "nahr") "disagree_biser2_nahr"
                 else "inconclusive"

  confidence <- if (concordance %in% c("agree_nahr", "agree_nhej")) "high"
                else if (grepl("_only_", concordance)) "medium"
                else if (grepl("disagree", concordance)) "low"  # flagged, needs check
                else "low"

  list(concordance = concordance, confidence = confidence)
}

# =============================================================================
# RUNNER — per-candidate
# =============================================================================

#' Run both angles for one candidate and return a structured result
#' @param chr chromosome
#' @param boundary_bps integer vector of breakpoints (typically left + right)
#' @param ref_fasta path to reference fasta (NULL to skip minimap2)
#' @param biser2_tsv path to BISER2 output TSV (NULL to skip BISER2)
#' @param paf_dir where to save the minimap2 PAF (NULL → temp + deleted)
#' @param run_minimap2 logical — if FALSE, skip Angle A entirely
run_sd_substrate <- function(chr, boundary_bps, ref_fasta = NULL,
                              biser2_tsv = NULL, paf_dir = NULL,
                              run_minimap2 = TRUE) {
  message("[sd_substrate] ", chr, ": ", length(boundary_bps),
          " breakpoints (", paste(boundary_bps, collapse = ","), ")")

  # ── Angle A: minimap2 ──────────────────────────────────────────
  mm_mech <- data.table()
  paf_path <- NA_character_
  if (run_minimap2 && !is.null(ref_fasta) && file.exists(ref_fasta)) {
    fa <- extract_breakpoint_flanks(ref_fasta, boundary_bps, chr)
    if (!is.null(fa)) {
      paf_save <- if (!is.null(paf_dir)) {
        dir.create(paf_dir, recursive = TRUE, showWarnings = FALSE)
        file.path(paf_dir,
                  sprintf("self_align_%s_%s.paf",
                          chr, paste(boundary_bps, collapse = "_")))
      } else NULL
      paf <- self_align_flanks(fa, save_paf_to = paf_save)
      if (!is.null(paf_save) && file.exists(paf_save)) paf_path <- paf_save
      message("[sd_substrate] minimap2: ", nrow(paf), " hits (",
              sum(paf$orientation == "inverted"), " inverted)")
      mm_mech <- classify_mechanism_minimap2(paf, boundary_bps)
    }
  } else {
    message("[sd_substrate] minimap2 skipped (", 
            if (!run_minimap2) "--no_minimap2"
            else if (is.null(ref_fasta)) "no --ref_fasta"
            else "ref FASTA missing", ")")
  }

  # ── Angle B: BISER2 lookup ─────────────────────────────────────
  b2_per_bp <- list()
  if (!is.null(biser2_tsv) && file.exists(biser2_tsv)) {
    sd_dt <- load_biser2(biser2_tsv)
    if (nrow(sd_dt) > 0) {
      for (i in seq_along(boundary_bps)) {
        flanking <- find_flanking_sds(sd_dt, chr, boundary_bps[i])
        b2_per_bp[[i]] <- classify_sd_mechanism_biser2(flanking)
      }
    } else {
      message("[sd_substrate] BISER2 file parsed but empty")
    }
  } else {
    message("[sd_substrate] BISER2 skipped (no --biser2_tsv or file missing)")
  }

  # ── Combine per-breakpoint ─────────────────────────────────────
  per_bp <- list()
  for (i in seq_along(boundary_bps)) {
    mm <- if (nrow(mm_mech) >= i) mm_mech$mechanism[i] else NA_character_
    b2 <- if (length(b2_per_bp) >= i) b2_per_bp[[i]]$biser2_nahr_class else NA_character_
    combo <- combine_angles(mm, b2)

    per_bp[[i]] <- list(
      bp                    = boundary_bps[i],
      minimap2_mechanism    = mm %||% NA_character_,
      minimap2_repeat_len   = if (nrow(mm_mech) >= i) mm_mech$repeat_length[i] else NA_integer_,
      minimap2_identity     = if (nrow(mm_mech) >= i) mm_mech$repeat_identity[i] else NA_real_,
      biser2_nahr_class     = b2 %||% NA_character_,
      biser2_orientation    = if (length(b2_per_bp) >= i) b2_per_bp[[i]]$biser2_orientation else NA_character_,
      biser2_identity       = if (length(b2_per_bp) >= i) b2_per_bp[[i]]$biser2_identity    else NA_real_,
      biser2_length         = if (length(b2_per_bp) >= i) b2_per_bp[[i]]$biser2_length      else NA_integer_,
      concordance           = combo$concordance,
      confidence            = combo$confidence
    )
  }

  # ── Candidate-level summary for schema fields ──────────────────
  # Left = first boundary, Right = last boundary
  first <- per_bp[[1]]
  last  <- per_bp[[length(per_bp)]]
  # Overall concordance: if all per-bp agree on same class, adopt; else "disagree"
  all_conc <- sapply(per_bp, function(x) x$concordance)
  overall_concordance <- if (all(all_conc == "agree_nahr")) "agree_nahr"
                         else if (all(all_conc == "agree_nhej")) "agree_nhej"
                         else if (any(grepl("disagree", all_conc))) "disagree"
                         else "unresolved"
  all_conf <- sapply(per_bp, function(x) x$confidence)
  overall_confidence <- if (any(all_conf == "low")) "low"
                        else if (all(all_conf == "high")) "high"
                        else "medium"

  list(
    block_type = "mechanism",
    # Per-side SD keys the schema expects
    sd_present_left        = first$minimap2_mechanism %in% c("NAHR_CANDIDATE", "NAHR_POSSIBLE") ||
                             first$biser2_nahr_class %in% c("strong", "weak"),
    sd_present_right       = last$minimap2_mechanism %in% c("NAHR_CANDIDATE", "NAHR_POSSIBLE") ||
                             last$biser2_nahr_class %in% c("strong", "weak"),
    sd_identity_left       = first$minimap2_identity %||% first$biser2_identity,
    sd_identity_right      = last$minimap2_identity  %||% last$biser2_identity,
    sd_length_left_bp      = first$minimap2_repeat_len %||% first$biser2_length,
    sd_length_right_bp     = last$minimap2_repeat_len  %||% last$biser2_length,
    sd_concordance         = overall_concordance,
    mechanism_confidence   = overall_confidence,
    # Structured detail (not in flat-key spec)
    minimap2_mechanism_per_bp = sapply(per_bp, function(x) x$minimap2_mechanism),
    biser2_nahr_class_per_bp  = sapply(per_bp, function(x) x$biser2_nahr_class),
    concordance_per_bp        = all_conc,
    paf_path                  = paf_path,
    # Provenance
    ran_minimap2 = run_minimap2 && nrow(mm_mech) > 0,
    ran_biser2   = length(b2_per_bp) > 0
  )
}

# Small helper because R has no %||% by default in base
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a[1]))) b else a

# =============================================================================
# CLI
# =============================================================================
parse_args <- function(argv) {
  opts <- list(
    candidate   = NA_character_,
    chrom       = NA_character_,
    left_bp     = NA_integer_,
    right_bp    = NA_integer_,
    ref_fasta   = NA_character_,
    biser2_tsv  = NA_character_,
    outdir      = NA_character_,
    paf_dir     = NA_character_,
    no_minimap2 = FALSE,
    registries_root = NA_character_
  )
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    v <- if (i < length(argv)) argv[i + 1] else NA_character_
    switch(a,
      "--candidate"       = { opts$candidate  <- v; i <- i + 2 },
      "--chrom"           = { opts$chrom      <- v; i <- i + 2 },
      "--left_bp"         = { opts$left_bp    <- as.integer(v); i <- i + 2 },
      "--right_bp"        = { opts$right_bp   <- as.integer(v); i <- i + 2 },
      "--ref_fasta"       = { opts$ref_fasta  <- v; i <- i + 2 },
      "--biser2_tsv"      = { opts$biser2_tsv <- v; i <- i + 2 },
      "--outdir"          = { opts$outdir     <- v; i <- i + 2 },
      "--paf_dir"         = { opts$paf_dir    <- v; i <- i + 2 },
      "--no_minimap2"     = { opts$no_minimap2 <- TRUE; i <- i + 1 },
      "--registries_root" = { opts$registries_root <- v; i <- i + 2 },
      "--help"            = { cat(readLines(sys.frame(1)$ofile)[1:100], sep = "\n"); quit(status = 0) },
      { i <- i + 1 }
    )
  }
  opts
}

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  opts <- parse_args(argv)

  # Required args
  for (k in c("candidate", "chrom", "left_bp", "right_bp", "outdir")) {
    if (is.na(opts[[k]])) {
      stop("[sd_substrate] missing required --", k)
    }
  }

  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  result <- run_sd_substrate(
    chr          = opts$chrom,
    boundary_bps = c(opts$left_bp, opts$right_bp),
    ref_fasta    = if (is.na(opts$ref_fasta)) NULL else opts$ref_fasta,
    biser2_tsv   = if (is.na(opts$biser2_tsv)) NULL else opts$biser2_tsv,
    paf_dir      = if (is.na(opts$paf_dir)) NULL else opts$paf_dir,
    run_minimap2 = !opts$no_minimap2
  )

  # Write structured block (jsonlite handles NA properly as null)
  out_struct <- file.path(opts$outdir, "structured")
  dir.create(out_struct, recursive = TRUE, showWarnings = FALSE)
  payload <- list(
    block_type    = "mechanism",
    candidate_id  = opts$candidate,
    source_script = "sd_substrate.R",
    data          = result
  )
  out_file <- file.path(out_struct, "mechanism_sd_substrate.json")
  writeLines(toJSON(payload, auto_unbox = TRUE, null = "null",
                    na = "null", pretty = TRUE),
             out_file)
  message("[sd_substrate] wrote ", out_file)

  # Registry write (optional)
  if (!is.na(opts$registries_root)) {
    reg_r <- file.path(opts$registries_root, "api", "R", "registry_loader.R")
    if (file.exists(reg_r)) {
      tryCatch({
        source(reg_r)
        reg <- load_registry()
        reg$evidence$write_block(
          candidate_id = opts$candidate,
          block_type   = "mechanism",
          data         = result,
          source_script = "sd_substrate.R"
        )
        message("[sd_substrate] wrote registry block for ", opts$candidate)
      }, error = function(e) {
        message("[sd_substrate] registry write failed: ", e$message)
      })
    }
  }

  message("[sd_substrate] ", opts$candidate, ": concordance=",
          result$sd_concordance, " confidence=", result$mechanism_confidence)
}

# Only run main if invoked via Rscript (not if sourced interactively)
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  # Detect whether this is Rscript vs interactive source()
  if (length(commandArgs(trailingOnly = FALSE)) > 0 &&
      any(grepl("--file=", commandArgs(trailingOnly = FALSE)))) {
    main()
  }
}
