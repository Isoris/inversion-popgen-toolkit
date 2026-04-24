#!/usr/bin/env Rscript
# =============================================================================
# cheat27_sd_nahr_substrate.R — BISER2 SD concordance check for NAHR
#
# ROLE: CONCORDANCE CHECK — not primary mechanism classifier.
#   Primary mechanism call comes from Cheat 14 (minimap2 self-alignment).
#   This cheat cross-checks against pre-computed BISER2 segmental duplications.
#
#   minimap2 (Cheat 14) + BISER2 agree → HIGH CONFIDENCE NAHR.
#   They disagree → FLAG for manual inspection.
#   BISER2 alone NEVER overrides minimap2.
#
# BIOLOGY:
#   NAHR requires inverted segmental duplications flanking breakpoints.
#   BISER2 already computed all SD pairs genome-wide. Cross-checking adds
#   confidence but minimap2 is the trusted tool.
#
#   Intrachromosomal SDs (Mac: 1,493 pairs, Gar: 1,414) are the NAHR
#   substrate. Interchromosomal SDs are irrelevant for inversions.
#
# INPUT:
#   biser2_file: BISER2 output TSV
#   boundary_catalog: from C01g
#   cheat14_result: from Cheat 14 minimap2 (optional, for concordance)
#
# OUTPUT per breakpoint:
#   biser2_has_sd:        TRUE/FALSE
#   biser2_orientation:   "inverted" / "direct" / "none"
#   biser2_identity:      % identity of the SD pair
#   biser2_length:        length of the SD
#   biser2_nahr_class:    "strong" / "weak" / "direct_only" / "none"
#   concordance:          "agree_nahr" / "agree_nhej" / "disagree" / "biser2_only" / "minimap2_only"
#   mechanism_confidence: "high" (both agree) / "moderate" (one only) / "check" (disagree)
#
# INTEGRATION:
#   C01g: annotate boundaries → cheat27_* columns
#   Cheat 14: compare → concordance column
#   Cheat 15: concordance boosts/downgrades recurrence prior
#
# REQUIRES: BISER2 output file, boundary catalog
# DIFFICULTY: Easy (pure BED interval overlap)
# DISPATCHER: No (no popgen computation)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
SD_FLANK_WINDOW  <- 50000L   # ±50 kb from breakpoint
SD_MIN_LEN_NAHR  <- 1000L    # minimum SD length for strong NAHR
SD_MIN_LEN_WEAK  <- 500L     # minimum for weak NAHR
SD_MIN_IDENT     <- 85       # minimum identity %
SD_STRONG_IDENT  <- 90       # identity for "strong" substrate

# ── Load BISER2 SDs ──────────────────────────────────────────────────────

load_biser2 <- function(biser2_file) {
  if (!file.exists(biser2_file)) {
    message("[cheat27] BISER2 file not found: ", biser2_file)
    return(data.table())
  }
  dt <- fread(biser2_file)
  # Standardize column names — BISER2 output format varies
  expected <- c("chr1", "start1", "end1", "chr2", "start2", "end2",
                "identity", "orientation")
  if (!all(expected[1:6] %in% names(dt))) {
    # Try positional assignment
    if (ncol(dt) >= 8) {
      setnames(dt, seq_len(min(8, ncol(dt))),
               expected[seq_len(min(8, ncol(dt)))])
    } else {
      message("[cheat27] Cannot parse BISER2 format: ", ncol(dt), " columns")
      return(data.table())
    }
  }
  dt[, sd_length := pmax(abs(end1 - start1), abs(end2 - start2))]
  dt
}

# ── Find SD pairs flanking a breakpoint ──────────────────────────────────

find_flanking_sds <- function(sd_dt, chr, boundary_bp, flank = SD_FLANK_WINDOW) {
  if (nrow(sd_dt) == 0) return(data.table())

  left_start  <- boundary_bp - flank
  left_end    <- boundary_bp
  right_start <- boundary_bp
  right_end   <- boundary_bp + flank

  # Intrachromosomal SDs where one copy is on each side of the breakpoint
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

# ── Classify mechanism from SD context ────────────────────────────────────

classify_sd_mechanism <- function(flanking_sds) {
  if (nrow(flanking_sds) == 0) {
    return(list(biser2_has_sd = FALSE, biser2_orientation = "none",
                biser2_identity = NA_real_, biser2_length = NA_integer_,
                biser2_nahr_class = "none"))
  }

  # Best inverted SD pair
  inv_sds <- flanking_sds[is_inverted == TRUE]
  if (nrow(inv_sds) > 0) {
    best <- inv_sds[which.max(sd_length)]
    ident <- best$identity[1]
    len <- best$sd_length[1]

    substrate <- if (len >= SD_MIN_LEN_NAHR && ident >= SD_STRONG_IDENT) "strong"
                 else if (len >= SD_MIN_LEN_WEAK && ident >= SD_MIN_IDENT) "weak"
                 else "subthreshold"

    return(list(
      biser2_has_sd = TRUE,
      biser2_orientation = "inverted",
      biser2_identity = round(ident, 1),
      biser2_length = as.integer(len),
      biser2_nahr_class = substrate
    ))
  }

  # Direct SDs only
  best_dir <- flanking_sds[which.max(sd_length)]
  list(
    biser2_has_sd = TRUE,
    biser2_orientation = "direct",
    biser2_identity = round(best_dir$identity[1], 1),
    biser2_length = as.integer(best_dir$sd_length[1]),
    biser2_nahr_class = "direct_only"
  )
}

# ── Runner ────────────────────────────────────────────────────────────────

run_cheat27 <- function(chr, boundary_bps, biser2_file_or_dt,
                         cheat14_mechanisms = NULL,
                         candidate_start = NULL, candidate_end = NULL) {
  sd_dt <- if (is.data.table(biser2_file_or_dt)) biser2_file_or_dt
           else load_biser2(biser2_file_or_dt)

  if (nrow(sd_dt) == 0) {
    message("[cheat27] No SD data")
    return(data.table(boundary_bp = boundary_bps,
                       biser2_has_sd = FALSE, biser2_nahr_class = "none",
                       concordance = "no_biser2_data",
                       mechanism_confidence = "moderate"))
  }

  results <- list()
  for (i in seq_along(boundary_bps)) {
    bp <- boundary_bps[i]
    flanking <- find_flanking_sds(sd_dt, chr, bp)
    mech <- classify_sd_mechanism(flanking)

    # Concordance with Cheat 14 minimap2 result
    c14_mech <- if (!is.null(cheat14_mechanisms) && length(cheat14_mechanisms) >= i)
      cheat14_mechanisms[i] else NA_character_

    biser2_says_nahr <- mech$biser2_nahr_class %in% c("strong", "weak")
    c14_says_nahr <- !is.na(c14_mech) && grepl("NAHR|inverted", c14_mech, ignore.case = TRUE)
    c14_says_nhej <- !is.na(c14_mech) && grepl("NHEJ|blunt|microhom", c14_mech, ignore.case = TRUE)

    concordance <- if (is.na(c14_mech)) {
      if (biser2_says_nahr) "biser2_only_nahr" else "no_cheat14"
    } else if (biser2_says_nahr && c14_says_nahr) {
      "agree_nahr"
    } else if (!biser2_says_nahr && c14_says_nhej) {
      "agree_nhej"
    } else if (biser2_says_nahr && !c14_says_nahr) {
      "disagree_biser2_nahr"
    } else if (!biser2_says_nahr && c14_says_nahr) {
      "disagree_minimap2_nahr"
    } else {
      "inconclusive"
    }

    confidence <- if (concordance %in% c("agree_nahr", "agree_nhej")) "high"
                  else if (grepl("disagree", concordance)) "check"
                  else "moderate"

    results[[i]] <- c(list(boundary_bp = bp), mech,
                       list(concordance = concordance,
                            mechanism_confidence = confidence))
  }

  dt <- rbindlist(results, fill = TRUE)
  n_agree <- sum(dt$concordance %in% c("agree_nahr", "agree_nhej"))
  n_disagree <- sum(grepl("disagree", dt$concordance))
  message("[cheat27] ", chr, ": ", length(boundary_bps), " breakpoints, ",
          n_agree, " concordant, ", n_disagree, " disagreements")

  # Candidate-level summary
  if (!is.null(candidate_start)) {
    attr(dt, "candidate_summary") <- list(
      n_biser2_sd = sum(dt$biser2_has_sd),
      n_biser2_nahr = sum(dt$biser2_nahr_class %in% c("strong", "weak")),
      n_concordant = n_agree,
      n_disagree = n_disagree,
      # Mechanism from CONCORDANCE, not BISER2 alone
      concordance_mechanism = if (n_agree > 0 && any(dt$concordance == "agree_nahr")) "NAHR_confirmed"
                              else if (n_agree > 0) "NHEJ_confirmed"
                              else "unresolved"
    )
  }
  dt
}
