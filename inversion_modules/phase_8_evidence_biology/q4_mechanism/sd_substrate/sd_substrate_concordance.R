#!/usr/bin/env Rscript
# =============================================================================
# sd_substrate_concordance.R  —  joint verdict from Angle A + Angle B
# =============================================================================
#
# ROLE:
#   Reads the two blocks written by:
#     sd_substrate_minimap2.R  (Angle A — TRUSTED minimap2 de novo)
#     sd_substrate_biser2.R    (Angle B — BISER2 catalog lookup)
#   and writes a third block with the joint concordance verdict.
#
#   Angle A owns the primary NAHR call. Angle B is a catalog cross-check.
#   When they agree the confidence is high; when they disagree, Angle A
#   wins for the `mechanism_confidence` key but the disagreement is
#   surfaced so manuscript-level review can look at it.
#
#   Lives at:
#     q4_mechanism/sd_substrate/sd_substrate_concordance.R
#
# =============================================================================
# CONCORDANCE MATRIX
# =============================================================================
#
#   Angle A says:            Angle B says:    Concordance:            Confidence:
#   NAHR_CANDIDATE           strong|weak      agree_nahr              high
#   NAHR_CANDIDATE           subthreshold     minimap2_stronger_nahr  medium
#   NAHR_CANDIDATE           direct_only|none disagree_biser2_missing high  (A wins)
#   NHEJ_CANDIDATE           none|direct_only agree_nhej              high
#   NHEJ_CANDIDATE           strong|weak      disagree_minimap2_miss  medium (review)
#   NAHR_POSSIBLE            strong|weak      agree_nahr_weak         medium
#   NAHR_POSSIBLE            none             minimap2_only_weak      low
#   COMPLEX_ARCHITECTURE     any              complex_flagged         low
#   (A missing)              strong|weak      biser2_only_nahr        medium
#   (A missing)              none|direct_only biser2_only_nhej        medium
#   (B missing)              NAHR_*           minimap2_only_nahr      medium
#   (B missing)              NHEJ_*           minimap2_only_nhej      medium
#   (both missing)                            no_data                 low
#
# Confidence is per-side first, then the worst of left and right becomes
# the candidate-level confidence.
#
# =============================================================================
# SCHEMA BLOCK WRITTEN: sd_substrate_concordance
# =============================================================================
# Flat keys:
#   q4_sd_concordance_left          string (enum above)
#   q4_sd_concordance_right         string
#   q4_sd_concordance_overall       string (agree_nahr / agree_nhej /
#                                           disagree / *_only / complex_flagged /
#                                           no_data)
#   q4_sd_confidence_left           high | medium | low
#   q4_sd_confidence_right          high | medium | low
#   q4_sd_confidence_overall        worst of {left, right}
#   q4_sd_primary_call_source       minimap2 | biser2 | none  (which angle
#                                   drove the q4_sd_concordance_overall call)
#   q4_sd_primary_call              NAHR | NHEJ | COMPLEX | UNCERTAIN
#
# This block ALSO mirrors into the existing mechanism.schema.json keys
# for backward compatibility with STEP_04_assign_structural_class.py:
#   q4_sd_concordance               = q4_sd_concordance_overall
#   q4_mechanism_confidence         = q4_sd_confidence_overall
#
# =============================================================================
# CLI
# =============================================================================
#   Rscript sd_substrate_concordance.R \
#     --candidate LG28_cand_1 \
#     --outdir /per-candidate/output/root
#     [--registries_root ...]
#
# Reads <outdir>/<candidate>/structured/sd_substrate_minimap2.json
#   and <outdir>/<candidate>/structured/sd_substrate_biser2.json
# Writes <outdir>/<candidate>/structured/sd_substrate_concordance.json
#
# REQUIRES: nothing external
# DIFFICULTY: easy (rule table on two JSON blocks)
# DISPATCHER: no
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

log_msg <- function(...) message("[sd_substrate_cc] ", ...)

# --------------------------------------------------------------------------
# Concordance rules (per side)
# --------------------------------------------------------------------------

classify_concordance_side <- function(mm2_mech, b2_class) {
  # mm2_mech ∈ {NAHR_CANDIDATE, NAHR_POSSIBLE, NHEJ_CANDIDATE,
  #             COMPLEX_ARCHITECTURE, NA}
  # b2_class ∈ {strong, weak, subthreshold, direct_only, none, NA}

  mm2_missing <- is.null(mm2_mech) || is.na(mm2_mech) || mm2_mech == ""
  b2_missing  <- is.null(b2_class) || is.na(b2_class) || b2_class == ""

  if (mm2_missing && b2_missing) {
    return(list(concordance = "no_data", confidence = "low", call = "UNCERTAIN"))
  }
  if (mm2_missing) {
    if (b2_class %in% c("strong", "weak")) {
      return(list(concordance = "biser2_only_nahr", confidence = "medium",
                  call = "NAHR"))
    }
    return(list(concordance = "biser2_only_nhej", confidence = "medium",
                call = "NHEJ"))
  }
  if (b2_missing) {
    if (mm2_mech %in% c("NAHR_CANDIDATE", "NAHR_POSSIBLE")) {
      return(list(concordance = "minimap2_only_nahr", confidence = "medium",
                  call = "NAHR"))
    }
    if (mm2_mech == "COMPLEX_ARCHITECTURE") {
      return(list(concordance = "complex_flagged", confidence = "low",
                  call = "COMPLEX"))
    }
    return(list(concordance = "minimap2_only_nhej", confidence = "medium",
                call = "NHEJ"))
  }

  # Both present — full matrix
  if (mm2_mech == "COMPLEX_ARCHITECTURE") {
    return(list(concordance = "complex_flagged", confidence = "low",
                call = "COMPLEX"))
  }

  if (mm2_mech == "NAHR_CANDIDATE") {
    if (b2_class %in% c("strong", "weak")) {
      return(list(concordance = "agree_nahr", confidence = "high", call = "NAHR"))
    }
    if (b2_class == "subthreshold") {
      return(list(concordance = "minimap2_stronger_nahr", confidence = "medium",
                  call = "NAHR"))
    }
    # b2 says direct_only or none → BISER2 missed it; trust minimap2
    return(list(concordance = "disagree_biser2_missing", confidence = "high",
                call = "NAHR"))
  }

  if (mm2_mech == "NAHR_POSSIBLE") {
    if (b2_class %in% c("strong", "weak")) {
      return(list(concordance = "agree_nahr_weak", confidence = "medium",
                  call = "NAHR"))
    }
    if (b2_class == "subthreshold") {
      return(list(concordance = "agree_nahr_weak", confidence = "low",
                  call = "NAHR"))
    }
    return(list(concordance = "minimap2_only_weak", confidence = "low",
                call = "NAHR"))
  }

  if (mm2_mech == "NHEJ_CANDIDATE") {
    if (b2_class %in% c("none", "direct_only")) {
      return(list(concordance = "agree_nhej", confidence = "high", call = "NHEJ"))
    }
    if (b2_class == "subthreshold") {
      return(list(concordance = "agree_nhej", confidence = "medium", call = "NHEJ"))
    }
    # minimap2 says NHEJ, BISER2 says strong/weak NAHR — review needed
    return(list(concordance = "disagree_minimap2_missing", confidence = "medium",
                call = "UNCERTAIN"))
  }

  list(concordance = "no_data", confidence = "low", call = "UNCERTAIN")
}

# Candidate-level combiner: worst of per-side confidences
combine_sides <- function(left, right) {
  conf_rank <- c(low = 1L, medium = 2L, high = 3L)
  worst <- names(conf_rank)[min(conf_rank[c(left$confidence, right$confidence)])]

  # Overall concordance: if both sides agree on same class, adopt; else "mixed"
  if (left$concordance == right$concordance) {
    overall_conc <- left$concordance
  } else if (grepl("agree_nahr", left$concordance) && grepl("agree_nahr", right$concordance)) {
    overall_conc <- "agree_nahr"
  } else if (left$concordance == "agree_nhej" && right$concordance == "agree_nhej") {
    overall_conc <- "agree_nhej"
  } else if (grepl("disagree", left$concordance) || grepl("disagree", right$concordance)) {
    overall_conc <- "disagree"
  } else if (left$concordance == "complex_flagged" || right$concordance == "complex_flagged") {
    overall_conc <- "complex_flagged"
  } else {
    overall_conc <- "mixed"
  }

  # Overall call: if both sides say same call, use it; else UNCERTAIN
  overall_call <- if (left$call == right$call) left$call else "UNCERTAIN"

  # Primary source: minimap2 wins when both angles are present
  primary_source <- if (grepl("minimap2_only", overall_conc) ||
                        grepl("minimap2_only", left$concordance) ||
                        grepl("minimap2_only", right$concordance)) {
    if (grepl("biser2_only", overall_conc)) "none" else "minimap2"
  } else if (grepl("biser2_only", overall_conc)) {
    "biser2"
  } else if (grepl("^no_data$", overall_conc)) {
    "none"
  } else {
    "minimap2"  # default when both present: minimap2 is trusted
  }

  list(overall_conc = overall_conc, overall_conf = worst,
       overall_call = overall_call, primary_source = primary_source)
}

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

parse_args <- function(argv) {
  opts <- list(candidate = NA_character_, outdir = NA_character_,
               registries_root = NA_character_)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]; v <- if (i < length(argv)) argv[i + 1] else NA_character_
    switch(a,
      "--candidate"       = { opts$candidate <- v; i <- i + 2 },
      "--outdir"          = { opts$outdir <- v; i <- i + 2 },
      "--registries_root" = { opts$registries_root <- v; i <- i + 2 },
      { i <- i + 1 }
    )
  }
  opts
}

read_block <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(fromJSON(path, simplifyVector = TRUE),
           error = function(e) { log_msg("failed to parse ", path, ": ", e$message); NULL })
}

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  opts <- parse_args(argv)
  for (k in c("candidate", "outdir")) {
    if (is.na(opts[[k]])) stop("[sd_substrate_concordance] missing required --", k)
  }

  cand_dir <- file.path(opts$outdir, opts$candidate)
  str_dir  <- file.path(cand_dir, "structured")
  if (!dir.exists(str_dir)) {
    stop("[sd_substrate_concordance] structured dir not found (run angles first): ",
         str_dir)
  }

  mm2_json <- file.path(str_dir, "sd_substrate_minimap2.json")
  b2_json  <- file.path(str_dir, "sd_substrate_biser2.json")

  mm2_block <- read_block(mm2_json)
  b2_block  <- read_block(b2_json)

  if (is.null(mm2_block) && is.null(b2_block)) {
    stop("[sd_substrate_concordance] neither angle's block found — nothing to combine")
  }

  # Extract per-side verdicts
  mm2_left  <- if (!is.null(mm2_block)) mm2_block$data$q4a_mm2_mechanism_left  else NA_character_
  mm2_right <- if (!is.null(mm2_block)) mm2_block$data$q4a_mm2_mechanism_right else NA_character_

  # BISER2 is currently candidate-level (nahr_class is one value for the whole
  # straddle). We apply the same class to both sides — if a single SD pair
  # straddles both breakpoints, both sides get the same verdict.
  b2_class  <- if (!is.null(b2_block)) b2_block$data$q4b_biser2_nahr_class else NA_character_

  left  <- classify_concordance_side(mm2_left,  b2_class)
  right <- classify_concordance_side(mm2_right, b2_class)
  combined <- combine_sides(left, right)

  # Structured block
  block <- list(
    block_type    = "sd_substrate_concordance",
    candidate_id  = opts$candidate,
    source_script = "sd_substrate_concordance.R",
    data = list(
      q4_sd_concordance_left      = left$concordance,
      q4_sd_concordance_right     = right$concordance,
      q4_sd_concordance_overall   = combined$overall_conc,
      q4_sd_confidence_left       = left$confidence,
      q4_sd_confidence_right      = right$confidence,
      q4_sd_confidence_overall    = combined$overall_conf,
      q4_sd_primary_call_source   = combined$primary_source,
      q4_sd_primary_call          = combined$overall_call,
      # Backward-compat mirrors into mechanism.schema.json keys
      q4_sd_concordance           = combined$overall_conc,
      q4_mechanism_confidence     = combined$overall_conf,
      # Provenance — path pointers
      q4_sd_source_minimap2_block = if (!is.null(mm2_block)) normalizePath(mm2_json, mustWork = FALSE) else NA_character_,
      q4_sd_source_biser2_block   = if (!is.null(b2_block))  normalizePath(b2_json, mustWork = FALSE)  else NA_character_
    )
  )

  out_json <- file.path(str_dir, "sd_substrate_concordance.json")
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
          block_type    = "sd_substrate_concordance",
          data          = block$data,
          source_script = "sd_substrate_concordance.R"
        )
        log_msg("registered block for ", opts$candidate)
      }, error = function(e) {
        log_msg("registry write failed: ", e$message)
      })
    }
  }

  log_msg(opts$candidate,
          " — overall concordance=", combined$overall_conc,
          " confidence=", combined$overall_conf,
          " call=", combined$overall_call,
          " (", combined$primary_source, " primary)")
}

invoked_as_script <- any(grepl("--file=", commandArgs(trailingOnly = FALSE)))
if (invoked_as_script) main()
