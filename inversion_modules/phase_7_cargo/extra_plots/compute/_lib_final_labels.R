#!/usr/bin/env Rscript
# =============================================================================
# _lib_final_labels.R — Phase 4 v7 label loader for phase_7_cargo
# =============================================================================
#
# Helper for BREEDING_A_broodstock_compatibility.R,
# BREEDING_C_founder_haplotype_tracker.R,
# BREEDING_D_recombination_atlas.R
#
# Reads the per-candidate final_label.json files written by
# phase_4_postprocessing/4g_final_classification/assign_structural_class_v7.py
# and returns:
#   - a data.table of (candidate_id, q_overall_structural_class, weakest_component)
#   - a convenience filter function for selecting usable candidates
#
# BREEDING scripts should exclude candidates labelled
# `complex_rearrangement_out_of_scope` and `family_confounded_locus` from
# broodstock compatibility scoring, founder haplotype tracking, and
# recombination atlas analysis — those are unresolved or spurious and
# using them for breeding decisions would be misleading.
#
# USAGE (at the top of each BREEDING_* script, after the config block):
#
#   source(file.path(Sys.getenv("PHASE7_LIB"), "_lib_final_labels.R"))
#   # or: source("phase_7_cargo/extra_plots/compute/_lib_final_labels.R")
#
#   V7_FINAL_DIR <- Sys.getenv("V7_FINAL_DIR",
#     file.path(BASE, "phase4_v7_blocks", "final"))
#
#   if (dir.exists(V7_FINAL_DIR)) {
#     labels    <- load_final_labels(V7_FINAL_DIR)
#     usable_cids <- filter_usable_candidates(labels)
#     cands     <- cands[candidate_id %in% usable_cids]
#     cat(sprintf("[breeding] filtered to %d usable candidates by v7 labels\n",
#                 length(usable_cids)))
#   } else {
#     cat("[breeding] V7_FINAL_DIR not found; running on ALL candidates (no label filter)\n")
#   }
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# -----------------------------------------------------------------------------
# load_final_labels — read every final_label.json under v7_final_dir
# -----------------------------------------------------------------------------
load_final_labels <- function(v7_final_dir) {
  if (!dir.exists(v7_final_dir)) {
    warning("v7_final_dir not found: ", v7_final_dir)
    return(data.table(candidate_id = character(),
                      q_overall_structural_class = character(),
                      weakest_component = character()))
  }

  label_files <- list.files(v7_final_dir,
                             pattern = "final_label\\.json$",
                             recursive = TRUE,
                             full.names = TRUE)

  if (length(label_files) == 0) {
    return(data.table(candidate_id = character(),
                      q_overall_structural_class = character(),
                      weakest_component = character()))
  }

  rows <- lapply(label_files, function(p) {
    j <- tryCatch(jsonlite::fromJSON(p, simplifyVector = TRUE),
                  error = function(e) NULL)
    if (is.null(j)) return(NULL)
    data.table(
      candidate_id               = j$candidate_id              %||% NA_character_,
      q_overall_structural_class = j$q_overall_structural_class %||% NA_character_,
      weakest_component          = j$weakest_component          %||% NA_character_
    )
  })
  rbindlist(Filter(Negate(is.null), rows), fill = TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -----------------------------------------------------------------------------
# filter_usable_candidates — return candidate_ids usable for breeding logic
# -----------------------------------------------------------------------------
#
# Excludes out-of-scope and family-confounded candidates. The "usable"
# set is anything matching ^supported_(balanced_inversion|shared|nested)
# which covers:
#   - supported_balanced_inversion_simple
#   - supported_balanced_inversion_NAHR_like_hypothesis
#   - supported_balanced_inversion_NAHR_like_supported_by_assembly
#   - supported_balanced_inversion_NHEJ_like_hypothesis
#   - supported_balanced_inversion_NHEJ_like_supported_by_assembly
#   - supported_balanced_inversion_with_substrate_mechanism_unresolved
#   - supported_balanced_inversion_with_edge_recombinants
#   - supported_shared_between_species_inversion
#   - supported_nested_inversion
#
filter_usable_candidates <- function(labels_dt,
                                       pattern = "^supported_(balanced_inversion|shared|nested)") {
  if (nrow(labels_dt) == 0) return(character())
  keep <- grepl(pattern, labels_dt$q_overall_structural_class)
  labels_dt$candidate_id[keep]
}

# -----------------------------------------------------------------------------
# label_summary — print a short histogram of labels
# -----------------------------------------------------------------------------
label_summary <- function(labels_dt) {
  if (nrow(labels_dt) == 0) {
    cat("[final_labels] no labels loaded\n")
    return(invisible(NULL))
  }
  tbl <- sort(table(labels_dt$q_overall_structural_class), decreasing = TRUE)
  cat(sprintf("[final_labels] %d candidates, %d distinct labels:\n",
              nrow(labels_dt), length(tbl)))
  for (i in seq_along(tbl)) {
    cat(sprintf("  %5d  %s\n", tbl[i], names(tbl)[i]))
  }
  invisible(tbl)
}
