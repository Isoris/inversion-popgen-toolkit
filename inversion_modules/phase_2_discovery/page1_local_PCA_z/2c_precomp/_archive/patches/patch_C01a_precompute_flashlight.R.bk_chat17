#!/usr/bin/env Rscript
# =============================================================================
# PATCH: Flashlight annotation for STEP_C01a_snake1_precompute.R
#
# INSERT this block AFTER the inv_like_dt is computed (after the
# compute_inv_likeness_all() call) and BEFORE the dosage het merge.
#
# What it does:
#   For each window in inv_like_dt, checks whether a DELLY/Manta INV call
#   overlaps that window. If yes, annotates with:
#     sv_inv_overlap      — 0/1 flag
#     sv_inv_confidence   — confidence level of the best overlapping INV
#     sv_inv_af           — allele frequency of the overlapping INV
#     sv_het_del_count    — number of large het-DELs in this window
#     sv_n_anchors        — number of anchor samples (MEDIUM+ confidence)
#
# These columns are purely ANNOTATION — they do NOT alter inv_likeness,
# structure_likeness, or family_likeness. Downstream scripts (merge,
# scoring, decomposition) consume them for informed-prior decisions.
#
# REQUIRES: flashlight_loader.R sourced (provides load_flashlight, etc.)
# =============================================================================

# ─── Source flashlight loader ────────────────────────────────────────
fl_loader_path <- Sys.getenv("FLASHLIGHT_LOADER", "")
if (!nzchar(fl_loader_path)) {
  # Auto-detect from codebase
  candidates <- c(
    file.path(dirname(outdir), "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path),
    file.path(dirname(dirname(outdir)), "inversion_codebase_v8.5.9", "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path),
    file.path(Sys.getenv("BASE", "."), "inversion_localpca_v7", "inversion_codebase_v8.5.9", "utils", "flashlight_loader_v2.R")
    if (!file.exists(fl_loader_path)) fl_loader_path <- sub("_v2", "", fl_loader_path)
  )
  for (cand in candidates) {
    if (file.exists(cand)) { fl_loader_path <- cand; break }
  }
}

has_flashlight <- FALSE
if (nzchar(fl_loader_path) && file.exists(fl_loader_path)) {
  tryCatch({
    source(fl_loader_path)
    has_flashlight <- TRUE
    message("[PRECOMP] Flashlight loader sourced: ", fl_loader_path)
  }, error = function(e) {
    message("[PRECOMP] Flashlight loader failed: ", e$message, " — continuing without")
  })
} else {
  message("[PRECOMP] No flashlight_loader.R found — SV annotation columns will be NA")
}

# ─── Annotate inv_like_dt with flashlight columns ────────────────────
if (has_flashlight && nrow(inv_like_dt) > 0) {
  message("[PRECOMP] ── FLASHLIGHT ANNOTATION ──")

  # Pre-initialize columns
  inv_like_dt[, `:=`(
    sv_inv_overlap    = 0L,
    sv_inv_confidence = NA_character_,
    sv_inv_af         = NA_real_,
    sv_het_del_count  = 0L,
    sv_n_anchors      = 0L
  )]

  n_annotated <- 0L
  for (chr in unique(inv_like_dt$chrom)) {
    fl <- load_flashlight(chr)
    if (is.null(fl) || nrow(fl$inv_calls) == 0) next

    chr_idx <- which(inv_like_dt$chrom == chr)
    if (length(chr_idx) == 0) next

    for (wi in chr_idx) {
      w_start <- inv_like_dt$start_bp[wi]
      w_end   <- inv_like_dt$end_bp[wi]

      # Cheat 1: Does an SV INV call overlap this window?
      overlapping <- fl$inv_calls[bp1 <= w_end & bp2 >= w_start]
      if (nrow(overlapping) > 0) {
        # Take the highest-confidence overlapping INV
        conf_order <- c("VERY_HIGH", "HIGH", "MEDIUM", "LOW", "MINIMAL", "UNKNOWN")
        best_idx <- 1L
        if ("confidence_level" %in% names(overlapping)) {
          ranks <- match(overlapping$confidence_level, conf_order)
          ranks[is.na(ranks)] <- length(conf_order)
          best_idx <- which.min(ranks)
        }

        inv_like_dt[wi, `:=`(
          sv_inv_overlap    = 1L,
          sv_inv_confidence = overlapping$confidence_level[best_idx],
          sv_inv_af         = overlapping$af[best_idx]
        )]

        # Count anchor samples for this region
        anchors <- fl$sample_inv_states[
          inv_id %in% overlapping$inv_id &
          sv_genotype != "MISSING" &
          sv_confidence %in% c("HIGH", "MEDIUM")
        ]
        inv_like_dt[wi, sv_n_anchors := nrow(anchors)]
        n_annotated <- n_annotated + 1L
      }

      # Cheat 2: Het-DELs overlapping this window
      if (nrow(fl$breakpoint_dels) > 0) {
        bp_dels_here <- fl$breakpoint_dels[
          bp_pos >= w_start & bp_pos <= w_end
        ]
        if (nrow(bp_dels_here) > 0) {
          inv_like_dt[wi, sv_het_del_count := nrow(bp_dels_here)]
        }
      }

      # Also count internal het-DELs in this window
      if (nrow(fl$internal_dels) > 0) {
        int_dels_here <- fl$internal_dels[
          del_start >= w_start & del_end <= w_end
        ]
        if (nrow(int_dels_here) > 0) {
          inv_like_dt[wi, sv_het_del_count := sv_het_del_count + nrow(int_dels_here)]
        }
      }
    }
  }

  n_with_sv <- sum(inv_like_dt$sv_inv_overlap == 1, na.rm = TRUE)
  n_with_del <- sum(inv_like_dt$sv_het_del_count > 0, na.rm = TRUE)
  message("[PRECOMP] Flashlight annotation complete:")
  message("  Windows with SV INV overlap: ", n_with_sv, " / ", nrow(inv_like_dt))
  message("  Windows with het-DEL: ", n_with_del, " / ", nrow(inv_like_dt))
  message("  Anchor assignments: ", sum(inv_like_dt$sv_n_anchors, na.rm = TRUE))
} else {
  # Add empty columns so downstream scripts don't break
  if (nrow(inv_like_dt) > 0) {
    inv_like_dt[, `:=`(
      sv_inv_overlap    = 0L,
      sv_inv_confidence = NA_character_,
      sv_inv_af         = NA_real_,
      sv_het_del_count  = 0L,
      sv_n_anchors      = 0L
    )]
  }
}

# NOTE: These columns propagate into per-chromosome precomp RDS files
# via the existing merge(dt, inv_like_dt) block in precompute_one_chr().
# No additional changes needed — the merge is column-agnostic.
