#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01g_boundary_catalog.R  (v9.3.4 — Engine B Fst for Cheat 4)
#
# UNIFIED BOUNDARY CATALOG — collects boundaries from ALL pipeline sources,
# merges by proximity, runs cheats 4/10/11 on each, produces concordance.
#
# FIVE BOUNDARY SOURCES:
#   1. PHASE_01C block_detect → boundary_catalog_<chr>.tsv.gz
#      (clear_hard / clear_soft / diffuse / inner_hard_* / inner_soft_*)
#      Also: blue_cross_verdicts (assembly errors, internal transitions)
#      Also: block_concordance (which blocks share sample grouping)
#
#   2. phase_2/2d staircase → scoring_table_<chr>.tsv
#      (block start/end from staircase vote peaks)
#
#   3. phase_2/2c C01b_1 seeded regions → seeded_regions_summary_<chr>.tsv.gz
#      (region start/end = where the seed grew to)
#      Legacy filename snake1_core_regions_<chr>.tsv.gz is accepted as
#      fallback for pre-rename output dirs.
#
#   4. phase_2/2c C00 SV prior → sv_prior_<chr>.rds
#      (DELLY/Manta INV bp1/bp2 at base-pair precision)
#      Legacy filename sv_flashlight_<chr>.rds is accepted as fallback.
#
#   5. Blue-cross inner boundaries from 01C (already in source 1)
#      (assembly errors, recombination zones, nested breakpoints)
#
# THREE CHEATS ON EACH BOUNDARY:
#   Cheat 4:  Fst/Hobs step function (requires Engine B or fallback)
#   Cheat 10: Read depth anomaly at boundary (requires BAM/mosdepth)
#   Cheat 11: Clipped read pileup (requires BAM)
#
# OUTPUT:
#   boundary_catalog_unified.tsv.gz  — one row per boundary, all sources merged
#   boundary_concordance.tsv.gz      — per-boundary cheat concordance + verdict
#   boundary_summary.tsv             — counts by type and concordance level
#
# Usage:
#   Rscript STEP_C01g_boundary_catalog.R \
#     --precomp <precomp_dir> \
#     --landscape <01C_output_dir> \
#     --staircase <inv_detect_output_dir> \
#     --cores <snake1_output_dir> \
#     --sv_prior <sv_prior_dir> \
#     --outdir <output_dir> \
#     [--chrom C_gar_LG01]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ── Source load_bridge.R (provides smap, reg, get_Q, get_region_stats) ───────
.bridge_file <- Sys.getenv("LOAD_BRIDGE", "")
if (!nzchar(.bridge_file)) {
  for (.bp in c("utils/load_bridge.R", "../utils/load_bridge.R",
                 file.path(Sys.getenv("BASE", ""), "inversion_codebase_v8.5/utils/load_bridge.R"))) {
    if (file.exists(.bp)) { .bridge_file <- .bp; break }
  }
}
.bridge_available <- FALSE
if (nzchar(.bridge_file) && file.exists(.bridge_file)) {
  tryCatch({
    source(.bridge_file)
    .bridge_available <- TRUE
    message("[C01g] load_bridge.R sourced — Engine B available for Cheat 4")
  }, error = function(e) {
    message("[C01g] WARNING: load_bridge.R failed: ", conditionMessage(e))
  })
} else {
  message("[C01g] load_bridge.R not found — Cheat 4 will use PCA proxy fallback")
}

args <- commandArgs(trailingOnly = TRUE)
precomp_dir <- NULL; landscape_dir <- NULL; staircase_dir <- NULL
cores_dir <- NULL; sv_prior_dir <- NULL; outdir <- "boundary_catalog"
chrom_filter <- NULL; samples_file <- NULL
bam_manifest <- NULL; mosdepth_dir <- NULL
repeatmasker_file <- NULL

# BUGFIX 2026-04-17: --ref_fasta was parsed but never read anywhere.
# ARCHIVED 2026-04-17: --scores (Cheat 17 fossil detection) removed.
# Fossil detection created a circular dependency with C01d (C01g reads
# C01d's candidate_scores.tsv.gz, C01d reads C01g's boundary_catalog),
# requiring a two-pass run. The fossil-vs-active distinction was a
# diagnostic annotation not in the manuscript main path, so removing it
# is simpler than orchestrating the two passes. Both flags are accepted
# (silently ignored) for back-compat with existing launchers.
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--precomp" && i < length(args))     { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--landscape" && i < length(args))  { landscape_dir <- args[i+1]; i <- i+2 }
  else if (a == "--staircase" && i < length(args))  { staircase_dir <- args[i+1]; i <- i+2 }
  else if (a == "--cores" && i < length(args))      { cores_dir <- args[i+1]; i <- i+2 }
  else if (a == "--sv_prior" && i < length(args)) { sv_prior_dir <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))     { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--chrom" && i < length(args))      { chrom_filter <- args[i+1]; i <- i+2 }
  else if (a == "--samples" && i < length(args))    { samples_file <- args[i+1]; i <- i+2 }
  else if (a == "--bam_manifest" && i < length(args)) { bam_manifest <- args[i+1]; i <- i+2 }
  else if (a == "--mosdepth_dir" && i < length(args)) { mosdepth_dir <- args[i+1]; i <- i+2 }
  else if (a == "--ref_fasta" && i < length(args))    { i <- i+2 }  # accepted but unused
  else if (a == "--repeatmasker" && i < length(args)) { repeatmasker_file <- args[i+1]; i <- i+2 }
  else if (a == "--scores" && i < length(args))       { i <- i+2 }  # accepted but unused (fossil detection archived)
  else { i <- i+1 }
}

if (is.null(precomp_dir)) stop("--precomp required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Source inversion config for paths ──
# The config provides BASE, SAMPLES_IND, FLASHLIGHT_DIR, LANDSCAPE_DIR, etc.
# R reads shell config via Sys.getenv or by parsing the file.
BASE <- Sys.getenv("BASE", "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04")

# Markdup BAMs and mosdepth from DELLY pipeline
MARKDUP_DIR <- file.path(BASE, "delly_sv/00_markdup")
SAMPLES_IND <- Sys.getenv("SAMPLES_IND",
  file.path(BASE, "het_roh/01_inputs_check/samples.ind"))

# Auto-detect if not provided via args
if (is.null(mosdepth_dir) && dir.exists(MARKDUP_DIR)) mosdepth_dir <- MARKDUP_DIR
if (is.null(sv_prior_dir)) {
  fl_auto <- file.path(BASE, "inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/sv_prior")
  if (dir.exists(fl_auto)) sv_prior_dir <- fl_auto
}
if (is.null(landscape_dir)) {
  ls_auto <- file.path(BASE, "inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/landscape")
  if (dir.exists(ls_auto)) landscape_dir <- ls_auto
}
if (is.null(samples_file) && file.exists(SAMPLES_IND)) samples_file <- SAMPLES_IND

message("═══════════════════════════════════════════════════════════")
message("[C01g] Unified Boundary Catalog v9.3.2")
message("═══════════════════════════════════════════════════════════")
message("  Precomp:    ", precomp_dir %||% "NONE")
message("  Landscape:  ", landscape_dir %||% "NONE")
message("  Staircase:  ", staircase_dir %||% "NONE")
message("  Cores:      ", cores_dir %||% "NONE")
message("  Flashlight: ", sv_prior_dir %||% "NONE")
message("  Output:     ", outdir)

# =============================================================================
# LOAD PRECOMP
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

precomp_list <- list()
for (f in rds_files) {
  obj <- readRDS(f)
  precomp_list[[obj$chrom]] <- obj
}

chroms <- names(precomp_list)
if (!is.null(chrom_filter)) chroms <- intersect(chroms, chrom_filter)
message("[C01g] Chromosomes: ", length(chroms))

# Sample names from precomp
sample_names <- NULL
for (chr_tmp in chroms) {
  pc1_cols <- grep("^PC_1_", names(precomp_list[[chr_tmp]]$dt), value = TRUE)
  if (length(pc1_cols) > 0) { sample_names <- sub("^PC_1_", "", pc1_cols); break }
}
n_samples <- length(sample_names)
message("[C01g] Samples: ", n_samples)

# Optionally load real sample names for mapping
real_names <- NULL
if (!is.null(samples_file) && file.exists(samples_file)) {
  real_names <- trimws(readLines(samples_file))
  real_names <- real_names[nzchar(real_names)]
  message("[C01g] Real sample names loaded: ", length(real_names))
}

# =============================================================================
# MERGE PROXIMITY THRESHOLD
# =============================================================================

MERGE_DIST_BP <- 50000L  # boundaries within 50 kb are merged

# Type priority for dedup (lower = keep)
TYPE_PRIORITY <- c(
  STRUCTURAL_SHARP = 1, OUTER_HARD = 2, STRUCTURAL_MODERATE = 3,
  INNER_HARD = 4, SYSTEM_JUNCTION = 5, HOTCORE_EDGE = 6,
  STRUCTURAL_DIFFUSE = 7, OUTER_SOFT = 8, INTERNAL_GAP = 9,
  INTERNAL_TRANSITION = 10, INNER_SOFT = 11,
  BACKGROUND_FADE = 12, COMPOSITE = 13
)

# =============================================================================
# COLLECT BOUNDARIES FROM ALL SOURCES
# =============================================================================

all_bounds <- list()
bid <- 0L

for (chr in chroms) {
  pc <- precomp_list[[chr]]
  dt <- pc$dt
  sim_mat <- pc$sim_mat
  n_w <- nrow(dt)
  if (n_w < 50) next

  message("\n[C01g] ═══ ", chr, " ═══")

  # ── SOURCE 1: PHASE_01C landscape boundaries ───────────────────────
  n_01c <- 0L
  if (!is.null(landscape_dir)) {
    f <- file.path(landscape_dir, paste0("boundary_catalog_", chr, ".tsv.gz"))
    if (file.exists(f)) {
      lnd <- fread(f)
      for (li in seq_len(nrow(lnd))) {
        lb <- lnd[li]
        bp <- lb$pos_bp %||% lb$boundary_bp
        if (is.null(bp) || !is.finite(bp)) next

        # Map 01C types to unified taxonomy
        lt <- lb$boundary_type %||% "unknown"
        mapped <- if (grepl("clear_hard", lt)) "STRUCTURAL_SHARP"
          else if (grepl("inner_hard_different", lt)) "SYSTEM_JUNCTION"
          else if (grepl("inner_hard_same", lt)) "HOTCORE_EDGE"
          else if (grepl("inner_hard_assembly|inner_hard_suspect", lt)) "INNER_HARD"
          else if (grepl("inner_soft", lt)) "INTERNAL_TRANSITION"
          else if (grepl("clear_soft", lt)) "STRUCTURAL_MODERATE"
          else if (grepl("diffuse", lt)) "STRUCTURAL_DIFFUSE"
          else if (grepl("chromosome_edge", lt)) "BACKGROUND_FADE"
          else "OUTER_SOFT"

        bid <- bid + 1L
        all_bounds[[bid]] <- data.table(
          boundary_id = bid, chrom = chr, boundary_bp = as.integer(bp),
          boundary_type = mapped,
          source = "01C_landscape", side = lb$side %||% "classified",
          drop_magnitude = lb$drop_magnitude %||% NA_real_,
          transition_width = lb$transition_width %||% NA_integer_,
          block_left = lb$block %||% NA_character_,
          block_right = NA_character_,
          candidate_id = NA_character_
        )
        n_01c <- n_01c + 1L
      }
    }

    # Also load blue-cross inner boundaries
    f_bx <- file.path(landscape_dir, paste0("blue_cross_verdicts_", chr, ".tsv.gz"))
    if (file.exists(f_bx)) {
      bx <- fread(f_bx)
      for (bi in seq_len(nrow(bx))) {
        b <- bx[bi]
        bp <- as.integer((b$start_bp + b$end_bp) / 2)
        it <- b$inner_type %||% "inner_ambiguous"

        mapped <- if (grepl("same_system", it)) "HOTCORE_EDGE"
          else if (grepl("different_system", it)) "SYSTEM_JUNCTION"
          else if (grepl("assembly|suspect", it)) "INNER_HARD"
          else if (grepl("soft|ambiguous", it)) "INTERNAL_TRANSITION"
          else "INTERNAL_GAP"

        bid <- bid + 1L
        all_bounds[[bid]] <- data.table(
          boundary_id = bid, chrom = chr, boundary_bp = bp,
          boundary_type = mapped,
          source = "01C_blue_cross", side = "inner",
          drop_magnitude = b$drop %||% NA_real_,
          transition_width = b$n_cross_windows %||% NA_integer_,
          block_left = b$block %||% NA_character_,
          block_right = NA_character_,
          candidate_id = NA_character_
        )
        n_01c <- n_01c + 1L
      }
    }
    message("  01C landscape: ", n_01c, " boundaries")
  }

  # ── SOURCE 2: Staircase blocks (left/right edges) ──────────────────
  n_stair <- 0L
  if (!is.null(staircase_dir)) {
    f <- file.path(staircase_dir, paste0("scoring_table_", chr, ".tsv"))
    if (!file.exists(f)) f <- file.path(staircase_dir, paste0("scoring_table_", chr, ".tsv.gz"))
    if (file.exists(f)) {
      st <- fread(f)
      # Staircase blocks have start_bin / end_bin or start_bp / end_bp
      bp_start_col <- intersect(c("start_bp", "block_start_bp"), names(st))
      bp_end_col <- intersect(c("end_bp", "block_end_bp"), names(st))

      if (length(bp_start_col) > 0 && length(bp_end_col) > 0) {
        for (si in seq_len(nrow(st))) {
          s <- st[si]
          s_bp <- s[[bp_start_col[1]]]
          e_bp <- s[[bp_end_col[1]]]
          cid <- s$block_id %||% s$candidate_id %||% paste0("stair_", si)

          # Left boundary
          bid <- bid + 1L
          all_bounds[[bid]] <- data.table(
            boundary_id = bid, chrom = chr, boundary_bp = as.integer(s_bp),
            boundary_type = "OUTER_HARD",
            source = "staircase", side = "left",
            drop_magnitude = NA_real_, transition_width = NA_integer_,
            block_left = "background", block_right = as.character(cid),
            candidate_id = as.character(cid)
          )

          # Right boundary
          bid <- bid + 1L
          all_bounds[[bid]] <- data.table(
            boundary_id = bid, chrom = chr, boundary_bp = as.integer(e_bp),
            boundary_type = "OUTER_HARD",
            source = "staircase", side = "right",
            drop_magnitude = NA_real_, transition_width = NA_integer_,
            block_left = as.character(cid), block_right = "background",
            candidate_id = as.character(cid)
          )
          n_stair <- n_stair + 2L
        }
      }
    }
    message("  Staircase: ", n_stair, " boundaries")
  }

  # ── SOURCE 3: Seeded-region edges (formerly snake cores) ───────────
  n_snake <- 0L
  if (!is.null(cores_dir)) {
    # New filename from STEP_C01b_1_seeded_regions.R
    f <- file.path(cores_dir, paste0("seeded_regions_summary_", chr, ".tsv.gz"))
    if (!file.exists(f)) {
      # Fallback for pre-rename output directories
      f <- file.path(cores_dir, paste0("snake1_core_regions_", chr, ".tsv.gz"))
    }
    if (file.exists(f)) {
      cores <- fread(f)
      # Normalise column names (new input uses region_id; legacy used region_id)
      if ("region_id" %in% names(cores) && !"region_id" %in% names(cores)) {
        cores[, region_id := region_id]
      }
      for (ci in seq_len(nrow(cores))) {
        cr <- cores[ci]
        cid <- cr$region_id %||% paste0("region_", ci)

        # Left edge
        bid <- bid + 1L
        all_bounds[[bid]] <- data.table(
          boundary_id = bid, chrom = chr, boundary_bp = as.integer(cr$start_bp),
          boundary_type = "OUTER_SOFT",  # seeded-region edges are softer than staircase
          source = "seeded_region", side = "left",
          drop_magnitude = NA_real_, transition_width = NA_integer_,
          block_left = "background", block_right = as.character(cid),
          candidate_id = as.character(cid)
        )

        # Right edge
        bid <- bid + 1L
        all_bounds[[bid]] <- data.table(
          boundary_id = bid, chrom = chr, boundary_bp = as.integer(cr$end_bp),
          boundary_type = "OUTER_SOFT",
          source = "seeded_region", side = "right",
          drop_magnitude = NA_real_, transition_width = NA_integer_,
          block_left = as.character(cid), block_right = "background",
          candidate_id = as.character(cid)
        )
        n_snake <- n_snake + 2L
      }
    }
    message("  Seeded regions: ", n_snake, " boundaries")
  }

  # ── SOURCE 4: SV prior breakpoints (formerly sv_prior) ───────────
  # The CLI flag is still --sv_prior for back-compat with launchers;
  # the underlying dir is what STEP_C00_build_sv_prior.R writes.
  n_sv <- 0L
  if (!is.null(sv_prior_dir)) {
    # New filename from STEP_C00_build_sv_prior.R
    fl_file <- file.path(sv_prior_dir, paste0("sv_prior_", chr, ".rds"))
    if (!file.exists(fl_file)) {
      # Fallback for pre-rename output directories
      fl_file <- file.path(sv_prior_dir, paste0("sv_flashlight_", chr, ".rds"))
    }
    if (file.exists(fl_file)) {
      fl <- tryCatch(readRDS(fl_file), error = function(e) NULL)
      if (!is.null(fl) && "inv_calls" %in% names(fl) && nrow(fl$inv_calls) > 0) {
        for (ii in seq_len(nrow(fl$inv_calls))) {
          inv <- fl$inv_calls[ii]
          for (bp_pos in c(inv$bp1, inv$bp2)) {
            if (!is.finite(bp_pos)) next

            conf <- inv$confidence_level %||% "MEDIUM"
            bid <- bid + 1L
            all_bounds[[bid]] <- data.table(
              boundary_id = bid, chrom = chr, boundary_bp = as.integer(bp_pos),
              boundary_type = if (conf %in% c("HIGH", "VERY_HIGH"))
                "STRUCTURAL_SHARP" else "STRUCTURAL_MODERATE",
              source = paste0("sv_", inv$caller %||% "unknown"),
              side = if (bp_pos == inv$bp1) "left" else "right",
              drop_magnitude = NA_real_, transition_width = NA_integer_,
              block_left = NA_character_, block_right = NA_character_,
              candidate_id = inv$inv_id %||% NA_character_
            )
            n_sv <- n_sv + 1L
          }
        }
      }

      # Cheat 8: BND-triangulated inversions (additional SV breakpoints)
      if (!is.null(fl) && "bnd_triangulated" %in% names(fl) &&
          nrow(fl$bnd_triangulated) > 0) {
        bnd_tri <- fl$bnd_triangulated
        for (bi in seq_len(nrow(bnd_tri))) {
          b <- bnd_tri[bi]
          for (bp_pos in c(b$bp1, b$bp2)) {
            if (!is.finite(bp_pos)) next

            conf <- if (b$carrier_jaccard > 0.6) "HIGH" else "MEDIUM"
            bid <- bid + 1L
            all_bounds[[bid]] <- data.table(
              boundary_id = bid, chrom = chr, boundary_bp = as.integer(bp_pos),
              boundary_type = if (conf == "HIGH")
                "STRUCTURAL_SHARP" else "STRUCTURAL_MODERATE",
              source = "sv_bnd_triangulation",
              side = if (bp_pos == b$bp1) "left" else "right",
              drop_magnitude = NA_real_, transition_width = NA_integer_,
              block_left = NA_character_, block_right = NA_character_,
              candidate_id = b$bnd_left_id %||% NA_character_
            )
            n_sv <- n_sv + 1L
          }
        }
        message("    + Cheat 8 BND: ", nrow(bnd_tri) * 2, " breakpoints")
      }
    }
    message("  Flashlight SV: ", n_sv, " boundaries")
  }
}

# =============================================================================
# COMBINE AND DEDUP
# =============================================================================

if (length(all_bounds) == 0) {
  message("[C01g] No boundaries found from any source. Exiting.")
  quit(save = "no", status = 0)
}

catalog <- rbindlist(all_bounds, fill = TRUE)
message("\n[C01g] Raw boundaries: ", nrow(catalog))

# Sort by chromosome and position
catalog <- catalog[order(chrom, boundary_bp)]

# Assign type priority
catalog[, type_rank := TYPE_PRIORITY[boundary_type]]
catalog[is.na(type_rank), type_rank := 99L]

# Group boundaries within MERGE_DIST_BP clusters
catalog[, cluster := cumsum(c(TRUE,
  diff(boundary_bp) > MERGE_DIST_BP | chrom != shift(chrom, 1, "")))]

# Within each cluster, keep highest-priority but collect all sources
deduped <- catalog[, {
  best_idx <- which.min(type_rank)
  sources_all <- unique(source)
  cands_all <- unique(na.omit(candidate_id))
  # BUGFIX 2026-04-17 (chat 7, FIX 35): produce matched_core_id as the
  # first non-NA candidate_id in the cluster. The registry-write block
  # at L1407 guards on `"matched_core_id" %in% names(confirmed)`, which
  # was always FALSE because no prior step set this column. When the
  # registry helper lib lands (flagged in chat 7 Finding 8/9), the
  # guard will now succeed and the per-boundary blocks will write.
  # Note: candidate_ids (semicolon-joined, all ids in cluster) is kept
  # for traceability; matched_core_id is the canonical single id for
  # the boundary's candidate-association (picked from best_idx row when
  # available, otherwise the first cluster candidate).
  best_cid <- if (!is.na(candidate_id[best_idx]) && nzchar(candidate_id[best_idx])) {
    candidate_id[best_idx]
  } else if (length(cands_all) > 0) {
    cands_all[1]
  } else {
    NA_character_
  }
  .(boundary_bp = boundary_bp[best_idx],
    boundary_type = boundary_type[best_idx],
    side = side[best_idx],
    drop_magnitude = max(drop_magnitude, na.rm = TRUE),
    transition_width = min(transition_width, na.rm = TRUE),
    block_left = block_left[best_idx],
    block_right = block_right[best_idx],
    sources = paste(sources_all, collapse = ";"),
    candidate_ids = paste(cands_all, collapse = ";"),
    matched_core_id = as.character(best_cid),
    n_sources = length(sources_all),
    n_within_cluster = .N
  )
}, by = .(chrom, cluster)]
deduped[, cluster := NULL]
deduped[, boundary_id := .I]

# Fix infinite values from max/min on empty
deduped[is.infinite(drop_magnitude), drop_magnitude := NA_real_]
deduped[is.infinite(transition_width), transition_width := NA_integer_]

message("[C01g] After dedup (", MERGE_DIST_BP / 1000, "kb): ", nrow(deduped), " boundaries")
message("  Sources: ", paste(names(table(catalog$source)), table(catalog$source),
                              sep = "=", collapse = ", "))
message("  Types: ", paste(names(table(deduped$boundary_type)),
                            table(deduped$boundary_type), sep = "=", collapse = ", "))

# =============================================================================
# REFINE CLASSIFICATION FROM SIM_MAT
# =============================================================================
# For each boundary, compute sim_mat hardness and color-based type

message("[C01g] Refining boundary classification from sim_mat...")

deduped[, hardness := NA_real_]
deduped[, left_color := NA_character_]
deduped[, right_color := NA_character_]
deduped[, cross_color := NA_character_]

for (chr in unique(deduped$chrom)) {
  pc <- precomp_list[[chr]]
  if (is.null(pc) || is.null(pc$sim_mat)) next
  sim_mat <- pc$sim_mat; dt_chr <- pc$dt
  n_w <- nrow(dt_chr)
  mid_bps <- (dt_chr$start_bp + dt_chr$end_bp) / 2
  bg_level <- median(sim_mat[sim_mat > 0 & sim_mat < 1], na.rm = TRUE)
  if (!is.finite(bg_level)) bg_level <- 0.5

  classify_level <- function(sim_val) {
    if (!is.finite(sim_val)) return("unknown")
    if (sim_val < bg_level * 0.7) "blue"
    else if (sim_val < bg_level * 1.1) "white"
    else if (sim_val < bg_level * 1.5) "orange"
    else "red"
  }

  chr_idx <- which(deduped$chrom == chr)
  for (bi in chr_idx) {
    bp <- deduped$boundary_bp[bi]
    bnd_w <- which.min(abs(mid_bps - bp))
    if (bnd_w <= 5 || bnd_w >= n_w - 4) next

    left_idx  <- max(1, bnd_w - 5):(bnd_w - 1)
    right_idx <- (bnd_w + 1):min(n_w, bnd_w + 5)

    left_sim  <- mean(sim_mat[left_idx, left_idx], na.rm = TRUE)
    right_sim <- mean(sim_mat[right_idx, right_idx], na.rm = TRUE)
    cross_sim <- mean(sim_mat[left_idx, right_idx], na.rm = TRUE)

    within_mean <- mean(c(left_sim, right_sim), na.rm = TRUE)
    hard <- within_mean - cross_sim
    if (!is.finite(hard)) next

    deduped[bi, hardness := round(hard, 4)]
    deduped[bi, left_color := classify_level(left_sim)]
    deduped[bi, right_color := classify_level(right_sim)]
    deduped[bi, cross_color := classify_level(cross_sim)]

    # Refine type based on colors (only if not already SV-confirmed)
    if (!grepl("sv_", deduped$sources[bi])) {
      lc <- classify_level(left_sim); rc <- classify_level(right_sim)
      cc <- classify_level(cross_sim)
      if (hard > 0.15) {
        if ((lc == "blue" && rc == "red") || (lc == "red" && rc == "blue"))
          deduped[bi, boundary_type := "STRUCTURAL_SHARP"]
        else if ((lc == "blue" && rc == "orange") || (lc == "orange" && rc == "blue"))
          deduped[bi, boundary_type := "STRUCTURAL_MODERATE"]
        else if ((lc == "orange" && rc == "red") || (lc == "red" && rc == "orange"))
          deduped[bi, boundary_type := "HOTCORE_EDGE"]
        else if (lc %in% c("red", "orange") && rc %in% c("red", "orange") && cc == "blue")
          deduped[bi, boundary_type := "SYSTEM_JUNCTION"]
      } else if (hard > 0.05) {
        if ((lc == "blue" && rc %in% c("orange", "red")) ||
            (rc == "blue" && lc %in% c("orange", "red")))
          deduped[bi, boundary_type := "STRUCTURAL_DIFFUSE"]
        else if (lc %in% c("orange", "red") && rc %in% c("orange", "red") && cc == "blue")
          deduped[bi, boundary_type := "INTERNAL_GAP"]
      } else if (hard < 0.02) {
        if (lc == "blue" && rc == "blue")
          deduped[bi, boundary_type := "BACKGROUND_FADE"]
      }
    }
  }
}

message("[C01g] Refined types: ", paste(names(table(deduped$boundary_type)),
  table(deduped$boundary_type), sep = "=", collapse = ", "))

# =============================================================================
# CHEAT 4: BOUNDARY SHARPNESS (Fst/Hobs step)
# =============================================================================
# Compute Fst between band1 vs band3 on LEFT side vs RIGHT side.
# Sharp Fst step = real structural boundary. Diffuse = family LD.
#
# ENGINE B MODE (preferred): calls get_region_stats(what=c("Fst","Hobs"))
#   for left and right flanks with band1 vs band3 groups in CGA space.
#   This gives real Hudson Fst from BEAGLE dosage.
#
# FALLBACK MODE: uses PC1 between-band variance / total variance as proxy.

message("\n[C01g] Running Cheat 4 (boundary sharpness)...")

HALF_FLANK <- 20L
deduped[, cheat4_sharpness := NA_real_]
deduped[, cheat4_class := NA_character_]
deduped[, cheat4_fst_step := NA_real_]
deduped[, cheat4_fst_left := NA_real_]
deduped[, cheat4_fst_right := NA_real_]
deduped[, cheat4_hobs_het := NA_real_]
deduped[, cheat4_method := NA_character_]

# Helper: extract single Fst from dispatcher result
extract_fst <- function(s) {
  if (is.null(s$Fst) || length(s$Fst) == 0) return(NA_real_)
  as.numeric(s$Fst[[1]])
}

# Helper: extract Hobs from dispatcher result
extract_hobs <- function(s, subset_name = NULL) {
  if (is.null(s$Hobs) || length(s$Hobs) == 0) return(NA_real_)
  if (!is.null(subset_name) && subset_name %in% names(s$Hobs))
    return(s$Hobs[[subset_name]]$mean_Hobs)
  # Take first available subset
  as.numeric(s$Hobs[[1]]$mean_Hobs)
}

# PCA proxy fallback function (unchanged from v9.3.2)
fst_proxy <- function(mat, b1_idx, b3_idx) {
  avg <- colMeans(mat, na.rm = TRUE)
  b1_mean <- mean(avg[b1_idx], na.rm = TRUE)
  b3_mean <- mean(avg[b3_idx], na.rm = TRUE)
  total_var <- var(avg, na.rm = TRUE)
  if (!is.finite(total_var) || total_var == 0) return(NA_real_)
  between <- ((b1_mean - mean(avg))^2 * length(b1_idx) +
               (b3_mean - mean(avg))^2 * length(b3_idx)) /
              (length(b1_idx) + length(b3_idx))
  pmin(1, between / total_var)
}

for (chr in unique(deduped$chrom)) {
  pc <- precomp_list[[chr]]
  if (is.null(pc)) next
  dt_chr <- pc$dt; sim_mat <- pc$sim_mat
  n_w <- nrow(dt_chr)
  mid_bps <- (dt_chr$start_bp + dt_chr$end_bp) / 2
  pc1_cols <- intersect(paste0("PC_1_", sample_names), names(dt_chr))
  if (length(pc1_cols) < 20) next

  chr_idx <- which(deduped$chrom == chr)
  n_engine_b <- 0L; n_proxy <- 0L

  for (bi in chr_idx) {
    bp <- deduped$boundary_bp[bi]
    bnd_w <- which.min(abs(mid_bps - bp))
    if (bnd_w <= HALF_FLANK || bnd_w >= n_w - HALF_FLANK) next

    left_w  <- (bnd_w - HALF_FLANK):(bnd_w - 1)
    right_w <- (bnd_w + 1):(bnd_w + HALF_FLANK)
    all_w <- c(left_w, right_w)

    # Get band assignments from region around boundary
    pc1_mat <- as.matrix(dt_chr[all_w, ..pc1_cols])
    avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
    valid <- is.finite(avg_pc1)
    if (sum(valid) < 20) next

    km <- tryCatch(kmeans(avg_pc1[valid], centers = 3, nstart = 5), error = function(e) NULL)
    if (is.null(km)) next
    co <- order(km$centers[, 1])
    bands <- integer(sum(valid))
    for (b in 1:3) bands[km$cluster == co[b]] <- b
    band_snames <- sub("^PC_1_", "", names(avg_pc1)[valid])
    names(bands) <- band_snames

    b1_idx <- which(bands == 1); b3_idx <- which(bands == 3)
    if (length(b1_idx) < 5 || length(b3_idx) < 5) next

    # Flank bp ranges
    left_start  <- dt_chr$start_bp[min(left_w)]
    left_end    <- dt_chr$end_bp[max(left_w)]
    right_start <- dt_chr$start_bp[min(right_w)]
    right_end   <- dt_chr$end_bp[max(right_w)]

    fst_left <- NA_real_; fst_right <- NA_real_; hobs_het <- NA_real_
    method_used <- "proxy"

    # ── Engine B path ──────────────────────────────────────────────────
    if (.bridge_available && exists("get_region_stats", mode = "function")) {
      # Map band names to CGA if needed
      b1_names <- band_snames[b1_idx]
      b3_names <- band_snames[b3_idx]
      b2_names <- band_snames[which(bands == 2)]

      if (!is.null(smap) && grepl("^Ind[0-9]", b1_names[1])) {
        b1_cga <- smap$to_real_vec(b1_names)
        b3_cga <- smap$to_real_vec(b3_names)
        b2_cga <- smap$to_real_vec(b2_names)
      } else if (!is.null(real_names) && grepl("^Ind[0-9]", b1_names[1])) {
        ind_to_real <- setNames(real_names, paste0("Ind", seq_along(real_names) - 1))
        b1_cga <- ind_to_real[b1_names]; b1_cga <- b1_cga[!is.na(b1_cga)]
        b3_cga <- ind_to_real[b3_names]; b3_cga <- b3_cga[!is.na(b3_cga)]
        b2_cga <- ind_to_real[b2_names]; b2_cga <- b2_cga[!is.na(b2_cga)]
      } else {
        b1_cga <- b1_names; b3_cga <- b3_names; b2_cga <- b2_names
      }

      if (length(b1_cga) >= 5 && length(b3_cga) >= 5) {
        groups_arg <- list(b1 = b1_cga, b3 = b3_cga)

        fst_left <- tryCatch({
          s <- get_region_stats(chr, left_start, left_end,
                                 what = c("Fst", "Hobs"), groups = groups_arg)
          fl <- extract_fst(s)
          # Also get Hobs for HET band if available
          if (length(b2_cga) >= 3) {
            hobs_het <- tryCatch({
              sh <- get_region_stats(chr, left_start, right_end,
                                      what = "Hobs", groups = list(het = b2_cga))
              extract_hobs(sh)
            }, error = function(e) NA_real_)
          }
          fl
        }, error = function(e) NA_real_)

        fst_right <- tryCatch({
          s <- get_region_stats(chr, right_start, right_end,
                                 what = "Fst", groups = groups_arg)
          extract_fst(s)
        }, error = function(e) NA_real_)

        if (is.finite(fst_left) && is.finite(fst_right)) {
          method_used <- "engine_b"
          n_engine_b <- n_engine_b + 1L
        }
      }
    }

    # ── Fallback: PCA proxy ──────────────────────────────────────────
    if (!is.finite(fst_left) || !is.finite(fst_right)) {
      left_pc1 <- as.matrix(dt_chr[left_w, ..pc1_cols[valid]])
      right_pc1 <- as.matrix(dt_chr[right_w, ..pc1_cols[valid]])
      fst_left  <- fst_proxy(left_pc1, b1_idx, b3_idx)
      fst_right <- fst_proxy(right_pc1, b1_idx, b3_idx)
      method_used <- "proxy"
      n_proxy <- n_proxy + 1L
    }

    if (!is.finite(fst_left) || !is.finite(fst_right)) next

    fst_step <- abs(fst_left - fst_right)
    sharpness <- pmin(1, fst_step / 0.15)

    deduped[bi, cheat4_sharpness := round(sharpness, 4)]
    deduped[bi, cheat4_fst_step := round(fst_step, 4)]
    deduped[bi, cheat4_fst_left := round(fst_left, 4)]
    deduped[bi, cheat4_fst_right := round(fst_right, 4)]
    deduped[bi, cheat4_hobs_het := round(hobs_het %||% NA_real_, 4)]
    deduped[bi, cheat4_method := method_used]
    deduped[bi, cheat4_class := if (sharpness >= 0.7) "sharp"
                                 else if (sharpness >= 0.3) "moderate"
                                 else "diffuse"]
  }
  n_done <- sum(!is.na(deduped$cheat4_sharpness[chr_idx]))
  message("  ", chr, ": ", n_done, "/", length(chr_idx), " boundaries scored",
          " (engine_b=", n_engine_b, " proxy=", n_proxy, ")")
}

# =============================================================================
# CHEAT 10: READ DEPTH ANOMALY (real mosdepth + sim_mat fallback)
# =============================================================================
# Real version: reads mosdepth per-sample BED files, computes depth dip/spike
# at each boundary, compares HET vs HOM samples.
# Fallback: uses sim_mat boundary isolation score when mosdepth unavailable.
#
# DATA: mosdepth per-sample: {sample}.markdup.mosdepth.regions.bed.gz
#   (from delly_sv/00_markdup/ or MODULE_4D_PAV/mosdepth_output/)

message("\n[C01g] Running Cheat 10 (depth anomaly)...")

# Helper: find mosdepth BED for a sample
find_mosdepth_bed <- function(sample_id) {
  # Primary: mosdepth_dir (from --mosdepth_dir or auto-detected MARKDUP_DIR)
  if (!is.null(mosdepth_dir)) {
    for (pattern in c(
      paste0(sample_id, ".markdup.mosdepth.regions.bed.gz"),
      paste0(sample_id, ".mosdepth.regions.bed.gz"),
      paste0(sample_id, ".regions.bed.gz")
    )) {
      f <- file.path(mosdepth_dir, pattern)
      if (file.exists(f)) return(f)
    }
  }
  # Fallback: search known locations
  for (d in c(file.path(BASE, "delly_sv/00_markdup"),
              file.path(BASE, "MODULE_4D_PAV/mosdepth_output"),
              file.path(BASE, "MODULE_4H_ALL_Manta/00_markdup"))) {
    for (pattern in c(
      paste0(sample_id, ".markdup.mosdepth.regions.bed.gz"),
      paste0(sample_id, ".mosdepth.regions.bed.gz")
    )) {
      f <- file.path(d, pattern)
      if (file.exists(f)) return(f)
    }
  }
  NULL
}

# Check if mosdepth is available
mosdepth_available <- FALSE
test_sample <- if (!is.null(real_names) && length(real_names) > 0) real_names[1] else sample_names[1]
if (!is.null(find_mosdepth_bed(test_sample))) {
  mosdepth_available <- TRUE
  message("  Mosdepth files FOUND — using real depth analysis")
} else {
  message("  Mosdepth files NOT found — using sim_mat isolation proxy")
}

deduped[, cheat10_depth_dip := NA_real_]
deduped[, cheat10_depth_ratio := NA_real_]
deduped[, cheat10_cv_het := NA_real_]
deduped[, cheat10_cv_hom := NA_real_]
deduped[, cheat10_score := NA_real_]
deduped[, cheat10_class := NA_character_]
deduped[, cheat10_source := NA_character_]

for (chr in unique(deduped$chrom)) {
  pc <- precomp_list[[chr]]
  if (is.null(pc)) next
  dt_chr <- pc$dt; sim_mat <- pc$sim_mat
  n_w <- nrow(dt_chr)
  mid_bps <- (dt_chr$start_bp + dt_chr$end_bp) / 2
  pc1_cols <- intersect(paste0("PC_1_", sample_names), names(dt_chr))

  # Get band assignments for per-band depth analysis
  chr_bands <- NULL
  if (length(pc1_cols) >= 20) {
    avg_all <- colMeans(as.matrix(dt_chr[, ..pc1_cols]), na.rm = TRUE)
    valid <- is.finite(avg_all)
    if (sum(valid) >= 20) {
      km <- tryCatch(kmeans(avg_all[valid], centers = 3, nstart = 5), error = function(e) NULL)
      if (!is.null(km)) {
        co <- order(km$centers[, 1])
        chr_bands <- integer(sum(valid))
        for (b in 1:3) chr_bands[km$cluster == co[b]] <- b
        snames <- sub("^PC_1_", "", names(avg_all)[valid])
        # Map to real names if needed
        if (!is.null(real_names) && grepl("^Ind", snames[1])) {
          ind_to_real <- setNames(real_names, paste0("Ind", seq_along(real_names) - 1))
          names(chr_bands) <- ind_to_real[snames]
        } else {
          names(chr_bands) <- snames
        }
        chr_bands <- chr_bands[!is.na(names(chr_bands))]
      }
    }
  }

  chr_idx <- which(deduped$chrom == chr)

  if (mosdepth_available && !is.null(chr_bands)) {
    # ── REAL DEPTH ANALYSIS ──
    # Subsample 30 samples stratified by band
    depth_samples <- names(chr_bands)
    if (length(depth_samples) > 30) {
      depth_samples <- c(
        sample(names(chr_bands)[chr_bands == 1], min(10, sum(chr_bands == 1))),
        sample(names(chr_bands)[chr_bands == 2], min(10, sum(chr_bands == 2))),
        sample(names(chr_bands)[chr_bands == 3], min(10, sum(chr_bands == 3)))
      )
    }

    for (bi in chr_idx) {
      bp <- deduped$boundary_bp[bi]
      bnd_w <- which.min(abs(mid_bps - bp))
      if (bnd_w <= 10 || bnd_w >= n_w - 9) next

      # Load depth ±200kb each side
      left_start <- dt_chr$start_bp[max(1, bnd_w - 80)]
      left_end <- dt_chr$end_bp[max(1, bnd_w - 1)]
      right_start <- dt_chr$start_bp[min(n_w, bnd_w + 1)]
      right_end <- dt_chr$end_bp[min(n_w, bnd_w + 80)]

      depth_left <- list(); depth_right <- list(); depth_dip_vals <- list()
      for (sid in depth_samples) {
        bed <- find_mosdepth_bed(sid)
        if (is.null(bed)) next

        # Query left side
        cmd_l <- if (file.exists(paste0(bed, ".tbi"))) {
          sprintf("tabix %s %s:%d-%d 2>/dev/null", bed, chr, left_start, left_end)
        } else {
          sprintf("zcat %s | awk '$1==\"%s\" && $2>=%d && $3<=%d'", bed, chr, left_start, left_end)
        }
        dl <- tryCatch(fread(cmd = cmd_l, header = FALSE, sep = "\t"), error = function(e) NULL)

        # Query right side
        cmd_r <- if (file.exists(paste0(bed, ".tbi"))) {
          sprintf("tabix %s %s:%d-%d 2>/dev/null", bed, chr, right_start, right_end)
        } else {
          sprintf("zcat %s | awk '$1==\"%s\" && $2>=%d && $3<=%d'", bed, chr, right_start, right_end)
        }
        dr <- tryCatch(fread(cmd = cmd_r, header = FALSE, sep = "\t"), error = function(e) NULL)

        # Query boundary ±5kb
        cmd_d <- if (file.exists(paste0(bed, ".tbi"))) {
          sprintf("tabix %s %s:%d-%d 2>/dev/null", bed, chr, bp - 5000, bp + 5000)
        } else {
          sprintf("zcat %s | awk '$1==\"%s\" && $2>=%d && $3<=%d'", bed, chr, bp - 5000, bp + 5000)
        }
        dd <- tryCatch(fread(cmd = cmd_d, header = FALSE, sep = "\t"), error = function(e) NULL)

        ml <- if (!is.null(dl) && nrow(dl) > 0) mean(dl[[4]], na.rm = TRUE) else NA_real_
        mr <- if (!is.null(dr) && nrow(dr) > 0) mean(dr[[4]], na.rm = TRUE) else NA_real_
        md <- if (!is.null(dd) && nrow(dd) > 0) mean(dd[[4]], na.rm = TRUE) else NA_real_

        band_id <- chr_bands[sid]
        depth_left[[sid]] <- list(depth = ml, band = band_id)
        depth_right[[sid]] <- list(depth = mr, band = band_id)
        depth_dip_vals[[sid]] <- list(depth = md, band = band_id)
      }

      # Compute metrics
      left_depths <- sapply(depth_left, function(x) x$depth)
      right_depths <- sapply(depth_right, function(x) x$depth)
      dip_depths <- sapply(depth_dip_vals, function(x) x$depth)
      bands_vec <- sapply(depth_left, function(x) x$band)

      mean_left <- mean(left_depths, na.rm = TRUE)
      mean_right <- mean(right_depths, na.rm = TRUE)
      mean_dip <- mean(dip_depths, na.rm = TRUE)
      context_mean <- mean(c(mean_left, mean_right), na.rm = TRUE)

      depth_ratio <- if (max(mean_left, mean_right, na.rm = TRUE) > 0)
        min(mean_left, mean_right, na.rm = TRUE) / max(mean_left, mean_right, na.rm = TRUE)
      else NA_real_

      depth_dip <- if (is.finite(mean_dip) && context_mean > 0)
        1 - mean_dip / context_mean else NA_real_

      # Per-band CV
      het_d <- dip_depths[bands_vec == 2]; het_d <- het_d[is.finite(het_d)]
      hom_d <- dip_depths[bands_vec %in% c(1, 3)]; hom_d <- hom_d[is.finite(hom_d)]
      cv_het <- if (length(het_d) >= 3 && mean(het_d) > 0) sd(het_d) / mean(het_d) else NA_real_
      cv_hom <- if (length(hom_d) >= 3 && mean(hom_d) > 0) sd(hom_d) / mean(hom_d) else NA_real_

      # Score
      score <- 0; nc <- 0
      if (is.finite(depth_dip) && abs(depth_dip) > 0.05) {
        score <- score + pmin(1, abs(depth_dip) / 0.3) * 0.4; nc <- nc + 1 }
      if (is.finite(depth_ratio) && depth_ratio < 0.9) {
        score <- score + (1 - depth_ratio) * 0.3; nc <- nc + 1 }
      if (is.finite(cv_het) && is.finite(cv_hom) && cv_het > cv_hom * 1.2) {
        score <- score + pmin(1, (cv_het - cv_hom) / 0.1) * 0.3; nc <- nc + 1 }
      if (nc > 0) score <- pmin(1, score / (nc * 0.35))

      deduped[bi, cheat10_depth_dip := round(depth_dip, 4)]
      deduped[bi, cheat10_depth_ratio := round(depth_ratio, 4)]
      deduped[bi, cheat10_cv_het := round(cv_het, 4)]
      deduped[bi, cheat10_cv_hom := round(cv_hom, 4)]
      deduped[bi, cheat10_score := round(score, 4)]
      deduped[bi, cheat10_source := "mosdepth"]
      deduped[bi, cheat10_class := if (score > 0.5) "strong_dip"
                                    else if (score > 0.2) "moderate_dip"
                                    else "no_dip"]
    }
  } else {
    # ── SIM_MAT ISOLATION PROXY (fallback) ──
    for (bi in chr_idx) {
      bp <- deduped$boundary_bp[bi]
      bnd_w <- which.min(abs(mid_bps - bp))
      if (bnd_w <= 10 || bnd_w >= n_w - 9) next

      left_nb  <- max(1, bnd_w - 10):(bnd_w - 1)
      right_nb <- (bnd_w + 1):min(n_w, bnd_w + 10)

      sim_to_left  <- mean(sim_mat[bnd_w, left_nb], na.rm = TRUE)
      sim_to_right <- mean(sim_mat[bnd_w, right_nb], na.rm = TRUE)
      left_internal  <- mean(sim_mat[left_nb, left_nb], na.rm = TRUE)
      right_internal <- mean(sim_mat[right_nb, right_nb], na.rm = TRUE)

      avg_to_bnd <- mean(c(sim_to_left, sim_to_right), na.rm = TRUE)
      avg_internal <- mean(c(left_internal, right_internal), na.rm = TRUE)
      isolation <- avg_internal - avg_to_bnd
      if (!is.finite(isolation)) next

      score <- pmin(1, pmax(0, isolation / 0.2))
      deduped[bi, cheat10_score := round(score, 4)]
      deduped[bi, cheat10_source := "sim_mat_proxy"]
      deduped[bi, cheat10_class := if (isolation > 0.15) "strong_dip"
                                    else if (isolation > 0.05) "moderate_dip"
                                    else "no_dip"]
    }
  }
  n_scored <- sum(!is.na(deduped$cheat10_score[chr_idx]))
  message("  ", chr, ": ", n_scored, "/", length(chr_idx),
          " boundaries (", if (mosdepth_available) "mosdepth" else "sim_mat proxy", ")")
}

# =============================================================================
# CHEAT 11: CLIPPED READ PILEUP (real BAM + band concordance fallback)
# =============================================================================
# Real version: samtools view on ±2kb around boundary, parse CIGAR for S/H,
# count per sample, test bimodality and per-band enrichment.
# Fallback: cross-boundary band concordance (from PC1 loadings).
#
# DATA: markdup BAMs from bam_manifest or delly_sv/00_markdup/

message("\n[C01g] Running Cheat 11 (clipped reads / band concordance)...")

# Helper: find BAM for a sample
find_bam_file <- function(sample_id) {
  # From manifest (column 1 = sample, column 3 = BAM path — Quentin's standard)
  if (!is.null(bam_manifest) && file.exists(bam_manifest)) {
    manifest <- tryCatch(fread(bam_manifest, header = FALSE), error = function(e) NULL)
    if (!is.null(manifest) && ncol(manifest) >= 3) {
      idx <- which(manifest[[1]] == sample_id)
      if (length(idx) > 0) {
        bam <- manifest[[3]][idx[1]]
        if (file.exists(bam)) return(bam)
      }
    }
  }

  # Standard locations: delly_sv/00_markdup and Manta
  for (d in c(MARKDUP_DIR,
              file.path(BASE, "MODULE_4H_ALL_Manta/00_markdup"))) {
    for (pattern in c(paste0(sample_id, ".markdup.bam"),
                       paste0(sample_id, ".bam"))) {
      f <- file.path(d, pattern)
      if (file.exists(f)) return(f)
    }
  }
  NULL
}

# Count clips in a narrow region for one sample
count_clips_narrow <- function(bam_path, chr, start_bp, end_bp) {
  if (!file.exists(bam_path)) return(list(n_total = NA, n_clip = NA, clip_frac = NA))
  region <- paste0(chr, ":", as.integer(start_bp), "-", as.integer(end_bp))
  cmd <- paste0(
    "samtools view -F 0x904 ", bam_path, " ", region, " 2>/dev/null | ",
    "awk '{total++; cigar=$6; ",
    "while (match(cigar, /([0-9]+)[SH]/, a)) {",
    "  if (a[1]+0>=5) {sc++; break}; ",
    "  cigar=substr(cigar, RSTART+RLENGTH)}} ",
    "END{print total+0, sc+0}'"
  )
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) NULL)
  if (is.null(out) || length(out) == 0) return(list(n_total = NA, n_clip = NA, clip_frac = NA))
  parts <- as.integer(strsplit(trimws(out[1]), " ")[[1]])
  if (length(parts) < 2) return(list(n_total = NA, n_clip = NA, clip_frac = NA))
  list(n_total = parts[1], n_clip = parts[2],
       clip_frac = if (parts[1] > 0) parts[2] / parts[1] else NA)
}

# Check if BAMs are available
bam_available <- FALSE
test_bam <- find_bam_file(test_sample)
if (!is.null(test_bam)) {
  bam_available <- TRUE
  message("  BAM files FOUND — using real clip analysis")
} else {
  message("  BAM files NOT found — using band concordance proxy")
}

deduped[, cheat11_clip_score := NA_real_]
deduped[, cheat11_clip_enrichment := NA_real_]
deduped[, cheat11_clip_bimodal := NA]
deduped[, cheat11_concordance := NA_real_]
deduped[, cheat11_jaccard := NA_real_]
deduped[, cheat11_class := NA_character_]
deduped[, cheat11_source := NA_character_]

for (chr in unique(deduped$chrom)) {
  pc <- precomp_list[[chr]]
  if (is.null(pc)) next
  dt_chr <- pc$dt; n_w <- nrow(dt_chr)
  mid_bps <- (dt_chr$start_bp + dt_chr$end_bp) / 2
  pc1_cols <- intersect(paste0("PC_1_", sample_names), names(dt_chr))

  chr_idx <- which(deduped$chrom == chr)

  if (bam_available) {
    # ── REAL CLIP ANALYSIS ──
    # Subsample 40 samples stratified by band
    clip_samples <- if (!is.null(real_names)) real_names else sample_names
    if (length(clip_samples) > 40) clip_samples <- sample(clip_samples, 40)

    for (bi in chr_idx) {
      bp <- deduped$boundary_bp[bi]

      # Narrow ±2kb clip count per sample
      clip_fracs <- numeric(0)
      ctx_fracs <- numeric(0)
      for (sid in clip_samples) {
        bam <- find_bam_file(sid)
        if (is.null(bam)) next
        narrow <- count_clips_narrow(bam, chr, bp - 2000, bp + 2000)
        context <- count_clips_narrow(bam, chr, bp - 50000, bp + 50000)
        if (is.finite(narrow$clip_frac)) clip_fracs <- c(clip_fracs, narrow$clip_frac)
        if (is.finite(context$clip_frac)) ctx_fracs <- c(ctx_fracs, context$clip_frac)
      }

      if (length(clip_fracs) < 10) next

      mean_clip <- mean(clip_fracs, na.rm = TRUE)
      mean_ctx <- mean(ctx_fracs, na.rm = TRUE)
      enrichment <- if (mean_ctx > 0) mean_clip / mean_ctx else NA_real_

      # Bimodality check
      cv_clip <- if (mean_clip > 0) sd(clip_fracs) / mean_clip else 0
      bimodal <- cv_clip > 0.5

      # Score
      score <- 0
      if (is.finite(enrichment) && enrichment > 1.5)
        score <- score + pmin(1, (enrichment - 1) / 3) * 0.5
      if (bimodal) score <- score + 0.3
      score <- pmin(1, score)

      deduped[bi, cheat11_clip_score := round(score, 4)]
      deduped[bi, cheat11_clip_enrichment := round(enrichment, 2)]
      deduped[bi, cheat11_clip_bimodal := bimodal]
      deduped[bi, cheat11_source := "bam_clips"]
      deduped[bi, cheat11_class := if (score > 0.5) "strong_clip_pileup"
                                    else if (score > 0.2) "moderate_clips"
                                    else "no_clip_signal"]
    }
  }

  # ── BAND CONCORDANCE (always computed — complement to clips) ──
  if (length(pc1_cols) >= 20) {
    get_side_bands <- function(w_idx) {
      mat <- as.matrix(dt_chr[w_idx, ..pc1_cols])
      avg <- colMeans(mat, na.rm = TRUE)
      valid <- is.finite(avg)
      if (sum(valid) < 20) return(NULL)
      km <- tryCatch(kmeans(avg[valid], centers = 3, nstart = 5), error = function(e) NULL)
      if (is.null(km)) return(NULL)
      co <- order(km$centers[, 1])
      bands <- character(sum(valid))
      for (b in 1:3) bands[km$cluster == co[b]] <- paste0("band", b)
      names(bands) <- sub("^PC_1_", "", names(avg)[valid])
      bands
    }

    for (bi in chr_idx) {
      bp <- deduped$boundary_bp[bi]
      bnd_w <- which.min(abs(mid_bps - bp))
      if (bnd_w <= 15 || bnd_w >= n_w - 14) next

      left_bands  <- get_side_bands((bnd_w - 15):(bnd_w - 2))
      right_bands <- get_side_bands((bnd_w + 2):(bnd_w + 15))
      if (is.null(left_bands) || is.null(right_bands)) next

      shared <- intersect(names(left_bands), names(right_bands))
      if (length(shared) < 15) next

      lb <- left_bands[shared]; rb <- right_bands[shared]
      concordance <- mean(lb == rb)

      jac_vals <- numeric(3)
      for (b in paste0("band", 1:3)) {
        l_set <- names(lb)[lb == b]; r_set <- names(rb)[rb == b]
        inter <- length(intersect(l_set, r_set))
        union <- length(union(l_set, r_set))
        jac_vals[as.integer(sub("band", "", b))] <- if (union > 0) inter / union else 0
      }
      jaccard <- mean(jac_vals)

      deduped[bi, cheat11_concordance := round(concordance, 4)]
      deduped[bi, cheat11_jaccard := round(jaccard, 4)]

      # Only set class from concordance if BAM wasn't available
      if (!bam_available || is.na(deduped$cheat11_clip_score[bi])) {
        deduped[bi, cheat11_source := "band_concordance"]
        deduped[bi, cheat11_class := if (concordance > 0.85 && jaccard > 0.7) "same_system"
                                      else if (concordance < 0.50) "different_systems"
                                      else if (concordance < 0.65) "system_transition"
                                      else "partial_overlap"]
      }
    }
  }
  n_scored <- sum(!is.na(deduped$cheat11_class[chr_idx]))
  message("  ", chr, ": ", n_scored, "/", length(chr_idx),
          " boundaries (", if (bam_available) "BAM clips" else "band concordance", ")")
}

# =============================================================================
# CHEAT 17: ARCHIVED 2026-04-17
# =============================================================================
# The fossil-breakpoint classification block used to live here. It created
# a circular dependency: Cheat 17 reads C01d's candidate_scores.tsv.gz via
# --scores to find boundaries OUTSIDE any candidate, while C01d reads this
# script's boundary_catalog_unified.tsv.gz for D11. Getting both fully
# populated required a two-pass run (C01g → C01d → C01g), which the
# orchestrator didn't wire up. In single-pass runs the block did nothing.
#
# The fossil-vs-active distinction was not used in the manuscript main
# path — it was a diagnostic annotation. Downstream code no longer
# references cheat17_class; boundary_activity is now computed
# unconditionally as "active" (all kept boundaries are treated as
# currently-segregating). If the fossil distinction becomes needed
# again, the archived code is preserved in
# _archive_superseded/cheat17_fossil_detection/ and can be brought back
# as a post-hoc annotation step reading the already-written
# boundary_catalog_unified.tsv.gz.
#
# Initialize cheat17 columns as NA for schema-stability with any
# downstream readers that still expect them.
deduped[, cheat17_class := NA_character_]
deduped[, cheat17_inv_likeness := NA_real_]

# =============================================================================
# CHEAT 21: TE ENRICHMENT/DEPLETION AT BREAKPOINTS
# =============================================================================
# Compare TE density at breakpoints vs random genome positions.
# Enriched → repeat-mediated mechanism (NAHR substrate).
# Depleted → structural constraint (breakpoint avoids TEs).

message("\n[C01g] Running Cheat 21 (TE enrichment)...")

deduped[, cheat21_te_density := NA_real_]
deduped[, cheat21_te_count := NA_integer_]
deduped[, cheat21_verdict := NA_character_]

# Auto-detect RepeatMasker file
if (is.null(repeatmasker_file)) {
  for (rp in c(
    file.path(BASE, "genome_annotation/repeatmasker/catfish_ref.fa.out"),
    file.path(BASE, "genome_annotation/repeatmasker/catfish_ref.repeatmasker.bed"),
    Sys.glob(file.path(BASE, "genome_annotation/repeatmasker/*.out"))[1],
    Sys.glob(file.path(BASE, "genome_annotation/repeatmasker/*.bed"))[1]
  )) {
    if (!is.na(rp) && file.exists(rp)) { repeatmasker_file <- rp; break }
  }
}

if (!is.null(repeatmasker_file) && file.exists(repeatmasker_file)) {
  message("[C01g] RepeatMasker: ", repeatmasker_file)

  # Load TE annotations
  rm_ext <- tools::file_ext(repeatmasker_file)
  te_dt <- tryCatch({
    if (rm_ext %in% c("bed", "gz")) {
      dt <- fread(repeatmasker_file)
      if (ncol(dt) >= 4 && !"chr" %in% names(dt))
        setnames(dt, seq_len(min(4, ncol(dt))),
                 c("chr", "start", "end", "te_class")[seq_len(min(4, ncol(dt)))])
      dt
    } else {
      # Parse .out format
      lines <- readLines(repeatmasker_file)
      lines <- lines[!grepl("^\\s*$|^\\s*SW", lines)]
      if (length(lines) > 3) lines <- lines[-(1:3)]
      raw <- fread(text = lines, fill = TRUE)
      if (ncol(raw) >= 11) {
        data.table(chr = raw[[5]], start = as.integer(raw[[6]]),
                    end = as.integer(raw[[7]]),
                    te_class = sub("/.*", "", raw[[11]]))
      } else data.table()
    }
  }, error = function(e) { message("[C01g] TE parse failed: ", e$message); data.table() })

  if (nrow(te_dt) > 0) {
    te_dt[, te_bp := end - start]
    TE_WINDOW <- 10000L  # ±10 kb

    for (chr in unique(deduped$chrom)) {
      chr_idx <- which(deduped$chrom == chr)
      if (length(chr_idx) == 0) next
      chr_te <- te_dt[chr == (..chr)]
      if (nrow(chr_te) == 0) {
        # BUGFIX 2026-04-17: was `grepl(chr, chr, fixed = TRUE)` — pattern
        # and target both the scalar `chr`, which is always TRUE and
        # returned the entire genome's TE set. Now targets te_dt$chr
        # against a simplified chr name AND against the full chr name
        # so naming mismatches like LG01 vs C_gar_LG01 are tolerated.
        short_chr <- sub("^[Cc]_?[A-Za-z]+_", "", chr)  # C_gar_LG01 → LG01
        chr_te <- te_dt[
          grepl(short_chr, te_dt$chr, fixed = TRUE) |
          grepl(chr, te_dt$chr, fixed = TRUE)
        ]
      }
      if (nrow(chr_te) == 0) next

      # Genome-wide TE density for this chr
      chr_len <- max(chr_te$end, na.rm = TRUE)
      genome_te_density <- sum(chr_te$te_bp) / chr_len

      for (bi in chr_idx) {
        bp <- deduped$boundary_bp[bi]
        local_te <- chr_te[start <= bp + TE_WINDOW & end >= bp - TE_WINDOW]
        local_bp <- sum(pmin(local_te$end, bp + TE_WINDOW) -
                         pmax(local_te$start, bp - TE_WINDOW))
        local_bp <- max(0, local_bp)
        local_density <- local_bp / (2 * TE_WINDOW)

        ratio <- if (genome_te_density > 0) local_density / genome_te_density else NA_real_

        verdict <- if (!is.finite(ratio)) "NO_DATA"
          else if (ratio > 1.5) "ENRICHED"
          else if (ratio < 0.5) "DEPLETED"
          else "NEUTRAL"

        deduped[bi, cheat21_te_density := round(local_density, 4)]
        deduped[bi, cheat21_te_count := nrow(local_te)]
        deduped[bi, cheat21_verdict := verdict]
      }
      n_enriched <- sum(deduped$cheat21_verdict[chr_idx] == "ENRICHED", na.rm = TRUE)
      n_depleted <- sum(deduped$cheat21_verdict[chr_idx] == "DEPLETED", na.rm = TRUE)
      message("  ", chr, ": ", n_enriched, " enriched, ", n_depleted, " depleted")
    }
  }
} else {
  message("[C01g] No RepeatMasker file found — Cheat 21 skipped")
}

# =============================================================================
# BUILD CONCORDANCE SUMMARY
# =============================================================================

message("\n[C01g] Building concordance summary...")

# ── STRUCTURAL CONCORDANCE (purely: is this a real boundary?) ──
# Cheats 4/10/11/21 + source counts. No fossil penalty — fossils
# are structurally real, just evolutionarily inactive.
deduped[, n_cheats_supporting := {
  n <- 0L
  if (!is.na(cheat4_sharpness) && cheat4_sharpness >= 0.3) n <- n + 1L
  if (!is.na(cheat10_score) && cheat10_score >= 0.3) n <- n + 1L
  if (!is.na(cheat11_clip_score) && cheat11_clip_score >= 0.3) n <- n + 1L
  if (is.na(cheat11_clip_score) && !is.na(cheat11_concordance) && cheat11_concordance > 0.70) n <- n + 1L
  # Source concordance counts too
  if (n_sources >= 2) n <- n + 1L
  if (grepl("sv_", sources)) n <- n + 1L  # SV evidence is strong
  # Cheat 21: TE enrichment = structural context (repeat-mediated substrate)
  if (!is.na(cheat21_verdict) && cheat21_verdict == "ENRICHED") n <- n + 1L
  n
}, by = boundary_id]

# ── BOUNDARY ACTIVITY ──
# ARCHIVED 2026-04-17: Fossil vs active distinction is no longer computed
# (see Cheat 17 archive note above). All kept boundaries are treated as
# currently-segregating. The column is still emitted as "active" so
# downstream readers/summary expecting it don't break.
deduped[, boundary_activity := "active"]

deduped[, boundary_verdict := fifelse(
  n_cheats_supporting >= 4, "confirmed_structural",
  fifelse(n_cheats_supporting == 3, "likely_structural",
  fifelse(n_cheats_supporting == 2, "probable_structural",
  fifelse(n_cheats_supporting == 1 &
            boundary_type %in% c("STRUCTURAL_SHARP", "STRUCTURAL_MODERATE"),
          "single_evidence_structural",
  fifelse(boundary_type == "SYSTEM_JUNCTION", "system_junction",
  fifelse(boundary_type %in% c("HOTCORE_EDGE", "INTERNAL_GAP",
                                 "INTERNAL_TRANSITION"), "internal_feature",
  fifelse(boundary_type == "BACKGROUND_FADE", "background",
          "unresolved"))))))
)]

# =============================================================================
# WRITE
# =============================================================================

message("\n[C01g] Writing outputs...")

fwrite(deduped, file.path(outdir, "boundary_catalog_unified.tsv.gz"), sep = "\t")

# Per-chromosome files for downstream scripts
for (chr in unique(deduped$chrom)) {
  fwrite(deduped[chrom == chr],
         file.path(outdir, paste0("boundary_catalog_", chr, ".tsv.gz")), sep = "\t")
}

# Summary table
summ <- deduped[, .(
  n_boundaries = .N,
  n_confirmed = sum(boundary_verdict == "confirmed_structural"),
  n_likely = sum(boundary_verdict == "likely_structural"),
  n_probable = sum(boundary_verdict == "probable_structural"),
  n_unresolved = sum(boundary_verdict == "unresolved"),
  n_internal = sum(boundary_verdict == "internal_feature"),
  n_background = sum(boundary_verdict == "background"),
  n_sv_confirmed = sum(grepl("sv_", sources)),
  mean_hardness = round(mean(hardness, na.rm = TRUE), 4),
  mean_sharpness = round(mean(cheat4_sharpness, na.rm = TRUE), 4),
  n_fossil = sum(boundary_activity == "fossil", na.rm = TRUE),
  n_active = sum(boundary_activity == "active", na.rm = TRUE),
  n_te_enriched = sum(cheat21_verdict == "ENRICHED", na.rm = TRUE)
), by = chrom]
fwrite(summ, file.path(outdir, "boundary_summary.tsv"), sep = "\t")

# Final summary
message("\n[C01g] ═══ BOUNDARY CATALOG SUMMARY ═══")
message("  Total boundaries: ", nrow(deduped))
tab <- table(deduped$boundary_verdict)
for (v in sort(names(tab))) message("  ", v, ": ", tab[v])
message("\n  Sources breakdown:")
src_tab <- table(unlist(strsplit(deduped$sources, ";")))
for (s in sort(names(src_tab))) message("    ", s, ": ", src_tab[s])
message("\n  Multi-source boundaries: ", sum(deduped$n_sources >= 2),
        " / ", nrow(deduped), " (", round(100 * mean(deduped$n_sources >= 2), 1), "%)")


# ═══════════════════════════════════════════════════════════════════════
# REGISTER IN EVIDENCE REGISTRY
# ═══════════════════════════════════════════════════════════════════════
tryCatch({
  for (.hf in c("utils/registry_key_helpers.R", "../utils/registry_key_helpers.R")) {
    if (file.exists(.hf)) { source(.hf); break }
  }

  if (.bridge_available && exists("reg") && !is.null(reg) && !is.null(reg$add_evidence) &&
      nrow(deduped) > 0) {
    message("[C01g] Registering boundaries in evidence registry...")
    n_keys_total <- 0L
    confirmed <- deduped[boundary_verdict %in% c("confirmed_structural", "likely_structural")]
    if (nrow(confirmed) > 0 && "matched_core_id" %in% names(confirmed)) {
      for (bi in seq_len(nrow(confirmed))) {
        bd <- confirmed[bi]
        core_id <- bd$matched_core_id
        if (is.na(core_id) || core_id == "") next
        cid <- paste0(bd$chrom, "_", core_id)
        side <- if (!is.null(bd$side) && !is.na(bd$side)) bd$side
                else if (bi %% 2 == 1) "left" else "right"

        if (exists("register_C01g_boundary", mode = "function")) {
          n_keys_total <- n_keys_total + register_C01g_boundary(bd, cid, side, outdir)
        }
        if (exists("store_C01g_boundary", mode = "function")) {
          store_C01g_boundary(bd, cid, side, outdir)
        }
      }
    }
    message("[C01g] Registered ", n_keys_total, " boundary evidence keys")
  }
}, error = function(e) message("[C01g] Registry wiring: ", conditionMessage(e)))

message("\n[DONE] -> ", outdir)
