#!/usr/bin/env Rscript
# =============================================================================
# STEP_C01i_b_multi_recomb.R — chat-12 rewrite: R(DAG) + G combination rule
# =============================================================================
# PRIOR STATE (removed):
#   - Inline cheat24 fallback ("classify_recombinant_event") that disagreed
#     with the real cheat24 on the gene_conversion posterior (Finding V,
#     chats 8/11). Deleted; if the real cheat24 isn't available we flag
#     but do not invent numbers.
#   - 100 kb mosaic-length threshold as a standalone gate for Signal 1
#     recombinant calls. Removed per chat-9 design (see HANDOFF §"Files
#     chat 12 SHOULD touch").
#   - Signal 1 (per-window PCA class switch) as a GATE. The PCA per-window
#     decomposition is a noisy 1D projection of regime membership — C01j
#     (sliding-window compatibility clustering) gives the regime track we
#     actually want. Signal 1 is KEPT but demoted to Tier-4 diagnostic:
#     recorded in the per-sample record for auditability, not used for
#     classification.
#
# NEW GATE (chat-12):
#   R AND G            → RECOMBINANT, HIGH
#   R AND !G           → RECOMBINANT, MEDIUM (regime_only_ghsl_insufficient)
#                          or recomb_disputed if G was conclusive NON-SPLIT
#   !R AND G           → recomb_ghsl_only
#   !R AND !G          → NOT_RECOMB (Tier-1 decompose class stands)
#   Orthogonally, GC tracts (from gene_conversion_detector) are recorded
#   as annotation per sample but do NOT gate the call.
#
# Event-class split:
#   Within RECOMBINANT, samples get subdivided:
#     RECOMBINANT_GC  if GC tracts > 0 and no long deviation (span_bp based)
#     RECOMBINANT_DCO if longest_deviation_bp > params$min_dco_bp
#     RECOMBINANT     otherwise (single-crossover or ambiguous)
#   These drive the per-candidate subgroup names registered by
#   STEP_C01i_d_seal (chat-13 wiring).
#
# INPUTS:
#   --candidates      TSV with candidate_id, chrom, start_bp, end_bp, tier
#   --regime_memb     regime_memberships.tsv.gz from C01j (required)
#   --decomp_dir      directory with per-candidate internal_dynamics block
#                      (for the baseline class per sample); the legacy
#                      per_window_class.rds used for the old Signal 1 is
#                      now optional
#   --ghsl_dir        directory with per-candidate GHSL output
#                      (from lib_ghsl_confirmation — ghsl_per_sample.tsv)
#   --gc_dir          OPTIONAL directory with per-candidate GC detector
#                      output (gene_conversion_tracts.json)
#   --outdir          fallback output root
#   --tier_max        default 3
#   --min_dev_frac    DAG R_fired gate on deviation_fraction (default 0.15)
#   --min_dev_bp      DAG R_fired gate on longest_deviation_bp (default 50000)
#   --min_dco_bp      longest deviation above this → RECOMBINANT_DCO (def 200000)
#
# OUTPUT: recombinant_map.json per candidate plus the new regime_sample_dag
#   block (for C01d Layer C scoring in chat-13 wiring).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

option_list <- list(
  make_option("--candidates",   type = "character", help = "Candidate table"),
  make_option("--regime_memb",  type = "character", default = NA_character_,
              help = "Path to regime_memberships.tsv.gz from C01j (required)"),
  make_option("--decomp_dir",   type = "character", default = "decomp_out",
              help = "Where STEP_C01i_decompose wrote per-candidate blocks"),
  make_option("--ghsl_dir",     type = "character", default = NA_character_,
              help = "Per-candidate GHSL output directory (optional)"),
  make_option("--gc_dir",       type = "character", default = NA_character_,
              help = "Per-candidate GC detector output directory (optional)"),
  make_option("--outdir",       type = "character", default = "recomb_out"),
  make_option("--tier_max",     type = "integer",   default = 3L),
  make_option("--min_dev_frac", type = "double",    default = 0.15),
  make_option("--min_dev_bp",   type = "integer",   default = 50000L),
  make_option("--min_dco_bp",   type = "integer",   default = 200000L)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$candidates) || !file.exists(opt$candidates)) {
  stop("--candidates required and must exist")
}
if (is.na(opt$regime_memb) || !file.exists(opt$regime_memb)) {
  stop("--regime_memb required (C01j regime_memberships.tsv.gz)")
}

# ── Helpers + library ───────────────────────────────────────────────────────
helper_paths <- c("lib_decompose_helpers.R", "R/lib_decompose_helpers.R",
                  file.path(dirname(sys.frame(1)$ofile %||% "."),
                             "lib_decompose_helpers.R"))
for (p in helper_paths) if (file.exists(p)) { source(p); break }

combo_paths <- c("lib_recomb_combination.R",
                 "R/lib_recomb_combination.R",
                 file.path(dirname(sys.frame(1)$ofile %||% "."),
                            "lib_recomb_combination.R"))
combo_found <- FALSE
for (p in combo_paths) if (file.exists(p)) {
  source(p); combo_found <- TRUE; break
}
if (!combo_found) stop("[multi_recomb] cannot find lib_recomb_combination.R")

paths <- resolve_decomp_paths(list(outdir = opt$outdir))
reg   <- try_load_registry()
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Source real cheat24 if available (Finding V: no inline fallback) ────────
.has_cheat24 <- FALSE
for (p in c("cheats/cheat24_recombinant_prior.R",
            "../cheats/cheat24_recombinant_prior.R",
            file.path(paths$base, "flashlight_v2", "cheats",
                       "cheat24_recombinant_prior.R"))) {
  if (file.exists(p)) {
    tryCatch({ source(p); .has_cheat24 <- TRUE; break },
             error = function(e) message("[multi_recomb] cheat24 source fail: ",
                                          conditionMessage(e)))
  }
}
if (!.has_cheat24) {
  message("[multi_recomb] cheat24 unavailable — recombinant_map will set ",
          "event_class=NA and posterior=NA for the affected samples. ",
          "The inline fallback from earlier drafts has been REMOVED ",
          "(Finding V): fabricating a posterior that disagrees with the ",
          "real cheat24 is worse than returning NA.")
}

# ─────────────────────────────────────────────────────────────────────────────
# Optional Tier-4 diagnostic: Signal 1 (per-window PCA class switch)
# ─────────────────────────────────────────────────────────────────────────────
# Kept for diagnostic purposes only — no longer gates the call. Reads the
# per_window_class.rds produced by STEP_C01i_decompose if present.
find_signal1_segment <- function(per_window_class_col, global_class,
                                    window_starts, window_ends,
                                    min_consecutive = 2L) {
  n <- length(per_window_class_col)
  if (n < min_consecutive + 1L) return(NULL)
  diff <- ifelse(!is.na(per_window_class_col) &
                  per_window_class_col != global_class, 1L, 0L)
  if (sum(diff) < min_consecutive) return(NULL)
  rle_d <- rle(diff)
  max_len <- 0L; max_start <- NA_integer_; idx <- 1L
  for (i in seq_along(rle_d$values)) {
    if (rle_d$values[i] == 1L && rle_d$lengths[i] >= min_consecutive) {
      if (rle_d$lengths[i] > max_len) {
        max_len   <- rle_d$lengths[i]
        max_start <- idx
      }
    }
    idx <- idx + rle_d$lengths[i]
  }
  if (is.na(max_start)) return(NULL)
  list(
    mosaic_start_bp  = as.integer(window_starts[max_start]),
    mosaic_end_bp    = as.integer(window_ends[max_start + max_len - 1L]),
    mosaic_length_bp = as.integer(window_ends[max_start + max_len - 1L] -
                                     window_starts[max_start]),
    n_switch_windows = as.integer(max_len)
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Load the C01j regime memberships ONCE, subset per candidate in the loop
# ─────────────────────────────────────────────────────────────────────────────
message("[multi_recomb] loading regime memberships from ", opt$regime_memb)
regime_all <- tryCatch(fread(opt$regime_memb), error = function(e) {
  stop("[multi_recomb] failed to read regime memberships: ",
       conditionMessage(e))
})
if (!nrow(regime_all)) stop("[multi_recomb] regime_memb is empty")

# Normalise column names so derive_R_from_regime finds them
if (!"sample" %in% names(regime_all) && "sample_id" %in% names(regime_all)) {
  setnames(regime_all, "sample_id", "sample")
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: load the GHSL per-sample call for a candidate, if available
# ─────────────────────────────────────────────────────────────────────────────
load_ghsl_for_cid <- function(cid) {
  if (is.na(opt$ghsl_dir) || !nzchar(opt$ghsl_dir)) {
    return(list(g_dt = NULL, resolution = "no_annot"))
  }
  f <- file.path(opt$ghsl_dir, cid, "ghsl_per_sample.tsv")
  if (!file.exists(f)) f <- file.path(opt$ghsl_dir, paste0(cid, "_ghsl_per_sample.tsv"))
  if (!file.exists(f)) return(list(g_dt = NULL, resolution = "no_annot"))
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(list(g_dt = NULL, resolution = "no_annot"))
  # Expect columns: sample_id, ghsl_call, (optional) resolution_flag
  if (!"sample_id" %in% names(dt) || !"ghsl_call" %in% names(dt)) {
    return(list(g_dt = NULL, resolution = "no_annot"))
  }
  res <- if ("resolution" %in% names(dt)) dt$resolution[1] else "sufficient"
  list(g_dt = dt[, .(sample_id, ghsl_call)], resolution = as.character(res))
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: load GC summary for a candidate if present
# ─────────────────────────────────────────────────────────────────────────────
load_gc_for_cid <- function(cid) {
  if (is.na(opt$gc_dir) || !nzchar(opt$gc_dir)) return(NULL)
  f <- file.path(opt$gc_dir, cid, "gene_conversion_tracts.json")
  if (!file.exists(f)) return(NULL)
  # Lazy-load jsonlite to avoid a hard dep when gc_dir is unused
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  blk <- tryCatch(jsonlite::fromJSON(f, simplifyVector = TRUE),
                  error = function(e) NULL)
  if (is.null(blk) || is.null(blk$per_sample_summary)) return(NULL)
  s <- as.data.table(blk$per_sample_summary)
  if (!"sample_id" %in% names(s) || !"n_tracts" %in% names(s)) return(NULL)
  s[, .(sample_id, n_tracts = as.integer(n_tracts))]
}

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
cand_dt <- fread(opt$candidates)
if (!"candidate_id" %in% names(cand_dt)) {
  cand_dt[, candidate_id := paste0(chrom, "_", interval_id)]
}
if (!"tier" %in% names(cand_dt)) cand_dt[, tier := 2L]
if (!"start_bp" %in% names(cand_dt))
  cand_dt[, start_bp := as.integer(start_mb * 1e6)]
if (!"end_bp" %in% names(cand_dt))
  cand_dt[, end_bp := as.integer(end_mb * 1e6)]
cand_dt <- cand_dt[tier <= opt$tier_max]

cat("[multi_recomb] processing ", nrow(cand_dt), " candidates at tier<=",
    opt$tier_max, "\n", sep = "")

for (ci in seq_len(nrow(cand_dt))) {
  cd  <- cand_dt[ci]
  cid <- cd$candidate_id
  chr <- cd$chrom
  s   <- as.integer(cd$start_bp)
  e   <- as.integer(cd$end_bp)
  cat("\n[multi_recomb] ─── ", cid, " (", chr, ":",
      round(s/1e6, 2), "-", round(e/1e6, 2), " Mb) ───\n", sep = "")

  # ── Restrict regime to this chrom/candidate interval ──
  rm_cid <- regime_all[chrom == chr]
  if (!nrow(rm_cid)) {
    cat("[multi_recomb]   SKIP: no regime windows on ", chr, "\n", sep = "")
    next
  }

  # ── R: DAG-based recombinant signal ──
  r_res <- derive_R_from_regime(
    regime_memb           = rm_cid,
    interval              = list(start_bp = s, end_bp = e),
    min_deviation_fraction= opt$min_dev_frac,
    min_deviation_bp      = opt$min_dev_bp
  )
  r_dt <- r_res$per_sample

  # ── G: GHSL ──
  ghsl <- load_ghsl_for_cid(cid)

  # ── C: GC tracts (annotation only; does not gate) ──
  gc_summary <- load_gc_for_cid(cid)

  # ── Combine (gate: R AND G) ──
  combo <- combine_cohort_recomb(
    r_dt         = r_dt,
    g_dt         = ghsl$g_dt,
    g_resolution = ghsl$resolution,
    gc_summary   = gc_summary,
    all_samples  = r_dt$sample_id
  )

  # ── Event-class split within RECOMBINANT (chat-12) ──
  # Join back to r_dt for per-sample longest_deviation_bp
  if (nrow(combo)) {
    combo <- merge(
      combo,
      r_dt[, .(sample_id, longest_deviation_bp, deviation_fraction,
               dominant_regime, n_edges, terminal_matches_start)],
      by = "sample_id", all.x = TRUE, sort = FALSE
    )
    combo[, event_class := fifelse(
      recomb_status == "RECOMBINANT" &
        !is.na(longest_deviation_bp) & longest_deviation_bp >= opt$min_dco_bp,
      "double_crossover",
      fifelse(recomb_status == "RECOMBINANT" & (n_gc_tracts %||% 0L) > 0L,
              "gene_conversion_embedded",
              fifelse(recomb_status == "RECOMBINANT", "single_crossover_or_ambiguous",
                      NA_character_))
    )]
    combo[, recomb_subgroup := fifelse(
      recomb_status == "RECOMBINANT" & event_class == "double_crossover",
      "RECOMBINANT_DCO",
      fifelse(recomb_status == "RECOMBINANT" & event_class == "gene_conversion_embedded",
              "RECOMBINANT_GC",
              fifelse(recomb_status == "RECOMBINANT", "RECOMBINANT",
                      NA_character_))
    )]
  } else {
    combo[, `:=`(event_class = character(), recomb_subgroup = character(),
                 longest_deviation_bp = integer(), deviation_fraction = numeric(),
                 dominant_regime = character(), n_edges = integer(),
                 terminal_matches_start = logical())]
  }

  # ── Posterior via real cheat24 (if present) — Finding V ──
  # We no longer fabricate a number inline. If cheat24 is present, call it
  # per RECOMBINANT sample; otherwise leave posterior = NA.
  combo[, posterior := NA_real_]
  if (.has_cheat24 && exists("classify_recombinant_event", mode = "function")) {
    rec_idx <- which(combo$recomb_status == "RECOMBINANT")
    for (k in rec_idx) {
      # cheat24's classify_recombinant_event contract (see
      # flashlight_v2/cheats/cheat24_recombinant_prior.R) takes
      # recomb_position, left_bp, right_bp, mosaic_length_bp. We don't
      # have a meaningful single-point switchpoint from the DAG, so
      # pass the interval midpoint and the longest deviation as the
      # mosaic length.
      res <- tryCatch(classify_recombinant_event(
        recomb_position = as.integer((s + e) / 2L),
        left_bp         = s,
        right_bp        = e,
        mosaic_length   = combo$longest_deviation_bp[k]
      ), error = function(e) NULL)
      if (!is.null(res)) {
        combo$posterior[k] <- as.numeric(res$posterior_probability %||% NA_real_)
      }
    }
  }

  # ── Tier-4 diagnostic: old Signal 1 (per-window PCA class switch) ──
  # Not a gate; recorded per-sample for auditing only.
  s1_rows <- list()
  pw_file <- file.path(opt$decomp_dir, cid, "per_window_class.rds")
  if (file.exists(pw_file)) {
    pw <- tryCatch(readRDS(pw_file), error = function(e) NULL)
    if (!is.null(pw)) {
      pw_class     <- pw$per_window_class
      global_class <- pw$class_assignment_global
      win_starts   <- pw$window_starts
      win_ends     <- pw$window_ends
      for (sid in pw$sample_names) {
        gc_class <- unname(global_class[sid])
        if (is.na(gc_class)) next
        s1 <- find_signal1_segment(pw_class[, sid], gc_class,
                                   win_starts, win_ends,
                                   min_consecutive = 2L)
        if (!is.null(s1)) {
          s1_rows[[length(s1_rows) + 1L]] <- data.table(
            sample_id        = sid,
            s1_mosaic_start_bp  = s1$mosaic_start_bp,
            s1_mosaic_end_bp    = s1$mosaic_end_bp,
            s1_mosaic_length_bp = s1$mosaic_length_bp,
            s1_n_switch_windows = s1$n_switch_windows
          )
        }
      }
    }
  }
  s1_dt <- if (length(s1_rows)) rbindlist(s1_rows) else
           data.table(sample_id = character())
  if (nrow(s1_dt) && nrow(combo)) {
    combo <- merge(combo, s1_dt, by = "sample_id", all.x = TRUE, sort = FALSE)
  }

  # ── Counts for logging ──
  n_recomb <- sum(combo$recomb_status == "RECOMBINANT", na.rm = TRUE)
  n_gc     <- sum(combo$recomb_subgroup == "RECOMBINANT_GC", na.rm = TRUE)
  n_dco    <- sum(combo$recomb_subgroup == "RECOMBINANT_DCO", na.rm = TRUE)
  n_disp   <- sum(combo$recomb_status == "recomb_disputed", na.rm = TRUE)
  n_gonly  <- sum(combo$recomb_status == "recomb_ghsl_only", na.rm = TRUE)
  cat("[multi_recomb]   ", n_recomb, " recombinants (GC=", n_gc,
      " DCO=", n_dco, ") | disputed=", n_disp, " ghsl_only=", n_gonly, "\n",
      sep = "")

  # ── Write recombinant_map block ──
  # list-column gc_tract_widths_bp from combine_cohort_recomb can't be
  # JSON-serialised as-is by all writers; drop before writing.
  recombinants_list <- lapply(seq_len(nrow(combo)), function(k) {
    row <- as.list(combo[k])
    row$gc_tract_widths_bp <- NULL
    row
  })

  block_data <- list(
    candidate_id      = cid,
    chrom             = chr,
    start_bp          = s,
    end_bp            = e,
    n_samples         = as.integer(nrow(combo)),
    n_recombinants    = as.integer(n_recomb),
    n_recombinant_gc  = as.integer(n_gc),
    n_recombinant_dco = as.integer(n_dco),
    n_disputed        = as.integer(n_disp),
    n_ghsl_only       = as.integer(n_gonly),
    cheat24_version   = if (.has_cheat24) "cheat24_recombinant_prior.R"
                        else "cheat24_unavailable_posterior_na",
    combination_rule  = "lib_recomb_combination.R (DAG R + G gate, chat-12)",
    gate_params       = list(
      min_deviation_fraction = opt$min_dev_frac,
      min_deviation_bp       = as.integer(opt$min_dev_bp),
      min_dco_bp             = as.integer(opt$min_dco_bp)
    ),
    recombinants      = recombinants_list
  )

  cand_outdir <- file.path(opt$outdir, cid)
  dir.create(cand_outdir, recursive = TRUE, showWarnings = FALSE)

  write_block_safe(
    reg            = reg,
    candidate_id   = cid,
    block_type     = "recombinant_map",
    data           = block_data,
    source_script  = "STEP_C01i_b_multi_recomb.R",
    outdir_fallback= cand_outdir
  )

  # ──────────────────────────────────────────────────────────────────────
  # Chat-16: register segment candidates + segment karyotype groups.
  #
  # Only runs if at least one recombinant sample has a tract (CO or DCO).
  # Creates child candidates in interval_registry with parent_id=cid, one
  # per distinct segment boundary set observed across all recombinant
  # samples. For each segment:
  #   - segment candidate inv_<cid>_seg_<LABEL> registered via add_candidate
  #   - karyotype groups inv_<cid>_seg_<LABEL>_<KARYO> registered via
  #     add_group for each sample's LOCAL arrangement in that segment
  #
  # This is what lets reg$compute$scan_pairwise() and scan_region_stat()
  # pick up the correct groups per window during a multi-window scan.
  # See DATABASE_DESIGN.md § "Segmented karyotypes for recombinant
  # candidates" for the full convention.
  #
  # Safe-by-design: if add_candidate / add_group are unavailable in reg
  # (e.g. called with a stub registry), the whole block is skipped with
  # a message — does not affect the recombinant_map block above.
  # ──────────────────────────────────────────────────────────────────────
  if (n_recomb > 0 && !is.null(reg$intervals) &&
      !is.null(reg$intervals$add_candidate) &&
      !is.null(reg$samples) && !is.null(reg$samples$add_group)) {
    tryCatch({
      # Collect every breakpoint across all recombinant samples.
      # combo has columns: sample_id, R_type (RECOMBINANT/RECOMBINANT_GC/
      # RECOMBINANT_DCO), regime_tract_start_bp, regime_tract_end_bp,
      # gc_tract_starts_bp, gc_tract_ends_bp (space-separated lists).
      rec_rows <- combo[combo$R_type %in%
                         c("RECOMBINANT", "RECOMBINANT_GC", "RECOMBINANT_DCO")]
      if (nrow(rec_rows) > 0) {
        # Union of all breakpoints within (s, e), deduplicated
        bps <- integer(0)
        for (i in seq_len(nrow(rec_rows))) {
          tss <- suppressWarnings(as.integer(rec_rows$regime_tract_start_bp[i]))
          tse <- suppressWarnings(as.integer(rec_rows$regime_tract_end_bp[i]))
          if (!is.na(tss) && tss > s && tss < e) bps <- c(bps, tss)
          if (!is.na(tse) && tse > s && tse < e) bps <- c(bps, tse)
          # GC inner tracts (if present)
          gs <- rec_rows$gc_tract_starts_bp[i]
          ge <- rec_rows$gc_tract_ends_bp[i]
          if (!is.null(gs) && !is.na(gs) && nzchar(as.character(gs))) {
            gss <- suppressWarnings(
              as.integer(strsplit(as.character(gs), "[ ,]+")[[1]]))
            gee <- suppressWarnings(
              as.integer(strsplit(as.character(ge), "[ ,]+")[[1]]))
            for (gj in seq_along(gss)) {
              if (!is.na(gss[gj]) && gss[gj] > s && gss[gj] < e) bps <- c(bps, gss[gj])
              if (!is.na(gee[gj]) && gee[gj] > s && gee[gj] < e) bps <- c(bps, gee[gj])
            }
          }
        }
        bps <- sort(unique(bps))

        # ──── Chat-16 SANITY GUARD — refuse noise-driven segmentation ────
        # Biology bound: even a triple-crossover gives only 4 segments
        # per candidate. More than that and something is very likely
        # wrong with the recombinant detector. Also refuse segments
        # smaller than 10 kb — below that the FST / ancestry signal is
        # dominated by noise at 9× coverage and the segment is
        # scientifically meaningless.
        MAX_N_SEGMENTS_SANE <- 6L         # >4 segments = suspicious
        MAX_N_SEGMENTS_HARD <- 10L        # >10 = refuse outright
        MIN_SEGMENT_BP      <- 10000L     # <10 kb = below signal floor
        candidate_span_bp   <- as.integer(e) - as.integer(s)

        n_proposed_segs <- length(bps) + 1L
        seg_sizes <- diff(c(as.integer(s), bps, as.integer(e)))
        n_tiny <- sum(seg_sizes < MIN_SEGMENT_BP)

        if (n_proposed_segs > MAX_N_SEGMENTS_HARD) {
          warning(sprintf(paste0(
            "[multi_recomb] %s: REFUSING segment registration — ",
            "%d segments proposed (> hard cap %d). This is almost ",
            "certainly noise in the recombinant detector. ",
            "recombinant_map block was written but NO segment candidates ",
            "or segment groups will be registered. Inspect the ",
            "breakpoint distribution before re-enabling."),
            cid, n_proposed_segs, MAX_N_SEGMENTS_HARD))
          return()   # exits the tryCatch body; parent function continues
        }
        if (n_proposed_segs > MAX_N_SEGMENTS_SANE || n_tiny > 0L) {
          warning(sprintf(paste0(
            "[multi_recomb] %s: SUSPICIOUS segmentation — %d segments ",
            "proposed (biology expects ≤4), %d tiny segments (<%d bp). ",
            "Registering anyway but flag for manual review. ",
            "Candidate span: %.2f Mb. Segment sizes (kb): %s"),
            cid, n_proposed_segs, n_tiny, MIN_SEGMENT_BP,
            candidate_span_bp / 1e6,
            paste(sprintf("%.1f", seg_sizes / 1000), collapse = ",")))
        }

        # Segment boundaries = [s, bps..., e]
        seg_bounds <- c(as.integer(s), bps, as.integer(e))
        if (length(seg_bounds) >= 3L) {
          # ≥ 2 segments exist
          n_segs <- length(seg_bounds) - 1L
          cat(sprintf("[multi_recomb] %s: registering %d segments from %d breakpoints (span %.2f Mb)\n",
                      cid, n_segs, length(bps), candidate_span_bp / 1e6))
          seg_ids <- character(n_segs)
          for (si in seq_len(n_segs)) {
            seg_s <- seg_bounds[si]
            seg_e <- seg_bounds[si + 1L]
            # Label: _seg_L for first, _seg_R for last, _seg_<N> for middles.
            # For 3 segments with 1 DCO: usually L, DCO1, R. But we can't
            # infer DCO vs CO purely from positions; leave as numbered
            # and let downstream scripts annotate via block_data lookups.
            seg_label <- if (n_segs == 3L && si == 2L) "DCO1"
                         else if (si == 1L) "L"
                         else if (si == n_segs) "R"
                         else sprintf("MID%d", si - 1L)
            seg_cid <- sprintf("%s_seg_%s", cid, seg_label)
            seg_ids[si] <- seg_cid
            # Register segment candidate (idempotent via add_candidate's
            # "skips if exists" behavior — chat-11 interval API contract).
            tryCatch(reg$intervals$add_candidate(
                       candidate_id = seg_cid,
                       chrom        = chr,
                       start_bp     = seg_s,
                       end_bp       = seg_e,
                       scale        = if (seg_label == "DCO1") "seg_dco" else "seg",
                       parent_id    = cid),
                     error = function(e)
                       message("[multi_recomb] add_candidate ", seg_cid,
                               " failed: ", conditionMessage(e)))
          }

          # Now build per-segment karyotype groups.
          # For each sample, determine its local karyotype in each segment:
          #   - If sample is NOT a recombinant (R_type in HOM_REF/HET/HOM_INV):
          #     local karyotype == cohort karyotype in ALL segments
          #   - If sample IS a recombinant: invert the karyotype inside
          #     the tract boundaries (seg_DCO for DCO carriers is HOM_REF
          #     if parent is HOM_INV, and vice versa), keep the parent
          #     karyotype in flanking segments
          cohort_karyo <- combo[, .(sample_id, R_type,
                                     regime_tract_start_bp,
                                     regime_tract_end_bp)]
          # Map R_type → parent karyotype label. Recombinant rows use
          # their parent's most-common arrangement; here we assume
          # HOM_INV as parent for RECOMBINANT_* (chat-12 convention).
          # Non-recombinant rows pass through their R_type directly.
          .infer_parent_karyo <- function(r_type) {
            if (r_type %in% c("RECOMBINANT","RECOMBINANT_GC","RECOMBINANT_DCO"))
              return("HOM_INV")
            r_type
          }
          .invert_karyo <- function(k) {
            if (k == "HOM_INV") "HOM_REF" else if (k == "HOM_REF") "HOM_INV" else k
          }
          # Build membership map: seg_cid -> list(karyotype -> sample_ids)
          seg_members <- setNames(
            lapply(seq_len(n_segs), function(x) list()),
            seg_ids)
          for (i in seq_len(nrow(cohort_karyo))) {
            sid <- cohort_karyo$sample_id[i]
            rtype <- cohort_karyo$R_type[i]
            tss <- suppressWarnings(as.integer(cohort_karyo$regime_tract_start_bp[i]))
            tse <- suppressWarnings(as.integer(cohort_karyo$regime_tract_end_bp[i]))
            parent_k <- .infer_parent_karyo(rtype)
            for (si in seq_len(n_segs)) {
              seg_s <- seg_bounds[si]; seg_e <- seg_bounds[si + 1L]
              seg_mid <- (seg_s + seg_e) %/% 2L
              # Determine local karyotype for this sample in this segment
              local_k <- if (rtype %in% c("RECOMBINANT","RECOMBINANT_GC","RECOMBINANT_DCO")) {
                # Is the segment midpoint inside this sample's recombinant tract?
                inside <- !is.na(tss) && !is.na(tse) &&
                           seg_mid >= tss && seg_mid <= tse
                if (inside) .invert_karyo(parent_k) else parent_k
              } else {
                parent_k  # non-recombinant: same everywhere
              }
              if (is.null(seg_members[[seg_ids[si]]][[local_k]])) {
                seg_members[[seg_ids[si]]][[local_k]] <- character(0)
              }
              seg_members[[seg_ids[si]]][[local_k]] <-
                c(seg_members[[seg_ids[si]]][[local_k]], sid)
            }
          }
          # Register the groups (skip empty ones; add_group is idempotent)
          for (si in seq_len(n_segs)) {
            seg_cid <- seg_ids[si]
            for (karyo in names(seg_members[[seg_cid]])) {
              members <- unique(seg_members[[seg_cid]][[karyo]])
              if (!length(members)) next
              gid <- sprintf("inv_%s_%s", seg_cid, karyo)
              tryCatch(reg$samples$add_group(
                         group_id    = gid,
                         sample_ids  = members,
                         chrom       = chr,
                         inv_id      = seg_cid,
                         description = sprintf(
                           "Segment karyotype from multi_recomb (%s, parent %s)",
                           seg_cid, cid)),
                       error = function(e)
                         message("[multi_recomb] add_group ", gid,
                                 " failed: ", conditionMessage(e)))
            }
          }
          cat(sprintf("[multi_recomb] %s: %d segment candidates + groups registered\n",
                      cid, n_segs))
        }
      }
    }, error = function(e) {
      message("[multi_recomb] segment registration for ", cid,
              " failed (non-fatal): ", conditionMessage(e))
    })
  }

  # ── Write the DAG evidence block (regime_sample_dag schema) ──
  # Trim per_sample to schema fields before writing
  ps <- r_res$per_sample[, .(
    sample_id, n_windows, n_distinct_nodes, n_edges,
    dominant_regime, dominant_run_windows,
    deviation_fraction, longest_deviation_windows, longest_deviation_bp,
    terminal_matches_start, R_fired,
    dag_nodes_compact
    # dag_node_weights intentionally omitted: list-columns don't round-trip
    # through the registry's JSON writer. The compact string is sufficient
    # for the diagnostic plot; the weights can be recomputed on-demand.
  )]
  dag_block <- list(
    candidate_id = cid,
    per_sample   = ps,
    cohort       = r_res$cohort,
    params       = r_res$params
  )
  write_block_safe(
    reg             = reg,
    candidate_id    = cid,
    block_type      = "regime_sample_dag",
    data            = dag_block,
    source_script   = "STEP_C01i_b_multi_recomb.R",
    outdir_fallback = cand_outdir
  )
}

cat("\n[multi_recomb] done\n")
