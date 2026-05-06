#!/usr/bin/env Rscript
# =============================================================================
# STEP_RD_TE_density_full.R — Comprehensive TE density layer for the Inversion Atlas
# =============================================================================
#
# Reads BOTH the EDTA TEanno.gff3 (homology + structural, all 26 classes) and
# the EDTA intact.gff3 (structural-only, with parent-child TSD records), and
# emits a single per-chromosome repeat_density_v2 JSON containing:
#
#   ── Per-class density arrays (78 total) ──────────────────────────────────
#   For each TE class C ∈ {LTR_retrotransposon, Tc1_Mariner_TIR_transposon,
#                          Mutator_TIR_transposon, ...} (26 classes):
#     C            : fraction of bp covered per 5-kb window (all elements)
#     C__young     : same, filtered to identity >= 0.95 (recently inserted)
#     C__old       : same, filtered to identity <  0.95 (ancient remnants)
#
#   ── Aggregate arrays (5 total) ───────────────────────────────────────────
#   all_TE              : sum of fraction across all classes (saturates → 1)
#   young_TE_all        : sum of young fractions
#   old_TE_all          : sum of old fractions
#   insertion_count     : raw count of TE records per window (any class)
#   intact_element_count: count of `repeat_region` parents from intact.gff3
#
#   ── TSD arrays (2 total) ─────────────────────────────────────────────────
#   target_site_duplication        : explicit-coord TSDs from intact.gff3
#                                    (the same layer as STEP_RD_TSD_density.R
#                                    output — 4k Gar / 1.8k Mac records)
#   target_site_duplication_all    : explicit + parsed-from-attribute TSDs
#                                    of TIR transposons (TSD=motif_motif_pct
#                                    attribute → infer position from element
#                                    flanks). Denser layer.
#
# Total: 85 density arrays per chromosome JSON. File size ~6-8 MB per chrom.
#
# Usage:
#   Rscript STEP_RD_TE_density_full.R \
#     --teanno_gff3 /path/to/EDTA.TEanno.gff3 \
#     --intact_gff3 /path/to/EDTA.intact.gff3 \
#     --fai         /path/to/genome.fa.fai \
#     --window_bp   5000 \
#     --step_bp     5000 \
#     --out_dir     /out_dir \
#     [--chrom      C_gar_LG01]   # if omitted, processes every chromosome
#     [--species    "Clarias gariepinus"]
#     [--identity_young_threshold 0.95]
#
# Runtime: ~3-5 minutes per species on the login node (bigger input than
# intact-only). Memory: ~2-3 GB peak for TEanno read.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# CLI parser
# -----------------------------------------------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    teanno_gff3 = NA_character_,
    intact_gff3 = NA_character_,
    fai         = NA_character_,
    window_bp   = 5000L,
    step_bp     = 5000L,
    out_dir     = ".",
    chrom       = NA_character_,
    species     = "Clarias gariepinus",
    identity_young_threshold = 0.95
  )
  i <- 1
  while (i <= length(args)) {
    flag <- args[i]; val <- if (i+1 <= length(args)) args[i+1] else NA
    name <- sub("^--", "", flag)
    if (is.null(defaults[[name]])) stop("unknown flag: ", flag)
    if (is.numeric(defaults[[name]]) && !is.na(val)) {
      val <- if (is.integer(defaults[[name]])) as.integer(val) else as.numeric(val)
    }
    defaults[[name]] <- val
    i <- i + 2
  }
  for (req in c("teanno_gff3", "intact_gff3", "fai")) {
    if (is.na(defaults[[req]])) stop("--", req, " is required")
  }
  defaults
}

opts <- parse_args()
cat("[STEP_RD_TE_FULL] starting\n")
cat("  teanno_gff3 :", opts$teanno_gff3, "\n")
cat("  intact_gff3 :", opts$intact_gff3, "\n")
cat("  fai         :", opts$fai, "\n")
cat("  window_bp   :", opts$window_bp, "\n")
cat("  step_bp     :", opts$step_bp, "\n")
cat("  out_dir     :", opts$out_dir, "\n")
cat("  species     :", opts$species, "\n")
cat("  young_thr   :", opts$identity_young_threshold, "\n")

# -----------------------------------------------------------------------------
# Load chromosome sizes
# -----------------------------------------------------------------------------
fai <- fread(opts$fai, header = FALSE, sep = "\t",
             select = 1:2, col.names = c("chrom", "size_bp"))
if (!is.na(opts$chrom)) {
  fai <- fai[chrom == opts$chrom]
  if (nrow(fai) == 0) stop("--chrom not found in .fai: ", opts$chrom)
}
cat("[STEP_RD_TE_FULL] chromosomes to process:", nrow(fai), "\n")

# -----------------------------------------------------------------------------
# Read TEanno.gff3 — full annotation, ~289 MB Gar / 271 MB Mac
# We need columns: 1 (seqid), 3 (type), 4 (start), 5 (end), 9 (attributes)
# -----------------------------------------------------------------------------
cat("[STEP_RD_TE_FULL] reading TEanno.gff3 ...\n")
ta <- fread(opts$teanno_gff3, sep = "\t", header = FALSE,
            comment.char = "#",
            select = c(1, 3, 4, 5, 9),
            col.names = c("chrom", "type", "start_bp", "end_bp", "attrs"))
cat("[STEP_RD_TE_FULL] TEanno records:", nrow(ta), "across", uniqueN(ta$chrom), "chromosomes\n")

# Extract identity from attribute strings.
# Attribute format: "ID=...;Name=...;classification=...;sequence_ontology=...;identity=0.993;method=homology"
# We use a fast regex extraction.
ta[, identity := suppressWarnings(as.numeric(
  sub(".*?identity=([0-9.]+).*", "\\1", attrs, perl = TRUE)
))]
# If the regex didn't match (no `identity=` field), as.numeric on the original
# returns NA. Sanity:
cat("[STEP_RD_TE_FULL] records with identity field:",
    sum(!is.na(ta$identity)), "/", nrow(ta), "\n")

# -----------------------------------------------------------------------------
# Read intact.gff3 — for explicit-coord TSDs and intact-element counts
# -----------------------------------------------------------------------------
cat("[STEP_RD_TE_FULL] reading intact.gff3 ...\n")
ig <- fread(opts$intact_gff3, sep = "\t", header = FALSE,
            comment.char = "#",
            select = c(1, 3, 4, 5),
            col.names = c("chrom", "type", "start_bp", "end_bp"))
cat("[STEP_RD_TE_FULL] intact records:", nrow(ig), "\n")

# -----------------------------------------------------------------------------
# TE classes — the 26 (or however many) sequence-ontology types in TEanno.
# We treat anything that's not a sub-feature (target_site_duplication,
# long_terminal_repeat, repeat_region) as a "real" TE class.
# -----------------------------------------------------------------------------
SUBFEATURE_TYPES <- c(
  "target_site_duplication",  # TSD (sub-feature in intact)
  "long_terminal_repeat",     # LTR halves of LTRRT (sub-feature)
  "repeat_region"             # parent of intact LTR elements (sub-feature)
)
te_classes <- sort(unique(ta$type))
te_classes <- te_classes[!te_classes %in% SUBFEATURE_TYPES]
cat("[STEP_RD_TE_FULL] TE classes detected:", length(te_classes), "\n")
for (c in te_classes) cat("  -", c, "\n")

# -----------------------------------------------------------------------------
# Helper: compute per-window bp-coverage fraction.
# Given intervals (start_bp, end_bp) on one chromosome and a window grid,
# returns the FRACTION of each window's bp that's covered by any interval.
#
# Vectorized via data.table foverlaps. For ~17k records × 78 class-strats per
# chrom this needs to be sub-second; the naïve double-loop would be hours.
# -----------------------------------------------------------------------------
compute_coverage_fraction <- function(starts, ends, window_bp, step_bp, n_windows) {
  if (length(starts) == 0) return(numeric(n_windows))

  # Build interval data.table (1-based, inclusive end — foverlaps convention).
  # GFF coords are 1-based inclusive already. Window grid is 0-based start,
  # exclusive end, so we adjust to 1-based inclusive for the overlap join.
  ivl <- data.table(start = as.integer(starts), end = as.integer(ends))
  ivl <- ivl[end >= start]
  if (nrow(ivl) == 0) return(numeric(n_windows))
  setkey(ivl, start, end)

  win_starts <- (seq_len(n_windows) - 1L) * step_bp + 1L  # 1-based
  win_ends   <- pmin(win_starts + window_bp - 1L,
                     win_starts[length(win_starts)] + window_bp - 1L)
  win <- data.table(widx = seq_len(n_windows),
                    start = win_starts, end = win_ends)
  setkey(win, start, end)

  # foverlaps returns one row per (interval × window) overlap pair.
  # Compute the bp of each overlap, then sum per widx.
  ov <- foverlaps(ivl, win, type = "any", nomatch = NULL)
  if (nrow(ov) == 0) return(numeric(n_windows))
  # Intersection bp = min(end, i.end) - max(start, i.start) + 1
  ov[, ov_bp := pmin(end, i.end) - pmax(start, i.start) + 1L]
  ov[ov_bp < 0L, ov_bp := 0L]
  cov_per_window <- ov[, .(bp = sum(ov_bp)), by = widx]

  cov_bp <- numeric(n_windows)
  cov_bp[cov_per_window$widx] <- cov_per_window$bp
  pmin(cov_bp / window_bp, 1.0)
}

# Helper: count records per window (each record = 1 event at midpoint)
count_per_window <- function(starts, ends, step_bp, n_windows) {
  if (length(starts) == 0) return(integer(n_windows))
  mid <- (starts + ends) %/% 2L
  widx <- (mid %/% step_bp) + 1L
  widx <- widx[widx >= 1L & widx <= n_windows]
  tabulate(widx, nbins = n_windows)
}

# Build window grid for one chromosome
build_window_grid <- function(chrom_size_bp, window_bp, step_bp) {
  starts <- seq(0L, chrom_size_bp - 1L, by = step_bp)
  ends   <- pmin(starts + window_bp, as.integer(chrom_size_bp))
  centers_mb <- (starts + ends) / 2 / 1e6
  list(
    n_windows         = length(starts),
    window_start_bp   = as.integer(starts),
    window_end_bp     = as.integer(ends),
    window_centers_mb = centers_mb
  )
}

# -----------------------------------------------------------------------------
# Per-chromosome JSON build
# -----------------------------------------------------------------------------
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite not installed; install.packages('jsonlite') in your assembly env")
}

YOUNG_THR <- as.numeric(opts$identity_young_threshold)

for (i in seq_len(nrow(fai))) {
  this_chrom <- fai$chrom[i]
  size_bp    <- fai$size_bp[i]
  if (size_bp < opts$window_bp) {
    cat("  skip", this_chrom, "(size", size_bp, "<", opts$window_bp, ")\n")
    next
  }

  grid <- build_window_grid(size_bp, opts$window_bp, opts$step_bp)
  ta_chr <- ta[chrom == this_chrom]
  ig_chr <- ig[chrom == this_chrom]

  by_class <- list()
  classes_emitted <- character(0)

  # (1) Per-class densities (all + young + old) ─────────────────────────────
  for (cls in te_classes) {
    rows_all <- ta_chr[type == cls]
    if (nrow(rows_all) == 0) next   # skip classes absent on this chromosome

    dens_all <- compute_coverage_fraction(
      rows_all$start_bp, rows_all$end_bp,
      opts$window_bp, opts$step_bp, grid$n_windows
    )

    # Young / old split based on identity
    rows_with_id <- rows_all[!is.na(identity)]
    rows_young   <- rows_with_id[identity >= YOUNG_THR]
    rows_old     <- rows_with_id[identity <  YOUNG_THR]

    dens_young <- compute_coverage_fraction(
      rows_young$start_bp, rows_young$end_bp,
      opts$window_bp, opts$step_bp, grid$n_windows
    )
    dens_old <- compute_coverage_fraction(
      rows_old$start_bp, rows_old$end_bp,
      opts$window_bp, opts$step_bp, grid$n_windows
    )

    by_class[[cls]] <- list(
      densities    = dens_all,
      max_density  = max(dens_all),
      n_records    = nrow(rows_all)
    )
    classes_emitted <- c(classes_emitted, cls)

    if (nrow(rows_young) > 0) {
      young_key <- paste0(cls, "__young")
      by_class[[young_key]] <- list(
        densities    = dens_young,
        max_density  = max(dens_young),
        n_records    = nrow(rows_young)
      )
      classes_emitted <- c(classes_emitted, young_key)
    }
    if (nrow(rows_old) > 0) {
      old_key <- paste0(cls, "__old")
      by_class[[old_key]] <- list(
        densities    = dens_old,
        max_density  = max(dens_old),
        n_records    = nrow(rows_old)
      )
      classes_emitted <- c(classes_emitted, old_key)
    }
  }

  # (2) Aggregate densities (all_TE, young_TE_all, old_TE_all) ──────────────
  # Sum across all classes (NOT sub-features — they're already excluded).
  all_te_rows <- ta_chr[!type %in% SUBFEATURE_TYPES]
  young_rows  <- all_te_rows[!is.na(identity) & identity >= YOUNG_THR]
  old_rows    <- all_te_rows[!is.na(identity) & identity <  YOUNG_THR]

  by_class[["all_TE"]] <- list(
    densities   = compute_coverage_fraction(all_te_rows$start_bp, all_te_rows$end_bp,
                                            opts$window_bp, opts$step_bp, grid$n_windows),
    max_density = NA, n_records = nrow(all_te_rows)
  )
  by_class[["all_TE"]]$max_density <- max(by_class[["all_TE"]]$densities)
  classes_emitted <- c(classes_emitted, "all_TE")

  if (nrow(young_rows) > 0) {
    by_class[["young_TE_all"]] <- list(
      densities   = compute_coverage_fraction(young_rows$start_bp, young_rows$end_bp,
                                              opts$window_bp, opts$step_bp, grid$n_windows),
      max_density = NA, n_records = nrow(young_rows)
    )
    by_class[["young_TE_all"]]$max_density <- max(by_class[["young_TE_all"]]$densities)
    classes_emitted <- c(classes_emitted, "young_TE_all")
  }
  if (nrow(old_rows) > 0) {
    by_class[["old_TE_all"]] <- list(
      densities   = compute_coverage_fraction(old_rows$start_bp, old_rows$end_bp,
                                              opts$window_bp, opts$step_bp, grid$n_windows),
      max_density = NA, n_records = nrow(old_rows)
    )
    by_class[["old_TE_all"]]$max_density <- max(by_class[["old_TE_all"]]$densities)
    classes_emitted <- c(classes_emitted, "old_TE_all")
  }

  # (3) Insertion count (raw count of TE records per window) ────────────────
  ic <- count_per_window(all_te_rows$start_bp, all_te_rows$end_bp,
                         opts$step_bp, grid$n_windows)
  by_class[["insertion_count"]] <- list(
    densities   = as.numeric(ic),
    max_density = if (length(ic) > 0) max(ic) else 0,
    n_records   = nrow(all_te_rows)
  )
  classes_emitted <- c(classes_emitted, "insertion_count")

  # (4) TSD layers ──────────────────────────────────────────────────────────
  # 4a. target_site_duplication: explicit coords from intact.gff3 (turn-1 layer)
  tsd_intact <- ig_chr[type == "target_site_duplication"]
  by_class[["target_site_duplication"]] <- list(
    densities   = as.numeric(count_per_window(tsd_intact$start_bp, tsd_intact$end_bp,
                                              opts$step_bp, grid$n_windows)),
    max_density = NA, n_records = nrow(tsd_intact)
  )
  by_class[["target_site_duplication"]]$max_density <- if (grid$n_windows > 0)
    max(by_class[["target_site_duplication"]]$densities) else 0
  classes_emitted <- c(classes_emitted, "target_site_duplication")

  # 4b. target_site_duplication_all: explicit + TIR-attribute-derived TSDs.
  # TIR records in TEanno.gff3 carry TSD as an attribute string like
  # "TSD=GTGTACTG_GTGTACTG_100.0" — we treat each such record as evidence of
  # one TSD pair flanking the element. We can't recover exact TSD coordinates
  # without re-parsing the FASTA, but we use the element's start_bp - tsd_len
  # and end_bp + tsd_len as approximations. The TSD motif length comes from
  # the attribute (length of motif1 before the underscore).
  has_tsd_attr <- grepl(";TSD=[^;]+;", ta_chr$attrs, perl = TRUE) |
                  grepl(";TSD=[^;]+$", ta_chr$attrs, perl = TRUE)
  tir_with_tsd <- ta_chr[has_tsd_attr]
  # Each TIR-with-TSD row contributes 2 events (left flank + right flank).
  if (nrow(tir_with_tsd) > 0) {
    # Extract motif length: TSD=GTGTACTG_GTGTACTG_100.0 → motif "GTGTACTG" → 8 bp
    motif1 <- sub(".*;TSD=([^_]+)_.*", "\\1", tir_with_tsd$attrs, perl = TRUE)
    tsd_len <- nchar(motif1)
    tsd_len[tsd_len > 30 | tsd_len < 2] <- 8L   # sanity-clamp on bad parses
    # Left flank ≈ element start - tsd_len; right flank ≈ element end + tsd_len
    left_pos  <- pmax(tir_with_tsd$start_bp - tsd_len, 1L)
    right_pos <- tir_with_tsd$end_bp + tsd_len
    all_tsd_pos <- c(tsd_intact$start_bp, left_pos, right_pos)
    all_tsd_pos2 <- c(tsd_intact$end_bp, left_pos, right_pos)
  } else {
    all_tsd_pos  <- tsd_intact$start_bp
    all_tsd_pos2 <- tsd_intact$end_bp
  }
  by_class[["target_site_duplication_all"]] <- list(
    densities   = as.numeric(count_per_window(all_tsd_pos, all_tsd_pos2,
                                              opts$step_bp, grid$n_windows)),
    max_density = NA, n_records = length(all_tsd_pos)
  )
  by_class[["target_site_duplication_all"]]$max_density <-
    if (grid$n_windows > 0) max(by_class[["target_site_duplication_all"]]$densities) else 0
  classes_emitted <- c(classes_emitted, "target_site_duplication_all")

  # (5) Intact-element count from intact.gff3 (repeat_region parents)
  intact_parents <- ig_chr[type == "repeat_region"]
  by_class[["intact_element_count"]] <- list(
    densities   = as.numeric(count_per_window(intact_parents$start_bp,
                                              intact_parents$end_bp,
                                              opts$step_bp, grid$n_windows)),
    max_density = NA, n_records = nrow(intact_parents)
  )
  by_class[["intact_element_count"]]$max_density <-
    if (grid$n_windows > 0) max(by_class[["intact_element_count"]]$densities) else 0
  classes_emitted <- c(classes_emitted, "intact_element_count")

  # -------------------------------------------------------------------------
  # Build the JSON object
  # -------------------------------------------------------------------------
  out <- list(
    version          = 2L,
    schema_version   = 2L,
    binning_source   = "scrubber_windows",
    species          = opts$species,
    precomp_chrom    = this_chrom,
    n_classes        = length(classes_emitted),
    classes          = as.list(classes_emitted),
    default_class    = "all_TE",
    loess_span       = 0.3,
    generated_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    layer_provenance = list(
      script         = "STEP_RD_TE_density_full.R",
      teanno_source  = basename(opts$teanno_gff3),
      intact_source  = basename(opts$intact_gff3),
      window_bp      = opts$window_bp,
      step_bp        = opts$step_bp,
      identity_young_threshold = YOUNG_THR,
      method = "TEanno + intact GFF3 binned to scrubber 5-kb grid; per-class densities are bp-coverage fractions; aggregates sum across classes; counts are raw event counts."
    ),
    chromosomes = list(
      list(
        chrom             = this_chrom,
        n_windows         = grid$n_windows,
        window_centers_mb = grid$window_centers_mb,
        window_start_bp   = grid$window_start_bp,
        window_end_bp     = grid$window_end_bp,
        by_class          = by_class
      )
    )
  )

  out_path <- file.path(opts$out_dir,
                        paste0(this_chrom, "_repeat_density_TEfull.json"))
  jsonlite::write_json(out, out_path, auto_unbox = TRUE,
                       null = "null", na = "null", digits = 6)

  size_kb <- round(file.info(out_path)$size / 1024, 0)
  n_classes_ok <- length(classes_emitted)
  cat(sprintf("  wrote %s · n_windows=%d · classes=%d · size=%dKB\n",
              basename(out_path), grid$n_windows, n_classes_ok, size_kb))
}

cat("[STEP_RD_TE_FULL] done.\n")
