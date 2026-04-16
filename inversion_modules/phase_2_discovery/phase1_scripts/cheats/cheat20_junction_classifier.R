#!/usr/bin/env Rscript
# =============================================================================
# cheat20_junction_classifier.R — Breakpoint junction microhomology classifier
#
# BIOLOGY:
#   Junction sequences reveal formation mechanism: blunt → NHEJ,
#   microhomology 1-9 bp → MMEJ/alt-EJ, insertion → template switching.
#
# INPUT:  junction TSVs from cheat20_extract_junctions.py, boundary types
# OUTPUT: per-breakpoint junction class, null comparison
# REQUIRES: cheat20_extract_junctions.py (pysam), markdup BAMs
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
PYTHON_SCRIPT <- file.path(dirname(sys.frame(1)$ofile %||% "."),
                            "cheat20_extract_junctions.py")
N_RANDOM_POS  <- 100L    # random positions for null comparison

# ── Run Python extractor on one BAM at one breakpoint ─────────────────

extract_junction_one <- function(bam_path, chr, bp, ref_fasta,
                                  out_tsv = tempfile(fileext = ".tsv")) {
  cmd <- paste("python3", shQuote(PYTHON_SCRIPT),
               "--bam", shQuote(bam_path),
               "--chr", chr,
               "--bp", bp,
               "--ref", shQuote(ref_fasta),
               "--out", shQuote(out_tsv),
               "2>/dev/null")
  ret <- system(cmd, intern = FALSE)
  if (ret != 0 || !file.exists(out_tsv) || file.size(out_tsv) == 0)
    return(data.table())
  tryCatch(fread(out_tsv), error = function(e) data.table())
}

# ── Aggregate across samples for one breakpoint ──────────────────────

classify_breakpoint_junction <- function(bam_paths, chr, bp, ref_fasta) {
  all_results <- list()
  for (bam in bam_paths) {
    if (!file.exists(bam)) next
    res <- extract_junction_one(bam, chr, bp, ref_fasta)
    if (nrow(res) > 0) all_results[[length(all_results)+1]] <- res
  }
  if (length(all_results) == 0)
    return(data.table(breakpoint_bp = bp, junction_class = "NO_EVIDENCE",
                       microhomology_length = 0L, n_reads_total = 0L,
                       confidence = 0))

  combined <- rbindlist(all_results, fill = TRUE)
  # Majority vote on junction class
  class_tab <- table(combined$junction_class)
  majority  <- names(which.max(class_tab))
  total_reads <- sum(combined$n_supporting_reads, na.rm = TRUE)
  mean_microhom <- round(mean(combined$microhomology_length, na.rm = TRUE), 1)
  conf <- round(mean(combined$confidence, na.rm = TRUE), 3)

  data.table(breakpoint_bp = bp, junction_class = majority,
             microhomology_length = mean_microhom,
             n_reads_total = total_reads, confidence = conf)
}

# ── Compare junction classes across boundary types ────────────────────

compare_junction_classes <- function(junction_results, boundary_types) {
  if (nrow(junction_results) == 0 || length(boundary_types) == 0)
    return(data.table())

  junction_results[, boundary_type := boundary_types[as.character(breakpoint_bp)]]
  junction_results[!is.na(boundary_type), .N, by = .(junction_class, boundary_type)]
}

# ── Null comparison: random positions ─────────────────────────────────

junction_vs_null <- function(junction_results, random_results) {
  if (nrow(junction_results) == 0 || nrow(random_results) == 0)
    return(list(enrichment = NA_real_, ks_p = NA_real_))

  bp_microhom  <- junction_results$microhomology_length
  rnd_microhom <- random_results$microhomology_length

  # KS test: are distributions different?
  ks <- tryCatch(ks.test(bp_microhom, rnd_microhom), error = function(e) NULL)
  ks_p <- if (!is.null(ks)) ks$p.value else NA_real_

  # Enrichment: mean breakpoint / mean random
  enrich <- if (mean(rnd_microhom, na.rm = TRUE) > 0)
    round(mean(bp_microhom, na.rm = TRUE) / mean(rnd_microhom, na.rm = TRUE), 3)
  else NA_real_

  list(enrichment = enrich, ks_p = ks_p,
       bp_mean_microhom = round(mean(bp_microhom, na.rm = TRUE), 2),
       random_mean_microhom = round(mean(rnd_microhom, na.rm = TRUE), 2))
}

# ── Search mode ────────────────────────────────────────────────────────

search_junction_class <- function(chr, zone_start, zone_end, ...) {
  data.table(method = "junction_classifier", best_bp = NA_integer_,
             score = 0, is_precise = FALSE,
             detail = "requires_bam_not_search_compatible")
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat20 <- function(chr, boundary_bp, bam_paths = NULL,
                         ref_fasta = NULL, boundary_types = NULL) {
  message("[cheat20] ", chr, ": ", length(boundary_bp), " breakpoints, ",
          length(bam_paths), " BAMs")

  if (is.null(bam_paths) || length(bam_paths) == 0 || is.null(ref_fasta)) {
    message("[cheat20] BAMs or reference not provided")
    return(list(junction_results = data.table(),
                comparison = data.table(),
                null_test = list(),
                search_result = search_junction_class(chr, 0, 0)))
  }

  # Classify each breakpoint
  results <- list()
  for (bp in boundary_bp) {
    message("[cheat20]   Processing bp=", bp)
    res <- classify_breakpoint_junction(bam_paths, chr, bp, ref_fasta)
    results[[length(results)+1]] <- res
  }
  junction_dt <- rbindlist(results)

  message("[cheat20] Junction classes: ",
          paste(table(junction_dt$junction_class), collapse = ", "))

  # Compare across boundary types
  comparison <- if (!is.null(boundary_types))
    compare_junction_classes(junction_dt, boundary_types) else data.table()

  list(junction_results = junction_dt, comparison = comparison,
       null_test = list(),  # null comparison requires random positions
       search_result = search_junction_class(chr, 0, 0))
}
