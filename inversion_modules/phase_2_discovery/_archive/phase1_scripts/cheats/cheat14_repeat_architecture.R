#!/usr/bin/env Rscript
# =============================================================================
# cheat14_repeat_architecture.R — Repeat/SD architecture at breakpoints
#
# BIOLOGY:
#   NAHR (non-allelic homologous recombination) between inverted repeats
#   → recurrent inversion, can arise on multiple haplotype backgrounds.
#   NHEJ (non-homologous end joining) from two DSBs → unique origin.
#   Flanking repeat architecture reveals which mechanism created each INV.
#
# INPUT:  ref FASTA, boundary catalog, optional RepeatMasker .out
# OUTPUT: per-breakpoint mechanism class, TE context, dotplots
# REQUIRES: samtools, minimap2 on PATH; RepeatMasker optional
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
FLANK_BP         <- 50000L
NAHR_MIN_LEN     <- 500L
NAHR_MIN_IDENT   <- 85
NAHR_POSS_MIN    <- 100L
TE_WINDOW_BP     <- 10000L

# ── Flank extraction via samtools ──────────────────────────────────────

extract_breakpoint_flanks <- function(ref_fasta, boundary_bp, chr,
                                       flank_size = FLANK_BP,
                                       out_fasta = tempfile(fileext = ".fa")) {
  if (!file.exists(ref_fasta)) {
    message("[cheat14] Reference FASTA not found: ", ref_fasta); return(NULL)
  }
  regions <- character()
  for (bp in boundary_bp) {
    ls <- max(1, bp - flank_size); le <- bp
    rs <- bp; re <- bp + flank_size
    regions <- c(regions,
      paste0(chr, ":", ls, "-", le),
      paste0(chr, ":", rs, "-", re))
  }
  cmd <- paste0("samtools faidx ", shQuote(ref_fasta), " ",
                 paste(regions, collapse = " "),
                 " > ", shQuote(out_fasta), " 2>/dev/null")
  if (system(cmd) != 0 || file.size(out_fasta) == 0) {
    message("[cheat14] samtools faidx failed"); return(NULL)
  }
  out_fasta
}

# ── Self-alignment with minimap2 ──────────────────────────────────────

self_align_flanks <- function(flank_fasta, method = "minimap2") {
  if (is.null(flank_fasta) || !file.exists(flank_fasta)) return(data.table())
  paf <- tempfile(fileext = ".paf")
  cmd <- paste0("minimap2 -X -c --eqx ", shQuote(flank_fasta), " ",
                 shQuote(flank_fasta), " > ", shQuote(paf), " 2>/dev/null")
  if (system(cmd) != 0 || !file.exists(paf) || file.size(paf) == 0) {
    message("[cheat14] minimap2 self-alignment failed"); return(data.table())
  }
  dt <- tryCatch(
    fread(paf, header = FALSE, fill = TRUE, select = 1:12,
          col.names = c("qname","qlen","qstart","qend","strand",
                         "tname","tlen","tstart","tend","nmatch","alen","mapq")),
    error = function(e) data.table())
  unlink(paf)
  if (nrow(dt) == 0) return(dt)
  # remove self-hits and short
  dt <- dt[qname != tname | abs(qstart - tstart) > 100]
  dt <- dt[alen >= 50]
  dt[, identity := round(nmatch / alen * 100, 1)]
  dt[, orientation := fifelse(strand == "+", "direct", "inverted")]
  dt
}

# ── Mechanism classification ──────────────────────────────────────────

classify_mechanism <- function(paf_results, boundary_bp) {
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
  # Summarise per breakpoint pair
  inv_hits <- paf_results[orientation == "inverted"]
  results <- list()
  for (i in seq_along(boundary_bp)) {
    label <- paste0("bp_", i)
    # Look for inverted hits involving this breakpoint's flanks
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

# ── TE context annotation ────────────────────────────────────────────

annotate_te_context <- function(boundary_bp, chr, repeatmasker_bed,
                                 window = TE_WINDOW_BP) {
  if (is.null(repeatmasker_bed) || !file.exists(repeatmasker_bed))
    return(data.table(bp = boundary_bp, te_density = NA_real_,
                       te_at_bp = FALSE, te_class_dominant = NA_character_))
  rm_dt <- tryCatch(fread(repeatmasker_bed), error = function(e) data.table())
  if (nrow(rm_dt) == 0)
    return(data.table(bp = boundary_bp, te_density = NA_real_,
                       te_at_bp = FALSE, te_class_dominant = NA_character_))
  # Standardize columns
  expected <- c("chr","start","end","te_class")
  if (ncol(rm_dt) >= 4) setnames(rm_dt, seq_len(min(4, ncol(rm_dt))),
                                   expected[seq_len(min(4, ncol(rm_dt)))])
  rm_dt <- rm_dt[chr == (..chr)]
  results <- list()
  for (i in seq_along(boundary_bp)) {
    bp <- boundary_bp[i]
    local <- rm_dt[start <= bp + window & end >= bp - window]
    te_bp_cov <- sum(pmin(local$end, bp + window) -
                      pmax(local$start, bp - window))
    density <- te_bp_cov / (2 * window)
    at_bp <- any(local$start <= bp & local$end >= bp)
    dom_class <- if (nrow(local) > 0 && "te_class" %in% names(local))
      local[, .N, by = te_class][which.max(N)]$te_class else NA_character_
    results[[i]] <- data.table(bp = bp, te_density = round(density, 4),
                                te_at_bp = at_bp, te_class_dominant = dom_class)
  }
  rbindlist(results)
}

# ── Search mode ────────────────────────────────────────────────────────

search_repeat_architecture <- function(chr, zone_start, zone_end,
                                        ref_fasta = NULL, ...) {
  empty <- data.table(method = "repeat_architecture", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_ref")
  if (is.null(ref_fasta) || !file.exists(ref_fasta)) return(empty)
  mid <- as.integer((zone_start + zone_end) / 2)
  fa <- extract_breakpoint_flanks(ref_fasta, mid, chr, flank_size = FLANK_BP)
  if (is.null(fa)) return(empty)
  paf <- self_align_flanks(fa)
  if (nrow(paf) == 0)
    return(data.table(method = "repeat_architecture", best_bp = mid,
                       score = 0.2, is_precise = FALSE, detail = "no_repeats"))
  inv_hits <- paf[orientation == "inverted"]
  if (nrow(inv_hits) == 0)
    return(data.table(method = "repeat_architecture", best_bp = mid,
                       score = 0.2, is_precise = FALSE, detail = "direct_only"))
  best <- inv_hits[which.max(alen)]
  sc <- pmin(1, best$alen / 1000) * pmin(1, best$identity / 90)
  data.table(method = "repeat_architecture", best_bp = mid,
             score = round(sc, 3), is_precise = FALSE,
             detail = paste0("inv_repeat_", best$alen, "bp_", best$identity, "%"))
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat14 <- function(chr, boundary_bp, ref_fasta,
                         repeatmasker_bed = NULL) {
  message("[cheat14] ", chr, ": ", length(boundary_bp), " breakpoints")

  fa <- extract_breakpoint_flanks(ref_fasta, boundary_bp, chr)
  if (is.null(fa)) {
    message("[cheat14] Flank extraction failed")
    return(list(mechanism = data.table(), te_context = data.table(),
                paf = data.table()))
  }
  paf <- self_align_flanks(fa)
  message("[cheat14] Self-alignment: ", nrow(paf), " hits (",
          sum(paf$orientation == "inverted"), " inverted)")

  mech <- classify_mechanism(paf, boundary_bp)
  message("[cheat14] Mechanism: ",
          paste(mech$mechanism, collapse = ", "))

  te <- annotate_te_context(boundary_bp, chr, repeatmasker_bed)
  list(mechanism = mech, te_context = te, paf = paf)
}
