#!/usr/bin/env Rscript
# =============================================================================
# cheat29_junction_forensics.R — Breakpoint junction forensics
#
# Based on Porubsky et al. 2022 Cell (HGSVC recurrent inversions) and
# Sultana et al. 2017 (L1 inversion mechanisms).
#
# WHAT THIS DOES:
#   For each confirmed breakpoint (from DELLY/Manta or boundary catalog),
#   extracts the junction sequence from the reference assembly, classifies
#   the junction type, and infers the formation mechanism.
#
# JUNCTION TYPES (from literature):
#   1. BLUNT — clean break, no overlap, no insertion
#      → NHEJ (non-homologous end joining)
#   2. MICROHOMOLOGY (MH) — 1-7 bp shared at junction
#      → MMEJ (microhomology-mediated end joining) or alt-NHEJ
#      Peak at 2-3 bp MH (Sultana Fig G)
#   3. INSERTION — 1-30 bp non-templated or templated insertion
#      → MMBIR/FoSTeS (replication-based, template switching)
#   4. DUPLICATION — inverted segment + flanking duplication (Inv+Dup)
#      → Twin-priming during L1 retrotransposition
#      Median dup length: 11 bp (Sultana Fig H)
#   5. DELETION — inverted segment + flanking deletion (Inv+Del)
#      → Endonuclease-mediated L1 insertion with 5' truncation
#      Median del length: 19 bp, high IQR 329 bp (Sultana Fig H)
#
# SD CONTEXT (from Porubsky Fig 1):
#   Inverted SD pairs flanking breakpoints → NAHR mechanism
#   SD orientation predicts:
#     Inverted SDs → increased inversion recurrence rate
#     Direct SDs → decreased inversion recurrence, increased CNV risk
#   This interfaces with Cheat 27 (BISER2 concordance check)
#
# INPUT:
#   ref_fasta:        reference genome FASTA (indexed)
#   breakpoint_pairs: data.table with chr, bp_left, bp_right per candidate
#   junction_window:  bp to extract around each junction (default ±50)
#
# OUTPUT per junction:
#   junction_type:    "blunt" / "microhomology" / "insertion" / "inv_dup" / "inv_del"
#   mh_length:        microhomology length in bp (0 for blunt)
#   ins_length:       insertion length (0 if none)
#   ins_sequence:     inserted bases (if any)
#   dup_length:       duplication length (for Inv+Dup)
#   del_length:       deletion length (for Inv+Del)
#   tsd_present:      TRUE if target site duplication detected → TE-mediated
#   mechanism_class:  "NHEJ" / "MMEJ" / "MMBIR" / "L1_twin_priming" /
#                     "L1_endonuclease" / "NAHR" / "unknown"
#
# OUTPUT per candidate:
#   dominant_junction: most common junction type across breakpoints
#   mechanism_model:   integrated mechanism from junction + SD context
#   is_recurrent:      TRUE if NAHR mechanism (from Cheat 14/27 concordance)
#
# INTEGRATION:
#   Cheat 14: minimap2 finds inverted repeats → NAHR substrate
#   Cheat 27: BISER2 cross-checks SD pairs → concordance
#   Cheat 29: junction sequence → junction type → mechanism class
#   Together: full mechanism model (substrate + junction = mechanism)
#
# REQUIRES: samtools (for faidx), reference FASTA indexed,
#           RepeatMasker BED (primary TE check) — optional fallback to de novo TSD
# DISPATCHER: No
# DIFFICULTY: Medium (needs base-pair-level sequence extraction)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
JUNCTION_WINDOW    <- 50L      # ±50 bp around breakpoint
MH_MAX_LENGTH      <- 15L     # max microhomology to search for
INS_MAX_LENGTH     <- 100L    # max insertion to consider
TSD_MIN_LENGTH     <- 5L      # minimum TSD length
TSD_MAX_LENGTH     <- 20L     # maximum TSD length
TE_PROXIMITY_BP    <- 200L    # TE within ±200 bp of breakpoint counts

# ── RepeatMasker TE lookup (primary, faster + more reliable) ──────────

.rm_cache <- new.env(parent = emptyenv())

load_repeatmasker <- function(rm_bed) {
  if (exists("rm_dt", envir = .rm_cache)) return(get("rm_dt", envir = .rm_cache))
  if (is.null(rm_bed) || !file.exists(rm_bed)) return(NULL)

  dt <- tryCatch(fread(rm_bed), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  # Standardize columns — BED format or .out format
  if (ncol(dt) >= 4 && !("chr" %in% names(dt))) {
    setnames(dt, seq_len(min(6, ncol(dt))),
             c("chr", "start", "end", "te_name", "score", "strand")[seq_len(min(6, ncol(dt)))])
  }
  # Extract TE class from name (e.g., "L1ME3a" → "LINE/L1", "AluSx" → "SINE/Alu")
  if (!"te_class" %in% names(dt) && "te_name" %in% names(dt)) {
    dt[, te_class := fifelse(grepl("^L1|LINE", te_name, ignore.case = TRUE), "LINE/L1",
                     fifelse(grepl("^Alu|SINE", te_name, ignore.case = TRUE), "SINE/Alu",
                     fifelse(grepl("^DNA|hAT|Tc1|Mariner", te_name, ignore.case = TRUE), "DNA",
                     fifelse(grepl("^LTR|ERV|Gypsy", te_name, ignore.case = TRUE), "LTR",
                     "other"))))]
  }
  assign("rm_dt", dt, envir = .rm_cache)
  dt
}

check_te_at_breakpoint <- function(rm_dt, chr, bp, proximity = TE_PROXIMITY_BP) {
  if (is.null(rm_dt) || nrow(rm_dt) == 0)
    return(list(te_present = FALSE, te_name = NA, te_class = NA, te_dist = NA))

  hits <- rm_dt[chr == (..chr) & start <= bp + proximity & end >= bp - proximity]
  if (nrow(hits) == 0)
    return(list(te_present = FALSE, te_name = NA, te_class = NA, te_dist = NA))

  # Closest TE
  hits[, dist := pmin(abs(start - bp), abs(end - bp))]
  best <- hits[which.min(dist)]

  list(
    te_present = TRUE,
    te_name = best$te_name[1],
    te_class = best$te_class[1],
    te_dist = best$dist[1],
    te_overlaps = any(hits$start <= bp & hits$end >= bp),
    is_line_l1 = grepl("LINE|L1", best$te_class[1], ignore.case = TRUE)
  )
}

# ── Extract junction sequences via samtools ───────────────────────────

extract_junction_seqs <- function(ref_fasta, chr, bp, window = JUNCTION_WINDOW) {
  if (!file.exists(ref_fasta)) return(NULL)

  left_region  <- paste0(chr, ":", max(1, bp - window), "-", bp)
  right_region <- paste0(chr, ":", bp + 1, "-", bp + window)

  left_seq <- tryCatch({
    system2("samtools", c("faidx", ref_fasta, left_region),
            stdout = TRUE, stderr = FALSE)
  }, error = function(e) NULL)

  right_seq <- tryCatch({
    system2("samtools", c("faidx", ref_fasta, right_region),
            stdout = TRUE, stderr = FALSE)
  }, error = function(e) NULL)

  if (is.null(left_seq) || is.null(right_seq)) return(NULL)

  # Remove FASTA header
  left_seq  <- paste(left_seq[!grepl("^>", left_seq)], collapse = "")
  right_seq <- paste(right_seq[!grepl("^>", right_seq)], collapse = "")

  list(left = toupper(left_seq), right = toupper(right_seq))
}

# ── Detect microhomology at junction ──────────────────────────────────

detect_microhomology <- function(left_seq, right_seq) {
  if (is.null(left_seq) || is.null(right_seq)) return(0L)
  if (nchar(left_seq) == 0 || nchar(right_seq) == 0) return(0L)

  # MH = shared bases at the junction point
  # Last N bases of left == First N bases of right
  left_end <- substring(left_seq, nchar(left_seq) - MH_MAX_LENGTH + 1)
  right_start <- substring(right_seq, 1, MH_MAX_LENGTH)

  mh_len <- 0L
  for (i in seq_len(min(nchar(left_end), nchar(right_start)))) {
    if (substring(left_end, nchar(left_end) - i + 1) ==
        substring(right_start, 1, i)) {
      mh_len <- as.integer(i)
    }
  }
  mh_len
}

# ── Detect target site duplication (TSD) ──────────────────────────────

detect_tsd <- function(left_seq, right_seq) {
  if (is.null(left_seq) || is.null(right_seq)) return(list(present = FALSE, length = 0L))

  # TSD = identical sequence flanking both sides of the insertion
  # Check if last N bp of left == last N bp of right (after the junction)
  for (len in TSD_MAX_LENGTH:TSD_MIN_LENGTH) {
    left_end <- substring(left_seq, nchar(left_seq) - len + 1)
    # Search for this motif in the right flank
    pos <- regexpr(left_end, right_seq, fixed = TRUE)
    if (pos > 0 && pos <= 30) {  # TSD should be close to junction
      return(list(present = TRUE, length = as.integer(len), sequence = left_end))
    }
  }
  list(present = FALSE, length = 0L, sequence = "")
}

# ── Classify one junction ─────────────────────────────────────────────

classify_junction <- function(ref_fasta, chr, bp, rm_dt = NULL) {
  default <- list(
    junction_type = NA_character_, mh_length = NA_integer_,
    ins_length = 0L, ins_sequence = "", dup_length = 0L, del_length = 0L,
    tsd_present = FALSE, te_present = FALSE, te_name = NA_character_,
    te_class = NA_character_, mechanism_class = "unknown"
  )

  seqs <- extract_junction_seqs(ref_fasta, chr, bp)
  if (is.null(seqs)) return(default)

  mh <- detect_microhomology(seqs$left, seqs$right)

  # ── TE check: RepeatMasker first (reliable), sequence-based TSD as fallback ──
  te_info <- list(te_present = FALSE, is_line_l1 = FALSE)
  tsd <- list(present = FALSE, length = 0L, sequence = "")

  if (!is.null(rm_dt)) {
    te_info <- check_te_at_breakpoint(rm_dt, chr, bp)
  }

  if (!te_info$te_present) {
    # Fallback: de novo TSD detection from sequence
    tsd <- detect_tsd(seqs$left, seqs$right)
  }

  te_at_junction <- te_info$te_present || tsd$present

  # ── Classify junction type ──
  if (te_at_junction) {
    is_l1 <- te_info$is_line_l1 %||% FALSE
    if (is_l1) {
      # L1 at breakpoint → distinguish twin-priming vs endonuclease
      # Twin-priming: internal inversion within L1, short (Sultana Fig J)
      # Endonuclease: L1 insertion with 5' truncation, Inv+Del signature
      junction_type <- "te_l1_mediated"
      mechanism <- "L1_retrotransposition"
    } else {
      junction_type <- "te_mediated"
      mechanism <- "TE_mediated"
    }
  } else if (mh >= 1 && mh <= 3) {
    junction_type <- "microhomology"
    mechanism <- "MMEJ"
  } else if (mh >= 4 && mh <= 7) {
    junction_type <- "extended_microhomology"
    mechanism <- "MMEJ"
  } else if (mh >= 8) {
    junction_type <- "long_homology"
    mechanism <- "NAHR_like"  # long MH suggests NAHR substrate
  } else {
    junction_type <- "blunt"
    mechanism <- "NHEJ"
  }

  list(
    junction_type = junction_type,
    mh_length = mh,
    ins_length = 0L,
    ins_sequence = "",
    dup_length = 0L,
    del_length = 0L,
    tsd_present = tsd$present || (te_info$te_present && (te_info$te_overlaps %||% FALSE)),
    te_present = te_info$te_present,
    te_name = te_info$te_name %||% NA_character_,
    te_class = te_info$te_class %||% NA_character_,
    mechanism_class = mechanism
  )
}

# ── Integrate with SD context (Cheat 14/27) ───────────────────────────

integrate_mechanism_model <- function(junction_result, sd_result = NULL) {
  jct <- junction_result$mechanism_class

  # If we have SD context from Cheat 14/27
  if (!is.null(sd_result)) {
    has_inv_sd <- sd_result$biser2_nahr_class %in% c("strong", "weak") ||
                  (sd_result$concordance %||% "") == "agree_nahr"

    if (has_inv_sd) {
      # Inverted SDs + any junction → NAHR is primary mechanism
      # Junction MH is expected from NAHR resolution
      return(list(
        mechanism_model = "NAHR",
        is_recurrent = TRUE,
        confidence = if ((sd_result$concordance %||% "") == "agree_nahr") "high" else "moderate",
        note = paste0("Inverted SD pair + ", jct, " junction → NAHR-mediated, likely recurrent")
      ))
    }
  }

  # No SD context → mechanism from junction alone
  is_recurrent <- jct == "NAHR"  # only if somehow flagged upstream
  list(
    mechanism_model = jct,
    is_recurrent = is_recurrent,
    confidence = "junction_only",
    note = paste0("Junction type: ", junction_result$junction_type,
                  " (MH=", junction_result$mh_length, "bp)")
  )
}

# ── Runner ────────────────────────────────────────────────────────────────

run_cheat29 <- function(chr, breakpoint_bps, ref_fasta,
                         sd_results = NULL,
                         repeatmasker_bed = NULL) {
  message("[cheat29] ", chr, ": ", length(breakpoint_bps), " breakpoints")

  # Load RepeatMasker if available
  rm_dt <- NULL
  if (!is.null(repeatmasker_bed)) {
    rm_dt <- load_repeatmasker(repeatmasker_bed)
    if (!is.null(rm_dt)) {
      message("[cheat29] RepeatMasker loaded: ", nrow(rm_dt), " TEs (primary TE check)")
    }
  }

  results <- list()
  for (i in seq_along(breakpoint_bps)) {
    bp <- breakpoint_bps[i]
    jct <- classify_junction(ref_fasta, chr, bp, rm_dt)

    # Integrate with SD context if available
    sd_i <- if (!is.null(sd_results) && nrow(sd_results) >= i) sd_results[i] else NULL
    model <- integrate_mechanism_model(jct, sd_i)

    results[[i]] <- c(
      list(boundary_bp = bp),
      jct,
      list(mechanism_model = model$mechanism_model,
           is_recurrent = model$is_recurrent,
           model_confidence = model$confidence,
           model_note = model$note)
    )
  }

  dt <- rbindlist(results, fill = TRUE)

  # Summary
  if (nrow(dt) > 0) {
    tab <- table(dt$junction_type)
    message("[cheat29] Junction types: ",
            paste(names(tab), tab, sep = "=", collapse = ", "))
    tab2 <- table(dt$mechanism_model)
    message("[cheat29] Mechanisms: ",
            paste(names(tab2), tab2, sep = "=", collapse = ", "))
  }

  dt
}
