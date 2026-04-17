#!/usr/bin/env Rscript
# =============================================================================
# cheat28_tandem_repeat_context.R — STR/VNTR breakpoint context + theta filter
#
# TWO USES:
#
# USE 1 — MECHANISM CLASSIFICATION (feeds Cheat 14/27):
#   STR/VNTR expansions at breakpoints indicate replication-based mechanisms
#   (FoSTeS/MMBIR) rather than NAHR or NHEJ. This is a THIRD mechanism class:
#     NAHR: inverted SD pair flanking breakpoint (Cheat 27)
#     NHEJ: blunt/microhomology junction, no flanking repeats (Cheat 14/20)
#     MMBIR/FoSTeS: tandem repeat expansion/insertion at junction
#
#   Mac haplotype: 570,453 TRF loci (8.67% genome), mean tract 136 bp
#   Gar haplotype: 559,925 TRF loci (8.34% genome), mean tract 144 bp
#   → ~1.07 SSR/kb (Mac), dense enough for breakpoint-level annotation
#
# USE 2 — THETA NOISE FILTER (feeds Cheat 12):
#   STR/VNTR length polymorphism inside inversions inflates per-sample θ_P.
#   Windows with high microsatellite density may show elevated theta NOT
#   because of heterokaryotype divergence but because of STR mutation rate.
#   Flagging these windows lets Cheat 12 downweight STR-driven noise.
#
#   Dinucleotide CA/TG + GT/AC = 53% of all SSRs → main noise source.
#   Tri/tetranucleotides are less polymorphic, less problematic.
#
# INPUT:
#   trf_bed:      TRF output in BED format (chr, start, end, period, copies, ...)
#   gmata_bed:    GMATA SSR output (chr, start, end, motif, motif_len, copies)
#   boundary_catalog: from C01g
#   window_dt:  from C01a (for theta filter annotation)
#
# OUTPUT per breakpoint (USE 1):
#   n_tr_at_bp:       count of tandem repeats within ±5 kb
#   tr_density_ratio: density at bp / genome-wide density
#   dominant_tr_class: "microsatellite" / "minisatellite" / "none"
#   dominant_motif:    e.g. "CA" or "AAT"
#   mmbir_signal:     TRUE if TR enrichment + insertion at junction
#
# OUTPUT per window (USE 2):
#   ssr_density:      SSRs per kb in this window
#   ssr_bp_fraction:  fraction of window covered by SSRs
#   dinuc_fraction:   fraction of SSRs that are dinucleotide (noisiest)
#   theta_noise_flag: TRUE if ssr_density > 2× genome-wide → downweight
#
# INTEGRATION:
#   C01g: annotate boundaries with TR context → cheat28_* columns
#   C01a: annotate window_dt with ssr_density → theta_noise_flag
#   Cheat 12: downweight windows where theta_noise_flag == TRUE
#   Cheat 14: add MMBIR mechanism class from TR context
#
# REQUIRES: TRF BED and/or GMATA BED
# DISPATCHER: No
# DIFFICULTY: Easy (BED interval overlap)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
BP_FLANK_TR       <- 5000L    # ±5 kb for breakpoint TR lookup
WINDOW_SIZE_TR    <- 50000L   # 50 kb window for density
MICROSATELLITE_MAX <- 6L      # period ≤ 6 bp = microsatellite
MINISATELLITE_MAX  <- 60L     # period 7-60 bp = minisatellite
NOISE_MULTIPLIER   <- 2.0     # flag if density > 2× genome-wide

# ── Load TRF output ──────────────────────────────────────────────────────

load_trf <- function(trf_file) {
  if (!file.exists(trf_file)) {
    message("[cheat28] TRF file not found: ", trf_file)
    return(data.table())
  }
  dt <- fread(trf_file)
  # TRF BED typically: chr, start, end, period, copies, consensus_size,
  #                     pct_matches, pct_indels, score, A, C, G, T, entropy, motif
  if (ncol(dt) >= 5) {
    if (!"chr" %in% names(dt)) {
      setnames(dt, seq_len(min(5, ncol(dt))),
               c("chr", "start", "end", "period", "copies")[seq_len(min(5, ncol(dt)))])
    }
  }
  dt[, tr_class := fifelse(period <= MICROSATELLITE_MAX, "microsatellite",
                    fifelse(period <= MINISATELLITE_MAX, "minisatellite", "macrosatellite"))]
  dt[, tr_length := end - start]
  dt
}

# ── Load GMATA SSR output ─────────────────────────────────────────────────

load_gmata <- function(gmata_file) {
  if (!file.exists(gmata_file)) {
    message("[cheat28] GMATA file not found: ", gmata_file)
    return(data.table())
  }
  dt <- fread(gmata_file)
  if (ncol(dt) >= 4 && !"chr" %in% names(dt)) {
    setnames(dt, seq_len(min(6, ncol(dt))),
             c("chr", "start", "end", "motif", "motif_len", "copies")[seq_len(min(6, ncol(dt)))])
  }
  if ("motif_len" %in% names(dt)) {
    dt[, is_dinucleotide := motif_len == 2]
  } else if ("motif" %in% names(dt)) {
    dt[, motif_len := nchar(motif)]
    dt[, is_dinucleotide := motif_len == 2]
  }
  dt
}

# ── USE 1: Breakpoint TR context ──────────────────────────────────────────

annotate_breakpoint_tr <- function(tr_dt, chr, boundary_bp,
                                    flank = BP_FLANK_TR) {
  if (nrow(tr_dt) == 0) {
    return(data.table(boundary_bp = boundary_bp,
                       n_tr_at_bp = 0L, tr_density_ratio = NA_real_,
                       dominant_tr_class = "none", mmbir_signal = FALSE))
  }

  # Genome-wide density for this chromosome
  chr_tr <- tr_dt[chr == (..chr)]
  if (nrow(chr_tr) == 0) chr_tr <- tr_dt[grepl(chr, chr)]  # fuzzy match
  genome_density <- nrow(chr_tr) / max(chr_tr$end, na.rm = TRUE) * 1000  # per kb

  results <- list()
  for (bp in boundary_bp) {
    local <- chr_tr[start <= bp + flank & end >= bp - flank]
    n_local <- nrow(local)
    local_density <- n_local / (2 * flank / 1000)  # per kb

    ratio <- if (genome_density > 0) local_density / genome_density else NA_real_

    dom_class <- "none"
    if (n_local > 0) {
      tab <- table(local$tr_class)
      dom_class <- names(which.max(tab))
    }

    # MMBIR signal: enriched TR + presence of larger tandem repeats
    mmbir <- n_local > 5 && ratio > 2.0 &&
      any(local$tr_class == "minisatellite" & local$tr_length > 100)

    results[[length(results) + 1]] <- data.table(
      boundary_bp = bp,
      n_tr_at_bp = n_local,
      tr_density_ratio = round(ratio, 2),
      dominant_tr_class = dom_class,
      mmbir_signal = mmbir
    )
  }
  rbindlist(results)
}

# ── USE 2: Per-window theta noise filter ──────────────────────────────────

compute_ssr_density_track <- function(ssr_dt, chr, window_starts, window_ends) {
  if (nrow(ssr_dt) == 0 || length(window_starts) == 0) {
    return(data.table(
      start_bp = window_starts, end_bp = window_ends,
      ssr_density = NA_real_, ssr_bp_fraction = NA_real_,
      dinuc_fraction = NA_real_, theta_noise_flag = NA
    ))
  }

  chr_ssr <- ssr_dt[chr == (..chr)]
  if (nrow(chr_ssr) == 0) chr_ssr <- ssr_dt[grepl(chr, chr)]

  # Genome-wide SSR density for this chr
  genome_ssr_density <- nrow(chr_ssr) / max(chr_ssr$end, na.rm = TRUE) * 1000

  results <- list()
  for (i in seq_along(window_starts)) {
    ws <- window_starts[i]; we <- window_ends[i]
    win_len <- (we - ws) / 1000  # kb

    local <- chr_ssr[start <= we & end >= ws]
    n_ssr <- nrow(local)
    density <- n_ssr / win_len
    bp_covered <- sum(pmin(local$end, we) - pmax(local$start, ws))
    bp_frac <- bp_covered / (we - ws)

    dinuc_frac <- NA_real_
    if (n_ssr > 0 && "is_dinucleotide" %in% names(local)) {
      dinuc_frac <- mean(local$is_dinucleotide, na.rm = TRUE)
    }

    noise_flag <- !is.na(density) && genome_ssr_density > 0 &&
      density > NOISE_MULTIPLIER * genome_ssr_density

    results[[i]] <- list(
      start_bp = ws, end_bp = we,
      ssr_density = round(density, 3),
      ssr_bp_fraction = round(bp_frac, 4),
      dinuc_fraction = round(dinuc_frac, 3),
      theta_noise_flag = noise_flag
    )
  }
  rbindlist(results)
}

# ── Runner ────────────────────────────────────────────────────────────────

run_cheat28 <- function(chr, boundary_bps = NULL,
                         window_starts = NULL, window_ends = NULL,
                         trf_file = NULL, gmata_file = NULL) {
  tr_dt <- data.table()
  ssr_dt <- data.table()

  if (!is.null(trf_file)) tr_dt <- load_trf(trf_file)
  if (!is.null(gmata_file)) ssr_dt <- load_gmata(gmata_file)

  result <- list()

  # USE 1: breakpoint annotation
  if (!is.null(boundary_bps) && length(boundary_bps) > 0 && nrow(tr_dt) > 0) {
    result$breakpoint_tr <- annotate_breakpoint_tr(tr_dt, chr, boundary_bps)
    n_mmbir <- sum(result$breakpoint_tr$mmbir_signal)
    message("[cheat28] ", chr, ": ", length(boundary_bps), " breakpoints, ",
            n_mmbir, " with MMBIR signal")
  }

  # USE 2: theta noise track
  if (!is.null(window_starts) && length(window_starts) > 0 && nrow(ssr_dt) > 0) {
    result$ssr_track <- compute_ssr_density_track(ssr_dt, chr,
                                                   window_starts, window_ends)
    n_noisy <- sum(result$ssr_track$theta_noise_flag, na.rm = TRUE)
    message("[cheat28] ", chr, ": ", length(window_starts), " windows, ",
            n_noisy, " flagged as theta-noisy (SSR-dense)")
  }

  result
}
