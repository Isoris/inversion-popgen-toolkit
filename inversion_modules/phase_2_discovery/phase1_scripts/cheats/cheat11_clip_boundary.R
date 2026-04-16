#!/usr/bin/env Rscript
# =============================================================================
# cheat11_clip_boundary.R — Soft/hard clipped read pileup at boundaries
#
# At real inversion breakpoints, reads that span the junction get
# soft-clipped or hard-clipped by the aligner. This creates a pileup
# of clipped reads at the exact breakpoint position.
#
# Signal:
#   - High clip count at position → structural breakpoint (INV, DEL, etc.)
#   - Clip count bimodal across samples → breakpoint is polymorphic
#   - HET samples: clip count intermediate (one haplotype clips, one doesn't)
#   - HOM_INV samples: low clip count (both haplotypes have the rearrangement,
#     so the aligner maps to the wrong location instead of clipping)
#   - HOM_REF samples: low clip count (no rearrangement to clip against)
#
# Speed strategy:
#   OPTION A (fast): Pre-compute a genome-wide clip-depth track per sample
#     using samtools in one pass. Then querying boundaries = table lookup.
#     Cost: ~10 min/sample × 226 = ~37 hours (SLURM array, parallelizable)
#
#   OPTION B (targeted): For each boundary, query the BAM with samtools view
#     in a ±2 kb window. Parse CIGAR for S/H operations.
#     Cost: ~0.5 sec/sample/boundary × 226 × 200 boundaries = ~6 hours
#
#   We implement OPTION B (targeted) here. OPTION A is a preprocessing
#   step that can be added later for speed.
#
# Multi-scale: ±80/±160/±320 windows each side for CONTEXT,
#   but clip counting is always in a NARROW window (±2 kb) around boundary.
#   The multi-scale context provides baseline clip rate for normalization.
#
# REQUIRES: samtools in PATH, load_bridge.R sourced
# =============================================================================

suppressPackageStartupMessages(library(data.table))

CLIP_WINDOW    <- 2000L    # ±2 kb around boundary for clip counting
CONTEXT_WINDOW <- 50000L   # ±50 kb for baseline clip rate
MAX_SAMPLES    <- 40L      # subsample for speed (stratified by band)
MIN_CLIP_LEN   <- 5L       # minimum clip length to count (ignore tiny clips)

#' Count clipped reads in a region for one sample
#' @param bam_path Path to markdup BAM
#' @param chr Chromosome
#' @param start_bp, end_bp Region bounds
#' @return list(n_softclip, n_hardclip, n_total_reads, clip_fraction)
count_clips_region <- function(bam_path, chr, start_bp, end_bp) {
  result <- list(n_softclip = NA_integer_, n_hardclip = NA_integer_,
                 n_total_reads = NA_integer_, clip_fraction = NA_real_)

  if (!file.exists(bam_path)) return(result)

  region <- paste0(chr, ":", as.integer(start_bp), "-", as.integer(end_bp))

  # Use samtools view + awk to count clips from CIGAR
  # CIGAR field is $6. S = soft clip, H = hard clip.
  # Count reads with S≥5 or H≥5 in CIGAR.
  cmd <- paste0(
    "samtools view -F 0x904 ", bam_path, " ", region, " 2>/dev/null | ",
    "awk '{",
    "  total++; ",
    "  cigar=$6; ",
    "  s=0; h=0; ",
    "  while (match(cigar, /([0-9]+)([SH])/, a)) {",
    "    if (a[2]==\"S\" && a[1]+0>=", MIN_CLIP_LEN, ") s++; ",
    "    if (a[2]==\"H\" && a[1]+0>=", MIN_CLIP_LEN, ") h++; ",
    "    cigar=substr(cigar, RSTART+RLENGTH); ",
    "  }",
    "  if (s>0) sc++; if (h>0) hc++; ",
    "} END {",
    "  print total+0, sc+0, hc+0",
    "}'"
  )

  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) NULL)
  if (is.null(out) || length(out) == 0) return(result)

  parts <- as.integer(strsplit(trimws(out[1]), " ")[[1]])
  if (length(parts) >= 3) {
    result$n_total_reads <- parts[1]
    result$n_softclip    <- parts[2]
    result$n_hardclip    <- parts[3]
    if (parts[1] > 0) {
      result$clip_fraction <- (parts[2] + parts[3]) / parts[1]
    }
  }
  result
}

#' Find BAM file for a sample
find_bam <- function(sample_id, base_dir = NULL) {
  if (is.null(base_dir)) base_dir <- Sys.getenv("BASE",
    "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04")
  bam <- file.path(base_dir, "delly_sv/00_markdup",
                     paste0(sample_id, ".markdup.bam"))
  if (file.exists(bam)) bam else NULL
}

#' Compute clip boundary score at one position
#' @param chr Chromosome
#' @param boundary_bp Boundary position
#' @param sample_ids CGA names (already subsampled)
#' @param pc1_bands Named integer vector (optional, for per-band analysis)
#' @param boundary_type "hard" or "soft"
#' @return data.table row
compute_clip_boundary <- function(chr, boundary_bp, sample_ids,
                                   pc1_bands = NULL, boundary_type = "hard") {
  # Narrow window: clip counting
  narrow_start <- boundary_bp - CLIP_WINDOW
  narrow_end   <- boundary_bp + CLIP_WINDOW

  # Context window: baseline
  ctx_start <- boundary_bp - CONTEXT_WINDOW
  ctx_end   <- boundary_bp + CONTEXT_WINDOW

  per_sample <- list()
  for (sid in sample_ids) {
    bam <- find_bam(sid)
    if (is.null(bam)) next

    narrow <- count_clips_region(bam, chr, narrow_start, narrow_end)
    context <- count_clips_region(bam, chr, ctx_start, ctx_end)

    # Enrichment: is clip fraction at boundary HIGHER than context?
    enrichment <- NA_real_
    if (is.finite(narrow$clip_fraction) && is.finite(context$clip_fraction) &&
        context$clip_fraction > 0) {
      enrichment <- narrow$clip_fraction / context$clip_fraction
    }

    per_sample[[length(per_sample) + 1]] <- data.table(
      sample_id = sid,
      n_reads_narrow = narrow$n_total_reads,
      n_softclip = narrow$n_softclip,
      n_hardclip = narrow$n_hardclip,
      clip_frac_narrow = round(narrow$clip_fraction, 4),
      clip_frac_context = round(context$clip_fraction, 4),
      clip_enrichment = round(enrichment, 2)
    )
  }

  if (length(per_sample) == 0) {
    return(data.table(boundary_bp = boundary_bp, boundary_type = boundary_type,
                       mean_clip_frac = NA_real_, clip_enrichment = NA_real_,
                       clip_bimodal = NA, clip_score = NA_real_))
  }

  ps_dt <- rbindlist(per_sample)

  # Summary metrics
  mean_clip_frac <- mean(ps_dt$clip_frac_narrow, na.rm = TRUE)
  mean_enrichment <- mean(ps_dt$clip_enrichment, na.rm = TRUE)

  # Bimodality: is clip fraction bimodal across samples?
  # (some samples clip a lot = breakpoint carriers; others don't = non-carriers)
  clip_bimodal <- FALSE
  valid_cf <- ps_dt$clip_frac_narrow[is.finite(ps_dt$clip_frac_narrow)]
  if (length(valid_cf) >= 15) {
    cv <- sd(valid_cf) / max(mean(valid_cf), 0.001)
    clip_bimodal <- cv > 0.5  # high CV suggests bimodality
  }

  # Per-band clip analysis (if bands available)
  clip_het_excess <- NA_real_
  if (!is.null(pc1_bands)) {
    het_ids <- intersect(names(pc1_bands)[pc1_bands == 2], ps_dt$sample_id)
    hom_ids <- intersect(names(pc1_bands)[pc1_bands %in% c(1, 3)], ps_dt$sample_id)
    if (length(het_ids) >= 3 && length(hom_ids) >= 3) {
      het_clip <- mean(ps_dt[sample_id %in% het_ids, clip_frac_narrow], na.rm = TRUE)
      hom_clip <- mean(ps_dt[sample_id %in% hom_ids, clip_frac_narrow], na.rm = TRUE)
      if (is.finite(het_clip) && is.finite(hom_clip) && hom_clip > 0) {
        clip_het_excess <- het_clip / hom_clip
      }
    }
  }

  # Composite score
  score <- 0; n_c <- 0
  if (is.finite(mean_enrichment) && mean_enrichment > 1.5) {
    score <- score + pmin(1, (mean_enrichment - 1) / 3) * 0.4; n_c <- n_c + 1
  }
  if (clip_bimodal) { score <- score + 0.3; n_c <- n_c + 1 }
  if (is.finite(clip_het_excess) && clip_het_excess > 1.3) {
    score <- score + pmin(1, (clip_het_excess - 1) / 2) * 0.3; n_c <- n_c + 1
  }
  if (n_c > 0) score <- pmin(1, score)

  data.table(
    boundary_bp = boundary_bp,
    boundary_type = boundary_type,
    n_samples_queried = nrow(ps_dt),
    mean_clip_frac = round(mean_clip_frac, 4),
    mean_clip_enrichment = round(mean_enrichment, 2),
    clip_bimodal = clip_bimodal,
    clip_het_excess = round(clip_het_excess, 2),
    total_softclips = sum(ps_dt$n_softclip, na.rm = TRUE),
    total_hardclips = sum(ps_dt$n_hardclip, na.rm = TRUE),
    clip_score = round(score, 4)
  )
}

#' Scan all boundaries of a candidate for clip evidence
#' @param chr Chromosome
#' @param boundaries_dt data.table with boundary_bp, boundary_type
#' @param sample_ids CGA sample names
#' @param pc1_bands Named integer vector (optional)
#' @return data.table with per-boundary clip scores
scan_clip_boundaries <- function(chr, boundaries_dt, sample_ids,
                                  pc1_bands = NULL) {
  if (nrow(boundaries_dt) == 0) return(data.table())

  # Subsample for speed
  sub_ids <- sample_ids
  if (length(sample_ids) > MAX_SAMPLES) {
    if (!is.null(pc1_bands)) {
      sub_ids <- c(
        sample(names(pc1_bands)[pc1_bands == 1], min(13, sum(pc1_bands == 1))),
        sample(names(pc1_bands)[pc1_bands == 2], min(14, sum(pc1_bands == 2))),
        sample(names(pc1_bands)[pc1_bands == 3], min(13, sum(pc1_bands == 3)))
      )
    } else {
      sub_ids <- sample(sample_ids, MAX_SAMPLES)
    }
  }

  results <- lapply(seq_len(nrow(boundaries_dt)), function(bi) {
    message("[cheat11] Boundary ", bi, "/", nrow(boundaries_dt),
            " (", boundaries_dt$boundary_type[bi], ") @ ",
            boundaries_dt$boundary_bp[bi])
    compute_clip_boundary(
      chr, boundaries_dt$boundary_bp[bi], sub_ids,
      pc1_bands, boundaries_dt$boundary_type[bi]
    )
  })
  rbindlist(results, fill = TRUE)
}

#' Pre-generate genome-wide clip-depth track for one sample (OPTION A)
#' Writes a BED-like TSV: chrom, start, end, n_reads, n_softclip, n_hardclip
#' Run as SLURM array: one job per sample
#' @param bam_path BAM file
#' @param output_path Output TSV path
#' @param bin_size Bin size in bp (default 1000)
precompute_clip_track <- function(bam_path, output_path, bin_size = 1000L) {
  if (!file.exists(bam_path)) stop("BAM not found: ", bam_path)

  # Use samtools + awk to bin clip counts genome-wide
  cmd <- paste0(
    "samtools view -F 0x904 ", bam_path, " | ",
    "awk -v BS=", bin_size, " -v MINC=", MIN_CLIP_LEN, " '",
    "BEGIN{OFS=\"\\t\"} {",
    "  chr=$3; pos=int($4/BS)*BS; ",
    "  key=chr\"\\t\"pos; ",
    "  total[key]++; ",
    "  cigar=$6; ",
    "  while (match(cigar, /([0-9]+)([SH])/, a)) {",
    "    if (a[2]==\"S\" && a[1]+0>=MINC) sc[key]++; ",
    "    if (a[2]==\"H\" && a[1]+0>=MINC) hc[key]++; ",
    "    cigar=substr(cigar, RSTART+RLENGTH); ",
    "  }",
    "} END {",
    "  for (k in total) print k, k+BS, total[k], sc[k]+0, hc[k]+0",
    "}' | sort -k1,1 -k2,2n | gzip > ", output_path
  )
  message("[cheat11] Precomputing clip track: ", basename(bam_path))
  system(cmd)
  message("[cheat11] Done: ", output_path)
}
