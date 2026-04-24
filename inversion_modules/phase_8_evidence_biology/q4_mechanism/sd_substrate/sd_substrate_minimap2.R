#!/usr/bin/env Rscript
# =============================================================================
# sd_substrate_minimap2.R  —  Angle A of SD substrate detection
# =============================================================================
#
# ROLE:
#   De novo minimap2 self-alignment on ±50 kb flanks around the candidate's
#   two breakpoints. Finds inverted repeats directly from the reference
#   FASTA. This is the TRUSTED angle for NAHR substrate calls (per user
#   preference) — angle B (BISER2 lookup) is a catalog cross-check.
#
#   Lives at:
#     q4_mechanism/sd_substrate/sd_substrate_minimap2.R
#
#   Companion scripts (same directory):
#     sd_substrate_biser2.R      — angle B (catalog lookup)
#     sd_substrate_concordance.R — reads A + B, writes joint verdict
#
# =============================================================================
# minimap2 INVOCATION (audit trail saved to self_align.cmd.txt)
# =============================================================================
#   minimap2 -ax asm10 -X --eqx --secondary=yes -N 50 -p 0.1 \
#       flanks.fa flanks.fa | samtools sort -O bam -o self_align.bam
#   samtools index self_align.bam
#
# Why these flags:
#   -ax asm10    assembly-to-assembly preset, ≤10% divergence. Catches SDs
#                at ≥90% identity (covers the biologically meaningful NAHR
#                substrate range 85-99%). Can be overridden with --mm2_preset.
#   -X           skip self/dual mappings (each query vs itself trivially).
#                Essential for self-alignment, otherwise every sequence
#                produces a trivial 100% hit to itself.
#   --eqx        use = / X CIGAR ops (exact identity computable from CIGAR).
#   --secondary=yes  keep secondary alignments. NECESSARY for SD detection
#                because an SD pair produces one primary + one secondary
#                alignment (A→B primary means B→A comes back as secondary).
#                Without this we lose half of every SD pair.
#   -N 50        report up to 50 secondary alignments per query. Default 5
#                is too few for repeat-rich regions (some SD-dense windows
#                have >20 partial-identity hits).
#   -p 0.1       secondary-to-primary score ratio threshold. Default 0.8
#                is too strict — drops partial-identity SD copies that
#                have score < 80% of the primary. 0.1 retains any SD with
#                ≥10% of primary score, letting the post-filter handle
#                quality.
#
# Post-filters applied AFTER minimap2:
#   - alignment length ≥ 50 bp (remove noise)
#   - remove hits where query and target regions overlap (true self-hits
#     that survived -X)
#   - classify orientation: + = direct, - = inverted
#   - inverted hits are the ones that matter for NAHR
#
# =============================================================================
# OUTPUTS (all in <outdir>/<candidate_id>/sd_substrate/minimap2/)
# =============================================================================
#   flanks.fa              ±50 kb flanks extracted with samtools faidx
#   flanks.fa.fai          index for the flanks
#   self_align.bam         sorted, indexed BAM (load in IGV/JBrowse)
#   self_align.bam.bai
#   self_align.cmd.txt     exact minimap2 + samtools invocation (audit)
#   hits.tsv               parsed alignment table, human-readable
#                          cols: qname tname strand qstart qend tstart tend
#                                alen identity orientation
#   inverted_hits.bed      BED9 with itemRgb. track-colored:
#                            red (255,0,0) = strand '-' (inverted SD)
#                            green (0,128,0) = strand '+' (direct repeat)
#   region.bed             BED6 marking the two breakpoints in blue
#                          (0,0,255) for navigation in IGV.
#
# And also:
#   <outdir>/<candidate_id>/structured/sd_substrate_minimap2.json
#     Tier-2 block with schema sd_substrate_minimap2
#
# =============================================================================
# SCHEMA BLOCK WRITTEN: sd_substrate_minimap2
# =============================================================================
# Flat keys extracted (via keys_extracted spec):
#   q4a_mm2_ran                    bool
#   q4a_mm2_mechanism_left         NAHR_CANDIDATE | NAHR_POSSIBLE |
#                                   COMPLEX_ARCHITECTURE | NHEJ_CANDIDATE | NA
#   q4a_mm2_mechanism_right        (analogous)
#   q4a_mm2_best_inv_length_left   int, bp of best inverted hit at left side
#   q4a_mm2_best_inv_length_right  int, at right side
#   q4a_mm2_best_inv_identity_left  float, %
#   q4a_mm2_best_inv_identity_right float, %
#   q4a_mm2_n_inverted_hits        int, total inverted hits found
#   q4a_mm2_n_direct_hits          int, total direct hits
#   q4a_mm2_preset_used            string, e.g. "asm10"
#   q4a_mm2_bam_path               string, path to the BAM (for audit)
#
# =============================================================================
# CLI
# =============================================================================
#   Rscript sd_substrate_minimap2.R \
#     --candidate LG28_cand_1 \
#     --chrom C_gar_LG28 \
#     --left_bp 15115243 \
#     --right_bp 18005891 \
#     --ref_fasta /path/to/ref.fasta \
#     --outdir /per-candidate/output/root
#     [--mm2_preset asm10]       # default; can use asm5 (stricter) or asm20
#     [--flank_bp 50000]         # default
#     [--overwrite]              # re-run even if BAM exists
#     [--registries_root ...]    # for auto-register the block
#
# REQUIRES: minimap2, samtools on PATH (both needed; we fail fast if absent).
# DIFFICULTY: medium (subprocess coordination + BAM/BED emit)
# DISPATCHER: no (pure sequence analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# --------------------------------------------------------------------------
# Tunables (override via CLI)
# --------------------------------------------------------------------------
DEFAULT_FLANK_BP   <- 50000L
DEFAULT_MM2_PRESET <- "asm10"   # was asm5 in pre-pass-19; asm10 widens to ≥90% id SDs
NAHR_MIN_LEN       <- 500L
NAHR_MIN_IDENT     <- 85
NAHR_POSS_MIN      <- 100L
MIN_ALEN           <- 50L       # post-filter: drop tiny hits
SECONDARY_TO_PRIMARY_RATIO <- 0.1  # minimap2 -p flag
MAX_SECONDARIES    <- 50L           # minimap2 -N flag

# --------------------------------------------------------------------------
# Utilities
# --------------------------------------------------------------------------

check_tool <- function(name) {
  p <- Sys.which(name)
  if (p == "" || !file.exists(p)) {
    stop("[sd_substrate_minimap2] required tool not found on PATH: ", name)
  }
  unname(p)
}

log_msg <- function(...) {
  message("[sd_substrate_mm2] ", ...)
}

# --------------------------------------------------------------------------
# Flank extraction
# --------------------------------------------------------------------------

extract_flanks <- function(ref_fasta, chr, bps, flank_bp, out_fasta) {
  # bps is length-2 (left, right)
  regions <- character()
  for (i in seq_along(bps)) {
    bp <- bps[i]
    ls <- max(1L, bp - flank_bp); le <- bp
    rs <- bp;                     re <- bp + flank_bp
    # Name regions explicitly so we can tell left/right later
    # samtools faidx reads coords from the region string; we can't rename
    # via samtools alone, so we'll run it and then rewrite headers below.
    regions <- c(regions,
                 paste0(chr, ":", ls, "-", le),
                 paste0(chr, ":", rs, "-", re))
  }
  raw_fa <- paste0(out_fasta, ".raw")
  cmd <- paste0("samtools faidx ", shQuote(ref_fasta), " ",
                paste(regions, collapse = " "),
                " > ", shQuote(raw_fa), " 2>/dev/null")
  rc <- system(cmd)
  if (rc != 0 || !file.exists(raw_fa) || file.size(raw_fa) == 0) {
    stop("[sd_substrate_minimap2] samtools faidx failed for ", ref_fasta, ": ", chr)
  }

  # Rewrite headers so downstream code can tell which flank is which.
  # New names: bp1_left, bp1_right, bp2_left, bp2_right
  new_names <- c(
    paste0("bp", 1, "_upstream"),   paste0("bp", 1, "_downstream"),
    paste0("bp", 2, "_upstream"),   paste0("bp", 2, "_downstream")
  )
  lines <- readLines(raw_fa)
  header_idx <- grep("^>", lines)
  if (length(header_idx) != length(new_names)) {
    stop("[sd_substrate_minimap2] unexpected number of flanks extracted: ",
         length(header_idx), " vs expected ", length(new_names))
  }
  lines[header_idx] <- paste0(">", new_names,
                              "  original=", sub("^>", "", lines[header_idx]))
  writeLines(lines, out_fasta)
  unlink(raw_fa)

  # Index the flanks fasta
  rc2 <- system(paste0("samtools faidx ", shQuote(out_fasta)))
  if (rc2 != 0) stop("[sd_substrate_minimap2] samtools faidx index failed")

  invisible(out_fasta)
}

# --------------------------------------------------------------------------
# minimap2 invocation
# --------------------------------------------------------------------------

run_minimap2 <- function(flanks_fa, bam_out, cmd_log, preset) {
  mm2  <- check_tool("minimap2")
  stl  <- check_tool("samtools")
  mm2_cmd <- paste(
    shQuote(mm2),
    "-ax", preset,
    "-X",
    "--eqx",
    "--secondary=yes",
    "-N", MAX_SECONDARIES,
    "-p", SECONDARY_TO_PRIMARY_RATIO,
    shQuote(flanks_fa), shQuote(flanks_fa)
  )
  pipe_cmd <- paste(mm2_cmd, "2>/dev/null |",
                    shQuote(stl), "sort", "-O", "bam",
                    "-o", shQuote(bam_out), "-")
  writeLines(c(
    "# minimap2 self-alignment command (pass-19 audit trail)",
    paste0("# Generated: ", as.character(Sys.time())),
    "",
    mm2_cmd,
    "  |",
    paste(shQuote(stl), "sort -O bam -o", shQuote(bam_out), "-"),
    "",
    paste(shQuote(stl), "index", shQuote(bam_out))
  ), cmd_log)

  rc <- system(pipe_cmd)
  if (rc != 0 || !file.exists(bam_out) || file.size(bam_out) == 0) {
    stop("[sd_substrate_minimap2] minimap2|samtools sort failed (rc=", rc, ")")
  }

  rc2 <- system(paste(shQuote(stl), "index", shQuote(bam_out)))
  if (rc2 != 0) stop("[sd_substrate_minimap2] samtools index failed")

  invisible(bam_out)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a[1]))) b else a

# --------------------------------------------------------------------------
# BAM → hits.tsv (parse via samtools view)
# --------------------------------------------------------------------------

parse_bam_to_hits <- function(bam, hits_tsv) {
  stl <- check_tool("samtools")
  # Get SAM records; skip header with -F 4 (mapped only) + -F 256? No — we
  # WANT secondary (256). We want mapped and (primary OR secondary).
  # -F 4 drops unmapped. That's enough.
  view_cmd <- paste(shQuote(stl), "view", "-F", "4", shQuote(bam))
  con <- pipe(view_cmd, "r")
  sam_lines <- readLines(con)
  close(con)

  if (length(sam_lines) == 0) {
    log_msg("no alignments in BAM")
    hits <- data.table(
      qname = character(), tname = character(), strand = character(),
      qstart = integer(), qend = integer(), tstart = integer(), tend = integer(),
      alen = integer(), identity = numeric(), orientation = character(),
      is_secondary = logical()
    )
    fwrite(hits, hits_tsv, sep = "\t")
    return(hits)
  }

  parse_sam_line <- function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    qname <- f[1]
    flag  <- as.integer(f[2])
    tname <- f[3]
    tstart <- as.integer(f[4])  # 1-based leftmost
    cigar <- f[6]
    seq_str <- f[10]
    is_secondary <- bitwAnd(flag, 256L) > 0L
    is_reverse   <- bitwAnd(flag, 16L) > 0L

    # Decode CIGAR for ref length + match/mismatch counts.
    # With --eqx we have =, X, I, D, S, H, N, M ops.
    op_matches <- gregexpr("[0-9]+[MIDNSHP=X]", cigar, perl = TRUE)[[1]]
    if (op_matches[1] == -1) return(NULL)
    ops_len <- attr(op_matches, "match.length")
    ops <- substring(cigar, op_matches, op_matches + ops_len - 1)
    ops_n   <- as.integer(sub("[MIDNSHP=X]$", "", ops))
    ops_op  <- substring(ops, nchar(ops), nchar(ops))

    # Reference-consuming ops: M = X D N
    ref_consume <- ops_op %in% c("M", "=", "X", "D", "N")
    tlen_on_ref <- sum(ops_n[ref_consume])

    # Query-consuming ops for alignment length: M = X I S H (but H is hard-clip
    # in primary, we should exclude; S and H stripped from SEQ anyway)
    query_align_consume <- ops_op %in% c("M", "=", "X", "I")
    alen <- sum(ops_n[query_align_consume])

    # Match / mismatch counts for identity
    n_match <- sum(ops_n[ops_op == "="])
    n_mism  <- sum(ops_n[ops_op == "X"])
    n_indel <- sum(ops_n[ops_op %in% c("I", "D")])
    identity <- if ((n_match + n_mism + n_indel) > 0) {
      round(100 * n_match / (n_match + n_mism + n_indel), 2)
    } else 0

    # Query coords. Soft-clips at either end are not in alen.
    left_clip  <- if (length(ops_op) > 0 && ops_op[1] %in% c("S", "H")) ops_n[1] else 0L
    right_clip <- if (length(ops_op) > 0 && tail(ops_op, 1) %in% c("S", "H")) tail(ops_n, 1) else 0L
    qstart <- left_clip + 1L
    qend   <- left_clip + alen

    list(
      qname = qname, tname = tname,
      strand = if (is_reverse) "-" else "+",
      qstart = qstart, qend = qend,
      tstart = tstart, tend = tstart + tlen_on_ref - 1L,
      alen = alen, identity = identity,
      orientation = if (is_reverse) "inverted" else "direct",
      is_secondary = is_secondary
    )
  }

  parsed <- lapply(sam_lines, parse_sam_line)
  parsed <- parsed[!sapply(parsed, is.null)]
  hits <- rbindlist(parsed)

  # Post-filters
  if (nrow(hits) > 0) {
    hits <- hits[alen >= MIN_ALEN]
    # Drop symmetric self-hits that survived -X (qname == tname with overlapping coords)
    hits <- hits[qname != tname | (qend <= tstart | tend <= qstart)]
  }

  fwrite(hits, hits_tsv, sep = "\t")
  log_msg("parsed ", nrow(hits), " alignments (",
          sum(hits$orientation == "inverted"), " inverted, ",
          sum(hits$orientation == "direct"), " direct)")
  hits
}

# --------------------------------------------------------------------------
# Classification per breakpoint
# --------------------------------------------------------------------------

classify_side <- function(hits, bp_tag) {
  # bp_tag is "bp1" or "bp2" (left or right breakpoint)
  # An inverted SD relevant to this breakpoint shows up as an inverted
  # alignment where one side of the pair is upstream of bp_tag and the
  # other is downstream.
  if (nrow(hits) == 0) {
    return(list(mechanism = "NHEJ_CANDIDATE", best_len = 0L, best_id = 0,
                n_inv = 0L, n_dir = 0L))
  }
  rel <- hits[grepl(bp_tag, qname) | grepl(bp_tag, tname)]
  inv <- rel[orientation == "inverted"]
  dir <- rel[orientation == "direct"]
  if (nrow(inv) == 0) {
    return(list(mechanism = "NHEJ_CANDIDATE", best_len = 0L, best_id = 0,
                n_inv = 0L, n_dir = nrow(dir)))
  }
  best <- inv[which.max(alen)]
  mech <- if (best$alen[1] >= NAHR_MIN_LEN && best$identity[1] >= NAHR_MIN_IDENT) {
    "NAHR_CANDIDATE"
  } else if (best$alen[1] >= NAHR_POSS_MIN) {
    "NAHR_POSSIBLE"
  } else {
    "COMPLEX_ARCHITECTURE"
  }
  list(mechanism = mech,
       best_len = as.integer(best$alen[1]),
       best_id  = as.numeric(best$identity[1]),
       n_inv    = nrow(inv),
       n_dir    = nrow(dir))
}

# --------------------------------------------------------------------------
# BED writers — track colored by strand
# --------------------------------------------------------------------------

# When extracting flanks, the FASTA coordinates reset to 1 for each flank.
# To convert back to genomic coords we need the flank's origin.
#
# Flank names: bp1_upstream, bp1_downstream, bp2_upstream, bp2_downstream
# upstream   covers (bp - flank, bp)   so genomic = bp - flank + local_coord
# downstream covers (bp, bp + flank)   so genomic = bp + local_coord

flank_to_genomic <- function(flank_name, chr, left_bp, right_bp, flank_bp) {
  if (flank_name == "bp1_upstream") {
    return(list(chr = chr, offset = max(1L, left_bp - flank_bp) - 1L))
  }
  if (flank_name == "bp1_downstream") {
    return(list(chr = chr, offset = left_bp - 1L))
  }
  if (flank_name == "bp2_upstream") {
    return(list(chr = chr, offset = max(1L, right_bp - flank_bp) - 1L))
  }
  if (flank_name == "bp2_downstream") {
    return(list(chr = chr, offset = right_bp - 1L))
  }
  list(chr = NA_character_, offset = NA_integer_)
}

write_inverted_hits_bed <- function(hits, chr, left_bp, right_bp, flank_bp, bed_out) {
  # BED9 with itemRgb — needs `track name=... itemRgb=On` header for UCSC/IGV
  lines <- c(
    "track name=\"sd_substrate_minimap2\" description=\"minimap2 self-alignment hits; red=inverted green=direct\" itemRgb=On"
  )
  if (nrow(hits) > 0) {
    for (i in seq_len(nrow(hits))) {
      row <- hits[i]
      # We emit one BED record for the QUERY side of each hit — that places
      # the repeat where the query's flank maps in the genome.
      q_info <- flank_to_genomic(row$qname, chr, left_bp, right_bp, flank_bp)
      if (is.na(q_info$chr)) next
      start <- q_info$offset + row$qstart - 1L  # BED is 0-based
      end   <- q_info$offset + row$qend
      color <- if (row$orientation == "inverted") "255,0,0" else "0,128,0"
      name  <- sprintf("%s>%s_%s_%.1f%%", row$qname, row$tname,
                       if (row$orientation == "inverted") "inv" else "dir",
                       row$identity)
      score <- as.integer(pmin(1000, round(row$identity * 10)))
      strand <- row$strand
      # BED9: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb
      lines <- c(lines, sprintf("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s",
                                q_info$chr, start, end, name, score, strand,
                                start, end, color))
    }
  }
  writeLines(lines, bed_out)
}

write_region_bed <- function(chr, left_bp, right_bp, flank_bp, bed_out) {
  # BED6 marking the two breakpoints + ±flank zones for IGV nav
  lines <- c(
    "track name=\"sd_substrate_region\" description=\"candidate breakpoints (blue) ±flank\" itemRgb=On",
    sprintf("%s\t%d\t%d\t%s\t0\t.\t%d\t%d\t0,0,255",
            chr, max(0L, left_bp - flank_bp), left_bp + flank_bp,
            "left_bp_window", max(0L, left_bp - flank_bp), left_bp + flank_bp),
    sprintf("%s\t%d\t%d\t%s\t0\t.\t%d\t%d\t0,0,255",
            chr, max(0L, right_bp - flank_bp), right_bp + flank_bp,
            "right_bp_window", max(0L, right_bp - flank_bp), right_bp + flank_bp),
    sprintf("%s\t%d\t%d\t%s\t1000\t.\t%d\t%d\t0,0,255",
            chr, max(0L, left_bp - 1L), left_bp,
            "left_breakpoint", max(0L, left_bp - 1L), left_bp),
    sprintf("%s\t%d\t%d\t%s\t1000\t.\t%d\t%d\t0,0,255",
            chr, max(0L, right_bp - 1L), right_bp,
            "right_breakpoint", max(0L, right_bp - 1L), right_bp)
  )
  writeLines(lines, bed_out)
}

# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

parse_args <- function(argv) {
  opts <- list(
    candidate = NA_character_,
    chrom = NA_character_,
    left_bp = NA_integer_,
    right_bp = NA_integer_,
    ref_fasta = NA_character_,
    outdir = NA_character_,
    mm2_preset = DEFAULT_MM2_PRESET,
    flank_bp = DEFAULT_FLANK_BP,
    overwrite = FALSE,
    registries_root = NA_character_
  )
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]; v <- if (i < length(argv)) argv[i + 1] else NA_character_
    switch(a,
      "--candidate"       = { opts$candidate <- v; i <- i + 2 },
      "--chrom"           = { opts$chrom <- v; i <- i + 2 },
      "--left_bp"         = { opts$left_bp <- as.integer(v); i <- i + 2 },
      "--right_bp"        = { opts$right_bp <- as.integer(v); i <- i + 2 },
      "--ref_fasta"       = { opts$ref_fasta <- v; i <- i + 2 },
      "--outdir"          = { opts$outdir <- v; i <- i + 2 },
      "--mm2_preset"      = { opts$mm2_preset <- v; i <- i + 2 },
      "--flank_bp"        = { opts$flank_bp <- as.integer(v); i <- i + 2 },
      "--overwrite"       = { opts$overwrite <- TRUE; i <- i + 1 },
      "--registries_root" = { opts$registries_root <- v; i <- i + 2 },
      { i <- i + 1 }
    )
  }
  opts
}

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  opts <- parse_args(argv)
  for (k in c("candidate", "chrom", "left_bp", "right_bp", "ref_fasta", "outdir")) {
    if (is.na(opts[[k]])) stop("[sd_substrate_minimap2] missing required --", k)
  }

  if (!opts$mm2_preset %in% c("asm5", "asm10", "asm20")) {
    log_msg("warning: --mm2_preset '", opts$mm2_preset,
            "' is not asm5/asm10/asm20; passing through anyway")
  }

  # Per-candidate layout
  cand_dir <- file.path(opts$outdir, opts$candidate)
  mm2_dir  <- file.path(cand_dir, "sd_substrate", "minimap2")
  str_dir  <- file.path(cand_dir, "structured")
  dir.create(mm2_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(str_dir, recursive = TRUE, showWarnings = FALSE)

  flanks_fa  <- file.path(mm2_dir, "flanks.fa")
  bam_out    <- file.path(mm2_dir, "self_align.bam")
  cmd_log    <- file.path(mm2_dir, "self_align.cmd.txt")
  hits_tsv   <- file.path(mm2_dir, "hits.tsv")
  inv_bed    <- file.path(mm2_dir, "inverted_hits.bed")
  region_bed <- file.path(mm2_dir, "region.bed")

  # Skip if BAM already exists (unless --overwrite)
  if (file.exists(bam_out) && !opts$overwrite) {
    log_msg("BAM already exists, skipping minimap2 (use --overwrite to re-run): ", bam_out)
  } else {
    log_msg(opts$candidate, " — extracting ±", opts$flank_bp, " bp flanks for ",
            opts$chrom, ":", opts$left_bp, " and ", opts$chrom, ":", opts$right_bp)
    extract_flanks(opts$ref_fasta, opts$chrom,
                   c(opts$left_bp, opts$right_bp),
                   opts$flank_bp, flanks_fa)

    log_msg(opts$candidate, " — running minimap2 -ax ", opts$mm2_preset,
            " (self-alignment)")
    run_minimap2(flanks_fa, bam_out, cmd_log, opts$mm2_preset)
  }

  # Region BED for navigation
  write_region_bed(opts$chrom, opts$left_bp, opts$right_bp,
                   opts$flank_bp, region_bed)

  # Parse BAM → hits.tsv (always re-parses)
  hits <- parse_bam_to_hits(bam_out, hits_tsv)

  # BED overlay colored by strand
  write_inverted_hits_bed(hits, opts$chrom, opts$left_bp, opts$right_bp,
                          opts$flank_bp, inv_bed)

  # Per-breakpoint classification
  left  <- classify_side(hits, "bp1")
  right <- classify_side(hits, "bp2")

  # Structured block
  block <- list(
    block_type    = "sd_substrate_minimap2",
    candidate_id  = opts$candidate,
    source_script = "sd_substrate_minimap2.R",
    data = list(
      q4a_mm2_ran                     = nrow(hits) >= 0,  # ran successfully
      q4a_mm2_preset_used             = opts$mm2_preset,
      q4a_mm2_flank_bp                = opts$flank_bp,
      q4a_mm2_n_inverted_hits         = sum(hits$orientation == "inverted"),
      q4a_mm2_n_direct_hits           = sum(hits$orientation == "direct"),
      q4a_mm2_mechanism_left          = left$mechanism,
      q4a_mm2_mechanism_right         = right$mechanism,
      q4a_mm2_best_inv_length_left    = left$best_len,
      q4a_mm2_best_inv_length_right   = right$best_len,
      q4a_mm2_best_inv_identity_left  = left$best_id,
      q4a_mm2_best_inv_identity_right = right$best_id,
      q4a_mm2_bam_path                = normalizePath(bam_out, mustWork = FALSE),
      q4a_mm2_bed_path                = normalizePath(inv_bed, mustWork = FALSE),
      q4a_mm2_cmd_log                 = normalizePath(cmd_log, mustWork = FALSE)
    )
  )
  out_json <- file.path(str_dir, "sd_substrate_minimap2.json")
  writeLines(toJSON(block, auto_unbox = TRUE, null = "null", na = "null",
                    pretty = TRUE), out_json)
  log_msg("wrote ", out_json)

  # Optional registry write
  if (!is.na(opts$registries_root)) {
    reg_r <- file.path(opts$registries_root, "api", "R", "registry_loader.R")
    if (file.exists(reg_r)) {
      tryCatch({
        source(reg_r)
        reg <- load_registry()
        reg$evidence$write_block(
          candidate_id  = opts$candidate,
          block_type    = "sd_substrate_minimap2",
          data          = block$data,
          source_script = "sd_substrate_minimap2.R"
        )
        log_msg("registered block for ", opts$candidate)
      }, error = function(e) {
        log_msg("registry write failed: ", e$message)
      })
    }
  }

  log_msg(opts$candidate,
          " — mechanism_left=", left$mechanism,
          " mechanism_right=", right$mechanism,
          " n_inv=", sum(hits$orientation == "inverted"))
}

# Only run main if invoked as Rscript
invoked_as_script <- any(grepl("--file=",
                               commandArgs(trailingOnly = FALSE)))
if (invoked_as_script) main()
