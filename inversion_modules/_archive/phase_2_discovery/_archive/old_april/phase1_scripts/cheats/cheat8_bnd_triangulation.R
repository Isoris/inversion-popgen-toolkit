#!/usr/bin/env Rscript
# =============================================================================
# cheat8_bnd_triangulation.R — BND paired breakpoint triangulation
#
# DELLY BND records with CT=3to3 and CT=5to5 define inversion boundaries
# as PAIRED positions. Even without a direct INV call, paired BNDs
# triangulate the inversion boundaries.
#
# Extends the flashlight RDS with triangulated inversions from BND pairs.
# Already partially implemented in MODULE_5A2 STEP06, but not wired as
# a flashlight prior.
#
# REQUIRES: flashlight_loader.R sourced, DELLY BND VCF available
# =============================================================================

suppressPackageStartupMessages(library(data.table))

MIN_INV_SIZE    <- 5000L
MAX_INV_SIZE    <- 20e6
MIN_CARRIER_OVL <- 0.30  # minimum Jaccard overlap between left+right BND carriers
MIN_CARRIERS    <- 3L

#' Parse BND VCF and extract paired breakpoints
#' @param bnd_vcf_path Path to DELLY BND VCF (catalog_226.BND.vcf.gz)
#' @param sample_ids CGA sample names in VCF order
#' @return data.table of candidate inversion boundaries from BND pairs
parse_bnd_pairs <- function(bnd_vcf_path, sample_ids = NULL) {
  if (!file.exists(bnd_vcf_path)) {
    message("[cheat8] BND VCF not found: ", bnd_vcf_path)
    return(data.table())
  }

  message("[cheat8] Parsing BND VCF: ", bnd_vcf_path)

  con <- if (grepl("\\.gz$", bnd_vcf_path)) gzfile(bnd_vcf_path, "rt") else file(bnd_vcf_path, "rt")
  on.exit(close(con))

  vcf_samples <- character(0)
  bnds <- list()
  n <- 0L

  while (TRUE) {
    line <- readLines(con, n = 1L)
    if (length(line) == 0) break
    if (startsWith(line, "##")) next

    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (startsWith(line, "#CHROM")) {
      vcf_samples <- fields[10:length(fields)]
      next
    }

    # Parse INFO for CT, CHR2, POS2
    info_str <- fields[8]
    info_pairs <- strsplit(info_str, ";", fixed = TRUE)[[1]]
    info <- list()
    for (ip in info_pairs) {
      if (grepl("=", ip, fixed = TRUE)) {
        kv <- strsplit(ip, "=", fixed = TRUE)[[1]]
        info[[kv[1]]] <- kv[2]
      }
    }

    ct   <- info[["CT"]]
    chr2 <- info[["CHR2"]]
    pos2 <- as.integer(info[["POS2"]] %||% info[["END"]] %||% "0")
    if (is.null(ct) || is.na(pos2)) next

    chrom <- fields[1]
    pos1  <- as.integer(fields[2])

    # Only process same-chromosome BNDs with inversion orientations
    if (chrom != chr2) next
    if (!ct %in% c("3to3", "5to5")) next

    # Parse per-sample GT
    format_fields <- strsplit(fields[9], ":", fixed = TRUE)[[1]]
    gt_idx <- match("GT", format_fields)
    carriers <- character(0)
    if (!is.na(gt_idx)) {
      for (si in seq_along(vcf_samples)) {
        sdata <- strsplit(fields[9 + si], ":", fixed = TRUE)[[1]]
        gt <- if (gt_idx <= length(sdata)) sdata[gt_idx] else "./."
        if (gt %in% c("0/1", "0|1", "1|0", "1/1", "1|1")) {
          carriers <- c(carriers, vcf_samples[si])
        }
      }
    }

    n <- n + 1L
    bnds[[n]] <- data.table(
      bnd_id = fields[3], chrom = chrom, pos = pos1,
      ct = ct, chr2 = chr2, pos2 = pos2,
      n_carriers = length(carriers),
      carriers = list(carriers)
    )
  }

  if (n == 0) return(data.table())
  bnd_dt <- rbindlist(bnds)
  message("[cheat8] Parsed ", nrow(bnd_dt), " same-chr inversion-oriented BNDs")

  # Split by orientation
  left_bps  <- bnd_dt[ct == "3to3"]   # left inversion breakpoints
  right_bps <- bnd_dt[ct == "5to5"]   # right inversion breakpoints

  message("[cheat8] Left (3to3): ", nrow(left_bps), " | Right (5to5): ", nrow(right_bps))

  # Pair matching
  pairs <- list()
  for (li in seq_len(nrow(left_bps))) {
    lb <- left_bps[li]
    # Find matching right breakpoints on same chromosome
    candidates <- right_bps[
      chrom == lb$chrom &
      pos > lb$pos &
      (pos - lb$pos) >= MIN_INV_SIZE &
      (pos - lb$pos) <= MAX_INV_SIZE
    ]
    if (nrow(candidates) == 0) next

    for (ri in seq_len(nrow(candidates))) {
      rb <- candidates[ri]
      # Check carrier overlap
      l_carr <- lb$carriers[[1]]
      r_carr <- rb$carriers[[1]]
      shared <- length(intersect(l_carr, r_carr))
      union_n <- length(union(l_carr, r_carr))
      jaccard <- if (union_n > 0) shared / union_n else 0

      if (jaccard < MIN_CARRIER_OVL || shared < MIN_CARRIERS) next

      pairs[[length(pairs) + 1]] <- data.table(
        chrom = lb$chrom,
        bp_left = lb$pos, bp_right = rb$pos,
        svlen = rb$pos - lb$pos,
        bnd_left_id = lb$bnd_id, bnd_right_id = rb$bnd_id,
        n_carriers_left = lb$n_carriers,
        n_carriers_right = rb$n_carriers,
        n_shared_carriers = shared,
        carrier_jaccard = round(jaccard, 3),
        shared_carriers = list(intersect(l_carr, r_carr)),
        source = "BND_triangulation"
      )
    }
  }

  if (length(pairs) == 0) {
    message("[cheat8] No matching BND pairs found")
    return(data.table())
  }

  pairs_dt <- rbindlist(pairs)
  message("[cheat8] Triangulated ", nrow(pairs_dt), " candidate inversions from BND pairs")
  pairs_dt
}

#' Add BND-triangulated inversions to the flashlight RDS
#' @param fl Flashlight object (from load_flashlight)
#' @param bnd_pairs data.table from parse_bnd_pairs
#' @return Updated flashlight with bnd_triangulated field
add_bnd_to_flashlight <- function(fl, bnd_pairs) {
  if (is.null(fl)) return(fl)
  fl$bnd_triangulated <- bnd_pairs
  fl
}
