#!/usr/bin/env Rscript
# =============================================================================
# q09_gap_characterization.R
# =============================================================================
# Detect contiguous zero-SNP / low-density regions on one chromosome from
# the precomp and emit a features table + BED.
#
# A "gap" is a run of precomp windows where either n_snps == 0, or where the
# bp-span the window had to cover to collect the fixed SNP count is above a
# liberal threshold (3x the 99.5th percentile or 80 kb, whichever is larger).
# Adjacent gap windows within 50 kb of each other are merged into a region.
#
# For each gap region we compute:
#   - span_kb
#   - frac_N              (from reference FASTA; counts any of N, n)
#   - frac_softmasked     (lowercase letters — typically repeat-masked)
#   - frac_uppercase      (1 - frac_N - frac_softmasked)
#   - frac_repeat_bed     (bp overlap with REPEAT_BED, if provided)
#   - cov_mean            (across-sample mean mosdepth coverage)
#   - cov_n_low           (n samples with < 0.5 x global-mean coverage here)
#
# Output TSV columns:
#   chrom  start_bp  end_bp  span_kb  n_windows
#   frac_N  frac_softmasked  frac_uppercase
#   frac_repeat_bed  cov_mean  cov_n_low  flag
#
# Output BED: chrom  start_bp  end_bp  name(=gap_id)  score(=span_kb)  strand(=.)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)[1]
  if (is.na(i) || i >= length(args)) return(default)
  args[i + 1]
}

PRECOMP    <- get_arg("--precomp")
CHROM      <- get_arg("--chrom")
OUT_TSV    <- get_arg("--out_tsv")
OUT_BED    <- get_arg("--out_bed")
REFERENCE  <- get_arg("--reference",   NULL)
REPEAT_BED <- get_arg("--repeat_bed",  NULL)
MOSDEPTH_DIR <- get_arg("--mosdepth_dir", NULL)

stopifnot(!is.null(PRECOMP), file.exists(PRECOMP),
          !is.null(CHROM), !is.null(OUT_TSV), !is.null(OUT_BED))

message("[q09] ", CHROM, ": loading precomp")
pc <- readRDS(PRECOMP)
dt <- as.data.table(pc$dt)
dt[, span_kb := (end_bp - start_bp) / 1e3]

# -- Detect gap windows --------------------------------------------------------
dt[, is_gap := FALSE]
if ("n_snps" %in% names(dt)) dt[n_snps == 0, is_gap := TRUE]
span_hi <- stats::quantile(dt$span_kb, 0.995, na.rm = TRUE)
dt[span_kb > max(span_hi * 3, 80), is_gap := TRUE]

g <- dt[is_gap == TRUE][order(start_bp)]
if (nrow(g) == 0) {
  message("[q09] ", CHROM, ": no gap windows detected")
  # Still emit empty headers so the ALL stitch doesn't error
  fwrite(data.table(chrom=character(), start_bp=integer(), end_bp=integer(),
                    span_kb=numeric(), n_windows=integer(),
                    frac_N=numeric(), frac_softmasked=numeric(),
                    frac_uppercase=numeric(), frac_repeat_bed=numeric(),
                    cov_mean=numeric(), cov_n_low=integer(),
                    flag=character()),
         OUT_TSV, sep = "\t")
  writeLines(character(), OUT_BED)
  quit(save = "no", status = 0)
}

# Merge consecutive gap windows within 50 kb into regions
g[, run_grp := cumsum(c(TRUE, diff(start_bp) > 50e3))]
regions <- g[, .(chrom     = CHROM,
                 start_bp  = min(start_bp),
                 end_bp    = max(end_bp),
                 n_windows = .N),
             by = run_grp]
regions[, span_kb := (end_bp - start_bp) / 1e3]
regions[, run_grp := NULL]
message("[q09] ", CHROM, ": ", nrow(regions), " gap regions")

# -- N-content / softmask --------------------------------------------------
regions[, frac_N           := NA_real_]
regions[, frac_softmasked  := NA_real_]
regions[, frac_uppercase   := NA_real_]

extract_fasta_substr <- function(ref_path, chrom, start_bp, end_bp) {
  # Uses samtools faidx if available; falls back to reading the whole chrom
  # once via Biostrings if not (slower but self-contained).
  cmd <- sprintf("samtools faidx %s %s:%d-%d 2>/dev/null | tail -n +2 | tr -d '\\n'",
                 shQuote(ref_path), chrom, start_bp + 1L, end_bp)
  tryCatch({
    out <- system(cmd, intern = TRUE)
    if (length(out) == 0) return(NA_character_)
    paste(out, collapse = "")
  }, error = function(e) NA_character_)
}

if (!is.null(REFERENCE) && file.exists(REFERENCE)) {
  message("[q09] computing N-content and softmask for ", nrow(regions), " regions")
  for (i in seq_len(nrow(regions))) {
    seq <- extract_fasta_substr(REFERENCE, CHROM,
                                regions$start_bp[i], regions$end_bp[i])
    if (is.na(seq) || nchar(seq) == 0) next
    nN   <- sum(strsplit(seq, "")[[1]] %in% c("N", "n"))
    nLow <- sum(grepl("[a-z]", strsplit(seq, "")[[1]]))
    total <- nchar(seq)
    regions$frac_N[i]          <- nN   / total
    regions$frac_softmasked[i] <- nLow / total
    regions$frac_uppercase[i]  <- 1 - regions$frac_N[i] - regions$frac_softmasked[i]
  }
}

# -- Repeat BED overlap ----------------------------------------------------
regions[, frac_repeat_bed := NA_real_]
if (!is.null(REPEAT_BED) && file.exists(REPEAT_BED)) {
  message("[q09] computing repeat BED overlap")
  rb <- tryCatch(fread(REPEAT_BED, header = FALSE,
                       col.names = c("chrom", "start", "end"),
                       select = 1:3),
                 error = function(e) NULL)
  if (!is.null(rb)) {
    rb <- rb[chrom == CHROM]
    setorder(rb, start)
    for (i in seq_len(nrow(regions))) {
      a <- regions$start_bp[i]; b <- regions$end_bp[i]
      hits <- rb[end > a & start < b]
      if (nrow(hits) == 0) { regions$frac_repeat_bed[i] <- 0; next }
      ovl <- sum(pmin(hits$end, b) - pmax(hits$start, a))
      regions$frac_repeat_bed[i] <- ovl / (b - a)
    }
  }
}

# -- Mosdepth coverage stats -----------------------------------------------
regions[, cov_mean  := NA_real_]
regions[, cov_n_low := NA_integer_]
if (!is.null(MOSDEPTH_DIR) && dir.exists(MOSDEPTH_DIR)) {
  mos_files <- list.files(MOSDEPTH_DIR, pattern = "regions\\.bed\\.gz$",
                          recursive = TRUE, full.names = TRUE)
  if (length(mos_files) > 0) {
    message("[q09] scanning mosdepth across ", length(mos_files), " samples")
    # For each region, stream mosdepth files and collect per-sample mean cov
    for (i in seq_len(nrow(regions))) {
      a <- regions$start_bp[i]; b <- regions$end_bp[i]
      per_sample <- rep(NA_real_, length(mos_files))
      global_cov <- rep(NA_real_, length(mos_files))
      for (j in seq_along(mos_files)) {
        tb <- tryCatch(fread(cmd = sprintf("zcat %s", shQuote(mos_files[j])),
                             col.names = c("chrom", "start", "end", "cov")),
                       error = function(e) NULL)
        if (is.null(tb) || !nrow(tb)) next
        tb_c <- tb[chrom == CHROM]
        if (!nrow(tb_c)) next
        global_cov[j] <- mean(tb_c$cov, na.rm = TRUE)
        hits <- tb_c[end > a & start < b]
        if (nrow(hits) > 0) {
          # Weighted mean coverage across overlapping bins
          bin_bp <- pmin(hits$end, b) - pmax(hits$start, a)
          per_sample[j] <- sum(hits$cov * bin_bp, na.rm = TRUE) / sum(bin_bp)
        }
      }
      regions$cov_mean[i]  <- mean(per_sample, na.rm = TRUE)
      low_threshold <- mean(global_cov, na.rm = TRUE) * 0.5
      regions$cov_n_low[i] <- sum(per_sample < low_threshold, na.rm = TRUE)
    }
  }
}

# -- Flag / classify ---------------------------------------------------------
regions[, flag := {
  f <- character(.N)
  for (i in seq_len(.N)) {
    if (!is.na(frac_N[i]) && frac_N[i] > 0.5) {
      f[i] <- "ASSEMBLY_GAP"
    } else if (!is.na(frac_softmasked[i]) && frac_softmasked[i] > 0.5) {
      f[i] <- "REPEAT_RICH"
    } else if (!is.na(cov_mean[i]) && !is.na(cov_n_low[i]) && cov_n_low[i] > 50) {
      f[i] <- "LOW_MAPPABILITY"
    } else {
      f[i] <- "DROPOUT"
    }
  }
  f
}]

# -- Write outputs ------------------------------------------------------------
setcolorder(regions, c("chrom", "start_bp", "end_bp", "span_kb", "n_windows",
                       "frac_N", "frac_softmasked", "frac_uppercase",
                       "frac_repeat_bed", "cov_mean", "cov_n_low", "flag"))
fwrite(regions, OUT_TSV, sep = "\t", na = "NA")

bed <- regions[, .(chrom, start_bp, end_bp,
                   name = sprintf("%s_gap%d", CHROM, seq_len(.N)),
                   score = round(span_kb, 1),
                   strand = ".")]
fwrite(bed, OUT_BED, sep = "\t", col.names = FALSE)

message("[q09] wrote ", OUT_TSV)
message("[q09] wrote ", OUT_BED)
