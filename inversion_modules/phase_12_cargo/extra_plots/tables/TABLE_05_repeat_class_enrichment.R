#!/usr/bin/env Rscript

# =============================================================================
# TABLE_05_repeat_class_enrichment.R
#
# Image 6C — repeat-class breakdown of inversion vs genome:
#   class (LINE/SINE/LTR/DNA/SimpleRepeat/LowComplexity), n_in_inversions,
#   pct_inv_seq, pct_genome_seq, fold_enrichment
#
# Requires REPEAT_BED with a 4th column = repeat class. Falls back gracefully
# if the BED is unavailable or has no class column.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

REPEAT_BED      <- Sys.getenv("REPEAT_BED", "")
SNAKE_CAND_FILE <- Sys.getenv("SNAKE_CAND_FILE")
REF_FAI         <- Sys.getenv("REF_FAI")
EXTRAS_TBL_DIR  <- Sys.getenv("EXTRAS_TBL_DIR")

if (!nzchar(REPEAT_BED) || !file.exists(REPEAT_BED)) {
  message("[TABLE_05] [skip] REPEAT_BED not available — set REPEAT_BED env var")
  quit(status = 0)
}
if (!file.exists(SNAKE_CAND_FILE) || !file.exists(REF_FAI)) {
  message("[TABLE_05] [skip] missing inputs"); quit(status = 0)
}

# Read repeats — try column 4 as class, else mark "unknown"
rd <- if (endsWith(REPEAT_BED, ".gz")) {
  fread(cmd = paste0("zcat ", shQuote(REPEAT_BED)), header = FALSE)
} else {
  fread(REPEAT_BED, header = FALSE)
}
if (ncol(rd) < 3) {
  message("[TABLE_05] [skip] REPEAT_BED has <3 columns"); quit(status = 0)
}
setnames(rd, 1:3, c("chrom", "start", "end"))
if (ncol(rd) >= 4) {
  setnames(rd, 4, "rep_class")
} else {
  rd[, rep_class := "unknown"]
}
rd[, start := as.integer(start)]; rd[, end := as.integer(end)]
rd[, len := end - start]

# Normalize class names — RepeatMasker can emit "LINE/L1" etc; collapse to top level
rd[, class_top := sub("/.*", "", rep_class)]
# Group "Simple_repeat", "simple_repeat" → SimpleRepeat
rd[, class_top := fcase(
  grepl("^LINE",  class_top, ignore.case = TRUE), "LINE",
  grepl("^SINE",  class_top, ignore.case = TRUE), "SINE",
  grepl("^LTR",   class_top, ignore.case = TRUE), "LTR",
  grepl("^DNA",   class_top, ignore.case = TRUE), "DNA",
  grepl("simple", class_top, ignore.case = TRUE), "SimpleRepeat",
  grepl("low.?complex", class_top, ignore.case = TRUE), "LowComplexity",
  default = class_top)]

# Genome length
fai <- fread(REF_FAI, header = FALSE,
             col.names = c("chrom", "len", "off", "x1", "x2"))
GENOME_LEN <- sum(as.numeric(fai$len))

# Inversion intervals
cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))[
  , .(chrom, start = as.integer(start_bp), end = as.integer(end_bp))]
INV_LEN <- sum(as.numeric(cands$end - cands$start))

# Overlap each repeat with inversion intervals → bp inside inversions
setkey(cands, chrom, start, end)
ov <- foverlaps(rd[, .(chrom, start, end, len, class_top)],
                 cands, type = "any")
ov[, ov_start := pmax(start, i.start)]
ov[, ov_end   := pmin(end, i.end)]
ov[, ov_len   := pmax(0, ov_end - ov_start)]
in_inv <- ov[!is.na(i.start), .(bp_inside = sum(ov_len)), by = class_top]

# Total bp per class genome-wide
total_bp <- rd[, .(bp_genome = sum(len)), by = class_top]

agg <- merge(total_bp, in_inv, by = "class_top", all.x = TRUE)
agg[is.na(bp_inside), bp_inside := 0]
agg[, pct_genome := round(100 * bp_genome / GENOME_LEN, 3)]
agg[, pct_inv := round(100 * bp_inside / INV_LEN, 3)]
agg[, fold_enrichment := round((bp_inside / INV_LEN) / (bp_genome / GENOME_LEN), 3)]
agg <- agg[order(-fold_enrichment)]
setnames(agg, "class_top", "repeat_class")

fwrite(agg, file.path(EXTRAS_TBL_DIR, "TABLE_05_repeat_class_enrichment.tsv"),
       sep = "\t")
message("[TABLE_05] Done")
print(agg)
