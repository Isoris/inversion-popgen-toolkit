#!/usr/bin/env Rscript

# =============================================================================
# PLOT_10_context_enrichment.R
#
# Variant fold-enrichment in genomic contexts. Image 3A.
#   contexts: genes (CDS/exon/intron), repeats, inversion intervals, ROH
#   variant classes: SNP, INDEL, DEL/DUP/INV/BND, PAV
# Fold enrichment = (variants overlapping context / total variants) /
#                    (context bp / genome bp)
#
# Skips contexts that don't exist on disk.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

VARIANT_MASTER  <- Sys.getenv("VARIANT_MASTER")
GENE_BED        <- Sys.getenv("GENE_BED")
REPEAT_BED      <- Sys.getenv("REPEAT_BED", "")
ROH_BED         <- Sys.getenv("ROH_BED", "")
SNAKE_CAND_FILE <- Sys.getenv("SNAKE_CAND_FILE")
REF_FAI         <- Sys.getenv("REF_FAI")
EXTRAS_FIG_DIR  <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR  <- Sys.getenv("EXTRAS_TBL_DIR")

if (!file.exists(VARIANT_MASTER) || !file.exists(REF_FAI)) {
  message("[PLOT_10] [skip] missing VARIANT_MASTER or REF_FAI"); quit(status = 0)
}

v <- fread(VARIANT_MASTER)
chr_col <- if ("chr" %in% names(v)) "chr" else "chrom"
classify <- function(v) {
  cls <- rep("OTHER", nrow(v))
  ref_len <- if ("ref" %in% names(v)) nchar(v$ref) else NA
  alt_len <- if ("alt" %in% names(v)) nchar(v$alt) else NA
  if (!is.na(ref_len[1])) {
    cls[ref_len == 1 & alt_len == 1] <- "SNP"
    cls[ref_len != alt_len & pmax(ref_len, alt_len) <= 50] <- "INDEL"
  }
  if ("sv_type" %in% names(v)) {
    sv <- toupper(v$sv_type)
    for (k in c("DEL","DUP","INV","BND","INS")) cls[sv == k] <- k
  }
  if ("variant_class" %in% names(v))
    cls[toupper(v$variant_class) == "PAV"] <- "PAV"
  cls
}
v[, class := classify(.SD)]

# Genome length
fai <- fread(REF_FAI, header = FALSE,
             col.names = c("chrom", "len", "off", "x1", "x2"))
GENOME_LEN <- sum(as.numeric(fai$len))

# Read context BEDs (try gz if uncompressed missing)
read_bed <- function(path) {
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) return(NULL)
  rd <- if (endsWith(path, ".gz")) {
    fread(cmd = paste0("zcat ", shQuote(path)), header = FALSE)
  } else {
    fread(path, header = FALSE)
  }
  if (ncol(rd) < 3) return(NULL)
  setnames(rd, 1:3, c("chrom", "start", "end"))
  rd[, .(chrom, start = as.integer(start), end = as.integer(end))]
}

# Inversion candidate intervals from SNAKE_CAND_FILE
read_cands_as_bed <- function() {
  if (!file.exists(SNAKE_CAND_FILE)) return(NULL)
  c1 <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))
  c1[, .(chrom, start = as.integer(start_bp), end = as.integer(end_bp))]
}

contexts <- list(
  Genes      = read_bed(GENE_BED),
  Repeats    = read_bed(REPEAT_BED),
  ROH        = read_bed(ROH_BED),
  Inversions = read_cands_as_bed())
contexts <- contexts[!vapply(contexts, is.null, logical(1))]
if (length(contexts) == 0) {
  message("[PLOT_10] [skip] no context BEDs available"); quit(status = 0)
}
message("[PLOT_10] Contexts available: ", paste(names(contexts), collapse = ", "))

# Compute total context bp (with simple overlap merge for repeats)
context_bp <- vapply(contexts, function(b) sum(b$end - b$start),
                      numeric(1))

# For each variant, mark whether it falls in each context
in_context <- function(varv, b) {
  setkey(b, chrom, start, end)
  out <- foverlaps(varv[, .(chrom = get(chr_col), start = pos, end = pos)],
                    b, type = "within", which = TRUE, mult = "first")
  !is.na(out$yid)
}
v[, idx := .I]
contexts_long <- list()
for (ctx_name in names(contexts)) {
  v[, paste0("in_", ctx_name) := in_context(.SD, contexts[[ctx_name]])]
}

# Aggregate
rows <- list()
for (ctx_name in names(contexts)) {
  col <- paste0("in_", ctx_name)
  for (cls in unique(v$class)) {
    n_in <- sum(v[class == cls, get(col)])
    n_tot <- sum(v$class == cls)
    if (n_tot == 0) next
    obs_frac <- n_in / n_tot
    exp_frac <- context_bp[ctx_name] / GENOME_LEN
    rows[[length(rows) + 1]] <- data.table(
      context = ctx_name,
      class = cls,
      n_overlapping = n_in,
      n_total = n_tot,
      obs_frac = obs_frac,
      exp_frac = exp_frac,
      fold_enrichment = obs_frac / max(exp_frac, 1e-9))
  }
}
agg <- rbindlist(rows)
fwrite(agg, file.path(EXTRAS_TBL_DIR, "PLOT_10_context_enrichment.tsv"), sep = "\t")

CLASS_PAL <- c(SNP = "#1f77b4", INDEL = "#2ca02c",
                DEL = "#d62728", DUP = "#ff7f0e",
                INV = "#9467bd", INS = "#e377c2",
                BND = "#8c564b", PAV = "#7f7f7f", OTHER = "grey80")

p <- ggplot(agg, aes(x = context, y = fold_enrichment, fill = class)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = CLASS_PAL) +
  scale_y_continuous(trans = "log2",
                       breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
                       labels = c("0.25x", "0.5x", "1x", "2x", "4x", "8x", "16x")) +
  labs(title = "Variant fold-enrichment by genomic context",
       y = "Fold enrichment (obs/expected, log2)",
       x = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_10_context_enrichment.pdf"),
       p, width = 9, height = 6, device = cairo_pdf)
message("[PLOT_10] Done")
