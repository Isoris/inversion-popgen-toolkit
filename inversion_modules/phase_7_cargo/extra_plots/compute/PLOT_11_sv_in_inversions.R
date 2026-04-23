#!/usr/bin/env Rscript

# =============================================================================
# PLOT_11_sv_in_inversions.R
#
# log2 fold-enrichment of large SVs (DEL/DUP/INV/BND/INS) inside vs outside
# inversion candidate intervals. Image 3D.
#
# Reads DELLY + Manta VCFs and the inversion candidate table; bins SVs by
# whether they overlap any candidate; computes density inside vs outside.
#
# Skips per-VCF if the file isn't present.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

DELLY_DEL <- Sys.getenv("DELLY_DEL_VCF")
DELLY_DUP <- Sys.getenv("DELLY_DUP_VCF")
DELLY_INV <- Sys.getenv("DELLY_INV_VCF")
DELLY_BND <- Sys.getenv("DELLY_BND_VCF")
MANTA     <- Sys.getenv("MANTA_VCF")
SNAKE     <- Sys.getenv("SNAKE_CAND_FILE")
REF_FAI   <- Sys.getenv("REF_FAI")
EXTRAS_FIG_DIR <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR <- Sys.getenv("EXTRAS_TBL_DIR")

if (!file.exists(SNAKE) || !file.exists(REF_FAI)) {
  message("[PLOT_11] [skip] missing inputs"); quit(status = 0)
}

# Read SV positions from a VCF (CHROM, POS, INFO/SVTYPE, INFO/END if present)
read_sv_vcf <- function(path, sv_type_default = NA) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  cmd <- paste0("zcat ", shQuote(path),
                 " | grep -v '^#' | awk -F'\\t' 'BEGIN{OFS=\"\\t\"} {print $1, $2, $8}'")
  rd <- tryCatch(fread(cmd = cmd, header = FALSE,
                        col.names = c("chrom", "pos", "info"),
                        sep = "\t"),
                 error = function(e) NULL)
  if (is.null(rd) || nrow(rd) == 0) return(NULL)
  rd[, svtype := {
    m <- regmatches(info, regexpr("SVTYPE=[A-Z]+", info))
    sub("SVTYPE=", "", m)
  }]
  rd[svtype == "" | is.na(svtype), svtype := sv_type_default]
  rd[, end := {
    m <- regmatches(info, regexpr("END=[0-9]+", info))
    suppressWarnings(as.integer(sub("END=", "", m)))
  }]
  rd[is.na(end), end := pos + 1L]
  rd[, .(chrom, start = pos, end, svtype)]
}

svs <- rbindlist(list(
  read_sv_vcf(DELLY_DEL, "DEL"),
  read_sv_vcf(DELLY_DUP, "DUP"),
  read_sv_vcf(DELLY_INV, "INV"),
  read_sv_vcf(DELLY_BND, "BND"),
  read_sv_vcf(MANTA, NA)),
  fill = TRUE)

if (nrow(svs) == 0) {
  message("[PLOT_11] [skip] no SV calls available"); quit(status = 0)
}
message("[PLOT_11] Loaded ", nrow(svs), " SV records")

cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE)))[
  , .(chrom, start = as.integer(start_bp), end = as.integer(end_bp))]
fai <- fread(REF_FAI, header = FALSE,
             col.names = c("chrom", "len", "off", "x1", "x2"))
GENOME_LEN <- sum(as.numeric(fai$len))
INV_LEN <- sum(as.numeric(cands$end - cands$start))
NONINV_LEN <- GENOME_LEN - INV_LEN
if (INV_LEN <= 0 || NONINV_LEN <= 0) {
  message("[PLOT_11] [skip] no candidate intervals"); quit(status = 0)
}

# Mark SVs overlapping any candidate
setkey(cands, chrom, start, end)
svs[, idx := .I]
ov <- foverlaps(svs[, .(chrom, start, end, idx)], cands,
                 type = "any", which = TRUE, mult = "first")
svs[, in_inv := !is.na(ov$yid)]

agg <- svs[!is.na(svtype), .(
  n_inside = sum(in_inv),
  n_outside = sum(!in_inv)),
  by = svtype]
agg[, density_inside := n_inside / (INV_LEN / 1e6)]      # per Mb
agg[, density_outside := n_outside / (NONINV_LEN / 1e6)]
agg[, log2_fc := log2((density_inside + 1e-6) / (density_outside + 1e-6))]
fwrite(agg, file.path(EXTRAS_TBL_DIR, "PLOT_11_sv_in_inversions.tsv"), sep = "\t")

agg[, svtype := factor(svtype,
                         levels = agg[order(-log2_fc), svtype])]
SV_PAL <- c(DEL = "#d62728", DUP = "#ff7f0e", INV = "#9467bd",
             BND = "#8c564b", INS = "#e377c2")

p <- ggplot(agg, aes(x = log2_fc, y = svtype, fill = svtype)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_text(aes(label = sprintf("%.1f×",
                                  ifelse(log2_fc >= 0,
                                          2^log2_fc, 1 / 2^log2_fc))),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = SV_PAL, drop = FALSE, guide = "none") +
  labs(title = "Large SV enrichment inside inversions",
       subtitle = paste0("Inversion span: ",
                          format(round(INV_LEN/1e6), big.mark = ","),
                          " Mb of ", format(round(GENOME_LEN/1e6), big.mark = ","),
                          " Mb genome"),
       x = "log2(density inside / density outside)",
       y = NULL) +
  xlim(min(0, min(agg$log2_fc) - 0.5),
       max(0, max(agg$log2_fc) + 1)) +
  theme_minimal()

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_11_sv_in_inversions.pdf"),
       p, width = 8, height = 4 + 0.3 * nrow(agg), device = cairo_pdf)
message("[PLOT_11] Done")
