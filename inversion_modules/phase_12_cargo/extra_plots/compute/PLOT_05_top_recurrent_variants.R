#!/usr/bin/env Rscript

# =============================================================================
# PLOT_05_top_recurrent_variants.R
#
# Top-N variants by carrier count. Emits:
#   - ${EXTRAS_TBL_DIR}/PLOT_05_top_recurrent_variants.tsv
#   - ${EXTRAS_FIG_DIR}/PLOT_05_top_recurrent_variants.pdf  (horizontal barchart)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

VARIANT_MASTER <- Sys.getenv("VARIANT_MASTER")
EXTRAS_FIG_DIR <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR <- Sys.getenv("EXTRAS_TBL_DIR")
TOPN <- as.integer(Sys.getenv("EXTRAS_TOPN_RECURRENT", "25"))

if (!file.exists(VARIANT_MASTER)) {
  message("[PLOT_05] [skip]"); quit(status = 0)
}
v <- fread(VARIANT_MASTER)
nc_col <- intersect(c("n_carriers", "carrier_count", "AC_samples"), names(v))
if (length(nc_col) > 0) {
  v[, n_carriers := get(nc_col[1])]
} else if ("samples_with_alt" %in% names(v)) {
  v[, n_carriers := lengths(strsplit(samples_with_alt, ";"))]
} else {
  message("[PLOT_05] [skip] no carrier counts"); quit(status = 0)
}

# Build var_id and class label
v[, var_id := if ("var_key" %in% names(v)) var_key else paste0(
  if ("chr" %in% names(.SD)) chr else chrom, ":", pos, ":",
  if ("ref" %in% names(.SD)) ref else "?", ":",
  if ("alt" %in% names(.SD)) alt else "?")]

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

top <- v[order(-n_carriers)][seq_len(min(TOPN, .N))]
fwrite(top[, .(var_id, class, n_carriers,
               annotation = if ("snpeff_annotation" %in% names(top))
                              snpeff_annotation else "")],
       file.path(EXTRAS_TBL_DIR, "PLOT_05_top_recurrent_variants.tsv"), sep = "\t")

CLASS_PAL <- c(SNP = "#1f77b4", INDEL = "#2ca02c",
                DEL = "#d62728", DUP = "#ff7f0e",
                INV = "#9467bd", INS = "#e377c2",
                BND = "#8c564b", PAV = "#7f7f7f", OTHER = "grey80")
top[, var_id := factor(var_id, levels = rev(top$var_id))]

p <- ggplot(top, aes(x = n_carriers, y = var_id, fill = class)) +
  geom_col() +
  scale_fill_manual(values = CLASS_PAL) +
  labs(title = paste0("Top ", TOPN, " most-recurrent variants"),
       x = "Carrier count",
       y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7, family = "mono"))

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_05_top_recurrent_variants.pdf"),
       p, width = 10, height = max(4, 0.25 * TOPN + 1), device = cairo_pdf)
message("[PLOT_05] Done")
