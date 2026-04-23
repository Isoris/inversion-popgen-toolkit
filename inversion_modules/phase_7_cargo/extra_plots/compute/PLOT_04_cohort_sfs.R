#!/usr/bin/env Rscript

# =============================================================================
# PLOT_04_cohort_sfs.R
#
# Cohort-wide variant frequency spectrum: histogram of carrier counts
# (number of samples carrying the variant), stratified by variant class.
# Image 1A / image 2A.
#
# Output: ${EXTRAS_FIG_DIR}/PLOT_04_cohort_sfs.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

VARIANT_MASTER <- Sys.getenv("VARIANT_MASTER")
EXTRAS_FIG_DIR <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR <- Sys.getenv("EXTRAS_TBL_DIR")

if (!file.exists(VARIANT_MASTER)) {
  message("[PLOT_04] [skip] VARIANT_MASTER missing"); quit(status = 0)
}
v <- fread(VARIANT_MASTER)

# Find a carrier-count column or derive one from samples_with_alt
nc_col <- intersect(c("n_carriers", "carrier_count", "AC_samples"), names(v))
if (length(nc_col) > 0) {
  v[, n_carriers := get(nc_col[1])]
} else if ("samples_with_alt" %in% names(v)) {
  v[, n_carriers := lengths(strsplit(samples_with_alt, ";"))]
} else {
  message("[PLOT_04] [skip] no carrier-count column"); quit(status = 0)
}

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

# Bin into 1, 2-3, 4-5, 6-10, 11-20, 21-50, >50
binner <- function(n) {
  fcase(n == 1, "1",
        n <= 3, "2-3",
        n <= 5, "4-5",
        n <= 10, "6-10",
        n <= 20, "11-20",
        n <= 50, "21-50",
        default = ">50")
}
v[, carrier_bin := factor(binner(n_carriers),
                            levels = c("1","2-3","4-5","6-10","11-20","21-50",">50"))]

agg <- v[, .N, by = .(class, carrier_bin)]
fwrite(agg, file.path(EXTRAS_TBL_DIR, "PLOT_04_cohort_sfs_table.tsv"), sep = "\t")

CLASS_PAL <- c(SNP = "#1f77b4", INDEL = "#2ca02c",
                DEL = "#d62728", DUP = "#ff7f0e",
                INV = "#9467bd", INS = "#e377c2",
                BND = "#8c564b", PAV = "#7f7f7f", OTHER = "grey80")
agg[, class := factor(class, levels = names(CLASS_PAL))]

p <- ggplot(agg, aes(x = carrier_bin, y = N, fill = class)) +
  geom_col(position = "stack") +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  scale_fill_manual(values = CLASS_PAL, drop = FALSE) +
  labs(title = "Cohort variant frequency spectrum",
       subtitle = paste0("Total variants: ", format(nrow(v), big.mark = ",")),
       x = "Number of carrier samples",
       y = "Number of variants (log10)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_04_cohort_sfs.pdf"),
       p, width = 9, height = 6, device = cairo_pdf)
message("[PLOT_04] Done")
