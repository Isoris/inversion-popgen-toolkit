#!/usr/bin/env Rscript

# =============================================================================
# TABLE_01_top_recurrent_with_group_AF.R
#
# Image 1B's "Top 12 Recurrent Variants" — the version that adds a Max-AF-by-
# group column and a single best-annotation column. Output is suitable for
# direct import into the manuscript.
#
# Columns:
#   variant_id, chrom, pos, type, n_carriers, carrier_freq,
#   max_af_group, max_af_value, annotation, gene
#
# Output:
#   ${EXTRAS_TBL_DIR}/TABLE_01_top_recurrent_with_group_AF.tsv
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

# Resolve script-relative path
.this <- {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca)
  if (length(m)) sub("^--file=", "", ca[m]) else "."
}
source(file.path(dirname(normalizePath(.this, mustWork = FALSE)),
                  "..", "compute", "_lib_group_carriership.R"))

VARIANT_MASTER  <- Sys.getenv("VARIANT_MASTER")
SAMPLE_LIST     <- Sys.getenv("SAMPLE_LIST")
NGSADMIX_Q_FILE <- Sys.getenv("NGSADMIX_Q_FILE")
EXTRAS_TBL_DIR  <- Sys.getenv("EXTRAS_TBL_DIR")
CANONICAL_K     <- as.integer(Sys.getenv("CANONICAL_K", "8"))
TOPN <- as.integer(Sys.getenv("EXTRAS_TOPN_RECURRENT", "25"))

if (!file.exists(VARIANT_MASTER)) {
  message("[TABLE_01] [skip] VARIANT_MASTER missing"); quit(status = 0)
}

# Tiny helper for nullable fallback
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
groups_from <- if (any(args == "--groups-from")) args[which(args == "--groups-from") + 1] else NULL

sample_to_group <- resolve_groups(groups_from, SAMPLE_LIST,
                                   NGSADMIX_Q_FILE, CANONICAL_K)

v <- fread(VARIANT_MASTER)
nc_col <- intersect(c("n_carriers", "carrier_count", "AC_samples"), names(v))
if (length(nc_col) > 0) {
  v[, n_carriers := get(nc_col[1])]
} else if ("samples_with_alt" %in% names(v)) {
  v[, n_carriers := lengths(strsplit(samples_with_alt, ";"))]
} else {
  message("[TABLE_01] [skip] no carrier counts"); quit(status = 0)
}

v[, class := classify_variants(.SD)]
v[, var_id := if ("var_key" %in% names(.SD)) var_key else paste0(
  if ("chr" %in% names(.SD)) chr else chrom, ":", pos, ":",
  if ("ref" %in% names(.SD)) ref else "?", ":",
  if ("alt" %in% names(.SD)) alt else "?")]
chr_col <- if ("chr" %in% names(v)) "chr" else "chrom"

# Top N by carrier count
top <- v[order(-n_carriers)][seq_len(min(TOPN, .N))]

# Add per-group max AF
if (!is.null(sample_to_group)) {
  message("[TABLE_01] computing per-group AF for top variants...")
  long_top <- build_carriership_long(VARIANT_MASTER, sample_to_group)
  if (!is.null(long_top) && nrow(long_top) > 0) {
    long_top <- long_top[var_id %in% top$var_id]
    # AF = (carriers / n_group) — under HWE-ish, AF ≈ carrier_freq for SNPs
    # We report carrier-frequency-by-group as proxy
    grp_max <- long_top[, {
      i <- which.max(freq_in_group)
      .(max_af_group = group[i], max_af_value = freq_in_group[i])
    }, by = var_id]
    top <- merge(top, grp_max, by = "var_id", all.x = TRUE)
  }
} else {
  top[, `:=`(max_af_group = "", max_af_value = NA_real_)]
}

# Total cohort sample count for carrier_freq
n_total_samples <- if (!is.null(sample_to_group)) length(sample_to_group) else
                   length(unique(unlist(strsplit(top$samples_with_alt %||% "", ";"))))
if (n_total_samples == 0) n_total_samples <- nrow(v) / 1000  # fallback heuristic
top[, carrier_freq := round(n_carriers / max(n_total_samples, 1), 4)]

# Build output
out <- top[, .(
  variant_id  = var_id,
  chrom       = get(chr_col),
  pos         = if ("pos" %in% names(top)) pos else NA,
  type        = class,
  n_carriers,
  carrier_freq,
  max_af_group = if ("max_af_group" %in% names(top)) max_af_group else "",
  max_af_value = if ("max_af_value" %in% names(top)) round(max_af_value, 4) else NA_real_,
  annotation  = if ("snpeff_annotation" %in% names(top)) snpeff_annotation else "",
  gene        = if ("gene_id" %in% names(top)) gene_id else
                if ("bcsq_gene" %in% names(top)) bcsq_gene else ""
)]

fwrite(out, file.path(EXTRAS_TBL_DIR, "TABLE_01_top_recurrent_with_group_AF.tsv"),
       sep = "\t")
message("[TABLE_01] Done — ", nrow(out), " rows → TABLE_01_top_recurrent_with_group_AF.tsv")
