#!/usr/bin/env Rscript

# =============================================================================
# TABLE_02_top_private_singletons.R
#
# Image 2C — top singleton variants (carrier_count == 1) and which sample/group
# they belong to. Useful for QC review of "private founder" candidates.
#
# Output:
#   ${EXTRAS_TBL_DIR}/TABLE_02_top_private_singletons.tsv
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

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
TOPN <- as.integer(Sys.getenv("EXTRAS_TOPN_SINGLETONS", "200"))

if (!file.exists(VARIANT_MASTER)) {
  message("[TABLE_02] [skip]"); quit(status = 0)
}

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
  message("[TABLE_02] [skip] no carrier counts"); quit(status = 0)
}

singletons <- v[n_carriers == 1]
if (nrow(singletons) == 0) {
  message("[TABLE_02] no singletons in dataset"); quit(status = 0)
}
message("[TABLE_02] ", nrow(singletons), " total singletons")

singletons[, class := classify_variants(.SD)]
chr_col <- if ("chr" %in% names(singletons)) "chr" else "chrom"

# Sample assignment — needs samples_with_alt or per-sample dosage
sample_for_singleton <- function(rows) {
  if ("samples_with_alt" %in% names(rows)) {
    return(sub(";.*", "", rows$samples_with_alt))
  }
  if (!is.null(sample_to_group)) {
    sample_cols <- intersect(names(sample_to_group), names(rows))
    if (length(sample_cols) > 0) {
      # For each row, find the column with dosage > 0
      M <- as.matrix(rows[, sample_cols, with = FALSE])
      idx <- max.col(M > 0, ties.method = "first")
      return(sample_cols[idx])
    }
  }
  rep(NA_character_, nrow(rows))
}
singletons[, sample := sample_for_singleton(.SD)]
singletons[, group := if (!is.null(sample_to_group)) sample_to_group[sample] else NA_character_]

# Take top-N (or all if fewer)
out_rows <- singletons[order(get(chr_col), pos)][seq_len(min(TOPN, .N))]

out <- out_rows[, .(
  variant_id = if ("var_key" %in% names(out_rows)) var_key else
               paste0(get(chr_col), ":", pos, ":",
                      if ("ref" %in% names(out_rows)) ref else "?", ":",
                      if ("alt" %in% names(out_rows)) alt else "?"),
  chrom = get(chr_col),
  pos,
  type = class,
  sample,
  group,
  annotation = if ("snpeff_annotation" %in% names(out_rows)) snpeff_annotation else "",
  gene       = if ("gene_id" %in% names(out_rows)) gene_id else
               if ("bcsq_gene" %in% names(out_rows)) bcsq_gene else ""
)]

fwrite(out, file.path(EXTRAS_TBL_DIR, "TABLE_02_top_private_singletons.tsv"),
       sep = "\t")
message("[TABLE_02] Done — ", nrow(out), " singletons → TABLE_02_top_private_singletons.tsv")
