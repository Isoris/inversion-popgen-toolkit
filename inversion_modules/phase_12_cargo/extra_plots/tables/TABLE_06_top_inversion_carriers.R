#!/usr/bin/env Rscript

# =============================================================================
# TABLE_06_top_inversion_carriers.R
#
# Image 7B — clean ranked table of inversions by carrier frequency.
#   inversion_id, chrom, start_bp, end_bp, length_mb,
#   n_HOM_REF, n_HET, n_HOM_INV, carrier_freq, allele_freq, annotation
#
# annotation: takes the most-affected gene from cargo inventory if available;
# otherwise marks "intergenic" if no genes overlap, else lists the first 1-3
# overlapping genes.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

SNAKE_CAND_FILE     <- Sys.getenv("SNAKE_CAND_FILE")
SAMPLE_REGISTRY     <- Sys.getenv("SAMPLE_REGISTRY")
CARGO_INVENTORY_DIR <- Sys.getenv("CARGO_INVENTORY_DIR")
EXTRAS_TBL_DIR      <- Sys.getenv("EXTRAS_TBL_DIR")
TOPN <- as.integer(Sys.getenv("EXTRAS_TOPN_INV", "30"))

if (!file.exists(SNAKE_CAND_FILE)) {
  message("[TABLE_06] [skip]"); quit(status = 0)
}
cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))

read_grp_size <- function(gid) {
  fp <- file.path(SAMPLE_REGISTRY, "groups", paste0(gid, ".txt"))
  if (!file.exists(fp)) return(0L)
  length(trimws(readLines(fp)))
}

inventory_lookup <- function(cid) {
  fp <- file.path(CARGO_INVENTORY_DIR, cid, "genes.tsv")
  if (!file.exists(fp)) return(NA_character_)
  g <- fread(fp)
  if (nrow(g) == 0) return("intergenic")
  # Prefer named genes
  named <- g[nzchar(preferred_name)]
  if (nrow(named) >= 1) {
    return(paste(head(named$preferred_name, 3), collapse = ","))
  }
  paste(head(g$gene_id, 3), collapse = ",")
}

rows <- list()
for (i in seq_len(nrow(cands))) {
  cid <- cands$candidate_id[i]
  n_ref <- read_grp_size(paste0("inv_", cid, "_HOM_REF"))
  n_het <- read_grp_size(paste0("inv_", cid, "_HET"))
  n_inv <- read_grp_size(paste0("inv_", cid, "_HOM_INV"))
  n_total <- n_ref + n_het + n_inv
  if (n_total == 0) next
  rows[[length(rows) + 1]] <- data.table(
    inversion_id = cid,
    chrom        = cands$chrom[i],
    start_bp     = cands$start_bp[i],
    end_bp       = cands$end_bp[i],
    length_mb    = round((cands$end_bp[i] - cands$start_bp[i]) / 1e6, 3),
    n_HOM_REF    = n_ref,
    n_HET        = n_het,
    n_HOM_INV    = n_inv,
    carrier_freq = round((n_het + n_inv) / n_total, 4),
    allele_freq  = round((n_het + 2 * n_inv) / (2 * n_total), 4),
    annotation   = inventory_lookup(cid))
}
if (length(rows) == 0) {
  message("[TABLE_06] no karyotype groups registered"); quit(status = 0)
}
dt <- rbindlist(rows)[order(-allele_freq)]
top <- dt[seq_len(min(TOPN, .N))]

# Also write the full ranked table — useful for supplementary
fwrite(dt,  file.path(EXTRAS_TBL_DIR, "TABLE_06_all_inversions_ranked.tsv"),  sep = "\t")
fwrite(top, file.path(EXTRAS_TBL_DIR, "TABLE_06_top_inversion_carriers.tsv"), sep = "\t")

message("[TABLE_06] Done — ", nrow(dt), " inversions ranked, top ", nrow(top),
        " in TABLE_06_top_inversion_carriers.tsv")
