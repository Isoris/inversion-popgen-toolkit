#!/usr/bin/env Rscript

# =============================================================================
# BREEDING_C_founder_haplotype_tracker.R
#
# For each inversion locus, decompose the current cohort's HOM_INV haplotypes
# by their inferred founder lineage (proxy: NGSadmix K-cluster argmax).
# Identify inversions where the HOM_INV state is overrepresented in a single
# founder lineage — those are bottleneck-risk loci, where the entire inverted
# allele in the hatchery traces back to one founder pair.
#
# Caveat: without explicit pedigree founders, "founder lineage" is a proxy.
# We use the canonical-K NGSadmix Q matrix and assign each sample to its
# argmax cluster. If you have actual pedigree data, pass --groups-from with
# columns sample / group_id where group_id is the founder lineage label.
#
# Outputs:
#   ${EXTRAS_TBL_DIR}/BREEDING_C_founder_inversion_matrix.tsv
#       inversion × founder_lineage matrix of HOM_INV carrier counts
#   ${EXTRAS_TBL_DIR}/BREEDING_C_bottleneck_flags.tsv
#       per-inversion bottleneck score:
#         max_lineage_share  = (max carriers in single lineage) / (total carriers)
#         shannon_diversity  = entropy of carriers across lineages (bits)
#         flag               = "BOTTLENECK" | "concentrated" | "diverse"
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

.this <- {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca)
  if (length(m)) sub("^--file=", "", ca[m]) else "."
}
source(file.path(dirname(normalizePath(.this, mustWork = FALSE)),
                  "..", "compute", "_lib_group_carriership.R"))

SNAKE_CAND_FILE <- Sys.getenv("SNAKE_CAND_FILE")
SAMPLE_REGISTRY <- Sys.getenv("SAMPLE_REGISTRY")
SAMPLE_LIST     <- Sys.getenv("SAMPLE_LIST")
NGSADMIX_Q_FILE <- Sys.getenv("NGSADMIX_Q_FILE")
EXTRAS_TBL_DIR  <- Sys.getenv("EXTRAS_TBL_DIR")
CANONICAL_K     <- as.integer(Sys.getenv("CANONICAL_K", "8"))

args <- commandArgs(trailingOnly = TRUE)
groups_from <- if (any(args == "--groups-from")) args[which(args == "--groups-from") + 1] else NULL
sample_to_group <- resolve_groups(groups_from, SAMPLE_LIST,
                                   NGSADMIX_Q_FILE, CANONICAL_K)
if (is.null(sample_to_group)) {
  message("[BREEDING_C] [skip] no founder-lineage groups available")
  quit(status = 0)
}
if (!file.exists(SNAKE_CAND_FILE)) {
  message("[BREEDING_C] [skip]"); quit(status = 0)
}

cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))
lineages <- sort(unique(unname(sample_to_group)))
n_lineage_total <- table(sample_to_group)

read_grp <- function(gid) {
  fp <- file.path(SAMPLE_REGISTRY, "groups", paste0(gid, ".txt"))
  if (!file.exists(fp)) return(character())
  trimws(readLines(fp))
}

# Per-candidate: count HOM_INV carriers per founder lineage
rows <- list()
mat_rows <- list()
for (i in seq_len(nrow(cands))) {
  cid <- cands$candidate_id[i]
  hom_inv <- read_grp(paste0("inv_", cid, "_HOM_INV"))
  het     <- read_grp(paste0("inv_", cid, "_HET"))
  inv_carriers <- c(hom_inv, het)  # any inv allele
  if (length(inv_carriers) == 0) next
  by_lin <- table(sample_to_group[inv_carriers])
  by_lin <- by_lin[lineages]; by_lin[is.na(by_lin)] <- 0L
  total <- sum(by_lin)

  # Diversity stats
  p <- by_lin / total
  shannon_bits <- -sum(p[p > 0] * log2(p[p > 0]))
  max_lin_share <- max(p)
  top_lineage <- names(by_lin)[which.max(by_lin)]

  # Compare to lineage-size baseline (what we'd expect if INV alleles were
  # randomly distributed in proportion to lineage size)
  expected <- as.numeric(n_lineage_total[lineages]) / sum(as.numeric(n_lineage_total))
  observed <- as.numeric(by_lin / total)
  # Normalized chi-square against null
  chi <- sum((observed - expected)^2 / pmax(expected, 1e-6)) * total
  pval <- tryCatch(pchisq(chi, df = length(lineages) - 1, lower.tail = FALSE),
                   error = function(e) NA_real_)

  flag <- fcase(max_lin_share >= 0.8 & total >= 5, "BOTTLENECK",
                  max_lin_share >= 0.5,             "concentrated",
                  default                           = "diverse")

  rows[[length(rows) + 1]] <- data.table(
    candidate_id = cid,
    chrom = cands$chrom[i],
    n_inv_carriers = total,
    top_lineage = top_lineage,
    max_lineage_share = round(max_lin_share, 3),
    shannon_diversity_bits = round(shannon_bits, 3),
    chisq_vs_uniform = round(chi, 2),
    chisq_p = if (is.na(pval)) NA_real_ else round(pval, 6),
    flag = flag)

  mat_row <- as.list(by_lin)
  names(mat_row) <- paste0("L_", lineages)
  mat_rows[[length(mat_rows) + 1]] <- c(list(candidate_id = cid), mat_row)
}
if (length(rows) == 0) {
  message("[BREEDING_C] no candidates with INV carriers"); quit(status = 0)
}

flags <- rbindlist(rows)
flags <- flags[order(-max_lineage_share, -n_inv_carriers)]
fwrite(flags, file.path(EXTRAS_TBL_DIR, "BREEDING_C_bottleneck_flags.tsv"),
       sep = "\t")

mat <- rbindlist(mat_rows, fill = TRUE)
fwrite(mat, file.path(EXTRAS_TBL_DIR, "BREEDING_C_founder_inversion_matrix.tsv"),
       sep = "\t")

message("[BREEDING_C] Done")
message("  ", flags[flag == "BOTTLENECK", .N], " loci flagged BOTTLENECK (≥80% in single lineage)")
message("  ", flags[flag == "concentrated", .N], " loci concentrated (≥50% in single lineage)")
message("  ", flags[flag == "diverse", .N], " loci diverse (<50%)")
print(head(flags, 10))
