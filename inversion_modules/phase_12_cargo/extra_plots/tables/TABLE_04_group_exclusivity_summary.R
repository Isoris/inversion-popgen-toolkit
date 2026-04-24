#!/usr/bin/env Rscript

# =============================================================================
# TABLE_04_group_exclusivity_summary.R
#
# Image 2's "Group Exclusivity Statistics" — one row per group:
#   group, n_samples, n_private_variants, n_exclusive_variants,
#   exclusivity_rate, top_variant_types
#
# A "private" variant is one observed in samples from only that group.
# An "exclusive" variant is one whose carrier-count is at least
# EXTRAS_EXCLUSIVITY_THRESHOLD (default 0.9) within that group.
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
EXC_THRESH      <- as.numeric(Sys.getenv("EXTRAS_EXCLUSIVITY_THRESHOLD", "0.9"))

if (!file.exists(VARIANT_MASTER)) {
  message("[TABLE_04] [skip]"); quit(status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
groups_from <- if (any(args == "--groups-from")) args[which(args == "--groups-from") + 1] else NULL
sample_to_group <- resolve_groups(groups_from, SAMPLE_LIST,
                                   NGSADMIX_Q_FILE, CANONICAL_K)
if (is.null(sample_to_group)) {
  message("[TABLE_04] [skip] no groups"); quit(status = 0)
}

long <- build_carriership_long(VARIANT_MASTER, sample_to_group)
if (is.null(long) || nrow(long) == 0) {
  message("[TABLE_04] [skip] no carriership"); quit(status = 0)
}

# Per-variant: which groups carry it, and what's the max group fraction
per_var <- long[, .(n_groups = uniqueN(group),
                     primary_group = group[which.max(n_carriers_in_group)],
                     primary_freq = max(freq_in_group),
                     tot_carriers = sum(n_carriers_in_group),
                     primary_carriers = max(n_carriers_in_group)),
                 by = .(var_id, class)]
per_var[, exclusivity := primary_carriers / tot_carriers]
per_var[, is_private := n_groups == 1]
per_var[, is_exclusive := exclusivity >= EXC_THRESH]

# Group sizes
group_sizes <- table(sample_to_group)

# Per-group summary
rows <- list()
for (g in sort(unique(unname(sample_to_group)))) {
  pg <- per_var[primary_group == g]
  prv <- pg[is_private == TRUE]
  exc <- pg[is_exclusive == TRUE]
  type_counts <- if (nrow(prv) > 0) {
    sort(table(prv$class), decreasing = TRUE)
  } else {
    integer(0)
  }
  top_types <- if (length(type_counts) == 0) "-" else
               paste(paste0(names(type_counts)[seq_len(min(3, length(type_counts)))],
                             "(", type_counts[seq_len(min(3, length(type_counts)))], ")"),
                      collapse = ", ")
  # Enrichment: -log10(p) for hypergeometric "more private than expected"
  N_total <- nrow(per_var)
  K_private <- per_var[is_private == TRUE, .N]
  n_in_g <- pg[, .N]
  k_obs <- nrow(prv)
  pval <- tryCatch(phyper(k_obs - 1, K_private, N_total - K_private, n_in_g,
                            lower.tail = FALSE),
                   error = function(e) NA_real_)
  rows[[length(rows) + 1]] <- data.table(
    group = g,
    n_samples = as.integer(group_sizes[g]),
    n_variants_with_primary_in_group = nrow(pg),
    n_private_variants = nrow(prv),
    n_exclusive_variants = nrow(exc),
    exclusivity_rate = if (nrow(pg) > 0) round(nrow(exc) / nrow(pg), 4) else NA_real_,
    enrichment_neg_log10p = if (!is.na(pval) && pval > 0) round(-log10(pval), 2) else NA_real_,
    top_variant_types = top_types)
}
out <- rbindlist(rows)
out <- out[order(-n_exclusive_variants)]

fwrite(out, file.path(EXTRAS_TBL_DIR, "TABLE_04_group_exclusivity_summary.tsv"),
       sep = "\t")
message("[TABLE_04] Done")
print(out)
