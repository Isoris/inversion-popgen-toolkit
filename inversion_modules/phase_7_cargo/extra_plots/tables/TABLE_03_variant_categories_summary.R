#!/usr/bin/env Rscript

# =============================================================================
# TABLE_03_variant_categories_summary.R
#
# Image 2's "Variant Categories Summary" — 4-row Private / Rare-Shared /
# Common-Shared / High-frequency table. Computed from cohort SFS bins.
#
# Output:
#   ${EXTRAS_TBL_DIR}/TABLE_03_variant_categories_summary.tsv
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

VARIANT_MASTER <- Sys.getenv("VARIANT_MASTER")
SAMPLE_LIST    <- Sys.getenv("SAMPLE_LIST")
EXTRAS_TBL_DIR <- Sys.getenv("EXTRAS_TBL_DIR")

if (!file.exists(VARIANT_MASTER)) {
  message("[TABLE_03] [skip]"); quit(status = 0)
}
v <- fread(VARIANT_MASTER)
nc_col <- intersect(c("n_carriers", "carrier_count", "AC_samples"), names(v))
if (length(nc_col) > 0) {
  v[, n_carriers := get(nc_col[1])]
} else if ("samples_with_alt" %in% names(v)) {
  v[, n_carriers := lengths(strsplit(samples_with_alt, ";"))]
} else {
  message("[TABLE_03] [skip] no carrier counts"); quit(status = 0)
}

n_total_samples <- if (file.exists(SAMPLE_LIST)) {
  length(readLines(SAMPLE_LIST))
} else {
  max(v$n_carriers) * 2  # rough fallback
}

# Categories — image 2 thresholds adapted to the actual cohort size
high_threshold <- floor(0.5 * n_total_samples)

categorize <- function(n) {
  fcase(
    n == 1,                                  "Private",
    n >= 2 & n <= 5,                          "Rare Shared",
    n >= 6 & n <= 50,                         "Common Shared",
    n > high_threshold,                       "High-frequency",
    default                                   = "Mid-frequency")
}
v[, category := categorize(n_carriers)]

agg <- v[, .(count = .N), by = category]
total <- nrow(v)
agg[, pct_of_total := round(100 * count / total, 1)]

# Add definition + key insight columns
defs <- data.table(
  category = c("Private", "Rare Shared", "Common Shared",
                "Mid-frequency", "High-frequency"),
  definition = c(
    "Observed in only 1 sample",
    "Observed in 2-5 samples",
    "Observed in 6-50 samples",
    sprintf("Observed in 51 to %d samples", high_threshold),
    sprintf("Observed in >%d samples (>50%% of cohort)", high_threshold)),
  key_insight = c(
    "Strong founder effects and recent private mutations",
    "Low-frequency shared variation within groups",
    "Shared across multiple individuals/groups",
    "Polymorphisms segregating across the cohort",
    "Ancestral or balancing variants"))

out <- merge(defs, agg, by = "category", all.x = TRUE)
out <- out[, .(category, definition, count = ifelse(is.na(count), 0L, count),
                pct_of_total = ifelse(is.na(pct_of_total), 0, pct_of_total),
                key_insight)]
# Order: Private, Rare, Common, Mid, High
out <- out[match(c("Private","Rare Shared","Common Shared","Mid-frequency","High-frequency"),
                  category)]

fwrite(out, file.path(EXTRAS_TBL_DIR, "TABLE_03_variant_categories_summary.tsv"),
       sep = "\t")
message("[TABLE_03] Done")
print(out[, .(category, count, pct_of_total)])
