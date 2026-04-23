#!/usr/bin/env Rscript

# =============================================================================
# PLOT_06_private_vs_shared_by_group.R
#
# Per-group stacked bar of variants segregating only within that group
# (private) vs shared with at least one other group (shared). Image 1C.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

# Resolve script directory (works under Rscript and source())
.this_script <- {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca)
  if (length(m)) sub("^--file=", "", ca[m]) else
    if (!is.null(sys.frames()) && length(sys.frames()) > 0 &&
        !is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else "."
}
source(file.path(dirname(normalizePath(.this_script, mustWork = FALSE)),
                  "_lib_group_carriership.R"))

VARIANT_MASTER  <- Sys.getenv("VARIANT_MASTER")
SAMPLE_LIST     <- Sys.getenv("SAMPLE_LIST")
NGSADMIX_Q_FILE <- Sys.getenv("NGSADMIX_Q_FILE")
EXTRAS_FIG_DIR  <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR  <- Sys.getenv("EXTRAS_TBL_DIR")
CANONICAL_K     <- as.integer(Sys.getenv("CANONICAL_K", "8"))

args <- commandArgs(trailingOnly = TRUE)
groups_from <- if (any(args == "--groups-from")) args[which(args == "--groups-from") + 1] else NULL

if (!file.exists(VARIANT_MASTER)) { message("[PLOT_06] [skip]"); quit(status = 0) }
sample_to_group <- resolve_groups(groups_from, SAMPLE_LIST,
                                   NGSADMIX_Q_FILE, CANONICAL_K)
if (is.null(sample_to_group)) {
  message("[PLOT_06] [skip] no groups available"); quit(status = 0)
}

long <- build_carriership_long(VARIANT_MASTER, sample_to_group)
if (is.null(long) || nrow(long) == 0) {
  message("[PLOT_06] [skip] no carriership data"); quit(status = 0)
}

# Variant is "private to group X" if X is the only group in which it appears
n_groups_per_var <- long[, .(n_groups = uniqueN(group)), by = var_id]
private_vars <- n_groups_per_var[n_groups == 1, var_id]
long[, status := ifelse(var_id %in% private_vars, "private", "shared")]

agg <- long[, .(n = uniqueN(var_id)), by = .(group, status)]
fwrite(agg, file.path(EXTRAS_TBL_DIR, "PLOT_06_private_vs_shared.tsv"), sep = "\t")

p <- ggplot(agg, aes(x = group, y = n, fill = status)) +
  geom_col() +
  scale_fill_manual(values = c(private = "#1f3a5f", shared = "#a3c4f3")) +
  labs(title = "Private vs shared variants by group",
       x = NULL, y = "Variant count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_06_private_vs_shared_by_group.pdf"),
       p, width = max(7, 0.5 * length(unique(agg$group))), height = 5,
       device = cairo_pdf)
message("[PLOT_06] Done")
