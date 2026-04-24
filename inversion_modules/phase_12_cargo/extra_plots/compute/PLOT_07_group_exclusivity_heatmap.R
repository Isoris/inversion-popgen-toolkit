#!/usr/bin/env Rscript

# =============================================================================
# PLOT_07_group_exclusivity_heatmap.R
#
# Heatmap of carrier frequency for the top-N group-specific variants, one row
# per variant, one column per group. Image 1D / image 2D.
#
# A variant is "group-specific" if its carrier frequency in the primary group
# is at least EXTRAS_EXCLUSIVITY_THRESHOLD (default 0.9) times its total
# carrier count.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

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
EXC_THRESH      <- as.numeric(Sys.getenv("EXTRAS_EXCLUSIVITY_THRESHOLD", "0.9"))
TOPN            <- 30L

args <- commandArgs(trailingOnly = TRUE)
groups_from <- if (any(args == "--groups-from")) args[which(args == "--groups-from") + 1] else NULL

if (!file.exists(VARIANT_MASTER)) { message("[PLOT_07] [skip]"); quit(status = 0) }
sample_to_group <- resolve_groups(groups_from, SAMPLE_LIST,
                                   NGSADMIX_Q_FILE, CANONICAL_K)
if (is.null(sample_to_group)) { message("[PLOT_07] [skip]"); quit(status = 0) }

long <- build_carriership_long(VARIANT_MASTER, sample_to_group)
if (is.null(long) || nrow(long) == 0) {
  message("[PLOT_07] [skip] no carriership data"); quit(status = 0)
}

# Compute exclusivity = max group fraction / total carriers across groups
tot <- long[, .(tot_carriers = sum(n_carriers_in_group),
                top_group_carriers = max(n_carriers_in_group),
                top_group = group[which.max(n_carriers_in_group)]),
            by = var_id]
tot[, exclusivity := top_group_carriers / tot_carriers]

# Pick top exclusive variants per group
exc <- tot[exclusivity >= EXC_THRESH & tot_carriers >= 3]
top_per_group <- exc[order(-tot_carriers), head(.SD, TOPN), by = top_group]
keep_vars <- top_per_group$var_id
if (length(keep_vars) == 0) {
  message("[PLOT_07] no group-exclusive variants pass threshold"); quit(status = 0)
}

mat <- dcast(long[var_id %in% keep_vars],
             var_id ~ group, value.var = "freq_in_group", fill = 0)
fwrite(mat, file.path(EXTRAS_TBL_DIR, "PLOT_07_group_exclusivity.tsv"), sep = "\t")

# Order rows by primary group then total
ord <- top_per_group[order(top_group, -tot_carriers), var_id]
mat[, var_id := factor(var_id, levels = ord)]
long_mat <- melt(mat, id.vars = "var_id",
                  variable.name = "group", value.name = "freq")

p <- ggplot(long_mat, aes(x = group, y = var_id, fill = freq)) +
  geom_tile() +
  scale_fill_viridis_c(name = "carrier freq",
                        limits = c(0, 1), option = "magma", direction = -1) +
  labs(title = paste0("Group-specific variants  (exclusivity ≥ ",
                        EXC_THRESH, ")"),
       subtitle = paste0("Top ", TOPN, " per group, ",
                          length(keep_vars), " total"),
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 4, family = "mono"),
        panel.grid = element_blank())

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_07_group_exclusivity_heatmap.pdf"),
       p, width = max(7, 0.4 * length(unique(long_mat$group))),
       height = max(8, 0.13 * length(keep_vars)),
       device = cairo_pdf)
message("[PLOT_07] Done")
