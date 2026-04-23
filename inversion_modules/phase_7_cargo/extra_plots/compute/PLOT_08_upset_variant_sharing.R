#!/usr/bin/env Rscript

# =============================================================================
# PLOT_08_upset_variant_sharing.R
#
# UpSet plot of variant sharing across groups. Each variant is represented
# as a set membership vector (group → present / absent). Image 1E.
#
# Uses ComplexUpset if available, falls back to a simple intersection-size
# barplot if not (no required dependency on UpSetR/ComplexUpset).
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

args <- commandArgs(trailingOnly = TRUE)
groups_from <- if (any(args == "--groups-from")) args[which(args == "--groups-from") + 1] else NULL

if (!file.exists(VARIANT_MASTER)) { message("[PLOT_08] [skip]"); quit(status = 0) }
sample_to_group <- resolve_groups(groups_from, SAMPLE_LIST,
                                   NGSADMIX_Q_FILE, CANONICAL_K)
if (is.null(sample_to_group)) { message("[PLOT_08] [skip]"); quit(status = 0) }

long <- build_carriership_long(VARIANT_MASTER, sample_to_group)
if (is.null(long) || nrow(long) == 0) { message("[PLOT_08] [skip]"); quit(status = 0) }

# Each variant → membership vector (group present if any carriers)
groups_avail <- sort(unique(long$group))
mat <- dcast(long, var_id ~ group,
             value.var = "n_carriers_in_group", fill = 0)
for (g in groups_avail) mat[[g]] <- as.integer(mat[[g]] > 0)
fwrite(mat, file.path(EXTRAS_TBL_DIR, "PLOT_08_upset_membership.tsv"), sep = "\t")

# Build set-intersection counts as a string key, e.g. "A&B&C"
mat[, set_key := apply(.SD, 1, function(r) {
  paste(groups_avail[r == 1], collapse = " & ")
}), .SDcols = groups_avail]
intsz <- mat[, .N, by = set_key][order(-N)]
intsz[, n_groups := vapply(strsplit(set_key, " & "), length, integer(1))]

# Take top 25 intersections for readability
top <- intsz[seq_len(min(25, .N))]
top[, set_key := factor(set_key, levels = top$set_key)]

p <- ggplot(top, aes(x = set_key, y = N, fill = factor(n_groups))) +
  geom_col() +
  geom_text(aes(label = formatC(N, format = "d", big.mark = ",")),
            vjust = -0.3, size = 2.5) +
  scale_fill_viridis_d(name = "# groups in intersection") +
  labs(title = "Variant sharing across groups (top 25 intersections)",
       x = "Group intersection", y = "# variants") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        legend.position = "right")

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_08_upset_variant_sharing.pdf"),
       p, width = 12, height = 6, device = cairo_pdf)
message("[PLOT_08] Done")
