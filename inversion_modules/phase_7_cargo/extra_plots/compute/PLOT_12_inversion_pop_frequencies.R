#!/usr/bin/env Rscript

# =============================================================================
# PLOT_12_inversion_pop_frequencies.R
#
# Distribution of inversion carrier frequencies across the cohort. Image 7A.
#   Panel A: histogram of HOM_INV+HET frequency across all candidates
#   Panel B: top-N inversions by carrier frequency (table-style barchart)
#
# Reads sample_registry groups inv_<cid>_HOM_INV / HET / HOM_REF.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

SNAKE_CAND_FILE <- Sys.getenv("SNAKE_CAND_FILE")
SAMPLE_REGISTRY <- Sys.getenv("SAMPLE_REGISTRY")
EXTRAS_FIG_DIR  <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR  <- Sys.getenv("EXTRAS_TBL_DIR")
TOPN <- 25L

if (!file.exists(SNAKE_CAND_FILE)) { message("[PLOT_12] [skip]"); quit(status = 0) }
cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))

read_grp <- function(gid) {
  fp <- file.path(SAMPLE_REGISTRY, "groups", paste0(gid, ".txt"))
  if (!file.exists(fp)) return(0L)
  length(trimws(readLines(fp)))
}

rows <- list()
for (cid in cands$candidate_id) {
  n_ref <- read_grp(paste0("inv_", cid, "_HOM_REF"))
  n_het <- read_grp(paste0("inv_", cid, "_HET"))
  n_inv <- read_grp(paste0("inv_", cid, "_HOM_INV"))
  n_total <- n_ref + n_het + n_inv
  if (n_total == 0) next
  rows[[length(rows) + 1]] <- data.table(
    candidate_id = cid,
    chrom = cands$chrom[cands$candidate_id == cid],
    n_HOM_REF = n_ref, n_HET = n_het, n_HOM_INV = n_inv,
    n_total = n_total,
    carrier_freq = (n_het + n_inv) / n_total,    # at least one INV haplotype
    inv_allele_freq = (n_het + 2 * n_inv) / (2 * n_total))
}
if (length(rows) == 0) {
  message("[PLOT_12] no karyotype groups registered yet"); quit(status = 0)
}
dt <- rbindlist(rows)
fwrite(dt, file.path(EXTRAS_TBL_DIR, "PLOT_12_inversion_pop_freq_table.tsv"),
       sep = "\t")

# Panel A: distribution
pA <- ggplot(dt, aes(x = inv_allele_freq)) +
  geom_histogram(bins = 30, fill = "#9467bd", color = "white", alpha = 0.85) +
  geom_vline(xintercept = c(0.05, 0.5), color = "grey40",
              linetype = "dashed") +
  annotate("text", x = 0.05, y = Inf, label = "rare (5%)",
            vjust = 1.5, hjust = -0.05, size = 3, color = "grey30") +
  annotate("text", x = 0.5, y = Inf, label = "balanced (50%)",
            vjust = 1.5, hjust = -0.05, size = 3, color = "grey30") +
  labs(title = "A. Inversion allele frequency distribution",
       subtitle = paste0(nrow(dt), " candidates"),
       x = "Inversion allele frequency",
       y = "# candidates") +
  theme_minimal()

# Panel B: top N
top <- dt[order(-inv_allele_freq)][seq_len(min(TOPN, .N))]
top[, label := paste0(candidate_id, "  (", chrom, ")")]
top[, label := factor(label, levels = rev(top$label))]
top_long <- melt(top, id.vars = "label",
                  measure.vars = c("n_HOM_REF", "n_HET", "n_HOM_INV"),
                  variable.name = "kary", value.name = "n")
top_long[, kary := factor(kary, levels = c("n_HOM_REF","n_HET","n_HOM_INV"))]

pB <- ggplot(top_long, aes(x = n, y = label, fill = kary)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(n_HOM_REF = "#1f77b4",
                                  n_HET = "#7f7f7f",
                                  n_HOM_INV = "#d62728"),
                       labels = c("HOM_REF","HET","HOM_INV"),
                       name = "karyotype") +
  labs(title = paste0("B. Top ", TOPN, " inversions by allele frequency"),
       x = "Sample count",
       y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7, family = "mono"),
        legend.position = "bottom")

fig <- pA | pB
fig <- fig + plot_annotation(title = "Inversion population frequencies")
ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_12_inversion_pop_frequencies.pdf"),
       fig, width = 14, height = max(6, 0.3 * TOPN + 2),
       device = cairo_pdf)
message("[PLOT_12] Done")
