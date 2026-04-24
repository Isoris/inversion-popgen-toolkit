#!/usr/bin/env Rscript

# =============================================================================
# PLOT_02_inversion_sharing_across_groups.R
#
# Heatmap of HOM_INV carrier frequency for every candidate inversion × every
# ancestry / family group. Answers the manuscript question: "do inversions
# segregate within hatchery families/lineages, or across them?"
#
# A row is a candidate. A column is a group. Cell colour = fraction of group
# members that are HOM_INV. Right-side annotation: chromosome, span (Mb).
# Bottom annotation: group size.
#
# Group definition: by default uses NGSadmix Q at canonical K — assigns each
# sample to its argmax cluster (hard-K assignment). If a custom group-set is
# preferred (e.g. families), pass --groups-from <path/to/sample_groups.tsv>
# with columns sample / group_id.
#
# Inputs:
#   ${SNAKE_CAND_FILE}
#   ${SAMPLE_REGISTRY}/groups/inv_<cid>_HOM_INV.txt   (per-candidate carriers)
#   ${SAMPLE_REGISTRY}/groups/inv_<cid>_HOM_REF.txt
#   ${NGSADMIX_Q_FILE}    OR   --groups-from <tsv>
#
# Output:
#   ${EXTRAS_FIG_DIR}/PLOT_02_inversion_sharing_across_groups.pdf
#   ${EXTRAS_TBL_DIR}/PLOT_02_inversion_sharing_table.tsv
#
# Usage:
#   Rscript PLOT_02_inversion_sharing_across_groups.R [--groups-from <tsv>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
groups_from <- NULL
if (any(args == "--groups-from")) {
  i <- which(args == "--groups-from")
  groups_from <- args[i + 1]
}

SNAKE_CAND_FILE  <- Sys.getenv("SNAKE_CAND_FILE")
SAMPLE_REGISTRY  <- Sys.getenv("SAMPLE_REGISTRY")
NGSADMIX_Q_FILE  <- Sys.getenv("NGSADMIX_Q_FILE")
EXTRAS_FIG_DIR   <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR   <- Sys.getenv("EXTRAS_TBL_DIR")
SAMPLE_LIST      <- Sys.getenv("SAMPLE_LIST")
CANONICAL_K      <- as.integer(Sys.getenv("CANONICAL_K", "8"))

# ── Resolve groups ───────────────────────────────────────────────────────────
sample_to_group <- NULL
group_label <- "ancestry_K"

if (!is.null(groups_from) && file.exists(groups_from)) {
  message("[PLOT_02] Loading groups from ", groups_from)
  gd <- fread(groups_from)
  stopifnot(all(c("sample", "group_id") %in% names(gd)))
  sample_to_group <- setNames(gd$group_id, gd$sample)
  group_label <- "custom"
} else if (file.exists(NGSADMIX_Q_FILE)) {
  message("[PLOT_02] Loading NGSadmix Q at K=", CANONICAL_K, ": ", NGSADMIX_Q_FILE)
  Q <- as.matrix(fread(NGSADMIX_Q_FILE, header = FALSE))
  if (!file.exists(SAMPLE_LIST)) stop("SAMPLE_LIST not found: ", SAMPLE_LIST)
  snames <- as.character(fread(SAMPLE_LIST, header = FALSE)[[1]])
  if (length(snames) != nrow(Q)) {
    message("  [warn] sample list (", length(snames),
            ") != Q rows (", nrow(Q), ") — using first min")
    n <- min(length(snames), nrow(Q))
    snames <- snames[seq_len(n)]; Q <- Q[seq_len(n), , drop = FALSE]
  }
  argmax_k <- apply(Q, 1, which.max)
  sample_to_group <- setNames(paste0("K", CANONICAL_K, "_Q", argmax_k), snames)
  group_label <- paste0("ancestry_K", CANONICAL_K)
} else {
  stop("No groups available — neither --groups-from nor NGSADMIX_Q_FILE")
}

groups_in_use <- sort(unique(unname(sample_to_group)))
group_sizes <- table(sample_to_group)
message("[PLOT_02] ", length(groups_in_use), " groups: ",
        paste(groups_in_use, collapse = " "))

# ── Build inversion × group HOM_INV-frequency matrix ─────────────────────────
cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))

read_group_members <- function(gid) {
  fp <- file.path(SAMPLE_REGISTRY, "groups", paste0(gid, ".txt"))
  if (!file.exists(fp)) return(character())
  trimws(readLines(fp))
}

rows <- list()
for (i in seq_len(nrow(cands))) {
  cid <- cands$candidate_id[i]
  hom_inv <- read_group_members(paste0("inv_", cid, "_HOM_INV"))
  if (length(hom_inv) == 0) next
  # For each ancestry group, fraction that is HOM_INV
  for (g in groups_in_use) {
    g_samples <- names(sample_to_group)[sample_to_group == g]
    n_carriers <- sum(g_samples %in% hom_inv)
    rows[[length(rows) + 1]] <- data.table(
      candidate_id = cid,
      chrom = cands$chrom[i],
      length_mb = (cands$end_bp[i] - cands$start_bp[i]) / 1e6,
      group = g,
      n_group = length(g_samples),
      n_HOM_INV_in_group = n_carriers,
      freq_HOM_INV_in_group = if (length(g_samples) > 0) n_carriers / length(g_samples) else 0)
  }
}
if (length(rows) == 0) {
  message("[PLOT_02] No inversion carrier groups found in sample_registry — skip")
  quit(status = 0)
}
dt <- rbindlist(rows)
fwrite(dt, file.path(EXTRAS_TBL_DIR, "PLOT_02_inversion_sharing_table.tsv"),
       sep = "\t")

# ── Order rows by chrom + start, columns by group label ──────────────────────
dt[, candidate_id := factor(candidate_id,
                              levels = unique(cands$candidate_id))]
dt[, group := factor(group, levels = groups_in_use)]

# ── Plot ─────────────────────────────────────────────────────────────────────

p <- ggplot(dt, aes(x = group, y = candidate_id, fill = freq_HOM_INV_in_group)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(name = "frac HOM_INV\nin group",
                        limits = c(0, 1), option = "magma") +
  labs(title = "Inversion sharing across groups",
       subtitle = paste0("Groups defined by ", group_label,
                          " (", length(groups_in_use), " groups, ",
                          length(unique(dt$candidate_id)), " candidates)"),
       x = paste0("Group  (", group_label, ")"),
       y = "Candidate inversion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))

# Choose plot height proportional to # candidates
h <- max(5, min(20, 1 + 0.18 * length(unique(dt$candidate_id))))
ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_02_inversion_sharing_across_groups.pdf"),
       p, width = max(7, 1 + 0.4 * length(groups_in_use)),
       height = h, device = cairo_pdf)
message("[PLOT_02] Done — wrote PLOT_02_inversion_sharing_across_groups.pdf")
