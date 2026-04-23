#!/usr/bin/env Rscript
# =============================================================================
# fig4e_cumulative_burden_per_group.R — cumulative burden stacked by group
# =============================================================================
#
# Purpose
# -------
# Burden-accumulation plot à la the Nature paper's Afghani/Jordanian/…
# panel. X = ROH blocks sorted by a chosen order; Y = cumulative count
# of unique deleterious-carrying entities discovered up to block i,
# partitioned by which group first contributed each entity.
#
# Default mode: 3a — entity = gene, source = genes_hit column of
# deleterious_burden_per_ROH.tsv.
#
# This script is parameterized so 3b and 3c are flag flips:
#
#   --entity gene           (default) — count unique genes
#   --entity variant        — count unique Class A/B variants per sample
#                             in each ROH (requires an additional input:
#                             deleterious_variant_genotype_matrix.tsv)
#   --entity gene_inv       — same logic but x-axis is INVERSIONS not ROH
#                             (requires an inversion catalog BED)
#
# 3a (default) and 3b share the same input file. 3c swaps to a different
# region source. All three share the attribution logic below.
#
# Attribution semantics
# ---------------------
# The Nature figure stacks groups so their y-sum equals the cumulative
# total of *distinct* entities (not the sum of per-group counts, which
# would double-count shared genes). This means each gene must be
# attributed to ONE group — the one whose samples FIRST hit it in the
# sorted order. Default: attribution by first-hitting group. Flag
# `--attribution per_group` switches to per-group independent cumulatives
# (stacks sum to more than the unique total; useful as a diagnostic).
#
# Inputs
# ------
#   --burden_roh <tsv>        deleterious_burden_per_ROH.tsv
#                             Required columns: sample_id, chr, roh_start,
#                             roh_end, roh_length, genes_hit
#                             (genes_hit = ';'-joined gene_id list; "none"
#                             if empty)
#
#   --sample_groups <tsv>     sample_id <tab> group_id
#                             The ancestry / lineage assignment for each
#                             sample. Order of appearance defines the
#                             stack order (first-seen on the bottom).
#
#   --out <file.pdf>          Output PDF. PNG + TSV alongside.
#
# Optional
# --------
#   --sort_by length|n_del|position    ROH-block x-axis sort order.
#                                      Default: length (ascending).
#   --attribution first|per_group      See above. Default: first.
#   --palette <csv>                    Group color overrides.
#   --title <str>                      Plot title.
#
# Outputs
# -------
#   <out>.pdf / .png
#   <out>.cumulative.tsv   block_rank, block_id, group, cumulative_count
#
# Edge cases
# ----------
# - Samples in burden file but not in sample_groups: assigned to a
#   synthetic "Unassigned" group. Not silently dropped.
# - Samples in sample_groups but not in burden file: ignored.
# - "none" in genes_hit: skipped (no genes contribute).
# - Duplicate gene_ids within one ROH's genes_hit: deduped.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ── CLI ──
args <- commandArgs(trailingOnly = TRUE)
burden_roh_path <- NULL
sample_groups_path <- NULL
out_path <- "cumulative_burden.pdf"
entity_mode <- "gene"
sort_by <- "length"
attribution <- "first"
palette_override <- NULL
plot_title <- NULL

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--burden_roh"    = { burden_roh_path    <- args[i+1]; i <- i+2 },
    "--sample_groups" = { sample_groups_path <- args[i+1]; i <- i+2 },
    "--out"           = { out_path           <- args[i+1]; i <- i+2 },
    "--entity"        = { entity_mode        <- args[i+1]; i <- i+2 },
    "--sort_by"       = { sort_by            <- args[i+1]; i <- i+2 },
    "--attribution"   = { attribution        <- args[i+1]; i <- i+2 },
    "--palette"       = { palette_override   <- args[i+1]; i <- i+2 },
    "--title"         = { plot_title         <- args[i+1]; i <- i+2 },
    { message("[cumburden] unknown arg: ", args[i]); i <- i+1 }
  )
}
stopifnot(!is.null(burden_roh_path), !is.null(sample_groups_path))

if (!entity_mode %in% c("gene", "variant", "gene_inv")) {
  stop("[cumburden] --entity must be gene / variant / gene_inv")
}
if (!sort_by %in% c("length", "n_del", "position")) {
  stop("[cumburden] --sort_by must be length / n_del / position")
}
if (!attribution %in% c("first", "per_group")) {
  stop("[cumburden] --attribution must be first / per_group")
}

# NOTE: for this first version, only entity_mode == "gene" is fully wired.
# variant and gene_inv modes stub out below with messages explaining what
# additional inputs they'd need. Parameterization-ready but not
# parameter-complete.
if (entity_mode != "gene") {
  stop(sprintf("[cumburden] entity_mode='%s' not yet implemented. ",
               entity_mode),
       "See header docstring for the additional inputs required.")
}

# ── Load inputs ──
burden <- fread(burden_roh_path)
req <- c("sample_id", "chr", "roh_start", "roh_end", "roh_length", "genes_hit")
missing <- setdiff(req, names(burden))
if (length(missing) > 0) stop("[cumburden] burden_roh missing columns: ",
                              paste(missing, collapse = ", "))

groups <- fread(sample_groups_path, header = FALSE,
                col.names = c("sample_id", "group_id"))

# ── Per-ROH sort order ──
# Collapse ROHs to unique blocks across the cohort. The x-axis is the
# set of ROH blocks (not sample × block instances), so two samples
# overlapping the same ROH both contribute to the same x-bin. We define
# a "block" by (chr, roh_start, roh_end) tuple — exact coord match. If
# two samples have similar-but-not-identical boundaries, they count as
# two blocks (the heuristic can be loosened later via bedtools-merge
# if needed).
burden[, block_id := paste0(chr, ":", roh_start, "-", roh_end)]

block_info <- unique(burden[, .(block_id, chr, roh_start, roh_end,
                                roh_length)])
# n_del per block = sum across samples overlapping the block of their
# n_del_variants (if present). Fall back to gene count if not.
if ("n_del_variants" %in% names(burden)) {
  n_del_per_block <- burden[, .(n_del_block = sum(n_del_variants, na.rm = TRUE)),
                             by = block_id]
  block_info <- merge(block_info, n_del_per_block, by = "block_id", all.x = TRUE)
  block_info[is.na(n_del_block), n_del_block := 0L]
} else {
  block_info[, n_del_block := 0L]
}

block_info[, rank := switch(sort_by,
                            length   = frank(roh_length, ties.method = "first"),
                            n_del    = frank(n_del_block, ties.method = "first"),
                            position = frank(paste0(chr, "_",
                                                    sprintf("%012d", roh_start)),
                                             ties.method = "first"))]
setorder(block_info, rank)

# ── Expand burden to per-(block × sample × gene) tuples ──
# genes_hit is ";"-joined (e.g. "GENE1;GENE2;GENE3"). "none" / empty /
# NA means no deleterious genes in this ROH — skip.
burden_sub <- burden[genes_hit != "none" & !is.na(genes_hit) &
                       nzchar(genes_hit)]
if (nrow(burden_sub) == 0) {
  stop("[cumburden] no ROH blocks with deleterious genes found in input")
}

# Explode genes_hit to one row per (sample, block, gene). Use rep() +
# strsplit() with lengths so the semantics are obvious and the final
# data.table has the right row count (not dependent on data.table's `by`
# recycling rules).
gene_lists <- strsplit(burden_sub$genes_hit, ";", fixed = TRUE)
n_per_row <- lengths(gene_lists)
burden_expanded <- data.table(
  sample_id = rep(burden_sub$sample_id, n_per_row),
  block_id  = rep(burden_sub$block_id,  n_per_row),
  gene_id   = unlist(gene_lists)
)

# Attach group
burden_expanded <- merge(burden_expanded,
                         groups[, .(sample_id, group_id)],
                         by = "sample_id", all.x = TRUE)
burden_expanded[is.na(group_id), group_id := "Unassigned"]

# De-dup identical gene hits within a (block, sample)
burden_expanded <- unique(burden_expanded,
                          by = c("block_id", "sample_id", "gene_id"))

# Attach block rank
burden_expanded <- merge(burden_expanded,
                         block_info[, .(block_id, rank)],
                         by = "block_id")

# ── First-hit attribution ──
# For each gene_id, find the (rank, group_id) of the first time any
# sample in a given group hit it. "First" = lowest rank. Ties: pick the
# most-frequent group at that rank (keeps the stack layer stable rather
# than arbitrary).
if (attribution == "first") {
  # For each gene, find the lowest rank at which any sample hit it.
  # Ties at the same rank go to the group that appears FIRST in the
  # sample_groups file order (reproducible across runs rather than
  # row-order-dependent).
  group_order <- unique(groups$group_id)   # preserves input order
  if (!("Unassigned" %in% group_order)) {
    group_order <- c(group_order, "Unassigned")
  }
  burden_expanded[, group_rank := match(group_id, group_order)]
  setorder(burden_expanded, gene_id, rank, group_rank)
  first_hit <- unique(burden_expanded, by = "gene_id")

  # Count genes attributed per (rank, group)
  attr_counts <- first_hit[, .N, by = .(rank, group_id)]
  setnames(attr_counts, "N", "n_new_genes")

  # Fill in all (rank, group) combos so the area plot doesn't have gaps
  all_ranks <- seq_len(max(block_info$rank))
  all_groups <- unique(c(groups$group_id, "Unassigned"))
  grid <- CJ(rank = all_ranks, group_id = all_groups)
  attr_counts <- merge(grid, attr_counts, by = c("rank", "group_id"),
                       all.x = TRUE)
  attr_counts[is.na(n_new_genes), n_new_genes := 0L]

  # Cumulative per group
  setorder(attr_counts, group_id, rank)
  attr_counts[, cum_genes := cumsum(n_new_genes), by = group_id]

} else {
  # per_group: each group's cumulative = unique genes seen by that group
  # up to block rank i. Groups stacks will sum to more than total.
  setorder(burden_expanded, rank)
  per_group_cum <- burden_expanded[, {
    # At each rank, new gene_ids not seen before in this group
    seen <- character(0)
    new_counts <- integer(length(unique(rank)))
    ranks_u <- sort(unique(rank))
    for (r_i in seq_along(ranks_u)) {
      genes_here <- unique(gene_id[rank == ranks_u[r_i]])
      new <- setdiff(genes_here, seen)
      new_counts[r_i] <- length(new)
      seen <- c(seen, new)
    }
    .(rank = ranks_u, n_new_genes = new_counts,
      cum_genes = cumsum(new_counts))
  }, by = group_id]

  # Fill to all ranks for ggplot's geom_area continuity
  all_ranks <- seq_len(max(block_info$rank))
  all_groups <- unique(per_group_cum$group_id)
  grid <- CJ(rank = all_ranks, group_id = all_groups)
  per_group_cum <- merge(grid, per_group_cum,
                         by = c("rank", "group_id"), all.x = TRUE)
  setorder(per_group_cum, group_id, rank)
  # Forward-fill cum_genes within group (ranks where the group saw no
  # new genes keep the previous cum value)
  per_group_cum[, cum_genes := nafill(cum_genes, type = "locf"),
                by = group_id]
  per_group_cum[is.na(cum_genes), cum_genes := 0L]
  attr_counts <- per_group_cum[, .(rank, group_id, n_new_genes, cum_genes)]
}

# ── Palette ──
all_groups <- unique(attr_counts$group_id)
default_pal_base <- c("#B2DF8A", "#F4A6A6", "#F9A65A", "#C23B22",
                      "#B5E2EA", "#2A8A8C", "#C2A4D6", "#888888")
if (!is.null(palette_override)) {
  kv <- strsplit(palette_override, ",", fixed = TRUE)[[1]]
  # "GroupA=#RRGGBB" pairs
  pal <- setNames(vapply(kv, function(s) strsplit(s, "=", fixed = TRUE)[[1]][2],
                         character(1)),
                  vapply(kv, function(s) strsplit(s, "=", fixed = TRUE)[[1]][1],
                         character(1)))
  # Fill missing groups from default cycle
  missing_g <- setdiff(all_groups, names(pal))
  if (length(missing_g) > 0) {
    cycle <- default_pal_base[((seq_along(missing_g) - 1) %%
                                 length(default_pal_base)) + 1]
    pal <- c(pal, setNames(cycle, missing_g))
  }
} else {
  pal <- setNames(default_pal_base[((seq_along(all_groups) - 1) %%
                                      length(default_pal_base)) + 1],
                  all_groups)
}

# ── Plot ──
if (is.null(plot_title)) {
  plot_title <- sprintf("Cumulative deleterious-gene discovery by group  (N blocks = %d)",
                        max(block_info$rank))
}

p <- ggplot(attr_counts, aes(x = rank, y = cum_genes, fill = group_id)) +
  geom_area(position = if (attribution == "first") "stack" else "identity",
            color = "white", linewidth = 0.1,
            alpha = if (attribution == "first") 1 else 0.6) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = pretty(range(attr_counts$rank), n = 6)) +
  labs(title = plot_title,
       subtitle = sprintf("ROH blocks sorted by %s; attribution = %s",
                          sort_by, attribution),
       x = "ROH blocks (sorted)",
       y = "Cumulative gene count (del-carrying)",
       fill = "Group") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        legend.position = "right",
        legend.title    = element_text(size = 9),
        legend.key.size = grid::unit(0.4, "cm"))

# ── Outputs ──
pw <- 7.5; ph <- 5
ggsave(out_path, p, width = pw, height = ph, device = cairo_pdf)
png_path <- paste0(tools::file_path_sans_ext(out_path), ".png")
ggsave(png_path, p, width = pw, height = ph, dpi = 200)

tsv_path <- paste0(tools::file_path_sans_ext(out_path), ".cumulative.tsv")
fwrite(attr_counts, tsv_path, sep = "\t")

message("[cumburden] Wrote ", out_path)
message("[cumburden] Wrote ", png_path)
message("[cumburden] Wrote ", tsv_path)
