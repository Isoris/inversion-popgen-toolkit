#!/usr/bin/env Rscript
# =============================================================================
# fig4e_class_along_inversion.R — stacked-proportion plot along an inversion
# =============================================================================
#
# Purpose
# -------
# For one inversion candidate, show how the partition of the cohort into
# karyotype classes (INV_INV / INV_nonINV / MID) shifts across genomic
# position inside the inversion ± flanks. The result is a stacked-area
# plot with x = bp position and y = proportion of samples in each class
# at each window.
#
# What this isn't
# ---------------
# This is NOT "ancestry-along-inversion" in the ADMIXTURE sense — the
# catfish repo's Module 2B runs genome-wide NGSadmix only, not per-window
# local ancestry, so a Q-matrix per window does not exist. An upstream
# local-ancestry compute (ELAI / RFMix / windowed NGSadmix) would be
# required to produce it.
#
# What this IS: the in-repo analog using GHSL v6 karyotype calls, which
# answers the biologically related question "does the cohort partition
# shift across this inversion" — if it does, the region is composite;
# if it stays flat, the inversion is a single coherent block.
#
# Inputs
# ------
#   --candidate_bed <file>     BED: chrom, start, end (1 row)
#   --karyo_calls <tsv.gz>     GHSL v6 snake3v6_karyotype_calls.tsv.gz
#                              Columns: sample_id, window_start, window_end
#                              (per-chrom window indices), n_windows, call,
#                              mean_rank
#   --window_coords <tsv>      Full per-chrom window coordinates: chrom,
#                              global_window_id, start_bp, end_bp
#   --out <file.pdf>           Output PDF. PNG alongside.
#
# Optional
# --------
#   --flank_kb <int>           Flank around the inversion (default 500)
#   --title <str>              Plot title (default: chrom + coords)
#   --smooth <int>             Smoothing window in # of bins (default 0)
#
# Output
# ------
#   <out>.pdf                  Stacked area plot.
#   <out>.png                  Same, 200 dpi.
#   <out>.class_proportions.tsv   Underlying per-window proportions.
#
# Known limits
# ------------
# - Samples whose karyotype calls don't touch a given window contribute
#   nothing to that window's denominator. If call coverage varies across
#   the region, the "proportion" is really "proportion of called samples"
#   per window, not "proportion of the cohort." A sanity line is drawn
#   showing the n_called_samples per window.
# - MID samples here = samples with no LOW / HIGH run at this window per
#   the STEP_C04b classifier's rules. These are often uninformative in
#   short-read data and should not be over-interpreted.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ── CLI ──
args <- commandArgs(trailingOnly = TRUE)
candidate_bed_path <- NULL
karyo_calls_path   <- NULL
window_coords_path <- NULL
out_path           <- "class_along_inversion.pdf"
flank_kb           <- 500L
plot_title         <- NULL
smooth_bins        <- 0L

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--candidate_bed" = { candidate_bed_path <- args[i+1]; i <- i+2 },
    "--karyo_calls"   = { karyo_calls_path   <- args[i+1]; i <- i+2 },
    "--window_coords" = { window_coords_path <- args[i+1]; i <- i+2 },
    "--out"           = { out_path           <- args[i+1]; i <- i+2 },
    "--flank_kb"      = { flank_kb           <- as.integer(args[i+1]); i <- i+2 },
    "--title"         = { plot_title         <- args[i+1]; i <- i+2 },
    "--smooth"        = { smooth_bins        <- as.integer(args[i+1]); i <- i+2 },
    { message("[class_along] unknown arg: ", args[i]); i <- i+1 }
  )
}
stopifnot(!is.null(candidate_bed_path),
          !is.null(karyo_calls_path),
          !is.null(window_coords_path))

# ── Load inputs ──
cand <- fread(candidate_bed_path, header = FALSE,
              col.names = c("chrom", "start", "end"))[1, ]
ext_start <- max(0, cand$start - flank_kb * 1000)
ext_end   <- cand$end + flank_kb * 1000

win_all <- fread(window_coords_path)
stopifnot(all(c("chrom", "global_window_id", "start_bp", "end_bp") %in%
              names(win_all)))
win_all <- win_all[chrom == cand$chrom]
setorder(win_all, start_bp)
win_all[, chrom_idx := .I]

win_dt <- win_all[end_bp >= ext_start & start_bp <= ext_end]
if (nrow(win_dt) == 0) stop("[class_along] no windows overlap the interval")
n_win <- nrow(win_dt)
region_chrom_lo <- min(win_dt$chrom_idx)
region_chrom_hi <- max(win_dt$chrom_idx)

karyo <- fread(karyo_calls_path)
stopifnot(all(c("sample_id", "window_start", "window_end", "call") %in%
              names(karyo)))
karyo <- karyo[!(window_end < region_chrom_lo | window_start > region_chrom_hi)]
if (nrow(karyo) == 0) stop("[class_along] no karyotype runs overlap the interval")

# ── Tally per-window call counts ──
# For each window position (local index 1..n_win), count how many samples
# have an INV_INV run there vs an INV_nonINV run. Samples without any
# overlapping run contribute to "MID / no call" per window.

count_mat <- matrix(0L, nrow = 3, ncol = n_win,
                    dimnames = list(c("INV_INV", "INV_nonINV", "MID"), NULL))

# Map each run's per-chrom window range to region-local indices
for (ri in seq_len(nrow(karyo))) {
  local_hits <- win_dt[chrom_idx >= karyo$window_start[ri] &
                         chrom_idx <= karyo$window_end[ri], chrom_idx]
  if (length(local_hits) == 0) next
  # Convert chrom_idx → region-local 1..n_win
  local_idx <- match(local_hits, win_dt$chrom_idx)
  row_key <- karyo$call[ri]
  if (!row_key %in% rownames(count_mat)) next
  count_mat[row_key, local_idx] <- count_mat[row_key, local_idx] + 1L
}

# Samples with any run anywhere in this region — the denominator for
# proportions needs to distinguish "sample present but no call here"
# (→ MID) from "sample absent from the cohort" (→ not counted).
samples_in_region <- unique(karyo$sample_id)
n_samples <- length(samples_in_region)
# MID at window w = n_samples - INV_INV[w] - INV_nonINV[w]
count_mat["MID", ] <- pmax(0L, n_samples -
                             count_mat["INV_INV", ] -
                             count_mat["INV_nonINV", ])

# ── Convert to proportions ──
total <- colSums(count_mat)
total[total == 0] <- 1L
prop_mat <- sweep(count_mat, 2, total, "/")

# Optional smoothing (rolling mean across bins)
if (smooth_bins > 1) {
  sm <- function(x, w) {
    if (length(x) <= w) return(x)
    stats::filter(x, rep(1/w, w), sides = 2) |> as.numeric() |>
      (\(y) { y[is.na(y)] <- x[is.na(y)]; y })()
  }
  prop_mat <- t(apply(prop_mat, 1, sm, w = smooth_bins))
}

# ── Build plot dataframe ──
df <- data.table(
  mid_bp     = (win_dt$start_bp + win_dt$end_bp) / 2,
  INV_INV    = prop_mat["INV_INV", ],
  MID        = prop_mat["MID", ],
  INV_nonINV = prop_mat["INV_nonINV", ]
)
df_long <- melt(df, id.vars = "mid_bp",
                variable.name = "class", value.name = "prop")
# Enforce stack order: INV_INV bottom, MID middle, INV_nonINV top
df_long[, class := factor(class, levels = c("INV_nonINV", "MID", "INV_INV"))]

class_pal <- c(INV_INV = "#2166AC", MID = "#BBBBBB", INV_nonINV = "#D6604D")

if (is.null(plot_title)) {
  plot_title <- sprintf("%s:%d-%d  (N=%d samples, flank ±%d kb)",
                        cand$chrom, cand$start, cand$end, n_samples, flank_kb)
}

p <- ggplot(df_long, aes(x = mid_bp / 1e6, y = prop, fill = class)) +
  geom_area(color = "white", linewidth = 0.1) +
  geom_vline(xintercept = c(cand$start, cand$end) / 1e6,
             linetype = "dashed", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = class_pal) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = plot_title,
       x = sprintf("%s position (Mb)", cand$chrom),
       y = "Proportion of samples",
       fill = "Karyotype") +
  theme_classic(base_size = 11) +
  theme(plot.title   = element_text(face = "bold", size = 11),
        legend.position = "top",
        legend.title    = element_text(size = 9),
        legend.key.size = grid::unit(0.4, "cm"))

# ── Write outputs ──
pw <- 10; ph <- 3.5
ggsave(out_path, p, width = pw, height = ph, device = cairo_pdf)
png_path <- paste0(tools::file_path_sans_ext(out_path), ".png")
ggsave(png_path, p, width = pw, height = ph, dpi = 200)

tsv_path <- paste0(tools::file_path_sans_ext(out_path), ".class_proportions.tsv")
out_df <- data.table(
  chrom    = cand$chrom,
  start_bp = win_dt$start_bp,
  end_bp   = win_dt$end_bp,
  mid_bp   = (win_dt$start_bp + win_dt$end_bp) / 2,
  INV_INV    = prop_mat["INV_INV", ],
  MID        = prop_mat["MID", ],
  INV_nonINV = prop_mat["INV_nonINV", ],
  n_called   = total
)
fwrite(out_df, tsv_path, sep = "\t")

message("[class_along] Wrote ", out_path)
message("[class_along] Wrote ", png_path)
message("[class_along] Wrote ", tsv_path)
