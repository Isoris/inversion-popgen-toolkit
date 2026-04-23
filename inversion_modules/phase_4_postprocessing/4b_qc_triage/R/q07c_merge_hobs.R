#!/usr/bin/env Rscript
# =============================================================================
# q07c_merge_hobs.R — merge per-group hobs_windower output into one wide TSV
# =============================================================================
# Input: three <chrom>.<group>.win<scale>.tsv files from hobs_windower
#        (one each for HOM1, HET, HOM2 at a primary scale).
#
# Output wide TSV columns (one row per window):
#   chrom, start_bp, end_bp, n_sites,
#   Hobs_HOM1, Hexp_HOM1, HoverE_HOM1,
#   Hobs_HET,  Hexp_HET,  HoverE_HET,
#   Hobs_HOM2, Hexp_HOM2, HoverE_HOM2
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)[1]
  if (is.na(i) || i >= length(args)) return(default)
  args[i + 1]
}
CHROM    <- get_arg("--chrom")
WIN_DIR  <- get_arg("--win_dir")
SCALE    <- get_arg("--scale", "10kb")
OUT      <- get_arg("--out")

stopifnot(!is.null(CHROM), !is.null(WIN_DIR), !is.null(OUT))

# hobs_windower emits <prefix>.win<scale>.tsv (no .gz). Read each, ratio per
# window, then merge on (start_bp, end_bp).
read_group <- function(group) {
  f <- file.path(WIN_DIR, sprintf("%s.%s.win%s.tsv", CHROM, group, SCALE))
  if (!file.exists(f)) {
    message("[q07c_merge] missing ", f)
    return(NULL)
  }
  d <- fread(f)
  # hobs_windower emits columns: chrom, start, end, n_sites,
  #   Hobs_mean, Hobs_median, Hobs_sd,
  #   Hexp_mean, Hexp_median, Hexp_sd,
  #   F_mean, F_median, F_sd,
  #   n_outliers
  # We keep the mean of Hobs/Hexp as the primary signal plus n_sites for QC.
  stopifnot(all(c("start", "end", "Hobs_mean", "Hexp_mean") %in% names(d)))
  d[, HoverE := fifelse(Hexp_mean > 1e-12, Hobs_mean / Hexp_mean, NA_real_)]
  out <- d[, .(start_bp = start, end_bp = end, n_sites,
               Hobs    = Hobs_mean,
               Hexp    = Hexp_mean,
               HoverE  = HoverE)]
  setnames(out, c("Hobs", "Hexp", "HoverE"),
           paste0(c("Hobs_", "Hexp_", "HoverE_"), group))
  setnames(out, "n_sites", paste0("n_sites_", group))
  out
}

h1 <- read_group("HOM1")
he <- read_group("HET")
h2 <- read_group("HOM2")

# Merge — use HET as the reference window grid since Het is usually the most
# informative group. Missing groups become NA columns.
ref <- he
if (is.null(ref)) ref <- h1
if (is.null(ref)) ref <- h2
if (is.null(ref)) {
  message("[q07c_merge] no group output present, emitting empty TSV")
  fwrite(data.table(chrom = character(), start_bp = integer(), end_bp = integer()),
         OUT, sep = "\t")
  quit(save = "no", status = 0)
}

merged <- ref
for (g in list(h1, he, h2)) {
  if (is.null(g)) next
  merged <- merge(merged, g, by = c("start_bp", "end_bp"), all = TRUE, suffixes = c("", ".dup"))
}
# Drop any accidental .dup columns from the merges
dup_cols <- grep("\\.dup$", names(merged), value = TRUE)
if (length(dup_cols) > 0) merged[, (dup_cols) := NULL]
merged[, chrom := CHROM]
setcolorder(merged, c("chrom", "start_bp", "end_bp", grep("n_sites_|Hobs_|Hexp_|HoverE_", names(merged), value = TRUE)))
setorder(merged, start_bp)

fwrite(merged, OUT, sep = "\t")
message("[q07c_merge] wrote ", OUT, " (", nrow(merged), " windows)")
