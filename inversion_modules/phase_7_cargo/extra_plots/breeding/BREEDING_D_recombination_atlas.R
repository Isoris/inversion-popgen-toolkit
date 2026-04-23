#!/usr/bin/env Rscript

# =============================================================================
# BREEDING_D_recombination_atlas.R
#
# Genome-wide map of where breeding-by-recombination is achievable vs locked.
# Combines three signals:
#
#   1. Inversion intervals (from candidate table) — by construction, these
#      are recombination-suppressed loci between arrangements.
#   2. ROH co-occurrence — regions where many samples share runs of
#      homozygosity. Long ROH = cryptic linkage block, often the result of
#      historical recombination suppression and/or recent ancestry sharing.
#   3. Inter-arrangement FST per site (from cargo, if available) — high FST
#      between HOM_REF and HOM_INV inside an inversion = recombination
#      suppression is real and effective. Low FST = the inversion is leaky.
#
# Output:
#   ${EXTRAS_FIG_DIR}/BREEDING_D_recombination_atlas.pdf
#       2-track manhattan: inversion span (top) + ROH density (bottom)
#   ${EXTRAS_TBL_DIR}/BREEDING_D_locked_regions.tsv
#       per-window summary: chrom, window, in_inversion, mean_FST_inv, ROH_density
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

SNAKE_CAND_FILE     <- Sys.getenv("SNAKE_CAND_FILE")
ROH_BED             <- Sys.getenv("ROH_BED", "")
REF_FAI             <- Sys.getenv("REF_FAI")
RESULTS_REGISTRY_DIR <- Sys.getenv("RESULTS_REGISTRY_DIR")
CARGO_BURDEN_DIR    <- Sys.getenv("CARGO_BURDEN_DIR")
EXTRAS_FIG_DIR      <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR      <- Sys.getenv("EXTRAS_TBL_DIR")
WIN <- 1e6

if (!file.exists(SNAKE_CAND_FILE) || !file.exists(REF_FAI)) {
  message("[BREEDING_D] [skip] missing inputs"); quit(status = 0)
}

cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))
fai <- fread(REF_FAI, header = FALSE,
             col.names = c("chrom", "len", "off", "x1", "x2"))
chr_order <- fai$chrom
chr_off   <- setNames(c(0, cumsum(fai$len[-nrow(fai)])), chr_order)

# ── Genome windows ───────────────────────────────────────────────────────────
windows <- rbindlist(lapply(seq_len(nrow(fai)), function(i) {
  starts <- seq(0, fai$len[i] - 1, by = WIN)
  data.table(chrom = fai$chrom[i],
              start = starts,
              end = pmin(starts + WIN, fai$len[i]))
}))

# ── Inversion overlap per window ─────────────────────────────────────────────
windows[, in_inversion := FALSE]
cands_clean <- cands[, .(chrom, start = as.integer(start_bp),
                          end = as.integer(end_bp))]
setkey(cands_clean, chrom, start, end)
ov <- foverlaps(windows[, .(chrom, start, end)], cands_clean, type = "any",
                 which = TRUE, mult = "first")
windows[, in_inversion := !is.na(ov$yid)]

# ── ROH density per window ───────────────────────────────────────────────────
windows[, roh_density := 0]
if (nzchar(ROH_BED) && file.exists(ROH_BED)) {
  roh <- if (endsWith(ROH_BED, ".gz")) {
    fread(cmd = paste0("zcat ", shQuote(ROH_BED)), header = FALSE)
  } else {
    fread(ROH_BED, header = FALSE)
  }
  if (ncol(roh) >= 3) {
    setnames(roh, 1:3, c("chrom", "start", "end"))
    roh[, start := as.integer(start)]; roh[, end := as.integer(end)]
    setkey(roh, chrom, start, end)
    # Count ROH segments overlapping each window
    ov2 <- foverlaps(windows[, .(chrom, start, end, win_idx = .I)], roh,
                      type = "any", nomatch = NA)
    roh_counts <- ov2[!is.na(start), .N, by = win_idx]
    windows[roh_counts$win_idx, roh_density := roh_counts$N]
  }
} else {
  message("  [info] no ROH_BED — track will be empty")
}

# ── Optional: per-inversion mean FST from cargo ──────────────────────────────
windows[, mean_fst_inv := NA_real_]
fst_files <- character()
if (nzchar(RESULTS_REGISTRY_DIR) && dir.exists(RESULTS_REGISTRY_DIR)) {
  fst_files <- list.files(file.path(RESULTS_REGISTRY_DIR,
                                       "pairwise", "cargo_fst_per_site"),
                            pattern = "\\.tsv(\\.gz)?$", recursive = TRUE,
                            full.names = TRUE)
}
if (length(fst_files) == 0 && nzchar(CARGO_BURDEN_DIR)) {
  fst_files <- list.files(CARGO_BURDEN_DIR, pattern = "fst_per_site\\.tsv$",
                            recursive = TRUE, full.names = TRUE)
}
if (length(fst_files) > 0) {
  message("[BREEDING_D] reading ", length(fst_files), " FST files")
  fst_all <- rbindlist(lapply(fst_files, function(f) {
    d <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(d) || !"var_key" %in% names(d)) return(NULL)
    d[, c("chrom", "pos") := tstrsplit(var_key, ":", keep = 1:2, type.convert = TRUE)]
    d[, .(chrom, pos = as.integer(pos),
          fst = suppressWarnings(as.numeric(fst_hudson)))]
  }), fill = TRUE)
  fst_all <- fst_all[!is.na(fst)]
  fst_all[, win_start := floor(pos / WIN) * WIN]
  fst_win <- fst_all[, .(mean_fst = mean(fst, na.rm = TRUE)),
                     by = .(chrom, start = win_start)]
  windows <- merge(windows, fst_win, by = c("chrom", "start"), all.x = TRUE)
  windows[!is.na(mean_fst), mean_fst_inv := mean_fst]
  windows[, mean_fst := NULL]
}

# ── Save table ───────────────────────────────────────────────────────────────
windows[, abs_x := chr_off[chrom] + start + WIN/2]
fwrite(windows[, .(chrom, start, end, in_inversion, roh_density, mean_fst_inv)],
       file.path(EXTRAS_TBL_DIR, "BREEDING_D_locked_regions.tsv"), sep = "\t")

# ── Plot ─────────────────────────────────────────────────────────────────────
chr_mid <- data.table(
  chrom = chr_order,
  abs_mid = chr_off + fai$len / 2)
windows[, chrom := factor(chrom, levels = chr_order)]
windows[, chr_band := match(chrom, chr_order) %% 2]

# Top: inversion + FST
pA <- ggplot(windows, aes(x = abs_x)) +
  geom_rect(data = windows[in_inversion == TRUE],
             aes(xmin = abs_x - WIN/2, xmax = abs_x + WIN/2,
                  ymin = -Inf, ymax = Inf),
             fill = "#9467bd", alpha = 0.18, inherit.aes = FALSE) +
  {if ("mean_fst_inv" %in% names(windows) && any(!is.na(windows$mean_fst_inv)))
    geom_point(aes(y = mean_fst_inv, color = factor(chr_band)),
               size = 0.7, alpha = 0.8, na.rm = TRUE)
   else geom_blank()} +
  scale_color_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"), guide = "none") +
  scale_x_continuous(breaks = chr_mid$abs_mid, labels = chr_mid$chrom,
                       expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(0, 1), name = "mean inv-vs-arrangement FST") +
  labs(title = "A. Inversion span (purple shading) + per-site FST",
       subtitle = "Purple = recombination-locked region; high FST = effective lockdown") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        plot.title = element_text(face = "bold"))

# Bottom: ROH density
pB <- ggplot(windows, aes(x = abs_x, y = roh_density,
                            color = factor(chr_band))) +
  geom_rect(data = windows[in_inversion == TRUE],
             aes(xmin = abs_x - WIN/2, xmax = abs_x + WIN/2,
                  ymin = -Inf, ymax = Inf),
             fill = "#9467bd", alpha = 0.18, inherit.aes = FALSE) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"), guide = "none") +
  scale_x_continuous(breaks = chr_mid$abs_mid, labels = chr_mid$chrom,
                       expand = c(0.005, 0.005)) +
  labs(title = "B. ROH density (segments per Mb window)",
       subtitle = "Hotspots = cryptic linkage blocks where breeding-by-recombination is slow",
       x = NULL, y = "ROH segments / window") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        plot.title = element_text(face = "bold"))

fig <- pA / pB +
  plot_annotation(title = "Recombination-suppressed atlas",
                   subtitle = "Where breeding-by-recombination is achievable (clear) vs locked (purple/dense)")
ggsave(file.path(EXTRAS_FIG_DIR, "BREEDING_D_recombination_atlas.pdf"),
       fig, width = 18, height = 9, device = cairo_pdf)

# Summary numbers
inv_bp <- sum(windows[in_inversion == TRUE, end - start])
total_bp <- sum(windows[, end - start])
message("[BREEDING_D] Done")
message(sprintf("  Inversion-locked: %.1f Mb (%.1f%% of genome)",
                 inv_bp / 1e6, 100 * inv_bp / total_bp))
message(sprintf("  ROH-dense (>5 segments/window): %d windows",
                 windows[roh_density > 5, .N]))
