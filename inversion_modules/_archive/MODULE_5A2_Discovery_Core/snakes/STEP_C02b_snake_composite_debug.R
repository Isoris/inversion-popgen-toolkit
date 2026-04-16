#!/usr/bin/env Rscript

# =============================================================================
# STEP_C02b_snake_composite_debug.R  (v1.0)
#
# COMPOSITE DEBUG VIEWS per merge-A region.
#
# For each merge-A candidate region (from Snake 1), produces a multi-panel
# composite PDF showing what Snake 2 sees at every level:
#
#   Panel 1 — LOCAL PCA (PC1 × PC2) facet wrap: one facet per window,
#             colored by 3-band assignment (tailA/middle/tailB).
#             Strip background = region color.
#
#   Panel 2 — DOSAGE PROFILE HEATMAP: samples × markers (markers in genomic
#             order), discretized 0/1/2, sorted by band then by within-band
#             similarity. Shows the raw "012 string" for visual inspection.
#
#   Panel 3 — fKNN SIMILARITY HEATMAP: sample × sample, colored by fuzzy kNN
#             weight overlap. Shows community structure within and across bands.
#             Averaged across windows in the region.
#
#   Panel 4 — BAND ASSIGNMENT RIBBON: windows (x) × samples (y), colored by
#             band membership per window. Shows where samples switch bands
#             (composite signal).
#
#   Panel 5 — SCORE CURVES: per-window combined score, middle stability,
#             profile concordance, and inv-likeness stacked.
#
# INPUTS:
#   <step10_outprefix>.mds.rds           — window grid + PC scores
#   <snake1_dir>/snake_regions.tsv.gz    — merge_A region definitions
#   <snake1_dir>/snake_windows.tsv.gz    — per-window snake membership
#   <snake1_dir>/snake_inv_likeness.tsv.gz — inv-likeness scores (optional)
#   <snake2_dir>/snake2_track.tsv.gz     — per-window community scores
#   <snake2_dir>/snake2_band_assignments.tsv.gz — per-window band data
#   <dosage_dir>/<chr>.dosage.tsv.gz     — raw dosage for profile heatmaps
#   <dosage_dir>/<chr>.sites.tsv.gz      — marker positions
#
# OUTPUTS:
#   <outdir>/composite_<chr>_<region_id>.pdf — one per qualifying region
#   <outdir>/composite_summary.tsv           — which regions were plotted
#
# Usage:
#   Rscript STEP_C02b_snake_composite_debug.R \
#     <step10_outprefix> <snake1_dir> <snake2_dir> <dosage_dir> <outdir> \
#     [--min_windows 10] [--samples_ind <path>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript STEP_C02b_snake_composite_debug.R ",
       "<step10_outprefix> <snake1_dir> <snake2_dir> <dosage_dir> <outdir> ",
       "[--min_windows 10] [--samples_ind <path>]")
}

step10_prefix <- args[1]
snake1_dir    <- args[2]
snake2_dir    <- args[3]
dosage_dir    <- args[4]
outdir        <- args[5]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

MIN_WINDOWS   <- 10L
samples_ind   <- NULL

i <- 6L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--min_windows" && i < length(args)) {
    MIN_WINDOWS <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--samples_ind" && i < length(args)) {
    samples_ind <- args[i + 1]; i <- i + 2L
  } else { i <- i + 1L }
}

# =============================================================================
# LOAD
# =============================================================================

safe_load <- function(path) {
  if (!file.exists(path)) { message("  [skip] ", path); return(NULL) }
  tryCatch(fread(path), error = function(e) { message("  [fail] ", path); NULL })
}

message("[COMPOSITE] Loading data...")

mds_obj <- readRDS(paste0(step10_prefix, ".mds.rds"))
per_chr <- mds_obj$per_chr

reg_dt   <- safe_load(file.path(snake1_dir, "snake_regions.tsv.gz"))
win_dt   <- safe_load(file.path(snake1_dir, "snake_windows.tsv.gz"))
inv_dt   <- safe_load(file.path(snake1_dir, "snake_inv_likeness.tsv.gz"))
s2_track <- safe_load(file.path(snake2_dir, "snake2_track.tsv.gz"))
s2_bands <- safe_load(file.path(snake2_dir, "snake2_band_assignments.tsv.gz"))

sample_names <- NULL
if (!is.null(samples_ind) && file.exists(samples_ind)) {
  sample_names <- fread(samples_ind, header = FALSE)[[1]]
  sample_names <- sample_names[nchar(sample_names) > 0]
}

if (is.null(reg_dt) || nrow(reg_dt) == 0) stop("No snake regions found")
if (is.null(s2_track) || nrow(s2_track) == 0) stop("No snake2 track found")

# =============================================================================
# PALETTE
# =============================================================================

band_pal <- c(tailA = "#2563eb", middle = "#f59e0b", tailB = "#dc2626", unknown = "#94a3b8")
reg_pal  <- c("#dc2626", "#2563eb", "#16a34a", "#d97706", "#7c3aed",
              "#059669", "#e11d48", "#0284c7", "#65a30d", "#c026d3")

# Dosage heatmap palette (0=blue, 1=yellow, 2=red, NA=grey)
dos_pal <- c("0" = "#2563eb", "1" = "#f59e0b", "2" = "#dc2626", "-1" = "#e2e8f0")

# =============================================================================
# PER-REGION COMPOSITE
# =============================================================================

# Focus on merge_A regions
merge_a <- reg_dt[grepl("merge_A", family)]
if (nrow(merge_a) == 0) merge_a <- reg_dt
merge_a <- merge_a[n_windows >= MIN_WINDOWS]

message("[COMPOSITE] ", nrow(merge_a), " regions with >= ", MIN_WINDOWS, " windows")

summary_rows <- list()

for (ri in seq_len(nrow(merge_a))) {
  r <- merge_a[ri]
  chr <- r$chrom
  chr_obj <- per_chr[[chr]]
  if (is.null(chr_obj)) next

  dt_chr <- as.data.table(chr_obj$out_dt)[order(start_bp)]
  dt_reg <- dt_chr[start_bp >= r$start_bp & end_bp <= r$end_bp]
  if (nrow(dt_reg) < MIN_WINDOWS) next

  wids <- dt_reg$global_window_id
  n_win <- nrow(dt_reg)
  reg_color <- reg_pal[((ri - 1L) %% length(reg_pal)) + 1L]
  reg_label <- paste0(chr, " ", r$family, " #", r$snake_id,
                      " (", round(r$start_bp/1e6, 2), "–", round(r$end_bp/1e6, 2), " Mb)")
  message("[COMPOSITE] ", reg_label, " — ", n_win, " windows")

  plots <- list()

  # ── PANEL 1: Local PCA facet wrap ────────────────────────────────────
  pc1_cols <- grep("^PC_1_", names(dt_reg), value = TRUE)
  pc2_cols <- grep("^PC_2_", names(dt_reg), value = TRUE)

  if (length(pc1_cols) > 0) {
    facet_rows <- list()
    for (wi in seq_len(n_win)) {
      scores1 <- as.numeric(dt_reg[wi, ..pc1_cols])
      scores2 <- if (length(pc2_cols) > 0) as.numeric(dt_reg[wi, ..pc2_cols]) else rep(0, length(scores1))
      valid <- is.finite(scores1)
      if (sum(valid) < 20) next

      km <- tryCatch(kmeans(scores1[valid], centers = 3, nstart = 3), error = function(e) NULL)
      if (is.null(km)) next
      co <- order(km$centers[, 1])
      rm <- setNames(c("tailA", "middle", "tailB"), as.character(co))
      roles <- rm[as.character(km$cluster)]

      pos_mb <- round((dt_reg$start_bp[wi] + dt_reg$end_bp[wi]) / 2e6, 2)
      facet_rows[[length(facet_rows) + 1]] <- data.table(
        wlabel = paste0(pos_mb, " Mb"), PC1 = scores1[valid], PC2 = scores2[valid],
        band = roles, worder = wi
      )
    }
    if (length(facet_rows) > 0) {
      fdt <- rbindlist(facet_rows)
      fdt[, wlabel := factor(wlabel, levels = unique(wlabel[order(worder)]))]
      nc <- min(6L, length(unique(fdt$wlabel)))
      plots$pca <- ggplot(fdt, aes(PC1, PC2, color = band)) +
        geom_point(size = 0.3, alpha = 0.5) +
        scale_color_manual(values = band_pal) +
        facet_wrap(~ wlabel, ncol = nc, scales = "free") +
        labs(title = paste0("Panel 1: Local PCA — ", reg_label)) +
        theme_minimal(base_size = 6) +
        theme(strip.background = element_rect(fill = reg_color, color = NA),
              strip.text = element_text(color = "white", size = 4, face = "bold"),
              legend.position = "bottom")
    }
  }

  # ── PANEL 2: Dosage profile heatmap ──────────────────────────────────
  dos_file <- file.path(dosage_dir, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(dosage_dir, paste0(chr, ".sites.tsv.gz"))

  if (file.exists(dos_file) && file.exists(sites_file)) {
    dos <- fread(dos_file)
    sites <- fread(sites_file)
    dos_cols <- setdiff(names(dos), "marker")
    if (!is.null(sample_names) && length(sample_names) == length(dos_cols)) {
      setnames(dos, old = dos_cols, new = sample_names)
      dos_cols <- sample_names
    }

    mk <- which(sites$pos >= r$start_bp & sites$pos <= r$end_bp)
    if (length(mk) >= 20) {
      Xr <- as.matrix(dos[mk, ..dos_cols])
      storage.mode(Xr) <- "double"
      Gr <- round(Xr); Gr[is.na(Gr)] <- -1L
      storage.mode(Gr) <- "integer"

      # Sort samples: first by band assignment (from first valid window), then by PC1
      first_band <- s2_bands[chrom == chr & global_window_id %in% wids]
      sample_band <- rep("unknown", length(dos_cols))
      if (nrow(first_band) > 0) {
        fb <- first_band[1]
        mid_idx <- as.integer(strsplit(fb$middle_idx_csv, ",")[[1]])
        tA_idx <- as.integer(strsplit(fb$tailA_idx_csv, ",")[[1]])
        tB_idx <- as.integer(strsplit(fb$tailB_idx_csv, ",")[[1]])
        sample_band[mid_idx] <- "middle"
        sample_band[tA_idx] <- "tailA"
        sample_band[tB_idx] <- "tailB"
      }

      # PCA for ordering within bands
      pc <- tryCatch(prcomp(t(Xr[, is.finite(colMeans(Xr))]), center = TRUE)$x[, 1],
                     error = function(e) seq_along(dos_cols))
      samp_order <- order(sample_band, pc)

      # Build long-format heatmap data
      Gr_sorted <- Gr[, samp_order]
      hm_long <- data.table(
        marker_idx = rep(seq_len(nrow(Gr_sorted)), ncol(Gr_sorted)),
        sample_idx = rep(seq_len(ncol(Gr_sorted)), each = nrow(Gr_sorted)),
        geno = as.character(as.vector(Gr_sorted))
      )

      plots$dosage <- ggplot(hm_long, aes(x = marker_idx, y = sample_idx, fill = geno)) +
        geom_raster() +
        scale_fill_manual(values = dos_pal, na.value = "#e2e8f0") +
        labs(title = paste0("Panel 2: Dosage 012 heatmap — ", nrow(Gr_sorted), " markers × ",
                            ncol(Gr_sorted), " samples"),
             x = "Marker (genomic order)", y = "Sample (sorted by band + PC1)") +
        theme_minimal(base_size = 7) +
        theme(legend.position = "bottom", axis.text = element_blank(),
              panel.grid = element_blank())
    }
  }

  # ── PANEL 3: Sample × sample fKNN similarity (averaged across windows) ──
  # Build average distance matrix from Snake 2's correlation distances
  # We reconstruct from the track data: for each window pair, use combined_score
  # as a proxy for community similarity

  # Simpler approach: use the band assignments to build a sample × window
  # membership matrix, then compute sample-pair co-occurrence Jaccard
  reg_bands <- s2_bands[chrom == chr & global_window_id %in% wids]
  if (nrow(reg_bands) >= 3 && length(dos_cols) > 0) {
    n_samp <- length(dos_cols)
    co_mid <- matrix(0, n_samp, n_samp)
    co_any <- matrix(0, n_samp, n_samp)
    n_counted <- 0L

    for (bi in seq_len(nrow(reg_bands))) {
      mid_idx <- as.integer(strsplit(reg_bands$middle_idx_csv[bi], ",")[[1]])
      tA_idx <- as.integer(strsplit(reg_bands$tailA_idx_csv[bi], ",")[[1]])
      tB_idx <- as.integer(strsplit(reg_bands$tailB_idx_csv[bi], ",")[[1]])

      # Co-occurrence: are samples i,j in the same band in this window?
      all_bands <- rep(0L, n_samp)
      all_bands[mid_idx] <- 1L; all_bands[tA_idx] <- 2L; all_bands[tB_idx] <- 3L

      valid <- which(all_bands > 0)
      if (length(valid) < 10) next
      n_counted <- n_counted + 1L

      for (ii in seq_len(length(valid) - 1L)) {
        for (jj in (ii + 1L):length(valid)) {
          si <- valid[ii]; sj <- valid[jj]
          if (all_bands[si] == all_bands[sj]) {
            co_any[si, sj] <- co_any[si, sj] + 1
            co_any[sj, si] <- co_any[sj, si] + 1
          }
          if (all_bands[si] == 1L && all_bands[sj] == 1L) {
            co_mid[si, sj] <- co_mid[si, sj] + 1
            co_mid[sj, si] <- co_mid[sj, si] + 1
          }
        }
      }
    }

    if (n_counted > 0) {
      co_frac <- co_any / n_counted
      diag(co_frac) <- 1

      # Use the sample ordering from panel 2
      co_sorted <- co_frac[samp_order, samp_order]
      co_long <- data.table(
        s1 = rep(seq_len(n_samp), n_samp),
        s2 = rep(seq_len(n_samp), each = n_samp),
        similarity = as.vector(co_sorted)
      )

      plots$fknn <- ggplot(co_long, aes(x = s1, y = s2, fill = similarity)) +
        geom_raster() +
        scale_fill_gradient2(low = "#1e3a5f", mid = "#f59e0b", high = "#dc2626", midpoint = 0.5) +
        labs(title = paste0("Panel 3: Band co-occurrence — ", n_counted, " windows averaged"),
             x = "Sample", y = "Sample") +
        theme_minimal(base_size = 7) +
        theme(legend.position = "bottom", axis.text = element_blank(), panel.grid = element_blank())
    }
  }

  # ── PANEL 4: Band assignment ribbon ──────────────────────────────────
  if (nrow(reg_bands) >= 3 && length(dos_cols) > 0) {
    ribbon_rows <- list()
    for (bi in seq_len(nrow(reg_bands))) {
      mid_idx <- as.integer(strsplit(reg_bands$middle_idx_csv[bi], ",")[[1]])
      tA_idx <- as.integer(strsplit(reg_bands$tailA_idx_csv[bi], ",")[[1]])
      tB_idx <- as.integer(strsplit(reg_bands$tailB_idx_csv[bi], ",")[[1]])

      all_bands <- rep("unknown", n_samp)
      all_bands[mid_idx] <- "middle"
      all_bands[tA_idx] <- "tailA"
      all_bands[tB_idx] <- "tailB"

      pos_mb <- round((reg_bands$start_bp[bi] + reg_bands$end_bp[bi]) / 2e6, 3)
      ribbon_rows[[bi]] <- data.table(
        window_pos = pos_mb, sample_idx = seq_len(n_samp)[samp_order],
        band = all_bands[samp_order]
      )
    }

    rdt <- rbindlist(ribbon_rows)
    plots$ribbon <- ggplot(rdt, aes(x = window_pos, y = sample_idx, fill = band)) +
      geom_raster() +
      scale_fill_manual(values = band_pal) +
      labs(title = "Panel 4: Band membership ribbon (samples × windows)",
           x = "Position (Mb)", y = "Sample") +
      theme_minimal(base_size = 7) +
      theme(legend.position = "bottom", axis.text.y = element_blank(), panel.grid = element_blank())
  }

  # ── PANEL 5: Score curves ────────────────────────────────────────────
  reg_track <- s2_track[chrom == chr & global_window_id %in% wids]
  if (nrow(reg_track) > 0) {
    reg_track[, pos_mb := (start_bp + end_bp) / 2e6]

    curve_long <- rbindlist(list(
      reg_track[, .(pos_mb, value = combined_score, curve = "combined")],
      reg_track[, .(pos_mb, value = mean_middle_stability, curve = "middle_stability")],
      if ("profile_concordance" %in% names(reg_track))
        reg_track[, .(pos_mb, value = profile_concordance, curve = "profile_concordance")]
      else NULL
    ))
    curve_long <- curve_long[is.finite(value)]

    # Add inv_likeness if available
    if (!is.null(inv_dt)) {
      il_reg <- inv_dt[global_window_id %in% wids]
      if (nrow(il_reg) > 0) {
        il_reg <- merge(il_reg, dt_reg[, .(global_window_id, start_bp, end_bp)], by = "global_window_id")
        il_reg[, pos_mb := (start_bp + end_bp) / 2e6]
        curve_long <- rbind(curve_long,
                            il_reg[, .(pos_mb, value = inv_likeness, curve = "inv_likeness")])
      }
    }

    if (nrow(curve_long) > 0) {
      plots$curves <- ggplot(curve_long, aes(x = pos_mb, y = value, color = curve)) +
        geom_line(linewidth = 0.5) + geom_point(size = 0.5) +
        scale_color_manual(values = c(combined = "#334155", middle_stability = "#f59e0b",
                                      profile_concordance = "#2563eb", inv_likeness = "#16a34a")) +
        facet_wrap(~ curve, ncol = 1, scales = "free_y") +
        labs(title = "Panel 5: Score curves", x = "Position (Mb)", y = "Score") +
        theme_minimal(base_size = 7) +
        theme(legend.position = "none", strip.text = element_text(face = "bold", size = 6))
    }
  }

  # ── ASSEMBLE COMPOSITE PDF ───────────────────────────────────────────
  if (length(plots) == 0) next

  fname <- paste0("composite_", chr, "_", r$family, "_", r$snake_id, ".pdf")
  fpath <- file.path(outdir, fname)

  n_panels <- length(plots)
  panel_h <- c(pca = 5, dosage = 4, fknn = 4, ribbon = 3, curves = 4)
  total_h <- sum(panel_h[names(plots)])

  pdf(fpath, width = 14, height = min(total_h, 30))
  for (p in plots) {
    print(p)
  }
  dev.off()
  message("  → ", fpath, " (", n_panels, " panels)")

  summary_rows[[length(summary_rows) + 1]] <- data.table(
    chrom = chr, family = r$family, snake_id = r$snake_id,
    start_mb = round(r$start_bp / 1e6, 2), end_mb = round(r$end_bp / 1e6, 2),
    n_windows = n_win, n_panels = n_panels, file = fname
  )
}

# Write summary
if (length(summary_rows) > 0) {
  summary_dt <- rbindlist(summary_rows)
  fwrite(summary_dt, file.path(outdir, "composite_summary.tsv"), sep = "\t")
  message("\n[DONE] ", nrow(summary_dt), " composite debug PDFs → ", outdir)
} else {
  message("\n[DONE] No regions qualified for composite debug views")
}
