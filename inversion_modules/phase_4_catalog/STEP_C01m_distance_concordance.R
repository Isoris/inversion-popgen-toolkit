#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01m_distance_concordance.R
#
# MULTI-SCALE SAMPLE CONCORDANCE ANALYSIS
#
# For each chromosome, computes sample × sample co-grouping matrices at
# multiple genomic distances (d = 80, 160, 320, 640 windows). Reveals
# which sample pairs maintain group co-membership across distance:
#
#   View 1: sim_mat (window × window) — already exists in precomp
#   View 2: sample × sample concordance at each distance d (THIS SCRIPT)
#   View 3: contrast matrix (concordance_d2 - concordance_d1) — highlights
#           pairs that PERSIST (inversion) vs DECAY (family LD)
#
# Biology:
#   - Family LD: sample pairs co-group at short distance, lose co-grouping
#     at longer distance due to recombination
#   - Inversion: sample pairs co-group at ALL distances within the inversion
#     because recombination is suppressed
#   - The distance at which a pair loses co-grouping = LD decay length
#   - Inversion carriers have "infinite" co-grouping within the inversion
#
# Inputs:
#   precomp_dir    — directory with <chr>.precomp.rds files
#   outdir         — output directory
#   [--chrom chr]  — process single chromosome
#   [--distances 80,160,320,640]  — window distances to test
#   [--k 3]        — k for band assignment
#
# Outputs:
#   concordance_<chr>_d<N>.tsv.gz  — 226×226 sample concordance at distance N
#   contrast_<chr>_d<A>_vs_d<B>.tsv.gz — difference matrices
#   persistence_<chr>.tsv.gz       — per-sample-pair persistence profile
#   plots/<chr>_D1_concordance_multiscale.pdf  — 4-panel concordance heatmaps
#   plots/<chr>_D2_contrast.pdf                — contrast matrices
#   plots/<chr>_D3_persistence_profile.pdf     — decay curves
#   plots/<chr>_D4_simmat_multiscale.pdf       — sim_mat at 1x/2x/4x/8x
#   plots/<chr>_D5_inversion_carrier_network.pdf — network of persistent pairs
#   plots/<chr>_D6_sample_decay_heatmap.pdf    — sample × distance decay
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ── chat-13 wiring: registry + write_block_safe helpers ─────────────────────
.try_source <- function(path) {
  if (!is.null(path) && nzchar(path) && file.exists(path)) {
    tryCatch({ source(path); TRUE }, error = function(e) {
      message("[C01m] could not source ", path, ": ", conditionMessage(e))
      FALSE
    })
  } else FALSE
}
.registries_root <- Sys.getenv("REGISTRIES", "")
.pipe_base       <- Sys.getenv("BASE", ".")
.reg_loader <- file.path(ifelse(nzchar(.registries_root),
                                  .registries_root,
                                  file.path(.pipe_base,
                                             "inversion-popgen-toolkit",
                                             "registries")),
                          "api", "R", "registry_loader.R")
.helpers    <- file.path(.pipe_base, "inversion-popgen-toolkit",
                           "inversion_modules", "phase_7_karyotype_groups",
                           "proposal", "lib_decompose_helpers.R")
.try_source(.reg_loader)
.try_source(.helpers)

reg <- if (exists("load_registry", mode = "function")) {
  tryCatch(load_registry(), error = function(e) {
    message("[C01m] load_registry() failed: ", conditionMessage(e),
            " — falling back to JSON sidecar writes")
    NULL
  })
} else NULL

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript STEP_C01m_distance_concordance.R <precomp_dir> <outdir> [--chrom chr] [--distances 80,160,320,640] [--scores candidate_scores.tsv.gz]")

precomp_dir <- args[1]
outdir <- args[2]
chrom_filter <- NULL
DISTANCES <- c(80L, 160L, 320L, 640L)
K_BANDS <- 3L
scores_file <- NULL  # chat-13: optional, enables per-candidate block writes

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--chrom" && i < length(args)) { chrom_filter <- args[i+1]; i <- i+2 }
  else if (a == "--distances" && i < length(args)) {
    DISTANCES <- as.integer(strsplit(args[i+1], ",")[[1]]); i <- i+2
  }
  else if (a == "--k" && i < length(args)) { K_BANDS <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--scores" && i < length(args)) { scores_file <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

message("[C01m] Distance Concordance Analysis")
message("[C01m] Distances: ", paste(DISTANCES, collapse = ", "), " windows")

# =============================================================================
# CORE: compute sample concordance at distance d
# =============================================================================

# For each pair of windows (w, w+d), assign samples to k bands using PC1.
# Use chromosome-wide fixed bands (NOT per-window k-means) for stability.
# Score: for each sample pair (i,j), what fraction of (w, w+d) pairs
# have i and j in the same band at BOTH windows?

compute_concordance_at_distance <- function(dt, pc1_cols, chr_bands, d, sample_names) {
  n_windows <- nrow(dt)
  n_samples <- length(sample_names)
  n_pairs <- n_windows - d
  if (n_pairs < 10) return(NULL)

  # For each window, get per-sample band assignment using fixed chr_bands
  # chr_bands assigns each sample to band 1/2/3 based on chromosome-wide average
  # BUT we want to know: at THIS window, is sample i still in its expected band?
  # Use per-window PC1 + nearest-center to chr_bands centers

  # Get chromosome-wide band centers from chr_bands
  pc1_mat <- as.matrix(dt[, ..pc1_cols])
  avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
  band_centers <- tapply(avg_pc1, chr_bands, mean, na.rm = TRUE)

  # Per-window band assignment: for each window, assign each sample to
  # nearest band center based on that window's PC1 values
  per_window_bands <- matrix(0L, nrow = n_windows, ncol = n_samples)
  for (wi in seq_len(n_windows)) {
    w_pc1 <- as.numeric(dt[wi, ..pc1_cols])
    for (si in seq_len(n_samples)) {
      if (!is.finite(w_pc1[si])) { per_window_bands[wi, si] <- NA; next }
      dists <- abs(w_pc1[si] - band_centers)
      per_window_bands[wi, si] <- which.min(dists)
    }
  }

  # Concordance matrix: for each pair (w, w+d), check co-grouping
  co_count <- matrix(0L, n_samples, n_samples)
  valid_count <- matrix(0L, n_samples, n_samples)

  for (wi in seq_len(n_pairs)) {
    b_w <- per_window_bands[wi, ]
    b_wd <- per_window_bands[wi + d, ]

    for (si in seq_len(n_samples)) {
      if (is.na(b_w[si]) || is.na(b_wd[si])) next
      for (sj in si:n_samples) {
        if (is.na(b_w[sj]) || is.na(b_wd[sj])) next
        valid_count[si, sj] <- valid_count[si, sj] + 1L
        valid_count[sj, si] <- valid_count[si, sj]
        # Co-grouped: same band at window w AND same band at window w+d
        if (b_w[si] == b_w[sj] && b_wd[si] == b_wd[sj]) {
          co_count[si, sj] <- co_count[si, sj] + 1L
          co_count[sj, si] <- co_count[si, sj]
        }
      }
    }

    if (wi %% 500 == 0) message("    window pair ", wi, "/", n_pairs)
  }

  # Normalize
  concordance <- co_count / pmax(valid_count, 1L)
  diag(concordance) <- 1
  rownames(concordance) <- sample_names
  colnames(concordance) <- sample_names
  concordance
}

# Faster vectorized version using band indicator matrices
compute_concordance_fast <- function(dt, pc1_cols, chr_bands, d, sample_names) {
  n_windows <- nrow(dt)
  n_samples <- length(sample_names)
  n_pairs <- n_windows - d
  if (n_pairs < 10) return(NULL)

  pc1_mat <- as.matrix(dt[, ..pc1_cols])
  avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
  band_centers <- as.numeric(tapply(avg_pc1, chr_bands, mean, na.rm = TRUE))

  # Per-window band assignment (vectorized)
  per_window_bands <- matrix(0L, nrow = n_windows, ncol = n_samples)
  for (wi in seq_len(n_windows)) {
    w_pc1 <- as.numeric(pc1_mat[wi, ])
    for (si in seq_len(n_samples)) {
      if (!is.finite(w_pc1[si])) { per_window_bands[wi, si] <- NA_integer_; next }
      per_window_bands[wi, si] <- which.min(abs(w_pc1[si] - band_centers))
    }
  }

  # Sample subset for speed (use all if <= 226, subsample if somehow more)
  use_samples <- seq_len(n_samples)

  # For each distance pair, build same-band indicator and accumulate
  concordance <- matrix(0, n_samples, n_samples)
  n_valid <- 0L

  # Process in chunks for memory
  chunk_size <- min(200L, n_pairs)
  for (chunk_start in seq(1, n_pairs, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1, n_pairs)

    for (wi in chunk_start:chunk_end) {
      b_w <- per_window_bands[wi, ]
      b_wd <- per_window_bands[wi + d, ]

      # Which samples are in the same band at window w?
      # This is an outer comparison: same_w[i,j] = (b_w[i] == b_w[j])
      valid_w <- !is.na(b_w)
      valid_wd <- !is.na(b_wd)
      valid_both <- valid_w & valid_wd

      if (sum(valid_both) < 10) next

      idx <- which(valid_both)
      bw <- b_w[idx]
      bwd <- b_wd[idx]

      # same band at w AND same band at w+d
      for (ii in seq_along(idx)) {
        for (jj in ii:length(idx)) {
          if (bw[ii] == bw[jj] && bwd[ii] == bwd[jj]) {
            concordance[idx[ii], idx[jj]] <- concordance[idx[ii], idx[jj]] + 1
            if (ii != jj) concordance[idx[jj], idx[ii]] <- concordance[idx[jj], idx[ii]] + 1
          }
        }
      }
      n_valid <- n_valid + 1L
    }
  }

  concordance <- concordance / pmax(n_valid, 1L)
  diag(concordance) <- 1
  rownames(concordance) <- sample_names
  colnames(concordance) <- sample_names
  concordance
}

# =============================================================================
# MAIN LOOP
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (!is.null(chrom_filter)) {
  rds_files <- rds_files[grepl(chrom_filter, rds_files)]
}

# chat-13: load candidate scores if provided, for per-candidate block writes
cand_dt <- NULL
if (!is.null(scores_file) && file.exists(scores_file)) {
  cand_dt <- tryCatch(fread(scores_file), error = function(e) {
    message("[C01m] could not read --scores ", scores_file, ": ",
            conditionMessage(e)); NULL
  })
  if (!is.null(cand_dt) && !"candidate_id" %in% names(cand_dt)) {
    if (all(c("chrom", "interval_id") %in% names(cand_dt))) {
      cand_dt[, candidate_id := paste0(chrom, "_", interval_id)]
    }
  }
  if (!is.null(cand_dt) && !"start_bp" %in% names(cand_dt) && "start_mb" %in% names(cand_dt)) {
    cand_dt[, start_bp := as.integer(start_mb * 1e6)]
  }
  if (!is.null(cand_dt) && !"end_bp" %in% names(cand_dt) && "end_mb" %in% names(cand_dt)) {
    cand_dt[, end_bp := as.integer(end_mb * 1e6)]
  }
  message("[C01m] loaded ", nrow(cand_dt %||% data.table()),
          " candidates for per-candidate block writes")
} else {
  message("[C01m] no --scores; skipping per-candidate block writes ",
          "(chromosome-wide outputs still produced)")
}

for (rds_file in rds_files) {
  pc <- readRDS(rds_file)
  chr <- pc$chrom
  dt <- pc$dt
  sim_mat <- pc$sim_mat

  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  if (length(pc1_cols) < 20) { message("[SKIP] ", chr); next }
  sample_names <- sub("^PC_1_", "", pc1_cols)
  n_samples <- length(sample_names)
  n_windows <- nrow(dt)

  message("\n[C01m] === ", chr, " (", n_windows, " windows, ", n_samples, " samples) ===")

  # Chromosome-wide fixed k=3 bands
  pc1_mat <- as.matrix(dt[, ..pc1_cols])
  avg_pc1 <- colMeans(pc1_mat, na.rm = TRUE)
  valid <- is.finite(avg_pc1)
  km <- tryCatch(kmeans(avg_pc1[valid], centers = K_BANDS, nstart = 10),
                  error = function(e) NULL)
  if (is.null(km)) { message("  k-means failed"); next }
  co <- order(km$centers[, 1])
  chr_bands <- integer(length(avg_pc1))
  chr_bands[valid] <- match(km$cluster, co)
  names(chr_bands) <- sample_names
  message("  Fixed bands: ", paste(table(chr_bands[valid]), collapse = "/"))

  # ── Compute concordance at each distance ──
  conc_list <- list()
  for (d in DISTANCES) {
    if (d >= n_windows / 2) {
      message("  d=", d, ": skipping (> half chromosome)")
      next
    }
    message("  Computing concordance at d=", d, " windows ...")
    t0 <- proc.time()
    conc <- compute_concordance_fast(dt, pc1_cols, chr_bands, d, sample_names)
    elapsed <- round((proc.time() - t0)[3], 1)
    if (!is.null(conc)) {
      conc_list[[as.character(d)]] <- conc
      message("    Done (", elapsed, "s)")

      # Save
      conc_dt <- as.data.table(conc)
      conc_dt[, sample := sample_names]
      fwrite(conc_dt, file.path(outdir, paste0("concordance_", chr, "_d", d, ".tsv.gz")),
             sep = "\t")
    }
  }

  if (length(conc_list) < 2) { message("  Too few distances computed"); next }

  # ── Compute contrast matrices ──
  d_names <- names(conc_list)
  contrast_list <- list()
  for (ci in seq_len(length(d_names) - 1)) {
    d_short <- d_names[ci]
    d_long <- d_names[ci + 1]
    contrast <- conc_list[[d_long]] - conc_list[[d_short]]
    contrast_list[[paste0(d_short, "_vs_", d_long)]] <- contrast

    c_dt <- as.data.table(contrast)
    c_dt[, sample := sample_names]
    fwrite(c_dt, file.path(outdir, paste0("contrast_", chr, "_d", d_short,
                                            "_vs_d", d_long, ".tsv.gz")), sep = "\t")
  }

  # ── Persistence profile: per sample pair, concordance at each distance ──
  persist_rows <- list()
  for (si in seq_len(n_samples)) {
    for (sj in (si+1):min(n_samples, si + 225)) {  # cap for memory
      if (sj > n_samples) break
      vals <- vapply(d_names, function(d) conc_list[[d]][si, sj], numeric(1))
      persist_rows[[length(persist_rows) + 1]] <- data.table(
        sample_i = sample_names[si], sample_j = sample_names[sj],
        band_i = chr_bands[si], band_j = chr_bands[sj],
        same_band = chr_bands[si] == chr_bands[sj],
        as.list(setNames(round(vals, 4), paste0("conc_d", d_names)))
      )
    }
  }
  persist_dt <- rbindlist(persist_rows, fill = TRUE)
  fwrite(persist_dt, file.path(outdir, paste0("persistence_", chr, ".tsv.gz")), sep = "\t")

  # ════════════════════════════════════════════════════════════════
  # PLOTS
  # ════════════════════════════════════════════════════════════════

  # ── D1: Multi-scale concordance heatmaps (4-panel) ──
  message("  Plotting D1: concordance heatmaps...")
  pdf(file.path(outdir, "plots", paste0(chr, "_D1_concordance_multiscale.pdf")),
      width = 16, height = 14)

  # Subsample for display
  sub_idx <- seq(1, n_samples, by = max(1, n_samples %/% 100))
  panels <- list()
  for (d_name in d_names) {
    sub_conc <- conc_list[[d_name]][sub_idx, sub_idx]
    m_dt <- data.table(
      x = rep(seq_along(sub_idx), each = length(sub_idx)),
      y = rep(seq_along(sub_idx), times = length(sub_idx)),
      value = as.vector(sub_conc),
      distance = paste0("d = ", d_name, " windows")
    )
    panels[[d_name]] <- m_dt
  }
  pdt <- rbindlist(panels)
  pdt[, distance := factor(distance, levels = unique(distance))]

  p1 <- ggplot(pdt, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("navy", "steelblue", "white", "orange", "red3"),
      values = c(0, 0.3, 0.5, 0.7, 1), name = "Concordance") +
    facet_wrap(~ distance, ncol = 2) +
    coord_fixed() +
    labs(title = paste0(chr, " — Sample Concordance at Multiple Distances"),
         subtitle = paste0("226 samples (subsampled for display) | k=", K_BANDS, " fixed bands"),
         x = "Sample", y = "Sample",
         caption = paste0("Concordance[i,j] = fraction of (window, window+d) pairs where samples i,j are in same band at both\n",
                         "Red blocks that PERSIST across distances = inversion carriers (recombination suppressed)\n",
                         "Red blocks that FADE = family LD (recombination breaks co-grouping with distance)")) +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(face = "bold"),
          plot.caption = element_text(size = 6, hjust = 0, color = "grey50"))
  print(p1)
  dev.off()

  # ── D2: Contrast matrices ──
  message("  Plotting D2: contrast matrices...")
  pdf(file.path(outdir, "plots", paste0(chr, "_D2_contrast.pdf")),
      width = 16, height = 6 * length(contrast_list))

  for (cname in names(contrast_list)) {
    sub_con <- contrast_list[[cname]][sub_idx, sub_idx]
    m_dt <- data.table(
      x = rep(seq_along(sub_idx), each = length(sub_idx)),
      y = rep(seq_along(sub_idx), times = length(sub_idx)),
      value = as.vector(sub_con)
    )

    p2 <- ggplot(m_dt, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "red3", mid = "white", high = "blue3",
                            midpoint = 0, name = "Δ concordance") +
      coord_fixed() +
      labs(title = paste0(chr, " — Concordance Contrast: d=", cname),
           subtitle = "Red = LOST co-grouping (family decay) | Blue = GAINED (unlikely, check) | White = stable",
           x = "Sample", y = "Sample",
           caption = "Contrast = concordance(longer_d) - concordance(shorter_d)\nNegative (red) = pair lost co-grouping = family LD decayed\nNear-zero (white) = pair maintained co-grouping = structural (inversion)") +
      theme_minimal(base_size = 9) +
      theme(plot.caption = element_text(size = 6, hjust = 0, color = "grey50"))
    print(p2)
  }
  dev.off()

  # ── D3: Persistence decay curves ──
  message("  Plotting D3: persistence profiles...")

  # Classify pairs: same-band vs different-band
  conc_cols <- grep("^conc_d", names(persist_dt), value = TRUE)
  decay_summary <- persist_dt[, lapply(.SD, mean, na.rm = TRUE),
                               by = same_band, .SDcols = conc_cols]
  decay_long <- melt(decay_summary, id.vars = "same_band",
                      variable.name = "distance", value.name = "mean_concordance")
  decay_long[, d_val := as.integer(sub("conc_d", "", distance))]
  decay_long[, group := ifelse(same_band, "Same chr-wide band", "Different bands")]

  p3 <- ggplot(decay_long, aes(x = d_val, y = mean_concordance, color = group)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Same chr-wide band" = "red3",
                                    "Different bands" = "steelblue")) +
    labs(title = paste0(chr, " — Concordance Decay with Distance"),
         subtitle = "Same-band pairs should decay (family LD). If they don't = inversion.",
         x = "Distance (windows)", y = "Mean concordance",
         color = NULL,
         caption = paste0("Same-band: sample pairs assigned to same chr-wide k=3 band\n",
                         "If red line stays flat = structural linkage (inversion)\n",
                         "If red line drops = family LD (normal recombination decay)")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top",
          plot.caption = element_text(size = 7, hjust = 0, color = "grey50"))

  ggsave(file.path(outdir, "plots", paste0(chr, "_D3_persistence_profile.pdf")),
         p3, width = 10, height = 7)
  ggsave(file.path(outdir, "plots", paste0(chr, "_D3_persistence_profile.png")),
         p3, width = 10, height = 7, dpi = 300)

  # ── D4: Sim_mat at multiple subsampling rates ──
  message("  Plotting D4: multi-scale sim_mat...")
  if (!is.null(sim_mat) && nrow(sim_mat) >= 100) {
    scales <- c(1, 2, 4, 8)
    pdf(file.path(outdir, "plots", paste0(chr, "_D4_simmat_multiscale.pdf")),
        width = 16, height = 14)

    sim_panels <- list()
    for (sc in scales) {
      idx <- seq(1, nrow(sim_mat), by = sc)
      idx <- idx[seq_len(min(500, length(idx)))]
      sub_sim <- sim_mat[idx, idx]
      pos_mb <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6

      m_dt <- data.table(
        x = rep(pos_mb, each = length(idx)),
        y = rep(pos_mb, times = length(idx)),
        value = as.vector(sub_sim),
        scale = paste0(sc, "× step (", length(idx), " windows)")
      )
      sim_panels[[as.character(sc)]] <- m_dt
    }
    spdt <- rbindlist(sim_panels)
    spdt[, scale := factor(scale, levels = unique(scale))]
    sim_med <- median(spdt$value, na.rm = TRUE)

    p4 <- ggplot(spdt[is.finite(value)], aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_gradientn(
        colours = c("navy", "blue3", "steelblue", "gold", "orange", "red3"),
        values = scales::rescale(c(0, sim_med*0.6, sim_med, sim_med*1.1, sim_med*1.3, 1)),
        name = "Similarity") +
      facet_wrap(~ scale, ncol = 2) +
      coord_fixed() +
      labs(title = paste0(chr, " — Similarity Matrix at Multiple Scales"),
           subtitle = "Blocks that persist across scales = real structure. Blocks that vanish = short-range LD.",
           x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
           caption = "1× = original resolution | 2× = every 2nd window | 4× = every 4th | 8× = every 8th\nInversions persist. Family blocks dissolve by 4-8×.") +
      theme_minimal(base_size = 8) +
      theme(strip.text = element_text(face = "bold", size = 9),
            plot.caption = element_text(size = 6, hjust = 0, color = "grey50"))
    print(p4)
    dev.off()
  }

  # ── D5: Persistent pair network ──
  message("  Plotting D5: persistence network...")
  # Find pairs that maintain high concordance at the longest distance
  longest_d <- tail(d_names, 1)
  longest_conc <- conc_list[[longest_d]]

  # Threshold: pairs with concordance > 0.5 at longest distance
  high_pairs <- which(longest_conc > 0.5 & upper.tri(longest_conc), arr.ind = TRUE)
  if (nrow(high_pairs) > 10) {
    net_dt <- data.table(
      sample_i = sample_names[high_pairs[, 1]],
      sample_j = sample_names[high_pairs[, 2]],
      concordance = longest_conc[high_pairs],
      band_i = chr_bands[high_pairs[, 1]],
      band_j = chr_bands[high_pairs[, 2]],
      same_band = chr_bands[high_pairs[, 1]] == chr_bands[high_pairs[, 2]]
    )

    # Per-sample count of persistent connections
    samp_counts <- c(net_dt[, .N, by = sample_i][, setNames(N, sample_i)],
                      net_dt[, .N, by = sample_j][, setNames(N, sample_j)])
    samp_tab <- tapply(samp_counts, names(samp_counts), sum)
    samp_rank <- data.table(sample = names(samp_tab),
                             n_persistent_pairs = as.integer(samp_tab),
                             band = chr_bands[names(samp_tab)])
    samp_rank <- samp_rank[order(-n_persistent_pairs)]

    p5 <- ggplot(samp_rank[seq_len(min(60, nrow(samp_rank)))],
                  aes(x = reorder(sample, n_persistent_pairs),
                      y = n_persistent_pairs, fill = factor(band))) +
      geom_col() +
      scale_fill_manual(values = c("1" = "#3b82f6", "2" = "#22c55e", "3" = "#ef4444"),
                         name = "Chr-wide band") +
      coord_flip() +
      labs(title = paste0(chr, " — Samples with Most Persistent Co-grouping (d=", longest_d, ")"),
           subtitle = paste0("Pairs maintaining concordance > 0.5 at ", longest_d, " windows apart"),
           x = NULL, y = "Number of persistent pair connections",
           caption = "Samples with many persistent connections at long distance = likely inversion carriers\nIf concentrated in one band = inversion genotype group") +
      theme_minimal(base_size = 8) +
      theme(plot.caption = element_text(size = 6, hjust = 0, color = "grey50"))

    ggsave(file.path(outdir, "plots", paste0(chr, "_D5_persistent_pairs.pdf")),
           p5, width = 10, height = 8)
    ggsave(file.path(outdir, "plots", paste0(chr, "_D5_persistent_pairs.png")),
           p5, width = 10, height = 8, dpi = 300)
  }

  # ── D6: Sample × distance decay heatmap ──
  message("  Plotting D6: sample decay heatmap...")
  # For each sample, compute its mean concordance with same-band samples at each distance
  decay_per_sample <- list()
  for (si in seq_len(n_samples)) {
    same_band_idx <- which(chr_bands == chr_bands[si] & seq_len(n_samples) != si)
    if (length(same_band_idx) < 3) next

    row <- list(sample = sample_names[si], band = chr_bands[si])
    for (d_name in d_names) {
      vals <- conc_list[[d_name]][si, same_band_idx]
      row[[paste0("d", d_name)]] <- round(mean(vals, na.rm = TRUE), 4)
    }
    decay_per_sample[[length(decay_per_sample) + 1]] <- row
  }
  decay_samp_dt <- rbindlist(decay_per_sample, fill = TRUE)

  if (nrow(decay_samp_dt) > 10) {
    # Compute decay rate: slope of concordance vs log(distance)
    d_cols <- grep("^d[0-9]", names(decay_samp_dt), value = TRUE)
    d_vals <- as.numeric(sub("^d", "", d_cols))
    decay_samp_dt[, decay_rate := {
      vals <- as.numeric(.SD)
      if (all(is.finite(vals)) && length(vals) >= 2) {
        coef(lm(vals ~ log(d_vals)))[2]
      } else NA_real_
    }, .SDcols = d_cols, by = seq_len(nrow(decay_samp_dt))]

    # Sort by decay rate (least decay = most persistent = inversion carrier?)
    decay_samp_dt <- decay_samp_dt[order(decay_rate, decreasing = TRUE)]
    decay_samp_dt[, sample_ord := seq_len(.N)]

    decay_long <- melt(decay_samp_dt, id.vars = c("sample", "band", "sample_ord", "decay_rate"),
                        measure.vars = d_cols, variable.name = "dist", value.name = "concordance")
    decay_long[, d_num := as.integer(sub("^d", "", dist))]

    p6 <- ggplot(decay_long, aes(x = factor(d_num), y = reorder(sample, -sample_ord),
                                   fill = concordance)) +
      geom_tile() +
      scale_fill_gradientn(
        colours = c("navy", "steelblue", "white", "orange", "red3"),
        values = c(0, 0.25, 0.45, 0.65, 1), name = "Concordance") +
      labs(title = paste0(chr, " — Per-Sample Concordance Decay"),
           subtitle = "Sorted by decay rate (top = least decay = possible inversion carrier)",
           x = "Distance (windows)", y = "Sample",
           caption = "Each row = one sample. Color = mean concordance with same-band samples at that distance.\nFlat rows (stays orange/red) = inversion. Fading rows (red→blue) = family LD.") +
      theme_minimal(base_size = 7) +
      theme(axis.text.y = element_text(size = 3),
            plot.caption = element_text(size = 6, hjust = 0, color = "grey50"))

    ggsave(file.path(outdir, "plots", paste0(chr, "_D6_sample_decay_heatmap.pdf")),
           p6, width = 10, height = 16)
    ggsave(file.path(outdir, "plots", paste0(chr, "_D6_sample_decay_heatmap.png")),
           p6, width = 10, height = 16, dpi = 200)
  }

  message("[C01m] ", chr, " complete: ", length(conc_list), " distances, ",
          length(contrast_list), " contrasts, 6 plot types")

  # ═══════════════════════════════════════════════════════════════════
  # chat-13 WIRING: per-candidate distance_concordance blocks
  # ═══════════════════════════════════════════════════════════════════
  # For each candidate on this chromosome, subset the per-sample-pair
  # concordance to window-pairs (w, w+d) that BOTH fall inside the
  # candidate interval, compute aggregate stats at each tested distance,
  # and emit the block per the schema at
  # registries/schemas/structured_block_schemas/distance_concordance.schema.json.
  #
  # The chromosome-wide conc_list[[d]] is a full sample×sample matrix —
  # aggregated over ALL window-pairs on the chromosome. For a per-
  # candidate block we need interval-restricted aggregates. We recompute
  # the aggregate inline using the per-window PC1 bands, which is cheap
  # relative to the chromosome-wide compute that already ran.
  if (!is.null(cand_dt) && exists("write_block_safe", mode = "function")) {
    chr_cands <- cand_dt[chrom == chr]
    if (nrow(chr_cands) > 0) {
      # Precompute per-window band assignment: re-use the chromosome-wide
      # chr_bands as the "canonical" band per sample, then for each window
      # check nearest-centre assignment (same logic as
      # compute_concordance_at_distance used internally).
      pc1_mat_full <- as.matrix(dt[, ..pc1_cols])
      # Nearest-to-chr-band-centre assignment per window per sample
      # (NA when PC1 is non-finite)
      chr_centers <- tapply(avg_pc1[valid], chr_bands[valid],
                             FUN = mean, default = NA_real_)
      chr_centers <- chr_centers[order(as.integer(names(chr_centers)))]
      window_band <- matrix(NA_integer_, nrow(pc1_mat_full),
                              ncol(pc1_mat_full))
      for (wi in seq_len(nrow(pc1_mat_full))) {
        vals_w <- pc1_mat_full[wi, ]
        finite_w <- is.finite(vals_w)
        if (any(finite_w)) {
          dists <- abs(outer(vals_w[finite_w], chr_centers, `-`))
          window_band[wi, finite_w] <- max.col(-dists, ties.method = "first")
        }
      }

      for (ci2 in seq_len(nrow(chr_cands))) {
        cdr <- chr_cands[ci2]
        cid2 <- cdr$candidate_id
        s_bp <- as.integer(cdr$start_bp)
        e_bp <- as.integer(cdr$end_bp)
        win_in <- which(dt$start_bp >= s_bp & dt$end_bp <= e_bp)
        if (length(win_in) < 20) {
          message("[C01m]   ", cid2, ": skip (only ", length(win_in),
                  " windows inside)")
          next
        }

        # Per-distance aggregates inside the candidate interval.
        per_dist <- list()
        shortest_d_conc <- NULL
        for (d in DISTANCES) {
          # keep pairs whose (w, w+d) both fall inside win_in
          pairs_in <- intersect(win_in, win_in - d)
          if (length(pairs_in) < 5) next
          # For each window-pair (w, w+d), for each sample-pair (i,j),
          # co-occurrence means window_band[w, i] == window_band[w, j]
          # AND window_band[w+d, i] == window_band[w+d, j], AND both
          # bands at w equal at w+d (persistent co-grouping).
          # Aggregate over window-pairs: fraction of pairs where (i,j)
          # co-group at both endpoints.
          nsamp <- ncol(pc1_mat_full)
          conc <- matrix(0, nsamp, nsamp)
          n_valid <- 0L
          for (w in pairs_in) {
            bw  <- window_band[w, ]
            bwd <- window_band[w + d, ]
            if (any(!is.na(bw)) && any(!is.na(bwd))) {
              # co-grouping at w AND w+d
              same_w  <- outer(bw,  bw,  `==`)
              same_wd <- outer(bwd, bwd, `==`)
              both <- same_w & same_wd
              both[is.na(both)] <- FALSE
              conc <- conc + both
              n_valid <- n_valid + 1L
            }
          }
          if (n_valid == 0L) next
          conc <- conc / n_valid
          # Aggregate stats over upper triangle (excluding diagonal)
          ut <- upper.tri(conc)
          vals <- conc[ut]
          mean_conc <- mean(vals, na.rm = TRUE)
          n_high <- sum(vals > 0.8, na.rm = TRUE)
          if (is.null(shortest_d_conc)) shortest_d_conc <- conc
          # Count decaying pairs: pairs whose concordance at this d is
          # > 0.3 LOWER than at the shortest d already computed.
          n_decay <- sum((shortest_d_conc[ut] - vals) > 0.3, na.rm = TRUE)

          per_dist[[length(per_dist) + 1L]] <- list(
            distance_windows = as.integer(d),
            mean_concordance = as.numeric(round(mean_conc, 4)),
            n_pairs_high     = as.integer(n_high),
            n_pairs_decaying = as.integer(n_decay)
          )
        }

        if (length(per_dist) == 0L) {
          message("[C01m]   ", cid2, ": no distance aggregates usable; skip")
          next
        }

        # Summary score: fraction of the shortest-d high-concordance
        # pairs that REMAIN high-concordance at the longest tested d.
        # 1 = strong inversion signal (no decay); 0 = pure family LD.
        shortest_vals <- shortest_d_conc[upper.tri(shortest_d_conc)]
        longest <- per_dist[[length(per_dist)]]
        longest_mean <- longest$mean_concordance
        shortest_mean <- per_dist[[1]]$mean_concordance
        inv_vs_fam <- if (is.finite(shortest_mean) && shortest_mean > 0)
                        pmax(0, pmin(1, longest_mean / shortest_mean))
                      else NA_real_

        # Persistent carriers: samples whose mean co-grouping at the
        # longest distance with their own band-mates stays high (> 0.6).
        persistent <- character(0)
        if (!is.null(shortest_d_conc)) {
          # Use the longest-d matrix for persistence: we need to
          # recompute one more time, but for simplicity use the
          # chromosome-wide conc_list[[longest]] if available.
          longest_d_name <- as.character(longest$distance_windows)
          if (!is.null(conc_list[[longest_d_name]])) {
            long_mat <- conc_list[[longest_d_name]]
            # sample-level persistence: max off-diagonal value
            for (si in seq_len(nrow(long_mat))) {
              off <- long_mat[si, -si]
              if (any(is.finite(off)) && max(off, na.rm = TRUE) > 0.6) {
                persistent <- c(persistent, sample_names[si])
              }
            }
          }
        }

        # decaying pair count = use last d's n_decay
        total_decay <- per_dist[[length(per_dist)]]$n_pairs_decaying

        dc_block <- list(
          candidate_id      = as.character(cid2),
          chrom             = chr,
          distances_tested  = as.integer(vapply(per_dist,
                                                  function(x) x$distance_windows,
                                                  integer(1))),
          per_distance_summary      = per_dist,
          persistent_carriers        = as.character(persistent),
          decaying_pairs_count       = as.integer(total_decay),
          inversion_vs_family_score  = as.numeric(round(inv_vs_fam %||% NA_real_, 4))
        )

        cand_outdir_m <- file.path(outdir, "blocks", as.character(cid2))
        dir.create(cand_outdir_m, recursive = TRUE, showWarnings = FALSE)
        write_block_safe(
          reg             = reg,
          candidate_id    = as.character(cid2),
          block_type      = "distance_concordance",
          data            = dc_block,
          source_script   = "STEP_C01m_distance_concordance.R",
          outdir_fallback = cand_outdir_m
        )
      }
    }
  }
}

message("\n[DONE] Distance concordance → ", outdir)
