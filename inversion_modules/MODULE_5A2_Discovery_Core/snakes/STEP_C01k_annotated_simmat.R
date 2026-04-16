#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01k_annotated_simmat.R  (v8.4)
#
# ANNOTATED SIMILARITY MATRIX  --  the final integrative figure.
#
# Takes the raw sim_mat heatmap and overlays everything the pipeline found:
#   - Triangle outlines (from C01c) with type labels
#   - System span bars (from C01i co-segregation blocks)
#   - Regime state strip (from C01j compatibility engine)
#   - Tier + verdict annotations (from C01d + C01f)
#   - Family vs inversion labels (from C01f hypothesis tests)
#   - Recombinant zone markers (from C01h)
#   - Bridge arcs between connected triangles (from C01c bridges)
#
# One figure per chromosome that tells the complete story.
#
# Inputs:
#   --precomp <dir>             -- sim_mat + window coordinates
#   --triangles <dir>           -- intervals, bridges, sub-regimes
#   --scores <file>             -- candidate scoring with tiers
#   [--hyp_dir <dir>]           -- hypothesis verdicts
#   [--coseg_dir <dir>]         -- multi-inversion systems
#   [--regime_dir <dir>]        -- regime states
#   [--recomb_dir <dir>]        -- recombinant calls
#   --outdir <dir>
#   [--chrom <chr>]             -- single chromosome or "all"
#
# Output:
#   <chr>_annotated_simmat.png  -- per-chromosome annotated heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
precomp_dir <- NULL; triangle_dir <- NULL; scores_file <- NULL
hyp_dir <- NULL; coseg_dir <- NULL; regime_dir <- NULL; recomb_dir <- NULL
outdir <- "annotated_simmat"; chrom_filter <- NULL

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--precomp" && i < length(args))     { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--triangles" && i < length(args))  { triangle_dir <- args[i+1]; i <- i+2 }
  else if (a == "--scores" && i < length(args))     { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--hyp_dir" && i < length(args))    { hyp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--coseg_dir" && i < length(args))  { coseg_dir <- args[i+1]; i <- i+2 }
  else if (a == "--regime_dir" && i < length(args)) { regime_dir <- args[i+1]; i <- i+2 }
  else if (a == "--recomb_dir" && i < length(args)) { recomb_dir <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))     { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--chrom" && i < length(args))      { chrom_filter <- args[i+1]; i <- i+2 }
  else { i <- i+1 }
}

if (is.null(precomp_dir)) stop("--precomp required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DPI <- 400

# =============================================================================
# LOAD ALL UPSTREAM DATA
# =============================================================================

message("[C01k] Loading upstream data...")

# Triangle intervals
iv_dt <- data.table()
if (!is.null(triangle_dir)) {
  f <- file.path(triangle_dir, "triangle_intervals.tsv.gz")
  if (file.exists(f)) iv_dt <- fread(f)
}

# Bridges
bridge_dt <- data.table()
if (!is.null(triangle_dir)) {
  f <- file.path(triangle_dir, "triangle_bridges.tsv.gz")
  if (file.exists(f)) bridge_dt <- fread(f)
}

# Sub-regimes
subreg_dt <- data.table()
if (!is.null(triangle_dir)) {
  f <- file.path(triangle_dir, "triangle_subregimes.tsv.gz")
  if (file.exists(f)) subreg_dt <- fread(f)
}

# Candidate scores
cand_dt <- data.table()
if (!is.null(scores_file) && file.exists(scores_file)) cand_dt <- fread(scores_file)

# Hypothesis verdicts
hyp_dt <- data.table()
if (!is.null(hyp_dir)) {
  f <- file.path(hyp_dir, "hypothesis_verdicts.tsv")
  if (file.exists(f)) hyp_dt <- fread(f)
}

# Co-segregation systems
sys_dt <- data.table()
if (!is.null(coseg_dir)) {
  f <- file.path(coseg_dir, "marker_coseg_blocks.tsv.gz")
  if (file.exists(f)) sys_dt <- fread(f)
}

# Regime states
regime_dt <- data.table()
if (!is.null(regime_dir)) {
  f <- file.path(regime_dir, "regime_state_labels.tsv.gz")
  if (file.exists(f)) regime_dt <- fread(f)
}

# Recombinant calls
recomb_dt <- data.table()
if (!is.null(recomb_dir)) {
  f <- file.path(recomb_dir, "recombinant_summary.tsv")
  if (file.exists(f)) recomb_dt <- fread(f)
}

# =============================================================================
# PER-CHROMOSOME ANNOTATED FIGURE
# =============================================================================

rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))

for (f in rds_files) {
  pc <- readRDS(f)
  chr <- pc$chrom
  if (!is.null(chrom_filter) && chrom_filter != "all" && chr != chrom_filter) next

  dt <- pc$dt; sim_mat <- pc$sim_mat; n_w <- nrow(sim_mat)
  if (n_w < 50) next
  message("\n[C01k] === ", chr, " (", n_w, " windows) ===")

  # Subsample sim_mat for plotting
  step <- max(1L, n_w %/% 500)
  idx <- seq(1, n_w, by = step)
  sub_sim <- sim_mat[idx, idx]
  ns <- length(idx)

  # Window positions in Mb
  win_mb <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6

  # Build heatmap data
  hm_rows <- list()
  for (i in seq_len(ns)) {
    for (j in i:ns) {
      hm_rows[[length(hm_rows) + 1]] <- data.table(
        x = win_mb[i], y = win_mb[j], sim = sub_sim[i, j])
    }
  }
  hm_dt <- rbindlist(hm_rows)

  # --- Base heatmap ---
  p <- ggplot(hm_dt, aes(x = x, y = y, fill = sim)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("navy", "blue3", "purple3", "orange", "gold", "lightyellow"),
      name = "Similarity") +
    coord_fixed()

  # --- Overlay 1: Triangle outlines with type coloring ---
  chr_iv <- iv_dt[chrom == chr]
  if (nrow(chr_iv) > 0) {
    type_colors <- c("strong_triangle" = "red", "moderate_triangle" = "darkorange",
                      "patchy_signal" = "gold", "diffuse_zone" = "grey70",
                      "weak_zone" = "grey50", "transition" = "grey30")

    for (ri in seq_len(nrow(chr_iv))) {
      iv <- chr_iv[ri]
      p <- p +
        annotate("rect", xmin = iv$start_mb, xmax = iv$end_mb,
                 ymin = iv$start_mb, ymax = iv$end_mb,
                 fill = NA,
                 color = type_colors[iv$interval_type] %||% "white",
                 linewidth = 0.8, linetype = "solid")
    }
  }

  # --- Overlay 2: Tier + verdict labels inside triangles ---
  chr_cand <- cand_dt[chrom == chr]
  if (nrow(chr_cand) > 0) {
    tier_colors <- c("1" = "red3", "2" = "darkorange", "3" = "steelblue", "4" = "grey60")

    for (ri in seq_len(nrow(chr_cand))) {
      cd <- chr_cand[ri]
      mid_mb <- (cd$start_mb + cd$end_mb) / 2

      # Get verdict if available
      verdict_label <- ""
      if (nrow(hyp_dt) > 0) {
        hv <- hyp_dt[chrom == chr & interval_id == cd$interval_id]
        if (nrow(hv) > 0) verdict_label <- paste0("\n", hv$verdict[1])
      }

      label <- paste0("T", cd$tier, " ", cd$pattern %||% "", verdict_label)

      p <- p + annotate("text", x = mid_mb, y = cd$start_mb - (cd$end_mb - cd$start_mb) * 0.1,
                          label = label, size = 1.8, fontface = "bold",
                          color = tier_colors[as.character(cd$tier)] %||% "grey40")
    }
  }

  # --- Overlay 3: System span bars (from C01i) along bottom ---
  chr_sys <- sys_dt[chrom == chr]
  if (nrow(chr_sys) > 0) {
    sys_colors <- c("steelblue", "darkorange", "green4", "red3",
                     "gold3", "purple3", "brown", "cyan4")
    y_base <- min(win_mb) - (max(win_mb) - min(win_mb)) * 0.02

    for (si in seq_len(nrow(chr_sys))) {
      s <- chr_sys[si]
      col <- sys_colors[pmin(s$system_id, length(sys_colors))]
      s_start <- s$span_start / 1e6; s_end <- s$span_end / 1e6
      y_pos <- y_base - si * (max(win_mb) - min(win_mb)) * 0.015

      p <- p + annotate("segment", x = s_start, xend = s_end,
                          y = y_pos, yend = y_pos,
                          color = col, linewidth = 3, alpha = 0.8) +
        annotate("text", x = (s_start + s_end) / 2, y = y_pos,
                  label = paste0("S", s$system_id, " ", round(s$frequency * 100), "%"),
                  size = 1.5, color = "white", fontface = "bold")
    }
  }

  # --- Overlay 4: Bridge arcs between connected intervals ---
  chr_bridges <- bridge_dt[chrom == chr & bridge_type != "no_bridge"]
  if (nrow(chr_bridges) > 0 && nrow(chr_iv) > 0) {
    bridge_colors <- c("strong_bridge" = "red", "moderate_bridge" = "orange",
                        "weak_bridge" = "gold")
    for (bi in seq_len(nrow(chr_bridges))) {
      br <- chr_bridges[bi]
      iv_a <- chr_iv[interval_id == br$id_a]
      iv_b <- chr_iv[interval_id == br$id_b]
      if (nrow(iv_a) == 0 || nrow(iv_b) == 0) next

      mid_a <- (iv_a$start_mb + iv_a$end_mb) / 2
      mid_b <- (iv_b$start_mb + iv_b$end_mb) / 2

      # Draw off-diagonal marker
      p <- p + annotate("rect",
                          xmin = iv_a$start_mb, xmax = iv_a$end_mb,
                          ymin = iv_b$start_mb, ymax = iv_b$end_mb,
                          fill = NA,
                          color = bridge_colors[br$bridge_type] %||% "gold",
                          linewidth = 0.5, linetype = "dashed")
    }
  }

  # --- Overlay 5: Regime state strip along the diagonal ---
  chr_reg <- regime_dt[chrom == chr]
  if (nrow(chr_reg) > 0) {
    state_colors <- c("clean_inversion" = "red3", "structured_moderate" = "darkorange",
                       "structured_complex" = "purple3", "weak_signal" = "gold3",
                       "background_soup" = "grey70", "transition" = "grey90")

    for (ri in seq_len(nrow(chr_reg))) {
      rg <- chr_reg[ri]
      p <- p + annotate("point", x = rg$pos_mid_mb, y = rg$pos_mid_mb,
                          color = state_colors[rg$state] %||% "grey50",
                          size = 0.3)
    }
  }

  # --- Overlay 6: Recombinant zone markers ---
  chr_recomb <- recomb_dt[chrom == chr]
  if (nrow(chr_recomb) > 0) {
    for (ri in seq_len(nrow(chr_recomb))) {
      rc <- chr_recomb[ri]
      if (rc$n_simple_recomb > 0 || rc$n_double_xo > 0) {
        p <- p + annotate("text",
                            x = rc$start_mb, y = rc$end_mb,
                            label = paste0("R:", rc$n_simple_recomb,
                                          " DXO:", rc$n_double_xo),
                            size = 1.5, color = "red3", fontface = "italic")
      }
    }
  }

  # --- Final labels ---
  p <- p + labs(
    title = paste0(chr, " -- Annotated Similarity Matrix"),
    subtitle = paste0(n_w, " windows | ",
                      nrow(chr_iv), " triangles | ",
                      nrow(chr_cand), " candidates | ",
                      nrow(chr_sys), " systems"),
    x = paste0(chr, " (Mb)"), y = paste0(chr, " (Mb)"),
    caption = paste0(
      "Squares: triangle outlines (red=strong, orange=moderate, grey=weak)\n",
      "Labels: Tier + pattern + hypothesis verdict\n",
      "Bottom bars: co-segregation systems (colored, with frequency)\n",
      "Dashed off-diagonal: bridge connections between intervals\n",
      "Diagonal dots: regime state (red=clean inversion, purple=complex, grey=soup)\n",
      "R/DXO: recombinant / double-crossover counts"
    )) +
    theme_minimal(base_size = 8) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          plot.caption = element_text(size = 5, color = "grey60", hjust = 0),
          legend.key.size = unit(0.3, "cm"))

  # Save
  fig_height <- max(10, 10 + nrow(chr_sys) * 0.3)
  tryCatch(
    ggsave(file.path(outdir, paste0(chr, "_annotated_simmat.png")),
           p, width = 12, height = fig_height, dpi = DPI),
    error = function(e) message("  [PLOT] ", e$message))

  message("  -> ", file.path(outdir, paste0(chr, "_annotated_simmat.png")))
}

message("\n[DONE] -> ", outdir)
