#!/usr/bin/env Rscript
# =============================================================================
# plot_sample_regime_dag.R — per-sample DAG diagnostic for multi_recomb
# =============================================================================
# Tiny diagnostic plot: for each sample in a candidate, render the run-length-
# encoded regime track as a horizontal strip of rectangles, width proportional
# to run length in windows, colour-coded by regime label. R_fired samples get
# a coloured border; dominant regime is the background colour.
#
# Consumes the regime_sample_dag block written by STEP_C01i_b_multi_recomb.R
# (per_sample with dag_nodes_compact + dag_node_weights, or can derive from
# the raw regime_memberships.tsv.gz). Prefer consuming the block since the
# weights-list may already have been dropped for JSON compatibility.
#
# Base R graphics only — no ggplot2, no ggraph. Handoff preference:
# "prefer plain R — fewer dependencies".
#
# USAGE:
#   Rscript plot_sample_regime_dag.R \
#     --regime_memb  path/to/regime_memberships.tsv.gz \
#     --candidates   path/to/candidates.tsv \
#     --outdir       plots/dag/ \
#     [--min_dev_frac 0.15] [--min_dev_bp 50000]
#
# Produces one PDF per candidate: <outdir>/<cid>_dag.pdf
# Grid layout: up to 6 columns × N rows. Samples sorted by deviation_fraction
# descending (most-deviating at top).
#
# CAN ALSO BE CALLED AS A LIBRARY:
#   source("plot_sample_regime_dag.R")
#   plot_candidate_dag(per_sample_dt, cid = "LG12_17", out_pdf = "foo.pdf")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ── Palette: deterministic colour per regime label ──────────────────────────
# ColorBrewer Set1-ish, extended. Stable across plots because we sort the
# label set and map by index.
.regime_palette <- function(labels) {
  pal <- c("#4477AA", "#EE6677", "#228833", "#CCBB44",
           "#66CCEE", "#AA3377", "#BBBBBB", "#F28E2B",
           "#59A14F", "#E15759", "#AF7AA1", "#9C755F")
  u  <- sort(unique(labels))
  # Extend palette by recycling if more regimes than colours
  cols <- pal[((seq_along(u) - 1L) %% length(pal)) + 1L]
  setNames(cols, u)
}

# ── One sample's tiny DAG strip ─────────────────────────────────────────────
# Draw into the current plot device at y == y_center, spanning x in [0, 1].
# run_labels, run_weights: same length, describe the rle.
# dominant_label:          string to highlight as "main colour"
# R_fired:                 TRUE → draw a yellow box around the whole strip
# sample_id:               left-margin label
# deviation_fraction:      right-margin tiny number
.draw_one_dag <- function(run_labels, run_weights,
                            dominant_label, R_fired,
                            sample_id, deviation_fraction,
                            y_center, palette,
                            strip_height = 0.8,
                            x_left = 0, x_right = 1) {
  total_w <- sum(run_weights)
  if (total_w <= 0L) return(invisible())
  widths <- (x_right - x_left) * run_weights / total_w
  x_starts <- x_left + cumsum(c(0, widths[-length(widths)]))
  y_bot <- y_center - strip_height / 2
  y_top <- y_center + strip_height / 2

  # Background box (for R_fired highlight) — drawn first so rectangles sit on top
  if (isTRUE(R_fired)) {
    rect(x_left - 0.01, y_bot - 0.08, x_right + 0.01, y_top + 0.08,
         col = "#FFF3B0", border = NA)
  }

  # Per-run rectangles
  for (i in seq_along(run_labels)) {
    col <- palette[[run_labels[i]]]
    if (is.null(col) || is.na(col)) col <- "#DDDDDD"
    # Dominant runs: normal fill; non-dominant runs: hatched with same fill
    # but a dark border to emphasise they are deviations.
    is_dom <- run_labels[i] == dominant_label
    rect(x_starts[i], y_bot, x_starts[i] + widths[i], y_top,
         col = col,
         border = if (is_dom) col else "black",
         lwd    = if (is_dom) 0 else 1)
    # Label runs with >= 2 windows
    if (run_weights[i] >= 2L && widths[i] > 0.03) {
      text(x_starts[i] + widths[i] / 2, y_center,
           run_labels[i], cex = 0.55,
           col = if (is_dom) "white" else "black")
    }
  }

  # Sample label
  text(x_left - 0.02, y_center, sample_id, cex = 0.55, adj = c(1, 0.5))

  # Deviation fraction on the right
  text(x_right + 0.02, y_center,
       sprintf("%.2f", deviation_fraction),
       cex = 0.5, adj = c(0, 0.5),
       col = if (isTRUE(R_fired)) "#CC3300" else "#666666")
}

# ── Main public entry: one candidate → one PDF page with a grid of strips ──
# per_sample_dt expected columns:
#   sample_id, dag_nodes_compact (A->B->A), dag_node_weights (list<int>),
#   dominant_regime, deviation_fraction, R_fired
# If dag_node_weights is absent (it was dropped for JSON round-trip in the
# registry-block version), pass the raw regime_memb so it can be reconstructed.
plot_candidate_dag <- function(per_sample_dt,
                                  cid,
                                  out_pdf,
                                  regime_memb = NULL,
                                  interval    = NULL,
                                  cols_per_row = 2L,
                                  strip_height = 0.7) {

  if (!nrow(per_sample_dt)) {
    message("[dag_plot] ", cid, ": no samples, skip"); return(invisible())
  }

  # Reconstruct weights from regime_memb if missing
  has_weights <- "dag_node_weights" %in% names(per_sample_dt) &&
                 is.list(per_sample_dt$dag_node_weights)
  if (!has_weights) {
    if (is.null(regime_memb) || is.null(interval)) {
      stop("[dag_plot] per_sample_dt lacks dag_node_weights and ",
           "regime_memb/interval not supplied for reconstruction")
    }
    rm <- copy(regime_memb)
    if (!"sample" %in% names(rm) && "sample_id" %in% names(rm))
      setnames(rm, "sample_id", "sample")
    pos_col <- if ("pos_mid_mb" %in% names(rm)) "pos_mid_mb"
               else if ("pos_mb" %in% names(rm)) "pos_mb"
               else stop("[dag_plot] regime_memb has no pos column")
    s_mb <- interval$start_mb %||% (interval$start_bp / 1e6)
    e_mb <- interval$end_mb   %||% (interval$end_bp   / 1e6)
    rm_in <- rm[get(pos_col) >= s_mb & get(pos_col) <= e_mb]
    setorderv(rm_in, c("sample", pos_col))

    weights_by_sample <- rm_in[, {
      r <- rle(as.character(group_id))
      list(labels = list(r$values), weights = list(as.integer(r$lengths)))
    }, by = sample]
    setnames(weights_by_sample, "sample", "sample_id")

    per_sample_dt <- merge(per_sample_dt, weights_by_sample,
                           by = "sample_id", all.x = TRUE, sort = FALSE)
    per_sample_dt[, dag_node_weights := weights]
    # also reconstruct labels from dag_nodes_compact or labels list
    if ("labels" %in% names(per_sample_dt)) {
      per_sample_dt[, dag_labels := labels]
    }
  }

  # Extract per-sample labels vector from dag_nodes_compact
  if (!"dag_labels" %in% names(per_sample_dt)) {
    per_sample_dt[, dag_labels := lapply(dag_nodes_compact, function(s) {
      if (is.na(s) || !nzchar(s)) character(0)
      else strsplit(s, "->", fixed = TRUE)[[1]]
    })]
  }

  # Colour palette across all regimes seen
  all_labels <- unique(unlist(per_sample_dt$dag_labels))
  all_labels <- all_labels[!is.na(all_labels) & nzchar(all_labels)]
  palette <- .regime_palette(all_labels)

  # Sort samples: R_fired first (red-ringed ones at top), then by
  # deviation_fraction descending
  setorder(per_sample_dt, -R_fired, -deviation_fraction)

  # Grid geometry
  n <- nrow(per_sample_dt)
  ncol <- cols_per_row
  nrow <- ceiling(n / ncol)
  strip_pitch <- 1.2  # vertical distance between strip centres

  pdf_w <- 3 + 4 * ncol
  pdf_h <- 1 + strip_pitch * nrow / 2
  pdf(out_pdf, width = pdf_w, height = pdf_h)
  on.exit(dev.off(), add = TRUE)

  par(mar = c(2, 0.5, 2.5, 0.5), xpd = NA)
  plot.new()
  plot.window(xlim = c(-0.15, 1.15 + (ncol - 1L) * 1.20),
              ylim = c(-(nrow * strip_pitch + 0.5), 1))
  title(sprintf("Per-sample regime DAG — candidate %s  (n=%d)", cid, n),
        cex.main = 0.9)

  # Legend at the top
  lx <- 0; ly <- 0.3
  text(lx - 0.02, ly, "regimes:", cex = 0.6, adj = c(1, 0.5))
  for (i in seq_along(palette)) {
    rect(lx + (i - 1L) * 0.07, ly - 0.15,
         lx + (i - 1L) * 0.07 + 0.05, ly + 0.15,
         col = palette[i], border = NA)
    text(lx + (i - 1L) * 0.07 + 0.025, ly - 0.32,
         names(palette)[i], cex = 0.5)
  }
  text(0.55, ly, "(yellow background = R_fired; dev_frac on right)",
       cex = 0.5, adj = c(0, 0.5), col = "#666666")

  # Draw each sample
  for (k in seq_len(n)) {
    row_i <- ((k - 1L) %/% ncol)
    col_i <- ((k - 1L) %%  ncol)
    x_off <- col_i * 1.20
    y_c   <- -(row_i + 0.5) * strip_pitch

    lbls <- per_sample_dt$dag_labels[[k]]
    wts  <- per_sample_dt$dag_node_weights[[k]]
    if (is.null(lbls) || is.null(wts) || length(lbls) == 0L || length(wts) == 0L) {
      text(x_off + 0.5, y_c, "(no track)", cex = 0.5, col = "#CC0000")
      next
    }

    .draw_one_dag(
      run_labels         = lbls,
      run_weights        = wts,
      dominant_label     = per_sample_dt$dominant_regime[k],
      R_fired            = per_sample_dt$R_fired[k],
      sample_id          = per_sample_dt$sample_id[k],
      deviation_fraction = per_sample_dt$deviation_fraction[k],
      y_center           = y_c,
      palette            = palette,
      strip_height       = strip_height,
      x_left             = x_off,
      x_right            = x_off + 1
    )
  }

  invisible(NULL)
}

# ─────────────────────────────────────────────────────────────────────────────
# CLI driver — only runs when invoked directly via Rscript, NOT when sourced.
# Detection: commandArgs(FALSE) contains --file=<this_script> only in the
# Rscript entry-point case; when sourced from another script the --file
# argument points at the caller, not this file.
# ─────────────────────────────────────────────────────────────────────────────
.this_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile %||% "", mustWork = FALSE),
  error = function(e) ""
)
.is_direct_run <- {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg  <- sub("^--file=", "", args_full[grep("^--file=", args_full)])
  if (!length(file_arg)) FALSE
  else identical(normalizePath(file_arg, mustWork = FALSE),
                 normalizePath(.this_file,  mustWork = FALSE)) ||
       identical(basename(file_arg), "plot_sample_regime_dag.R")
}

if (isTRUE(.is_direct_run)) {

  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option("--regime_memb", type = "character", help = "C01j regime_memberships.tsv.gz"),
    make_option("--candidates",  type = "character", help = "Candidate table TSV"),
    make_option("--outdir",      type = "character", default = "dag_plots/"),
    make_option("--min_dev_frac",type = "double",    default = 0.15),
    make_option("--min_dev_bp",  type = "integer",   default = 50000L),
    make_option("--tier_max",    type = "integer",   default = 3L),
    make_option("--cols",        type = "integer",   default = 2L)
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$regime_memb) || !file.exists(opt$regime_memb))
    stop("--regime_memb required and must exist")
  if (is.null(opt$candidates) || !file.exists(opt$candidates))
    stop("--candidates required and must exist")

  # Bring in derive_R_from_regime from sibling lib
  combo_paths <- c("lib_recomb_combination.R",
                   file.path(dirname(sys.frame(1)$ofile %||% "."),
                              "lib_recomb_combination.R"),
                   "../4b_group_proposal/lib_recomb_combination.R")
  loaded <- FALSE
  for (p in combo_paths) if (file.exists(p)) { source(p); loaded <- TRUE; break }
  if (!loaded) stop("[dag_plot] cannot find lib_recomb_combination.R")

  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  rm_all <- fread(opt$regime_memb)
  if (!"sample" %in% names(rm_all) && "sample_id" %in% names(rm_all))
    setnames(rm_all, "sample_id", "sample")

  cand <- fread(opt$candidates)
  if (!"candidate_id" %in% names(cand))
    cand[, candidate_id := paste0(chrom, "_", interval_id)]
  if (!"start_bp" %in% names(cand)) cand[, start_bp := as.integer(start_mb * 1e6)]
  if (!"end_bp"   %in% names(cand)) cand[, end_bp   := as.integer(end_mb   * 1e6)]
  if (!"tier" %in% names(cand)) cand[, tier := 2L]
  cand <- cand[tier <= opt$tier_max]

  for (ci in seq_len(nrow(cand))) {
    cid <- cand$candidate_id[ci]
    chr <- cand$chrom[ci]
    s   <- as.integer(cand$start_bp[ci])
    e   <- as.integer(cand$end_bp[ci])
    rm_cid <- rm_all[chrom == chr]
    if (!nrow(rm_cid)) next
    r <- derive_R_from_regime(rm_cid, list(start_bp = s, end_bp = e),
                              min_deviation_fraction = opt$min_dev_frac,
                              min_deviation_bp       = opt$min_dev_bp)
    if (!nrow(r$per_sample)) next
    out_pdf <- file.path(opt$outdir, paste0(cid, "_dag.pdf"))
    cat("[dag_plot] ", cid, " → ", out_pdf, "  (",
        nrow(r$per_sample), " samples)\n", sep = "")
    plot_candidate_dag(r$per_sample, cid = cid, out_pdf = out_pdf,
                       regime_memb = rm_cid,
                       interval    = list(start_bp = s, end_bp = e),
                       cols_per_row= opt$cols)
  }
}
