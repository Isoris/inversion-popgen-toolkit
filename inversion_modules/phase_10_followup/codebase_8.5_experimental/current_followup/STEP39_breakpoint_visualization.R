#!/usr/bin/env Rscript

# =============================================================================
# STEP39_breakpoint_visualization.R  (v1.0)
#
# BREAKPOINT SUPPORT VISUALIZATION — Two plot families:
#
#   FAMILY A — Population-level chromosome overview
#     All breakpoint clusters on a chromosome: density ridges, median lines,
#     paired arc/ribbon links, width-encoded by sample count and read support.
#
#   FAMILY B — Focal candidate highlighted view
#     Same as Family A but with population background greyed out (α=0.3)
#     and one focal candidate remaining colored and emphasized.
#
# VISUAL ENCODING RULES:
#   Ribbon width     = number of unique samples supporting the paired cluster
#   Line thickness   = total read/support evidence (SR + PE)
#   Alpha            = support confidence (more reads = more opaque)
#   Color (pop view) = svtype_mode (INV=blue, BND=orange, COMBINED=purple)
#   Color (focal)    = focal candidate colored, everything else grey
#
# INPUTS:
#   <config.R>  — standard followup config
#   [--sv_dir <path>]   — directory with genome-wide DELLY/Manta parsed TSVs
#   [--cid <int>]       — focal candidate for highlighted view (optional)
#
# OUTPUTS:
#   breakpoint_paired_clusters.tsv.gz       — genome-wide paired clusters
#   breakpoint_density_plot_table.tsv.gz    — per-call density-ready table
#   breakpoint_link_plot_table.tsv.gz       — paired link/ribbon table
#   chromosome_breakpoint_overview.pdf      — Family A population plots
#   candidate_<ID>_breakpoint_highlighted.pdf — Family B per candidate
#
# Usage:
#   Rscript STEP39_breakpoint_visualization.R <config.R> \
#     --sv_dir <path> [--cid all] [--boundary_bp 50000] [--cluster_bp 5000]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"

sv_dir       <- NULL
cid_filter   <- NA_integer_
BOUNDARY_BP  <- 50000L
CLUSTER_BP   <- 5000L
PAIR_TOL_BP  <- 20000L  # tolerance for pairing left+right breakpoints
MIN_SAMPLES  <- 2L

i <- 2L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--sv_dir" && i < length(args)) {
    sv_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--cid" && i < length(args)) {
    v <- args[i + 1]
    if (v != "all") cid_filter <- as.integer(v)
    i <- i + 2L
  } else if (a == "--boundary_bp" && i < length(args)) {
    BOUNDARY_BP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--cluster_bp" && i < length(args)) {
    CLUSTER_BP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--pair_tol_bp" && i < length(args)) {
    PAIR_TOL_BP <- as.integer(args[i + 1]); i <- i + 2L
  } else {
    i <- i + 1L
  }
}

source(config_file)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand_focal <- cand[candidate_id == cid_filter] else cand_focal <- cand

VIZ_DIR <- file.path(INV_ROOT, "22_breakpoint_visualization")
ensure_dir(VIZ_DIR)

message("[STEP39] Breakpoint visualization module")
message("[STEP39] Candidates: ", nrow(cand), " (focal: ", nrow(cand_focal), ")")

# =============================================================================
# STEP 1 — LOAD ALL SV CALLS GENOME-WIDE
# =============================================================================

load_all_sv_calls <- function(followup_dir, chroms) {
  all_calls <- list()
  for (chr in chroms) {
    for (cid in cand$candidate_id) {
      cdir_new <- file.path(followup_dir, paste0(chr, "_cand", cid, "_",
                           sprintf("%.2fMb", as.numeric(cand[candidate_id == cid]$start_bp) / 1e6)))
      cdir_old <- file.path(followup_dir, paste0(chr, ".candidate_", cid))
      cdir <- if (dir.exists(cdir_new)) cdir_new else cdir_old
      f <- file.path(cdir, "candidate_bp_sv_calls.tsv.gz")
      if (file.exists(f)) {
        dt <- tryCatch(fread(f), error = function(e) NULL)
        if (!is.null(dt) && nrow(dt) > 0) {
          dt[, source_candidate_id := cid]
          all_calls[[length(all_calls) + 1]] <- dt
        }
      }
    }
  }
  if (length(all_calls) == 0) return(data.table())
  rbindlist(all_calls, fill = TRUE)
}

all_sv <- load_all_sv_calls(FOLLOWUP_DIR, unique(cand$chrom))
message("[STEP39] Total SV calls loaded: ", nrow(all_sv))

if (nrow(all_sv) == 0) {
  message("[STEP39] No SV calls found — generating empty outputs")
  # Write empty files and exit
  fwrite(data.table(), file.path(VIZ_DIR, "breakpoint_paired_clusters.tsv.gz"), sep = "\t")
  fwrite(data.table(), file.path(VIZ_DIR, "breakpoint_density_plot_table.tsv.gz"), sep = "\t")
  fwrite(data.table(), file.path(VIZ_DIR, "breakpoint_link_plot_table.tsv.gz"), sep = "\t")
  message("[DONE] STEP39 — no data to plot")
  quit(save = "no", status = 0)
}

# Ensure CT field (connection type) exists
if (!"ct" %in% names(all_sv)) {
  # Infer CT from svtype: INV calls are 5to5 or 3to3, BND is variable
  all_sv[, ct := fifelse(toupper(svtype) == "INV", "5to5", "unknown")]
}

# =============================================================================
# STEP 2 — GENOME-WIDE INTERVAL CLUSTERING
# =============================================================================

# For INV calls: pos1 = left breakpoint, pos2 = right breakpoint
# For BND calls: pos1 = this position, pos2 = mate position (if available)
# Cluster both ends independently, then pair them

cluster_positions <- function(positions, tol_bp) {
  if (length(positions) == 0) return(integer(0))
  ord <- order(positions)
  pos_sorted <- positions[ord]
  cids <- integer(length(positions))
  cid <- 1L
  cids[ord[1]] <- cid
  for (k in seq_along(ord)[-1]) {
    if (pos_sorted[k] - pos_sorted[k - 1] > tol_bp) cid <- cid + 1L
    cids[ord[k]] <- cid
  }
  cids
}

# Build paired interval clusters
build_paired_clusters <- function(sv_dt) {
  # Only work with calls that have both pos1 and pos2
  paired <- sv_dt[!is.na(pos1) & !is.na(pos2) & pos1 < pos2]
  if (nrow(paired) == 0) return(data.table())

  result_list <- list()
  pcid <- 0L

  for (chr in unique(paired$chrom1)) {
    chr_dt <- paired[chrom1 == chr]

    # Cluster by first breakpoint
    chr_dt[, bp1_cluster := cluster_positions(pos1, CLUSTER_BP)]
    # Within each bp1 cluster, cluster by second breakpoint
    for (b1c in unique(chr_dt$bp1_cluster)) {
      sub <- chr_dt[bp1_cluster == b1c]
      sub[, bp2_cluster := cluster_positions(pos2, CLUSTER_BP)]

      for (b2c in unique(sub$bp2_cluster)) {
        pair_sub <- sub[bp2_cluster == b2c]
        n_samp <- uniqueN(pair_sub$sample)
        if (n_samp < MIN_SAMPLES) next

        pcid <- pcid + 1L
        total_sr <- sum(pair_sub$support_sr, na.rm = TRUE)
        total_pe <- sum(pair_sub$support_pe, na.rm = TRUE)

        result_list[[pcid]] <- data.table(
          paired_cluster_id = pcid,
          chrom = chr,
          bp1_median = as.integer(median(pair_sub$pos1)),
          bp2_median = as.integer(median(pair_sub$pos2)),
          bp1_start = min(pair_sub$pos1),
          bp1_end = max(pair_sub$pos1),
          bp2_start = min(pair_sub$pos2),
          bp2_end = max(pair_sub$pos2),
          bp1_span = max(pair_sub$pos1) - min(pair_sub$pos1),
          bp2_span = max(pair_sub$pos2) - min(pair_sub$pos2),
          interval_size = as.integer(median(pair_sub$pos2)) - as.integer(median(pair_sub$pos1)),
          n_calls = nrow(pair_sub),
          n_unique_samples = n_samp,
          n_inv = sum(toupper(pair_sub$svtype) == "INV"),
          n_bnd = sum(toupper(pair_sub$svtype) == "BND"),
          total_sr = total_sr,
          total_pe = total_pe,
          total_support = total_sr + total_pe,
          svtype_mode = names(sort(table(toupper(pair_sub$svtype)), decreasing = TRUE))[1],
          ct_mode = names(sort(table(pair_sub$ct), decreasing = TRUE))[1],
          callers = paste(sort(unique(pair_sub$caller)), collapse = ","),
          samples = paste(sort(unique(pair_sub$sample)), collapse = ";")
        )
      }
    }
  }

  if (length(result_list) == 0) return(data.table())
  rbindlist(result_list)
}

paired_clusters <- build_paired_clusters(all_sv)
message("[STEP39] Paired interval clusters: ", nrow(paired_clusters))

# =============================================================================
# STEP 3 — MAP CLUSTERS TO CANDIDATES (for focal highlighting)
# =============================================================================

if (nrow(paired_clusters) > 0 && nrow(cand) > 0) {
  paired_clusters[, focal_candidate_id := NA_integer_]
  for (ci in seq_len(nrow(cand))) {
    cr <- cand[ci]
    mask <- paired_clusters$chrom == cr$chrom &
            abs(paired_clusters$bp1_median - as.numeric(cr$start_bp)) < BOUNDARY_BP &
            abs(paired_clusters$bp2_median - as.numeric(cr$end_bp)) < BOUNDARY_BP
    paired_clusters[mask & is.na(focal_candidate_id),
                    focal_candidate_id := cr$candidate_id]
  }
  paired_clusters[, is_focal := !is.na(focal_candidate_id)]
}

# =============================================================================
# STEP 4 — BUILD PLOT-READY TABLES
# =============================================================================

# Density table: one row per SV call with cluster membership
density_dt <- copy(all_sv)
density_dt[, svtype_upper := toupper(svtype)]
if (!"ct" %in% names(density_dt)) density_dt[, ct := "unknown"]

# Add focal flag from paired clusters
density_dt[, focal_candidate_id := NA_integer_]
if (nrow(paired_clusters) > 0) {
  for (ci in seq_len(nrow(cand))) {
    cr <- cand[ci]
    mask <- density_dt$chrom1 == cr$chrom &
            ((abs(density_dt$pos1 - as.numeric(cr$start_bp)) < BOUNDARY_BP) |
             (abs(density_dt$pos1 - as.numeric(cr$end_bp)) < BOUNDARY_BP))
    density_dt[mask & is.na(focal_candidate_id),
               focal_candidate_id := cr$candidate_id]
  }
}
density_dt[, is_focal := !is.na(focal_candidate_id)]

# Link table: one row per paired cluster
link_dt <- copy(paired_clusters)

# =============================================================================
# WRITE PLOT-READY TABLES
# =============================================================================

f1 <- file.path(VIZ_DIR, "breakpoint_paired_clusters.tsv.gz")
f2 <- file.path(VIZ_DIR, "breakpoint_density_plot_table.tsv.gz")
f3 <- file.path(VIZ_DIR, "breakpoint_link_plot_table.tsv.gz")

fwrite(paired_clusters, f1, sep = "\t")
fwrite(density_dt, f2, sep = "\t")
fwrite(link_dt, f3, sep = "\t")

message("[STEP39] Plot tables written")
message("  ", f1)
message("  ", f2)
message("  ", f3)

# =============================================================================
# PLOTTING HELPERS
# =============================================================================

SVTYPE_COLORS <- c("INV" = "#2563eb", "BND" = "#ea580c", "COMBINED" = "#7c3aed")
GREY_BG      <- "#b0b0b0"
FOCAL_ALPHA   <- 0.85
BG_ALPHA      <- 0.25

# Arc/ribbon helper: compute Bezier-like arc points for paired breakpoints
arc_points <- function(x_start, x_end, y_base = 0, n_pts = 50, arc_height = NULL) {
  if (is.null(arc_height)) arc_height <- (x_end - x_start) * 0.15
  t <- seq(0, pi, length.out = n_pts)
  x <- x_start + (x_end - x_start) * (1 - cos(t)) / 2
  y <- y_base + arc_height * sin(t)
  data.table(x = x, y = y)
}

# Build arc polygons for ribbons (width varies)
build_ribbon_arc <- function(x_start, x_end, width_samples, max_samples,
                             y_base = 0, max_height = NULL) {
  if (is.null(max_height)) max_height <- (x_end - x_start) * 0.2
  # Ribbon half-width proportional to sample count
  hw <- max(0.002, width_samples / max(max_samples, 1) * 0.03)

  n_pts <- 60
  t <- seq(0, pi, length.out = n_pts)
  x_mid <- (x_start + x_end) / 2
  half_span <- (x_end - x_start) / 2
  x <- x_mid - half_span * cos(t)
  y_center <- y_base + max_height * sin(t)

  data.table(
    x = c(x, rev(x)),
    y = c(y_center + hw * max_height, rev(y_center - hw * max_height))
  )
}

# =============================================================================
# FAMILY A — POPULATION OVERVIEW PLOT (per chromosome)
# =============================================================================

build_population_plot <- function(chr, density_dt, link_dt, cand_dt) {
  chr_density <- density_dt[chrom1 == chr]
  chr_links   <- link_dt[chrom == chr]

  if (nrow(chr_density) == 0 && nrow(chr_links) == 0) return(NULL)

  # Determine x range
  all_pos <- c(chr_density$pos1, chr_links$bp1_median, chr_links$bp2_median)
  x_min <- min(all_pos, na.rm = TRUE) - 100000
  x_max <- max(all_pos, na.rm = TRUE) + 100000

  # Candidate region shading
  chr_cands <- cand_dt[chrom == chr]

  p <- ggplot()

  # Candidate region background shading
  if (nrow(chr_cands) > 0) {
    for (ci in seq_len(nrow(chr_cands))) {
      p <- p + annotate("rect",
        xmin = as.numeric(chr_cands$start_bp[ci]),
        xmax = as.numeric(chr_cands$end_bp[ci]),
        ymin = -Inf, ymax = Inf, fill = "#ede9fe", alpha = 0.4)
    }
  }

  # ── Arc/ribbon links ────────────────────────────────────────────────
  if (nrow(chr_links) > 0) {
    max_samp <- max(chr_links$n_unique_samples, na.rm = TRUE)
    for (li in seq_len(nrow(chr_links))) {
      lk <- chr_links[li]
      arc_dt <- build_ribbon_arc(
        lk$bp1_median, lk$bp2_median,
        lk$n_unique_samples, max_samp,
        y_base = 0, max_height = max(0.5, lk$interval_size * 0.00001)
      )
      arc_dt[, group_id := li]

      # Support-based alpha
      support_alpha <- min(0.8, 0.2 + lk$total_support / max(100, max(chr_links$total_support)))
      sv_color <- SVTYPE_COLORS[lk$svtype_mode] %||% "#7c3aed"

      p <- p + geom_polygon(data = arc_dt, aes(x = x, y = y),
                             fill = sv_color, alpha = support_alpha,
                             color = sv_color, linewidth = 0.2)
    }
  }

  # ── Breakpoint density (jittered strip) ─────────────────────────────
  if (nrow(chr_density) > 0) {
    chr_density[, y_jitter := runif(.N, -0.15, -0.05)]
    # Line thickness encoded as point size proportional to support
    chr_density[, support_total := support_sr + support_pe]
    chr_density[, pt_size := pmin(3, 0.8 + support_total / max(1, max(support_total, na.rm = TRUE)) * 2)]

    p <- p + geom_point(data = chr_density,
                         aes(x = pos1, y = y_jitter,
                             color = svtype_upper, size = pt_size),
                         alpha = 0.6, shape = 16) +
      scale_size_identity()
  }

  # ── Median breakpoint lines ─────────────────────────────────────────
  if (nrow(chr_links) > 0) {
    # Left breakpoint medians
    p <- p + geom_segment(data = chr_links,
                           aes(x = bp1_median, xend = bp1_median,
                               y = -0.18, yend = -0.01),
                           color = "#1e3a5f", linewidth = 0.8, alpha = 0.7) +
      # Right breakpoint medians
      geom_segment(data = chr_links,
                   aes(x = bp2_median, xend = bp2_median,
                       y = -0.18, yend = -0.01),
                   color = "#5f1e1e", linewidth = 0.8, alpha = 0.7)

    # Sample count labels
    p <- p + geom_text(data = chr_links,
                        aes(x = (bp1_median + bp2_median) / 2,
                            y = max(0.5, interval_size * 0.000008),
                            label = paste0("n=", n_unique_samples)),
                        size = 2.5, color = "#334155", vjust = -0.3)
  }

  p <- p +
    scale_color_manual(values = SVTYPE_COLORS, name = "SV Type") +
    scale_x_continuous(labels = function(x) sprintf("%.1f", x / 1e6)) +
    labs(x = paste0(chr, " position (Mb)"), y = NULL,
         title = paste0(chr, " — Breakpoint Support Overview"),
         subtitle = paste0(nrow(chr_links), " paired clusters, ",
                           nrow(chr_density), " SV calls")) +
    theme_inversion(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
          panel.grid.major.y = element_blank())

  p
}

# =============================================================================
# FAMILY B — FOCAL CANDIDATE HIGHLIGHTED PLOT
# =============================================================================

build_focal_plot <- function(focal_cid, chr, density_dt, link_dt, cand_dt) {
  focal_cand <- cand_dt[candidate_id == focal_cid]
  if (nrow(focal_cand) == 0) return(NULL)

  c_start <- as.numeric(focal_cand$start_bp)
  c_end   <- as.numeric(focal_cand$end_bp)
  flank   <- max(200000, (c_end - c_start) * 0.5)
  x_min   <- c_start - flank
  x_max   <- c_end + flank

  chr_density <- density_dt[chrom1 == chr & pos1 >= x_min & pos1 <= x_max]
  chr_links   <- link_dt[chrom == chr &
                          bp1_median >= x_min & bp2_median <= x_max]

  if (nrow(chr_density) == 0 && nrow(chr_links) == 0) return(NULL)

  # Classify focal vs background
  chr_density[, plot_focal := (!is.na(focal_candidate_id) & focal_candidate_id == focal_cid)]
  chr_links[, plot_focal := (!is.na(focal_candidate_id) & focal_candidate_id == focal_cid)]

  p <- ggplot()

  # Candidate region shading
  p <- p + annotate("rect", xmin = c_start, xmax = c_end,
                     ymin = -Inf, ymax = Inf, fill = "#ede9fe", alpha = 0.4)
  p <- p + geom_vline(xintercept = c(c_start, c_end),
                       linetype = "dashed", color = "#7c3aed", linewidth = 0.5)

  # ── BACKGROUND arcs (grey, low alpha) ──────────────────────────────
  bg_links <- chr_links[plot_focal == FALSE]
  if (nrow(bg_links) > 0) {
    max_samp <- max(chr_links$n_unique_samples, na.rm = TRUE)
    for (li in seq_len(nrow(bg_links))) {
      lk <- bg_links[li]
      arc_dt <- build_ribbon_arc(lk$bp1_median, lk$bp2_median,
                                  lk$n_unique_samples, max_samp,
                                  y_base = 0,
                                  max_height = max(0.3, lk$interval_size * 0.000008))
      p <- p + geom_polygon(data = arc_dt, aes(x = x, y = y),
                             fill = GREY_BG, alpha = BG_ALPHA,
                             color = GREY_BG, linewidth = 0.1)
    }
  }

  # ── FOCAL arcs (colored, full alpha) ───────────────────────────────
  focal_links <- chr_links[plot_focal == TRUE]
  if (nrow(focal_links) > 0) {
    max_samp <- max(chr_links$n_unique_samples, na.rm = TRUE)
    for (li in seq_len(nrow(focal_links))) {
      lk <- focal_links[li]
      arc_dt <- build_ribbon_arc(lk$bp1_median, lk$bp2_median,
                                  lk$n_unique_samples, max_samp,
                                  y_base = 0,
                                  max_height = max(0.5, lk$interval_size * 0.00001))
      sv_color <- SVTYPE_COLORS[lk$svtype_mode] %||% "#7c3aed"
      support_alpha <- min(FOCAL_ALPHA, 0.4 + lk$total_support /
                             max(100, max(chr_links$total_support)))

      p <- p + geom_polygon(data = arc_dt, aes(x = x, y = y),
                             fill = sv_color, alpha = support_alpha,
                             color = sv_color, linewidth = 0.3)
    }
  }

  # ── Background density points (grey) ───────────────────────────────
  bg_pts <- chr_density[plot_focal == FALSE]
  if (nrow(bg_pts) > 0) {
    bg_pts[, y_jitter := runif(.N, -0.15, -0.05)]
    p <- p + geom_point(data = bg_pts, aes(x = pos1, y = y_jitter),
                         color = GREY_BG, alpha = BG_ALPHA, size = 1, shape = 16)
  }

  # ── Focal density points (colored) ─────────────────────────────────
  focal_pts <- chr_density[plot_focal == TRUE]
  if (nrow(focal_pts) > 0) {
    focal_pts[, y_jitter := runif(.N, -0.15, -0.05)]
    focal_pts[, support_total := support_sr + support_pe]
    focal_pts[, pt_size := pmin(3, 0.8 + support_total /
                                  max(1, max(support_total, na.rm = TRUE)) * 2)]
    p <- p + geom_point(data = focal_pts,
                         aes(x = pos1, y = y_jitter,
                             color = svtype_upper, size = pt_size),
                         alpha = FOCAL_ALPHA, shape = 16) +
      scale_size_identity()
  }

  # ── Focal median lines ─────────────────────────────────────────────
  if (nrow(focal_links) > 0) {
    p <- p +
      geom_segment(data = focal_links,
                   aes(x = bp1_median, xend = bp1_median, y = -0.18, yend = -0.01),
                   color = "#1e3a5f", linewidth = 1.0, alpha = 0.9) +
      geom_segment(data = focal_links,
                   aes(x = bp2_median, xend = bp2_median, y = -0.18, yend = -0.01),
                   color = "#5f1e1e", linewidth = 1.0, alpha = 0.9) +
      geom_text(data = focal_links,
                aes(x = (bp1_median + bp2_median) / 2,
                    y = max(0.4, interval_size * 0.000008),
                    label = paste0("n=", n_unique_samples, " (",
                                   svtype_mode, ")")),
                size = 3, color = "#1a1a2e", fontface = "bold", vjust = -0.3)
  }

  # Boundary annotations
  p <- p +
    annotate("text", x = c_start, y = -0.22,
             label = paste0("L: ", sprintf("%.2f", c_start / 1e6)),
             size = 2.8, color = "#7c3aed", hjust = 0.5) +
    annotate("text", x = c_end, y = -0.22,
             label = paste0("R: ", sprintf("%.2f", c_end / 1e6)),
             size = 2.8, color = "#7c3aed", hjust = 0.5)

  p <- p +
    scale_color_manual(values = SVTYPE_COLORS, name = "SV Type") +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6),
                       limits = c(x_min, x_max)) +
    labs(x = paste0(chr, " position (Mb)"), y = NULL,
         title = paste0("Candidate ", focal_cid, " — Breakpoint Support (Highlighted)"),
         subtitle = paste0(chr, ":", sprintf("%.2f", c_start / 1e6), "–",
                           sprintf("%.2f", c_end / 1e6), " Mb | ",
                           nrow(focal_links), " focal clusters, ",
                           nrow(focal_pts), " focal calls | ",
                           "grey = population background")) +
    theme_inversion(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
          panel.grid.major.y = element_blank())

  p
}

# =============================================================================
# GENERATE PLOTS
# =============================================================================

# ── Family A: Population overview per chromosome ─────────────────────
message("[STEP39] Generating population overview plots...")

pop_plots <- list()
for (chr in unique(cand$chrom)) {
  p <- build_population_plot(chr, density_dt, link_dt, cand)
  if (!is.null(p)) pop_plots[[chr]] <- p
}

if (length(pop_plots) > 0) {
  pop_pdf <- file.path(VIZ_DIR, "chromosome_breakpoint_overview.pdf")
  pdf(pop_pdf, width = 12, height = 5)
  for (chr in names(pop_plots)) {
    print(pop_plots[[chr]])
  }
  dev.off()
  message("[STEP39] Population overview: ", pop_pdf, " (", length(pop_plots), " pages)")
}

# ── Family B: Focal candidate highlighted ────────────────────────────
message("[STEP39] Generating focal candidate highlighted plots...")

for (ci in seq_len(nrow(cand_focal))) {
  focal_cid <- cand_focal$candidate_id[ci]
  focal_chr <- cand_focal$chrom[ci]

  p <- build_focal_plot(focal_cid, focal_chr, density_dt, link_dt, cand)
  if (!is.null(p)) {
    focal_pdf <- file.path(VIZ_DIR,
                           paste0("candidate_", focal_cid, "_breakpoint_highlighted.pdf"))
    ggsave(focal_pdf, p, width = 12, height = 5, device = "pdf")
    message("[STEP39] Focal plot: ", focal_pdf)
  }
}

# =============================================================================
# SUMMARY
# =============================================================================

message("\n[DONE] STEP39 breakpoint visualization complete")
message("  Output dir: ", VIZ_DIR)
message("  Paired clusters: ", nrow(paired_clusters))
if (nrow(paired_clusters) > 0) {
  message("  Focal clusters: ", sum(paired_clusters$is_focal, na.rm = TRUE))
  message("  Background clusters: ", sum(!paired_clusters$is_focal, na.rm = TRUE))
  sv_tab <- paired_clusters[, .N, by = svtype_mode]
  for (i in seq_len(nrow(sv_tab))) {
    message("    ", sv_tab$svtype_mode[i], ": ", sv_tab$N[i], " clusters")
  }
}
