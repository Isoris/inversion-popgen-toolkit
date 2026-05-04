#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01j_stream_graph.R  (v1.0)
#
# STREAM GRAPH of per-window compatibility cluster composition across a
# chromosome or candidate interval. Consumes C01j outputs directly.
#
# WHAT IT SHOWS:
#   For each window along the chromosome, the clusters that exist and how
#   many samples are in each. Plotted as a stacked area ("stream") graph so
#   that:
#     - Band thickness at position x   = number of samples in that cluster
#     - Total height at position x     = total samples (nearly constant)
#     - Band color                     = cluster's dominant dosage class
#                                          (REF_like / HET_like / INV_like)
#     - Color saturation / opacity     = cluster's thickness score
#                                          (from regime_windows.tsv.gz)
#
# BIOLOGICAL INTERPRETATION:
#   - Clean inversion region:  3 thick stable bands (REF/HET/INV)
#   - Background soup:         many thin noisy bands
#   - Composite / nested:      bands split or fuse mid-region
#   - Recombinant sample:      visible as an individual trajectory
#                              (from_cluster -> to_cluster) that crosses bands
#
# INPUTS (required):
#   --memberships <regime_memberships.tsv.gz>  — from C01j
#   --windows     <regime_windows.tsv.gz>      — from C01j
#
# INPUTS (optional):
#   --transitions <regime_transitions.tsv>     — from C01j, annotate verticals
#   --segments    <regime_segments.tsv[.gz]>   — from C01j, state ribbon below
#   --chrom       <chr>                        — filter if multi-chrom file
#   --xrange      <mb_start,mb_end>            — zoom in (e.g. "14,19")
#   --out         <out_path_base>              — without extension
#                                                 (writes .pdf and .png)
#   --sort_mode   "dosage" | "size"            — within same dosage class,
#                                                 sort bands by mean size
#                                                 (default "dosage")
#   --show_trajectory <sample>                 — overlay one sample's
#                                                 per-window band as a line
#
# OUTPUT:
#   <out>.pdf
#   <out>.png
#   <out>_bands.tsv   — the band table driving the plot (for auditability)
#
# Usage examples:
#   Rscript STEP_C01j_stream_graph.R \
#     --memberships results/c01j/regime_memberships.tsv.gz \
#     --windows     results/c01j/regime_windows.tsv.gz \
#     --transitions results/c01j/regime_transitions.tsv \
#     --segments    results/c01j/regime_segments.tsv.gz \
#     --chrom       C_gar_LG28 \
#     --out         results/c01j/stream_LG28
#
# REGISTRY WIRING (chat-18, minimal):
#   06 is a standalone visualization and stays CLI-driven. However, if it
#   can locate utils/registry_bridge.R it will source it opportunistically
#   so you can pass --candidate <cid> instead of --memberships / --windows
#   (in which case the tool looks up the C01j output files via the registry
#   convention: evidence_registry/per_candidate/<cid>/raw/regime_*).
#
#   If --candidate is NOT given, nothing registry-dependent runs — the tool
#   works exactly as before.
#
#   After the figure is written, the PDF path is registered via
#     reg$evidence$add_evidence(cid, "figure_stream_graph_path", ...)
#   if the bridge was sourced. Otherwise, the figure is simply written to
#   the --out path and nothing else happens.
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i >= length(args)) return(default)
  args[i + 1]
}

MEMB_FILE        <- get_arg("--memberships")
WIN_FILE         <- get_arg("--windows")
TRANS_FILE       <- get_arg("--transitions", NULL)
SEG_FILE         <- get_arg("--segments", NULL)
CHROM_FILTER     <- get_arg("--chrom", NULL)
XRANGE_STR       <- get_arg("--xrange", NULL)
OUT_BASE         <- get_arg("--out", "regime_stream_graph")
SORT_MODE        <- get_arg("--sort_mode", "dosage")
SHOW_TRAJECTORY  <- get_arg("--show_trajectory", NULL)
CANDIDATE_ID     <- get_arg("--candidate", NULL)

# ---- Optional registry bridge -----------------------------------------------
# If --candidate is given, we MUST have the bridge (to look up C01j paths).
# If --candidate is NOT given, the bridge is sourced if findable (so we can
# register the output PDF path), but its absence is not an error.
.have_reg <- FALSE
.bridge <- Sys.getenv("REGISTRY_BRIDGE", "utils/registry_bridge.R")
if (!file.exists(.bridge)) {
  for (p in c("utils/registry_bridge.R",
              "../utils/registry_bridge.R",
              file.path(Sys.getenv("BASE", ""), "utils/registry_bridge.R"))) {
    if (file.exists(p)) { .bridge <- p; break }
  }
}
if (file.exists(.bridge)) {
  Sys.setenv(CURRENT_SCRIPT = "06_regime_stream_graph.R")
  tryCatch({
    source(.bridge)
    .have_reg <- exists("reg") && !is.null(reg)
  }, error = function(e) {
    message("[stream] bridge source failed (", conditionMessage(e),
            ") — continuing in standalone mode")
  })
}

# If --candidate given, resolve C01j output paths from the evidence-registry
# convention. Overrides any explicit --memberships / --windows.
if (!is.null(CANDIDATE_ID)) {
  if (!.have_reg) {
    stop("[stream] --candidate ", CANDIDATE_ID,
         " requires utils/registry_bridge.R to be sourceable")
  }
  cand_raw <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                         "evidence_registry", "per_candidate", CANDIDATE_ID,
                         "raw")
  # Convention: C01j writes regime_memberships.tsv.gz, regime_windows.tsv.gz,
  # regime_transitions.tsv, regime_segments.tsv.gz under per-candidate raw/.
  resolved_memb <- file.path(cand_raw, "regime_memberships.tsv.gz")
  resolved_win  <- file.path(cand_raw, "regime_windows.tsv.gz")
  resolved_tr   <- file.path(cand_raw, "regime_transitions.tsv")
  resolved_seg  <- file.path(cand_raw, "regime_segments.tsv.gz")

  if (!file.exists(resolved_memb) || !file.exists(resolved_win)) {
    stop("[stream] C01j outputs not found under ", cand_raw,
         ". Expected regime_memberships.tsv.gz + regime_windows.tsv.gz. ",
         "Has C01j been run and its outputs placed in the evidence raw/ dir ",
         "for candidate ", CANDIDATE_ID, "?")
  }
  MEMB_FILE <- resolved_memb
  WIN_FILE  <- resolved_win
  if (is.null(TRANS_FILE) && file.exists(resolved_tr))  TRANS_FILE <- resolved_tr
  if (is.null(SEG_FILE)   && file.exists(resolved_seg)) SEG_FILE   <- resolved_seg
  # Default output path under evidence per-candidate figures/ if --out was
  # not explicitly set
  if (identical(OUT_BASE, "regime_stream_graph")) {
    fig_dir <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                          "evidence_registry", "per_candidate", CANDIDATE_ID,
                          "figures")
    if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
    OUT_BASE <- file.path(fig_dir, paste0("regime_stream_", CANDIDATE_ID))
  }
  # Pull chromosome from the candidate record if --chrom wasn't given
  if (is.null(CHROM_FILTER)) {
    cand_row <- tryCatch(reg$intervals$get_candidate(CANDIDATE_ID),
                          error = function(e) NULL)
    if (!is.null(cand_row) && !is.null(cand_row$chrom)) {
      CHROM_FILTER <- cand_row$chrom
    }
  }
  message("[stream] --candidate ", CANDIDATE_ID, " resolved to:")
  message("[stream]   MEMB_FILE = ", MEMB_FILE)
  message("[stream]   WIN_FILE  = ", WIN_FILE)
  message("[stream]   OUT_BASE  = ", OUT_BASE)
}

if (is.null(MEMB_FILE) || is.null(WIN_FILE)) {
  cat("ERROR: --memberships and --windows are required (or --candidate <cid>)\n\n")
  cat("Usage:\n")
  cat("  Rscript 06_regime_stream_graph.R --memberships <f> --windows <f> [--out <base>]\n")
  cat("  Rscript 06_regime_stream_graph.R --candidate LG28_1       (registry mode)\n")
  quit(status = 2)
}
if (!file.exists(MEMB_FILE)) stop("memberships file not found: ", MEMB_FILE)
if (!file.exists(WIN_FILE))  stop("windows file not found: ",     WIN_FILE)

ensure_dir <- function(p) {
  dp <- dirname(p)
  if (!dir.exists(dp)) dir.create(dp, recursive = TRUE)
  invisible(p)
}
ensure_dir(OUT_BASE)

# -----------------------------------------------------------------------------
# Load C01j outputs
# -----------------------------------------------------------------------------
message("[stream] reading: ", MEMB_FILE)
memb <- fread(MEMB_FILE)

message("[stream] reading: ", WIN_FILE)
win  <- fread(WIN_FILE)

required_memb_cols <- c("window", "sample", "group_id", "dosage_class", "pos_mid_mb")
missing_cols <- setdiff(required_memb_cols, names(memb))
if (length(missing_cols) > 0)
  stop("memberships file missing required columns: ",
       paste(missing_cols, collapse = ", "))

# Filter to chromosome if requested and column exists
if (!is.null(CHROM_FILTER) && "chrom" %in% names(memb)) {
  memb <- memb[chrom == CHROM_FILTER]
  if ("chrom" %in% names(win)) win <- win[chrom == CHROM_FILTER]
}
if (nrow(memb) == 0) stop("no rows after chrom filter")

# Zoom
if (!is.null(XRANGE_STR)) {
  parts <- as.numeric(strsplit(XRANGE_STR, ",")[[1]])
  if (length(parts) != 2 || any(!is.finite(parts)))
    stop("--xrange must be 'mb_start,mb_end', e.g. '14,19'")
  memb <- memb[pos_mid_mb >= parts[1] & pos_mid_mb <= parts[2]]
  win  <- if ("pos_mid_mb" %in% names(win))
    win[pos_mid_mb >= parts[1] & pos_mid_mb <= parts[2]] else win
  message("[stream] zoomed to ", parts[1], "-", parts[2], " Mb")
}

message("[stream] memberships: ", nrow(memb), " rows (",
        length(unique(memb$window)), " windows x ",
        length(unique(memb$sample)), " samples)")

# -----------------------------------------------------------------------------
# Cluster matching across windows (Hungarian assignment)
# -----------------------------------------------------------------------------
# C01j's group_id is arbitrary per window. To draw a coherent stream graph,
# we need global cluster identities. Strategy:
#   1. Walk windows in order
#   2. For each adjacent pair (w_{i-1}, w_i), build a confusion matrix of
#      sample -> (group_{i-1}, group_i)
#   3. Match groups to maximize overlap via greedy assignment
#      (Hungarian would be optimal; greedy is O(k^2) and sufficient for
#       k <= ~8 clusters)
#   4. Unmatched new clusters get fresh global IDs
#
# Output: per-(window,group_id) -> global_cluster_id
# -----------------------------------------------------------------------------
message("[stream] matching clusters across windows ...")

win_ids <- sort(unique(memb$window))
# Wide: row = sample, col = window, value = group_id
memb_dt <- copy(memb)
setkey(memb_dt, window, sample)

# Store assignments
n_windows <- length(win_ids)
# Each window has a vector: group_id -> global_id
global_map <- vector("list", n_windows)
names(global_map) <- as.character(win_ids)

# Initialize window 1 with identity mapping
w0 <- win_ids[1]
sub0 <- memb_dt[window == w0]
uniq_g0 <- sort(unique(sub0$group_id))
global_map[[as.character(w0)]] <- setNames(seq_along(uniq_g0), as.character(uniq_g0))
next_global_id <- length(uniq_g0) + 1L

# Sequentially match
for (k in 2:n_windows) {
  w_prev <- win_ids[k - 1]
  w_curr <- win_ids[k]
  sub_prev <- memb_dt[window == w_prev]
  sub_curr <- memb_dt[window == w_curr]
  # Align on sample
  both <- merge(sub_prev[, .(sample, prev_g = group_id)],
                sub_curr[, .(sample, curr_g = group_id)],
                by = "sample")
  if (nrow(both) == 0) {
    # Disjoint — assign all curr groups fresh IDs
    uniq_g <- sort(unique(sub_curr$group_id))
    global_map[[as.character(w_curr)]] <- setNames(
      next_global_id:(next_global_id + length(uniq_g) - 1L),
      as.character(uniq_g)
    )
    next_global_id <- next_global_id + length(uniq_g)
    next
  }
  # Confusion matrix
  conf <- table(both$prev_g, both$curr_g)
  prev_ids <- as.integer(rownames(conf))
  curr_ids <- as.integer(colnames(conf))

  # Greedy assignment: pick max overlaps one at a time, break ties by size
  conf_mat <- as.matrix(conf)
  assigned_prev <- rep(FALSE, length(prev_ids))
  assigned_curr <- rep(FALSE, length(curr_ids))
  mapping <- integer(0)
  names_mapping <- character(0)

  tmp_mat <- conf_mat
  while (TRUE) {
    if (all(tmp_mat == 0) || all(assigned_prev) || all(assigned_curr)) break
    ij <- which(tmp_mat == max(tmp_mat), arr.ind = TRUE)[1, ]
    i <- ij[1]; j <- ij[2]
    if (tmp_mat[i, j] == 0) break
    prev_g <- prev_ids[i]; curr_g <- curr_ids[j]
    prev_global <- global_map[[as.character(w_prev)]][as.character(prev_g)]
    mapping <- c(mapping, prev_global)
    names_mapping <- c(names_mapping, as.character(curr_g))
    assigned_prev[i] <- TRUE
    assigned_curr[j] <- TRUE
    tmp_mat[i, ] <- 0
    tmp_mat[, j] <- 0
  }
  # Unmatched curr groups get fresh global IDs
  unmatched_idx <- which(!assigned_curr)
  if (length(unmatched_idx) > 0) {
    for (j in unmatched_idx) {
      curr_g <- curr_ids[j]
      mapping <- c(mapping, next_global_id)
      names_mapping <- c(names_mapping, as.character(curr_g))
      next_global_id <- next_global_id + 1L
    }
  }
  names(mapping) <- names_mapping
  global_map[[as.character(w_curr)]] <- mapping
}

# Attach global_id to memberships
memb_dt[, global_id := {
  m <- global_map[[as.character(window[1])]]
  as.integer(m[as.character(group_id)])
}, by = window]

n_global <- length(unique(memb_dt$global_id))
message("[stream] ", n_global, " global clusters across ", n_windows, " windows")

# -----------------------------------------------------------------------------
# Per-cluster properties: dominant dosage class, mean size, mean thickness
# -----------------------------------------------------------------------------
cluster_info <- memb_dt[, {
  # Dominant dosage class for this global cluster
  dc_tab <- table(dosage_class)
  dom_dc <- names(dc_tab)[which.max(dc_tab)]
  list(
    dom_dosage_class = dom_dc,
    n_obs = .N,
    n_windows_present = uniqueN(window),
    mean_size = .N / uniqueN(window)
  )
}, by = global_id]

# Pull thickness from regime_windows.tsv.gz if available.
# win may have per-window columns like mean_thickness; we can't split to cluster
# level from that alone, so we approximate: thickness per cluster = mean_thickness
# of the window weighted by cluster's size in that window.
if ("mean_thickness" %in% names(win)) {
  # Average mean_thickness across windows, weighted by cluster presence
  win_small <- win[, .(window, mean_thickness)]
  memb_win <- merge(memb_dt[, .(window, global_id)], win_small, by = "window",
                     all.x = TRUE)
  thick_per_cluster <- memb_win[, .(
    mean_thickness_approx = round(mean(mean_thickness, na.rm = TRUE), 2)
  ), by = global_id]
  cluster_info <- merge(cluster_info, thick_per_cluster, by = "global_id", all.x = TRUE)
}

# -----------------------------------------------------------------------------
# Build band table: per window × per cluster -> count
# -----------------------------------------------------------------------------
band_dt <- memb_dt[, .(n_samples = .N), by = .(window, pos_mid_mb, global_id)]

# Merge in dominant dosage class for coloring
band_dt <- merge(band_dt,
                 cluster_info[, .(global_id, dom_dosage_class, mean_size)],
                 by = "global_id", all.x = TRUE)

# Sort: within each window, order clusters so bands stay vertically stable
# Sort key: (dosage class REF_like -> HET_like -> INV_like -> other),
#           then by mean_size desc
dc_order <- c("REF_like" = 1, "HET_like" = 2, "INV_like" = 3)
band_dt[, dc_rank := dc_order[dom_dosage_class]]
band_dt[is.na(dc_rank), dc_rank := 4L]

if (SORT_MODE == "size") {
  band_dt[, sort_key := -mean_size]
} else {
  band_dt[, sort_key := dc_rank * 1e6 - mean_size]
}
setorder(band_dt, window, sort_key)

# Compute stacked y range (center baseline so that the graph stays balanced)
band_dt[, cum_up := cumsum(n_samples), by = window]
band_dt[, total_w := sum(n_samples), by = window]
band_dt[, ymax := cum_up - total_w / 2]
band_dt[, ymin := ymax - n_samples]

# Save the band table for auditability
ensure_dir(paste0(OUT_BASE, "_bands.tsv"))
fwrite(band_dt[order(window, sort_key)], paste0(OUT_BASE, "_bands.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# Color palette
# -----------------------------------------------------------------------------
# Dominant-dosage-class base colors; clusters within the same class get varying
# hues/shades distinguished by global_id modulo some palette step.
base_colors <- c(
  "REF_like" = "#1f5d8e",   # blue
  "HET_like" = "#8d8d8d",   # grey
  "INV_like" = "#b03a2e"    # red
)

# Assign a shade to each global cluster based on mean_size (larger = more
# saturated). This helps distinguish multiple REF_like bands in soup regions.
cluster_info[, color_base := base_colors[dom_dosage_class]]
cluster_info[is.na(color_base), color_base := "#c49a6c"]  # "other" tan

# Shade: map mean_size to [0.55, 1.0] scaling
if (nrow(cluster_info) > 1) {
  ms_range <- range(cluster_info$mean_size)
  if (diff(ms_range) > 0) {
    cluster_info[, shade_factor := 0.55 + 0.45 *
                   (mean_size - ms_range[1]) / (ms_range[2] - ms_range[1])]
  } else {
    cluster_info[, shade_factor := 0.9]
  }
} else {
  cluster_info[, shade_factor := 0.9]
}

mix_with_white <- function(hex, f) {
  # f = 1 -> original color, f < 1 -> lighter
  rgb_col <- col2rgb(hex) / 255
  mixed <- rgb_col * f + (1 - f)
  grDevices::rgb(mixed[1], mixed[2], mixed[3])
}
cluster_info[, color := mapply(mix_with_white, color_base, shade_factor)]

# Attach to bands
band_dt <- merge(band_dt,
                 cluster_info[, .(global_id, color)],
                 by = "global_id", all.x = TRUE)

# -----------------------------------------------------------------------------
# Build polygons: each cluster is a polygon across windows
# (geom_ribbon needs continuous x -> use pos_mid_mb)
# -----------------------------------------------------------------------------
# We also need to handle gaps — if a cluster doesn't exist at window w, the
# ribbon should close. Simplest: compute ymin/ymax per (window, global_id);
# for missing (window, global_id), insert NA so ribbon breaks.

all_win <- unique(band_dt[, .(window, pos_mid_mb)])
setorder(all_win, window)
all_clusters <- cluster_info$global_id

full_grid <- CJ(window = all_win$window, global_id = all_clusters)
full_grid <- merge(full_grid, all_win, by = "window")
full_grid <- merge(full_grid,
                    band_dt[, .(window, global_id, n_samples, ymin, ymax, color)],
                    by = c("window", "global_id"), all.x = TRUE)
# For clusters absent at a window, the ribbon breaks (NA creates break in geom_ribbon)

# Colors need to be consistent per cluster (fill = factor(global_id), manual scale)
full_grid[, cluster_label := paste0("C", global_id)]
color_map <- setNames(cluster_info$color, paste0("C", cluster_info$global_id))

# -----------------------------------------------------------------------------
# Plot: stream graph
# -----------------------------------------------------------------------------
message("[stream] rendering ...")

xrange <- range(all_win$pos_mid_mb, na.rm = TRUE)
yrange <- range(c(band_dt$ymin, band_dt$ymax), na.rm = TRUE)
y_pad <- diff(yrange) * 0.05

main_title <- if (!is.null(CHROM_FILTER))
  sprintf("Regime stream graph  —  %s", CHROM_FILTER) else
  "Regime stream graph"
main_subtitle <- sprintf("%d global clusters across %d windows  |  %d samples total  |  C01j output",
                          n_global, n_windows, length(unique(memb_dt$sample)))

p <- ggplot(full_grid,
             aes(x = pos_mid_mb, ymin = ymin, ymax = ymax,
                  fill = cluster_label, group = cluster_label)) +
  geom_ribbon(alpha = 0.92, color = NA) +
  scale_fill_manual(values = color_map, guide = "none") +
  scale_x_continuous(limits = xrange, expand = expansion(0.005),
                       breaks = scales::pretty_breaks(n = 8)) +
  scale_y_continuous(limits = c(yrange[1] - y_pad, yrange[2] + y_pad),
                       expand = c(0, 0)) +
  labs(title = main_title, subtitle = main_subtitle,
       x = "Position (Mb)",
       y = "Samples (stacked, centered)") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "#555", size = 9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "#eeeeee", size = 0.2)
  )

# -----------------------------------------------------------------------------
# Optional overlay: regime transitions as vertical lines
# -----------------------------------------------------------------------------
if (!is.null(TRANS_FILE) && file.exists(TRANS_FILE)) {
  tr <- tryCatch(fread(TRANS_FILE), error = function(e) NULL)
  if (!is.null(tr) && nrow(tr) > 0) {
    if (!is.null(CHROM_FILTER) && "chrom" %in% names(tr))
      tr <- tr[chrom == CHROM_FILTER]
    pos_col <- if ("pos_mb" %in% names(tr)) "pos_mb"
                else if ("position_mb" %in% names(tr)) "position_mb"
                else if ("start_mb" %in% names(tr)) "start_mb"
                else if ("pos_bp" %in% names(tr)) "pos_bp"
                else NULL
    if (!is.null(pos_col) && nrow(tr) > 0) {
      tr_plot <- copy(tr)
      if (pos_col == "pos_bp") tr_plot[, pos_mb := get(pos_col) / 1e6]
      else setnames(tr_plot, pos_col, "pos_mb")
      tr_plot <- tr_plot[pos_mb >= xrange[1] & pos_mb <= xrange[2]]

      ttype_col <- if ("transition_type" %in% names(tr_plot)) "transition_type"
                    else if ("type" %in% names(tr_plot)) "type" else NULL
      if (!is.null(ttype_col) && nrow(tr_plot) > 0) {
        # Color by transition type
        tt_colors <- c(
          "enter_structured"   = "#005f73",
          "exit_structured"    = "#bb3e03",
          "complexity_increase" = "#ca6702",
          "complexity_decrease" = "#0a9396",
          "signal_weakening"   = "#ee9b00",
          "signal_strengthening" = "#94d2bd"
        )
        p <- p +
          geom_vline(data = tr_plot,
                      aes(xintercept = pos_mb, color = get(ttype_col)),
                      linetype = "dashed", size = 0.5, alpha = 0.7) +
          scale_color_manual(values = tt_colors,
                              name = "Transition",
                              na.value = "#888888") +
          theme(legend.position = "bottom",
                 legend.key.size = unit(0.35, "cm"),
                 legend.text = element_text(size = 7))
      }
    }
  }
}

# -----------------------------------------------------------------------------
# Optional: trajectory line for a specific sample
# -----------------------------------------------------------------------------
if (!is.null(SHOW_TRAJECTORY)) {
  traj <- memb_dt[sample == SHOW_TRAJECTORY]
  if (nrow(traj) > 0) {
    # For each (window, sample), find the y-center of the cluster it belongs to
    traj <- merge(traj,
                   band_dt[, .(window, global_id, ymin, ymax)],
                   by = c("window", "global_id"), all.x = TRUE)
    traj[, y_center := (ymin + ymax) / 2]
    traj <- traj[is.finite(y_center)]
    if (nrow(traj) > 0) {
      p <- p +
        geom_line(data = traj, aes(x = pos_mid_mb, y = y_center),
                   inherit.aes = FALSE,
                   size = 0.6, color = "black") +
        geom_point(data = traj, aes(x = pos_mid_mb, y = y_center),
                    inherit.aes = FALSE,
                    size = 1, color = "black") +
        labs(caption = paste0("Trajectory overlay: ", SHOW_TRAJECTORY))
    }
  }
}

# -----------------------------------------------------------------------------
# Optional: state ribbon below the stream graph (from regime_segments)
# -----------------------------------------------------------------------------
p_segments <- NULL
if (!is.null(SEG_FILE) && file.exists(SEG_FILE)) {
  seg <- tryCatch(fread(SEG_FILE), error = function(e) NULL)
  if (!is.null(seg) && nrow(seg) > 0 && "state" %in% names(seg)) {
    if (!is.null(CHROM_FILTER) && "chrom" %in% names(seg))
      seg <- seg[chrom == CHROM_FILTER]

    s_start_col <- if ("start_mb" %in% names(seg)) "start_mb"
                    else if ("start_bp" %in% names(seg)) "start_bp"
                    else NULL
    s_end_col <- if ("end_mb" %in% names(seg)) "end_mb"
                    else if ("end_bp" %in% names(seg)) "end_bp"
                    else NULL
    if (!is.null(s_start_col) && !is.null(s_end_col) && nrow(seg) > 0) {
      seg_plot <- copy(seg)
      if (s_start_col == "start_bp") {
        seg_plot[, start_mb := start_bp / 1e6]
        seg_plot[, end_mb   := end_bp / 1e6]
      } else {
        # already mb
      }
      seg_plot <- seg_plot[end_mb >= xrange[1] & start_mb <= xrange[2]]

      state_colors <- c(
        "clean_inversion"     = "#005f73",
        "structured_moderate" = "#0a9396",
        "structured_complex"  = "#94d2bd",
        "weak_signal"         = "#e9d8a6",
        "transition"          = "#ee9b00",
        "background_soup"     = "#bb3e03"
      )

      p_segments <- ggplot(seg_plot) +
        geom_rect(aes(xmin = start_mb, xmax = end_mb,
                       ymin = 0, ymax = 1, fill = state), color = NA) +
        scale_fill_manual(values = state_colors, drop = FALSE,
                           name = "State") +
        scale_x_continuous(limits = xrange, expand = expansion(0.005)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        labs(y = NULL, x = NULL) +
        theme_classic(base_size = 9) +
        theme(axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.line = element_blank(),
               legend.position = "bottom",
               legend.key.size = unit(0.3, "cm"),
               legend.text = element_text(size = 7),
               plot.margin = margin(0, 5, 5, 5))
    }
  }
}

# Combine p + p_segments if patchwork available
has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))
if (!is.null(p_segments) && has_patchwork) {
  combined <- p / p_segments + patchwork::plot_layout(heights = c(5, 0.6))
} else {
  combined <- p
}

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
pdf_path <- paste0(OUT_BASE, ".pdf")
png_path <- paste0(OUT_BASE, ".png")
width_in  <- max(10, min(18, 2 + 0.15 * n_windows))
height_in <- if (!is.null(p_segments)) 6 else 5

tryCatch(
  ggsave(pdf_path, combined, width = width_in, height = height_in,
          device = cairo_pdf),
  error = function(e) ggsave(pdf_path, combined, width = width_in, height = height_in)
)
tryCatch(
  ggsave(png_path, combined, width = width_in, height = height_in, dpi = 150),
  error = function(e) message("[stream] PNG failed: ", conditionMessage(e))
)

message("[stream] done.")
message("  PDF  -> ", pdf_path)
message("  PNG  -> ", png_path)
message("  TSV  -> ", paste0(OUT_BASE, "_bands.tsv"))

# Register the figure path in evidence_registry if we ran in registry mode
if (.have_reg && !is.null(CANDIDATE_ID)) {
  tryCatch({
    reg$evidence$add_evidence(CANDIDATE_ID, "figure_stream_graph_pdf_path",
                                pdf_path)
    reg$evidence$add_evidence(CANDIDATE_ID, "figure_stream_graph_png_path",
                                png_path)
    reg$evidence$add_evidence(CANDIDATE_ID, "figure_stream_graph_bands_tsv_path",
                                paste0(OUT_BASE, "_bands.tsv"))
    message("[stream] registered figure paths under candidate ", CANDIDATE_ID)
  }, error = function(e) {
    message("[stream] could not register figure paths: ", conditionMessage(e))
  })
}
