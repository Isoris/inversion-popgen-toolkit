# ============================================================================
# C01a_diag_matrices_NN_PATCH.R
# ============================================================================
#
# Adds NN-smoothed sim_mat heatmaps to the C01a diag plots.
#
# The existing diag only plots the raw nn0 sim_mat (build_sim_heatmap) and
# cropped-band views of it (build_local_sim_heatmap_80/_160). Those are not
# the same thing as the NN-smoothed sim_mats that D02/D09 actually consume.
# NN smoothing averages MDS coords across k-nearest-neighbors in MDS space
# and then rebuilds the sim_mat — blocks that are noise disappear at coarser
# scales while real inversions persist. That's the signal we want to SEE.
#
# This patch adds:
#   18_nn_smoothed_heatmap_nn20    — one heatmap per scale, full matrix
#   18_nn_smoothed_heatmap_nn40
#   18_nn_smoothed_heatmap_nn80
#   18_nn_smoothed_heatmap_nn160
#   18_nn_smoothed_heatmap_nn320
#   19_nn_smoothed_grid            — all scales side-by-side in one panel
#
# File convention it looks for:
#   <precomp_dir>/sim_mats/<chr>.sim_mat_nn<k>.rds
#   <precomp_dir>/<chr>.sim_mat_nn<k>.rds
#   <precomp_dir>/sim_mat_nn<k>/<chr>.rds
# First pattern that resolves wins.
#
# If no sim_mat_nn files found, these builders return NULL (silent skip).
#
# USAGE
# -----
#   A) Drop into STEP_C01a_diag_matrices.R — append this whole file to the
#      bottom, then add the new entries to plot_builders and plot_dims.
#
#   B) Source it after STEP_C01a_diag_matrices.R to patch at runtime:
#         source("STEP_C01a_diag_common.R")
#         source("STEP_C01a_diag_matrices.R")
#         source("C01a_diag_matrices_NN_PATCH.R")   # <- injects new builders
#         render_plots(plot_builders, plot_dims, chroms, precomp_list, outdir)
#
# DEPENDENCIES (already present in the diag environment)
#   data.table, ggplot2, scales, patchwork (optional for grid)
# ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    # Grid composite will fall back to individual plots
  }
})


# ----------------------------------------------------------------------------
# Locate and load a single sim_mat_nn<k>.rds for a given chromosome
# ----------------------------------------------------------------------------

.find_nn_simmat_file <- function(precomp_dir, chr, k) {
  candidates <- c(
    file.path(precomp_dir, "sim_mats",
              sprintf("%s.sim_mat_nn%d.rds", chr, k)),
    file.path(precomp_dir, "sim_mats",
              sprintf("%s_sim_mat_nn%d.rds", chr, k)),
    file.path(precomp_dir,
              sprintf("%s.sim_mat_nn%d.rds", chr, k)),
    file.path(precomp_dir,
              sprintf("%s_sim_mat_nn%d.rds", chr, k)),
    file.path(precomp_dir, sprintf("sim_mat_nn%d", k),
              sprintf("%s.rds", chr)),
    file.path(precomp_dir, sprintf("sim_mat_nn%d", k),
              sprintf("%s_sim_mat.rds", chr))
  )
  for (f in candidates) if (file.exists(f)) return(f)
  NA_character_
}


load_nn_simmat <- function(precomp_dir, chr, k) {
  f <- .find_nn_simmat_file(precomp_dir, chr, k)
  if (is.na(f)) return(NULL)
  obj <- tryCatch(readRDS(f), error = function(e) {
    message("[diag nn", k, "] failed to read ", f, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(obj)) return(NULL)

  # Some pipelines saved a list with $sim_mat, others saved the matrix directly
  if (is.list(obj) && !is.null(obj$sim_mat)) return(obj$sim_mat)
  if (is.matrix(obj)) return(obj)
  if (is.numeric(obj) && length(dim(obj)) == 2) return(obj)
  message("[diag nn", k, "] unrecognised RDS schema at ", f)
  NULL
}


# ----------------------------------------------------------------------------
# Core NN heatmap builder (one scale)
# ----------------------------------------------------------------------------

build_nn_smoothed_heatmap <- function(chr, k,
                                       precomp_dir = NULL,
                                       axis_mode = "mb") {
  # Get precomp_dir from env if not passed
  if (is.null(precomp_dir)) {
    precomp_dir <- if (exists("cli") && !is.null(cli$precomp_dir)) cli$precomp_dir
                   else if (exists("PRECOMP_DIR"))  PRECOMP_DIR
                   else stop("precomp_dir not resolvable")
  }

  # nn0 is in the precomp RDS itself; other k come from sim_mats/
  if (k == 0L) {
    sim_mat <- precomp_list[[chr]]$sim_mat
  } else {
    sim_mat <- load_nn_simmat(precomp_dir, chr, k)
    if (is.null(sim_mat)) return(NULL)
  }

  dt <- precomp_list[[chr]]$dt
  n  <- nrow(sim_mat)
  if (is.null(n) || n < 20) return(NULL)

  step <- max(1L, n %/% 800L)     # downsample for visualization
  idx  <- seq(1L, n, by = step)

  pdt <- if (exists("build_simmat_tiledata", mode = "function")) {
    build_simmat_tiledata(sim_mat, dt, step, idx, axis_mode)
  } else {
    # Minimal inline fallback
    if (axis_mode == "mb" && "start_bp" %in% names(dt)) {
      x_vals <- (dt$start_bp[idx] + dt$end_bp[idx]) / 2e6
    } else {
      x_vals <- idx
    }
    data.table(x = rep(x_vals, length(idx)),
               y = rep(x_vals, each = length(idx)),
               sim = as.vector(sim_mat[idx, idx]))
  }

  med_sim <- median(pdt$sim, na.rm = TRUE)

  x_label <- if (axis_mode == "mb") paste0(chr, " (Mb)") else paste0(chr, " (window)")

  title_tag <- if (k == 0L) "nn0 (raw)" else paste0("nn", k, " (MDS-NN smoothed)")

  p <- ggplot(pdt, aes(x = x, y = y, fill = sim)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("#0C1E3C", "#2563EB", "#F8F9FA", "#E8913A", "#C0392B"),
      values  = scales::rescale(c(0, med_sim * 0.7, med_sim,
                                   med_sim * 1.15, 1)),
      name    = "Similarity"
    ) +
    coord_fixed() +
    labs(x = x_label, y = x_label,
         title    = paste0(chr, " — ", title_tag),
         subtitle = paste0("N = ", n, " windows | median = ",
                           round(med_sim, 3)),
         caption  = "White = median | Red = high (block persists at this scale) | Blue = low") +
    (if (exists("THEME_BASE")) THEME_BASE else theme_minimal(base_size = 9))

  p
}


# ----------------------------------------------------------------------------
# Per-scale wrappers (one entry per builder to register)
# ----------------------------------------------------------------------------

build_nn_heatmap_nn20  <- function(chr) build_nn_smoothed_heatmap(chr, 20L)
build_nn_heatmap_nn40  <- function(chr) build_nn_smoothed_heatmap(chr, 40L)
build_nn_heatmap_nn80  <- function(chr) build_nn_smoothed_heatmap(chr, 80L)
build_nn_heatmap_nn160 <- function(chr) build_nn_smoothed_heatmap(chr, 160L)
build_nn_heatmap_nn320 <- function(chr) build_nn_smoothed_heatmap(chr, 320L)


# ----------------------------------------------------------------------------
# Grid composite — all scales side by side in one figure
# ----------------------------------------------------------------------------

build_nn_smoothed_grid <- function(chr,
                                    scales = c(0L, 20L, 40L, 80L, 160L, 320L),
                                    precomp_dir = NULL) {

  if (is.null(precomp_dir)) {
    precomp_dir <- if (exists("cli") && !is.null(cli$precomp_dir)) cli$precomp_dir
                   else if (exists("PRECOMP_DIR")) PRECOMP_DIR
                   else stop("precomp_dir not resolvable")
  }

  plots <- list()
  for (k in scales) {
    p <- tryCatch(build_nn_smoothed_heatmap(chr, k, precomp_dir = precomp_dir),
                  error = function(e) {
                    message("[diag grid nn", k, "] ", conditionMessage(e))
                    NULL
                  })
    if (!is.null(p)) plots[[paste0("nn", k)]] <- p
  }

  if (length(plots) == 0) return(NULL)

  # Strip legends except the last one for the grid
  n_plots <- length(plots)
  for (i in seq_len(n_plots - 1L)) {
    plots[[i]] <- plots[[i]] + theme(legend.position = "none")
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    g <- Reduce(`+`, plots)
    ncols <- if (n_plots <= 3L) n_plots else 3L
    g <- g + patchwork::plot_layout(ncol = ncols, guides = "collect") +
         patchwork::plot_annotation(
           title = paste0(chr, " — NN-smoothed sim_mat across scales"),
           subtitle = paste0("Same data, different k. Persistent blocks = real ",
                             "inversions. Dissolving blocks = family LD or noise.")
         )
    return(g)
  }

  # Fallback: just the first plot (render_plots handles one ggplot)
  plots[[1]]
}


# ----------------------------------------------------------------------------
# Registration snippet — append to the existing plot_builders / plot_dims
# lists. If sourced AFTER STEP_C01a_diag_matrices.R, the existing lists are
# already defined; we append in place.
# ----------------------------------------------------------------------------

.NN_BUILDERS <- list(
  "18_nn_heatmap_nn20"   = build_nn_heatmap_nn20,
  "18_nn_heatmap_nn40"   = build_nn_heatmap_nn40,
  "18_nn_heatmap_nn80"   = build_nn_heatmap_nn80,
  "18_nn_heatmap_nn160"  = build_nn_heatmap_nn160,
  "18_nn_heatmap_nn320"  = build_nn_heatmap_nn320,
  "19_nn_smoothed_grid"  = build_nn_smoothed_grid
)

.NN_DIMS <- list(
  "18_nn_heatmap_nn20"   = c(w = 9,  h = 8),
  "18_nn_heatmap_nn40"   = c(w = 9,  h = 8),
  "18_nn_heatmap_nn80"   = c(w = 9,  h = 8),
  "18_nn_heatmap_nn160"  = c(w = 9,  h = 8),
  "18_nn_heatmap_nn320"  = c(w = 9,  h = 8),
  "19_nn_smoothed_grid"  = c(w = 24, h = 16)
)


if (exists("plot_builders") && is.list(plot_builders)) {
  plot_builders <- c(plot_builders, .NN_BUILDERS)
  message("[NN patch] Added ", length(.NN_BUILDERS),
          " builders (total now ", length(plot_builders), ")")
}

if (exists("plot_dims") && is.list(plot_dims)) {
  plot_dims <- c(plot_dims, .NN_DIMS)
}


# ----------------------------------------------------------------------------
# If the caller wants to restrict to available scales only (skip missing
# files silently instead of producing empty plots), they can call:
#   filter_nn_builders_to_available(precomp_dir, chroms)
# after sourcing.
# ----------------------------------------------------------------------------

filter_nn_builders_to_available <- function(precomp_dir, chroms) {
  missing_scales <- integer(0)
  for (k in c(20L, 40L, 80L, 160L, 320L)) {
    found <- any(vapply(chroms, function(ch)
      !is.na(.find_nn_simmat_file(precomp_dir, ch, k)),
      logical(1)))
    if (!found) missing_scales <- c(missing_scales, k)
  }
  for (k in missing_scales) {
    key <- paste0("18_nn_heatmap_nn", k)
    plot_builders[[key]] <<- NULL
    plot_dims[[key]]     <<- NULL
    message("[NN patch] No sim_mat_nn", k, " found — dropping ", key)
  }
  invisible(missing_scales)
}
