#!/usr/bin/env Rscript

# =============================================================================
# 05_four_encoding_diagnostic.R  (v1.1 — registry-wired Mode B)
#
# FOUR-ENCODING PARALLEL SAMPLE COMPARISON
#
# Sample × sample distance matrices computed under four different encodings of
# the same dosage data. The principle (from your earlier design): don't commit
# to one allele-polarity convention — compute under all plausible conventions
# and test whether the sample clustering is robust to the choice.
#
# The four encodings:
#
#   1. MINOR dosage        raw BEAGLE dosage                (minor-allele count)
#   2. MAJOR dosage        2 - minor                         (major-allele count)
#   3. 012 discrete        hard-called genotype             (threshold on dosage)
#   4. POLARIZED dosage    STEP29 L1+L2 polarity flip       (homozygote-anchored)
#
# For each encoding:
#   - compute pairwise sample Manhattan distance (markers in genomic order)
#   - cluster via hclust (ward.D2)
#   - project to 2D via classical MDS
#   - render a heatmap + dendrogram + MDS plot
#
# Final:
#   - side-by-side 4-panel comparison
#   - inter-encoding agreement table: ARI between each pair of clusterings
#     (if all four encodings agree on sample clustering, structure is robust;
#      if some disagree, the odd-encoding-out diagnoses which allele-coding
#      conventions break the signal)
#
# TWO MODES:
#   MODE A (standalone CLI):  Rscript 05_... --dosage ... --sites ... --groups ...
#      Unchanged from v1.0; useful for one-off diagnostic runs on any matrix.
#
#   MODE B (registry-driven): Rscript 05_... [cid=all]
#      Pulls dosage/sites from $DOSAGE_DIR, karyotype groups via
#      reg$samples$get_groups_for_candidate(cid), writes encoding_robustness
#      block via reg$evidence$write_block. This is what the pipeline wrapper
#      calls.
#
# Auto-detection: if the first non-flag arg ends in .R, Mode B (legacy config).
#                  else if the first non-flag arg is a simple token (no slashes,
#                    no --), Mode B with that token as cid_filter.
#                  else Mode A.
#
# OUTPUTS (Mode B, via registry):
#   evidence_registry/per_candidate/<cid>/structured/encoding_robustness.json
#   evidence_registry/per_candidate/<cid>/figures/encoding_heatmaps.pdf
#   evidence_registry/per_candidate/<cid>/figures/encoding_mds.pdf
#
# OUTPUTS (Mode A, to --outdir):
#   <label>_distance_{minor,major,012,polarized}.tsv.gz
#   <label>_clusters.tsv
#   <label>_agreement.tsv
#   <label>_heatmaps_panel.pdf + .png
#   <label>_mds_panel.pdf + .png
#   <label>_dendrograms.pdf
#
# Usage (Mode A):
#   Rscript 05_four_encoding_diagnostic.R \
#     --dosage    results/dosage/chrX.dosage.tsv.gz \
#     --sites     results/dosage/chrX.sites.tsv.gz \
#     --polarity  ...  --groups ... --outdir ... --label cand1
#
# Usage (Mode B):
#   Rscript 05_four_encoding_diagnostic.R LG28_1
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

has_patchwork  <- suppressWarnings(require(patchwork,  quietly = TRUE))
has_viridis    <- suppressWarnings(require(viridisLite, quietly = TRUE))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -----------------------------------------------------------------------------
# Argument parsing — MODE A (CLI flags) vs MODE B (registry-driven by cid)
#
# Decision rule:
#   MODE A if any arg starts with "--"
#   MODE B otherwise. If a single non-flag arg is given, treat it as cid.
#   Legacy MODE B (config.R + cid) also supported for back-compat.
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || i >= length(args)) return(default)
  args[i + 1]
}

# Decide mode
is_mode_a <- any(grepl("^--", args))
# Legacy: first arg is a .R file = config-driven Mode B
is_legacy_mode_b <- (length(args) >= 1 && !is_mode_a &&
                       file.exists(args[1]) && grepl("\\.(R|r)$", args[1]))

if (is_mode_a) {
  # ------- MODE A: standalone -------
  CONFIG_MODE   <- FALSE
  REGISTRY_MODE <- FALSE
  DOSAGE_FILE   <- get_arg("--dosage")
  SITES_FILE    <- get_arg("--sites",      NULL)
  POLARITY_FILE <- get_arg("--polarity",   NULL)
  GROUPS_FILE   <- get_arg("--groups",     NULL)
  OUTDIR        <- get_arg("--outdir",     "four_encoding_output")
  LABEL         <- get_arg("--label",      "four_encoding")
  GT_THRESH     <- as.numeric(get_arg("--gt_thresh",   "0.5"))
  MIN_MARKERS   <- as.integer(get_arg("--min_markers", "20"))
  MIN_SAMPLES   <- as.integer(get_arg("--min_samples", "10"))

  if (is.null(DOSAGE_FILE)) {
    cat("ERROR: --dosage is required (or pass cid for Mode B)\n")
    cat("See script header for usage.\n")
    quit(status = 2)
  }
  if (!file.exists(DOSAGE_FILE)) stop("dosage file not found: ", DOSAGE_FILE)

} else {
  # ------- MODE B: registry-driven -------
  # Source the bridge to get reg / smap / BRIDGE_PATHS
  Sys.setenv(CURRENT_SCRIPT = "05_four_encoding_diagnostic.R")
  .bridge <- Sys.getenv("REGISTRY_BRIDGE", "utils/registry_bridge.R")
  if (!file.exists(.bridge)) {
    for (p in c("utils/registry_bridge.R",
                "../utils/registry_bridge.R",
                file.path(Sys.getenv("BASE", ""), "utils/registry_bridge.R"))) {
      if (file.exists(p)) { .bridge <- p; break }
    }
  }
  if (!file.exists(.bridge)) {
    stop("[05_four_encoding] cannot locate utils/registry_bridge.R")
  }
  source(.bridge)

  CONFIG_MODE   <- TRUE
  REGISTRY_MODE <- TRUE

  # Legacy config.R still accepted (first arg as file), else parse as cid
  if (is_legacy_mode_b) {
    config_file <- args[1]
    cid_filter  <- if (length(args) >= 2 && args[2] != "all") args[2] else NA_character_
    source(config_file)
  } else {
    cid_filter  <- if (length(args) >= 1 && args[1] != "all") args[1] else NA_character_
  }

  # Paths: dosage still on disk (too large for registry); everything else
  # comes from reg.
  if (!exists("DOSAGE_DIR")) {
    DOSAGE_DIR <- Sys.getenv("DOSAGE_DIR", "")
    if (!nzchar(DOSAGE_DIR)) {
      DOSAGE_DIR <- file.path(BRIDGE_PATHS$BASE, "popstruct_thin",
                               "04_beagle_byRF_majmin")
    }
  }

  GT_THRESH    <- if (exists("BP05_GT_THRESH"))   BP05_GT_THRESH   else 0.5
  MIN_MARKERS  <- if (exists("BP05_MIN_MARKERS")) BP05_MIN_MARKERS else 20L
  MIN_SAMPLES  <- if (exists("BP05_MIN_SAMPLES")) BP05_MIN_SAMPLES else 10L

  # Build list of candidates from interval_registry
  cand_path <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                          "interval_registry", "candidate_intervals.tsv")
  if (!file.exists(cand_path)) {
    stop("[05_four_encoding] candidate_intervals.tsv not found at ", cand_path)
  }
  cand_table_loaded <- fread(cand_path)
  if (!is.na(cid_filter)) {
    cand_table_loaded <- cand_table_loaded[candidate_id == cid_filter]
  }
  message("[05_four_encoding] MODE B (registry): ", nrow(cand_table_loaded),
          " candidates")
}

ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  invisible(p)
}

# =============================================================================
# MAIN PROCESSING — wrapped in a function so it can be called either in
# standalone Mode A (once) or config-driven Mode B (per candidate in a loop).
# =============================================================================
run_four_encoding <- function(DOSAGE_FILE, SITES_FILE, POLARITY_FILE,
                                GROUPS_FILE, OUTDIR, LABEL,
                                GT_THRESH, MIN_MARKERS, MIN_SAMPLES,
                                region_start_bp = NA, region_end_bp = NA) {
  ensure_dir(OUTDIR)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
message("[4enc] loading dosage: ", DOSAGE_FILE)
dos <- fread(DOSAGE_FILE)
marker_col <- names(dos)[1]
sample_cols <- setdiff(names(dos), marker_col)

X <- as.matrix(dos[, ..sample_cols])
storage.mode(X) <- "double"
markers <- dos[[marker_col]]

n_m <- nrow(X); n_s <- ncol(X)
message("[4enc]   ", n_m, " markers × ", n_s, " samples")

if (n_s < MIN_SAMPLES) stop("too few samples: ", n_s, " < ", MIN_SAMPLES)

# Sites (optional — used for genomic ordering and for region filtering)
if (!is.null(SITES_FILE) && file.exists(SITES_FILE)) {
  sites <- fread(SITES_FILE)
  if ("marker" %in% names(sites) && "pos" %in% names(sites)) {
    # Reorder X so markers are in genomic order
    ord <- match(markers, sites$marker)
    if (all(!is.na(ord))) {
      genomic_pos <- sites$pos[ord]
      ord2 <- order(genomic_pos)
      X <- X[ord2, , drop = FALSE]
      markers <- markers[ord2]
      genomic_pos <- genomic_pos[ord2]
      message("[4enc]   markers reordered to genomic order")

      # Region filter if explicit bounds given (Mode B)
      if (is.finite(region_start_bp) && is.finite(region_end_bp)) {
        keep <- which(genomic_pos >= region_start_bp & genomic_pos <= region_end_bp)
        if (length(keep) < MIN_MARKERS) {
          message("[4enc]   only ", length(keep), " markers in region ",
                  region_start_bp, "-", region_end_bp, " — skip")
          return(invisible(NULL))
        }
        X <- X[keep, , drop = FALSE]
        markers <- markers[keep]
        genomic_pos <- genomic_pos[keep]
        n_m <- nrow(X)
        message("[4enc]   filtered to region: ", n_m, " markers")
      }
    }
  }
}

if (n_m < MIN_MARKERS) stop("too few markers: ", n_m, " < ", MIN_MARKERS)

# Polarity (optional — enables the 4th encoding)
polarity_sign <- NULL
if (!is.null(POLARITY_FILE) && file.exists(POLARITY_FILE)) {
  pol <- tryCatch(fread(POLARITY_FILE), error = function(e) NULL)
  if (!is.null(pol)) {
    # Auto-detect sign column
    sign_col <- NULL
    for (cand_col in c("l2_sign", "L2_sign", "l1_sign", "L1_sign",
                         "polarity", "sign", "flip")) {
      if (cand_col %in% names(pol)) { sign_col <- cand_col; break }
    }
    if (!is.null(sign_col) && "marker" %in% names(pol)) {
      ord_pol <- match(markers, pol$marker)
      if (sum(!is.na(ord_pol)) > n_m * 0.5) {
        polarity_sign <- as.numeric(pol[[sign_col]][ord_pol])
        polarity_sign[!is.finite(polarity_sign) | polarity_sign == 0] <- 1
        message("[4enc]   loaded polarity from column '", sign_col,
                "' (", sum(polarity_sign == -1), " flips / ", n_m, " markers)")
      } else {
        message("[4enc]   polarity file present but marker IDs don't match; skipping")
      }
    } else {
      message("[4enc]   polarity file present but no recognizable sign column; skipping")
    }
  }
}
has_polarity <- !is.null(polarity_sign)
n_encodings <- if (has_polarity) 4L else 3L
message("[4enc]   running ", n_encodings, " encodings")

# Groups (optional annotation)
group_vec <- rep("all", n_s)
if (!is.null(GROUPS_FILE) && file.exists(GROUPS_FILE)) {
  grp <- tryCatch(fread(GROUPS_FILE), error = function(e) NULL)
  if (!is.null(grp) && "sample" %in% names(grp)) {
    grp_col <- NULL
    for (cand_col in c("coarse_group_refined", "coarse_group", "group", "invgt"))
      if (cand_col %in% names(grp)) { grp_col <- cand_col; break }
    if (!is.null(grp_col)) {
      ord <- match(sample_cols, grp$sample)
      if (sum(!is.na(ord)) > 0) {
        group_vec <- as.character(grp[[grp_col]][ord])
        group_vec[is.na(group_vec)] <- "unknown"
        message("[4enc]   loaded groups from '", grp_col, "'")
      }
    }
  }
}
names(group_vec) <- sample_cols

# -----------------------------------------------------------------------------
# Build the four encoded matrices
# -----------------------------------------------------------------------------
build_encoded_matrices <- function(X, gt_thresh = 0.5, polarity_sign = NULL) {
  # Minor-allele dosage: identity (X is already minor-dosage convention)
  E_minor <- X

  # Major-allele dosage: 2 - minor
  E_major <- 2 - X

  # 012 discrete: threshold
  E_012 <- X
  E_012[X <  gt_thresh]             <- 0
  E_012[X >= gt_thresh & X < 2 - gt_thresh] <- 1
  E_012[X >= 2 - gt_thresh]         <- 2

  encodings <- list(
    minor = E_minor,
    major = E_major,
    `012` = E_012
  )

  if (!is.null(polarity_sign)) {
    # Apply L2 polarity: flip selected markers.
    # Convention: after flip, all markers canonicalized so homozygote-2 is "high".
    # Flip is: x_flipped = 2 - x when sign = -1, else x.
    E_pol <- X
    flip_mask <- polarity_sign == -1
    if (any(flip_mask)) E_pol[flip_mask, ] <- 2 - E_pol[flip_mask, ]
    encodings$polarized <- E_pol
  }

  encodings
}

enc_list <- build_encoded_matrices(X, gt_thresh = GT_THRESH,
                                     polarity_sign = polarity_sign)

# -----------------------------------------------------------------------------
# Compute pairwise sample Manhattan distance for each encoding
# -----------------------------------------------------------------------------
compute_sample_distance <- function(E, method = "manhattan") {
  # E: markers × samples
  # Returns: samples × samples dist matrix
  # Normalize by number of non-NA markers for fair comparison
  tE <- t(E)  # samples × markers
  # Manhattan distance (sum of absolute differences per pair)
  d <- as.matrix(dist(tE, method = method))
  # Normalize by markers (average per-marker distance)
  d / ncol(tE)
}

message("[4enc] computing sample × sample distances under each encoding ...")
dist_list <- lapply(names(enc_list), function(nm) {
  message("[4enc]   ", nm, " ...")
  compute_sample_distance(enc_list[[nm]], method = "manhattan")
})
names(dist_list) <- names(enc_list)

# -----------------------------------------------------------------------------
# Cluster each distance matrix
# -----------------------------------------------------------------------------
message("[4enc] hclust + cut ...")
hc_list <- lapply(dist_list, function(d) {
  tryCatch(hclust(as.dist(d), method = "ward.D2"), error = function(e) NULL)
})
# Cut at k = 3 (HOMO_1 / HET / HOMO_2) as baseline
K_CUT <- 3L
cluster_list <- lapply(hc_list, function(hc) {
  if (is.null(hc)) return(rep(NA_integer_, n_s))
  cutree(hc, k = K_CUT)
})
# Apply sample-column order to all clusters
for (nm in names(cluster_list))
  names(cluster_list[[nm]]) <- sample_cols

# -----------------------------------------------------------------------------
# Inter-encoding agreement: ARI per pair
# -----------------------------------------------------------------------------
compute_ari <- function(cl1, cl2) {
  ok <- !is.na(cl1) & !is.na(cl2)
  if (sum(ok) < 2) return(NA_real_)
  tab <- table(cl1[ok], cl2[ok])
  n <- sum(tab)
  sc <- sum(choose(tab, 2))
  sa <- sum(choose(rowSums(tab), 2))
  sb <- sum(choose(colSums(tab), 2))
  e <- sa * sb / choose(n, 2)
  mx <- 0.5 * (sa + sb)
  if (mx == e) 1 else (sc - e) / (mx - e)
}

enc_names <- names(cluster_list)
n_enc <- length(enc_names)
ari_mat <- matrix(NA_real_, n_enc, n_enc,
                   dimnames = list(enc_names, enc_names))
for (i in seq_along(enc_names)) {
  for (j in seq_along(enc_names)) {
    ari_mat[i, j] <- if (i == j) 1 else compute_ari(cluster_list[[i]],
                                                       cluster_list[[j]])
  }
}
message("[4enc] inter-encoding ARI:")
print(round(ari_mat, 3))

# -----------------------------------------------------------------------------
# MDS per encoding
# -----------------------------------------------------------------------------
mds_list <- lapply(dist_list, function(d) {
  fit <- tryCatch(cmdscale(as.dist(d), k = 2, eig = FALSE),
                   error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  colnames(fit) <- c("MDS1", "MDS2")
  data.table(sample = rownames(fit), MDS1 = fit[,1], MDS2 = fit[,2])
})

# -----------------------------------------------------------------------------
# Write tabular outputs
# -----------------------------------------------------------------------------
message("[4enc] writing distance matrices + cluster tables ...")

write_dist_tsv <- function(d, path) {
  dt <- as.data.table(d)
  dt <- cbind(data.table(sample = rownames(d)), dt)
  fwrite(dt, path, sep = "\t")
}

for (nm in names(dist_list)) {
  rownames(dist_list[[nm]]) <- sample_cols
  colnames(dist_list[[nm]]) <- sample_cols
  write_dist_tsv(dist_list[[nm]],
                  file.path(OUTDIR, paste0(LABEL, "_distance_", nm, ".tsv.gz")))
}

clusters_dt <- data.table(sample = sample_cols, group = group_vec)
for (nm in names(cluster_list))
  clusters_dt[, (paste0("cluster_", nm)) := cluster_list[[nm]]]
fwrite(clusters_dt, file.path(OUTDIR, paste0(LABEL, "_clusters.tsv")), sep = "\t")

agreement_dt <- data.table(
  encoding_A = rep(enc_names, each  = n_enc),
  encoding_B = rep(enc_names, times = n_enc),
  ARI        = as.vector(t(ari_mat))
)
fwrite(agreement_dt, file.path(OUTDIR, paste0(LABEL, "_agreement.tsv")), sep = "\t")

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
make_heatmap <- function(d, title_str) {
  # Reorder by hclust order for visual block structure
  hc <- tryCatch(hclust(as.dist(d), method = "ward.D2"), error = function(e) NULL)
  ord <- if (!is.null(hc)) hc$order else seq_len(nrow(d))
  d_ord <- d[ord, ord]
  dt <- as.data.table(d_ord)
  dt[, sample_row := factor(rownames(d_ord), levels = rownames(d_ord))]
  dt_m <- melt(dt, id.vars = "sample_row",
                variable.name = "sample_col", value.name = "distance")
  dt_m[, sample_col := factor(sample_col, levels = colnames(d_ord))]

  fill_scale <- if (has_viridis)
    scale_fill_viridis_c(option = "magma", direction = -1, name = "Distance") else
    scale_fill_gradient(low = "#ffffcc", high = "#800026", name = "Distance")

  ggplot(dt_m, aes(x = sample_col, y = sample_row, fill = distance)) +
    geom_raster() +
    fill_scale +
    labs(title = title_str, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      panel.grid  = element_blank(),
      plot.title  = element_text(size = 10, face = "bold"),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width  = unit(0.3, "cm")
    ) +
    coord_fixed()
}

make_mds <- function(mds_dt, clusters, title_str, group_map = NULL) {
  if (is.null(mds_dt)) return(NULL)
  pd <- copy(mds_dt)
  pd[, cluster := factor(clusters[match(sample, names(clusters))])]
  if (!is.null(group_map)) {
    pd[, group := group_map[match(sample, names(group_map))]]
  }
  cluster_colors <- c("1" = "#1f5d8e", "2" = "#8d8d8d", "3" = "#b03a2e",
                       "4" = "#6a4c93", "5" = "#2a9d8f", "6" = "#e76f51")
  p <- ggplot(pd, aes(x = MDS1, y = MDS2))
  if (!is.null(group_map)) {
    p <- p + geom_point(aes(color = cluster, shape = group), size = 2, alpha = 0.85) +
      scale_shape_manual(values = c(16, 17, 15, 4, 8, 3), name = "Group")
  } else {
    p <- p + geom_point(aes(color = cluster), size = 2, alpha = 0.85)
  }
  p + scale_color_manual(values = cluster_colors, drop = TRUE, name = "Cluster") +
    labs(title = title_str) +
    theme_classic(base_size = 9) +
    theme(
      plot.title  = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 7)
    )
}

# Panel A: heatmaps
message("[4enc] rendering heatmap panel ...")
heatmap_plots <- lapply(names(dist_list), function(nm) {
  make_heatmap(dist_list[[nm]],
                 paste0("Distance — encoding: ", nm))
})
names(heatmap_plots) <- names(dist_list)

# Panel B: MDS
message("[4enc] rendering MDS panel ...")
mds_plots <- lapply(names(mds_list), function(nm) {
  make_mds(mds_list[[nm]], cluster_list[[nm]],
             paste0("MDS — encoding: ", nm),
             group_map = group_vec)
})
names(mds_plots) <- names(mds_list)

# Combine via patchwork if available
if (has_patchwork) {
  n_plots <- length(heatmap_plots)
  ncol_panel <- if (n_plots <= 2) n_plots else 2

  hm_panel <- Reduce(`+`, heatmap_plots) +
    patchwork::plot_layout(ncol = ncol_panel, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(LABEL, " — 4-encoding sample × sample distance heatmaps"),
      subtitle = paste0(n_m, " markers × ", n_s, " samples  |  ",
                         "ARI between encodings in agreement.tsv"),
      theme = theme(plot.title = element_text(face = "bold", size = 12),
                     plot.subtitle = element_text(color = "#555", size = 9))
    )

  mds_panel <- Reduce(`+`, Filter(Negate(is.null), mds_plots)) +
    patchwork::plot_layout(ncol = ncol_panel, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(LABEL, " — 4-encoding MDS projections"),
      subtitle = "Clusters from ward.D2 at k=3  |  Shapes = external group assignments",
      theme = theme(plot.title = element_text(face = "bold", size = 12),
                     plot.subtitle = element_text(color = "#555", size = 9))
    )
} else {
  # Fall back to stacking plots into one page each
  hm_panel  <- heatmap_plots[[1]]
  mds_panel <- mds_plots[[1]]
  message("[4enc]   patchwork not available — saving only first panel of each")
}

panel_w <- if (n_encodings == 4) 12 else 10
panel_h <- if (n_encodings == 4) 11 else 6

hm_pdf <- file.path(OUTDIR, paste0(LABEL, "_heatmaps_panel.pdf"))
hm_png <- file.path(OUTDIR, paste0(LABEL, "_heatmaps_panel.png"))
mds_pdf <- file.path(OUTDIR, paste0(LABEL, "_mds_panel.pdf"))
mds_png <- file.path(OUTDIR, paste0(LABEL, "_mds_panel.png"))

tryCatch(ggsave(hm_pdf, hm_panel, width = panel_w, height = panel_h,
                  device = cairo_pdf),
          error = function(e) ggsave(hm_pdf, hm_panel, width = panel_w, height = panel_h))
tryCatch(ggsave(hm_png, hm_panel, width = panel_w, height = panel_h, dpi = 150),
          error = function(e) message("[4enc]   heatmap PNG failed: ", conditionMessage(e)))

tryCatch(ggsave(mds_pdf, mds_panel, width = panel_w, height = panel_h,
                  device = cairo_pdf),
          error = function(e) ggsave(mds_pdf, mds_panel, width = panel_w, height = panel_h))
tryCatch(ggsave(mds_png, mds_panel, width = panel_w, height = panel_h, dpi = 150),
          error = function(e) message("[4enc]   MDS PNG failed: ", conditionMessage(e)))

# Dendrograms: stack into one tall PDF
dend_pdf <- file.path(OUTDIR, paste0(LABEL, "_dendrograms.pdf"))
pdf(dend_pdf, width = 12, height = 3 * n_encodings)
par(mfrow = c(n_encodings, 1), mar = c(2, 4, 3, 1))
for (nm in names(hc_list)) {
  if (is.null(hc_list[[nm]])) next
  plot(hc_list[[nm]], main = paste0("Dendrogram — encoding: ", nm),
       xlab = "", sub = "", cex = 0.4)
  abline(h = mean(hc_list[[nm]]$height[seq(length(hc_list[[nm]]$height) - K_CUT + 1,
                                             length(hc_list[[nm]]$height))]),
         col = "red", lty = 2)
}
dev.off()

# -----------------------------------------------------------------------------
# Summary to stdout
# -----------------------------------------------------------------------------
message("[4enc] ============================================================")
message("[4enc] DONE  label=", LABEL)
message("[4enc]   outputs in: ", OUTDIR)
message("[4enc]   encodings run: ", paste(enc_names, collapse = ", "))
message("[4enc]   min ARI across encoding pairs: ", round(min(ari_mat[upper.tri(ari_mat)],
                                                                 na.rm = TRUE), 3))
message("[4enc]   max ARI across encoding pairs: ", round(max(ari_mat[upper.tri(ari_mat)],
                                                                 na.rm = TRUE), 3))
if (has_polarity) {
  message("[4enc]   polarity was applied — if polarized agrees with minor (ARI ≈ 1),")
  message("[4enc]   STEP29 is redundant. If polarized disagrees strongly, it's")
  message("[4enc]   capturing something the raw encodings miss.")
}
message("[4enc] ============================================================")

  invisible(list(
    encodings = names(enc_list),
    ari = ari_mat,
    clusters = cluster_list
  ))
}   # end run_four_encoding

# =============================================================================
# DRIVER — dispatch to Mode A (standalone) or Mode B (per candidate)
# =============================================================================
if (CONFIG_MODE) {
  # Mode B: iterate candidates from interval_registry; use reg$samples /
  # reg$evidence to look up groups, polarity, and write back the robustness
  # block. Falls back to per-candidate followup dir only when a required
  # artifact isn't registered yet.
  ensure_dir <- function(p) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE); invisible(p)
  }

  for (ci in seq_len(nrow(cand_table_loaded))) {
    row <- cand_table_loaded[ci]
    cid <- as.character(row$candidate_id)
    chr <- as.character(row$chrom)
    c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)

    # Dosage stays on disk (too large for registry)
    dos_file   <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
    sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
    if (!file.exists(dos_file)) {
      message("[05_four_encoding] cid=", cid, " SKIP — no dosage file for ", chr)
      next
    }

    # Groups: resolve via registry. Write a minimal (sample, group) TSV to a
    # temp path so run_four_encoding (which expects a file) can consume it.
    grs <- tryCatch(reg$samples$get_groups_for_candidate(cid),
                    error = function(e) NULL)
    grp_file <- NULL
    if (!is.null(grs) &&
        length(grs$HOM_REF) > 0L && length(grs$HET) > 0L &&
        length(grs$HOM_INV) > 0L) {
      grp_dt <- rbindlist(list(
        data.table(sample = grs$HOM_REF, group = "HOM_REF"),
        data.table(sample = grs$HET,     group = "HET"),
        data.table(sample = grs$HOM_INV, group = "HOM_INV")
      ))
      cand_raw_dir <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                                 "evidence_registry", "per_candidate", cid,
                                 "raw")
      ensure_dir(cand_raw_dir)
      grp_file <- file.path(cand_raw_dir, "groups_for_encoding.tsv")
      fwrite(grp_dt, grp_file, sep = "\t")
    } else {
      message("[05_four_encoding] cid=", cid, " no karyotype groups — running unannotated")
    }

    # Polarity: optional. Try to find a STEP29 output in evidence raw/. If not
    # there, skip the 4th (polarized) encoding — Mode B runs 3 encodings and
    # reports min_ari across those 3.
    pol_file <- NULL
    cand_raw_pol <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                                "evidence_registry", "per_candidate", cid,
                                "raw", "step29_polarity.tsv")
    if (file.exists(cand_raw_pol)) pol_file <- cand_raw_pol

    # Output dir: under evidence_registry per-candidate figures/
    outdir_cid <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                             "evidence_registry", "per_candidate", cid,
                             "figures", "encoding")
    ensure_dir(outdir_cid)
    label_cid <- paste0(cid, "_encoding")

    message("[05_four_encoding] cid=", cid, " ", chr, " ", c_start, "-", c_end)
    result <- tryCatch(
      run_four_encoding(
        DOSAGE_FILE     = dos_file,
        SITES_FILE      = sites_file,
        POLARITY_FILE   = pol_file,
        GROUPS_FILE     = grp_file,
        OUTDIR          = outdir_cid,
        LABEL           = label_cid,
        GT_THRESH       = GT_THRESH,
        MIN_MARKERS     = MIN_MARKERS,
        MIN_SAMPLES     = MIN_SAMPLES,
        region_start_bp = c_start,
        region_end_bp   = c_end
      ),
      error = function(e) {
        message("[05_four_encoding] cid=", cid, " ERROR: ", conditionMessage(e))
        NULL
      }
    )

    # Write the encoding_robustness block. If run_four_encoding didn't return
    # a structured summary (current v1.0 returns invisibly), read the ARI TSV
    # it wrote and summarize from that.
    ari_tsv <- file.path(outdir_cid, paste0(label_cid, "_agreement.tsv"))
    cluster_tsv <- file.path(outdir_cid, paste0(label_cid, "_clusters.tsv"))
    if (file.exists(ari_tsv)) {
      ari_dt <- tryCatch(fread(ari_tsv), error = function(e) NULL)
      if (!is.null(ari_dt) && nrow(ari_dt) > 0) {
        # Expect long-format: (encoding_a, encoding_b, ari). Pair keys use "_vs_"
        pair_key <- function(a, b) paste0(a, "_vs_", b)
        ari_map <- list()
        for (r in seq_len(nrow(ari_dt))) {
          a <- as.character(ari_dt$encoding_a[r])
          b <- as.character(ari_dt$encoding_b[r])
          if (a == b) next
          ari_map[[pair_key(a, b)]] <- as.numeric(ari_dt$ari[r])
        }
        ari_vals <- unlist(ari_map)
        min_ari  <- if (length(ari_vals)) min(ari_vals, na.rm = TRUE) else NA_real_
        mean_ari <- if (length(ari_vals)) mean(ari_vals, na.rm = TRUE) else NA_real_
        verdict <- if (!is.finite(min_ari)) NA_character_ else
                    if (min_ari >= 0.9) "robust" else
                    if (min_ari >= 0.5) "partial" else "ambiguous"

        # which encoding pairs fall below 0.9 (i.e., "disagreeing")
        disagreeing <- names(ari_map)[vapply(ari_map,
                                              function(v) isTRUE(v < 0.9),
                                              logical(1))]

        reg$evidence$write_block(cid, "encoding_robustness", list(
          status             = "ok",
          n_markers_used     = if ("n_markers" %in% names(ari_dt)) ari_dt$n_markers[1] else NA_integer_,
          n_samples_used     = if ("n_samples" %in% names(ari_dt)) ari_dt$n_samples[1] else NA_integer_,
          encodings_used     = unique(c(as.character(ari_dt$encoding_a),
                                         as.character(ari_dt$encoding_b))),
          ari_matrix         = ari_map,
          min_ari            = min_ari,
          mean_ari           = mean_ari,
          verdict            = verdict,
          clusters_disagreeing = if (length(disagreeing)) disagreeing else character(0),
          heatmap_pdf_path   = file.path(outdir_cid, paste0(label_cid, "_heatmaps_panel.pdf")),
          mds_pdf_path       = file.path(outdir_cid, paste0(label_cid, "_mds_panel.pdf"))
        ))
      } else {
        reg$evidence$write_block(cid, "encoding_robustness", list(
          status = "ok",
          reason = "ari_tsv empty — check run_four_encoding output",
          heatmap_pdf_path = file.path(outdir_cid, paste0(label_cid, "_heatmaps_panel.pdf")),
          mds_pdf_path     = file.path(outdir_cid, paste0(label_cid, "_mds_panel.pdf"))
        ))
      }
    } else {
      reg$evidence$write_block(cid, "encoding_robustness", list(
        status = "insufficient_markers",
        reason = paste0("ARI tsv not produced at ", ari_tsv,
                        " — run_four_encoding likely skipped this candidate"),
        encodings_used = character(0)
      ))
    }
  }
  message("[05_four_encoding] DONE (Mode B — registry)")
} else {
  # Mode A: single standalone run
  run_four_encoding(
    DOSAGE_FILE   = DOSAGE_FILE,
    SITES_FILE    = SITES_FILE,
    POLARITY_FILE = POLARITY_FILE,
    GROUPS_FILE   = GROUPS_FILE,
    OUTDIR        = OUTDIR,
    LABEL         = LABEL,
    GT_THRESH     = GT_THRESH,
    MIN_MARKERS   = MIN_MARKERS,
    MIN_SAMPLES   = MIN_SAMPLES
  )
}
