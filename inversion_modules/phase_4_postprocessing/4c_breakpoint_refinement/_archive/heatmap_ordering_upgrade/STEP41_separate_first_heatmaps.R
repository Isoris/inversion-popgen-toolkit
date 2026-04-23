#!/usr/bin/env Rscript

# =============================================================================
# STEP41_separate_first_heatmaps.R  v1
#
# Per-subgroup harmonized heatmaps implementing the "separate first, correct
# second" principle for composite candidate inversions.
#
# For each candidate classified as COMPOSITE_INTERNAL, COMPOSITE_OVERLAP, or
# UNRESOLVED, this script:
#
#   1. Reads candidate_sample_subgroups.tsv (from STEP40)
#   2. For each subgroup, computes LOCAL polarity harmonization using only
#      that subgroup's samples (so a single global model cannot average
#      incompatible systems together)
#   3. Renders a stack of heatmaps:
#        - raw global heatmap (reference)
#        - globally-harmonized heatmap (for comparison)
#        - per-subgroup locally-harmonized heatmaps
#        - marker-family-colored top annotation (which markers support which
#          sample partition)
#
# Output files per candidate:
#   heatmap_separate_first_{raw,global_harm,subgroup_N}.png/.pdf
#   candidate_separate_first_audit.tsv
#
# Usage:
#   Rscript STEP41_separate_first_heatmaps.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

has_CH <- suppressWarnings(require(ComplexHeatmap, quietly = TRUE))
has_circlize <- suppressWarnings(require(circlize, quietly = TRUE))
if (!has_CH || !has_circlize) {
  message("[WARN] ComplexHeatmap/circlize unavailable — skipping heatmaps")
  quit("no", status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_
if (file.exists(config_file)) source(config_file)

if (!exists("FOLLOWUP_DIR"))    FOLLOWUP_DIR    <- "results/followup"
if (!exists("DOSAGE_DIR"))      DOSAGE_DIR      <- "results/dosage"
if (!exists("CANDIDATE_TABLE")) CANDIDATE_TABLE <- "results/candidates.tsv"
if (!exists("PLOTS_DIR"))       PLOTS_DIR       <- "results/plots"
if (!exists("ensure_dir"))      ensure_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); invisible(p) }
if (!exists("GROUP_COLORS"))    GROUP_COLORS <- c(HOMO_1 = "#2166AC", HET = "#F4A582", HOMO_2 = "#B2182B")

# ── Parameters ─────────────────────────────────────────────────────────────
MAX_MARKERS_DISPLAY <- 300L
MIN_SUBGROUP_SIZE   <- 8L

# =============================================================================
# Local polarity harmonization WITHIN a set of samples
#
# This is the "correct second" step, applied AFTER samples have been separated.
# For each marker, we compute the sign of (Hom2_local_mean - Hom1_local_mean)
# but using only the current subgroup's members of each karyotype.
# If a subgroup lacks a karyotype (e.g., all Hom1 and Het, no Hom2), we fall
# back to a PC1-sign-based orientation computed within the subgroup.
# =============================================================================
local_harmonize <- function(X_sub, groups_sub, min_per_grp = 3L) {
  # X_sub: markers × samples (subset)
  # groups_sub: same length as ncol(X_sub), values in HOMO_1/HET/HOMO_2
  n_markers <- nrow(X_sub)
  h1 <- which(groups_sub == "HOMO_1")
  h2 <- which(groups_sub == "HOMO_2")

  flip <- rep(FALSE, n_markers)
  method <- "none"

  if (length(h1) >= min_per_grp && length(h2) >= min_per_grp) {
    # Group-difference sign
    mh1 <- rowMeans(X_sub[, h1, drop = FALSE], na.rm = TRUE)
    mh2 <- rowMeans(X_sub[, h2, drop = FALSE], na.rm = TRUE)
    delta <- mh2 - mh1
    delta[!is.finite(delta)] <- 0
    # Dominant sign in this subgroup
    abs_d <- abs(delta)
    dom_sign <- sign(median(delta[abs_d > 0.1], na.rm = TRUE))
    if (!is.finite(dom_sign) || dom_sign == 0) dom_sign <- 1
    flip <- (sign(delta) != dom_sign) & (abs_d > 0.1)
    method <- "local_group_difference"
  } else {
    # Fallback: PC1-sign polarity within subgroup
    M <- X_sub
    for (j in seq_len(ncol(M))) {
      bad <- !is.finite(M[, j])
      if (any(bad)) M[bad, j] <- mean(M[, j], na.rm = TRUE)
    }
    M[!is.finite(M)] <- 0
    v <- apply(M, 1, var, na.rm = TRUE)
    kk <- which(is.finite(v) & v > 1e-10)
    if (length(kk) >= 5) {
      pc <- tryCatch(prcomp(t(M[kk, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1),
                     error = function(e) NULL)
      if (!is.null(pc)) {
        # loadings sign; flip markers with negative loading
        loadings <- numeric(n_markers)
        loadings[kk] <- pc$rotation[, 1]
        dom_sign <- sign(median(loadings[abs(loadings) > quantile(abs(loadings), 0.5)]))
        if (!is.finite(dom_sign) || dom_sign == 0) dom_sign <- 1
        flip <- sign(loadings) != dom_sign & abs(loadings) > 0
        method <- "local_pc1_sign"
      }
    }
  }

  # Apply flip
  X_harm <- X_sub
  if (any(flip)) X_harm[flip, ] <- 2 - X_harm[flip, ]
  list(X = X_harm, flip = flip, method = method,
       n_flipped = sum(flip))
}

# =============================================================================
# Render helper
# =============================================================================
render_panel <- function(X_display, row_info, col_info, title, filepath,
                          show_marker_family = TRUE) {
  # X_display: samples × markers (rows = samples)
  if (nrow(X_display) < 3 || ncol(X_display) < 3) return(invisible())

  # Left row annotation: group + subgroup
  ann_list <- list(Group = row_info$coarse_group)
  col_list <- list(Group = GROUP_COLORS[intersect(names(GROUP_COLORS),
                                                   unique(row_info$coarse_group))])
  if ("subgroup_id" %in% names(row_info) && length(unique(row_info$subgroup_id)) > 1) {
    ann_list$Subgroup <- as.character(row_info$subgroup_id)
    subg_uniq <- as.character(sort(unique(row_info$subgroup_id)))
    pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
    col_list$Subgroup <- setNames(pal[seq_along(subg_uniq)], subg_uniq)
  }
  ha_left <- ComplexHeatmap::rowAnnotation(
    df = as.data.frame(ann_list), col = col_list,
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 7)
  )

  # Top column annotation: marker family + flip status
  ha_top <- NULL
  if (show_marker_family && "marker_family" %in% names(col_info)) {
    col_ann <- list()
    col_colors <- list()
    if ("marker_family" %in% names(col_info) && length(unique(col_info$marker_family)) > 1) {
      col_ann$Family <- as.character(col_info$marker_family)
      fam_uniq <- as.character(sort(unique(col_info$marker_family)))
      pal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
      col_colors$Family <- setNames(pal[seq_along(fam_uniq)], fam_uniq)
    }
    if ("flipped" %in% names(col_info)) {
      col_ann$Flipped <- as.character(col_info$flipped)
      col_colors$Flipped <- c("FALSE" = "#CCCCCC", "TRUE" = "#D73027")
    }
    if (length(col_ann) > 0) {
      ha_top <- ComplexHeatmap::columnAnnotation(
        df = as.data.frame(col_ann), col = col_colors,
        show_annotation_name = TRUE,
        annotation_name_gp = grid::gpar(fontsize = 7)
      )
    }
  }

  w <- max(8, ncol(X_display) * 0.03 + 5)
  h <- max(6, nrow(X_display) * 0.04 + 3)

  ht <- ComplexHeatmap::Heatmap(
    X_display, name = "Dosage",
    col = circlize::colorRamp2(c(0, 1, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    left_annotation = ha_left, top_annotation = ha_top,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    use_raster = (nrow(X_display) * ncol(X_display) > 50000)
  )

  pdf(sub("\\.png$", ".pdf", filepath), width = w, height = h)
  ComplexHeatmap::draw(ht); dev.off()
  png(filepath, width = w, height = h, units = "in", res = 200)
  ComplexHeatmap::draw(ht); dev.off()

  message("  Saved: ", basename(filepath))
}

# =============================================================================
# Main loop
# =============================================================================
cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id; chr <- row$chrom
  c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)
  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  plot_dir <- file.path(PLOTS_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(plot_dir)

  arch_file <- file.path(cand_dir, "candidate_architecture_class.tsv")
  if (!file.exists(arch_file)) {
    message("[SKIP] cid=", cid, " — no architecture class (run STEP40 first)")
    next
  }
  arch <- fread(arch_file)

  # Only run separate-first rendering for composite / unresolved classes
  if (!arch$architecture_class %in% c("COMPOSITE_INTERNAL", "COMPOSITE_OVERLAP", "UNRESOLVED")) {
    message("[SKIP] cid=", cid, " — class=", arch$architecture_class, " (no separation needed)")
    next
  }

  subg_file <- file.path(cand_dir, "candidate_sample_subgroups.tsv")
  fam_file  <- file.path(cand_dir, "candidate_marker_families.tsv")
  rot_file  <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  if (!all(file.exists(c(subg_file, fam_file, rot_file)))) {
    message("[SKIP] cid=", cid, " — missing subgroup or family files")
    next
  }

  subg <- fread(subg_file)
  fam  <- fread(fam_file)
  rot  <- fread(rot_file)

  # Load dosage
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) next
  dos <- fread(dos_file); sites <- fread(sites_file)
  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 30) next
  dos_reg <- dos[keep]; sites_reg <- sites[keep]

  sc <- setdiff(names(dos_reg), "marker")
  if (all(grepl("^Ind", sc)) && length(sc) == nrow(rot))
    setnames(dos_reg, old = sc, new = rot$sample)
  sample_cols <- intersect(subg$sample, setdiff(names(dos_reg), "marker"))
  if (length(sample_cols) < 10) next

  X <- as.matrix(dos_reg[, ..sample_cols]); storage.mode(X) <- "double"

  subg <- subg[match(sample_cols, sample)]
  groups <- subg$coarse_group

  message(sprintf("[INFO] cid=%d %s — class=%s, %d subgroups, %d markers, %d samples",
                  cid, chr, arch$architecture_class,
                  length(unique(subg$subgroup_id)),
                  nrow(X), ncol(X)))

  # Select top-variance markers for display
  mv <- apply(X, 1, var, na.rm = TRUE)
  n_show <- min(MAX_MARKERS_DISPLAY, nrow(X))
  top_idx <- sort(order(mv, decreasing = TRUE)[seq_len(n_show)])
  X_show <- X[top_idx, , drop = FALSE]

  # Map marker-family labels to display markers
  fam_vec <- fam[match(top_idx, marker_idx), marker_family]
  fam_vec[is.na(fam_vec)] <- 1L

  # ── Sample ordering: subgroup primary, then u-primary within each combo ──
  # Uses the shared helper so STEP41 ordering is consistent with STEP36.
  # Within each (subgroup_id, coarse_group) combination, we sort by u then
  # by a configurable secondary axis (default: hclust, matching STEP36 V3).
  # Set STEP41_SECONDARY_ORDER = "none" in the config to reproduce the
  # pre-upgrade behaviour exactly (setorder(subgroup, group, u)).
  STEP41_SECONDARY <- if (exists("STEP41_SECONDARY_ORDER")) STEP41_SECONDARY_ORDER else "hclust"

  # Resolve the helper
  if (!exists("order_samples_with_u_primary")) {
    helper_path <- file.path(dirname(sys.frame(1)$ofile %||% "."),
                              "within_group_ordering.R")
    if (file.exists(helper_path)) source(helper_path)
  }

  sub_order_dt <- data.table(
    sample = sample_cols,
    coarse_group = groups,
    subgroup_id = subg$subgroup_id,
    u = rot[match(sample_cols, sample), u],
    v = rot[match(sample_cols, sample), v]
  )

  # Primary: subgroup_id. Within each subgroup, apply the u-primary ordering
  # helper on that subgroup's samples separately.
  ordered_samples <- character(0)
  for (sg in sort(unique(sub_order_dt$subgroup_id))) {
    sg_rows <- sub_order_dt[subgroup_id == sg]
    if (nrow(sg_rows) < 1) next

    if (exists("order_samples_with_u_primary") && nrow(sg_rows) >= 5) {
      # Dosage restricted to this subgroup's samples
      sg_cols_all <- intersect(sg_rows$sample, colnames(X))
      X_sg <- X[, sg_cols_all, drop = FALSE]
      sg_ri_for_helper <- sg_rows[match(sg_cols_all, sample)]
      sg_res <- tryCatch(
        order_samples_with_u_primary(
          samples      = sg_cols_all,
          groups       = sg_ri_for_helper$coarse_group,
          u            = sg_ri_for_helper$u,
          v            = sg_ri_for_helper$v,
          dosage       = X_sg,
          secondary    = STEP41_SECONDARY,
          group_levels = c("HOMO_1", "HET", "HOMO_2"),
          min_group_n  = 5L,
          top_var_n    = 500L
        ),
        error = function(e) {
          message("  [WARN] helper failed for subgroup ", sg, ": ",
                  conditionMessage(e))
          NULL
        }
      )
      if (!is.null(sg_res)) {
        ordered_samples <- c(ordered_samples, sg_res$ordered_samples)
        next
      }
    }
    # Fallback: the old method
    setorder(sg_rows, coarse_group, u)
    ordered_samples <- c(ordered_samples, sg_rows$sample)
  }

  sub_order_dt <- sub_order_dt[match(ordered_samples, sample)]
  order_idx <- match(ordered_samples, sample_cols)
  X_ordered <- t(X_show[, order_idx, drop = FALSE])
  rownames(X_ordered) <- ordered_samples

  # ── Panel 1: RAW heatmap ──────────────────────────────────────────────
  col_info_raw <- data.table(marker_idx = top_idx, marker_family = fam_vec,
                             flipped = rep(FALSE, length(top_idx)))
  render_panel(
    X_ordered, sub_order_dt, col_info_raw,
    sprintf("RAW | cid=%d %s:%s-%s | %s",
            cid, chr, format(c_start, big.mark = ","),
            format(c_end, big.mark = ","), arch$label),
    file.path(plot_dir, "heatmap_separate_first_01_raw.png")
  )

  # ── Panel 2: GLOBAL harmonization (the naive approach, for comparison) ──
  # Use the existing candidate_marker_polarity.tsv if available; else compute.
  pol_file <- file.path(cand_dir, "candidate_marker_polarity.tsv")
  global_flip_full <- rep(FALSE, nrow(X))
  if (file.exists(pol_file)) {
    pol <- fread(pol_file)
    if (nrow(pol) == nrow(X) && "final_flip_decision" %in% names(pol)) {
      global_flip_full <- as.logical(pol$final_flip_decision)
    }
  }
  if (!any(global_flip_full) || all(is.na(global_flip_full))) {
    # Fallback: compute a simple global flip vector
    h1 <- which(groups == "HOMO_1"); h2 <- which(groups == "HOMO_2")
    if (length(h1) >= 3 && length(h2) >= 3) {
      delta <- rowMeans(X[, h2, drop = FALSE], na.rm = TRUE) -
               rowMeans(X[, h1, drop = FALSE], na.rm = TRUE)
      delta[!is.finite(delta)] <- 0
      dom <- sign(median(delta[abs(delta) > 0.1], na.rm = TRUE))
      if (!is.finite(dom) || dom == 0) dom <- 1
      global_flip_full <- (sign(delta) != dom) & (abs(delta) > 0.1)
    }
  }
  X_global_harm <- X
  rev_idx <- which(global_flip_full)
  if (length(rev_idx) > 0) X_global_harm[rev_idx, ] <- 2 - X_global_harm[rev_idx, ]
  X_gh_show <- t(X_global_harm[top_idx, order_idx, drop = FALSE])
  rownames(X_gh_show) <- ordered_samples
  col_info_gh <- data.table(marker_idx = top_idx, marker_family = fam_vec,
                            flipped = global_flip_full[top_idx])
  render_panel(
    X_gh_show, sub_order_dt, col_info_gh,
    sprintf("GLOBAL HARMONIZATION (reference) | cid=%d | %d/%d markers flipped",
            cid, sum(global_flip_full), length(global_flip_full)),
    file.path(plot_dir, "heatmap_separate_first_02_global_harm.png")
  )

  # ── Panel 3+: PER-SUBGROUP LOCAL HARMONIZATION ────────────────────────
  # Key principle: for each subgroup, compute polarity using ONLY that
  # subgroup's samples. Different subgroups may flip different markers.
  subgroup_ids <- sort(unique(subg$subgroup_id))
  audit_rows <- list()

  # Combined "separate-first" full heatmap: per-marker flip varies by subgroup.
  # We build a sample-indexed flip matrix: for each sample, use its subgroup's
  # local flip vector on every marker.
  flip_by_sample <- matrix(FALSE, nrow = nrow(X), ncol = ncol(X))
  colnames(flip_by_sample) <- sample_cols
  per_subgroup_flip <- list()

  for (sg in subgroup_ids) {
    samples_in_sg <- subg[subgroup_id == sg, sample]
    if (length(samples_in_sg) < MIN_SUBGROUP_SIZE) {
      message(sprintf("  [subgroup %d] only %d samples — skipped local harmonization",
                      sg, length(samples_in_sg)))
      next
    }
    sg_cols <- match(samples_in_sg, sample_cols)
    groups_sg <- groups[sg_cols]
    X_sg <- X[, sg_cols, drop = FALSE]
    lh <- local_harmonize(X_sg, groups_sg)
    per_subgroup_flip[[as.character(sg)]] <- lh$flip

    # Fill flip_by_sample
    for (s in sg_cols) flip_by_sample[lh$flip, s] <- TRUE

    # Disagreement vs global flip
    disagree_with_global <- sum(lh$flip != global_flip_full)

    audit_rows[[length(audit_rows) + 1]] <- data.table(
      candidate_id = cid, subgroup_id = sg,
      n_samples = length(samples_in_sg),
      n_hom1 = sum(groups_sg == "HOMO_1"),
      n_het  = sum(groups_sg == "HET"),
      n_hom2 = sum(groups_sg == "HOMO_2"),
      harmonization_method = lh$method,
      n_markers_flipped = lh$n_flipped,
      n_disagree_with_global = disagree_with_global,
      frac_disagree_with_global = round(disagree_with_global / nrow(X), 4)
    )

    # Render per-subgroup panel
    sg_display_samples <- intersect(ordered_samples, samples_in_sg)
    if (length(sg_display_samples) < 3) next
    sg_idx_in_order <- match(sg_display_samples, ordered_samples)

    # Use local flip on X_show
    X_sg_local_harm <- X_show
    rev_i <- which(lh$flip[top_idx])
    if (length(rev_i) > 0) X_sg_local_harm[rev_i, ] <- 2 - X_sg_local_harm[rev_i, ]
    X_sg_show <- t(X_sg_local_harm[, match(sg_display_samples, sample_cols), drop = FALSE])
    rownames(X_sg_show) <- sg_display_samples

    sg_row_info <- sub_order_dt[sample %in% sg_display_samples]
    sg_row_info <- sg_row_info[match(sg_display_samples, sample)]
    col_info_sg <- data.table(marker_idx = top_idx, marker_family = fam_vec,
                              flipped = lh$flip[top_idx])

    render_panel(
      X_sg_show, sg_row_info, col_info_sg,
      sprintf("SEPARATE-FIRST subgroup %d | cid=%d | %s | n=%d | flipped=%d (%s)",
              sg, cid, arch$label, length(sg_display_samples), lh$n_flipped, lh$method),
      file.path(plot_dir, sprintf("heatmap_separate_first_03_subgroup_%d.png", sg))
    )
  }

  # ── Panel: combined separate-first heatmap (each sample uses its
  # subgroup's local polarity). This is the key novel visualization.
  X_sepfirst <- X
  for (s_j in seq_len(ncol(X_sepfirst))) {
    fvec <- flip_by_sample[, s_j]
    if (any(fvec)) X_sepfirst[fvec, s_j] <- 2 - X_sepfirst[fvec, s_j]
  }
  X_sf_show <- t(X_sepfirst[top_idx, order_idx, drop = FALSE])
  rownames(X_sf_show) <- ordered_samples
  # "Flipped" on the top annotation here is ambiguous (varies by sample) — use
  # the "fraction of samples where this marker is flipped" as a continuous
  # signal, but for binary annotation we mark >50% flipped as TRUE.
  frac_flipped <- rowMeans(flip_by_sample)
  col_info_sf <- data.table(marker_idx = top_idx, marker_family = fam_vec,
                            flipped = frac_flipped[top_idx] > 0.5)
  render_panel(
    X_sf_show, sub_order_dt, col_info_sf,
    sprintf("SEPARATE-FIRST (combined) | cid=%d %s | each subgroup polarity-harmonized separately",
            cid, arch$label),
    file.path(plot_dir, "heatmap_separate_first_04_combined.png")
  )

  # ── Write audit table ─────────────────────────────────────────────────
  if (length(audit_rows) > 0) {
    audit_dt <- rbindlist(audit_rows)
    fwrite(audit_dt, file.path(cand_dir, "candidate_separate_first_audit.tsv"), sep = "\t")
    message(sprintf("  Audit: %d subgroups harmonized, max disagreement with global = %.1f%%",
                    nrow(audit_dt),
                    100 * max(audit_dt$frac_disagree_with_global, na.rm = TRUE)))
  }
}

message("[DONE] STEP41 complete")
