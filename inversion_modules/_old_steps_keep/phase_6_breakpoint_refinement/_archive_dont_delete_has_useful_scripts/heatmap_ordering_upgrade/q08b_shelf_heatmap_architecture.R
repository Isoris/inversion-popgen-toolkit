#!/usr/bin/env Rscript
# =============================================================================
# q08b_shelf_heatmap_architecture.R   v3
#
# Upgrade of q08_shelf_heatmap.R v2 that integrates with the new candidate
# architecture classifier (STEP40) and the separate-first heatmap workflow
# (STEP41).
#
# Backward compatible: if neither classifier output nor subgroup file exists,
# this script behaves exactly like v2 — same 9 panels, same CLI, same output
# paths. No disruption to any existing run_chrom.sh invocation.
#
# New when STEP40/STEP41 outputs are present:
#
#   1. Title overlays
#      Every panel title gets the architecture label appended, e.g.
#         "A5: 1 top-var / 500-SNP window  [COMPOSITE_INTERNAL-B]"
#      Classifier feature values shown in a small caption line under each
#      panel: simple_axis, boundary, marker_family_k, monotonic_fraction.
#
#   2. Subgroup row annotation
#      ComplexHeatmap left annotation gets a second column, "subgroup", next
#      to the existing "invgt" column. Samples within the same subgroup are
#      visually contiguous (subgroup is the primary sort key, invgt class
#      the secondary, PC1 within that).
#
#   3. Three new panels (D1/D2/D3) — separate-first locally-harmonized
#      versions of A5/B5/C5 (the three mid-scale schemes). Each subgroup
#      gets its own polarity vector computed from ONLY its own samples'
#      Hom1/Hom2 means, and the dosages are flipped per sample accordingly.
#      These panels are only rendered for candidates the classifier labels
#      as COMPOSITE_INTERNAL, COMPOSITE_OVERLAP, or UNRESOLVED. For
#      SIMPLE_STRONG / SIMPLE_WEAK candidates, D1/D2/D3 are skipped (a
#      one-line stub page replaces them, since a global flip is already
#      correct).
#
#   4. Classifier catalog sidecar
#      The console dump at the end lists the assigned architecture class,
#      tier, reason, and the subgroup count + sizes.
#
# New CLI options (all optional — absence falls back to v2 behaviour):
#   --arch_class    path to candidate_architecture_class.tsv (STEP40 output)
#   --subgroup_file path to candidate_sample_subgroups.tsv   (STEP40 output)
#   --polarity_file path to candidate_marker_polarity.tsv    (STEP29 output,
#                   used for the "global flip" reference panel annotation)
#
# Produces an extra TSV:
#   ${OUT_DIR}/panel_architecture_meta.tsv
#       one row per config, columns: id, label, n_markers, ok,
#       architecture_class, tier, reason, simple_axis, boundary,
#       marker_family_k, monotonic_fraction, n_subgroups, n_markers_flipped
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
BEAGLE      <- get_arg("--beagle")
POS_FILE    <- get_arg("--pos")
PRECOMP     <- get_arg("--precomp")
CHROM       <- get_arg("--chrom")
SHELF_A     <- as.numeric(get_arg("--shelf_start_mb"))
SHELF_B     <- as.numeric(get_arg("--shelf_end_mb"))
BP1_MB      <- as.numeric(get_arg("--breakpoint1_mb", NA))
BP2_MB      <- as.numeric(get_arg("--breakpoint2_mb", NA))
SAMPLE_LIST <- get_arg("--sample_list")
INVGT_ASSIGN_FILE <- get_arg("--invgt_assign", NULL)
ARCH_FILE   <- get_arg("--arch_class", NULL)       # NEW
SUBG_FILE   <- get_arg("--subgroup_file", NULL)    # NEW
POLARITY_FILE <- get_arg("--polarity_file", NULL)  # NEW
OUT         <- get_arg("--out")
OUT_DIR     <- get_arg("--out_dir", NULL)
invisible(get_arg("--snp_stride", NULL))           # back-compat

stopifnot(file.exists(BEAGLE), file.exists(POS_FILE), file.exists(PRECOMP))
if (is.null(OUT_DIR)) OUT_DIR <- sub("\\.pdf$", "", OUT)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

has_CH <- suppressWarnings(require(ComplexHeatmap, quietly = TRUE))
has_circlize <- suppressWarnings(require(circlize, quietly = TRUE))
has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))
render_mode <- if (has_CH && has_circlize) "complexheatmap" else "ggplot"
message("[q08b] Render mode: ", render_mode)

# =============================================================================
# 1. Parse positions -> find shelf SNPs
# =============================================================================
message("[q08b] Loading positions: ", POS_FILE)
pos <- fread(POS_FILE, header = FALSE)
pos_col <- if (ncol(pos) >= 2) 2L else 1L
pos_vec <- as.integer(pos[[pos_col]])
shelf_mask <- pos_vec >= SHELF_A * 1e6 & pos_vec <= SHELF_B * 1e6
shelf_row_idx <- which(shelf_mask)
shelf_pos <- pos_vec[shelf_row_idx]
n_shelf <- length(shelf_row_idx)
if (n_shelf < 30) stop("[q08b] Only ", n_shelf, " SNPs in shelf — too few for any heatmap")
message("[q08b] Shelf SNPs: ", n_shelf)

# =============================================================================
# 2. Read sample master list, stream BEAGLE for full shelf dosage matrix
# =============================================================================
sample_master <- readLines(SAMPLE_LIST)
n_samples <- length(sample_master)
message("[q08b] Master samples: ", n_samples)

message("[q08b] Streaming BEAGLE for ", n_shelf, " shelf SNPs...")
dosage_mat <- matrix(NA_real_, nrow = n_shelf, ncol = n_samples)
rownames(dosage_mat) <- as.character(shelf_pos)
colnames(dosage_mat) <- sample_master

con <- gzfile(BEAGLE, "r")
header <- readLines(con, n = 1)
header_fields <- strsplit(header, "\t", fixed = TRUE)[[1]]
n_header_samples <- (length(header_fields) - 3) / 3
if (n_header_samples != n_samples) {
  message("[q08b] WARN: header says ", n_header_samples,
          " samples, master list says ", n_samples, ". Using header count.")
  n_samples <- as.integer(n_header_samples)
  dosage_mat <- matrix(NA_real_, nrow = n_shelf, ncol = n_samples)
  rownames(dosage_mat) <- as.character(shelf_pos)
  colnames(dosage_mat) <- if (length(sample_master) >= n_samples)
                           sample_master[seq_len(n_samples)]
                         else paste0("Ind", seq_len(n_samples) - 1L)
}

# env-backed hash (keeps the v2 fix from LG28 debugging)
target_env <- new.env(hash = TRUE, size = length(shelf_row_idx) * 2L)
for (i in seq_along(shelf_row_idx)) {
  assign(as.character(shelf_row_idx[i]), i, envir = target_env)
}

data_row_idx <- 0L; hits <- 0L; chunk_size <- 20000L
max_target <- max(shelf_row_idx)
repeat {
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0L) break
  for (ln in lines) {
    data_row_idx <- data_row_idx + 1L
    key <- as.character(data_row_idx)
    hi <- target_env[[key]]
    if (!is.null(hi)) {
      flds <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      gl_block <- as.numeric(flds[-(1:3)])
      gl_mat <- matrix(gl_block, nrow = n_samples, byrow = TRUE)
      dosage_mat[hi, ] <- gl_mat[, 2] + 2 * gl_mat[, 3]
      hits <- hits + 1L
      if (hits %% 2000L == 0L) message("  ", hits, " / ", n_shelf)
    }
    if (data_row_idx > max_target) break
  }
  if (data_row_idx > max_target) break
}
close(con)
message("[q08b] Extracted ", hits, " / ", n_shelf, " shelf rows")

snp_var <- apply(dosage_mat, 1L, var, na.rm = TRUE)
snp_var[!is.finite(snp_var)] <- 0

# =============================================================================
# 3. Sample ordering and class assignment (incl. optional subgroup)
# =============================================================================
sample_order <- NULL
sample_class <- NULL       # Hom1/Het/Hom2
sample_subgroup <- NULL    # NEW: integer subgroup id from STEP40

if (!is.null(INVGT_ASSIGN_FILE) && file.exists(INVGT_ASSIGN_FILE)) {
  message("[q08b] Using Q07 invgt assignments: ", INVGT_ASSIGN_FILE)
  inv <- fread(INVGT_ASSIGN_FILE)
  if ("sample_master" %in% names(inv) && "pc1_mean" %in% names(inv)) {
    inv[, invgt := factor(invgt, levels = c("Hom1", "Het", "Hom2"))]
    setorder(inv, invgt, pc1_mean)
    sample_order <- inv$sample_master
    sample_class <- setNames(as.character(inv$invgt), inv$sample_master)
  }
}
if (is.null(sample_order)) {
  message("[q08b] Falling back to PC1-mean ordering from precomp")
  pc <- readRDS(PRECOMP)
  dt <- as.data.table(pc$dt)
  dt[, mb := (start_bp + end_bp) / 2 / 1e6]
  shelf_mask_p <- dt$mb >= SHELF_A & dt$mb <= SHELF_B
  pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
  if (length(pc1_cols) > 0 && sum(shelf_mask_p) > 0) {
    pc1_means <- colMeans(as.matrix(dt[shelf_mask_p, ..pc1_cols]), na.rm = TRUE)
    pc_samples <- sub("^PC_1_", "", pc1_cols)
    if (all(grepl("^Ind[0-9]+$", pc_samples))) {
      ind_idx <- as.integer(sub("^Ind", "", pc_samples))
      if (length(sample_master) >= max(ind_idx) + 1L) {
        names(pc1_means) <- sample_master[ind_idx + 1L]
      } else names(pc1_means) <- pc_samples
    } else names(pc1_means) <- pc_samples
    sample_order <- names(sort(pc1_means))
  } else sample_order <- colnames(dosage_mat)
}
sample_order <- intersect(sample_order, colnames(dosage_mat))

# Resolve the shared ordering helper (shared with STEP36 and STEP41)
Q08B_SECONDARY <- if (exists("Q08B_SECONDARY_ORDER")) Q08B_SECONDARY_ORDER else "hclust"
if (!exists("order_samples_with_u_primary")) {
  .helper_candidates <- c(
    file.path(dirname(sys.frame(1)$ofile %||% "."), "within_group_ordering.R"),
    "within_group_ordering.R",
    file.path(dirname(sys.frame(1)$ofile %||% "."), "R", "within_group_ordering.R")
  )
  for (.h in .helper_candidates) {
    if (file.exists(.h)) { source(.h); break }
  }
}

# NEW: subgroup lookup from STEP40 output
if (!is.null(SUBG_FILE) && file.exists(SUBG_FILE)) {
  message("[q08b] Reading subgroup assignments: ", SUBG_FILE)
  subg <- tryCatch(fread(SUBG_FILE), error = function(e) NULL)
  if (!is.null(subg) && "sample" %in% names(subg) && "subgroup_id" %in% names(subg)) {
    sample_subgroup <- setNames(as.integer(subg$subgroup_id), subg$sample)

    # Try to extract u / v from the rotated PCA if available, and apply the
    # same u-primary + secondary-axis ordering used by STEP36 V3 / STEP41.
    # If the helper isn't available, fall back to the old subgroup→invgt→PC1 sort.
    u_by_sample <- NULL; v_by_sample <- NULL
    if (!is.null(CANDIDATE_DIR <- if (exists("CANDIDATE_DIR")) CANDIDATE_DIR else NULL)) {
      rot_file_q <- file.path(CANDIDATE_DIR, "candidate_pca_rotated.tsv")
      if (file.exists(rot_file_q)) {
        rot_q <- tryCatch(fread(rot_file_q), error = function(e) NULL)
        if (!is.null(rot_q) && "sample" %in% names(rot_q) && "u" %in% names(rot_q)) {
          u_by_sample <- setNames(as.numeric(rot_q$u), rot_q$sample)
          if ("v" %in% names(rot_q))
            v_by_sample <- setNames(as.numeric(rot_q$v), rot_q$sample)
        }
      }
    }

    if (exists("order_samples_with_u_primary") && !is.null(u_by_sample)) {
      # Use the helper. Primary grouping is still subgroup_id; within each
      # subgroup we apply the group→u→secondary ordering.
      message("[q08b]   Using u-primary ordering helper (secondary='",
              Q08B_SECONDARY, "')")
      new_order <- character(0)
      for (sg in sort(unique(sample_subgroup[sample_order]))) {
        sg_samples <- sample_order[sample_subgroup[sample_order] == sg]
        if (length(sg_samples) < 5) {
          new_order <- c(new_order, sg_samples); next
        }
        sg_groups <- sample_class[sg_samples]
        sg_groups[is.na(sg_groups)] <- "UNKNOWN"
        # Map Hom1/Het/Hom2 labels to the helper's expected HOMO_1/HET/HOMO_2
        sg_groups_mapped <- fifelse(sg_groups == "Hom1", "HOMO_1",
                               fifelse(sg_groups == "Hom2", "HOMO_2",
                                 fifelse(sg_groups == "Het", "HET", sg_groups)))
        sg_u <- u_by_sample[sg_samples]
        sg_v <- if (!is.null(v_by_sample)) v_by_sample[sg_samples] else NULL
        sg_dos <- dosage_mat[, sg_samples, drop = FALSE]
        res <- tryCatch(
          order_samples_with_u_primary(
            samples      = sg_samples,
            groups       = sg_groups_mapped,
            u            = sg_u,
            v            = sg_v,
            dosage       = sg_dos,
            secondary    = Q08B_SECONDARY,
            group_levels = c("HOMO_1", "HET", "HOMO_2"),
            min_group_n  = 5L,
            top_var_n    = 500L
          ), error = function(e) NULL)
        if (!is.null(res)) {
          new_order <- c(new_order, res$ordered_samples)
        } else {
          new_order <- c(new_order, sg_samples)
        }
      }
      sample_order <- intersect(new_order, colnames(dosage_mat))
      message("[q08b]   ", length(unique(sample_subgroup[sample_order])),
              " subgroups visible (u-primary order)")
    } else {
      # Fallback: old subgroup → invgt → PC1 rank order
      ord_dt <- data.table(
        sample = sample_order,
        subgroup = sample_subgroup[sample_order] %||% rep(1L, length(sample_order)),
        invgt    = sample_class[sample_order]   %||% rep(NA_character_, length(sample_order)),
        rank_pc1 = seq_along(sample_order)
      )
      ord_dt[is.na(subgroup), subgroup := 0L]
      setorder(ord_dt, subgroup, invgt, rank_pc1)
      sample_order <- ord_dt$sample
      message("[q08b]   ", length(unique(ord_dt$subgroup)),
              " subgroups visible (fallback PC1 order — helper unavailable)")
    }
  }
}
message("[q08b] Samples ordered: ", length(sample_order),
        if (!is.null(sample_class)) " (invgt)" else "",
        if (!is.null(sample_subgroup)) " (+subgroup)" else "")

# =============================================================================
# 3b. Optional architecture classifier metadata
# =============================================================================
arch_info <- NULL
if (!is.null(ARCH_FILE) && file.exists(ARCH_FILE)) {
  message("[q08b] Reading architecture class: ", ARCH_FILE)
  arch_dt <- tryCatch(fread(ARCH_FILE), error = function(e) NULL)
  if (!is.null(arch_dt) && nrow(arch_dt) >= 1) {
    # If multi-row file (catalog), select matching chrom+interval
    if ("chrom" %in% names(arch_dt)) {
      m <- arch_dt[chrom == CHROM]
      if (nrow(m) >= 1) arch_dt <- m[1]
    }
    arch_info <- as.list(arch_dt[1])
    message(sprintf("[q08b]   class = %s  tier = %s  reason = %s",
                    arch_info$architecture_class %||% "NA",
                    arch_info$tier %||% "NA",
                    arch_info$reason %||% ""))
  }
}

# Optional global polarity vector (from STEP29) — only used for D panel
# captions showing how much the local subgroup flips diverge from global.
global_flip <- NULL
if (!is.null(POLARITY_FILE) && file.exists(POLARITY_FILE)) {
  pol_dt <- tryCatch(fread(POLARITY_FILE), error = function(e) NULL)
  if (!is.null(pol_dt) && "final_flip_decision" %in% names(pol_dt)
      && "pos" %in% names(pol_dt)) {
    gf <- setNames(as.logical(pol_dt$final_flip_decision), as.character(pol_dt$pos))
    global_flip <- gf
    message("[q08b] Loaded global polarity: ", sum(gf, na.rm = TRUE), " flipped / ", length(gf))
  }
}

# =============================================================================
# 4. Scheme definitions (unchanged from v2)
# =============================================================================
scheme_A <- function(window_snps) {
  n_bins <- ceiling(n_shelf / window_snps)
  keep <- integer(0)
  for (b in seq_len(n_bins)) {
    lo <- (b - 1L) * window_snps + 1L
    hi <- min(n_shelf, b * window_snps)
    if (hi < lo) next
    seg_var <- snp_var[lo:hi]
    if (all(seg_var == 0)) next
    keep <- c(keep, lo + which.max(seg_var) - 1L)
  }
  sort(unique(keep))
}
scheme_B <- function(top_n) {
  ord <- order(snp_var, decreasing = TRUE)
  keep <- ord[seq_len(min(top_n, length(ord)))]
  keep <- keep[snp_var[keep] > 0]
  sort(unique(keep))
}
scheme_C <- function(bp_bin) {
  bin_id <- floor(shelf_pos / bp_bin)
  keep <- integer(0)
  for (b in unique(bin_id)) {
    idxs <- which(bin_id == b)
    if (!length(idxs)) next
    seg_var <- snp_var[idxs]
    if (all(seg_var == 0)) next
    keep <- c(keep, idxs[which.max(seg_var)])
  }
  sort(unique(keep))
}

configs <- list(
  list(id = "A1",  label = "1 top-var / 100-SNP window",  idx_fn = function() scheme_A(100),  kind = "raw"),
  list(id = "A5",  label = "1 top-var / 500-SNP window",  idx_fn = function() scheme_A(500),  kind = "raw"),
  list(id = "A10", label = "1 top-var / 1000-SNP window", idx_fn = function() scheme_A(1000), kind = "raw"),
  list(id = "B1",  label = "top 30 var markers",          idx_fn = function() scheme_B(30),   kind = "raw"),
  list(id = "B5",  label = "top 100 var markers",         idx_fn = function() scheme_B(100),  kind = "raw"),
  list(id = "B10", label = "top 500 var markers",         idx_fn = function() scheme_B(500),  kind = "raw"),
  list(id = "C1",  label = "1 top-var / 5 kb",            idx_fn = function() scheme_C(5e3),  kind = "raw"),
  list(id = "C5",  label = "1 top-var / 25 kb",           idx_fn = function() scheme_C(25e3), kind = "raw"),
  list(id = "C10", label = "1 top-var / 100 kb",          idx_fn = function() scheme_C(100e3),kind = "raw"),
  # NEW D panels: separate-first locally-harmonized versions of A5, B5, C5
  list(id = "D1",  label = "SEP-FIRST 1 top-var / 500-SNP window", idx_fn = function() scheme_A(500), kind = "sepfirst"),
  list(id = "D2",  label = "SEP-FIRST top 100 var markers",        idx_fn = function() scheme_B(100), kind = "sepfirst"),
  list(id = "D3",  label = "SEP-FIRST 1 top-var / 25 kb",          idx_fn = function() scheme_C(25e3),kind = "sepfirst")
)

# =============================================================================
# 4b. Separate-first local harmonization
# =============================================================================
# For each sample subgroup independently, compute a per-marker flip vector
# using only that subgroup's Hom1/Hom2 means (falls back to local PC1 sign
# if one karyotype is absent in the subgroup). Returns a samples × markers
# logical flip matrix + per-subgroup diagnostics.
#
# IMPORTANT: this is the "correct second" step. The "separate first" has
# already been done by STEP40 and stored in SUBG_FILE.
# =============================================================================
compute_sep_first_flip <- function(X_mark_samp, sample_ids, marker_positions) {
  n_m <- nrow(X_mark_samp); n_s <- ncol(X_mark_samp)
  flip_by_sample <- matrix(FALSE, nrow = n_m, ncol = n_s)
  diag_list <- list()

  if (is.null(sample_subgroup)) {
    # No subgroup info at all — compute one global flip and apply uniformly.
    h1 <- which(sample_class[sample_ids] == "Hom1")
    h2 <- which(sample_class[sample_ids] == "Hom2")
    if (length(h1) >= 3 && length(h2) >= 3) {
      mh1 <- rowMeans(X_mark_samp[, h1, drop = FALSE], na.rm = TRUE)
      mh2 <- rowMeans(X_mark_samp[, h2, drop = FALSE], na.rm = TRUE)
      delta <- mh2 - mh1; delta[!is.finite(delta)] <- 0
      dom <- sign(median(delta[abs(delta) > 0.1], na.rm = TRUE))
      if (!is.finite(dom) || dom == 0) dom <- 1
      flip_v <- (sign(delta) != dom) & (abs(delta) > 0.1)
      for (j in seq_len(n_s)) flip_by_sample[, j] <- flip_v
      diag_list[["all"]] <- list(n_samples = n_s, n_flipped = sum(flip_v),
                                  method = "global_group_difference")
    }
    return(list(flip_by_sample = flip_by_sample, diag = diag_list))
  }

  subg_ids <- sample_subgroup[sample_ids]
  subg_ids[is.na(subg_ids)] <- 0L

  for (sg in unique(subg_ids)) {
    idx_in_sg <- which(subg_ids == sg)
    if (length(idx_in_sg) < 5) {
      diag_list[[as.character(sg)]] <- list(n_samples = length(idx_in_sg),
                                             n_flipped = 0,
                                             method = "skip_too_few")
      next
    }
    X_sg <- X_mark_samp[, idx_in_sg, drop = FALSE]
    cls_sg <- sample_class[sample_ids[idx_in_sg]]
    h1 <- which(cls_sg == "Hom1"); h2 <- which(cls_sg == "Hom2")
    flip_v <- rep(FALSE, n_m); method <- "none"
    if (length(h1) >= 3 && length(h2) >= 3) {
      mh1 <- rowMeans(X_sg[, h1, drop = FALSE], na.rm = TRUE)
      mh2 <- rowMeans(X_sg[, h2, drop = FALSE], na.rm = TRUE)
      delta <- mh2 - mh1; delta[!is.finite(delta)] <- 0
      dom <- sign(median(delta[abs(delta) > 0.1], na.rm = TRUE))
      if (!is.finite(dom) || dom == 0) dom <- 1
      flip_v <- (sign(delta) != dom) & (abs(delta) > 0.1)
      method <- "local_group_difference"
    } else {
      # PC1-sign fallback within subgroup
      M <- X_sg
      for (j in seq_len(ncol(M))) {
        bad <- !is.finite(M[, j])
        if (any(bad)) M[bad, j] <- mean(M[, j], na.rm = TRUE)
      }
      M[!is.finite(M)] <- 0
      v <- apply(M, 1, var, na.rm = TRUE)
      kk <- which(is.finite(v) & v > 1e-10)
      if (length(kk) >= 5) {
        pc <- tryCatch(prcomp(t(M[kk, , drop = FALSE]), center = TRUE,
                              scale. = FALSE, rank. = 1),
                       error = function(e) NULL)
        if (!is.null(pc)) {
          loadings <- numeric(n_m); loadings[kk] <- pc$rotation[, 1]
          big <- abs(loadings) > quantile(abs(loadings), 0.5)
          dom <- sign(median(loadings[big]))
          if (!is.finite(dom) || dom == 0) dom <- 1
          flip_v <- (sign(loadings) != dom) & (abs(loadings) > 0)
          method <- "local_pc1_sign"
        }
      }
    }
    for (j in idx_in_sg) flip_by_sample[, j] <- flip_v
    diag_list[[as.character(sg)]] <- list(n_samples = length(idx_in_sg),
                                           n_flipped = sum(flip_v),
                                           method = method)
  }
  list(flip_by_sample = flip_by_sample, diag = diag_list)
}

# =============================================================================
# 5. Render helpers (extended with subgroup annotation + caption)
# =============================================================================
dosage_cols <- c("#2166AC", "#F7F7F7", "#B2182B")
class_cols  <- c("Hom1" = "#1f4e79", "Het" = "#c0504d", "Hom2" = "#2c7a39")
subg_pal    <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                  "#E6AB02", "#A6761D", "#666666")

build_title <- function(base_title, config_kind) {
  suffix <- ""
  if (!is.null(arch_info)) {
    cls <- arch_info$architecture_class %||% "?"
    tr  <- arch_info$tier %||% "?"
    suffix <- sprintf("  [%s-%s]", cls, tr)
  }
  paste0(base_title, suffix)
}
build_caption <- function(config_kind, n_markers_flipped = NA_integer_) {
  if (is.null(arch_info)) return("")
  parts <- c()
  for (fn in c("simple_axis_score", "boundary_consistency",
               "monotonic_fraction", "marker_family_k")) {
    v <- arch_info[[fn]]
    if (!is.null(v) && is.finite(suppressWarnings(as.numeric(v)))) {
      parts <- c(parts, sprintf("%s=%s", sub("_score|_fraction|_consistency", "", fn),
                                format(round(as.numeric(v), 3))))
    }
  }
  cap <- paste(parts, collapse = "  ")
  if (config_kind == "sepfirst" && is.finite(n_markers_flipped)) {
    cap <- paste0(cap, sprintf("  |  n_flipped(max_subgrp)=%d", n_markers_flipped))
  }
  cap
}

render_ch <- function(X, positions, title, caption, pdf_path, png_path) {
  col_fun <- circlize::colorRamp2(c(0, 1, 2), dosage_cols)

  ha_left <- NULL
  left_annots <- list()
  left_cols <- list()
  if (!is.null(sample_class)) {
    class_v <- sample_class[rownames(X)]; class_v[is.na(class_v)] <- "unknown"
    col_map <- c(class_cols, "unknown" = "#CCCCCC")
    col_map <- col_map[intersect(names(col_map), unique(class_v))]
    left_annots$invgt <- class_v
    left_cols$invgt <- col_map
  }
  if (!is.null(sample_subgroup)) {
    sg_v <- sample_subgroup[rownames(X)]; sg_v[is.na(sg_v)] <- 0L
    sg_v_chr <- as.character(sg_v)
    uniq <- sort(unique(sg_v))
    col_map <- setNames(subg_pal[((seq_along(uniq) - 1L) %% length(subg_pal)) + 1L],
                         as.character(uniq))
    left_annots$subgroup <- sg_v_chr
    left_cols$subgroup <- col_map
  }
  if (length(left_annots) > 0) {
    ha_left <- ComplexHeatmap::rowAnnotation(
      df = as.data.frame(left_annots, stringsAsFactors = FALSE),
      col = left_cols,
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 7),
      annotation_width = grid::unit(rep(4, length(left_annots)), "mm")
    )
  }

  pos_mb <- positions / 1e6
  annot_list <- list(
    pos_mb = ComplexHeatmap::anno_points(
      pos_mb, pch = 16, size = grid::unit(1, "mm"),
      gp = grid::gpar(col = "#555e69"),
      axis_param = list(gp = grid::gpar(fontsize = 6))
    )
  )
  bp_keep <- c(); bp_labels <- character()
  for (bp in c(BP1_MB, BP2_MB)) {
    if (is.finite(bp)) {
      j <- which.min(abs(pos_mb - bp))
      if (length(j) == 1L) {
        bp_keep <- c(bp_keep, j)
        bp_labels <- c(bp_labels, sprintf("BP %.3f Mb", bp))
      }
    }
  }
  if (length(bp_keep) > 0) {
    annot_list$breakpoint <- ComplexHeatmap::anno_mark(
      at = bp_keep, labels = bp_labels, which = "column", side = "top",
      labels_gp = grid::gpar(fontsize = 7, col = "#8a2b30"),
      lines_gp  = grid::gpar(col = "#e0555c", lwd = 1)
    )
  }
  ha_top <- do.call(ComplexHeatmap::columnAnnotation,
                    c(annot_list,
                      list(height = grid::unit(6, "mm"),
                           annotation_name_gp = grid::gpar(fontsize = 7))))

  n_cell <- nrow(X) * ncol(X)
  w <- max(6, min(16, 3 + ncol(X) * 0.06))
  h <- max(5, min(18, 3 + nrow(X) * 0.035))

  ht <- ComplexHeatmap::Heatmap(
    X, name = "Dosage", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    left_annotation = ha_left,
    top_annotation = ha_top,
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    use_raster = n_cell > 50000,
    raster_quality = 3,
    na_col = "#ffffff"
  )

  draw_with_caption <- function() {
    ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 2, 2, 2), "mm"))
    if (nzchar(caption)) {
      grid::grid.text(caption, x = 0.01, y = 0.01, just = c("left", "bottom"),
                      gp = grid::gpar(fontsize = 7, col = "#444444"))
    }
  }
  tryCatch({
    pdf(pdf_path, width = w, height = h); draw_with_caption(); dev.off()
    png(png_path, width = w, height = h, units = "in", res = 150)
    draw_with_caption(); dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("  [WARN] CH render failed for ", basename(pdf_path), ": ",
            conditionMessage(e))
  })
}

render_ggplot <- function(X, positions, title, caption, pdf_path, png_path) {
  df <- as.data.table(reshape2::melt(X, varnames = c("sample", "snp"),
                                      value.name = "dosage"))
  df[, sample := factor(sample, levels = rownames(X))]
  df[, snp := factor(snp, levels = colnames(X))]
  df[, gt := fcase(
    is.na(dosage), "NA",
    dosage < 0.5, "0/0",
    dosage < 1.5, "0/1",
    default = "1/1"
  )]
  gt_pal <- c("0/0" = dosage_cols[1], "0/1" = "#d0d0d0",
              "1/1" = dosage_cols[3], "NA" = "#ffffff")
  df[, gt := factor(gt, levels = names(gt_pal))]

  p <- ggplot(df, aes(snp, sample, fill = gt)) +
    geom_raster() +
    scale_fill_manual(values = gt_pal, name = "gt", na.value = "#ffffff") +
    labs(title = title,
         subtitle = if (nzchar(caption)) caption else NULL,
         x = NULL, y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.key.width = unit(0.25, "cm"),
          legend.key.height = unit(0.35, "cm"),
          plot.subtitle = element_text(size = 7, colour = "#444444"))

  # Subgroup+invgt strip
  strip_plots <- list()
  if (!is.null(sample_class)) {
    class_v <- sample_class[rownames(X)]
    strip_df <- data.table(sample = factor(rownames(X), levels = rownames(X)),
                            class = class_v)
    strip_plots$invgt <- ggplot(strip_df, aes(x = 1, y = sample, fill = class)) +
      geom_raster() + scale_fill_manual(values = class_cols, na.value = "#CCCCCC") +
      theme_void() + theme(legend.position = "none")
  }
  if (!is.null(sample_subgroup)) {
    sg_v <- sample_subgroup[rownames(X)]
    sg_df <- data.table(sample = factor(rownames(X), levels = rownames(X)),
                        subgroup = factor(sg_v))
    uniq <- sort(unique(as.integer(as.character(sg_df$subgroup))))
    col_map <- setNames(subg_pal[((seq_along(uniq) - 1L) %% length(subg_pal)) + 1L],
                         as.character(uniq))
    strip_plots$subgroup <- ggplot(sg_df, aes(x = 1, y = sample, fill = subgroup)) +
      geom_raster() + scale_fill_manual(values = col_map, na.value = "#CCCCCC") +
      theme_void() + theme(legend.position = "none")
  }
  if (length(strip_plots) > 0 && has_patchwork) {
    strips <- Reduce(`|`, strip_plots)
    p <- (strips | p) + patchwork::plot_layout(widths = c(0.03 * length(strip_plots), 1))
  }

  w <- max(6, min(16, 3 + ncol(X) * 0.06))
  h <- max(5, min(18, 3 + nrow(X) * 0.035))
  ggsave(pdf_path, p, width = w, height = h, device = cairo_pdf)
  ggsave(png_path, p, width = w, height = h, dpi = 150)
}

render_empty <- function(title, pdf_path, png_path, reason) {
  p <- ggplot() + labs(title = title, subtitle = paste0("(", reason, ")")) +
    theme_void() + theme(plot.title = element_text(size = 9))
  ggsave(pdf_path, p, width = 6, height = 5, device = cairo_pdf)
  ggsave(png_path, p, width = 6, height = 5, dpi = 100)
}

# =============================================================================
# 6. Generate panels (A1..C10 raw + D1..D3 separate-first)
# =============================================================================
panel_files <- character(length(configs))
meta_rows   <- list()

# Determine if this candidate needs D-panels
run_sepfirst <- TRUE
if (!is.null(arch_info)) {
  cls <- arch_info$architecture_class %||% "UNRESOLVED"
  run_sepfirst <- cls %in% c("COMPOSITE_INTERNAL", "COMPOSITE_OVERLAP", "UNRESOLVED")
  message("[q08b] D-panels (separate-first) enabled: ", run_sepfirst,
          "  (class=", cls, ")")
}

for (ci in seq_along(configs)) {
  cf <- configs[[ci]]
  message(sprintf("[q08b] Rendering %s: %s (%s)", cf$id, cf$label, cf$kind))

  pdf_path <- file.path(OUT_DIR, paste0("panel_", cf$id, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0("panel_", cf$id, ".png"))
  panel_files[ci] <- pdf_path

  meta <- list(id = cf$id, label = cf$label, kind = cf$kind,
               n_markers = 0L, ok = FALSE,
               architecture_class = arch_info$architecture_class %||% NA_character_,
               tier = arch_info$tier %||% NA_character_,
               reason = arch_info$reason %||% NA_character_,
               n_subgroups = if (!is.null(sample_subgroup))
                 length(unique(sample_subgroup[sample_order])) else 1L,
               n_markers_flipped_max = NA_integer_)

  if (cf$kind == "sepfirst" && !run_sepfirst) {
    render_empty(build_title(paste0(cf$id, ": ", cf$label), cf$kind),
                  pdf_path, png_path,
                  "skipped: architecture class is simple (no subgroup decomposition needed)")
    meta$ok <- FALSE
    meta_rows[[length(meta_rows) + 1]] <- meta
    next
  }

  idx <- tryCatch(cf$idx_fn(), error = function(e) integer(0))
  n_markers <- length(idx); meta$n_markers <- n_markers

  if (n_markers < 3) {
    render_empty(build_title(paste0(cf$id, ": ", cf$label), cf$kind),
                  pdf_path, png_path,
                  paste0("only ", n_markers, " markers available"))
    meta_rows[[length(meta_rows) + 1]] <- meta
    next
  }

  X_shelf <- dosage_mat[idx, , drop = FALSE]              # SNPs × samples
  panel_positions <- shelf_pos[idx]

  # For sepfirst panels: compute per-sample flip matrix on this panel's markers
  if (cf$kind == "sepfirst") {
    sf <- compute_sep_first_flip(X_shelf, colnames(X_shelf), panel_positions)
    max_flipped <- if (length(sf$diag) > 0)
      max(vapply(sf$diag, function(d) as.integer(d$n_flipped %||% 0L), integer(1)))
      else 0L
    meta$n_markers_flipped_max <- as.integer(max_flipped)

    # Apply flip per sample
    X_harm <- X_shelf
    fbs <- sf$flip_by_sample
    for (s_j in seq_len(ncol(X_harm))) {
      fv <- fbs[, s_j]
      if (any(fv)) X_harm[fv, s_j] <- 2 - X_harm[fv, s_j]
    }
    X_display <- t(X_harm)[sample_order, , drop = FALSE]
  } else {
    X_display <- t(X_shelf)[sample_order, , drop = FALSE]
  }
  colnames(X_display) <- as.character(panel_positions)

  base_title <- sprintf("%s: %s  (%d SNPs × %d samples)",
                        cf$id, cf$label, n_markers, nrow(X_display))
  title   <- build_title(base_title, cf$kind)
  caption <- build_caption(cf$kind, meta$n_markers_flipped_max)

  if (render_mode == "complexheatmap") {
    render_ch(X_display, panel_positions, title, caption, pdf_path, png_path)
  } else if (requireNamespace("reshape2", quietly = TRUE)) {
    render_ggplot(X_display, panel_positions, title, caption, pdf_path, png_path)
  } else {
    render_empty(title, pdf_path, png_path, "reshape2 missing for ggplot path")
  }
  meta$ok <- TRUE
  meta_rows[[length(meta_rows) + 1]] <- meta
}

meta_df <- rbindlist(lapply(meta_rows, function(m) as.data.table(m)), fill = TRUE)

# =============================================================================
# 7. Combined PDF (same strategy as v2)
# =============================================================================
tryCatch({
  qpdf_bin <- Sys.which("qpdf")
  if (nzchar(qpdf_bin)) {
    existing <- panel_files[file.exists(panel_files)]
    if (length(existing) > 1) {
      cmd <- sprintf("%s --empty --pages %s -- %s",
                     qpdf_bin, paste(shQuote(existing), collapse = " "),
                     shQuote(OUT))
      ret <- system(cmd, intern = FALSE)
      if (ret == 0 && file.exists(OUT))
        message("[q08b] Combined PDF via qpdf -> ", OUT)
      else stop("qpdf returned ", ret)
    } else stop("no individual panels found")
  } else stop("qpdf not found")
}, error = function(e) {
  message("[q08b] qpdf unavailable (", conditionMessage(e),
          "); writing index PDF")
  pdf(OUT, width = 11, height = 8, onefile = TRUE)
  on.exit(try(dev.off(), silent = TRUE))
  for (pf in panel_files) {
    if (!file.exists(pf)) next
    grid::grid.newpage()
    grid::grid.text(paste0("see individual file: ", basename(pf)),
                    gp = grid::gpar(fontsize = 12))
  }
  on.exit(NULL); try(dev.off(), silent = TRUE)
  message("[q08b] Index composite -> ", OUT)
})

# =============================================================================
# 8. Summary TSVs + console dump
# =============================================================================
summary_tsv <- file.path(OUT_DIR, "panel_summary.tsv")
fwrite(meta_df[, .(id, label, n_markers, ok)], summary_tsv, sep = "\t", quote = FALSE)

meta_tsv <- file.path(OUT_DIR, "panel_architecture_meta.tsv")
fwrite(meta_df, meta_tsv, sep = "\t", quote = FALSE)

cat("\n========= Q08b SHELF HEATMAPS (architecture-aware) =========\n")
cat(sprintf("Chromosome : %s\n", CHROM))
cat(sprintf("Shelf      : %.2f - %.2f Mb (%d SNPs total)\n", SHELF_A, SHELF_B, n_shelf))
cat(sprintf("Samples    : %d  %s%s\n", length(sample_order),
            if (!is.null(sample_class)) "(invgt)" else "(unannotated)",
            if (!is.null(sample_subgroup)) " (+subgroup)" else ""))
cat(sprintf("Render     : %s\n", render_mode))
if (!is.null(arch_info)) {
  cat(sprintf("Architecture: %s-%s\n",
              arch_info$architecture_class %||% "?",
              arch_info$tier %||% "?"))
  if (!is.null(arch_info$reason)) cat(sprintf("  reason   : %s\n", arch_info$reason))
  cat(sprintf("  D-panels : %s\n", if (run_sepfirst) "rendered" else "skipped (simple class)"))
}
cat("Panels:\n")
for (i in seq_len(nrow(meta_df))) {
  cat(sprintf("  %-5s [%-8s] %-38s n=%5d  %s\n",
              meta_df$id[i], meta_df$kind[i],
              substr(meta_df$label[i], 1, 38),
              meta_df$n_markers[i],
              if (meta_df$ok[i]) "OK" else "EMPTY"))
}
cat(sprintf("Summary    : %s\n", summary_tsv))
cat(sprintf("Arch meta  : %s\n", meta_tsv))
cat(sprintf("Combined   : %s\n", OUT))
cat(sprintf("Panels dir : %s\n", OUT_DIR))
cat("============================================================\n")
message("[q08b] DONE")
