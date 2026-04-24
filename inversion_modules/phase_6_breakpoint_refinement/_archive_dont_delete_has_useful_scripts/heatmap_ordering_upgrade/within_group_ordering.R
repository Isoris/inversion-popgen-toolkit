# =============================================================================
# within_group_ordering.R   v2 — u-axis aware
#
# Revised after reading the existing STEP20b / STEP21 anchor geometry pipeline:
#
#   - STEP21 already computes a rotated (u, v) coordinate system per candidate,
#     where u is the projection onto the Hom1 → Hom2 centroid axis. This
#     ordering is biologically meaningful: within Hom1, u measures how far
#     into the Hom1 arrangement a sample sits, placing samples nearer to Het
#     at one end and samples deep in Hom1 at the other.
#
#   - For simple inversions this is already the correct ordering. Keeping it
#     makes the three karyotype bands stack cleanly and preserves the visual
#     consistency with all previously published heatmaps.
#
#   - For composite inversions like LG28, u captures only the main inversion
#     axis. Any secondary substructure inside Hom1 (sub-haplotypes, founder
#     groups, gene-conversion recipients) is perpendicular to u and therefore
#     invisible to u-sorting alone. This is what makes LG28 look scrambled.
#
# Fix: two-level ordering
#
#   Level 1: group → u   (same as before, preserves the main V-geometry)
#   Level 2: within each group, tie-break on a secondary axis that captures
#            residual variance after the u signal is removed.
#
# Secondary methods:
#   "v"             → just use v from the rotated PCA (fastest, no clustering)
#   "residual_pc1"  → PCA on the group's dosage matrix after regressing out
#                      u; order by PC1 of the residual. Reveals sub-haplotypes.
#   "hclust"        → hierarchical clustering on the group's dosage profile
#                      (correlation distance, average linkage). Dendrogram
#                      leaf order as the score. Most flexible.
#   "error_profile" → cluster on sign-of-deviation from group median.
#                      Directly implements the user's intuition that samples
#                      sharing an error pattern should end up adjacent.
#   "none"          → purely u-based ordering (old behaviour).
#
# u always dominates. The secondary score only breaks ties, so on simple
# candidates the heatmap is nearly identical to what STEP28/35 produce today.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

order_samples_with_u_primary <- function(samples, groups, u,
                                           v = NULL,
                                           dosage = NULL,
                                           secondary = c("none", "v", "residual_pc1",
                                                         "hclust", "error_profile"),
                                           group_levels = c("HOMO_1", "HET", "HOMO_2"),
                                           min_group_n = 5L,
                                           top_var_n = 500L) {
  secondary <- match.arg(secondary)
  stopifnot(length(samples) == length(groups), length(samples) == length(u))
  if (!is.null(v)) stopifnot(length(v) == length(samples))

  # Restrict dosage to these samples (in order)
  if (!is.null(dosage) && !is.null(colnames(dosage))) {
    keep <- samples %in% colnames(dosage)
    if (!all(keep)) {
      samples <- samples[keep]; groups <- groups[keep]
      u <- u[keep]; if (!is.null(v)) v <- v[keep]
    }
    dosage <- dosage[, samples, drop = FALSE]
  }

  need_dosage <- secondary %in% c("residual_pc1", "hclust", "error_profile")
  D_top <- NULL
  if (need_dosage) {
    if (is.null(dosage) || nrow(dosage) < 10) {
      message("[within_group_ordering] secondary='", secondary,
              "' needs dosage; falling back to 'v' / 'none'")
      secondary <- if (!is.null(v)) "v" else "none"
      need_dosage <- FALSE
    } else {
      mv <- apply(dosage, 1L, var, na.rm = TRUE)
      mv[!is.finite(mv)] <- 0
      top_idx <- order(mv, decreasing = TRUE)[seq_len(min(top_var_n, nrow(dosage)))]
      top_idx <- top_idx[mv[top_idx] > 0]
      if (length(top_idx) < 10L) {
        secondary <- if (!is.null(v)) "v" else "none"
        need_dosage <- FALSE
      } else {
        D_top <- dosage[top_idx, , drop = FALSE]
        storage.mode(D_top) <- "double"
      }
    }
  }

  ordered_samples <- character(0)
  audit_rows <- list()

  for (g in group_levels) {
    idx_g <- which(groups == g)
    if (!length(idx_g)) next
    g_samples <- samples[idx_g]
    g_u <- u[idx_g]

    sec_score <- rep(NA_real_, length(g_samples))
    method_used <- "u_only"

    if (length(g_samples) >= min_group_n) {
      if (secondary == "v" && !is.null(v)) {
        sec_score <- v[idx_g]
        method_used <- "u_then_v"

      } else if (secondary == "residual_pc1" && !is.null(D_top)) {
        Dg <- D_top[, g_samples, drop = FALSE]
        for (i in seq_len(nrow(Dg))) {
          bad <- !is.finite(Dg[i, ])
          if (any(bad)) Dg[i, bad] <- mean(Dg[i, !bad])
        }
        Dg[!is.finite(Dg)] <- 0
        # Regress out u from each marker row
        u_c <- g_u - mean(g_u)
        u_ss <- sum(u_c^2)
        if (u_ss > 1e-10) {
          for (i in seq_len(nrow(Dg))) {
            row_c <- Dg[i, ] - mean(Dg[i, ])
            beta <- sum(row_c * u_c) / u_ss
            Dg[i, ] <- row_c - beta * u_c
          }
        }
        sv <- apply(Dg, 1L, var)
        kk <- which(is.finite(sv) & sv > 0)
        if (length(kk) >= 2L) {
          pc <- tryCatch(
            prcomp(t(Dg[kk, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 1),
            error = function(e) NULL
          )
          if (!is.null(pc)) {
            sec_score <- as.numeric(pc$x[, 1])
            method_used <- "u_then_residual_pc1"
          }
        }

      } else if (secondary == "hclust" && !is.null(D_top)) {
        Dg <- D_top[, g_samples, drop = FALSE]
        for (i in seq_len(nrow(Dg))) {
          bad <- !is.finite(Dg[i, ])
          if (any(bad)) Dg[i, bad] <- mean(Dg[i, !bad])
        }
        Dg[!is.finite(Dg)] <- 0
        sv <- apply(Dg, 2L, var)
        d <- if (any(!is.finite(sv) | sv == 0)) {
          dist(t(Dg), method = "euclidean")
        } else {
          as.dist(1 - suppressWarnings(cor(Dg, use = "pairwise.complete.obs")))
        }
        hc <- tryCatch(hclust(d, method = "average"), error = function(e) NULL)
        if (!is.null(hc)) {
          leaf_rank <- seq_along(hc$order)
          names(leaf_rank) <- g_samples[hc$order]
          sec_score <- leaf_rank[g_samples]
          method_used <- "u_then_hclust_leaf"
        }

      } else if (secondary == "error_profile" && !is.null(D_top)) {
        Dg <- D_top[, g_samples, drop = FALSE]
        for (i in seq_len(nrow(Dg))) {
          bad <- !is.finite(Dg[i, ])
          if (any(bad)) Dg[i, bad] <- mean(Dg[i, !bad])
        }
        Dg[!is.finite(Dg)] <- 0
        grp_med <- apply(Dg, 1L, median)
        dev_sign <- sign(Dg - grp_med)
        dev_sign[!is.finite(dev_sign)] <- 0
        d <- dist(t(dev_sign), method = "manhattan")
        hc <- tryCatch(hclust(d, method = "average"), error = function(e) NULL)
        if (!is.null(hc)) {
          leaf_rank <- seq_along(hc$order)
          names(leaf_rank) <- g_samples[hc$order]
          sec_score <- leaf_rank[g_samples]
          method_used <- "u_then_error_profile"
        }
      }
    } else {
      method_used <- "group_too_small_u_only"
    }

    # Primary: u. Secondary: sec_score (if available).
    ord_dt <- data.table(sample = g_samples, u = g_u, sec = sec_score)
    if (all(is.na(ord_dt$sec))) {
      setorder(ord_dt, u)
    } else {
      setorder(ord_dt, u, sec)
    }

    ordered_samples <- c(ordered_samples, ord_dt$sample)
    audit_rows[[length(audit_rows) + 1]] <- data.table(
      sample = ord_dt$sample, group = g,
      u = round(ord_dt$u, 6),
      secondary_score = round(ord_dt$sec, 6),
      original_u_rank = order(order(g_u))[match(ord_dt$sample, g_samples)],
      final_rank = seq_along(ord_dt$sample) +
                   (length(ordered_samples) - nrow(ord_dt)),
      method_used = method_used
    )
  }

  leftover <- setdiff(samples, ordered_samples)
  if (length(leftover)) {
    ordered_samples <- c(ordered_samples, leftover)
    audit_rows[[length(audit_rows) + 1]] <- data.table(
      sample = leftover,
      group = groups[match(leftover, samples)],
      u = u[match(leftover, samples)],
      secondary_score = NA_real_,
      original_u_rank = NA_integer_,
      final_rank = seq_along(leftover) +
                   (length(ordered_samples) - length(leftover)),
      method_used = "ungrouped_passthrough"
    )
  }

  list(
    ordered_samples = ordered_samples,
    within_group_audit = rbindlist(audit_rows)
  )
}
