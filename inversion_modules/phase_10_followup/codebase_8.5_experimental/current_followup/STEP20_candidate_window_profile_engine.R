#!/usr/bin/env Rscript

# =============================================================================
# STEP20_candidate_window_profile_engine.R  (v4)
#
# Core engine: sample-structure-first inversion candidate analysis.
#
# DESIGN PRINCIPLES:
#   1. Returns to marker-level dosage after STEP09 PCA compression.
#      STEP09 compresses each window into covariance/eigen summaries.
#      This script explicitly bypasses that compression by loading the
#      raw dosage matrix from STEP08 for each candidate region.
#
#   2. Raw ordered per-sample marker vectors are FIRST-CLASS OBJECTS.
#      They are preserved in three encodings (minor, major, combined 012)
#      before any PCA, MDS, polarity, or heatmap transformation.
#      All profile distances and reference scores are computed from
#      these raw vectors, not from PCA summaries.
#
#   3. Clair3 is a parallel optional enrichment layer, not a dependency.
#      MODE 1 (dosage-only) works now from discovery-side markers alone.
#      MODE 2 (Clair3) plugs additional signatures into the same framework
#      without changing downstream logic or output formats.
#
#   4. Three encodings are maintained independently:
#      - MINOR: dosage as expected minor allele count (0→2)
#      - MAJOR: 2 - minor dosage (expected major allele count)
#      - COMBINED_012: discretized 0/1/2 from minor dosage
#      Even though minor and major are scalar inverses, they behave
#      differently operationally after discretization, polarity
#      harmonization, and profile distance computation.
#
# OUTPUTS PER CANDIDATE:
#
#   Layer 1 — Window Profile:
#     candidate_retained_windows.tsv
#     candidate_window_profile_table.tsv
#
#   Layer 2 — Raw Ordered Marker Vectors (first-class, not collapsed):
#     candidate_raw_minor_vectors.tsv.gz       (markers × samples, continuous)
#     candidate_raw_major_vectors.tsv.gz       (markers × samples, continuous)
#     candidate_012_marker_matrix.tsv.gz       (markers × samples, discrete 0/1/2)
#     candidate_sample_profile_strings.tsv     (per-sample: 012 string, AB string)
#
#   Layer 3 — Profile Distance (computed from raw vectors):
#     candidate_pairwise_profile_distance_minor.tsv.gz
#     candidate_pairwise_profile_distance_major.tsv.gz
#     candidate_pairwise_profile_distance_012.tsv.gz
#
#   Layer 4 — Reference/Consensus Comparison:
#     candidate_reference_profile_scores.tsv
#     candidate_best_profile_assignment.tsv
#
#   Layer 5 — Window Grouping by Induced Sample Structure:
#     candidate_window_similarity_matrix.tsv
#     candidate_window_group_assignments.tsv
#     candidate_window_group_summary.tsv
#
#   Layer 6 — Signature Summary:
#     candidate_marker_informativeness.tsv
#
# Usage:
#   Rscript STEP20_candidate_window_profile_engine.R <config.R> [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- if (length(args) >= 2 && args[2] != "all") as.integer(args[2]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

# ── Clair3 availability check ───────────────────────────────────────────────
CLAIR3_ROOT <- file.path(PROJECT_ROOT, "clair3_indel_discovery", "postprocess_results")
clair3_available <- dir.exists(CLAIR3_ROOT)
message("[INFO] Signature mode: ", if (clair3_available) "MODE 2 (dosage + Clair3)" else "MODE 1 (dosage-only)")

# ── Read inputs ──────────────────────────────────────────────────────────────
cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]
message("[INFO] Processing ", nrow(cand), " candidates")

# Load STEP09 for window coordinates
step09_rds <- list.files(STEP09_DIR, pattern = "\\.window_pca\\.rds$",
                         full.names = TRUE, recursive = TRUE)
step09_by_chr <- list()
for (f in step09_rds) {
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.null(obj)) step09_by_chr[[as.character(obj$chrom)]] <- obj
}

# ═══════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

discretize_012 <- function(x, t1 = 0.5, t2 = 1.5) {
  # Discretize continuous dosage to 0/1/2
  out <- rep(NA_integer_, length(x))
  out[!is.na(x) & x < t1] <- 0L
  out[!is.na(x) & x >= t1 & x <= t2] <- 1L
  out[!is.na(x) & x > t2] <- 2L
  out
}

encode_string_012 <- function(x) {
  paste(ifelse(is.na(x), "9", as.character(x)), collapse = "")
}

encode_string_AB <- function(x) {
  map <- c("A", "H", "B")
  paste(ifelse(is.na(x), "N", map[x + 1L]), collapse = "")
}

# Per-sample continuous dosage distance (Manhattan on raw dosage)
manhattan_dist_continuous <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  n_v <- sum(ok)
  if (n_v == 0) return(list(dist = NA_real_, n_valid = 0L))
  list(dist = sum(abs(a[ok] - b[ok])) / n_v, n_valid = n_v)
}

# Per-sample discrete distance (Hamming on 0/1/2)
hamming_012 <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  n_v <- sum(ok)
  if (n_v == 0) return(list(dist = NA_real_, n_match = 0L, n_mismatch = 0L, n_missing = length(a), n_valid = 0L))
  mm <- sum(a[ok] != b[ok])
  list(dist = mm / n_v, n_match = n_v - mm, n_mismatch = mm,
       n_missing = length(a) - n_v, n_valid = n_v)
}

# Consensus profile (majority vote per marker)
consensus_012 <- function(mat) {
  apply(mat, 1, function(row) {
    row <- row[!is.na(row)]
    if (length(row) == 0) return(NA_integer_)
    as.integer(which.max(tabulate(row + 1L, nbins = 3L)) - 1L)
  })
}

# Adjusted Rand Index
adj_rand <- function(a, b) {
  if (length(a) != length(b) || length(a) < 2) return(NA_real_)
  tab <- table(a, b); n <- sum(tab)
  if (n < 2) return(NA_real_)
  sc <- sum(choose(tab, 2))
  sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  e <- sa * sb / choose(n, 2); mx <- 0.5 * (sa + sb)
  if (mx == e) return(1)
  (sc - e) / (mx - e)
}

# Load STEP12 regional PCA for a candidate
load_step12 <- function(cid) {
  ff <- list.files(STEP12_DIR, pattern = paste0("candidate_", cid, "\\.regional_pca_samples"),
                   full.names = TRUE)
  if (length(ff) == 0) return(NULL)
  tryCatch(fread(ff[1]), error = function(e) NULL)
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════════

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id; chr <- row$chrom
  c_start <- as.numeric(row$start_bp); c_end <- as.numeric(row$end_bp)

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  ensure_dir(cand_dir)

  message("\n[INFO] ════════════════════════════════════════════════════════")
  message("[INFO] Candidate ", cid, " (", chr, ":",
          format(c_start, big.mark = ","), "-", format(c_end, big.mark = ","), ")")

  # ─────────────────────────────────────────────────────────────────────────
  # LOAD RAW MARKER-LEVEL DOSAGE (bypassing STEP09 PCA compression)
  # ─────────────────────────────────────────────────────────────────────────
  dos_file <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[SKIP] Missing dosage/sites for ", chr); next
  }

  dos <- fread(dos_file); sites <- fread(sites_file)
  sample_cols <- setdiff(names(dos), "marker")

  # Map Ind-style columns to real sample names
  step12 <- load_step12(cid)
  if (!is.null(step12) && all(grepl("^Ind", sample_cols)) && "sample" %in% names(step12)) {
    real_names <- step12$sample
    if (length(real_names) == length(sample_cols)) {
      setnames(dos, old = sample_cols, new = real_names)
      sample_cols <- real_names
    }
  }

  # Filter to candidate region
  keep <- which(sites$pos >= c_start & sites$pos <= c_end)
  if (length(keep) < 20) { message("[SKIP] <20 markers"); next }
  sites_reg <- sites[keep]; dos_reg <- dos[keep]

  # ─────────────────────────────────────────────────────────────────────────
  # BUILD THREE RAW ORDERED MARKER VECTOR ENCODINGS
  # These are first-class objects: preserved before any PCA/MDS/polarity.
  # ─────────────────────────────────────────────────────────────────────────

  # X_minor: raw continuous minor dosage (markers × samples)
  X_minor <- as.matrix(dos_reg[, ..sample_cols])
  storage.mode(X_minor) <- "double"
  n_markers <- nrow(X_minor); n_samples <- ncol(X_minor)

  # X_major: raw continuous major dosage (2 - minor)
  X_major <- 2.0 - X_minor

  # X_012: discretized from minor dosage
  X_012 <- matrix(NA_integer_, n_markers, n_samples)
  for (mi in seq_len(n_markers)) {
    X_012[mi, ] <- discretize_012(X_minor[mi, ])
  }

  message("[INFO] Markers: ", n_markers, "  Samples: ", n_samples)

  # Write raw ordered marker vectors
  write_marker_matrix <- function(X, fname, sample_names, markers) {
    dt <- data.table(marker = markers)
    for (si in seq_len(ncol(X))) {
      set(dt, j = sample_names[si], value = round(X[, si], 6))
    }
    fwrite(dt, file.path(cand_dir, fname), sep = "\t")
  }

  message("[INFO] Writing raw ordered marker vectors (3 encodings)")
  write_marker_matrix(X_minor, "candidate_raw_minor_vectors.tsv.gz", sample_cols, sites_reg$marker)
  write_marker_matrix(X_major, "candidate_raw_major_vectors.tsv.gz", sample_cols, sites_reg$marker)

  # 012 matrix
  dt_012 <- data.table(marker = sites_reg$marker)
  for (si in seq_len(n_samples)) set(dt_012, j = sample_cols[si], value = X_012[, si])
  fwrite(dt_012, file.path(cand_dir, "candidate_012_marker_matrix.tsv.gz"), sep = "\t")

  # ─────────────────────────────────────────────────────────────────────────
  # PER-SAMPLE PROFILE STRINGS
  # ─────────────────────────────────────────────────────────────────────────

  prof_dt <- data.table(
    sample = sample_cols,
    candidate_id = cid,
    profile_012 = sapply(seq_len(n_samples), function(si) encode_string_012(X_012[, si])),
    profile_AB  = sapply(seq_len(n_samples), function(si) encode_string_AB(X_012[, si])),
    n_markers   = n_markers,
    n_missing   = colSums(is.na(X_012)),
    n_hom0      = colSums(X_012 == 0L, na.rm = TRUE),
    n_het       = colSums(X_012 == 1L, na.rm = TRUE),
    n_hom2      = colSums(X_012 == 2L, na.rm = TRUE)
  )

  # Attach coarse group if available
  coarse_groups <- NULL
  if (!is.null(step12) && "group_label" %in% names(step12)) {
    coarse_groups <- setNames(as.character(step12$group_label), step12$sample)
    prof_dt[, coarse_group := coarse_groups[sample]]
  }

  fwrite(prof_dt, file.path(cand_dir, "candidate_sample_profile_strings.tsv"), sep = "\t")

  # ─────────────────────────────────────────────────────────────────────────
  # PAIRWISE PROFILE DISTANCES — computed independently for each encoding
  # ─────────────────────────────────────────────────────────────────────────

  message("[INFO] Computing pairwise profile distances (3 encodings)")

  compute_pairwise <- function(X, dist_fn, label) {
    n <- ncol(X)
    pairs <- vector("list", n * (n - 1) / 2)
    pi <- 0L
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        pi <- pi + 1L
        d <- dist_fn(X[, i], X[, j])
        pairs[[pi]] <- data.table(
          sample1 = sample_cols[i], sample2 = sample_cols[j],
          dist = round(d$dist, 6), n_valid = d$n_valid
        )
      }
    }
    dt <- rbindlist(pairs)
    dt[, encoding := label]
    dt[, candidate_id := cid]
    dt
  }

  # Minor encoding: Manhattan distance on continuous dosage
  pd_minor <- compute_pairwise(X_minor, manhattan_dist_continuous, "minor_continuous")
  fwrite(pd_minor, file.path(cand_dir, "candidate_pairwise_profile_distance_minor.tsv.gz"), sep = "\t")

  # Major encoding: Manhattan distance on continuous dosage
  pd_major <- compute_pairwise(X_major, manhattan_dist_continuous, "major_continuous")
  fwrite(pd_major, file.path(cand_dir, "candidate_pairwise_profile_distance_major.tsv.gz"), sep = "\t")

  # 012 encoding: Hamming distance on discrete genotype
  compute_hamming_pairwise <- function(X_int) {
    n <- ncol(X_int)
    pairs <- vector("list", n * (n - 1) / 2)
    pi <- 0L
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        pi <- pi + 1L
        h <- hamming_012(X_int[, i], X_int[, j])
        pairs[[pi]] <- data.table(
          sample1 = sample_cols[i], sample2 = sample_cols[j],
          hamming_dist = round(h$dist, 6),
          n_match = h$n_match, n_mismatch = h$n_mismatch,
          n_missing = h$n_missing, n_valid = h$n_valid
        )
      }
    }
    dt <- rbindlist(pairs)
    dt[, encoding := "discrete_012"]
    dt[, candidate_id := cid]
    dt
  }

  pd_012 <- compute_hamming_pairwise(X_012)
  fwrite(pd_012, file.path(cand_dir, "candidate_pairwise_profile_distance_012.tsv.gz"), sep = "\t")

  # ─────────────────────────────────────────────────────────────────────────
  # REFERENCE / CONSENSUS COMPARISON
  # ─────────────────────────────────────────────────────────────────────────

  message("[INFO] Computing reference profile scores")

  if (!is.null(coarse_groups)) {
    groups <- unique(na.omit(coarse_groups[sample_cols]))
    consensuses_012 <- list()
    consensuses_minor <- list()

    for (g in groups) {
      g_idx <- which(coarse_groups[sample_cols] == g)
      if (length(g_idx) >= 3) {
        consensuses_012[[g]] <- consensus_012(X_012[, g_idx, drop = FALSE])
        consensuses_minor[[g]] <- rowMeans(X_minor[, g_idx, drop = FALSE], na.rm = TRUE)
      }
    }

    if (length(consensuses_012) >= 2) {
      ref_list <- list()
      for (si in seq_len(n_samples)) {
        for (g in names(consensuses_012)) {
          # 012 Hamming comparison
          h <- hamming_012(X_012[, si], consensuses_012[[g]])
          # Minor continuous comparison
          mc <- manhattan_dist_continuous(X_minor[, si], consensuses_minor[[g]])

          ref_list[[length(ref_list) + 1]] <- data.table(
            sample = sample_cols[si],
            reference = g,
            hamming_dist = round(h$dist, 6),
            hamming_n_match = h$n_match,
            hamming_n_mismatch = h$n_mismatch,
            manhattan_dist_minor = round(mc$dist, 6),
            n_valid = h$n_valid,
            match_fraction = round(h$n_match / max(h$n_valid, 1), 4),
            candidate_id = cid
          )
        }
      }

      ref_dt <- rbindlist(ref_list)
      fwrite(ref_dt, file.path(cand_dir, "candidate_reference_profile_scores.tsv"), sep = "\t")

      # Best assignment (using 012 Hamming as primary)
      best_dt <- ref_dt[, {
        ord <- order(hamming_dist)
        b <- .SD[ord[1]]
        s <- if (length(ord) >= 2) .SD[ord[2]] else b
        list(
          best_reference = b$reference,
          best_hamming = b$hamming_dist,
          best_match_frac = b$match_fraction,
          second_reference = s$reference,
          second_hamming = s$hamming_dist,
          margin = round(s$hamming_dist - b$hamming_dist, 6),
          best_manhattan_minor = b$manhattan_dist_minor
        )
      }, by = sample]
      best_dt[, candidate_id := cid]

      if (!is.null(coarse_groups)) {
        best_dt[, coarse_group := coarse_groups[sample]]
      }

      # Quality tiers from margin
      best_dt[, quality_tier := fifelse(
        margin >= 0.15, "core",
        fifelse(margin >= 0.05, "peripheral", "ambiguous")
      )]

      # Agreement check: does best_reference match coarse_group?
      best_dt[, agrees_with_coarse := (best_reference == coarse_group)]

      fwrite(best_dt, file.path(cand_dir, "candidate_best_profile_assignment.tsv"), sep = "\t")
    }
  }

  # ─────────────────────────────────────────────────────────────────────────
  # LAYER 1 — WINDOW PROFILE (from STEP09 window coordinates + raw dosage)
  # ─────────────────────────────────────────────────────────────────────────

  s09 <- step09_by_chr[[chr]]
  retained_windows <- data.table()

  if (!is.null(s09)) {
    wmeta <- s09$window_meta
    woverlap <- wmeta[end_bp >= c_start & start_bp <= c_end]

    if (nrow(woverlap) > 0) {
      win_profiles <- list()
      window_group_vectors <- list()

      for (wi in seq_len(nrow(woverlap))) {
        wrow <- woverlap[wi]
        w_start <- wrow$start_bp; w_end <- wrow$end_bp
        w_idx <- which(sites_reg$pos >= w_start & sites_reg$pos <= w_end)
        if (length(w_idx) < 10) next

        # Per-window PCA on RAW DOSAGE (not STEP09 summary)
        Xw <- X_minor[w_idx, , drop = FALSE]
        pca_w <- tryCatch({
          pc <- prcomp(t(Xw), center = TRUE, scale. = FALSE)
          list(pc1 = pc$x[, 1], var1 = pc$sdev[1]^2 / sum(pc$sdev^2))
        }, error = function(e) NULL)
        if (is.null(pca_w)) next

        # Adaptive k on PC1
        km2 <- tryCatch(kmeans(matrix(pca_w$pc1, ncol = 1), centers = 2, nstart = 20), error = function(e) NULL)
        km3 <- tryCatch(kmeans(matrix(pca_w$pc1, ncol = 1), centers = 3, nstart = 20), error = function(e) NULL)

        best_k <- 3L
        km_best <- km3
        if (!is.null(km3) && !is.null(km2)) {
          if (km3$tot.withinss / km2$tot.withinss > 0.85) {
            best_k <- 2L; km_best <- km2
          }
        } else if (!is.null(km2)) {
          best_k <- 2L; km_best <- km2
        }
        if (is.null(km_best)) next

        # Order groups by PC1 centroid
        gc <- tapply(pca_w$pc1, km_best$cluster, mean)
        go <- order(gc)
        ordered_grp <- match(km_best$cluster, go)

        win_profiles[[length(win_profiles) + 1]] <- data.table(
          window_idx = wi,
          window_start = w_start, window_end = w_end,
          window_mid = (w_start + w_end) / 2,
          n_markers = length(w_idx), best_k = best_k,
          pc1_var_frac = round(pca_w$var1, 4),
          grp1_n = sum(ordered_grp == 1),
          grp2_n = sum(ordered_grp == 2),
          grp3_n = if (best_k >= 3) sum(ordered_grp == 3) else NA_integer_
        )

        window_group_vectors[[as.character(length(win_profiles))]] <- ordered_grp
      }

      if (length(win_profiles) > 0) {
        retained_windows <- rbindlist(win_profiles)
        retained_windows[, candidate_id := cid]
      }
    }
  }

  fwrite(retained_windows, file.path(cand_dir, "candidate_retained_windows.tsv"), sep = "\t")
  message("[INFO] Retained windows: ", nrow(retained_windows))

  # ─────────────────────────────────────────────────────────────────────────
  # WINDOW GROUPING BY INDUCED SAMPLE STRUCTURE
  # ─────────────────────────────────────────────────────────────────────────

  if (exists("window_group_vectors") && length(window_group_vectors) >= 2) {
    message("[INFO] Computing window similarity by induced sample structure")
    wnames <- names(window_group_vectors)
    nw <- length(wnames)

    sim_mat <- matrix(NA_real_, nw, nw, dimnames = list(wnames, wnames))
    for (i in seq_len(nw)) {
      sim_mat[i, i] <- 1.0
      if (i < nw) {
        for (j in (i + 1):nw) {
          ari <- adj_rand(window_group_vectors[[wnames[i]]], window_group_vectors[[wnames[j]]])
          sim_mat[i, j] <- ari; sim_mat[j, i] <- ari
        }
      }
    }

    sim_dt <- as.data.table(sim_mat, keep.rownames = "window_idx")
    sim_dt[, candidate_id := cid]
    fwrite(sim_dt, file.path(cand_dir, "candidate_window_similarity_matrix.tsv"), sep = "\t")

    # Cluster windows
    if (nw >= 3) {
      sc <- sim_mat; sc[is.na(sc)] <- 0
      hc <- hclust(as.dist(1 - pmax(sc, 0)), method = "average")
      wclust <- cutree(hc, h = 0.5)
    } else {
      wclust <- setNames(rep(1L, nw), wnames)
    }

    wg_dt <- data.table(
      window_idx = as.integer(wnames),
      window_group = wclust,
      mean_ari = sapply(seq_along(wnames), function(i) {
        g <- wclust[i]; sg <- which(wclust == g)
        if (length(sg) <= 1) return(1.0)
        mean(sim_mat[i, sg], na.rm = TRUE)
      })
    )
    wg_dt[, candidate_id := cid]

    max_grp <- wg_dt[, .N, by = window_group][which.max(N), window_group]
    wg_dt[, window_system := fifelse(window_group == max_grp, "core",
                              fifelse(mean_ari >= 0.5, "coherent_minor", "weak"))]

    fwrite(wg_dt, file.path(cand_dir, "candidate_window_group_assignments.tsv"), sep = "\t")

    ws_dt <- wg_dt[, .(n_windows = .N, mean_ari = round(mean(mean_ari, na.rm = TRUE), 4)),
                    by = .(window_group, window_system)]
    ws_dt[, candidate_id := cid]
    fwrite(ws_dt, file.path(cand_dir, "candidate_window_group_summary.tsv"), sep = "\t")
  }

  # ─────────────────────────────────────────────────────────────────────────
  # MARKER INFORMATIVENESS
  # ─────────────────────────────────────────────────────────────────────────

  if (!is.null(coarse_groups)) {
    groups <- unique(na.omit(coarse_groups[sample_cols]))
    if (length(groups) >= 2) {
      mk_info <- data.table(marker = sites_reg$marker, chrom = chr, pos = sites_reg$pos)
      for (g in groups) {
        g_idx <- which(coarse_groups[sample_cols] == g)
        if (length(g_idx) >= 2) {
          gn <- gsub("[^a-zA-Z0-9]", "_", g)
          mk_info[, paste0("mean_minor_", gn) := rowMeans(X_minor[, g_idx, drop = FALSE], na.rm = TRUE)]
          mk_info[, paste0("mean_major_", gn) := rowMeans(X_major[, g_idx, drop = FALSE], na.rm = TRUE)]
          mk_info[, paste0("sd_", gn) := apply(X_minor[, g_idx, drop = FALSE], 1, sd, na.rm = TRUE)]
        }
      }

      minor_cols <- grep("^mean_minor_", names(mk_info), value = TRUE)
      if (length(minor_cols) >= 2) {
        mm <- as.matrix(mk_info[, ..minor_cols])
        mk_info[, max_minor_diff := apply(mm, 1, function(r) diff(range(r, na.rm = TRUE)))]
        mk_info[, informative := max_minor_diff >= 0.3]
      }
      mk_info[, candidate_id := cid]
      fwrite(mk_info, file.path(cand_dir, "candidate_marker_informativeness.tsv"), sep = "\t")
    }
  }

  # ─────────────────────────────────────────────────────────────────────────
  # WINDOW PROFILE TABLE (consolidated summary)
  # ─────────────────────────────────────────────────────────────────────────

  if (nrow(retained_windows) > 0) {
    pt <- copy(retained_windows)
    if (exists("mk_info") && "informative" %in% names(mk_info)) {
      pt[, n_informative := sum(mk_info$informative, na.rm = TRUE)]
      pt[, frac_informative := round(sum(mk_info$informative, na.rm = TRUE) / nrow(mk_info), 4)]
    }
    pt[, signature_mode := if (clair3_available) "dosage+clair3" else "dosage_only"]
    fwrite(pt, file.path(cand_dir, "candidate_window_profile_table.tsv"), sep = "\t")
  }

  message("[INFO] Candidate ", cid, " complete")
}

message("\n[DONE] STEP20 v4 complete")
