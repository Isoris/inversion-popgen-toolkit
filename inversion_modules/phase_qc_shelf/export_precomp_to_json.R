#!/usr/bin/env Rscript
# =============================================================================
# export_precomp_to_json.R  (v3 — with theta + ancestry track support)
#
# Convert a precomp RDS into a compact JSON bundle for the PCA scrubber.
# v3 optionally embeds:
#   - per-sample per-window θπ from aggregated pestPG tracks (Q05)
#   - per-window ancestry delta12/entropy/ena from Engine B cache (Q06)
#
# Usage:
#   Rscript export_precomp_to_json.R \
#     --precomp   <LGNN.precomp.rds> \
#     --sim_mat   <optional separate sim_mat RDS for v9.1> \
#     --samples   <tab file: ind<TAB>cga<TAB>ancestry> \
#     --theta     <theta.<CHR>.<SCALE>.tsv.gz from Q05> \
#     --ancestry  <ancestry_window.<CHR>.tsv from Q06> \
#     --out       <output .json>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ---- CLI ---------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) == 0) return(default); args[i + 1]
}
PRECOMP <- get_arg("--precomp")
SIMMAT  <- get_arg("--sim_mat", NULL)
SAMPLES <- get_arg("--samples", NULL)
THETA   <- get_arg("--theta", NULL)
ANCESTRY <- get_arg("--ancestry", NULL)
ANCESTRY_SAMPLES <- get_arg("--ancestry_samples", NULL)  # maxQ wide matrix
OUT     <- get_arg("--out", "precomp_scrubber.json")
THUMB_N <- as.integer(get_arg("--thumb_n", "200"))
stopifnot(!is.null(PRECOMP), file.exists(PRECOMP))

message("[export] Loading precomp: ", PRECOMP)
pc <- readRDS(PRECOMP)
dt <- as.data.table(pc$dt)
chrom <- pc$chrom %||% dt$chrom[1]
n_win <- pc$n_windows %||% nrow(dt)

# ---- PC columns --------------------------------------------------------------
pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
pc2_cols <- grep("^PC_2_", names(dt), value = TRUE)
if (length(pc1_cols) == 0) stop("No PC_1_* columns in precomp$dt")

sample_ids <- sub("^PC_1_", "", pc1_cols)
n_samp <- length(sample_ids)
message("[export] Windows: ", n_win, "   Samples: ", n_samp)

# ---- Sample metadata ---------------------------------------------------------
sample_meta <- data.table(ind = sample_ids, cga = sample_ids, ancestry = "unknown")
if (!is.null(SAMPLES) && file.exists(SAMPLES)) {
  message("[export] Loading sample metadata: ", SAMPLES)
  smeta <- fread(SAMPLES, header = TRUE)
  setnames(smeta, tolower(names(smeta)))
  if ("ind" %in% names(smeta)) {
    sample_meta <- merge(sample_meta[, .(ind)], smeta, by = "ind",
                         all.x = TRUE, sort = FALSE)
  } else if ("cga" %in% names(smeta)) {
    # Precomp uses CGA names directly
    sample_meta[, cga := ind]
    sample_meta <- merge(sample_meta, smeta, by = "cga", all.x = TRUE, sort = FALSE)
  }
  if (!"ancestry" %in% names(sample_meta)) sample_meta[, ancestry := "unknown"]
  if (!"cga" %in% names(sample_meta))      sample_meta[, cga := ind]
  sample_meta[is.na(cga),      cga      := ind]
  sample_meta[is.na(ancestry), ancestry := "unknown"]
}

# ---- Z column ----------------------------------------------------------------
z_col <- NULL
for (cand in c("robust_z", "z_robust", "z", "mds_z_robust", "mds_z1_robust")) {
  if (cand %in% names(dt)) { z_col <- cand; break }
}
if (is.null(z_col)) {
  message("[export] No Z column found, setting to 0")
  dt[, .export_z := 0]; z_col <- ".export_z"
}

# ---- Theta alignment ---------------------------------------------------------
# Output: per-sample vector of theta_pi_persite aligned to precomp windows
theta_by_sample <- NULL
theta_range <- NULL
if (!is.null(THETA) && file.exists(THETA)) {
  message("[export] Loading theta track: ", THETA)
  theta <- fread(THETA)
  theta <- theta[chrom == chrom]
  samples_in_theta <- intersect(sample_meta$cga, unique(theta$sample))
  message("[export] Theta samples matched to precomp samples: ",
          length(samples_in_theta), " / ", n_samp)

  # Precomp window centers
  win_centers <- as.integer((dt$start_bp + dt$end_bp) / 2)

  # Per-sample alignment: for each precomp window, find nearest theta window
  theta_by_sample <- matrix(NA_real_, nrow = n_win, ncol = n_samp)
  colnames(theta_by_sample) <- sample_meta$cga

  for (s in samples_in_theta) {
    ts <- theta[sample == s]
    if (nrow(ts) < 2) next
    setorder(ts, win_center_bp)
    # Max step distance (twice the median spacing)
    spacing <- median(diff(ts$win_center_bp), na.rm = TRUE)
    max_gap <- 2 * spacing

    # Nearest join via findInterval
    idx <- findInterval(win_centers, ts$win_center_bp, all.inside = TRUE)
    left  <- pmax(1L, idx)
    right <- pmin(nrow(ts), idx + 1L)
    dL <- abs(win_centers - ts$win_center_bp[left])
    dR <- abs(win_centers - ts$win_center_bp[right])
    nearest <- ifelse(dL <= dR, left, right)
    dmin <- pmin(dL, dR)

    vals <- ts$theta_pi_persite[nearest]
    vals[dmin > max_gap] <- NA_real_
    col_idx <- match(s, sample_meta$cga)
    theta_by_sample[, col_idx] <- vals
  }

  theta_range <- range(theta_by_sample, na.rm = TRUE, finite = TRUE)
  message(sprintf("[export] Theta range: [%.3e, %.3e]",
                  theta_range[1], theta_range[2]))
}

# ---- Ancestry alignment (per-window delta12/entropy/ena) --------------------
anc_by_window <- NULL   # list with delta12, entropy, ena vectors aligned to precomp windows
anc_range <- NULL
if (!is.null(ANCESTRY) && file.exists(ANCESTRY)) {
  message("[export] Loading ancestry window track: ", ANCESTRY)
  anc_raw <- fread(ANCESTRY)
  anc_raw <- anc_raw[chrom == chrom]
  if (nrow(anc_raw) > 1) {
    setorder(anc_raw, window_start_bp)
    # Engine B windows are coarser than precomp. Nearest-window join.
    anc_centers <- as.integer((anc_raw$window_start_bp + anc_raw$window_end_bp) / 2)
    win_centers <- as.integer((dt$start_bp + dt$end_bp) / 2)
    spacing <- median(diff(anc_centers), na.rm = TRUE)
    max_gap <- 2 * spacing

    idx <- findInterval(win_centers, anc_centers, all.inside = TRUE)
    left  <- pmax(1L, idx)
    right <- pmin(length(anc_centers), idx + 1L)
    dL <- abs(win_centers - anc_centers[left])
    dR <- abs(win_centers - anc_centers[right])
    nearest <- ifelse(dL <= dR, left, right)
    dmin <- pmin(dL, dR)

    dlt <- anc_raw$delta12[nearest]
    ent <- anc_raw$entropy[nearest]
    ena <- anc_raw$ena[nearest]
    mx_label <- anc_raw$maxQ_label[nearest]
    cvd <- if ("cv_delta12_across_samples" %in% names(anc_raw)) {
      anc_raw$cv_delta12_across_samples[nearest]
    } else rep(NA_real_, length(nearest))

    bad <- dmin > max_gap
    dlt[bad] <- NA_real_; ent[bad] <- NA_real_; ena[bad] <- NA_real_
    mx_label[bad] <- NA_character_; cvd[bad] <- NA_real_

    anc_by_window <- list(
      delta12   = round(dlt, 5),
      entropy   = round(ent, 5),
      ena       = round(ena, 4),
      cv_delta12 = round(cvd, 4),
      maxQ_label = mx_label
    )
    anc_range <- list(
      delta12 = range(dlt, na.rm = TRUE, finite = TRUE),
      entropy = range(ent, na.rm = TRUE, finite = TRUE)
    )
    message(sprintf("[export] Ancestry delta12 range: [%.3f, %.3f]",
                    anc_range$delta12[1], anc_range$delta12[2]))
  }
}

# ---- Per-sample maxQ alignment (for PCA coloring) ---------------------------
# Wide matrix from Q06: window_mid_bp + one column per sample (CGA), values = K1..K8
maxq_by_sample <- NULL   # matrix n_win x n_samp, integer codes of K labels (1..K) or NA
maxq_levels <- NULL
if (!is.null(ANCESTRY_SAMPLES) && file.exists(ANCESTRY_SAMPLES)) {
  message("[export] Loading per-sample maxQ wide matrix: ", ANCESTRY_SAMPLES)
  wide <- fread(ANCESTRY_SAMPLES)
  if ("window_mid_bp" %in% names(wide)) {
    anc_win_centers <- as.integer(wide$window_mid_bp)
    sample_cols <- setdiff(names(wide), "window_mid_bp")
    # Intersect with precomp samples
    common <- intersect(sample_meta$cga, sample_cols)
    message("[export] Per-sample maxQ available for ", length(common),
            " / ", n_samp, " samples")

    if (length(common) > 0) {
      # Align by nearest-window
      win_centers <- as.integer((dt$start_bp + dt$end_bp) / 2)
      setorder(wide, window_mid_bp)
      spacing <- median(diff(wide$window_mid_bp), na.rm = TRUE)
      max_gap <- 2 * spacing

      idx <- findInterval(win_centers, wide$window_mid_bp, all.inside = TRUE)
      left  <- pmax(1L, idx)
      right <- pmin(nrow(wide), idx + 1L)
      dL <- abs(win_centers - wide$window_mid_bp[left])
      dR <- abs(win_centers - wide$window_mid_bp[right])
      nearest <- ifelse(dL <= dR, left, right)
      dmin <- pmin(dL, dR)

      # Collect all unique labels for indexing
      all_labels <- unique(unlist(lapply(common, function(s) unique(wide[[s]]))))
      all_labels <- all_labels[!is.na(all_labels) & nzchar(all_labels)]
      # Sort K1, K2, ..., K10 naturally
      all_labels <- all_labels[order(as.integer(sub("^K", "", all_labels)))]
      maxq_levels <- all_labels
      label_to_idx <- setNames(seq_along(all_labels), all_labels)

      maxq_by_sample <- matrix(NA_integer_, nrow = n_win, ncol = n_samp)
      colnames(maxq_by_sample) <- sample_meta$cga

      for (s in common) {
        col_vals <- wide[[s]][nearest]
        col_vals[dmin > max_gap] <- NA_character_
        col_idx <- match(s, sample_meta$cga)
        maxq_by_sample[, col_idx] <- label_to_idx[col_vals]
      }

      message("[export] Embedded per-sample maxQ. Levels: ",
              paste(all_labels, collapse = ", "))
    }
  } else {
    message("[export] Wide matrix missing window_mid_bp column, skipping")
  }
}

# ---- Build windows array -----------------------------------------------------
message("[export] Building per-window PC arrays ...")
pc1_mat <- as.matrix(dt[, ..pc1_cols])
pc2_mat <- as.matrix(dt[, ..pc2_cols])
round4 <- function(x) round(x, 4)

windows <- vector("list", n_win)
for (i in seq_len(n_win)) {
  w <- list(
    idx       = i - 1L,
    start_bp  = as.integer(dt$start_bp[i]),
    end_bp    = as.integer(dt$end_bp[i]),
    center_mb = round(((dt$start_bp[i] + dt$end_bp[i]) / 2) / 1e6, 4),
    z         = round4(dt[[z_col]][i]),
    lam1      = round4(dt$lam_1[i] %||% NA_real_),
    lam2      = round4(dt$lam_2[i] %||% NA_real_),
    pc1       = round4(pc1_mat[i, ]),
    pc2       = round4(pc2_mat[i, ])
  )
  if (!is.null(theta_by_sample)) {
    w$theta <- signif(theta_by_sample[i, ], 4)
  }
  if (!is.null(anc_by_window)) {
    w$anc_delta12   <- anc_by_window$delta12[i]
    w$anc_entropy   <- anc_by_window$entropy[i]
    w$anc_ena       <- anc_by_window$ena[i]
    w$anc_cv_delta12 <- anc_by_window$cv_delta12[i]
    w$anc_maxQ_label <- anc_by_window$maxQ_label[i]
  }
  if (!is.null(maxq_by_sample)) {
    # integer-encoded per sample; null for samples without data
    w$maxq <- as.integer(maxq_by_sample[i, ])
  }
  windows[[i]] <- w
  if (i %% 500 == 0) message("  window ", i, " / ", n_win)
}

# ---- Sim_mat thumbnail -------------------------------------------------------
sim_thumb <- NULL; sim_thumb_n <- 0
sim_source <- pc$sim_mat
if (is.null(sim_source) && !is.null(SIMMAT) && file.exists(SIMMAT)) {
  message("[export] Loading separate sim_mat: ", SIMMAT)
  sim_source <- readRDS(SIMMAT)
}
if (!is.null(sim_source) && is.matrix(sim_source)) {
  message("[export] Building sim_mat thumbnail (", THUMB_N, "x", THUMB_N, ")")
  N <- nrow(sim_source)
  step <- max(1, floor(N / THUMB_N))
  idx <- seq(1, N, by = step)[seq_len(min(THUMB_N, ceiling(N / step)))]
  thumb <- sim_source[idx, idx, drop = FALSE]
  sim_thumb <- round(as.vector(t(thumb)), 3)
  sim_thumb_n <- length(idx)
}

# ---- Assemble & write --------------------------------------------------------
out_list <- list(
  chrom       = chrom,
  n_windows   = n_win,
  n_samples   = n_samp,
  samples     = lapply(seq_len(nrow(sample_meta)), function(i) {
    list(ind = sample_meta$ind[i], cga = sample_meta$cga[i],
         ancestry = sample_meta$ancestry[i])
  }),
  windows     = windows,
  sim_thumb   = sim_thumb,
  sim_thumb_n = sim_thumb_n,
  theta_range = if (!is.null(theta_range)) as.numeric(theta_range) else NULL,
  has_theta   = !is.null(theta_by_sample),
  has_ancestry = !is.null(anc_by_window),
  has_sample_ancestry = !is.null(maxq_by_sample),
  maxq_levels = maxq_levels,
  anc_delta12_range = if (!is.null(anc_range)) as.numeric(anc_range$delta12) else NULL,
  anc_entropy_range = if (!is.null(anc_range)) as.numeric(anc_range$entropy) else NULL
)

message("[export] Writing JSON: ", OUT)
write_json(out_list, OUT, auto_unbox = TRUE, digits = 6, na = "null")
fs <- file.info(OUT)$size
message(sprintf("[export] Done. Output size: %.1f MB", fs / 1e6))
