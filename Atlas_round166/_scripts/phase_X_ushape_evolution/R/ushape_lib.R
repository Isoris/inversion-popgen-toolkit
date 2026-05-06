# =============================================================================
# R/ushape_lib.R
# =============================================================================
# Pure stat functions for the U-shape evolution module. Sourceable from CLI
# scripts AND from the popstats-server R subprocess. Side-effect free.
#
# Two entry surfaces:
#
#   compute_window_stats_from_dosage(dosage_dt, candidate, groups, params)
#       -> data.table of per-window stats
#       Used by STEP_U01_cli.R when running offline against dosage TSVs.
#
#   compute_window_stats_from_popstats_json(popstats_payload, candidate, params)
#       -> data.table of per-window stats
#       Used by the server endpoint: popstats_server already produces FST/dXY/
#       theta-pi per window; we just rezone, normalize, and score.
#
# The downstream scoring/classification functions consume EITHER source
# identically — they only see the window data.table.
#
# Conventions:
#   - All coordinate inputs 1-based, inclusive on both sides (matches the
#     project's existing TSVs and popstats_server contracts).
#   - Group labels in input: HOMO_1, HET, HOMO_2 (extras allowed but ignored
#     unless explicitly opted in).
#   - Hudson FST is the default. Weir-Cockerham is documented as a stub —
#     it requires per-individual genotype counts that the server output
#     doesn't currently expose.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# -----------------------------------------------------------------------------
# zoning
# -----------------------------------------------------------------------------
# Given a candidate interval [start, end] and flank length, classify a window
# midpoint into one of: left_flank, left_edge, center, right_edge, right_flank.
.zone_of_midpoint <- function(mid, cand_start, cand_end,
                              edge_fraction = 0.20,
                              center_fraction = 0.60) {
  L <- cand_end - cand_start + 1L
  edge_bp <- L * edge_fraction

  if (mid < cand_start)                                return("left_flank")
  if (mid > cand_end)                                  return("right_flank")
  if (mid <= cand_start + edge_bp)                     return("left_edge")
  if (mid >= cand_end   - edge_bp)                     return("right_edge")
  return("center")
}

assign_zones <- function(dt, cand_start, cand_end,
                         edge_fraction = 0.20,
                         center_fraction = 0.60) {
  stopifnot("start_bp" %in% names(dt), "end_bp" %in% names(dt))
  mids <- (dt$start_bp + dt$end_bp) / 2
  dt[, zone := vapply(mids, .zone_of_midpoint, character(1L),
                      cand_start = cand_start, cand_end = cand_end,
                      edge_fraction = edge_fraction,
                      center_fraction = center_fraction)]
  dt[, relative_position := pmin(1, pmax(0, (mids - cand_start) / max(1, cand_end - cand_start)))]
  dt[]
}


# -----------------------------------------------------------------------------
# per-SNP stats (used only by the offline dosage path)
# -----------------------------------------------------------------------------
.snp_stats_from_dosage <- function(dosage_row, samples_g1, samples_g2) {
  d1 <- as.numeric(dosage_row[samples_g1])
  d2 <- as.numeric(dosage_row[samples_g2])
  n1 <- sum(!is.na(d1))
  n2 <- sum(!is.na(d2))
  if (n1 < 1L || n2 < 1L) return(NULL)
  p1 <- mean(d1, na.rm = TRUE) / 2
  p2 <- mean(d2, na.rm = TRUE) / 2
  pi1 <- 2 * p1 * (1 - p1)
  pi2 <- 2 * p2 * (1 - p2)
  dxy <- p1 * (1 - p2) + p2 * (1 - p1)
  list(p1 = p1, p2 = p2, pi1 = pi1, pi2 = pi2, dxy = dxy,
       afd = abs(p1 - p2), n1 = n1, n2 = n2)
}

# Hudson FST estimator (Bhatia et al. 2013 ratio-of-averages),
# computed per window from per-SNP numerator/denominator.
.hudson_fst_terms <- function(p1, p2, n1, n2) {
  # numerator: (p1 - p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)
  # denominator: p1*(1-p2) + p2*(1-p1)
  num <- (p1 - p2)^2 -
    ifelse(n1 > 1, p1 * (1 - p1) / (n1 - 1), 0) -
    ifelse(n2 > 1, p2 * (1 - p2) / (n2 - 1), 0)
  den <- p1 * (1 - p2) + p2 * (1 - p1)
  list(num = num, den = den)
}


# -----------------------------------------------------------------------------
# offline dosage path
# -----------------------------------------------------------------------------
# dosage_dt: data.table with columns chrom, pos, then one column per sample
# candidate: list(chrom, start_bp, end_bp, left_bp, right_bp)
# groups: list(HOMO_1 = c(...), HET = c(...), HOMO_2 = c(...))
# params: list(n_windows_inside, max_flank_bp, min_flank_bp, min_window_bp,
#              edge_fraction, center_fraction)
compute_window_stats_from_dosage <- function(dosage_dt, candidate, groups, params) {
  chrom_dat <- dosage_dt[chrom == candidate$chrom]
  if (nrow(chrom_dat) == 0L) {
    warning("no dosage rows for chrom ", candidate$chrom)
    return(data.table())
  }

  L <- candidate$end_bp - candidate$start_bp + 1L
  flank_bp <- max(params$min_flank_bp,
                  min(params$max_flank_bp, ceiling(L * 0.5)))
  win_bp <- max(params$min_window_bp,
                ceiling(L / params$n_windows_inside))

  # build window scaffold across [start - flank, end + flank]
  scan_lo <- max(1L, candidate$start_bp - flank_bp)
  scan_hi <- candidate$end_bp + flank_bp

  win_starts <- seq(scan_lo, scan_hi, by = win_bp)
  windows <- data.table(
    start_bp = win_starts,
    end_bp   = pmin(scan_hi, win_starts + win_bp - 1L)
  )

  g1 <- intersect(groups$HOMO_1, names(chrom_dat))
  g2 <- intersect(groups$HOMO_2, names(chrom_dat))
  if (length(g1) < params$min_group_n || length(g2) < params$min_group_n) {
    warning("insufficient HOMO_1/HOMO_2 samples for ", candidate$candidate_id)
    return(data.table())
  }

  # iterate windows (vectorisable but kept readable)
  out <- vector("list", nrow(windows))
  for (i in seq_len(nrow(windows))) {
    w_lo <- windows$start_bp[i]; w_hi <- windows$end_bp[i]
    snps <- chrom_dat[pos >= w_lo & pos <= w_hi]
    if (nrow(snps) == 0L) {
      out[[i]] <- data.table(
        start_bp = w_lo, end_bp = w_hi, n_snps = 0L,
        pi_homo1 = NA_real_, pi_homo2 = NA_real_,
        dxy_homo1_homo2 = NA_real_, fst_homo1_homo2 = NA_real_,
        allele_freq_delta = NA_real_,
        private_homo1 = 0L, private_homo2 = 0L,
        fixed_diff = 0L, shared_poly = 0L
      )
      next
    }
    pi1 <- pi2 <- dxy <- afd <- numeric(nrow(snps))
    fst_num <- fst_den <- numeric(nrow(snps))
    priv1 <- priv2 <- fixed <- shared <- 0L
    for (s in seq_len(nrow(snps))) {
      ss <- .snp_stats_from_dosage(as.list(snps[s]), g1, g2)
      if (is.null(ss)) {
        pi1[s] <- pi2[s] <- dxy[s] <- afd[s] <- NA_real_
        fst_num[s] <- fst_den[s] <- NA_real_
        next
      }
      pi1[s] <- ss$pi1; pi2[s] <- ss$pi2; dxy[s] <- ss$dxy; afd[s] <- ss$afd
      ft <- .hudson_fst_terms(ss$p1, ss$p2, ss$n1, ss$n2)
      fst_num[s] <- ft$num; fst_den[s] <- ft$den
      if (ss$p1 > 0 && ss$p2 == 0)              priv1 <- priv1 + 1L
      if (ss$p2 > 0 && ss$p1 == 0)              priv2 <- priv2 + 1L
      if ((ss$p1 == 1 && ss$p2 == 0) ||
          (ss$p1 == 0 && ss$p2 == 1))           fixed <- fixed + 1L
      if (ss$p1 > 0 && ss$p1 < 1 &&
          ss$p2 > 0 && ss$p2 < 1)               shared <- shared + 1L
    }
    den_sum <- sum(fst_den, na.rm = TRUE)
    fst_w <- if (is.finite(den_sum) && den_sum > 0) sum(fst_num, na.rm = TRUE) / den_sum else NA_real_
    out[[i]] <- data.table(
      start_bp = w_lo, end_bp = w_hi, n_snps = nrow(snps),
      pi_homo1 = mean(pi1, na.rm = TRUE),
      pi_homo2 = mean(pi2, na.rm = TRUE),
      dxy_homo1_homo2 = mean(dxy, na.rm = TRUE),
      fst_homo1_homo2 = pmax(0, pmin(1, fst_w)),
      allele_freq_delta = mean(afd, na.rm = TRUE),
      private_homo1 = priv1, private_homo2 = priv2,
      fixed_diff = fixed, shared_poly = shared
    )
  }
  res <- rbindlist(out)
  res <- assign_zones(res, candidate$start_bp, candidate$end_bp,
                      edge_fraction = params$edge_fraction,
                      center_fraction = params$center_fraction)
  res[, candidate_id := candidate$candidate_id]
  res[, chrom := candidate$chrom]
  res[, window_id := seq_len(.N)]
  setcolorder(res, c("candidate_id", "chrom", "window_id",
                     "start_bp", "end_bp", "relative_position", "zone",
                     "n_snps",
                     "pi_homo1", "pi_homo2",
                     "dxy_homo1_homo2", "fst_homo1_homo2",
                     "allele_freq_delta",
                     "private_homo1", "private_homo2",
                     "fixed_diff", "shared_poly"))
  res
}


# -----------------------------------------------------------------------------
# server-side path: rezone an existing popstats_groupwise payload
# -----------------------------------------------------------------------------
# popstats_payload: parsed JSON from POST /api/popstats/groupwise with
#                   groups = {HOMO_1: [...], HOMO_2: [...]}.
#                   payload$columns names the metric columns;
#                   payload$windows is a list-of-lists or matrix of values.
#
# We expect at minimum columns: start_bp, end_bp, n_snps,
# fst_HOMO_1_vs_HOMO_2, dxy_HOMO_1_vs_HOMO_2,
# theta_pi_HOMO_1, theta_pi_HOMO_2.
# (Server emits these names by convention; a small adapter handles
# other naming styles.)
#
# This function does NOT compute private_homo* / fixed_diff because those
# require per-SNP iteration which the popstats engine does not expose.
# The server endpoint adds an optional second call to a small helper
# (compute_arrangement_snp_classes) that scans the BEAGLE block once.
# For now, those four columns come back as NA when called from server.
compute_window_stats_from_popstats_json <- function(popstats_payload,
                                                    candidate, params,
                                                    snp_class_dt = NULL) {
  cols <- popstats_payload$columns
  rows <- popstats_payload$windows
  if (is.null(rows) || length(rows) == 0L) return(data.table())

  # rows may be list-of-lists or already a data.frame
  if (is.list(rows) && !is.data.frame(rows)) {
    dt <- rbindlist(lapply(rows, function(r) as.list(setNames(r, cols))),
                    fill = TRUE)
  } else {
    dt <- as.data.table(rows)
    if (!is.null(cols) && length(cols) == ncol(dt)) setnames(dt, cols)
  }

  # column-name normalisation (be tolerant)
  rename_if <- function(dt, target, candidates) {
    for (cand in candidates) {
      if (cand %in% names(dt)) { setnames(dt, cand, target); return(dt) }
    }
    dt[, (target) := NA_real_]
    dt
  }
  dt <- rename_if(dt, "fst_homo1_homo2",
                  c("fst_HOMO_1_vs_HOMO_2", "fst", "FST"))
  dt <- rename_if(dt, "dxy_homo1_homo2",
                  c("dxy_HOMO_1_vs_HOMO_2", "dxy", "dXY"))
  dt <- rename_if(dt, "pi_homo1",
                  c("theta_pi_HOMO_1", "pi_HOMO_1", "pi1"))
  dt <- rename_if(dt, "pi_homo2",
                  c("theta_pi_HOMO_2", "pi_HOMO_2", "pi2"))
  if (!"n_snps" %in% names(dt)) dt[, n_snps := NA_integer_]
  if (!"start_bp" %in% names(dt)) stop("popstats payload lacks start_bp")
  if (!"end_bp"   %in% names(dt)) stop("popstats payload lacks end_bp")

  dt[, allele_freq_delta := NA_real_]   # not provided by region_popstats
  dt[, private_homo1 := NA_integer_]
  dt[, private_homo2 := NA_integer_]
  dt[, fixed_diff := NA_integer_]
  dt[, shared_poly := NA_integer_]

  # if the server has run the per-SNP class scan, merge those counts in
  if (!is.null(snp_class_dt)) {
    setkey(snp_class_dt, start_bp, end_bp)
    dt[, c("private_homo1", "private_homo2", "fixed_diff", "shared_poly") :=
        snp_class_dt[.(start_bp, end_bp), .(private_homo1, private_homo2,
                                            fixed_diff, shared_poly), on = .(start_bp, end_bp)]]
  }

  dt <- assign_zones(dt, candidate$start_bp, candidate$end_bp,
                     edge_fraction = params$edge_fraction,
                     center_fraction = params$center_fraction)
  dt[, candidate_id := candidate$candidate_id]
  dt[, chrom := candidate$chrom]
  dt[, window_id := seq_len(.N)]
  setcolorder(dt, c("candidate_id", "chrom", "window_id",
                    "start_bp", "end_bp", "relative_position", "zone",
                    "n_snps",
                    "pi_homo1", "pi_homo2",
                    "dxy_homo1_homo2", "fst_homo1_homo2",
                    "allele_freq_delta",
                    "private_homo1", "private_homo2",
                    "fixed_diff", "shared_poly"))
  dt
}


# -----------------------------------------------------------------------------
# zone summaries
# -----------------------------------------------------------------------------
zone_means <- function(window_dt, value_col) {
  v <- window_dt[[value_col]]
  zones <- window_dt$zone
  list(
    inside      = mean(v[zones %in% c("left_edge","center","right_edge")], na.rm = TRUE),
    flank       = mean(v[zones %in% c("left_flank","right_flank")],         na.rm = TRUE),
    left_flank  = mean(v[zones == "left_flank"],  na.rm = TRUE),
    left_edge   = mean(v[zones == "left_edge"],   na.rm = TRUE),
    center      = mean(v[zones == "center"],      na.rm = TRUE),
    right_edge  = mean(v[zones == "right_edge"],  na.rm = TRUE),
    right_flank = mean(v[zones == "right_flank"], na.rm = TRUE)
  )
}

candidate_raw_summary <- function(window_dt, candidate, groups) {
  d <- zone_means(window_dt, "dxy_homo1_homo2")
  f <- zone_means(window_dt, "fst_homo1_homo2")
  p1 <- zone_means(window_dt, "pi_homo1")
  p2 <- zone_means(window_dt, "pi_homo2")
  af <- zone_means(window_dt, "allele_freq_delta")
  data.table(
    candidate_id = candidate$candidate_id,
    chrom = candidate$chrom,
    start_bp = candidate$start_bp, end_bp = candidate$end_bp,
    length_bp = candidate$end_bp - candidate$start_bp + 1L,
    n_homo1 = length(groups$HOMO_1),
    n_homo2 = length(groups$HOMO_2),
    n_het   = length(groups$HET %||% character(0)),
    dxy_inside_mean      = d$inside,
    dxy_flank_mean       = d$flank,
    dxy_left_edge_mean   = d$left_edge,
    dxy_right_edge_mean  = d$right_edge,
    dxy_center_mean      = d$center,
    fst_inside_mean      = f$inside,
    fst_flank_mean       = f$flank,
    fst_left_edge_mean   = f$left_edge,
    fst_right_edge_mean  = f$right_edge,
    fst_center_mean      = f$center,
    pi_homo1_inside_mean = p1$inside,
    pi_homo2_inside_mean = p2$inside,
    pi_homo1_flank_mean  = p1$flank,
    pi_homo2_flank_mean  = p2$flank,
    allele_freq_delta_inside_mean = af$inside,
    private_snp_homo1_count = sum(window_dt$private_homo1, na.rm = TRUE),
    private_snp_homo2_count = sum(window_dt$private_homo2, na.rm = TRUE),
    fixed_diff_count        = sum(window_dt$fixed_diff,    na.rm = TRUE),
    shared_snp_count        = sum(window_dt$shared_poly,   na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# flank-resampling matched background
# -----------------------------------------------------------------------------
# Bootstraps the flank windows to build a local null distribution for the
# inside summary statistics. Returns a one-row data.table of z-scores.
#
# Used by:
#   - STEP_U03 (offline batch)
#   - the server endpoint (inline, same call, same numbers)
#
# Why flank-only: cohort-scale chromosome-wide windows aren't on disk and
# are expensive to recompute. Flank-resampling is honest: it captures the
# LOCAL variance instead of pretending the flank mean is a point estimate
# with zero variance. It is NOT a chromosome-wide null — flag this in any
# JSON output as bg_mode = "flank_resample".
flank_resample_zscores <- function(window_dt, n_random = 1000L, seed = 42L) {
  set.seed(seed)
  flank  <- window_dt[zone %in% c("left_flank", "right_flank")]
  inside <- window_dt[zone %in% c("left_edge", "center", "right_edge")]
  edge   <- window_dt[zone %in% c("left_edge", "right_edge")]
  centr  <- window_dt[zone == "center"]

  flank_n <- nrow(flank); inside_n <- nrow(inside)
  edge_n  <- nrow(edge);  cent_n   <- nrow(centr)

  if (inside_n < 5L || flank_n < 5L) {
    return(data.table(
      z_dxy_inside = NA_real_, z_fst_inside = NA_real_,
      z_dxy_u_score = NA_real_, z_dxy_internal_peak = NA_real_,
      bg_n_pool = flank_n, bg_mode = "skip_low_n"
    ))
  }

  obs_inside_dxy <- mean(inside$dxy_homo1_homo2, na.rm = TRUE)
  obs_inside_fst <- mean(inside$fst_homo1_homo2, na.rm = TRUE)
  obs_dxy_u      <- mean(edge$dxy_homo1_homo2, na.rm = TRUE) /
                    max(1e-12, mean(centr$dxy_homo1_homo2, na.rm = TRUE))
  obs_dxy_peak   <- suppressWarnings(max(centr$dxy_homo1_homo2, na.rm = TRUE)) /
                    max(1e-12, mean(edge$dxy_homo1_homo2, na.rm = TRUE))

  rep_inside_dxy <- numeric(n_random); rep_inside_fst <- numeric(n_random)
  rep_dxy_u      <- numeric(n_random); rep_dxy_peak   <- numeric(n_random)

  for (b in seq_len(n_random)) {
    samp_in <- flank[sample.int(flank_n, inside_n, replace = TRUE)]
    rep_inside_dxy[b] <- mean(samp_in$dxy_homo1_homo2, na.rm = TRUE)
    rep_inside_fst[b] <- mean(samp_in$fst_homo1_homo2, na.rm = TRUE)
    if (edge_n >= 1L && cent_n >= 1L) {
      samp_e <- flank[sample.int(flank_n, edge_n, replace = TRUE)]
      samp_c <- flank[sample.int(flank_n, cent_n, replace = TRUE)]
      rep_dxy_u[b] <- mean(samp_e$dxy_homo1_homo2, na.rm = TRUE) /
                       max(1e-12, mean(samp_c$dxy_homo1_homo2, na.rm = TRUE))
      rep_dxy_peak[b] <- suppressWarnings(max(samp_c$dxy_homo1_homo2, na.rm = TRUE)) /
                          max(1e-12, mean(samp_e$dxy_homo1_homo2, na.rm = TRUE))
    } else {
      rep_dxy_u[b] <- NA_real_; rep_dxy_peak[b] <- NA_real_
    }
  }

  z <- function(o, r) {
    rs <- sd(r, na.rm = TRUE)
    if (!is.finite(rs) || rs < 1e-12) return(NA_real_)
    (o - mean(r, na.rm = TRUE)) / rs
  }
  data.table(
    z_dxy_inside        = z(obs_inside_dxy, rep_inside_dxy),
    z_fst_inside        = z(obs_inside_fst, rep_inside_fst),
    z_dxy_u_score       = z(obs_dxy_u,      rep_dxy_u),
    z_dxy_internal_peak = z(obs_dxy_peak,   rep_dxy_peak),
    bg_n_pool           = flank_n,
    bg_mode             = "flank_resample"
  )
}


`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
