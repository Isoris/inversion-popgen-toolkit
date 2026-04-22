#!/usr/bin/env Rscript

# =============================================================================
# 02_ancestral_fragments.R  (v1.1 — registry-wired)
#
# PER-CARRIER ANCESTRAL FRAGMENT SCAN — second stage.
#
# For each inversion carrier (HET + HOMO_2 or HOMO_1 depending on which is
# minor), walk outward from the core block one SNP at a time, testing whether
# THIS SAMPLE's dosage still correlates with the core consensus vector. Stop
# when correlation drops. The stopping positions define this sample's
# personal inversion-ancestry fragment.
#
# The DISTRIBUTION of per-carrier fragment boundaries across all carriers,
# aggregated on each side, gives:
#   - modal position  = population-level breakpoint estimate
#   - spread (bootstrap) = CI on the breakpoint
#   - tail mass = recombinant ancestry fraction
#
# REGISTRY WIRING (chat-18):
#   - Candidate coords from reg$intervals$get_candidate(cid)
#   - Karyotype groups from reg$samples$get_groups_for_candidate(cid)
#   - Block coords from reg$evidence$read_block(cid, "dosage_blocks")
#     (must run 01_dosage_signal.R first)
#   - Per-marker informative TSV read from raw/dosage_informative_markers.tsv.gz
#   - Output: per-candidate summary block via reg$evidence$write_block
#   - Per-sample fragments TSV stays in raw/ (large)
#
# OUTPUTS:
#   evidence_registry/per_candidate/<cid>/structured/ancestral_fragments_summary.json
#   evidence_registry/per_candidate/<cid>/raw/ancestral_fragments_per_sample.tsv.gz
#
# References:
#   compute_ancestral_fragments — ported from STEP_C01i
#
# Usage:
#   Rscript 02_ancestral_fragments.R [cid=all]
#   Rscript 02_ancestral_fragments.R LG28_1
#   Rscript 02_ancestral_fragments.R --config my_overrides.R LG28_1
# =============================================================================

# --- Entry: one line sources registry + ancestry stack -----------------------
Sys.setenv(CURRENT_SCRIPT = "02_ancestral_fragments.R")
.bridge <- Sys.getenv("REGISTRY_BRIDGE", "utils/registry_bridge.R")
if (!file.exists(.bridge)) {
  for (p in c("utils/registry_bridge.R",
              "../utils/registry_bridge.R",
              file.path(Sys.getenv("BASE", ""), "utils/registry_bridge.R"))) {
    if (file.exists(p)) { .bridge <- p; break }
  }
}
if (!file.exists(.bridge)) {
  stop("[02_ancestral_fragments] cannot locate utils/registry_bridge.R. ",
       "Set $REGISTRY_BRIDGE or $BASE and retry.")
}
source(.bridge)

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
config_file <- ""
cid_filter  <- NA_character_
.ai <- 1L
while (.ai <= length(args)) {
  a <- args[.ai]
  if (a == "--config" && .ai < length(args)) {
    config_file <- args[.ai + 1L]; .ai <- .ai + 2L
  } else if (a != "all") {
    cid_filter <- a; .ai <- .ai + 1L
  } else {
    .ai <- .ai + 1L
  }
}
if (nzchar(config_file) && file.exists(config_file)) source(config_file)

if (!exists("DOSAGE_DIR")) {
  DOSAGE_DIR <- Sys.getenv("DOSAGE_DIR", "")
  if (!nzchar(DOSAGE_DIR)) {
    DOSAGE_DIR <- file.path(BRIDGE_PATHS$BASE, "popstruct_thin",
                             "04_beagle_byRF_majmin")
  }
}
if (!exists("ensure_dir"))
  ensure_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); invisible(p) }

BP02_PARAMS <- list(
  rho_frag             = 0.5,     # per-sample correlation threshold
  max_scan_snps        = 5000L,   # max SNPs to scan past block edge per carrier
  min_valid_samples    = 10L,     # min non-NA samples at a SNP to test
  bootstrap_reps       = 1000L,   # bootstrap iterations for CI
  ci_level             = 0.95,    # 95% CI
  kde_bandwidth_method = "nrd0",  # Silverman's rule for density()
  min_carriers         = 5L       # need at least this many carriers to compute CI
)
if (exists("BP02_PARAMS_OVERRIDE") && is.list(BP02_PARAMS_OVERRIDE)) {
  for (k in names(BP02_PARAMS_OVERRIDE)) BP02_PARAMS[[k]] <- BP02_PARAMS_OVERRIDE[[k]]
}

# =============================================================================
# Per-carrier scan — port of C01i compute_ancestral_fragments
# =============================================================================
scan_carrier_fragments <- function(X, core_idx, carrier_idx,
                                     block_left, block_right,
                                     positions,
                                     rho_frag = 0.5, max_scan = 5000L,
                                     min_valid = 10L) {
  if (length(carrier_idx) == 0L) return(data.table())

  core_mat <- X[core_idx, , drop = FALSE]
  core_consensus <- colMeans(core_mat, na.rm = TRUE)

  n_total <- nrow(X)
  frag_rows <- vector("list", length(carrier_idx))

  test_cor_at <- function(mi) {
    m_gt <- X[mi, ]
    valid <- is.finite(m_gt) & is.finite(core_consensus)
    if (sum(valid) < min_valid) return(NA_real_)
    r <- suppressWarnings(cor(m_gt[valid], core_consensus[valid]))
    if (!is.finite(r)) NA_real_ else abs(r)
  }

  for (k in seq_along(carrier_idx)) {
    ci_col <- carrier_idx[k]
    # Left-side scan
    left_boundary <- block_left
    left_cor_at <- NA_real_
    if (block_left > 1L) {
      for (mi in seq(block_left - 1L, max(1L, block_left - max_scan), by = -1L)) {
        # Only proceed if this sample has a call at this marker
        if (!is.finite(X[mi, ci_col])) next
        r <- test_cor_at(mi)
        if (is.na(r) || r < rho_frag) { left_cor_at <- r; break }
        left_boundary <- mi
      }
    }

    # Right-side scan
    right_boundary <- block_right
    right_cor_at <- NA_real_
    if (block_right < n_total) {
      for (mi in seq(block_right + 1L, min(n_total, block_right + max_scan), by = 1L)) {
        if (!is.finite(X[mi, ci_col])) next
        r <- test_cor_at(mi)
        if (is.na(r) || r < rho_frag) { right_cor_at <- r; break }
        right_boundary <- mi
      }
    }

    frag_rows[[k]] <- list(
      sample_col_idx      = ci_col,
      frag_left_snp       = left_boundary,
      frag_right_snp      = right_boundary,
      frag_start_bp       = positions[left_boundary],
      frag_end_bp         = positions[right_boundary],
      frag_length_bp      = positions[right_boundary] - positions[left_boundary],
      extended_left_snps  = block_left  - left_boundary,
      extended_right_snps = right_boundary - block_right,
      left_cor_at_boundary  = left_cor_at,
      right_cor_at_boundary = right_cor_at
    )
  }
  rbindlist(lapply(frag_rows, as.data.table))
}

# =============================================================================
# Fragment-distribution summary: KDE mode + bootstrap CI
# =============================================================================
summarize_fragment_side <- function(boundary_positions,
                                      bootstrap_reps = 1000L,
                                      ci_level = 0.95) {
  boundary_positions <- boundary_positions[is.finite(boundary_positions)]
  n <- length(boundary_positions)
  if (n < BP02_PARAMS$min_carriers) {
    return(list(
      n                  = n,
      mode_bp            = NA_real_,
      ci_low             = NA_real_,
      ci_high            = NA_real_,
      ci_width_kb        = NA_real_,
      median_bp          = if (n > 0) median(boundary_positions) else NA_real_,
      mean_bp            = if (n > 0) mean(boundary_positions)   else NA_real_,
      mad_kb             = if (n > 0) mad(boundary_positions)/1000 else NA_real_,
      reason             = sprintf("only %d carriers (need >= %d)", n, BP02_PARAMS$min_carriers)
    ))
  }

  # KDE mode estimation. Use a density() call with Silverman's rule bandwidth
  # but computed on the main mass of the distribution (exclude extreme 10%
  # tails for bandwidth) to avoid the recombinant tail biasing bandwidth.
  q_lo <- quantile(boundary_positions, 0.05)
  q_hi <- quantile(boundary_positions, 0.95)
  core_mass <- boundary_positions[boundary_positions >= q_lo & boundary_positions <= q_hi]
  # Silverman bandwidth on core mass
  bw <- if (length(core_mass) >= 5) {
    tryCatch(bw.nrd0(core_mass), error = function(e) NA_real_)
  } else NA_real_
  if (!is.finite(bw) || bw <= 0) bw <- max(1, sd(boundary_positions, na.rm = TRUE) / 3)

  kde_fit <- function(x) {
    dens <- tryCatch(
      density(x, bw = bw, n = 2048),
      error = function(e) NULL
    )
    if (is.null(dens)) return(NA_real_)
    dens$x[which.max(dens$y)]
  }

  point_mode <- kde_fit(boundary_positions)

  # Bootstrap CI
  boot_modes <- numeric(bootstrap_reps)
  for (b in seq_len(bootstrap_reps)) {
    bs <- sample(boundary_positions, n, replace = TRUE)
    boot_modes[b] <- kde_fit(bs)
  }
  boot_modes <- boot_modes[is.finite(boot_modes)]
  if (length(boot_modes) < bootstrap_reps * 0.5) {
    return(list(
      n = n, mode_bp = point_mode,
      ci_low = NA_real_, ci_high = NA_real_, ci_width_kb = NA_real_,
      median_bp = median(boundary_positions),
      mean_bp = mean(boundary_positions),
      mad_kb = mad(boundary_positions) / 1000,
      reason = "bootstrap unstable"
    ))
  }
  alpha <- 1 - ci_level
  ci_low  <- quantile(boot_modes, alpha / 2)
  ci_high <- quantile(boot_modes, 1 - alpha / 2)

  list(
    n           = n,
    mode_bp     = as.numeric(point_mode),
    ci_low      = as.numeric(ci_low),
    ci_high     = as.numeric(ci_high),
    ci_width_kb = as.numeric((ci_high - ci_low) / 1000),
    median_bp   = median(boundary_positions),
    mean_bp     = mean(boundary_positions),
    mad_kb      = mad(boundary_positions) / 1000,
    reason      = "ok"
  )
}

# =============================================================================
# Process one candidate
# =============================================================================
process_candidate <- function(cid, chr) {
  # Read 01's block from the registry
  blocks <- tryCatch(reg$evidence$read_block(cid, "dosage_blocks"),
                     error = function(e) NULL)
  if (is.null(blocks) || is.null(blocks$data)) {
    message("[02_ancestral_fragments] cid=", cid,
            " SKIP — no dosage_blocks block (run 01_dosage_signal.R first)")
    return(invisible(NULL))
  }
  bdata <- blocks$data
  if (is.null(bdata$status) || bdata$status != "ok") {
    message("[02_ancestral_fragments] cid=", cid,
            " SKIP — dosage_blocks status=", bdata$status %||% "NA")
    return(invisible(NULL))
  }

  # Karyotype groups via registry
  grs <- tryCatch(reg$samples$get_groups_for_candidate(cid),
                  error = function(e) NULL)
  if (is.null(grs) ||
      length(grs$HOM_REF) == 0L ||
      length(grs$HET)     == 0L ||
      length(grs$HOM_INV) == 0L) {
    message("[02_ancestral_fragments] cid=", cid,
            " SKIP — karyotype groups not registered")
    return(invisible(NULL))
  }
  rot <- rbindlist(list(
    data.table(sample = grs$HOM_REF, coarse_group_refined = "HOMO_1"),
    data.table(sample = grs$HET,     coarse_group_refined = "HET"),
    data.table(sample = grs$HOM_INV, coarse_group_refined = "HOMO_2")
  ))

  dos_file   <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[02_ancestral_fragments] cid=", cid,
            " SKIP — dosage files missing under ", DOSAGE_DIR)
    return(invisible(NULL))
  }
  dos   <- fread(dos_file)
  sites <- fread(sites_file)

  # Same scan region as script 01 (flanked)
  c_start <- as.integer(bdata$input_start_bp)
  c_end   <- as.integer(bdata$input_end_bp)
  span <- c_end - c_start
  flank <- max(500000L, as.integer(span * 0.3))
  scan_start <- max(1L, c_start - flank)
  scan_end   <- c_end + flank

  keep <- which(sites$pos >= scan_start & sites$pos <= scan_end)
  dos_reg   <- dos[keep]
  sites_reg <- sites[keep]

  sc <- setdiff(names(dos_reg), "marker")
  if (all(grepl("^Ind", sc)) && length(sc) == nrow(rot)) {
    setnames(dos_reg, old = sc, new = rot$sample)
  }
  sample_cols <- intersect(rot$sample, setdiff(names(dos_reg), "marker"))

  X <- as.matrix(dos_reg[, ..sample_cols]); storage.mode(X) <- "double"
  groups <- rot[match(sample_cols, sample), coarse_group_refined]

  # Block indices within the SCAN region
  block_left  <- as.integer(bdata$ext_left_marker_idx)
  block_right <- as.integer(bdata$ext_right_marker_idx)

  # Load per-marker table from raw/ to derive core markers
  markers_tsv <- bdata$markers_tsv_path %||% file.path(
    BRIDGE_PATHS$REGISTRIES_ROOT, "data", "evidence_registry",
    "per_candidate", cid, "raw", "dosage_informative_markers.tsv.gz")

  if (!file.exists(markers_tsv)) {
    message("[02_ancestral_fragments] cid=", cid,
            " SKIP — markers TSV missing at ", markers_tsv)
    return(invisible(NULL))
  }
  info_dt <- fread(markers_tsv)
  info_idx_in_block <- info_dt[marker_idx >= as.integer(bdata$core_left_marker_idx) &
                                 marker_idx <= as.integer(bdata$core_right_marker_idx) &
                                 informative == TRUE, marker_idx]
  if (length(info_idx_in_block) < 5L) {
    info_idx_in_block <- seq(as.integer(bdata$core_left_marker_idx),
                              as.integer(bdata$core_right_marker_idx))
  }
  core_idx <- info_idx_in_block

  # Carriers: minority-karyotype + HET
  n_h1 <- sum(groups == "HOMO_1"); n_h2 <- sum(groups == "HOMO_2")
  minor_group <- if (n_h2 <= n_h1) "HOMO_2" else "HOMO_1"
  carrier_mask <- groups %in% c("HET", minor_group)
  carrier_idx <- which(carrier_mask)
  n_carriers <- length(carrier_idx)
  n_het     <- sum(groups[carrier_idx] == "HET")
  n_hom_inv <- sum(groups[carrier_idx] == minor_group)

  message("[02_ancestral_fragments] cid=", cid, " ", chr,
          "  block=", block_left, "-", block_right,
          "  carriers=", n_carriers, " (", minor_group, "+HET)",
          "  core_markers=", length(core_idx))

  if (n_carriers < BP02_PARAMS$min_carriers) {
    message("[02_ancestral_fragments]   too few carriers; writing status block")
    reg$evidence$write_block(cid, "ancestral_fragments_summary", list(
      status = "insufficient_carriers",
      reason = paste0("need >= ", BP02_PARAMS$min_carriers,
                      " carriers; got ", n_carriers),
      n_carriers = n_carriers,
      n_het = n_het,
      n_hom_inv = n_hom_inv,
      rho_frag = BP02_PARAMS$rho_frag,
      bootstrap_reps = BP02_PARAMS$bootstrap_reps,
      ci_level = BP02_PARAMS$ci_level
    ))
    return(invisible(NULL))
  }

  # --- Scan ----------------------------------------------------------------
  frag_dt <- scan_carrier_fragments(
    X, core_idx, carrier_idx,
    block_left = block_left, block_right = block_right,
    positions = sites_reg$pos,
    rho_frag = BP02_PARAMS$rho_frag,
    max_scan = BP02_PARAMS$max_scan_snps,
    min_valid = BP02_PARAMS$min_valid_samples
  )

  # Decorate with sample info
  frag_dt[, sample := sample_cols[sample_col_idx]]
  frag_dt[, coarse_group := groups[sample_col_idx]]
  frag_dt[, candidate_id := cid]
  frag_dt[, chrom := chr]
  setcolorder(frag_dt, c("candidate_id", "chrom", "sample", "coarse_group",
                          "frag_start_bp", "frag_end_bp", "frag_length_bp",
                          "extended_left_snps", "extended_right_snps",
                          "left_cor_at_boundary", "right_cor_at_boundary",
                          "frag_left_snp", "frag_right_snp", "sample_col_idx"))
  frag_dt[, left_cor_at_boundary := round(left_cor_at_boundary, 4)]
  frag_dt[, right_cor_at_boundary := round(right_cor_at_boundary, 4)]

  # Per-sample fragments TSV → raw/
  cand_raw_dir <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                             "evidence_registry", "per_candidate", cid, "raw")
  ensure_dir(cand_raw_dir)
  fragments_tsv <- file.path(cand_raw_dir, "ancestral_fragments_per_sample.tsv.gz")
  fwrite(frag_dt, fragments_tsv, sep = "\t", compress = "gzip")

  # --- Summarize distributions -------------------------------------------
  left_summary  <- summarize_fragment_side(frag_dt$frag_start_bp,
                                            bootstrap_reps = BP02_PARAMS$bootstrap_reps,
                                            ci_level = BP02_PARAMS$ci_level)
  right_summary <- summarize_fragment_side(frag_dt$frag_end_bp,
                                            bootstrap_reps = BP02_PARAMS$bootstrap_reps,
                                            ci_level = BP02_PARAMS$ci_level)

  # Decide overall status
  left_ok  <- isTRUE(left_summary$reason == "ok")
  right_ok <- isTRUE(right_summary$reason == "ok")
  overall_status <- if (left_ok && right_ok) "ok" else "bootstrap_failed"

  reg$evidence$write_block(cid, "ancestral_fragments_summary", list(
    status             = overall_status,
    reason             = if (overall_status != "ok") {
                            paste0("left=", left_summary$reason,
                                   "; right=", right_summary$reason)
                         } else NA_character_,
    n_carriers         = n_carriers,
    n_het              = n_het,
    n_hom_inv          = n_hom_inv,
    frag_left_bp_mode  = if (is.finite(left_summary$mode_bp))
                           as.integer(round(left_summary$mode_bp)) else NA_integer_,
    frag_right_bp_mode = if (is.finite(right_summary$mode_bp))
                           as.integer(round(right_summary$mode_bp)) else NA_integer_,
    frag_left_ci_low   = if (is.finite(left_summary$ci_low))
                           as.integer(round(left_summary$ci_low)) else NA_integer_,
    frag_left_ci_high  = if (is.finite(left_summary$ci_high))
                           as.integer(round(left_summary$ci_high)) else NA_integer_,
    frag_right_ci_low  = if (is.finite(right_summary$ci_low))
                           as.integer(round(right_summary$ci_low)) else NA_integer_,
    frag_right_ci_high = if (is.finite(right_summary$ci_high))
                           as.integer(round(right_summary$ci_high)) else NA_integer_,
    frag_left_ci_width_kb  = round(left_summary$ci_width_kb, 2),
    frag_right_ci_width_kb = round(right_summary$ci_width_kb, 2),
    frag_left_mad_kb       = round(left_summary$mad_kb, 2),
    frag_right_mad_kb      = round(right_summary$mad_kb, 2),
    rho_frag               = BP02_PARAMS$rho_frag,
    bootstrap_reps         = BP02_PARAMS$bootstrap_reps,
    ci_level               = BP02_PARAMS$ci_level,
    fragments_tsv_path     = fragments_tsv
  ))

  message("[02_ancestral_fragments]   left_mode=",
          format(round(left_summary$mode_bp),  big.mark = ","),
          " (CI ", round(left_summary$ci_width_kb),  " kb)",
          "  right_mode=",
          format(round(right_summary$mode_bp), big.mark = ","),
          " (CI ", round(right_summary$ci_width_kb), " kb)")

  invisible(list(fragments = frag_dt,
                  left = left_summary, right = right_summary))
}

# =============================================================================
# Driver
# =============================================================================
main <- function() {
  cand_path <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                          "interval_registry", "candidate_intervals.tsv")
  if (!file.exists(cand_path)) {
    stop("[02_ancestral_fragments] candidate_intervals.tsv not found at ",
         cand_path)
  }
  cand <- fread(cand_path)
  if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]
  if (nrow(cand) == 0) {
    message("[02_ancestral_fragments] No candidates to process",
            if (!is.na(cid_filter)) paste0(" (filter: ", cid_filter, ")") else "")
    return(invisible(NULL))
  }
  message("[02_ancestral_fragments] Processing ", nrow(cand),
          " candidate(s) from interval_registry")
  for (ci in seq_len(nrow(cand))) {
    row <- cand[ci]
    tryCatch(
      process_candidate(as.character(row$candidate_id),
                         as.character(row$chrom)),
      error = function(e) {
        message("[02_ancestral_fragments] cid=", row$candidate_id,
                " ERROR: ", conditionMessage(e))
      }
    )
  }
  message("[02_ancestral_fragments] DONE")
}

main()
