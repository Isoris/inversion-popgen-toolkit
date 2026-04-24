#!/usr/bin/env Rscript

# =============================================================================
# STEP37_breakpoint_support_from_sv.R  (v1.0)
#
# BREAKPOINT SUPPORT MODULE for inversion candidates using existing DELLY2
# and Manta VCF outputs.
#
# DESIGN: Simple, statistically interpretable breakpoint-support evidence.
#   - Does NOT require exact breakpoint-base precision
#   - Does NOT require all samples to agree exactly
#   - Uses breakpoint recurrence + sample-group association
#   - Support layer, not a breakpoint assembler
#
# For each inversion candidate region:
#   1. Parse/standardize SV calls from DELLY2 + Manta
#   2. Collect INV + BND calls near candidate boundary zones
#   3. Cluster nearby calls into breakpoint-support zones
#   4. Build sample × breakpoint-cluster support matrix
#   5. Test association with provisional sample groups (Fisher exact)
#   6. Assess left/right boundary concordance within samples
#
# INPUTS:
#   <config.R>                        — standard followup config
#   --delly_dir <path>                — directory with DELLY2 VCFs or parsed TSVs
#   --manta_dir <path>                — directory with Manta VCFs or parsed TSVs
#   --group_source <path_or_auto>     — sample group assignments (auto = from STEP21)
#   --boundary_bp <int>               — half-width of boundary zone (default 50000)
#   [cid=all]
#
# OUTPUTS (per candidate):
#   candidate_bp_sv_calls.tsv.gz       — standardized SV calls near boundaries
#   candidate_bp_clusters.tsv          — breakpoint-support clusters
#   candidate_bp_sample_matrix.tsv.gz  — sample × cluster support matrix
#   candidate_bp_group_tests.tsv       — Fisher exact tests per cluster × group
#   candidate_bp_concordance.tsv       — left/right boundary concordance per sample
#   candidate_bp_summary.tsv           — one-row candidate-level summary
#
# Usage:
#   Rscript STEP37_breakpoint_support_from_sv.R <config.R> [cid=all] \
#     --delly_dir <path> --manta_dir <path> \
#     [--boundary_bp 50000] [--cluster_bp 5000] [--min_support 2]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
cid_filter  <- NA_integer_

# Named argument parsing
delly_dir    <- NULL
manta_dir    <- NULL
group_source <- "auto"
BOUNDARY_BP  <- 50000L   # half-width of boundary zone
CLUSTER_BP   <- 5000L    # max distance to merge calls into one cluster
MIN_SUPPORT  <- 2L       # minimum calls in a cluster to keep
MIN_SAMPLES  <- 2L       # minimum unique samples in a cluster

i <- 2L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--delly_dir" && i < length(args)) {
    delly_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--manta_dir" && i < length(args)) {
    manta_dir <- args[i + 1]; i <- i + 2L
  } else if (a == "--group_source" && i < length(args)) {
    group_source <- args[i + 1]; i <- i + 2L
  } else if (a == "--boundary_bp" && i < length(args)) {
    BOUNDARY_BP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--cluster_bp" && i < length(args)) {
    CLUSTER_BP <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a == "--min_support" && i < length(args)) {
    MIN_SUPPORT <- as.integer(args[i + 1]); i <- i + 2L
  } else if (a != "all" && !grepl("^--", a)) {
    cid_filter <- suppressWarnings(as.integer(a))
    i <- i + 1L
  } else if (a == "all") {
    i <- i + 1L
  } else {
    i <- i + 1L
  }
}

source(config_file)
ensure_dir(FOLLOWUP_DIR)

cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

message("[STEP37] Breakpoint support module — ", nrow(cand), " candidates")
message("[STEP37] DELLY dir: ", delly_dir %||% "(not provided)")
message("[STEP37] Manta dir: ", manta_dir %||% "(not provided)")
message("[STEP37] Boundary zone: ±", BOUNDARY_BP, " bp")
message("[STEP37] Cluster tolerance: ", CLUSTER_BP, " bp")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# STEP 1 — PARSE AND STANDARDIZE SV CALLS
# =============================================================================

# Standard schema for all SV calls
# caller | sample | chrom1 | pos1 | chrom2 | pos2 | svtype | qual | filter |
# genotype | support_sr | support_pe | svlen

parse_delly_tsv <- function(dir_path, chroms_of_interest) {
  # Expects pre-parsed TSV files: <svtype>.filtered_calls.tsv.gz or similar
  # Format: chrom, pos, end, svtype, sample, genotype, qual, filter, PE, SR
  if (is.null(dir_path) || !dir.exists(dir_path)) return(data.table())

  files <- list.files(dir_path, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE,
                      recursive = TRUE)
  # Filter to INV and BND types
  inv_files <- files[grepl("INV|inv|BND|bnd", files, ignore.case = TRUE)]
  if (length(inv_files) == 0) {
    # Try generic files
    inv_files <- files
  }

  rows <- list()
  for (f in inv_files) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next

    # Standardize column names (handle various DELLY output formats)
    nms <- tolower(names(dt))
    names(dt) <- nms

    # Map common column names
    chrom_col <- intersect(nms, c("chrom", "chr", "chrom1", "#chrom"))[1]
    pos_col   <- intersect(nms, c("pos", "pos1", "start"))[1]
    end_col   <- intersect(nms, c("end", "pos2", "stop"))[1]
    type_col  <- intersect(nms, c("svtype", "type", "sv_type"))[1]
    samp_col  <- intersect(nms, c("sample", "sample_id", "sampleid"))[1]
    gt_col    <- intersect(nms, c("genotype", "gt", "geno"))[1]
    qual_col  <- intersect(nms, c("qual", "quality"))[1]
    filt_col  <- intersect(nms, c("filter", "filt"))[1]
    pe_col    <- intersect(nms, c("pe", "paired_end", "pe_support"))[1]
    sr_col    <- intersect(nms, c("sr", "split_read", "sr_support"))[1]

    if (is.na(chrom_col) || is.na(pos_col) || is.na(samp_col)) next

    std <- data.table(
      caller    = "DELLY",
      sample    = as.character(dt[[samp_col]]),
      chrom1    = as.character(dt[[chrom_col]]),
      pos1      = as.integer(dt[[pos_col]]),
      chrom2    = as.character(if (!is.na(end_col)) dt[[chrom_col]] else NA),
      pos2      = as.integer(if (!is.na(end_col)) dt[[end_col]] else NA),
      svtype    = as.character(if (!is.na(type_col)) dt[[type_col]] else "UNKNOWN"),
      qual      = as.numeric(if (!is.na(qual_col)) dt[[qual_col]] else NA),
      filter    = as.character(if (!is.na(filt_col)) dt[[filt_col]] else "UNKNOWN"),
      genotype  = as.character(if (!is.na(gt_col)) dt[[gt_col]] else NA),
      support_sr = as.integer(if (!is.na(sr_col)) dt[[sr_col]] else 0L),
      support_pe = as.integer(if (!is.na(pe_col)) dt[[pe_col]] else 0L)
    )

    # Filter to INV and BND only
    std <- std[toupper(svtype) %in% c("INV", "BND")]
    # Filter to chromosomes of interest
    if (length(chroms_of_interest) > 0) {
      std <- std[chrom1 %in% chroms_of_interest]
    }

    rows[[length(rows) + 1]] <- std
  }

  if (length(rows) == 0) return(data.table())
  rbindlist(rows, fill = TRUE)
}

parse_manta_tsv <- function(dir_path, chroms_of_interest) {
  # Same logic as DELLY but for Manta output format
  if (is.null(dir_path) || !dir.exists(dir_path)) return(data.table())

  files <- list.files(dir_path, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE,
                      recursive = TRUE)
  inv_files <- files[grepl("INV|inv|BND|bnd|sv|SV", files, ignore.case = TRUE)]
  if (length(inv_files) == 0) inv_files <- files

  rows <- list()
  for (f in inv_files) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next

    nms <- tolower(names(dt))
    names(dt) <- nms

    chrom_col <- intersect(nms, c("chrom", "chr", "chrom1", "#chrom"))[1]
    pos_col   <- intersect(nms, c("pos", "pos1", "start"))[1]
    end_col   <- intersect(nms, c("end", "pos2", "stop"))[1]
    type_col  <- intersect(nms, c("svtype", "type", "sv_type"))[1]
    samp_col  <- intersect(nms, c("sample", "sample_id"))[1]
    gt_col    <- intersect(nms, c("genotype", "gt"))[1]
    qual_col  <- intersect(nms, c("qual", "quality"))[1]
    filt_col  <- intersect(nms, c("filter", "filt"))[1]
    pe_col    <- intersect(nms, c("pe", "pr", "paired_end"))[1]
    sr_col    <- intersect(nms, c("sr", "split_read"))[1]

    if (is.na(chrom_col) || is.na(pos_col) || is.na(samp_col)) next

    std <- data.table(
      caller    = "Manta",
      sample    = as.character(dt[[samp_col]]),
      chrom1    = as.character(dt[[chrom_col]]),
      pos1      = as.integer(dt[[pos_col]]),
      chrom2    = as.character(if (!is.na(end_col)) dt[[chrom_col]] else NA),
      pos2      = as.integer(if (!is.na(end_col)) dt[[end_col]] else NA),
      svtype    = as.character(if (!is.na(type_col)) dt[[type_col]] else "UNKNOWN"),
      qual      = as.numeric(if (!is.na(qual_col)) dt[[qual_col]] else NA),
      filter    = as.character(if (!is.na(filt_col)) dt[[filt_col]] else "UNKNOWN"),
      genotype  = as.character(if (!is.na(gt_col)) dt[[gt_col]] else NA),
      support_sr = as.integer(if (!is.na(sr_col)) dt[[sr_col]] else 0L),
      support_pe = as.integer(if (!is.na(pe_col)) dt[[pe_col]] else 0L)
    )

    std <- std[toupper(svtype) %in% c("INV", "BND")]
    if (length(chroms_of_interest) > 0) {
      std <- std[chrom1 %in% chroms_of_interest]
    }
    rows[[length(rows) + 1]] <- std
  }

  if (length(rows) == 0) return(data.table())
  rbindlist(rows, fill = TRUE)
}

# =============================================================================
# STEP 2 — CLUSTER BREAKPOINT-SUPPORT CALLS
# =============================================================================

cluster_breakpoints <- function(calls_dt, cluster_bp = CLUSTER_BP) {
  # calls_dt must have: pos1, boundary_side, svtype_class
  # Returns calls_dt with cluster_id added
  if (nrow(calls_dt) == 0) return(calls_dt[, cluster_id := integer(0)])

  calls_dt <- calls_dt[order(boundary_side, pos1)]
  calls_dt[, cluster_id := NA_integer_]
  cid <- 0L

  for (side in c("left", "right")) {
    sub_idx <- which(calls_dt$boundary_side == side)
    if (length(sub_idx) == 0) next

    positions <- calls_dt$pos1[sub_idx]
    cid <- cid + 1L
    calls_dt$cluster_id[sub_idx[1]] <- cid

    for (k in seq_along(sub_idx)[-1]) {
      if (positions[k] - positions[k - 1] > cluster_bp) {
        cid <- cid + 1L
      }
      calls_dt$cluster_id[sub_idx[k]] <- cid
    }
  }

  calls_dt
}

summarize_clusters <- function(calls_dt) {
  # Produce per-cluster summary
  calls_dt[, .(
    boundary_side     = boundary_side[1],
    cluster_start     = min(pos1),
    cluster_end       = max(pos1),
    cluster_median    = as.integer(median(pos1)),
    cluster_span_bp   = max(pos1) - min(pos1),
    n_calls           = .N,
    n_unique_samples  = uniqueN(sample),
    n_inv             = sum(toupper(svtype) == "INV"),
    n_bnd             = sum(toupper(svtype) == "BND"),
    callers           = paste(sort(unique(caller)), collapse = ","),
    svtype_mode       = names(sort(table(toupper(svtype)), decreasing = TRUE))[1],
    mean_support_sr   = round(mean(support_sr, na.rm = TRUE), 1),
    mean_support_pe   = round(mean(support_pe, na.rm = TRUE), 1),
    samples           = paste(sort(unique(sample)), collapse = ";")
  ), by = cluster_id]
}

# =============================================================================
# STEP 3 — BUILD SAMPLE × CLUSTER SUPPORT MATRIX
# =============================================================================

build_support_matrix <- function(calls_dt, all_samples, mode = "binary") {
  # Returns a matrix: rows = samples, cols = cluster_ids
  clusters <- sort(unique(calls_dt$cluster_id))
  mat <- matrix(0, nrow = length(all_samples), ncol = length(clusters),
                dimnames = list(all_samples, paste0("C", clusters)))

  for (ci in seq_along(clusters)) {
    cl <- clusters[ci]
    sub <- calls_dt[cluster_id == cl]

    for (s in unique(sub$sample)) {
      si <- match(s, all_samples)
      if (is.na(si)) next

      if (mode == "binary") {
        mat[si, ci] <- 1L
      } else if (mode == "support") {
        s_calls <- sub[sample == s]
        mat[si, ci] <- sum(s_calls$support_sr + s_calls$support_pe, na.rm = TRUE)
      } else if (mode == "caller_count") {
        mat[si, ci] <- uniqueN(sub[sample == s]$caller)
      }
    }
  }

  mat
}

# =============================================================================
# STEP 4 — GROUP ASSOCIATION TESTING (Fisher exact)
# =============================================================================

test_group_association <- function(support_mat, group_assignments, cluster_summary) {
  # group_assignments: named vector, names = sample, values = group label
  # Returns one row per cluster × group combination
  results <- list()
  clusters <- colnames(support_mat)
  groups <- sort(unique(group_assignments))

  for (ci in seq_along(clusters)) {
    cl_name <- clusters[ci]
    has_support <- support_mat[, ci] > 0

    for (grp in groups) {
      in_group <- names(group_assignments)[group_assignments == grp]
      not_in_group <- names(group_assignments)[group_assignments != grp]

      # 2×2 table:
      #                 has_support  no_support
      # in_group            a            b
      # not_in_group        c            d
      a <- sum(has_support[names(has_support) %in% in_group])
      b <- sum(!has_support[names(has_support) %in% in_group])
      c <- sum(has_support[names(has_support) %in% not_in_group])
      d <- sum(!has_support[names(has_support) %in% not_in_group])

      # Fisher exact test
      ft <- tryCatch(
        fisher.test(matrix(c(a, c, b, d), nrow = 2)),
        error = function(e) list(p.value = NA, estimate = NA)
      )

      # Enrichment: proportion with support in-group vs out-of-group
      prop_in  <- if ((a + b) > 0) a / (a + b) else 0
      prop_out <- if ((c + d) > 0) c / (c + d) else 0

      cl_info <- cluster_summary[cluster_id == as.integer(gsub("^C", "", cl_name))]

      results[[length(results) + 1]] <- data.table(
        cluster_id     = cl_name,
        boundary_side  = if (nrow(cl_info) > 0) cl_info$boundary_side[1] else NA,
        group          = grp,
        n_in_group     = a + b,
        n_with_support_in_group = a,
        n_not_in_group = c + d,
        n_with_support_out_group = c,
        prop_in_group  = round(prop_in, 4),
        prop_out_group = round(prop_out, 4),
        enrichment     = round(prop_in - prop_out, 4),
        odds_ratio     = round(as.numeric(ft$estimate), 4),
        fisher_p       = signif(ft$p.value, 4),
        fisher_signif  = ft$p.value < 0.05
      )
    }
  }

  if (length(results) == 0) return(data.table())
  rbindlist(results)
}

# =============================================================================
# STEP 5 — LEFT/RIGHT CONCORDANCE PER SAMPLE
# =============================================================================

compute_concordance <- function(calls_dt, all_samples) {
  # For each sample: does it have support near BOTH left and right boundaries?
  left_samples  <- unique(calls_dt[boundary_side == "left"]$sample)
  right_samples <- unique(calls_dt[boundary_side == "right"]$sample)

  data.table(
    sample = all_samples,
    has_left_support  = all_samples %in% left_samples,
    has_right_support = all_samples %in% right_samples,
    concordant_both   = all_samples %in% left_samples & all_samples %in% right_samples,
    left_only         = all_samples %in% left_samples & !(all_samples %in% right_samples),
    right_only        = !(all_samples %in% left_samples) & all_samples %in% right_samples,
    neither           = !(all_samples %in% left_samples) & !(all_samples %in% right_samples)
  )
}

# =============================================================================
# LOAD SV CALLS
# =============================================================================

chroms_needed <- unique(cand$chrom)
message("[STEP37] Loading and standardizing SV calls...")

sv_delly <- parse_delly_tsv(delly_dir, chroms_needed)
sv_manta <- parse_manta_tsv(manta_dir, chroms_needed)
sv_all   <- rbindlist(list(sv_delly, sv_manta), fill = TRUE)

if (nrow(sv_all) > 0) {
  sv_all[, support_total := support_sr + support_pe]
}

message("[STEP37] DELLY calls: ", nrow(sv_delly),
        " (INV: ", sum(toupper(sv_delly$svtype) == "INV", na.rm = TRUE),
        ", BND: ", sum(toupper(sv_delly$svtype) == "BND", na.rm = TRUE), ")")
message("[STEP37] Manta calls: ", nrow(sv_manta),
        " (INV: ", sum(toupper(sv_manta$svtype) == "INV", na.rm = TRUE),
        ", BND: ", sum(toupper(sv_manta$svtype) == "BND", na.rm = TRUE), ")")

# =============================================================================
# LOAD SAMPLE NAMES
# =============================================================================

sample_names <- tryCatch({
  s <- fread(SAMPLE_IND_FILE, header = FALSE)[[1]]
  s[nchar(s) > 0]
}, error = function(e) NULL)
if (is.null(sample_names)) stop("Cannot read sample names from: ", SAMPLE_IND_FILE)

# =============================================================================
# MAIN LOOP — PER CANDIDATE
# =============================================================================

all_summaries <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom
  c_start <- as.numeric(row$start_bp)
  c_end   <- as.numeric(row$end_bp)

  cand_prefix <- paste0(chr, "_cand", cid, "_", sprintf("%.2fMb", c_start / 1e6))
  cand_dir <- file.path(FOLLOWUP_DIR, cand_prefix)
  ensure_dir(cand_dir)

  message("\n[STEP37] ═══════ Candidate ", cid, " (", chr, ":",
          round(c_start / 1e6, 2), "–", round(c_end / 1e6, 2), " Mb) ═══════")

  # ── Define boundary zones ──────────────────────────────────────────
  left_zone_start  <- c_start - BOUNDARY_BP
  left_zone_end    <- c_start + BOUNDARY_BP
  right_zone_start <- c_end - BOUNDARY_BP
  right_zone_end   <- c_end + BOUNDARY_BP

  # ── Collect SV calls near boundaries ───────────────────────────────
  if (nrow(sv_all) == 0) {
    message("[STEP37] No SV calls available — skipping")
    all_summaries[[ci]] <- data.table(
      candidate_id = cid, chrom = chr, n_sv_calls = 0L,
      n_clusters = 0L, n_left_clusters = 0L, n_right_clusters = 0L,
      has_breakpoint_support = FALSE
    )
    next
  }

  # Match calls to boundary zones
  sv_chr <- sv_all[chrom1 == chr]
  if (nrow(sv_chr) == 0) {
    message("[STEP37] No SV calls on ", chr)
    all_summaries[[ci]] <- data.table(
      candidate_id = cid, chrom = chr, n_sv_calls = 0L,
      n_clusters = 0L, n_left_clusters = 0L, n_right_clusters = 0L,
      has_breakpoint_support = FALSE
    )
    next
  }

  # Assign boundary side
  sv_chr[, boundary_side := NA_character_]
  sv_chr[pos1 >= left_zone_start & pos1 <= left_zone_end, boundary_side := "left"]
  sv_chr[pos1 >= right_zone_start & pos1 <= right_zone_end, boundary_side := "right"]

  # Also check pos2 for INV calls (the other end of the inversion)
  sv_chr[!is.na(pos2) & pos2 >= left_zone_start & pos2 <= left_zone_end &
           is.na(boundary_side), boundary_side := "left"]
  sv_chr[!is.na(pos2) & pos2 >= right_zone_start & pos2 <= right_zone_end &
           is.na(boundary_side), boundary_side := "right"]

  # Keep only calls assigned to a boundary
  calls <- sv_chr[!is.na(boundary_side)]
  calls[, candidate_id := cid]

  message("[STEP37] Calls near boundaries: ", nrow(calls),
          " (left: ", sum(calls$boundary_side == "left"),
          ", right: ", sum(calls$boundary_side == "right"), ")")

  if (nrow(calls) == 0) {
    all_summaries[[ci]] <- data.table(
      candidate_id = cid, chrom = chr, n_sv_calls = 0L,
      n_clusters = 0L, n_left_clusters = 0L, n_right_clusters = 0L,
      has_breakpoint_support = FALSE
    )
    next
  }

  # ── Write standardized calls ───────────────────────────────────────
  fwrite(calls, file.path(cand_dir, "candidate_bp_sv_calls.tsv.gz"), sep = "\t")

  # ── Cluster calls ──────────────────────────────────────────────────
  # Run 3 modes: INV-only, BND-only, combined
  modes <- list(
    INV_only  = calls[toupper(svtype) == "INV"],
    BND_only  = calls[toupper(svtype) == "BND"],
    combined  = calls
  )

  cluster_results <- list()
  for (mode_name in names(modes)) {
    mode_calls <- modes[[mode_name]]
    if (nrow(mode_calls) < MIN_SUPPORT) next

    mode_calls <- cluster_breakpoints(copy(mode_calls), cluster_bp = CLUSTER_BP)
    cl_summary <- summarize_clusters(mode_calls)

    # Filter clusters by minimum support
    cl_keep <- cl_summary[n_calls >= MIN_SUPPORT & n_unique_samples >= MIN_SAMPLES]
    if (nrow(cl_keep) == 0) next

    mode_calls <- mode_calls[cluster_id %in% cl_keep$cluster_id]
    cl_keep[, mode := mode_name]
    cluster_results[[mode_name]] <- list(calls = mode_calls, summary = cl_keep)
  }

  # Use combined mode as primary (fall back to INV_only or BND_only)
  primary_mode <- if ("combined" %in% names(cluster_results)) "combined" else
    if ("INV_only" %in% names(cluster_results)) "INV_only" else
    if ("BND_only" %in% names(cluster_results)) "BND_only" else NULL

  if (is.null(primary_mode)) {
    message("[STEP37] No clusters passed filters")
    all_summaries[[ci]] <- data.table(
      candidate_id = cid, chrom = chr, n_sv_calls = nrow(calls),
      n_clusters = 0L, n_left_clusters = 0L, n_right_clusters = 0L,
      has_breakpoint_support = FALSE
    )
    next
  }

  primary <- cluster_results[[primary_mode]]
  cl_summary <- primary$summary
  cl_calls   <- primary$calls

  # Append all mode summaries
  all_cl_summaries <- rbindlist(lapply(cluster_results, function(x) x$summary), fill = TRUE)
  all_cl_summaries[, candidate_id := cid]
  fwrite(all_cl_summaries, file.path(cand_dir, "candidate_bp_clusters.tsv"), sep = "\t")

  message("[STEP37] Clusters (", primary_mode, "): ", nrow(cl_summary),
          " (left: ", sum(cl_summary$boundary_side == "left"),
          ", right: ", sum(cl_summary$boundary_side == "right"), ")")

  # ── Build sample × cluster support matrix ──────────────────────────
  mat_binary  <- build_support_matrix(cl_calls, sample_names, mode = "binary")
  mat_support <- build_support_matrix(cl_calls, sample_names, mode = "support")

  # Write binary matrix
  mat_dt <- as.data.table(mat_binary)
  mat_dt[, sample := sample_names]
  setcolorder(mat_dt, c("sample", setdiff(names(mat_dt), "sample")))
  mat_dt[, candidate_id := cid]
  fwrite(mat_dt, file.path(cand_dir, "candidate_bp_sample_matrix.tsv.gz"), sep = "\t")

  # ── Load group assignments ─────────────────────────────────────────
  group_dt <- NULL

  if (group_source == "auto") {
    # Try STEP21 output first
    g_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
    if (file.exists(g_file)) {
      gdt <- tryCatch(fread(g_file), error = function(e) NULL)
      if (!is.null(gdt) && "group" %in% names(gdt) && "sample" %in% names(gdt)) {
        group_dt <- gdt[, .(sample, group)]
      }
    }
    # Fall back to STEP12 groups
    if (is.null(group_dt)) {
      g_file2 <- file.path(cand_dir, "candidate_group_assignments.tsv")
      if (file.exists(g_file2)) {
        gdt2 <- tryCatch(fread(g_file2), error = function(e) NULL)
        if (!is.null(gdt2)) {
          grp_col <- intersect(names(gdt2), c("group", "cluster", "assignment",
                                                "state", "inversion_group"))
          if (length(grp_col) > 0) {
            samp_col <- intersect(names(gdt2), c("sample", "sample_id"))
            if (length(samp_col) > 0) {
              group_dt <- gdt2[, .SD, .SDcols = c(samp_col[1], grp_col[1])]
              setnames(group_dt, c("sample", "group"))
            }
          }
        }
      }
    }
  } else if (file.exists(group_source)) {
    group_dt <- tryCatch(fread(group_source), error = function(e) NULL)
    if (!is.null(group_dt)) {
      grp_col <- intersect(names(group_dt), c("group", "cluster", "assignment",
                                                "provisional_group"))
      samp_col <- intersect(names(group_dt), c("sample", "sample_id"))
      if (length(grp_col) > 0 && length(samp_col) > 0) {
        group_dt <- group_dt[, .SD, .SDcols = c(samp_col[1], grp_col[1])]
        setnames(group_dt, c("sample", "group"))
      } else {
        group_dt <- NULL
      }
    }
  }

  # ── Group association tests ────────────────────────────────────────
  if (!is.null(group_dt) && nrow(group_dt) > 0) {
    # Build named vector
    group_vec <- setNames(as.character(group_dt$group), group_dt$sample)
    # Keep only samples in our matrix
    group_vec <- group_vec[names(group_vec) %in% sample_names]

    if (length(group_vec) >= 10 && length(unique(group_vec)) >= 2) {
      test_dt <- test_group_association(mat_binary, group_vec, cl_summary)
      test_dt[, candidate_id := cid]
      fwrite(test_dt, file.path(cand_dir, "candidate_bp_group_tests.tsv"), sep = "\t")

      n_sig <- sum(test_dt$fisher_signif, na.rm = TRUE)
      message("[STEP37] Group association tests: ", nrow(test_dt),
              " (", n_sig, " significant at p<0.05)")
    } else {
      message("[STEP37] Too few samples/groups for association testing")
      test_dt <- data.table()
    }
  } else {
    message("[STEP37] No group assignments found — skipping association tests")
    test_dt <- data.table()
  }

  # ── Left/right concordance ─────────────────────────────────────────
  conc <- compute_concordance(cl_calls, sample_names)
  if (!is.null(group_dt) && nrow(group_dt) > 0) {
    conc <- merge(conc, group_dt, by = "sample", all.x = TRUE)
  }
  conc[, candidate_id := cid]
  fwrite(conc, file.path(cand_dir, "candidate_bp_concordance.tsv"), sep = "\t")

  n_concordant <- sum(conc$concordant_both)
  n_left_only  <- sum(conc$left_only)
  n_right_only <- sum(conc$right_only)
  message("[STEP37] Concordance: both=", n_concordant,
          " left_only=", n_left_only,
          " right_only=", n_right_only)

  # ── Candidate-level summary ────────────────────────────────────────
  best_enrichment <- NA_real_
  best_enrichment_group <- NA_character_
  best_fisher_p <- NA_real_
  if (nrow(test_dt) > 0 && any(test_dt$fisher_signif, na.rm = TRUE)) {
    best_row <- test_dt[which.min(fisher_p)]
    best_enrichment <- best_row$enrichment
    best_enrichment_group <- best_row$group
    best_fisher_p <- best_row$fisher_p
  }

  all_summaries[[ci]] <- data.table(
    candidate_id           = cid,
    chrom                  = chr,
    n_sv_calls             = nrow(calls),
    n_inv_calls            = sum(toupper(calls$svtype) == "INV"),
    n_bnd_calls            = sum(toupper(calls$svtype) == "BND"),
    n_clusters             = nrow(cl_summary),
    n_left_clusters        = sum(cl_summary$boundary_side == "left"),
    n_right_clusters       = sum(cl_summary$boundary_side == "right"),
    n_samples_with_support = sum(rowSums(mat_binary) > 0),
    n_concordant_both      = n_concordant,
    n_left_only            = n_left_only,
    n_right_only           = n_right_only,
    primary_mode           = primary_mode,
    best_enrichment        = round(best_enrichment, 4),
    best_enrichment_group  = best_enrichment_group,
    best_fisher_p          = best_fisher_p,
    has_breakpoint_support = nrow(cl_summary) > 0 && n_concordant >= MIN_SAMPLES
  )

  fwrite(all_summaries[[ci]],
         file.path(cand_dir, "candidate_bp_summary.tsv"), sep = "\t")
}

# =============================================================================
# GLOBAL SUMMARY
# =============================================================================

global_summary <- if (length(all_summaries) > 0) rbindlist(all_summaries, fill = TRUE) else {
  data.table(candidate_id = integer(), has_breakpoint_support = logical())
}

f_global <- file.path(FOLLOWUP_DIR, "breakpoint_support_global_summary.tsv")
fwrite(global_summary, f_global, sep = "\t")

message("\n[DONE] STEP37 breakpoint support module complete")
message("  Global summary: ", f_global)
message("  Candidates with support: ",
        sum(global_summary$has_breakpoint_support, na.rm = TRUE),
        " / ", nrow(global_summary))

if (nrow(global_summary) > 0) {
  message("\n[STEP37] Per-candidate breakdown:")
  for (i in seq_len(nrow(global_summary))) {
    r <- global_summary[i]
    flag <- if (r$has_breakpoint_support) "✓" else "·"
    message("  ", flag, " Candidate ", r$candidate_id, " (", r$chrom, "): ",
            r$n_sv_calls, " calls, ", r$n_clusters, " clusters, ",
            r$n_concordant_both, " concordant samples",
            if (!is.na(r$best_fisher_p) && r$best_fisher_p < 0.05)
              paste0(" [enriched in ", r$best_enrichment_group,
                     ", p=", signif(r$best_fisher_p, 3), "]") else "")
  }
}
