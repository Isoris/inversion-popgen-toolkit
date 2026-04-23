#!/usr/bin/env Rscript
# =============================================================================
# q10_register.R
# =============================================================================
# Bridge: takes phase_qc_shelf outputs for ONE chromosome and writes them into
# the four registries (sample / interval / results / evidence).
#
# Design principles:
#   - Idempotent: re-running with same inputs yields same registry state.
#     If cid already exists, update the evidence bundle and refresh the result
#     rows for this method, leaving other methods' rows untouched.
#   - Fault-tolerant: any registry API loader that fails to source is fatal
#     (can't silently no-op — the user explicitly asked to wire everything).
#   - Provenance: every write records this script's name, timestamp, sha256 of
#     key inputs (precomp, invgt_file) so any downstream query knows where the
#     numbers came from.
#
# CID convention:
#   <prefix>_<CHR>_<START_MB_int>_<END_MB_int>
#   e.g.   inv_C_gar_LG28_15_18
#
# Sample group convention (per chat-16 design):
#   inv_<cid>_HOM1   inv_<cid>_HET   inv_<cid>_HOM2
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)[1]
  if (is.na(i) || i >= length(args)) return(default)
  args[i + 1]
}

CHR              <- get_arg("--chrom")
REGISTRIES_ROOT  <- get_arg("--registries_root")
REGISTRY_API_R   <- get_arg("--registry_api_r")
PRECOMP          <- get_arg("--precomp")
INVGT_FILE       <- get_arg("--invgt_file")
POPSTATS_FILE    <- get_arg("--popstats_file")
HOBS_FILE        <- get_arg("--hobs_file", NA_character_)
GAP_FILE         <- get_arg("--gap_file", NA)
SAMPLE_GROUP     <- get_arg("--sample_group", "all_226")
METHOD_TAG       <- get_arg("--method_tag", "phase_qc_shelf")
CID_PREFIX       <- get_arg("--cid_prefix", "inv")
FORCE_OVERWRITE  <- as.integer(get_arg("--force_overwrite", "0")) == 1L
SHELF_A          <- as.numeric(get_arg("--shelf_start_mb", NA))
SHELF_B          <- as.numeric(get_arg("--shelf_end_mb",   NA))
BP1_MB           <- as.numeric(get_arg("--breakpoint1_mb", NA))
BP2_MB           <- as.numeric(get_arg("--breakpoint2_mb", NA))
EVIDENCE_FILES   <- get_arg("--evidence_files", "")

# Resolve breakpoints — prefer explicit BP, fall back to shelf bounds
bp_start_mb <- if (is.finite(BP1_MB)) BP1_MB else SHELF_A
bp_end_mb   <- if (is.finite(BP2_MB)) BP2_MB else SHELF_B
if (!is.finite(bp_start_mb) || !is.finite(bp_end_mb)) {
  stop("q10: breakpoint / shelf coordinates required")
}
start_bp <- as.integer(round(bp_start_mb * 1e6))
end_bp   <- as.integer(round(bp_end_mb   * 1e6))

# cid = <prefix>_<CHR>_<START_MB_int>_<END_MB_int>
cid <- sprintf("%s_%s_%d_%d", CID_PREFIX, CHR,
               as.integer(floor(bp_start_mb)),
               as.integer(floor(bp_end_mb)))

message("[q10] ", CHR, " -> cid=", cid)
message("[q10]   breakpoints: ", bp_start_mb, " - ", bp_end_mb, " Mb")
message("[q10]   sample_group: ", SAMPLE_GROUP)

# -----------------------------------------------------------------------------
# Load registry APIs. Each registry has a loader in REGISTRY_API_R/
# (sample_registry.R, interval_registry.R, evidence_registry.R,
# results_registry.R). If some aren't present yet (e.g. results_registry was
# built in chat-16 but not all chats had it), we degrade gracefully.
# -----------------------------------------------------------------------------
api_files <- c(
  sample_registry   = file.path(REGISTRY_API_R, "sample_registry.R"),
  interval_registry = file.path(REGISTRY_API_R, "interval_registry.R"),
  evidence_registry = file.path(REGISTRY_API_R, "evidence_registry.R"),
  results_registry  = file.path(REGISTRY_API_R, "results_registry.R")
)
apis_loaded <- character()
for (name in names(api_files)) {
  f <- api_files[[name]]
  if (file.exists(f)) {
    tryCatch({ source(f); apis_loaded <- c(apis_loaded, name) },
             error = function(e) message("[q10]   WARN: failed to source ", f, ": ", conditionMessage(e)))
  } else {
    message("[q10]   NOTE: API file not found: ", f, " (will skip this registry)")
  }
}
if (length(apis_loaded) == 0L) {
  stop("[q10] no registry APIs could be loaded; check REGISTRY_API_R=", REGISTRY_API_R)
}
message("[q10]   loaded APIs: ", paste(apis_loaded, collapse = ", "))

# Data dir convention: <REGISTRIES>/data/<registry_name>/
data_dir <- function(registry) {
  file.path(REGISTRIES_ROOT, "data", registry)
}

# -----------------------------------------------------------------------------
# Provenance helper
# -----------------------------------------------------------------------------
file_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}
ts_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

provenance <- list(
  script       = "phase_qc_shelf/R/q10_register.R",
  method       = METHOD_TAG,
  written_at   = ts_now(),
  precomp_sha  = file_sha256(PRECOMP),
  invgt_sha    = file_sha256(INVGT_FILE),
  popstats_sha = file_sha256(POPSTATS_FILE)
)

# -----------------------------------------------------------------------------
# 1. INTERVAL REGISTRY
# -----------------------------------------------------------------------------
if ("interval_registry" %in% apis_loaded) {
  message("[q10] 1/4: interval_registry")
  interval_reg <- load_interval_registry(data_dir("interval_registry"))

  interval_row <- list(
    cid          = cid,
    chrom        = CHR,
    start_bp     = start_bp,
    end_bp       = end_bp,
    span_bp      = end_bp - start_bp,
    method       = METHOD_TAG,
    confidence   = "high",
    breakpoint1_mb = bp_start_mb,
    breakpoint2_mb = bp_end_mb,
    shelf_start_mb = SHELF_A,
    shelf_end_mb   = SHELF_B,
    provenance   = provenance
  )

  if (interval_reg$has(cid) && !FORCE_OVERWRITE) {
    message("[q10]   cid exists, updating in place")
    interval_reg$update(cid, interval_row)
  } else {
    interval_reg$add(interval_row)
  }
  message("[q10]   registered interval: ", cid,
          " (", sprintf("%.3f-%.3f Mb", bp_start_mb, bp_end_mb), ")")
}

# -----------------------------------------------------------------------------
# 2. SAMPLE REGISTRY: karyotype groups Hom1 / Het / Hom2
# -----------------------------------------------------------------------------
if ("sample_registry" %in% apis_loaded) {
  message("[q10] 2/4: sample_registry")
  sample_reg <- load_sample_registry(data_dir("sample_registry"))

  inv <- fread(INVGT_FILE)
  # Expect columns sample_id + invgt (or similar). Accept a few aliases.
  sid_col  <- intersect(c("sample_id", "sample", "ind"),        names(inv))[1]
  grp_col  <- intersect(c("invgt", "group", "karyotype", "assignment"), names(inv))[1]
  if (is.null(sid_col) || is.null(grp_col)) {
    stop("[q10] invgt_file missing expected cols: found ", paste(names(inv), collapse=","))
  }
  setnames(inv, c(sid_col, grp_col), c("sample_id", "invgt"))
  inv[, invgt := toupper(as.character(invgt))]
  # Normalize labels to HOM1 / HET / HOM2
  inv[invgt %in% c("HOM1", "HOM_REF", "HOMOZYGOUS_REF", "REF"),       invgt := "HOM1"]
  inv[invgt %in% c("HET",  "HETEROZYGOUS"),                            invgt := "HET"]
  inv[invgt %in% c("HOM2", "HOM_ALT", "HOMOZYGOUS_ALT", "INV"),        invgt := "HOM2"]

  groups_registered <- 0L
  for (k in c("HOM1", "HET", "HOM2")) {
    ids <- inv[invgt == k, sample_id]
    if (length(ids) == 0L) next
    gid <- sprintf("inv_%s_%s", cid, k)
    if (sample_reg$has_group(gid) && !FORCE_OVERWRITE) {
      message("[q10]   ", gid, ": already registered (", length(ids), " samples), skipping")
      next
    }
    sample_reg$add_group(
      group_id     = gid,
      sample_ids   = ids,
      chrom        = CHR,
      inv_id       = cid,
      subgroup     = k,
      dimension    = "karyotype",
      parent_group = SAMPLE_GROUP,
      description  = sprintf("%s %s (%s method, %.3f-%.3f Mb)",
                             cid, k, METHOD_TAG, bp_start_mb, bp_end_mb)
    )
    groups_registered <- groups_registered + 1L
    message("[q10]   ", gid, ": n=", length(ids))
  }
  message("[q10]   registered ", groups_registered, " karyotype groups")
}

# -----------------------------------------------------------------------------
# 3. RESULTS REGISTRY: pairwise Fst Hom1-vs-Hom2 + summary stats
# -----------------------------------------------------------------------------
if ("results_registry" %in% apis_loaded) {
  message("[q10] 3/4: results_registry")
  results_reg <- load_results_registry(data_dir("results_registry"))

  ps <- fread(POPSTATS_FILE)
  setnames(ps, tolower(names(ps)))
  n_used_col <- intersect(c("n_sites_used", "n_used", "nsites"), names(ps))[1]
  if (!is.null(n_used_col)) ps <- ps[get(n_used_col) >= 20]
  start_col <- intersect(c("start_bp", "window_start", "start"), names(ps))[1]
  end_col   <- intersect(c("end_bp",   "window_end",   "end"),   names(ps))[1]
  ps[, mb := (get(start_col) + get(end_col)) / 2 / 1e6]

  inside  <- ps[mb >= bp_start_mb & mb <= bp_end_mb]
  outside <- ps[mb <  bp_start_mb | mb >  bp_end_mb]

  fst_col <- grep("^(fst|hudson_fst)_hom1.*hom2|hom2.*hom1", names(ps),
                  value = TRUE, ignore.case = TRUE)[1]
  pi_cols <- grep("^theta_pi_(hom1|het|hom2)$", names(ps), value = TRUE)

  summary_stats <- list(
    n_windows_inside  = nrow(inside),
    n_windows_outside = nrow(outside),
    fst_hom1_hom2_mean_inside  = if (!is.na(fst_col)) mean(inside[[fst_col]],  na.rm = TRUE) else NA_real_,
    fst_hom1_hom2_mean_outside = if (!is.na(fst_col)) mean(outside[[fst_col]], na.rm = TRUE) else NA_real_,
    fst_enrichment_fold        = NA_real_
  )
  if (is.finite(summary_stats$fst_hom1_hom2_mean_outside) &&
      summary_stats$fst_hom1_hom2_mean_outside > 0) {
    summary_stats$fst_enrichment_fold <-
      summary_stats$fst_hom1_hom2_mean_inside /
      summary_stats$fst_hom1_hom2_mean_outside
  }
  for (pcol in pi_cols) {
    summary_stats[[paste0("pi_inside_",  pcol)]] <- mean(inside[[pcol]],  na.rm = TRUE)
    summary_stats[[paste0("pi_outside_", pcol)]] <- mean(outside[[pcol]], na.rm = TRUE)
  }

  # 3a. Pairwise Fst profile — full per-window table for the chromosome.
  # Written under kind="pairwise", stat="fst", who="hom1_vs_hom2", where=cid.
  fst_profile_file <- file.path(tempdir(), sprintf("fst_%s.tsv", cid))
  if (!is.na(fst_col)) {
    out_df <- ps[, .(chrom = CHR,
                     start_bp = get(start_col),
                     end_bp   = get(end_col),
                     fst      = get(fst_col))]
    fwrite(out_df, fst_profile_file, sep = "\t")
    results_reg$put_pairwise(
      kind    = "pairwise",
      stat    = "fst",
      who     = "hom1_vs_hom2",
      where   = cid,
      chrom   = CHR,
      file    = fst_profile_file,
      sample_group = SAMPLE_GROUP,
      method  = METHOD_TAG,
      provenance = provenance
    )
    message("[q10]   registered Fst profile (", nrow(out_df), " windows)")
  }

  # 3b. Summary stats row — kind="interval_summary", stat="fst_summary"
  results_reg$put_interval_summary(
    cid          = cid,
    stat         = "fst_summary",
    values       = summary_stats,
    chrom        = CHR,
    sample_group = SAMPLE_GROUP,
    method       = METHOD_TAG,
    provenance   = provenance
  )
  message("[q10]   summary: Fst_in=", sprintf("%.3f", summary_stats$fst_hom1_hom2_mean_inside),
          " Fst_out=", sprintf("%.3f", summary_stats$fst_hom1_hom2_mean_outside),
          " fold=",    sprintf("%.1fx", summary_stats$fst_enrichment_fold))

  # 3c. Engine H: per-group HoverE profile + summary (Merot-style HWE)
  if (!is.na(HOBS_FILE) && file.exists(HOBS_FILE)) {
    hh <- fread(HOBS_FILE)
    setnames(hh, tolower(names(hh)))
    start_h <- intersect(c("start_bp", "start"), names(hh))[1]
    end_h   <- intersect(c("end_bp",   "end"),   names(hh))[1]
    hh[, mb := (get(start_h) + get(end_h)) / 2 / 1e6]
    h_inside  <- hh[mb >= bp_start_mb & mb <= bp_end_mb]
    h_outside <- hh[mb <  bp_start_mb | mb >  bp_end_mb]
    hobs_summary <- list()
    for (g in c("hom1", "het", "hom2")) {
      col <- paste0("hovere_", g)
      if (col %in% names(hh)) {
        hobs_summary[[paste0("hovere_", g, "_inside")]]  <- mean(h_inside[[col]],  na.rm = TRUE)
        hobs_summary[[paste0("hovere_", g, "_outside")]] <- mean(h_outside[[col]], na.rm = TRUE)
      }
    }
    # Full profile
    hobs_profile_file <- file.path(tempdir(), sprintf("hobs_%s.tsv", cid))
    fwrite(hh, hobs_profile_file, sep = "\t")
    results_reg$put_pairwise(
      kind    = "group_hwe",
      stat    = "hobs_hexp",
      who     = "hom1_het_hom2",
      where   = cid,
      chrom   = CHR,
      file    = hobs_profile_file,
      sample_group = SAMPLE_GROUP,
      method  = METHOD_TAG,
      provenance = provenance
    )
    results_reg$put_interval_summary(
      cid          = cid,
      stat         = "hobs_summary",
      values       = hobs_summary,
      chrom        = CHR,
      sample_group = SAMPLE_GROUP,
      method       = METHOD_TAG,
      provenance   = provenance
    )
    message("[q10]   Engine H summary: HoverE_Het_in=",
            sprintf("%.2f", hobs_summary$hovere_het_inside %||% NA),
            " HoverE_Hom1_in=",
            sprintf("%.2f", hobs_summary$hovere_hom1_inside %||% NA),
            " HoverE_Hom2_in=",
            sprintf("%.2f", hobs_summary$hovere_hom2_inside %||% NA))
  }
}

# -----------------------------------------------------------------------------
# 4. EVIDENCE REGISTRY: figure bundle + tables
# -----------------------------------------------------------------------------
if ("evidence_registry" %in% apis_loaded) {
  message("[q10] 4/4: evidence_registry")
  evidence_reg <- load_evidence_registry(data_dir("evidence_registry"))

  ev_files <- if (nzchar(EVIDENCE_FILES))
                strsplit(EVIDENCE_FILES, ":", fixed = TRUE)[[1]] else character(0)
  ev_files <- ev_files[file.exists(ev_files)]

  # Write a summary.json containing the interval definition + stats
  sum_json <- list(
    cid                = cid,
    chrom              = CHR,
    breakpoint1_mb     = bp_start_mb,
    breakpoint2_mb     = bp_end_mb,
    method             = METHOD_TAG,
    shelf_start_mb     = SHELF_A,
    shelf_end_mb       = SHELF_B,
    sample_group       = SAMPLE_GROUP,
    karyotype_groups   = sprintf("inv_%s_{HOM1,HET,HOM2}", cid),
    provenance         = provenance
  )
  if (exists("summary_stats")) sum_json$summary_stats <- summary_stats
  sum_path <- file.path(tempdir(), sprintf("summary_%s.json", cid))
  writeLines(toJSON(sum_json, pretty = TRUE, auto_unbox = TRUE), sum_path)

  ev_files <- c(ev_files, sum_path)

  evidence_reg$add_candidate(
    cid         = cid,
    chrom       = CHR,
    start_bp    = start_bp,
    end_bp      = end_bp,
    method      = METHOD_TAG,
    files       = ev_files,
    overwrite   = FORCE_OVERWRITE,
    provenance  = provenance
  )
  message("[q10]   registered ", length(ev_files), " evidence files")
}

message("[q10] DONE: cid=", cid)
