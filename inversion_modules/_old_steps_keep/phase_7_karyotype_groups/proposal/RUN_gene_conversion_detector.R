#!/usr/bin/env Rscript
# =============================================================================
# RUN_gene_conversion_detector.R — per-candidate driver for
#   gene_conversion_detector.R
# =============================================================================
# Wires the orphaned gene-conversion detector into the batch pipeline.
# Before this driver, gene_conversion_detector.R was built and tested
# (phase_9_classification/tests/test_gc_detector.R) but never called
# from production. C01i_b_multi_recomb.R had a load_gc_for_cid() helper
# and a --gc_dir CLI option expecting this driver's output; the driver
# itself was missing.
#
# =============================================================================
# EXECUTION ORDER
# =============================================================================
# This script runs AFTER STEP_C01i_decompose.R (needs its
# per_window_class.rds for baseline classes per sample) and BEFORE
# STEP_C01i_b_multi_recomb.R (which reads --gc_dir pointing at this
# driver's output). Typical chain:
#
#   C01d → C01i_decompose → RUN_gene_conversion_detector → C01i_b
#        → C01i_c_nested_composition → C01i_d_seal → C01f
#
# It does NOT require C01i_d_seal to have registered groups — it pulls
# the baseline class directly from decompose's RDS. This keeps the
# driver independent of the seal step and safe to re-run on subsets.
#
# =============================================================================
# INPUTS
# =============================================================================
#   --candidates   TSV with candidate_id, chrom, start_bp, end_bp, tier
#   --decomp_dir   STEP_C01i_decompose output root; per_window_class.rds
#                  lives at <decomp_dir>/<cid>/per_window_class.rds
#   --outdir       Where to write per-candidate gene_conversion_tracts.json.
#                  Pass this path as --gc_dir to STEP_C01i_b_multi_recomb.
#   --dosage_dir   Per-chromosome dosage files. Default: env DOSAGE_DIR
#                  or ${BASE}/popstruct_thin/04_beagle_byRF_majmin.
#                  Expected layout: <dosage_dir>/<chr>.dosage.tsv.gz
#                  (markers × samples) + <chr>.sites.tsv.gz (pos column).
#                  Same convention as phase_6/01_dosage_signal.R.
#   --tier_max     Only process candidates with tier <= this (default 3).
#   --flank_kb     Extra flank kb to include when slicing SNPs. Default 0:
#                  GC tracts are by definition INSIDE the inversion
#                  interval, so we don't need flank. Tunable for sensitivity
#                  studies.
#
# Detector-parameter passthroughs (all optional; defaults match
# detect_cohort_gene_conversion_v2 in gene_conversion_detector.R):
#   --min_run, --max_tolerance, --max_span_bp, --max_flagged_snps,
#   --min_delta_af, --max_het_fraction, --min_maf, --min_call_rate,
#   --depth_z_abs, --min_per_class
#
# =============================================================================
# OUTPUTS
# =============================================================================
#   <outdir>/<cid>/gene_conversion_tracts.json    full detector block
#   <outdir>/<cid>/gene_conversion_tracts.rds     same, as RDS (faster load)
#   <outdir>/_summary.tsv                         per-candidate driver summary
#
# Plus, if the registry is available:
#   reg$evidence$write_block(cid, "gene_conversion_tracts", block)
# which emits the 3 flat keys q2_gc_total_tracts, q2_gc_total_samples_with_gc,
# q2_gc_n_snps_pass_qc_diagnostic per the schema's keys_extracted rule.
#
# =============================================================================
# REGISTRY_CONTRACT
#   BLOCKS_WRITTEN:
#     - gene_conversion_tracts: registries/schemas/structured_block_schemas/gene_conversion_tracts.schema.json
#       keys: q2_gc_total_tracts, q2_gc_total_samples_with_gc,
#             q2_gc_n_snps_pass_qc_diagnostic,
#             q2_gc_n_samples_scannable,
#             q2_gc_mean_n_diagnostic_scanned_per_sample
#       status: WIRED
#       note: Per-candidate write_block call wrapped in tryCatch — the JSON
#             output is the primary handoff (C01i_b reads it via
#             load_gc_for_cid), the registry write is additive so the 5 flat
#             schema keys surface in keys.tsv. A registry-unavailable run
#             still produces valid JSON and is safe to re-join later.
#             Driver written 2026-04-24 (chat D) to unblock the previously
#             orphaned detector. Keys q2_gc_n_samples_scannable and
#             q2_gc_mean_n_diagnostic_scanned_per_sample added same day
#             with the diagnostic_snps block-field addition — audit the
#             per-sample scan coverage, useful for distinguishing
#             "no tracts because no opportunity" from "no tracts because
#             nothing was there".
#   KEYS_IN: none
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ── Parse args ───────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--candidates", type = "character", help = "Candidate table TSV"),
  make_option("--decomp_dir", type = "character", default = "decomp_out",
              help = "STEP_C01i_decompose output root"),
  make_option("--outdir",     type = "character", default = "gc_out",
              help = "Per-candidate GC detector output root"),
  make_option("--dosage_dir", type = "character", default = NA_character_,
              help = paste0("Per-chromosome dosage dir. If unset, falls back ",
                            "to env DOSAGE_DIR, then ${BASE}/popstruct_thin/",
                            "04_beagle_byRF_majmin (matching 01_dosage_signal.R).")),
  make_option("--tier_max",   type = "integer", default = 3L),
  make_option("--flank_kb",   type = "integer", default = 0L,
              help = "Flank kb around candidate interval when slicing SNPs."),
  # Detector-parameter passthroughs
  make_option("--min_run",          type = "integer", default = 2L),
  make_option("--max_tolerance",    type = "integer", default = 1L),
  make_option("--max_span_bp",      type = "integer", default = 20000L),
  make_option("--max_flagged_snps", type = "integer", default = 10L),
  make_option("--min_delta_af",     type = "double",  default = 0.5),
  make_option("--max_het_fraction", type = "double",  default = 0.70),
  make_option("--min_maf",          type = "double",  default = 0.02),
  make_option("--min_call_rate",    type = "double",  default = 0.80),
  make_option("--depth_z_abs",      type = "double",  default = 3.0),
  make_option("--min_per_class",    type = "integer", default = 3L)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$candidates) || !file.exists(opt$candidates)) {
  stop("--candidates must be an existing file")
}

# ── Resolve DOSAGE_DIR ──────────────────────────────────────────────────────
# Same resolution order as 01_dosage_signal.R: --dosage_dir > env > BASE path.
resolve_dosage_dir <- function(arg_val) {
  if (!is.na(arg_val) && nzchar(arg_val)) return(arg_val)
  env_val <- Sys.getenv("DOSAGE_DIR", "")
  if (nzchar(env_val)) return(env_val)
  base <- Sys.getenv("BASE",
    "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04")
  file.path(base, "popstruct_thin", "04_beagle_byRF_majmin")
}
DOSAGE_DIR <- resolve_dosage_dir(opt$dosage_dir)
cat("[gc_driver] DOSAGE_DIR = ", DOSAGE_DIR, "\n", sep = "")

# ── Source the detector ─────────────────────────────────────────────────────
# Try local, script-relative, and repo-root locations so the driver works
# whether launched from proposal/, repo root, or a SLURM job.
detector_paths <- c(
  "gene_conversion_detector.R",
  file.path(dirname(sys.frame(1)$ofile %||% "."),
            "gene_conversion_detector.R"),
  file.path(Sys.getenv("BASE", "."),
            "inversion-popgen-toolkit",
            "inversion_modules", "phase_7_karyotype_groups",
            "proposal", "gene_conversion_detector.R")
)
det_found <- FALSE
for (p in detector_paths) {
  if (file.exists(p)) {
    source(p); det_found <- TRUE; break
  }
}
if (!det_found) {
  stop("[gc_driver] cannot locate gene_conversion_detector.R — tried: ",
       paste(detector_paths, collapse = ", "))
}
if (!exists("detect_cohort_gene_conversion_v2", mode = "function")) {
  stop("[gc_driver] detect_cohort_gene_conversion_v2() not defined after ",
       "sourcing the detector — maybe a stale copy of the file?")
}

# ── Optional registry ────────────────────────────────────────────────────────
# Follow the bridge-source pattern used by C01i_decompose. If the registry
# isn't wired, the driver still writes JSON/RDS and runs fine.
reg <- NULL
bridge_paths <- c(
  "utils/registry_bridge.R",
  "../../../utils/registry_bridge.R",
  file.path(Sys.getenv("BASE", "."),
            "inversion-popgen-toolkit", "utils", "registry_bridge.R")
)
for (p in bridge_paths) {
  if (file.exists(p)) {
    tryCatch({ source(p) }, error = function(e)
              message("[gc_driver] registry_bridge source failed: ",
                      conditionMessage(e)))
    break
  }
}
if (exists("load_registry", mode = "function")) {
  reg <- tryCatch(load_registry(),
                  error = function(e) {
                    message("[gc_driver] load_registry() failed: ",
                            conditionMessage(e), " — JSON-only mode.")
                    NULL
                  })
}

# ── Candidate table ─────────────────────────────────────────────────────────
cand_dt <- fread(opt$candidates)
if (!"candidate_id" %in% names(cand_dt)) {
  if ("interval_id" %in% names(cand_dt) && "chrom" %in% names(cand_dt)) {
    cand_dt[, candidate_id := paste0(chrom, "_", interval_id)]
  } else {
    stop("[gc_driver] candidates table needs candidate_id or chrom+interval_id")
  }
}
if (!"tier" %in% names(cand_dt))     cand_dt[, tier := 2L]
if (!"start_bp" %in% names(cand_dt)) cand_dt[, start_bp := as.integer(start_mb * 1e6)]
if (!"end_bp"   %in% names(cand_dt)) cand_dt[, end_bp   := as.integer(end_mb * 1e6)]
cand_dt <- cand_dt[tier <= opt$tier_max]

cat("[gc_driver] processing ", nrow(cand_dt), " candidates at tier<=",
    opt$tier_max, "\n", sep = "")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Per-chrom dosage/sites cache — most runs have multiple candidates per chrom
# and re-reading the same 400+ MB TSV is the dominant cost. Cache by chrom.
.dos_cache <- new.env(parent = emptyenv())

load_chr_dosage <- function(chr) {
  key <- chr
  if (exists(key, envir = .dos_cache, inherits = FALSE)) {
    return(get(key, envir = .dos_cache, inherits = FALSE))
  }
  dos_file   <- file.path(DOSAGE_DIR, paste0(chr, ".dosage.tsv.gz"))
  sites_file <- file.path(DOSAGE_DIR, paste0(chr, ".sites.tsv.gz"))
  if (!file.exists(dos_file) || !file.exists(sites_file)) {
    message("[gc_driver]   dosage files missing for ", chr, " under ",
            DOSAGE_DIR)
    return(NULL)
  }
  dos   <- fread(dos_file)
  sites <- fread(sites_file)
  obj <- list(dos = dos, sites = sites)
  assign(key, obj, envir = .dos_cache)
  obj
}

# ── Per-candidate loop ───────────────────────────────────────────────────────
summary_rows <- list()

for (ci in seq_len(nrow(cand_dt))) {
  cd  <- cand_dt[ci]
  cid <- cd$candidate_id
  chr <- cd$chrom
  s   <- as.integer(cd$start_bp)
  e   <- as.integer(cd$end_bp)
  cat("[gc_driver] (", ci, "/", nrow(cand_dt), ") ", cid, " ",
      chr, ":", s, "-", e, "\n", sep = "")

  # 1. Load baseline class from decompose's per_window_class.rds
  pw_file <- file.path(opt$decomp_dir, cid, "per_window_class.rds")
  if (!file.exists(pw_file)) {
    message("[gc_driver]   SKIP — no per_window_class.rds at ", pw_file,
            " (did STEP_C01i_decompose run for this cid?)")
    next
  }
  pw <- tryCatch(readRDS(pw_file), error = function(err) NULL)
  if (is.null(pw) || is.null(pw$class_assignment_global)) {
    message("[gc_driver]   SKIP — decompose RDS lacks class_assignment_global")
    next
  }
  baseline_class <- pw$class_assignment_global  # named character vector

  # 2. Load per-chrom dosage + sites
  chrdat <- load_chr_dosage(chr)
  if (is.null(chrdat)) {
    message("[gc_driver]   SKIP — no dosage for ", chr)
    next
  }

  # 3. Slice to candidate interval (+ optional flank)
  flank_bp <- opt$flank_kb * 1000L
  scan_start <- max(1L, s - flank_bp)
  scan_end   <- e + flank_bp
  keep <- which(chrdat$sites$pos >= scan_start &
                chrdat$sites$pos <= scan_end)
  if (length(keep) < 20L) {
    message("[gc_driver]   SKIP — only ", length(keep),
            " SNPs in interval ", chr, ":", scan_start, "-", scan_end)
    next
  }
  dos_sl   <- chrdat$dos[keep]
  sites_sl <- chrdat$sites[keep]

  # 4. Build samples × SNPs dosage matrix. dos has one row per marker and
  #    one column per sample; the detector wants samples × SNPs.
  sample_cols <- setdiff(names(dos_sl), "marker")
  if (length(sample_cols) == 0L) {
    message("[gc_driver]   SKIP — no sample columns in dosage file")
    next
  }
  # Intersect with samples that HAVE a baseline class (must be non-NA,
  # HOM_REF/HET/HOM_INV — AMBIGUOUS tolerated; NA dropped).
  class_samples <- intersect(sample_cols, names(baseline_class))
  class_samples <- class_samples[
    !is.na(baseline_class[class_samples]) &
    baseline_class[class_samples] %in% c("HOM_REF", "HET", "HOM_INV",
                                           "AMBIGUOUS")
  ]
  if (length(class_samples) < 10L) {
    message("[gc_driver]   SKIP — only ", length(class_samples),
            " samples have both dosage and a baseline class (<10 threshold)")
    next
  }
  # Build the matrix: transpose and keep only class_samples.
  mat <- t(as.matrix(dos_sl[, ..class_samples]))
  storage.mode(mat) <- "double"
  rownames(mat) <- class_samples
  snp_pos <- as.integer(sites_sl$pos)

  # 5. Run the detector
  block <- tryCatch(
    detect_cohort_gene_conversion_v2(
      dosage_mat       = mat,
      snp_pos          = snp_pos,
      baseline_class_by_sample = baseline_class[class_samples],
      candidate_id     = cid,
      min_run          = opt$min_run,
      max_tolerance    = opt$max_tolerance,
      max_span_bp      = opt$max_span_bp,
      max_flagged_snps = opt$max_flagged_snps,
      min_delta_af     = opt$min_delta_af,
      max_het_fraction = opt$max_het_fraction,
      min_maf          = opt$min_maf,
      min_call_rate    = opt$min_call_rate,
      depth_z_abs      = opt$depth_z_abs,
      min_per_class    = opt$min_per_class
    ),
    error = function(err) {
      message("[gc_driver]   detector error for ", cid, ": ",
              conditionMessage(err))
      NULL
    }
  )
  if (is.null(block)) next

  # 6. Persist
  cand_out <- file.path(opt$outdir, cid)
  dir.create(cand_out, recursive = TRUE, showWarnings = FALSE)
  # JSON for C01i_b's load_gc_for_cid
  json_path <- file.path(cand_out, "gene_conversion_tracts.json")
  tryCatch(
    writeLines(toJSON(block, auto_unbox = TRUE, pretty = TRUE,
                       na = "null", dataframe = "columns"),
                json_path),
    error = function(err) message("[gc_driver]   JSON write failed: ",
                                   conditionMessage(err))
  )
  # RDS for quick re-reads
  rds_path  <- file.path(cand_out, "gene_conversion_tracts.rds")
  tryCatch(saveRDS(block, rds_path),
            error = function(err) message("[gc_driver]   RDS write failed: ",
                                           conditionMessage(err)))

  # 7. Registry write (additive; schema keys_extracted surface 3 q2_gc_* keys)
  if (!is.null(reg) && !is.null(reg$evidence) &&
      !is.null(reg$evidence$write_block)) {
    Sys.setenv(CURRENT_SCRIPT = "RUN_gene_conversion_detector.R")
    tryCatch(
      reg$evidence$write_block(candidate_id = cid,
                                block_type   = "gene_conversion_tracts",
                                data         = block),
      error = function(err) message("[gc_driver]   registry write failed: ",
                                     conditionMessage(err))
    )
  }

  # 8. Row for the summary TSV
  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    candidate_id              = cid,
    chrom                     = chr,
    start_bp                  = s,
    end_bp                    = e,
    n_samples_scanned         = length(class_samples),
    n_snps_in_interval        = length(keep),
    n_snps_pass_qc            = block$snp_qc$n_snps_pass_qc %||% NA_integer_,
    n_snps_pass_qc_diagnostic = block$snp_qc$n_snps_pass_qc_diagnostic %||% NA_integer_,
    total_tracts              = block$total_tracts %||% 0L,
    total_samples_with_gc     = block$total_samples_with_gc %||% 0L
  )
}

# ── Driver summary ──────────────────────────────────────────────────────────
if (length(summary_rows) > 0L) {
  summ <- rbindlist(summary_rows, fill = TRUE)
  summ_path <- file.path(opt$outdir, "_summary.tsv")
  fwrite(summ, summ_path, sep = "\t")
  cat("\n[gc_driver] summary written to ", summ_path, "\n", sep = "")
  cat("[gc_driver] tract count distribution:\n")
  print(summary(summ$total_tracts))
} else {
  cat("\n[gc_driver] no candidates produced GC output.\n")
}

cat("\n[gc_driver] done.\n")
