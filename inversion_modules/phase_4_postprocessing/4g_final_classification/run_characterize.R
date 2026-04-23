#!/usr/bin/env Rscript
# =============================================================================
# run_characterize.R — Driver for characterize_candidate() (Phase 4e)
# =============================================================================
#
# Chat 10 (2026-04-17) — FIX 45. Resolves Finding AF.
#
# This is the driver script for the Phase 4e characterization step.
# characterize_candidate.R is a *function library*; on its own it has no
# entry point. This driver reads per-candidate registry keys, calls
# characterize_candidate(keys) for each, and writes per-candidate +
# cohort-level outputs.
#
# CLI CONTRACT (mirrors compute_candidate_status.R):
#
#   Rscript run_characterize.R <registry_dir> <outdir>
#
# where <registry_dir> contains one of:
#     - per-candidate *.keys.tsv files, or
#     - merged evidence_registry.tsv, or
#     - candidate_scores.tsv.gz (fallback; produces Q1-only keys)
#
# OUTPUTS (written to <outdir>):
#
#   1. characterization.tsv
#        One row per candidate. Columns:
#          candidate_id, group_level, group_reason,
#          Q1_status, Q1_label, Q1_reason,
#          ... Q2 .. Q7 ...,
#          n_answered, n_partial, n_measured, n_contradicted, n_empty,
#          characterization_string
#
#   2. per_candidate/<cid>/characterization.txt
#        Human-readable report from format_characterization().
#
#   3. characterization_summary.txt
#        Genome-wide summary: status-cut counts, per-Q breakdowns,
#        evolutionary label distribution, top CONTRADICTED candidates.
#
#   4. candidate_full_report.tsv    (ONLY IF candidate_status.tsv exists)
#        Merge of candidate_status.tsv × characterization.tsv on
#        candidate_id. "One row, all three axes + per-Q convergence"
#        deliverable.
#
# INTEGRATION:
#   - `compute_candidate_status.R` computes the 3 axes (tier / completion /
#     evolutionary class). Run it first.
#   - `run_characterize.R` (this script) computes per-Q convergence.
#     Run it after, same registry_dir, same outdir. The combined report
#     merges on candidate_id.
#
# ERROR MODEL:
#   - Candidate with NO keys in registry  → row with all Qs = EMPTY
#     (status EMPTY, label NA, reason "no keys loaded"). Stderr warning.
#     Pipeline MUST tolerate this because aspirational keys are missing
#     for most candidates in the current pipeline state.
#   - characterize_candidate() throws an error → stderr warning, row
#     skipped, processing continues.
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Null-coalesce operator, matching characterize_candidate.R's semantics ──
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a[1]))) b else a
}

# =============================================================================
# LOCATE SIBLING/GATE SCRIPTS — explicit sourcing for SLURM robustness
# =============================================================================
# FIX 44 (chat 9) already added a path-search block inside
# characterize_candidate.R. We source it here FIRST and explicitly, so the
# gate functions are in scope regardless of cwd at SLURM run time.

get_script_dir <- function() {
  # Try several mechanisms in order; return "." if all fail.
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  # When sourced: sys.frame(1)$ofile may be set
  frames <- sys.frames()
  for (fr in rev(frames)) {
    ofile <- tryCatch(fr$ofile, error = function(e) NULL)
    if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile)))
  }
  "."
}

script_dir <- get_script_dir()

# Candidate paths for gate.R (most-specific first). GROUP_VALIDATION_GATE
# env var takes priority for explicit override.
gate_candidates <- c(
  Sys.getenv("GROUP_VALIDATION_GATE", ""),
  file.path(script_dir, "..", "4e_group_validation", "group_validation_gate.R"),
  "../4e_group_validation/group_validation_gate.R",
  "4e_group_validation/group_validation_gate.R",
  file.path(script_dir, "group_validation_gate.R")   # co-located fallback
)
gate_file <- NA_character_
for (p in gate_candidates) {
  if (nzchar(p) && file.exists(p)) { gate_file <- p; break }
}
if (is.na(gate_file)) {
  stop("[run_characterize] Could not locate group_validation_gate.R. ",
       "Set GROUP_VALIDATION_GATE env var or place it in the expected ",
       "location. Searched: ", paste(gate_candidates, collapse = "; "))
}
source(gate_file)

# characterize_candidate.R: same-directory sibling.
char_candidates <- c(
  file.path(script_dir, "characterize_candidate.R"),
  "characterize_candidate.R",
  "4g_final_classification/characterize_candidate.R"
)
char_file <- NA_character_
for (p in char_candidates) {
  if (nzchar(p) && file.exists(p)) { char_file <- p; break }
}
if (is.na(char_file)) {
  stop("[run_characterize] Could not locate characterize_candidate.R. ",
       "Expected as sibling of run_characterize.R.")
}
source(char_file)

# Sanity check: required functions are now in scope.
for (fn in c("characterize_candidate", "format_characterization",
             "assess_group_validation", "check_group_gate")) {
  if (!exists(fn, mode = "function")) {
    stop("[run_characterize] Required function '", fn, "' not in scope ",
         "after sourcing. Check gate.R and characterize_candidate.R.")
  }
}

# =============================================================================
# CLI
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript run_characterize.R <registry_dir> <outdir>\n")
  cat("  registry_dir: directory containing per-candidate *.keys.tsv files,\n")
  cat("                OR merged evidence_registry.tsv,\n")
  cat("                OR candidate_scores.tsv.gz (fallback)\n")
  cat("  outdir:       output directory (will be created if missing)\n")
  quit(save = "no", status = 1)
}
registry_dir <- args[1]
outdir       <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
per_cand_dir <- file.path(outdir, "per_candidate")
dir.create(per_cand_dir, showWarnings = FALSE, recursive = TRUE)

cat("[CHAR] registry_dir = ", registry_dir, "\n", sep = "")
cat("[CHAR] outdir       = ", outdir, "\n", sep = "")
cat("[CHAR] gate.R       = ", gate_file, "\n", sep = "")
cat("[CHAR] characterize = ", char_file, "\n", sep = "")

# =============================================================================
# LOADER — three-tier fallback, mirroring compute_candidate_status.R L727–827
# =============================================================================
# Tier 1: per-candidate *.keys.tsv (highest fidelity — all registered keys)
# Tier 2: merged evidence_registry.tsv (same keys, flat format)
# Tier 3: candidate_scores.tsv.gz (Q1-only reconstruction — most candidates
#         will have sparse Qs and trigger many EMPTY statuses — that's fine)

key_files <- list.files(registry_dir, pattern = "\\.keys\\.tsv$",
                        full.names = TRUE, recursive = TRUE)
merged_registry <- file.path(registry_dir, "evidence_registry.tsv")
score_file <- NULL

reg_dt <- NULL   # Tier 2 source (data.table with candidate_id, key, value)
sdt    <- NULL   # Tier 3 source (data.table with per-candidate row)

if (length(key_files) > 0) {
  cat("[CHAR] Tier 1: found ", length(key_files),
      " per-candidate key files\n", sep = "")
  candidates <- gsub("\\.keys\\.tsv$", "", basename(key_files))
} else if (file.exists(merged_registry)) {
  cat("[CHAR] Tier 2: reading merged registry ", merged_registry, "\n", sep = "")
  reg_dt <- fread(merged_registry)
  if (!all(c("candidate_id", "key", "value") %in% names(reg_dt))) {
    stop("[run_characterize] merged registry ", merged_registry,
         " missing required columns (candidate_id, key, value).")
  }
  candidates <- unique(reg_dt$candidate_id)
} else {
  # Tier 3 fallback.
  for (f in c(
    file.path(registry_dir, "candidate_scores.tsv.gz"),
    file.path(registry_dir, "..", "candidate_scores", "candidate_scores.tsv.gz"),
    Sys.glob(file.path(registry_dir, "..", "*", "candidate_scores.tsv.gz"))
  )) {
    if (length(f) == 1 && file.exists(f)) { score_file <- f; break }
  }
  if (is.null(score_file)) {
    stop("[run_characterize] No registry data found in ", registry_dir,
         ". Expected *.keys.tsv files, evidence_registry.tsv, or ",
         "candidate_scores.tsv.gz.")
  }
  cat("[CHAR] Tier 3: reconstructing from ", score_file, "\n", sep = "")
  sdt <- fread(score_file)
  if (!"candidate_id" %in% names(sdt)) {
    if (all(c("chrom", "interval_id") %in% names(sdt))) {
      sdt[, candidate_id := paste0(chrom, "_", interval_id)]
    } else {
      stop("[run_characterize] candidate_scores.tsv.gz missing candidate_id ",
           "and cannot construct it from chrom/interval_id.")
    }
  }
  candidates <- sdt$candidate_id
}

# Q1 reconstruction map for Tier 3 fallback (same mapping as
# compute_candidate_status.R L801–819).
score_col_to_key <- c(
  d1_block_strength       = "q1_d01_block_strength",
  d2_block_shape          = "q1_d02_block_shape",
  d3_nn_persistence       = "q1_d03_nn_persistence",
  d4_decay_flatness       = "q1_d04_decay_flatness",
  d5_interior_quality     = "q1_d05_interior_quality",
  d6_consensus            = "q1_d06_consensus",
  d7_sv_breakpoint        = "q1_d07_sv_breakpoint",
  d8_peel_or_hyp          = "q1_d08_peel_or_hyp",
  d9_pca_clusters         = "q1_d09_pca_clusters",
  d10_partition           = "q1_d10_partition",
  d11_boundary_concordance= "q1_d11_boundary_concordance",
  d12_snake_concordance   = "q1_d12_snake_concordance",
  final_score             = "q1_composite_score",
  dim_positive            = "q1_dim_positive",
  shape_class             = "q1_shape_class",
  span_mb                 = "q1_span_kb",
  tier                    = "q7_tier"
)

load_keys_for_candidate <- function(cid) {
  keys <- list()

  # Tier 1: per-candidate file
  kf <- file.path(registry_dir, paste0(cid, ".keys.tsv"))
  if (!file.exists(kf)) {
    # also try nested layout
    nested <- list.files(registry_dir,
                          pattern = paste0("^", cid, "\\.keys\\.tsv$"),
                          full.names = TRUE, recursive = TRUE)
    if (length(nested) > 0) kf <- nested[1]
  }
  if (file.exists(kf)) {
    kdt <- tryCatch(fread(kf), error = function(e) NULL)
    if (!is.null(kdt) && all(c("key", "value") %in% names(kdt))) {
      for (ki in seq_len(nrow(kdt))) {
        keys[[ kdt$key[ki] ]] <- kdt$value[ki]
      }
    }
  }

  # Tier 2: merged registry supplement (never overwrite Tier 1 — Tier 1
  # is the highest-fidelity source).
  if (!is.null(reg_dt)) {
    cand_rows <- reg_dt[candidate_id == cid]
    if (nrow(cand_rows) > 0) {
      for (ki in seq_len(nrow(cand_rows))) {
        k <- cand_rows$key[ki]
        if (is.null(keys[[k]])) keys[[k]] <- cand_rows$value[ki]
      }
    }
  }

  # Tier 3: scores data.table supplement
  if (!is.null(sdt) && cid %in% sdt$candidate_id) {
    row <- sdt[candidate_id == cid][1]
    for (col in names(score_col_to_key)) {
      if (col %in% names(row)) {
        val <- row[[col]]
        if (!is.null(val) && !is.na(val)) {
          target_key <- score_col_to_key[[col]]
          if (is.null(keys[[target_key]])) {
            if (col == "span_mb") val <- as.numeric(val) * 1000  # Mb → kb
            keys[[target_key]] <- as.character(val)
          }
        }
      }
    }
  }

  keys
}

# =============================================================================
# EMPTY-ROW BUILDER: for candidates with no registry data, so processing
# never crashes and the output TSV row count matches candidate count.
# =============================================================================
empty_row <- function(cid, reason = "no keys loaded for candidate") {
  empty_q <- function() list(status = "EMPTY", label = NA,
                              evidence_for = character(),
                              evidence_against = character(),
                              reason = reason)
  list(
    per_question = setNames(replicate(7, empty_q(), simplify = FALSE),
                             paste0("Q", 1:7)),
    group_validation = list(level = "NONE",
                             reason = "no keys loaded"),
    summary = list(
      answered = 0L, partial = 0L, measured = 0L,
      contradicted = 0L, empty = 7L,
      group_level = "NONE",
      characterization_string = paste0(
        "0/7 ANSWERED, 0 PARTIAL, 0 MEASURED, 0 CONTRADICTED, 7 EMPTY ",
        "[groups=NONE] (no keys)"
      )
    )
  )
}

# =============================================================================
# MAIN LOOP
# =============================================================================
cat("[CHAR] Processing ", length(candidates), " candidates...\n", sep = "")

char_rows     <- vector("list", length(candidates))
# Tracking for the summary report
n_empty_load  <- 0L
n_err         <- 0L
err_log       <- character()

for (ci in seq_along(candidates)) {
  cid <- candidates[ci]

  # ── Load keys ──
  keys <- tryCatch(load_keys_for_candidate(cid),
                   error = function(e) { err_log <<- c(err_log,
                                          sprintf("[LOAD ERROR] %s: %s", cid, e$message));
                                         list() })

  if (length(keys) == 0) {
    message("[CHAR] WARN: no keys loaded for ", cid,
            " — emitting all-EMPTY row")
    char_result <- empty_row(cid)
    n_empty_load <- n_empty_load + 1L
  } else {
    # ── Characterize ──
    char_result <- tryCatch(
      characterize_candidate(keys),
      error = function(e) {
        msg <- sprintf("[CHAR ERROR] %s: %s", cid, conditionMessage(e))
        message(msg)
        err_log <<- c(err_log, msg)
        n_err <<- n_err + 1L
        empty_row(cid, reason = paste0("characterize_candidate() error: ",
                                        conditionMessage(e)))
      }
    )
  }

  # ── Build TSV row ──
  pq <- char_result$per_question
  gv <- char_result$group_validation %||% list(level = "NONE", reason = "")
  sm <- char_result$summary %||% list()

  row <- list(candidate_id = cid,
              group_level  = gv$level %||% "NONE",
              group_reason = gv$reason %||% "")
  for (q in paste0("Q", 1:7)) {
    r <- pq[[q]] %||% list(status = "EMPTY", label = NA, reason = "")
    row[[paste0(q, "_status")]] <- r$status %||% "EMPTY"
    row[[paste0(q, "_label")]]  <- if (is.null(r$label) || length(r$label) == 0 ||
                                        is.na(r$label[1])) NA_character_
                                    else as.character(r$label)
    row[[paste0(q, "_reason")]] <- r$reason %||% ""
  }
  row$n_answered               <- sm$answered      %||% 0L
  row$n_partial                <- sm$partial       %||% 0L
  row$n_measured               <- sm$measured      %||% 0L
  row$n_contradicted           <- sm$contradicted  %||% 0L
  row$n_empty                  <- sm$empty         %||% 7L
  row$characterization_string  <- sm$characterization_string %||% ""

  char_rows[[ci]] <- as.data.table(row)

  # ── Per-candidate human-readable report ──
  cand_outdir <- file.path(per_cand_dir, cid)
  dir.create(cand_outdir, showWarnings = FALSE, recursive = TRUE)
  txt <- tryCatch(format_characterization(cid, char_result),
                   error = function(e) paste0(
                     "ERROR formatting ", cid, ": ", conditionMessage(e)))
  writeLines(txt, file.path(cand_outdir, "characterization.txt"))

  if (ci %% 20 == 0 || ci == length(candidates)) {
    cat("  [", ci, "/", length(candidates), "]\n", sep = "")
  }
}

# =============================================================================
# OUTPUT 1: characterization.tsv
# =============================================================================
char_dt <- rbindlist(char_rows, fill = TRUE)
fwrite(char_dt, file.path(outdir, "characterization.tsv"), sep = "\t")
cat("[CHAR] wrote ", nrow(char_dt), " rows -> ",
    file.path(outdir, "characterization.tsv"), "\n", sep = "")

# =============================================================================
# OUTPUT 2: per-candidate characterization.txt already written in main loop
# =============================================================================

# =============================================================================
# OUTPUT 3: characterization_summary.txt — genome-wide overview
# =============================================================================
status_cuts <- c("ANSWERED", "PARTIAL", "MEASURED", "CONTRADICTED", "EMPTY")
group_cuts  <- c("NONE", "SUSPECT", "UNCERTAIN", "SUPPORTED", "VALIDATED")

count_status <- function(q) {
  col <- paste0(q, "_status")
  if (!col %in% names(char_dt)) return(setNames(integer(length(status_cuts)),
                                                 status_cuts))
  tab <- table(factor(char_dt[[col]], levels = status_cuts))
  as.integer(tab)
}

# Per-Q status matrix (7 questions × 5 status cuts)
per_q_matrix <- do.call(rbind, lapply(paste0("Q", 1:7), count_status))
rownames(per_q_matrix) <- paste0("Q", 1:7)
colnames(per_q_matrix) <- status_cuts

# Group-level distribution
group_tab <- table(factor(char_dt$group_level, levels = group_cuts))

# Top CONTRADICTED candidates (rows with any Q status = CONTRADICTED)
contradicted_mask <- apply(
  as.matrix(char_dt[, paste0("Q", 1:7, "_status"), with = FALSE]),
  1, function(r) any(r == "CONTRADICTED", na.rm = TRUE)
)
contradicted_cids <- char_dt$candidate_id[contradicted_mask]
n_contradicted_total <- length(contradicted_cids)

# Label distribution across all candidates × Qs (compact)
all_labels <- unlist(lapply(paste0("Q", 1:7), function(q) {
  col <- paste0(q, "_label")
  if (col %in% names(char_dt)) paste0(q, ":", char_dt[[col]]) else character(0)
}))
all_labels <- all_labels[!grepl(":NA$", all_labels) & !is.na(all_labels)]
label_tab <- sort(table(all_labels), decreasing = TRUE)

summary_lines <- c(
  paste(rep("=", 80), collapse = ""),
  "PHASE 4e CHARACTERIZATION SUMMARY",
  paste(rep("=", 80), collapse = ""),
  sprintf("  Total candidates processed:    %d", nrow(char_dt)),
  sprintf("  Candidates with no keys (empty): %d", n_empty_load),
  sprintf("  Candidates with errors:         %d", n_err),
  "",
  paste(rep("-", 80), collapse = ""),
  "PER-QUESTION STATUS DISTRIBUTION (counts across candidates)",
  paste(rep("-", 80), collapse = ""),
  sprintf("  %-4s %12s %8s %8s %12s %8s",
          "Q", status_cuts[1], status_cuts[2], status_cuts[3],
          status_cuts[4], status_cuts[5])
)
for (q in paste0("Q", 1:7)) {
  v <- per_q_matrix[q, ]
  summary_lines <- c(summary_lines,
    sprintf("  %-4s %12d %8d %8d %12d %8d",
            q, v["ANSWERED"], v["PARTIAL"], v["MEASURED"],
            v["CONTRADICTED"], v["EMPTY"]))
}

summary_lines <- c(summary_lines,
  "",
  paste(rep("-", 80), collapse = ""),
  "GROUP VALIDATION DISTRIBUTION (across candidates)",
  paste(rep("-", 80), collapse = "")
)
for (g in group_cuts) {
  summary_lines <- c(summary_lines,
    sprintf("  %-12s %6d", g, as.integer(group_tab[g] %||% 0)))
}

summary_lines <- c(summary_lines,
  "",
  paste(rep("-", 80), collapse = ""),
  sprintf("CANDIDATES WITH AT LEAST ONE CONTRADICTION: %d", n_contradicted_total),
  paste(rep("-", 80), collapse = "")
)
if (n_contradicted_total > 0) {
  head_n <- min(n_contradicted_total, 25L)
  for (i in seq_len(head_n)) {
    cid <- contradicted_cids[i]
    cs  <- char_dt[candidate_id == cid, characterization_string][1]
    summary_lines <- c(summary_lines,
      sprintf("  %-30s  %s", cid, cs))
  }
  if (n_contradicted_total > head_n) {
    summary_lines <- c(summary_lines,
      sprintf("  ... (+%d more CONTRADICTED candidates — see characterization.tsv)",
              n_contradicted_total - head_n))
  }
}

summary_lines <- c(summary_lines,
  "",
  paste(rep("-", 80), collapse = ""),
  "TOP EVOLUTIONARY LABELS ASSIGNED (across all Q × candidates)",
  paste(rep("-", 80), collapse = "")
)
top_n <- min(length(label_tab), 15L)
if (top_n > 0) {
  for (i in seq_len(top_n)) {
    summary_lines <- c(summary_lines,
      sprintf("  %-45s %6d", names(label_tab)[i], as.integer(label_tab[i])))
  }
} else {
  summary_lines <- c(summary_lines, "  (no labels assigned yet)")
}

if (length(err_log) > 0) {
  summary_lines <- c(summary_lines, "",
    paste(rep("-", 80), collapse = ""),
    sprintf("ERROR LOG (%d entries)", length(err_log)),
    paste(rep("-", 80), collapse = ""),
    head(err_log, 50L))
  if (length(err_log) > 50L) {
    summary_lines <- c(summary_lines,
      sprintf("  ... (+%d more errors)", length(err_log) - 50L))
  }
}

summary_lines <- c(summary_lines, "", paste(rep("=", 80), collapse = ""))
writeLines(summary_lines, file.path(outdir, "characterization_summary.txt"))
cat("[CHAR] wrote summary -> ",
    file.path(outdir, "characterization_summary.txt"), "\n", sep = "")

# =============================================================================
# OUTPUT 4: candidate_full_report.tsv (IF candidate_status.tsv is present)
# =============================================================================
status_file <- file.path(outdir, "candidate_status.tsv")
if (file.exists(status_file)) {
  status_dt <- tryCatch(fread(status_file), error = function(e) NULL)
  if (!is.null(status_dt) && "candidate_id" %in% names(status_dt)) {
    full_dt <- merge(status_dt, char_dt, by = "candidate_id",
                      all = TRUE, suffixes = c("", ".char"))
    fwrite(full_dt, file.path(outdir, "candidate_full_report.tsv"), sep = "\t")
    cat("[CHAR] wrote ", nrow(full_dt), " rows -> ",
        file.path(outdir, "candidate_full_report.tsv"),
        " (merged with candidate_status.tsv)\n", sep = "")
  } else {
    cat("[CHAR] candidate_status.tsv present but unreadable / missing ",
        "candidate_id column — skipping merge\n", sep = "")
  }
} else {
  cat("[CHAR] candidate_status.tsv not found in outdir — skipping ",
      "candidate_full_report.tsv merge. Run compute_candidate_status.R ",
      "first if you want the combined report.\n", sep = "")
}

cat("\n[DONE] Phase 4e characterization complete -> ", outdir, "\n", sep = "")
