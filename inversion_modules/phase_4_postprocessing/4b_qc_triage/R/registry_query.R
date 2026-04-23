#!/usr/bin/env Rscript
# =============================================================================
# registry_query.R — command-line interrogation of the four registries
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)[1]
  if (is.na(i) || i >= length(args)) return(default)
  args[i + 1]
}

REGISTRIES_ROOT <- get_arg("--registries_root")
REGISTRY_API_R  <- get_arg("--registry_api_r")
CMD             <- get_arg("--cmd", "summary")

# Positional args after the named ones
pos <- args[!grepl("^--", args) & c(FALSE, grepl("^--", args[-length(args)]))]
pos <- setdiff(pos, args[seq(1, length(args), 2)[seq(1, length(args), 2) <= length(args)]])
# Simpler: extract positional args = any arg not immediately after a flag
flag_ix <- which(grepl("^--", args))
flag_val_ix <- flag_ix + 1
consumed <- sort(unique(c(flag_ix, flag_val_ix[flag_val_ix <= length(args)])))
pos_ix <- setdiff(seq_along(args), consumed)
pos <- args[pos_ix]

data_dir <- function(registry) file.path(REGISTRIES_ROOT, "data", registry)

for (nm in c("interval_registry", "sample_registry", "results_registry", "evidence_registry")) {
  f <- file.path(REGISTRY_API_R, paste0(nm, ".R"))
  if (file.exists(f)) source(f)
}

safe_load <- function(loader, reg_name) {
  tryCatch(loader(data_dir(reg_name)), error = function(e) {
    cat(sprintf("[query] %s: not available (%s)\n", reg_name, conditionMessage(e)))
    NULL
  })
}

intervals <- safe_load(load_interval_registry, "interval_registry")
samples   <- safe_load(load_sample_registry,   "sample_registry")
results   <- safe_load(load_results_registry,  "results_registry")
evidence  <- safe_load(load_evidence_registry, "evidence_registry")

# ---- Commands ---------------------------------------------------------------

cmd_summary <- function() {
  cat("\n=== REGISTRY SUMMARY ===\n\n")
  if (!is.null(intervals)) {
    ivs <- intervals$list()
    cat(sprintf("Intervals: %d candidates across %d chromosomes\n",
                nrow(ivs), length(unique(ivs$chrom))))
    if (nrow(ivs) > 0) {
      by_chr <- ivs[, .N, by = chrom][order(-N)]
      cat("  by chromosome:\n")
      for (i in seq_len(min(10, nrow(by_chr)))) {
        cat(sprintf("    %s: %d\n", by_chr$chrom[i], by_chr$N[i]))
      }
      if (nrow(by_chr) > 10) cat(sprintf("    (+ %d more)\n", nrow(by_chr) - 10))
    }
  }
  if (!is.null(samples)) {
    g <- samples$list_groups()
    cat(sprintf("\nSample groups: %d groups, %d samples in master\n",
                nrow(g), nrow(samples$get_master())))
    dims <- g[, .N, by = dimension][order(-N)]
    for (i in seq_len(nrow(dims))) {
      cat(sprintf("  %s: %d groups\n", dims$dimension[i], dims$N[i]))
    }
  }
  if (!is.null(results)) {
    m <- results$list_manifest()
    cat(sprintf("\nResults manifest: %d rows\n", nrow(m)))
    if (nrow(m) > 0) {
      by_kind <- m[, .N, by = .(kind, stat)][order(-N)]
      for (i in seq_len(min(10, nrow(by_kind)))) {
        cat(sprintf("  kind=%s stat=%s : %d\n",
                    by_kind$kind[i], by_kind$stat[i], by_kind$N[i]))
      }
    }
  }
  if (!is.null(evidence)) {
    ec <- evidence$list_candidates()
    cat(sprintf("\nEvidence: %d candidate bundles\n", nrow(ec)))
  }
  cat("\n")
}

cmd_list_candidates <- function() {
  if (is.null(intervals)) { cat("interval_registry unavailable\n"); return() }
  filter_chrom <- get_arg("--chrom", NULL)
  filter_method <- get_arg("--method", NULL)
  ivs <- intervals$list()
  if (!is.null(filter_chrom))  ivs <- ivs[chrom == filter_chrom]
  if (!is.null(filter_method)) ivs <- ivs[method == filter_method]
  if (nrow(ivs) == 0) { cat("no candidates match\n"); return() }
  cat(sprintf("\n=== %d CANDIDATES ===\n\n", nrow(ivs)))
  print(ivs[, .(cid, chrom, start_bp, end_bp,
                span_mb = round((end_bp - start_bp) / 1e6, 3),
                method, confidence)])
}

cmd_describe <- function() {
  if (length(pos) == 0L) { cat("usage: describe <cid>\n"); return() }
  cid <- pos[1]
  cat(sprintf("\n=== CANDIDATE: %s ===\n\n", cid))
  if (!is.null(intervals) && intervals$has(cid)) {
    iv <- intervals$get(cid)
    cat("INTERVAL:\n")
    for (k in names(iv)) cat(sprintf("  %-20s %s\n", k, paste(iv[[k]], collapse=",")))
  }
  if (!is.null(samples)) {
    subs <- samples$get_subgroups(cid)
    cat(sprintf("\nKARYOTYPE GROUPS (%d):\n", nrow(subs)))
    if (nrow(subs) > 0) {
      for (i in seq_len(nrow(subs))) {
        cat(sprintf("  %-35s karyo=%-5s dim=%-10s n=%d\n",
                    subs$group_id[i], subs$karyotype[i],
                    subs$dimension[i], subs$n[i]))
      }
    }
  }
  if (!is.null(results)) {
    rows <- results$ask(where = cid)
    cat(sprintf("\nRESULTS (%d rows):\n", nrow(rows)))
    if (nrow(rows) > 0) {
      for (i in seq_len(nrow(rows))) {
        cat(sprintf("  kind=%-10s stat=%-15s who=%-20s method=%s\n",
                    rows$kind[i], rows$stat[i], rows$who[i], rows$method[i]))
      }
    }
  }
  if (!is.null(evidence)) {
    ev <- evidence$get_candidate(cid)
    if (!is.null(ev)) {
      cat(sprintf("\nEVIDENCE FILES (%d):\n", length(ev$files)))
      for (f in ev$files) cat(sprintf("  %s\n", f))
    }
  }
  cat("\n")
}

cmd_karyotypes <- function() {
  if (length(pos) == 0L) { cat("usage: karyotypes <cid>\n"); return() }
  if (is.null(samples)) { cat("sample_registry unavailable\n"); return() }
  cid <- pos[1]
  subs <- samples$get_subgroups(cid)
  if (nrow(subs) == 0) { cat("no karyotype groups for ", cid, "\n"); return() }
  for (i in seq_len(nrow(subs))) {
    ids <- samples$get_group(subs$group_id[i])
    cat(sprintf("\n== %s (karyo=%s, n=%d) ==\n", subs$group_id[i],
                subs$karyotype[i], length(ids)))
    cat(paste(ids, collapse = "\n"), "\n")
  }
}

cmd_fst <- function() {
  if (length(pos) == 0L) { cat("usage: fst <cid>\n"); return() }
  if (is.null(results)) { cat("results_registry unavailable\n"); return() }
  cid <- pos[1]
  rows <- results$ask(where = cid, kind = "pairwise", stat = "fst")
  if (nrow(rows) == 0) { cat("no Fst rows for ", cid, "\n"); return() }
  for (i in seq_len(nrow(rows))) {
    f <- rows$file[i]
    if (!file.exists(f)) next
    dt <- fread(f)
    cat(sprintf("\n== %s (who=%s, %d windows) ==\n",
                rows$kind[i], rows$who[i], nrow(dt)))
    # Print summary
    if ("fst" %in% names(dt)) {
      cat(sprintf("  Fst mean=%.3f  median=%.3f  max=%.3f  n=%d\n",
                  mean(dt$fst, na.rm=TRUE), median(dt$fst, na.rm=TRUE),
                  max(dt$fst, na.rm=TRUE), nrow(dt)))
    }
  }
  # Also show summary stats
  sum_rows <- results$ask(where = cid, stat = "fst_summary")
  if (nrow(sum_rows) > 0) {
    cat("\n== SUMMARY STATS ==\n")
    for (i in seq_len(nrow(sum_rows))) {
      f <- sum_rows$file[i]
      if (file.exists(f)) {
        v <- fromJSON(f)
        for (k in names(v)) cat(sprintf("  %-35s %s\n", k, paste(v[[k]], collapse=",")))
      }
    }
  }
}

# Dispatch
switch(CMD,
  summary         = cmd_summary(),
  list_candidates = cmd_list_candidates(),
  describe        = cmd_describe(),
  karyotypes      = cmd_karyotypes(),
  fst             = cmd_fst(),
  {
    cat("Unknown command: ", CMD, "\n")
    cat("Available: summary, list_candidates, describe, karyotypes, fst\n")
  }
)
