#!/usr/bin/env Rscript

# =============================================================================
# export_relatedness_to_json_v1.R
#
# Reads ngsRelate pairwise output (or its NAToRA-format 3-column derivative)
# and emits a per-cohort relatedness JSON for the scrubber. The JSON carries:
#
#   - hub_id_1st : integer per-sample, connected-component id at theta >= 0.177
#                  (1st-degree threshold: parent-offspring, full sibs)
#   - hub_id_2nd : connected-component id at theta >= 0.0884
#                  (2nd-degree: half sibs, grandparent, aunt-uncle)
#   - hub_id_3rd : connected-component id at theta >= 0.0442
#                  (3rd-degree: first cousins; in hatchery cohorts this may
#                   reflect background kinship, not literal first-cousinship)
#   - pairs      : raw {a, b, theta} array for interactive threshold override
#                  in the scrubber (capped at top-N pairs by theta to keep JSON
#                  size reasonable)
#
# This is a COHORT-LEVEL artifact, not per-chromosome — relatedness doesn't
# vary along the genome (it's a per-sample-pair scalar). The scrubber loads
# it once and uses it for every chromosome.
#
# Input formats (auto-detected):
#
#   A) NAToRA 3-column format (catfish_226_for_natora.txt):
#      id1<TAB>id2<TAB>theta    no header, sample names directly
#
#   B) Full ngsRelate .res format (catfish_226_relatedness.res):
#      header line + many columns. ida/idb in cols 3-4, theta in col 18.
#      (when -z flag was passed to ngsRelate, ida/idb are sample names;
#      otherwise they're 0-based numeric indices and we need a sample list).
#
# Usage:
#   Rscript export_relatedness_to_json_v1.R \
#     --pairs    catfish_226_for_natora.txt \
#     --samples  /path/to/sample_list.txt \
#     --out      relatedness.json
#
#   Optional:
#     --max_pairs 50000   (cap raw pairs in JSON; default 50k)
#     --t1 0.177          (override 1st-degree threshold)
#     --t2 0.0884         (override 2nd-degree threshold)
#     --t3 0.0442         (override 3rd-degree threshold)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# =============================================================================
# Parse args
# =============================================================================
usage <- function() {
  cat(paste(
    "Usage: Rscript export_relatedness_to_json_v1.R --pairs <pairs.txt> --samples <samples.txt> --out <relatedness.json>",
    "  --pairs      ngsRelate pairwise output (.res or 3-col NAToRA format)",
    "  --samples    sample list (one ID per line, in cohort order)",
    "  --out        output JSON path",
    "  [--max_pairs N]  cap raw pairs in JSON (default 50000, top-N by theta)",
    "  [--t1 0.177]     1st-degree threshold (default 0.177)",
    "  [--t2 0.0884]    2nd-degree threshold (default 0.0884)",
    "  [--t3 0.0442]    3rd-degree threshold (default 0.0442)",
    sep = "\n"), "\n")
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)
PAIRS    <- NULL
SAMPLES  <- NULL
OUT      <- NULL
MAX_PAIRS <- 50000L
T1 <- 0.177
T2 <- 0.0884
T3 <- 0.0442

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if      (a == "--pairs"     && i < length(args)) { PAIRS    <- args[i + 1]; i <- i + 2L }
  else if (a == "--samples"   && i < length(args)) { SAMPLES  <- args[i + 1]; i <- i + 2L }
  else if (a == "--out"       && i < length(args)) { OUT      <- args[i + 1]; i <- i + 2L }
  else if (a == "--max_pairs" && i < length(args)) { MAX_PAIRS <- as.integer(args[i + 1]); i <- i + 2L }
  else if (a == "--t1"        && i < length(args)) { T1 <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--t2"        && i < length(args)) { T2 <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "--t3"        && i < length(args)) { T3 <- as.numeric(args[i + 1]); i <- i + 2L }
  else if (a == "-h" || a == "--help") { usage() }
  else { i <- i + 1L }
}

if (is.null(PAIRS) || is.null(SAMPLES) || is.null(OUT)) usage()
if (!file.exists(PAIRS))   stop("Pairs file not found: ", PAIRS)
if (!file.exists(SAMPLES)) stop("Samples file not found: ", SAMPLES)

# =============================================================================
# Load samples (defines cohort order; hub_id arrays are aligned to this order)
# =============================================================================
cat("[REL_EXPORT] Loading samples: ", SAMPLES, "\n", sep = "")
sample_ids <- as.character(fread(SAMPLES, header = FALSE)[[1]])
sample_ids <- sample_ids[nchar(sample_ids) > 0]
n_samples  <- length(sample_ids)
cat("[REL_EXPORT]   n_samples = ", n_samples, "\n", sep = "")
if (n_samples < 2) stop("Need at least 2 samples")

sample_to_idx <- setNames(seq_along(sample_ids) - 1L, sample_ids)   # 0-based

# =============================================================================
# Load pairs (auto-detect NAToRA 3-col vs ngsRelate .res)
# =============================================================================
cat("[REL_EXPORT] Loading pairs: ", PAIRS, "\n", sep = "")
# Sniff first line to decide format
con <- file(PAIRS, open = "r")
first_line <- readLines(con, n = 1L)
close(con)

raw <- tryCatch(fread(PAIRS), error = function(e) NULL)
if (is.null(raw) || nrow(raw) == 0) stop("Could not read pairs file")
cat("[REL_EXPORT]   raw rows: ", nrow(raw), ", cols: ",
    paste(names(raw), collapse = ", "), "\n", sep = "")

# Format detection
pairs_dt <- NULL
if (ncol(raw) == 3 && !any(grepl("^[A-Za-z_]", first_line))) {
  # NAToRA format: 3 cols, first line starts with sample id
  setnames(raw, c("id1", "id2", "theta"))
  raw[, theta := suppressWarnings(as.numeric(theta))]
  pairs_dt <- raw[is.finite(theta)]
  cat("[REL_EXPORT]   format: NAToRA 3-column\n")
} else if ("theta" %in% names(raw)) {
  # ngsRelate .res with proper header (theta column present)
  # ida/idb columns hold sample names if -z was used during ngsRelate, else
  # they're numeric indices. We try names first.
  id1_col <- if ("ida" %in% names(raw)) "ida" else NULL
  id2_col <- if ("idb" %in% names(raw)) "idb" else NULL
  if (is.null(id1_col) || is.null(id2_col)) {
    # Numeric index columns: a / b
    raw[, id1 := sample_ids[a + 1L]]
    raw[, id2 := sample_ids[b + 1L]]
  } else {
    raw[, id1 := raw[[id1_col]]]
    raw[, id2 := raw[[id2_col]]]
  }
  pairs_dt <- raw[, .(id1, id2, theta = as.numeric(theta))]
  pairs_dt <- pairs_dt[is.finite(theta)]
  cat("[REL_EXPORT]   format: ngsRelate .res with header\n")
} else if ("rab" %in% names(raw)) {
  # ngsRelate variant where 'rab' is the kinship metric and 'theta' missing
  # (this happens when ngsRelate is run without sample IDs, just indices)
  raw[, id1 := sample_ids[a + 1L]]
  raw[, id2 := sample_ids[b + 1L]]
  pairs_dt <- raw[, .(id1, id2, theta = as.numeric(rab))]
  pairs_dt <- pairs_dt[is.finite(theta)]
  cat("[REL_EXPORT]   format: ngsRelate .res indexed (using rab as theta)\n")
} else {
  stop("Could not detect pairs file format. Expected NAToRA 3-col or ngsRelate .res.")
}

cat("[REL_EXPORT]   parsed ", nrow(pairs_dt), " valid pairs\n", sep = "")

# Validate that id1 / id2 match samples list
n_unmatched_id1 <- sum(!(pairs_dt$id1 %in% sample_ids))
n_unmatched_id2 <- sum(!(pairs_dt$id2 %in% sample_ids))
if (n_unmatched_id1 > 0 || n_unmatched_id2 > 0) {
  cat("[REL_EXPORT] WARNING: ", n_unmatched_id1 + n_unmatched_id2,
      " pair endpoints don't match the sample list — those pairs will be dropped.\n", sep = "")
  pairs_dt <- pairs_dt[id1 %in% sample_ids & id2 %in% sample_ids]
}
if (nrow(pairs_dt) == 0) stop("No pairs survive sample-id matching. Check that the samples file matches the pairs file.")

# =============================================================================
# Distribution summary (for the user — sanity check before clustering)
# =============================================================================
cat("[REL_EXPORT] Theta distribution:\n")
n_first  <- sum(pairs_dt$theta >= T1)
n_second <- sum(pairs_dt$theta >= T2 & pairs_dt$theta < T1)
n_third  <- sum(pairs_dt$theta >= T3 & pairs_dt$theta < T2)
n_below  <- sum(pairs_dt$theta < T3)
cat(sprintf("[REL_EXPORT]   1st-degree pairs (theta >= %.4f):       %d\n", T1, n_first))
cat(sprintf("[REL_EXPORT]   2nd-degree pairs (%.4f <= theta < %.4f): %d\n", T2, T1, n_second))
cat(sprintf("[REL_EXPORT]   3rd-degree pairs (%.4f <= theta < %.4f): %d\n", T3, T2, n_third))
cat(sprintf("[REL_EXPORT]   below threshold (theta < %.4f):           %d\n", T3, n_third + n_below - n_third))

# =============================================================================
# Connected-components (union-find) at three thresholds
# =============================================================================
# Hub IDs are integer per-sample. Singletons (samples with no edges above
# threshold) get their own unique ID.
make_uf <- function(n) {
  parent <- as.integer(seq_len(n) - 1L)   # 0-based parents
  rank_  <- integer(n)
  list(parent = parent, rank_ = rank_)
}
uf_find <- function(uf, x) {
  while (uf$parent[x + 1L] != x) {
    uf$parent[x + 1L] <- uf$parent[uf$parent[x + 1L] + 1L]
    x <- uf$parent[x + 1L]
  }
  x
}
uf_union <- function(uf, a, b) {
  ra <- uf_find(uf, a)
  rb <- uf_find(uf, b)
  if (ra == rb) return(uf)
  # union by rank
  if (uf$rank_[ra + 1L] < uf$rank_[rb + 1L]) {
    uf$parent[ra + 1L] <- rb
  } else if (uf$rank_[ra + 1L] > uf$rank_[rb + 1L]) {
    uf$parent[rb + 1L] <- ra
  } else {
    uf$parent[rb + 1L] <- ra
    uf$rank_[ra + 1L]  <- uf$rank_[ra + 1L] + 1L
  }
  uf
}

build_hubs_at <- function(threshold) {
  uf <- make_uf(n_samples)
  edges <- pairs_dt[theta >= threshold]
  if (nrow(edges) > 0) {
    e1 <- sample_to_idx[edges$id1]
    e2 <- sample_to_idx[edges$id2]
    for (k in seq_len(nrow(edges))) {
      uf <- uf_union(uf, e1[k], e2[k])
    }
  }
  # Resolve roots and remap to dense 0-based hub IDs
  roots <- vapply(seq_len(n_samples) - 1L, function(x) uf_find(uf, x), integer(1))
  remap <- integer(0)
  hub_id <- integer(n_samples)
  next_id <- 0L
  for (i in seq_len(n_samples)) {
    r <- roots[i]
    key <- as.character(r)
    if (is.na(remap[key])) {
      remap[key] <- next_id
      next_id <- next_id + 1L
    }
    hub_id[i] <- remap[[key]]
  }
  list(hub_id = hub_id, n_hubs = next_id)
}

cat("[REL_EXPORT] Computing connected components at three thresholds...\n")
hubs1 <- build_hubs_at(T1)
hubs2 <- build_hubs_at(T2)
hubs3 <- build_hubs_at(T3)
cat("[REL_EXPORT]   1st-degree hubs: ", hubs1$n_hubs,
    " (", n_samples - hubs1$n_hubs, " samples in non-singleton hubs)\n", sep = "")
cat("[REL_EXPORT]   2nd-degree hubs: ", hubs2$n_hubs,
    " (", n_samples - hubs2$n_hubs, " samples in non-singleton hubs)\n", sep = "")
cat("[REL_EXPORT]   3rd-degree hubs: ", hubs3$n_hubs,
    " (", n_samples - hubs3$n_hubs, " samples in non-singleton hubs)\n", sep = "")

# Hub-size summary for the user
hub_size_summary <- function(hubs, label) {
  sizes <- table(hubs$hub_id)
  multi <- sizes[sizes >= 2]
  cat(sprintf("[REL_EXPORT]   %s: %d singletons, %d multi-hubs (sizes: min=%d, median=%.1f, max=%d)\n",
              label, sum(sizes == 1), length(multi),
              if (length(multi)) min(multi) else 0,
              if (length(multi)) median(multi) else 0,
              if (length(multi)) max(multi) else 0))
}
hub_size_summary(hubs1, "1st")
hub_size_summary(hubs2, "2nd")
hub_size_summary(hubs3, "3rd")

# =============================================================================
# Pack pairs (cap by theta to keep JSON size bounded)
# =============================================================================
# We sort pairs by theta descending and keep top MAX_PAIRS. Most relatedness-
# downstream views care most about high-theta edges. The scrubber UI for
# threshold sliders works on this capped set; below the cap, lower-theta
# pairs are mostly noise anyway.
setorder(pairs_dt, -theta)
n_pairs_keep <- min(MAX_PAIRS, nrow(pairs_dt))
pairs_kept <- pairs_dt[seq_len(n_pairs_keep)]
cat("[REL_EXPORT] Keeping top ", n_pairs_keep, " pairs by theta\n", sep = "")

# Map pair endpoints to sample indices for compactness
pairs_a <- sample_to_idx[pairs_kept$id1]
pairs_b <- sample_to_idx[pairs_kept$id2]
pairs_theta <- round(pairs_kept$theta, 5)

# =============================================================================
# Assemble JSON
# =============================================================================
# To match the scrubber's enrichment merge pattern (layer name -> object),
# everything goes under a top-level `relatedness:` key, with the layer
# metadata at the root and the data inside.

relatedness_layer <- list(
  n_samples      = as.integer(n_samples),
  samples        = sample_ids,
  thresholds     = list(t1 = T1, t2 = T2, t3 = T3),
  hub_id_1st     = as.integer(hubs1$hub_id),
  hub_id_2nd     = as.integer(hubs2$hub_id),
  hub_id_3rd     = as.integer(hubs3$hub_id),
  n_hubs         = list(n1 = as.integer(hubs1$n_hubs),
                         n2 = as.integer(hubs2$n_hubs),
                         n3 = as.integer(hubs3$n_hubs)),
  pair_count     = list(total = as.integer(nrow(pairs_dt)),
                        kept  = as.integer(n_pairs_keep),
                        n1    = as.integer(n_first),
                        n2    = as.integer(n_second),
                        n3    = as.integer(n_third)),
  pairs          = list(
    a     = as.integer(pairs_a),
    b     = as.integer(pairs_b),
    theta = pairs_theta
  )
)

out <- list(
  schema_version  = 2,
  `_layers_present` = list("relatedness"),
  relatedness     = relatedness_layer,
  `_generated_at` = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  `_generator`    = "export_relatedness_to_json_v1.R"
)

cat("[REL_EXPORT] Serializing JSON...\n")
t_ser <- proc.time()
json_str <- jsonlite::toJSON(out, auto_unbox = TRUE, na = "null",
                             digits = NA, pretty = FALSE)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(json_str, OUT)
cat("[REL_EXPORT]   serialized in ", round((proc.time() - t_ser)[3], 1), "s\n", sep = "")

f_size_kb <- round(file.info(OUT)$size / 1024, 1)
cat("\n[REL_EXPORT] === DONE ===\n")
cat("[REL_EXPORT] Output: ", OUT, "\n", sep = "")
cat("[REL_EXPORT] Size:   ", f_size_kb, " KB\n", sep = "")
cat("[REL_EXPORT] schema_version: 2\n")
cat("[REL_EXPORT] n_samples:      ", n_samples, "\n", sep = "")
cat("[REL_EXPORT] hubs:           1st=", hubs1$n_hubs,
    "  2nd=", hubs2$n_hubs, "  3rd=", hubs3$n_hubs, "\n", sep = "")
cat("[REL_EXPORT] pairs kept:     ", n_pairs_keep, " (of ", nrow(pairs_dt), " total)\n", sep = "")
cat("\n")
cat("[REL_EXPORT] Drop this JSON into the scrubber as a cohort-level enrichment.\n")
cat("[REL_EXPORT] The scrubber registers it as the 'relatedness' layer; the\n")
cat("[REL_EXPORT] family-hubs contingency view (next session) reads hub_id_1st/2nd/3rd\n")
cat("[REL_EXPORT] and renders cluster x hub partitions per L2 envelope.\n")
