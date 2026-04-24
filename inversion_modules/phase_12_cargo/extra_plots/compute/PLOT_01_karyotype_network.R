#!/usr/bin/env Rscript

# =============================================================================
# PLOT_01_karyotype_network.R
#
# Per-candidate sample-similarity network, coloured by karyotype.
#
# Replaces the "founder partition" mockup (image 2 panel F) with what's
# actually meaningful in this dataset: a force-directed network where each
# node is a sample and each edge is a high-similarity sample pair, computed
# from the per-candidate dosage matrix restricted to the inversion span.
# Nodes are coloured by HOM_REF / HET / HOM_INV from sample_karyotypes.tsv.
#
# Why this is meaningful:
#   - The within-candidate IBS naturally separates samples into the three
#     karyotype clusters, so the figure visually verifies what STEP_C01e/i
#     produces numerically.
#   - For candidates with sub-cluster structure (sub1/sub2/__noise from the
#     C01i d_seal step), nodes are additionally outlined by sub-cluster.
#
# Inputs:
#   ${DOSAGE_DIR}/<chr>.dosage.tsv.gz
#   ${DOSAGE_DIR}/<chr>.sites.tsv.gz
#   ${INVDIR}/<cand_dir>/data/sample_karyotypes.tsv
#   ${SNAKE_CAND_FILE}
#
# Output:
#   ${EXTRAS_FIG_DIR}/PLOT_01_karyotype_network__<cid>.pdf
#
# Usage:
#   Rscript PLOT_01_karyotype_network.R [candidate_id|all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(ggraph)
})

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[1] else "all"

DOSAGE_DIR      <- Sys.getenv("DOSAGE_DIR")
INVDIR          <- Sys.getenv("INVDIR")
SNAKE_CAND_FILE <- Sys.getenv("SNAKE_CAND_FILE")
SAMPLE_REGISTRY <- Sys.getenv("SAMPLE_REGISTRY")
EXTRAS_FIG_DIR  <- Sys.getenv("EXTRAS_FIG_DIR")
SAMPLE_LIST     <- Sys.getenv("SAMPLE_LIST")

stopifnot(nzchar(DOSAGE_DIR), nzchar(INVDIR),
          nzchar(SNAKE_CAND_FILE), nzchar(EXTRAS_FIG_DIR))

KARY_PAL <- c(REF = "#1f77b4", HET = "#7f7f7f", INV = "#d62728",
              HOM_REF = "#1f77b4", HOM_INV = "#d62728",
              UNKNOWN = "grey85")

cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))

# ── Read sample list (defines column order in dosage matrix) ──
sample_names <- NULL
if (nzchar(SAMPLE_LIST) && file.exists(SAMPLE_LIST)) {
  sample_names <- as.character(fread(SAMPLE_LIST, header = FALSE)[[1]])
}

# ── Helpers ──────────────────────────────────────────────────────────────────

read_dosage_in_window <- function(chrom, start_bp, end_bp) {
  dos_f   <- file.path(DOSAGE_DIR, paste0(chrom, ".dosage.tsv.gz"))
  sites_f <- file.path(DOSAGE_DIR, paste0(chrom, ".sites.tsv.gz"))
  if (!file.exists(dos_f) || !file.exists(sites_f)) {
    message("  [skip] missing dosage / sites for ", chrom)
    return(NULL)
  }
  sites <- fread(cmd = paste0("zcat ", shQuote(sites_f)),
                  header = FALSE,
                  col.names = c("chrom", "pos"))
  in_win <- which(sites$chrom == chrom & sites$pos >= start_bp & sites$pos < end_bp)
  if (length(in_win) < 5) {
    message("  [skip] <5 sites in window")
    return(NULL)
  }
  con <- gzfile(dos_f, "r")
  on.exit(close(con), add = TRUE)
  m <- as.matrix(fread(cmd = paste0("zcat ", shQuote(dos_f)),
                        header = FALSE))
  m <- m[in_win, , drop = FALSE]
  if (!is.null(sample_names) && ncol(m) == length(sample_names)) {
    colnames(m) <- sample_names
  } else {
    colnames(m) <- paste0("S", seq_len(ncol(m)))
  }
  m
}

ibs_distance <- function(dos_mat) {
  # Pairwise Euclidean distance on dosage, normalised by number of valid sites.
  # Works on small samples × variable-sites matrices — typical inversion has
  # 50–500 polymorphic sites, 50–80 samples in a karyotype, so this is fast.
  X <- t(dos_mat)
  X[is.na(X)] <- 1  # impute mean dosage 1 for missing
  D <- as.matrix(dist(X, method = "manhattan"))
  D <- D / ncol(X)  # per-site mean Manhattan distance
  D
}

build_network <- function(D, k_neighbors = 4) {
  # k-nearest-neighbours network: each node connects to its k closest others.
  # Symmetrise (an edge exists if either endpoint listed the other as a neighbour).
  n <- nrow(D)
  if (n < 3) return(NULL)
  diag(D) <- Inf
  k <- min(k_neighbors, n - 1)
  adj <- matrix(0L, n, n, dimnames = dimnames(D))
  for (i in seq_len(n)) {
    nn <- order(D[i, ])[seq_len(k)]
    adj[i, nn] <- 1L
  }
  adj <- pmax(adj, t(adj))
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  # Edge weights = inverse distance (close pairs draw together in layout)
  E(g)$weight <- vapply(seq_len(ecount(g)), function(e) {
    ends <- ends(g, e)
    1 / (D[ends[1], ends[2]] + 0.01)
  }, numeric(1))
  g
}

# ── Per-candidate driver ─────────────────────────────────────────────────────

process_one <- function(cid) {
  cand <- cands[candidate_id == cid]
  if (nrow(cand) == 0) {
    message("[PLOT_01] candidate not found: ", cid); return(invisible())
  }
  chr <- cand$chrom
  s   <- as.integer(cand$start_bp)
  e   <- as.integer(cand$end_bp)

  # Karyotype labels
  cand_dirs <- list.files(INVDIR, pattern = paste0("candidate_", cid, "$"),
                           recursive = TRUE, include.dirs = TRUE,
                           full.names = TRUE)
  kfile <- NULL
  for (cd in cand_dirs) {
    fp <- file.path(cd, "data", "sample_karyotypes.tsv")
    if (file.exists(fp)) { kfile <- fp; break }
  }
  if (is.null(kfile)) {
    message("[PLOT_01] no karyotype file for ", cid); return(invisible())
  }
  kdt <- fread(kfile)
  if (!"sample" %in% names(kdt) || !"karyotype" %in% names(kdt)) {
    message("[PLOT_01] karyotype file missing required cols for ", cid)
    return(invisible())
  }

  # Dosage in window
  message("[PLOT_01] ", cid, " — reading dosage on ", chr, ":",
          format(s, big.mark = ","), "-", format(e, big.mark = ","))
  m <- read_dosage_in_window(chr, s, e)
  if (is.null(m)) return(invisible())

  # Restrict to samples that have karyotype calls
  keep <- intersect(colnames(m), kdt$sample)
  if (length(keep) < 6) {
    message("[PLOT_01] ", cid, ": <6 samples with karyotype + dosage — skip")
    return(invisible())
  }
  m <- m[, keep, drop = FALSE]
  kdt_sub <- kdt[sample %in% keep][match(keep, sample)]

  D <- ibs_distance(m)
  g <- build_network(D, k_neighbors = 4)
  if (is.null(g)) {
    message("[PLOT_01] ", cid, ": graph too small"); return(invisible())
  }

  V(g)$karyotype <- kdt_sub$karyotype
  V(g)$sample    <- kdt_sub$sample

  set.seed(42 + nchar(cid))
  layout_xy <- create_layout(g, layout = "fr", weights = E(g)$weight)

  p <- ggraph(layout_xy) +
    geom_edge_link(alpha = 0.18, edge_width = 0.4, color = "grey50") +
    geom_node_point(aes(color = karyotype), size = 3.5, alpha = 0.9) +
    scale_color_manual(values = KARY_PAL, na.value = "grey85", drop = FALSE) +
    labs(title = paste0("Sample IBS network — ", cid),
         subtitle = paste0(chr, ":", format(s, big.mark = ","), "–",
                            format(e, big.mark = ","),
                            "   |   ", length(keep), " samples,  ",
                            nrow(m), " polymorphic sites"),
         color = "karyotype") +
    theme_void() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  outpdf <- file.path(EXTRAS_FIG_DIR,
                       paste0("PLOT_01_karyotype_network__", cid, ".pdf"))
  ggsave(outpdf, p, width = 8, height = 8, device = cairo_pdf)
  message("  wrote ", outpdf)
}

# ── Main ─────────────────────────────────────────────────────────────────────

if (target == "all") {
  cids <- cands$candidate_id
} else {
  cids <- target
}
message("[PLOT_01] ", length(cids), " candidate(s)")
for (cid in cids) {
  tryCatch(process_one(cid),
           error = function(e) message("[PLOT_01] ERROR on ", cid, ": ", conditionMessage(e)))
}
message("[PLOT_01] Done")
