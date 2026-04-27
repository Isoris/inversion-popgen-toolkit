#!/usr/bin/env Rscript
# ============================================================================
# STEP_D16_tree_node_pca_facets.R  v3.2
# ----------------------------------------------------------------------------
# v3.2 vs v3.1
# ------------
# - New per-node "halves coherence" metric: split the node windows into
#   first half / second half by chromosome position, average per-window PC1
#   independently in each half, Pearson-correlate the two 226-sample
#   vectors. A coherent inversion block holds coherence ~ 1.0 because every
#   window inside sorts samples the same way. An over-merged composite
#   (the tree drew the rectangle too wide and lumped several adjacent
#   blocks into one node) drops to coherence < 0.6.
# - New badge [INCOHERENT] when coherence < 0.60. This fires when the
#   "wrong square" failure mode applies — node geometry is wrong, not
#   biology.
# - New 4th panel in _summary.pdf: coherence vs nW with the threshold line.
# - node_summary.tsv now carries `coherence` and `incoherent` columns.
#
# This addresses the case where bands look "many" not because the
# population biology is high-K, but because the tree merged adjacent
# distinct blocks into one rectangle. Distinct from FAM-LD? (high
# fam_purity) which is the rare-allele family-LD failure mode.
#
# v3.1 vs v3
# ----------
# - New flag --bamlist <txt>: maps Ind0..IndN -> CGAxxx using one-per-line
#   bam list (same order as the precomp's PC_1_Ind* columns).
# - New flag --family <tsv>: two-column CGA<TAB>family_id (e.g. natora's
#   *_familyList.txt). When given, emits a SECOND PDF (_family.pdf) with
#   the rug colored by family ID instead of K-cluster, and adds a per-node
#   "family purity" metric to the facet title.
# - fam_purity: the size-weighted mean across K-clusters of "fraction of
#   samples in this cluster that come from its single most-common family."
#     ~1.0  -> bands ARE families => family LD, NOT inversion (FAM-LD? badge)
#     ~0.2  -> bands span many families => real inversion-like signal
# - The _node_summary.tsv now carries fam_purity and family_ld columns.
#
# v3 vs v2
# --------
# v3 replaced fixed k=3 with adaptive density-gated silhouette k=2..7 and
# added the over-merge diagnostic plot (best_k vs nW with COMPOSITE? flags).
# Default layout is now 6x4 (24 facets per page).
#
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# One PC1 strip-plot per tree node, all on faceted multi-page PDFs.
# This is the "three stripes" diagnostic: a real inversion shows three
# clean clusters along PC1 (HOM_REF / HET / HOM_INV); family LD shows a
# smear; no signal shows a blob.
#
# Why this exists
# ---------------
# D13d gives you the heatmap with tree-node outlines. That tells you WHERE
# the algorithm thinks blocks are. It does NOT confirm that those blocks
# are real inversions vs family LD vs noise. The sim_mat numbers
# (e.g. 0.99 inside, 0.77 across) are objective but abstract — you can't
# look at a number and decide whether to trust it.
#
# The local-PCA strip plot, on the other hand, shows you the SHAPE of the
# population structure inside the block. Three stripes = inversion.
# Smear = family LD. Single tight blob = no signal.
# This is the eye-test that the numbers can't replace.
#
# Why no PC2
# ----------
# The slim precomp on disk only carries PC_1_* per-sample columns. PC2 was
# never stored at the per-window level in this codebase (only the PC2
# eigenvalue lam_2 is, in the dt). So we plot PC1 along x, with a tiny
# vertical jitter to spread overlapping samples — three stripes show as
# three vertical clouds along x. Same diagnostic value as PC1 vs PC2.
#
# Inputs
# ------
#   --precomp        precomp.slim.rds  (provides $dt with PC_1_* columns)
#   --tree           nn_tree_<chr>.tsv (D09 output)
#   --samples        TSV with columns: ind  cga  ancestry   (optional)
#                    Used to color points by ancestry, like STEP27.
#                    If not given, colors come from k=3 k-means on PC1.
#   --tree_filter    "inversion_only" | "candidate_or_better" | "weak_or_better" | "all"
#                    (default: candidate_or_better)
#   --min_width      drop nodes narrower than this many windows (default 10)
#   --outdir         output directory
#   --chr            chromosome label (for filenames)
#   --ncol           facets per row (default 5)
#   --nrow           facets per page row (default 8 → 40 nodes per page)
#   --page_width_in  PDF width in inches (default 14)
#   --page_height_in PDF height in inches (default 18)
#
# Output
# ------
#   <outdir>/<chr>_tree_pca_facets_<filter>_bands.pdf       — k=3 band coloring
#   <outdir>/<chr>_tree_pca_facets_<filter>_ancestry.pdf    — ancestry coloring (if --samples)
#
# Per-facet content
# -----------------
#   x = averaged PC1 across all precomp windows in the node
#   y = vertical jitter (Uniform[-0.5, +0.5]; no biological meaning,
#       just so points don't all stack on one horizontal line)
#   color = k-means k=3 band on PC1 (Band1=blue, Band2=grey, Band3=amber)
#           Bands ordered along PC1: lowest centroid = Band1, highest = Band3.
#           In a real inversion these correspond to HOM_REF, HET, HOM_INV
#           (or vice-versa — the orientation of PC1 is arbitrary).
#   facet label = node_id, classification, span Mb, n windows
#
# Reading the figure
# ------------------
#   Three clean vertical clouds along x  -> real inversion
#   One single cloud (smear)             -> family LD or nothing biological
#   Two clumps                            -> two-population structure, not inv
#   Dots smeared along entire x range
#     with no visible cluster gaps        -> noisy signal — needs deeper look
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# ---- CLI --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

precomp_f    <- get_arg("--precomp")
tree_f       <- get_arg("--tree")
samples_f    <- get_arg("--samples")
bamlist_f    <- get_arg("--bamlist")    # one CGA ID per line, precomp order
family_f     <- get_arg("--family")     # CGA_id<TAB>family_id (natora style)
pairs_f      <- get_arg("--pairs")      # raw pairwise relatedness file
                                         # 3 cols: id1 id2 theta
theta_cutoff <- as.numeric(get_arg("--theta_cutoff", "0.177"))
                                         # default = 1st-degree threshold
                                         # used for connected-component
                                         # family detection from --pairs
outdir       <- get_arg("--outdir", ".")
chr_label    <- get_arg("--chr", "chr")
tree_filter  <- get_arg("--tree_filter", "candidate_or_better")
min_width    <- as.integer(get_arg("--min_width", "10"))
ncol_f       <- as.integer(get_arg("--ncol", "6"))
nrow_f       <- as.integer(get_arg("--nrow", "4"))
page_w       <- as.numeric(get_arg("--page_width_in", "16"))
page_h       <- as.numeric(get_arg("--page_height_in", "12"))
free_scales  <- as.logical(get_arg("--free_scales", "FALSE"))
k_max        <- as.integer(get_arg("--k_max", "7"))     # adaptive k tries 2..k_max
k_min        <- as.integer(get_arg("--k_min", "2"))

if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("--precomp is required and must exist")
if (is.na(tree_f) || !file.exists(tree_f))
  stop("--tree is required and must exist")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

panels_per_page <- ncol_f * nrow_f

cat("[D16] precomp:    ", precomp_f, "\n")
cat("[D16] tree:       ", tree_f, "\n")
cat("[D16] samples:    ", samples_f %||% "(none)", "\n")
cat("[D16] tree_filter:", tree_filter, "\n")
cat("[D16] facets:     ", ncol_f, "x", nrow_f, "= ",
    panels_per_page, " per page\n", sep = "")
cat("[D16] adaptive k: range ", k_min, "..", k_max, " (silhouette-selected)\n", sep = "")

# ---- Load precomp -----------------------------------------------------------

pc <- readRDS(precomp_f)
if (is.null(pc$dt)) stop("[D16] precomp lacks $dt")
dt <- as.data.table(pc$dt)
N <- nrow(dt)

# Detect PC columns. Only PC1 is required — the slim precomp on disk
# typically does NOT carry per-sample PC2 columns (only the PC2 eigenvalue
# lam_2 lives in dt). We plot PC1 with a small vertical jitter instead of
# a true PC1 vs PC2 scatter — same diagnostic content (three stripes).
pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
pc2_cols <- grep("^PC_2_", names(dt), value = TRUE)
if (length(pc1_cols) == 0)
  stop("[D16] precomp$dt has no PC_1_* columns. ",
       "Did you load the precomp.slim or the full precomp?")

have_pc2 <- length(pc2_cols) == length(pc1_cols)
if (!have_pc2) {
  cat("[D16] no PC_2_* columns in precomp$dt (slim precomp) — ",
      "falling back to PC1-strip layout (jittered y)\n", sep = "")
}

sample_ids <- sub("^PC_1_", "", pc1_cols)
n_samp <- length(sample_ids)
cat("[D16] N windows: ", N, " | N samples: ", n_samp, "\n", sep = "")

# Extract PC matrices once (windows x samples). PC1 is mandatory; PC2 only
# if available. ~ 4302 * 226 * 8 bytes = ~7.8 MB per matrix.
pc1_mat <- as.matrix(dt[, ..pc1_cols])
storage.mode(pc1_mat) <- "double"
if (have_pc2) {
  pc2_mat <- as.matrix(dt[, ..pc2_cols])
  storage.mode(pc2_mat) <- "double"
} else {
  pc2_mat <- NULL
}

# ---- Load samples metadata (optional) ---------------------------------------

sample_meta <- data.table(ind = sample_ids, cga = sample_ids, ancestry = "unknown")
if (!is.na(samples_f) && file.exists(samples_f)) {
  cat("[D16] loading sample metadata\n")
  smeta <- fread(samples_f, header = TRUE)
  setnames(smeta, tolower(names(smeta)))
  if ("ind" %in% names(smeta)) {
    sample_meta <- merge(sample_meta[, .(ind)], smeta,
                         by = "ind", all.x = TRUE, sort = FALSE)
    if (!"ancestry" %in% names(sample_meta)) sample_meta[, ancestry := "unknown"]
    if (!"cga" %in% names(sample_meta)) sample_meta[, cga := ind]
    sample_meta[is.na(ancestry), ancestry := "unknown"]
    sample_meta[is.na(cga), cga := ind]
  } else {
    warning("[D16] samples file lacks 'ind' column — using IDs as-is")
  }
}

# ---- Bamlist remap (Ind0..IndN -> CGAxxx) ---------------------------------
# The precomp's PC_1_Ind* columns use sequential indices in bamfile order.
# The bamlist file is one CGA ID per line, in the SAME order. Row 1 of
# bamlist = Ind0, row 2 = Ind1, etc.
have_bamlist <- !is.na(bamlist_f) && file.exists(bamlist_f)
if (have_bamlist) {
  cat("[D16] loading bamlist for Ind -> CGA remap\n")
  bam_ids <- readLines(bamlist_f)
  bam_ids <- bam_ids[nzchar(trimws(bam_ids))]
  if (length(bam_ids) != n_samp) {
    warning(sprintf(
      "[D16] bamlist has %d entries but precomp has %d samples — using bamlist length min",
      length(bam_ids), n_samp))
  }
  ncommon <- min(length(bam_ids), n_samp)
  ind_to_cga <- setNames(bam_ids[seq_len(ncommon)],
                          paste0("Ind", seq_len(ncommon) - 1L))
  # Update sample_meta and sample_ids order: keep the precomp Ind ordering
  # but rename to the CGA equivalents for display.
  sample_meta[, cga := ind_to_cga[ind]]
  sample_meta[is.na(cga), cga := ind]
} else {
  cat("[D16] no --bamlist given; samples will be labeled with Ind*\n")
}

# ---- Family info (optional) ------------------------------------------------
# Two ways to get family IDs, in priority order:
#
# 1. --pairs <file> + --theta_cutoff <theta>: build a graph from the raw
#    pairwise relatedness file (cols: id1 id2 theta), keep edges with
#    theta >= theta_cutoff, find connected components. Each component =
#    one family. This reproduces the network figure's components when
#    theta_cutoff = 0.177. RECOMMENDED — matches the visible figure
#    structure exactly (e.g. the 16-member C8 hub appears as one family,
#    not as 13 singletons + 3 retained as natora's pruning logic produces).
#
# 2. --family <file>: a two-column TSV (CGA_id<TAB>family_id) like
#    natora's *_familyList.txt. Use this if you want natora's pruning-
#    aware labels rather than the raw graph components.
#
# Both produce a `family_id` column on sample_meta. Unassigned/unmatched
# samples get family_id = -1L.

build_components_from_pairs <- function(pairs_path, theta_thr, all_ids) {
  # Read pairs (3 cols: id1 id2 theta), keep edges with theta >= cutoff,
  # find connected components by union-find. Returns a data.table with
  # columns (cga, family_id). Singletons get unique family_ids too.
  cat("[D16]   reading pairs file: ", pairs_path, "\n", sep = "")
  # File may or may not have a header row depending on how it was written.
  # Sniff first byte: if column 3 of row 1 is numeric, no header. Otherwise
  # treat row 1 as header. fread's auto detection sometimes gets this wrong,
  # so be explicit.
  first_line <- readLines(pairs_path, n = 1L)
  parts <- strsplit(first_line, "[\t ]+")[[1]]
  has_header <- !suppressWarnings(is.finite(as.numeric(parts[length(parts)])))
  pairs_dt <- fread(pairs_path,
                    header = has_header,
                    col.names = c("id1", "id2", "theta"))
  # Coerce theta to numeric in case fread got it as character
  pairs_dt[, theta := as.numeric(theta)]
  cat("[D16]   total pairs in file: ", nrow(pairs_dt),
      " (header_detected=", has_header, ")\n", sep = "")

  # Always print theta distribution so the user can verify the third
  # column actually IS theta (not some other relatedness metric like rab).
  qs <- quantile(pairs_dt$theta, c(0.5, 0.9, 0.95, 0.99, 1.0), na.rm = TRUE)
  cat("[D16]   theta distribution: ",
      "median=", round(qs[1], 4),
      " p90=", round(qs[2], 4),
      " p95=", round(qs[3], 4),
      " p99=", round(qs[4], 4),
      " max=", round(qs[5], 4),
      "\n", sep = "")
  cat("[D16]   (sanity check: 1st-degree pairs should have theta ~0.25,",
      " 2nd-degree ~0.125, unrelated ~0)\n")

  # NB: filter the data.table by an explicit external variable name.
  # Earlier the function param was also called `theta`, which made
  # `pairs_dt[theta >= !!theta]` ambiguous in data.table NSE — it ended
  # up keeping 5400/25425 rows instead of the ~400 with theta >= 0.177.
  thr <- theta_thr
  pairs_dt <- pairs_dt[theta >= thr]
  cat("[D16]   pairs above theta_cutoff=", thr, ": ", nrow(pairs_dt), "\n", sep = "")

  # Universe = union of all_ids and any IDs in the kept pairs (some samples
  # in all_ids may have NO edges and become singleton components).
  all_ids_in_pairs <- unique(c(pairs_dt$id1, pairs_dt$id2))
  universe <- unique(c(all_ids, all_ids_in_pairs))
  N <- length(universe)
  id2idx <- setNames(seq_len(N), universe)

  # Integer-indexed parent vector — fast union-find with path compression
  parent <- seq_len(N)

  find_root <- function(x) {
    while (parent[x] != x) {
      parent[parent[x]] <<- parent[parent[x]]   # path halving
      x <- parent[x]
    }
    x
  }
  union_xy <- function(a, b) {
    ra <- find_root(a); rb <- find_root(b)
    if (ra != rb) parent[ra] <<- rb
  }

  if (nrow(pairs_dt) > 0) {
    a_idx <- id2idx[pairs_dt$id1]
    b_idx <- id2idx[pairs_dt$id2]
    valid <- !is.na(a_idx) & !is.na(b_idx)
    a_idx <- a_idx[valid]; b_idx <- b_idx[valid]
    for (i in seq_along(a_idx)) {
      union_xy(a_idx[i], b_idx[i])
    }
  }

  # Read out roots
  roots <- vapply(seq_len(N), find_root, integer(1))
  # Stable family_id ordering: family 1 = biggest component
  comp_size <- table(roots)
  ordered_roots <- as.integer(names(sort(comp_size, decreasing = TRUE)))
  fam_id_map <- setNames(seq_along(ordered_roots), ordered_roots)
  family_id <- fam_id_map[as.character(roots)]

  out <- data.table(cga = universe, family_id = as.integer(family_id))
  out <- out[cga %in% all_ids]
  out
}

have_pairs  <- !is.na(pairs_f)  && file.exists(pairs_f)
have_family <- !is.na(family_f) && file.exists(family_f)

if (have_pairs) {
  cat("[D16] loading family info from pairs file (theta >= ", theta_cutoff, ")\n",
      sep = "")
  fam_dt <- build_components_from_pairs(
    pairs_f, theta_cutoff,
    all_ids = sample_meta$cga
  )
  if (have_bamlist) {
    sample_meta <- merge(sample_meta, fam_dt, by = "cga", all.x = TRUE,
                          sort = FALSE)
  } else {
    setnames(fam_dt, "cga", "ind")
    sample_meta <- merge(sample_meta, fam_dt, by = "ind", all.x = TRUE,
                          sort = FALSE)
  }
  sample_meta[is.na(family_id), family_id := -1L]
  fam_size_quick <- sample_meta[family_id != -1L,
                                .N, by = family_id][order(-N)]
  n_hub_quick <- sum(fam_size_quick$N >= 4L)
  cat("[D16]   total components: ", nrow(fam_size_quick),
      "  |  hubs (>=4 samples): ", n_hub_quick, "\n", sep = "")
  if (n_hub_quick > 0) {
    cat("[D16]   biggest components: ",
        paste0("F", fam_size_quick$family_id[1:min(5, n_hub_quick)],
               "(n=", fam_size_quick$N[1:min(5, n_hub_quick)], ")",
               collapse = ", "),
        "\n", sep = "")
  }
  have_family <- TRUE   # downstream logic checks this flag
} else if (have_family) {
  cat("[D16] loading family info from --family file (natora-style)\n")
  fam_dt <- fread(family_f, header = FALSE,
                  col.names = c("cga", "family_id"))
  if (have_bamlist) {
    sample_meta <- merge(sample_meta, fam_dt, by = "cga", all.x = TRUE,
                          sort = FALSE)
  } else {
    setnames(fam_dt, "cga", "ind")
    sample_meta <- merge(sample_meta, fam_dt, by = "ind", all.x = TRUE,
                          sort = FALSE)
  }
  sample_meta[is.na(family_id), family_id := -1L]
  cat("[D16]   ", length(unique(sample_meta$family_id)),
      " distinct family IDs (incl. unmatched=-1)\n", sep = "")
} else {
  sample_meta[, family_id := -1L]
}

# Make sure sample_meta is in precomp order (Ind0, Ind1, ...) for downstream
# vectorized operations that assume positional alignment.
sample_meta <- sample_meta[match(sample_ids, ind)]

# ---- Load tree --------------------------------------------------------------

tree <- fread(tree_f)
required <- c("node_id", "start", "end", "start_mb", "end_mb", "classification")
miss <- setdiff(required, names(tree))
if (length(miss) > 0)
  stop("[D16] tree TSV missing columns: ", paste(miss, collapse = ", "))

keep_classes <- switch(tree_filter,
  "inversion_only"      = c("INVERSION"),
  "candidate_or_better" = c("INVERSION", "CANDIDATE"),
  "weak_or_better"      = c("INVERSION", "CANDIDATE", "WEAK_CANDIDATE"),
  "all"                 = c("INVERSION", "CANDIDATE", "WEAK_CANDIDATE", "FAMILY_LD"),
  stop("[D16] unknown --tree_filter: ", tree_filter)
)

tree <- tree[classification %in% keep_classes]
tree[, width := end - start + 1L]
tree <- tree[width >= min_width]

# Sort by chromosome position so PDF pages walk left-to-right along chromosome
setorder(tree, start_mb)

cat("[D16] tree nodes after filter (", tree_filter, ", min_width=",
    min_width, "): ", nrow(tree), "\n", sep = "")

if (nrow(tree) == 0) {
  cat("[D16] no tree nodes pass the filter — nothing to plot. exiting.\n")
  quit(status = 0)
}

# ---- Build per-node averaged PC scatter --------------------------------------
# For each node, average PC1 / PC2 across all precomp windows whose row
# index falls inside [start, end]. Then k=3 k-means on PC1 to assign bands.

CLASS_COLOR <- c(
  INVERSION       = "#7f1d1d",
  CANDIDATE       = "#b45309",
  WEAK_CANDIDATE  = "#0d9488",
  FAMILY_LD       = "#94a3b8"
)

# K-palette: cluster 1 = lowest PC1, cluster best_k = highest PC1.
# Diverging palette so high-K and low-K are visually distinguishable.
BAND_COLOR <- c(
  K1 = "#1e3a8a",   # deep blue (lowest PC1)
  K2 = "#4fa3ff",   # light blue
  K3 = "#67e8f9",   # cyan
  K4 = "#b8b8b8",   # grey (middle)
  K5 = "#facc15",   # yellow
  K6 = "#f5a524",   # amber
  K7 = "#dc2626",   # red (highest PC1)
  `?` = "#475569"
)

# ---- Diagnostic helpers -----------------------------------------------------

# Pearson bimodality coefficient. BC > 0.555 -> non-unimodal.
# Uses sample skewness and excess kurtosis with bias correction.
# References: SAS/STAT BIMODALITY chapter; Pfister et al. 2013.
bimodality_coef <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4L) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  z <- (x - m) / s
  skew <- mean(z^3) * sqrt(n * (n - 1)) / (n - 2)
  kurt <- (n - 1) / ((n - 2) * (n - 3)) *
          ((n + 1) * (mean(z^4) - 3) + 6)
  bc_corr <- 3 * (n - 1)^2 / ((n - 2) * (n - 3))
  (skew^2 + 1) / (kurt + bc_corr)
}

# Count peaks in a 1D density estimate. Coarse — just counts local maxima
# in the density curve, ignoring tiny ripples below max_density * peak_floor.
# bw_adj < 1 narrows the kernel (more peaks visible); bw_adj > 1 broadens.
density_peak_count <- function(x, peak_floor = 0.05, n = 256L, bw_adj = 1.0) {
  x <- x[is.finite(x)]
  if (length(x) < 5L) return(NA_integer_)
  d <- tryCatch(stats::density(x, n = n, adjust = bw_adj), error = function(e) NULL)
  if (is.null(d)) return(NA_integer_)
  y <- d$y
  thresh <- max(y) * peak_floor
  # Local maxima with neighbor comparison
  is_peak <- c(FALSE,
               y[-c(1, length(y))] > y[-c(length(y) - 1, length(y))] &
               y[-c(1, length(y))] > y[-c(1, 2)],
               FALSE)
  sum(is_peak & y >= thresh)
}

# Adaptive k via density-gated silhouette in 1D.
# Logic:
#   1. Wide-bandwidth (Silverman default) density must have >= 2 peaks
#      to even consider multimodality. Single-mode noise stays k=1.
#   2. Fine-bandwidth density (adjust=0.7) sets the candidate k. This
#      sees finer structure but is sanity-checked by step 1.
#   3. Silhouette of kmeans(centers = candidate_k) must exceed
#      sil_threshold (default 0.45) to ratify the split.
# If any gate fails, return k=1 (= "essentially unimodal").
# Cluster ids are reordered 1..best_k by centroid PC1 ascending so the
# color palette is stable across facets (cluster 1 always = lowest PC1).
adaptive_k <- function(x, k_min = 2L, k_max = 7L, sil_threshold = 0.45) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 6L) return(list(best_k = 1L, best_sil = NA_real_,
                           cluster = rep(1L, n)))
  npk_wide <- density_peak_count(x, bw_adj = 1.0)
  npk_fine <- density_peak_count(x, bw_adj = 0.7)
  # Gate 1: must be multimodal at default bandwidth
  if (is.na(npk_wide) || npk_wide < 2L) {
    return(list(best_k = 1L, best_sil = NA_real_, cluster = rep(1L, n)))
  }
  cand_k <- min(max(npk_fine, npk_wide, k_min), k_max)
  if (cand_k >= n) cand_k <- max(2L, n - 1L)

  km <- tryCatch(kmeans(x, centers = cand_k, nstart = 10L, iter.max = 50),
                 error = function(e) NULL)
  if (is.null(km)) {
    return(list(best_k = 1L, best_sil = NA_real_, cluster = rep(1L, n)))
  }
  cl <- km$cluster

  # Mean silhouette in 1D
  d_full <- as.matrix(dist(x))
  sil_vec <- numeric(n)
  cluster_ids <- sort(unique(cl))
  for (i in seq_len(n)) {
    own <- which(cl == cl[i])
    own <- setdiff(own, i)
    if (length(own) == 0L) { sil_vec[i] <- 0; next }
    a_i <- mean(d_full[i, own])
    other_ids <- setdiff(cluster_ids, cl[i])
    b_candidates <- vapply(other_ids, function(oc) {
      mean(d_full[i, cl == oc])
    }, numeric(1))
    b_i <- min(b_candidates)
    sil_vec[i] <- (b_i - a_i) / max(a_i, b_i)
  }
  mean_sil <- mean(sil_vec)

  # Gate 3: silhouette must clear threshold
  if (is.na(mean_sil) || mean_sil < sil_threshold) {
    return(list(best_k = 1L, best_sil = mean_sil, cluster = rep(1L, n)))
  }

  # Re-label clusters by centroid PC1 ascending
  centroids <- km$centers[, 1]
  perm <- order(centroids)
  new_cl <- match(km$cluster, perm)
  list(best_k = cand_k, best_sil = mean_sil, cluster = new_cl)
}

# Stable per-sample y-jitter (same value for the same sample across all
# panels — easier on the eye than fresh randomness per facet). When PC2
# is available it overrides this. Range chosen to spread 226 samples
# without collisions.
set.seed(42L)
y_jitter <- runif(n_samp, min = -0.5, max = 0.5)

panels <- vector("list", nrow(tree))
for (ni in seq_len(nrow(tree))) {
  s <- tree$start[ni]
  e <- tree$end[ni]
  s <- as.integer(max(1L, s))
  e <- as.integer(min(N, e))
  if (s > e) next

  rows <- s:e
  # Average across windows in node — colMeans on the rows of the matrix.
  # Some windows may be all-NA for a sample (depth dropouts); na.rm=TRUE
  # handles it. If a sample is NA in EVERY row of this node, we propagate
  # NA so the dot is dropped.
  avg_pc1 <- colMeans(pc1_mat[rows, , drop = FALSE], na.rm = TRUE)
  avg_pc1[is.nan(avg_pc1)] <- NA_real_

  if (have_pc2) {
    avg_pc2 <- colMeans(pc2_mat[rows, , drop = FALSE], na.rm = TRUE)
    avg_pc2[is.nan(avg_pc2)] <- NA_real_
  } else {
    # Jittered y: no biological meaning, just spreads out points.
    avg_pc2 <- y_jitter
  }

  valid <- is.finite(avg_pc1) & is.finite(avg_pc2)
  if (sum(valid) < 10L) next   # need a minimum

  # Adaptive k via density-gated silhouette — pick the k in [k_min..k_max]
  # that best explains the averaged PC1 distribution. Cluster ids 1..best_k
  # are ordered by centroid PC1 ascending so colors stay consistent across
  # facets (cluster 1 = lowest PC1, cluster best_k = highest PC1).
  ak <- adaptive_k(avg_pc1[valid],
                   k_min = k_min,
                   k_max = k_max,
                   sil_threshold = 0.45)
  band_id <- rep(NA_integer_, n_samp)
  band_id[valid] <- ak$cluster
  band <- ifelse(is.na(band_id), "?", paste0("K", band_id))

  cls    <- as.character(tree$classification[ni])
  nid    <- tree$node_id[ni]
  span_mb <- round(tree$end_mb[ni] - tree$start_mb[ni], 2)
  nwin    <- tree$width[ni]

  bc       <- bimodality_coef(avg_pc1[valid])
  n_peaks  <- density_peak_count(avg_pc1[valid])
  bc_lab   <- if (is.finite(bc)) sprintf("BC=%.2f", bc) else "BC=NA"
  pk_lab   <- if (!is.na(n_peaks)) sprintf("peaks=%d", n_peaks) else "peaks=?"
  k_lab    <- sprintf("best_k=%d", ak$best_k)
  sil_lab  <- if (is.finite(ak$best_sil)) sprintf("sil=%.2f", ak$best_sil) else "sil=NA"

  # Halves coherence: split node windows into halves by chromosome position,
  # compute averaged PC1 of each half independently across the 226 samples,
  # then Pearson-correlate the two vectors.
  #
  # A coherent biological block (one inversion or one consistent signal):
  # both halves sort samples the same way along PC1 -> coherence ~ 1.0.
  # An over-merged composite (the tree drew the rectangle too wide and now
  # encloses two independent structural signals): halves sort samples
  # differently -> coherence drops to 0.5 or lower.
  #
  # This is the right diagnostic when the issue is "wrong rectangle
  # geometry" rather than "family LD inflating bands". It fires when the
  # node averages across multiple distinct sortings.
  coherence <- NA_real_
  if (length(rows) >= 6L) {
    half_size <- floor(length(rows) / 2L)
    rows_a <- rows[1:half_size]
    rows_b <- rows[(length(rows) - half_size + 1L):length(rows)]
    pc1_a <- colMeans(pc1_mat[rows_a, , drop = FALSE], na.rm = TRUE)
    pc1_b <- colMeans(pc1_mat[rows_b, , drop = FALSE], na.rm = TRUE)
    valid_ab <- is.finite(pc1_a) & is.finite(pc1_b)
    if (sum(valid_ab) >= 5L &&
        sd(pc1_a[valid_ab]) > 0 && sd(pc1_b[valid_ab]) > 0) {
      coherence <- cor(pc1_a[valid_ab], pc1_b[valid_ab])
    }
  }
  coh_lab <- if (is.finite(coherence)) sprintf("coh=%.2f", coherence) else "coh=NA"

  # Family purity per node: for each cluster, compute the fraction of
  # samples that come from the most-common family in that cluster. Then
  # average those fractions weighted by cluster size.
  #
  # Reading guide:
  #   purity ~ 1.0  -> each cluster IS one family ==> family LD, NOT inversion
  #   purity ~ 0.2  -> clusters mix many families  ==> inversion-like signal
  #   purity NA     -> no family info loaded
  fam_purity <- NA_real_
  if (have_family && ak$best_k >= 2L) {
    fams <- sample_meta$family_id[valid]
    cls_ids <- ak$cluster
    purity_vec <- numeric(0)
    size_vec   <- numeric(0)
    for (k in seq_len(ak$best_k)) {
      in_k <- which(cls_ids == k)
      if (length(in_k) == 0) next
      fam_in_k <- fams[in_k]
      # Ignore unmatched (-1) when computing purity within cluster
      fam_in_k <- fam_in_k[fam_in_k != -1L]
      if (length(fam_in_k) == 0) next
      tab <- table(fam_in_k)
      max_count <- max(tab)
      purity_vec <- c(purity_vec, max_count / length(fam_in_k))
      size_vec   <- c(size_vec, length(fam_in_k))
    }
    if (length(purity_vec) > 0) {
      fam_purity <- weighted.mean(purity_vec, w = size_vec)
    }
  }
  fam_lab <- if (is.finite(fam_purity)) sprintf("fam_purity=%.2f", fam_purity) else ""

  # COMPOSITE? badge: best_k > 4 OR peaks > 3 = candidate for over-merge
  is_composite <- (!is.na(ak$best_k) && ak$best_k > 4L) ||
                   (!is.na(n_peaks) && n_peaks > 3L)
  comp_lab <- if (is_composite) " [COMPOSITE?]" else ""

  # FAMILY-LD? badge: high purity = bands ARE families = NOT inversion
  is_familyld <- is.finite(fam_purity) && fam_purity >= 0.70
  fld_lab <- if (is_familyld) " [FAM-LD?]" else ""

  # INCOHERENT? badge: halves of the node disagree on sample sorting.
  # The rectangle the tree drew on the heatmap is enclosing more than one
  # block. This is the "wrong square" diagnosis: the node should be split.
  # Threshold 0.6 chosen so that genuine biological noise inside one
  # inversion still passes (correlations ~0.85+ within real blocks) while
  # composite cases (correlations ~0.3-0.5) get flagged.
  is_incoherent <- is.finite(coherence) && coherence < 0.60
  inc_lab <- if (is_incoherent) " [INCOHERENT]" else ""

  facet_label <- sprintf(
    "n%s | %s | %.2f-%.2f Mb | %dw | %s | %s | %s | %s | %s%s%s%s%s",
    as.character(nid), cls,
    tree$start_mb[ni], tree$end_mb[ni], nwin,
    k_lab, sil_lab, pk_lab, bc_lab, coh_lab,
    if (nzchar(fam_lab)) paste0(" | ", fam_lab) else "",
    comp_lab, fld_lab, inc_lab
  )

  panels[[ni]] <- data.table(
    sample      = sample_ids,
    PC1         = avg_pc1,
    PC2         = avg_pc2,
    band        = band,
    facet_label = facet_label,
    classification = cls,
    node_id     = nid,
    facet_order = ni,        # preserves chromosome-position order
    nwin_node   = nwin,
    best_k      = ak$best_k,
    best_sil    = ak$best_sil,
    n_peaks     = n_peaks,
    bc          = bc,
    coherence   = coherence,
    composite   = is_composite,
    fam_purity  = fam_purity,
    family_ld   = is_familyld,
    incoherent  = is_incoherent
  )
}

# Drop NULLs from skipped nodes
panels <- panels[!vapply(panels, is.null, logical(1))]
if (length(panels) == 0) {
  cat("[D16] no panels built — exiting\n")
  quit(status = 0)
}

panel_all <- rbindlist(panels, fill = TRUE)

# Merge ancestry, cga, family_id from sample_meta
panel_all <- merge(panel_all,
                   sample_meta[, .(ind, cga, ancestry, family_id)],
                   by.x = "sample", by.y = "ind", all.x = TRUE, sort = FALSE)
panel_all[is.na(ancestry), ancestry := "unknown"]
panel_all[is.na(family_id), family_id := -1L]
panel_all[is.na(cga), cga := sample]

# Family display label: prefix with "F" so it's clearly categorical, and
# tag unmatched as "F?"
panel_all[, family_label := ifelse(family_id == -1L, "F?",
                                    paste0("F", family_id))]

# ---- Pagination -------------------------------------------------------------

facet_ord <- unique(panel_all[, .(facet_label, facet_order)])
setorder(facet_ord, facet_order)
facet_ord[, page := rep(seq_len(ceiling(.N / panels_per_page)),
                        each = panels_per_page)[seq_len(.N)]]

panel_all <- merge(panel_all,
                   facet_ord[, .(facet_label, page)],
                   by = "facet_label", sort = FALSE)

# Preserve facet ordering in ggplot
panel_all[, facet_label := factor(facet_label,
                                  levels = facet_ord$facet_label)]

n_pages <- max(facet_ord$page)
cat("[D16] panels: ", nrow(facet_ord), " | pages: ", n_pages,
    " (", panels_per_page, " per page)\n", sep = "")

# ---- Plot ------------------------------------------------------------------
# Two color schemes per page:
#   - by k=3 band  (the eye-test)
#   - by ancestry  (only meaningful if --samples was passed)
# We always emit the band version. The ancestry version is emitted only
# if at least one sample has a non-"unknown" ancestry.

emit_ancestry <- !all(panel_all$ancestry == "unknown")

theme_facet <- function() {
  theme_minimal(base_size = 7) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#f8fafc", color = NA),
      panel.grid       = element_blank(),
      panel.border     = element_rect(fill = NA, color = "#cbd5e1", linewidth = 0.2),
      strip.background = element_rect(fill = "#0e1116", color = NA),
      strip.text       = element_text(size = 5.5, color = "#e7edf3", margin = margin(2, 2, 2, 2)),
      legend.position  = "bottom",
      legend.key.size  = unit(0.3, "cm"),
      legend.title     = element_text(size = 7),
      legend.text      = element_text(size = 6),
      axis.text        = element_text(size = 4, color = "grey50"),
      axis.title       = element_text(size = 7),
      panel.spacing    = unit(0.15, "lines"),
      plot.title       = element_text(size = 10, face = "bold"),
      plot.subtitle    = element_text(size = 7, color = "grey40")
    )
}

render_one_page <- function(page_data, pg, n_pages, color_var, color_pal,
                            color_label) {
  scale_arg <- if (free_scales) "free_x" else "fixed"

  # The rug sits below the density (negative y). Density y is on a positive
  # scale; we don't show its axis (it's a relative value for shape only).
  # Rug y is in [-0.20, -0.04], density in [0, +max_density_y].
  rug_y_top <- -0.04
  rug_y_bot <- -0.20

  # Pre-compute rug y position per row using a stable per-sample jitter so
  # each sample sits at the same y across facets — easier to track.
  set.seed(7L)
  per_sample_yjit <- runif(length(unique(page_data$sample)),
                            min = rug_y_bot, max = rug_y_top)
  yjit_map <- setNames(per_sample_yjit, unique(page_data$sample))
  page_data <- copy(page_data)
  page_data[, rug_y := yjit_map[sample]]

  ggplot(page_data, aes(x = PC1)) +
    # 1D density — the eye reads peak count from this curve. Neutral
    # grey fill so the curve doesn't bias toward any clustering interpretation.
    geom_density(aes(y = after_stat(scaled)),
                 fill = "#94a3b8", color = "#475569",
                 alpha = 0.45, linewidth = 0.3,
                 bw = "nrd0", n = 256L) +
    # Sample rug below the density. Color = kmeans band (orientation only).
    geom_point(aes(y = rug_y, color = .data[[color_var]]),
               size = 0.45, alpha = 0.85) +
    scale_color_manual(values = color_pal, na.value = "grey60",
                        name = color_label,
                        guide = guide_legend(nrow = 2, byrow = TRUE,
                                              override.aes = list(size = 2.5))) +
    facet_wrap(~ facet_label, ncol = ncol_f, scales = scale_arg) +
    theme_facet() +
    labs(
      title    = sprintf("%s | tree-node PC1 density facets | %s | page %d/%d",
                        chr_label, color_label, pg, n_pages),
      subtitle = sprintf("filter=%s | min_width=%dw | %d nodes total | facets=%dx%d | %s | density curve = shape, rug = samples",
                        tree_filter, min_width,
                        nrow(facet_ord), ncol_f, nrow_f,
                        if (free_scales) "free x scales" else "fixed x scale"),
      x = "Local PC1 (averaged across node windows)",
      y = NULL
    ) +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      legend.position = "bottom",
      legend.box      = "horizontal",
      legend.direction = "horizontal",
      legend.margin   = margin(0, 0, 0, 0),
      legend.spacing.x = unit(0.15, "cm")
    )
}

write_pdf <- function(color_var, color_pal, color_label, tag) {
  pdf_path <- file.path(outdir,
    sprintf("%s_tree_pca_facets_%s_%s.pdf", chr_label, tree_filter, tag))
  cat("[D16] writing ", pdf_path, "\n", sep = "")
  pdf(pdf_path, width = page_w, height = page_h)
  on.exit(dev.off(), add = TRUE)
  for (pg in seq_len(n_pages)) {
    pdat <- panel_all[page == pg]
    p <- render_one_page(pdat, pg, n_pages, color_var, color_pal,
                          color_label = color_label)
    print(p)
  }
  dev.off()
  on.exit()
  cat("[D16]   done: ", pdf_path, "\n", sep = "")
}

write_pdf("band", BAND_COLOR, "adaptive K cluster on PC1", "bands")

if (emit_ancestry) {
  # Build ancestry palette dynamically from observed values
  anc_levels <- sort(unique(panel_all$ancestry))
  # Use a simple categorical palette
  base_pal <- c("#7f1d1d", "#b45309", "#0d9488", "#4338ca", "#be185d",
                "#7c2d12", "#15803d", "#9a3412", "#1e40af", "#831843")
  if (length(anc_levels) > length(base_pal)) {
    base_pal <- rep(base_pal, length.out = length(anc_levels))
  }
  anc_pal <- setNames(base_pal[seq_along(anc_levels)], anc_levels)
  anc_pal["unknown"] <- "#94a3b8"

  write_pdf("ancestry", anc_pal, "ancestry / coarse group", "ancestry")
}

# ---- Family-coloring PDF (only if --family was given) -----------------------
# Each sample is colored by its family ID. A real inversion will mix many
# family colors within each K-band on the rug. A family-LD region will
# show coherent family colors per K-band (the family-purity in the title
# also tells you this numerically).
#
# Color strategy: large hub families get distinct colors; small families
# (singletons or pairs) all get a neutral grey "small_family" color so
# they do not visually dominate. This matches the relatedness network
# figure logic.
if (have_family) {
  # Count DISTINCT SAMPLES per family — not panel rows. panel_all has 1 row
  # per (sample, node), so naive table(panel_all$family_label) inflates every
  # family by n_nodes. The earlier version had this bug and put every family
  # in the legend regardless of size.
  fam_size <- panel_all[, .(n_samp = uniqueN(sample)), by = family_label]
  hub_threshold <- 4L
  hub_families <- fam_size[n_samp >= hub_threshold & family_label != "F?",
                            family_label]
  small_families <- setdiff(unique(panel_all$family_label),
                             c(hub_families, "F?"))
  cat("[D16] family sizes: ",
      "n_singletons=", sum(fam_size$n_samp == 1L), " ",
      "n_pairs=", sum(fam_size$n_samp == 2L), " ",
      "n_triplets=", sum(fam_size$n_samp == 3L), " ",
      "n_hubs(>=", hub_threshold, ")=", length(hub_families),
      "\n", sep = "")

  # Big distinct categorical palette (Tol qualitative + extras)
  fam_pal_base <- c(
    "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77",
    "#CC6677", "#AA4499", "#882255", "#999933", "#661100",
    "#6699CC", "#888888", "#ed7953", "#0d9488", "#7c3aed",
    "#dc2626", "#16a34a", "#ca8a04", "#0284c7", "#9333ea",
    "#be123c", "#a16207", "#15803d", "#1d4ed8", "#7e22ce",
    "#9f1239", "#854d0e", "#166534", "#1e40af", "#6b21a8"
  )
  if (length(hub_families) > length(fam_pal_base)) {
    fam_pal_base <- rep(fam_pal_base, length.out = length(hub_families))
  }
  fam_pal <- setNames(fam_pal_base[seq_along(hub_families)], hub_families)
  fam_pal["small_fam"] <- "#cbd5e1"
  fam_pal["F?"]        <- "#94a3b8"

  # Collapse all small families to a single label "small_fam" so the legend
  # only lists the hub families plus one bucket. Without this the legend
  # explodes to hundreds of entries (~210 singletons in the first-degree
  # natora file) and dominates the page.
  panel_all[, family_label_collapsed := ifelse(
    family_label %in% hub_families, family_label,
    ifelse(family_label == "F?", "F?", "small_fam"))]

  n_small <- length(small_families)
  n_singletons_in_data <- sum(fam_size$n_samp == 1L)
  cat("[D16] family-color PDF: ", length(hub_families), " hub families (>= ",
      hub_threshold, " samples), ",
      n_small, " small families collapsed to 'small_fam'\n", sep = "")

  write_pdf("family_label_collapsed", fam_pal,
            paste0("family (n_hub=", length(hub_families),
                    ", small=", n_small, ")"),
            "family")
}

# ---- Over-merge diagnostic --------------------------------------------------
# One row per node: nwin_node vs best_k, colored by classification, sized
# by silhouette score, and marked if COMPOSITE? was flagged. If wide nodes
# systematically have higher best_k than narrow ones (positive trend on
# the scatter), that's the over-merge signature: D09 merged independent
# structural signals into one node, and averaging their PC1s collapsed
# the signal into a higher-K mixture.

node_summary <- unique(panel_all[, .(node_id, classification, nwin_node,
                                      best_k, best_sil, n_peaks, bc,
                                      coherence, composite,
                                      fam_purity, family_ld, incoherent)])

# Write a TSV alongside the PDFs for downstream inspection
summary_tsv <- file.path(outdir,
  sprintf("%s_tree_pca_facets_%s_node_summary.tsv", chr_label, tree_filter))
fwrite(node_summary, summary_tsv, sep = "\t")
cat("[D16] node summary written: ", summary_tsv, "\n", sep = "")

# Spearman correlation between nwin_node and best_k as a one-number
# over-merge indicator. rho > 0 with low p = positive trend = over-merge.
ns_for_cor <- node_summary[is.finite(best_k) & is.finite(nwin_node)]
if (nrow(ns_for_cor) >= 5) {
  ct <- suppressWarnings(
    cor.test(ns_for_cor$nwin_node, ns_for_cor$best_k, method = "spearman")
  )
  rho_lab <- sprintf("Spearman rho(nW, best_k) = %.3f  (p = %.3g, n = %d)",
                     ct$estimate, ct$p.value, nrow(ns_for_cor))
} else {
  rho_lab <- sprintf("Spearman rho not computed (only %d nodes)", nrow(ns_for_cor))
}
cat("[D16] ", rho_lab, "\n", sep = "")

# Composite count
n_composite <- sum(node_summary$composite, na.rm = TRUE)
cat("[D16] COMPOSITE? flagged: ", n_composite, "/", nrow(node_summary),
    " nodes\n", sep = "")

# Incoherent count
n_incoherent <- sum(node_summary$incoherent, na.rm = TRUE)
cat("[D16] INCOHERENT  flagged: ", n_incoherent, "/", nrow(node_summary),
    " nodes (halves coherence < 0.60)\n", sep = "")

# Family-LD count (if family info loaded)
if (any(!is.na(node_summary$family_ld))) {
  n_familyld <- sum(node_summary$family_ld, na.rm = TRUE)
  cat("[D16] FAM-LD?    flagged: ", n_familyld, "/", nrow(node_summary),
      " nodes (fam_purity >= 0.70)\n", sep = "")
}

summary_pdf <- file.path(outdir,
  sprintf("%s_tree_pca_facets_%s_summary.pdf", chr_label, tree_filter))
pdf(summary_pdf, width = 14, height = 9)

p_main <- ggplot(node_summary,
                 aes(x = nwin_node, y = best_k,
                     color = classification, size = best_sil,
                     shape = composite)) +
  geom_jitter(width = 0.0, height = 0.15, alpha = 0.85) +
  scale_color_manual(values = CLASS_COLOR, na.value = "grey60") +
  scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 16),
                     name = "COMPOSITE?",
                     labels = c(`TRUE` = "yes", `FALSE` = "no")) +
  scale_size_continuous(range = c(1.2, 4.5),
                        name = "silhouette",
                        breaks = c(0.1, 0.3, 0.5, 0.7)) +
  scale_y_continuous(breaks = seq(1, k_max, 1),
                     limits = c(0.7, k_max + 0.5)) +
  scale_x_log10() +
  labs(
    title    = sprintf("%s | Over-merge diagnostic | best_k vs nW per tree node",
                       chr_label),
    subtitle = sprintf("filter=%s  |  %s  |  COMPOSITE? n=%d/%d  (best_k>4 OR peaks>3)",
                       tree_filter, rho_lab, n_composite, nrow(node_summary)),
    x = "n windows in node (log10)",
    y = "best_k (silhouette-selected)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#f8fafc", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e2e8f0"),
    legend.position  = "right",
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 10, color = "grey30")
  )
print(p_main)

# Second panel: best_k distribution per classification
p_box <- ggplot(node_summary,
                aes(x = classification, y = best_k, fill = classification)) +
  geom_boxplot(alpha = 0.6, outlier.size = 1) +
  scale_fill_manual(values = CLASS_COLOR, guide = "none") +
  scale_y_continuous(breaks = seq(1, k_max, 1)) +
  labs(
    title    = sprintf("%s | best_k distribution by tree classification", chr_label),
    subtitle = "If INVERSION nodes cluster at best_k=3 and CANDIDATE nodes at higher k, the tree is splitting things correctly.",
    x = NULL, y = "best_k"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#f8fafc", color = NA),
    panel.grid.minor = element_blank()
  )
print(p_box)

# Third panel: fam_purity vs best_k — directly answers "are bands families?"
# Only emitted if family info was loaded.
ns_with_purity <- node_summary[is.finite(fam_purity)]
if (nrow(ns_with_purity) >= 5) {
  p_purity <- ggplot(ns_with_purity,
                     aes(x = best_k, y = fam_purity,
                         color = classification, size = nwin_node,
                         shape = family_ld)) +
    geom_jitter(width = 0.18, height = 0.0, alpha = 0.85) +
    geom_hline(yintercept = 0.70, linetype = "dashed",
               color = "#dc2626", linewidth = 0.4) +
    annotate("text", x = max(ns_with_purity$best_k), y = 0.72,
             label = "FAM-LD threshold (0.70)", hjust = 1, size = 3,
             color = "#dc2626") +
    scale_color_manual(values = CLASS_COLOR, na.value = "grey60") +
    scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 16),
                       name = "FAM-LD?",
                       labels = c(`TRUE` = "yes", `FALSE` = "no")) +
    scale_size_continuous(range = c(1.2, 4.5),
                          trans = "log10",
                          name = "n windows") +
    scale_x_continuous(breaks = seq(1, k_max, 1)) +
    scale_y_continuous(limits = c(0, 1.05),
                       breaks = seq(0, 1, 0.2)) +
    labs(
      title    = sprintf("%s | Family LD diagnostic | fam_purity vs best_k",
                         chr_label),
      subtitle = sprintf("Nodes above 0.70 are likely family LD, not biological inversions. n_FAM-LD=%d/%d",
                         sum(node_summary$family_ld, na.rm = TRUE),
                         nrow(node_summary)),
      x = "best_k", y = "family purity (size-weighted)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#f8fafc", color = NA),
      panel.grid.minor = element_blank()
    )
  print(p_purity)
}

# Fourth panel: coherence vs nwin_node — the "wrong square" diagnostic.
# This is the most direct test of the over-merge hypothesis you raised:
# wide nodes that enclose multiple distinct sortings of the samples will
# show low halves coherence. Real inversion blocks (even very wide ones)
# should hold coherence near 1.0 because every window inside them sorts
# samples the same way.
ns_with_coh <- node_summary[is.finite(coherence)]
if (nrow(ns_with_coh) >= 5) {
  # Spearman rho between nwin and coherence — should be NEGATIVE for over-merge
  ct_coh <- suppressWarnings(
    cor.test(ns_with_coh$nwin_node, ns_with_coh$coherence, method = "spearman")
  )
  rho_coh_lab <- sprintf("Spearman rho(nW, coherence) = %.3f  (p = %.3g, n = %d)",
                          ct_coh$estimate, ct_coh$p.value, nrow(ns_with_coh))

  p_coh <- ggplot(ns_with_coh,
                  aes(x = nwin_node, y = coherence,
                      color = classification, size = best_sil,
                      shape = incoherent)) +
    geom_jitter(width = 0.0, height = 0.0, alpha = 0.85) +
    geom_hline(yintercept = 0.60, linetype = "dashed",
               color = "#dc2626", linewidth = 0.4) +
    annotate("text", x = max(ns_with_coh$nwin_node), y = 0.62,
             label = "INCOHERENT threshold (0.60)", hjust = 1, size = 3,
             color = "#dc2626") +
    scale_color_manual(values = CLASS_COLOR, na.value = "grey60") +
    scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = 16),
                       name = "INCOHERENT?",
                       labels = c(`TRUE` = "yes", `FALSE` = "no")) +
    scale_size_continuous(range = c(1.2, 4.5), name = "silhouette",
                          breaks = c(0.1, 0.3, 0.5, 0.7)) +
    scale_x_log10() +
    scale_y_continuous(limits = c(NA, 1.05),
                       breaks = seq(-1, 1, 0.25)) +
    labs(
      title    = sprintf("%s | Wrong-square diagnostic | halves coherence vs nW",
                        chr_label),
      subtitle = sprintf("Low coherence -> the rectangle encloses >1 block. n_INCOHERENT=%d/%d  |  %s",
                         n_incoherent, nrow(node_summary), rho_coh_lab),
      x = "n windows in node (log10)",
      y = "halves coherence (Pearson r of avg PC1 between halves)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#f8fafc", color = NA),
      panel.grid.minor = element_blank()
    )
  print(p_coh)
}

dev.off()
cat("[D16] over-merge summary written: ", summary_pdf, "\n", sep = "")

cat("[D16] done. v3.2 outputs:\n")
cat("       _bands.pdf       facets colored by adaptive K cluster on PC1\n")
cat("       _family.pdf      facets colored by family ID (if --family given)\n")
cat("       _ancestry.pdf    facets colored by ancestry (if --samples given)\n")
cat("       _summary.pdf     four panels: best_k vs nW, best_k boxplot,\n")
cat("                        fam_purity vs best_k, coherence vs nW\n")
cat("       _node_summary.tsv per-node table with all metrics\n")
cat("\n")
cat("       Reading the facets (in priority order):\n")
cat("       coh<0.60                       -> [INCOHERENT] wrong square,\n")
cat("                                          rectangle drawn too wide,\n")
cat("                                          node should be split\n")
cat("       coh>=0.85, best_k=3, sil>0.5   -> classic clean inversion\n")
cat("       coh>=0.85, best_k=2, fam>0.7   -> family LD on a coherent block\n")
cat("       coh>=0.85, best_k>=5 [COMP?]   -> COMPOSITE-real (multiple\n")
cat("                                          inversion alleles at same locus)\n")
cat("       best_k=1                       -> unimodal, no structure\n")
