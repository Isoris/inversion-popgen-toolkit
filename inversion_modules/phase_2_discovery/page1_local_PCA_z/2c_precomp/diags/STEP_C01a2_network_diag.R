#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01a2_network_diag.R  (v8.4)
#
# SAMPLE NETWORK DIAGNOSTICS from precomputed sim_mat.
#
# For each chromosome, identifies blocks of highly similar windows
# (k-NN clusters), then within each block extracts the SAMPLE-level
# PC1 scores and builds a sample similarity network. This shows which
# samples cluster together in each inversion candidate region.
#
# Helps distinguish:
#   - Real inversions: samples from MULTIPLE families cluster together
#   - Family artifacts: one family dominates a cluster
#
# Inputs:
#   <precomp_dir>/*.precomp.rds  -- from C01a
#   [optional] <ancestry_file>   -- NGSadmix Q-matrix (best_seed K=8)
#   [optional] <relatedness_file> -- ngsRelate pairwise theta
#
# Outputs per chromosome:
#   01_window_cluster_map.pdf    -- which windows belong to which NN cluster
#   02_sample_network_facet.pdf  -- one facet per cluster: sample network
#   03_sample_network_colored.pdf -- same but nodes colored by cluster membership
#   04_ancestry_check.pdf        -- if Q-matrix provided: ancestry diversity per cluster
#
# Usage:
#   Rscript STEP_C01a2_network_diag.R <precomp_dir> <outdir> [--ancestry <qopt>] \
#     [--relatedness <pairs.tsv>] [--samples <list.txt>] [--chrom <chr>]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript STEP_C01a2_network_diag.R <precomp_dir> <outdir> ",
       "[--ancestry <qopt>] [--relatedness <pairs.tsv>] [--samples <list>] [--chrom <chr>]")
}

precomp_dir   <- args[1]
outdir        <- args[2]
ancestry_file <- NULL
relate_file   <- NULL
samples_file  <- NULL
chrom_filter  <- NULL

i <- 3L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--ancestry" && i < length(args))    { ancestry_file <- args[i+1]; i <- i + 2L }
  else if (a == "--relatedness" && i < length(args)) { relate_file <- args[i+1]; i <- i + 2L }
  else if (a == "--samples" && i < length(args)) { samples_file <- args[i+1]; i <- i + 2L }
  else if (a == "--chrom" && i < length(args))   { chrom_filter <- args[i+1]; i <- i + 2L }
  else { i <- i + 1L }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DPI <- 350
THEME_BASE <- theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#555555"),
        plot.caption = element_text(size = 6, color = "#888888", hjust = 0),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))


# Cluster palette (named R colors for HPC compatibility -- no hex codes)
CLUSTER_PAL <- c(
  "red3", "royalblue", "green4", "darkorange", "purple3",
  "cyan4", "deeppink3", "dodgerblue", "olivedrab", "magenta3",
  "turquoise4", "goldenrod", "slateblue", "forestgreen", "violetred",
  "mediumpurple", "chocolate", "cadetblue", "darkviolet", "sienna",
  "steelblue", "tomato", "seagreen", "plum3", "tan3",
  "palevioletred", "darkslategray", "coral", "aquamarine4", "rosybrown"
)

# =============================================================================
# LOAD DATA
# =============================================================================

message("[C01a2] Loading precomputed data...")
rds_files <- sort(list.files(precomp_dir, pattern = "\\.precomp\\.rds$", full.names = TRUE))
if (length(rds_files) == 0) stop("No .precomp.rds files in: ", precomp_dir)

precomp_list <- list()
for (f in rds_files) {
  obj <- readRDS(f)
  precomp_list[[obj$chrom]] <- obj
}
chroms <- names(precomp_list)
if (!is.null(chrom_filter) && chrom_filter != "all") chroms <- intersect(chroms, chrom_filter)
message("[C01a2] ", length(chroms), " chromosomes")

# Sample names from precomp (may be Ind0, Ind1... from ANGSD)
sample_names <- NULL
for (chr_tmp in chroms) {
  pc1_cols <- grep("^PC_1_", names(precomp_list[[chr_tmp]]$dt), value = TRUE)
  if (length(pc1_cols) > 0) { sample_names <- sub("^PC_1_", "", pc1_cols); break }
}
if (is.null(sample_names)) stop("No PC_1_ columns found in precomp data")
n_samples <- length(sample_names)
message("[C01a2] Precomp sample names: ", sample_names[1], " ... ", sample_names[n_samples],
        " (", n_samples, " samples)")

# Build mapping from precomp names (Ind0...) to real sample names (CGA009...)
# If samples_file is provided AND precomp uses Ind0/Ind1 naming, create the map
real_sample_names <- sample_names  # default: same as precomp
name_map <- NULL  # precomp_name -> real_name
if (!is.null(samples_file) && file.exists(samples_file)) {
  slist <- as.character(fread(samples_file, header = FALSE)[[1]])
  if (length(slist) == n_samples && grepl("^Ind[0-9]", sample_names[1])) {
    # Precomp uses Ind0, Ind1... -> map to real names from sample list
    name_map <- setNames(slist, sample_names)
    real_sample_names <- slist
    message("[C01a2] Sample name mapping: ", sample_names[1], " -> ", slist[1],
            " ... ", sample_names[n_samples], " -> ", slist[n_samples])
  } else if (length(slist) == n_samples) {
    real_sample_names <- slist
  }
}
message("[C01a2] Real sample names: ", real_sample_names[1], " ... ",
        real_sample_names[n_samples])

# Optional: ancestry Q-matrix
q_mat <- NULL
if (!is.null(ancestry_file) && file.exists(ancestry_file)) {
  q_raw <- as.matrix(fread(ancestry_file, header = FALSE))
  if (nrow(q_raw) == n_samples) {
    rownames(q_raw) <- real_sample_names
    colnames(q_raw) <- paste0("K", seq_len(ncol(q_raw)))
    q_mat <- q_raw
    message("[C01a2] Ancestry Q-matrix loaded: ", nrow(q_mat), " samples x K=", ncol(q_mat))
  } else {
    message("[C01a2] Q-matrix rows (", nrow(q_raw), ") != samples (", n_samples, ") -- skipping")
  }
}

# Optional: relatedness
relate_dt <- NULL
if (!is.null(relate_file) && file.exists(relate_file)) {
  # Try multiple formats:
  # 1. NATOra format: sample1 sample2 theta (3 cols, no header)
  # 2. ngsRelate .res format: many columns with header including "a" "b" "rab"
  raw_dt <- tryCatch(fread(relate_file), error = function(e) NULL)
  if (!is.null(raw_dt) && nrow(raw_dt) > 0) {
    if (all(c("sample1", "sample2", "theta") %in% names(raw_dt))) {
      # Already in expected format
      relate_dt <- raw_dt[, .(sample1, sample2, theta)]
    } else if (ncol(raw_dt) == 3 && !any(grepl("[a-zA-Z]", raw_dt[[3]][1:min(5, nrow(raw_dt))]))) {
      # 3 columns, third is numeric -> NATOra format
      setnames(raw_dt, c("sample1", "sample2", "theta"))
      raw_dt[, theta := as.numeric(theta)]
      relate_dt <- raw_dt
    } else if ("rab" %in% names(raw_dt)) {
      # ngsRelate .res format with "a" "b" "rab" columns (indices)
      # Need sample list to convert indices to names
      if (!is.null(samples_file) && file.exists(samples_file)) {
        slist <- as.character(fread(samples_file, header = FALSE)[[1]])
        raw_dt[, sample1 := slist[a + 1L]]  # ngsRelate uses 0-based
        raw_dt[, sample2 := slist[b + 1L]]
        raw_dt[, theta := rab]
        relate_dt <- raw_dt[, .(sample1, sample2, theta)]
      }
    } else if (ncol(raw_dt) >= 3) {
      # Fallback: assume first 3 columns
      setnames(raw_dt, 1:3, c("sample1", "sample2", "theta"))
      raw_dt[, theta := as.numeric(theta)]
      relate_dt <- raw_dt[, .(sample1, sample2, theta)]
    }
    if (!is.null(relate_dt)) {
      relate_dt <- relate_dt[is.finite(theta)]
      # Match against real sample names (CGA009 etc.)
      relate_dt <- relate_dt[sample1 %in% real_sample_names & sample2 %in% real_sample_names]
      message("[C01a2] Relatedness pairs loaded: ", nrow(relate_dt))
      if (nrow(relate_dt) == 0) relate_dt <- NULL
    }
  }
}

# =============================================================================
# HELPER: Find window clusters from sim_mat using NN connectivity
# =============================================================================

find_window_clusters <- function(sim_mat, dt, min_cluster_size = 5L, 
                                  nn_k = 5L, sim_threshold_quantile = 0.75) {
  n_w <- nrow(sim_mat)
  if (n_w < 20) return(NULL)
  
  # For each window, find k most similar windows
  # Build adjacency: connect windows whose mutual similarity > threshold
  nn_sims <- numeric(n_w)
  for (wi in seq_len(n_w)) {
    sims <- sim_mat[wi, ]
    sims[wi] <- -Inf
    top_k <- order(sims, decreasing = TRUE)[seq_len(nn_k)]
    nn_sims[wi] <- mean(sims[top_k])
  }
  
  # Threshold: windows with high mean k-NN similarity
  sim_thresh <- quantile(nn_sims, sim_threshold_quantile, na.rm = TRUE)
  
  # Connected components of high-similarity windows
  # Two windows are connected if they are mutual k-NN AND both above threshold
  adj <- matrix(FALSE, n_w, n_w)
  for (wi in seq_len(n_w)) {
    if (nn_sims[wi] < sim_thresh) next
    sims <- sim_mat[wi, ]
    sims[wi] <- -Inf
    top_k <- order(sims, decreasing = TRUE)[seq_len(nn_k)]
    for (ki in top_k) {
      if (nn_sims[ki] >= sim_thresh && abs(wi - ki) <= 100) {
        # Only connect if within 100 windows (local clustering, not long-range)
        adj[wi, ki] <- TRUE
        adj[ki, wi] <- TRUE
      }
    }
  }
  
  # Find connected components
  visited <- rep(FALSE, n_w)
  cluster_id <- rep(0L, n_w)
  cid <- 0L
  
  for (wi in seq_len(n_w)) {
    if (visited[wi] || !any(adj[wi, ])) next
    cid <- cid + 1L
    # BFS
    queue <- wi
    while (length(queue) > 0) {
      curr <- queue[1]; queue <- queue[-1]
      if (visited[curr]) next
      visited[curr] <- TRUE
      cluster_id[curr] <- cid
      neighbors <- which(adj[curr, ] & !visited)
      queue <- c(queue, neighbors)
    }
  }
  
  # Filter small clusters
  cluster_tab <- table(cluster_id[cluster_id > 0])
  valid_clusters <- as.integer(names(cluster_tab[cluster_tab >= min_cluster_size]))
  cluster_id[!cluster_id %in% valid_clusters] <- 0L
  
  # Renumber
  if (length(valid_clusters) > 0) {
    new_ids <- setNames(seq_along(valid_clusters), as.character(valid_clusters))
    cluster_id[cluster_id > 0] <- new_ids[as.character(cluster_id[cluster_id > 0])]
  }
  
  data.table(
    window_idx = seq_len(n_w),
    pos_mb = (dt$start_bp + dt$end_bp) / 2e6,
    start_bp = dt$start_bp,
    end_bp = dt$end_bp,
    cluster_id = cluster_id,
    mean_knn_sim = nn_sims
  )
}

# =============================================================================
# HELPER: Build sample network within a window cluster
# =============================================================================

build_sample_network <- function(dt, window_idxs, sample_names, name_map = NULL) {
  pc1_cols <- paste0("PC_1_", sample_names)
  available <- intersect(pc1_cols, names(dt))
  if (length(available) < 20) return(NULL)
  
  mat <- as.matrix(dt[window_idxs, ..available])
  if (nrow(mat) > 1) {
    avg_loadings <- colMeans(mat, na.rm = TRUE)
  } else {
    avg_loadings <- mat[1, ]
  }
  
  valid <- is.finite(avg_loadings)
  if (sum(valid) < 20) return(NULL)
  
  vals <- avg_loadings[valid]
  snames <- sub("^PC_1_", "", names(vals))
  
  # Map Ind0 -> CGA009 etc. if mapping exists
  if (!is.null(name_map)) {
    mapped <- name_map[snames]
    snames[!is.na(mapped)] <- mapped[!is.na(mapped)]
  }
  
  # k=3 clustering
  km <- tryCatch(kmeans(vals, centers = 3, nstart = 5), error = function(e) NULL)
  if (is.null(km)) return(NULL)
  
  center_order <- order(km$centers[, 1])
  band_labels <- character(length(vals))
  for (bi in seq_along(center_order)) {
    band_labels[km$cluster == center_order[bi]] <- paste0("band", bi)
  }
  
  data.table(
    sample = snames,
    pc1_loading = vals,
    band = band_labels
  )
}

# =============================================================================
# MAIN LOOP
# =============================================================================

message("[C01a2] Processing chromosomes...")

for (chr in chroms) {
  pc <- precomp_list[[chr]]
  if (is.null(pc) || pc$n_windows < 20) next
  
  message("\n[C01a2] === ", chr, " ===")
  dt <- pc$dt
  sim_mat <- pc$sim_mat
  n_w <- nrow(sim_mat)
  
  chr_outdir <- file.path(outdir, chr)
  dir.create(chr_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # -- Step 1: Find window clusters ------------------------------------
  wc <- find_window_clusters(sim_mat, dt)
  if (is.null(wc) || max(wc$cluster_id) == 0) {
    message("  No significant clusters found")
    next
  }
  
  n_clusters <- max(wc$cluster_id)
  message("  Found ", n_clusters, " window clusters")
  
  # Cap at 30 largest clusters to avoid palette overflow and giant plots
  if (n_clusters > 30) {
    cluster_sizes <- wc[cluster_id > 0, .N, by = cluster_id][order(-N)]
    keep_ids <- cluster_sizes$cluster_id[seq_len(min(30, nrow(cluster_sizes)))]
    wc[!cluster_id %in% keep_ids, cluster_id := 0L]
    # Renumber
    old_ids <- sort(unique(wc$cluster_id[wc$cluster_id > 0]))
    new_map <- setNames(seq_along(old_ids), as.character(old_ids))
    wc[cluster_id > 0, cluster_id := new_map[as.character(cluster_id)]]
    n_clusters <- max(wc$cluster_id)
    message("  Capped to ", n_clusters, " largest clusters")
  }
  
  # -- Plot 01: Window cluster map -------------------------------------
  # Use continuous color scale (cluster_id as numeric) to avoid RGB issues
  # with large discrete palettes on HPC
  wc[, cluster_color := fifelse(cluster_id > 0, as.numeric(cluster_id), NA_real_)]
  
  p1 <- ggplot(wc, aes(x = pos_mb, y = mean_knn_sim)) +
    geom_point(data = wc[cluster_id == 0], color = "grey80", size = 0.3, alpha = 0.3) +
    geom_point(data = wc[cluster_id > 0], aes(color = cluster_color),
               size = 0.6, alpha = 0.7) +
    scale_color_viridis_c(
                          name = "Cluster ID") +
    labs(x = paste0(chr, " (Mb)"), y = "Mean k-NN similarity",
         title = paste0(chr, " -- Window Clusters (k-NN connectivity)"),
         subtitle = paste0(n_w, " windows, ", n_clusters, " clusters found"),
         caption = "Colored = cluster members | Grey = background") +
    THEME_BASE

  tryCatch({
    use_cairo <- capabilities("cairo")
    f1 <- file.path(chr_outdir, "01_window_cluster_map.png")
    ggsave(f1, p1, width = 14, height = 5, dpi = DPI)
    message("  -> ", f1)
  }, error = function(e) message("  [FAIL] 01: ", e$message))
  
  # -- Plot 02 + 03: Sample networks per cluster (faceted) -------------
  cluster_sample_rows <- list()
  
  for (ci in seq_len(n_clusters)) {
    widxs <- wc[cluster_id == ci]$window_idx
    if (length(widxs) < 3) next
    
    sn <- build_sample_network(dt, widxs, sample_names, name_map)
    if (is.null(sn)) next
    
    cluster_start <- round(min(wc[cluster_id == ci]$pos_mb), 2)
    cluster_end   <- round(max(wc[cluster_id == ci]$pos_mb), 2)
    
    sn[, cluster := paste0("C", ci, " (", cluster_start, "-", cluster_end, " Mb, ",
                            length(widxs), " win)")]
    sn[, cluster_id := ci]
    cluster_sample_rows[[ci]] <- sn
  }
  
  if (length(cluster_sample_rows) == 0) {
    message("  No clusters with enough windows for sample networks")
    next
  }
  
  all_samples <- rbindlist(cluster_sample_rows)
  
  # Plot 02: Sample PC1 distributions per cluster (faceted strip plots)
  p2 <- ggplot(all_samples, aes(x = pc1_loading, fill = band)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("band1" = "steelblue", "band2" = "gold3", "band3" = "firebrick")) +
    facet_wrap(~ cluster, scales = "free", ncol = min(3, n_clusters)) +
    labs(x = "Average PC1 loading", y = "Sample count",
         title = paste0(chr, " -- Sample PC1 Distribution per Cluster"),
         subtitle = paste0(n_clusters, " clusters, ", n_samples, " samples per cluster\n",
                          "3 bands: band1 (blue) = low PC1, band3 (red) = high PC1"),
         caption = "Trimodal = likely inversion (band1=AA, band2=AB, band3=BB)\nUnimodal = no inversion signal in this cluster") +
    THEME_BASE +
    theme(strip.text = element_text(size = 7, face = "bold"))
  
  tryCatch({
    f2 <- file.path(chr_outdir, "02_sample_pc1_per_cluster.png")
    ggsave(f2, p2, width = 14, height = min(30, max(5, ceiling(n_clusters / 3) * 4)), dpi = DPI)
    message("  -> ", f2)
  }, error = function(e) message("  [FAIL] 02: ", e$message))
  
  # Plot 03: Sample dot plot colored by CLUSTER membership
  # Shows which samples appear in which band across clusters
  # If same sample is band1 in C1 and band1 in C2 -> consistent inversion
  
  # Pivot: for each sample, what band are they in for each cluster?
  sample_bands <- dcast(all_samples, sample ~ cluster_id, value.var = "band")
  
  # For the dot plot, show samples sorted by their band in C1
  c1_data <- all_samples[cluster_id == 1]
  if (nrow(c1_data) > 0) {
    c1_order <- c1_data[order(pc1_loading)]$sample
    all_samples[, sample := factor(sample, levels = c1_order)]
  }
  
  p3 <- ggplot(all_samples, aes(x = pc1_loading, y = sample, color = as.numeric(cluster_id))) +
    geom_point(size = 1.0, alpha = 0.7) +
    scale_color_viridis_c(
                          name = "Cluster ID") +
    facet_wrap(~ cluster, scales = "free_x", ncol = min(4, n_clusters)) +
    labs(x = "Average PC1 loading", y = NULL,
         title = paste0(chr, " -- Sample Positions Across Clusters"),
         subtitle = "Same sample should appear at similar relative position if same inversion",
         caption = "Consistent band membership across clusters = same inversion system\nScattered = different inversion systems or noise") +
    THEME_BASE +
    theme(axis.text.y = element_blank(),
          strip.text = element_text(size = 7, face = "bold"))
  
  tryCatch({
    f3 <- file.path(chr_outdir, "03_sample_positions_across_clusters.png")
    ggsave(f3, p3, width = 14, height = min(30, max(5, n_samples / 30)), dpi = DPI)
    message("  -> ", f3)
  }, error = function(e) message("  [FAIL] 03: ", e$message))
  
  # -- Plot 04: Ancestry check (if Q-matrix provided) ------------------
  if (!is.null(q_mat)) {
    # For each cluster, compute ancestry diversity of each band
    anc_rows <- list()
    for (ci in seq_len(n_clusters)) {
      cs <- all_samples[cluster_id == ci]
      for (bi in paste0("band", 1:3)) {
        band_samples <- cs[band == bi]$sample
        band_samples <- intersect(band_samples, rownames(q_mat))
        if (length(band_samples) < 3) next
        
        q_sub <- q_mat[band_samples, , drop = FALSE]
        # Diversity: how many K components are represented?
        # Use effective number of components (1/sum(p^2) where p = mean proportion)
        mean_q <- colMeans(q_sub)
        eff_k <- 1 / sum(mean_q^2)
        # Also: dominant component fraction
        dom_frac <- max(mean_q)
        
        cluster_start <- round(min(wc[cluster_id == ci]$pos_mb), 2)
        cluster_end   <- round(max(wc[cluster_id == ci]$pos_mb), 2)
        
        anc_rows[[length(anc_rows) + 1]] <- data.table(
          cluster = paste0("C", ci, " (", cluster_start, "-", cluster_end, " Mb)"),
          cluster_id = ci,
          band = bi,
          n_samples = length(band_samples),
          effective_K = round(eff_k, 2),
          dominant_fraction = round(dom_frac, 3)
        )
      }
    }
    
    if (length(anc_rows) > 0) {
      anc_dt <- rbindlist(anc_rows)
      
      p4 <- ggplot(anc_dt, aes(x = band, y = effective_K, fill = band)) +
        geom_col(alpha = 0.8) +
        geom_text(aes(label = paste0("n=", n_samples)), vjust = -0.3, size = 2.5) +
        scale_fill_manual(values = c("band1" = "steelblue", "band2" = "gold3", "band3" = "firebrick")) +
        facet_wrap(~ cluster, scales = "free_x", ncol = min(4, n_clusters)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red3") +
        labs(x = "Band", y = "Effective K (ancestry diversity)",
             title = paste0(chr, " -- Ancestry Diversity per Band per Cluster"),
             subtitle = paste0("Effective K > 2 = multiple families (real inversion)\n",
                              "Effective K ~ 1 = one family dominates (artifact suspect)"),
             caption = paste0("Red dashed = K=1 (single ancestry component)\n",
                             "Higher = more diverse ancestry = more likely real inversion\n",
                             "n = number of samples in that band")) +
        THEME_BASE +
        theme(strip.text = element_text(size = 7, face = "bold"))
      
      tryCatch({
        f4 <- file.path(chr_outdir, "04_ancestry_diversity.png")
        ggsave(f4, p4, width = 14, height = min(30, max(5, ceiling(n_clusters / 4) * 4)), dpi = DPI)
        message("  -> ", f4)
      }, error = function(e) message("  [FAIL] 04: ", e$message))
    }
  }
  
  # -- Plot 05: Relatedness within bands (if provided) -----------------
  if (!is.null(relate_dt) && nrow(relate_dt) > 0) {
    rel_rows <- list()
    for (ci in seq_len(n_clusters)) {
      cs <- all_samples[cluster_id == ci]
      for (bi in paste0("band", 1:3)) {
        band_samples <- cs[band == bi]$sample
        if (length(band_samples) < 3) next
        
        # Within-band relatedness
        within_pairs <- relate_dt[sample1 %in% band_samples & sample2 %in% band_samples]
        if (nrow(within_pairs) > 0) {
          rel_rows[[length(rel_rows) + 1]] <- data.table(
            cluster_id = ci, band = bi, comparison = "within-band",
            mean_theta = round(mean(within_pairs$theta, na.rm = TRUE), 4),
            n_pairs = nrow(within_pairs)
          )
        }
        
        # Between-band relatedness (band vs all other bands)
        other_samples <- cs[band != bi]$sample
        between_pairs <- relate_dt[(sample1 %in% band_samples & sample2 %in% other_samples) |
                                    (sample2 %in% band_samples & sample1 %in% other_samples)]
        if (nrow(between_pairs) > 0) {
          rel_rows[[length(rel_rows) + 1]] <- data.table(
            cluster_id = ci, band = bi, comparison = "between-band",
            mean_theta = round(mean(between_pairs$theta, na.rm = TRUE), 4),
            n_pairs = nrow(between_pairs)
          )
        }
      }
    }
    
    if (length(rel_rows) > 0) {
      rel_dt_summary <- rbindlist(rel_rows)
      
      p5 <- ggplot(rel_dt_summary, aes(x = band, y = mean_theta, fill = comparison)) +
        geom_col(position = "dodge", alpha = 0.8) +
        scale_fill_manual(values = c("within-band" = "steelblue", "between-band" = "grey60")) +
        facet_wrap(~ paste0("C", cluster_id), ncol = min(4, n_clusters)) +
        labs(x = "Band", y = "Mean theta (relatedness)",
             title = paste0(chr, " -- Within vs Between Band Relatedness"),
             subtitle = paste0("If within-band >> between-band -> family artifact\n",
                              "If within-band ~ between-band -> real inversion (unrelated fish cluster)"),
             caption = "Blue = within-band pairs, Grey = between-band pairs\nHigh within-band theta = siblings clustering together (suspect)") +
        THEME_BASE +
        theme(strip.text = element_text(size = 7, face = "bold"))
      
      tryCatch({
        f5 <- file.path(chr_outdir, "05_relatedness_within_vs_between.png")
        ggsave(f5, p5, width = 14, height = min(30, max(5, ceiling(n_clusters / 4) * 4)), dpi = DPI)
        message("  -> ", f5)
      }, error = function(e) message("  [FAIL] 05: ", e$message))
    }
  }
  
  # Save cluster table
  fwrite(wc, file.path(chr_outdir, "window_clusters.tsv.gz"), sep = "\t")
  if (nrow(all_samples) > 0) {
    fwrite(all_samples, file.path(chr_outdir, "sample_band_assignments.tsv.gz"), sep = "\t")
  }
}

message("\n[DONE] STEP_C01a2_network_diag complete -> ", outdir)
