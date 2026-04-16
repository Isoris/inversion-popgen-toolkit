#!/usr/bin/env Rscript

# =============================================================================
# STEP23_candidate_ancestry_association.R
#
# Ancestry / Q interpretation layer for inversion candidates.
#
# For each candidate:
#   1. Build clean sample-level ancestry table from existing NGSadmix
#   2. Merge with STEP21 coarse groups and STEP22 subclusters
#   3. Compute per-group ancestry composition, entropy, purity
#   4. Test association: Cramér's V, chi-square, correlation(v, Q)
#   5. Generate interpretation labels
#
# Inputs:
#   - NGSadmix best_seed_by_K / Q tables
#   - STEP21 rotated PCA tables
#   - STEP22 subcluster tables
#
# Outputs per candidate:
#   - candidate_ancestry_association.tsv
#   - candidate_cluster_ancestry_association.tsv
#
# Usage:
#   Rscript STEP23_candidate_ancestry_association.R <config.R> [K=best] [cid=all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) >= 1) args[1] else "config_inversion_followup.R"
K_choice    <- if (length(args) >= 2 && args[2] != "best") as.integer(args[2]) else NA_integer_
cid_filter  <- if (length(args) >= 3 && args[3] != "all") as.integer(args[3]) else NA_integer_

source(config_file)
ensure_dir(FOLLOWUP_DIR)

# ── Build sample ancestry table ──────────────────────────────────────────────
message("[INFO] Building sample ancestry table...")

# Determine best K
if (is.na(K_choice) && file.exists(BEST_SEED_FILE)) {
  best_seed <- fread(BEST_SEED_FILE)
  # Assume columns like K, seed, loglike — pick K with best loglike
  if ("K" %in% names(best_seed)) {
    if ("loglike" %in% names(best_seed)) {
      K_choice <- best_seed[which.max(loglike), K]
    } else {
      K_choice <- max(best_seed$K)
    }
  }
}
if (is.na(K_choice)) K_choice <- 3L
message("[INFO] Using K = ", K_choice)

# Try to load pre-built sample_main_ancestry table
ancestry_dt <- NULL
if (file.exists(SAMPLE_ANCESTRY)) {
  anc_raw <- fread(SAMPLE_ANCESTRY)
  if ("K" %in% names(anc_raw)) {
    ancestry_dt <- anc_raw[K == K_choice]
  } else {
    ancestry_dt <- anc_raw
  }
}

# If not available, try to parse Q files directly
if (is.null(ancestry_dt) || nrow(ancestry_dt) == 0) {
  message("[INFO] sample_main_ancestry not usable; trying Q matrix files...")

  # Look for Q file patterns
  q_pattern <- paste0("K", K_choice)
  q_files <- list.files(Q_RUNS_DIR, pattern = q_pattern, full.names = TRUE, recursive = TRUE)
  q_files <- q_files[grepl("\\.qopt$|\\.Q$|qopt\\.gz$", q_files)]

  if (length(q_files) > 0) {
    # Pick the one from best seed if possible
    q_file <- q_files[1]
    if (file.exists(BEST_SEED_FILE)) {
      bs <- fread(BEST_SEED_FILE)
      if ("seed" %in% names(bs) && K_choice %in% bs$K) {
        best_s <- bs[K == K_choice, seed]
        seed_match <- q_files[grepl(as.character(best_s), q_files)]
        if (length(seed_match) > 0) q_file <- seed_match[1]
      }
    }

    q_mat <- fread(q_file, header = FALSE)
    q_cols <- paste0("Q", seq_len(ncol(q_mat)))
    setnames(q_mat, q_cols)

    # Need sample order
    if (file.exists(SAMPLE_IND_FILE)) {
      samp_order <- fread(SAMPLE_IND_FILE, header = FALSE)[[1]]
      samp_order <- samp_order[nchar(samp_order) > 0]
      if (length(samp_order) == nrow(q_mat)) {
        q_mat[, sample := samp_order]
      } else {
        q_mat[, sample := paste0("Ind", seq_len(.N) - 1)]
      }
    } else {
      q_mat[, sample := paste0("Ind", seq_len(.N) - 1)]
    }

    # Derive maxQ label
    q_only <- as.matrix(q_mat[, ..q_cols])
    q_mat[, maxQ_col := apply(q_only, 1, which.max)]
    q_mat[, maxQ_label := paste0("Pop", maxQ_col)]
    q_mat[, qmax := apply(q_only, 1, max)]
    q_mat[, K := K_choice]

    ancestry_dt <- q_mat
    message("[INFO] Loaded Q matrix: ", nrow(ancestry_dt), " samples, K=", K_choice)
  }
}

if (is.null(ancestry_dt) || nrow(ancestry_dt) == 0) {
  message("[WARN] No ancestry data available — writing empty ancestry outputs")
  # Still write empty template files so downstream scripts don't break
  template <- data.table(candidate_id = integer(), level_type = character(),
                         level_label = character(), n_samples = integer())
  fwrite(template, file.path(FOLLOWUP_DIR, "all_candidates_ancestry_association.tsv.gz"),
         sep = "\t")
  message("[DONE] STEP23 complete (no ancestry data)")
  quit("no", status = 0)
}

# Ensure required columns
q_cols <- grep("^Q[0-9]", names(ancestry_dt), value = TRUE)
if (!("maxQ_label" %in% names(ancestry_dt))) {
  if (length(q_cols) > 0) {
    q_only <- as.matrix(ancestry_dt[, ..q_cols])
    ancestry_dt[, maxQ_label := paste0("Pop", apply(q_only, 1, which.max))]
    ancestry_dt[, qmax := apply(q_only, 1, max)]
  }
}

message("[INFO] Ancestry table: ", nrow(ancestry_dt), " samples, ",
        length(q_cols), " Q columns, K=", K_choice)

# ── Cramér's V ───────────────────────────────────────────────────────────────
cramers_v <- function(x, y) {
  tab <- table(x, y)
  if (min(dim(tab)) < 2) return(NA_real_)
  n <- sum(tab)
  chi2 <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 2000))
  k <- min(nrow(tab), ncol(tab))
  v <- sqrt(chi2$statistic / (n * (k - 1)))
  as.numeric(v)
}

# ── Entropy helper ───────────────────────────────────────────────────────────
shannon_entropy <- function(probs) {
  probs <- probs[probs > 0]
  -sum(probs * log2(probs))
}

# ── Read candidate table ─────────────────────────────────────────────────────
cand <- fread(CANDIDATE_TABLE)
if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]

# ── Main loop ────────────────────────────────────────────────────────────────
all_assoc <- list()
all_cluster_assoc <- list()

for (ci in seq_len(nrow(cand))) {
  row <- cand[ci]
  cid <- row$candidate_id
  chr <- row$chrom

  cand_dir <- file.path(FOLLOWUP_DIR, paste0(chr, ".candidate_", cid))
  rot_file <- file.path(cand_dir, "candidate_pca_rotated.tsv")
  sub_file <- file.path(cand_dir, "candidate_subclusters.tsv")

  if (!file.exists(rot_file)) next

  rot <- fread(rot_file)
  sub <- if (file.exists(sub_file)) fread(sub_file) else NULL

  # Merge ancestry
  merged <- merge(rot, ancestry_dt[, c("sample", "maxQ_label", "qmax", q_cols), with = FALSE],
                  by = "sample", all.x = TRUE)

  if (sum(!is.na(merged$maxQ_label)) < 3) {
    message("[SKIP] Too few ancestry-matched samples for candidate ", cid)
    next
  }

  # ── Per-group ancestry analysis ────────────────────────────────────────
  for (level_col in c("coarse_group_refined", "subcluster_label")) {
    if (level_col == "subcluster_label" && is.null(sub)) next

    if (level_col == "subcluster_label") {
      merged2 <- merge(merged, sub[, .(sample, subcluster_label)], by = "sample", all.x = TRUE)
      group_var <- "subcluster_label"
    } else {
      merged2 <- copy(merged)
      group_var <- "coarse_group_refined"
    }

    levels_present <- unique(merged2[[group_var]])
    level_type <- if (level_col == "coarse_group_refined") "coarse_group" else "subcluster"

    for (lev in levels_present) {
      idx <- merged2[[group_var]] == lev
      sub_m <- merged2[idx & !is.na(maxQ_label)]
      n_lev <- nrow(sub_m)
      if (n_lev < 2) next

      # Ancestry counts
      anc_tab <- table(sub_m$maxQ_label)
      dominant <- names(anc_tab)[which.max(anc_tab)]
      dominant_frac <- max(anc_tab) / sum(anc_tab)

      # Entropy
      probs <- as.numeric(anc_tab) / sum(anc_tab)
      ent <- shannon_entropy(probs)

      # Mean Q values
      q_means <- colMeans(sub_m[, ..q_cols], na.rm = TRUE)

      assoc_row <- data.table(
        candidate_id = cid,
        level_type = level_type,
        level_label = lev,
        n_samples = n_lev,
        dominant_maxQ = dominant,
        dominant_maxQ_fraction = round(dominant_frac, 4),
        ancestry_entropy = round(ent, 4),
        n_ancestry_labels = length(anc_tab)
      )
      # Add mean Q columns
      for (qc in q_cols) {
        assoc_row[[paste0("mean_", qc)]] <- round(q_means[[qc]], 4)
      }

      all_assoc[[length(all_assoc) + 1]] <- assoc_row
    }

    # Cramér's V: maxQ vs group
    valid <- merged2[!is.na(maxQ_label) & !is.na(get(group_var))]
    if (nrow(valid) >= 5 && length(unique(valid$maxQ_label)) >= 2 &&
        length(unique(valid[[group_var]])) >= 2) {
      cv <- cramers_v(valid$maxQ_label, valid[[group_var]])

      # Chi-square p-value
      chi_p <- tryCatch({
        ct <- chisq.test(table(valid$maxQ_label, valid[[group_var]]),
                         simulate.p.value = TRUE, B = 5000)
        ct$p.value
      }, error = function(e) NA_real_)

      # Correlation of v with Q components
      v_q_cors <- sapply(q_cols, function(qc) {
        suppressWarnings(cor(valid$v, valid[[qc]], use = "complete.obs"))
      })

      assoc_summary <- data.table(
        candidate_id = cid,
        level_type = level_type,
        cramers_v = round(cv, 4),
        chi_sq_pvalue = chi_p,
        n_samples_tested = nrow(valid),
        strongest_v_Q_cor = round(max(abs(v_q_cors), na.rm = TRUE), 4),
        strongest_v_Q_component = q_cols[which.max(abs(v_q_cors))]
      )

      # Interpretation label
      interp <- if (!is.na(cv) && cv > 0.4) {
        "ancestry_substructure"
      } else if (!is.na(cv) && cv > 0.2) {
        "moderate_ancestry_association"
      } else {
        "no_strong_ancestry_association"
      }
      assoc_summary[, interpretation := interp]

      all_cluster_assoc[[length(all_cluster_assoc) + 1]] <- assoc_summary
    }
  }

  # ── Save per-candidate ─────────────────────────────────────────────────
  if (length(all_assoc) > 0) {
    cand_assoc <- rbindlist(all_assoc[sapply(all_assoc, function(x) x$candidate_id[1] == cid)],
                            fill = TRUE)
    if (nrow(cand_assoc) > 0) {
      fwrite(cand_assoc, file.path(cand_dir, "candidate_ancestry_association.tsv"),
             sep = "\t")
    }
  }

  message("[INFO] Candidate ", cid, " ancestry analysis done")
}

# ── Global outputs ───────────────────────────────────────────────────────────
if (length(all_assoc) > 0) {
  fwrite(rbindlist(all_assoc, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_ancestry_association.tsv.gz"),
         sep = "\t")
}
if (length(all_cluster_assoc) > 0) {
  fwrite(rbindlist(all_cluster_assoc, fill = TRUE),
         file.path(FOLLOWUP_DIR, "all_candidates_ancestry_cramers_summary.tsv.gz"),
         sep = "\t")
}

message("[DONE] STEP23 complete")
