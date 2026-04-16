#!/usr/bin/env Rscript

# =============================================================================
# config_inversion_followup.R
#
# Central configuration for the stripe-aware inversion candidate interpretation
# framework (MODULE_5B).  Source this file at the top of every STEP20+ script.
#
# All paths point to the HPC layout under PROJECT_ROOT.
# =============================================================================

# ── Project root ─────────────────────────────────────────────────────────────
PROJECT_ROOT <- "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"

# ── Input roots ──────────────────────────────────────────────────────────────
INV_ROOT          <- file.path(PROJECT_ROOT, "inversion_localpca_v7")
CANDIDATE_TABLE   <- {
  snake_cand <- file.path(INV_ROOT, "06_mds_candidates",
                          "snake_regions_multiscale",
                          "snake_candidate_regions.tsv.gz")
  control_cand <- file.path(INV_ROOT, "06_mds_candidates",
                            "inversion_localpca.candidate_regions.tsv.gz")
  if (file.exists(snake_cand)) snake_cand else control_cand
}
DOSAGE_DIR        <- file.path(INV_ROOT, "04_dosage_by_chr")
STEP09_DIR        <- file.path(INV_ROOT, "05_local_pca")
STEP10_DIR        <- file.path(INV_ROOT, "06_mds_candidates")
STEP12_DIR        <- file.path(INV_ROOT, "08_regional_pca")
THETA_OVERLAP_DIR <- file.path(INV_ROOT, "07_het_overlap")
GROUP_TABLE_DIR   <- file.path(INV_ROOT, "11_candidate_group_tables")
COHERENCE_DIR     <- file.path(INV_ROOT, "12_block_coherence")
SCORING_DIR       <- file.path(INV_ROOT, "13_scoring_table")
COMPOSITE_DIR     <- file.path(INV_ROOT, "15_composite_figures")

# ── Ancestry / Q paths ──────────────────────────────────────────────────────
ANCESTRY_ROOT     <- file.path(PROJECT_ROOT, "popstruct_thin", "05_ngsadmix_global")
BEST_SEED_FILE    <- file.path(ANCESTRY_ROOT, "best_seed_by_K.tsv")
SAMPLE_ANCESTRY   <- file.path(ANCESTRY_ROOT, "sample_main_ancestry_by_K.tsv")
Q_RUNS_DIR        <- file.path(ANCESTRY_ROOT, "runs_thin500")

# ── Reference / annotation ──────────────────────────────────────────────────
REF_FASTA         <- file.path(PROJECT_ROOT, "00-samples", "fClaHyb_Gar_LG.fa")
REF_FAI           <- file.path(PROJECT_ROOT, "00-samples", "fClaHyb_Gar_LG.fa.fai")
GFF_FILE          <- file.path(PROJECT_ROOT, "00-samples",
                               "fClaHyb_Gar_LG.from_CGAR.gff3_polished.sorted.gff3")

# ── Sample metadata ─────────────────────────────────────────────────────────
SAMPLE_IND_FILE   <- file.path(PROJECT_ROOT, "het_roh", "01_inputs_check", "samples.ind")
CHROM_SIZES_FILE  <- file.path(PROJECT_ROOT, "het_roh", "01_inputs_check", "chrom_sizes.tsv")

# ── New output roots (additive, never overwrite old dirs) ────────────────────
FOLLOWUP_DIR      <- file.path(INV_ROOT, "20_candidate_followup")
PLOTS_DIR         <- file.path(INV_ROOT, "21_candidate_plots")

# ── Color palettes ───────────────────────────────────────────────────────────
GROUP_COLORS <- c(
  "HOMO_1"    = "#2166AC",
  "HET"       = "#D6604D",
  "HOMO_2"    = "#1B7837",
  "AMBIGUOUS" = "#878787",
  "Homo_1"    = "#2166AC",
  "Het"       = "#D6604D",
  "Homo_2"    = "#1B7837"
)

STATE_COLORS <- c(
  "FULL_A"  = "#2166AC",
  "HALF"    = "#D6604D",
  "FULL_B"  = "#1B7837",
  "NOISE"   = "#BDBDBD",
  "COMPLEX" = "#9970AB"
)

CONFIDENCE_COLORS <- c(
  "high"   = "#1A9850",
  "medium" = "#F46D43",
  "low"    = "#D73027"
)

TOPOLOGY_COLORS <- c(
  "clean_3band"       = "#2166AC",
  "tilted"            = "#F4A582",
  "split_het"         = "#D6604D",
  "split_homo"        = "#B2182B",
  "gradient"          = "#FDB863",
  "curved"            = "#E08214",
  "weak"              = "#D9D9D9",
  "ambiguous"         = "#BDBDBD",
  "multimodal"        = "#9970AB",
  "diffuse"           = "#E0E0E0"
)

GEOMETRY_LABELS <- c("compact", "split_discrete", "continuous_gradient",
                     "curved", "diffuse")

# ── Theme helper ─────────────────────────────────────────────────────────────
theme_inversion <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, color = "grey30"),
      strip.background = ggplot2::element_rect(fill = "#F0F0F0"),
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ── Utility functions ────────────────────────────────────────────────────────

#' Format bp as human-readable Mb string
fmt_mb <- function(bp) sprintf("%.2f Mb", bp / 1e6)

#' Build informative candidate header string
candidate_header <- function(cid, chrom, start_bp, end_bp,
                             n_snps = NA, n_samples = NA,
                             n_windows = NA, group_counts = NULL,
                             interpretation = NA) {
  span <- end_bp - start_bp
  parts <- c(
    paste0("Candidate ", cid),
    chrom,
    paste0(fmt_mb(start_bp), " – ", fmt_mb(end_bp)),
    paste0("span ", fmt_mb(span))
  )
  if (!is.na(n_snps))    parts <- c(parts, paste0(format(n_snps, big.mark = ","), " SNPs"))
  if (!is.na(n_samples)) parts <- c(parts, paste0(n_samples, " samples"))
  if (!is.na(n_windows)) parts <- c(parts, paste0(n_windows, " windows"))
  if (!is.null(group_counts)) {
    gc_str <- paste0(names(group_counts), "=", group_counts, collapse = " ")
    parts <- c(parts, gc_str)
  }
  if (!is.na(interpretation)) parts <- c(parts, interpretation)
  paste(parts, collapse = " | ")
}

#' Ensure output directory exists
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# ── Subtitle wrapping (prevents overflow) ────────────────────────────────────
#' Wrap a long subtitle into multiple lines at ~max_chars per line
wrap_subtitle <- function(text, max_chars = 90) {
  if (nchar(text) <= max_chars) return(text)
  words <- strsplit(text, " \\| ")[[1]]
  lines <- character(0); cur <- ""
  for (w in words) {
    test <- if (nchar(cur) == 0) w else paste0(cur, " | ", w)
    if (nchar(test) > max_chars && nchar(cur) > 0) {
      lines <- c(lines, cur)
      cur <- w
    } else {
      cur <- test
    }
  }
  if (nchar(cur) > 0) lines <- c(lines, cur)
  paste(lines, collapse = "\n")
}

# ── Safe file reader ─────────────────────────────────────────────────────────
safe_fread <- function(path) {
  tryCatch(fread(path), error = function(e) NULL)
}

# ── Safe color/shape scale helpers ───────────────────────────────────────────
safe_group_colors <- function(groups, fallback = "#BDBDBD") {
  groups <- unique(as.character(groups))
  cols <- setNames(rep(fallback, length(groups)), groups)
  hit <- intersect(groups, names(GROUP_COLORS))
  cols[hit] <- GROUP_COLORS[hit]
  cols
}

SHAPE_MAP <- c("HOMO_1" = 16, "HET" = 17, "HOMO_2" = 15, "AMBIGUOUS" = 4)

safe_shape_values <- function(groups) {
  groups <- unique(as.character(groups))
  shp <- setNames(rep(16L, length(groups)), groups)
  hit <- intersect(groups, names(SHAPE_MAP))
  shp[hit] <- SHAPE_MAP[hit]
  shp
}

message("[CONFIG] Inversion followup config loaded (v7.4)")
message("[CONFIG] PROJECT_ROOT = ", PROJECT_ROOT)
message("[CONFIG] INV_ROOT     = ", INV_ROOT)
