#!/usr/bin/env Rscript

# =============================================================================
# PLOT_03_variant_burden_dashboard.R
#
# Three-panel cohort QC figure (image 5):
#   A: per-sample total variant burden, stacked by class (sorted by total)
#   B: per-group total variant burden, boxplot
#   C: per-chromosome variant burden by class, heatmap
#
# Inputs:
#   ${VARIANT_MASTER}                — variant_master_scored.tsv
#                                       must have: chr, pos, snpeff_annotation
#                                       and per-sample carrier columns OR a
#                                       'samples_with_alt' semicolon-list
#   --groups-from <tsv>             — sample group_id (default: ancestry K)
#
# Output:
#   ${EXTRAS_FIG_DIR}/PLOT_03_variant_burden_dashboard.pdf
#   ${EXTRAS_TBL_DIR}/PLOT_03_per_sample_burden.tsv
#   ${EXTRAS_TBL_DIR}/PLOT_03_per_chrom_burden.tsv
#
# Class derivation:
#   SNP: snpeff_annotation grep "missense|synonymous|stop_|splice_" + length(ref)==1 && length(alt)==1
#   INDEL: length(ref) != length(alt) and both <= 50bp
#   DEL/DUP/INV/BND/INS: from sv_type column if present
#   PAV: from variant_class column if present (catches PAV)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

VARIANT_MASTER  <- Sys.getenv("VARIANT_MASTER")
SAMPLE_LIST     <- Sys.getenv("SAMPLE_LIST")
NGSADMIX_Q_FILE <- Sys.getenv("NGSADMIX_Q_FILE")
EXTRAS_FIG_DIR  <- Sys.getenv("EXTRAS_FIG_DIR")
EXTRAS_TBL_DIR  <- Sys.getenv("EXTRAS_TBL_DIR")
CANONICAL_K     <- as.integer(Sys.getenv("CANONICAL_K", "8"))

if (!file.exists(VARIANT_MASTER)) {
  message("[PLOT_03] [skip] VARIANT_MASTER not found: ", VARIANT_MASTER)
  quit(status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
groups_from <- NULL
if (any(args == "--groups-from")) {
  groups_from <- args[which(args == "--groups-from") + 1]
}

message("[PLOT_03] Reading variant master: ", VARIANT_MASTER)
v <- fread(VARIANT_MASTER)

# ── Classify variants ────────────────────────────────────────────────────────
classify <- function(v) {
  cls <- rep("OTHER", nrow(v))
  has_sv <- "sv_type" %in% names(v)
  has_vc <- "variant_class" %in% names(v)
  ref_len <- if ("ref" %in% names(v)) nchar(v$ref) else NA
  alt_len <- if ("alt" %in% names(v)) nchar(v$alt) else NA
  if (!is.na(ref_len[1])) {
    is_snp <- ref_len == 1 & alt_len == 1
    is_indel <- ref_len != alt_len & pmax(ref_len, alt_len) <= 50
    cls[is_snp]   <- "SNP"
    cls[is_indel] <- "INDEL"
  }
  if (has_sv) {
    sv <- toupper(v$sv_type)
    cls[sv == "DEL"] <- "DEL"
    cls[sv == "DUP"] <- "DUP"
    cls[sv == "INV"] <- "INV"
    cls[sv == "BND"] <- "BND"
    cls[sv == "INS"] <- "INS"
  }
  if (has_vc) cls[toupper(v$variant_class) == "PAV"] <- "PAV"
  cls
}
v[, class := classify(.SD)]
message("[PLOT_03] Variant classes: ",
        paste(names(table(v$class)), table(v$class), sep = "=", collapse = ", "))

# ── Sample list ──────────────────────────────────────────────────────────────
samples <- if (file.exists(SAMPLE_LIST)) {
  as.character(fread(SAMPLE_LIST, header = FALSE)[[1]])
} else stop("SAMPLE_LIST not found")

# ── Resolve per-sample carriership ───────────────────────────────────────────
# Try in this order:
#   1. samples_with_alt column (semicolon list)
#   2. per-sample dosage columns matching sample names
#   3. carriers column (count) — give up on per-sample, only do per-chrom

per_sample_counts <- NULL
if ("samples_with_alt" %in% names(v)) {
  message("[PLOT_03] Using samples_with_alt column")
  long <- v[, .(sample = unlist(strsplit(samples_with_alt, ";"))),
            by = .(class, chr = if ("chr" %in% names(v)) chr else chrom)]
  long <- long[nzchar(sample)]
  per_sample_counts <- long[, .N, by = .(sample, class)]
} else {
  sample_cols <- intersect(samples, names(v))
  if (length(sample_cols) >= 5) {
    message("[PLOT_03] Using per-sample dosage columns (n=", length(sample_cols), ")")
    # Vectorised melt instead of per-sample loop — much faster on a wide table.
    # Use chunked melts to keep memory bounded on a 226-sample × 384k-variant table.
    chunk_n <- 30L
    rows <- list()
    for (chunk_start in seq(1, length(sample_cols), by = chunk_n)) {
      chunk_end <- min(chunk_start + chunk_n - 1L, length(sample_cols))
      cols_chunk <- sample_cols[chunk_start:chunk_end]
      mlt <- melt(v[, c("class", cols_chunk), with = FALSE],
                   id.vars = "class",
                   variable.name = "sample",
                   value.name = "dosage")
      rows[[length(rows) + 1]] <- mlt[dosage > 0, .N, by = .(sample, class)]
    }
    per_sample_counts <- rbindlist(rows)
    per_sample_counts[, sample := as.character(sample)]
  } else {
    message("[PLOT_03] [warn] no per-sample data — Panel A/B will be skipped")
  }
}

# ── Resolve per-sample groups ────────────────────────────────────────────────
sample_to_group <- NULL
if (!is.null(groups_from) && file.exists(groups_from)) {
  gd <- fread(groups_from)
  sample_to_group <- setNames(gd$group_id, gd$sample)
} else if (file.exists(NGSADMIX_Q_FILE)) {
  Q <- as.matrix(fread(NGSADMIX_Q_FILE, header = FALSE))
  n <- min(length(samples), nrow(Q))
  argmax_k <- apply(Q[seq_len(n), , drop = FALSE], 1, which.max)
  sample_to_group <- setNames(paste0("K", CANONICAL_K, "_Q", argmax_k),
                                samples[seq_len(n)])
}

CLASS_PAL <- c(SNP = "#1f77b4", INDEL = "#2ca02c",
                DEL = "#d62728", DUP = "#ff7f0e",
                INV = "#9467bd", INS = "#e377c2",
                BND = "#8c564b", PAV = "#7f7f7f", OTHER = "grey80")

# ── Panel A: per-sample stacked bars ─────────────────────────────────────────
pA <- NULL
if (!is.null(per_sample_counts)) {
  totals <- per_sample_counts[, .(tot = sum(N)), by = sample][order(tot)]
  per_sample_counts[, sample := factor(sample, levels = totals$sample)]
  per_sample_counts[, class := factor(class, levels = names(CLASS_PAL))]
  fwrite(per_sample_counts,
         file.path(EXTRAS_TBL_DIR, "PLOT_03_per_sample_burden.tsv"),
         sep = "\t")
  pA <- ggplot(per_sample_counts, aes(x = sample, y = N, fill = class)) +
    geom_col() +
    scale_fill_manual(values = CLASS_PAL, drop = FALSE) +
    labs(title = "A. Per-sample variant burden",
         x = NULL, y = "Variant count") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(face = "bold"))
}

# ── Panel B: per-group boxplot ───────────────────────────────────────────────
pB <- NULL
if (!is.null(per_sample_counts) && !is.null(sample_to_group)) {
  totals <- per_sample_counts[, .(N = sum(N)), by = sample]
  totals[, group := sample_to_group[sample]]
  totals <- totals[!is.na(group)]
  pB <- ggplot(totals, aes(x = group, y = N, fill = group)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    labs(title = "B. Per-group total variant burden",
         x = NULL, y = "Variants per sample") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(face = "bold"))
}

# ── Panel C: per-chromosome heatmap ──────────────────────────────────────────
chr_col <- if ("chr" %in% names(v)) "chr" else if ("chrom" %in% names(v)) "chrom" else NULL
if (is.null(chr_col)) stop("variant master has no chr/chrom column")
per_chrom <- v[, .N, by = c(chr_col, "class")]
setnames(per_chrom, chr_col, "chrom")
per_chrom[, chrom := factor(chrom, levels = unique(chrom))]
per_chrom[, class := factor(class, levels = names(CLASS_PAL))]
fwrite(per_chrom, file.path(EXTRAS_TBL_DIR, "PLOT_03_per_chrom_burden.tsv"),
       sep = "\t")
pC <- ggplot(per_chrom, aes(x = class, y = chrom, fill = N)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = formatC(N, format = "d", big.mark = ",")),
            size = 2.5, color = "black") +
  scale_fill_viridis_c(name = "Count", trans = "log10",
                       option = "viridis", direction = -1) +
  labs(title = "C. Per-chromosome variant burden",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))

# ── Compose ──────────────────────────────────────────────────────────────────
panels <- list(pA, pB, pC)
panels <- panels[!vapply(panels, is.null, logical(1))]
fig <- if (length(panels) == 3) {
  (pA / pB) | pC
} else if (length(panels) == 2) {
  panels[[1]] | panels[[2]]
} else {
  panels[[1]]
}
fig <- fig + plot_annotation(title = "Variant burden dashboard")

ggsave(file.path(EXTRAS_FIG_DIR, "PLOT_03_variant_burden_dashboard.pdf"),
       fig, width = 14, height = 10, device = cairo_pdf)
message("[PLOT_03] Done")
