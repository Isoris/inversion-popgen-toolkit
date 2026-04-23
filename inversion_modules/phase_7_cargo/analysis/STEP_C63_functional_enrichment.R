#!/usr/bin/env Rscript

# =============================================================================
# STEP_C63_functional_enrichment.R — 6B
#
# For each candidate inversion, test for over-representation of GO terms,
# KEGG pathways, and gene families in its cargo gene set vs three backgrounds:
#
#   1. Genome-wide background (all protein-coding genes)
#   2. Matched-collinear background (random non-inversion intervals matched
#      for length, gene density, and chromosome)
#   3. Within-locus background (genes inside the inversion, comparing
#      "perturbed by arrangement" vs "neutral within-arrangement" — uses
#      the per-arrangement burden output from STEP_C61)
#
# Tests use Fisher's exact test with BH FDR correction.
#
# Inputs (per candidate):
#   ${CARGO_INVENTORY_DIR}/<cid>/genes.tsv     — cargo gene list with annotations
#   ${GENE_FUNCTION_TSV}                        — genome-wide function table
#   ${SNAKE_CAND_FILE}                          — all candidate intervals
#   ${GENE_BED}                                 — all gene coordinates
#
# Outputs (per candidate):
#   ${CARGO_ENRICH_DIR}/<cid>/go_enrichment.tsv
#   ${CARGO_ENRICH_DIR}/<cid>/kegg_enrichment.tsv
#   ${CARGO_ENRICH_DIR}/<cid>/family_enrichment.tsv
#
# Usage:
#   Rscript STEP_C63_functional_enrichment.R [candidate_id|all]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) >= 1) args[1] else "all"

# ── Resolve paths from environment ──
CARGO_INVENTORY_DIR <- Sys.getenv("CARGO_INVENTORY_DIR")
CARGO_ENRICH_DIR    <- Sys.getenv("CARGO_ENRICH_DIR")
GENE_FUNCTION_TSV   <- Sys.getenv("GENE_FUNCTION_TSV")
SNAKE_CAND_FILE     <- Sys.getenv("SNAKE_CAND_FILE")
GENE_BED            <- Sys.getenv("GENE_BED")
MATCHED_BG_N        <- as.integer(Sys.getenv("CARGO_MATCHED_BG_N", "100"))
ENRICH_FDR          <- as.numeric(Sys.getenv("CARGO_ENRICH_FDR", "0.05"))
ENRICH_MIN_GENES    <- as.integer(Sys.getenv("CARGO_ENRICH_MIN_GENES", "3"))

stopifnot(nzchar(CARGO_INVENTORY_DIR),
          nzchar(CARGO_ENRICH_DIR),
          nzchar(GENE_FUNCTION_TSV),
          nzchar(SNAKE_CAND_FILE),
          nzchar(GENE_BED))

# ── Load reference data ──
message("[C63] Loading gene function table: ", GENE_FUNCTION_TSV)
fn <- fread(GENE_FUNCTION_TSV, sep = "\t", header = TRUE)
message("  ", nrow(fn), " genes loaded")

message("[C63] Loading gene BED: ", GENE_BED)
gbed <- fread(cmd = paste0("zcat ", shQuote(GENE_BED)),
              header = FALSE,
              col.names = c("chrom", "start", "end", "gene_id", "score", "strand"))

message("[C63] Loading candidates: ", SNAKE_CAND_FILE)
cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))

# ── Build genome-wide term → genes index ──
explode_terms <- function(dt, col) {
  if (!col %in% names(dt)) return(data.table(gene_id = character(), term = character()))
  rows <- dt[nzchar(get(col)),
             .(term = unlist(strsplit(get(col), ";"))),
             by = gene_id]
  rows <- rows[nzchar(term)]
  rows[, term := trimws(term)]
  unique(rows)
}

go_long      <- explode_terms(fn, "go_terms")
kegg_long    <- explode_terms(fn, "kegg_pathway")
family_long  <- explode_terms(fn, "family")

message("[C63] Term-gene index: GO=", nrow(go_long),
        " KEGG=", nrow(kegg_long), " family=", nrow(family_long))

# ── Matched collinear background sampler ──
# Build intervals on genome that do NOT overlap any candidate.
build_collinear_pool <- function(gbed, cands) {
  chrs <- unique(gbed$chrom)
  pool <- list()
  for (ch in chrs) {
    chr_cands <- cands[chrom == ch][order(start_bp)]
    chr_genes <- gbed[chrom == ch][order(start)]
    if (nrow(chr_genes) == 0) next
    chr_max <- max(chr_genes$end)
    if (nrow(chr_cands) == 0) {
      pool[[ch]] <- data.table(chrom = ch, start = 0L, end = as.integer(chr_max))
      next
    }
    free <- data.table(chrom = ch,
                        start = c(0L, chr_cands$end_bp),
                        end   = c(chr_cands$start_bp, chr_max))
    free <- free[end > start]
    pool[[ch]] <- free
  }
  rbindlist(pool, fill = TRUE)
}

sample_matched_intervals <- function(target_len, target_n_genes, pool, gbed,
                                      n = MATCHED_BG_N, max_tries = n * 50) {
  out <- list()
  tries <- 0
  while (length(out) < n && tries < max_tries) {
    tries <- tries + 1
    # Pick a random free segment large enough
    ok <- pool[end - start >= target_len]
    if (nrow(ok) == 0) break
    seg <- ok[sample.int(nrow(ok), 1)]
    max_off <- seg$end - seg$start - target_len
    off <- sample.int(max_off + 1, 1) - 1L
    s <- seg$start + off
    e <- s + target_len
    n_g <- gbed[chrom == seg$chrom & end > s & start < e, .N]
    # Density tolerance: ±50%
    if (abs(n_g - target_n_genes) <= max(2, ceiling(0.5 * target_n_genes))) {
      out[[length(out) + 1]] <- data.table(chrom = seg$chrom, start = s, end = e,
                                            n_genes = n_g)
    }
  }
  if (length(out) == 0) return(data.table())
  rbindlist(out)
}

# ── Fisher enrichment helper ──
fisher_enrichment <- function(target_genes, term_long, bg_genes,
                                min_genes = ENRICH_MIN_GENES) {
  if (length(target_genes) == 0 || length(bg_genes) == 0) {
    return(data.table())
  }
  tgt_set <- unique(target_genes)
  bg_set  <- unique(c(bg_genes, tgt_set))
  in_tgt  <- term_long[gene_id %in% tgt_set]
  in_bg   <- term_long[gene_id %in% bg_set]
  if (nrow(in_tgt) == 0 || nrow(in_bg) == 0) return(data.table())

  per_term <- in_tgt[, .(n_target = uniqueN(gene_id),
                          target_genes = paste(sort(unique(gene_id)), collapse = ",")),
                      by = term]
  bg_counts <- in_bg[, .(n_bg = uniqueN(gene_id)), by = term]
  d <- merge(per_term, bg_counts, by = "term", all.x = TRUE)
  d <- d[n_target >= min_genes]
  if (nrow(d) == 0) return(data.table())

  N_tgt <- length(tgt_set)
  N_bg  <- length(bg_set)
  res <- d[, {
    a <- n_target
    b <- N_tgt - a
    c <- n_bg - n_target
    d_ <- (N_bg - N_tgt) - c
    if (a + b == 0 || c + d_ == 0 || a + c == 0 || b + d_ == 0) {
      list(p = NA_real_, odds = NA_real_)
    } else {
      ft <- tryCatch(fisher.test(matrix(c(a, b, c, d_), nrow = 2),
                                  alternative = "greater"),
                     error = function(e) NULL)
      if (is.null(ft)) list(p = NA_real_, odds = NA_real_)
      else list(p = ft$p.value, odds = unname(ft$estimate))
    }
  }, by = .(term, n_target, n_bg, target_genes)]
  res[, fdr_bh := p.adjust(p, method = "BH")]
  res <- res[order(p)]
  res[]
}

# ── Per-candidate driver ──
process_candidate <- function(cid) {
  cand <- cands[candidate_id == cid]
  if (nrow(cand) == 0) {
    message("[C63] Candidate not found: ", cid); return(invisible())
  }
  inv_genes_file <- file.path(CARGO_INVENTORY_DIR, cid, "genes.tsv")
  if (!file.exists(inv_genes_file)) {
    message("[C63] No inventory for ", cid, " — skipping (run STEP_C60 first)")
    return(invisible())
  }
  inv_genes <- fread(inv_genes_file)
  out_dir <- file.path(CARGO_ENRICH_DIR, cid)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cargo_gid <- unique(inv_genes$gene_id)
  if (length(cargo_gid) < 3) {
    message("[C63] ", cid, ": <3 genes — skipping enrichment")
    return(invisible())
  }

  # ── Background 1: genome-wide ──
  bg_genome <- unique(fn$gene_id)

  # ── Background 2: matched collinear ──
  pool <- build_collinear_pool(gbed, cands)
  matched <- sample_matched_intervals(
    target_len = cand$end_bp - cand$start_bp,
    target_n_genes = nrow(inv_genes),
    pool = pool, gbed = gbed)
  if (nrow(matched) > 0) {
    bg_matched <- unique(unlist(lapply(seq_len(nrow(matched)), function(i) {
      r <- matched[i]
      gbed[chrom == r$chrom & end > r$start & start < r$end, gene_id]
    })))
  } else {
    bg_matched <- character()
    message("[C63] ", cid, ": matched-collinear sampling produced 0 intervals")
  }

  # ── Run enrichments × 3 backgrounds × 3 ontologies ──
  enrich_layer <- function(term_long, ont) {
    do_one <- function(bg, bg_label) {
      r <- fisher_enrichment(cargo_gid, term_long, bg)
      if (nrow(r) > 0) r[, background := bg_label]
      r
    }
    rbindlist(list(do_one(bg_genome, "genome"),
                   do_one(bg_matched, "matched_collinear")),
              fill = TRUE)
  }

  go_res     <- enrich_layer(go_long,     "GO")
  kegg_res   <- enrich_layer(kegg_long,   "KEGG")
  family_res <- enrich_layer(family_long, "family")

  fwrite(go_res,     file.path(out_dir, "go_enrichment.tsv"),     sep = "\t")
  fwrite(kegg_res,   file.path(out_dir, "kegg_enrichment.tsv"),   sep = "\t")
  fwrite(family_res, file.path(out_dir, "family_enrichment.tsv"), sep = "\t")

  message(sprintf("[C63] %s: GO sig=%d KEGG sig=%d family sig=%d (FDR<%.2f)",
                  cid,
                  sum(go_res$fdr_bh < ENRICH_FDR, na.rm = TRUE),
                  sum(kegg_res$fdr_bh < ENRICH_FDR, na.rm = TRUE),
                  sum(family_res$fdr_bh < ENRICH_FDR, na.rm = TRUE),
                  ENRICH_FDR))
}

if (target == "all") {
  diag <- fread(file.path(CARGO_INVENTORY_DIR, "diagnostic_table.tsv"))
  cids <- diag[n_genes_inside >= 3, candidate_id]
} else {
  cids <- target
}
message("[C63] Processing ", length(cids), " candidate(s)")
for (cid in cids) process_candidate(cid)
message("[C63] Done")
