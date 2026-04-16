#!/usr/bin/env Rscript
# =============================================================================
# cheat21_te_context.R — TE enrichment/depletion test
#
# BIOLOGY:
#   Don't assume TEs cause inversions — test it. Some species show TE
#   enrichment at breakpoints (repeat-mediated), others depletion.
#   The catfish-specific answer matters for mechanism interpretation.
#
# INPUT:  RepeatMasker .out/.bed, boundary catalog, genome FAI
# OUTPUT: per-breakpoint TE density, enrichment vs random, per-TE-class
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
TE_WINDOW_BP   <- 10000L
N_RANDOM       <- 10000L
EXT_FLANK_BP   <- 50000L

# ── Load RepeatMasker annotation ──────────────────────────────────────

load_repeat_annotation <- function(repeatmasker_file) {
  if (!file.exists(repeatmasker_file)) {
    message("[cheat21] RepeatMasker file not found: ", repeatmasker_file)
    return(data.table())
  }

  ext <- tools::file_ext(repeatmasker_file)
  if (ext %in% c("bed", "gz")) {
    dt <- tryCatch(fread(repeatmasker_file), error = function(e) data.table())
    if (ncol(dt) >= 4) {
      setnames(dt, seq_len(min(4, ncol(dt))),
               c("chr", "start", "end", "te_class")[seq_len(min(4, ncol(dt)))])
    }
  } else {
    # Parse .out format (skip header lines, whitespace-delimited)
    dt <- tryCatch({
      lines <- readLines(repeatmasker_file)
      lines <- lines[!grepl("^\\s*$|^\\s*SW", lines)]
      lines <- lines[-(1:min(3, length(lines)))]
      raw <- fread(text = lines, fill = TRUE)
      if (ncol(raw) >= 11) {
        data.table(chr = raw[[5]], start = as.integer(raw[[6]]),
                    end = as.integer(raw[[7]]),
                    te_class = sub("/.*", "", raw[[11]]),
                    te_family = raw[[11]], te_name = raw[[10]])
      } else data.table()
    }, error = function(e) data.table())
  }

  if (nrow(dt) == 0) return(dt)
  dt[, te_bp := end - start]
  dt
}

# ── Compute TE density around positions ───────────────────────────────

compute_te_density <- function(positions, chr_name, repeat_ann,
                                window = TE_WINDOW_BP) {
  if (nrow(repeat_ann) == 0 || length(positions) == 0)
    return(data.table(position = positions, te_density = NA_real_,
                       te_count = 0L))
  chr_te <- repeat_ann[chr == chr_name]
  if (nrow(chr_te) == 0)
    return(data.table(position = positions, te_density = NA_real_,
                       te_count = 0L))

  results <- list()
  for (i in seq_along(positions)) {
    pos <- positions[i]
    local <- chr_te[start <= pos + window & end >= pos - window]
    overlap_bp <- sum(pmin(local$end, pos + window) -
                       pmax(local$start, pos - window))
    overlap_bp <- max(0, overlap_bp)
    density <- overlap_bp / (2 * window)

    class_counts <- if (nrow(local) > 0 && "te_class" %in% names(local))
      as.list(local[, .N, by = te_class][, setNames(N, te_class)])
    else list()

    results[[i]] <- data.table(position = pos,
                                te_density = round(density, 4),
                                te_count = nrow(local))
  }
  rbindlist(results)
}

# ── Permutation enrichment test ───────────────────────────────────────

te_enrichment_test <- function(breakpoint_positions, chr_name,
                                repeat_ann, genome_fai,
                                n_random = N_RANDOM) {
  bp_density <- compute_te_density(breakpoint_positions, chr_name, repeat_ann)
  if (all(is.na(bp_density$te_density)))
    return(list(enrichment_ratio = NA_real_, enrichment_p = NA_real_,
                verdict = "NO_DATA"))

  # Load chromosome sizes
  fai <- tryCatch(fread(genome_fai, header = FALSE,
                         col.names = c("chr","size","offset","bases","bytes")),
                   error = function(e) data.table())
  chr_size <- if (nrow(fai) > 0) fai[chr == chr_name]$size[1] else 1e8

  # Random positions
  random_pos <- sample(seq(10000, chr_size - 10000), min(n_random, chr_size/100))
  rnd_density <- compute_te_density(random_pos, chr_name, repeat_ann)

  bp_mean  <- mean(bp_density$te_density, na.rm = TRUE)
  rnd_mean <- mean(rnd_density$te_density, na.rm = TRUE)
  ratio    <- if (rnd_mean > 0) round(bp_mean / rnd_mean, 3) else NA_real_

  # Permutation p-value
  rnd_vals <- rnd_density$te_density[!is.na(rnd_density$te_density)]
  if (length(rnd_vals) > 0 && !is.na(bp_mean)) {
    p_val <- mean(rnd_vals >= bp_mean)
  } else {
    p_val <- NA_real_
  }

  verdict <- if (is.na(ratio)) "NO_DATA"
    else if (ratio > 1.5 && p_val < 0.05) "ENRICHED"
    else if (ratio < 0.67 && p_val > 0.95) "DEPLETED"
    else "NEUTRAL"

  list(enrichment_ratio = ratio, enrichment_p = round(p_val, 4),
       bp_mean_density = round(bp_mean, 4),
       random_mean_density = round(rnd_mean, 4),
       verdict = verdict)
}

# ── Four-way comparison: breakpoint/interior/external/background ─────

te_breakpoints_vs_interior <- function(boundary_bp, candidate_start,
                                        candidate_end, chr_name,
                                        repeat_ann, genome_fai,
                                        ext_flank = EXT_FLANK_BP) {
  regions <- list(
    breakpoint_flanks = boundary_bp,
    interior = seq(candidate_start + 10000, candidate_end - 10000,
                    length.out = min(20, (candidate_end - candidate_start) / 10000)),
    external = c(seq(candidate_start - ext_flank, candidate_start - 10000,
                      length.out = 5),
                  seq(candidate_end + 10000, candidate_end + ext_flank,
                      length.out = 5))
  )

  results <- list()
  for (rname in names(regions)) {
    dens <- compute_te_density(regions[[rname]], chr_name, repeat_ann)
    results[[rname]] <- data.table(region = rname,
                                    mean_density = round(mean(dens$te_density, na.rm = TRUE), 4),
                                    n_positions = nrow(dens))
  }
  rbindlist(results)
}

# ── Search mode ────────────────────────────────────────────────────────

search_te_context <- function(chr, zone_start, zone_end,
                               repeat_ann = NULL, ...) {
  empty <- data.table(method = "te_context", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(repeat_ann) || nrow(repeat_ann) == 0) return(empty)
  mid <- as.integer((zone_start + zone_end) / 2)
  dens <- compute_te_density(mid, chr, repeat_ann)
  sc <- if (!is.na(dens$te_density)) pmin(1, dens$te_density * 2) else 0
  data.table(method = "te_context", best_bp = mid,
             score = round(sc, 3), is_precise = FALSE,
             detail = paste0("te_density=", dens$te_density))
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat21 <- function(chr, boundary_bp, candidate_start, candidate_end,
                         repeatmasker_file, genome_fai) {
  message("[cheat21] ", chr, ": ", length(boundary_bp), " breakpoints")

  repeat_ann <- load_repeat_annotation(repeatmasker_file)
  message("[cheat21] Loaded ", nrow(repeat_ann), " TE annotations")

  if (nrow(repeat_ann) == 0) {
    message("[cheat21] No TE data — skipping")
    return(list(enrichment = list(verdict = "NO_DATA"),
                four_way = data.table(),
                search_result = search_te_context(chr, candidate_start,
                  candidate_end)))
  }

  enrich <- te_enrichment_test(boundary_bp, chr, repeat_ann, genome_fai)
  message("[cheat21] Enrichment: ", enrich$enrichment_ratio,
          "× (p=", enrich$enrichment_p, ") → ", enrich$verdict)

  four_way <- te_breakpoints_vs_interior(boundary_bp, candidate_start,
                                          candidate_end, chr, repeat_ann,
                                          genome_fai)
  message("[cheat21] Four-way:\n",
          paste(capture.output(print(four_way)), collapse = "\n"))

  list(enrichment = enrich, four_way = four_way,
       search_result = search_te_context(chr, candidate_start, candidate_end,
                                          repeat_ann))
}
