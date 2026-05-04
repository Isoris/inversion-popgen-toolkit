# ============================================================================
# STEP_D06_sv_overlap.R — SV breakpoint overlap per candidate
# ============================================================================
# For each candidate, find the best-matching SV call from DELLY/Manta.
# Report overlap percentage and distance of breakpoints to triangle edges.
#
# SV file format expected (TSV):
#   chr  start  end  sv_type  sv_id  caller  [optional: qual, genotype_info]
#
# Returns data.frame with columns:
#   candidate_id, best_sv_id, best_sv_caller, sv_type,
#   sv_overlap_pct, sv_left_dist_kb, sv_right_dist_kb,
#   sv_span_ratio, n_sv_hits
# ============================================================================

compute_sv_overlap <- function(candidates, sv_file, chr,
                               max_dist_kb=500, window_size_bp=50000) {
  # Read SV calls
  sv <- read.delim(sv_file, stringsAsFactors=FALSE)

  # Standardize column names
  names(sv) <- tolower(names(sv))
  if (!"chr" %in% names(sv) && "chrom" %in% names(sv)) names(sv)[names(sv)=="chrom"] <- "chr"
  if (!"sv_id" %in% names(sv) && "id" %in% names(sv)) names(sv)[names(sv)=="id"] <- "sv_id"
  if (!"sv_type" %in% names(sv) && "svtype" %in% names(sv)) names(sv)[names(sv)=="svtype"] <- "sv_type"
  if (!"sv_id" %in% names(sv)) sv$sv_id <- paste0("SV_", seq_len(nrow(sv)))
  if (!"caller" %in% names(sv)) sv$caller <- "unknown"

  # Filter to chromosome
  # Handle various chr naming conventions
  chr_clean <- gsub("C_gar_", "", chr)
  sv_chr <- sv[sv$chr == chr | sv$chr == chr_clean |
               sv$chr == paste0("C_gar_", sv$chr), ]

  # Also try matching without prefix
  if (nrow(sv_chr) == 0) {
    sv_chr <- sv[grepl(chr_clean, sv$chr, fixed=TRUE), ]
  }

  results <- list()

  for (i in seq_len(nrow(candidates))) {
    cid <- candidates$candidate_id[i]
    c_start_bp <- (candidates$start[i] - 1) * window_size_bp
    c_end_bp   <- candidates$end[i] * window_size_bp
    c_width_bp <- c_end_bp - c_start_bp

    if (nrow(sv_chr) == 0) {
      results[[i]] <- make_empty_sv_row(cid)
      next
    }

    # Find all SVs within max_dist_kb of the candidate
    max_dist_bp <- max_dist_kb * 1000
    nearby <- sv_chr[
      sv_chr$end >= (c_start_bp - max_dist_bp) &
      sv_chr$start <= (c_end_bp + max_dist_bp), ]

    if (nrow(nearby) == 0) {
      results[[i]] <- make_empty_sv_row(cid)
      next
    }

    # Score each nearby SV
    scores <- numeric(nrow(nearby))
    for (j in seq_len(nrow(nearby))) {
      sv_start <- nearby$start[j]
      sv_end   <- nearby$end[j]
      sv_width <- sv_end - sv_start

      # Overlap
      ov_start <- max(c_start_bp, sv_start)
      ov_end   <- min(c_end_bp, sv_end)
      overlap  <- max(0, ov_end - ov_start)

      # Overlap as fraction of candidate
      overlap_pct <- overlap / c_width_bp * 100

      # Span ratio: how well do the widths match?
      span_ratio <- if (c_width_bp > 0) sv_width / c_width_bp else 0

      # Distance of breakpoints to candidate edges
      left_dist  <- abs(sv_start - c_start_bp) / 1000  # kb
      right_dist <- abs(sv_end - c_end_bp) / 1000      # kb

      # Combined score: high overlap + close breakpoints + similar span
      scores[j] <- overlap_pct -
                   (left_dist + right_dist) * 0.1 +
                   (1 - abs(1 - span_ratio)) * 20
    }

    # Best match
    best <- which.max(scores)
    sv_start <- nearby$start[best]
    sv_end   <- nearby$end[best]

    ov_start <- max(c_start_bp, sv_start)
    ov_end   <- min(c_end_bp, sv_end)
    overlap_pct <- max(0, ov_end - ov_start) / c_width_bp * 100

    results[[i]] <- data.frame(
      candidate_id    = cid,
      best_sv_id      = nearby$sv_id[best],
      best_sv_caller  = if ("caller" %in% names(nearby)) nearby$caller[best] else "unknown",
      sv_type         = if ("sv_type" %in% names(nearby)) nearby$sv_type[best] else "unknown",
      sv_overlap_pct  = round(overlap_pct, 1),
      sv_left_dist_kb = round(abs(sv_start - c_start_bp) / 1000, 1),
      sv_right_dist_kb = round(abs(sv_end - c_end_bp) / 1000, 1),
      sv_span_ratio   = round(if (c_width_bp > 0) (sv_end - sv_start) / c_width_bp else 0, 3),
      n_sv_hits       = nrow(nearby),
      stringsAsFactors = FALSE
    )
  }

  return(do.call(rbind, results))
}


# ---- Helper: Empty SV row ----
make_empty_sv_row <- function(cid) {
  data.frame(
    candidate_id     = cid,
    best_sv_id       = NA_character_,
    best_sv_caller   = NA_character_,
    sv_type          = NA_character_,
    sv_overlap_pct   = NA_real_,
    sv_left_dist_kb  = NA_real_,
    sv_right_dist_kb = NA_real_,
    sv_span_ratio    = NA_real_,
    n_sv_hits        = 0L,
    stringsAsFactors = FALSE
  )
}
