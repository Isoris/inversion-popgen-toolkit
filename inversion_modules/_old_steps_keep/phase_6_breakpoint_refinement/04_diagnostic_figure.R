#!/usr/bin/env Rscript

# =============================================================================
# 04_diagnostic_figure.R  (v1.1 — registry-wired)
#
# MULTI-TRACK DIAGNOSTIC FIGURE per candidate, stacked on common x-axis
# (genomic position in Mb).
#
# TRACKS (top to bottom):
#
#   T1  Smoothed |delta| = |h2_mean - h1_mean| along x
#       Shows where the population-level inversion signal is strongest.
#
#   T2  Per-carrier ancestral-fragment left-boundary rug + density
#   T3  Per-carrier ancestral-fragment right-boundary rug + density
#       Per-sample votes + KDE. The peak of each density IS the population
#       breakpoint estimate on that side. Bootstrap CI shown as shaded band.
#
#   T4  C01j regime state ribbon (if available)
#   T5  C01l Delta_12 per segment (if available)
#   T6  STEP41 switch pile-up rug (if available)
#   T7  STEP37 SV clusters (if available, low alpha)
#
# Final consensus breakpoints drawn as red vertical lines with CI ribbons.
#
# REGISTRY WIRING (chat-18):
#   - Input data read from registry blocks:
#       dosage_blocks, ancestral_fragments_summary, boundary_refined_{side}
#     plus per-sample TSV from raw/ (fragments), per-marker TSV from raw/.
#   - Output PDF/PNG goes to evidence_registry/per_candidate/<cid>/figures/
#   - Path registered via reg$evidence$add_evidence(cid, "figure_breakpoint_diagnostic_path", <path>)
#
# Usage:
#   Rscript 04_diagnostic_figure.R [cid=all]
# =============================================================================

# --- Entry: registry + ancestry stack --------------------------------------
Sys.setenv(CURRENT_SCRIPT = "04_diagnostic_figure.R")
.bridge <- Sys.getenv("REGISTRY_BRIDGE", "utils/registry_bridge.R")
if (!file.exists(.bridge)) {
  for (p in c("utils/registry_bridge.R",
              "../utils/registry_bridge.R",
              file.path(Sys.getenv("BASE", ""), "utils/registry_bridge.R"))) {
    if (file.exists(p)) { .bridge <- p; break }
  }
}
if (!file.exists(.bridge)) {
  stop("[04_diagnostic_figure] cannot locate utils/registry_bridge.R.")
}
source(.bridge)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
config_file <- ""
cid_filter  <- NA_character_
.ai <- 1L
while (.ai <= length(args)) {
  a <- args[.ai]
  if (a == "--config" && .ai < length(args)) {
    config_file <- args[.ai + 1L]; .ai <- .ai + 2L
  } else if (a != "all") {
    cid_filter <- a; .ai <- .ai + 1L
  } else {
    .ai <- .ai + 1L
  }
}
if (nzchar(config_file) && file.exists(config_file)) source(config_file)

# Optional upstream disk-based sources (same as script 03)
if (!exists("C01J_DIR"))   C01J_DIR   <- Sys.getenv("C01J_DIR",   NA_character_)
if (!exists("C01L_DIR"))   C01L_DIR   <- Sys.getenv("C01L_DIR",   NA_character_)
if (!exists("STEP37_DIR")) STEP37_DIR <- Sys.getenv("STEP37_DIR", NA_character_)
for (nm in c("C01J_DIR","C01L_DIR","STEP37_DIR")) {
  v <- get(nm)
  if (is.na(v) || !nzchar(v)) assign(nm, NULL)
}

if (!exists("ensure_dir"))
  ensure_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); invisible(p) }

# =============================================================================
# Per-candidate render
# =============================================================================
render_candidate <- function(cid, chr, c_start, c_end) {
  # Output location lives under evidence_registry/<cid>/figures/
  cand_root <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                          "evidence_registry", "per_candidate", cid)
  plot_dir  <- file.path(cand_root, "figures")
  ensure_dir(plot_dir)

  # --- Required inputs via registry blocks + raw/ TSVs ----------------------
  blocks_blk <- tryCatch(reg$evidence$read_block(cid, "dosage_blocks"),
                         error = function(e) NULL)
  blocks <- if (!is.null(blocks_blk)) blocks_blk$data else NULL

  frag_blk <- tryCatch(reg$evidence$read_block(cid, "ancestral_fragments_summary"),
                       error = function(e) NULL)
  frag_summary <- if (!is.null(frag_blk)) frag_blk$data else NULL

  ref_left  <- tryCatch(reg$evidence$read_block(cid, "boundary_refined_left"),
                        error = function(e) NULL)
  ref_right <- tryCatch(reg$evidence$read_block(cid, "boundary_refined_right"),
                        error = function(e) NULL)

  # Per-marker TSV from raw/
  markers_tsv <- if (!is.null(blocks) && !is.null(blocks$markers_tsv_path))
                   blocks$markers_tsv_path
                 else file.path(cand_root, "raw", "dosage_informative_markers.tsv.gz")
  if (!file.exists(markers_tsv)) {
    message("[04_diagnostic_figure] cid=", cid,
            " SKIP — no markers TSV at ", markers_tsv)
    return(invisible(NULL))
  }
  markers <- fread(markers_tsv)

  # Per-sample fragments TSV from raw/
  fragments_tsv <- if (!is.null(frag_summary) && !is.null(frag_summary$fragments_tsv_path))
                     frag_summary$fragments_tsv_path
                   else file.path(cand_root, "raw", "ancestral_fragments_per_sample.tsv.gz")
  frag <- if (file.exists(fragments_tsv)) tryCatch(fread(fragments_tsv),
                                                    error = function(e) NULL) else NULL

  # Consensus breakpoints from refined blocks
  final_left    <- if (!is.null(ref_left$data$final_bp))   as.numeric(ref_left$data$final_bp)  else NA_real_
  final_right   <- if (!is.null(ref_right$data$final_bp))  as.numeric(ref_right$data$final_bp) else NA_real_
  left_ci_low   <- if (!is.null(ref_left$data$ci_low))     as.numeric(ref_left$data$ci_low)    else NA_real_
  left_ci_high  <- if (!is.null(ref_left$data$ci_high))    as.numeric(ref_left$data$ci_high)   else NA_real_
  right_ci_low  <- if (!is.null(ref_right$data$ci_low))    as.numeric(ref_right$data$ci_low)   else NA_real_
  right_ci_high <- if (!is.null(ref_right$data$ci_high))   as.numeric(ref_right$data$ci_high)  else NA_real_

  # Determine x-range from markers
  xrange_bp <- range(markers$pos, na.rm = TRUE)
  xrange_mb <- xrange_bp / 1e6

  # -----------------------------------------------------------------------
  # T1: smoothed |delta| track
  # -----------------------------------------------------------------------
  markers[, pos_mb := pos / 1e6]
  # Smooth with a moving-average kernel across ~50 kb physical bp
  order_idx <- order(markers$pos)
  m_sorted <- markers[order_idx]
  bin_bp <- 50000L
  m_sorted[, bin := (pos %/% bin_bp) * bin_bp]
  smoothed <- m_sorted[, .(pos_mb = mean(pos_mb),
                            smoothed_abs_delta = mean(abs_delta, na.rm = TRUE)),
                        by = bin]

  add_bp_lines <- function(p) {
    if (is.finite(final_left)) {
      p <- p + geom_vline(xintercept = final_left / 1e6,
                           linetype = "solid", color = "#c00000", size = 0.6)
      if (is.finite(left_ci_low) && is.finite(left_ci_high)) {
        p <- p + annotate("rect",
                            xmin = left_ci_low / 1e6, xmax = left_ci_high / 1e6,
                            ymin = -Inf, ymax = Inf,
                            fill = "#c00000", alpha = 0.15)
      }
    }
    if (is.finite(final_right)) {
      p <- p + geom_vline(xintercept = final_right / 1e6,
                           linetype = "solid", color = "#c00000", size = 0.6)
      if (is.finite(right_ci_low) && is.finite(right_ci_high)) {
        p <- p + annotate("rect",
                            xmin = right_ci_low / 1e6, xmax = right_ci_high / 1e6,
                            ymin = -Inf, ymax = Inf,
                            fill = "#c00000", alpha = 0.15)
      }
    }
    p
  }

  p_T1 <- ggplot(smoothed, aes(x = pos_mb, y = smoothed_abs_delta)) +
    geom_area(fill = "#4a6fa5", color = "#233e5c", size = 0.3, alpha = 0.7) +
    labs(title = paste0("Candidate ", cid, "  ", chr,
                          "   input shelf: ",
                          round(c_start / 1e6, 3), " - ",
                          round(c_end / 1e6, 3), " Mb"),
         subtitle = "T1: smoothed |dosage delta|  =  |h2_mean - h1_mean|  (50 kb bins)",
         y = "|delta|") +
    scale_x_continuous(limits = xrange_mb, expand = expansion(0.01)) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, color = "#444"),
          axis.title.x = element_blank())
  p_T1 <- add_bp_lines(p_T1)

  # -----------------------------------------------------------------------
  # T2: left-boundary rug + density
  # -----------------------------------------------------------------------
  p_T2 <- NULL
  if (!is.null(frag) && nrow(frag) > 5) {
    frag_dt <- copy(frag)
    frag_dt[, frag_start_mb := frag_start_bp / 1e6]
    frag_dt[, frag_end_mb   := frag_end_bp   / 1e6]

    p_T2 <- ggplot(frag_dt, aes(x = frag_start_mb)) +
      geom_density(fill = "#2a9d8f", color = "#1a4a3e", size = 0.4, alpha = 0.6) +
      geom_rug(sides = "b", alpha = 0.5, length = unit(0.08, "npc"), color = "#1a4a3e") +
      labs(subtitle = paste0("T2: per-carrier LEFT fragment boundaries  (n=", nrow(frag_dt),
                              ")  rug + KDE density")) +
      scale_x_continuous(limits = xrange_mb, expand = expansion(0.01)) +
      theme_classic(base_size = 9) +
      theme(plot.subtitle = element_text(size = 8, color = "#444"),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_T2 <- add_bp_lines(p_T2)
  }

  # -----------------------------------------------------------------------
  # T3: right-boundary rug + density
  # -----------------------------------------------------------------------
  p_T3 <- NULL
  if (!is.null(frag) && nrow(frag) > 5) {
    p_T3 <- ggplot(frag, aes(x = frag_end_bp / 1e6)) +
      geom_density(fill = "#e76f51", color = "#b05437", size = 0.4, alpha = 0.6) +
      geom_rug(sides = "b", alpha = 0.5, length = unit(0.08, "npc"), color = "#b05437") +
      labs(subtitle = paste0("T3: per-carrier RIGHT fragment boundaries  (n=", nrow(frag),
                              ")  rug + KDE density")) +
      scale_x_continuous(limits = xrange_mb, expand = expansion(0.01)) +
      theme_classic(base_size = 9) +
      theme(plot.subtitle = element_text(size = 8, color = "#444"),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_T3 <- add_bp_lines(p_T3)
  }

  # -----------------------------------------------------------------------
  # T4: regime ribbon (C01j)
  # -----------------------------------------------------------------------
  p_T4 <- NULL
  if (!is.null(C01J_DIR)) {
    f_candidates <- c(
      file.path(C01J_DIR, paste0("regime_segments_", chr, "_cid", cid, ".tsv")),
      file.path(C01J_DIR, paste0("regime_segments_", chr, ".tsv")),
      file.path(C01J_DIR, "regime_segments.tsv.gz"),
      file.path(C01J_DIR, "regime_segments.tsv")
    )
    f <- f_candidates[file.exists(f_candidates)][1]
    if (!is.na(f)) {
      seg <- tryCatch(fread(f), error = function(e) NULL)
      if (!is.null(seg) && nrow(seg) > 0 && "state" %in% names(seg)) {
        if ("chrom" %in% names(seg)) seg <- seg[chrom == chr]
        # Coerce start/end to mb if needed
        if ("start_mb" %in% names(seg) && "end_mb" %in% names(seg)) {
          seg_plot <- seg[start_mb >= xrange_mb[1] & end_mb <= xrange_mb[2]]
        } else if ("start_bp" %in% names(seg) && "end_bp" %in% names(seg)) {
          seg_plot <- seg[start_bp >= xrange_bp[1] & end_bp <= xrange_bp[2]]
          seg_plot[, start_mb := start_bp / 1e6]
          seg_plot[, end_mb   := end_bp   / 1e6]
        } else seg_plot <- data.table()

        if (nrow(seg_plot) > 0) {
          state_colors <- c(
            "clean_inversion"     = "#005f73",
            "structured_moderate" = "#0a9396",
            "structured_complex"  = "#94d2bd",
            "weak_signal"         = "#e9d8a6",
            "transition"          = "#ee9b00",
            "background_soup"     = "#bb3e03"
          )
          p_T4 <- ggplot(seg_plot) +
            geom_rect(aes(xmin = start_mb, xmax = end_mb,
                          ymin = 0, ymax = 1, fill = state), color = NA) +
            scale_fill_manual(values = state_colors, drop = FALSE, name = NULL) +
            scale_x_continuous(limits = xrange_mb, expand = expansion(0.01)) +
            scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            labs(subtitle = "T4: C01j regime state ribbon") +
            theme_classic(base_size = 9) +
            theme(plot.subtitle = element_text(size = 8, color = "#444"),
                  axis.title = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.position = "bottom",
                  legend.key.size = unit(0.3, "cm"),
                  legend.text = element_text(size = 6))
          p_T4 <- add_bp_lines(p_T4)
        }
      }
    }
  }

  # -----------------------------------------------------------------------
  # T5: C01l Delta_12 bars per segment (if available)
  # -----------------------------------------------------------------------
  p_T5 <- NULL
  if (!is.null(C01L_DIR)) {
    f_candidates <- c(
      file.path(C01L_DIR, "segment_summary.tsv.gz"),
      file.path(C01L_DIR, "segment_summary.tsv")
    )
    f <- f_candidates[file.exists(f_candidates)][1]
    if (!is.na(f)) {
      seg <- tryCatch(fread(f), error = function(e) NULL)
      if (!is.null(seg) && nrow(seg) > 0) {
        if ("candidate_id" %in% names(seg)) seg <- seg[candidate_id == cid]
        delta_col <- if ("delta_12" %in% names(seg)) "delta_12"
                      else if ("Delta_12" %in% names(seg)) "Delta_12" else NA
        if (nrow(seg) > 0 && !is.na(delta_col) && "segment" %in% names(seg)) {
          # Approximate x position per segment using c_start/c_end
          # (left_flank: before c_start; inv_core: middle; right_flank: after c_end)
          seg_levels <- c("left_flank", "inv_left_half", "inv_core",
                           "inv_right_half", "right_flank")
          seg[, segment := factor(segment, levels = seg_levels)]
          # Use segment centers on the x axis
          inv_span <- c_end - c_start
          seg_xcenters <- data.table(
            segment = seg_levels,
            x_mb = c(c_start - inv_span * 0.35, c_start + inv_span / 6,
                       c_start + inv_span * 0.5, c_end - inv_span / 6,
                       c_end + inv_span * 0.35) / 1e6
          )
          seg <- merge(seg, seg_xcenters, by = "segment", all.x = TRUE)

          p_T5 <- ggplot(seg, aes(x = x_mb, y = get(delta_col))) +
            geom_col(fill = "#6a4c93", color = "#4a2c7a", width = inv_span / 6e6) +
            geom_text(aes(label = segment), y = 0, vjust = 1.2, size = 2.5, color = "#444") +
            scale_x_continuous(limits = xrange_mb, expand = expansion(0.01)) +
            labs(subtitle = paste0("T5: C01l Delta_12 per segment (high flanks = validated boundary)"),
                 y = "Delta_12") +
            theme_classic(base_size = 9) +
            theme(plot.subtitle = element_text(size = 8, color = "#444"),
                  axis.title.x = element_blank())
          p_T5 <- add_bp_lines(p_T5)
        }
      }
    }
  }

  # -----------------------------------------------------------------------
  # T6: STEP41 switch pile-ups (if available via $STEP41_DIR)
  # -----------------------------------------------------------------------
  p_T6 <- NULL
  step41_dir <- Sys.getenv("STEP41_DIR", "")
  if (nzchar(step41_dir)) {
    sw_file_candidates <- c(
      file.path(step41_dir, cid, "candidate_switching_events.tsv.gz"),
      file.path(step41_dir, cid, "candidate_switching_events.tsv"),
      file.path(step41_dir, paste0(chr, ".candidate_", cid),
                "candidate_switching_events.tsv.gz"),
      file.path(step41_dir, paste0(chr, ".candidate_", cid),
                "candidate_switching_events.tsv")
    )
    sw_file <- sw_file_candidates[file.exists(sw_file_candidates)][1]
    if (!is.na(sw_file)) {
      sw <- tryCatch(fread(sw_file), error = function(e) NULL)
      if (!is.null(sw) && nrow(sw) >= 3 && "switch_bp" %in% names(sw)) {
        p_T6 <- ggplot(sw, aes(x = switch_bp / 1e6)) +
          geom_density(fill = "#8d99ae", color = "#2b2d42", size = 0.3, alpha = 0.5) +
          geom_rug(sides = "b", alpha = 0.4, length = unit(0.05, "npc"), color = "#2b2d42") +
          labs(subtitle = paste0("T6: STEP41 sample switching events (pile-ups = transition zones, n=",
                                  nrow(sw), ")")) +
          scale_x_continuous(limits = xrange_mb, expand = expansion(0.01)) +
          theme_classic(base_size = 9) +
          theme(plot.subtitle = element_text(size = 8, color = "#444"),
                axis.title.x = element_blank(), axis.title.y = element_blank(),
                axis.text.y = element_blank(), axis.ticks.y = element_blank())
        p_T6 <- add_bp_lines(p_T6)
      }
    }
  }

  # -----------------------------------------------------------------------
  # T7: STEP37 SV clusters (bottom layer, low alpha, user: "if you have energy")
  # -----------------------------------------------------------------------
  p_T7 <- NULL
  if (!is.null(STEP37_DIR)) {
    f_candidates <- c(
      file.path(STEP37_DIR, "breakpoint_support_per_candidate.tsv"),
      file.path(STEP37_DIR, paste0("bp_clusters_", chr, ".tsv"))
    )
    f <- f_candidates[file.exists(f_candidates)][1]
    if (!is.na(f)) {
      sv <- tryCatch(fread(f), error = function(e) NULL)
      if (!is.null(sv) && nrow(sv) > 0) {
        if ("candidate_id" %in% names(sv)) sv <- sv[candidate_id == cid]
        if (nrow(sv) > 0 && "cluster_median" %in% names(sv)) {
          sv[, pos_mb := cluster_median / 1e6]
          side_col <- if ("boundary_side" %in% names(sv)) "boundary_side" else "side"
          # jitter y slightly
          sv[, y := runif(.N, 0.3, 0.7)]
          p_T7 <- ggplot(sv, aes(x = pos_mb, y = y)) +
            geom_point(alpha = 0.35, size = 2, color = "#8a2b30") +
            scale_x_continuous(limits = xrange_mb, expand = expansion(0.01)) +
            scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            labs(subtitle = paste0("T7: STEP37 SV breakpoint clusters [BOTTOM LAYER — weight 0.5 — DELLY/Manta shown for reference only]"),
                 x = paste0("position on ", chr, " (Mb)")) +
            theme_classic(base_size = 9) +
            theme(plot.subtitle = element_text(size = 8, color = "#8a2b30",
                                                face = "italic"),
                  axis.title.y = element_blank(),
                  axis.text.y = element_blank(), axis.ticks.y = element_blank())
          p_T7 <- add_bp_lines(p_T7)
        }
      }
    }
  }

  # -----------------------------------------------------------------------
  # Assemble
  # -----------------------------------------------------------------------
  panels <- list(T1 = p_T1, T2 = p_T2, T3 = p_T3, T4 = p_T4,
                  T5 = p_T5, T6 = p_T6, T7 = p_T7)
  panels <- panels[!sapply(panels, is.null)]
  if (length(panels) == 0) {
    message("[04_diagnostic_figure] cid=", cid, " SKIP — no panels to render")
    return(invisible(NULL))
  }
  bottom_name <- names(panels)[length(panels)]
  if (bottom_name != "T7") {
    panels[[bottom_name]] <- panels[[bottom_name]] +
      labs(x = paste0("position on ", chr, " (Mb)")) +
      theme(axis.title.x = element_text())
  }

  if (has_patchwork && length(panels) > 1) {
    combined <- Reduce(`/`, panels) +
      patchwork::plot_layout(heights = c(3, rep(1.5, length(panels) - 1)))
  } else {
    combined <- panels[[1]]
  }

  pdf_path <- file.path(plot_dir, "candidate_breakpoint_diagnostic.pdf")
  png_path <- file.path(plot_dir, "candidate_breakpoint_diagnostic.png")
  total_h <- max(7, 1.4 + 1.8 * length(panels))
  total_w <- 11

  tryCatch(ggsave(pdf_path, combined, width = total_w, height = total_h,
                   device = cairo_pdf),
           error = function(e) ggsave(pdf_path, combined, width = total_w, height = total_h))
  tryCatch(ggsave(png_path, combined, width = total_w, height = total_h,
                   dpi = 150),
           error = function(e) message("[04_diagnostic_figure]   PNG failed: ",
                                        conditionMessage(e)))

  # Register figure paths as flat evidence keys so downstream gallery /
  # atlas scripts can discover them via reg$evidence$get_keys().
  tryCatch({
    reg$evidence$add_evidence(cid, "figure_breakpoint_diagnostic_pdf_path", pdf_path)
    reg$evidence$add_evidence(cid, "figure_breakpoint_diagnostic_png_path", png_path)
  }, error = function(e) {
    message("[04_diagnostic_figure]   could not register figure paths: ",
            conditionMessage(e))
  })

  message("[04_diagnostic_figure] cid=", cid, " -> ", basename(pdf_path),
          " (", length(panels), " tracks)")
  invisible(NULL)
}

# =============================================================================
# Driver
# =============================================================================
main <- function() {
  cand_path <- file.path(BRIDGE_PATHS$REGISTRIES_ROOT, "data",
                          "interval_registry", "candidate_intervals.tsv")
  if (!file.exists(cand_path)) {
    stop("[04_diagnostic_figure] candidate_intervals.tsv not found at ", cand_path)
  }
  cand <- fread(cand_path)
  if (!is.na(cid_filter)) cand <- cand[candidate_id == cid_filter]
  if (nrow(cand) == 0) {
    message("[04_diagnostic_figure] No candidates to process",
            if (!is.na(cid_filter)) paste0(" (filter: ", cid_filter, ")") else "")
    return(invisible(NULL))
  }
  message("[04_diagnostic_figure] Rendering ", nrow(cand), " candidate(s)")
  for (ci in seq_len(nrow(cand))) {
    row <- cand[ci]
    tryCatch(
      render_candidate(as.character(row$candidate_id),
                        as.character(row$chrom),
                        as.numeric(row$start_bp),
                        as.numeric(row$end_bp)),
      error = function(e) {
        message("[04_diagnostic_figure] cid=", row$candidate_id,
                " ERROR: ", conditionMessage(e))
      }
    )
  }
  message("[04_diagnostic_figure] DONE")
}

main()

