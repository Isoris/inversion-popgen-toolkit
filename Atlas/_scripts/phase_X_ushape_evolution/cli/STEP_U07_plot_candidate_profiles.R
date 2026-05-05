#!/usr/bin/env Rscript
# =============================================================================
# cli/STEP_U07_plot_candidate_profiles.R
# =============================================================================
# One PDF per candidate showing dXY / FST / pi_HOMO1 / pi_HOMO2 across the
# inversion + flanks, with zone shading and class label in the title.
# Plus a combined multipage PDF "all_candidates.pdf".
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(jsonlite); library(ggplot2); library(data.table)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--ushape_json", type = "character"),
  make_option("--out_dir",     type = "character")
)))

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
J <- fromJSON(opts$ushape_json, simplifyVector = FALSE)

zone_fill <- c(left_flank  = "#F2F2F2", left_edge   = "#FDE0DD",
               center      = "#FFFFFF", right_edge  = "#FDE0DD",
               right_flank = "#F2F2F2")
zone_order <- c("left_flank","left_edge","center","right_edge","right_flank")

plot_candidate <- function(cand) {
  if (length(cand$windows) == 0L) return(NULL)
  W <- rbindlist(lapply(cand$windows, as.data.table), fill = TRUE)
  W[, zone := factor(zone, levels = zone_order)]

  # zone backdrop polygons
  zb <- W[, .(start = min(start_bp), end = max(end_bp)), by = zone]
  zb[, fill := zone_fill[as.character(zone)]]

  sc <- cand$shape_scores
  cls <- cand$classification
  subt <- sprintf("U=%s  inside/flank=%s  internal_peak=%s  oldness=%s",
                  format(round(sc$dxy_u_score %||% NA, 2), nsmall = 2),
                  format(round(sc$dxy_inside_flank_ratio %||% NA, 2), nsmall = 2),
                  format(round(sc$dxy_internal_peak_score %||% NA, 2), nsmall = 2),
                  format(round(sc$oldness_score %||% NA, 2), nsmall = 2))

  long <- melt(W,
               id.vars = c("start_bp","end_bp","zone","relative_position"),
               measure.vars = c("dxy_homo1_homo2","fst_homo1_homo2",
                                "pi_homo1","pi_homo2"),
               variable.name = "metric", value.name = "value")
  long[, metric := factor(metric,
       levels = c("dxy_homo1_homo2","fst_homo1_homo2","pi_homo1","pi_homo2"))]

  ggplot() +
    geom_rect(data = zb, aes(xmin = start, xmax = end,
                              ymin = -Inf, ymax = Inf, fill = zone),
              alpha = 0.6, inherit.aes = FALSE) +
    scale_fill_manual(values = zone_fill, guide = "none") +
    geom_line(data = long, aes(x = (start_bp + end_bp)/2,
                                y = value, colour = metric)) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1L) +
    labs(title = sprintf("%s — %s (%s)",
                          cand$candidate_id, cls$primary_class, cls$confidence),
         subtitle = subt, x = "position (bp)", y = NULL, colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

pdfs <- character(0)
for (cand in J$candidates) {
  p <- plot_candidate(cand)
  if (is.null(p)) next
  out_pdf <- file.path(opts$out_dir, paste0(cand$candidate_id, ".pdf"))
  out_png <- file.path(opts$out_dir, paste0(cand$candidate_id, ".png"))
  ggsave(out_pdf, p, width = 9, height = 7)
  ggsave(out_png, p, width = 9, height = 7, dpi = 150)
  pdfs <- c(pdfs, out_pdf)
}

# combined multipage
if (length(pdfs) > 0L && requireNamespace("qpdf", quietly = TRUE)) {
  qpdf::pdf_combine(pdfs, output = file.path(opts$out_dir, "all_candidates.pdf"))
  message("[U07] combined multipage written")
} else {
  message("[U07] qpdf not installed; skipping multipage combine")
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
