#!/usr/bin/env Rscript
# =============================================================================
# cheat19_diversity_gradient.R — Within-inversion diversity gradient (U-shape)
#
# BIOLOGY:
#   Inside a real inversion, nucleotide diversity follows a U-shape:
#   LOW at breakpoints (recombination strongly suppressed), gradually
#   INCREASING toward center (more exchange possible). The depth of the U
#   correlates with inversion age — older inversions have deeper U-shapes.
#
# INPUT:  per-sample θ_P from MODULE_3, decomposition classes, boundaries
# OUTPUT: U-shape parameters, per-class comparison, gradient plot data
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ── Parameters ──────────────────────────────────────────────────────────
N_BINS_DEFAULT   <- 20L
MIN_SAMPLES_FIT  <- 5L
MIN_BINS_FIT     <- 8L
U_PVAL_THRESHOLD <- 0.05

# ── Compute diversity by position ─────────────────────────────────────

compute_diversity_by_position <- function(theta_z, decomp_classes,
                                           candidate_start, candidate_end,
                                           left_bp, right_bp,
                                           n_bins = N_BINS_DEFAULT) {
  if (nrow(theta_z) == 0 || length(decomp_classes) == 0) return(data.table())

  # Add genotype class
  theta_z[, geno_class := decomp_classes[sample_id]]
  theta_z <- theta_z[!is.na(geno_class)]
  if (nrow(theta_z) == 0) return(data.table())

  # Relative position: 0 = left bp, 1 = right bp
  span <- right_bp - left_bp
  if (span <= 0) return(data.table())

  theta_z[, rel_pos := (wincenter - left_bp) / span]
  theta_z <- theta_z[rel_pos >= -0.1 & rel_pos <= 1.1]  # allow small padding

  # Bin by position
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  theta_z[, bin := findInterval(rel_pos, bin_breaks, rightmost.closed = TRUE)]
  theta_z[bin == 0, bin := 1L]
  theta_z[bin > n_bins, bin := n_bins]
  theta_z[, bin_center := (bin_breaks[bin] + bin_breaks[pmin(bin + 1L, n_bins + 1L)]) / 2]

  # Per-bin per-class mean diversity
  result <- theta_z[, .(
    mean_diversity = mean(tP_z, na.rm = TRUE),
    sd_diversity = sd(tP_z, na.rm = TRUE),
    n_samples = uniqueN(sample_id),
    n_windows = .N
  ), by = .(geno_class, bin, bin_center)]

  theta_z[, c("geno_class", "rel_pos", "bin", "bin_center") := NULL]
  result
}

# ── Fit U-shape quadratic ─────────────────────────────────────────────

fit_u_shape <- function(div_by_pos, genotype_class = "HOM_INV") {
  empty <- list(u_coefficient = NA_real_, u_pvalue = NA_real_,
                u_depth = NA_real_, is_u_shaped = FALSE,
                interpretation = "insufficient_data",
                model = NULL)

  class_dt <- div_by_pos[geno_class == genotype_class]
  if (nrow(class_dt) < MIN_BINS_FIT) return(empty)

  x <- class_dt$bin_center
  y <- class_dt$mean_diversity

  # Model: diversity ~ a * x * (1 - x) + c
  # This is equivalent to: y ~ I(x * (1 - x))
  # a > 0 → center elevated → U-shape confirmed
  quad_term <- x * (1 - x)

  fit_quad <- tryCatch(lm(y ~ quad_term), error = function(e) NULL)
  fit_null <- tryCatch(lm(y ~ 1), error = function(e) NULL)

  if (is.null(fit_quad) || is.null(fit_null)) return(empty)

  a_coeff <- coef(fit_quad)["quad_term"]
  a_se    <- tryCatch(summary(fit_quad)$coefficients["quad_term", "Std. Error"],
                       error = function(e) NA_real_)

  # F-test: quadratic vs null
  f_test <- tryCatch(anova(fit_null, fit_quad), error = function(e) NULL)
  p_val  <- if (!is.null(f_test)) f_test$`Pr(>F)`[2] else NA_real_

  # U-depth: diversity at center (x=0.5) minus mean at breakpoints (x~0 and x~1)
  pred_center <- predict(fit_quad, newdata = data.frame(quad_term = 0.25))
  pred_edge   <- predict(fit_quad, newdata = data.frame(quad_term = 0))
  u_depth     <- pred_center - pred_edge

  is_u <- !is.na(a_coeff) && a_coeff > 0 &&
           !is.na(p_val) && p_val < U_PVAL_THRESHOLD

  interpretation <- if (is.na(a_coeff)) "insufficient_data"
    else if (a_coeff > 0 && !is.na(p_val) && p_val < 0.01) "strong_U_old_inversion"
    else if (a_coeff > 0 && !is.na(p_val) && p_val < 0.05) "weak_U_young_inversion"
    else if (abs(a_coeff) < a_se) "flat_very_young_or_false"
    else if (a_coeff < 0) "inverted_U_suspicious"
    else "ambiguous"

  list(u_coefficient = round(a_coeff, 4),
       u_pvalue = p_val,
       u_depth = round(u_depth, 4),
       is_u_shaped = is_u,
       interpretation = interpretation,
       model = fit_quad)
}

# ── Compare U-shape across genotype classes ───────────────────────────

compare_classes <- function(div_by_pos) {
  classes <- unique(div_by_pos$geno_class)
  results <- list()
  for (cl in classes) {
    ushape <- fit_u_shape(div_by_pos, genotype_class = cl)
    results[[cl]] <- data.table(
      geno_class = cl,
      u_coefficient = ushape$u_coefficient,
      u_pvalue = ushape$u_pvalue,
      u_depth = ushape$u_depth,
      is_u_shaped = ushape$is_u_shaped,
      interpretation = ushape$interpretation)
  }
  if (length(results) > 0) rbindlist(results) else data.table()
}

# ── Search mode ────────────────────────────────────────────────────────

search_diversity_gradient <- function(chr, zone_start, zone_end,
                                       theta_z = NULL,
                                       decomp_classes = NULL,
                                       left_bp = NULL, right_bp = NULL,
                                       ...) {
  empty <- data.table(method = "diversity_gradient", best_bp = NA_integer_,
                       score = 0, is_precise = FALSE, detail = "no_data")
  if (is.null(theta_z) || nrow(theta_z) == 0) return(empty)
  if (is.null(left_bp)) left_bp <- zone_start
  if (is.null(right_bp)) right_bp <- zone_end

  div <- compute_diversity_by_position(theta_z, decomp_classes,
                                        zone_start, zone_end,
                                        left_bp, right_bp)
  if (nrow(div) == 0) return(empty)

  ushape <- fit_u_shape(div, "HOM_INV")
  sc <- if (ushape$is_u_shaped) pmin(1, abs(ushape$u_coefficient) * 2) else 0

  data.table(method = "diversity_gradient",
             best_bp = as.integer((zone_start + zone_end) / 2),
             score = round(sc, 3), is_precise = FALSE,
             detail = paste0("u_coeff=", ushape$u_coefficient,
                              ",p=", round(ushape$u_pvalue, 4)))
}

# ── Convenience runner ─────────────────────────────────────────────────

run_cheat19 <- function(chr, candidate_start, candidate_end,
                         left_bp, right_bp, theta_dir,
                         sample_ids, decomp_classes) {
  message("[cheat19] ", chr, ":", round(candidate_start/1e6,1), "-",
          round(candidate_end/1e6,1), " Mb")

  # Load theta z-scores (reuse cheat12 infrastructure if available)
  theta_z <- if (exists("compute_theta_zscores", mode = "function")) {
    compute_theta_zscores(theta_dir, sample_ids, chr,
                           candidate_start, candidate_end)
  } else {
    message("[cheat19] compute_theta_zscores not found — source cheat12 first")
    data.table()
  }

  if (nrow(theta_z) == 0) {
    message("[cheat19] No theta data")
    return(list(div_by_pos = data.table(), u_shape = list(),
                class_comparison = data.table(),
                search_result = data.table(method = "diversity_gradient",
                  best_bp = NA_integer_, score = 0, is_precise = FALSE,
                  detail = "no_data")))
  }

  message("[cheat19] Loaded ", uniqueN(theta_z$sample_id), " samples × ",
          uniqueN(theta_z$wincenter), " windows")

  div <- compute_diversity_by_position(theta_z, decomp_classes,
                                        candidate_start, candidate_end,
                                        left_bp, right_bp)

  u_inv <- fit_u_shape(div, "HOM_INV")
  message("[cheat19] HOM_INV U-shape: coeff=", u_inv$u_coefficient,
          " p=", round(u_inv$u_pvalue, 4),
          " depth=", u_inv$u_depth,
          " → ", u_inv$interpretation)

  class_comp <- compare_classes(div)
  if (nrow(class_comp) > 0)
    message("[cheat19] Class comparison:\n",
            paste(capture.output(print(class_comp)), collapse = "\n"))

  list(div_by_pos = div, u_shape = u_inv,
       class_comparison = class_comp,
       search_result = search_diversity_gradient(chr, candidate_start,
                        candidate_end, theta_z, decomp_classes,
                        left_bp, right_bp))
}
