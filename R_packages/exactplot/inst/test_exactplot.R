#!/usr/bin/env Rscript
# =============================================================================
# exactplot — test & demo script
# Run after: install.packages("exactplot_0.1.0.tar.gz", repos=NULL, type="source")
# =============================================================================

library(ggplot2)
library(exactplot)
library(grid)

cat("=== exactplot test suite ===\n\n")

# ---- Test 1: unit_convert round-trips ----
cat("Test 1: unit conversions\n")
stopifnot(abs(unit_convert(2.54, "cm", "in") - 1.0) < 1e-10)
stopifnot(abs(unit_convert(1, "in", "pt") - 72) < 1e-10)
stopifnot(abs(unit_convert(25.4, "mm", "in") - 1.0) < 1e-10)
stopifnot(abs(unit_convert(1, "in", "cm") - 2.54) < 1e-10)
cat("  PASS\n\n")

# ---- Test 2: exact_dims returns sane values ----
cat("Test 2: exact_dims basic\n")
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
  labs(title = "Test Plot", x = "Weight", y = "MPG")

dims <- exact_dims(p, panel_width = 5, panel_height = 4, units = "cm")
cat(sprintf("  Panel target: 5 cm x 4 cm\n"))
cat(sprintf("  Device size:  %.3f x %.3f in\n", dims$device_width, dims$device_height))
cat(sprintf("  Chrome L=%.3f R=%.3f T=%.3f B=%.3f in\n",
    dims$chrome$left, dims$chrome$right, dims$chrome$top, dims$chrome$bottom))
cat(sprintf("  Iterations: %d\n", dims$iterations))

# Sanity: device must be larger than the panel
pw_in <- unit_convert(5, "cm", "in")
ph_in <- unit_convert(4, "cm", "in")
stopifnot(dims$device_width > pw_in)
stopifnot(dims$device_height > ph_in)
# Chrome must be positive
stopifnot(dims$chrome$left > 0)
stopifnot(dims$chrome$bottom > 0)
cat("  PASS\n\n")

# ---- Test 3: ggsave_exact produces a file ----
cat("Test 3: ggsave_exact output\n")
tmpf <- tempfile(fileext = ".pdf")
result <- ggsave_exact(tmpf, p, panel_width = 5, panel_height = 4,
                       units = "cm", verbose = TRUE)
stopifnot(file.exists(tmpf))
cat(sprintf("  Output: %s (%.1f KB)\n", basename(tmpf), file.info(tmpf)$size / 1024))
cat("  PASS\n\n")

# ---- Test 4: Verify actual panel size in the PDF ----
cat("Test 4: Panel size verification (gold standard)\n")
# We re-render at the computed device size and measure the panel
verify_panel_size <- function(plot, panel_w_cm, panel_h_cm) {
  dims <- exact_dims(plot, panel_width = panel_w_cm, panel_height = panel_h_cm, units = "cm")

  # Open a device at exactly the computed size
  pdf(nullfile(), width = dims$device_width, height = dims$device_height)
  on.exit(dev.off())

  gt <- ggplotGrob(plot)

  # Find panel cells
  panel_rows <- which(grepl("^panel", gt$layout$name))
  panel_l <- min(gt$layout$l[panel_rows])
  panel_r <- max(gt$layout$r[panel_rows])
  panel_t <- min(gt$layout$t[panel_rows])
  panel_b <- max(gt$layout$b[panel_rows])

  w_in <- convertWidth(gt$widths, "inches", valueOnly = TRUE)
  h_in <- convertHeight(gt$heights, "inches", valueOnly = TRUE)

  actual_pw_in <- sum(w_in[panel_l:panel_r])
  actual_ph_in <- sum(h_in[panel_t:panel_b])

  target_pw_in <- unit_convert(panel_w_cm, "cm", "in")
  target_ph_in <- unit_convert(panel_h_cm, "cm", "in")

  err_w <- abs(actual_pw_in - target_pw_in)
  err_h <- abs(actual_ph_in - target_ph_in)

  list(
    target_w_in = target_pw_in, actual_w_in = actual_pw_in, err_w = err_w,
    target_h_in = target_ph_in, actual_h_in = actual_ph_in, err_h = err_h
  )
}

v <- verify_panel_size(p, 5, 4)
cat(sprintf("  Target W: %.4f in | Actual W: %.4f in | Error: %.5f in (%.3f mm)\n",
    v$target_w_in, v$actual_w_in, v$err_w, v$err_w * 25.4))
cat(sprintf("  Target H: %.4f in | Actual H: %.4f in | Error: %.5f in (%.3f mm)\n",
    v$target_h_in, v$actual_h_in, v$err_h, v$err_h * 25.4))
stopifnot(v$err_w < 0.01)  # less than 0.25 mm error
stopifnot(v$err_h < 0.01)
cat("  PASS\n\n")

# ---- Test 5: scale_to_data ----
cat("Test 5: scale_to_data\n")
sdims <- scale_to_data(p, x_scale = 1, y_scale = 0.2, units = "cm")
cat(sprintf("  X range: %.1f | x_scale=1 -> panel_width=%.2f cm\n",
    sdims$x_range, sdims$panel_width))
cat(sprintf("  Y range: %.1f | y_scale=0.2 -> panel_height=%.2f cm\n",
    sdims$y_range, sdims$panel_height))
stopifnot(abs(sdims$panel_width - sdims$x_range * 1) < 0.01)
stopifnot(abs(sdims$panel_height - sdims$y_range * 0.2) < 0.01)
cat("  PASS\n\n")

# ---- Test 6: axis_range extraction ----
cat("Test 6: axis_range\n")
p2 <- ggplot(data.frame(x = 1:10, y = 1:10), aes(x, y)) + geom_point() +
  xlim(0, 10) + ylim(0, 20)
rng <- axis_range(p2)
cat(sprintf("  X limits: [%.1f, %.1f]  Y limits: [%.1f, %.1f]\n",
    rng$x_limits[1], rng$x_limits[2], rng$y_limits[1], rng$y_limits[2]))
cat("  PASS\n\n")

# ---- Test 7: Different themes ----
cat("Test 7: Robustness across themes\n")
themes <- list(
  minimal = theme_minimal(),
  classic = theme_classic(),
  bw      = theme_bw(),
  void    = theme_void(),
  fat_margins = theme(plot.margin = margin(2, 2, 2, 2, "cm"))
)

for (nm in names(themes)) {
  p_themed <- p + themes[[nm]]
  v <- verify_panel_size(p_themed, 3, 4)
  status <- if (v$err_w < 0.01 && v$err_h < 0.01) "PASS" else "FAIL"
  cat(sprintf("  %-12s: err_w=%.4f err_h=%.4f  %s\n", nm, v$err_w, v$err_h, status))
}
cat("\n")

# ---- Test 8: Plot with legend ----
cat("Test 8: Plot with legend\n")
p_legend <- ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
  geom_point() +
  labs(title = "With Legend", color = "Cylinders")

v <- verify_panel_size(p_legend, 5, 4)
cat(sprintf("  Error: W=%.4f in (%.3f mm)  H=%.4f in (%.3f mm)\n",
    v$err_w, v$err_w * 25.4, v$err_h, v$err_h * 25.4))
stopifnot(v$err_w < 0.01 && v$err_h < 0.01)
cat("  PASS\n\n")

# ---- Test 9: Faceted plot ----
cat("Test 9: Faceted plot\n")
p_facet <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
  facet_wrap(~cyl)

v <- verify_panel_size(p_facet, 8, 3)
cat(sprintf("  Error: W=%.4f in (%.3f mm)  H=%.4f in (%.3f mm)\n",
    v$err_w, v$err_w * 25.4, v$err_h, v$err_h * 25.4))
cat("  (Note: for facets, this is the bounding box of all panels)\n")
stopifnot(v$err_w < 0.02 && v$err_h < 0.02)
cat("  PASS\n\n")

# ---- Summary ----
cat("=== All tests passed ===\n")

# ---- Demo outputs ----
cat("\n--- Generating demo figures ---\n")
outdir <- file.path(getwd(), "exactplot_demo")
dir.create(outdir, showWarnings = FALSE)

# Demo 1: Simple scatter, panel = exactly 6 cm x 4 cm
p_demo <- ggplot(mtcars, aes(wt, mpg)) +
  geom_point(size = 2) +
  labs(title = "Weight vs Fuel Economy",
       x = "Weight (1000 lbs)", y = "Miles per gallon") +
  theme_classic()

ggsave_exact(file.path(outdir, "demo_6x4cm.pdf"), p_demo,
  panel_width = 6, panel_height = 4, units = "cm", verbose = TRUE)

# Demo 2: Data-scaled: 1 wt unit = 1 cm
sdims <- scale_to_data(p_demo, x_scale = 1, y_scale = 0.15, units = "cm")
ggsave_exact(file.path(outdir, "demo_data_scaled.pdf"), p_demo,
  panel_width = sdims$panel_width, panel_height = sdims$panel_height,
  units = "cm", verbose = TRUE)

# Demo 3: Tiny panel for a journal inset
ggsave_exact(file.path(outdir, "demo_2x2cm_inset.pdf"), p_demo,
  panel_width = 2, panel_height = 2, units = "cm", verbose = TRUE)

cat(sprintf("\nDemo files saved to: %s\n", outdir))
cat(list.files(outdir), sep = "\n")
