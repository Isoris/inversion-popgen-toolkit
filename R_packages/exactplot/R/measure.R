#' Unit conversion utilities
#'
#' Convert between physical units used for plot dimensions.
#'
#' @param value Numeric value to convert.
#' @param from Source unit: one of "cm", "in", "mm", "pt".
#' @param to Target unit: one of "cm", "in", "mm", "pt".
#' @return Numeric value in the target unit.
#' @export
#' @examples
#' unit_convert(2.54, "cm", "in")  # 1
#' unit_convert(1, "in", "pt")     # 72
unit_convert <- function(value, from = "cm", to = "in") {
  # Everything goes through inches as the canonical internal unit
  to_inches <- c(cm = 1 / 2.54, mm = 1 / 25.4, "in" = 1, pt = 1 / 72)
  from_inches <- c(cm = 2.54, mm = 25.4, "in" = 1, pt = 72)

  from <- match.arg(from, names(to_inches))
  to   <- match.arg(to, names(from_inches))

  value * to_inches[[from]] * from_inches[[to]]
}


#' Measure non-panel grob sizes from a built ggplot
#'
#' Renders a ggplot into a null device at a trial size, then measures
#' how much space is consumed by everything that is NOT the panel
#' (axis labels, titles, legends, margins, etc.).
#'
#' @param plot A ggplot2 object.
#' @param trial_width Trial device width in inches (used for text wrapping).
#' @param trial_height Trial device height in inches.
#' @return A list with components: left, right, top, bottom (all in inches),
#'   representing the non-panel chrome on each side.
#' @keywords internal
measure_chrome <- function(plot, trial_width = 7, trial_height = 7) {

  # Open a null device so we can build grobs without side-effects
  cur_dev <- grDevices::dev.cur()
  grDevices::pdf(nullfile(), width = trial_width, height = trial_height)
  on.exit({
    grDevices::dev.off()
    if (cur_dev > 1L) try(grDevices::dev.set(cur_dev), silent = TRUE)
  }, add = TRUE)

  # Build the plot into a gtable
  gt <- ggplot2::ggplotGrob(plot)

  # Find the panel cell(s). In single-panel plots there's one; in facets

  # there are multiple. We want the bounding box of ALL panels.
  panel_rows <- which(grepl("^panel", gt$layout$name))

  if (length(panel_rows) == 0L) {
    stop("exactplot: could not find any panel grobs in the plot. ",
         "Is this a valid ggplot2 object?")
  }

  panel_t <- min(gt$layout$t[panel_rows])
  panel_b <- max(gt$layout$b[panel_rows])
  panel_l <- min(gt$layout$l[panel_rows])
  panel_r <- max(gt$layout$r[panel_rows])

  # Convert gtable widths/heights to inches
  w_in <- grid::convertWidth(gt$widths, "inches", valueOnly = TRUE)
  h_in <- grid::convertHeight(gt$heights, "inches", valueOnly = TRUE)

  # Sum the chrome: everything outside the panel cell range
  chrome_left   <- sum(w_in[seq_len(panel_l - 1L)])
  chrome_right  <- sum(w_in[seq(from = panel_r + 1L, to = length(w_in))])
  chrome_top    <- sum(h_in[seq_len(panel_t - 1L)])
  chrome_bottom <- sum(h_in[seq(from = panel_b + 1L, to = length(h_in))])

  list(
    left   = chrome_left,
    right  = chrome_right,
    top    = chrome_top,
    bottom = chrome_bottom,
    panel_width  = sum(w_in[panel_l:panel_r]),
    panel_height = sum(h_in[panel_t:panel_b])
  )
}


#' Compute exact device dimensions for a desired panel size
#'
#' Two-pass algorithm:
#' 1. Render at a trial size to measure chrome.
#' 2. Compute device = panel_target + chrome.
#' 3. Re-render at that device size to re-measure (chrome can shift
#'    slightly when the device size changes due to text wrapping).
#' 4. Iterate up to `max_iter` times until convergence.
#'
#' @param plot A ggplot2 object.
#' @param panel_width Desired panel width (numeric).
#' @param panel_height Desired panel height (numeric).
#' @param units Unit for panel_width/panel_height: "cm", "in", "mm", "pt".
#' @param max_iter Maximum refinement iterations (default 5).
#' @param tol Convergence tolerance in inches (default 0.005 ≈ 0.13 mm).
#' @return A list with device_width, device_height (in inches), and the
#'   final chrome measurements.
#' @export
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' dims <- exact_dims(p, panel_width = 3, panel_height = 4, units = "cm")
#' dims$device_width   # total device width in inches
#' dims$device_height  # total device height in inches
exact_dims <- function(plot,
                       panel_width,
                       panel_height,
                       units = "cm",
                       max_iter = 5L,
                       tol = 0.005) {

  # Convert target to inches
  pw_in <- unit_convert(panel_width, from = units, to = "in")
  ph_in <- unit_convert(panel_height, from = units, to = "in")

  # Initial trial: generous device so text doesn't get clipped
  dev_w <- pw_in + 3

  dev_h <- ph_in + 3

  for (i in seq_len(max_iter)) {
    chrome <- measure_chrome(plot, trial_width = dev_w, trial_height = dev_h)

    new_dev_w <- pw_in + chrome$left + chrome$right
    new_dev_h <- ph_in + chrome$top  + chrome$bottom

    # Check convergence
    if (abs(new_dev_w - dev_w) < tol && abs(new_dev_h - dev_h) < tol) {
      dev_w <- new_dev_w
      dev_h <- new_dev_h
      break
    }

    dev_w <- new_dev_w
    dev_h <- new_dev_h
  }

  list(
    device_width  = dev_w,
    device_height = dev_h,
    panel_width   = pw_in,
    panel_height  = ph_in,
    chrome        = chrome,
    iterations    = i,
    units_requested = units
  )
}
