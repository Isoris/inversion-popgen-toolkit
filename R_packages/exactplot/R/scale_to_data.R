#' Define panel size by axis range and physical scale
#'
#' Instead of saying "3 cm wide", say "1 unit of X = 0.5 cm" and let the
#' axis range determine the panel width automatically. This is useful for
#' maps, karyoplots, or any figure where the data-to-physical ratio matters.
#'
#' @param plot A ggplot2 object.
#' @param x_scale Physical size per data unit on the x-axis (e.g., 0.5 means
#'   1 data unit = 0.5 cm). Set to NULL to leave width unconstrained and use
#'   `panel_width` instead.
#' @param y_scale Physical size per data unit on the y-axis. Set to NULL to
#'   leave height unconstrained and use `panel_height` instead.
#' @param panel_width Fixed panel width (used when x_scale is NULL).
#' @param panel_height Fixed panel height (used when y_scale is NULL).
#' @param units Physical units: "cm", "in", "mm", "pt".
#' @return A list with `panel_width` and `panel_height` in the given units.
#' @export
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#'
#' # 1 unit of weight = 1 cm, 1 unit of mpg = 0.2 cm
#' dims <- scale_to_data(p, x_scale = 1, y_scale = 0.2, units = "cm")
#' # Now use: ggsave_exact("fig.pdf", p,
#' #   panel_width = dims$panel_width, panel_height = dims$panel_height)
scale_to_data <- function(plot,
                          x_scale = NULL,
                          y_scale = NULL,
                          panel_width = NULL,
                          panel_height = NULL,
                          units = "cm") {

  rng <- axis_range(plot)

  if (!is.null(x_scale)) {
    pw <- rng$x_range * x_scale
  } else if (!is.null(panel_width)) {
    pw <- panel_width
  } else {
    stop("Provide either `x_scale` or `panel_width`.")
  }

  if (!is.null(y_scale)) {
    ph <- rng$y_range * y_scale
  } else if (!is.null(panel_height)) {
    ph <- panel_height
  } else {
    stop("Provide either `y_scale` or `panel_height`.")
  }

  list(
    panel_width  = pw,
    panel_height = ph,
    units        = units,
    x_range      = rng$x_range,
    y_range      = rng$y_range,
    x_limits     = rng$x_limits,
    y_limits     = rng$y_limits
  )
}


#' Extract the effective axis ranges from a ggplot
#'
#' Builds the plot to determine the final axis limits after scale
#' transformations, expansions, and coord adjustments.
#'
#' @param plot A ggplot2 object.
#' @return A list with x_limits (length-2 numeric), y_limits,
#'   x_range, and y_range.
#' @export
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' axis_range(p)
axis_range <- function(plot) {

  # Build the plot to resolve scales
  built <- ggplot2::ggplot_build(plot)

  # For standard Cartesian coords, panel_params holds the resolved ranges
  pp <- built$layout$panel_params

  if (length(pp) == 0L) {
    stop("exactplot: could not extract panel parameters. ",
         "Is this a valid ggplot2 object with data?")
  }

  # Use the first panel (for facets, they should share ranges if
  # scales are fixed, which is the common case for exact sizing)
  p1 <- pp[[1]]

  # ggplot2 >= 3.5 uses $x$range, older uses $x.range
  x_lim <- tryCatch(
    p1$x$limits %||% p1$x.range,
    error = function(e) NULL
  )
  y_lim <- tryCatch(
    p1$y$limits %||% p1$y.range,
    error = function(e) NULL
  )

  # Fallback: try continuous_range

  if (is.null(x_lim)) {
    x_lim <- tryCatch(
      p1$x$continuous_range %||% p1$x.range,
      error = function(e) NULL
    )
  }
  if (is.null(y_lim)) {
    y_lim <- tryCatch(
      p1$y$continuous_range %||% p1$y.range,
      error = function(e) NULL
    )
  }

  # Final attempt via panel_scales
  if (is.null(x_lim)) {
    x_lim <- tryCatch({
      sc <- built$layout$panel_scales_x[[1]]
      sc$range$range
    }, error = function(e) {
      stop("exactplot: could not determine x-axis range.")
    })
  }
  if (is.null(y_lim)) {
    y_lim <- tryCatch({
      sc <- built$layout$panel_scales_y[[1]]
      sc$range$range
    }, error = function(e) {
      stop("exactplot: could not determine y-axis range.")
    })
  }

  list(
    x_limits = x_lim,
    y_limits = y_lim,
    x_range  = diff(x_lim),
    y_range  = diff(y_lim)
  )
}

# Internal null-coalescing operator (same as rlang's %||%)
`%||%` <- function(a, b) if (!is.null(a)) a else b
