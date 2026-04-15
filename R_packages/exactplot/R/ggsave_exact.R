#' Save a ggplot with exact plot-area dimensions
#'
#' Drop-in replacement for [ggplot2::ggsave()] where `width` and `height`
#' refer to the **panel area** (the region bounded by the axes), not the
#' full device including margins, axis text, titles, and legends.
#'
#' The function uses an iterative two-pass approach: it renders the plot
#' to measure the "chrome" (everything outside the panel), then inflates
#' the device dimensions so the panel lands at exactly the requested size.
#'
#' @param filename Output file path. Extension determines format
#'   (pdf, png, tiff, svg, etc.).
#' @param plot A ggplot2 object. Defaults to [ggplot2::last_plot()].
#' @param panel_width Desired width of the plot panel area.
#' @param panel_height Desired height of the plot panel area.
#' @param units Units for panel_width and panel_height:
#'   `"cm"` (default), `"in"`, `"mm"`, or `"pt"`.
#' @param dpi Resolution in dots per inch (for raster formats).
#' @param bg Background colour (default `"white"`).
#' @param device Graphics device. If NULL, inferred from filename extension.
#'   Can be a string ("pdf", "png", etc.) or a function.
#' @param max_iter Maximum iterations for the convergence loop.
#' @param tol Convergence tolerance in inches.
#' @param verbose If TRUE, print the computed device dimensions and chrome.
#' @param ... Additional arguments passed to the graphics device.
#' @return Invisibly returns the `exact_dims` result list.
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
#'   labs(title = "Weight vs MPG", x = "Weight (1000 lbs)", y = "Miles/gallon")
#'
#' # Panel area will be exactly 3 cm x 4 cm
#' ggsave_exact("figure1.pdf", p, panel_width = 3, panel_height = 4, units = "cm")
#'
#' # Same but in inches
#' ggsave_exact("figure1.png", p, panel_width = 2, panel_height = 1.5,
#'              units = "in", dpi = 300)
#'
#' # With verbose output to see what happened
#' ggsave_exact("figure1.tiff", p, panel_width = 80, panel_height = 60,
#'              units = "mm", dpi = 600, verbose = TRUE)
#' }
ggsave_exact <- function(filename,
                         plot = ggplot2::last_plot(),
                         panel_width,
                         panel_height,
                         units = "cm",
                         dpi = 300,
                         bg = "white",
                         device = NULL,
                         max_iter = 5L,
                         tol = 0.005,
                         verbose = FALSE,
                         ...) {

  if (!inherits(plot, "ggplot") && !inherits(plot, "gg")) {
    stop("`plot` must be a ggplot2 object.")
  }

  # Compute the exact device dimensions
  dims <- exact_dims(
    plot         = plot,
    panel_width  = panel_width,
    panel_height = panel_height,
    units        = units,
    max_iter     = max_iter,
    tol          = tol
  )

  if (verbose) {
    message(sprintf(
      "exactplot: panel target = %.3f x %.3f in | device = %.3f x %.3f in | iterations = %d",
      dims$panel_width, dims$panel_height,
      dims$device_width, dims$device_height,
      dims$iterations
    ))
    message(sprintf(
      "  chrome: L=%.3f R=%.3f T=%.3f B=%.3f in",
      dims$chrome$left, dims$chrome$right,
      dims$chrome$top, dims$chrome$bottom
    ))
  }

  # Delegate to ggplot2::ggsave with the inflated device dimensions
  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = dims$device_width,
    height   = dims$device_height,
    units    = "in",
    dpi      = dpi,
    bg       = bg,
    device   = device,
    ...
  )

  invisible(dims)
}
