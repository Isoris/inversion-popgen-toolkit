#' exactplot: Export ggplot2 Figures with Exact Plot-Area Dimensions
#'
#' @description
#' When you say "3 cm x 4 cm", you mean the **plot area** — the rectangle
#' bounded by the axes where your data lives. Not the full figure with
#' margins, axis labels, titles, and legends tacked on.
#'
#' `exactplot` fixes this by measuring the non-panel "chrome" and inflating
#' the device dimensions so your panel lands at exactly the size you asked for.
#'
#' @section Main functions:
#' \describe{
#'   \item{[ggsave_exact()]}{Drop-in replacement for `ggsave()` — just
#'     specify `panel_width` and `panel_height` instead of `width`/`height`.}
#'   \item{[exact_dims()]}{Compute the required device dimensions without
#'     saving (useful if you want to open a device manually).}
#'   \item{[scale_to_data()]}{Define panel size by data-coordinate scale
#'     (e.g., "1 unit = 0.5 cm") instead of absolute dimensions.}
#'   \item{[axis_range()]}{Extract the effective axis limits from a built plot.}
#'   \item{[unit_convert()]}{Convert between cm, in, mm, pt.}
#' }
#'
#' @docType package
#' @name exactplot-package
"_PACKAGE"
