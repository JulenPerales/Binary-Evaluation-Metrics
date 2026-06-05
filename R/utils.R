# ---------------------------------------------------------------------------
# Package ggplot2 theme ----
# ---------------------------------------------------------------------------

#' binmetrics ggplot2 Theme
#'
#' A clean, publication-ready ggplot2 theme used as the default for all
#' **binmetrics** plots.
#'
#' @param base_size Base font size (default 12).
#' @param base_family Base font family.
#' @return A `ggplot2::theme` object.
#' @export
#' @examples
#' ggplot2::ggplot() + binm_theme()
binm_theme <- function(base_size = 12, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.background    = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background   = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border       = ggplot2::element_rect(colour = "black", fill = NA),
      panel.grid.major   = ggplot2::element_line(colour = "grey90", linetype = "dashed"),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text          = ggplot2::element_text(size = ggplot2::rel(0.85)),
      axis.title         = ggplot2::element_text(face = "bold"),
      legend.background  = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key         = ggplot2::element_rect(fill = "white", colour = NA),
      legend.title       = ggplot2::element_text(face = "bold"),
      legend.position    = "right",
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle      = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      strip.background   = ggplot2::element_rect(fill = "grey95"),
      strip.text         = ggplot2::element_text(face = "bold")
    )
}

# ---------------------------------------------------------------------------
# Colour palette ----
# ---------------------------------------------------------------------------

.BINM_COLOURS <- list(
  classifier = "#298d8d",
  upper      = "grey20",
  lower      = "grey60",
  random     = "#950000",
  tipping    = c("darkgreen", "grey85", "darkred")
)

# ---------------------------------------------------------------------------
# Package-level documentation ----
# ---------------------------------------------------------------------------

#' binmetrics: Binary Classification Evaluation Metrics
#'
#' Tools for computing, comparing, and visualising binary classification
#' evaluation metrics across the full range of decision thresholds.
#'
#' @section Main workflows:
#' 1. **Quick metric at a single threshold**: Use [confusion_matrix()] +
#'    individual metric functions ([fom()], [oa()], [gss()], …) or
#'    [all_metrics()].
#' 2. **Full threshold analysis**: Build a [toc_data()] (or use
#'    [toc_from_scores()]) and call [toc_metric_curve()], [afom()],
#'    [auc_roc()], etc.
#' 3. **Object-oriented API**: Create a [BinaryClassifier] and call its
#'    `$summary()`, `$plot_toc()`, `$plot_fom()`, `$plot_all_metrics()`
#'    methods.
#' 4. **Interactive exploration**: [launch_threshold_explorer()] and
#'    [launch_bounds_explorer()] open Shiny applications.
#'
#' @docType package
#' @name binmetrics-package
#' @aliases binmetrics
"_PACKAGE"
