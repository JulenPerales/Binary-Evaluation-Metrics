# ---------------------------------------------------------------------------
# Shared legend data ----
# ---------------------------------------------------------------------------

.curve_scale_colour <- function(classifier_name = "Classifier") {
  list(
    ggplot2::scale_colour_manual(
      name   = "",
      values = c(
        "Classifier"        = .BINM_COLOURS$classifier,
        "Upper Bound"       = .BINM_COLOURS$upper,
        "Lower Bound"       = .BINM_COLOURS$lower,
        "Random Classifier" = .BINM_COLOURS$random
      ),
      breaks = c("Classifier", "Random Classifier", "Upper Bound", "Lower Bound")
    ),
    ggplot2::scale_linetype_manual(
      name   = "",
      values = c(
        "Classifier"        = "solid",
        "Upper Bound"       = "solid",
        "Lower Bound"       = "longdash",
        "Random Classifier" = "dotted"
      ),
      breaks = c("Classifier", "Random Classifier", "Upper Bound", "Lower Bound")
    )
  )
}

# ---------------------------------------------------------------------------
# plot_toc ----
# ---------------------------------------------------------------------------

#' Plot the Threshold Operating Characteristic (TOC) Curve
#'
#' Displays Hits (TP count) versus the number of predicted positives `k`,
#' together with the upper bound, lower bound and random-classifier reference.
#'
#' @param toc A `toc_data` object.
#' @param classifier_name Character. Legend label for the classifier curve.
#' @param percent Logical. Scale axes as percentages of the total (default
#'   `FALSE` — raw counts).
#' @param annotate_auc Logical. Overlay the AUC value (default `TRUE`).
#' @return A `ggplot2` object.
#' @export
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(200), runif(200) > 0.5)
#' plot_toc(td)
plot_toc <- function(toc, classifier_name = "Classifier",
                     percent = FALSE, annotate_auc = TRUE) {
  stopifnot(inherits(toc, "toc_data"))
  P  <- attr(toc, "n_positive")
  N  <- attr(toc, "n_total")
  sc <- if (percent) N else 1
  sy <- if (percent) P else 1

  df <- tibble::tibble(
    k          = toc$k / sc,
    Classifier = toc$hits / sy,
    `Upper Bound`       = toc$hits_upper / sy,
    `Lower Bound`       = toc$hits_lower / sy,
    `Random Classifier` = toc$hits_random / sy
  ) |>
    tidyr::pivot_longer(-k, names_to = "curve", values_to = "hits")

  xlab <- if (percent) "Predicted Positives (fraction of N)" else "Predicted Positives (k)"
  ylab <- if (percent) "Hits (fraction of P)"                else "Hits"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = hits,
                                         colour = curve, linetype = curve)) +
    ggplot2::geom_line(linewidth = 1.1) +
    .curve_scale_colour(classifier_name) +
    ggplot2::labs(
      x        = xlab,
      y        = ylab,
      title    = "TOC Curve",
      subtitle = sprintf("P = %d | Q = %d | prev = %.3f",
                         P, attr(toc, "n_negative"), P / N)
    ) +
    binm_theme()

  if (annotate_auc) {
    auc_val <- auc_roc(toc)
    p <- p + ggplot2::annotate(
      "text",
      x = Inf, y = -Inf,
      label  = sprintf("AUC = %.3f", auc_val),
      hjust  = 1.1, vjust = -0.5,
      size   = 3.5,
      colour = "grey30"
    )
  }
  p
}

# ---------------------------------------------------------------------------
# plot_roc ----
# ---------------------------------------------------------------------------

#' Plot the ROC Curve
#'
#' Receiver Operating Characteristic curve: TPR (Sensitivity) versus FPR
#' (1 - Specificity).
#'
#' @param toc A `toc_data` object.
#' @param classifier_name Character. Legend label.
#' @param annotate_auc Logical. Overlay AUC value.
#' @return A `ggplot2` object.
#' @export
#' @examples
#' set.seed(2)
#' td <- toc_from_scores(runif(400), runif(400) > 0.45)
#' plot_roc(td)
plot_roc <- function(toc, classifier_name = "Classifier", annotate_auc = TRUE) {
  stopifnot(inherits(toc, "toc_data"))

  df_clf <- tibble::tibble(
    fpr   = toc$fpr,
    tpr   = toc$tpr,
    curve = classifier_name
  )
  df_rnd <- tibble::tibble(
    fpr   = c(0, 1),
    tpr   = c(0, 1),
    curve = "Random Classifier"
  )
  df <- dplyr::bind_rows(df_clf, df_rnd)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = fpr, y = tpr,
                                         colour = curve, linetype = curve)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 1),
                                 labels = scales::percent) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                                 labels = scales::percent) +
    ggplot2::scale_colour_manual(
      name   = "",
      values = c(
        setNames(.BINM_COLOURS$classifier, classifier_name),
        "Random Classifier" = .BINM_COLOURS$random
      )
    ) +
    ggplot2::scale_linetype_manual(
      name   = "",
      values = c(setNames("solid", classifier_name),
                 "Random Classifier" = "dotted")
    ) +
    ggplot2::labs(x = "False Positive Rate", y = "True Positive Rate",
                  title = "ROC Curve") +
    binm_theme()

  if (annotate_auc) {
    p <- p + ggplot2::annotate(
      "text",
      x = 0.98, y = 0.02,
      label  = sprintf("AUC = %.3f", auc_roc(toc)),
      hjust  = 1, vjust = 0,
      size   = 3.5, colour = "grey30"
    )
  }
  p
}

# ---------------------------------------------------------------------------
# plot_fom ----
# ---------------------------------------------------------------------------

#' Plot the FOM Curve with Bounds
#'
#' Displays the FOM (CSI) as a function of threshold `k`, together with the
#' upper bound, lower bound and random-classifier baseline.  Annotates the
#' area ratio (AFOM relative to upper/lower area).
#'
#' @param toc A `toc_data` object.
#' @param classifier_name Character. Legend label.
#' @param percent_x Logical. Scale x-axis as percentage of N (default `TRUE`).
#' @return A `ggplot2` object.
#' @export
#' @examples
#' set.seed(3)
#' td <- toc_from_scores(runif(300), runif(300) > 0.5)
#' plot_fom(td)
plot_fom <- function(toc, classifier_name = "Classifier", percent_x = TRUE) {
  stopifnot(inherits(toc, "toc_data"))
  N <- attr(toc, "n_total")
  P <- attr(toc, "n_positive")
  Q <- attr(toc, "n_negative")

  x_vals <- if (percent_x) toc$k / N else toc$k

  df <- tibble::tibble(
    k                   = x_vals,
    Classifier          = toc$fom_curve,
    `Upper Bound`       = toc$fom_upper,
    `Lower Bound`       = toc$fom_lower,
    `Random Classifier` = toc$fom_random
  ) |>
    tidyr::pivot_longer(-k, names_to = "curve", values_to = "fom")

  ar  <- area_ratio(x_vals, toc$fom_curve, toc$fom_upper, toc$fom_lower)
  dau <- dfom(toc)

  xlab <- if (percent_x) "Predicted Positives (% of N)" else "Predicted Positives (k)"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = fom,
                                         colour = curve, linetype = curve)) +
    ggplot2::geom_line(linewidth = 1.1) +
    .curve_scale_colour(classifier_name) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0), limits = c(0, max(x_vals)),
      labels = if (percent_x) scales::percent else scales::label_number()
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0), limits = c(0, 1),
      labels = scales::percent, breaks = seq(0, 1, 0.2)
    ) +
    ggplot2::labs(
      x        = xlab,
      y        = "FOM / CSI",
      title    = "FOM Curve with Bounds",
      subtitle = sprintf(
        "P = %d | Q = %d | prev = %.3f", P, Q, P / N
      )
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = sprintf(
        "AFOM = %.3f\nDFOM = %.3f\nMFOM = %.3f",
        afom(toc), dau, mfom(toc)
      ),
      hjust = 1.05, vjust = 1.2,
      size  = 3.2, colour = "grey30"
    ) +
    binm_theme()
  p
}

# ---------------------------------------------------------------------------
# plot_metric_curve ----
# ---------------------------------------------------------------------------

#' Plot Any Metric Curve with Bounds
#'
#' Generic plotting function that accepts the output of [toc_metric_curve()]
#' and displays the classifier curve alongside upper, lower and random bounds.
#'
#' @param mc A tibble returned by [toc_metric_curve()].
#' @param ylab Character. Y-axis label (defaults to `mc$metric_name[1]`).
#' @param percent_x Logical. Scale x-axis as percentage (default `TRUE`).
#' @return A `ggplot2` object.
#' @export
#' @examples
#' set.seed(4)
#' td  <- toc_from_scores(runif(300), runif(300) > 0.5)
#' mc  <- toc_metric_curve(td, gss, "GSS")
#' plot_metric_curve(mc, "GSS")
plot_metric_curve <- function(mc, ylab = mc$metric_name[1], percent_x = TRUE) {
  df <- tibble::tibble(
    k                   = mc$k_pct,
    Classifier          = mc$curve,
    `Upper Bound`       = mc$upper,
    `Lower Bound`       = mc$lower,
    `Random Classifier` = mc$random
  ) |>
    tidyr::pivot_longer(-k, names_to = "curve", values_to = "value")

  type <- mc$type[1]
  base_col <- if (type == "Skill") "random" else "lower"

  y_base  <- if (type == "Skill") mc$random else mc$lower
  ar <- area_ratio(mc$k_pct, mc$curve, mc$upper, y_base)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = value,
                                         colour = curve, linetype = curve)) +
    ggplot2::geom_line(linewidth = 1.1) +
    .curve_scale_colour("Classifier") +
    ggplot2::scale_x_continuous(
      expand = c(0, 0), limits = c(0, 1),
      labels = scales::percent
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      x        = "Predicted Positives (% of N)",
      y        = ylab,
      title    = sprintf("%s curve with bounds", ylab),
      subtitle = sprintf("Type: %s | Area ratio = %.3f", type, ar)
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label  = sprintf("Area ratio\n= %.3f", ar),
      hjust  = 1.05, vjust = 1.2,
      size   = 3.2, colour = "grey30"
    ) +
    binm_theme()
  p
}

# ---------------------------------------------------------------------------
# plot_precision_recall ----
# ---------------------------------------------------------------------------

#' Plot the Precision-Recall Curve
#'
#' @param toc A `toc_data` object.
#' @param classifier_name Character. Legend label.
#' @return A `ggplot2` object.
#' @export
#' @examples
#' set.seed(5)
#' td <- toc_from_scores(runif(400), runif(400) > 0.55)
#' plot_precision_recall(td)
plot_precision_recall <- function(toc, classifier_name = "Classifier") {
  stopifnot(inherits(toc, "toc_data"))
  P  <- attr(toc, "n_positive")
  N  <- attr(toc, "n_total")
  pr <- P / N

  recall    <- toc$hits / P
  precision <- ifelse(toc$k == 0, 1, toc$hits / toc$k)

  df <- tibble::tibble(
    recall    = recall,
    Classifier = precision,
    Baseline   = pr
  ) |>
    tidyr::pivot_longer(-recall, names_to = "curve", values_to = "precision")

  aucprc <- auc_prc(toc)
  skill  <- prc_skill(toc)

  ggplot2::ggplot(df, ggplot2::aes(x = recall, y = precision,
                                    colour = curve, linetype = curve)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::scale_colour_manual(
      name   = "",
      values = c(
        setNames(.BINM_COLOURS$classifier, classifier_name),
        "Classifier" = .BINM_COLOURS$classifier,
        "Baseline"   = .BINM_COLOURS$random
      )
    ) +
    ggplot2::scale_linetype_manual(
      name   = "",
      values = c(
        setNames("solid", classifier_name),
        "Classifier" = "solid",
        "Baseline"   = "dotted"
      )
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 1),
                                 labels = scales::percent) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                                 labels = scales::percent) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x        = "Recall (TPR)",
      y        = "Precision (PPV)",
      title    = "Precision-Recall Curve"
    ) +
    ggplot2::annotate(
      "text", x = 0.02, y = 0.02,
      label  = sprintf("AUCPRC = %.3f\nPRC Skill = %.3f", aucprc, skill),
      hjust  = 0, vjust = 0,
      size   = 3.2, colour = "grey30"
    ) +
    binm_theme()
}

# ---------------------------------------------------------------------------
# compare_classifiers ----
# ---------------------------------------------------------------------------

#' Compare Multiple Classifiers on a Single Plot
#'
#' Overlays the TOC or ROC curves of several classifiers, each passed as a
#' `toc_data` object.
#'
#' @param toc_list A named list of `toc_data` objects.
#' @param type Character. One of `"toc"`, `"roc"`, or `"fom"`.
#' @return A `ggplot2` object.
#' @export
#' @examples
#' set.seed(10)
#' td_a <- toc_from_scores(runif(300), runif(300) > 0.4)
#' td_b <- toc_from_scores(runif(300), runif(300) > 0.5)
#' compare_classifiers(list(A = td_a, B = td_b), type = "roc")
compare_classifiers <- function(toc_list, type = c("roc", "toc", "fom")) {
  type <- match.arg(type)
  stopifnot(
    is.list(toc_list),
    length(toc_list) >= 2,
    all(sapply(toc_list, inherits, "toc_data"))
  )
  nms <- names(toc_list)
  if (is.null(nms)) nms <- paste0("Clf_", seq_along(toc_list))

  build_df <- function(td, nm) {
    switch(type,
      roc = tibble::tibble(x = td$fpr,          y = td$tpr,        clf = nm),
      toc = tibble::tibble(x = td$k,             y = td$hits,       clf = nm),
      fom = tibble::tibble(x = td$k / attr(td, "n_total"),
                           y = td$fom_curve,     clf = nm)
    )
  }

  df  <- dplyr::bind_rows(mapply(build_df, toc_list, nms, SIMPLIFY = FALSE))
  pal <- scales::hue_pal()(length(nms))

  axis_labels <- switch(type,
    roc = list("False Positive Rate", "True Positive Rate"),
    toc = list("Predicted Positives (k)", "Hits"),
    fom = list("Predicted Positives (% of N)", "FOM / CSI")
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = clf)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_colour_manual(name = "", values = setNames(pal, nms)) +
    ggplot2::labs(x = axis_labels[[1]], y = axis_labels[[2]],
                  title = sprintf("Classifier Comparison — %s", toupper(type))) +
    binm_theme()

  if (type == "roc") {
    p <- p +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                            colour = "grey50") +
      ggplot2::coord_equal() +
      ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1))
  }
  p
}
