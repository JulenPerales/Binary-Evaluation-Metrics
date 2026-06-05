#' BinaryClassifier R6 Class
#'
#' An object-oriented interface to the full binary classifier evaluation
#' workflow. Wraps a `toc_data` tibble and exposes metric computation,
#' visualisation, and integrated analysis through a coherent API.
#'
#' @section Construction:
#' ```r
#' # From a score + observed vector (most common use case)
#' clf <- BinaryClassifier$new(scores, observed)
#'
#' # From pre-computed hit counts (e.g. from a TOC table)
#' clf <- BinaryClassifier$from_hits(hits_vector, n_positive, n_negative)
#' ```
#'
#' @section Key methods:
#' * `$at_threshold(k)` — confusion matrix at a specific rank threshold
#' * `$compute_metric_curve(metric_fn)` — any metric across all thresholds
#' * `$plot_toc()` — TOC curve
#' * `$plot_roc()` — ROC curve
#' * `$plot_fom()` — FOM curve with bounds
#' * `$plot_metric(metric_fn, ylab)` — any metric curve with bounds
#' * `$plot_all_metrics()` — patchwork panel of all standard metrics
#'
#' @section Key active bindings:
#' * `$toc` — the underlying `toc_data` tibble
#' * `$prevalence` — P / N
#' * `$mfom` — maximum FOM across all thresholds
#' * `$afom` — integrated FOM
#' * `$dfom` — FOM skill vs. random baseline
#' * `$auc` — area under ROC curve
#' * `$dauc` — AUC skill score
#'
#' @export
#' @examples
#' set.seed(42)
#' n     <- 1000
#' score <- runif(n)
#' obs   <- score > 0.5 + rnorm(n, 0, 0.25)
#' clf   <- BinaryClassifier$new(score, obs)
#' clf
#' clf$summary()
#' clf$plot_toc()
#' clf$plot_fom()
BinaryClassifier <- R6::R6Class(
  "BinaryClassifier",
  cloneable = TRUE,

  # ── Private state ──────────────────────────────────────────────────────
  private = list(
    .toc  = NULL,
    .name = NULL
  ),

  # ── Active bindings ────────────────────────────────────────────────────
  active = list(
    toc        = function() private$.toc,
    name       = function(v) {
      if (missing(v)) private$.name else { private$.name <- v; invisible(self) }
    },
    n_positive = function() attr(private$.toc, "n_positive"),
    n_negative = function() attr(private$.toc, "n_negative"),
    n_total    = function() attr(private$.toc, "n_total"),
    prevalence = function() self$n_positive / self$n_total,
    mfom       = function() mfom(private$.toc),
    afom       = function() afom(private$.toc),
    aufom      = function() aufom(private$.toc),
    dfom       = function() dfom(private$.toc),
    auc        = function() auc_roc(private$.toc),
    dauc       = function() dauc(private$.toc),
    aucprc     = function() auc_prc(private$.toc)
  ),

  # ── Public methods ─────────────────────────────────────────────────────
  public = list(

    #' @description
    #' Create a `BinaryClassifier` from score and observed vectors, or from a
    #' pre-built `toc_data` object.
    #' @param scores Numeric vector of classifier scores, **or** a `toc_data`
    #'   object (in which case `observed` should be omitted).
    #' @param observed Logical / integer (0/1) vector of true classes.
    #' @param name Optional character label for the classifier.
    initialize = function(scores, observed = NULL, name = "Classifier") {
      if (inherits(scores, "toc_data")) {
        private$.toc <- scores
      } else {
        stopifnot(!is.null(observed))
        private$.toc <- toc_from_scores(scores, observed)
      }
      private$.name <- name
      invisible(self)
    },

    #' @description
    #' Evaluate the confusion matrix at a single threshold rank `k`.
    #' @param k Integer. Number of predicted positives (0 to N).
    #' @return A `confusion_matrix` object.
    at_threshold = function(k) {
      td <- private$.toc
      stopifnot(k >= 0, k <= self$n_total)
      row <- td[td$k == k, , drop = FALSE]
      confusion_matrix(row$hits, row$fa, row$misses, row$cr)
    },

    #' @description
    #' Compute any metric across all decision thresholds.
    #' @param metric_fn A metric function (e.g. [oa()], [gss()]).
    #' @param metric_name Optional label.
    #' @return A tibble from [toc_metric_curve()].
    compute_metric_curve = function(metric_fn,
                                    metric_name = deparse(substitute(metric_fn))) {
      toc_metric_curve(private$.toc, metric_fn, metric_name)
    },

    #' @description Print a concise summary.
    print = function(...) {
      cat(sprintf(
        "BinaryClassifier <%s>\n  N = %d | P = %d | Q = %d | prev = %.4f\n",
        private$.name, self$n_total, self$n_positive,
        self$n_negative, self$prevalence
      ))
      invisible(self)
    },

    #' @description Full metric summary printed to the console.
    summary = function() {
      self$print()
      cat("\n  ── Point Metrics at Optimal Threshold ──────────────────\n")
      cm <- self$at_threshold(mfom_threshold(private$.toc))
      metrics_vec <- all_metrics(cm)
      nms <- names(metrics_vec)
      for (i in seq_along(metrics_vec)) {
        cat(sprintf("  %-12s  %7.4f\n", nms[i], metrics_vec[i]))
      }
      cat("\n  ── Integrated Metrics ───────────────────────────────────\n")
      int <- integrated_metrics(private$.toc)
      for (i in seq_along(int)) {
        cat(sprintf("  %-12s  %7.4f\n", names(int)[i], int[i]))
      }
      invisible(self)
    },

    # ── Visualisation methods ──────────────────────────────────────────

    #' @description TOC curve (Hits vs k).
    #' @param ... Additional arguments passed to [plot_toc()].
    plot_toc = function(...) {
      plot_toc(private$.toc, classifier_name = private$.name, ...)
    },

    #' @description ROC curve (TPR vs FPR).
    #' @param ... Additional arguments passed to [plot_roc()].
    plot_roc = function(...) {
      plot_roc(private$.toc, classifier_name = private$.name, ...)
    },

    #' @description FOM curve with bounds.
    #' @param ... Additional arguments passed to [plot_fom()].
    plot_fom = function(...) {
      plot_fom(private$.toc, classifier_name = private$.name, ...)
    },

    #' @description Any metric curve with bounds.
    #' @param metric_fn A metric function from this package.
    #' @param ylab Character y-axis label.
    #' @param ... Passed to [plot_metric_curve()].
    plot_metric = function(metric_fn,
                           ylab = deparse(substitute(metric_fn)),
                           ...) {
      mc <- self$compute_metric_curve(metric_fn, ylab)
      plot_metric_curve(mc, ylab = ylab, ...)
    },

    #' @description Grid of all standard metric curves.
    #' @return A `patchwork` composite plot.
    plot_all_metrics = function() {
      td <- private$.toc
      P  <- self$n_positive
      Q  <- self$n_negative
      N  <- self$n_total

      p_oa   <- plot_metric_curve(toc_metric_curve(td, oa,    "OA"),    "OA")
      p_ba   <- plot_metric_curve(toc_metric_curve(td, ba,    "BA"),    "BA")
      p_mcc  <- plot_metric_curve(toc_metric_curve(td, mcc,   "MCC"),   "MCC")
      p_kap  <- plot_metric_curve(toc_metric_curve(td, kappa_score, "Kappa"), "Kappa")
      p_fom  <- plot_metric_curve(toc_metric_curve(td, fom,   "FOM"),   "FOM / CSI")
      p_gss  <- plot_metric_curve(toc_metric_curve(td, gss,   "GSS"),   "GSS")
      p_pss  <- plot_metric_curve(toc_metric_curve(td, pss,   "PSS"),   "PSS")
      p_hss  <- plot_metric_curve(toc_metric_curve(td, hss,   "HSS"),   "HSS")

      patchwork::wrap_plots(
        p_oa, p_ba, p_mcc, p_kap,
        p_fom, p_gss, p_pss, p_hss,
        ncol = 4
      ) +
        patchwork::plot_annotation(
          title    = sprintf("Metric profiles — %s", private$.name),
          subtitle = sprintf(
            "N = %d | P = %d | Q = %d | prev = %.3f",
            N, P, Q, P / N
          )
        )
    }
  )
)

# ── Alternative constructor from raw hit vector ─────────────────────────────

#' Create a BinaryClassifier from Pre-computed Hit Counts
#'
#' Wraps [toc_data()] and [BinaryClassifier] for the case where you already
#' have the cumulative hit counts at each threshold (e.g., from published TOC
#' tables or simulation output).
#'
#' @param hits_vector Numeric vector (length N+1). See [toc_data()].
#' @param n_positive Integer. Total positives P.
#' @param n_negative Integer. Total negatives Q.
#' @param name Optional character label.
#' @return A `BinaryClassifier` object.
#' @export
#' @examples
#' P <- 5; Q <- 10
#' H <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5)
#' clf <- BinaryClassifier_from_hits(H, P, Q, name = "My Classifier")
#' clf$summary()
BinaryClassifier_from_hits <- function(hits_vector, n_positive, n_negative,
                                        name = "Classifier") {
  td <- toc_data(hits_vector, n_positive, n_negative)
  BinaryClassifier$new(td, name = name)
}
