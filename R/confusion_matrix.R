#' Create a Confusion Matrix Object
#'
#' Constructs a `confusion_matrix` object from the four cells of a binary
#' classification contingency table, following the meteorological convention:
#' **Hits** (TP), **False Alarms** (FP), **Misses** (FN), and **Correct
#' Rejections** (TN).
#'
#' @param hits Non-negative number. True positives (Hits / TP).
#' @param fa Non-negative number. False positives (False Alarms / FP).
#' @param misses Non-negative number. False negatives (Misses / FN).
#' @param cr Non-negative number. True negatives (Correct Rejections / TN).
#'
#' @return An S3 object of class `confusion_matrix` with elements:
#'   `hits`, `fa`, `misses`, `cr`, `n_positive`, `n_negative`, `n_total`,
#'   `prevalence`.
#'
#' @export
#' @examples
#' cm <- confusion_matrix(hits = 28, fa = 72, misses = 23, cr = 2680)
#' print(cm)
#' summary(cm)
confusion_matrix <- function(hits, fa, misses, cr) {
  stopifnot(
    is.numeric(hits),   is.numeric(fa),
    is.numeric(misses), is.numeric(cr),
    hits >= 0, fa >= 0, misses >= 0, cr >= 0
  )
  n_positive <- hits + misses
  n_negative <- fa + cr
  n_total    <- n_positive + n_negative
  if (n_total == 0) stop("Confusion matrix total must be > 0.")
  structure(
    list(
      hits       = hits,
      fa         = fa,
      misses     = misses,
      cr         = cr,
      n_positive = n_positive,
      n_negative = n_negative,
      n_total    = n_total,
      prevalence = n_positive / n_total
    ),
    class = "confusion_matrix"
  )
}

#' @export
print.confusion_matrix <- function(x, ...) {
  cat("  Confusion Matrix\n")
  cat("  ───────────────────────────────────\n")
  cat(sprintf("  %16s  %8s  %8s\n", "", "Obs (+)", "Obs (-)"))
  cat(sprintf("  %16s  %8g  %8g\n", "Predicted (+)", x$hits,   x$fa))
  cat(sprintf("  %16s  %8g  %8g\n", "Predicted (-)", x$misses, x$cr))
  cat("  ───────────────────────────────────\n")
  cat(sprintf("  N = %g  |  P = %g  |  Q = %g  |  prev = %.4f\n",
              x$n_total, x$n_positive, x$n_negative, x$prevalence))
  invisible(x)
}

#' @export
summary.confusion_matrix <- function(object, ...) {
  x <- object
  print(x)
  cat("\n  ── Agreement Metrics ────────────────\n")
  cat(sprintf("  %-8s  %6.4f  Overall Accuracy\n",     "OA",    oa(x)))
  cat(sprintf("  %-8s  %6.4f  Balanced Accuracy\n",    "BA",    ba(x)))
  cat(sprintf("  %-8s  %6.4f  Matthews Correlation\n", "MCC",   mcc(x)))
  cat(sprintf("  %-8s  %6.4f  Cohen's Kappa\n",        "Kappa", kappa_score(x)))
  cat(sprintf("  %-8s  %6.4f  F1 Score\n",             "F1",    f1_score(x)))
  cat("\n  ── Skill Metrics ─────────────────────\n")
  cat(sprintf("  %-8s  %6.4f  Figure of Merit / CSI\n","FOM",   fom(x)))
  cat(sprintf("  %-8s  %6.4f  Gilbert Skill Score\n",  "GSS",   gss(x)))
  cat(sprintf("  %-8s  %6.4f  Peirce Skill Score\n",   "PSS",   pss(x)))
  cat(sprintf("  %-8s  %6.4f  Heidke Skill Score\n",   "HSS",   hss(x)))
  invisible(x)
}

#' @export
as.data.frame.confusion_matrix <- function(x, ...) {
  tibble::tibble(
    hits       = x$hits,
    fa         = x$fa,
    misses     = x$misses,
    cr         = x$cr,
    n_positive = x$n_positive,
    n_negative = x$n_negative,
    n_total    = x$n_total,
    prevalence = x$prevalence
  )
}

#' @export
`[.confusion_matrix` <- function(x, i, ...) {
  unclass(x)[i]
}
