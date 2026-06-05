# ---------------------------------------------------------------------------
# TOC bounds (hit counts) ----
# ---------------------------------------------------------------------------

#' Upper Bound Hits at Threshold k
#'
#' For a given threshold rank `k` (number of predicted positives), the
#' theoretically maximum achievable number of Hits is `min(k, P)`, reached
#' when the highest-ranked cases are all true positives.
#'
#' @param k Integer vector. Number of predicted positives (0 to N).
#' @param n_positive Integer. Total number of positives (P).
#' @return Numeric vector of upper-bound Hit counts.
#' @export
#' @seealso [hits_lower()], [hits_random()]
#' @examples
#' hits_upper(0:10, n_positive = 5)
hits_upper <- function(k, n_positive) {
  pmin(k, n_positive)
}

#' Lower Bound Hits at Threshold k
#'
#' The theoretical minimum number of Hits at threshold `k` is
#' `max(0, k - Q)` where `Q = n_negative`. This occurs when all true
#' negatives are ranked first and negatives spill into the positive predictions
#' only when forced.
#'
#' @param k Integer vector. Number of predicted positives (0 to N).
#' @param n_positive Integer. Total number of positives (P).
#' @param n_negative Integer. Total number of negatives (Q).
#' @return Numeric vector of lower-bound Hit counts.
#' @export
#' @examples
#' hits_lower(0:10, n_positive = 5, n_negative = 5)
hits_lower <- function(k, n_positive, n_negative) {
  pmax(0, k - n_negative)
}

#' Random Classifier Hits at Threshold k
#'
#' A classifier with no discrimination selects cases uniformly at random,
#' yielding `k * prevalence` expected Hits.
#'
#' @param k Integer vector. Number of predicted positives (0 to N).
#' @param prevalence Numeric scalar in (0, 1). Fraction of positives.
#' @return Numeric vector of expected Hit counts for a random classifier.
#' @export
#' @examples
#' hits_random(0:10, prevalence = 0.3)
hits_random <- function(k, prevalence) {
  k * prevalence
}

# ---------------------------------------------------------------------------
# FOM bounds (metric values) ----
# ---------------------------------------------------------------------------

#' Upper Bound FOM at Threshold k
#'
#' The maximum achievable FOM at threshold `k` is achieved when all predicted
#' positives that can be Hits are Hits:
#' - `k <= P`: `FOM_upper = k / P`
#' - `k >  P`: `FOM_upper = P / k`
#'
#' Equivalently: `min(k, P) / (P + max(0, k - P))`.
#'
#' @param k Integer vector. Number of predicted positives.
#' @param n_positive Integer. Total positives (P).
#' @param n_negative Integer. Total negatives (Q). Not used in this bound but
#'   included for interface consistency.
#' @return Numeric vector of upper-bound FOM values in \[0, 1\].
#' @export
#' @examples
#' fom_upper(0:100, n_positive = 30, n_negative = 70)
fom_upper <- function(k, n_positive, n_negative) {
  H  <- hits_upper(k, n_positive)
  FA <- pmax(0, k - n_positive)
  M  <- n_positive - H
  denom <- H + FA + M
  ifelse(denom == 0, 0, H / denom)
}

#' Lower Bound FOM at Threshold k
#'
#' The minimum achievable FOM at threshold `k` is:
#' - `k <= Q`: `FOM_lower = 0`
#' - `k >  Q`: `FOM_lower = (k - Q) / N`
#'
#' @param k Integer vector. Number of predicted positives.
#' @param n_positive Integer. Total positives (P).
#' @param n_negative Integer. Total negatives (Q).
#' @return Numeric vector of lower-bound FOM values in \[0, 1\].
#' @export
#' @examples
#' fom_lower(0:100, n_positive = 30, n_negative = 70)
fom_lower <- function(k, n_positive, n_negative) {
  H  <- hits_lower(k, n_positive, n_negative)
  FA <- k - H
  M  <- n_positive - H
  denom <- H + FA + M
  ifelse(denom == 0, 0, H / denom)
}

#' Random Classifier FOM at Threshold k
#'
#' FOM of a classifier with no discrimination (`AUC = 0.5`), which assigns
#' probabilities uniformly.
#' `FOM_random = (k * prevalence) / (k + P*(1 - prevalence))
#'             = k*P/N / (k + P - k*P/N)`
#'
#' @param k Integer vector. Number of predicted positives.
#' @param n_positive Integer. Total positives (P).
#' @param n_negative Integer. Total negatives (Q).
#' @return Numeric vector of random-classifier FOM values.
#' @export
#' @examples
#' fom_random(0:100, n_positive = 30, n_negative = 70)
fom_random <- function(k, n_positive, n_negative) {
  N <- n_positive + n_negative
  H <- hits_random(k, n_positive / N)   # k * prevalence
  FA <- k - H
  M  <- n_positive - H
  denom <- H + FA + M
  ifelse(denom == 0, 0, H / denom)
}

# ---------------------------------------------------------------------------
# Apply any metric function to bound hit counts ----
# ---------------------------------------------------------------------------

#' Evaluate a Metric Along a Bound Curve
#'
#' Given a vector of threshold values `k` and a bound hit-count function,
#' computes the metric value at each threshold.  Useful for overlaying any
#' metric (OA, BA, GSS, …) alongside its bounds.
#'
#' @param k Integer vector. Threshold values (0 to N).
#' @param hits_fn Function `f(k, n_positive, n_negative)` returning Hit counts.
#' @param metric_fn Function from this package: `fom`, `oa`, `gss`, etc.
#'   Must accept a `confusion_matrix` object.
#' @param n_positive Integer. Total positives (P).
#' @param n_negative Integer. Total negatives (Q).
#' @return Numeric vector of metric values.
#' @export
#' @examples
#' k <- 0:100
#' fom_upper_vals <- metric_along_bound(k, hits_upper, fom, 30, 70)
metric_along_bound <- function(k, hits_fn, metric_fn, n_positive, n_negative) {
  H <- hits_fn(k, n_positive, n_negative)
  mapply(function(h, ki) {
    FA <- ki - h
    M  <- n_positive - h
    CR <- n_negative - FA
    tryCatch(
      metric_fn(confusion_matrix(h, FA, M, CR)),
      error = function(e) NA_real_
    )
  }, H, k)
}
