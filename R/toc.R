#' Build a TOC (Threshold Operating Characteristic) Dataset
#'
#' The TOC dataset is the central data structure of **binmetrics**. It tabulates
#' the confusion matrix at every possible threshold `k` (from 0 to N), and
#' precomputes the upper, lower and random bounds for both the Hit count and the
#' FOM metric.
#'
#' @section Parametrisation:
#' A classifier is represented by the vector of cumulative **Hit counts** at
#' each threshold rank `k = 0, 1, …, N`:
#' - `k = 0`: no positives predicted → Hits = 0
#' - `k = N`: everything predicted positive → Hits = P
#'
#' @param hits_vector Numeric vector of length `N + 1`. Element `k+1` is the
#'   number of Hits when the `k` highest-ranked cases are labelled positive.
#'   Must be non-decreasing and bounded by `[hits_lower(k), hits_upper(k)]`.
#' @param n_positive Integer. Total number of positives P.
#' @param n_negative Integer. Total number of negatives Q.
#'
#' @return A `toc_data` tibble with one row per threshold value `k`, containing:
#'   \itemize{
#'     \item `k` — threshold rank (predicted positives count)
#'     \item `hits`, `fa`, `misses`, `cr` — confusion matrix cells
#'     \item `hits_upper`, `hits_lower`, `hits_random` — bound hit counts
#'     \item `fom_curve`, `fom_upper`, `fom_lower`, `fom_random` — FOM values
#'     \item `tpr`, `fpr` — True and False Positive Rates (for ROC)
#'   }
#'
#' @export
#' @seealso [toc_from_scores()], [toc_metric_curve()]
#' @examples
#' # Small synthetic example: P=5, Q=10
#' hits_vec <- c(0, 0, 1, 2, 2, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5)
#' td <- toc_data(hits_vec, n_positive = 5, n_negative = 10)
#' td
toc_data <- function(hits_vector, n_positive, n_negative) {
  stopifnot(
    is.numeric(hits_vector),
    length(hits_vector) == n_positive + n_negative + 1,
    all(diff(hits_vector) >= -1e-9),   # non-decreasing
    hits_vector[1] == 0,
    tail(hits_vector, 1) == n_positive
  )

  N <- n_positive + n_negative
  k <- 0:N

  H  <- hits_vector
  FA <- k - H
  M  <- n_positive - H
  CR <- n_negative - FA

  tibble::tibble(
    k           = k,
    hits        = H,
    fa          = FA,
    misses      = M,
    cr          = CR,
    # TOC bounds (hit counts)
    hits_upper  = hits_upper(k, n_positive),
    hits_lower  = hits_lower(k, n_positive, n_negative),
    hits_random = hits_random(k, n_positive / N),
    # FOM metric curves
    fom_curve   = .safe_fom(H,                          k, n_positive),
    fom_upper   = .safe_fom(hits_upper(k, n_positive),  k, n_positive),
    fom_lower   = .safe_fom(hits_lower(k, n_positive, n_negative), k, n_positive),
    fom_random  = .safe_fom(hits_random(k, n_positive / N),        k, n_positive),
    # ROC coordinates
    tpr         = H / n_positive,
    fpr         = FA / n_negative
  ) |>
    structure(
      class      = c("toc_data", "tbl_df", "tbl", "data.frame"),
      n_positive = n_positive,
      n_negative = n_negative,
      n_total    = N
    )
}

#' Build TOC Dataset from Score Vector
#'
#' Convenience constructor that converts a vector of continuous classifier
#' scores and a binary observed vector into a `toc_data` object. Cases are
#' ranked from highest to lowest score, and the cumulative Hit count is
#' computed at every rank.
#'
#' @param scores Numeric vector of classifier scores (higher = more likely
#'   positive). Need not be in \[0, 1\].
#' @param observed Logical or integer (0/1) vector indicating observed class
#'   (`TRUE` / `1` = positive event).
#' @return A `toc_data` tibble.
#' @export
#' @examples
#' set.seed(42)
#' scores   <- runif(100)
#' observed <- scores > 0.4 + rnorm(100, 0, 0.2)
#' td <- toc_from_scores(scores, observed)
toc_from_scores <- function(scores, observed) {
  stopifnot(
    is.numeric(scores),
    length(scores) == length(observed)
  )
  observed   <- as.logical(observed)
  n_positive <- sum(observed)
  n_negative <- sum(!observed)
  N          <- length(scores)

  ord        <- order(scores, decreasing = TRUE)
  hits_vec   <- c(0, cumsum(observed[ord]))

  toc_data(hits_vec, n_positive, n_negative)
}

#' @export
print.toc_data <- function(x, ...) {
  P  <- attr(x, "n_positive")
  Q  <- attr(x, "n_negative")
  N  <- attr(x, "n_total")
  pr <- P / N
  cat(sprintf(
    "TOC Dataset: N = %d | P = %d | Q = %d | prevalence = %.4f\n",
    N, P, Q, pr
  ))
  NextMethod()
}

# ---------------------------------------------------------------------------
# Derive a metric curve from a toc_data object ----
# ---------------------------------------------------------------------------

#' Compute Any Metric Curve from TOC Data
#'
#' Applies a metric function row-by-row to a `toc_data` object, returning a
#' tibble with the metric value alongside the corresponding upper, lower and
#' random bounds for the same metric.
#'
#' @param toc A `toc_data` object from [toc_data()] or [toc_from_scores()].
#' @param metric_fn A metric function from this package (e.g. [oa()], [gss()],
#'   [ba()]).  Must accept a `confusion_matrix` object.
#' @param metric_name Character. Label used for the y-axis in plots.
#' @return A tibble with columns `k`, `k_pct`, `curve`, `upper`, `lower`,
#'   `random`, and classification (`"Agreement"` or `"Skill"`).
#' @export
#' @examples
#' td  <- toc_from_scores(runif(200), runif(200) > 0.6)
#' oa_curve <- toc_metric_curve(td, oa, "OA")
#' gss_curve <- toc_metric_curve(td, gss, "GSS")
toc_metric_curve <- function(toc, metric_fn, metric_name = deparse(substitute(metric_fn))) {
  stopifnot(inherits(toc, "toc_data"))
  P <- attr(toc, "n_positive")
  Q <- attr(toc, "n_negative")
  N <- attr(toc, "n_total")

  apply_metric <- function(H, FA, M, CR) {
    tryCatch(
      metric_fn(confusion_matrix(H, FA, M, CR)),
      error = function(e) NA_real_
    )
  }

  out <- tibble::tibble(
    k       = toc$k,
    k_pct   = toc$k / N,
    curve   = mapply(apply_metric, toc$hits,        toc$fa,                    toc$misses,        toc$cr),
    upper   = mapply(apply_metric, toc$hits_upper,  toc$k - toc$hits_upper,   P - toc$hits_upper,  Q - (toc$k - toc$hits_upper)),
    lower   = mapply(apply_metric, toc$hits_lower,  toc$k - toc$hits_lower,   P - toc$hits_lower,  Q - (toc$k - toc$hits_lower)),
    random  = mapply(apply_metric, toc$hits_random, toc$k - toc$hits_random,  P - toc$hits_random, Q - (toc$k - toc$hits_random))
  )

  random_sum <- sum(out$random, na.rm = TRUE)
  type <- if (random_sum < 0.5) "Skill" else "Agreement"

  out$type        <- type
  out$metric_name <- metric_name
  out
}

# ---------------------------------------------------------------------------
# Internal helpers ----
# ---------------------------------------------------------------------------

.safe_fom <- function(H, k, P) {
  denom <- H + (k - H) + (P - H)
  ifelse(denom == 0, 0, H / denom)
}
