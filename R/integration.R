# ---------------------------------------------------------------------------
# Area ratios and integrated metrics ----
# ---------------------------------------------------------------------------

#' Area Ratio Between Two Curves
#'
#' Computes the normalised area between a classifier curve and a reference
#' baseline, divided by the area between the upper bound and the same baseline.
#' Uses the trapezoidal rule for integration.
#'
#' `area_ratio = integral(curve - base) / integral(upper - base)`
#'
#' Negative areas are clamped to zero before integration so that the ratio
#' reflects performance above the baseline only.
#'
#' @param x Numeric vector. Common x-axis values.
#' @param y_curve Numeric vector. Classifier curve values.
#' @param y_upper Numeric vector. Upper bound values.
#' @param y_base Numeric vector. Baseline/lower reference values.
#' @return Numeric scalar in \[0, 1\] (can exceed 1 only if `y_curve > y_upper`
#'   at some points, which should not occur for valid classifiers).
#' @export
#' @examples
#' x   <- seq(0, 1, length.out = 101)
#' area_ratio(x, y_curve = x^0.5, y_upper = rep(1, 101), y_base = x)
area_ratio <- function(x, y_curve, y_upper, y_base) {
  # Drop NAs consistently across all vectors before integrating
  ok  <- complete.cases(x, y_curve, y_upper, y_base)
  if (!any(ok)) return(NA_real_)
  num <- pracma::trapz(x[ok], pmax(0, y_curve[ok] - y_base[ok]))
  den <- pracma::trapz(x[ok], pmax(0, y_upper[ok] - y_base[ok]))
  if (isTRUE(den == 0) || is.na(den)) return(NA_real_)
  num / den
}

# ---------------------------------------------------------------------------
# AFOM: Integrated FOM ----
# ---------------------------------------------------------------------------

#' Integrated FOM (AFOM)
#'
#' Area under the FOM curve, normalised by the area between the FOM upper
#' bound and the FOM lower bound (integrated over all thresholds k).
#'
#' `AFOM = integral(FOM_curve - FOM_lower) / integral(FOM_upper - FOM_lower)`
#'
#' @param toc A `toc_data` object from [toc_data()] or [toc_from_scores()].
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @seealso [aufom()], [dfom()]
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' afom(td)
afom <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  area_ratio(
    x       = toc$k,
    y_curve = toc$fom_curve,
    y_upper = toc$fom_upper,
    y_base  = toc$fom_lower
  )
}

#' Integrated FOM of a Uniform (Random) Classifier (AUFOM)
#'
#' AFOM value that a random classifier with AUC = 0.5 would achieve.  Serves
#' as the chance baseline for [dfom()].
#'
#' `AUFOM = integral(FOM_random - FOM_lower) / integral(FOM_upper - FOM_lower)`
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[0, 1\].
#' @export
aufom <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  area_ratio(
    x       = toc$k,
    y_curve = toc$fom_random,
    y_upper = toc$fom_upper,
    y_base  = toc$fom_lower
  )
}

#' FOM Skill: Normalised Improvement of AFOM over Random Baseline (DFOM)
#'
#' Standardises AFOM relative to the random baseline:
#' `DFOM = (AFOM - AUFOM) / (1 - AUFOM)`.
#' A value of 0 corresponds to a random classifier; 1 corresponds to a perfect
#' classifier.
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar.
#' @export
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' dfom(td)
dfom <- function(toc) {
  a  <- afom(toc)
  au <- aufom(toc)
  if (is.na(a) || is.na(au)) return(NA_real_)
  if (1 - au == 0) return(NA_real_)
  (a - au) / (1 - au)
}

# ---------------------------------------------------------------------------
# AUC (Area Under ROC curve) ----
# ---------------------------------------------------------------------------

#' Area Under the ROC Curve (AUC)
#'
#' Computed via trapezoidal integration of TPR over FPR.
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' set.seed(42)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' auc_roc(td)
auc_roc <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  ord <- order(toc$fpr)
  pracma::trapz(toc$fpr[ord], toc$tpr[ord])
}

#' AUC Skill Score (DAUC)
#'
#' Normalised deviation of AUC from the random classifier baseline (0.5):
#' `DAUC = 2 * (AUC - 0.5)`.
#' Ranges from -1 (worst) through 0 (random) to 1 (perfect).
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[-1, 1\].
#' @export
dauc <- function(toc) {
  2 * (auc_roc(toc) - 0.5)
}

# ---------------------------------------------------------------------------
# MFOM: Maximum FOM ----
# ---------------------------------------------------------------------------

#' Maximum FOM Across All Thresholds (MFOM)
#'
#' The best FOM value achievable by optimally choosing the decision threshold.
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[0, 1\].
#' @export
mfom <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  max(toc$fom_curve, na.rm = TRUE)
}

#' Threshold at Which MFOM is Achieved
#'
#' @param toc A `toc_data` object.
#' @return Integer scalar. The value of `k` that maximises FOM.
#' @export
mfom_threshold <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  toc$k[which.max(toc$fom_curve)]
}

# ---------------------------------------------------------------------------
# Area Under Precision-Recall Curve (AUCPRC) ----
# ---------------------------------------------------------------------------

#' Area Under the Precision-Recall Curve (AUCPRC)
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[0, 1\].
#' @export
auc_prc <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  P   <- attr(toc, "n_positive")
  recall    <- toc$hits / P
  precision <- ifelse(toc$k == 0, 1, toc$hits / toc$k)
  ord       <- order(recall)
  pracma::trapz(recall[ord], precision[ord])
}

#' Precision-Recall Skill Score
#'
#' `PRC_Skill = (AUCPRC - prevalence) / (1 - prevalence)`
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar.
#' @export
prc_skill <- function(toc) {
  P  <- attr(toc, "n_positive")
  N  <- attr(toc, "n_total")
  pr <- P / N
  (auc_prc(toc) - pr) / (1 - pr)
}

# ---------------------------------------------------------------------------
# Convenience: integrated summary ----
# ---------------------------------------------------------------------------

#' Integrated Metric Summary
#'
#' Returns a named numeric vector with all integrated metrics for a classifier.
#'
#' @param toc A `toc_data` object.
#' @return Named numeric vector.
#' @export
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' integrated_metrics(td)
integrated_metrics <- function(toc) {
  c(
    MFOM    = mfom(toc),
    AFOM    = afom(toc),
    AUFOM   = aufom(toc),
    DFOM    = dfom(toc),
    AUC     = auc_roc(toc),
    DAUC    = dauc(toc),
    AUCPRC  = auc_prc(toc),
    PRC_Skill = prc_skill(toc)
  )
}
