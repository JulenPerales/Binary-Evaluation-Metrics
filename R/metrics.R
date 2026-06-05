# ---------------------------------------------------------------------------
# Internal helper ----
# ---------------------------------------------------------------------------

#' @keywords internal
.parse_cm <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  if (inherits(hits, "confusion_matrix")) return(hits)
  confusion_matrix(hits = hits, fa = fa, misses = misses, cr = cr)
}

# ---------------------------------------------------------------------------
# Agreement metrics ----
# ---------------------------------------------------------------------------

#' Overall Accuracy (OA)
#'
#' Proportion of correctly classified cases.
#' `OA = (TP + TN) / N`
#'
#' @param hits A `confusion_matrix` object, **or** the TP (Hits) count.
#' @param fa FP (False Alarms) count; ignored when `hits` is a `confusion_matrix`.
#' @param misses FN (Misses) count.
#' @param cr TN (Correct Rejections) count.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' oa(confusion_matrix(28, 72, 23, 2680))
#' oa(28, 72, 23, 2680)
oa <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  (cm$hits + cm$cr) / cm$n_total
}

#' Balanced Accuracy (BA)
#'
#' Arithmetic mean of sensitivity and specificity. Unaffected by class
#' imbalance.
#' `BA = 0.5 * (TPR + TNR)`
#'
#' @inheritParams oa
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' ba(confusion_matrix(28, 72, 23, 2680))
ba <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  0.5 * (cm$hits / cm$n_positive + cm$cr / cm$n_negative)
}

#' Matthews Correlation Coefficient (MCC)
#'
#' Balanced metric that accounts for all four cells of the confusion matrix.
#' Equivalent to the Pearson correlation between observed and predicted.
#'
#' @inheritParams oa
#' @return Numeric scalar in \[-1, 1\].
#' @export
#' @examples
#' mcc(confusion_matrix(28, 72, 23, 2680))
mcc <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm  <- .parse_cm(hits, fa, misses, cr)
  num <- cm$hits * cm$cr - cm$fa * cm$misses
  den <- sqrt(
    (cm$hits + cm$fa) * (cm$hits + cm$misses) *
    (cm$cr  + cm$fa) * (cm$cr  + cm$misses)
  )
  if (den == 0) return(NA_real_)
  num / den
}

#' Cohen's Kappa
#'
#' Chance-corrected agreement metric.
#' `Kappa = (po - pe) / (1 - pe)`.
#'
#' @inheritParams oa
#' @return Numeric scalar.
#' @export
#' @examples
#' kappa_score(confusion_matrix(28, 72, 23, 2680))
kappa_score <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  N  <- cm$n_total
  po <- (cm$hits + cm$cr) / N
  pe <- ((cm$hits + cm$fa) * (cm$hits + cm$misses) +
         (cm$misses + cm$cr) * (cm$fa + cm$cr)) / N^2
  if (1 - pe == 0) return(NA_real_)
  (po - pe) / (1 - pe)
}

#' F1 Score
#'
#' Harmonic mean of precision and recall.
#' `F1 = 2*TP / (2*TP + FP + FN) = 2*TP / (Predicted+ + Observed+)`.
#'
#' @param hits A `confusion_matrix` object or TP count.
#' @param fa FP count.
#' @param misses FN count.
#' @param cr TN count. **Not required for F1** — omit when passing raw counts.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' f1_score(confusion_matrix(28, 72, 23, 2680))
#' f1_score(28, 72, 23)
f1_score <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  if (inherits(hits, "confusion_matrix")) {
    cm    <- hits
    denom <- cm$n_positive + cm$hits + cm$fa
    return(ifelse(denom == 0, NA_real_, 2 * cm$hits / denom))
  }
  stopifnot(!is.null(fa), !is.null(misses))
  denom <- (hits + misses) + hits + fa   # P + TP + FP
  ifelse(denom == 0, NA_real_, 2 * hits / denom)
}

#' F-beta Score
#'
#' Generalised F score weighting recall `beta` times more than precision.
#' `F_beta = (1 + beta^2) * TP / ((1 + beta^2)*TP + beta^2*FN + FP)`.
#'
#' @param hits A `confusion_matrix` object or TP count.
#' @param fa FP count.
#' @param misses FN count.
#' @param cr TN count. **Not required for F-beta** — omit when passing raw counts.
#' @param beta Positive numeric weight. `beta = 1` gives F1; `beta = 2`
#'   weights recall twice as much.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' f_beta(confusion_matrix(28, 72, 23, 2680), beta = 2)
#' f_beta(28, 72, 23, beta = 2)
f_beta <- function(hits, fa = NULL, misses = NULL, cr = NULL, beta = 1) {
  if (inherits(hits, "confusion_matrix")) {
    cm    <- hits
    b2    <- beta^2
    num   <- (1 + b2) * cm$hits
    denom <- num + b2 * cm$misses + cm$fa
    return(ifelse(denom == 0, NA_real_, num / denom))
  }
  stopifnot(!is.null(fa), !is.null(misses))
  b2    <- beta^2
  num   <- (1 + b2) * hits
  denom <- num + b2 * misses + fa
  ifelse(denom == 0, NA_real_, num / denom)
}

# ---------------------------------------------------------------------------
# Skill metrics ----
# ---------------------------------------------------------------------------

#' Figure of Merit (FOM) / Critical Success Index (CSI)
#'
#' Fraction of events that were correctly forecast, out of all that were either
#' observed or forecast (or both). Also known as the Jaccard similarity
#' coefficient.
#' `FOM = TP / (TP + FP + FN)`
#'
#' @param hits A `confusion_matrix` object, **or** the TP (Hits) count.
#' @param fa FP (False Alarms) count; ignored when `hits` is a `confusion_matrix`.
#' @param misses FN (Misses) count.
#' @param cr TN (Correct Rejections) count. **Not required for FOM** — omit
#'   when passing raw counts.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' fom(confusion_matrix(28, 72, 23, 2680))
#' fom(28, 72, 23)   # cr not needed for FOM
fom <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  if (inherits(hits, "confusion_matrix")) {
    cm    <- hits
    denom <- cm$hits + cm$fa + cm$misses
    return(ifelse(denom == 0, NA_real_, cm$hits / denom))
  }
  stopifnot(!is.null(fa), !is.null(misses))
  denom <- hits + fa + misses
  ifelse(denom == 0, NA_real_, hits / denom)
}

#' Gilbert Skill Score (GSS) / Equitable Threat Score (ETS)
#'
#' Chance-corrected version of FOM. Accounts for the fraction of Hits expected
#' by random chance (`Hr`).
#' `GSS = (TP - Hr) / (TP + FP + FN - Hr)`, where
#' `Hr = (TP + FP)(TP + FN) / N`.
#'
#' @inheritParams oa
#' @return Numeric scalar in \[-1/3, 1\].
#' @export
#' @examples
#' gss(confusion_matrix(28, 72, 23, 2680))
gss <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  k  <- cm$hits + cm$fa       # predicted positives
  Hr <- k * cm$n_positive / cm$n_total
  den <- k + cm$misses - Hr
  if (den == 0) return(NA_real_)
  (cm$hits - Hr) / den
}

#' Peirce Skill Score (PSS) / Hanssen-Kuipers Discriminant
#'
#' Difference between True Positive Rate and False Positive Rate.
#' `PSS = TPR - FPR = TP/(TP+FN) - FP/(FP+TN)`.
#'
#' @inheritParams oa
#' @return Numeric scalar in \[-1, 1\].
#' @export
#' @examples
#' pss(confusion_matrix(28, 72, 23, 2680))
pss <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  cm$hits / cm$n_positive - cm$fa / cm$n_negative
}

#' Heidke Skill Score (HSS)
#'
#' `HSS = 2*(TP*TN - FP*FN) / ((TP+FP)(FP+TN) + (TP+FN)(FN+TN))`
#'
#' @inheritParams oa
#' @return Numeric scalar in \(-\infty, 1\].
#' @export
#' @examples
#' hss(confusion_matrix(28, 72, 23, 2680))
hss <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm  <- .parse_cm(hits, fa, misses, cr)
  num <- 2 * (cm$hits * cm$cr - cm$fa * cm$misses)
  den <- (cm$hits + cm$fa) * (cm$fa + cm$cr) +
         (cm$hits + cm$misses) * (cm$misses + cm$cr)
  if (den == 0) return(NA_real_)
  num / den
}

#' Equitable Threat Score (ETS)
#'
#' Alias for [gss()]. Both names refer to the same metric.
#'
#' @inheritParams oa
#' @return Numeric scalar.
#' @export
ets_score <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  gss(hits, fa, misses, cr)
}

# ---------------------------------------------------------------------------
# Auxiliary rates ----
# ---------------------------------------------------------------------------

#' True Positive Rate (Sensitivity / Recall / POD)
#'
#' `TPR = TP / (TP + FN)`
#'
#' @param hits A `confusion_matrix` object or TP count.
#' @param fa FP count (ignored for this metric unless using a `confusion_matrix`).
#' @param misses FN count.
#' @param cr TN count. Not required.
#' @return Numeric scalar in \[0, 1\].
#' @export
sensitivity <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  if (inherits(hits, "confusion_matrix")) {
    cm <- hits; return(cm$hits / cm$n_positive)
  }
  stopifnot(!is.null(misses))
  hits / (hits + misses)
}

#' True Negative Rate (Specificity)
#'
#' `TNR = TN / (TN + FP)`
#'
#' @inheritParams oa
#' @return Numeric scalar in \[0, 1\].
#' @export
specificity <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  cm$cr / cm$n_negative
}

#' False Positive Rate (POFD)
#'
#' `FPR = FP / (FP + TN)`
#'
#' @inheritParams oa
#' @return Numeric scalar in \[0, 1\].
#' @export
false_positive_rate <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  cm$fa / cm$n_negative
}

#' False Alarm Ratio (FAR)
#'
#' Fraction of positive predictions that were wrong.
#' `FAR = FP / (TP + FP)`.
#'
#' @param hits A `confusion_matrix` object or TP count.
#' @param fa FP count.
#' @param misses FN count. Not required.
#' @param cr TN count. Not required.
#' @return Numeric scalar in \[0, 1\].
#' @export
false_alarm_ratio <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  if (inherits(hits, "confusion_matrix")) {
    cm    <- hits
    denom <- cm$hits + cm$fa
    return(ifelse(denom == 0, NA_real_, cm$fa / denom))
  }
  stopifnot(!is.null(fa))
  denom <- hits + fa
  ifelse(denom == 0, NA_real_, fa / denom)
}

# ---------------------------------------------------------------------------
# Convenience: compute all metrics at once ----
# ---------------------------------------------------------------------------

#' Compute All Metrics from a Confusion Matrix
#'
#' Returns a named numeric vector with all agreement and skill metrics.
#'
#' @param cm A `confusion_matrix` object.
#' @return A named numeric vector.
#' @export
#' @examples
#' all_metrics(confusion_matrix(28, 72, 23, 2680))
all_metrics <- function(cm) {
  stopifnot(inherits(cm, "confusion_matrix"))
  c(
    OA          = oa(cm),
    BA          = ba(cm),
    MCC         = mcc(cm),
    Kappa       = kappa_score(cm),
    F1          = f1_score(cm),
    FOM         = fom(cm),
    GSS         = gss(cm),
    PSS         = pss(cm),
    HSS         = hss(cm),
    Sensitivity = sensitivity(cm),
    Specificity = specificity(cm),
    FPR         = false_positive_rate(cm),
    FAR         = false_alarm_ratio(cm)
  )
}
