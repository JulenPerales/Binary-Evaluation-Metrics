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
#' Balanced **Skill** metric that accounts for all four cells of the confusion
#' matrix. Equivalent to the Pearson correlation between observed and predicted.
#' The dissertation explicitly classifies MCC as a Skill metric.
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
#' For binary classification, Cohen's Kappa equals the Heidke Skill Score (HSS).
#' The dissertation (CSI Framework) criticizes the use of Kappa in rare event
#' evaluation because it conflates agreement and skill concepts.
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
# Agreement metrics (continued): CSI and related ----
# ---------------------------------------------------------------------------

#' Critical Success Index (CSI) / Figure of Merit (FOM) / Threat Score
#'
#' An **Agreement** metric. Fraction of events that were correctly forecast,
#' out of all that were either observed or forecast (or both). Also known as
#' the Jaccard similarity coefficient or Threat Score.
#' `CSI = TP / (TP + FP + FN)`
#'
#' This is the central metric of the CSI Framework (dissertation). Note that
#' despite its name, CSI is classified as an **Agreement** metric (not Skill)
#' because it does not correct for chance.
#'
#' @param hits A `confusion_matrix` object, **or** the TP (Hits) count.
#' @param fa FP (False Alarms) count; ignored when `hits` is a `confusion_matrix`.
#' @param misses FN (Misses) count.
#' @param cr TN (Correct Rejections) count. **Not required for CSI** — omit
#'   when passing raw counts.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' csi(confusion_matrix(28, 72, 23, 2680))
#' csi(28, 72, 23)   # cr not needed for CSI
csi <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  if (inherits(hits, "confusion_matrix")) {
    cm    <- hits
    denom <- cm$hits + cm$fa + cm$misses
    return(ifelse(denom == 0, NA_real_, cm$hits / denom))
  }
  stopifnot(!is.null(fa), !is.null(misses))
  denom <- hits + fa + misses
  ifelse(denom == 0, NA_real_, hits / denom)
}

#' Figure of Merit (FOM) — Alias for CSI
#'
#' Alias for [csi()]. Kept for backward compatibility.
#' `FOM = CSI = TP / (TP + FP + FN)`
#'
#' @inheritParams csi
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' fom(confusion_matrix(28, 72, 23, 2680))
#' fom(28, 72, 23)
fom <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  csi(hits, fa, misses, cr)
}

#' User's Accuracy (UA) / Precision / Positive Predictive Value
#'
#' An **Agreement** metric. Fraction of predicted positives that are true
#' positives. Corresponds to User's Accuracy in the remote sensing literature.
#' `UA = TP / (TP + FP)`
#'
#' @param hits A `confusion_matrix` object, **or** the TP (Hits) count.
#' @param fa FP (False Alarms) count; ignored when `hits` is a `confusion_matrix`.
#' @param misses FN (Misses) count. Not required.
#' @param cr TN (Correct Rejections) count. Not required.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' ua(confusion_matrix(28, 72, 23, 2680))
#' ua(28, 72)
ua <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  if (inherits(hits, "confusion_matrix")) {
    cm    <- hits
    denom <- cm$hits + cm$fa
    return(ifelse(denom == 0, NA_real_, cm$hits / denom))
  }
  stopifnot(!is.null(fa))
  denom <- hits + fa
  ifelse(denom == 0, NA_real_, hits / denom)
}

#' False Alarm Rate (FAR) / False Positive Rate / POFD
#'
#' An **Agreement** metric. Fraction of observed negatives that were
#' incorrectly predicted as positive. Also called False Positive Rate (FPR)
#' or Probability of False Detection (POFD).
#' `FAR = FP / (FP + TN) = F / (C + F)`
#'
#' Note: This is distinct from the False Alarm *Ratio* (= FP/(TP+FP)),
#' which is available as [false_alarm_ratio()].
#'
#' @param hits A `confusion_matrix` object, **or** the TP (Hits) count.
#' @param fa FP (False Alarms) count.
#' @param misses FN (Misses) count. Not required.
#' @param cr TN (Correct Rejections) count.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' far(confusion_matrix(28, 72, 23, 2680))
far <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  cm$fa / cm$n_negative
}

#' False Omission Rate (FOR)
#'
#' An **Agreement** metric. Fraction of predicted negatives that are
#' actually positive (false omissions).
#' `FOR = FN / (FN + TN) = M / (C + M)`
#'
#' The function is named `for_metric` because `for` is a reserved keyword in R.
#'
#' @param hits A `confusion_matrix` object, **or** the TP (Hits) count.
#' @param fa FP (False Alarms) count.
#' @param misses FN (Misses) count.
#' @param cr TN (Correct Rejections) count.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' for_metric(confusion_matrix(28, 72, 23, 2680))
for_metric <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm    <- .parse_cm(hits, fa, misses, cr)
  denom <- cm$misses + cm$cr
  ifelse(denom == 0, NA_real_, cm$misses / denom)
}

# ---------------------------------------------------------------------------
# Skill metrics ----
# ---------------------------------------------------------------------------

#' Revised Gilbert Skill Score (GSS) / Equitable Threat Score (ETS)
#'
#' A **Skill** metric. Chance-corrected version of CSI. The dissertation
#' (CSI Framework) identifies that the standard GSS formula is conceptually
#' wrong because its denominator for the random baseline uses H+F+M (the
#' actual union), whereas CSI_random should use H_rand+F_rand+M_rand (the
#' random union). This implementation uses the revised formula:
#'
#' `GSS = (CSI - CSI_random) / (1 - CSI_random)`
#'
#' where `CSI_random = H_rand / (H_rand + F_rand + M_rand)` and
#' - `H_rand = k * P / N`
#' - `F_rand = k * Q / N`
#' - `M_rand = P * (N - k) / N`
#' (`k = H + F`, `P = n_positive`, `Q = n_negative`, `N = n_total`).
#'
#' @inheritParams oa
#' @return Numeric scalar in \[-1, 1\].
#' @export
#' @examples
#' gss(confusion_matrix(28, 72, 23, 2680))
gss <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  k  <- cm$hits + cm$fa          # H + F (predicted positives)
  P  <- cm$n_positive
  E  <- cm$n_total
  H_rand <- k * P / E
  F_rand <- k * (E - P) / E
  M_rand <- P * (E - k) / E
  denom_random <- H_rand + F_rand + M_rand
  csi_actual   <- cm$hits / (cm$hits + cm$fa + cm$misses)
  csi_random   <- ifelse(denom_random == 0, NA_real_, H_rand / denom_random)
  denom_actual <- cm$hits + cm$fa + cm$misses
  if (denom_actual == 0) return(NA_real_)
  denom_skill  <- 1 - csi_random
  ifelse(is.na(denom_skill) | denom_skill == 0, NA_real_,
         (csi_actual - csi_random) / denom_skill)
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

#' Producer's Accuracy (PA) / Sensitivity / Recall / POD
#'
#' An **Agreement** metric. Fraction of observed positives that were correctly
#' predicted. Corresponds to Producer's Accuracy in the remote sensing
#' literature and to Sensitivity / Recall / Probability of Detection (POD) in
#' other fields.
#' `PA = TP / (TP + FN)`
#'
#' @param hits A `confusion_matrix` object or TP count.
#' @param fa FP count (ignored for this metric unless using a `confusion_matrix`).
#' @param misses FN count.
#' @param cr TN count. Not required.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' pa(confusion_matrix(28, 72, 23, 2680))
#' pa(28, fa = NULL, misses = 23)
pa <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  if (inherits(hits, "confusion_matrix")) {
    cm <- hits; return(cm$hits / cm$n_positive)
  }
  stopifnot(!is.null(misses))
  hits / (hits + misses)
}

#' Sensitivity / Recall / POD — Alias for PA
#'
#' Alias for [pa()]. Kept for backward compatibility.
#' `Sensitivity = PA = TP / (TP + FN)`
#'
#' @inheritParams pa
#' @return Numeric scalar in \[0, 1\].
#' @export
sensitivity <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  pa(hits, fa, misses, cr)
}

#' SP / Specificity / True Negative Rate
#'
#' An **Agreement** metric. Fraction of observed negatives that were correctly
#' predicted as negative.
#' `SP = TN / (TN + FP)`
#'
#' @inheritParams oa
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' sp(confusion_matrix(28, 72, 23, 2680))
sp <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  cm <- .parse_cm(hits, fa, misses, cr)
  cm$cr / cm$n_negative
}

#' Specificity — Alias for SP
#'
#' Alias for [sp()]. Kept for backward compatibility.
#' `Specificity = SP = TN / (TN + FP)`
#'
#' @inheritParams oa
#' @return Numeric scalar in \[0, 1\].
#' @export
specificity <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  sp(hits, fa, misses, cr)
}

#' False Positive Rate (POFD) — Alias for FAR
#'
#' Alias for [far()]. Kept for backward compatibility.
#' `FPR = FP / (FP + TN)`
#'
#' @inheritParams oa
#' @return Numeric scalar in \[0, 1\].
#' @export
false_positive_rate <- function(hits, fa = NULL, misses = NULL, cr = NULL) {
  far(hits, fa, misses, cr)
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
#' Returns a tibble with all agreement and skill metrics grouped according to
#' the CSI Framework (dissertation Chapter 1). Skill metrics carry two extra
#' columns: `baseline_agreement` (the value of the corresponding agreement
#' metric that a random classifier would achieve at the same prediction rate)
#' and `corresponding_agreement_metric` (the name of that agreement metric).
#'
#' **Agreement metrics** (do not correct for chance): OA, PA, SP, FAR, UA,
#' FOR, CSI, F1.
#' **Skill metrics** (chance-corrected): GSS -> CSI, HSS -> OA, PSS -> PA,
#' MCC -> OA.
#'
#' @param cm A `confusion_matrix` object.
#' @return A tibble with columns `metric`, `type`, `value`,
#'   `baseline_agreement`, `corresponding_agreement_metric`.
#' @export
#' @examples
#' all_metrics(confusion_matrix(28, 72, 23, 2680))
all_metrics <- function(cm) {
  stopifnot(inherits(cm, "confusion_matrix"))
  k <- cm$hits + cm$fa
  P <- cm$n_positive; N <- cm$n_total; Q <- cm$n_negative

  # Random baseline CSI, given the actual prediction count k (for GSS)
  denom_csi_rand <- N * (k + P) - k * P
  csi_rand <- if (denom_csi_rand == 0) NA_real_ else k * P / denom_csi_rand

  # Random baseline OA, given k (for HSS and MCC)
  oa_rand <- (k * P + (N - k) * Q) / N^2

  # Random baseline PA = prediction rate k/N (for PSS)
  pa_rand <- if (N == 0) NA_real_ else k / N

  tibble::tibble(
    metric = c("OA", "PA", "SP", "FAR", "UA", "FOR", "CSI", "F1",
               "GSS", "HSS", "PSS", "MCC"),
    type   = c(rep("Agreement", 8), rep("Skill", 4)),
    value  = c(
      oa(cm), pa(cm), sp(cm), far(cm), ua(cm), for_metric(cm), csi(cm), f1_score(cm),
      gss(cm), hss(cm), pss(cm), mcc(cm)
    ),
    baseline_agreement = c(
      rep(NA_real_, 8),
      csi_rand,  # GSS: CSI for random classifier
      oa_rand,   # HSS: OA for random classifier
      pa_rand,   # PSS: PA for random classifier
      oa_rand    # MCC: OA for random classifier
    ),
    corresponding_agreement_metric = c(
      rep(NA_character_, 8),
      "CSI", "OA", "PA", "OA"
    )
  )
}
