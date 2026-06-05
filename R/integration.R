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
# AUCSI: Area Under the CSI Curve (= AFOM) ----
# ---------------------------------------------------------------------------

#' Area Under the CSI Curve (AUCSI)
#'
#' Area under the CSI curve normalised by the feasible area (between the CSI
#' upper bound and the CSI lower bound), integrated over all thresholds k.
#' This is the primary integrated CSI metric from the dissertation.
#'
#' `AUCSI = integral(CSI_curve - CSI_lower) / integral(CSI_upper - CSI_lower)`
#'
#' @param toc A `toc_data` object from [toc_data()] or [toc_from_scores()].
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @seealso [aucsi_baseline()], [aucsis()]
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' aucsi(td)
aucsi <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  area_ratio(
    x       = toc$k,
    y_curve = toc$csi,
    y_upper = toc$csi_upper,
    y_base  = toc$csi_lower
  )
}

#' Integrated FOM (AFOM) — Alias for AUCSI
#'
#' Alias for [aucsi()]. Kept for backward compatibility.
#'
#' @inheritParams aucsi
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @seealso [aucsi()]
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' afom(td)
afom <- function(toc) {
  aucsi(toc)
}

#' AUCSI Baseline (Random Classifier) — AUCSI_baseline
#'
#' AUCSI value that a random classifier with AUC = 0.5 would achieve. Serves
#' as the chance baseline for [aucsis()].
#'
#' `AUCSI_baseline = integral(CSI_baseline - CSI_lower) / integral(CSI_upper - CSI_lower)`
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' aucsi_baseline(td)
aucsi_baseline <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  area_ratio(
    x       = toc$k,
    y_curve = toc$csi_baseline,
    y_upper = toc$csi_upper,
    y_base  = toc$csi_lower
  )
}

#' Integrated FOM of a Uniform (Random) Classifier (AUFOM) — Alias for AUCSI_baseline
#'
#' Alias for [aucsi_baseline()]. Kept for backward compatibility.
#'
#' @inheritParams aucsi_baseline
#' @return Numeric scalar in \[0, 1\].
#' @export
aufom <- function(toc) {
  aucsi_baseline(toc)
}

#' AUCSI Skill Score (AUCSIS)
#'
#' Standardises AUCSI relative to the random classifier baseline:
#' `AUCSIS = (AUCSI - AUCSI_baseline) / (1 - AUCSI_baseline)`.
#' A value of 0 corresponds to a random classifier; 1 corresponds to a perfect
#' classifier.
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar.
#' @export
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' aucsis(td)
aucsis <- function(toc) {
  a  <- aucsi(toc)
  au <- aucsi_baseline(toc)
  if (is.na(a) || is.na(au)) return(NA_real_)
  if (1 - au == 0) return(NA_real_)
  (a - au) / (1 - au)
}

#' FOM Skill: Normalised Improvement of AFOM over Random Baseline (DFOM) — Alias for AUCSIS
#'
#' Alias for [aucsis()]. Kept for backward compatibility.
#'
#' @inheritParams aucsis
#' @return Numeric scalar.
#' @export
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' dfom(td)
dfom <- function(toc) {
  aucsis(toc)
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
#' auc(td)
auc <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  ord <- order(toc$fpr)
  pracma::trapz(toc$fpr[ord], toc$tpr[ord])
}

#' Area Under the ROC Curve (AUC-ROC) — Alias for AUC
#'
#' Alias for [auc()]. Kept for backward compatibility.
#'
#' @inheritParams auc
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' set.seed(42)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' auc_roc(td)
auc_roc <- function(toc) {
  auc(toc)
}

#' AUC Skill Score (AUCS)
#'
#' Normalised deviation of AUC from the random classifier baseline (0.5):
#' `AUCS = 2 * (AUC - 0.5) = 2*AUC - 1`.
#' Ranges from -1 (worst) through 0 (random) to 1 (perfect).
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[-1, 1\].
#' @export
#' @examples
#' set.seed(42)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' aucs(td)
aucs <- function(toc) {
  2 * (auc(toc) - 0.5)
}

#' AUC Skill Score (DAUC) — Alias for AUCS
#'
#' Alias for [aucs()]. Kept for backward compatibility.
#'
#' @inheritParams aucs
#' @return Numeric scalar in \[-1, 1\].
#' @export
dauc <- function(toc) {
  aucs(toc)
}

# ---------------------------------------------------------------------------
# MaxCSI: Maximum CSI ----
# ---------------------------------------------------------------------------

#' Maximum CSI Across All Thresholds (MaxCSI)
#'
#' The best CSI value achievable by optimally choosing the decision threshold.
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' mcsi(td)
mcsi <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  max(toc$csi, na.rm = TRUE)
}

#' Maximum FOM Across All Thresholds (MFOM) — Alias for MaxCSI
#'
#' Alias for [mcsi()]. Kept for backward compatibility.
#'
#' @inheritParams mcsi
#' @return Numeric scalar in \[0, 1\].
#' @export
mfom <- function(toc) {
  mcsi(toc)
}

#' Threshold at Which MaxCSI is Achieved
#'
#' @param toc A `toc_data` object.
#' @return Integer scalar. The value of `k` that maximises CSI.
#' @export
#' @examples
#' set.seed(1)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' mcsi_threshold(td)
mcsi_threshold <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  toc$k[which.max(toc$csi)]
}

#' Threshold at Which MFOM is Achieved — Alias for mcsi_threshold
#'
#' Alias for [mcsi_threshold()]. Kept for backward compatibility.
#'
#' @inheritParams mcsi_threshold
#' @return Integer scalar.
#' @export
mfom_threshold <- function(toc) {
  mcsi_threshold(toc)
}

# ---------------------------------------------------------------------------
# AUPRC: Area Under Precision-Recall Curve ----
# ---------------------------------------------------------------------------

#' Area Under the Precision-Recall Curve (AUPRC)
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar in \[0, 1\].
#' @export
#' @examples
#' set.seed(42)
#' td <- toc_from_scores(runif(500), runif(500) > 0.4)
#' auprc(td)
auprc <- function(toc) {
  stopifnot(inherits(toc, "toc_data"))
  P   <- attr(toc, "n_positive")
  recall    <- toc$hits / P
  precision <- ifelse(toc$k == 0, 1, toc$hits / toc$k)
  ord       <- order(recall)
  pracma::trapz(recall[ord], precision[ord])
}

#' Area Under the Precision-Recall Curve (AUC-PRC) — Alias for AUPRC
#'
#' Alias for [auprc()]. Kept for backward compatibility.
#'
#' @inheritParams auprc
#' @return Numeric scalar in \[0, 1\].
#' @export
auc_prc <- function(toc) {
  auprc(toc)
}

#' Precision-Recall Skill Score
#'
#' `PRC_Skill = (AUPRC - prevalence) / (1 - prevalence)`
#'
#' @param toc A `toc_data` object.
#' @return Numeric scalar.
#' @export
prc_skill <- function(toc) {
  P  <- attr(toc, "n_positive")
  N  <- attr(toc, "n_total")
  pr <- P / N
  (auprc(toc) - pr) / (1 - pr)
}

# ---------------------------------------------------------------------------
# Convenience: integrated summary ----
# ---------------------------------------------------------------------------

#' Integrated Metric Summary
#'
#' Returns a named numeric vector with all integrated metrics for a classifier,
#' using the dissertation (CSI Framework) naming convention.
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
    MaxCSI    = mcsi(toc),
    AUCSI     = aucsi(toc),
    AUCSI_baseline = aucsi_baseline(toc),
    AUCSIS    = aucsis(toc),
    AUC       = auc(toc),
    AUCS      = aucs(toc),
    AUPRC     = auprc(toc),
    PRC_Skill = prc_skill(toc)
  )
}
