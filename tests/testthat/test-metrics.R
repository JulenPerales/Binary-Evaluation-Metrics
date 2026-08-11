test_that("confusion_matrix constructs correctly", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  expect_equal(cm$hits, 28)
  expect_equal(cm$n_positive, 51)
  expect_equal(cm$n_negative, 2752)
  expect_equal(cm$n_total, 2803)
  expect_equal(cm$prevalence, 51 / 2803, tolerance = 1e-9)
})

test_that("confusion_matrix rejects invalid inputs", {
  expect_error(confusion_matrix(-1, 0, 0, 10))
  expect_error(confusion_matrix(0, 0, 0, 0))
})

test_that("OA is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  expect_equal(oa(cm), (28 + 2680) / 2803, tolerance = 1e-9)
  expect_equal(oa(28, 72, 23, 2680), oa(cm))
})

test_that("BA is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  tpr <- 28 / 51
  tnr <- 2680 / 2752
  expect_equal(ba(cm), 0.5 * (tpr + tnr), tolerance = 1e-9)
})

test_that("MCC is correct for Finley example", {
  # Finley example from the R scripts
  cm <- confusion_matrix(28, 72, 23, 2680)
  H <- 28; FA <- 72; M <- 23; CR <- 2680
  num <- H * CR - FA * M
  den <- sqrt((H + FA) * (H + M) * (CR + FA) * (CR + M))
  expect_equal(mcc(cm), num / den, tolerance = 1e-9)
})

test_that("FOM equals CSI formula", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  expect_equal(fom(cm), 28 / (28 + 72 + 23), tolerance = 1e-9)
})

test_that("GSS (Revised) matches dissertation formula", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  H <- 28; FA <- 72; M <- 23; CR <- 2680
  k <- H + FA; P <- H + M; E <- H + FA + M + CR
  H_rand <- k * P / E
  F_rand <- k * (E - P) / E
  M_rand <- P * (E - k) / E
  csi_actual <- H / (H + FA + M)
  csi_random <- H_rand / (H_rand + F_rand + M_rand)
  expected   <- (csi_actual - csi_random) / (1 - csi_random)
  expect_equal(gss(cm), expected, tolerance = 1e-9)
  expect_equal(ets_score(cm), gss(cm))
})

test_that("PSS equals TPR - FPR", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  expect_equal(pss(cm), 28 / 51 - 72 / 2752, tolerance = 1e-9)
})

test_that("HSS formula is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  H <- 28; FA <- 72; M <- 23; CR <- 2680
  num <- 2 * (H * CR - FA * M)
  den <- (H + FA) * (FA + CR) + (H + M) * (M + CR)
  expect_equal(hss(cm), num / den, tolerance = 1e-9)
})

test_that("kappa_score is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  N  <- 2803
  po <- (28 + 2680) / N
  pe <- ((28 + 72) * (28 + 23) + (23 + 2680) * (72 + 2680)) / N^2
  expect_equal(kappa_score(cm), (po - pe) / (1 - pe), tolerance = 1e-6)
})

test_that("all_metrics returns a tibble with CSI Framework structure", {
  cm  <- confusion_matrix(28, 72, 23, 2680)
  res <- all_metrics(cm)
  expect_s3_class(res, "data.frame")
  # Core columns present
  expect_true(all(c("metric", "type", "value",
                     "baseline_agreement", "corresponding_agreement_metric") %in% names(res)))
  # All expected metrics present
  expect_true(all(c("OA","PA","SP","FAR","UA","FOR","CSI","F1",
                     "GSS","HSS","PSS","MCC") %in% res$metric))
  # Agreement vs Skill classification
  agreement <- res[res$type == "Agreement", ]
  skill     <- res[res$type == "Skill", ]
  expect_equal(nrow(agreement), 8)
  expect_equal(nrow(skill), 4)
  # Skill metrics have baseline_agreement and corresponding metric
  expect_true(all(!is.na(skill$baseline_agreement)))
  expect_true(all(!is.na(skill$corresponding_agreement_metric)))
  expect_equal(skill$corresponding_agreement_metric,
               c("CSI", "OA", "PA", "OA"))
  # Agreement metrics have NA baselines
  expect_true(all(is.na(agreement$baseline_agreement)))
  # baseline_agreement is last numeric column (before corresponding_agreement_metric)
  col_idx_base  <- which(names(res) == "baseline_agreement")
  col_idx_corr  <- which(names(res) == "corresponding_agreement_metric")
  expect_equal(col_idx_corr, col_idx_base + 1L)
})

test_that("perfect classifier has FOM = OA = 1", {
  cm <- confusion_matrix(100, 0, 0, 200)
  expect_equal(fom(cm), 1)
  expect_equal(oa(cm),  1)
  expect_equal(ba(cm),  1)
})

test_that("metrics accept raw arguments", {
  # FOM, F1, sensitivity do not require cr
  expect_equal(fom(28, 72, 23),       28 / (28 + 72 + 23))
  expect_equal(f1_score(28, 72, 23),  2 * 28 / ((28 + 23) + 28 + 72))
  expect_equal(sensitivity(28, fa = NULL, misses = 23), 28 / (28 + 23))
  # Metrics requiring all 4 cells
  expect_equal(oa(28, 72, 23, 2680),  (28 + 2680) / 2803, tolerance = 1e-9)
  expect_equal(gss(28, 72, 23, 2680), gss(confusion_matrix(28, 72, 23, 2680)))
})

test_that("UA (User's Accuracy / Precision) is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  expect_equal(ua(cm), 28 / (28 + 72), tolerance = 1e-9)
  expect_equal(ua(28, 72), 28 / (28 + 72), tolerance = 1e-9)
})

test_that("PA (Producer's Accuracy / Sensitivity) is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  expect_equal(pa(cm), 28 / (28 + 23), tolerance = 1e-9)
  expect_equal(sensitivity(cm), pa(cm))
})

test_that("SP (Specificity) is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  expect_equal(sp(cm), 2680 / (2680 + 72), tolerance = 1e-9)
  expect_equal(specificity(cm), sp(cm))
})

test_that("FAR (False Alarm Rate = FPR) is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  # FAR = F / (C + F) = FP / (TN + FP)
  expect_equal(far(cm), 72 / (2680 + 72), tolerance = 1e-9)
  expect_equal(false_positive_rate(cm), far(cm))
})

test_that("FOR (False Omission Rate) is correct", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  # FOR = M / (C + M) = FN / (TN + FN)
  expect_equal(for_metric(cm), 23 / (2680 + 23), tolerance = 1e-9)
})

test_that("CSI is primary and fom() is alias", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  expect_equal(csi(cm), 28 / (28 + 72 + 23), tolerance = 1e-9)
  expect_equal(fom(cm), csi(cm))
})

test_that("Revised GSS is skill: 0 for random, 1 for perfect", {
  # Perfect classifier: H=P, F=0, M=0
  cm_perfect <- confusion_matrix(50, 0, 0, 50)
  expect_equal(gss(cm_perfect), 1, tolerance = 1e-9)
  # For the random case, CSI_actual == CSI_random, so GSS = 0
  # Random: H_rand = k*P/N for k=P, H=H_rand=P^2/N
  # Use k=P=50, N=100, so H_rand=25
  cm_random <- confusion_matrix(25, 25, 25, 25)
  expect_equal(gss(cm_random), 0, tolerance = 1e-6)
})
