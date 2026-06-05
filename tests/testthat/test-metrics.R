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

test_that("GSS equals ETS", {
  cm <- confusion_matrix(28, 72, 23, 2680)
  Hr  <- (28 + 72) * (28 + 23) / 2803
  expected <- (28 - Hr) / (28 + 72 + 23 - Hr)
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

test_that("all_metrics returns named vector", {
  cm  <- confusion_matrix(28, 72, 23, 2680)
  res <- all_metrics(cm)
  expect_true(is.numeric(res))
  expect_named(res)
  expect_true("FOM" %in% names(res))
  expect_true("OA"  %in% names(res))
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
