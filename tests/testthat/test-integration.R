test_that("AUC of perfect classifier is 1", {
  P <- 10; Q <- 10; N <- P + Q
  # Perfect: all positives ranked first
  H <- c(rep(0, 1), seq(1, P), rep(P, Q))
  td <- toc_data(H, P, Q)
  expect_equal(auc_roc(td), 1, tolerance = 1e-6)
})

test_that("AUC of random classifier is ~0.5", {
  P <- 100; Q <- 100; N <- P + Q
  H <- seq(0, P, length.out = N + 1)
  td <- toc_data(H, P, Q)
  expect_equal(auc_roc(td), 0.5, tolerance = 1e-6)
})

test_that("DAUC of perfect classifier is 1", {
  P <- 10; Q <- 10
  H <- c(rep(0, 1), seq(1, P), rep(P, Q))
  td <- toc_data(H, P, Q)
  expect_equal(dauc(td), 1, tolerance = 1e-6)
})

test_that("AFOM of upper-bound classifier is 1", {
  P <- 20; Q <- 30; N <- P + Q
  H <- hits_upper(0:N, P)
  td <- toc_data(H, P, Q)
  expect_equal(afom(td), 1, tolerance = 1e-6)
})

test_that("MFOM returns value in [0,1]", {
  set.seed(7)
  td <- toc_from_scores(runif(200), runif(200) > 0.4)
  m  <- mfom(td)
  expect_true(m >= 0 && m <= 1)
})

test_that("area_ratio respects [0,1] range for valid inputs", {
  x  <- seq(0, 1, length.out = 101)
  ar <- area_ratio(x, y_curve = x * 0.8, y_upper = x, y_base = rep(0, 101))
  expect_true(ar >= 0 && ar <= 1)
})

test_that("integrated_metrics returns named vector", {
  set.seed(8)
  td  <- toc_from_scores(runif(100), runif(100) > 0.5)
  res <- integrated_metrics(td)
  expect_true(is.numeric(res))
  expect_true("MFOM"  %in% names(res))
  expect_true("AUC"   %in% names(res))
  expect_true("AFOM"  %in% names(res))
})
