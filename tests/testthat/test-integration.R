test_that("AUC of perfect classifier is 1", {
  P <- 10; Q <- 10; N <- P + Q
  # Perfect: all positives ranked first
  H <- c(rep(0, 1), seq(1, P), rep(P, Q))
  td <- toc_data(H, P, Q)
  expect_equal(auc(td), 1, tolerance = 1e-6)
  expect_equal(auc_roc(td), auc(td))  # backward compat alias
})

test_that("AUC of random classifier is ~0.5", {
  P <- 100; Q <- 100; N <- P + Q
  H <- seq(0, P, length.out = N + 1)
  td <- toc_data(H, P, Q)
  expect_equal(auc(td), 0.5, tolerance = 1e-6)
})

test_that("AUCS of perfect classifier is 1", {
  P <- 10; Q <- 10
  H <- c(rep(0, 1), seq(1, P), rep(P, Q))
  td <- toc_data(H, P, Q)
  expect_equal(aucs(td), 1, tolerance = 1e-6)
  expect_equal(dauc(td), aucs(td))  # backward compat alias
})

test_that("AUCSI of upper-bound classifier is 1", {
  P <- 20; Q <- 30; N <- P + Q
  H <- hits_upper(0:N, P)
  td <- toc_data(H, P, Q)
  expect_equal(aucsi(td), 1, tolerance = 1e-6)
  expect_equal(afom(td), aucsi(td))  # backward compat alias
})

test_that("MaxCSI returns value in [0,1]", {
  set.seed(7)
  td <- toc_from_scores(runif(200), runif(200) > 0.4)
  m  <- mcsi(td)
  expect_true(m >= 0 && m <= 1)
  expect_equal(mfom(td), mcsi(td))  # backward compat alias
})

test_that("AUPRC works and auc_prc is alias", {
  set.seed(9)
  td <- toc_from_scores(runif(200), runif(200) > 0.4)
  expect_equal(auprc(td), auc_prc(td))
  expect_true(auprc(td) >= 0 && auprc(td) <= 1)
})

test_that("area_ratio respects [0,1] range for valid inputs", {
  x  <- seq(0, 1, length.out = 101)
  ar <- area_ratio(x, y_curve = x * 0.8, y_upper = x, y_base = rep(0, 101))
  expect_true(ar >= 0 && ar <= 1)
})

test_that("integrated_metrics returns named vector with dissertation names", {
  set.seed(8)
  td  <- toc_from_scores(runif(100), runif(100) > 0.5)
  res <- integrated_metrics(td)
  expect_true(is.numeric(res))
  expect_true("MaxCSI"         %in% names(res))
  expect_true("AUCSI"          %in% names(res))
  expect_true("AUCSI_baseline" %in% names(res))
  expect_true("AUCSIS"         %in% names(res))
  expect_true("AUC"            %in% names(res))
  expect_true("AUCS"           %in% names(res))
  expect_true("AUPRC"          %in% names(res))
})
