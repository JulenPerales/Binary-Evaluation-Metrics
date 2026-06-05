test_that("toc_data constructs correctly", {
  # Perfect classifier: P=5, Q=10
  P  <- 5; Q <- 10; N <- P + Q
  H  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5)
  td <- toc_data(H, P, Q)
  expect_s3_class(td, "toc_data")
  expect_equal(nrow(td), N + 1)
  expect_equal(attr(td, "n_positive"), P)
  expect_equal(attr(td, "n_negative"), Q)
})

test_that("toc_data rejects invalid hits_vector", {
  expect_error(toc_data(c(1, 0, 1), 1, 1))   # not non-decreasing (diff = -1)
  expect_error(toc_data(c(0, 1),    1, 1))    # wrong length (need N+1 = 3)
  expect_error(toc_data(c(0, 1, 3), 1, 1))   # last element 3 != n_positive 1
})

test_that("toc_from_scores produces valid toc_data", {
  set.seed(1)
  scores <- runif(100)
  obs    <- scores > 0.5
  td     <- toc_from_scores(scores, obs)
  expect_s3_class(td, "toc_data")
  P      <- sum(obs)
  expect_equal(attr(td, "n_positive"), P)
  expect_equal(nrow(td), 101)
  expect_equal(td$hits[1], 0)
  expect_equal(tail(td$hits, 1), P)
})

test_that("TOC bounds are correct", {
  set.seed(2)
  td <- toc_from_scores(runif(50), runif(50) > 0.4)
  P  <- attr(td, "n_positive")
  Q  <- attr(td, "n_negative")
  expect_true(all(td$hits_upper == pmin(td$k, P)))
  expect_true(all(td$hits_lower == pmax(0, td$k - Q)))
})

test_that("CSI upper >= CSI curve >= CSI lower", {
  set.seed(3)
  td <- toc_from_scores(runif(100), runif(100) > 0.45)
  # Allow small floating point tolerance
  expect_true(all(td$csi_upper >= td$csi - 1e-9))
  expect_true(all(td$csi       >= td$csi_lower - 1e-9))
  # Backward-compatible columns still exist
  expect_true("fom_curve"  %in% names(td))
  expect_true("fom_upper"  %in% names(td))
  expect_true("fom_lower"  %in% names(td))
  expect_true("fom_random" %in% names(td))
})

test_that("toc_metric_curve runs for all metrics", {
  set.seed(4)
  td <- toc_from_scores(runif(80), runif(80) > 0.4)
  for (fn in list(oa, ba, mcc, fom, gss, pss, hss)) {
    mc <- toc_metric_curve(td, fn)
    expect_true(is.data.frame(mc))
    expect_true("curve" %in% names(mc))
  }
})
