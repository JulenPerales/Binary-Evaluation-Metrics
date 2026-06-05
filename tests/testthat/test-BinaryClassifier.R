test_that("BinaryClassifier$new from scores", {
  set.seed(1)
  scores <- runif(100)
  obs    <- scores > 0.4
  clf    <- BinaryClassifier$new(scores, obs)
  expect_s3_class(clf, "R6")
  expect_equal(clf$n_total, 100)
  expect_equal(clf$n_positive, sum(obs))
  expect_equal(clf$prevalence, sum(obs) / 100)
})

test_that("BinaryClassifier$new from toc_data object", {
  P <- 5; Q <- 10
  H <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5)
  td  <- toc_data(H, P, Q)
  clf <- BinaryClassifier$new(td, name = "From TOC")
  expect_equal(clf$n_positive, P)
  expect_equal(clf$n_total, P + Q)
})

test_that("BinaryClassifier_from_hits constructor", {
  P <- 5; Q <- 10
  # Perfect classifier: all positives ranked first
  H <- c(0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
  clf <- BinaryClassifier_from_hits(H, P, Q, name = "Perfect")
  expect_equal(clf$n_positive, P)
  expect_equal(clf$n_negative, Q)
  expect_equal(clf$mfom, 1, tolerance = 1e-9)   # perfect classifier MFOM = 1
  expect_equal(clf$auc,  1, tolerance = 1e-6)   # perfect AUC = 1
})

test_that("at_threshold returns confusion_matrix", {
  set.seed(2)
  scores <- runif(50)
  obs    <- scores > 0.5
  clf    <- BinaryClassifier$new(scores, obs)
  cm     <- clf$at_threshold(10)
  expect_s3_class(cm, "confusion_matrix")
  expect_equal(cm$hits + cm$fa, 10)  # k=10 means 10 predicted positives
})

test_that("integrated metrics are in valid range", {
  set.seed(3)
  scores <- runif(200)
  obs    <- scores > 0.4
  clf    <- BinaryClassifier$new(scores, obs)
  expect_true(clf$auc   >= 0 && clf$auc   <= 1)
  expect_true(clf$afom  >= 0 && clf$afom  <= 1)
  expect_true(clf$mfom  >= 0 && clf$mfom  <= 1)
  expect_true(clf$dauc  >= -1 && clf$dauc  <= 1)
})

test_that("compute_metric_curve returns expected structure", {
  set.seed(4)
  clf <- BinaryClassifier$new(runif(80), runif(80) > 0.45)
  mc  <- clf$compute_metric_curve(oa, "OA")
  expect_true(is.data.frame(mc))
  expect_true(all(c("k", "k_pct", "curve", "upper", "lower", "random") %in% names(mc)))
  expect_equal(nrow(mc), clf$n_total + 1)
})
