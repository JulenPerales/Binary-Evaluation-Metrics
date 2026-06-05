set.seed(2024)
N  <- 500
pr <- 0.2
observed     <- sample(c(TRUE, FALSE), N, replace = TRUE, prob = c(pr, 1 - pr))
noise_a      <- rnorm(N, sd = 0.5)
noise_b      <- rnorm(N, sd = 0.6)
scores_a     <- ifelse(observed, 0.65 + noise_a, 0.35 + noise_a)
scores_b     <- ifelse(observed, 0.70 + noise_b, 0.30 + noise_b)
scores_naive <- rep(0.5, N)
example_classifiers <- list(
  observed     = observed,
  scores_A     = scores_a,
  scores_B     = scores_b,
  scores_naive = scores_naive,
  N            = N,
  prevalence   = mean(observed)
)
save(example_classifiers, file = "data/example_classifiers.rda", compress = "xz")
cat("Saved data/example_classifiers.rda\n")
