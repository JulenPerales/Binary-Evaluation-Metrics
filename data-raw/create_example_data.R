# Script to generate the package example datasets.
# Run this once from the package root to regenerate data/example_classifiers.rda

set.seed(2024)
N   <- 500
pr  <- 0.2

# True labels
observed <- sample(c(TRUE, FALSE), N, replace = TRUE, prob = c(pr, 1 - pr))

# Classifier A: moderately good (AUC ~ 0.82)
noise_a     <- rnorm(N, sd = 0.5)
scores_a    <- ifelse(observed, 0.65 + noise_a, 0.35 + noise_a)

# Classifier B: different profile — same AUC, shifted tipping point
noise_b     <- rnorm(N, sd = 0.6)
scores_b    <- ifelse(observed, 0.70 + noise_b, 0.30 + noise_b)

# Naive model: always predicts a constant score
scores_naive <- rep(0.5, N)

example_classifiers <- list(
  observed = observed,
  scores_A = scores_a,
  scores_B = scores_b,
  scores_naive = scores_naive,
  N        = N,
  prevalence = mean(observed)
)

usethis::use_data(example_classifiers, overwrite = TRUE)
