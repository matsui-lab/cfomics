test_that("cf_fit_bcf basic execution", {
  skip_if_not_installed("bcf")

  set.seed(123)
  n <- 200
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  tau_true <- 2.0 + 0.5 * X[,1]
  Y <- tau_true * T + 0.5 * X[,1] + 0.3 * X[,2] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  # Build formula
  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  # Use fewer iterations for faster testing
  fit <- cf_fit(fml, data = data, method = "bcf",
                n_burn = 100, n_iter = 200)

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "bcf")

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  # True ATE is approximately mean(2.0 + 0.5 * X[,1]) ~ 2.0
  expect_true(ate > 1.0 && ate < 3.5)

  ite <- predict(fit, type = "ite")
  expect_equal(length(ite), n)
})

test_that("cf_fit_bcf returns posterior samples", {
  skip_if_not_installed("bcf")

  set.seed(456)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, 0.5)
  Y <- 2.0 * T + X[,1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "bcf",
                n_burn = 50, n_iter = 100)

  # Check that posterior samples are available
  samples <- predict(fit, type = "samples")
  expect_true("tau_posterior" %in% names(samples))
  expect_true("ate_posterior" %in% names(samples))
})

test_that("cf_fit_bcf provides credible intervals", {
  skip_if_not_installed("bcf")

  set.seed(789)
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, 0.5)
  Y <- 2.0 * T + X[,1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "bcf",
                n_burn = 50, n_iter = 100)

  summary_stats <- predict(fit, type = "summary")
  expect_true("ate_ci_lower" %in% names(summary_stats))
  expect_true("ate_ci_upper" %in% names(summary_stats))
  expect_true(summary_stats$ate_ci_lower < summary_stats$ate)
  expect_true(summary_stats$ate_ci_upper > summary_stats$ate)
})
