test_that("cf_fit_hdps basic execution", {
  skip_if_not_installed("glmnet")

  set.seed(123)
  n <- 200
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  beta_t <- c(rep(0.3, 5), rep(0, p - 5))
  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)
  Y <- 2.0 * T + 0.5 * X[,1] + 0.3 * X[,2] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  # Build formula with explicit covariate names
  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "hdps")

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "hdps")

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_true(ate > 1.0 && ate < 3.0)
})

test_that("cf_fit_hdps returns IPW weights and propensity scores", {
  skip_if_not_installed("glmnet")

  set.seed(456)
  n <- 150
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, 0.5)
  Y <- 2.0 * T + X[,1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "hdps")

  # Check that propensity scores are stored
  expect_true("propensity" %in% names(fit$fit$res))
  expect_length(fit$fit$res$propensity, n)
  expect_true(all(fit$fit$res$propensity > 0 & fit$fit$res$propensity < 1))

  # Check that weights are stored
  expect_true("weights" %in% names(fit$fit$res))
  expect_length(fit$fit$res$weights, n)
})

test_that("cf_fit_hdps respects trim parameter", {
  skip_if_not_installed("glmnet")

  set.seed(789)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  # Create strong confounding for extreme propensity scores
  beta_t <- c(rep(1.5, 3), rep(0, p - 3))
  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)
  Y <- 2.0 * T + X[,1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  # Fit with tight trimming
  fit <- cf_fit(fml, data = data, method = "hdps", trim = c(0.1, 0.9))

  # Propensity scores should be trimmed
  ps_hat <- fit$fit$res$propensity
  expect_true(all(ps_hat >= 0.1 & ps_hat <= 0.9))
})
