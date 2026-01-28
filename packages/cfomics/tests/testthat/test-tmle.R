test_that("cf_fit_tmle basic execution", {
  skip_if_not_installed("glmnet")

  set.seed(123)
  n <- 200
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.5 * X[,1] + 0.3 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- 2.0 * T + 0.5 * X[,1] + 0.3 * X[,2] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  # Build formula
  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "tmle")

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "tmle")

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_true(ate > 1.0 && ate < 3.0)
})

test_that("cf_fit_tmle returns confidence intervals", {
  skip_if_not_installed("glmnet")

  set.seed(456)
  n <- 150
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, 0.5)
  Y <- 2.0 * T + X[,1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "tmle")

  summary_stats <- predict(fit, type = "summary")
  expect_true("ate_ci_lower" %in% names(summary_stats))
  expect_true("ate_ci_upper" %in% names(summary_stats))
  expect_true(summary_stats$ate_ci_lower < summary_stats$ate)
  expect_true(summary_stats$ate_ci_upper > summary_stats$ate)
})

test_that("cf_fit_tmle returns valid ITE estimates", {
  skip_if_not_installed("glmnet")

  set.seed(789)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, 0.5)
  Y <- 2.0 * T + X[,1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "tmle")

  ite <- predict(fit, type = "ite")
  expect_length(ite, n)
  expect_type(ite, "double")

  y0 <- predict(fit, type = "y0")
  y1 <- predict(fit, type = "y1")
  expect_length(y0, n)
  expect_length(y1, n)
  # ITE should be y1 - y0
  expect_equal(ite, y1 - y0)
})
