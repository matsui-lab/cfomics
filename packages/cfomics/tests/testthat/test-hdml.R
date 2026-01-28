test_that("cf_fit_hdml basic execution", {
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

  # Build formula with explicit covariate names (cfomics doesn't expand ".")
  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "hdml")

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "hdml")

  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_true(ate > 1.0 && ate < 3.0)

  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("ate_ci_lower" %in% names(summary_stats))
  expect_true("ate_ci_upper" %in% names(summary_stats))
})

test_that("cf_fit_hdml supports different penalty types", {
  skip_if_not_installed("glmnet")

  set.seed(456)
  n <- 150
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  beta_t <- c(rep(0.3, 5), rep(0, p - 5))
  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)
  Y <- 2.0 * T + 0.5 * X[,1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  # Build formula with explicit covariate names
  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  # Test ridge
  fit_ridge <- cf_fit(fml, data = data, method = "hdml", penalty = "ridge")
  expect_s3_class(fit_ridge, "cf_model")

  # Test elastic_net
  fit_enet <- cf_fit(fml, data = data, method = "hdml", penalty = "elastic_net")
  expect_s3_class(fit_enet, "cf_model")
})

test_that("cf_fit_hdml returns valid ITE estimates", {
  skip_if_not_installed("glmnet")

  set.seed(789)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  T <- rbinom(n, 1, 0.5)
  Y <- 2.0 * T + X[,1] + rnorm(n)
  data <- data.frame(Y = Y, T = T, X)

  # Build formula with explicit covariate names
  cov_names <- paste0("X", 1:p)
  fml <- as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  fit <- cf_fit(fml, data = data, method = "hdml")

  ite <- predict(fit, type = "ite")
  expect_length(ite, n)
  expect_type(ite, "double")

  y0 <- predict(fit, type = "y0")
  y1 <- predict(fit, type = "y1")
  expect_length(y0, n)
  expect_length(y1, n)
  expect_equal(ite, y1 - y0)
})
