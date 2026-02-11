test_that("cf_fit_grf basic execution", {
  skip_if_not_installed("grf")
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)
  
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")
  
  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "grf")
  
  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_true(ate > 0.5 && ate < 2.5)
  
  ite <- predict(fit, type = "ite")
  expect_equal(length(ite), n)
  
  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("ate_ci_lower" %in% names(summary_stats))
  expect_true("ate_ci_upper" %in% names(summary_stats))
})

test_that("cf_fit_ipw basic execution", {
  skip_if_not_installed("ipw")
  skip_if_not_installed("survey")
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)
  
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "ipw")
  
  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "ipw")
  
  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  
  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("ate_ci_lower" %in% names(summary_stats))
  expect_true("ate_ci_upper" %in% names(summary_stats))
})

test_that("cf_fit_gformula basic execution", {
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)
  
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula", nsimul = 100L)
  
  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "gformula")
  
  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_true(ate > 0.5 && ate < 2.5)
  
  ite <- predict(fit, type = "ite")
  expect_equal(length(ite), n)
  
  y0 <- predict(fit, type = "y0")
  y1 <- predict(fit, type = "y1")
  expect_equal(length(y0), n)
  expect_equal(length(y1), n)
  
  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("ate_ci_lower" %in% names(summary_stats))
  expect_true("ate_ci_upper" %in% names(summary_stats))
})

test_that("IPW returns NA for ITE with informative message", {
  skip_if_not_installed("ipw")
  skip_if_not_installed("survey")

  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1, data = data, method = "ipw")

  # ITE should be NA (IPW cannot estimate individual effects)
  expect_true(all(is.na(fit$fit$res$ite)))

  # ATE should still work (not NA)
  expect_false(is.na(fit$fit$res$ate))
})

test_that("gformula prediction on new data", {
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)
  
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula", nsimul = 50L)
  
  newdata <- data.frame(X1 = rnorm(10), X2 = rnorm(10))
  
  ite_new <- predict_cf_gformula(fit, newdata = newdata, type = "ite")
  expect_equal(length(ite_new), 10)
})
