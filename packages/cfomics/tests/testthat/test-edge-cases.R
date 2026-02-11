# test-edge-cases.R - Edge case tests and input validation for cf_fit()

test_that("cf_fit handles single covariate", {
  set.seed(123)
  data <- data.frame(Y = rnorm(50), T = rbinom(50, 1, 0.5), X1 = rnorm(50))
  expect_no_error(cf_fit(Y ~ T | X1, data = data, method = "gformula"))
})

test_that("cf_fit handles small sample sizes (minimum 10)", {
  set.seed(123)
  # 20 observations should work
  data <- data.frame(Y = rnorm(20), T = c(rep(0, 10), rep(1, 10)), X1 = rnorm(20))
  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula")
  expect_s3_class(fit, "cfomics_result")
})

test_that("cf_fit rejects sample sizes below minimum", {
  set.seed(123)
  # 9 observations should fail (minimum is 10)
  data <- data.frame(Y = rnorm(9), T = c(rep(0, 5), rep(1, 4)), X1 = rnorm(9))
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "Insufficient sample size"
  )
})

test_that("cf_fit validates binary treatment", {
  set.seed(123)
  data <- data.frame(Y = rnorm(50), T = rnorm(50), X1 = rnorm(50))  # T not binary
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "binary"
  )
})

test_that("cf_fit validates numeric outcome", {
  set.seed(123)
  # Create a factor outcome - won't be converted to NA like character
  data <- data.frame(
    Y = factor(sample(c("low", "high"), 50, replace = TRUE)),
    T = rbinom(50, 1, 0.5),
    X1 = rnorm(50),
    stringsAsFactors = FALSE
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "numeric"
  )
})

test_that("cf_fit requires both treatment groups", {
  set.seed(123)
  # All treated
  data <- data.frame(Y = rnorm(50), T = rep(1, 50), X1 = rnorm(50))
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "treatment and control"
  )

  # All control
  data2 <- data.frame(Y = rnorm(50), T = rep(0, 50), X1 = rnorm(50))
  expect_error(
    cf_fit(Y ~ T | X1, data = data2, method = "gformula"),
    "treatment and control"
  )
})

test_that("cf_fit handles multiple covariates", {
  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n),
    X3 = rnorm(n),
    X4 = rnorm(n),
    X5 = rnorm(n)
  )
  fit <- cf_fit(Y ~ T | X1 + X2 + X3 + X4 + X5, data = data, method = "gformula")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$meta$p, 5)
})

test_that("cf_fit validates numeric covariates", {
  set.seed(123)
  # Create a factor covariate - won't be converted to NA
  data <- data.frame(
    Y = rnorm(50),
    T = rbinom(50, 1, 0.5),
    X1 = factor(sample(c("A", "B", "C"), 50, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "numeric covariates"
  )
})

test_that("cf_fit rejects NA values", {
  set.seed(123)
  # NA in outcome
  data <- data.frame(Y = c(NA, rnorm(49)), T = rbinom(50, 1, 0.5), X1 = rnorm(50))
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "NA values"
  )

  # NA in treatment
  data2 <- data.frame(Y = rnorm(50), T = c(NA, rbinom(49, 1, 0.5)), X1 = rnorm(50))
  expect_error(
    cf_fit(Y ~ T | X1, data = data2, method = "gformula"),
    "NA values"
  )

  # NA in covariate
  data3 <- data.frame(Y = rnorm(50), T = rbinom(50, 1, 0.5), X1 = c(NA, rnorm(49)))
  expect_error(
    cf_fit(Y ~ T | X1, data = data3, method = "gformula"),
    "NA values"
  )
})

test_that("cf_fit handles logical treatment (converts to integer)", {
  set.seed(123)
  data <- data.frame(
    Y = rnorm(50),
    T = sample(c(TRUE, FALSE), 50, replace = TRUE),
    X1 = rnorm(50)
  )
  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula")
  expect_s3_class(fit, "cfomics_result")
})

test_that("cf_fit handles highly imbalanced treatment", {
  set.seed(123)
  # 95% treated, 5% control (but still valid)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = c(rep(1, 95), rep(0, 5)),
    X1 = rnorm(n)
  )
  # Should work but may produce warnings about rank-deficient fit in bootstrap
  # Suppress expected warnings from lm() with imbalanced data
  fit <- suppressWarnings(cf_fit(Y ~ T | X1, data = data, method = "gformula"))
  expect_s3_class(fit, "cfomics_result")
})

test_that("cf_fit returns correct result structure", {
  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula")

  # Check class

  expect_s3_class(fit, "cfomics_result")
  expect_s3_class(fit, "cf_model")

  # Check structure

  expect_true("method" %in% names(fit))
  expect_true("fit" %in% names(fit))
  expect_true("meta" %in% names(fit))

  # Check meta fields
  expect_equal(fit$meta$n, n)
  expect_equal(fit$meta$p, 2)
  expect_equal(fit$meta$outcome_name, "Y")
  expect_equal(fit$meta$treatment_name, "T")
  expect_equal(fit$meta$covariate_names, c("X1", "X2"))

  # Check results
  expect_true("res" %in% names(fit$fit))
  expect_true("ate" %in% names(fit$fit$res))
  expect_true("ite" %in% names(fit$fit$res))
  expect_length(fit$fit$res$ite, n)
})

test_that("cf_fit handles empty data gracefully", {
  data <- data.frame(Y = numeric(0), T = integer(0), X1 = numeric(0))
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "no observations|at least 10"
  )
})

test_that("cf_fit validates variable name constraints", {
  set.seed(123)
  # Variable names with special characters should be rejected
  data <- data.frame(
    Y = rnorm(50),
    T = rbinom(50, 1, 0.5),
    `X.1` = rnorm(50),  # Contains period
    check.names = FALSE
  )
  expect_error(
    cf_fit(Y ~ T | `X.1`, data = data, method = "gformula"),
    "variable names"
  )
})
