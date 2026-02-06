# benchmarks/external/wrappers/tests/test_hdm.R
# Test hdm (high-dimensional metrics) wrapper

test_that("wrapper_hdm returns correct structure", {
  skip_if_not_installed("hdm")

  source(here::here("benchmarks/external/wrappers/wrapper_hdm.R"))

  set.seed(123)
  n <- 200
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.5 * X[,1] + 0.3 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_hdm(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1.5)
})

test_that("wrapper_hdm returns confidence intervals", {
  skip_if_not_installed("hdm")

  source(here::here("benchmarks/external/wrappers/wrapper_hdm.R"))

  set.seed(456)
  n <- 200
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.3 * X[,1] + 0.2 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- 0.5 * X[,1] + 1.5 * T + rnorm(n)

  result <- wrapper_hdm(X, T, Y)

  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(result$ci_lower < result$ate)
  expect_true(result$ci_upper > result$ate)
})

test_that("wrapper_hdm returns counterfactual placeholders", {
  skip_if_not_installed("hdm")

  source(here::here("benchmarks/external/wrappers/wrapper_hdm.R"))

  set.seed(789)
  n <- 150
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.4 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2.5 * T + rnorm(n)

  result <- wrapper_hdm(X, T, Y)

  # hdm doesn't provide counterfactual predictions, so they should be NA
  expect_true("y0_hat" %in% names(result))
  expect_true("y1_hat" %in% names(result))
  expect_length(result$y0_hat, n)
  expect_length(result$y1_hat, n)
  expect_true(all(is.na(result$y0_hat)))
  expect_true(all(is.na(result$y1_hat)))
})

test_that("wrapper_hdm handles unnamed covariate matrix", {
  skip_if_not_installed("hdm")

  source(here::here("benchmarks/external/wrappers/wrapper_hdm.R"))

  set.seed(101)
  n <- 100
  p <- 25
  X <- matrix(rnorm(n * p), n, p)
  # No column names
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  # Should work without error
  result <- wrapper_hdm(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_length(result$ite, n)
})

test_that("wrapper_hdm uses constant ITE", {
  skip_if_not_installed("hdm")

  source(here::here("benchmarks/external/wrappers/wrapper_hdm.R"))

  set.seed(202)
  n <- 150
  p <- 40
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_hdm(X, T, Y)

  # hdm only provides ATE, not heterogeneous effects

# ITE should be constant (all equal to ATE)
  expect_true(all(result$ite == result$ate))
})

test_that("wrapper_hdm works with high-dimensional data (p > n)", {
  skip_if_not_installed("hdm")

  source(here::here("benchmarks/external/wrappers/wrapper_hdm.R"))

  set.seed(303)
  n <- 100
  p <- 150  # More covariates than observations
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  # Only first few covariates are relevant
  ps <- plogis(0.5 * X[,1] + 0.3 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 0.5 * X[,2] + 2 * T + rnorm(n)

  result <- wrapper_hdm(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_length(result$ite, n)
  # Should still get reasonable estimate despite p > n
  expect_true(abs(result$ate - 2) < 2)
})

test_that("wrapper_hdm supports different methods", {
  skip_if_not_installed("hdm")

  source(here::here("benchmarks/external/wrappers/wrapper_hdm.R"))

  set.seed(404)
  n <- 200
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  # Test with partialling out method
  result_po <- wrapper_hdm(X, T, Y, method = "partialling out")

  expect_type(result_po, "list")
  expect_true("ate" %in% names(result_po))
  expect_length(result_po$ite, n)

  # Test with double selection method (default)
  result_ds <- wrapper_hdm(X, T, Y, method = "double selection")

  expect_type(result_ds, "list")
  expect_true("ate" %in% names(result_ds))
})
