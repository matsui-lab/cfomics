# benchmarks/external/wrappers/tests/test_weightit.R
# Test WeightIt wrapper

test_that("wrapper_weightit returns correct structure", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("marginaleffects")

  source(here::here("benchmarks/external/wrappers/wrapper_weightit.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_weightit(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(is.numeric(result$ate))
  expect_true(abs(result$ate - 2) < 1)  # Rough check
})

test_that("wrapper_weightit handles different weighting methods", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("marginaleffects")

  source(here::here("benchmarks/external/wrappers/wrapper_weightit.R"))

  set.seed(456)
  n <- 200
  X <- matrix(rnorm(n * 3), n, 3)
  colnames(X) <- paste0("X", 1:3)
  ps <- plogis(0.3 * X[,1] + 0.2 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- 0.5 * X[,1] + 1.5 * T + rnorm(n)

  # Test with default GLM method
  result <- wrapper_weightit(X, T, Y, method = "glm", estimand = "ATE")

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
})

test_that("wrapper_weightit returns counterfactual predictions", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("marginaleffects")

  source(here::here("benchmarks/external/wrappers/wrapper_weightit.R"))

  set.seed(789)
  n <- 150
  X <- matrix(rnorm(n * 4), n, 4)
  colnames(X) <- paste0("X", 1:4)
  ps <- plogis(0.4 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2.5 * T + rnorm(n)

  result <- wrapper_weightit(X, T, Y)

  # Check counterfactual predictions exist
  expect_true("y0_hat" %in% names(result))
  expect_true("y1_hat" %in% names(result))
  expect_length(result$y0_hat, n)
  expect_length(result$y1_hat, n)

  # ITE should be y1_hat - y0_hat
  expect_equal(result$ite, result$y1_hat - result$y0_hat)
})

test_that("wrapper_weightit returns confidence intervals", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("marginaleffects")

  source(here::here("benchmarks/external/wrappers/wrapper_weightit.R"))

  set.seed(101)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_weightit(X, T, Y)

  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(result$ci_lower < result$ate)
  expect_true(result$ci_upper > result$ate)
})

test_that("wrapper_weightit handles unnamed covariate matrix", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("marginaleffects")

  source(here::here("benchmarks/external/wrappers/wrapper_weightit.R"))

  set.seed(202)
  n <- 100
  X <- matrix(rnorm(n * 3), n, 3)
  # No column names
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  # Should work without error
  result <- wrapper_weightit(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_length(result$ite, n)
})

test_that("wrapper_weightit supports ATT estimand", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("marginaleffects")

  source(here::here("benchmarks/external/wrappers/wrapper_weightit.R"))

  set.seed(303)
  n <- 200
  X <- matrix(rnorm(n * 4), n, 4)
  colnames(X) <- paste0("X", 1:4)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_weightit(X, T, Y, estimand = "ATT")

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
})
