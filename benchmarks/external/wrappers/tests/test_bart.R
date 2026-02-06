# benchmarks/external/wrappers/tests/test_bart.R
# Test BART (Bayesian Additive Regression Trees) wrapper

test_that("wrapper_bart returns correct structure", {
  skip_if_not_installed("dbarts")

  source(here::here("benchmarks/external/wrappers/wrapper_bart.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_bart(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1)
})

test_that("wrapper_bart returns confidence intervals", {
  skip_if_not_installed("dbarts")

  source(here::here("benchmarks/external/wrappers/wrapper_bart.R"))

  set.seed(456)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.3 * X[,1] + 0.2 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- 0.5 * X[,1] + 1.5 * T + rnorm(n)

  result <- wrapper_bart(X, T, Y)

  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(result$ci_lower < result$ate)
  expect_true(result$ci_upper > result$ate)
})

test_that("wrapper_bart returns real counterfactual predictions (NOT NA)", {
  skip_if_not_installed("dbarts")

  source(here::here("benchmarks/external/wrappers/wrapper_bart.R"))

  set.seed(789)
  n <- 150
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.4 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2.5 * T + rnorm(n)

  result <- wrapper_bart(X, T, Y)

  # BART provides REAL counterfactual predictions (unlike hdm/DoubleML/tmle3)
  expect_true("y0_hat" %in% names(result))
  expect_true("y1_hat" %in% names(result))
  expect_length(result$y0_hat, n)
  expect_length(result$y1_hat, n)
  expect_false(any(is.na(result$y0_hat)))
  expect_false(any(is.na(result$y1_hat)))

  # ITE should be y1_hat - y0_hat
  expect_equal(result$ite, result$y1_hat - result$y0_hat)
})

test_that("wrapper_bart handles unnamed covariate matrix", {
  skip_if_not_installed("dbarts")

  source(here::here("benchmarks/external/wrappers/wrapper_bart.R"))

  set.seed(101)
  n <- 100
  X <- matrix(rnorm(n * 5), n, 5)
  # No column names
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  # Should work without error
  result <- wrapper_bart(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_length(result$ite, n)
})

test_that("wrapper_bart provides heterogeneous ITE (NOT constant)", {
  skip_if_not_installed("dbarts")

  source(here::here("benchmarks/external/wrappers/wrapper_bart.R"))

  set.seed(202)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  # Heterogeneous treatment effect: effect depends on X1
  Y <- X[,1] + (2 + 0.5 * X[,1]) * T + rnorm(n)

  result <- wrapper_bart(X, T, Y)

  # BART provides heterogeneous effects (unlike hdm/DoubleML/tmle3)
  # ITE should NOT be constant (should have variance)
  expect_true(sd(result$ite) > 0)
  # Check that ITE values are different (not all equal to ATE)
  expect_false(all(result$ite == result$ate))
})

test_that("wrapper_bart supports MCMC parameters", {
  skip_if_not_installed("dbarts")

  source(here::here("benchmarks/external/wrappers/wrapper_bart.R"))

  set.seed(303)
  n <- 150
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  # Test with custom MCMC parameters
  result <- wrapper_bart(
    X, T, Y,
    n_trees = 100L,
    n_burn = 100L,
    n_iter = 500L
  )

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_length(result$ite, n)
  # Should still get reasonable estimate
  expect_true(abs(result$ate - 2) < 1.5)
})

test_that("wrapper_bart CI contains true effect", {
  skip_if_not_installed("dbarts")

  source(here::here("benchmarks/external/wrappers/wrapper_bart.R"))

  set.seed(404)
  n <- 300
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  true_ate <- 2
  Y <- X[,1] + true_ate * T + rnorm(n)

  result <- wrapper_bart(X, T, Y)

  # The 95% CI should contain the true effect
  expect_true(result$ci_lower < true_ate)
  expect_true(result$ci_upper > true_ate)
})
