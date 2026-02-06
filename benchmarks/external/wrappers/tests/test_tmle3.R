# benchmarks/external/wrappers/tests/test_tmle3.R
# Test tmle3 wrapper

test_that("wrapper_tmle3 returns correct structure", {
  skip_if_not_installed("tmle3")
  skip_if_not_installed("sl3")

  source(here::here("benchmarks/external/wrappers/wrapper_tmle3.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_tmle3(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1)
})

test_that("wrapper_tmle3 returns confidence intervals", {
  skip_if_not_installed("tmle3")
  skip_if_not_installed("sl3")

  source(here::here("benchmarks/external/wrappers/wrapper_tmle3.R"))

  set.seed(456)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.3 * X[,1] + 0.2 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- 0.5 * X[,1] + 1.5 * T + rnorm(n)

  result <- wrapper_tmle3(X, T, Y)

  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(result$ci_lower < result$ate)
  expect_true(result$ci_upper > result$ate)
})

test_that("wrapper_tmle3 returns counterfactual placeholders", {
  skip_if_not_installed("tmle3")
  skip_if_not_installed("sl3")

  source(here::here("benchmarks/external/wrappers/wrapper_tmle3.R"))

  set.seed(789)
  n <- 150
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.4 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2.5 * T + rnorm(n)

  result <- wrapper_tmle3(X, T, Y)

  # tmle3 doesn't provide individual counterfactual predictions, so they should be NA
  expect_true("y0_hat" %in% names(result))
  expect_true("y1_hat" %in% names(result))
  expect_length(result$y0_hat, n)
  expect_length(result$y1_hat, n)
  expect_true(all(is.na(result$y0_hat)))
  expect_true(all(is.na(result$y1_hat)))
})

test_that("wrapper_tmle3 handles unnamed covariate matrix", {
  skip_if_not_installed("tmle3")
  skip_if_not_installed("sl3")

  source(here::here("benchmarks/external/wrappers/wrapper_tmle3.R"))

  set.seed(101)
  n <- 100
  X <- matrix(rnorm(n * 5), n, 5)
  # No column names
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  # Should work without error
  result <- wrapper_tmle3(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_length(result$ite, n)
})

test_that("wrapper_tmle3 uses constant ITE", {
  skip_if_not_installed("tmle3")
  skip_if_not_installed("sl3")

  source(here::here("benchmarks/external/wrappers/wrapper_tmle3.R"))

  set.seed(202)
  n <- 150
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_tmle3(X, T, Y)

  # TMLE provides population-level ATE, not heterogeneous individual effects
  # ITE should be constant (all equal to ATE)
  expect_true(all(result$ite == result$ate))
})

test_that("wrapper_tmle3 estimates reasonable ATE", {
  skip_if_not_installed("tmle3")
  skip_if_not_installed("sl3")

  source(here::here("benchmarks/external/wrappers/wrapper_tmle3.R"))

  set.seed(303)
  n <- 300
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  # True ATE = 3
  Y <- X[,1] + 3 * T + rnorm(n)

  result <- wrapper_tmle3(X, T, Y)

  # With larger sample, should be closer to true effect
  expect_true(abs(result$ate - 3) < 1)
  # CI should contain true value
  expect_true(result$ci_lower < 3 && result$ci_upper > 3)
})
