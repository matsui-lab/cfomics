# benchmarks/external/wrappers/tests/test_doubleml.R
# Test DoubleML wrapper

test_that("wrapper_doubleml returns correct structure", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")

  source(here::here("benchmarks/external/wrappers/wrapper_doubleml.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_doubleml(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1)
})

test_that("wrapper_doubleml returns confidence intervals", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")

  source(here::here("benchmarks/external/wrappers/wrapper_doubleml.R"))

  set.seed(456)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.3 * X[,1] + 0.2 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- 0.5 * X[,1] + 1.5 * T + rnorm(n)

  result <- wrapper_doubleml(X, T, Y)

  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(result$ci_lower < result$ate)
  expect_true(result$ci_upper > result$ate)
})

test_that("wrapper_doubleml returns counterfactual placeholders", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")

  source(here::here("benchmarks/external/wrappers/wrapper_doubleml.R"))

  set.seed(789)
  n <- 150
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.4 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2.5 * T + rnorm(n)

  result <- wrapper_doubleml(X, T, Y)

  # DoubleML PLR doesn't provide counterfactual predictions, so they should be NA
  expect_true("y0_hat" %in% names(result))
  expect_true("y1_hat" %in% names(result))
  expect_length(result$y0_hat, n)
  expect_length(result$y1_hat, n)
  expect_true(all(is.na(result$y0_hat)))
  expect_true(all(is.na(result$y1_hat)))
})

test_that("wrapper_doubleml handles unnamed covariate matrix", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")

  source(here::here("benchmarks/external/wrappers/wrapper_doubleml.R"))

  set.seed(101)
  n <- 100
  X <- matrix(rnorm(n * 5), n, 5)
  # No column names
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  # Should work without error
  result <- wrapper_doubleml(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_length(result$ite, n)
})

test_that("wrapper_doubleml uses constant ITE", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")
  skip_if_not_installed("ranger")

  source(here::here("benchmarks/external/wrappers/wrapper_doubleml.R"))

  set.seed(202)
  n <- 150
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_doubleml(X, T, Y)

  # DoubleML PLR only provides ATE, not heterogeneous effects
  # ITE should be constant (all equal to ATE)
  expect_true(all(result$ite == result$ate))
})

test_that("wrapper_doubleml supports different learners", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  source(here::here("benchmarks/external/wrappers/wrapper_doubleml.R"))

  set.seed(303)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  # Test with different learners (cv_glmnet for lasso)
  skip_if_not_installed("glmnet")
  result <- wrapper_doubleml(X, T, Y, ml_l = "regr.cv_glmnet", ml_m = "classif.cv_glmnet")

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_length(result$ite, n)
})
