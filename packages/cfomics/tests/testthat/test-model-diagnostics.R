test_that("cf_model_diagnostics returns valid diagnostics", {
  set.seed(42)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  ps <- plogis(0.5 * X$X1 - 0.3 * X$X2)
  T_var <- rbinom(n, 1, ps)
  Y <- 1.0 * X$X1 + 2.0 * T_var + rnorm(n)
  df <- data.frame(Y = Y, T = T_var, X)

  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 2.0,
      ite = rep(2.0, n),
      y0_hat = Y - 2.0 * T_var,
      y1_hat = Y - 2.0 * T_var + 2.0,
      summary = list(ate = 2.0, ate_ci_lower = 1.5, ate_ci_upper = 2.5)
    )),
    meta = list(formula = Y ~ T | X1 + X2, n = as.integer(n), p = 2L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = c("X1", "X2"))
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  diag <- cf_model_diagnostics(mock_result, df, Y ~ T | X1 + X2)

  expect_s3_class(diag, "cf_model_diagnostics")
  expect_true(is.numeric(diag$ps_cstat))
  expect_true(diag$ps_cstat > 0.5 && diag$ps_cstat < 1.0)
  expect_true(is.numeric(diag$extreme_weight_pct))
  expect_true(is.data.frame(diag$residual_summary))
})

test_that("cf_model_diagnostics handles methods without PS", {
  n <- 100
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 1.0, ite = rep(1.0, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 1),
      summary = list(ate = 1.0, ate_ci_lower = 0.5, ate_ci_upper = 1.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  diag <- cf_model_diagnostics(mock_result, df, Y ~ T | X1)

  # Should still work -- PS is estimated internally for diagnostics
  expect_s3_class(diag, "cf_model_diagnostics")
})

test_that("cf_model_diagnostics print works", {
  n <- 100
  mock_result <- list(
    method = "hdml",
    fit = list(res = list(
      ate = 1.0, ite = rep(1.0, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 1),
      summary = list(ate = 1.0, ate_ci_lower = 0.5, ate_ci_upper = 1.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  diag <- cf_model_diagnostics(mock_result, df, Y ~ T | X1)
  expect_output(print(diag), "Model Diagnostics")
})
