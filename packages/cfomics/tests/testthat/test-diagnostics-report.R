test_that("cf_diagnostics_report runs all diagnostics", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n)
  T_var <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- X1 + 2.0 * T_var + rnorm(n)
  df <- data.frame(Y = Y, T = T_var, X1 = X1)

  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 2.0, ite = rep(2.0, n),
      y0_hat = X1, y1_hat = X1 + 2.0,
      summary = list(ate = 2.0, ate_ci_lower = 1.5, ate_ci_upper = 2.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  report <- cf_diagnostics_report(mock_result, df, Y ~ T | X1)

  expect_s3_class(report, "cf_diagnostics_report")
  expect_s3_class(report$sensitivity, "cf_sensitivity")
  expect_s3_class(report$model, "cf_model_diagnostics")
  expect_s3_class(report$cate, "cf_cate_diagnostics")
  expect_s3_class(report$balance, "cf_balance")
  expect_true(is.data.frame(report$summary))
})

test_that("cf_diagnostics_report respects include parameter", {
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
  report <- cf_diagnostics_report(mock_result, df, Y ~ T | X1,
                                   include = "sensitivity")

  expect_s3_class(report$sensitivity, "cf_sensitivity")
  expect_null(report$model)
  expect_null(report$cate)
})

test_that("cf_diagnostics_report print works", {
  n <- 100
  mock_result <- list(
    method = "hdml",
    fit = list(res = list(
      ate = 1.5, ite = rep(1.5, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 1.5),
      summary = list(ate = 1.5, ate_ci_lower = 1.0, ate_ci_upper = 2.0)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  report <- cf_diagnostics_report(mock_result, df, Y ~ T | X1)
  expect_output(print(report), "Diagnostics Report")
})
