test_that("cf_cate_diagnostics works with heterogeneous ITE", {
  set.seed(42)
  n <- 500
  X1 <- rnorm(n)
  T_var <- rbinom(n, 1, 0.5)
  ite <- 1.0 + 2.0 * X1  # heterogeneous
  Y <- X1 + ite * T_var + rnorm(n)
  df <- data.frame(Y = Y, T = T_var, X1 = X1)

  mock_result <- list(
    method = "grf",
    fit = list(res = list(
      ate = mean(ite), ite = ite,
      y0_hat = X1, y1_hat = X1 + ite,
      summary = list(ate = mean(ite), ate_ci_lower = 0.5, ate_ci_upper = 1.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  cate <- cf_cate_diagnostics(mock_result, df, Y ~ T | X1)

  expect_s3_class(cate, "cf_cate_diagnostics")
  expect_true(!is.null(cate$gates))
  expect_equal(nrow(cate$gates), 4)  # quartiles
  expect_true(!is.null(cate$blp))
  expect_true(cate$heterogeneity_detected)
})

test_that("cf_cate_diagnostics skips GATES/BLP for constant ITE", {
  n <- 200
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 2.0, ite = rep(2.0, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 2),
      summary = list(ate = 2.0, ate_ci_lower = 1.5, ate_ci_upper = 2.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  cate <- cf_cate_diagnostics(mock_result, df, Y ~ T | X1)

  expect_s3_class(cate, "cf_cate_diagnostics")
  expect_true(is.null(cate$gates))
  expect_true(is.null(cate$blp))
  expect_false(cate$heterogeneity_detected)
  expect_true(cate$constant_ite)
})

test_that("cf_cate_diagnostics print works", {
  n <- 200
  mock_result <- list(
    method = "grf",
    fit = list(res = list(
      ate = 2.0, ite = rnorm(n, 2, 1),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 2),
      summary = list(ate = 2.0, ate_ci_lower = 1.5, ate_ci_upper = 2.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  cate <- cf_cate_diagnostics(mock_result, df, Y ~ T | X1)
  expect_output(print(cate), "CATE")
})
