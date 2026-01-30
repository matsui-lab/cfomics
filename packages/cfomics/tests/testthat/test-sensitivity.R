test_that("cf_sensitivity computes E-value for a cfomics_result", {
  # Create minimal cfomics_result
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 1.5,
      ite = rep(1.5, 100),
      y0_hat = rnorm(100, 5, 1),
      y1_hat = rnorm(100, 6.5, 1),
      summary = list(ate = 1.5, ate_ci_lower = 1.0, ate_ci_upper = 2.0)
    )),
    meta = list(formula = Y ~ T | X, n = 100L, p = 5L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = paste0("X", 1:5))
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  result <- cf_sensitivity(mock_result)

  expect_s3_class(result, "cf_sensitivity")
  expect_true(is.numeric(result$evalue_point))
  expect_true(is.numeric(result$evalue_ci))
  expect_true(result$evalue_point > 1)
  expect_true(result$evalue_ci > 1)
  # E-value for point should be >= E-value for CI lower bound
  expect_gte(result$evalue_point, result$evalue_ci)
  expect_true(is.character(result$interpretation))
})

test_that("cf_sensitivity returns E-value = 1 when ATE is 0", {
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 0,
      ite = rep(0, 100),
      y0_hat = rnorm(100, 5, 1),
      y1_hat = rnorm(100, 5, 1),
      summary = list(ate = 0, ate_ci_lower = -0.5, ate_ci_upper = 0.5)
    )),
    meta = list(formula = Y ~ T | X, n = 100L, p = 5L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = paste0("X", 1:5))
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  result <- cf_sensitivity(mock_result)
  expect_equal(result$evalue_point, 1)
})

test_that("cf_sensitivity print method works", {
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 1.5,
      ite = rep(1.5, 100),
      y0_hat = rnorm(100, 5, 1),
      y1_hat = rnorm(100, 6.5, 1),
      summary = list(ate = 1.5, ate_ci_lower = 1.0, ate_ci_upper = 2.0)
    )),
    meta = list(formula = Y ~ T | X, n = 100L, p = 5L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = paste0("X", 1:5))
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  sens <- cf_sensitivity(mock_result)
  expect_output(print(sens), "E-value")
})
