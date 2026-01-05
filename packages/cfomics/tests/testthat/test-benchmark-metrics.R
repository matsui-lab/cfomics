test_that("cf_benchmark_compute_metrics computes bias correctly", {
  truth <- list(
    ate_true = 1.0,
    ite_true = rep(1.0, 10)
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 1.2,
    ite_hat = rep(1.0, 10),
    summary_hat = NULL,
    truth = truth
  )
  
  expect_equal(metrics$bias_ate, 0.2)
  expect_equal(metrics$abs_bias_ate, 0.2)
  expect_equal(metrics$mse_ate, 0.04)
})

test_that("cf_benchmark_compute_metrics computes PEHE correctly", {
  truth <- list(
    ate_true = 1.0,
    ite_true = c(1.0, 2.0, 3.0, 4.0)
  )
  
  ite_hat <- c(1.5, 2.5, 3.5, 4.5)
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 2.5,
    ite_hat = ite_hat,
    summary_hat = NULL,
    truth = truth
  )
  
  expected_pehe <- sqrt(mean((ite_hat - truth$ite_true)^2))
  expect_equal(metrics$pehe, expected_pehe)
  expect_equal(metrics$pehe, 0.5)
})

test_that("cf_benchmark_compute_metrics handles perfect estimation", {
  truth <- list(
    ate_true = 2.0,
    ite_true = c(1.5, 2.0, 2.5)
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 2.0,
    ite_hat = c(1.5, 2.0, 2.5),
    summary_hat = NULL,
    truth = truth
  )
  
  expect_equal(metrics$bias_ate, 0)
  expect_equal(metrics$abs_bias_ate, 0)
  expect_equal(metrics$mse_ate, 0)
  expect_equal(metrics$pehe, 0)
})

test_that("cf_benchmark_compute_metrics computes coverage when CI provided", {
  truth <- list(
    ate_true = 1.0,
    ite_true = rep(1.0, 5)
  )
  
  summary_with_coverage <- list(
    ate_ci_lower = 0.8,
    ate_ci_upper = 1.2
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 1.0,
    ite_hat = rep(1.0, 5),
    summary_hat = summary_with_coverage,
    truth = truth
  )
  
  expect_equal(metrics$coverage_ate, 1)
  expect_equal(metrics$ci_len_ate, 0.4)
})

test_that("cf_benchmark_compute_metrics coverage is 0 when true ATE outside CI", {
  truth <- list(
    ate_true = 2.0,
    ite_true = rep(2.0, 5)
  )
  
  summary_no_coverage <- list(
    ate_ci_lower = 0.5,
    ate_ci_upper = 1.5
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 1.0,
    ite_hat = rep(1.0, 5),
    summary_hat = summary_no_coverage,
    truth = truth
  )
  
  expect_equal(metrics$coverage_ate, 0)
  expect_equal(metrics$ci_len_ate, 1.0)
})

test_that("cf_benchmark_compute_metrics returns NA for coverage when no CI", {
  truth <- list(
    ate_true = 1.0,
    ite_true = rep(1.0, 5)
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 1.0,
    ite_hat = rep(1.0, 5),
    summary_hat = NULL,
    truth = truth
  )
  
  expect_true(is.na(metrics$coverage_ate))
  expect_true(is.na(metrics$ci_len_ate))
})

test_that("cf_benchmark_compute_metrics handles NULL ate_hat", {
  truth <- list(
    ate_true = 1.0,
    ite_true = rep(1.0, 5)
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = NULL,
    ite_hat = rep(1.0, 5),
    summary_hat = NULL,
    truth = truth
  )
  
  expect_true(is.na(metrics$bias_ate))
  expect_true(is.na(metrics$abs_bias_ate))
  expect_true(is.na(metrics$mse_ate))
})

test_that("cf_benchmark_compute_metrics handles NULL ite_hat", {
  truth <- list(
    ate_true = 1.0,
    ite_true = rep(1.0, 5)
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 1.0,
    ite_hat = NULL,
    summary_hat = NULL,
    truth = truth
  )
  
  expect_true(is.na(metrics$pehe))
  expect_equal(metrics$bias_ate, 0)
})

test_that("cf_benchmark_compute_metrics handles mismatched ite lengths", {
  truth <- list(
    ate_true = 1.0,
    ite_true = rep(1.0, 10)
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 1.0,
    ite_hat = rep(1.0, 5),
    summary_hat = NULL,
    truth = truth
  )
  
  expect_true(is.na(metrics$pehe))
})

test_that("cf_benchmark_compute_metrics returns all numeric values", {
  truth <- list(
    ate_true = 1.5,
    ite_true = c(1.0, 1.5, 2.0)
  )
  
  summary_hat <- list(
    ate_ci_lower = 1.0,
    ate_ci_upper = 2.0
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 1.6,
    ite_hat = c(1.1, 1.6, 2.1),
    summary_hat = summary_hat,
    truth = truth
  )
  
  expect_type(metrics$bias_ate, "double")
  expect_type(metrics$abs_bias_ate, "double")
  expect_type(metrics$mse_ate, "double")
  expect_type(metrics$pehe, "double")
  expect_type(metrics$coverage_ate, "double")
  expect_type(metrics$ci_len_ate, "double")
})

test_that("cf_benchmark_compute_metrics handles negative bias", {
  truth <- list(
    ate_true = 2.0,
    ite_true = rep(2.0, 5)
  )
  
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = 1.5,
    ite_hat = rep(1.5, 5),
    summary_hat = NULL,
    truth = truth
  )
  
  expect_equal(metrics$bias_ate, -0.5)
  expect_equal(metrics$abs_bias_ate, 0.5)
  expect_equal(metrics$mse_ate, 0.25)
})
