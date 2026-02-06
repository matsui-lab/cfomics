test_that("cf_benchmark_run_once returns correct structure", {
  result <- cfomics:::cf_benchmark_run_once(
    scenario = "linear_homogeneous",
    method = "gformula",
    n = 50L,
    p = 5L,
    seed_dgp = 1L,
    formula = stats::as.formula("Y ~ T | X1 + X2"),
    effect_size = 1.0,
    noise_sd = 1.0,
    allow_python_env_install = FALSE,
    nsimul = 10L
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  
  expected_cols <- c(
    "scenario_id", "method", "replicate_id", "n", "p",
    "seed_dgp", "seed_method", "ate_true", "ate_hat",
    "bias_ate", "squared_error_ate", "pehe", "coverage_ate", "ci_len_ate",
    "time_fit_sec", "time_predict_sec", "status", "error_message"
  )
  
  for (col in expected_cols) {
    expect_true(col %in% names(result), info = sprintf("Missing column: %s", col))
  }
})

test_that("cf_benchmark_run_once gformula completes successfully", {
  result <- cfomics:::cf_benchmark_run_once(
    scenario = "linear_homogeneous",
    method = "gformula",
    n = 50L,
    p = 5L,
    seed_dgp = 42L,
    formula = stats::as.formula("Y ~ T | X1 + X2"),
    effect_size = 1.0,
    noise_sd = 1.0,
    allow_python_env_install = FALSE,
    nsimul = 10L
  )
  
  expect_equal(result$status, "ok")
  expect_true(is.numeric(result$ate_hat))
  expect_true(is.numeric(result$bias_ate))
  expect_true(is.numeric(result$pehe))
  expect_true(result$time_fit_sec > 0)
})

test_that("cf_benchmark_run_once skips Python methods when dependencies unavailable", {
  result <- cfomics:::cf_benchmark_run_once(
    scenario = "linear_homogeneous",
    method = "dowhy_gcm",
    n = 50L,
    p = 5L,
    seed_dgp = 1L,
    allow_python_env_install = FALSE
  )
  
  expect_true(result$status %in% c("ok", "skipped", "error"))
  
  if (result$status == "skipped") {
    expect_true(!is.na(result$error_message))
  }
})

test_that("cf_benchmark_run returns data.frame with correct dimensions", {
  skip_on_cran()
  
  results <- cfomics:::cf_benchmark_run(
    scenarios = c("linear_homogeneous"),
    methods = c("gformula"),
    n = 30L,
    p = 5L,
    n_reps = 2L,
    base_seed = 1L,
    covariates_k = 2L,
    allow_python_env_install = FALSE,
    nsimul = 5L
  )
  
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 2)
})

test_that("cf_benchmark_summarize produces correct summary", {
  skip_on_cran()
  
  results <- cfomics:::cf_benchmark_run(
    scenarios = c("linear_homogeneous"),
    methods = c("gformula"),
    n = 30L,
    p = 5L,
    n_reps = 3L,
    base_seed = 1L,
    covariates_k = 2L,
    allow_python_env_install = FALSE,
    nsimul = 5L
  )
  
  summary_df <- cfomics:::cf_benchmark_summarize(results)
  
  expect_s3_class(summary_df, "data.frame")
  expect_equal(nrow(summary_df), 1)
  
  expected_cols <- c(
    "scenario_id", "method", "rmse_ate", "mean_bias_ate", "mean_pehe",
    "coverage_rate", "median_time_fit_sec", "median_time_predict_sec",
    "n_ok", "n_skipped", "n_error"
  )
  
  for (col in expected_cols) {
    expect_true(col %in% names(summary_df), info = sprintf("Missing column: %s", col))
  }
  
  expect_equal(summary_df$n_ok, 3)
})

test_that("cf_benchmark_run handles multiple scenarios and methods", {
  skip_on_cran()
  
  results <- cfomics:::cf_benchmark_run(
    scenarios = c("linear_homogeneous", "heterogeneous_ite"),
    methods = c("gformula"),
    n = 30L,
    p = 5L,
    n_reps = 2L,
    base_seed = 1L,
    covariates_k = 2L,
    allow_python_env_install = FALSE,
    nsimul = 5L
  )
  
  expect_equal(nrow(results), 4)
  expect_equal(length(unique(results$scenario_id)), 2)
})

test_that("cf_benchmark_run with grf method", {
  skip_if_not_installed("grf")
  skip_on_cran()
  
  results <- cfomics:::cf_benchmark_run(
    scenarios = c("linear_homogeneous"),
    methods = c("grf"),
    n = 200L,
    p = 5L,
    n_reps = 2L,
    base_seed = 1L,
    covariates_k = 2L,
    allow_python_env_install = FALSE,
    num.trees = 100L
  )
  
  expect_equal(nrow(results), 2)
  expect_true(all(results$status == "ok"))
})

test_that("cf_benchmark_run with ipw method", {
  skip_if_not_installed("ipw")
  skip_if_not_installed("survey")
  skip_on_cran()
  
  results <- cfomics:::cf_benchmark_run(
    scenarios = c("linear_homogeneous"),
    methods = c("ipw"),
    n = 50L,
    p = 5L,
    n_reps = 2L,
    base_seed = 1L,
    covariates_k = 2L,
    allow_python_env_install = FALSE
  )
  
  expect_equal(nrow(results), 2)
  expect_true(all(results$status == "ok"))
})

test_that("cf_benchmark_plot handles missing ggplot2 gracefully", {
  summary_df <- data.frame(
    scenario_id = "linear_homogeneous",
    method = "gformula",
    rmse_ate = 0.1,
    mean_bias_ate = 0.05,
    mean_pehe = 0.2,
    coverage_rate = 0.95,
    median_time_fit_sec = 1.0,
    median_time_predict_sec = 0.1,
    n_ok = 10,
    n_skipped = 0,
    n_error = 0,
    stringsAsFactors = FALSE
  )
  
  result <- cfomics:::cf_benchmark_plot(summary_df, out_dir = NULL)
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_type(result, "list")
  } else {
    expect_null(result)
  }
})

test_that("covariates_k is bounded by p", {
  results <- cfomics:::cf_benchmark_run(
    scenarios = c("linear_homogeneous"),
    methods = c("gformula"),
    n = 30L,
    p = 3L,
    n_reps = 1L,
    base_seed = 1L,
    covariates_k = 10L,
    allow_python_env_install = FALSE,
    nsimul = 5L
  )
  
  expect_equal(nrow(results), 1)
  expect_equal(results$status, "ok")
})

test_that(".check_method_dependencies returns NULL for available R methods", {
  result <- cfomics:::.check_method_dependencies("gformula", FALSE)
  expect_null(result)
})

test_that(".check_method_dependencies returns skip reason for missing packages", {
  skip_if(requireNamespace("grf", quietly = TRUE))
  
  result <- cfomics:::.check_method_dependencies("grf", FALSE)
  expect_true(!is.null(result))
  expect_true(grepl("grf", result))
})
