# Tests for reporting functions (cf_table, cf_export, cf_session_info)

# Helper to create a fitted model for testing
create_test_fit <- function() {
  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  data$Y <- 1.5 * data$T + 0.5 * data$X1 + 0.3 * data$X2 + rnorm(n)

  cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula", nsimul = 50L)
}

# --- cf_table tests ---

test_that("cf_table returns data.frame by default", {
  fit <- create_test_fit()

  result <- cf_table(fit)

  expect_s3_class(result, "data.frame")
  expect_true("Metric" %in% names(result))
  expect_true("Value" %in% names(result))
})

test_that("cf_table includes required metrics", {
  fit <- create_test_fit()

  result <- cf_table(fit)
  metrics <- result$Metric

  expect_true("Method" %in% metrics)
  expect_true("Sample size (N)" %in% metrics)
  expect_true("Covariates (p)" %in% metrics)
  expect_true("ATE" %in% metrics)
})

test_that("cf_table includes confidence intervals when available", {
  fit <- create_test_fit()

  result <- cf_table(fit)
  metrics <- result$Metric

  expect_true("ATE 95% CI" %in% metrics)
})

test_that("cf_table respects digits parameter", {
  fit <- create_test_fit()

  result <- cf_table(fit, digits = 2)
  ate_value <- result$Value[result$Metric == "ATE"]

  # ATE should be formatted with 2 decimal places
  expect_true(grepl("^-?[0-9]+\\.[0-9]{1,2}$", as.character(ate_value)))
})

test_that("cf_table respects include_ite parameter", {
  fit <- create_test_fit()

  result_with_ite <- cf_table(fit, include_ite = TRUE)
  result_without_ite <- cf_table(fit, include_ite = FALSE)

  # With ITE should have more rows
  expect_gte(nrow(result_with_ite), nrow(result_without_ite))
})

test_that("cf_table rejects non-cf_model input", {
  expect_error(
    cf_table("not a model"),
    "cf_model"
  )

  expect_error(
    cf_table(list(a = 1)),
    "cf_model"
  )
})

test_that("cf_table returns markdown format", {
  skip_if_not_installed("knitr")

  fit <- create_test_fit()

  result <- cf_table(fit, format = "markdown")

  expect_true(is.character(result))
  expect_true(grepl("\\|", result[1]))  # markdown table uses |
})

test_that("cf_table returns latex format", {
  skip_if_not_installed("knitr")

  fit <- create_test_fit()

  result <- cf_table(fit, format = "latex")

  expect_true(is.character(result))
  expect_true(any(grepl("\\\\", result)))  # latex uses backslashes
})

test_that("cf_table returns html format", {
  skip_if_not_installed("knitr")

  fit <- create_test_fit()

  result <- cf_table(fit, format = "html")

  expect_true(is.character(result))
  expect_true(any(grepl("<table", result)))
})

# --- cf_export tests ---

test_that("cf_export creates CSV file", {
  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(c(temp_file, sub("\\.csv$", "_summary.csv", temp_file))), add = TRUE)

  expect_message(
    cf_export(fit, temp_file),
    "exported"
  )

  expect_true(file.exists(temp_file))

  # Check CSV content
  content <- utils::read.csv(temp_file)
  expect_true("ite" %in% names(content))
  expect_equal(nrow(content), fit$meta$n)
})

test_that("cf_export creates summary CSV alongside main file", {
  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".csv")
  summary_file <- sub("\\.csv$", "_summary.csv", temp_file)
  on.exit(unlink(c(temp_file, summary_file)), add = TRUE)

  cf_export(fit, temp_file)

  expect_true(file.exists(summary_file))

  summary_content <- utils::read.csv(summary_file)
  expect_true("Metric" %in% names(summary_content))
})

test_that("cf_export creates RDS file", {
  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  expect_message(
    cf_export(fit, temp_file),
    "exported"
  )

  expect_true(file.exists(temp_file))

  # Check RDS content
  content <- readRDS(temp_file)
  expect_true(is.list(content))
  expect_true("ate" %in% names(content))
  expect_true("method" %in% names(content))
})

test_that("cf_export creates JSON file", {
  skip_if_not_installed("jsonlite")

  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".json")
  on.exit(unlink(temp_file), add = TRUE)

  expect_message(
    cf_export(fit, temp_file),
    "exported"
  )

  expect_true(file.exists(temp_file))

  # Check JSON content
  content <- jsonlite::read_json(temp_file)
  expect_true("ate" %in% names(content))
  expect_true("method" %in% names(content))
})

test_that("cf_export auto-detects format from extension", {
  fit <- create_test_fit()

  # Test RDS detection
  temp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_rds), add = TRUE)

  cf_export(fit, temp_rds)
  content <- readRDS(temp_rds)
  expect_true(is.list(content))
})

test_that("cf_export rejects non-cf_model input", {
  temp_file <- tempfile(fileext = ".csv")

  expect_error(
    cf_export("not a model", temp_file),
    "cf_model"
  )
})

test_that("cf_export includes ITE when include_ite=TRUE", {
  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  cf_export(fit, temp_file, include_ite = TRUE)

  content <- readRDS(temp_file)
  expect_true("ite" %in% names(content))
  expect_equal(length(content$ite), fit$meta$n)
})

test_that("cf_export returns file path invisibly", {
  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  result <- cf_export(fit, temp_file)

  expect_equal(result, temp_file)
})

# --- cf_session_info tests ---

test_that("cf_session_info returns structured list", {
  info <- cf_session_info()

  expect_s3_class(info, "cf_session_info")
  expect_s3_class(info, "list")

  expect_true("timestamp" %in% names(info))
  expect_true("cfomics_version" %in% names(info))
  expect_true("r_version" %in% names(info))
  expect_true("platform" %in% names(info))
})

test_that("cf_session_info includes package versions", {
  info <- cf_session_info()

  expect_true("packages" %in% names(info))
  expect_true("reticulate" %in% names(info$packages))
})

test_that("cf_session_info includes Python info if available", {
  info <- cf_session_info()

  expect_true("python" %in% names(info))
})

test_that("cf_session_info includes fit info when provided", {
  fit <- create_test_fit()

  info <- cf_session_info(fit)

  expect_true("fit" %in% names(info))
  expect_equal(info$fit$method, "gformula")
  expect_equal(info$fit$n, fit$meta$n)
  expect_equal(info$fit$p, fit$meta$p)
})

test_that("print.cf_session_info produces output", {
  info <- cf_session_info()

  expect_output(print(info), "cfomics Session Information")
  expect_output(print(info), "R version")
  expect_output(print(info), "Platform")
})

test_that("print.cf_session_info shows fit info when available", {
  fit <- create_test_fit()
  info <- cf_session_info(fit)

  expect_output(print(info), "Fit Information")
  expect_output(print(info), "Method")
})

# --- Edge cases ---

test_that("cf_table works with minimal fit (no CI)", {
  # Create a fit object manually without CI
  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )
  data$Y <- data$T + rnorm(n)

  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula", nsimul = 10L)

  result <- cf_table(fit)

  expect_s3_class(result, "data.frame")
  expect_true("ATE" %in% result$Metric)
})

test_that("cf_export handles models without y0_hat/y1_hat", {
  skip_if_not_installed("ipw")
  skip_if_not_installed("survey")

  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )
  data$Y <- data$T + rnorm(n)

  fit <- cf_fit(Y ~ T | X1, data = data, method = "ipw")

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  # Should not error
  expect_message(
    cf_export(fit, temp_file),
    "exported"
  )
})

test_that("cf_table works with GRF method", {
  skip_if_not_installed("grf")

  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )
  data$Y <- 1.5 * data$T + rnorm(n)

  fit <- cf_fit(Y ~ T | X1, data = data, method = "grf")

  result <- cf_table(fit)

  expect_s3_class(result, "data.frame")
  expect_true(any(result$Value == "grf"))
})

# --- File overwrite protection tests ---

test_that("cf_export errors when file exists and overwrite=FALSE", {
  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  # Create the file first
  cf_export(fit, temp_file)

  # Should error when trying to overwrite without permission
expect_error(
    cf_export(fit, temp_file, overwrite = FALSE),
    "already exist"
  )
})

test_that("cf_export succeeds when file exists and overwrite=TRUE", {
  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  # Create the file first
  cf_export(fit, temp_file)

  # Should succeed with overwrite=TRUE
  expect_message(
    cf_export(fit, temp_file, overwrite = TRUE),
    "exported"
  )
})

test_that("cf_export checks both CSV files for overwrite protection", {
  fit <- create_test_fit()

  temp_file <- tempfile(fileext = ".csv")
  summary_file <- sub("\\.csv$", "_summary.csv", temp_file)
  on.exit(unlink(c(temp_file, summary_file)), add = TRUE)

  # Create the files first
  cf_export(fit, temp_file)

  # Should error when trying to overwrite without permission
  expect_error(
    cf_export(fit, temp_file, overwrite = FALSE),
    "already exist"
  )

  # Should succeed with overwrite=TRUE
  expect_message(
    cf_export(fit, temp_file, overwrite = TRUE),
    "exported"
  )
})
