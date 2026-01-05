# Tests for cf_metadata() function

test_that("cf_metadata requires cf_model object", {
  # Non cf_model object should error
  expect_error(
    cf_metadata(list(a = 1, b = 2)),
    "must be a cf_model"
  )

  expect_error(
    cf_metadata("not a model"),
    "must be a cf_model"
  )

  expect_error(
    cf_metadata(NULL),
    "must be a cf_model"
  )
})

test_that("cf_metadata returns NULL when no metadata available", {
  skip_if_not_installed("grf")

  # Create test data
  set.seed(42)
  n <- 100
  test_data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    T = rbinom(n, 1, 0.5),
    Y = rnorm(n)
  )

  # Fit model (R methods typically don't set metadata)
  fit <- cf_fit(Y ~ T | X1 + X2, data = test_data, method = "gformula")

  # Capture message about no metadata
  expect_message(
    result <- cf_metadata(fit),
    "No metadata available"
  )

  expect_null(result)
})

test_that("cf_metadata works with grf model", {
  skip_if_not_installed("grf")

  set.seed(42)
  n <- 100
  test_data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    T = rbinom(n, 1, 0.5),
    Y = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = test_data, method = "grf")

  # GRF doesn't set metadata, so this should return NULL with message
  expect_message(
    result <- cf_metadata(fit),
    "No metadata available"
  )

  expect_null(result)
})

test_that("cf_metadata works with ipw model", {
  skip_if_not_installed("ipw")
  skip_if_not_installed("survey")

  set.seed(42)
  n <- 100
  test_data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    T = rbinom(n, 1, 0.5),
    Y = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = test_data, method = "ipw")

  expect_message(
    result <- cf_metadata(fit),
    "No metadata available"
  )

  expect_null(result)
})

test_that("cf_metadata works with gformula model", {
  set.seed(42)
  n <- 100
  test_data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    T = rbinom(n, 1, 0.5),
    Y = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = test_data, method = "gformula")

  expect_message(
    result <- cf_metadata(fit),
    "No metadata available"
  )

  expect_null(result)
})

test_that("cf_metadata handles cfomics_result class", {
  set.seed(42)
  n <- 100
  test_data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    T = rbinom(n, 1, 0.5),
    Y = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = test_data, method = "gformula")

  # Check dual class
  expect_true(inherits(fit, "cf_model"))
  expect_true(inherits(fit, "cfomics_result"))

  # cf_metadata should work with either class
  expect_message(
    cf_metadata(fit),
    "No metadata available"
  )
})

# Test structure expectations when metadata IS present
test_that("cf_metadata returns metadata when set", {
  set.seed(42)
  n <- 100
  test_data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    T = rbinom(n, 1, 0.5),
    Y = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = test_data, method = "gformula")

  # Manually add metadata to test extraction
  fit$fit$res$metadata <- list(
    input_hash = "abc123",
    dag_structure = "Y <- T <- X1",
    random_seed = 42,
    processing_time = 1.5
  )

  result <- cf_metadata(fit)

  expect_true(is.list(result))
  expect_equal(result$input_hash, "abc123")
  expect_equal(result$random_seed, 42)
  expect_equal(result$processing_time, 1.5)
})

test_that("cf_metadata preserves complex metadata structures", {
  set.seed(42)
  n <- 100
  test_data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    T = rbinom(n, 1, 0.5),
    Y = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = test_data, method = "gformula")

  # Add nested metadata
  fit$fit$res$metadata <- list(
    input_hash = "test123",
    model_info = list(
      type = "gformula",
      n_samples = 100,
      covariates = c("X1", "X2")
    ),
    timing = list(
      start = Sys.time() - 10,
      end = Sys.time()
    )
  )

  result <- cf_metadata(fit)

  expect_true(is.list(result))
  expect_true(is.list(result$model_info))
  expect_equal(result$model_info$type, "gformula")
  expect_equal(result$model_info$n_samples, 100)
  expect_equal(result$model_info$covariates, c("X1", "X2"))
})
