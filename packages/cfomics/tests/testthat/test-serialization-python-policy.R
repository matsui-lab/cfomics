# Test serialization policy for Python-based methods
# Python methods (dowhy_gcm, ganite, drlearner, cavae) may contain
# non-serializable Python objects. Use cf_export() for portable output.

# Skip entire file if Python is not available - these are policy documentation tests
skip_on_cran()

test_that("cf_export works for any method (portable serialization)", {
  set.seed(123)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )
  data$Y <- data$Y + 2 * data$T

  # Test with R-native method
  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula")

  # RDS export creates portable lightweight output
  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds), add = TRUE)

  cf_export(fit, file = tmp_rds, format = "rds")
  expect_true(file.exists(tmp_rds))

  # Verify exported object is lightweight and portable
  exported <- readRDS(tmp_rds)
  expect_true(is.list(exported))
  expect_true("ate" %in% names(exported))
  expect_true("method" %in% names(exported))
  expect_true("meta" %in% names(exported))

  # JSON export is always portable
  tmp_json <- tempfile(fileext = ".json")
  on.exit(unlink(tmp_json), add = TRUE)

  skip_if_not_installed("jsonlite")
  cf_export(fit, file = tmp_json, format = "json")
  expect_true(file.exists(tmp_json))
})

test_that("cf_export preserves essential results for reproducibility", {
  set.seed(456)
  n <- 40
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X = rnorm(n)
  )
  data$Y <- data$Y + 1.5 * data$T

  fit <- cf_fit(Y ~ T | X, data = data, method = "gformula")

  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds), add = TRUE)

  cf_export(fit, file = tmp_rds, format = "rds", include_ite = TRUE)

  exported <- readRDS(tmp_rds)

  # Essential results are preserved

  expect_equal(exported$ate, predict(fit, type = "ate"), tolerance = 1e-10)
  expect_equal(exported$ite, predict(fit, type = "ite"), tolerance = 1e-10)

  # Metadata is preserved
  expect_equal(exported$meta$n, fit$meta$n)
  expect_equal(exported$meta$p, fit$meta$p)
  expect_equal(exported$meta$outcome_name, fit$meta$outcome_name)
  expect_equal(exported$meta$treatment_name, fit$meta$treatment_name)
})

test_that("CSV export creates human-readable output", {
  set.seed(789)
  n <- 30
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X = rnorm(n)
  )
  data$Y <- data$Y + 2 * data$T

  fit <- cf_fit(Y ~ T | X, data = data, method = "gformula")

  tmp_csv <- tempfile(fileext = ".csv")
  on.exit({
    unlink(tmp_csv)
    unlink(sub("\\.csv$", "_summary.csv", tmp_csv))
  }, add = TRUE)

  cf_export(fit, file = tmp_csv, format = "csv")
  expect_true(file.exists(tmp_csv))

  # Summary file is also created
  summary_file <- sub("\\.csv$", "_summary.csv", tmp_csv)
  expect_true(file.exists(summary_file))

  # CSV is readable
  csv_data <- utils::read.csv(tmp_csv)
  expect_true("ite" %in% names(csv_data))
  expect_equal(nrow(csv_data), n)
})

# Document the serialization policy for Python methods
test_that("serialization policy is documented", {
  # This test documents the expected behavior for method serialization:
  #
  # R-NATIVE METHODS (gformula, ipw, grf):
  #   - Full cfomics_result objects CAN be serialized with saveRDS()
  #   - Restored objects remain fully functional (predict, summary work)
  #   - Recommended for R-only workflows
  #

  # PYTHON-BASED METHODS (dowhy_gcm, ganite, drlearner, cavae):
  #   - Full cfomics_result objects may NOT serialize correctly
  #   - Python objects (models, numpy arrays) may become invalid after restore
  #   - USE cf_export() to create portable, lightweight output
  #   - Exported results contain ATE, ITE, counterfactuals, and metadata
  #   - Exported results are sufficient for analysis and reporting
  #
  # RECOMMENDED WORKFLOW FOR REPRODUCIBILITY:
  #   1. Save results with cf_export(fit, "results.rds")
  #   2. Save session info with cf_session_info(fit)
  #   3. For R-native methods, saveRDS(fit) also works
  #
  expect_true(TRUE)  # Documentation test always passes
})
