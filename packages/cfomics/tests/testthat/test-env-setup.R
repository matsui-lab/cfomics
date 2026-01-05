# Tests for environment setup and diagnostics (cf_check_env, cf_config)

test_that("cf_check_env returns valid structure", {
  result <- cf_check_env(verbose = FALSE)

  expect_true(is.list(result))
  expect_s3_class(result, "cf_env_check")

  # Check required top-level components
  expect_true("r_version" %in% names(result))
  expect_true("os" %in% names(result))
  expect_true("python" %in% names(result))
  expect_true("conda" %in% names(result))
  expect_true("packages" %in% names(result))
  expect_true("methods" %in% names(result))
  expect_true("resources" %in% names(result))
  expect_true("timestamp" %in% names(result))
})

test_that("cf_check_env r_version check works", {
  result <- cf_check_env(verbose = FALSE)

  expect_true(is.list(result$r_version))
  expect_true("version" %in% names(result$r_version))
  expect_true("ok" %in% names(result$r_version))
  expect_true("message" %in% names(result$r_version))

  # Version should be a character string
  expect_type(result$r_version$version, "character")
  expect_type(result$r_version$ok, "logical")
})

test_that("cf_check_env methods check works", {
  result <- cf_check_env(verbose = FALSE, methods = c("grf", "gformula"))

  expect_true(is.list(result$methods))
  expect_true("grf" %in% names(result$methods))
  expect_true("gformula" %in% names(result$methods))

  # Each method should have available, type, message
  for (m in names(result$methods)) {
    expect_true("available" %in% names(result$methods[[m]]))
    expect_true("type" %in% names(result$methods[[m]]))
    expect_true("message" %in% names(result$methods[[m]]))
  }
})

test_that("cf_check_env gformula is always available", {
  result <- cf_check_env(verbose = FALSE, methods = "gformula")

  # gformula is built-in, should always be available

  expect_true(result$methods$gformula$available)
  expect_equal(result$methods$gformula$type, "R")
})

test_that("cf_check_env resources check works", {
  result <- cf_check_env(verbose = FALSE)

  expect_true(is.list(result$resources))
  expect_true("memory" %in% names(result$resources))
  expect_true("cores" %in% names(result$resources))

  # Cores should be numeric or NA
  expect_true(is.numeric(result$resources$cores) || is.na(result$resources$cores))
})

test_that("cf_check_env packages check structure", {
  result <- cf_check_env(verbose = FALSE)

  expect_true(is.list(result$packages))
  expect_true("r" %in% names(result$packages))
  expect_true("python" %in% names(result$packages))

  # R packages should have required and suggested
  expect_true("required" %in% names(result$packages$r))
  expect_true("suggested" %in% names(result$packages$r))
})

test_that("cf_config returns defaults when no config file exists", {
  # Get default method
  default_method <- cf_config("methods.default")

  # Should be character or NULL
  expect_true(is.character(default_method) || is.null(default_method))

  # If not null, should be a valid method name
  if (!is.null(default_method)) {
    valid_methods <- c("grf", "ipw", "gformula", "dowhy_gcm", "drlearner", "ganite", "cavae")
    expect_true(default_method %in% valid_methods)
  }
})

test_that("cf_config returns NULL for non-existent keys", {
  result <- cf_config("nonexistent.deeply.nested.key")
  expect_null(result)
})

test_that("cf_config can get nested values", {
  # Test accessing a nested value
  python_path <- cf_config("python.path")

  # Should return the default "auto" or NULL if not configured
  expect_true(is.null(python_path) || is.character(python_path))
})

test_that("cf_config shows all config with no arguments", {
  # cf_config() prints output, so capture it
  output <- capture.output({
    result <- cf_config()
  })

  # Result should be a list with config structure
  expect_true(is.list(result))

  # Output should contain config information
  expect_true(length(output) > 0)
})

test_that("print.cf_env_check produces output", {
  result <- cf_check_env(verbose = FALSE)

  # Capture output
  output <- capture.output(print(result))

  # Should produce non-empty output

  expect_true(length(output) > 0)

  # Should contain key sections
  output_text <- paste(output, collapse = "\n")
  expect_true(grepl("cfomics", output_text, ignore.case = TRUE))
  expect_true(grepl("R version", output_text))
})
