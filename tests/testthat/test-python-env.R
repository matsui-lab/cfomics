# Tests for Python environment management functions
# (cf_install_python_env, cf_use_python_env, cf_list_python_envs, ensure_python_env, ensure_conda_env)

# Note: These tests focus on input validation and structure.
# Actual conda operations are mocked or skipped to avoid system-level changes.

test_that("cf_list_python_envs returns valid structure", {
  result <- cf_list_python_envs()

  expect_true(is.data.frame(result))
  expect_true("method" %in% names(result))
  expect_true("env_name" %in% names(result))
  expect_true("installed" %in% names(result))

  # Should have all known methods
  expected_methods <- c("unified", "unified_full", "dowhy", "ganite", "econml", "pyro")
  expect_true(all(expected_methods %in% result$method))

  # installed should be logical

  expect_type(result$installed, "logical")
})

test_that("cf_list_python_envs env_name matches expected pattern", {
  result <- cf_list_python_envs()

  # All env names should start with "cfomics_"
  expect_true(all(grepl("^cfomics_", result$env_name)))

  # Specific env names
  expect_equal(result$env_name[result$method == "unified"], "cfomics_unified")
  expect_equal(result$env_name[result$method == "dowhy"], "cfomics_dowhy")
  expect_equal(result$env_name[result$method == "ganite"], "cfomics_ganite")
})

test_that("cf_install_python_env validates method argument", {
  # Invalid method should error
  expect_error(
    cf_install_python_env(method = "invalid_method"),
    "should be one of"
  )
})

test_that("cf_use_python_env validates method argument", {
  # Invalid method should error
  expect_error(
    cf_use_python_env(method = "invalid_method"),
    "should be one of"
  )
})

test_that("setup_python_env validates method argument", {
  # Invalid method should error
  expect_error(
    setup_python_env(method = "invalid_method"),
    "should be one of"
  )
})

test_that("ensure_python_env works when Python is available", {
  skip_if_not(
    reticulate::py_available(initialize = FALSE) ||
      nzchar(Sys.getenv("RETICULATE_PYTHON")),
    "Python not available"
  )

  # Skip if Python is not built with shared library (common in containerized envs)
  result <- tryCatch({
    ensure_python_env()
  }, error = function(e) {
    if (grepl("shared library", e$message)) {
      skip("Python not built with shared library")
    }
    stop(e)
  })
  expect_true(result)
})

test_that("ensure_python_env errors gracefully without Python", {
  skip_if(
    reticulate::py_available(initialize = FALSE) ||
      nzchar(Sys.getenv("RETICULATE_PYTHON")),
    "Python is available, skipping 'no Python' test"
  )

  expect_error(
    ensure_python_env(),
    "No Python environment configured"
  )
})

test_that("ensure_conda_env handles invalid method gracefully", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))

  # When Python is already available, it should return immediately
  if (reticulate::py_available(initialize = FALSE)) {
    result <- ensure_conda_env(method = "nonexistent_method")
    expect_true(result)
  }
})

test_that("ensure_conda_env defaults to dowhy for unknown method", {
  # This tests the fallback behavior
  skip_if(reticulate::py_available(initialize = FALSE))

  # Should attempt to use dowhy env when given unknown method
  # This will error if env doesn't exist, which is expected
  expect_error(
    ensure_conda_env(method = "unknown"),
    "cfomics_dowhy"
  )
})

test_that(".cfomics_env_specs has correct structure", {
  # Access internal specs
  specs <- cfomics:::.cfomics_env_specs

  expect_true(is.list(specs))

  # Check each spec has required fields
  for (name in names(specs)) {
    spec <- specs[[name]]
    expect_true("name" %in% names(spec), info = paste("Missing 'name' in", name))
    expect_true("packages" %in% names(spec), info = paste("Missing 'packages' in", name))
    expect_true("pip_packages" %in% names(spec), info = paste("Missing 'pip_packages' in", name))
  }
})

test_that(".cfomics_env_specs unified contains expected packages", {
  specs <- cfomics:::.cfomics_env_specs

  # unified should have numpy, pandas
  expect_true("numpy" %in% specs$unified$packages)
  expect_true("pandas" %in% specs$unified$packages)

  # unified should have dowhy and econml via pip
  expect_true("dowhy==0.11.1" %in% specs$unified$pip_packages)
  expect_true("econml" %in% specs$unified$pip_packages)
})

test_that(".cfomics_env_specs unified_full contains deep learning packages", {
  specs <- cfomics:::.cfomics_env_specs

  # unified_full should have tensorflow, torch, pyro
  pip_pkgs <- specs$unified_full$pip_packages
  expect_true("tensorflow" %in% pip_pkgs)
  expect_true("torch" %in% pip_pkgs)
  expect_true("pyro-ppl" %in% pip_pkgs)
})

test_that(".cfomics_env_specs method-specific envs are correctly defined", {
  specs <- cfomics:::.cfomics_env_specs

  # dowhy
  expect_true("dowhy==0.11.1" %in% specs$dowhy$pip_packages)

  # ganite
  expect_true("tensorflow" %in% specs$ganite$pip_packages)

  # econml
  expect_true("econml" %in% specs$econml$pip_packages)

  # pyro
  expect_true("torch" %in% specs$pyro$pip_packages)
  expect_true("pyro-ppl" %in% specs$pyro$pip_packages)
})

# Tests for setup_python_env logic (without actual conda operations)
test_that("setup_python_env uses correct env names", {
  specs <- cfomics:::.cfomics_env_specs

  expect_equal(specs$unified$name, "cfomics_unified")
  expect_equal(specs$unified_full$name, "cfomics_unified_full")
  expect_equal(specs$dowhy$name, "cfomics_dowhy")
  expect_equal(specs$ganite$name, "cfomics_ganite")
  expect_equal(specs$econml$name, "cfomics_econml")
  expect_equal(specs$pyro$name, "cfomics_pyro")
})
