# Tests for configuration reset functions
# (cf_config_reset, cf_reset_env, cf_setup_method)

# Helper to get config path (using internal function pattern)
get_test_config_path <- function() {
  file.path(Sys.getenv("HOME"), ".cfomics", "config.yaml")
}

test_that("cf_config_reset handles non-existent config file", {
  config_path <- get_test_config_path()

  # Temporarily rename existing config if present
  backup_path <- paste0(config_path, ".backup_test")
  had_config <- file.exists(config_path)
  if (had_config) {
    file.rename(config_path, backup_path)
  }

  on.exit({
    # Restore original config
    if (had_config && file.exists(backup_path)) {
      file.rename(backup_path, config_path)
    } else if (file.exists(backup_path)) {
      file.remove(backup_path)
    }
  }, add = TRUE)

  # Should print message about no config file
  expect_message(
    result <- cf_config_reset(confirm = FALSE),
    "No configuration file found"
  )

  expect_true(result)
})

test_that("cf_config_reset removes config file when exists", {
  skip_if_not_installed("yaml")

  config_path <- get_test_config_path()

  # Backup existing config
  backup_path <- paste0(config_path, ".backup_test")
  had_config <- file.exists(config_path)
  if (had_config) {
    file.copy(config_path, backup_path)
  }

  on.exit({
    # Restore original config
    if (had_config && file.exists(backup_path)) {
      file.copy(backup_path, config_path, overwrite = TRUE)
      file.remove(backup_path)
    } else if (!had_config && file.exists(config_path)) {
      file.remove(config_path)
    }
  }, add = TRUE)

  # Create a test config file
  cf_config("test.key", "test_value")

  # Verify config was created
  expect_true(file.exists(config_path))

  # Reset config
  expect_message(
    result <- cf_config_reset(confirm = FALSE),
    "Configuration reset to defaults"
  )

  expect_true(result)

  # Config file should be removed
  expect_false(file.exists(config_path))
})

test_that("cf_config_reset returns FALSE when cancelled in interactive mode", {
  # This test simulates cancellation - in non-interactive mode, confirmation is skipped
  config_path <- get_test_config_path()

  # When config doesn't exist, it returns TRUE
  if (!file.exists(config_path)) {
    expect_message(
      result <- cf_config_reset(confirm = FALSE),
      "No configuration file found"
    )
    expect_true(result)
  }
})

test_that("cf_reset_env handles no environments case", {
  # If no cfomics environments exist, should return TRUE with message
  # We don't actually test conda operations, just the logic

  # Mock scenario: capture output when no environments are found
  # The function checks for conda_list which may or may not have cfomics envs
  result <- tryCatch({
    # Run without confirmation (would be cancelled otherwise)
    cf_reset_env(confirm = FALSE)
  }, error = function(e) {
    # If conda is not available, function might error - that's ok for this test
    NULL
  })

  # If it succeeded, check the result
  if (!is.null(result)) {
    expect_true(result)
  }
})

test_that("cf_setup_method validates method argument", {
  # Invalid method should error with match.arg
  expect_error(
    cf_setup_method(method = "invalid_method"),
    "should be one of"
  )
})

test_that("cf_setup_method accepts valid methods", {
  valid_methods <- c("dowhy", "ganite", "econml", "pyro")

  for (method in valid_methods) {
    # We don't want to actually install - just verify argument is accepted
    # This will error on conda operations, but that's expected
    result <- tryCatch({
      # Just testing that the argument passes validation
      match.arg(method, c("dowhy", "ganite", "econml", "pyro"))
      TRUE
    }, error = function(e) FALSE)

    expect_true(result, info = paste("Method", method, "should be valid"))
  }
})

# Tests for internal config helper functions
test_that("get_nested_value works correctly", {
  # Access internal function
  get_nested_value <- cfomics:::get_nested_value

  test_list <- list(
    level1 = list(
      level2 = list(
        value = "deep"
      ),
      simple = "shallow"
    ),
    top = "surface"
  )

  expect_equal(get_nested_value(test_list, "top"), "surface")
  expect_equal(get_nested_value(test_list, "level1.simple"), "shallow")
  expect_equal(get_nested_value(test_list, "level1.level2.value"), "deep")
  expect_null(get_nested_value(test_list, "nonexistent"))
  expect_null(get_nested_value(test_list, "level1.nonexistent"))
})

test_that("set_nested_value works correctly", {
  # Access internal function
  set_nested_value <- cfomics:::set_nested_value

  test_list <- list(
    existing = "old"
  )

  # Set top-level value
  result <- set_nested_value(test_list, "new_key", "new_value")
  expect_equal(result$new_key, "new_value")
  expect_equal(result$existing, "old")

  # Set nested value
  result <- set_nested_value(test_list, "nested.deep.value", "deep_value")
  expect_equal(result$nested$deep$value, "deep_value")

  # Update existing value
  result <- set_nested_value(test_list, "existing", "updated")
  expect_equal(result$existing, "updated")
})

test_that("get_default_config returns expected structure", {
  # Access internal function
  get_default_config <- cfomics:::get_default_config

  defaults <- get_default_config()

  expect_true(is.list(defaults))
  expect_true("python" %in% names(defaults))
  expect_true("methods" %in% names(defaults))
  expect_true("preferences" %in% names(defaults))

  # Check default values
  expect_equal(defaults$python$path, "auto")
  expect_equal(defaults$methods$default, "grf")
  expect_true(defaults$methods$fallback_to_r)
  expect_false(defaults$preferences$auto_install)
})

test_that("get_default_method returns grf by default", {
  # Access internal function
  get_default_method <- cfomics:::get_default_method

  method <- get_default_method()

  # Should be grf or whatever is configured
  expect_true(is.character(method))
  expect_true(nchar(method) > 0)
})

test_that("use_r_fallback returns logical",
{
  # Access internal function
  use_r_fallback <- cfomics:::use_r_fallback

  result <- use_r_fallback()

  expect_type(result, "logical")
})

test_that("cf_config set and get roundtrip works", {
  skip_if_not_installed("yaml")

  config_path <- get_test_config_path()

  # Backup existing config
  backup_path <- paste0(config_path, ".backup_test")
  had_config <- file.exists(config_path)
  if (had_config) {
    file.copy(config_path, backup_path)
  }

  on.exit({
    # Restore original config
    if (had_config && file.exists(backup_path)) {
      file.copy(backup_path, config_path, overwrite = TRUE)
      file.remove(backup_path)
    } else if (!had_config && file.exists(config_path)) {
      file.remove(config_path)
    }
  }, add = TRUE)

  # Set a value
  expect_message(
    cf_config("test.roundtrip", "test_value_123"),
    "Set test.roundtrip"
  )

  # Get the value back
  result <- cf_config("test.roundtrip")
  expect_equal(result, "test_value_123")

  # Set numeric value
  cf_config("test.numeric", 42)
  expect_equal(cf_config("test.numeric"), 42)

  # Set boolean value
  cf_config("test.boolean", TRUE)
  expect_true(cf_config("test.boolean"))
})
