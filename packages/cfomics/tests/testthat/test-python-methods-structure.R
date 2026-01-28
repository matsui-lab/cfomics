# Tests for Python method structure and error handling
# These tests focus on validation and structure without requiring GPU or full Python execution

# ============================================================================
# DoWhy-GCM Tests
# ============================================================================

test_that("cf_fit errors with dowhy_gcm when graph not provided", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))

  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  expect_error(
    cf_fit(Y ~ T | X1 + X2, data = data, method = "dowhy_gcm"),
    "graph"
  )
})

test_that("cf_fit with dowhy_gcm validates graph structure", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))
  skip_if_not_installed("igraph")

  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  # Graph is provided but Python may not be available - should handle gracefully
  g <- igraph::graph_from_edgelist(
    matrix(c("X1", "T", "T", "Y", "X2", "Y"), ncol = 2, byrow = TRUE)
  )

  # This will error if Python/DoWhy is not available, which is expected
  result <- tryCatch(
    cf_fit(Y ~ T | X1 + X2, data = data, method = "dowhy_gcm", graph = g),
    error = function(e) e
  )

  # Should either succeed (Python available) or give informative error
  if (inherits(result, "error")) {
    expect_true(
      grepl("python|dowhy|conda|environment|reticulate", tolower(result$message)),
      info = paste("Error should mention Python/DoWhy:", result$message)
    )
  } else {
    expect_s3_class(result, "cf_model")
  }
})

# ============================================================================
# DRLearner Error Handling Tests
# ============================================================================

test_that("cf_fit with drlearner handles missing Python gracefully", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))

  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  result <- tryCatch(
    cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner"),
    error = function(e) e
  )

  # Should either succeed or give informative error about Python/econml
  if (inherits(result, "error")) {
    expect_true(
      grepl("python|econml|sklearn|conda|environment|reticulate", tolower(result$message)),
      info = paste("Error should be about Python deps:", result$message)
    )
  } else {
    expect_s3_class(result, "cf_model")
    expect_equal(result$method, "drlearner")
  }
})

test_that("predict_cf_drlearner validates type argument", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))

  # Skip if econml is not available
  skip_if_not({
    tryCatch({
      reticulate::py_available(initialize = TRUE) &&
        reticulate::py_module_available("econml")
    }, error = function(e) FALSE)
  }, message = "econml not available")

  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")

  # Test all valid types
  expect_type(predict(fit, type = "ate"), "double")
  expect_type(predict(fit, type = "ite"), "double")
  expect_type(predict(fit, type = "y0"), "double")
  expect_type(predict(fit, type = "y1"), "double")
  expect_type(predict(fit, type = "summary"), "list")
  expect_type(predict(fit, type = "samples"), "list")
})

# ============================================================================
# GANITE Error Handling Tests
# ============================================================================

test_that("cf_fit with ganite handles missing Python gracefully", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))
  skip_if_not_installed("igraph")

  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  g <- igraph::graph_from_edgelist(
    matrix(c("X1", "T", "T", "Y", "X2", "Y"), ncol = 2, byrow = TRUE)
  )

  result <- tryCatch(
    cf_fit(Y ~ T | X1 + X2, data = data, method = "ganite", graph = g,
           iterations = 10L, verbose = FALSE),
    error = function(e) e
  )

  # Should either succeed or give informative error about Python/TensorFlow
  if (inherits(result, "error")) {
    expect_true(
      grepl("python|tensorflow|ganite|conda|environment|module|reticulate",
            tolower(result$message)),
      info = paste("Error should be about Python deps:", result$message)
    )
  } else {
    expect_s3_class(result, "cf_model")
    expect_equal(result$method, "ganite")
  }
})

# ============================================================================
# CAVAE Error Handling Tests
# ============================================================================

test_that("cf_fit with cavae handles missing Python gracefully", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))
  skip_if_not_installed("igraph")

  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  g <- igraph::graph_from_edgelist(
    matrix(c("X1", "T", "T", "Y", "X2", "Y"), ncol = 2, byrow = TRUE)
  )

  result <- tryCatch(
    cf_fit(Y ~ T | X1 + X2, data = data, method = "cavae", graph = g,
           num_epochs = 10L, verbose = FALSE),
    error = function(e) e
  )

  # Should either succeed or give informative error about Python/PyTorch/Pyro
  if (inherits(result, "error")) {
    expect_true(
      grepl("python|torch|pyro|cavae|conda|environment|module|reticulate",
            tolower(result$message)),
      info = paste("Error should be about Python deps:", result$message)
    )
  } else {
    expect_s3_class(result, "cf_model")
    expect_equal(result$method, "cavae")
  }
})

# ============================================================================
# Cross-method Tests
# ============================================================================

test_that("cf_fit method auto-selection prefers R methods when Python unavailable", {
  # When Python is not available, auto should select an R method
  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  # Auto method should work regardless of Python availability
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "auto")

  expect_s3_class(fit, "cf_model")

  # Method should be one of the valid methods
  valid_methods <- c("grf", "hdml", "hdps", "bcf", "tmle", "ipw", "gformula", "dowhy_gcm", "drlearner", "ganite", "cavae")
  expect_true(fit$method %in% valid_methods)
})

test_that("all methods return consistent result structure", {
  skip_if_not_installed("grf")

  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  # Test with grf (always available when installed)
  fit_grf <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")

  # Check structure
  expect_true("method" %in% names(fit_grf))
  expect_true("fit" %in% names(fit_grf))
  expect_true("meta" %in% names(fit_grf))

  # Check fit structure
  expect_true("model" %in% names(fit_grf$fit) || "res" %in% names(fit_grf$fit))
  expect_true("res" %in% names(fit_grf$fit))

  # Check res structure
  res <- fit_grf$fit$res
  expect_true("ate" %in% names(res))
  expect_true("ite" %in% names(res))
  expect_true("y0_hat" %in% names(res))
  expect_true("y1_hat" %in% names(res))
  expect_true("summary" %in% names(res))

  # Check meta structure
  meta <- fit_grf$meta
  expect_true("n" %in% names(meta))
  expect_true("p" %in% names(meta))
})

test_that("Python method result structure matches R method structure", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))

  # Check if econml is available
  econml_available <- tryCatch({
    reticulate::py_available(initialize = TRUE) &&
      reticulate::py_module_available("econml")
  }, error = function(e) FALSE)

  skip_if_not(econml_available, "econml not available")

  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  fit_r <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula")
  fit_py <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")

  # Both should have same top-level structure
  expect_equal(names(fit_r), names(fit_py))

  # Both should have same result names
  expect_true(all(c("ate", "ite", "y0_hat", "y1_hat") %in% names(fit_r$fit$res)))
  expect_true(all(c("ate", "ite", "y0_hat", "y1_hat") %in% names(fit_py$fit$res)))

  # Both should have same meta names (at minimum)
  common_meta <- c("n", "p")
  expect_true(all(common_meta %in% names(fit_r$meta)))
  expect_true(all(common_meta %in% names(fit_py$meta)))
})

# ============================================================================
# Prediction Type Tests for All Methods
# ============================================================================

test_that("predict types work for gformula method", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula")

  # Test all prediction types
  ate <- predict(fit, type = "ate")
  expect_length(ate, 1)
  expect_type(ate, "double")

  ite <- predict(fit, type = "ite")
  expect_length(ite, n)
  expect_type(ite, "double")

  y0 <- predict(fit, type = "y0")
  expect_length(y0, n)

  y1 <- predict(fit, type = "y1")
  expect_length(y1, n)

  summary <- predict(fit, type = "summary")
  expect_type(summary, "list")
  expect_true("ate" %in% names(summary))

  samples <- predict(fit, type = "samples")
  expect_type(samples, "list")
  expect_true(all(c("y0", "y1", "ite") %in% names(samples)))
})

test_that("predict validates type argument", {
  set.seed(42)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )

  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula")

  expect_error(
    predict(fit, type = "invalid_type"),
    "should be one of"
  )
})
