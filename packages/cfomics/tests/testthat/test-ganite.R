test_that("cf_fit with method='ganite' basic execution", {
  skip_if_not(
    requireNamespace("reticulate", quietly = TRUE),
    message = "reticulate not installed"
  )

  skip_if_not({
    tryCatch({
      reticulate::py_run_string("import tensorflow")
      TRUE
    }, error = function(e) FALSE)
  }, message = "Python module 'tensorflow' not importable")

  mod <- tryCatch(.get_ganite_module(), error = function(e) NULL)
  skip_if(is.null(mod), "GANITE module not found")

  # Synthetic Data
  set.seed(123)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- as.integer(X1 > 0)
  Y <- X1 + X2 + T + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  # Graph (T -> Y, X1 -> T, X1 -> Y, X2 -> Y)
  library(igraph)
  g <- graph_from_edgelist(
    matrix(c("X1", "T",
             "T", "Y",
             "X1", "Y",
             "X2", "Y"), ncol = 2, byrow = TRUE)
  )

  # Fit with GANITE method (use fewer iterations for testing)
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "ganite", graph = g,
                iterations = 100L, batch_size = 32L, verbose = FALSE)

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "ganite")

  # Predict ATE
  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  expect_length(ate, 1)

  # Predict ITE
  ite <- predict(fit, type = "ite")
  expect_type(ite, "double")
  expect_equal(length(ite), n)

  # Predict Y0 and Y1
  y0 <- predict(fit, type = "y0")
  y1 <- predict(fit, type = "y1")
  expect_equal(length(y0), n)
  expect_equal(length(y1), n)

  # Verify ITE = Y1 - Y0
  expect_equal(ite, y1 - y0, tolerance = 1e-6)
})

test_that("cf_fit_ganite returns correct structure", {
  skip_if_not(
    requireNamespace("reticulate", quietly = TRUE),
    message = "reticulate not installed"
  )

  skip_if_not({
    tryCatch({
      reticulate::py_run_string("import tensorflow")
      TRUE
    }, error = function(e) FALSE)
  }, message = "Python module 'tensorflow' not importable")

  mod <- tryCatch(.get_ganite_module(), error = function(e) NULL)
  skip_if(is.null(mod), "GANITE module not found")

  # Synthetic Data
  set.seed(456)
  n <- 50
  X <- matrix(rnorm(n * 3), ncol = 3)
  T <- rbinom(n, 1, 0.5)
  Y <- X[, 1] + 0.5 * X[, 2] + T + rnorm(n)

  # Fit directly using internal function
  fit <- cf_fit_ganite(
    X = X, T = T, Y = Y,
    dag_spec = NULL,
    covariate_names = c("X1", "X2", "X3"),
    random_state = 42L,
    iterations = 50L,
    batch_size = 16L,
    verbose = FALSE
  )

  expect_s3_class(fit, "cf_model_ganite")
  expect_equal(fit$backend, "ganite")
  expect_true(!is.null(fit$res$ate))
  expect_true(!is.null(fit$res$ite))
  expect_true(!is.null(fit$res$y0_hat))
  expect_true(!is.null(fit$res$y1_hat))
  expect_true(!is.null(fit$res$summary))
  expect_true(!is.null(fit$res$samples))
})

test_that("predict_cf_ganite returns correct types", {
  skip_if_not(
    requireNamespace("reticulate", quietly = TRUE),
    message = "reticulate not installed"
  )

  skip_if_not({
    tryCatch({
      reticulate::py_run_string("import tensorflow")
      TRUE
    }, error = function(e) FALSE)
  }, message = "Python module 'tensorflow' not importable")

  mod <- tryCatch(.get_ganite_module(), error = function(e) NULL)
  skip_if(is.null(mod), "GANITE module not found")

  # Synthetic Data
  set.seed(789)
  n <- 30
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- as.integer(X1 > 0)
  Y <- X1 + X2 + T + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  library(igraph)
  g <- graph_from_edgelist(
    matrix(c("X1", "T", "T", "Y", "X2", "Y"), ncol = 2, byrow = TRUE)
  )

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "ganite", graph = g,
                iterations = 50L, batch_size = 16L, verbose = FALSE)

  # Test summary type
  summary_result <- predict(fit, type = "summary")
  expect_type(summary_result, "list")
  expect_true("ate" %in% names(summary_result))
  expect_true("ite_mean" %in% names(summary_result))
  expect_true("ite_std" %in% names(summary_result))
  expect_true("ite_quantiles" %in% names(summary_result))

  # Test samples type
  samples_result <- predict(fit, type = "samples")
  expect_type(samples_result, "list")
  expect_true("y0" %in% names(samples_result))
  expect_true("y1" %in% names(samples_result))
  expect_true("ite" %in% names(samples_result))
})
