test_that("cf_fit_drlearner basic execution", {
  skip_if_not(
    requireNamespace("reticulate", quietly = TRUE),
    message = "reticulate not installed"
  )
  
  if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required = FALSE)
  }
  
  skip_if_not(
    reticulate::py_available(initialize = TRUE),
    message = "Python not available"
  )
  
  skip_if_not(
    reticulate::py_module_available("sklearn"),
    message = "Python module 'sklearn' not available"
  )
  
  econml_available <- tryCatch({
    reticulate::import("econml")
    TRUE
  }, error = function(e) FALSE)
  
  skip_if_not(econml_available, "econml not available")
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)
  
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")
  
  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "drlearner")
  
  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  
  ite <- predict(fit, type = "ite")
  expect_equal(length(ite), n)
  
  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("ite_mean" %in% names(summary_stats))
})

test_that("drlearner prediction on new data", {
  skip_if_not(
    requireNamespace("reticulate", quietly = TRUE),
    message = "reticulate not installed"
  )
  
  if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required = FALSE)
  }
  
  skip_if_not(
    reticulate::py_available(initialize = TRUE),
    message = "Python not available"
  )
  
  skip_if_not(
    reticulate::py_module_available("sklearn"),
    message = "Python module 'sklearn' not available"
  )
  
  econml_available <- tryCatch({
    reticulate::import("econml")
    TRUE
  }, error = function(e) FALSE)
  
  skip_if_not(econml_available, "econml not available")
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1.5 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)
  
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "drlearner")
  
  newdata <- data.frame(X1 = rnorm(10), X2 = rnorm(10))
  
  ite_new <- predict_cf_drlearner(fit, newdata = newdata, type = "ite")
  expect_equal(length(ite_new), 10)
})
