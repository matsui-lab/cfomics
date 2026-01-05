test_that("cf_fit_cavae basic execution", {
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
  
  pyro_available <- tryCatch({
    reticulate::import("pyro")
    TRUE
  }, error = function(e) FALSE)
  
  skip_if_not(pyro_available, "Pyro not available")
  
  torch_available <- tryCatch({
    reticulate::import("torch")
    TRUE
  }, error = function(e) FALSE)
  
  skip_if_not(torch_available, "PyTorch not available")
  
  # CAVAE requires GPU for practical use - skip on CPU-only environments
  # This is intentional: CAVAE training is very slow on CPU and may timeout
  torch <- reticulate::import("torch", convert = TRUE)
  cuda_available <- tryCatch({
    torch$cuda$is_available()
  }, error = function(e) FALSE)
  
  skip_if_not(cuda_available, "CUDA/GPU not available (CAVAE requires GPU)")
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- rbinom(n, 1, plogis(1.5 * T + 0.8 * X1))
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)
  
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "cavae",
                num_epochs = 5L, batch_size = 50L)
  
  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")
  expect_equal(fit$method, "cavae")
  
  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")
  
  ite <- predict(fit, type = "ite")
  expect_equal(length(ite), n)
  
  summary_stats <- predict(fit, type = "summary")
  expect_true("ate" %in% names(summary_stats))
  expect_true("ite_mean" %in% names(summary_stats))
})

test_that("cavae prediction on new data", {
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
  
  pyro_available <- tryCatch({
    reticulate::import("pyro")
    TRUE
  }, error = function(e) FALSE)
  
  skip_if_not(pyro_available, "Pyro not available")
  
  torch_available <- tryCatch({
    reticulate::import("torch")
    TRUE
  }, error = function(e) FALSE)
  
  skip_if_not(torch_available, "PyTorch not available")
  
  # CAVAE requires GPU for practical use - skip on CPU-only environments
  torch <- reticulate::import("torch", convert = TRUE)
  cuda_available <- tryCatch({
    torch$cuda$is_available()
  }, error = function(e) FALSE)
  
  skip_if_not(cuda_available, "CUDA/GPU not available (CAVAE requires GPU)")
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- rbinom(n, 1, plogis(1.5 * T + 0.8 * X1))
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)
  
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "cavae",
                num_epochs = 5L, batch_size = 50L)
  
  newdata <- data.frame(X1 = rnorm(10), X2 = rnorm(10))
  
  ite_new <- predict_cf_cavae(fit, newdata = newdata, type = "ite")
  expect_equal(length(ite_new), 10)
})
