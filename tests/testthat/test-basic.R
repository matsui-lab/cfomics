test_that("cf_fit_dowhy_gcm basic execution", {
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
    reticulate::py_module_available("dowhy"),
    message = "Python module 'dowhy' not available"
  )

  mod <- tryCatch(.get_dowhy_module(), error = function(e) NULL)
  skip_if(is.null(mod), "DoWhy module not found")

  # Synthetic Data
  set.seed(123)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- as.integer(X1 > 0)
  Y <- X1 + X2 + T + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  # Graph
  # T -> Y, X1 -> T, X1 -> Y, X2 -> Y
  library(igraph)
  g <- graph_from_edgelist(
      matrix(c("X1", "T",
               "T", "Y",
               "X1", "Y",
               "X2", "Y"), ncol=2, byrow=TRUE)
  )

  # Fit
  fit <- cf_dowhy(Y ~ T | X1 + X2, data = data, graph = g)

  expect_s3_class(fit, "cf_model")
  expect_s3_class(fit, "cfomics_result")

  # Predict
  ate <- predict(fit, type = "ate")
  expect_type(ate, "double")

  ite <- predict(fit, type = "ite")
  expect_equal(length(ite), n)
})

test_that("parse_cf_formula validation", {
  data <- data.frame(Y = rnorm(10), T = rbinom(10, 1, 0.5), X = rnorm(10), C = letters[1:10])

  # Factor/Character X should fail
  expect_error(cf_fit(Y ~ T | C, data = data, graph = igraph::make_empty_graph()),
               "numeric covariates")

  # Bad variable names
  data$bad.name <- rnorm(10)
  expect_error(cf_fit(Y ~ T | bad.name, data = data, graph = igraph::make_empty_graph()),
               "variable names")
})
