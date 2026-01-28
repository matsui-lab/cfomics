# Tests for input validation and error handling

test_that("cf_fit rejects non-numeric outcome", {
  data <- data.frame(
    Y = letters[1:20],
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "numeric"
  )
})

test_that("cf_fit rejects non-binary treatment", {
  data <- data.frame(
    Y = rnorm(20),
    T = 1:20,
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "binary"
  )
})

test_that("cf_fit rejects treatment with values other than 0/1", {
  data <- data.frame(
    Y = rnorm(20),
    T = sample(c(1, 2), 20, replace = TRUE),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "binary"
  )
})

test_that("cf_fit accepts logical treatment and converts to 0/1", {
  set.seed(42)
  data <- data.frame(
    Y = rnorm(50),
    T = sample(c(TRUE, FALSE), 50, replace = TRUE),
    X1 = rnorm(50)
  )

  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula", nsimul = 10L)
  expect_s3_class(fit, "cf_model")
})

test_that("cf_fit rejects non-numeric covariates", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X1 = letters[1:20]
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    "numeric"
  )
})

test_that("cf_fit rejects variable names with special characters", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    `x.1` = rnorm(20),
    check.names = FALSE
  )
  expect_error(
    cf_fit(Y ~ T | `x.1`, data = data, method = "gformula"),
    "letters|digits|underscore"
  )
})

test_that("cf_fit handles NA values in outcome", {
  data <- data.frame(
    Y = c(NA, rnorm(19)),
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  # NA in outcome should cause model.frame to fail or be handled
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    regexp = NULL
  )
})

test_that("cf_fit handles NA values in treatment", {
  data <- data.frame(
    Y = rnorm(20),
    T = c(NA, rbinom(19, 1, 0.5)),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    regexp = NULL
  )
})

test_that("cf_fit handles NA values in covariates", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X1 = c(NA, rnorm(19))
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    regexp = NULL
  )
})

test_that("cf_fit handles empty data", {
  data <- data.frame(
    Y = numeric(0),
    T = integer(0),
    X1 = numeric(0)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    regexp = NULL
  )
})

test_that("cf_fit handles single observation", {
  data <- data.frame(
    Y = 1.0,
    T = 1L,
    X1 = 0.5
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    regexp = NULL
  )
})

test_that("cf_fit handles single treatment group (all treated)", {
  data <- data.frame(
    Y = rnorm(20),
    T = rep(1L, 20),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    regexp = NULL
  )
})

test_that("cf_fit handles single treatment group (all control)", {
  data <- data.frame(
    Y = rnorm(20),
    T = rep(0L, 20),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "gformula"),
    regexp = NULL
  )
})

test_that("cf_fit rejects unknown method", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "unknown_method"),
    "arg"
  )
})

test_that("cf_dowhy requires graph argument", {
  data <- data.frame(
    Y = rnorm(20),
    T = rbinom(20, 1, 0.5),
    X1 = rnorm(20)
  )
  expect_error(
    cf_fit(Y ~ T | X1, data = data, method = "dowhy_gcm"),
    "graph"
  )
})

test_that("parse_cf_formula extracts correct components", {
  set.seed(123)
  data <- data.frame(
    Y = rnorm(50),
    T = rbinom(50, 1, 0.5),
    X1 = rnorm(50),
    X2 = rnorm(50)
  )

  result <- parse_cf_formula(Y ~ T | X1 + X2, data)

  expect_equal(result$outcome_name, "Y")
  expect_equal(result$treatment_name, "T")
  expect_equal(result$covariate_names, c("X1", "X2"))
  expect_equal(length(result$Y), 50)
  expect_equal(length(result$T), 50)
  expect_equal(ncol(result$X), 2)
  expect_true(is.numeric(result$Y))
  expect_true(is.integer(result$T))
  expect_true(is.matrix(result$X))
})

test_that("as_cf_data.default rejects unsupported types", {
  expect_error(
    as_cf_data("not a data.frame"),
    "does not support"
  )
  expect_error(
    as_cf_data(list(a = 1, b = 2)),
    "does not support"
  )
})

test_that("cf_fit handles data with many covariates (p >> n)", {
  set.seed(42)
  n <- 30
  p <- 100
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("X", 1:p)

  data <- as.data.frame(X)
  data$Y <- rnorm(n)
  data$T <- rbinom(n, 1, 0.5)

  # GRF can handle high-dimensional data
  skip_if_not_installed("grf")

  formula <- as.formula(paste("Y ~ T |", paste(paste0("X", 1:p), collapse = " + ")))

  fit <- cf_fit(formula, data = data, method = "grf")
  expect_s3_class(fit, "cf_model")
  expect_equal(fit$meta$p, p)
})

test_that("cf_fit handles very small sample (n=5)", {
  skip_if_not_installed("grf")

  set.seed(123)
  data <- data.frame(
    Y = rnorm(5),
    T = c(0L, 0L, 1L, 1L, 1L),
    X1 = rnorm(5)
  )

  # Should work but may have warnings
  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula", nsimul = 10L)
  expect_s3_class(fit, "cf_model")
})

test_that("cf_fit handles zero-variance covariate", {
  set.seed(42)
  data <- data.frame(
    Y = rnorm(50),
    T = rbinom(50, 1, 0.5),
    X1 = rep(1, 50),  # constant
    X2 = rnorm(50)
  )

  # Should still fit (method may or may not use the constant variable)
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula", nsimul = 10L)
  expect_s3_class(fit, "cf_model")
})

test_that("select_best_method returns valid method", {
  method <- select_best_method()
  expect_true(method %in% c("grf", "hdml", "hdps", "bcf", "tmle", "ipw", "gformula", "dowhy_gcm", "drlearner"))
})
