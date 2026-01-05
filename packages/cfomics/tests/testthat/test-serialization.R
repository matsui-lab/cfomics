# Test serialization (save/load) of cfomics results
# R-native methods should be fully serializable
# Python methods have documented limitations

test_that("R-native gformula results can be saved and restored", {
  # Create simple test data
  set.seed(123)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  data$Y <- data$Y + 2 * data$T + 0.5 * data$X1

  # Fit model
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula")

  # Save and restore
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  # Verify structure preserved
  expect_s3_class(fit2, "cfomics_result")
  expect_equal(fit2$method, "gformula")

  # Verify predictions work
  ate1 <- predict(fit, type = "ate")
  ate2 <- predict(fit2, type = "ate")
  expect_equal(ate1, ate2, tolerance = 1e-10)

  ite1 <- predict(fit, type = "ite")
  ite2 <- predict(fit2, type = "ite")
  expect_equal(ite1, ite2, tolerance = 1e-10)

  # Verify summary works
  expect_no_error(summary(fit2))
})

test_that("R-native ipw results can be saved and restored", {
  skip_if_not_installed("ipw")
  skip_if_not_installed("survey")

  set.seed(456)
  n <- 50
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )
  data$Y <- data$Y + 1.5 * data$T

  fit <- cf_fit(Y ~ T | X1, data = data, method = "ipw")

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  expect_s3_class(fit2, "cfomics_result")
  expect_equal(predict(fit, type = "ate"), predict(fit2, type = "ate"), tolerance = 1e-10)
})

test_that("R-native grf results can be saved and restored", {
  skip_if_not_installed("grf")

  set.seed(789)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  data$Y <- data$Y + 2 * data$T + 0.3 * data$X1

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  expect_s3_class(fit2, "cfomics_result")

  # ATE should be preserved
  ate1 <- predict(fit, type = "ate")
  ate2 <- predict(fit2, type = "ate")
  expect_equal(ate1, ate2, tolerance = 1e-10)

  # ITE should be preserved
  ite1 <- predict(fit, type = "ite")
  ite2 <- predict(fit2, type = "ite")
  expect_equal(ite1, ite2, tolerance = 1e-10)
})

test_that("cf_export produces loadable results", {
  set.seed(111)
  n <- 30
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X = rnorm(n)
  )
  data$Y <- data$Y + 2 * data$T

  fit <- cf_fit(Y ~ T | X, data = data, method = "gformula")

  # Test RDS export
  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_rds), add = TRUE)

  cf_export(fit, file = tmp_rds, format = "rds")
  expect_true(file.exists(tmp_rds))

  loaded <- readRDS(tmp_rds)
  expect_true(is.list(loaded))
  expect_true("ate" %in% names(loaded) || "summary" %in% names(loaded))

  # Test CSV export
  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv), add = TRUE)

  cf_export(fit, file = tmp_csv, format = "csv")
  expect_true(file.exists(tmp_csv))
})

test_that("meta information is preserved after serialization", {
  set.seed(222)
  n <- 40
  data <- data.frame(
    outcome = rnorm(n),
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 50, 10)
  )
  data$outcome <- data$outcome + 1.5 * data$treatment

  fit <- cf_fit(outcome ~ treatment | age, data = data, method = "gformula",
                random_state = 42)

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  fit2 <- readRDS(tmp)

  # Meta should be fully preserved
  expect_equal(fit2$meta$n, fit$meta$n)
  expect_equal(fit2$meta$p, fit$meta$p)
  expect_equal(fit2$meta$outcome_name, fit$meta$outcome_name)
  expect_equal(fit2$meta$treatment_name, fit$meta$treatment_name)
  expect_equal(fit2$meta$covariate_names, fit$meta$covariate_names)
  expect_equal(fit2$meta$random_state, fit$meta$random_state)
  expect_equal(fit2$meta$estimand, fit$meta$estimand)
})
