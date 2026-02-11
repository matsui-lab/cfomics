test_that("Results are reproducible with same random_state (grf)", {
  skip_if_not_installed("grf")

  set.seed(999)
  data <- data.frame(
    Y = rnorm(100),
    T = rbinom(100, 1, 0.5),
    X1 = rnorm(100)
  )

  fit1 <- cf_fit(Y ~ T | X1, data = data, method = "grf", random_state = 42)
  fit2 <- cf_fit(Y ~ T | X1, data = data, method = "grf", random_state = 42)

  expect_equal(fit1$fit$res$ate, fit2$fit$res$ate)
  expect_equal(fit1$fit$res$ite, fit2$fit$res$ite)
})

test_that("Results are reproducible with gformula (deterministic method)", {
  # gformula is deterministic given the same data, so results should be identical
  set.seed(123)
  data <- data.frame(
    Y = rnorm(100),
    T = rbinom(100, 1, 0.5),
    X1 = rnorm(100)
  )

  fit1 <- cf_fit(Y ~ T | X1, data = data, method = "gformula")
  fit2 <- cf_fit(Y ~ T | X1, data = data, method = "gformula")

  expect_equal(fit1$fit$res$ate, fit2$fit$res$ate)
  expect_equal(fit1$fit$res$ite, fit2$fit$res$ite)
})

test_that("meta contains software versions (grf)", {
  skip_if_not_installed("grf")

  set.seed(123)
  data <- data.frame(Y = rnorm(50), T = rbinom(50, 1, 0.5), X1 = rnorm(50))
  fit <- cf_fit(Y ~ T | X1, data = data, method = "grf")

  expect_true("software_versions" %in% names(fit$meta))
  expect_true("R" %in% names(fit$meta$software_versions))
  expect_true("cfomics" %in% names(fit$meta$software_versions))
  expect_true("timestamp" %in% names(fit$meta$software_versions))
})

test_that("meta contains software versions (gformula)", {
  set.seed(123)
  data <- data.frame(Y = rnorm(50), T = rbinom(50, 1, 0.5), X1 = rnorm(50))
  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula")

  expect_true("software_versions" %in% names(fit$meta))
  expect_true("R" %in% names(fit$meta$software_versions))
  expect_true("cfomics" %in% names(fit$meta$software_versions))
  expect_true("timestamp" %in% names(fit$meta$software_versions))

  # Verify R version format

  expect_match(fit$meta$software_versions$R, "^[0-9]+\\.[0-9]+")

  # Verify timestamp is a POSIXct object
  expect_s3_class(fit$meta$software_versions$timestamp, "POSIXct")
})

test_that("software_versions includes method-specific packages", {
  skip_if_not_installed("grf")

  set.seed(123)
  data <- data.frame(Y = rnorm(50), T = rbinom(50, 1, 0.5), X1 = rnorm(50))
  fit <- cf_fit(Y ~ T | X1, data = data, method = "grf")

  # Check that packages list includes grf version
  expect_true("packages" %in% names(fit$meta$software_versions))
  expect_true("grf" %in% names(fit$meta$software_versions$packages))
})
