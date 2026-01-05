# Tests for visualization functions
# Note: These functions require Python and matplotlib

# Helper to skip if Python visualization is not available
skip_if_no_viz <- function() {
  skip_if_no_python()
  skip_if_no_python_package("matplotlib")
  skip_if_no_python_package("numpy")
}

# --- cf_plot_ite tests ---

test_that("cf_plot_ite requires numeric input", {
  skip_if_no_viz()

  expect_error(
    cf_plot_ite("not numeric"),
    regexp = NULL
  )
})

test_that("cf_plot_ite works with numeric vector", {
  skip_if_no_viz()

  set.seed(42)
  ite <- rnorm(100, mean = 1.5, sd = 0.5)

  result <- cf_plot_ite(ite)

  # Should return invisibly
expect_true(is.list(result) || is.null(result))
})

test_that("cf_plot_ite can save to file", {
  skip_if_no_viz()

  set.seed(42)
  ite <- rnorm(100, mean = 1.5, sd = 0.5)

  temp_file <- tempfile(fileext = ".png")
  on.exit(unlink(temp_file), add = TRUE)

  cf_plot_ite(ite, save_path = temp_file)

  expect_true(file.exists(temp_file))
})

# --- cf_plot_outcome_shift tests ---

test_that("cf_plot_outcome_shift requires equal length vectors", {
  skip_if_no_viz()

  y0 <- rnorm(100)
  y1 <- rnorm(50)  # different length

  expect_error(
    cf_plot_outcome_shift(y0, y1),
    regexp = NULL
  )
})

test_that("cf_plot_outcome_shift works with matched vectors", {
  skip_if_no_viz()

  set.seed(42)
  y0 <- rnorm(100, mean = 5)
  y1 <- y0 + rnorm(100, mean = 1.5)

  result <- cf_plot_outcome_shift(y0, y1)

  expect_true(is.list(result) || is.null(result))
})

test_that("cf_plot_outcome_shift can save to file", {
  skip_if_no_viz()

  set.seed(42)
  y0 <- rnorm(100, mean = 5)
  y1 <- y0 + rnorm(100, mean = 1.5)

  temp_file <- tempfile(fileext = ".png")
  on.exit(unlink(temp_file), add = TRUE)

  cf_plot_outcome_shift(y0, y1, save_path = temp_file)

  expect_true(file.exists(temp_file))
})

# --- cf_plot_uplift tests ---

test_that("cf_plot_uplift works with numeric vector", {
  skip_if_no_viz()

  set.seed(42)
  ite <- rnorm(100, mean = 1.5, sd = 0.5)

  result <- cf_plot_uplift(ite)

  expect_true(is.list(result) || is.null(result))
})

test_that("cf_plot_uplift can save to file", {
  skip_if_no_viz()

  set.seed(42)
  ite <- rnorm(100, mean = 1.5, sd = 0.5)

  temp_file <- tempfile(fileext = ".png")
  on.exit(unlink(temp_file), add = TRUE)

  cf_plot_uplift(ite, save_path = temp_file)

  expect_true(file.exists(temp_file))
})

# --- cf_plot_ate_bootstrap tests ---

test_that("cf_plot_ate_bootstrap works with bootstrap samples", {
  skip_if_no_viz()

  set.seed(42)
  bootstrap_samples <- rnorm(500, mean = 1.5, sd = 0.2)

  result <- cf_plot_ate_bootstrap(bootstrap_samples)

  expect_true(is.list(result) || is.null(result))
})

test_that("cf_plot_ate_bootstrap with point estimate", {
  skip_if_no_viz()

  set.seed(42)
  bootstrap_samples <- rnorm(500, mean = 1.5, sd = 0.2)
  ate_point <- 1.5

  result <- cf_plot_ate_bootstrap(bootstrap_samples, ate_point = ate_point)

  expect_true(is.list(result) || is.null(result))
})

test_that("cf_plot_ate_bootstrap can save to file", {
  skip_if_no_viz()

  set.seed(42)
  bootstrap_samples <- rnorm(500, mean = 1.5, sd = 0.2)

  temp_file <- tempfile(fileext = ".png")
  on.exit(unlink(temp_file), add = TRUE)

  cf_plot_ate_bootstrap(bootstrap_samples, save_path = temp_file)

  expect_true(file.exists(temp_file))
})

# --- cf_plot_all tests ---

test_that("cf_plot_all requires cf_model object", {
  expect_error(
    cf_plot_all("not a model"),
    "cf_model"
  )

  expect_error(
    cf_plot_all(list(a = 1)),
    "cf_model"
  )
})

test_that("cf_plot_all works with cf_model", {
  skip_if_no_viz()

  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  data$Y <- 1.5 * data$T + 0.5 * data$X1 + rnorm(n)

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula", nsimul = 20L)

  result <- cf_plot_all(fit)

  expect_true(is.list(result))
  expect_true("ite" %in% names(result))
  expect_true("outcome_shift" %in% names(result))
  expect_true("uplift" %in% names(result))
})

test_that("cf_plot_all saves to directory", {
  skip_if_no_viz()

  set.seed(42)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )
  data$Y <- 1.5 * data$T + 0.5 * data$X1 + rnorm(n)

  fit <- cf_fit(Y ~ T | X1, data = data, method = "gformula", nsimul = 20L)

  temp_dir <- tempdir()
  output_dir <- file.path(temp_dir, "test_plots")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  cf_plot_all(fit, output_dir = output_dir, prefix = "test_")

  expect_true(dir.exists(output_dir))
  expect_true(file.exists(file.path(output_dir, "test_ite_distribution.png")))
  expect_true(file.exists(file.path(output_dir, "test_outcome_shift.png")))
  expect_true(file.exists(file.path(output_dir, "test_uplift_curve.png")))
})

# --- figsize parameter tests ---

test_that("cf_plot_ite accepts custom figsize", {
  skip_if_no_viz()

  set.seed(42)
  ite <- rnorm(100)

  # Should not error with custom figsize
  result <- cf_plot_ite(ite, figsize = c(8, 5))

  expect_true(is.list(result) || is.null(result))
})

test_that("cf_plot_outcome_shift accepts custom figsize", {
  skip_if_no_viz()

  set.seed(42)
  y0 <- rnorm(100)
  y1 <- y0 + 1

  result <- cf_plot_outcome_shift(y0, y1, figsize = c(12, 8))

  expect_true(is.list(result) || is.null(result))
})

# --- Integration with fitted models ---

test_that("visualization works end-to-end with grf", {
  skip_if_no_viz()
  skip_if_not_installed("grf")

  set.seed(42)
  n <- 150
  data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  data$T <- rbinom(n, 1, plogis(0.3 * data$X1))
  data$Y <- 1.5 * data$T + 0.8 * data$X1 + 0.3 * data$X2 + rnorm(n)

  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")

  # Get predictions
  ite <- predict(fit, type = "ite")
  y0 <- predict(fit, type = "y0")
  y1 <- predict(fit, type = "y1")

  # Plot each type
  expect_silent(cf_plot_ite(ite))
  expect_silent(cf_plot_outcome_shift(y0, y1))
  expect_silent(cf_plot_uplift(ite))
})
