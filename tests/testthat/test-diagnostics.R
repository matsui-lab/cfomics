# Tests for diagnostic functions (balance check, overlap check)

# --- cf_balance_check tests ---

test_that("cf_balance_check computes SMD correctly", {
  set.seed(42)
  n <- 200
  data <- data.frame(
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 50, 10),
    income = rnorm(n, 50000, 10000)
  )

  balance <- cf_balance_check(data, "treatment", c("age", "income"))

  expect_s3_class(balance, "cf_balance")
  expect_s3_class(balance, "data.frame")
  expect_equal(nrow(balance), 2)
  expect_true("smd" %in% names(balance))
  expect_true("balanced" %in% names(balance))
  expect_true("variable" %in% names(balance))
})

test_that("cf_balance_check detects imbalance", {
  set.seed(42)
  n <- 200

  # Create imbalanced data where treated have higher age
  data <- data.frame(
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 50, 10),
    income = rnorm(n, 50000, 10000)
  )
  # Make age highly imbalanced
  data$age[data$treatment == 1] <- data$age[data$treatment == 1] + 20

  balance <- cf_balance_check(data, "treatment", c("age", "income"), threshold = 0.1)

  # Age should be imbalanced
  age_row <- balance[balance$variable == "age", ]
  expect_false(age_row$balanced)
  expect_true(abs(age_row$smd) > 0.1)
})

test_that("cf_balance_check auto-detects numeric covariates", {
  set.seed(42)
  n <- 100
  data <- data.frame(
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 50, 10),
    income = rnorm(n, 50000, 10000),
    name = sample(letters, n, replace = TRUE)  # character column
  )

  balance <- cf_balance_check(data, "treatment")

  # Should only include numeric columns (age, income)
  expect_equal(nrow(balance), 2)
  expect_true(all(c("age", "income") %in% balance$variable))
  expect_false("name" %in% balance$variable)
})

test_that("cf_balance_check respects threshold parameter", {
  set.seed(42)
  n <- 200
  data <- data.frame(
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 50, 10)
  )
  # Add moderate imbalance
  data$age[data$treatment == 1] <- data$age[data$treatment == 1] + 2

  # With strict threshold, should be imbalanced
  balance_strict <- cf_balance_check(data, "treatment", "age", threshold = 0.05)

  # With lenient threshold, might be balanced

  balance_lenient <- cf_balance_check(data, "treatment", "age", threshold = 0.5)

  expect_true(attr(balance_strict, "threshold") == 0.05)
  expect_true(attr(balance_lenient, "threshold") == 0.5)
})

test_that("cf_balance_check rejects missing treatment variable", {
  data <- data.frame(
    age = rnorm(100),
    income = rnorm(100)
  )

  expect_error(
    cf_balance_check(data, "treatment", "age"),
    "not found"
  )
})

test_that("cf_balance_check rejects no numeric covariates", {
  data <- data.frame(
    treatment = rbinom(100, 1, 0.5),
    name = sample(letters, 100, replace = TRUE)
  )

  expect_error(
    cf_balance_check(data, "treatment"),
    "No numeric"
  )
})

test_that("print.cf_balance produces output", {
  set.seed(42)
  data <- data.frame(
    treatment = rbinom(100, 1, 0.5),
    age = rnorm(100, 50, 10),
    income = rnorm(100, 50000, 10000)
  )

  balance <- cf_balance_check(data, "treatment")

  expect_output(print(balance), "Covariate Balance")
  expect_output(print(balance), "Treatment")
  expect_output(print(balance), "age")
})

test_that("plot.cf_balance works with ggplot2", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  data <- data.frame(
    treatment = rbinom(100, 1, 0.5),
    age = rnorm(100, 50, 10),
    income = rnorm(100, 50000, 10000)
  )

  balance <- cf_balance_check(data, "treatment")
  p <- plot(balance)

  expect_s3_class(p, "ggplot")
})

test_that("plot.cf_balance works with base R", {
  set.seed(42)
  data <- data.frame(
    treatment = rbinom(100, 1, 0.5),
    age = rnorm(100, 50, 10),
    income = rnorm(100, 50000, 10000)
  )

  balance <- cf_balance_check(data, "treatment")

  # Temporarily hide ggplot2 if available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    skip("ggplot2 is available, base R plot not tested")
  }

  # Should not error
  expect_invisible(plot(balance))
})

# --- cf_overlap_check tests ---

test_that("cf_overlap_check computes propensity scores", {
  set.seed(42)
  n <- 200
  data <- data.frame(
    age = rnorm(n, 50, 10),
    income = rnorm(n, 50000, 10000)
  )
  # Treatment assignment related to covariates
  prob <- plogis(-2 + 0.02 * data$age + 0.00001 * data$income)
  data$treatment <- rbinom(n, 1, prob)

  overlap <- cf_overlap_check(data, "treatment", c("age", "income"))

  expect_s3_class(overlap, "cf_overlap")
  expect_s3_class(overlap, "list")
  expect_equal(length(overlap$propensity_scores), n)
  expect_true(all(overlap$propensity_scores >= 0 & overlap$propensity_scores <= 1))
})

test_that("cf_overlap_check summary structure is correct", {
  set.seed(42)
  n <- 200
  data <- data.frame(
    age = rnorm(n, 50, 10),
    treatment = rbinom(n, 1, 0.5)
  )

  overlap <- cf_overlap_check(data, "treatment", "age")

  expect_true("summary" %in% names(overlap))
  expect_equal(nrow(overlap$summary), 2)
  expect_true(all(c("group", "n", "min_ps", "max_ps", "mean_ps", "median_ps") %in% names(overlap$summary)))
})

test_that("cf_overlap_check detects extreme propensity scores", {
  set.seed(42)
  n <- 200

  # Create data with near-perfect separation
  data <- data.frame(
    x = c(rnorm(100, -5), rnorm(100, 5)),
    treatment = c(rep(0, 100), rep(1, 100))
  )

  overlap <- cf_overlap_check(data, "treatment", "x")

  # Should detect extreme PS values
  expect_true(overlap$violations$extreme_ps)
})

test_that("cf_overlap_check identifies overlap region", {
  set.seed(42)
  n <- 200
  data <- data.frame(
    age = rnorm(n, 50, 10),
    treatment = rbinom(n, 1, 0.5)
  )

  overlap <- cf_overlap_check(data, "treatment", "age")

  expect_true("overlap_region" %in% names(overlap))
  expect_true(overlap$overlap_region["lower"] <= overlap$overlap_region["upper"])
})

test_that("print.cf_overlap produces output", {
  set.seed(42)
  data <- data.frame(
    age = rnorm(100, 50, 10),
    treatment = rbinom(100, 1, 0.5)
  )

  overlap <- cf_overlap_check(data, "treatment", "age")

  expect_output(print(overlap), "Overlap")
  expect_output(print(overlap), "Control")
  expect_output(print(overlap), "Treated")
})

test_that("print.cf_overlap warns about violations", {
  set.seed(42)
  n <- 200

  # Create data with extreme separation
  data <- data.frame(
    x = c(rnorm(100, -10), rnorm(100, 10)),
    treatment = c(rep(0, 100), rep(1, 100))
  )

  overlap <- cf_overlap_check(data, "treatment", "x")

  expect_output(print(overlap), "WARNING|extreme")
})

test_that("plot.cf_overlap works with ggplot2", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  data <- data.frame(
    age = rnorm(100, 50, 10),
    treatment = rbinom(100, 1, 0.5)
  )

  overlap <- cf_overlap_check(data, "treatment", "age")
  p <- plot(overlap)

  expect_s3_class(p, "ggplot")
})

test_that("cf_overlap_check model is returned", {
  set.seed(42)
  data <- data.frame(
    age = rnorm(100, 50, 10),
    treatment = rbinom(100, 1, 0.5)
  )

  overlap <- cf_overlap_check(data, "treatment", "age")

  expect_true("model" %in% names(overlap))
  expect_s3_class(overlap$model, "glm")
})

# --- Integration with SummarizedExperiment ---

test_that("cf_balance_check works with SummarizedExperiment",
{
  skip_if_not_installed("SummarizedExperiment")

  set.seed(42)
  n <- 50
  counts <- matrix(rnorm(100 * n), nrow = 100, ncol = n)
  rownames(counts) <- paste0("gene", 1:100)
  colnames(counts) <- paste0("sample", 1:n)

  col_data <- data.frame(
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 50, 10),
    row.names = colnames(counts)
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )

  balance <- cf_balance_check(se, "treatment", "age")

  expect_s3_class(balance, "cf_balance")
  expect_equal(nrow(balance), 1)
})

test_that("cf_overlap_check works with SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")

  set.seed(42)
  n <- 50
  counts <- matrix(rnorm(100 * n), nrow = 100, ncol = n)
  rownames(counts) <- paste0("gene", 1:100)
  colnames(counts) <- paste0("sample", 1:n)

  col_data <- data.frame(
    treatment = rbinom(n, 1, 0.5),
    age = rnorm(n, 50, 10),
    row.names = colnames(counts)
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )

  overlap <- cf_overlap_check(se, "treatment", "age")

  expect_s3_class(overlap, "cf_overlap")
  expect_equal(length(overlap$propensity_scores), n)
})
