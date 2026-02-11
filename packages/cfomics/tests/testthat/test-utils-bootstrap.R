# Tests for bootstrap utilities

test_that(".bootstrap_ci returns correct structure", {
  set.seed(42)
  data <- data.frame(x = rnorm(100))
  result <- .bootstrap_ci(data, function(d) mean(d$x), n_boot = 100)

  expect_true(is.list(result))
  expect_named(result, c("estimate", "ci_lower", "ci_upper", "se", "boot_samples"))
  expect_true(is.numeric(result$estimate))
  expect_true(is.numeric(result$ci_lower))

  expect_true(is.numeric(result$ci_upper))
  expect_true(is.numeric(result$se))
  expect_length(result$boot_samples, 100)
})

test_that(".bootstrap_ci confidence interval contains estimate typically", {
  set.seed(42)
  data <- data.frame(x = rnorm(100, mean = 5))
  result <- .bootstrap_ci(data, function(d) mean(d$x), n_boot = 200)

  # CI should typically contain the estimate for well-behaved data
  expect_true(result$ci_lower <= result$estimate)
  expect_true(result$ci_upper >= result$estimate)
})

test_that(".bootstrap_ci is reproducible with seed", {
  data <- data.frame(x = rnorm(50))
  r1 <- .bootstrap_ci(data, function(d) mean(d$x), n_boot = 50, seed = 123)
  r2 <- .bootstrap_ci(data, function(d) mean(d$x), n_boot = 50, seed = 123)

  expect_equal(r1$estimate, r2$estimate)
  expect_equal(r1$ci_lower, r2$ci_lower)
  expect_equal(r1$ci_upper, r2$ci_upper)
  expect_equal(r1$se, r2$se)
  expect_equal(r1$boot_samples, r2$boot_samples)
})

test_that(".bootstrap_ci respects conf_level parameter", {
  set.seed(42)
  data <- data.frame(x = rnorm(200))

  result_90 <- .bootstrap_ci(data, function(d) mean(d$x), n_boot = 500,
                              conf_level = 0.90, seed = 456)
  result_99 <- .bootstrap_ci(data, function(d) mean(d$x), n_boot = 500,
                              conf_level = 0.99, seed = 456)

  # 99% CI should be wider than 90% CI
  width_90 <- result_90$ci_upper - result_90$ci_lower
  width_99 <- result_99$ci_upper - result_99$ci_lower
  expect_true(width_99 > width_90)
})

test_that(".bootstrap_ci validates data argument", {
  expect_error(
    .bootstrap_ci("not a dataframe", function(d) mean(d$x)),
    class = "cfomics_bootstrap_error"
  )

  expect_error(
    .bootstrap_ci(list(x = 1:10), function(d) mean(d$x)),
    class = "cfomics_bootstrap_error"
  )
})

test_that(".bootstrap_ci validates statistic_fn argument", {
  data <- data.frame(x = 1:10)

  expect_error(
    .bootstrap_ci(data, "not a function"),
    class = "cfomics_bootstrap_error"
  )

  expect_error(
    .bootstrap_ci(data, NULL),
    class = "cfomics_bootstrap_error"
  )
})

test_that(".bootstrap_ci validates n_boot argument", {
  data <- data.frame(x = 1:10)
  stat_fn <- function(d) mean(d$x)

  expect_error(
    .bootstrap_ci(data, stat_fn, n_boot = -1),
    class = "cfomics_bootstrap_error"
  )

  expect_error(
    .bootstrap_ci(data, stat_fn, n_boot = "100"),
    class = "cfomics_bootstrap_error"
  )

  expect_error(
    .bootstrap_ci(data, stat_fn, n_boot = c(100, 200)),
    class = "cfomics_bootstrap_error"
  )
})

test_that(".bootstrap_ci validates conf_level argument", {
  data <- data.frame(x = 1:10)
  stat_fn <- function(d) mean(d$x)

  expect_error(
    .bootstrap_ci(data, stat_fn, conf_level = 0),
    class = "cfomics_bootstrap_error"
  )

  expect_error(
    .bootstrap_ci(data, stat_fn, conf_level = 1),
    class = "cfomics_bootstrap_error"
  )

  expect_error(
    .bootstrap_ci(data, stat_fn, conf_level = 1.5),
    class = "cfomics_bootstrap_error"
  )

  expect_error(
    .bootstrap_ci(data, stat_fn, conf_level = "0.95"),
    class = "cfomics_bootstrap_error"
  )
})

test_that(".bootstrap_ci validates statistic_fn return value", {
  data <- data.frame(x = 1:10)

  # Function that returns vector
  expect_error(
    .bootstrap_ci(data, function(d) d$x),
    class = "cfomics_bootstrap_error"
  )

  # Function that returns character
  expect_error(
    .bootstrap_ci(data, function(d) "text"),
    class = "cfomics_bootstrap_error"
  )
})

test_that(".gformula_ate_statistic computes correct ATE", {
  # Create simple data where true ATE is known
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  # True effect of A is 2
  Y <- 1 + 0.5 * x + 2 * A + rnorm(n, sd = 0.1)
  data <- data.frame(Y = Y, A = A, x = x)

  outcome_formula <- Y ~ A + x

  ate <- .gformula_ate_statistic(data, outcome_formula)

  # ATE should be close to 2
  expect_true(abs(ate - 2) < 0.2)
})

test_that(".gformula_ate_statistic works with bootstrap_ci", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.5 * x + 1.5 * A + rnorm(n, sd = 0.5)
  data <- data.frame(Y = Y, A = A, x = x)

  outcome_formula <- Y ~ A + x

  # Use .bootstrap_ci with .gformula_ate_statistic
  result <- .bootstrap_ci(
    data,
    function(d) .gformula_ate_statistic(d, outcome_formula),
    n_boot = 100,
    seed = 123
  )

  expect_true(is.list(result))
  expect_true(result$ci_lower < result$ci_upper)
  # True ATE is 1.5, should be within reasonable CI
  expect_true(result$ci_lower < 2.5)
  expect_true(result$ci_upper > 0.5)
})
