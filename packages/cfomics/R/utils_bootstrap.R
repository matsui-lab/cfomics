# Bootstrap Utilities for cfomics
#
# Internal functions for bootstrap-based confidence interval estimation.
# These utilities support various causal inference methods that require
# uncertainty quantification.

#' Internal bootstrap utility for confidence interval estimation
#'
#' Performs bootstrap resampling to estimate confidence intervals for a
#' statistic computed from data. Uses the percentile method for CI calculation.
#'
#' @param data Data frame to bootstrap
#' @param statistic_fn Function taking a data frame and returning a single numeric value.
#'   The function should compute the statistic of interest on the resampled data.
#' @param n_boot Number of bootstrap replicates (default 1000)
#' @param conf_level Confidence level for the interval (default 0.95)
#' @param seed Random seed for reproducibility (optional)
#'
#' @return A list with components:
#' \describe{
#'   \item{estimate}{The statistic computed on the original data}
#'   \item{ci_lower}{Lower bound of the confidence interval}
#'   \item{ci_upper}{Upper bound of the confidence interval}
#'   \item{se}{Bootstrap standard error}
#'   \item{boot_samples}{Vector of bootstrap replicate values}
#' }
#'
#' @examples
#' \dontrun{
#' data <- data.frame(x = rnorm(100))
#' result <- .bootstrap_ci(data, function(d) mean(d$x), n_boot = 500)
#' }
#'
#' @keywords internal
#' @noRd
.bootstrap_ci <- function(data, statistic_fn, n_boot = 1000,
                          conf_level = 0.95, seed = NULL) {
  # Input validation
  if (!is.data.frame(data)) {
    rlang::abort(
      message = "`data` must be a data frame",
      class = "cfomics_bootstrap_error"
    )
  }

  if (!is.function(statistic_fn)) {
    rlang::abort(
      message = "`statistic_fn` must be a function",
      class = "cfomics_bootstrap_error"
    )
  }

  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1) {
    rlang::abort(
      message = "`n_boot` must be a positive integer",
      class = "cfomics_bootstrap_error"
    )
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    rlang::abort(
      message = "`conf_level` must be a number between 0 and 1 (exclusive)",
      class = "cfomics_bootstrap_error"
    )
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- nrow(data)

  # Compute estimate on original data
  estimate <- statistic_fn(data)

  if (!is.numeric(estimate) || length(estimate) != 1) {
    rlang::abort(
      message = "`statistic_fn` must return a single numeric value",
      class = "cfomics_bootstrap_error"
    )
  }

  # Bootstrap resampling
  boot_samples <- vapply(seq_len(n_boot), function(i) {
    idx <- sample(n, n, replace = TRUE)
    statistic_fn(data[idx, , drop = FALSE])
  }, numeric(1))

  # Calculate confidence interval using percentile method
  alpha <- 1 - conf_level
  ci <- stats::quantile(boot_samples, probs = c(alpha / 2, 1 - alpha / 2),
                        na.rm = TRUE)

  list(
    estimate = estimate,
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2]),
    se = stats::sd(boot_samples, na.rm = TRUE),
    boot_samples = boot_samples
  )
}


#' Bootstrap confidence interval for ATE using G-formula
#'
#' Internal function that computes a single bootstrap replicate of the ATE
#' using the G-formula approach. This is used by \code{cf_fit_gformula()}.
#'
#' @param data Data frame with columns Y (outcome), A (treatment), and covariates
#' @param outcome_formula Formula for the outcome model
#'
#' @return Numeric value representing the ATE for this bootstrap sample
#'
#' @keywords internal
#' @noRd
.gformula_ate_statistic <- function(data, outcome_formula) {
  # Fit outcome model
  model <- stats::lm(outcome_formula, data = data)

  # Predict under control (A=0)
  data_control <- data
  data_control$A <- 0
  y0 <- stats::predict(model, newdata = data_control)

  # Predict under treatment (A=1)
  data_treated <- data
  data_treated$A <- 1
  y1 <- stats::predict(model, newdata = data_treated)

  # Return ATE

  mean(y1 - y0)
}
