#' Sensitivity Analysis for Unmeasured Confounding
#'
#' Computes the E-value (VanderWeele & Ding, 2017) to quantify robustness
#' of causal effect estimates to unmeasured confounding.
#'
#' @param result A cfomics_result object from cf_fit()
#' @param type Character, type of sensitivity analysis. Currently only "evalue".
#' @param alpha Numeric, significance level for CI-based E-value (default 0.05)
#' @return A list of class "cf_sensitivity"
#' @export
#' @examples
#' \donttest{
#' # Create example data and fit model
#' set.seed(123)
#' n <- 100
#' demo_data <- data.frame(
#'   X1 = rnorm(n),
#'   X2 = rnorm(n),
#'   T = rbinom(n, 1, 0.5)
#' )
#' demo_data$Y <- 2 * demo_data$T + 0.5 * demo_data$X1 + rnorm(n)
#'
#' # Fit model and run sensitivity analysis
#' fit <- cf_fit(Y ~ T | X1 + X2, data = demo_data, method = "gformula")
#' sens <- cf_sensitivity(fit)
#' print(sens)
#' }
cf_sensitivity <- function(result, type = "evalue", alpha = 0.05) {
  if (!inherits(result, "cfomics_result")) {
    rlang::abort("Input must be a cfomics_result object",
                 class = "cfomics_invalid_input")
  }

  type <- match.arg(type, choices = "evalue")

  ate <- result$fit$res$ate
  ci_lower <- result$fit$res$summary$ate_ci_lower
  ci_upper <- result$fit$res$summary$ate_ci_upper

  # Compute SD of observed outcomes from counterfactuals
  y0 <- result$fit$res$y0_hat
  y1 <- result$fit$res$y1_hat
  sd_y <- stats::sd(c(y0, y1))
  if (is.na(sd_y) || sd_y == 0) sd_y <- 1

  # Convert ATE to approximate risk ratio
  rr_point <- exp(abs(ate) / sd_y)
  # Use the CI bound closer to null for conservative E-value
  if (ci_lower > 0) {
    rr_ci <- exp(ci_lower / sd_y)
  } else if (ci_upper < 0) {
    rr_ci <- exp(abs(ci_upper) / sd_y)
  } else {
    # CI crosses null
    rr_ci <- 1
  }

  # E-value formula
  evalue_fn <- function(rr) {
    if (rr <= 1) return(1)
    rr + sqrt(rr * (rr - 1))
  }

  ev_point <- evalue_fn(rr_point)
  ev_ci <- evalue_fn(rr_ci)

  # Interpretation
  interp <- if (ev_point > 2) {
    "Robust: an unmeasured confounder would need strong associations with both treatment and outcome to explain away this effect."
  } else if (ev_point > 1.5) {
    "Moderate robustness: an unmeasured confounder with moderate associations could potentially explain this effect."
  } else {
    "Weak robustness: even a weak unmeasured confounder could explain away this effect."
  }

  result_out <- list(
    evalue_point = ev_point,
    evalue_ci = ev_ci,
    ate = ate,
    ci = c(lower = ci_lower, upper = ci_upper),
    rr_point = rr_point,
    sd_y = sd_y,
    interpretation = interp,
    method = result$method
  )

  class(result_out) <- c("cf_sensitivity", "list")
  result_out
}

#' @export
print.cf_sensitivity <- function(x, ...) {
  cat("\nSensitivity Analysis (E-value)\n")
  cat("==============================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("ATE: %.3f [%.3f, %.3f]\n", x$ate, x$ci["lower"], x$ci["upper"]))
  cat(sprintf("E-value (point): %.2f\n", x$evalue_point))
  cat(sprintf("E-value (CI):    %.2f\n", x$evalue_ci))
  cat(sprintf("\n%s\n", x$interpretation))
  invisible(x)
}

#' @export
plot.cf_sensitivity <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for sensitivity plot")
    return(invisible(NULL))
  }

  # Bias plot: range of confounder-treatment and confounder-outcome RRs
  rr_seq <- seq(1, max(x$evalue_point * 1.5, 3), length.out = 100)
  bias_data <- data.frame(
    rr_outcome = rr_seq,
    rr_treatment_point = vapply(rr_seq, function(rr_uy) {
      if (rr_uy <= 1) return(Inf)
      x$rr_point * (rr_uy - 1 + x$rr_point) / (rr_uy * x$rr_point - rr_uy + 1)
    }, numeric(1))
  )
  bias_data$rr_treatment_point[bias_data$rr_treatment_point < 1] <- NA
  bias_data$rr_treatment_point[bias_data$rr_treatment_point > max(rr_seq) * 2] <- NA

  p <- ggplot2::ggplot(bias_data, ggplot2::aes(
    x = .data[["rr_outcome"]], y = .data[["rr_treatment_point"]]
  )) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_point(
      data = data.frame(x = x$evalue_point, y = x$evalue_point),
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      color = "red", size = 3
    ) +
    ggplot2::annotate("text", x = x$evalue_point, y = x$evalue_point,
                      label = sprintf("E = %.1f", x$evalue_point),
                      vjust = -1, color = "red") +
    ggplot2::labs(x = "Confounder-Outcome RR",
                  y = "Confounder-Treatment RR",
                  title = "Sensitivity Analysis: E-value Bias Plot",
                  subtitle = sprintf("ATE = %.3f, E-value = %.2f",
                                     x$ate, x$evalue_point)) +
    ggplot2::theme_minimal()
  print(p)
  invisible(p)
}
