#' Compute benchmark metrics for causal inference evaluation
#'
#' This function computes evaluation metrics comparing estimated treatment
#' effects against known true values from simulation data.
#'
#' @param ate_hat Numeric scalar, estimated average treatment effect
#' @param ite_hat Numeric vector, estimated individual treatment effects
#' @param summary_hat List (optional), summary statistics from predict(type="summary")
#'   containing ate_ci_lower and ate_ci_upper for coverage calculation
#' @param truth List containing true values:
#'   \itemize{
#'     \item ate_true: numeric scalar, true ATE
#'     \item ite_true: numeric vector, true ITEs
#'   }
#'
#' @return A named list with:
#'   \itemize{
#'     \item bias_ate: ATE bias (ate_hat - ate_true)
#'     \item abs_bias_ate: Absolute ATE bias
#'     \item squared_error_ate: Squared error of ATE (bias^2). Note: This is
#'       the squared error for a single estimate. True MSE requires averaging
#'       over replications at the analysis stage.
#'     \item pehe: Precision in Estimation of Heterogeneous Effect
#'       (sqrt(mean((ite_hat - ite_true)^2)))
#'     \item coverage_ate: 0/1 indicator if true ATE is within CI (NA if no CI)
#'     \item ci_len_ate: Length of confidence interval (NA if no CI)
#'   }
#'
#' @export
cf_benchmark_compute_metrics <- function(
    ate_hat,
    ite_hat,
    summary_hat = NULL,
    truth
) {
  ate_true <- truth$ate_true
  ite_true <- truth$ite_true
  
  if (is.null(ate_hat) || !is.numeric(ate_hat) || length(ate_hat) != 1) {
    ate_hat <- NA_real_
  }
  
  bias_ate <- ate_hat - ate_true
  abs_bias_ate <- abs(bias_ate)
  # Note: This is the squared error for a single estimate, not MSE.
  # True MSE requires averaging over replications: MSE = E[(estimate - true)^2]
  # Aggregation to MSE happens at the analysis/reporting stage.
  squared_error_ate <- bias_ate^2
  
  if (is.null(ite_hat) || !is.numeric(ite_hat) || length(ite_hat) != length(ite_true)) {
    pehe <- NA_real_
  } else {
    pehe <- sqrt(mean((ite_hat - ite_true)^2))
  }
  
  coverage_ate <- NA_real_
  ci_len_ate <- NA_real_
  
  if (!is.null(summary_hat)) {
    ci_lower <- summary_hat$ate_ci_lower
    ci_upper <- summary_hat$ate_ci_upper
    
    if (!is.null(ci_lower) && !is.null(ci_upper) &&
        is.numeric(ci_lower) && is.numeric(ci_upper)) {
      coverage_ate <- as.numeric(ate_true >= ci_lower && ate_true <= ci_upper)
      ci_len_ate <- ci_upper - ci_lower
    }
  }
  
  list(
    bias_ate = as.numeric(bias_ate),
    abs_bias_ate = as.numeric(abs_bias_ate),
    squared_error_ate = as.numeric(squared_error_ate),
    pehe = as.numeric(pehe),
    coverage_ate = as.numeric(coverage_ate),
    ci_len_ate = as.numeric(ci_len_ate)
  )
}

#' Compute computational metrics
#'
#' Measures execution time and peak memory usage for a function call.
#' This is useful for benchmarking the computational cost of different
#' causal inference methods.
#'
#' @param fit_fn Function to call for fitting
#' @param ... Arguments to pass to fit_fn
#'
#' @return A named list with:
#'   \itemize{
#'     \item result: The return value from fit_fn
#'     \item time_sec: Elapsed time in seconds
#'     \item peak_memory_mb: Approximate peak memory usage in MB
#'   }
#'
#' @export
cf_benchmark_compute_computational_metrics <- function(fit_fn, ...) {
  # Memory before
  gc(reset = TRUE)
  mem_before <- sum(gc()[, 2])

  # Time
  time_result <- system.time({
    result <- fit_fn(...)
  })

  # Memory after
  gc()
  mem_after <- sum(gc()[, 2])
  peak_memory_mb <- mem_after - mem_before

  list(
    result = result,
    time_sec = as.numeric(time_result["elapsed"]),
    peak_memory_mb = peak_memory_mb
  )
}
