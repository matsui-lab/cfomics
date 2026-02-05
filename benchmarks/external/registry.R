# benchmarks/external/registry.R
# Registry for external package wrappers

#' Run external method wrapper
#'
#' Unified interface for external causal inference packages.
#' All wrappers return the same structure for fair comparison.
#'
#' @param method Character, external method name
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param ... Additional method-specific arguments
#' @return List with: ate, ite, y0_hat, y1_hat, ci_lower, ci_upper, time_sec
run_external_method <- function(method, X, T, Y, ...) {
  wrapper_fn <- get_external_wrapper(method)
  if (is.null(wrapper_fn)) {
    stop(sprintf("Unknown external method: %s", method))
  }

  start_time <- Sys.time()
  result <- wrapper_fn(X = X, T = T, Y = Y, ...)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Validate result structure
  required_fields <- c("ate", "ite")
  missing <- setdiff(required_fields, names(result))
  if (length(missing) > 0) {
    stop(sprintf("Wrapper for %s missing required fields: %s",
                 method, paste(missing, collapse = ", ")))
  }

  # Add defaults for optional fields
  result$y0_hat <- result$y0_hat %||% rep(NA_real_, length(Y))
  result$y1_hat <- result$y1_hat %||% rep(NA_real_, length(Y))
  result$ci_lower <- result$ci_lower %||% NA_real_
  result$ci_upper <- result$ci_upper %||% NA_real_

  result$time_sec <- elapsed

  result
}

#' Get wrapper function for external method
#' @param method Character, method name
#' @return Function or NULL if not found
get_external_wrapper <- function(method) {
  wrappers <- list(
    "matchit" = wrapper_matchit,
    "weightit" = wrapper_weightit,
    "hdm" = wrapper_hdm,
    "doubleml" = wrapper_doubleml,
    "tmle3" = wrapper_tmle3,
    "bart" = wrapper_bart,
    "superlearner" = wrapper_superlearner
  )
  wrappers[[method]]
}

#' List available external methods
#' @return Character vector of method names
list_external_methods <- function() {
  c("matchit", "weightit", "hdm", "doubleml", "tmle3", "bart", "superlearner")
}

#' Check if external method dependencies are available
#' @param method Character, method name
#' @return NULL if available, character string with reason if not
check_external_dependencies <- function(method) {
  deps <- list(
    matchit = c("MatchIt", "marginaleffects"),
    weightit = c("WeightIt", "marginaleffects"),
    hdm = c("hdm"),
    doubleml = c("DoubleML", "mlr3", "mlr3learners"),
    tmle3 = c("tmle3", "sl3"),
    bart = c("dbarts"),
    superlearner = c("SuperLearner")
  )

  required <- deps[[method]]
  if (is.null(required)) return(sprintf("Unknown method: %s", method))

  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    return(sprintf("Missing packages: %s", paste(missing, collapse = ", ")))
  }
  NULL
}

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
