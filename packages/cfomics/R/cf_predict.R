#' Predict treatment effects from a cf_model
#'
#' @param object cf_model object
#' @param type One of "ate", "ite", "y0", "y1", "summary", or "samples"
#' @param ... Additional arguments
#' @return Depending on type:
#'   - "ate": numeric scalar (average treatment effect)
#'   - "ite": numeric vector (individual treatment effects)
#'   - "y0": numeric vector (counterfactual outcomes under T=0)
#'   - "y1": numeric vector (counterfactual outcomes under T=1)
#'   - "summary": list with ATE, confidence intervals, and ITE statistics
#'   - "samples": list with y0, y1, ite, and optionally ate_bootstrap samples
#' @export
#' @examples
#' # Fit a model first
#' set.seed(42)
#' demo_data <- data.frame(
#'   X = rnorm(50),
#'   T = rbinom(50, 1, 0.5)
#' )
#' demo_data$Y <- 1.5 * demo_data$T + 0.3 * demo_data$X + rnorm(50, sd = 0.5)
#' fit <- cf_fit(Y ~ T | X, data = demo_data, method = "gformula")
#'
#' # Average treatment effect (scalar)
#' predict(fit, type = "ate")
#'
#' # Individual treatment effects (vector)
#' ite <- predict(fit, type = "ite")
#' length(ite)  # Same as number of observations
#'
#' # Counterfactual outcomes
#' y0 <- predict(fit, type = "y0")  # Predicted Y if T=0
#' y1 <- predict(fit, type = "y1")  # Predicted Y if T=1
predict.cf_model <- function(object,
                             type = c("ate", "ite", "y0", "y1", "summary", "samples"),
                             ...) {
  # Validate input is a valid cfomics_result
 if (!inherits(object, c("cf_model", "cfomics_result"))) {
    rlang::abort(
      message = "object must be a cf_model from cf_fit()",
      class = "cfomics_invalid_input"
    )
  }

  # Check required structure
  if (is.null(object$method) || is.null(object$fit)) {
    rlang::abort(
      message = "Invalid cf_model structure: missing 'method' or 'fit' component",
      class = "cfomics_invalid_result"
    )
  }

  type <- match.arg(type)

  switch(
    object$method,
    "dowhy_gcm" = predict_cf_dowhy_gcm(object, type = type, ...),
    "ganite"    = predict_cf_ganite(object, type = type, ...),
    "drlearner" = predict_cf_drlearner(object, type = type, ...),
    "cavae"     = predict_cf_cavae(object, type = type, ...),
    "grf"       = predict_cf_grf(object, type = type, ...),
    "ipw"       = predict_cf_ipw(object, type = type, ...),
    "gformula"  = predict_cf_gformula(object, type = type, ...),
    "hdml"      = predict_cf_hdml(object, type = type, ...),
    "hdps"      = predict_cf_hdps(object, type = type, ...),
    "bcf"       = predict_cf_bcf(object, type = type, ...),
    stop(sprintf("Unknown method: %s", object$method))
  )
}

#' @rdname predict.cf_model
#' @export
predict.cfomics_result <- function(object,
                                   type = c("ate", "ite", "y0", "y1", "summary", "samples"),
                                   ...) {
  predict.cf_model(object, type = type, ...)
}

predict_cf_dowhy_gcm <- function(object,
                                 type = c("ate", "ite", "y0", "y1", "summary", "samples"),
                                 ...) {
  type <- match.arg(type)
  res  <- object$fit$res

  switch(
    type,
    "ate"     = res$ate,
    "ite"     = res$ite,
    "y0"      = res$y0_hat,
    "y1"      = res$y1_hat,
    "summary" = res$summary,
    "samples" = res$samples
  )
}

#' Get summary statistics from a cf_model
#'
#' @param object cf_model object
#' @param ... Additional arguments (unused)
#' @return A list containing summary statistics including ATE, confidence
#'   intervals (if bootstrap was enabled), and ITE distribution statistics.
#' @export
#' @examples
#' # Fit a model
#' set.seed(123)
#' demo_data <- data.frame(
#'   X = rnorm(50),
#'   T = rbinom(50, 1, 0.5)
#' )
#' demo_data$Y <- 2 * demo_data$T + 0.5 * demo_data$X + rnorm(50)
#' fit <- cf_fit(Y ~ T | X, data = demo_data, method = "gformula")
#'
#' # Print summary
#' summary(fit)
summary.cf_model <- function(object, ...) {
  # Validate input
  if (!inherits(object, c("cf_model", "cfomics_result"))) {
    rlang::abort(
      message = "object must be a cf_model from cf_fit()",
      class = "cfomics_invalid_input"
    )
  }

  if (is.null(object$fit) || is.null(object$fit$res)) {
    rlang::abort(
      message = "Invalid cf_model: missing fit results",
      class = "cfomics_invalid_result"
    )
  }

  res <- object$fit$res

  cat("Causal Inference Model Summary\n")
  cat("==============================\n")
  cat(sprintf("Method: %s\n", object$method))
  cat(sprintf("N samples: %d\n", object$meta$n))
  cat(sprintf("N covariates: %d\n", object$meta$p))
  cat("\n")
  
  cat("Treatment Effect Estimates:\n")
  cat(sprintf("  ATE: %.4f\n", res$ate))
  
  if (!is.null(res$summary)) {
    summary_stats <- res$summary
    
    if (!is.null(summary_stats$ate_ci_lower)) {
      cat(sprintf("  ATE 95%% CI: [%.4f, %.4f]\n", 
                  summary_stats$ate_ci_lower, summary_stats$ate_ci_upper))
    }
    
    if (!is.null(summary_stats$ite_mean)) {
      cat(sprintf("\nITE Distribution:\n"))
      cat(sprintf("  Mean: %.4f\n", summary_stats$ite_mean))
      cat(sprintf("  Std:  %.4f\n", summary_stats$ite_std))
      
      if (!is.null(summary_stats$ite_quantiles)) {
        q <- summary_stats$ite_quantiles
        cat(sprintf("  5th percentile:  %.4f\n", q$q05))
        cat(sprintf("  50th percentile: %.4f\n", q$q50))
        cat(sprintf("  95th percentile: %.4f\n", q$q95))
      }
    }
  }
  
  invisible(res$summary)
}

#' @rdname summary.cf_model
#' @export
summary.cfomics_result <- function(object, ...) {
  summary.cf_model(object, ...)
}

#' Get metadata from a cf_model
#'
#' @param object cf_model object
#' @return A list containing metadata about the model fit, including
#'   input hash, DAG structure, random seed, and processing information.
#' @export
cf_metadata <- function(object) {
  if (!inherits(object, "cf_model")) {
    stop("object must be a cf_model from cf_fit() or cf_dowhy()")
  }
  
  res <- object$fit$res
  
  if (is.null(res$metadata)) {
    message("No metadata available. Set return_metadata = TRUE when fitting the model.")
    return(NULL)
  }
  
  res$metadata
}
