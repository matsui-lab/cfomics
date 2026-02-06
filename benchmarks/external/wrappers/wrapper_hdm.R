# benchmarks/external/wrappers/wrapper_hdm.R
# hdm (high-dimensional metrics) wrapper for benchmark comparison

#' hdm wrapper
#'
#' Uses double selection Lasso via hdm package for high-dimensional data.
#' This method is specifically designed for p >> n settings where standard
#' regression would fail.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param method Method for inference: "double selection" (default) or
#'   "partialling out"
#' @param post Logical, whether to use post-Lasso (default TRUE)
#' @return List with ate, ite, y0_hat, y1_hat, ci_lower, ci_upper
wrapper_hdm <- function(X, T, Y, method = "double selection", post = TRUE) {
  if (!requireNamespace("hdm", quietly = TRUE)) {
    stop("Package 'hdm' required")
  }

  # hdm::rlassoEffect for treatment effect with double selection
  fit <- hdm::rlassoEffect(
    x = X,
    y = Y,
    d = T,
    method = method,
    post = post
  )

  # Extract ATE and standard error from fit object
  ate <- as.numeric(fit$alpha)
  se <- as.numeric(fit$se)

  # Compute confidence interval
  ci_lower <- ate - 1.96 * se
  ci_upper <- ate + 1.96 * se

  # ITE: hdm doesn't provide heterogeneous effects, use constant
  # This is a limitation of the method - it only estimates ATE
  n <- length(Y)
  ite <- rep(ate, n)

  list(
    ate = ate,
    ite = ite,
    y0_hat = rep(NA_real_, n),
    y1_hat = rep(NA_real_, n),
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}
