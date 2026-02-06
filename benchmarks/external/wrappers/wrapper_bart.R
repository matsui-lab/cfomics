# benchmarks/external/wrappers/wrapper_bart.R
# BART (Bayesian Additive Regression Trees) wrapper for benchmark comparison

#' BART wrapper
#'
#' Uses Bayesian Additive Regression Trees via dbarts package. BART is a
#' Bayesian "sum of trees" model that provides full posterior inference
#' for counterfactual predictions, enabling heterogeneous treatment effect
#' estimation with uncertainty quantification.
#'
#' Unlike methods that only estimate ATE (e.g., hdm, tmle3), BART provides:
#' - Real counterfactual predictions (y0_hat, y1_hat are NOT NA)
#' - Heterogeneous ITE (not constant)
#' - Full posterior-based uncertainty quantification
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param n_trees Number of trees (default: 200)
#' @param n_burn Burn-in iterations (default: 250)
#' @param n_iter Post-burn iterations (default: 1000)
#' @return List with ate, ite, y0_hat, y1_hat, ci_lower, ci_upper
wrapper_bart <- function(X, T, Y, n_trees = 200L, n_burn = 250L, n_iter = 1000L) {
  if (!requireNamespace("dbarts", quietly = TRUE)) {
    stop("Package 'dbarts' required")
  }

  n <- length(Y)

  # Combine X and T for BART
  X_with_T <- cbind(T = T, X)

  # Fit BART
  # keeptrees = TRUE is required for predict() on new data
  bart_fit <- dbarts::bart(
    x.train = X_with_T,
    y.train = Y,
    ntree = n_trees,
    nskip = n_burn,
    ndpost = n_iter,
    keeptrees = TRUE,
    verbose = FALSE
  )

  # Counterfactual predictions
  X_t0 <- cbind(T = 0, X)
  X_t1 <- cbind(T = 1, X)

  # Posterior predictions (matrix: n_iter x n)
  pred_t0 <- predict(bart_fit, newdata = X_t0)
  pred_t1 <- predict(bart_fit, newdata = X_t1)

  # Point estimates (posterior mean)
  y0_hat <- colMeans(pred_t0)
  y1_hat <- colMeans(pred_t1)
  ite <- y1_hat - y0_hat
  ate <- mean(ite)

  # CI from posterior
  ite_samples <- pred_t1 - pred_t0
  ate_samples <- rowMeans(ite_samples)
  ci_lower <- quantile(ate_samples, 0.025)
  ci_upper <- quantile(ate_samples, 0.975)

  list(
    ate = ate,
    ite = as.numeric(ite),
    y0_hat = as.numeric(y0_hat),
    y1_hat = as.numeric(y1_hat),
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper)
  )
}
