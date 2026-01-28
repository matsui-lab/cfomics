#' Fit Bayesian Causal Forests model
#'
#' Uses the bcf package to estimate heterogeneous treatment effects
#' via Bayesian Additive Regression Trees (BART).
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param n_burn Integer, number of burn-in iterations
#' @param n_iter Integer, number of MCMC iterations after burn-in
#' @param n_chains Integer, number of chains (currently only 1 supported)
#' @param pihat Optional propensity scores; if NULL, estimated internally
#' @param ... Additional arguments passed to bcf()
#' @return List with model fit and results
#' @keywords internal
cf_fit_bcf <- function(X, T, Y,
                       covariate_names = NULL,
                       n_burn = 1000L,
                       n_iter = 2000L,
                       n_chains = 1L,
                       pihat = NULL,
                       ...) {
  if (!requireNamespace("bcf", quietly = TRUE)) {
    stop("Package 'bcf' is required for method='bcf'. Please install it.")
  }

  X_mat <- as.matrix(X)
  W <- as.numeric(T)
  Y_vec <- as.numeric(Y)
  n <- length(Y_vec)

  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }

  # Estimate propensity score if not provided
  if (is.null(pihat)) {
    if (requireNamespace("glmnet", quietly = TRUE)) {
      ps_cv <- glmnet::cv.glmnet(X_mat, W, family = "binomial", alpha = 1)
      pihat <- as.numeric(predict(ps_cv, newx = X_mat, s = "lambda.min",
                                   type = "response"))
    } else {
      # Simple logistic regression fallback
      df_ps <- data.frame(W = W, X_mat)
      ps_model <- stats::glm(W ~ ., data = df_ps, family = stats::binomial)
      pihat <- stats::predict(ps_model, type = "response")
    }
    pihat <- pmax(pmin(pihat, 0.99), 0.01)
  }

  # Fit BCF
  bcf_fit <- bcf::bcf(
    y = Y_vec,
    z = W,
    x_control = X_mat,
    x_moderate = X_mat,
    pihat = pihat,
    nburn = as.integer(n_burn),
    nsim = as.integer(n_iter),
    ...
  )

  # Extract treatment effects
  # bcf returns tau samples of dimension (n_iter x n)
  tau_samples <- bcf_fit$tau
  ite <- colMeans(tau_samples)
  ate <- mean(ite)

  # Credible intervals from posterior
  ate_samples <- rowMeans(tau_samples)
  ate_ci_lower <- unname(stats::quantile(ate_samples, 0.025))
  ate_ci_upper <- unname(stats::quantile(ate_samples, 0.975))

  # Posterior predictive for counterfactuals
  # yhat contains the fitted values, we need to extract mu (prognostic)
  # In bcf, yhat = mu + z * tau, so mu = yhat - z * tau
  yhat_samples <- bcf_fit$yhat
  mu_samples <- yhat_samples - outer(rep(1, nrow(tau_samples)), W) * tau_samples
  y0_hat <- colMeans(mu_samples)
  y1_hat <- y0_hat + ite

  summary_stats <- list(
    ate = ate,
    ate_ci_lower = ate_ci_lower,
    ate_ci_upper = ate_ci_upper,
    ite_mean = mean(ite),
    ite_std = stats::sd(ite),
    ite_quantiles = list(
      q05 = unname(stats::quantile(ite, 0.05)),
      q50 = unname(stats::quantile(ite, 0.50)),
      q95 = unname(stats::quantile(ite, 0.95))
    )
  )

  # ITE credible intervals
  ite_ci_lower <- apply(tau_samples, 2, stats::quantile, 0.025)
  ite_ci_upper <- apply(tau_samples, 2, stats::quantile, 0.975)

  list(
    model = bcf_fit,
    res = list(
      ite = ite,
      ate = ate,
      y0_hat = y0_hat,
      y1_hat = y1_hat,
      summary = summary_stats,
      samples = list(
        y0 = y0_hat,
        y1 = y1_hat,
        ite = ite,
        tau_posterior = tau_samples,
        ate_posterior = ate_samples
      ),
      credible_intervals = list(
        ite_lower = ite_ci_lower,
        ite_upper = ite_ci_upper
      ),
      propensity = pihat,
      metadata = list(
        n_burn = n_burn,
        n_iter = n_iter
      )
    )
  )
}

#' Predict from BCF model
#'
#' @param object cf_model object with method="bcf"
#' @param newdata Optional new data for prediction (not yet supported)
#' @param type Type of prediction: "ate", "ite", "y0", "y1", "summary", "samples"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_bcf <- function(object, newdata = NULL, type = "ite", ...) {
  if (!is.null(newdata)) {
    stop("BCF prediction on new data not yet implemented")
  }

  res <- object$fit$res
  switch(
    type,
    "ate" = res$ate,
    "ite" = res$ite,
    "y0" = res$y0_hat,
    "y1" = res$y1_hat,
    "summary" = res$summary,
    "samples" = res$samples,
    stop(sprintf("Unknown type: %s", type))
  )
}
