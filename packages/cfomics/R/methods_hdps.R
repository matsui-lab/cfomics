#' Fit High-Dimensional Propensity Score model
#'
#' Uses LASSO-regularized propensity score estimation with IPW for
#' causal effect estimation in high-dimensional settings.
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param lambda Character or numeric, lambda selection: "cv" or "1se"
#' @param trim Numeric vector of length 2, propensity score trimming bounds
#' @param nfolds Integer, number of CV folds
#' @param ... Additional arguments
#' @return List with model fit and results
#' @keywords internal
cf_fit_hdps <- function(X, T, Y,
                        covariate_names = NULL,
                        lambda = "cv",
                        trim = c(0.01, 0.99),
                        nfolds = 10L,
                        ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for method='hdps'. Please install it.")
  }

  X_mat <- as.matrix(X)
  W <- as.numeric(T)
  Y_vec <- as.numeric(Y)
  n <- length(Y_vec)

  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }

  # Step 1: Estimate propensity score with LASSO logistic regression
  ps_cv <- glmnet::cv.glmnet(X_mat, W, family = "binomial",
                             alpha = 1, nfolds = nfolds)
  ps_lambda <- if (lambda == "1se") ps_cv$lambda.1se else ps_cv$lambda.min
  ps_hat <- as.numeric(predict(ps_cv, newx = X_mat, s = ps_lambda, type = "response"))

  # Trim propensity scores
  ps_hat <- pmax(pmin(ps_hat, trim[2]), trim[1])

  # Step 2: Compute IPW weights
  weights <- W / ps_hat + (1 - W) / (1 - ps_hat)

  # Step 3: Weighted mean difference (Horvitz-Thompson estimator)
  y1_weighted <- sum(W * Y_vec / ps_hat) / sum(W / ps_hat)
  y0_weighted <- sum((1 - W) * Y_vec / (1 - ps_hat)) / sum((1 - W) / (1 - ps_hat))
  ate <- y1_weighted - y0_weighted

  # Step 4: Variance estimation via influence function
  mu1 <- W * Y_vec / ps_hat
  mu0 <- (1 - W) * Y_vec / (1 - ps_hat)
  influence_fn <- mu1 - mu0 - ate
  ate_se <- sqrt(var(influence_fn) / n)

  ate_ci_lower <- ate - 1.96 * ate_se
  ate_ci_upper <- ate + 1.96 * ate_se

  # IPW doesn't estimate ITE, use ATE for all
  ite <- rep(ate, n)
  y0_hat <- Y_vec - W * ate
  y1_hat <- Y_vec + (1 - W) * ate

  summary_stats <- list(
    ate = ate,
    ate_ci_lower = ate_ci_lower,
    ate_ci_upper = ate_ci_upper,
    ite_mean = ate,
    ite_std = 0,
    ite_quantiles = list(
      q05 = ate,
      q50 = ate,
      q95 = ate
    )
  )

  list(
    model = list(
      ps_model = ps_cv,
      ps_lambda = ps_lambda
    ),
    res = list(
      ite = ite,
      ate = ate,
      y0_hat = y0_hat,
      y1_hat = y1_hat,
      summary = summary_stats,
      samples = list(
        y0 = y0_hat,
        y1 = y1_hat,
        ite = ite
      ),
      propensity = ps_hat,
      weights = weights,
      metadata = list(
        lambda = lambda,
        trim = trim,
        n_trimmed = sum(ps_hat == trim[1] | ps_hat == trim[2])
      )
    )
  )
}

#' Predict from HDPS model
#'
#' @param object cf_model object with method="hdps"
#' @param newdata Not supported for HDPS
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_hdps <- function(object, newdata = NULL, type = "ite", ...) {
  if (!is.null(newdata)) {
    stop("HDPS does not support prediction on new data")
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
