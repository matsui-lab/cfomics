#' Fit High-Dimensional Machine Learning (Debiased Lasso) model
#'
#' Uses regularized regression with AIPW/doubly-robust estimation for
#' causal effect estimation in high-dimensional settings.
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param penalty Character, penalty type: "lasso", "ridge", "elastic_net"
#' @param lambda Character or numeric, lambda selection: "cv", "1se", or numeric
#' @param alpha Numeric, elastic net mixing parameter (1=lasso, 0=ridge)
#' @param nfolds Integer, number of CV folds
#' @param ... Additional arguments
#' @return List with model fit and results
#' @keywords internal
cf_fit_hdml <- function(X, T, Y,
                        covariate_names = NULL,
                        penalty = "lasso",
                        lambda = "cv",
                        alpha = NULL,
                        nfolds = 10L,
                        ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for method='hdml'. Please install it.")
  }

  X_mat <- as.matrix(X)
  W <- as.numeric(T)

  Y_vec <- as.numeric(Y)
  n <- length(Y_vec)

  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }

  # Set alpha based on penalty type
  if (is.null(alpha)) {
    alpha <- switch(penalty,
      "lasso" = 1,
      "ridge" = 0,
      "elastic_net" = 0.5,
      1
    )
  }

  # Step 1: Estimate propensity score with regularized logistic regression
  ps_cv <- glmnet::cv.glmnet(X_mat, W, family = "binomial",
                             alpha = alpha, nfolds = nfolds)
  ps_lambda <- if (lambda == "1se") ps_cv$lambda.1se else ps_cv$lambda.min
  ps_hat <- as.numeric(predict(ps_cv, newx = X_mat, s = ps_lambda, type = "response"))
  ps_hat <- pmax(pmin(ps_hat, 0.99), 0.01)  # Trim extreme values

  # Step 2: Estimate outcome model with regularized regression
  outcome_cv <- glmnet::cv.glmnet(cbind(X_mat, W), Y_vec,
                                   alpha = alpha, nfolds = nfolds)
  outcome_lambda <- if (lambda == "1se") outcome_cv$lambda.1se else outcome_cv$lambda.min

  # Step 3: Predict counterfactual outcomes
  X_with_T0 <- cbind(X_mat, 0)
  X_with_T1 <- cbind(X_mat, 1)
  y0_hat <- as.numeric(predict(outcome_cv, newx = X_with_T0, s = outcome_lambda))
  y1_hat <- as.numeric(predict(outcome_cv, newx = X_with_T1, s = outcome_lambda))

  # Step 4: Compute AIPW/doubly-robust ATE estimator
  ite_dr <- y1_hat - y0_hat +
    W * (Y_vec - y1_hat) / ps_hat -
    (1 - W) * (Y_vec - y0_hat) / (1 - ps_hat)

  ate <- mean(ite_dr)

  # Step 5: Estimate variance using influence function
  influence_fn <- ite_dr - ate
  ate_se <- sqrt(var(influence_fn) / n)

  ate_ci_lower <- ate - 1.96 * ate_se
  ate_ci_upper <- ate + 1.96 * ate_se

  # ITE is approximated by y1_hat - y0_hat (model-based)
  ite <- y1_hat - y0_hat

  summary_stats <- list(
    ate = ate,
    ate_ci_lower = ate_ci_lower,
    ate_ci_upper = ate_ci_upper,
    ite_mean = mean(ite),
    ite_std = sd(ite),
    ite_quantiles = list(
      q05 = unname(quantile(ite, 0.05)),
      q50 = unname(quantile(ite, 0.50)),
      q95 = unname(quantile(ite, 0.95))
    )
  )

  list(
    model = list(
      ps_model = ps_cv,
      outcome_model = outcome_cv,
      ps_lambda = ps_lambda,
      outcome_lambda = outcome_lambda
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
      metadata = list(
        penalty = penalty,
        alpha = alpha,
        lambda = lambda
      )
    )
  )
}

#' Predict from HDML model
#'
#' @param object cf_model object with method="hdml"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_hdml <- function(object, newdata = NULL, type = "ite", ...) {
  if (is.null(newdata)) {
    res <- object$fit$res
    return(switch(
      type,
      "ate" = res$ate,
      "ite" = res$ite,
      "y0" = res$y0_hat,
      "y1" = res$y1_hat,
      "summary" = res$summary,
      "samples" = res$samples,
      stop(sprintf("Unknown type: %s", type))
    ))
  }

  # Prediction on new data
  outcome_model <- object$fit$model$outcome_model
  outcome_lambda <- object$fit$model$outcome_lambda

  X_new <- as.matrix(newdata)
  X_with_T0 <- cbind(X_new, 0)
  X_with_T1 <- cbind(X_new, 1)

  y0_hat <- as.numeric(predict(outcome_model, newx = X_with_T0, s = outcome_lambda))
  y1_hat <- as.numeric(predict(outcome_model, newx = X_with_T1, s = outcome_lambda))
  ite <- y1_hat - y0_hat

  switch(
    type,
    "ate" = mean(ite),
    "ite" = ite,
    "y0" = y0_hat,
    "y1" = y1_hat,
    stop(sprintf("type '%s' not supported for newdata prediction", type))
  )
}
