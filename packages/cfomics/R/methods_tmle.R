#' Fit Targeted Maximum Likelihood Estimation model
#'
#' Implements a basic TMLE estimator for ATE using glmnet for nuisance
#' parameter estimation. This is a simplified implementation; for full
#' Super Learner integration, consider using the tmle3 package directly.
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param sl_lib Character vector of Super Learner libraries (currently unused,
#'   for future compatibility with tmle3/sl3)
#' @param nfolds Integer, number of CV folds for nuisance estimation
#' @param ... Additional arguments
#' @return List with model fit and results
#' @keywords internal
cf_fit_tmle <- function(X, T, Y,
                        covariate_names = NULL,
                        sl_lib = NULL,
                        nfolds = 10L,
                        ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for method='tmle'. Please install it.")
  }

  X_mat <- as.matrix(X)
  W <- as.numeric(T)
  Y_vec <- as.numeric(Y)
  n <- length(Y_vec)

  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }

  # Step 1: Initial outcome model Q(A, W)
  # Fit outcome model using cross-validation
  Q_cv <- glmnet::cv.glmnet(cbind(X_mat, W), Y_vec, alpha = 0.5, nfolds = nfolds)
  Q_lambda <- Q_cv$lambda.min

  # Get initial predictions
  Q1_init <- as.numeric(predict(Q_cv, newx = cbind(X_mat, 1), s = Q_lambda))
  Q0_init <- as.numeric(predict(Q_cv, newx = cbind(X_mat, 0), s = Q_lambda))
  QA_init <- W * Q1_init + (1 - W) * Q0_init

  # Step 2: Propensity score model g(W)
  g_cv <- glmnet::cv.glmnet(X_mat, W, family = "binomial", alpha = 1, nfolds = nfolds)
  g_lambda <- g_cv$lambda.min
  g_hat <- as.numeric(predict(g_cv, newx = X_mat, s = g_lambda, type = "response"))
  g_hat <- pmax(pmin(g_hat, 0.99), 0.01)  # Bound propensity scores

  # Step 3: Clever covariate H(A, W)
  H1 <- 1 / g_hat
  H0 <- -1 / (1 - g_hat)
  HA <- W * H1 + (1 - W) * H0

  # Step 4: Targeting step - fit epsilon
  # Regress Y - Q_init on H with offset Q_init
  # Using logistic regression for bounded outcomes would be better,
  # but for simplicity we use linear here
  epsilon_model <- stats::lm(Y_vec ~ -1 + offset(QA_init) + HA)
  epsilon <- stats::coef(epsilon_model)["HA"]
  if (is.na(epsilon)) epsilon <- 0

  # Step 5: Update predictions
  Q1_star <- Q1_init + epsilon * H1
  Q0_star <- Q0_init + epsilon * H0

  # ATE estimate
  ate <- mean(Q1_star - Q0_star)

  # Step 6: Influence curve based variance
  IC <- (Q1_star - Q0_star - ate) +
        H1 * W * (Y_vec - Q1_star) +
        H0 * (1 - W) * (Y_vec - Q0_star)
  ate_se <- sqrt(stats::var(IC) / n)

  ate_ci_lower <- ate - 1.96 * ate_se
  ate_ci_upper <- ate + 1.96 * ate_se

  # ITE approximation (model-based, not true CATE)
  ite <- Q1_star - Q0_star

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

  list(
    model = list(
      Q_model = Q_cv,
      g_model = g_cv,
      epsilon = epsilon
    ),
    res = list(
      ite = ite,
      ate = ate,
      y0_hat = Q0_star,
      y1_hat = Q1_star,
      summary = summary_stats,
      samples = list(
        y0 = Q0_star,
        y1 = Q1_star,
        ite = ite
      ),
      propensity = g_hat,
      influence_curve = IC,
      metadata = list(
        epsilon = epsilon,
        sl_lib = sl_lib
      )
    )
  )
}

#' Predict from TMLE model
#'
#' @param object cf_model object with method="tmle"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1", "summary", "samples"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_tmle <- function(object, newdata = NULL, type = "ite", ...) {
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
  Q_model <- object$fit$model$Q_model
  g_model <- object$fit$model$g_model
  epsilon <- object$fit$model$epsilon

  X_new <- as.matrix(newdata)

  Q1_init <- as.numeric(predict(Q_model, newx = cbind(X_new, 1), s = "lambda.min"))
  Q0_init <- as.numeric(predict(Q_model, newx = cbind(X_new, 0), s = "lambda.min"))

  g_hat <- as.numeric(predict(g_model, newx = X_new, s = "lambda.min", type = "response"))
  g_hat <- pmax(pmin(g_hat, 0.99), 0.01)

  H1 <- 1 / g_hat
  H0 <- -1 / (1 - g_hat)

  Q1_star <- Q1_init + epsilon * H1
  Q0_star <- Q0_init + epsilon * H0
  ite <- Q1_star - Q0_star

  switch(
    type,
    "ate" = mean(ite),
    "ite" = ite,
    "y0" = Q0_star,
    "y1" = Q1_star,
    stop(sprintf("type '%s' not supported", type))
  )
}
