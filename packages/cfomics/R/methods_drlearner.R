#' Fit DRLearner model using econml
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary or multi-level)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param random_state Integer random seed
#' @param model_propensity Model for propensity score (default: "auto")
#' @param model_regression Model for outcome regression (default: "auto")
#' @param cv Integer, number of cross-validation folds (default: 5)
#' @param ... Additional arguments passed to DRLearner
#' @return List with model fit and results
#' @keywords internal
cf_fit_drlearner <- function(X, T, Y,
                             covariate_names = NULL,
                             random_state = 0L,
                             model_propensity = "auto",
                             model_regression = "auto",
                             cv = 5L,
                             ...) {
  # Ensure Python environment is properly configured
  cf_require_python("drlearner")
  
  econml <- reticulate::import("econml")
  np <- reticulate::import("numpy")
  DRLearner <- reticulate::import("econml.dr")$DRLearner
  
  if (!is.null(random_state)) {
    np$random$seed(as.integer(random_state))
  }
  
  X_py <- reticulate::r_to_py(as.matrix(X))
  T_py <- reticulate::r_to_py(as.numeric(T))
  Y_py <- reticulate::r_to_py(as.numeric(Y))
  
  dr <- DRLearner(cv = as.integer(cv))

  # Wrap Python fitting call with user-friendly error handling
  .wrap_python_call(
    dr$fit(Y_py, T_py, X = X_py),
    method_name = "DRLearner"
  )

  n <- nrow(X)
  n_treatments <- length(unique(T))

  # Wrap Python effect estimation with user-friendly error handling
  if (n_treatments == 2) {
    ite_py <- .wrap_python_call(
      dr$effect(X_py, T0 = 0L, T1 = 1L),
      method_name = "DRLearner"
    )
    ite <- reticulate::py_to_r(ite_py)
  } else {
    cme_py <- .wrap_python_call(
      dr$const_marginal_effect(X_py),
      method_name = "DRLearner"
    )
    ite <- reticulate::py_to_r(cme_py)
  }
  
  ate <- mean(ite)
  
  y0_hat <- Y - T * ite
  y1_hat <- Y + (1 - T) * ite
  
  ite_summary <- list(
    mean = mean(ite),
    std = sd(ite),
    quantile_05 = quantile(ite, 0.05),
    quantile_50 = quantile(ite, 0.50),
    quantile_95 = quantile(ite, 0.95)
  )
  
  summary_stats <- list(
    ate = ate,
    ate_ci_lower = NA,
    ate_ci_upper = NA,
    ite_mean = ite_summary$mean,
    ite_std = ite_summary$std,
    ite_quantiles = list(
      q05 = ite_summary$quantile_05,
      q50 = ite_summary$quantile_50,
      q95 = ite_summary$quantile_95
    )
  )
  
  list(
    model = dr,
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
      metadata = NULL
    )
  )
}

#' Predict from DRLearner model
#'
#' @param object cf_model object with method="drlearner"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_drlearner <- function(object, newdata = NULL, type = "ite", ...) {
  if (is.null(newdata)) {
    res <- object$fit$res
    return(switch(
      type,
      "ate" = res$ate,
      "ite" = res$ite,
      "y0" = res$y0_hat,
      "y1" = res$y1_hat,
      "summary" = res$summary,
      "samples" = res$samples
    ))
  }
  
  np <- reticulate::import("numpy")
  X_new <- reticulate::r_to_py(as.matrix(newdata))
  
  dr <- object$fit$model
  ite_py <- dr$effect(X_new, T0 = 0L, T1 = 1L)
  ite <- reticulate::py_to_r(ite_py)
  
  switch(
    type,
    "ate" = mean(ite),
    "ite" = ite,
    stop(sprintf("type '%s' not supported for newdata prediction", type))
  )
}
