#' Fit Generalized Random Forest (GRF) causal forest model
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param random_state Integer random seed
#' @param num.trees Integer, number of trees (default: 2000)
#' @param tune.parameters Character, which parameters to tune (default: "all")
#' @param estimate.variance Logical, whether to estimate variance (default: TRUE)
#' @param ... Additional arguments passed to causal_forest
#' @return List with model fit and results
#' @keywords internal
cf_fit_grf <- function(X, T, Y,
                       covariate_names = NULL,
                       random_state = 0L,
                       num.trees = 2000L,
                       tune.parameters = "all",
                       estimate.variance = TRUE,
                       ...) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' is required for method='grf'. Please install it.")
  }
  
  if (!is.null(random_state)) {
    set.seed(random_state)
  }
  
  X_mat <- as.matrix(X)
  W <- as.numeric(T)
  Y_vec <- as.numeric(Y)
  
  forest.W <- grf::regression_forest(X_mat, W, tune.parameters = tune.parameters)
  W.hat <- predict(forest.W)$predictions
  
  forest.Y <- grf::regression_forest(X_mat, Y_vec, tune.parameters = tune.parameters)
  Y.hat <- predict(forest.Y)$predictions
  
  tau.forest <- grf::causal_forest(
    X_mat, Y_vec, W,
    W.hat = W.hat,
    Y.hat = Y.hat,
    num.trees = as.integer(num.trees),
    tune.parameters = tune.parameters,
    ...
  )
  
  tau.pred <- predict(tau.forest, estimate.variance = estimate.variance)
  ite <- tau.pred$predictions
  
  ate <- grf::average_treatment_effect(tau.forest)
  ate_estimate <- ate["estimate"]
  ate_se <- ate["std.err"]
  
  y0_hat <- Y_vec - W * ite
  y1_hat <- Y_vec + (1 - W) * ite
  
  ite_summary <- list(
    mean = mean(ite),
    std = sd(ite),
    quantile_05 = quantile(ite, 0.05),
    quantile_50 = quantile(ite, 0.50),
    quantile_95 = quantile(ite, 0.95)
  )
  
  summary_stats <- list(
    ate = ate_estimate,
    ate_ci_lower = ate_estimate - 1.96 * ate_se,
    ate_ci_upper = ate_estimate + 1.96 * ate_se,
    ite_mean = ite_summary$mean,
    ite_std = ite_summary$std,
    ite_quantiles = list(
      q05 = ite_summary$quantile_05,
      q50 = ite_summary$quantile_50,
      q95 = ite_summary$quantile_95
    )
  )
  
  variance <- if (estimate.variance && "variance.estimates" %in% names(tau.pred)) {
    tau.pred$variance.estimates
  } else {
    NULL
  }
  
  list(
    model = tau.forest,
    res = list(
      ite = ite,
      ate = ate_estimate,
      y0_hat = y0_hat,
      y1_hat = y1_hat,
      summary = summary_stats,
      samples = list(
        y0 = y0_hat,
        y1 = y1_hat,
        ite = ite
      ),
      variance = variance,
      metadata = NULL
    )
  )
}

#' Fit Inverse Probability Weighting (IPW) model
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param stabilized Logical, whether to use stabilized weights (default: TRUE)
#' @param ... Additional arguments passed to ipwpoint
#' @return List with model fit and results
#' @keywords internal
cf_fit_ipw <- function(X, T, Y,
                       covariate_names = NULL,
                       stabilized = TRUE,
                       ...) {
  if (!requireNamespace("ipw", quietly = TRUE)) {
    stop("Package 'ipw' is required for method='ipw'. Please install it.")
  }
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package 'survey' is required for method='ipw'. Please install it.")
  }
  
  X_mat <- as.matrix(X)
  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }
  
  df <- as.data.frame(X_mat)
  colnames(df) <- covariate_names
  df$A <- as.numeric(T)
  df$Y <- as.numeric(Y)
  
  denom_formula_str <- paste("~", paste(covariate_names, collapse = " + "))
  
  ipw_fit <- eval(substitute(
    ipw::ipwpoint(
      exposure = A,
      family = "binomial",
      link = "logit",
      numerator = ~ 1,
      denominator = DENOM,
      data = df
    ),
    list(DENOM = as.formula(denom_formula_str))
  ))
  
  df$sw <- ipw_fit$ipw.weights
  
  design <- survey::svydesign(~1, weights = ~sw, data = df)
  msm <- survey::svyglm(Y ~ A, design = design)
  
  coefs <- coef(msm)
  ci <- confint(msm)
  
  ate <- coefs["A"]
  ate_ci_lower <- ci["A", 1]
  ate_ci_upper <- ci["A", 2]
  
  y0_hat <- df$Y - df$A * ate
  y1_hat <- df$Y + (1 - df$A) * ate


  ite <- rep(NA_real_, nrow(df))
  rlang::inform(
    c("IPW estimates population-level ATE only.",
      "i" = "Individual treatment effects (ITE) are not available with this method.",
      "i" = "Consider using 'grf' or 'drlearner' for individual-level estimates."),
    class = "cfomics_ipw_no_ite"
  )

  summary_stats <- list(
    ate = ate,
    ate_ci_lower = ate_ci_lower,
    ate_ci_upper = ate_ci_upper,
    ite_mean = NA_real_,
    ite_std = NA_real_,
    ite_quantiles = list(
      q05 = NA_real_,
      q50 = NA_real_,
      q95 = NA_real_
    )
  )
  
  list(
    model = list(ipw = ipw_fit, msm = msm),
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
      weights = df$sw,
      metadata = NULL
    )
  )
}

#' Fit G-formula model
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param nsimul Integer, number of simulations (default: 1000)
#' @param ... Additional arguments
#' @return List with model fit and results
#' @keywords internal
cf_fit_gformula <- function(X, T, Y,
                            covariate_names = NULL,
                            nsimul = 1000L,
                            ...) {
  X_mat <- as.matrix(X)
  if (is.null(covariate_names)) {
    covariate_names <- paste0("X", seq_len(ncol(X_mat)))
  }
  
  df <- as.data.frame(X_mat)
  colnames(df) <- covariate_names
  df$A <- as.numeric(T)
  df$Y <- as.numeric(Y)
  
  outcome_formula <- as.formula(paste("Y ~ A +", paste(covariate_names, collapse = " + ")))
  outcome_model <- lm(outcome_formula, data = df)
  
  df0 <- df
  df0$A <- 0
  y0_hat <- predict(outcome_model, newdata = df0)
  
  df1 <- df
  df1$A <- 1
  y1_hat <- predict(outcome_model, newdata = df1)
  
  ite <- y1_hat - y0_hat
  ate <- mean(ite)
  
  boot_ates <- numeric(nsimul)
  n <- nrow(df)
  
  for (b in seq_len(nsimul)) {
    idx <- sample(n, n, replace = TRUE)
    df_boot <- df[idx, ]
    
    model_boot <- lm(outcome_formula, data = df_boot)
    
    df0_boot <- df_boot
    df0_boot$A <- 0
    y0_boot <- predict(model_boot, newdata = df0_boot)
    
    df1_boot <- df_boot
    df1_boot$A <- 1
    y1_boot <- predict(model_boot, newdata = df1_boot)
    
    boot_ates[b] <- mean(y1_boot - y0_boot)
  }
  
  ate_ci_lower <- quantile(boot_ates, 0.025)
  ate_ci_upper <- quantile(boot_ates, 0.975)
  
  ite_summary <- list(
    mean = mean(ite),
    std = sd(ite),
    quantile_05 = quantile(ite, 0.05),
    quantile_50 = quantile(ite, 0.50),
    quantile_95 = quantile(ite, 0.95)
  )
  
  summary_stats <- list(
    ate = ate,
    ate_ci_lower = ate_ci_lower,
    ate_ci_upper = ate_ci_upper,
    ite_mean = ite_summary$mean,
    ite_std = ite_summary$std,
    ite_quantiles = list(
      q05 = ite_summary$quantile_05,
      q50 = ite_summary$quantile_50,
      q95 = ite_summary$quantile_95
    )
  )
  
  list(
    model = outcome_model,
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
        ate_bootstrap = boot_ates
      ),
      metadata = NULL
    )
  )
}

#' Predict from GRF model
#'
#' @param object cf_model object with method="grf"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param estimate.variance Logical, whether to estimate variance
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_grf <- function(object, newdata = NULL, type = "ite",
                           estimate.variance = TRUE, ...) {
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
  
  tau.forest <- object$fit$model
  X_new <- as.matrix(newdata)
  
  tau.pred <- predict(tau.forest, X_new, estimate.variance = estimate.variance)
  ite <- tau.pred$predictions
  
  switch(
    type,
    "ate" = mean(ite),
    "ite" = ite,
    stop(sprintf("type '%s' not supported for newdata prediction", type))
  )
}

#' Predict from IPW model
#'
#' @param object cf_model object with method="ipw"
#' @param newdata Not supported for IPW
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_ipw <- function(object, newdata = NULL, type = "ite", ...) {
  if (!is.null(newdata)) {
    stop("IPW does not support prediction on new data")
  }
  
  res <- object$fit$res
  switch(
    type,
    "ate" = res$ate,
    "ite" = res$ite,
    "y0" = res$y0_hat,
    "y1" = res$y1_hat,
    "summary" = res$summary,
    "samples" = res$samples
  )
}

#' Predict from G-formula model
#'
#' @param object cf_model object with method="gformula"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_gformula <- function(object, newdata = NULL, type = "ite", ...) {
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
  
  outcome_model <- object$fit$model
  covariate_names <- names(coef(outcome_model))[-c(1, 2)]
  
  df_new <- as.data.frame(as.matrix(newdata))
  colnames(df_new) <- covariate_names
  
  df0 <- df_new
  df0$A <- 0
  y0_hat <- predict(outcome_model, newdata = df0)
  
  df1 <- df_new
  df1$A <- 1
  y1_hat <- predict(outcome_model, newdata = df1)
  
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
