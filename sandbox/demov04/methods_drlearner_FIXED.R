#' Fit DRLearner model using econml (FIXED VERSION)
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary or multi-level)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param random_state Integer random seed
#' @param model_propensity Model for propensity score (default: "auto")
#' @param model_regression Model for outcome regression (default: "auto")
#' @param model_final Model for final CATE estimation (default: "auto")
#' @param cv Integer, number of cross-validation folds (default: 5)
#' @param min_propensity Minimum propensity score for clipping (default: 0.05)
#' @param standardize_y Whether to standardize Y before fitting (default: TRUE)
#' @param ... Additional arguments passed to DRLearner
#' @return List with model fit and results
#' @keywords internal
cf_fit_drlearner <- function(X, T, Y,
                             covariate_names = NULL,
                             random_state = 0L,
                             model_propensity = "auto",
                             model_regression = "auto",
                             model_final = "auto",
                             cv = 5L,
                             min_propensity = 0.05,
                             standardize_y = TRUE,
                             ...) {
  # Ensure Python environment is properly configured
  cf_require_python("drlearner")

  econml <- reticulate::import("econml")
  sklearn <- reticulate::import("sklearn")
  np <- reticulate::import("numpy")
  DRLearner <- reticulate::import("econml.dr")$DRLearner

  # ============================================================
  # FIX #4: Y の標準化
  # ============================================================
  Y_orig <- Y
  if (standardize_y) {
    Y_mean <- mean(Y)
    Y_sd <- sd(Y)
    Y <- (Y - Y_mean) / Y_sd
  } else {
    Y_mean <- 0
    Y_sd <- 1
  }

  # ============================================================
  # FIX #5: 乱数シード制御の改善
  # ============================================================
  if (!is.null(random_state)) {
    np$random$seed(as.integer(random_state))
    # Python の random モジュールもシード設定
    random <- reticulate::import("random")
    random$seed(as.integer(random_state))
  }

  # ============================================================
  # FIX #1: デフォルトモデルの設定と引数の渡し
  # ============================================================
  # Propensity model: デフォルトは LogisticRegressionCV（正則化付き）
  if (identical(model_propensity, "auto")) {
    LogisticRegressionCV <- sklearn$linear_model$LogisticRegressionCV
    model_prop_py <- LogisticRegressionCV(
      cv = as.integer(cv),
      random_state = as.integer(random_state),
      max_iter = 1000L
    )
  } else if (is.null(model_propensity)) {
    model_prop_py <- NULL  # EconML default
  } else {
    model_prop_py <- model_propensity
  }

  # Regression model: デフォルトは Ridge（正則化付き）
  if (identical(model_regression, "auto")) {
    RidgeCV <- sklearn$linear_model$RidgeCV
    model_reg_py <- RidgeCV(cv = as.integer(cv))
  } else if (is.null(model_regression)) {
    model_reg_py <- NULL  # EconML default
  } else {
    model_reg_py <- model_regression
  }

  # Final model: デフォルトは Lasso（外れ値にロバスト）
  if (identical(model_final, "auto")) {
    LassoCV <- sklearn$linear_model$LassoCV
    model_final_py <- LassoCV(
      cv = as.integer(cv),
      random_state = as.integer(random_state)
    )
  } else if (is.null(model_final)) {
    model_final_py <- NULL  # EconML default
  } else {
    model_final_py <- model_final
  }

  # ============================================================
  # FIX #1 & #2 & #5: DRLearner コンストラクタに引数を渡す
  # ============================================================
  X_py <- reticulate::r_to_py(as.matrix(X))
  T_py <- reticulate::r_to_py(as.numeric(T))
  Y_py <- reticulate::r_to_py(as.numeric(Y))

  dr <- DRLearner(
    model_propensity = model_prop_py,
    model_regression = model_reg_py,
    model_final = model_final_py,
    cv = as.integer(cv),
    min_propensity = min_propensity,  # FIX #2: 傾向スコアクリッピング
    random_state = as.integer(random_state)  # FIX #5: EconMLに直接渡す
  )

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

  # ============================================================
  # FIX #4: 標準化を元に戻す
  # ============================================================
  if (standardize_y) {
    ite <- ite * Y_sd  # ITE は差分なのでスケールのみ戻す
  }

  ate <- mean(ite)

  y0_hat <- Y_orig - T * ite
  y1_hat <- Y_orig + (1 - T) * ite

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
      metadata = list(
        standardize_y = standardize_y,
        Y_mean = Y_mean,
        Y_sd = Y_sd,
        min_propensity = min_propensity
      )
    )
  )
}

#' Predict from DRLearner model (FIXED VERSION)
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

  # 標準化を元に戻す
  metadata <- object$fit$res$metadata
  if (!is.null(metadata) && metadata$standardize_y) {
    ite <- ite * metadata$Y_sd
  }

  switch(
    type,
    "ate" = mean(ite),
    "ite" = ite,
    stop(sprintf("type '%s' not supported for newdata prediction", type))
  )
}
