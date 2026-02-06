# benchmarks/external/wrappers/wrapper_doubleml.R
# DoubleML wrapper for benchmark comparison

#' DoubleML wrapper
#'
#' Uses Double Machine Learning via DoubleML package. This method uses
#' cross-fitting to avoid regularization bias and provides valid inference
#' for treatment effects in partially linear regression models.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param ml_l Learner for outcome model (default: "regr.ranger")
#' @param ml_m Learner for propensity model (default: "classif.ranger")
#' @return List with ate, ite, y0_hat, y1_hat, ci_lower, ci_upper
wrapper_doubleml <- function(X, T, Y, ml_l = "regr.ranger", ml_m = "classif.ranger") {
  if (!requireNamespace("DoubleML", quietly = TRUE)) {
    stop("Package 'DoubleML' required")
  }
  if (!requireNamespace("mlr3", quietly = TRUE)) {
    stop("Package 'mlr3' required")
  }
  if (!requireNamespace("mlr3learners", quietly = TRUE)) {
    stop("Package 'mlr3learners' required")
  }

  # Prepare data
  df <- data.frame(Y = Y, T = T, X)
  cov_names <- colnames(X)
  if (is.null(cov_names)) cov_names <- paste0("X", seq_len(ncol(X)))
  colnames(df) <- c("Y", "T", cov_names)

  # Create DoubleML data object
  dml_data <- DoubleML::DoubleMLData$new(
    df,
    y_col = "Y",
    d_cols = "T",
    x_cols = cov_names
  )

  # Create learners
  # Suppress mlr3 logging to avoid verbose output
  lgr::get_logger("mlr3")$set_threshold("warn")

  ml_l_obj <- mlr3::lrn(ml_l)
  ml_m_obj <- mlr3::lrn(ml_m)

  # Fit PLR (Partially Linear Regression) model
  # Note: ml_l is the learner for the outcome model (previously called ml_g)
  dml_plr <- DoubleML::DoubleMLPLR$new(
    dml_data,
    ml_l = ml_l_obj,
    ml_m = ml_m_obj
  )
  dml_plr$fit()

  # Extract results
  ate <- dml_plr$coef
  se <- dml_plr$se
  ci_lower <- ate - 1.96 * se
  ci_upper <- ate + 1.96 * se

  # ITE: PLR gives constant effect (no heterogeneity estimation)
  # This is a limitation of the PLR model - it only estimates ATE
  n <- length(Y)
  ite <- rep(ate, n)

  list(
    ate = as.numeric(ate),
    ite = as.numeric(ite),
    y0_hat = rep(NA_real_, n),
    y1_hat = rep(NA_real_, n),
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper)
  )
}
