# benchmarks/external/wrappers/wrapper_weightit.R
# WeightIt wrapper for benchmark comparison

#' WeightIt wrapper
#'
#' Uses inverse probability weighting via WeightIt package.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param method Weighting method (default: "glm")
#' @param estimand Estimand (default: "ATE")
#' @return List with ate, ite, y0_hat, y1_hat, ci_lower, ci_upper
wrapper_weightit <- function(X, T, Y, method = "glm", estimand = "ATE") {
  if (!requireNamespace("WeightIt", quietly = TRUE)) {
    stop("Package 'WeightIt' required")
  }
  if (!requireNamespace("marginaleffects", quietly = TRUE)) {
    stop("Package 'marginaleffects' required")
  }

  # Prepare data
  df <- data.frame(Y = Y, T = T, X)
  cov_names <- colnames(X)
  if (is.null(cov_names)) cov_names <- paste0("X", seq_len(ncol(X)))
  colnames(df) <- c("Y", "T", cov_names)

  # Build formula for propensity score model
  ps_formula <- as.formula(paste("T ~", paste(cov_names, collapse = " + ")))

  # Get weights via WeightIt
  w_out <- WeightIt::weightit(
    ps_formula,
    data = df,
    method = method,
    estimand = estimand
  )

  df$weights <- w_out$weights

  # Fit weighted outcome model
  out_formula <- as.formula(paste("Y ~ T +", paste(cov_names, collapse = " + ")))
  fit <- lm(out_formula, data = df, weights = weights)

  # Get ATE with CI using marginaleffects
  ate_est <- marginaleffects::avg_comparisons(
    fit,
    variables = "T",
    vcov = "HC3",
    newdata = df,
    wts = "weights"
  )

  ate <- ate_est$estimate
  ci_lower <- ate_est$conf.low
  ci_upper <- ate_est$conf.high

  # ITE: predict Y(1) - Y(0) for each unit
  df0 <- df1 <- df
  df0$T <- 0
  df1$T <- 1
  y0_hat <- predict(fit, newdata = df0)
  y1_hat <- predict(fit, newdata = df1)
  ite <- y1_hat - y0_hat

  list(
    ate = ate,
    ite = as.numeric(ite),
    y0_hat = as.numeric(y0_hat),
    y1_hat = as.numeric(y1_hat),
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}
