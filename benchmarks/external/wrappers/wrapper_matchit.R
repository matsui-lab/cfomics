# benchmarks/external/wrappers/wrapper_matchit.R
# MatchIt wrapper for benchmark comparison

#' MatchIt wrapper
#'
#' Uses propensity score matching via MatchIt package.
#' By default uses "nearest" matching with ATT estimand (most widely available).
#' For ATE estimation with "full" method, the 'optmatch' package is required.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param method Matching method (default: "nearest").
#' @param estimand Estimand (default: "ATT"). For "nearest", use "ATT" or "ATC".
#'   For "full" or "subclass", can use "ATE".
#' @return List with ate, ite, y0_hat, y1_hat, ci_lower, ci_upper
wrapper_matchit <- function(X, T, Y, method = "nearest", estimand = "ATT") {
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' required")
  }
  if (!requireNamespace("marginaleffects", quietly = TRUE)) {
    stop("Package 'marginaleffects' required")
  }

  # Prepare data
  df <- data.frame(Y = Y, T = T, X)
  cov_names <- colnames(X)
  if (is.null(cov_names)) cov_names <- paste0("X", seq_len(ncol(X)))
  colnames(df) <- c("Y", "T", cov_names)

  # Build formula
  ps_formula <- as.formula(paste("T ~", paste(cov_names, collapse = " + ")))

  # Determine estimand based on method
  # For methods that don't support ATE, use ATT
  methods_supporting_ate <- c("full", "subclass", "cem", "exact", "cardinality")
  if (estimand == "ATE" && !(method %in% methods_supporting_ate)) {
    estimand <- "ATT"  # Fall back to ATT for nearest, genetic, optimal
  }

  # Check for optmatch package if using full matching
  if (method == "full" && !requireNamespace("optmatch", quietly = TRUE)) {
    warning("'optmatch' package required for full matching; falling back to 'nearest'")
    method <- "nearest"
    if (estimand == "ATE") estimand <- "ATT"
  }

  # Run matching
  m_out <- MatchIt::matchit(
    ps_formula,
    data = df,
    method = method,
    estimand = estimand
  )

  # Get matched data
  m_data <- MatchIt::match.data(m_out)

  # Fit outcome model on matched data
  out_formula <- as.formula(paste("Y ~ T +", paste(cov_names, collapse = " + ")))
  fit <- lm(out_formula, data = m_data, weights = weights)

  # Get ATE with CI using marginaleffects
  # Use default vcov (TRUE) for robust standard errors
  ate_est <- marginaleffects::avg_comparisons(
    fit,
    variables = "T",
    vcov = TRUE,
    newdata = m_data,
    wts = "weights"
  )

  ate <- ate_est$estimate
  ci_lower <- ate_est$conf.low
  ci_upper <- ate_est$conf.high

  # ITE: predict Y(1) - Y(0) for each unit
  # Use original data for predictions
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
