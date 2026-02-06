# benchmarks/external/wrappers/wrapper_tmle3.R
# tmle3 (tlverse) wrapper for benchmark comparison

#' tmle3 wrapper
#'
#' Uses Targeted Learning via tmle3/tlverse packages. TMLE is a doubly-robust,
#' semiparametric efficient estimator that combines machine learning with
#' efficient influence function-based targeting to estimate causal effects.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @return List with ate, ite, y0_hat, y1_hat, ci_lower, ci_upper
#'
#' @note tmle3 and sl3 are tlverse packages available from GitHub only:
#'   remotes::install_github("tlverse/sl3")
#'   remotes::install_github("tlverse/tmle3")
wrapper_tmle3 <- function(X, T, Y) {
  if (!requireNamespace("tmle3", quietly = TRUE)) {
    stop("Package 'tmle3' required. Install from GitHub: remotes::install_github('tlverse/tmle3')")
  }
  if (!requireNamespace("sl3", quietly = TRUE)) {
    stop("Package 'sl3' required. Install from GitHub: remotes::install_github('tlverse/sl3')")
  }

  # Prepare data
  # Note: tlverse convention uses 'A' for treatment, but we use 'A' internally
  df <- data.frame(Y = Y, A = T, X)
  cov_names <- colnames(X)
  if (is.null(cov_names)) cov_names <- paste0("X", seq_len(ncol(X)))
  colnames(df) <- c("Y", "A", cov_names)

  # Define node list for TMLE
  node_list <- list(
    W = cov_names,
    A = "A",
    Y = "Y"
  )

  # Super Learner stack with parametric and non-parametric learners
  # Note: For production use, consider adding Lrnr_ranger, Lrnr_xgboost
  lrnr_glm <- sl3::Lrnr_glm$new()
  lrnr_mean <- sl3::Lrnr_mean$new()

  # Add ridge regression for regularization
  lrnr_ridge <- tryCatch(
    sl3::Lrnr_glmnet$new(alpha = 0),
    error = function(e) {
      warning("glmnet not available, using GLM-only stack")
      NULL
    }
  )

  learners <- list(lrnr_glm, lrnr_mean)
  if (!is.null(lrnr_ridge)) {
    learners <- c(learners, list(lrnr_ridge))
  }

  sl <- sl3::Lrnr_sl$new(learners = learners)

  # Learner list: one learner for treatment model (A), one for outcome model (Y)
  learner_list <- list(A = sl, Y = sl)

  # Define specification for Average Treatment Effect (ATE)
  ate_spec <- tmle3::tmle_ATE(
    treatment_level = 1,
    control_level = 0
  )

  # Fit TMLE
  tmle_fit <- tmle3::tmle3(
    ate_spec,
    data = df,
    node_list = node_list,
    learner_list = learner_list
  )

  # Extract results from summary
  summ <- tmle_fit$summary
  ate <- summ$psi
  ci_lower <- summ$lower
  ci_upper <- summ$upper

  # ITE: TMLE gives population-level ATE, not individual treatment effects
  # We return constant ITE equal to the ATE for all observations
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
