#' Generate simulation data for benchmarking causal inference methods
#'
#' This function generates synthetic data with known true treatment effects
#' for benchmarking causal inference methods. The data includes covariates,
#' binary treatment, and outcome with known ATE and ITE.
#'
#' @param scenario Character, the data generating scenario. One of:
#'   \itemize{
#'     \item "linear_homogeneous": Constant treatment effect (tau = effect_size)
#'     \item "nonlinear_outcome": Nonlinear outcome functions mu0, mu1
#'     \item "heterogeneous_ite": Treatment effect varies with X1
#'     \item "strong_confounding": Strong confounding with peaked propensity
#'   }
#' @param n Integer, number of observations (default: 500)
#' @param p Integer, number of covariates (default: 20)
#' @param seed Integer, random seed for reproducibility (default: 1)
#' @param effect_size Numeric, base treatment effect size (default: 1.0)
#' @param noise_sd Numeric, standard deviation of outcome noise (default: 1.0)
#' @param return_graph Logical, whether to return an igraph DAG for dowhy_gcm
#'   (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item data: data.frame with columns Y, T, X1, ..., Xp
#'     \item truth: list with ate_true, ite_true, mu0_true, mu1_true, propensity_true
#'     \item graph: igraph object (if return_graph=TRUE) or NULL
#'     \item meta: list with scenario, n, p, seed, and parameter values
#'   }
#'
#' @export
cf_benchmark_generate_data <- function(
    scenario = c("linear_homogeneous", "nonlinear_outcome", "heterogeneous_ite", "strong_confounding"),
    n = 500L,
    p = 20L,
    seed = 1L,
    effect_size = 1.0,
    noise_sd = 1.0,
    return_graph = TRUE
) {
  scenario <- match.arg(scenario)
  n <- as.integer(n)
  p <- as.integer(p)
  seed <- as.integer(seed)
  
  set.seed(seed)
  
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  
  result <- switch(
    scenario,
    "linear_homogeneous" = .dgp_linear_homogeneous(X, effect_size, noise_sd),
    "nonlinear_outcome" = .dgp_nonlinear_outcome(X, effect_size, noise_sd),
    "heterogeneous_ite" = .dgp_heterogeneous_ite(X, effect_size, noise_sd),
    "strong_confounding" = .dgp_strong_confounding(X, effect_size, noise_sd)
  )
  
  data <- data.frame(
    Y = result$Y,
    T = result$T,
    X
  )
  
  graph <- NULL
  if (return_graph) {
    graph <- .create_benchmark_dag(p)
  }
  
  list(
    data = data,
    truth = list(
      ate_true = result$ate_true,
      ite_true = result$ite_true,
      mu0_true = result$mu0_true,
      mu1_true = result$mu1_true,
      propensity_true = result$propensity_true
    ),
    graph = graph,
    meta = list(
      scenario = scenario,
      n = n,
      p = p,
      seed = seed,
      effect_size = effect_size,
      noise_sd = noise_sd
    )
  )
}

.dgp_linear_homogeneous <- function(X, effect_size, noise_sd) {
  n <- nrow(X)
  X1 <- X[, 1]
  X2 <- X[, 2]
  
  propensity <- plogis(0.5 * X1)
  T <- rbinom(n, 1, propensity)
  
  tau <- rep(effect_size, n)
  
  mu0 <- 0.5 * X1 + 0.3 * X2
  mu1 <- mu0 + tau
  
  Y <- mu0 + T * tau + rnorm(n, sd = noise_sd)
  
  list(
    Y = Y,
    T = as.integer(T),
    ate_true = effect_size,
    ite_true = tau,
    mu0_true = mu0,
    mu1_true = mu1,
    propensity_true = propensity
  )
}

.dgp_nonlinear_outcome <- function(X, effect_size, noise_sd) {
  n <- nrow(X)
  X1 <- X[, 1]
  X2 <- X[, 2]
  X3 <- if (ncol(X) >= 3) X[, 3] else rnorm(n)
  
  propensity <- plogis(0.5 * X1 + 0.2 * X2)
  T <- rbinom(n, 1, propensity)
  
  tau <- rep(effect_size, n)
  
  mu0 <- sin(X1) + X2^2 + 0.5 * X1 * X2 + 0.3 * X3
  mu1 <- mu0 + tau
  
  Y <- mu0 + T * tau + rnorm(n, sd = noise_sd)
  
  list(
    Y = Y,
    T = as.integer(T),
    ate_true = effect_size,
    ite_true = tau,
    mu0_true = mu0,
    mu1_true = mu1,
    propensity_true = propensity
  )
}

.dgp_heterogeneous_ite <- function(X, effect_size, noise_sd) {
  n <- nrow(X)
  X1 <- X[, 1]
  X2 <- X[, 2]
  
  propensity <- plogis(0.5 * X1)
  T <- rbinom(n, 1, propensity)
  
  tau <- effect_size * (1 + 0.5 * X1)
  
  mu0 <- 0.5 * X1 + 0.3 * X2
  mu1 <- mu0 + tau
  
  Y <- mu0 + T * tau + rnorm(n, sd = noise_sd)
  
  ate_true <- mean(tau)
  
  list(
    Y = Y,
    T = as.integer(T),
    ate_true = ate_true,
    ite_true = tau,
    mu0_true = mu0,
    mu1_true = mu1,
    propensity_true = propensity
  )
}

.dgp_strong_confounding <- function(X, effect_size, noise_sd) {
  n <- nrow(X)
  X1 <- X[, 1]
  X2 <- X[, 2]
  
  propensity <- plogis(2.0 * X1 + 1.5 * X2)
  T <- rbinom(n, 1, propensity)
  
  tau <- rep(effect_size, n)
  
  mu0 <- 1.0 * X1 + 0.8 * X2
  mu1 <- mu0 + tau
  
  Y <- mu0 + T * tau + rnorm(n, sd = noise_sd)
  
  list(
    Y = Y,
    T = as.integer(T),
    ate_true = effect_size,
    ite_true = tau,
    mu0_true = mu0,
    mu1_true = mu1,
    propensity_true = propensity
  )
}

.create_benchmark_dag <- function(p) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    return(NULL)
  }

  edges <- rbind(
    c("X1", "T"),
    c("T", "Y"),
    c("X1", "Y"),
    c("X2", "Y")
  )

  igraph::graph_from_edgelist(edges, directed = TRUE)
}

# ============================================================================
# Standalone DGP functions for high-dimensional benchmarking
# ============================================================================

#' Generate baseline DGP data (S1)
#'
#' Generates data with sparse confounding where only the first 10 covariates
#' affect treatment and outcome.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, propensity_score
#' @export
dgp_baseline <- function(n = 500, p = 50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # Sparse confounding: first 10 variables
  n_conf <- min(10, p)
  beta_t <- c(rep(0.3, n_conf), rep(0, p - n_conf))
  beta_y <- c(rep(0.2, n_conf), rep(0, p - n_conf))

  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "baseline",
    dgp_params = list(n = n, p = p)
  )
}

#' Generate high-dimensional DGP data (S2-S3 dimension sweep)
#'
#' Generates data with specified n and p to test methods under different
#' n/p ratios.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, propensity_score
#' @export
dgp_dimension_sweep <- function(n = 200, p = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  n_conf <- min(10, p)
  beta_t <- c(rep(0.5, n_conf), rep(0, p - n_conf))
  beta_y <- c(rep(0.3, n_conf), rep(0, p - n_conf))

  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "dimension_sweep",
    dgp_params = list(n = n, p = p, n_p_ratio = n / p)
  )
}

#' Generate linear heterogeneous treatment effect DGP (S4a/e/f)
#'
#' Treatment effect varies linearly with X1 and X2.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param strength Numeric, heterogeneity strength
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_heterogeneous_linear <- function(n = 500, p = 500, strength = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau_base <- 2.0
  tau <- tau_base + strength * X[,1] + 0.3 * strength * X[,2]

  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = mean(tau),
    true_ite = tau,
    ite_sd = stats::sd(tau),
    propensity_score = as.numeric(ps),
    dgp_name = "heterogeneous_linear",
    dgp_params = list(n = n, p = p, strength = strength)
  )
}

#' Generate nonlinear heterogeneous treatment effect DGP (S4b)
#'
#' Treatment effect varies nonlinearly with covariates.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_heterogeneous_nonlinear <- function(n = 500, p = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0 + 1.0 * sin(pi * X[,1]) + 0.5 * X[,2]^2

  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = mean(tau),
    true_ite = tau,
    propensity_score = as.numeric(ps),
    dgp_name = "heterogeneous_nonlinear",
    dgp_params = list(n = n, p = p)
  )
}

#' Generate subgroup heterogeneous treatment effect DGP (S4c/g/h)
#'
#' Treatment effect varies by discrete subgroups.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param subgroup_props Numeric vector, proportion in each subgroup
#' @param subgroup_effects Numeric vector, treatment effect in each subgroup
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, subgroup
#' @export
dgp_heterogeneous_subgroup <- function(n = 500, p = 500,
                                       subgroup_props = c(1/3, 1/3, 1/3),
                                       subgroup_effects = c(0.5, 2.0, 4.0),
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  # Subgroup assignment based on X1
  breaks <- stats::qnorm(cumsum(subgroup_props)[1:(length(subgroup_props)-1)])
  subgroup <- cut(X[,1], breaks = c(-Inf, breaks, Inf), labels = FALSE)

  tau <- subgroup_effects[subgroup]

  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = sum(subgroup_props * subgroup_effects),
    true_ite = tau,
    subgroup = subgroup,
    subgroup_n = table(subgroup),
    propensity_score = as.numeric(ps),
    dgp_name = "heterogeneous_subgroup",
    dgp_params = list(n = n, p = p, subgroup_props = subgroup_props,
                      subgroup_effects = subgroup_effects)
  )
}

#' Generate qualitative interaction DGP (S4d)
#'
#' Treatment effect can be positive or negative depending on X1.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, prop_harmed
#' @export
dgp_heterogeneous_qualitative <- function(n = 500, p = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0 * X[,1]  # Negative for X1 < 0

  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = mean(tau),
    true_ite = tau,
    prop_harmed = mean(tau < 0),
    prop_benefited = mean(tau > 0),
    propensity_score = as.numeric(ps),
    dgp_name = "heterogeneous_qualitative",
    dgp_params = list(n = n, p = p)
  )
}

#' Generate nonlinear confounding DGP (S5)
#'
#' Confounding relationships include nonlinear terms.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param nonlinear_type Character: "quadratic", "trigonometric", "interaction",
#'   "combined", "threshold"
#' @param strength Numeric, nonlinearity strength (0=linear, 1=full nonlinear)
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_nonlinear_confounding <- function(n = 500, p = 500,
                                       nonlinear_type = "combined",
                                       strength = 1.0,
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  linear_t <- 0.3 * X[,1] + 0.3 * X[,2]
  linear_y <- 0.2 * X[,1] + 0.2 * X[,2]

  nonlinear_t <- nonlinear_y <- 0

  if (nonlinear_type %in% c("quadratic", "combined")) {
    nonlinear_t <- nonlinear_t + 0.5 * X[,1]^2
    nonlinear_y <- nonlinear_y + 0.3 * X[,1]^2
  }

  if (nonlinear_type %in% c("trigonometric", "combined")) {
    nonlinear_t <- nonlinear_t + 0.3 * sin(pi * X[,3])
    nonlinear_y <- nonlinear_y + 0.2 * sin(pi * X[,2])
  }

  if (nonlinear_type %in% c("interaction", "combined")) {
    nonlinear_t <- nonlinear_t + 0.2 * X[,1] * X[,2]
    nonlinear_y <- nonlinear_y + 0.1 * X[,1] * X[,3]
  }

  if (nonlinear_type == "threshold") {
    nonlinear_t <- 0.5 * as.numeric(X[,1] > 0)
    nonlinear_y <- 0.3 * as.numeric(X[,2] > 0.5)
  }

  ps <- stats::plogis(linear_t + strength * nonlinear_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  Y <- linear_y + strength * nonlinear_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "nonlinear_confounding",
    dgp_params = list(n = n, p = p, nonlinear_type = nonlinear_type, strength = strength)
  )
}

#' Generate dense confounding DGP (S6)
#'
#' Many covariates affect both treatment and outcome.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param n_confounders Integer, number of confounding variables
#' @param coef_scaling Character: "fixed", "sqrt", "linear"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_dense_confounding <- function(n = 500, p = 500,
                                   n_confounders = 100,
                                   coef_scaling = "fixed",
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  base_coef <- 0.1
  coef_size <- switch(coef_scaling,
    "fixed"  = base_coef,
    "sqrt"   = base_coef / sqrt(n_confounders),
    "linear" = base_coef * 10 / n_confounders,
    base_coef
  )

  n_conf <- min(n_confounders, p)
  beta_t <- c(rep(coef_size, n_conf), rep(0, p - n_conf))
  beta_y <- c(rep(coef_size, n_conf), rep(0, p - n_conf))

  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    confounding_info = list(
      n_confounders = n_conf,
      coef_size = coef_size,
      coef_scaling = coef_scaling
    ),
    dgp_name = "dense_confounding",
    dgp_params = list(n = n, p = p, n_confounders = n_confounders,
                      coef_scaling = coef_scaling)
  )
}

#' Generate weak overlap DGP (S7)
#'
#' Creates propensity scores with limited overlap between treatment groups.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param overlap_strength Character: "good", "moderate", "weak", "extreme"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, ps_summary
#' @export
dgp_weak_overlap <- function(n = 500, p = 500,
                             overlap_strength = "weak",
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  coef_scale <- switch(overlap_strength,
    "good"     = 0.3,
    "moderate" = 0.6,
    "weak"     = 1.2,
    "extreme"  = 2.0,
    1.2
  )

  beta_t <- c(rep(coef_scale, 10), rep(0, p - 10))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  ps_summary <- list(
    min = min(ps),
    max = max(ps),
    mean = mean(ps),
    prop_extreme = mean(ps < 0.05 | ps > 0.95)
  )

  tau <- 2.0
  beta_y <- c(rep(0.3, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    ps_summary = ps_summary,
    overlap_strength = overlap_strength,
    dgp_name = "weak_overlap",
    dgp_params = list(n = n, p = p, overlap_strength = overlap_strength)
  )
}

#' Generate covariate shift DGP (S8)
#'
#' Training and test distributions differ in covariate distribution.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param shift_type Character: "mean", "variance", "correlation", "combined"
#' @param shift_magnitude Numeric, magnitude of the shift
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_covariate_shift <- function(n = 500, p = 500,
                                 shift_type = "mean",
                                 shift_magnitude = 1.0,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate shifted covariates
  if (shift_type == "mean") {
    X <- matrix(stats::rnorm(n * p, mean = shift_magnitude), n, p)
  } else if (shift_type == "variance") {
    X <- matrix(stats::rnorm(n * p, sd = 1 + shift_magnitude), n, p)
  } else {
    X <- matrix(stats::rnorm(n * p), n, p)
  }
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "covariate_shift",
    dgp_params = list(n = n, p = p, shift_type = shift_type,
                      shift_magnitude = shift_magnitude)
  )
}

#' Generate correlated confounding DGP (S9)
#'
#' Confounders are correlated with each other.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param correlation_type Character: "block", "ar1", "factor"
#' @param correlation_strength Numeric, strength of correlation (0-1)
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_correlated_confounding <- function(n = 500, p = 500,
                                        correlation_type = "block",
                                        correlation_strength = 0.5,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate correlated covariates
  if (correlation_type == "ar1") {
    # AR(1) correlation structure
    rho <- correlation_strength
    Sigma <- rho^abs(outer(1:p, 1:p, "-"))
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  } else if (correlation_type == "block") {
    # Block correlation: first 10 variables correlated
    block_size <- min(10, p)
    Sigma <- diag(p)
    Sigma[1:block_size, 1:block_size] <- correlation_strength
    diag(Sigma) <- 1
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  } else {
    # Factor model
    n_factors <- 3
    loadings <- matrix(stats::rnorm(p * n_factors), p, n_factors)
    factors <- matrix(stats::rnorm(n * n_factors), n, n_factors)
    X <- factors %*% t(loadings) * sqrt(correlation_strength) +
         matrix(stats::rnorm(n * p), n, p) * sqrt(1 - correlation_strength)
  }
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "correlated_confounding",
    dgp_params = list(n = n, p = p, correlation_type = correlation_type,
                      correlation_strength = correlation_strength)
  )
}

#' Generate unobserved confounding DGP (S10)
#'
#' Includes unmeasured confounders that affect both treatment and outcome.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of observed covariates
#' @param n_unobserved Integer, number of unobserved confounders
#' @param unobserved_strength Numeric, strength of unobserved confounding
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, X_unobserved
#' @export
dgp_unobserved_confounding <- function(n = 500, p = 500,
                                        n_unobserved = 5,
                                        unobserved_strength = 0.5,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # Unobserved confounders
  U <- matrix(stats::rnorm(n * n_unobserved), n, n_unobserved)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  beta_t_u <- rep(unobserved_strength, n_unobserved)

  ps <- stats::plogis(X %*% beta_t + U %*% beta_t_u)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  beta_y_u <- rep(unobserved_strength, n_unobserved)

  Y <- X %*% beta_y + U %*% beta_y_u + tau * T + stats::rnorm(n)

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    X_unobserved = U,
    dgp_name = "unobserved_confounding",
    dgp_params = list(n = n, p = p, n_unobserved = n_unobserved,
                      unobserved_strength = unobserved_strength)
  )
}

#' Generate collider DGP (S11)
#'
#' Includes a collider variable that is affected by both treatment and outcome.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param collider_strength Numeric, strength of collider relationships
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, collider
#' @export
dgp_collider <- function(n = 500, p = 500,
                          collider_strength = 0.5,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + stats::rnorm(n)

  # Collider is affected by both T and Y
  collider <- collider_strength * T + collider_strength * Y + stats::rnorm(n)

  # Include collider as first column of X so methods that adjust for all X columns see it
  X <- cbind(collider, X)
  colnames(X) <- c("collider", paste0("X", 1:p))

  list(
    X = X,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    collider = collider,
    dgp_name = "collider",
    dgp_params = list(n = n, p = p, collider_strength = collider_strength)
  )
}

#' Generate nonlinear outcome DGP (S12)
#'
#' Generates data with nonlinear outcome model but linear propensity score.
#' This breaks gformula (OLS) while leaving hdps (IPW-only) unaffected.
#' hdml and tmle are partially affected since their outcome models are also linear.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param nonlinearity Character, "moderate" or "severe"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, propensity_score, dgp_name, dgp_params
#' @export
dgp_nonlinear_outcome <- function(n = 500, p = 500,
                                   nonlinearity = "moderate",
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # Linear propensity score (all methods can estimate this correctly)
  beta_t <- c(0.5, -0.3, 0.2, rep(0, p - 3))
  ps <- stats::plogis(X %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0

  # Nonlinear outcome surface
  if (nonlinearity == "moderate") {
    mu0 <- X[, 1]^2 + X[, 2] * X[, 3] + sin(X[, 1])
  } else {
    mu0 <- X[, 1]^3 / 3 + exp(X[, 2] / 2) + X[, 3] * X[, 4] * (X[, 1] > 0) +
      sin(2 * X[, 1]) * cos(X[, 2])
  }

  Y <- as.numeric(mu0 + tau * T + stats::rnorm(n))

  list(
    X = X,
    T = T,
    Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "nonlinear_outcome",
    dgp_params = list(n = n, p = p, nonlinearity = nonlinearity)
  )
}

#' Generate nonlinear propensity score DGP (S13)
#'
#' Generates data with nonlinear propensity score but linear outcome.
#' This breaks hdps (IPW with linear PS model) while gformula (no PS) is unaffected.
#' hdml and tmle are partially protected by double robustness.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param nonlinearity Character, "moderate" or "severe"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, propensity_score, dgp_name, dgp_params
#' @export
dgp_nonlinear_propensity <- function(n = 500, p = 500,
                                      nonlinearity = "moderate",
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # Nonlinear propensity score
  if (nonlinearity == "moderate") {
    logit_ps <- X[, 1]^2 - 0.5 + X[, 2] * X[, 3]
  } else {
    logit_ps <- sin(2 * X[, 1]) + X[, 2]^2 * (X[, 3] > 0) +
      abs(X[, 1] * X[, 2]) - 1
  }
  ps <- stats::plogis(logit_ps)
  ps <- pmin(pmax(ps, 0.05), 0.95)
  T <- stats::rbinom(n, 1, ps)

  # Linear outcome (gformula can handle this correctly)
  tau <- 2.0
  beta_y <- c(1.0, -0.8, 0.6, -0.4, 0.2, rep(0, p - 5))
  Y <- as.numeric(X %*% beta_y + tau * T + stats::rnorm(n))

  list(
    X = X,
    T = T,
    Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "nonlinear_propensity",
    dgp_params = list(n = n, p = p, nonlinearity = nonlinearity)
  )
}

#' Generate double nonlinear DGP (S14)
#'
#' Both outcome and propensity score are nonlinear.
#' All linear methods struggle. Doubly robust methods lose DR protection
#' since both nuisance models are misspecified.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, propensity_score, dgp_name, dgp_params
#' @export
dgp_double_nonlinear <- function(n = 500, p = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # Nonlinear propensity score
  logit_ps <- sin(2 * X[, 1]) + X[, 2]^2 * (X[, 3] > 0) +
    abs(X[, 1] * X[, 2]) - 1
  ps <- stats::plogis(logit_ps)
  ps <- pmin(pmax(ps, 0.05), 0.95)
  T <- stats::rbinom(n, 1, ps)

  # Nonlinear outcome
  tau <- 2.0
  mu0 <- X[, 1]^3 / 3 + exp(X[, 2] / 2) + X[, 3] * X[, 4] * (X[, 1] > 0) +
    sin(2 * X[, 1]) * cos(X[, 2])
  Y <- as.numeric(mu0 + tau * T + stats::rnorm(n))

  list(
    X = X,
    T = T,
    Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "double_nonlinear",
    dgp_params = list(n = n, p = p)
  )
}
