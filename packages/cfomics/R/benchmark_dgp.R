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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
    # Use empirical mean of ITEs for consistency with other heterogeneous DGPs
    # (theoretical ATE may differ from sample ATE due to finite sample variation)
    true_ate = mean(tau),
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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
    Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- c("collider", paste0("X", seq_len(p)))

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

#' Generate benchmark data by DGP name
#'
#' Dispatches to the appropriate DGP function based on a string name.
#'
#' @param dgp Character, the DGP name (e.g., "baseline", "nonlinear_outcome")
#' @param params List of parameters passed to the DGP function
#' @return List with X, T, Y, true_ate, true_ite, etc.
#' @export
generate_benchmark_data <- function(dgp, params = list()) {
  dgp_fn <- switch(dgp,
    "baseline"               = dgp_baseline,
    "dimension_sweep"        = dgp_dimension_sweep,
    "heterogeneous_linear"   = dgp_heterogeneous_linear,
    "heterogeneous_nonlinear" = dgp_heterogeneous_nonlinear,
    "heterogeneous_subgroup" = dgp_heterogeneous_subgroup,
    "heterogeneous_qualitative" = dgp_heterogeneous_qualitative,
    "nonlinear_confounding"  = dgp_nonlinear_confounding,
    "dense_confounding"      = dgp_dense_confounding,
    "weak_overlap"           = dgp_weak_overlap,
    "covariate_shift"        = dgp_covariate_shift,
    "correlated_confounding" = dgp_correlated_confounding,
    "unobserved_confounding" = dgp_unobserved_confounding,
    "collider"               = dgp_collider,
    "nonlinear_outcome"      = dgp_nonlinear_outcome,
    "nonlinear_propensity"   = dgp_nonlinear_propensity,
    "double_nonlinear"       = dgp_double_nonlinear,
    "missing_data"           = dgp_missing_data,
    "high_dimensional_omics" = dgp_high_dimensional_omics,
    rlang::abort(paste0("Unknown DGP: '", dgp, "'"), class = "cfomics_unknown_dgp")
  )
  do.call(dgp_fn, params)
}

#' Generate nonlinear outcome DGP (S12)
#'
#' Generates data with nonlinear confounding: both the outcome and propensity
#' score depend on the same nonlinear function h(X). This creates omitted-variable
#' bias for methods using linear outcome models (gformula).
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
  colnames(X) <- paste0("X", seq_len(p))

  tau <- 2.0

  # Nonlinear function h(X) drives both outcome and treatment
  if (nonlinearity == "moderate") {
    h <- X[, 1]^2 + X[, 2] * X[, 3] + sin(X[, 1])
  } else {
    h <- X[, 1]^3 / 3 + exp(X[, 2] / 2) + X[, 3] * X[, 4] * (X[, 1] > 0) +
      sin(2 * X[, 1]) * cos(X[, 2])
  }

  # PS depends on h(X) — creates nonlinear confounding
  ps_coef <- if (nonlinearity == "moderate") 1.2 else 0.8
  h_centered <- h - mean(h)
  ps <- stats::plogis(ps_coef * h_centered)
  ps <- pmin(pmax(ps, 0.05), 0.95)
  T <- stats::rbinom(n, 1, ps)

  # Outcome depends on same h(X)
  Y <- as.numeric(h + tau * T + stats::rnorm(n))

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
  colnames(X) <- paste0("X", seq_len(p))

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
  colnames(X) <- paste0("X", seq_len(p))

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

#' Generate missing data DGP (S15)
#'
#' Generates data with missing values under MCAR or MAR mechanisms.
#' Omics data often has missing values, and this scenario tests how
#' causal inference methods handle incomplete covariate data.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param missing_rate Numeric, proportion of missing values (0-1)
#' @param missing_type Character, "MCAR" or "MAR"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, missing_mask
#' @export
dgp_missing_data <- function(n = 500, p = 500,
                              missing_rate = 0.2,
                              missing_type = "MCAR",
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate complete data first
  X_complete <- matrix(stats::rnorm(n * p), n, p)
  colnames(X_complete) <- paste0("X", seq_len(p))

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X_complete %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X_complete %*% beta_y + tau * T + stats::rnorm(n)

  # Create missing mask
  if (missing_type == "MCAR") {
    # Missing Completely At Random: each cell has equal probability of missing
    missing_mask <- matrix(
      stats::rbinom(n * p, 1, missing_rate),
      n, p
    )
  } else {
    # MAR: missingness depends on observed values (X1 always observed)
    missing_prob <- stats::plogis(-1 + 0.5 * X_complete[, 1])
    missing_mask <- matrix(0, n, p)
    for (j in 2:p) {
      missing_mask[, j] <- stats::rbinom(n, 1, missing_prob * missing_rate * 2)
    }
  }

  # Apply missing values
  X <- X_complete
  X[missing_mask == 1] <- NA

  # Validate missing rate
  actual_missing_rate <- mean(is.na(X))
  if (abs(actual_missing_rate - missing_rate) > 0.1) {
    warning(sprintf("Actual missing rate (%.2f) differs from requested (%.2f)",
                    actual_missing_rate, missing_rate))
  }

  # Validate dimensions match

  stopifnot(
    identical(dim(X), dim(X_complete)),
    all(!is.na(X_complete))
  )

  list(
    X = X,
    X_complete = X_complete,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    missing_mask = missing_mask,
    missing_rate_actual = mean(missing_mask),
    dgp_name = "missing_data",
    dgp_params = list(n = n, p = p, missing_rate = missing_rate,
                      missing_type = missing_type)
  )
}

#' High-dimensional omics-like DGP (S16)
#'
#' Generates data with correlated features, sparse true effects,
#' and block correlation structure typical of gene expression data.
#' This DGP is designed to mimic the characteristics of real omics
#' data where genes within pathways are correlated and only a sparse
#' subset of features drives the treatment effect.
#'
#' @param n Integer, sample size
#' @param p Integer, number of features (default 500)
#' @param n_blocks Integer, number of correlation blocks (default 10)
#' @param block_cor Numeric, within-block correlation (default 0.5)
#' @param sparsity Numeric, proportion of non-zero effects (default 0.05)
#' @param seed Integer, random seed (optional)
#' @return List with:
#'   \itemize{
#'     \item data: data.frame with columns Y, T, X1, ..., Xp
#'     \item truth: list with ate_true, ite_true, beta, active_features
#'     \item structure: list with n_blocks, block_cor, sparsity
#'   }
#' @export
#' @examples
#' \donttest{
#' result <- dgp_high_dimensional_omics(n = 100, p = 200, seed = 42)
#' dim(result$data)
#' # Check sparse effects
#' sum(result$truth$beta != 0)
#' }
dgp_high_dimensional_omics <- function(
    n,
    p = 500,
    n_blocks = 10,
    block_cor = 0.5,
    sparsity = 0.05,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  if (!requireNamespace("MASS", quietly = TRUE)) {
    rlang::abort(
      "Package 'MASS' is required for block-correlated data generation",
      class = "cfomics_missing_package"
    )
  }

  # Generate block-correlated features
  block_size <- p %/% n_blocks
  remainder <- p %% n_blocks

  X_list <- lapply(seq_len(n_blocks), function(b) {
    bs <- if (b <= remainder) block_size + 1 else block_size
    Sigma <- matrix(block_cor, bs, bs)
    diag(Sigma) <- 1
    MASS::mvrnorm(n, mu = rep(0, bs), Sigma = Sigma)
  })
  X <- do.call(cbind, X_list)
  colnames(X) <- paste0("X", seq_len(ncol(X)))

  # Sparse true effects
  n_active <- max(1, ceiling(p * sparsity))
  active_idx <- sample(p, n_active)
  beta <- rep(0, p)
  beta[active_idx] <- stats::rnorm(n_active, mean = 0, sd = 0.5)

  # Treatment assignment (confounded by first few active features)
  n_confounders <- min(5, n_active)
  ps_linear <- X[, active_idx[seq_len(n_confounders)], drop = FALSE] %*%
    beta[active_idx[seq_len(n_confounders)]]
  ps <- stats::plogis(as.numeric(ps_linear))
  # Clip propensity scores to ensure overlap
  ps <- pmin(pmax(ps, 0.05), 0.95)
  T_var <- stats::rbinom(n, 1, ps)

  # Outcome with treatment effect heterogeneity
  base_effect <- 2
  het_effect <- X[, active_idx[1], drop = FALSE] * 0.5
  ite_true <- base_effect + as.numeric(het_effect)
  Y <- ite_true * T_var + as.numeric(X %*% beta) + stats::rnorm(n)

  data <- data.frame(Y = Y, T = T_var, X)

  list(
    data = data,
    truth = list(
      ate_true = mean(ite_true),
      ite_true = ite_true,
      beta = beta,
      active_features = active_idx,
      propensity_score = ps
    ),
    structure = list(
      n_blocks = n_blocks,
      block_cor = block_cor,
      sparsity = sparsity
    ),
    dgp_name = "high_dimensional_omics",
    dgp_params = list(n = n, p = p, n_blocks = n_blocks,
                      block_cor = block_cor, sparsity = sparsity)
  )
}
