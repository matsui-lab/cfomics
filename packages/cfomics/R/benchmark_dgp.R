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
