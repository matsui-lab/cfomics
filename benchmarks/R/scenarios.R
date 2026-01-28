# benchmarks/R/scenarios.R
# Scenario data generation - dispatches to cfomics DGP functions

generate_scenario_data <- function(scenario, rep, base_seed) {
  seed <- base_seed + rep

  dgp_fn <- get(scenario$dgp, envir = asNamespace("cfomics"))

  args <- c(
    list(n = scenario$n, p = scenario$p, seed = seed),
    scenario$params
  )

  dgp <- do.call(dgp_fn, args)

  # Build data.frame for cf_fit
  p <- ncol(dgp$X)
  df <- data.frame(Y = dgp$Y, T = dgp$T, dgp$X)

  list(
    data = df,
    truth = list(
      ate_true = dgp$true_ate,
      ite_true = dgp$true_ite
    ),
    scenario_id = scenario$id,
    n = scenario$n,
    p = p,
    rep = rep,
    seed = seed
  )
}
