test_that("cf_benchmark_generate_data returns correct structure", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 100L,
    p = 10L,
    seed = 123L
  )
  
  expect_type(result, "list")
  expect_true("data" %in% names(result))
  expect_true("truth" %in% names(result))
  expect_true("graph" %in% names(result))
  expect_true("meta" %in% names(result))
  
  expect_s3_class(result$data, "data.frame")
  expect_equal(nrow(result$data), 100)
  
  expect_true("Y" %in% names(result$data))
  expect_true("T" %in% names(result$data))
  expect_true("X1" %in% names(result$data))
  expect_true("X10" %in% names(result$data))
})

test_that("cf_benchmark_generate_data produces reproducible results with same seed", {
  result1 <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 50L,
    p = 5L,
    seed = 42L
  )
  
  result2 <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 50L,
    p = 5L,
    seed = 42L
  )
  
  expect_equal(result1$data$Y, result2$data$Y)
  expect_equal(result1$data$T, result2$data$T)
  expect_equal(result1$data$X1, result2$data$X1)
  expect_equal(result1$truth$ate_true, result2$truth$ate_true)
  expect_equal(result1$truth$ite_true, result2$truth$ite_true)
})

test_that("cf_benchmark_generate_data treatment is binary 0/1 integer", {
  for (scenario in c("linear_homogeneous", "nonlinear_outcome", 
                     "heterogeneous_ite", "strong_confounding")) {
    result <- cfomics:::cf_benchmark_generate_data(
      scenario = scenario,
      n = 100L,
      p = 5L,
      seed = 1L
    )
    
    T_vals <- result$data$T
    expect_type(T_vals, "integer")
    expect_true(all(T_vals %in% c(0L, 1L)))
  }
})

test_that("cf_benchmark_generate_data truth values have correct types", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "heterogeneous_ite",
    n = 100L,
    p = 5L,
    seed = 1L
  )
  
  expect_type(result$truth$ate_true, "double")
  expect_length(result$truth$ate_true, 1)
  
  expect_type(result$truth$ite_true, "double")
  expect_length(result$truth$ite_true, 100)
  
  expect_type(result$truth$mu0_true, "double")
  expect_length(result$truth$mu0_true, 100)
  
  expect_type(result$truth$mu1_true, "double")
  expect_length(result$truth$mu1_true, 100)
  
  expect_type(result$truth$propensity_true, "double")
  expect_length(result$truth$propensity_true, 100)
})

test_that("cf_benchmark_generate_data linear_homogeneous has constant ITE", {
  effect_size <- 2.5
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 100L,
    p = 5L,
    seed = 1L,
    effect_size = effect_size
  )
  
  expect_equal(result$truth$ate_true, effect_size)
  expect_true(all(result$truth$ite_true == effect_size))
})

test_that("cf_benchmark_generate_data heterogeneous_ite has varying ITE", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "heterogeneous_ite",
    n = 100L,
    p = 5L,
    seed = 1L,
    effect_size = 1.0
  )
  
  ite_sd <- sd(result$truth$ite_true)
  expect_true(ite_sd > 0)
  
  expect_equal(result$truth$ate_true, mean(result$truth$ite_true))
})

test_that("cf_benchmark_generate_data returns igraph when return_graph=TRUE", {
  skip_if_not_installed("igraph")
  
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 50L,
    p = 5L,
    seed = 1L,
    return_graph = TRUE
  )
  
  expect_s3_class(result$graph, "igraph")
  
  edges <- igraph::as_data_frame(result$graph, what = "edges")
  expect_true(any(edges$from == "X1" & edges$to == "T"))
  expect_true(any(edges$from == "T" & edges$to == "Y"))
  expect_true(any(edges$from == "X1" & edges$to == "Y"))
})

test_that("cf_benchmark_generate_data returns NULL graph when return_graph=FALSE", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 50L,
    p = 5L,
    seed = 1L,
    return_graph = FALSE
  )
  
  expect_null(result$graph)
})

test_that("cf_benchmark_generate_data meta contains correct values", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "strong_confounding",
    n = 200L,
    p = 15L,
    seed = 99L,
    effect_size = 3.0,
    noise_sd = 0.5
  )
  
  expect_equal(result$meta$scenario, "strong_confounding")
  expect_equal(result$meta$n, 200L)
  expect_equal(result$meta$p, 15L)
  expect_equal(result$meta$seed, 99L)
  expect_equal(result$meta$effect_size, 3.0)
  expect_equal(result$meta$noise_sd, 0.5)
})

test_that("cf_benchmark_generate_data all scenarios work", {
  scenarios <- c("linear_homogeneous", "nonlinear_outcome",
                 "heterogeneous_ite", "strong_confounding")

  for (scenario in scenarios) {
    result <- cfomics:::cf_benchmark_generate_data(
      scenario = scenario,
      n = 50L,
      p = 5L,
      seed = 1L
    )

    expect_equal(nrow(result$data), 50)
    expect_true(is.numeric(result$truth$ate_true))
    expect_length(result$truth$ite_true, 50)
  }
})

# ============================================================================
# Tests for standalone DGP functions
# ============================================================================

test_that("dgp_baseline generates valid data", {
  dgp <- dgp_baseline(n = 100, p = 50, seed = 123)

  expect_type(dgp, "list")
  expect_equal(nrow(dgp$X), 100)
  expect_equal(ncol(dgp$X), 50)
  expect_equal(length(dgp$T), 100)
  expect_equal(length(dgp$Y), 100)
  expect_true("true_ate" %in% names(dgp))
  expect_true("true_ite" %in% names(dgp))
  expect_true("propensity_score" %in% names(dgp))
  expect_equal(dgp$dgp_name, "baseline")
})

test_that("dgp_dimension_sweep generates valid high-dim data", {
  dgp <- dgp_dimension_sweep(n = 100, p = 500, seed = 123)

  expect_equal(nrow(dgp$X), 100)
  expect_equal(ncol(dgp$X), 500)
  expect_true(dgp$dgp_params$n_p_ratio < 1)
  expect_equal(dgp$dgp_name, "dimension_sweep")
})

test_that("dgp_heterogeneous_linear generates valid HTE data", {
  dgp <- dgp_heterogeneous_linear(n = 100, p = 50, strength = 1.0, seed = 123)

  expect_equal(length(unique(dgp$true_ite)), 100)  # ITE varies
  expect_true(sd(dgp$true_ite) > 0)
  expect_equal(dgp$dgp_name, "heterogeneous_linear")
})

test_that("dgp_heterogeneous_nonlinear generates valid data", {
  dgp <- dgp_heterogeneous_nonlinear(n = 100, p = 50, seed = 123)

  expect_true(sd(dgp$true_ite) > 0)
  expect_equal(dgp$dgp_name, "heterogeneous_nonlinear")
})

test_that("dgp_heterogeneous_subgroup generates subgroup data", {
  dgp <- dgp_heterogeneous_subgroup(n = 300, p = 50, seed = 123)

  expect_true("subgroup" %in% names(dgp))
  expect_equal(length(dgp$subgroup), 300)
  expect_equal(dgp$dgp_name, "heterogeneous_subgroup")
})

test_that("dgp_heterogeneous_qualitative generates qualitative interaction data", {
  dgp <- dgp_heterogeneous_qualitative(n = 200, p = 50, seed = 123)

  expect_true("prop_harmed" %in% names(dgp))
  expect_true("prop_benefited" %in% names(dgp))
  expect_true(dgp$prop_harmed > 0)  # Some should be harmed
  expect_true(dgp$prop_benefited > 0)  # Some should benefit
  expect_equal(dgp$dgp_name, "heterogeneous_qualitative")
})

test_that("dgp_nonlinear_confounding generates valid data", {
  types <- c("quadratic", "trigonometric", "interaction", "combined", "threshold")

  for (type in types) {
    dgp <- dgp_nonlinear_confounding(n = 100, p = 50, nonlinear_type = type, seed = 123)
    expect_equal(length(dgp$Y), 100)
    expect_equal(dgp$dgp_params$nonlinear_type, type)
  }
})

test_that("dgp_dense_confounding generates valid data", {
  dgp <- dgp_dense_confounding(n = 100, p = 200, n_confounders = 100, seed = 123)

  expect_equal(length(dgp$Y), 100)
  expect_true("confounding_info" %in% names(dgp))
  expect_equal(dgp$confounding_info$n_confounders, 100)
  expect_equal(dgp$dgp_name, "dense_confounding")
})

test_that("dgp_weak_overlap generates valid data with varying overlap", {
  strengths <- c("good", "moderate", "weak", "extreme")

  for (strength in strengths) {
    dgp <- dgp_weak_overlap(n = 200, p = 50, overlap_strength = strength, seed = 123)
    expect_equal(length(dgp$Y), 200)
    expect_true("ps_summary" %in% names(dgp))
    expect_equal(dgp$overlap_strength, strength)
  }
})

test_that("dgp_covariate_shift generates valid data", {
  dgp <- dgp_covariate_shift(n = 100, p = 50, shift_type = "mean", seed = 123)

  expect_equal(length(dgp$Y), 100)
  expect_equal(dgp$dgp_params$shift_type, "mean")
  expect_equal(dgp$dgp_name, "covariate_shift")
})

test_that("dgp_correlated_confounding generates valid data", {
  skip_if_not_installed("MASS")

  types <- c("block", "ar1", "factor")

  for (type in types) {
    dgp <- dgp_correlated_confounding(n = 100, p = 50, correlation_type = type, seed = 123)
    expect_equal(length(dgp$Y), 100)
    expect_equal(dgp$dgp_params$correlation_type, type)
  }
})

test_that("dgp_unobserved_confounding generates valid data with unobserved confounders", {
  dgp <- dgp_unobserved_confounding(n = 100, p = 50, n_unobserved = 5, seed = 123)

  expect_equal(length(dgp$Y), 100)
  expect_true("X_unobserved" %in% names(dgp))
  expect_equal(ncol(dgp$X_unobserved), 5)
  expect_equal(dgp$dgp_name, "unobserved_confounding")
})

test_that("dgp_collider generates valid data with collider", {
  dgp <- dgp_collider(n = 100, p = 50, seed = 123)

  expect_equal(length(dgp$Y), 100)
  expect_true("collider" %in% names(dgp))
  expect_equal(length(dgp$collider), 100)
  expect_equal(dgp$dgp_name, "collider")
})

test_that("dgp_collider includes collider in X matrix", {
  set.seed(42)
  result <- dgp_collider(n = 200, p = 10, collider_strength = 2.0)
  expect_equal(ncol(result$X), 11)
  collider_col <- result$X[, 1]
  cor_T <- abs(cor(collider_col, result$T))
  cor_Y <- abs(cor(collider_col, result$Y))
  expect_gt(cor_T, 0.1)
  expect_gt(cor_Y, 0.1)
})

test_that("dgp_collider with different strengths produces different X matrices", {
  set.seed(42)
  r1 <- dgp_collider(n = 500, p = 10, collider_strength = 0.5)
  set.seed(42)
  r2 <- dgp_collider(n = 500, p = 10, collider_strength = 2.0)
  expect_false(all(r1$X[, 1] == r2$X[, 1]))
})

test_that("dgp_nonlinear_outcome produces valid data", {
  set.seed(42)
  result <- dgp_nonlinear_outcome(n = 200, p = 10, nonlinearity = "moderate")
  expect_equal(length(result$Y), 200)
  expect_equal(ncol(result$X), 10)
  expect_equal(length(result$T), 200)
  expect_true(all(result$T %in% c(0, 1)))
  expect_true(is.numeric(result$true_ate))
  expect_equal(result$dgp_name, "nonlinear_outcome")
})

test_that("dgp_nonlinear_outcome moderate and severe produce different outcomes", {
  set.seed(42)
  r_mod <- dgp_nonlinear_outcome(n = 1000, p = 10, nonlinearity = "moderate")
  set.seed(42)
  r_sev <- dgp_nonlinear_outcome(n = 1000, p = 10, nonlinearity = "severe")
  # Different nonlinearity levels produce different Y
  expect_false(all(r_mod$Y == r_sev$Y))
  # Both have valid structure
  expect_equal(r_mod$dgp_name, "nonlinear_outcome")
  expect_equal(r_sev$dgp_name, "nonlinear_outcome")
})

test_that("dgp_nonlinear_outcome creates confounding that biases OLS", {
  set.seed(123)
  r <- dgp_nonlinear_outcome(n = 5000, p = 10, nonlinearity = "moderate")
  # OLS should be biased because it omits h(X) which drives both PS and Y
  fit <- lm(r$Y ~ r$T + r$X)
  ate_ols <- coef(fit)["r$T"]
  bias <- abs(ate_ols - r$true_ate)
  # With nonlinear confounding, OLS bias should be substantial (> 0.3)
  expect_gt(bias, 0.3)
})

test_that("dgp_nonlinear_propensity produces valid data", {
  set.seed(42)
  result <- dgp_nonlinear_propensity(n = 200, p = 10, nonlinearity = "moderate")
  expect_equal(length(result$Y), 200)
  expect_equal(ncol(result$X), 10)
  expect_true(all(result$T %in% c(0, 1)))
  expect_true(is.numeric(result$true_ate))
  expect_equal(result$dgp_name, "nonlinear_propensity")
})

test_that("dgp_nonlinear_propensity has linear outcome recoverable by OLS", {
  set.seed(42)
  r <- dgp_nonlinear_propensity(n = 5000, p = 10, nonlinearity = "moderate")
  fit <- lm(r$Y ~ r$T + r$X)
  ate_ols <- coef(fit)["r$T"]
  expect_lt(abs(ate_ols - r$true_ate), 0.3)
})

test_that("dgp_double_nonlinear produces valid data", {
  set.seed(42)
  result <- dgp_double_nonlinear(n = 200, p = 10)
  expect_equal(length(result$Y), 200)
  expect_equal(ncol(result$X), 10)
  expect_true(all(result$T %in% c(0, 1)))
  expect_equal(result$dgp_name, "double_nonlinear")
})

test_that("dgp_missing_data generates correct missing pattern", {
  result <- dgp_missing_data(n = 200, p = 10, missing_rate = 0.2, seed = 123)

  # Check structure
  expect_true("X_complete" %in% names(result))
  expect_true("X" %in% names(result))
  expect_equal(dim(result$X), dim(result$X_complete))

  # Check missing values exist in X but not in X_complete
  expect_true(any(is.na(result$X)))
  expect_false(any(is.na(result$X_complete)))

  # Check approximate missing rate (within tolerance)
  actual_rate <- mean(is.na(result$X))
  expect_true(abs(actual_rate - 0.2) < 0.1)
})
