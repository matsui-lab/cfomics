# benchmarks/config.R
# Centralized benchmark configuration

benchmark_config <- function() {
  # Compute absolute path to scenarios.R for use in parallel workers
  # (workers may have different working directories)
  scenarios_path <- NULL

  # Try to find scenarios.R from current working directory or common locations
  candidates <- c(
    file.path("benchmarks", "R", "scenarios.R"),
    file.path("..", "benchmarks", "R", "scenarios.R"),
    file.path(Sys.getenv("CFOMICS_ROOT", "."), "benchmarks", "R", "scenarios.R")
  )
  for (cand in candidates) {
    if (file.exists(cand)) {
      scenarios_path <- normalizePath(cand, mustWork = TRUE)
      break
    }
  }

  if (is.null(scenarios_path)) {
    warning("Could not locate scenarios.R - parallel workers may fail. ",
            "Set CFOMICS_ROOT environment variable to the project root.")
    # Default fallback (will be validated at runtime)
    scenarios_path <- file.path("benchmarks", "R", "scenarios.R")
  }

  list(
    # Absolute path to scenarios.R for parallel workers
    scenarios_path = scenarios_path,

    # Methods to benchmark
    methods = c("gformula", "hdml", "hdps", "bcf", "tmle"),

    # Number of replications per scenario setting
    n_reps = 50L,

    # Base random seed (each job gets seed = base_seed + job_index)
    base_seed = 20260128L,

    # Parallel workers (0 = auto-detect)
    n_workers = 0L,

    # BCF MCMC settings
    bcf_n_burn = 1000L,
    bcf_n_iter = 2000L,

    # Formula: how many covariates to include
    formula_k = 10L,

    # Output directory (normalized to absolute path)
    results_dir = normalizePath(file.path("benchmarks", "results"), mustWork = FALSE),

    # Scenario definitions
    scenarios = list(
      # S1: Baseline
      list(id = "S1_n500_p50",   dgp = "dgp_baseline", n = 500L,  p = 50L,  params = list()),
      list(id = "S1_n500_p100",  dgp = "dgp_baseline", n = 500L,  p = 100L, params = list()),
      list(id = "S1_n1000_p50",  dgp = "dgp_baseline", n = 1000L, p = 50L,  params = list()),
      list(id = "S1_n1000_p100", dgp = "dgp_baseline", n = 1000L, p = 100L, params = list()),
      list(id = "S1_n2000_p50",  dgp = "dgp_baseline", n = 2000L, p = 50L,  params = list()),
      list(id = "S1_n2000_p100", dgp = "dgp_baseline", n = 2000L, p = 100L, params = list()),

      # S2: Dimension sweep
      list(id = "S2_n500_p100",  dgp = "dgp_dimension_sweep", n = 500L,  p = 100L,  params = list()),
      list(id = "S2_n500_p500",  dgp = "dgp_dimension_sweep", n = 500L,  p = 500L,  params = list()),
      list(id = "S2_n500_p1000", dgp = "dgp_dimension_sweep", n = 500L,  p = 1000L, params = list()),
      list(id = "S2_n1000_p100", dgp = "dgp_dimension_sweep", n = 1000L, p = 100L,  params = list()),
      list(id = "S2_n1000_p500", dgp = "dgp_dimension_sweep", n = 1000L, p = 500L,  params = list()),
      list(id = "S2_n1000_p1000",dgp = "dgp_dimension_sweep", n = 1000L, p = 1000L, params = list()),

      # S3: Heterogeneous linear (strength sweep)
      list(id = "S3_str0",   dgp = "dgp_heterogeneous_linear", n = 1000L, p = 100L, params = list(strength = 0)),
      list(id = "S3_str0.5", dgp = "dgp_heterogeneous_linear", n = 1000L, p = 100L, params = list(strength = 0.5)),
      list(id = "S3_str1",   dgp = "dgp_heterogeneous_linear", n = 1000L, p = 100L, params = list(strength = 1.0)),
      list(id = "S3_str2",   dgp = "dgp_heterogeneous_linear", n = 1000L, p = 100L, params = list(strength = 2.0)),

      # S4a-c: Other heterogeneity types
      list(id = "S4a_nonlinear",    dgp = "dgp_heterogeneous_nonlinear",    n = 1000L, p = 100L, params = list()),
      list(id = "S4b_subgroup",     dgp = "dgp_heterogeneous_subgroup",     n = 1000L, p = 100L, params = list()),
      list(id = "S4c_qualitative",  dgp = "dgp_heterogeneous_qualitative",  n = 1000L, p = 100L, params = list()),

      # S5: Nonlinear confounding (type sweep)
      list(id = "S5_quadratic",    dgp = "dgp_nonlinear_confounding", n = 1000L, p = 100L, params = list(nonlinear_type = "quadratic")),
      list(id = "S5_trigonometric",dgp = "dgp_nonlinear_confounding", n = 1000L, p = 100L, params = list(nonlinear_type = "trigonometric")),
      list(id = "S5_interaction",  dgp = "dgp_nonlinear_confounding", n = 1000L, p = 100L, params = list(nonlinear_type = "interaction")),
      list(id = "S5_combined",     dgp = "dgp_nonlinear_confounding", n = 1000L, p = 100L, params = list(nonlinear_type = "combined")),
      list(id = "S5_threshold",    dgp = "dgp_nonlinear_confounding", n = 1000L, p = 100L, params = list(nonlinear_type = "threshold")),

      # S6: Dense confounding (n_confounders sweep)
      list(id = "S6_conf10",  dgp = "dgp_dense_confounding", n = 1000L, p = 500L, params = list(n_confounders = 10L)),
      list(id = "S6_conf50",  dgp = "dgp_dense_confounding", n = 1000L, p = 500L, params = list(n_confounders = 50L)),
      list(id = "S6_conf100", dgp = "dgp_dense_confounding", n = 1000L, p = 500L, params = list(n_confounders = 100L)),
      list(id = "S6_conf200", dgp = "dgp_dense_confounding", n = 1000L, p = 500L, params = list(n_confounders = 200L)),

      # S7: Weak overlap (strength sweep)
      list(id = "S7_good",     dgp = "dgp_weak_overlap", n = 1000L, p = 100L, params = list(overlap_strength = "good")),
      list(id = "S7_moderate", dgp = "dgp_weak_overlap", n = 1000L, p = 100L, params = list(overlap_strength = "moderate")),
      list(id = "S7_weak",     dgp = "dgp_weak_overlap", n = 1000L, p = 100L, params = list(overlap_strength = "weak")),
      list(id = "S7_extreme",  dgp = "dgp_weak_overlap", n = 1000L, p = 100L, params = list(overlap_strength = "extreme")),

      # S8: Covariate shift (type sweep)
      list(id = "S8_mean",     dgp = "dgp_covariate_shift", n = 1000L, p = 100L, params = list(shift_type = "mean")),
      list(id = "S8_variance", dgp = "dgp_covariate_shift", n = 1000L, p = 100L, params = list(shift_type = "variance")),

      # S9: Correlated confounding (type sweep)
      list(id = "S9_block",  dgp = "dgp_correlated_confounding", n = 1000L, p = 100L, params = list(correlation_type = "block")),
      list(id = "S9_ar1",    dgp = "dgp_correlated_confounding", n = 1000L, p = 100L, params = list(correlation_type = "ar1")),
      list(id = "S9_factor", dgp = "dgp_correlated_confounding", n = 1000L, p = 100L, params = list(correlation_type = "factor")),

      # S10: Unobserved confounding (strength sweep)
      list(id = "S10_str0",   dgp = "dgp_unobserved_confounding", n = 1000L, p = 100L, params = list(unobserved_strength = 0)),
      list(id = "S10_str0.5", dgp = "dgp_unobserved_confounding", n = 1000L, p = 100L, params = list(unobserved_strength = 0.5)),
      list(id = "S10_str1",   dgp = "dgp_unobserved_confounding", n = 1000L, p = 100L, params = list(unobserved_strength = 1.0)),
      list(id = "S10_str2",   dgp = "dgp_unobserved_confounding", n = 1000L, p = 100L, params = list(unobserved_strength = 2.0)),

      # S11: Collider (strength sweep)
      list(id = "S11_str0",   dgp = "dgp_collider", n = 1000L, p = 100L, params = list(collider_strength = 0)),
      list(id = "S11_str0.5", dgp = "dgp_collider", n = 1000L, p = 100L, params = list(collider_strength = 0.5)),
      list(id = "S11_str1",   dgp = "dgp_collider", n = 1000L, p = 100L, params = list(collider_strength = 1.0)),
      list(id = "S11_str2",   dgp = "dgp_collider", n = 1000L, p = 100L, params = list(collider_strength = 2.0))
    )
  )
}

# Quick config for smoke testing (3 scenarios, 2 reps)
benchmark_config_smoke <- function() {
  cfg <- benchmark_config()
  cfg$n_reps <- 2L
  cfg$scenarios <- cfg$scenarios[c(1, 13, 30)]  # S1_n500_p50, S3_str0, S7_good
  cfg$bcf_n_burn <- 100L
  cfg$bcf_n_iter <- 200L
  cfg
}
