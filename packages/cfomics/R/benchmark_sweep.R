#' Sweep Experiment Functions for Systematic Benchmarking
#'
#' These functions generate benchmark datasets across parameter ranges
#' for systematic evaluation of causal inference methods.
#'
#' @name benchmark_sweep
NULL

#' Run dimension sweep experiment
#'
#' Generates datasets with varying n and p to test method performance
#' under different sample size and dimensionality regimes.
#'
#' @param n_values Integer vector, sample sizes to test
#' @param p_values Integer vector, number of covariates to test
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each n-p combination
#' @export
cf_benchmark_dimension_sweep <- function(
    n_values = c(100, 200, 500, 1000),
    p_values = c(50, 100, 500, 1000),
    n_reps = 10,
    seed = 123
) {
  results <- list()
  idx <- 1

  for (n in n_values) {
    for (p in p_values) {
      for (rep in seq_len(n_reps)) {
        dgp <- dgp_dimension_sweep(
          n = n,
          p = p,
          seed = seed + idx
        )
        dgp$sweep_params <- list(
          n = n,
          p = p,
          n_p_ratio = n / p,
          rep = rep
        )
        results[[idx]] <- dgp
        idx <- idx + 1
      }
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "dimension",
            n_values = n_values,
            p_values = p_values)
}

#' Run heterogeneity sweep experiment
#'
#' Generates datasets with varying treatment effect heterogeneity strength.
#'
#' @param strength_values Numeric vector, heterogeneity strength values to test
#' @param n Integer, sample size
#' @param p Integer, number of covariates
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each strength value
#' @export
cf_benchmark_heterogeneity_sweep <- function(
    strength_values = c(0, 0.5, 1.0, 2.0, 5.0),
    n = 500,
    p = 100,
    n_reps = 10,
    seed = 123
) {
  results <- list()
  idx <- 1

  for (strength in strength_values) {
    for (rep in seq_len(n_reps)) {
      dgp <- dgp_heterogeneous_linear(
        n = n,
        p = p,
        strength = strength,
        seed = seed + idx
      )
      dgp$sweep_params <- list(
        strength = strength,
        rep = rep
      )
      results[[idx]] <- dgp
      idx <- idx + 1
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "heterogeneity",
            strength_values = strength_values)
}

#' Run nonlinearity sweep experiment
#'
#' Generates datasets with varying nonlinearity strength in confounding.
#'
#' @param strength_values Numeric vector, nonlinearity strength values (0=linear, 1=full)
#' @param nonlinear_type Character, type of nonlinearity
#' @param n Integer, sample size
#' @param p Integer, number of covariates
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each strength value
#' @export
cf_benchmark_nonlinearity_sweep <- function(
    strength_values = c(0, 0.25, 0.5, 0.75, 1.0),
    nonlinear_type = "combined",
    n = 500,
    p = 100,
    n_reps = 10,
    seed = 123
) {
  results <- list()
  idx <- 1

  for (strength in strength_values) {
    for (rep in seq_len(n_reps)) {
      dgp <- dgp_nonlinear_confounding(
        n = n,
        p = p,
        nonlinear_type = nonlinear_type,
        strength = strength,
        seed = seed + idx
      )
      dgp$sweep_params <- list(
        strength = strength,
        nonlinear_type = nonlinear_type,
        rep = rep
      )
      results[[idx]] <- dgp
      idx <- idx + 1
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "nonlinearity",
            strength_values = strength_values)
}

#' Run confounder density sweep experiment
#'
#' Generates datasets with varying number of confounders.
#'
#' @param n_confounders_values Integer vector, number of confounders to test
#' @param coef_scaling Character, coefficient scaling method
#' @param n Integer, sample size
#' @param p Integer, number of covariates (must be >= max(n_confounders_values))
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each confounder count
#' @export
cf_benchmark_density_sweep <- function(
    n_confounders_values = c(10, 50, 100, 200, 500),
    coef_scaling = "sqrt",
    n = 500,
    p = 500,
    n_reps = 10,
    seed = 123
) {
  results <- list()
  idx <- 1

  for (n_conf in n_confounders_values) {
    for (rep in seq_len(n_reps)) {
      dgp <- dgp_dense_confounding(
        n = n,
        p = p,
        n_confounders = n_conf,
        coef_scaling = coef_scaling,
        seed = seed + idx
      )
      dgp$sweep_params <- list(
        n_confounders = n_conf,
        coef_scaling = coef_scaling,
        rep = rep
      )
      results[[idx]] <- dgp
      idx <- idx + 1
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "density",
            n_confounders_values = n_confounders_values)
}

#' Run overlap sweep experiment
#'
#' Generates datasets with varying propensity score overlap.
#'
#' @param overlap_values Character vector, overlap strength categories
#' @param n Integer, sample size
#' @param p Integer, number of covariates
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each overlap setting
#' @export
cf_benchmark_overlap_sweep <- function(
    overlap_values = c("good", "moderate", "weak", "extreme"),
    n = 500,
    p = 100,
    n_reps = 10,
    seed = 123
) {
  results <- list()
  idx <- 1

  for (overlap in overlap_values) {
    for (rep in seq_len(n_reps)) {
      dgp <- dgp_weak_overlap(
        n = n,
        p = p,
        overlap_strength = overlap,
        seed = seed + idx
      )
      dgp$sweep_params <- list(
        overlap_strength = overlap,
        rep = rep
      )
      results[[idx]] <- dgp
      idx <- idx + 1
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "overlap",
            overlap_values = overlap_values)
}

#' Run covariate shift sweep experiment
#'
#' Generates datasets with varying covariate distribution shifts.
#'
#' @param shift_magnitudes Numeric vector, shift magnitude values to test
#' @param shift_type Character, type of shift ("mean", "variance")
#' @param n Integer, sample size
#' @param p Integer, number of covariates
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each shift magnitude
#' @export
cf_benchmark_covariate_shift_sweep <- function(
    shift_magnitudes = c(0, 0.5, 1.0, 2.0, 3.0),
    shift_type = "mean",
    n = 500,
    p = 100,
    n_reps = 10,
    seed = 123
) {
  results <- list()
  idx <- 1

  for (magnitude in shift_magnitudes) {
    for (rep in seq_len(n_reps)) {
      dgp <- dgp_covariate_shift(
        n = n,
        p = p,
        shift_type = shift_type,
        shift_magnitude = magnitude,
        seed = seed + idx
      )
      dgp$sweep_params <- list(
        shift_magnitude = magnitude,
        shift_type = shift_type,
        rep = rep
      )
      results[[idx]] <- dgp
      idx <- idx + 1
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "covariate_shift",
            shift_magnitudes = shift_magnitudes)
}

#' Run correlation sweep experiment
#'
#' Generates datasets with varying confounder correlation strength.
#'
#' @param correlation_values Numeric vector, correlation strength values (0-1)
#' @param correlation_type Character, correlation structure type
#' @param n Integer, sample size
#' @param p Integer, number of covariates
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each correlation strength
#' @export
cf_benchmark_correlation_sweep <- function(
    correlation_values = c(0, 0.25, 0.5, 0.75, 0.9),
    correlation_type = "block",
    n = 500,
    p = 100,
    n_reps = 10,
    seed = 123
) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for correlation sweep. Please install it.")
  }

  results <- list()
  idx <- 1

  for (corr in correlation_values) {
    for (rep in seq_len(n_reps)) {
      dgp <- dgp_correlated_confounding(
        n = n,
        p = p,
        correlation_type = correlation_type,
        correlation_strength = corr,
        seed = seed + idx
      )
      dgp$sweep_params <- list(
        correlation_strength = corr,
        correlation_type = correlation_type,
        rep = rep
      )
      results[[idx]] <- dgp
      idx <- idx + 1
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "correlation",
            correlation_values = correlation_values)
}

#' Run unmeasured confounding sweep experiment
#'
#' Generates datasets with varying strength of unobserved confounding.
#'
#' @param strength_values Numeric vector, unobserved confounding strength values
#' @param n_unobserved Integer, number of unobserved confounders
#' @param n Integer, sample size
#' @param p Integer, number of observed covariates
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each strength value
#' @export
cf_benchmark_unmeasured_sweep <- function(
    strength_values = c(0, 0.25, 0.5, 1.0, 2.0),
    n_unobserved = 5,
    n = 500,
    p = 100,
    n_reps = 10,
    seed = 123
) {
  results <- list()
  idx <- 1

  for (strength in strength_values) {
    for (rep in seq_len(n_reps)) {
      dgp <- dgp_unobserved_confounding(
        n = n,
        p = p,
        n_unobserved = n_unobserved,
        unobserved_strength = strength,
        seed = seed + idx
      )
      dgp$sweep_params <- list(
        unobserved_strength = strength,
        n_unobserved = n_unobserved,
        rep = rep
      )
      results[[idx]] <- dgp
      idx <- idx + 1
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "unmeasured",
            strength_values = strength_values)
}

#' Run collider sweep experiment
#'
#' Generates datasets with varying collider strength.
#'
#' @param strength_values Numeric vector, collider strength values
#' @param n Integer, sample size
#' @param p Integer, number of covariates
#' @param n_reps Integer, number of replications per setting
#' @param seed Integer, base random seed
#' @return List of DGP results for each strength value
#' @export
cf_benchmark_collider_sweep <- function(
    strength_values = c(0, 0.25, 0.5, 1.0, 2.0),
    n = 500,
    p = 100,
    n_reps = 10,
    seed = 123
) {
  results <- list()
  idx <- 1

  for (strength in strength_values) {
    for (rep in seq_len(n_reps)) {
      dgp <- dgp_collider(
        n = n,
        p = p,
        collider_strength = strength,
        seed = seed + idx
      )
      dgp$sweep_params <- list(
        collider_strength = strength,
        rep = rep
      )
      results[[idx]] <- dgp
      idx <- idx + 1
    }
  }

  structure(results, class = c("cf_sweep_result", "list"),
            sweep_type = "collider",
            strength_values = strength_values)
}

#' Run benchmark on sweep results
#'
#' Applies causal inference methods to sweep experiment results.
#'
#' @param sweep_result cf_sweep_result object from a sweep function
#' @param methods Character vector, methods to evaluate
#' @param formula_k Integer, number of covariates to include in formula
#' @param ... Additional arguments passed to cf_fit
#' @return data.frame with benchmark results
#' @export
cf_benchmark_run_sweep <- function(
    sweep_result,
    methods = c("gformula", "hdml"),
    formula_k = 10,
    ...
) {
  if (!inherits(sweep_result, "cf_sweep_result")) {
    stop("sweep_result must be a cf_sweep_result object from a sweep function")
  }

  results_list <- list()
  idx <- 1
  total <- length(sweep_result) * length(methods)

  message(sprintf("Running %d benchmark evaluations...", total))

  for (i in seq_along(sweep_result)) {
    dgp <- sweep_result[[i]]
    p <- ncol(dgp$X)
    k <- min(formula_k, p)

    cov_names <- paste0("X", 1:k)
    fml <- stats::as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

    data <- data.frame(
      Y = dgp$Y,
      T = dgp$T,
      dgp$X
    )

    for (method in methods) {
      result_row <- tryCatch({
        fit_time <- system.time({
          fit <- cf_fit(fml, data = data, method = method, ...)
        })

        ate_hat <- predict(fit, type = "ate")
        ite_hat <- predict(fit, type = "ite")

        data.frame(
          sweep_type = attr(sweep_result, "sweep_type"),
          method = method,
          ate_true = dgp$true_ate,
          ate_hat = ate_hat,
          bias = ate_hat - dgp$true_ate,
          mse_ate = (ate_hat - dgp$true_ate)^2,
          pehe = sqrt(mean((ite_hat - dgp$true_ite)^2)),
          time_sec = as.numeric(fit_time["elapsed"]),
          status = "ok",
          error = NA_character_,
          dgp$sweep_params,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(
          sweep_type = attr(sweep_result, "sweep_type"),
          method = method,
          ate_true = dgp$true_ate,
          ate_hat = NA_real_,
          bias = NA_real_,
          mse_ate = NA_real_,
          pehe = NA_real_,
          time_sec = NA_real_,
          status = "error",
          error = conditionMessage(e),
          dgp$sweep_params,
          stringsAsFactors = FALSE
        )
      })

      results_list[[idx]] <- result_row
      idx <- idx + 1
    }

    if (idx %% 20 == 0) {
      message(sprintf("  Progress: %d/%d", idx - 1, total))
    }
  }

  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  message("Sweep benchmark complete.")
  results_df
}

#' Summarize sweep benchmark results
#'
#' Aggregates sweep benchmark results by sweep parameter and method.
#'
#' @param results_df data.frame from cf_benchmark_run_sweep
#' @param group_by Character vector, columns to group by (in addition to method)
#' @return data.frame with summary statistics
#' @export
cf_benchmark_summarize_sweep <- function(results_df, group_by = NULL) {
  sweep_type <- unique(results_df$sweep_type)[1]

  # Determine grouping columns based on sweep type
  if (is.null(group_by)) {
    group_by <- switch(
      sweep_type,
      "dimension" = c("n", "p"),
      "heterogeneity" = "strength",
      "nonlinearity" = "strength",
      "density" = "n_confounders",
      "overlap" = "overlap_strength",
      "covariate_shift" = "shift_magnitude",
      "correlation" = "correlation_strength",
      "unmeasured" = "unobserved_strength",
      "collider" = "collider_strength",
      character(0)
    )
  }

  group_cols <- c(group_by, "method")

  # Create grouping key
  ok_df <- results_df[results_df$status == "ok", ]

  if (nrow(ok_df) == 0) {
    warning("No successful runs to summarize")
    return(NULL)
  }

  # Aggregate
  summary_list <- list()
  unique_groups <- unique(ok_df[, group_cols, drop = FALSE])

  for (i in seq_len(nrow(unique_groups))) {
    group_vals <- unique_groups[i, , drop = FALSE]

    mask <- rep(TRUE, nrow(ok_df))
    for (col in group_cols) {
      mask <- mask & (ok_df[[col]] == group_vals[[col]])
    }

    subset_df <- ok_df[mask, ]

    if (nrow(subset_df) > 0) {
      summary_row <- data.frame(
        group_vals,
        n_reps = nrow(subset_df),
        rmse_ate = sqrt(mean(subset_df$mse_ate, na.rm = TRUE)),
        mean_bias = mean(subset_df$bias, na.rm = TRUE),
        mean_pehe = mean(subset_df$pehe, na.rm = TRUE),
        median_time_sec = stats::median(subset_df$time_sec, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      summary_list[[i]] <- summary_row
    }
  }

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  summary_df
}
