#' Run a single benchmark experiment
#'
#' This function runs a single benchmark experiment for a given scenario and
#' method combination, returning a tidy one-row data.frame with results.
#'
#' @param scenario Character, the DGP scenario name
#' @param method Character, the causal inference method to use
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed_dgp Integer, random seed for data generation
#' @param seed_method Integer, random seed for method (default: seed_dgp)
#' @param formula Formula for cf_fit (default: Y ~ T | X1 + X2)
#' @param effect_size Numeric, base treatment effect size
#' @param noise_sd Numeric, standard deviation of outcome noise
#' @param allow_python_env_install Logical, whether to allow Python environment
#'   installation (default: FALSE). When FALSE, Python methods will be skipped
#'   if dependencies are not available.
#' @param ... Additional arguments passed to cf_fit
#'
#' @return A one-row data.frame with columns:
#'   scenario_id, method, replicate_id, n, p, seed_dgp, seed_method,
#'   ate_true, ate_hat, bias_ate, squared_error_ate, pehe, coverage_ate, ci_len_ate,
#'   time_fit_sec, time_predict_sec, status, error_message
#'
#' @export
cf_benchmark_run_once <- function(
    scenario,
    method,
    n,
    p,
    seed_dgp,
    seed_method = seed_dgp,
    formula = stats::as.formula("Y ~ T | X1 + X2"),
    effect_size = 1.0,
    noise_sd = 1.0,
    allow_python_env_install = FALSE,
    ...
) {
  result_row <- data.frame(
    scenario_id = scenario,
    method = method,
    replicate_id = seed_dgp,
    n = n,
    p = p,
    seed_dgp = seed_dgp,
    seed_method = seed_method,
    ate_true = NA_real_,
    ate_hat = NA_real_,
    bias_ate = NA_real_,
    squared_error_ate = NA_real_,
    pehe = NA_real_,
    coverage_ate = NA_real_,
    ci_len_ate = NA_real_,
    time_fit_sec = NA_real_,
    time_predict_sec = NA_real_,
    status = "pending",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  
  skip_reason <- .check_method_dependencies(method, allow_python_env_install)
  if (!is.null(skip_reason)) {
    result_row$status <- "skipped"
    result_row$error_message <- skip_reason
    return(result_row)
  }
  
  tryCatch({
    dgp_result <- cf_benchmark_generate_data(
      scenario = scenario,
      n = n,
      p = p,
      seed = seed_dgp,
      effect_size = effect_size,
      noise_sd = noise_sd,
      return_graph = (method == "dowhy_gcm")
    )
    
    result_row$ate_true <- dgp_result$truth$ate_true
    
    graph <- dgp_result$graph
    
    time_fit <- system.time({
      fit <- cf_fit(
        formula = formula,
        data = dgp_result$data,
        method = method,
        graph = graph,
        random_state = as.integer(seed_method),
        ...
      )
    })
    result_row$time_fit_sec <- as.numeric(time_fit["elapsed"])
    
    time_predict <- system.time({
      ate_hat <- predict(fit, type = "ate")
      ite_hat <- predict(fit, type = "ite")
      summary_hat <- tryCatch(
        predict(fit, type = "summary"),
        error = function(e) NULL
      )
    })
    result_row$time_predict_sec <- as.numeric(time_predict["elapsed"])
    
    result_row$ate_hat <- as.numeric(ate_hat)
    
    metrics <- cf_benchmark_compute_metrics(
      ate_hat = ate_hat,
      ite_hat = ite_hat,
      summary_hat = summary_hat,
      truth = dgp_result$truth
    )
    
    result_row$bias_ate <- metrics$bias_ate
    result_row$squared_error_ate <- metrics$squared_error_ate
    result_row$pehe <- metrics$pehe
    result_row$coverage_ate <- metrics$coverage_ate
    result_row$ci_len_ate <- metrics$ci_len_ate
    result_row$status <- "ok"
    
  }, error = function(e) {
    result_row$status <<- "error"
    result_row$error_message <<- conditionMessage(e)
  })
  
  result_row
}

#' Check if method dependencies are available
#'
#' @param method Character, the method name
#' @param allow_python_env_install Logical, whether to allow Python env install
#' @return NULL if dependencies are available, or a character string with skip reason
#' @keywords internal
.check_method_dependencies <- function(method, allow_python_env_install) {
  r_only_methods <- c("grf", "ipw", "gformula")
  
  if (method %in% r_only_methods) {
    if (method == "grf" && !requireNamespace("grf", quietly = TRUE)) {
      return("Package 'grf' not installed")
    }
    if (method == "ipw") {
      if (!requireNamespace("ipw", quietly = TRUE)) {
        return("Package 'ipw' not installed")
      }
      if (!requireNamespace("survey", quietly = TRUE)) {
        return("Package 'survey' not installed")
      }
    }
    return(NULL)
  }
  
  python_methods <- c("dowhy_gcm", "drlearner", "cavae", "ganite")
  
  if (method %in% python_methods) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      return("Package 'reticulate' not installed")
    }
    
    if (!allow_python_env_install) {
      py_available <- tryCatch(
        reticulate::py_available(initialize = TRUE),
        error = function(e) FALSE
      )
      
      if (!py_available) {
        return("Python not available (allow_python_env_install=FALSE)")
      }
      
      module_check <- .check_python_modules(method)
      if (!is.null(module_check)) {
        return(module_check)
      }
    }
  }
  
  NULL
}

#' Check if required Python modules are available for a method
#'
#' @param method Character, the method name
#' @return NULL if modules are available, or a character string with skip reason
#' @keywords internal
.check_python_modules <- function(method) {
  tryCatch({
    switch(
      method,
      "dowhy_gcm" = {
        if (!reticulate::py_module_available("dowhy")) {
          return("Python module 'dowhy' not available")
        }
      },
      "drlearner" = {
        if (!reticulate::py_module_available("sklearn")) {
          return("Python module 'sklearn' not available")
        }
        econml_ok <- tryCatch({
          reticulate::import("econml", convert = FALSE)
          TRUE
        }, error = function(e) FALSE)
        if (!econml_ok) {
          return("Python module 'econml' not available")
        }
      },
      "cavae" = {
        torch_ok <- tryCatch({
          torch <- reticulate::import("torch", convert = FALSE)
          TRUE
        }, error = function(e) FALSE)
        if (!torch_ok) {
          return("PyTorch not available")
        }
        if (!reticulate::py_module_available("pyro")) {
          return("Python module 'pyro' not available")
        }
      },
      "ganite" = {
        if (!reticulate::py_module_available("tensorflow")) {
          return("Python module 'tensorflow' not available")
        }
      }
    )
    NULL
  }, error = function(e) {
    paste0("Error checking Python modules: ", conditionMessage(e))
  })
}

#' Run batch benchmark experiments
#'
#' This function runs benchmark experiments across multiple scenarios, methods,
#' and replications, returning a tidy data.frame with all results.
#'
#' @param scenarios Character vector of scenario names
#' @param methods Character vector of method names
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param n_reps Integer, number of replications per scenario-method combination
#' @param base_seed Integer, base random seed (seeds will be base_seed + 1:n_reps)
#' @param effect_size Numeric, base treatment effect size
#' @param noise_sd Numeric, standard deviation of outcome noise
#' @param formula_template Character, template for formula (default: "Y ~ T | \%s")
#' @param covariates_k Integer, number of covariates to include in formula
#' @param allow_python_env_install Logical, whether to allow Python environment
#'   installation (default: FALSE)
#' @param ... Additional arguments passed to cf_fit
#'
#' @return A data.frame with one row per scenario-method-replication combination
#'
#' @export
cf_benchmark_run <- function(
    scenarios,
    methods,
    n = 500L,
    p = 20L,
    n_reps = 20L,
    base_seed = 1L,
    effect_size = 1.0,
    noise_sd = 1.0,
    formula_template = "Y ~ T | %s",
    covariates_k = 10L,
    allow_python_env_install = FALSE,
    ...
) {
  n <- as.integer(n)
  p <- as.integer(p)
  n_reps <- as.integer(n_reps)
  base_seed <- as.integer(base_seed)
  covariates_k <- as.integer(covariates_k)
  
  k <- min(covariates_k, p)
  covariate_names <- paste0("X", seq_len(k))
  covariate_str <- paste(covariate_names, collapse = " + ")
  formula <- stats::as.formula(sprintf(formula_template, covariate_str))
  
  results_list <- list()
  idx <- 1L
  
  total_runs <- length(scenarios) * length(methods) * n_reps
  message(sprintf("Running %d benchmark experiments...", total_runs))
  
  for (scenario in scenarios) {
    for (method in methods) {
      for (rep_i in seq_len(n_reps)) {
        seed_dgp <- base_seed + rep_i
        
        result <- cf_benchmark_run_once(
          scenario = scenario,
          method = method,
          n = n,
          p = p,
          seed_dgp = seed_dgp,
          seed_method = seed_dgp,
          formula = formula,
          effect_size = effect_size,
          noise_sd = noise_sd,
          allow_python_env_install = allow_python_env_install,
          ...
        )
        
        results_list[[idx]] <- result
        idx <- idx + 1L
        
        if (idx %% 10 == 0 || idx == total_runs) {
          message(sprintf("  Completed %d/%d runs", idx - 1, total_runs))
        }
      }
    }
  }
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  message("Benchmark complete.")
  results_df
}

#' Summarize benchmark results
#'
#' This function aggregates benchmark results by scenario and method,
#' computing summary statistics across replications.
#'
#' @param results_df data.frame from cf_benchmark_run()
#'
#' @return A data.frame with one row per scenario-method combination containing:
#'   rmse_ate, mean_bias_ate, mean_pehe, coverage_rate, median_time_fit_sec,
#'   median_time_predict_sec, n_ok, n_skipped, n_error
#'
#' @export
cf_benchmark_summarize <- function(results_df) {
  scenarios <- unique(results_df$scenario_id)
  methods <- unique(results_df$method)
  
  summary_list <- list()
  idx <- 1L
  
  for (scenario in scenarios) {
    for (method in methods) {
      subset_df <- results_df[
        results_df$scenario_id == scenario & results_df$method == method,
      ]
      
      ok_df <- subset_df[subset_df$status == "ok", ]
      
      n_ok <- sum(subset_df$status == "ok", na.rm = TRUE)
      n_skipped <- sum(subset_df$status == "skipped", na.rm = TRUE)
      n_error <- sum(subset_df$status == "error", na.rm = TRUE)
      
      if (n_ok > 0) {
        rmse_ate <- sqrt(mean(ok_df$squared_error_ate, na.rm = TRUE))
        mean_bias_ate <- mean(ok_df$bias_ate, na.rm = TRUE)
        mean_pehe <- mean(ok_df$pehe, na.rm = TRUE)
        coverage_rate <- mean(ok_df$coverage_ate, na.rm = TRUE)
        median_time_fit_sec <- median(ok_df$time_fit_sec, na.rm = TRUE)
        median_time_predict_sec <- median(ok_df$time_predict_sec, na.rm = TRUE)
      } else {
        rmse_ate <- NA_real_
        mean_bias_ate <- NA_real_
        mean_pehe <- NA_real_
        coverage_rate <- NA_real_
        median_time_fit_sec <- NA_real_
        median_time_predict_sec <- NA_real_
      }
      
      summary_list[[idx]] <- data.frame(
        scenario_id = scenario,
        method = method,
        rmse_ate = rmse_ate,
        mean_bias_ate = mean_bias_ate,
        mean_pehe = mean_pehe,
        coverage_rate = coverage_rate,
        median_time_fit_sec = median_time_fit_sec,
        median_time_predict_sec = median_time_predict_sec,
        n_ok = n_ok,
        n_skipped = n_skipped,
        n_error = n_error,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  summary_df
}
