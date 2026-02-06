# benchmarks/external/nn/run_nn_benchmark.R
# R interface to run neural network benchmarks via reticulate
#
# These models are for PAPER COMPARISON ONLY (not part of cfomics package).
# Implements standard causal inference NN architectures:
#   - TARNet: Two-head architecture for Y(0) and Y(1)
#   - CFRNet: TARNet + IPM regularization
#   - DragonNet: TARNet + propensity score co-learning

#' Run a single neural network benchmark job
#'
#' Executes a causal inference neural network model (TARNet, CFRNet, or DragonNet)
#' on a benchmark scenario and computes evaluation metrics.
#'
#' @param scenario List containing scenario specification (id, n, p, dgp, params)
#' @param method Character, one of "tarnet", "cfrnet", or "dragonnet"
#' @param rep Integer, replication number
#' @param cfg List with benchmark configuration including:
#'   - scenarios_path: Path to scenarios.R
#'   - base_seed: Base random seed
#' @param nn_config List with neural network hyperparameters:
#'   - hidden_dim: Hidden layer dimension (default: 200)
#'   - n_layers: Number of representation layers (default: 3)
#'   - n_epochs: Training epochs (default: 300)
#'   - lr: Learning rate (default: 1e-3)
#'   - batch_size: Mini-batch size (default: 64)
#'   - alpha: IPM weight for CFRNet (default: 1.0)
#'
#' @return Data frame with one row containing benchmark results
#'
#' @examples
#' \dontrun{
#' # Load benchmark configuration
#' source("benchmarks/config.R")
#' cfg <- get_benchmark_config()
#'
#' # Run TARNet on first scenario
#' result <- run_nn_benchmark(
#'   scenario = cfg$scenarios[[1]],
#'   method = "tarnet",
#'   rep = 1,
#'   cfg = cfg
#' )
#' }
run_nn_benchmark <- function(scenario, method, rep, cfg, nn_config = list()) {
  # Check reticulate availability

if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' required for NN benchmarks")
  }

  # Source scenarios for data generation
  source(cfg$scenarios_path, local = TRUE)
  gen <- generate_scenario_data(scenario, rep, cfg$base_seed)

  # Extract data components
  X <- as.matrix(gen$data[, grep("^X", names(gen$data))])
  T <- gen$data$T
  Y <- gen$data$Y

  # Locate Python models
  nn_dir <- file.path(dirname(cfg$scenarios_path), "..", "external", "nn")
  models_path <- file.path(nn_dir, "models.py")

  if (!file.exists(models_path)) {
    stop("Neural network models not found at: ", models_path)
  }

  # Source Python module
  reticulate::source_python(models_path)

  # Select model class
  model_class <- switch(method,
    "tarnet" = TARNet,
    "cfrnet" = CFRNet,
    "dragonnet" = DragonNet,
    stop("Unknown NN method: ", method)
  )

  # Model hyperparameters with defaults
  input_dim <- ncol(X)
  hidden_dim <- nn_config$hidden_dim %||% 200L
  n_layers <- nn_config$n_layers %||% 3L
  n_epochs <- nn_config$n_epochs %||% 300L
  lr <- nn_config$lr %||% 1e-3
  batch_size <- nn_config$batch_size %||% 64L
  alpha <- nn_config$alpha %||% 1.0

  # Ensure integer types for Python
  hidden_dim <- as.integer(hidden_dim)
  n_layers <- as.integer(n_layers)
  n_epochs <- as.integer(n_epochs)
  batch_size <- as.integer(batch_size)

  # Time the training and prediction
  start_time <- Sys.time()

  result <- tryCatch({
    # Instantiate model
    if (method == "cfrnet") {
      model <- model_class(input_dim, hidden_dim, n_layers, alpha)
    } else {
      model <- model_class(input_dim, hidden_dim, n_layers)
    }

    # Train model
    model <- train_model(
      model = model,
      X = X,
      T = T,
      Y = Y,
      n_epochs = n_epochs,
      lr = lr,
      batch_size = batch_size,
      alpha = alpha,
      verbose = FALSE
    )

    # Predict ITE
    pred <- predict_ite(model, X)

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    # Compute metrics using cfomics internal function
    metrics <- cfomics:::cf_benchmark_compute_metrics(
      ate_hat = pred$ate,
      ite_hat = pred$ite,
      summary_hat = NULL,  # NN models don't provide CIs
      truth = gen$truth
    )

    data.frame(
      scenario_id = gen$scenario_id,
      method = paste0("nn_", method),
      rep = gen$rep,
      n = gen$n,
      p = gen$p,
      seed = gen$seed,
      ate_true = gen$truth$ate_true,
      ate_hat = pred$ate,
      bias_ate = metrics$bias_ate,
      abs_bias_ate = metrics$abs_bias_ate,
      mse_ate = metrics$mse_ate,
      pehe = metrics$pehe,
      coverage_ate = NA_real_,
      ci_len_ate = NA_real_,
      time_sec = elapsed,
      status = "ok",
      error_msg = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    data.frame(
      scenario_id = gen$scenario_id,
      method = paste0("nn_", method),
      rep = gen$rep,
      n = gen$n,
      p = gen$p,
      seed = gen$seed,
      ate_true = gen$truth$ate_true,
      ate_hat = NA_real_,
      bias_ate = NA_real_,
      abs_bias_ate = NA_real_,
      mse_ate = NA_real_,
      pehe = NA_real_,
      coverage_ate = NA_real_,
      ci_len_ate = NA_real_,
      time_sec = elapsed,
      status = "error",
      error_msg = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  result
}


#' Run all neural network methods on all scenarios
#'
#' Convenience function to run TARNet, CFRNet, and DragonNet on all
#' configured benchmark scenarios.
#'
#' @param cfg Benchmark configuration list
#' @param nn_methods Character vector of NN methods to run
#'   (default: all three)
#' @param nn_config Hyperparameters for NN models
#' @param save_results Logical, save individual results to RDS files
#'
#' @return Data frame with all benchmark results
#'
#' @examples
#' \dontrun{
#' source("benchmarks/config.R")
#' cfg <- get_benchmark_config()
#' results <- run_all_nn_benchmarks(cfg)
#' }
run_all_nn_benchmarks <- function(
    cfg,
    nn_methods = c("tarnet", "cfrnet", "dragonnet"),
    nn_config = list(),
    save_results = TRUE
) {
  results <- list()
  idx <- 1L

  results_dir <- cfg$results_dir
  if (save_results) {
    dir.create(file.path(results_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  }

  total_jobs <- length(cfg$scenarios) * length(nn_methods) * cfg$n_reps
  message(sprintf("Running %d NN benchmark jobs...", total_jobs))

  for (scenario in cfg$scenarios) {
    for (method in nn_methods) {
      for (rep in seq_len(cfg$n_reps)) {
        message(sprintf("[%d/%d] %s / %s / rep %d",
                        idx, total_jobs, scenario$id, method, rep))

        res <- run_nn_benchmark(
          scenario = scenario,
          method = method,
          rep = rep,
          cfg = cfg,
          nn_config = nn_config
        )

        if (save_results) {
          fname <- sprintf("%s__nn_%s__rep%03d.rds", scenario$id, method, rep)
          saveRDS(res, file.path(results_dir, "raw", fname))
        }

        results[[idx]] <- res
        idx <- idx + 1L
      }
    }
  }

  do.call(rbind, results)
}


#' Check if Python environment is ready for NN benchmarks
#'
#' Verifies that reticulate and required Python packages are available.
#'
#' @return Logical, TRUE if ready
#'
#' @examples
#' \dontrun{
#' if (check_nn_environment()) {
#'   results <- run_all_nn_benchmarks(cfg)
#' }
#' }
check_nn_environment <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    message("reticulate package not installed")
    return(FALSE)
  }

  if (!reticulate::py_available(initialize = TRUE)) {
    message("Python not available via reticulate")
    return(FALSE)
  }

  required_modules <- c("torch", "numpy")
  missing <- character(0)

  for (mod in required_modules) {
    if (!reticulate::py_module_available(mod)) {
      missing <- c(missing, mod)
    }
  }

  if (length(missing) > 0) {
    message("Missing Python modules: ", paste(missing, collapse = ", "))
    message("Install with: pip install ", paste(missing, collapse = " "))
    return(FALSE)
  }

  TRUE
}


# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
