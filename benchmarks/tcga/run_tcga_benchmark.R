#!/usr/bin/env Rscript

# benchmarks/tcga/run_tcga_benchmark.R
# Run benchmarks on TCGA semi-synthetic data
#
# Usage:
#   Rscript benchmarks/tcga/run_tcga_benchmark.R              # Run BRCA benchmark
#   Rscript benchmarks/tcga/run_tcga_benchmark.R --smoke      # Quick smoke test
#   Rscript benchmarks/tcga/run_tcga_benchmark.R --all        # Run all projects
#   Rscript benchmarks/tcga/run_tcga_benchmark.R --project TCGA-LUAD

# Determine script directory robustly
get_script_dir <- function() {
  # Method 1: commandArgs (works with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }

  # Method 2: Check if benchmarks/tcga/preprocess_tcga.R exists relative to cwd
  if (file.exists("benchmarks/tcga/preprocess_tcga.R")) {
    return("benchmarks/tcga")
  }

  # Method 3: Check if preprocess_tcga.R exists in current directory
  if (file.exists("preprocess_tcga.R")) {
    return(".")
  }

  stop("Cannot determine script directory. Run from project root or benchmarks/tcga/ directory.")
}

script_dir <- get_script_dir()

# Source preprocessing functions
source(file.path(script_dir, "preprocess_tcga.R"))

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Run TCGA semi-synthetic benchmark
#'
#' @param project TCGA project identifier (e.g., "TCGA-BRCA")
#' @param methods Character vector of cfomics methods to benchmark
#' @param n_reps Number of replications (different random seeds)
#' @param dgp_params List of parameters for semi-synthetic data generation
#' @param n_genes Number of genes to use (after variance filtering)
#' @param output_dir Directory to save results
#' @param n_workers Number of parallel workers (0 = auto-detect)
#' @param formula_k Number of covariates to include in formula
#'
#' @return Data frame with benchmark results
run_tcga_benchmark <- function(
    project = "TCGA-BRCA",
    methods = c("gformula", "hdml", "hdps", "tmle"),
    n_reps = 20,
    dgp_params = list(),
    n_genes = 500,
    output_dir = NULL,
    n_workers = 0L,
    formula_k = 10L
) {
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- file.path(script_dir, "results")
  }

  # Load TCGA data
  se_path <- file.path(script_dir, "data", paste0(project, "_se.rds"))
  if (!file.exists(se_path)) {
    stop(sprintf(
      "TCGA data not found at %s.\nRun download_tcga.R first to download %s data.",
      se_path, project
    ))
  }

  message(sprintf("Loading %s data from %s...", project, se_path))
  se <- readRDS(se_path)

  # Preprocess
  message(sprintf("Preprocessing: selecting top %d genes by variance...", n_genes))

  tcga_data <- preprocess_tcga(se, n_genes = n_genes)

  message(sprintf("TCGA %s: %d samples, %d genes",
                  project, nrow(tcga_data$X), ncol(tcga_data$X)))

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  raw_dir <- file.path(output_dir, "raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  # Build job list
  jobs <- list()
  idx <- 1L
  for (rep in seq_len(n_reps)) {
    for (method in methods) {
      rds_path <- file.path(raw_dir, sprintf("%s__%s__rep%03d.rds", project, method, rep))
      jobs[[idx]] <- list(
        project = project,
        method = method,
        rep = rep,
        rds_path = rds_path
      )
      idx <- idx + 1L
    }
  }

  total_jobs <- length(jobs)
  message(sprintf("Total jobs: %d (%d methods x %d reps)",
                  total_jobs, length(methods), n_reps))

  # Filter completed jobs
  pending_jobs <- Filter(function(j) !file.exists(j$rds_path), jobs)
  message(sprintf("Completed: %d, Pending: %d",
                  total_jobs - length(pending_jobs), length(pending_jobs)))

  if (length(pending_jobs) == 0) {
    message("All jobs complete. Loading existing results...")
    results <- lapply(jobs, function(j) readRDS(j$rds_path))
    results_df <- do.call(rbind, results)
    return(results_df)
  }

  # Load cfomics package
  if (!requireNamespace("cfomics", quietly = TRUE)) {
    stop("Package 'cfomics' required. Install with: devtools::install('packages/cfomics')")
  }

  # Set up parallel execution
  if (n_workers <= 0L) {
    n_workers <- max(1L, parallel::detectCores() - 1L)
  }

  message(sprintf("Running %d jobs on %d workers...", length(pending_jobs), n_workers))

  # Run jobs (sequential for simplicity; can be parallelized with future/furrr if needed)
  results <- list()
  for (i in seq_along(pending_jobs)) {
    job <- pending_jobs[[i]]
    message(sprintf("  [%d/%d] %s / %s / rep%d",
                    i, length(pending_jobs), job$project, job$method, job$rep))

    result <- run_single_tcga_job(
      job = job,
      tcga_data = tcga_data,
      dgp_params = dgp_params,
      formula_k = formula_k
    )

    # Atomic write
    tmp_path <- paste0(job$rds_path, ".tmp")
    saveRDS(result, tmp_path)
    file.rename(tmp_path, job$rds_path)

    results[[i]] <- result
  }

  # Load all results (including previously completed)
  all_results <- lapply(jobs, function(j) readRDS(j$rds_path))
  results_df <- do.call(rbind, all_results)

  # Save combined results
  output_file <- file.path(output_dir, paste0(project, "_results.rds"))
  saveRDS(results_df, output_file)
  message(sprintf("Results saved to %s", output_file))

  # Print summary
  print_tcga_summary(results_df)

  results_df
}


#' Run a single TCGA benchmark job
#'
#' @param job List with project, method, rep, rds_path
#' @param tcga_data Preprocessed TCGA data from preprocess_tcga()
#' @param dgp_params Parameters for semi-synthetic data generation
#' @param formula_k Number of covariates in formula
#'
#' @return Data frame with single row of results
run_single_tcga_job <- function(job, tcga_data, dgp_params, formula_k) {
  # Use robust seed: base_seed * multiplier + rep offset for better randomness
  set.seed(20260206 * 1000 + job$rep)

  # Create semi-synthetic data
  sim_data <- create_semi_synthetic_data(tcga_data, dgp_params)

  result <- tryCatch({
    # Build formula with PC columns
    cov_names <- grep("^PC", names(sim_data$data), value = TRUE)
    if (length(cov_names) == 0) {
      # Fall back to X columns if no PC columns
      cov_names <- grep("^X", names(sim_data$data), value = TRUE)
    }
    k <- min(formula_k, length(cov_names))
    fml <- stats::as.formula(paste("Y ~ T |", paste(cov_names[1:k], collapse = " + ")))

    # Fit model
    fit_time <- system.time({
      fit <- cfomics::cf_fit(fml, data = sim_data$data, method = job$method)
    })

    # Get predictions
    ate_hat <- predict(fit, type = "ate")
    ite_hat <- predict(fit, type = "ite")
    summary_hat <- tryCatch(predict(fit, type = "summary"), error = function(e) NULL)

    # Compute metrics
    metrics <- cfomics:::cf_benchmark_compute_metrics(
      ate_hat = ate_hat,
      ite_hat = ite_hat,
      summary_hat = summary_hat,
      truth = sim_data$truth
    )

    data.frame(
      project = job$project,
      method = job$method,
      rep = job$rep,
      n = nrow(sim_data$data),
      p = length(cov_names),
      ate_true = sim_data$truth$ate_true,
      ate_hat = ate_hat,
      bias_ate = metrics$bias_ate,
      abs_bias_ate = metrics$abs_bias_ate,
      squared_error_ate = metrics$squared_error_ate,
      pehe = metrics$pehe,
      coverage_ate = metrics$coverage_ate,
      ci_len_ate = metrics$ci_len_ate,
      time_sec = as.numeric(fit_time["elapsed"]),
      status = "ok",
      error_msg = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      project = job$project,
      method = job$method,
      rep = job$rep,
      n = NA_integer_,
      p = NA_integer_,
      ate_true = NA_real_,
      ate_hat = NA_real_,
      bias_ate = NA_real_,
      abs_bias_ate = NA_real_,
      squared_error_ate = NA_real_,
      pehe = NA_real_,
      coverage_ate = NA_real_,
      ci_len_ate = NA_real_,
      time_sec = NA_real_,
      status = "error",
      error_msg = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  result
}


#' Print summary of TCGA benchmark results
#'
#' @param results_df Data frame with benchmark results
print_tcga_summary <- function(results_df) {
  message("\n=== TCGA Benchmark Summary ===")

  # Status counts
  status_counts <- table(results_df$status)
  message(sprintf("Total: %d, OK: %d, Errors: %d",
                  nrow(results_df),
                  status_counts["ok"] %||% 0,
                  status_counts["error"] %||% 0))

  # Per-method summary (for OK results only)
  ok_results <- results_df[results_df$status == "ok", ]
  if (nrow(ok_results) > 0) {
    message("\nPer-method results:")
    for (method in unique(ok_results$method)) {
      m_data <- ok_results[ok_results$method == method, ]
      message(sprintf("  %s: Bias=%.3f (%.3f), PEHE=%.3f, Time=%.1fs",
                      method,
                      mean(m_data$bias_ate, na.rm = TRUE),
                      sd(m_data$bias_ate, na.rm = TRUE),
                      mean(m_data$pehe, na.rm = TRUE),
                      mean(m_data$time_sec, na.rm = TRUE)))
    }
  }

  # Report errors if any
  err_results <- results_df[results_df$status == "error", ]
  if (nrow(err_results) > 0) {
    message("\nErrors:")
    for (i in seq_len(min(5, nrow(err_results)))) {
      r <- err_results[i, ]
      message(sprintf("  - %s / rep%d: %s", r$method, r$rep, r$error_msg))
    }
    if (nrow(err_results) > 5) {
      message(sprintf("  ... and %d more errors", nrow(err_results) - 5))
    }
  }
}


#' TCGA benchmark configuration (similar to simulation benchmark config)
tcga_benchmark_config <- function() {
  list(
    projects = c("TCGA-BRCA"),
    methods = c("gformula", "hdml", "hdps", "tmle"),
    n_reps = 20L,
    n_genes = 500L,
    formula_k = 10L,
    dgp_params = list(
      tau = 2.0,
      heterogeneous = FALSE,
      noise_sd = 1.0
    )
  )
}

#' Smoke test configuration for quick testing
tcga_benchmark_config_smoke <- function() {
  cfg <- tcga_benchmark_config()
  cfg$n_reps <- 2L
  cfg$methods <- c("gformula")
  cfg
}


# Main execution
if (sys.nframe() == 0) {
  library(cfomics)

  # Parse command-line arguments
  args <- commandArgs(trailingOnly = TRUE)

  # Handle --help flag
  if ("--help" %in% args) {
    cat("Usage: Rscript run_tcga_benchmark.R [options]\n")
    cat("\n")
    cat("Options:\n")
    cat("  --smoke           Run quick smoke test (2 reps, gformula only)\n")
    cat("  --all             Run all TCGA projects (BRCA, LUAD, COAD)\n")
    cat("  --project NAME    Run specific project (e.g., TCGA-BRCA)\n")
    cat("  --help            Show this help message and exit\n")
    cat("\n")
    cat("Examples:\n")
    cat("  Rscript benchmarks/tcga/run_tcga_benchmark.R              # Run BRCA\n")
    cat("  Rscript benchmarks/tcga/run_tcga_benchmark.R --smoke      # Quick test\n")
    cat("  Rscript benchmarks/tcga/run_tcga_benchmark.R --all        # All projects\n")
    cat("  Rscript benchmarks/tcga/run_tcga_benchmark.R --project TCGA-LUAD\n")
    quit(status = 0, save = "no")
  }

  smoke <- "--smoke" %in% args
  run_all <- "--all" %in% args

  # Parse --project argument
  project_idx <- which(args == "--project")
  custom_project <- NULL
  if (length(project_idx) > 0 && project_idx < length(args)) {
    custom_project <- args[project_idx + 1]
  }

  # Select configuration
  if (smoke) {
    message("Running smoke test...")
    cfg <- tcga_benchmark_config_smoke()
  } else {
    message("Running TCGA benchmark...")
    cfg <- tcga_benchmark_config()
  }

  # Determine projects to run
  if (!is.null(custom_project)) {
    projects <- custom_project
  } else if (run_all) {
    projects <- c("TCGA-BRCA", "TCGA-LUAD", "TCGA-COAD")
  } else {
    projects <- cfg$projects
  }

  # Display configuration
  message(sprintf("  Projects: %s", paste(projects, collapse = ", ")))
  message(sprintf("  Methods: %s", paste(cfg$methods, collapse = ", ")))
  message(sprintf("  Replications: %d", cfg$n_reps))
  message(sprintf("  Genes: %d", cfg$n_genes))
  message("")

  # Run benchmarks for each project
  all_results <- list()
  for (proj in projects) {
    message(sprintf("\n=== Running %s benchmark ===", proj))
    result <- tryCatch({
      run_tcga_benchmark(
        project = proj,
        methods = cfg$methods,
        n_reps = cfg$n_reps,
        dgp_params = cfg$dgp_params,
        n_genes = cfg$n_genes,
        formula_k = cfg$formula_k
      )
    }, error = function(e) {
      message(sprintf("Error running %s: %s", proj, e$message))
      NULL
    })
    if (!is.null(result)) {
      all_results[[proj]] <- result
    }
  }

  # Final summary
  message("\n=== All Benchmarks Complete ===")
  for (proj in names(all_results)) {
    r <- all_results[[proj]]
    n_ok <- sum(r$status == "ok")
    n_err <- sum(r$status == "error")
    message(sprintf("  %s: %d OK, %d errors", proj, n_ok, n_err))
  }
}
