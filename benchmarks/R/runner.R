# benchmarks/R/runner.R
# Parallel execution engine with interrupt/resume support

library(future)
library(furrr)
library(progressr)

# Build filename for a single job result
make_result_filename <- function(scenario_id, method, rep, results_dir) {
  fname <- sprintf("%s__%s__rep%03d.rds", scenario_id, method, rep)
  file.path(results_dir, "raw", fname)
}

# Run a single benchmark job
run_single_job <- function(scenario, method, rep, cfg) {
  suppressPackageStartupMessages(library(cfomics))

  gen <- generate_scenario_data(scenario, rep, cfg$base_seed)

  k <- min(cfg$formula_k, gen$p)
  cov_names <- paste0("X", seq_len(k))
  fml <- stats::as.formula(paste("Y ~ T |", paste(cov_names, collapse = " + ")))

  # Method-specific args
  extra_args <- list()
  if (method == "bcf") {
    extra_args <- list(n_burn = cfg$bcf_n_burn, n_iter = cfg$bcf_n_iter)
  }

  result <- tryCatch({
    fit_time <- system.time({
      fit <- do.call(cf_fit, c(list(formula = fml, data = gen$data, method = method), extra_args))
    })

    ate_hat <- predict(fit, type = "ate")
    ite_hat <- predict(fit, type = "ite")

    summary_hat <- tryCatch(predict(fit, type = "summary"), error = function(e) NULL)

    metrics <- cfomics:::cf_benchmark_compute_metrics(
      ate_hat = ate_hat,
      ite_hat = ite_hat,
      summary_hat = summary_hat,
      truth = gen$truth
    )

    data.frame(
      scenario_id = gen$scenario_id,
      method = method,
      rep = gen$rep,
      n = gen$n,
      p = gen$p,
      seed = gen$seed,
      ate_true = gen$truth$ate_true,
      ate_hat = ate_hat,
      bias_ate = metrics$bias_ate,
      abs_bias_ate = metrics$abs_bias_ate,
      mse_ate = metrics$mse_ate,
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
      scenario_id = gen$scenario_id,
      method = method,
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
      time_sec = NA_real_,
      status = "error",
      error_msg = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  result
}

# Main runner: build job list, skip completed, run in parallel
run_benchmark <- function(cfg, retry_errors = FALSE) {
  results_dir <- cfg$results_dir
  dir.create(file.path(results_dir, "raw"), recursive = TRUE, showWarnings = FALSE)

  # Build full job list
  jobs <- list()
  idx <- 1L
  for (scenario in cfg$scenarios) {
    for (method in cfg$methods) {
      for (rep in seq_len(cfg$n_reps)) {
        rds_path <- make_result_filename(scenario$id, method, rep, results_dir)
        jobs[[idx]] <- list(
          scenario = scenario,
          method = method,
          rep = rep,
          rds_path = rds_path
        )
        idx <- idx + 1L
      }
    }
  }

  total <- length(jobs)
  message(sprintf("Total jobs: %d", total))

  # Filter: skip completed (unless retry_errors)
  existing <- vapply(jobs, function(j) {
    if (!file.exists(j$rds_path)) return(FALSE)
    if (retry_errors) {
      res <- readRDS(j$rds_path)
      return(res$status == "ok")
    }
    TRUE
  }, logical(1))

  pending_jobs <- jobs[!existing]
  message(sprintf("Completed: %d, Pending: %d", sum(existing), length(pending_jobs)))

  if (length(pending_jobs) == 0) {
    message("All jobs complete.")
    return(list(
      total = total,
      completed = sum(existing),
      pending = 0L,
      executed = 0L,
      results = list()
    ))
  }

  # Set up parallel
  n_workers <- cfg$n_workers
  if (n_workers <= 0L) n_workers <- max(1L, parallel::detectCores() - 1L)
  plan(multisession, workers = n_workers)
  on.exit(plan(sequential), add = TRUE)

  message(sprintf("Running %d jobs on %d workers...", length(pending_jobs), n_workers))

  # Execute in parallel with progressr
  results <- progressr::with_progress({
    p <- progressr::progressor(along = pending_jobs)
    furrr::future_map(pending_jobs, function(job) {
      source("benchmarks/R/scenarios.R", local = TRUE)
      res <- run_single_job(job$scenario, job$method, job$rep, cfg)
      saveRDS(res, job$rds_path)
      p()
      res
    }, .options = furrr::furrr_options(seed = TRUE))
  })

  message(sprintf("Done. %d jobs executed.", length(results)))

  # Return structured summary
  list(
    total = total,
    completed = sum(existing),
    pending = length(pending_jobs),
    executed = length(results),
    results = results
  )
}
