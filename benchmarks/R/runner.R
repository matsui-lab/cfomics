# benchmarks/R/runner.R
# Parallel execution engine with interrupt/resume support

library(future)
library(furrr)
library(progressr)

# Null coalescing operator (used for cfg$external_methods)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Build filename for a single job result
make_result_filename <- function(scenario_id, method, rep, results_dir) {
  fname <- sprintf("%s__%s__rep%03d.rds", scenario_id, method, rep)
  file.path(results_dir, "raw", fname)
}

# Run a single benchmark job
# Note: cfomics package is loaded via furrr_options(packages = "cfomics")
# All cfomics functions use explicit namespace prefix cfomics:: for clarity
run_single_job <- function(scenario, method, rep, cfg) {
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
    gc(reset = TRUE)
    mem_before <- sum(gc()[, 2])

    fit_time <- system.time({
      fit <- do.call(cfomics::cf_fit, c(list(formula = fml, data = gen$data, method = method), extra_args))
    })

    gc()
    mem_after <- sum(gc()[, 2])
    peak_memory_mb <- max(0, mem_after - mem_before)

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
      squared_error_ate = metrics$squared_error_ate,
      pehe = metrics$pehe,
      coverage_ate = metrics$coverage_ate,
      ci_len_ate = metrics$ci_len_ate,
      time_sec = as.numeric(fit_time["elapsed"]),
      peak_memory_mb = peak_memory_mb,
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
      squared_error_ate = NA_real_,
      pehe = NA_real_,
      coverage_ate = NA_real_,
      ci_len_ate = NA_real_,
      time_sec = NA_real_,
      peak_memory_mb = NA_real_,
      status = "error",
      error_msg = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  result
}

# Run a single benchmark job for EXTERNAL method
run_single_job_external <- function(scenario, method, rep, cfg) {
  gen <- generate_scenario_data(scenario, rep, cfg$base_seed)

  # Source external wrappers
  external_dir <- file.path(dirname(cfg$scenarios_path), "..", "external")
  source(file.path(external_dir, "registry.R"), local = TRUE)
  source(file.path(external_dir, "wrappers", paste0("wrapper_", method, ".R")), local = TRUE)

  result <- tryCatch({
    gc(reset = TRUE)
    mem_before <- sum(gc()[, 2])

    ext_result <- run_external_method(
      method = method,
      X = gen$data[, grep("^X", names(gen$data))],
      T = gen$data$T,
      Y = gen$data$Y
    )

    gc()
    mem_after <- sum(gc()[, 2])
    peak_memory_mb <- max(0, mem_after - mem_before)

    # Compute metrics
    metrics <- cfomics:::cf_benchmark_compute_metrics(
      ate_hat = ext_result$ate,
      ite_hat = ext_result$ite,
      summary_hat = list(
        ate_ci_lower = ext_result$ci_lower,
        ate_ci_upper = ext_result$ci_upper
      ),
      truth = gen$truth
    )

    data.frame(
      scenario_id = gen$scenario_id,
      method = paste0("ext_", method),
      rep = gen$rep,
      n = gen$n,
      p = gen$p,
      seed = gen$seed,
      ate_true = gen$truth$ate_true,
      ate_hat = ext_result$ate,
      bias_ate = metrics$bias_ate,
      abs_bias_ate = metrics$abs_bias_ate,
      squared_error_ate = metrics$squared_error_ate,
      pehe = metrics$pehe,
      coverage_ate = metrics$coverage_ate,
      ci_len_ate = metrics$ci_len_ate,
      time_sec = ext_result$time_sec,
      peak_memory_mb = peak_memory_mb,
      status = "ok",
      error_msg = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      scenario_id = gen$scenario_id,
      method = paste0("ext_", method),
      rep = gen$rep,
      n = gen$n,
      p = gen$p,
      seed = gen$seed,
      ate_true = gen$truth$ate_true,
      ate_hat = NA_real_,
      bias_ate = NA_real_,
      abs_bias_ate = NA_real_,
      squared_error_ate = NA_real_,
      pehe = NA_real_,
      coverage_ate = NA_real_,
      ci_len_ate = NA_real_,
      time_sec = NA_real_,
      peak_memory_mb = NA_real_,
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

  # Build full job list (cfomics methods)
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
          rds_path = rds_path,
          is_external = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  # Add external methods to job list
  external_methods <- cfg$external_methods %||% character(0)
  for (scenario in cfg$scenarios) {
    for (method in external_methods) {
      for (rep in seq_len(cfg$n_reps)) {
        # External methods use "ext_" prefix in filename
        rds_path <- make_result_filename(scenario$id, paste0("ext_", method), rep, results_dir)
        jobs[[idx]] <- list(
          scenario = scenario,
          method = method,
          rep = rep,
          rds_path = rds_path,
          is_external = TRUE
        )
        idx <- idx + 1L
      }
    }
  }

  total <- length(jobs)
  message(sprintf("Total jobs: %d (cfomics: %d, external: %d)",
                  total, length(cfg$methods) * length(cfg$scenarios) * cfg$n_reps,
                  length(external_methods) * length(cfg$scenarios) * cfg$n_reps))

  # Filter: skip completed (unless retry_errors)
  existing <- vapply(jobs, function(j) {
    if (!file.exists(j$rds_path)) return(FALSE)
    if (retry_errors) {
      # Handle corrupted RDS files gracefully - treat as incomplete
      res <- tryCatch(readRDS(j$rds_path), error = function(e) NULL)
      if (is.null(res)) return(FALSE)
      return(isTRUE(res$status == "ok"))
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
  # n_workers <= 0 means auto-detect: use all cores minus 1 (the -1L reserves
  # one core for the main R process and OS operations to keep the system responsive)
  if (n_workers <= 0L) n_workers <- max(1L, parallel::detectCores() - 1L)
  # Capture existing plan to restore on exit (before modifying)
  old_plan <- plan()
  on.exit(plan(old_plan), add = TRUE)
  plan(multisession, workers = n_workers)

  message(sprintf("Running %d jobs on %d workers...", length(pending_jobs), n_workers))

  # Execute in parallel with progressr
  results <- progressr::with_progress({
    p <- progressr::progressor(along = pending_jobs)
    furrr::future_map(pending_jobs, function(job) {
      # Source scenarios.R using absolute path from config (workers may have different cwd)
      source(cfg$scenarios_path, local = TRUE)

      # Dispatch to appropriate runner based on method type
      if (isTRUE(job$is_external)) {
        res <- run_single_job_external(job$scenario, job$method, job$rep, cfg)
      } else {
        res <- run_single_job(job$scenario, job$method, job$rep, cfg)
      }

      # Atomic write: save to temp file first, then rename to avoid partial writes
      tmp_path <- paste0(job$rds_path, ".tmp")
      saveRDS(res, tmp_path)
      file.rename(tmp_path, job$rds_path)
      p()
      res
    }, .options = furrr::furrr_options(seed = TRUE, packages = "cfomics"))
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
