# Benchmark Simulation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a production-scale benchmark simulation (12,500 jobs) comparing 5 causal inference methods across 13 DGPs, with publication-quality bilingual reports.

**Architecture:** Standalone `benchmarks/` directory at repo root. `config.R` defines all parameters. `R/runner.R` dispatches parallel jobs via future+furrr, saving one RDS per job for interrupt/resume. `aggregate_results.R` collects into summary CSV. `R/reporting.R` generates figures and LaTeX tables. `report/benchmark_report.Rmd` renders parameterized PDF (en/ja).

**Tech Stack:** R, cfomics package, future+furrr (parallel), ggplot2+patchwork (viz), kableExtra (LaTeX), PMCMRplus (stats tests), rmarkdown+bookdown (reports)

---

### Task 1: Install Dependencies and Create Directory Structure

**Files:**
- Create: `benchmarks/` directory tree
- Create: `benchmarks/.gitignore`

**Step 1: Install missing packages**

```r
install.packages(c("future", "furrr", "PMCMRplus"))
```

Run: `Rscript -e 'install.packages(c("future", "furrr", "PMCMRplus"), repos="https://cran.r-project.org")'`
Expected: packages install successfully

**Step 2: Create directory structure**

```bash
mkdir -p benchmarks/R
mkdir -p benchmarks/report/sections
mkdir -p benchmarks/report/figures
mkdir -p benchmarks/report/tables
mkdir -p benchmarks/results/raw
mkdir -p benchmarks/results/summary
mkdir -p benchmarks/results/plots
mkdir -p benchmarks/results/tables
```

**Step 3: Create .gitignore for results**

Create `benchmarks/.gitignore`:
```
results/raw/
results/summary/
results/plots/
results/tables/
report/figures/
report/tables/
*.pdf
```

**Step 4: Commit**

```bash
git add benchmarks/.gitignore
git commit -m "chore: create benchmarks directory structure"
```

---

### Task 2: Create config.R

**Files:**
- Create: `benchmarks/config.R`

**Step 1: Write config.R**

```r
# benchmarks/config.R
# Centralized benchmark configuration

benchmark_config <- function() {
  list(
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

    # Output directory
    results_dir = file.path("benchmarks", "results"),

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
```

**Step 2: Verify config loads**

Run: `Rscript -e 'source("benchmarks/config.R"); cfg <- benchmark_config(); cat("Scenarios:", length(cfg$scenarios), "\nTotal jobs:", length(cfg$scenarios) * length(cfg$methods) * cfg$n_reps, "\n")'`
Expected: `Scenarios: 47` and `Total jobs: 11750`

**Step 3: Commit**

```bash
git add benchmarks/config.R
git commit -m "feat(bench): add centralized benchmark configuration"
```

---

### Task 3: Create R/scenarios.R

**Files:**
- Create: `benchmarks/R/scenarios.R`

**Step 1: Write scenarios.R**

This file provides `generate_scenario_data()` which takes a scenario config entry and rep index, calls the appropriate DGP function, and returns data in a format ready for cf_fit.

```r
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
```

**Step 2: Verify with smoke test**

Run:
```bash
Rscript -e '
  library(cfomics)
  source("benchmarks/config.R")
  source("benchmarks/R/scenarios.R")
  cfg <- benchmark_config()
  d <- generate_scenario_data(cfg$scenarios[[1]], rep = 1, base_seed = cfg$base_seed)
  cat("rows:", nrow(d$data), "cols:", ncol(d$data), "ate:", d$truth$ate_true, "\n")
'
```
Expected: `rows: 500 cols: 52 ate: <some number>`

**Step 3: Commit**

```bash
git add benchmarks/R/scenarios.R
git commit -m "feat(bench): add scenario data generation dispatcher"
```

---

### Task 4: Create R/runner.R

**Files:**
- Create: `benchmarks/R/runner.R`

**Step 1: Write runner.R**

```r
# benchmarks/R/runner.R
# Parallel execution engine with interrupt/resume support

library(future)
library(furrr)

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
    return(invisible(NULL))
  }

  # Set up parallel
  n_workers <- cfg$n_workers
  if (n_workers <= 0L) n_workers <- max(1L, parallel::detectCores() - 1L)
  plan(multisession, workers = n_workers)
  on.exit(plan(sequential), add = TRUE)

  message(sprintf("Running %d jobs on %d workers...", length(pending_jobs), n_workers))

  # Execute in parallel
  results <- furrr::future_map(pending_jobs, function(job) {
    source("benchmarks/R/scenarios.R", local = TRUE)
    res <- run_single_job(job$scenario, job$method, job$rep, cfg)
    saveRDS(res, job$rds_path)
    res
  }, .options = furrr::furrr_options(seed = TRUE),
     .progress = TRUE)

  message(sprintf("Done. %d jobs executed.", length(results)))
  invisible(results)
}
```

**Step 2: Smoke test with 2 reps, 1 scenario, 1 method**

Run:
```bash
Rscript -e '
  source("benchmarks/config.R")
  source("benchmarks/R/scenarios.R")
  source("benchmarks/R/runner.R")
  cfg <- benchmark_config_smoke()
  cfg$methods <- "gformula"
  cfg$scenarios <- cfg$scenarios[1]
  cfg$n_reps <- 2L
  run_benchmark(cfg)
  f <- list.files("benchmarks/results/raw", full.names = TRUE)
  cat("Result files:", length(f), "\n")
  cat(readRDS(f[1])$status, "\n")
'
```
Expected: `Result files: 2` and `ok`

**Step 3: Test interrupt/resume (re-run should skip)**

Run the same command again.
Expected: `Completed: 2, Pending: 0` and `All jobs complete.`

**Step 4: Commit**

```bash
git add benchmarks/R/runner.R
git commit -m "feat(bench): add parallel runner with interrupt/resume"
```

---

### Task 5: Create run_full_benchmark.R

**Files:**
- Create: `benchmarks/run_full_benchmark.R`

**Step 1: Write run_full_benchmark.R**

```r
#!/usr/bin/env Rscript
# benchmarks/run_full_benchmark.R
# Main entry point for production benchmark

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args
retry <- "--retry" %in% args

message("=== cfomics Full Benchmark ===")
message(sprintf("Mode: %s", if (smoke) "SMOKE TEST" else "PRODUCTION"))
message(sprintf("Time: %s", Sys.time()))

suppressPackageStartupMessages(library(cfomics))

source("benchmarks/config.R")
source("benchmarks/R/scenarios.R")
source("benchmarks/R/runner.R")

cfg <- if (smoke) benchmark_config_smoke() else benchmark_config()

message(sprintf("Scenarios: %d", length(cfg$scenarios)))
message(sprintf("Methods: %s", paste(cfg$methods, collapse = ", ")))
message(sprintf("Reps: %d", cfg$n_reps))
message(sprintf("Total jobs: %d", length(cfg$scenarios) * length(cfg$methods) * cfg$n_reps))

t0 <- proc.time()
run_benchmark(cfg, retry_errors = retry)
elapsed <- (proc.time() - t0)["elapsed"]

message(sprintf("\nTotal elapsed: %.1f minutes", elapsed / 60))
message("Done.")
```

**Step 2: Smoke test**

Run: `Rscript benchmarks/run_full_benchmark.R --smoke`
Expected: Runs ~30 jobs (3 scenarios x 5 methods x 2 reps) and completes within a few minutes.

**Step 3: Commit**

```bash
git add benchmarks/run_full_benchmark.R
git commit -m "feat(bench): add main benchmark entry point script"
```

---

### Task 6: Create aggregate_results.R

**Files:**
- Create: `benchmarks/aggregate_results.R`

**Step 1: Write aggregate_results.R**

```r
#!/usr/bin/env Rscript
# benchmarks/aggregate_results.R
# Collect raw RDS files into a single summary CSV

message("=== Aggregating Benchmark Results ===")

results_dir <- file.path("benchmarks", "results")
raw_dir <- file.path(results_dir, "raw")
summary_dir <- file.path(results_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Read all RDS files
rds_files <- list.files(raw_dir, pattern = "\\.rds$", full.names = TRUE)
message(sprintf("Found %d result files", length(rds_files)))

if (length(rds_files) == 0) {
  stop("No result files found. Run run_full_benchmark.R first.")
}

results_list <- lapply(rds_files, readRDS)
results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

# Add scenario group (S1, S2, ..., S11)
results_df$scenario_group <- sub("_.*", "", results_df$scenario_id)

# Save raw collected results
raw_csv <- file.path(summary_dir, "all_results.csv")
write.csv(results_df, raw_csv, row.names = FALSE)
message(sprintf("Saved raw results: %s (%d rows)", raw_csv, nrow(results_df)))

# Compute summary statistics per scenario_id x method
ok_df <- results_df[results_df$status == "ok", ]
n_errors <- sum(results_df$status == "error")
if (n_errors > 0) {
  message(sprintf("WARNING: %d jobs had errors", n_errors))
  err_df <- results_df[results_df$status == "error", ]
  err_summary <- file.path(summary_dir, "errors.csv")
  write.csv(err_df, err_summary, row.names = FALSE)
}

summary_df <- do.call(rbind, lapply(
  split(ok_df, list(ok_df$scenario_id, ok_df$method), drop = TRUE),
  function(df) {
    data.frame(
      scenario_id = df$scenario_id[1],
      scenario_group = df$scenario_group[1],
      method = df$method[1],
      n = df$n[1],
      p = df$p[1],
      n_reps = nrow(df),
      mean_bias = mean(df$bias_ate, na.rm = TRUE),
      sd_bias = sd(df$bias_ate, na.rm = TRUE),
      rmse_ate = sqrt(mean(df$mse_ate, na.rm = TRUE)),
      mean_pehe = mean(df$pehe, na.rm = TRUE),
      sd_pehe = sd(df$pehe, na.rm = TRUE),
      median_pehe = median(df$pehe, na.rm = TRUE),
      coverage = mean(df$coverage_ate, na.rm = TRUE),
      mean_ci_len = mean(df$ci_len_ate, na.rm = TRUE),
      median_time = median(df$time_sec, na.rm = TRUE),
      mean_time = mean(df$time_sec, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
rownames(summary_df) <- NULL

summary_csv <- file.path(summary_dir, "summary.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)
message(sprintf("Saved summary: %s (%d rows)", summary_csv, nrow(summary_df)))

message("Aggregation complete.")
```

**Step 2: Verify with smoke results**

Run: `Rscript benchmarks/aggregate_results.R`
Expected: Creates `benchmarks/results/summary/all_results.csv` and `summary.csv`

**Step 3: Commit**

```bash
git add benchmarks/aggregate_results.R
git commit -m "feat(bench): add result aggregation script"
```

---

### Task 7: Create R/metrics.R (Statistical Tests)

**Files:**
- Create: `benchmarks/R/metrics.R`

**Step 1: Write metrics.R**

```r
# benchmarks/R/metrics.R
# Statistical tests and advanced metrics for benchmark comparison

library(PMCMRplus)

# Friedman test + Nemenyi post-hoc per scenario group
run_friedman_tests <- function(results_df) {
  ok_df <- results_df[results_df$status == "ok", ]
  groups <- unique(ok_df$scenario_group)
  methods <- unique(ok_df$method)

  friedman_results <- list()

  for (grp in groups) {
    grp_df <- ok_df[ok_df$scenario_group == grp, ]

    # Build matrix: rows = reps (blocks), cols = methods, values = RMSE per rep
    # For Friedman we need matched data: same rep index across methods
    reps <- sort(unique(grp_df$rep))
    scenario_ids <- unique(grp_df$scenario_id)

    # Use MSE per rep as the performance measure
    mat_list <- list()
    for (sid in scenario_ids) {
      for (r in reps) {
        row <- list()
        for (m in methods) {
          val <- grp_df$mse_ate[grp_df$scenario_id == sid &
                                grp_df$rep == r &
                                grp_df$method == m]
          row[[m]] <- if (length(val) == 1) val else NA_real_
        }
        mat_list[[length(mat_list) + 1]] <- unlist(row)
      }
    }

    mat <- do.call(rbind, mat_list)
    mat <- mat[complete.cases(mat), , drop = FALSE]

    if (nrow(mat) < 3) next

    # Friedman test
    ft <- tryCatch({
      friedman.test(mat)
    }, error = function(e) NULL)

    # Nemenyi post-hoc
    nemenyi <- tryCatch({
      frdAllPairsNemenyiTest(mat)
    }, error = function(e) NULL)

    # Mean ranks
    ranks <- apply(mat, 1, rank)
    mean_ranks <- rowMeans(ranks)

    friedman_results[[grp]] <- list(
      scenario_group = grp,
      friedman_p = if (!is.null(ft)) ft$p.value else NA_real_,
      friedman_stat = if (!is.null(ft)) ft$statistic else NA_real_,
      mean_ranks = mean_ranks,
      nemenyi = nemenyi,
      n_blocks = nrow(mat)
    )
  }

  friedman_results
}

# Compute bias-variance decomposition
bias_variance_decomposition <- function(results_df) {
  ok_df <- results_df[results_df$status == "ok", ]

  do.call(rbind, lapply(
    split(ok_df, list(ok_df$scenario_id, ok_df$method), drop = TRUE),
    function(df) {
      bias <- mean(df$ate_hat - df$ate_true, na.rm = TRUE)
      variance <- var(df$ate_hat, na.rm = TRUE)
      data.frame(
        scenario_id = df$scenario_id[1],
        scenario_group = df$scenario_group[1],
        method = df$method[1],
        bias_sq = bias^2,
        variance = variance,
        mse = bias^2 + variance,
        stringsAsFactors = FALSE
      )
    }
  ))
}

# Compute overall method ranks across all scenarios
overall_method_ranks <- function(summary_df, metric = "rmse_ate") {
  # For each scenario, rank methods
  ranked <- do.call(rbind, lapply(
    split(summary_df, summary_df$scenario_id),
    function(df) {
      df$rank <- rank(df[[metric]])
      df[, c("scenario_id", "method", "rank")]
    }
  ))

  # Mean rank per method
  aggregate(rank ~ method, data = ranked, FUN = mean)
}
```

**Step 2: Commit**

```bash
git add benchmarks/R/metrics.R
git commit -m "feat(bench): add Friedman/Nemenyi tests and bias-variance decomposition"
```

---

### Task 8: Create R/reporting.R (Figures)

**Files:**
- Create: `benchmarks/R/reporting.R`

**Step 1: Write reporting.R**

```r
# benchmarks/R/reporting.R
# Publication-quality figure and table generation

library(ggplot2)
library(patchwork)
library(scales)
library(kableExtra)

# Theme for all plots
theme_benchmark <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank()
    )
}

# Method color palette (consistent across all figures)
method_colors <- c(
  "gformula" = "#E69F00",
  "hdml"     = "#56B4E9",
  "hdps"     = "#009E73",
  "bcf"      = "#F0E442",
  "tmle"     = "#CC79A7"
)

# Fig 2: Method comparison heatmap (RMSE)
plot_heatmap <- function(summary_df, metric = "rmse_ate", title = "RMSE (ATE)") {
  ggplot(summary_df, aes(x = method, y = scenario_id, fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", .data[[metric]])), size = 2.5) +
    scale_fill_viridis_c(option = "plasma", direction = -1) +
    labs(title = title, x = "Method", y = "Scenario", fill = metric) +
    theme_benchmark() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Fig 3: Bias boxplots
plot_bias_boxplot <- function(results_df, scenario_groups = NULL) {
  df <- results_df[results_df$status == "ok", ]
  if (!is.null(scenario_groups)) df <- df[df$scenario_group %in% scenario_groups, ]

  ggplot(df, aes(x = method, y = bias_ate, fill = method)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = method_colors) +
    facet_wrap(~scenario_id, scales = "free_y") +
    labs(title = "Bias Distribution", x = "Method", y = "Bias (ATE)") +
    theme_benchmark() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Fig 4: Coverage heatmap
plot_coverage_heatmap <- function(summary_df) {
  ggplot(summary_df, aes(x = method, y = scenario_id, fill = coverage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", coverage * 100)), size = 2.5) +
    scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                         midpoint = 0.95, limits = c(0, 1),
                         labels = percent) +
    labs(title = "95% CI Coverage", x = "Method", y = "Scenario", fill = "Coverage") +
    theme_benchmark() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Fig 6: Sweep curves
plot_sweep_curve <- function(summary_df, scenario_prefix, x_var, x_label,
                             metric = "rmse_ate", y_label = "RMSE (ATE)") {
  df <- summary_df[grepl(paste0("^", scenario_prefix), summary_df$scenario_id), ]
  if (nrow(df) == 0) return(NULL)

  # Extract sweep value from scenario_id
  df$sweep_val <- as.numeric(gsub(paste0(scenario_prefix, "_[a-z]*"), "", df$scenario_id))
  # Handle non-numeric sweep values
  if (all(is.na(df$sweep_val))) {
    df$sweep_val <- as.factor(gsub(paste0(scenario_prefix, "_"), "", df$scenario_id))
  }

  ggplot(df, aes(x = sweep_val, y = .data[[metric]], color = method, group = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = method_colors) +
    labs(title = paste("Sensitivity:", x_label), x = x_label, y = y_label) +
    theme_benchmark()
}

# Fig 7: Bias-Variance decomposition
plot_bias_variance <- function(bv_df, scenario_groups = NULL) {
  if (!is.null(scenario_groups)) bv_df <- bv_df[bv_df$scenario_group %in% scenario_groups, ]

  df_long <- data.frame(
    scenario_id = rep(bv_df$scenario_id, 2),
    method = rep(bv_df$method, 2),
    component = rep(c("Bias^2", "Variance"), each = nrow(bv_df)),
    value = c(bv_df$bias_sq, bv_df$variance),
    stringsAsFactors = FALSE
  )

  ggplot(df_long, aes(x = method, y = value, fill = component)) +
    geom_col(position = "stack") +
    facet_wrap(~scenario_id, scales = "free_y") +
    scale_fill_manual(values = c("Bias^2" = "#d73027", "Variance" = "#4575b4")) +
    labs(title = "Bias-Variance Decomposition", x = "Method", y = "MSE") +
    theme_benchmark() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Fig 8: Computation time
plot_time <- function(results_df) {
  df <- results_df[results_df$status == "ok", ]

  ggplot(df, aes(x = method, y = time_sec, fill = method)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_y_log10(labels = label_number(suffix = "s")) +
    scale_fill_manual(values = method_colors) +
    labs(title = "Computation Time", x = "Method", y = "Time (seconds, log scale)") +
    theme_benchmark()
}

# Save a plot as PDF + PNG
save_plot <- function(p, filename, out_dir, width = 10, height = 7) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(filename, ".pdf")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 300)
}

# Generate all figures
generate_all_figures <- function(results_df, summary_df, bv_df, out_dir) {
  message("Generating figures...")

  save_plot(plot_heatmap(summary_df, "rmse_ate", "RMSE (ATE)"),
            "fig2_rmse_heatmap", out_dir, width = 10, height = 14)

  save_plot(plot_heatmap(summary_df, "mean_pehe", "PEHE"),
            "fig2b_pehe_heatmap", out_dir, width = 10, height = 14)

  save_plot(plot_bias_boxplot(results_df),
            "fig3_bias_boxplot", out_dir, width = 16, height = 12)

  save_plot(plot_coverage_heatmap(summary_df),
            "fig4_coverage_heatmap", out_dir, width = 10, height = 14)

  # Fig 6: Sweep curves
  sweep_configs <- list(
    list(prefix = "S3",  x_var = "strength",           x_label = "Heterogeneity Strength"),
    list(prefix = "S7",  x_var = "overlap_strength",   x_label = "Overlap Strength"),
    list(prefix = "S10", x_var = "unobserved_strength", x_label = "Unobserved Confounding Strength"),
    list(prefix = "S11", x_var = "collider_strength",   x_label = "Collider Strength")
  )
  sweep_plots <- list()
  for (sc in sweep_configs) {
    p <- plot_sweep_curve(summary_df, sc$prefix, sc$x_var, sc$x_label)
    if (!is.null(p)) sweep_plots[[length(sweep_plots) + 1]] <- p
  }
  if (length(sweep_plots) > 0) {
    combined <- patchwork::wrap_plots(sweep_plots, ncol = 2)
    save_plot(combined, "fig6_sweep_curves", out_dir, width = 14, height = 10)
  }

  save_plot(plot_bias_variance(bv_df),
            "fig7_bias_variance", out_dir, width = 16, height = 12)

  save_plot(plot_time(results_df),
            "fig8_computation_time", out_dir, width = 8, height = 6)

  message("Figures saved to: ", out_dir)
}

# ============================================================================
# LaTeX Tables
# ============================================================================

# Tab 1: Main results table (bold=best, underline=2nd)
generate_main_table <- function(summary_df, out_dir) {
  metrics <- c("rmse_ate", "mean_pehe", "mean_bias")
  methods <- unique(summary_df$method)
  scenarios <- unique(summary_df$scenario_id)

  # For each scenario x metric, find best and 2nd best
  format_cell <- function(val, is_best, is_second) {
    s <- sprintf("%.3f", val)
    if (is_best) s <- paste0("\\textbf{", s, "}")
    if (is_second) s <- paste0("\\underline{", s, "}")
    s
  }

  for (metric in metrics) {
    tab_df <- do.call(rbind, lapply(scenarios, function(sid) {
      row <- summary_df[summary_df$scenario_id == sid, ]
      vals <- row[[metric]]
      sorted_idx <- order(abs(vals))
      best_method <- row$method[sorted_idx[1]]
      second_method <- if (length(sorted_idx) > 1) row$method[sorted_idx[2]] else NA

      cells <- vapply(seq_len(nrow(row)), function(i) {
        format_cell(vals[i],
                    row$method[i] == best_method,
                    !is.na(second_method) && row$method[i] == second_method)
      }, character(1))
      names(cells) <- row$method
      c(Scenario = sid, cells[methods])
    }))

    tab_df <- as.data.frame(tab_df, stringsAsFactors = FALSE)

    tex <- kableExtra::kbl(tab_df, format = "latex", booktabs = TRUE, escape = FALSE,
                           caption = paste("Benchmark Results:", metric)) |>
      kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(tex, file.path(out_dir, paste0("tab1_", metric, ".tex")))
  }
}

# Tab 2: Coverage table
generate_coverage_table <- function(summary_df, out_dir) {
  tab <- summary_df[, c("scenario_id", "method", "coverage", "mean_ci_len")]
  wide_cov <- reshape(tab[, c("scenario_id", "method", "coverage")],
                      idvar = "scenario_id", timevar = "method",
                      direction = "wide")

  tex <- kableExtra::kbl(wide_cov, format = "latex", booktabs = TRUE,
                         caption = "95\\% CI Coverage by Method and Scenario",
                         digits = 3) |>
    kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(tex, file.path(out_dir, "tab2_coverage.tex"))
}

# Tab 3: Computation time
generate_time_table <- function(summary_df, out_dir) {
  tab <- summary_df[, c("scenario_id", "method", "median_time")]
  wide <- reshape(tab, idvar = "scenario_id", timevar = "method", direction = "wide")

  tex <- kableExtra::kbl(wide, format = "latex", booktabs = TRUE,
                         caption = "Median Computation Time (seconds)",
                         digits = 2) |>
    kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(tex, file.path(out_dir, "tab3_time.tex"))
}

# Tab 4: Friedman test results
generate_friedman_table <- function(friedman_results, out_dir) {
  rows <- lapply(friedman_results, function(fr) {
    ranks_str <- paste(sprintf("%s=%.2f", names(fr$mean_ranks), fr$mean_ranks),
                       collapse = ", ")
    data.frame(
      Scenario = fr$scenario_group,
      Friedman_p = fr$friedman_p,
      n_blocks = fr$n_blocks,
      Mean_Ranks = ranks_str,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)

  tex <- kableExtra::kbl(tab, format = "latex", booktabs = TRUE,
                         caption = "Friedman Test Results and Mean Ranks",
                         digits = 4) |>
    kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(tex, file.path(out_dir, "tab4_friedman.tex"))
}

# Generate all tables
generate_all_tables <- function(summary_df, friedman_results, out_dir) {
  message("Generating LaTeX tables...")
  generate_main_table(summary_df, out_dir)
  generate_coverage_table(summary_df, out_dir)
  generate_time_table(summary_df, out_dir)
  generate_friedman_table(friedman_results, out_dir)
  message("Tables saved to: ", out_dir)
}
```

**Step 2: Commit**

```bash
git add benchmarks/R/reporting.R
git commit -m "feat(bench): add figure and LaTeX table generation"
```

---

### Task 9: Create Report Rmd (Main + Sections EN)

**Files:**
- Create: `benchmarks/report/benchmark_report.Rmd`
- Create: `benchmarks/report/sections/01_introduction_en.Rmd`
- Create: `benchmarks/report/sections/02_methods_en.Rmd`
- Create: `benchmarks/report/sections/03_simulation_design_en.Rmd`
- Create: `benchmarks/report/sections/04_results_en.Rmd`
- Create: `benchmarks/report/sections/05_discussion_en.Rmd`
- Create: `benchmarks/report/sections/06_appendix_en.Rmd`

**Step 1: Write benchmark_report.Rmd**

```rmd
---
title: "Benchmark of High-Dimensional Causal Inference Methods"
subtitle: "cfomics Simulation Study"
author: "cfomics development team"
date: "`r Sys.Date()`"
output:
  bookdown::pdf_document2:
    toc: true
    toc_depth: 3
    number_sections: true
    latex_engine: xelatex
    keep_tex: true
params:
  lang: "en"
  results_dir: "../results"
bibliography: bibliography.bib
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
                      fig.width = 10, fig.height = 7, dpi = 300,
                      out.width = "100%")

library(ggplot2)
library(patchwork)
library(scales)
library(kableExtra)

source("../../R/reporting.R")
source("../../R/metrics.R")

# Load data
results_df <- read.csv(file.path(params$results_dir, "summary", "all_results.csv"),
                       stringsAsFactors = FALSE)
summary_df <- read.csv(file.path(params$results_dir, "summary", "summary.csv"),
                       stringsAsFactors = FALSE)
bv_df <- bias_variance_decomposition(results_df)
friedman_results <- run_friedman_tests(results_df)

lang <- params$lang
```

```{r child-sections, child=paste0("sections/01_introduction_", lang, ".Rmd")}
```

```{r child-methods, child=paste0("sections/02_methods_", lang, ".Rmd")}
```

```{r child-design, child=paste0("sections/03_simulation_design_", lang, ".Rmd")}
```

```{r child-results, child=paste0("sections/04_results_", lang, ".Rmd")}
```

```{r child-discussion, child=paste0("sections/05_discussion_", lang, ".Rmd")}
```

```{r child-appendix, child=paste0("sections/06_appendix_", lang, ".Rmd")}
```
```

**Step 2: Write English section files**

Each section should contain the academic prose and embedded R code chunks for figures/tables. The content should be publication-quality.

`01_introduction_en.Rmd`: Background on high-dimensional causal inference, motivation for systematic comparison, research questions.

`02_methods_en.Rmd`: Mathematical description of each method (G-formula, HDML/AIPW, HDPS/IPW, BCF, TMLE), their estimands, assumptions, theoretical properties.

`03_simulation_design_en.Rmd`: DGP definitions with equations for S1-S11, rationale for each scenario, evaluation metrics (Bias, RMSE, PEHE, Coverage, CI Length), statistical testing framework (Friedman + Nemenyi).

`04_results_en.Rmd`: Embeds all figures (Fig 2-8) and tables (Tab 1-5) with captions and inline results discussion.

`05_discussion_en.Rmd`: Cross-scenario interpretation, method recommendation matrix, limitations, future work.

`06_appendix_en.Rmd`: Full parameter table, session info, error summary.

(Full content for each section file will be written during implementation - each is 200-500 lines of Rmd with embedded R chunks.)

**Step 3: Create empty bibliography.bib**

```bibtex
@article{friedman1937,
  title={The use of ranks to avoid the assumption of normality implicit in the analysis of variance},
  author={Friedman, Milton},
  journal={Journal of the American Statistical Association},
  volume={32},
  pages={675--701},
  year={1937}
}
```

**Step 4: Commit**

```bash
git add benchmarks/report/
git commit -m "feat(bench): add report Rmd framework with English sections"
```

---

### Task 10: Create Report Sections (Japanese)

**Files:**
- Create: `benchmarks/report/sections/01_introduction_ja.Rmd`
- Create: `benchmarks/report/sections/02_methods_ja.Rmd`
- Create: `benchmarks/report/sections/03_simulation_design_ja.Rmd`
- Create: `benchmarks/report/sections/04_results_ja.Rmd`
- Create: `benchmarks/report/sections/05_discussion_ja.Rmd`
- Create: `benchmarks/report/sections/06_appendix_ja.Rmd`

**Step 1: Write all 6 Japanese section files**

Mirror the English sections with identical R code chunks but Japanese prose. Same figures and tables are embedded (code chunks are language-independent).

**Step 2: Commit**

```bash
git add benchmarks/report/sections/*_ja.Rmd
git commit -m "feat(bench): add Japanese report sections"
```

---

### Task 11: Create README.md

**Files:**
- Create: `benchmarks/README.md`

**Step 1: Write README.md**

```markdown
# cfomics Benchmark Simulation

Production-scale benchmark comparing 5 causal inference methods across 13 DGP scenarios.

## Quick Start

```bash
# Install dependencies
Rscript -e 'install.packages(c("future", "furrr", "PMCMRplus"))'

# Smoke test (few minutes)
Rscript benchmarks/run_full_benchmark.R --smoke

# Full benchmark (3-6 hours on 8+ cores)
Rscript benchmarks/run_full_benchmark.R

# Aggregate results
Rscript benchmarks/aggregate_results.R

# Generate reports
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="en"))'
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="ja"))'
```

## Methods

gformula, hdml, hdps, bcf, tmle

## Scenarios

S1 (Baseline), S2 (Dimension), S3 (Heterogeneity), S4a-c (HTE types),
S5 (Nonlinearity), S6 (Dense confounding), S7 (Overlap), S8 (Covariate shift),
S9 (Correlation), S10 (Unmeasured confounding), S11 (Collider)

## Interrupt/Resume

Re-run `run_full_benchmark.R` to continue from where it stopped.
Add `--retry` to re-run failed jobs.
```

**Step 2: Commit**

```bash
git add benchmarks/README.md
git commit -m "docs(bench): add benchmark README"
```

---

### Task 12: End-to-End Smoke Test

**Step 1: Clean previous smoke results**

```bash
rm -rf benchmarks/results/raw/*
```

**Step 2: Run full smoke pipeline**

```bash
Rscript benchmarks/run_full_benchmark.R --smoke
Rscript benchmarks/aggregate_results.R
```

Expected: Smoke completes in 2-5 minutes with ~30 result files and summary CSV.

**Step 3: Verify results structure**

```bash
Rscript -e '
  df <- read.csv("benchmarks/results/summary/summary.csv")
  cat("Rows:", nrow(df), "\n")
  cat("Methods:", paste(unique(df$method), collapse=", "), "\n")
  cat("Scenarios:", length(unique(df$scenario_id)), "\n")
  cat("Columns:", paste(names(df), collapse=", "), "\n")
'
```

Expected: 15 rows (3 scenarios x 5 methods), all expected columns present.

**Step 4: Generate figures (test reporting.R)**

```bash
Rscript -e '
  source("benchmarks/R/reporting.R")
  source("benchmarks/R/metrics.R")
  results_df <- read.csv("benchmarks/results/summary/all_results.csv", stringsAsFactors=FALSE)
  summary_df <- read.csv("benchmarks/results/summary/summary.csv", stringsAsFactors=FALSE)
  bv_df <- bias_variance_decomposition(results_df)
  generate_all_figures(results_df, summary_df, bv_df, "benchmarks/results/plots")
  cat("Plots:", length(list.files("benchmarks/results/plots")), "\n")
'
```

Expected: Multiple PDF/PNG files created in plots directory.

**Step 5: Commit any fixes**

```bash
git add -A benchmarks/
git commit -m "fix(bench): fixes from smoke test"
```

---

### Task 13: Run Full Production Benchmark

**Step 1: Clean smoke results**

```bash
rm -rf benchmarks/results/raw/*
rm -rf benchmarks/results/summary/*
rm -rf benchmarks/results/plots/*
rm -rf benchmarks/results/tables/*
```

**Step 2: Start production run**

```bash
Rscript benchmarks/run_full_benchmark.R
```

Expected: ~11,750 jobs across 13 workers. Estimated 3-6 hours.

**Step 3: Monitor progress**

```bash
ls benchmarks/results/raw/ | wc -l
```

Check periodically. Target: 11,750 files.

**Step 4: Aggregate**

```bash
Rscript benchmarks/aggregate_results.R
```

**Step 5: Generate all figures and tables**

```bash
Rscript -e '
  source("benchmarks/R/reporting.R")
  source("benchmarks/R/metrics.R")
  results_df <- read.csv("benchmarks/results/summary/all_results.csv", stringsAsFactors=FALSE)
  summary_df <- read.csv("benchmarks/results/summary/summary.csv", stringsAsFactors=FALSE)
  bv_df <- bias_variance_decomposition(results_df)
  friedman_results <- run_friedman_tests(results_df)
  generate_all_figures(results_df, summary_df, bv_df, "benchmarks/results/plots")
  generate_all_tables(summary_df, friedman_results, "benchmarks/results/tables")
'
```

**Step 6: Render reports**

```bash
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="en"), output_dir="../results")'
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="ja"), output_dir="../results")'
```

Expected: Two PDF reports in `benchmarks/results/`.
