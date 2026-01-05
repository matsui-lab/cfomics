# PR-6: Sim/Bench Full Implementation

## 目的

既存の `benchmark_*.R` ファイルから必要なコードを `cfomicsSim` と `cfomicsBench` に移植し、論文用ベンチマークが実行可能な状態にする。

## 前提条件

- [ ] PR-5 がマージ済み（cfomicsSim/cfomicsBench skeleton）

## ブランチ

```bash
git checkout main
git pull origin main
git checkout -b feature/sim-bench-full
```

---

## 作業手順

### Step 1: 移植対象ファイルの確認

**packages/cfomics/ にある benchmark 関連ファイル**:
- `R/benchmark_dgp.R` → cfomicsSim
- `R/benchmark_metrics.R` → cfomicsBench
- `R/benchmark_runner.R` → cfomicsBench
- `R/benchmark_report.R` → cfomicsBench
- `inst/benchmarks/DGP_v1.0.md` → cfomicsSim/inst/extdata/
- `inst/benchmarks/run_benchmark.R` → cfomicsBench/inst/scripts/

### Step 2: benchmark_dgp.R の移植

既存の DGP 実装を `cfomicsSim` に統合:

```bash
# 内容を確認
cat packages/cfomics/R/benchmark_dgp.R
```

**統合方針**:
1. 既存の DGP 関数を `cfomicsSim/R/dgp.R` として移植
2. `cf_simulate()` から呼び出せるように連携
3. 既存のシナリオ定義を YAML に変換

**ファイル**: `packages/cfomicsSim/R/dgp.R`

```r
#' @title Data Generating Processes
#' @description Advanced DGP implementations for cfomicsSim
#' @name dgp

#' Generate data from DGP specification
#'
#' @param dgp_spec DGP specification list
#' @param n Number of samples
#' @param p Number of covariates
#' @param seed Random seed
#'
#' @return List with data and truth
#' @keywords internal
#' @noRd
generate_from_dgp <- function(dgp_spec, n, p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Dispatch based on DGP type
  dgp_type <- dgp_spec$type %||% "linear"

  switch(dgp_type,
    "linear" = .dgp_linear(n, p, dgp_spec),
    "nonlinear" = .dgp_nonlinear(n, p, dgp_spec),
    "heterogeneous" = .dgp_heterogeneous(n, p, dgp_spec),
    "omics" = .dgp_omics(n, p, dgp_spec),
    rlang::abort(paste("Unknown DGP type:", dgp_type))
  )
}

#' Linear DGP
#' @keywords internal
#' @noRd
.dgp_linear <- function(n, p, spec) {
  # Covariates
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  # Treatment assignment
  conf_vars <- spec$confounding_vars %||% 1
  conf_strength <- spec$confounding_strength %||% 0.5

  propensity_linear <- conf_strength * rowSums(X[, seq_len(min(conf_vars, p)), drop = FALSE])
  propensity <- plogis(propensity_linear)
  T <- rbinom(n, 1, propensity)

  # Outcome
  outcome_vars <- spec$outcome_vars %||% 2
  outcome_coefs <- spec$outcome_coefs %||% rep(0.5, outcome_vars)
  ate <- spec$ate %||% 2.0

  Y0_linear <- spec$intercept %||% 1
  for (i in seq_len(min(outcome_vars, p))) {
    Y0_linear <- Y0_linear + outcome_coefs[i] * X[, i]
  }
  noise_sd <- spec$noise_sd %||% 1.0
  Y0 <- Y0_linear + rnorm(n, sd = noise_sd)
  Y1 <- Y0 + ate

  Y <- ifelse(T == 1, Y1, Y0)

  data <- as.data.frame(X)
  data$T <- T
  data$Y <- Y

  truth <- list(
    ite = rep(ate, n),
    ate = ate,
    y0 = Y0,
    y1 = Y1,
    propensity = propensity
  )

  list(data = data, truth = truth)
}

#' Nonlinear DGP
#' @keywords internal
#' @noRd
.dgp_nonlinear <- function(n, p, spec) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  conf_strength <- spec$confounding_strength %||% 0.5
  propensity <- plogis(conf_strength * X[, 1] + 0.3 * X[, 1]^2)
  T <- rbinom(n, 1, propensity)

  ate <- spec$ate %||% 2.0
  noise_sd <- spec$noise_sd %||% 1.0

  Y0 <- 1 + sin(X[, 1] * pi / 2) + 0.5 * X[, 2]^2 + rnorm(n, sd = noise_sd)
  Y1 <- Y0 + ate

  Y <- ifelse(T == 1, Y1, Y0)

  data <- as.data.frame(X)
  data$T <- T
  data$Y <- Y

  truth <- list(
    ite = rep(ate, n),
    ate = ate,
    y0 = Y0,
    y1 = Y1,
    propensity = propensity
  )

  list(data = data, truth = truth)
}

#' Heterogeneous treatment effect DGP
#' @keywords internal
#' @noRd
.dgp_heterogeneous <- function(n, p, spec) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  conf_strength <- spec$confounding_strength %||% 0.5
  propensity <- plogis(conf_strength * X[, 1])
  T <- rbinom(n, 1, propensity)

  base_ate <- spec$ate %||% 2.0
  het_strength <- spec$heterogeneity_strength %||% 0.5
  noise_sd <- spec$noise_sd %||% 1.0

  # Heterogeneous ITE
  ite <- base_ate + het_strength * X[, 1] - het_strength * 0.5 * X[, 2]

  Y0 <- 1 + 0.5 * X[, 1] + 0.3 * X[, 2] + rnorm(n, sd = noise_sd)
  Y1 <- Y0 + ite

  Y <- ifelse(T == 1, Y1, Y0)

  data <- as.data.frame(X)
  data$T <- T
  data$Y <- Y

  truth <- list(
    ite = ite,
    ate = mean(ite),
    y0 = Y0,
    y1 = Y1,
    propensity = propensity
  )

  list(data = data, truth = truth)
}

#' Omics-like DGP (high-dimensional, sparse effects)
#' @keywords internal
#' @noRd
.dgp_omics <- function(n, p, spec) {
  # High-dimensional covariates with correlation
  rho <- spec$correlation %||% 0.3

  # AR(1) correlation structure
  Sigma <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      Sigma[i, j] <- rho^abs(i - j)
    }
  }

  # Generate correlated covariates
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- Z %*% L
  colnames(X) <- paste0("X", seq_len(p))

  # Sparse confounding (only first few variables matter)
  n_conf <- spec$n_confounders %||% 5
  conf_strength <- spec$confounding_strength %||% 0.3
  propensity_linear <- conf_strength * rowSums(X[, seq_len(min(n_conf, p)), drop = FALSE])
  propensity <- plogis(propensity_linear)
  T <- rbinom(n, 1, propensity)

  # Sparse outcome (only first few variables matter)
  n_outcome <- spec$n_outcome_vars %||% 10
  outcome_coefs <- c(rep(0.3, n_outcome), rep(0, p - n_outcome))
  noise_sd <- spec$noise_sd %||% 1.0

  Y0 <- as.vector(X %*% outcome_coefs) + rnorm(n, sd = noise_sd)

  ate <- spec$ate %||% 1.5
  Y1 <- Y0 + ate

  Y <- ifelse(T == 1, Y1, Y0)

  data <- as.data.frame(X)
  data$T <- T
  data$Y <- Y

  truth <- list(
    ite = rep(ate, n),
    ate = ate,
    y0 = Y0,
    y1 = Y1,
    propensity = propensity
  )

  list(data = data, truth = truth)
}
```

### Step 3: シナリオ定義の拡充

**ファイル**: `packages/cfomicsSim/inst/extdata/scenarios/S03_omics.yml`

```yaml
id: S03_omics
name: High-dimensional Omics DGP
description: Sparse effects in high-dimensional correlated data
type: omics
params:
  ate: 1.5
  n_confounders: 5
  n_outcome_vars: 10
  confounding_strength: 0.3
  correlation: 0.3
  noise_sd: 1.0
mechanisms:
  confounding: true
  heterogeneity: false
  high_dimensional: true
  sparse: true
```

### Step 4: benchmark_metrics.R の移植

**ファイル**: `packages/cfomicsBench/R/metrics.R`

```r
#' @title Evaluation Metrics
#' @description Metrics for evaluating causal inference methods
#' @name metrics

#' Calculate all evaluation metrics
#'
#' @param fit cfomics_result object
#' @param truth Truth list from simulation
#'
#' @return List of metric values
#' @export
cf_calculate_metrics <- function(fit, truth) {
  est <- fit$fit$res

  metrics <- list()

  # ATE metrics
  if (!is.null(est$ate) && !is.null(truth$ate)) {
    metrics$ate_estimate <- est$ate
    metrics$ate_true <- truth$ate
    metrics$ate_bias <- est$ate - truth$ate
    metrics$ate_error <- abs(est$ate - truth$ate)
    metrics$ate_relative_error <- abs(est$ate - truth$ate) / abs(truth$ate)
  }

  # ITE metrics
  if (!is.null(est$ite) && !is.null(truth$ite)) {
    diff <- est$ite - truth$ite
    metrics$ite_rmse <- sqrt(mean(diff^2))
    metrics$ite_mae <- mean(abs(diff))
    metrics$ite_bias <- mean(diff)

    # Correlation (if variance exists)
    if (sd(est$ite) > 0 && sd(truth$ite) > 0) {
      metrics$ite_correlation <- cor(est$ite, truth$ite)
      metrics$ite_spearman <- cor(est$ite, truth$ite, method = "spearman")
    }

    # Quantile metrics
    metrics$ite_q05_error <- abs(quantile(est$ite, 0.05) - quantile(truth$ite, 0.05))
    metrics$ite_q95_error <- abs(quantile(est$ite, 0.95) - quantile(truth$ite, 0.95))
  }

  # Counterfactual metrics
  if (!is.null(est$y0_hat) && !is.null(truth$y0)) {
    metrics$y0_rmse <- sqrt(mean((est$y0_hat - truth$y0)^2))
  }

  if (!is.null(est$y1_hat) && !is.null(truth$y1)) {
    metrics$y1_rmse <- sqrt(mean((est$y1_hat - truth$y1)^2))
  }

  # Coverage (if CI available)
  if (!is.null(est$summary$ate_ci_lower) && !is.null(est$summary$ate_ci_upper)) {
    metrics$ate_covered <- truth$ate >= est$summary$ate_ci_lower &&
                           truth$ate <= est$summary$ate_ci_upper
    metrics$ci_width <- est$summary$ate_ci_upper - est$summary$ate_ci_lower
  }

  metrics
}

#' Aggregate metrics across replications
#'
#' @param metrics_list List of metric results
#'
#' @return data.frame with mean, sd, and quantiles
#' @export
cf_aggregate_metrics <- function(metrics_list) {
  if (length(metrics_list) == 0) {
    return(data.frame())
  }

  # Get all metric names
  all_names <- unique(unlist(lapply(metrics_list, names)))

  # Aggregate each metric
  results <- lapply(all_names, function(nm) {
    values <- vapply(metrics_list, function(m) {
      v <- m[[nm]]
      if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
    }, numeric(1))

    values <- values[!is.na(values)]

    if (length(values) == 0) {
      return(data.frame(
        metric = nm,
        mean = NA_real_,
        sd = NA_real_,
        q05 = NA_real_,
        q50 = NA_real_,
        q95 = NA_real_,
        n = 0L
      ))
    }

    data.frame(
      metric = nm,
      mean = mean(values),
      sd = sd(values),
      q05 = quantile(values, 0.05),
      q50 = median(values),
      q95 = quantile(values, 0.95),
      n = length(values)
    )
  })

  do.call(rbind, results)
}
```

### Step 5: benchmark_runner.R の機能拡充

**ファイル**: `packages/cfomicsBench/R/benchmark_run.R` に追加

```r
#' Run benchmark with advanced options
#'
#' @param scenarios Character vector or list of scenario specs
#' @param methods Character vector of method names
#' @param n_reps Number of replications
#' @param n_samples Number of samples per simulation
#' @param n_covariates Number of covariates
#' @param seeds Seeds for reproducibility
#' @param parallel Logical. Use parallel processing?
#' @param n_cores Number of cores for parallel processing
#' @param save_fits Logical. Save fitted models?
#' @param output_dir Directory to save results
#' @param verbose Show progress?
#' @param ... Additional arguments to cf_fit
#'
#' @return cf_benchmark object
#' @export
cf_benchmark_run_full <- function(scenarios,
                                   methods,
                                   n_reps = 100,
                                   n_samples = 500,
                                   n_covariates = 10,
                                   seeds = NULL,
                                   parallel = FALSE,
                                   n_cores = 2,
                                   save_fits = FALSE,
                                   output_dir = NULL,
                                   verbose = TRUE,
                                   ...) {
  # Validate
  if (is.null(seeds)) {
    seeds <- seq_len(n_reps)
  }

  # Create output directory if needed
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Prepare experiment grid
  experiments <- expand.grid(
    scenario = scenarios,
    method = methods,
    seed = seeds,
    stringsAsFactors = FALSE
  )

  n_total <- nrow(experiments)

  if (verbose) {
    cli::cli_alert_info("Running {n_total} experiments")
    cli::cli_alert_info("Scenarios: {paste(scenarios, collapse = ', ')}")
    cli::cli_alert_info("Methods: {paste(methods, collapse = ', ')}")
  }

  # Run experiments
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    results <- .run_parallel(experiments, n_samples, n_covariates, n_cores, ...)
  } else {
    results <- .run_sequential(experiments, n_samples, n_covariates, verbose, ...)
  }

  # Build output
  output <- list(
    results = results,
    config = list(
      scenarios = scenarios,
      methods = methods,
      n_reps = n_reps,
      n_samples = n_samples,
      n_covariates = n_covariates,
      seeds = seeds,
      parallel = parallel
    ),
    provenance = list(
      r_version = R.version.string,
      cfomics_version = as.character(packageVersion("cfomics")),
      cfomicsSim_version = as.character(packageVersion("cfomicsSim")),
      cfomicsBench_version = as.character(packageVersion("cfomicsBench")),
      run_at = Sys.time(),
      hostname = Sys.info()["nodename"]
    )
  )

  # Save if requested
  if (!is.null(output_dir)) {
    saveRDS(output, file.path(output_dir, "benchmark_results.rds"))
    if (verbose) cli::cli_alert_success("Results saved to {output_dir}")
  }

  structure(output, class = "cf_benchmark")
}

#' Run experiments sequentially
#' @keywords internal
#' @noRd
.run_sequential <- function(experiments, n_samples, n_covariates, verbose, ...) {
  n_total <- nrow(experiments)

  if (verbose) {
    pb <- cli::cli_progress_bar("Running experiments", total = n_total)
  }

  results <- vector("list", n_total)

  for (i in seq_len(n_total)) {
    if (verbose) cli::cli_progress_update(id = pb)

    exp <- experiments[i, ]

    results[[i]] <- tryCatch({
      .run_single_experiment(
        scenario_id = exp$scenario,
        method = exp$method,
        n_samples = n_samples,
        n_covariates = n_covariates,
        seed = exp$seed,
        ...
      )
    }, error = function(e) {
      list(
        scenario = exp$scenario,
        method = exp$method,
        seed = exp$seed,
        error = conditionMessage(e),
        metrics = NULL
      )
    })
  }

  if (verbose) cli::cli_progress_done(id = pb)

  results
}

#' Run experiments in parallel
#' @keywords internal
#' @noRd
.run_parallel <- function(experiments, n_samples, n_covariates, n_cores, ...) {
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl))

  # Export necessary functions and packages
  parallel::clusterEvalQ(cl, {
    library(cfomics)
    library(cfomicsSim)
    library(cfomicsBench)
  })

  dots <- list(...)

  results <- parallel::parLapply(cl, seq_len(nrow(experiments)), function(i) {
    exp <- experiments[i, ]
    tryCatch({
      .run_single_experiment(
        scenario_id = exp$scenario,
        method = exp$method,
        n_samples = n_samples,
        n_covariates = n_covariates,
        seed = exp$seed
      )
    }, error = function(e) {
      list(
        scenario = exp$scenario,
        method = exp$method,
        seed = exp$seed,
        error = conditionMessage(e),
        metrics = NULL
      )
    })
  })

  results
}
```

### Step 6: レポート生成

**ファイル**: `packages/cfomicsBench/R/benchmark_report.R`

```r
#' Generate benchmark report
#'
#' @param x cf_benchmark object
#' @param output_file Output file path (HTML or PDF)
#' @param title Report title
#'
#' @return Path to generated report
#' @export
cf_benchmark_report <- function(x,
                                 output_file = "benchmark_report.html",
                                 title = "cfomics Benchmark Report") {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    rlang::abort("rmarkdown is required. Install with: install.packages('rmarkdown')")
  }

  # For now, generate a simple text report
  # Full Quarto/Rmd template can be added later

  summary <- cf_benchmark_summarize(x)

  report <- c(
    paste("#", title),
    "",
    paste("Generated:", Sys.time()),
    "",
    "## Configuration",
    "",
    paste("- Scenarios:", paste(x$config$scenarios, collapse = ", ")),
    paste("- Methods:", paste(x$config$methods, collapse = ", ")),
    paste("- Replications:", x$config$n_reps),
    paste("- Samples per simulation:", x$config$n_samples),
    "",
    "## Results Summary",
    "",
    knitr::kable(summary, format = "markdown")
  )

  writeLines(report, sub("\\.html$", ".md", output_file))

  cli::cli_alert_success("Report saved to {output_file}")

  invisible(output_file)
}

#' Export benchmark results to various formats
#'
#' @param x cf_benchmark object
#' @param file Output file path
#' @param format Format: "csv", "rds", "json"
#'
#' @return Invisible file path
#' @export
cf_benchmark_export <- function(x, file, format = c("csv", "rds", "json")) {
  format <- match.arg(format)

  summary <- cf_benchmark_summarize(x)

  switch(format,
    csv = write.csv(summary, file, row.names = FALSE),
    rds = saveRDS(x, file),
    json = {
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        rlang::abort("jsonlite is required for JSON export")
      }
      jsonlite::write_json(summary, file, pretty = TRUE)
    }
  )

  cli::cli_alert_success("Exported to {file}")
  invisible(file)
}
```

### Step 7: 既存 benchmark ファイルの削除

移植完了後、Core から削除:

```bash
git rm packages/cfomics/R/benchmark_dgp.R
git rm packages/cfomics/R/benchmark_metrics.R
git rm packages/cfomics/R/benchmark_runner.R
git rm packages/cfomics/R/benchmark_report.R
git rm -r packages/cfomics/inst/benchmarks/
```

### Step 8: テストの更新

**ファイル**: `packages/cfomicsBench/tests/testthat/test-metrics.R`

```r
test_that("cf_calculate_metrics works", {
  skip_if_not_installed("cfomics")
  skip_if_not_installed("cfomicsSim")

  sim <- cfomicsSim::cf_simulate("S00_linear", n = 100, seed = 42)
  fit <- cfomics::cf_fit(Y ~ T | ., data = sim$data, method = "gformula")
  truth <- cfomicsSim::cf_truth(sim)

  metrics <- cf_calculate_metrics(fit, truth)

  expect_type(metrics, "list")
  expect_true("ate_bias" %in% names(metrics))
  expect_true("ite_rmse" %in% names(metrics))
})

test_that("cf_aggregate_metrics works", {
  metrics_list <- list(
    list(ate_bias = 0.1, ite_rmse = 0.5),
    list(ate_bias = -0.1, ite_rmse = 0.6),
    list(ate_bias = 0.05, ite_rmse = 0.55)
  )

  agg <- cf_aggregate_metrics(metrics_list)

  expect_s3_class(agg, "data.frame")
  expect_true("mean" %in% names(agg))
  expect_true("sd" %in% names(agg))
})
```

---

## 検証手順

```bash
# cfomicsSim
cd packages/cfomicsSim
Rscript -e "devtools::document()"
R CMD build .
R CMD check cfomicsSim_*.tar.gz --as-cran

# cfomicsBench
cd packages/cfomicsBench
Rscript -e "devtools::document()"
R CMD build .
R CMD check cfomicsBench_*.tar.gz --as-cran

# Core (benchmark ファイル削除後)
cd packages/cfomics
R CMD build .
R CMD check cfomics_*.tar.gz --as-cran
```

---

## DoD チェックリスト

- [ ] DGP 実装が cfomicsSim に移植済み
- [ ] Metrics 実装が cfomicsBench に移植済み
- [ ] Runner 実装が cfomicsBench に移植済み
- [ ] Core から benchmark ファイルが削除済み
- [ ] 全パッケージの `R CMD check` が通る
- [ ] `cf_benchmark_run_full()` が動作する
- [ ] 結果のエクスポートが動作する

---

## コミット手順

```bash
git add packages/cfomicsSim/R/dgp.R packages/cfomicsSim/inst/
git commit -m "feat(sim): add advanced DGP implementations

- Add omics-like DGP with high-dimensional sparse effects
- Add configurable DGP parameters
- Add YAML scenario definitions

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"

git add packages/cfomicsBench/R/
git commit -m "feat(bench): add full benchmark implementation

- Add comprehensive metrics calculation
- Add parallel execution support
- Add report generation and export

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"

git rm packages/cfomics/R/benchmark_*.R packages/cfomics/inst/benchmarks/
git commit -m "refactor: remove benchmark code from core (moved to cfomicsBench)

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

*Last updated: 2026-01-05*
