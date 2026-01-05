# PR-5: cfomicsSim / cfomicsBench Skeleton

## 目的

`packages/cfomicsSim/` と `packages/cfomicsBench/` パッケージの骨組みを作成し、単独で `R CMD check` が通る状態にする。

## 前提条件

- [ ] PR-1 がマージ済み（packages/ 構造）
- [ ] PR-2 がマージ済み（registry システム）

## ブランチ

```bash
git checkout main
git pull origin main
git checkout -b feature/sim-bench-skeleton
```

---

## Part A: cfomicsSim

### Step A1: ディレクトリ構造の作成

```bash
mkdir -p packages/cfomicsSim/R
mkdir -p packages/cfomicsSim/inst/extdata/scenarios
mkdir -p packages/cfomicsSim/tests/testthat
mkdir -p packages/cfomicsSim/man
```

### Step A2: DESCRIPTION

**ファイル**: `packages/cfomicsSim/DESCRIPTION`

```
Package: cfomicsSim
Title: Simulation and Data Generating Processes for cfomics
Version: 0.1.0
Authors@R:
    person("Yusuke", "Matsui", , "matsui@example.com", role = c("aut", "cre"),
           comment = c(ORCID = "0000-0000-0000-0000"))
Description: Provides simulation scenarios and data generating processes (DGPs)
    for benchmarking causal inference methods in the cfomics ecosystem.
    Includes realistic omics-like data generation with known ground truth.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.1
Depends:
    R (>= 4.2.0)
Imports:
    cfomics,
    rlang,
    cli,
    stats,
    yaml
Suggests:
    testthat (>= 3.0.0)
Config/testthat/edition: 3
```

### Step A3: LICENSE

**ファイル**: `packages/cfomicsSim/LICENSE`

```
YEAR: 2024
COPYRIGHT HOLDER: Yusuke Matsui
```

### Step A4: パッケージドキュメント

**ファイル**: `packages/cfomicsSim/R/cfomicsSim-package.R`

```r
#' @keywords internal
"_PACKAGE"

#' @import cfomics
#' @importFrom rlang abort warn
#' @importFrom cli cli_alert_info
#' @importFrom stats rnorm rbinom plogis
#' @importFrom yaml read_yaml
NULL
```

### Step A5: シナリオ管理

**ファイル**: `packages/cfomicsSim/R/scenarios.R`

```r
#' List available simulation scenarios
#'
#' @return A data.frame with scenario information
#' @export
#'
#' @examples
#' cf_scenarios()
cf_scenarios <- function() {
  scenarios <- .get_builtin_scenarios()

  data.frame(
    id = vapply(scenarios, `[[`, character(1), "id"),
    name = vapply(scenarios, `[[`, character(1), "name"),
    description = vapply(scenarios, function(s) s$description %||% "", character(1)),
    stringsAsFactors = FALSE
  )
}

#' Get a specific scenario by ID
#'
#' @param id Character. Scenario identifier
#'
#' @return A scenario specification list
#' @export
#'
#' @examples
#' scenario <- cf_scenario_get("S00_linear")
#' str(scenario)
cf_scenario_get <- function(id) {
  scenarios <- .get_builtin_scenarios()
  scenario_ids <- vapply(scenarios, `[[`, character(1), "id")

  idx <- match(id, scenario_ids)

  if (is.na(idx)) {
    available <- paste(scenario_ids, collapse = ", ")
    rlang::abort(
      paste0("Scenario '", id, "' not found.\nAvailable: ", available),
      class = "cfomics_scenario_error"
    )
  }

  scenarios[[idx]]
}

#' Get built-in scenarios
#' @keywords internal
#' @noRd
.get_builtin_scenarios <- function() {
  list(
    list(
      id = "S00_linear",
      name = "Linear DGP",
      description = "Simple linear treatment effect with confounding",
      params = list(
        ate_true = 2.0,
        effect_type = "constant",
        confounding_strength = 0.5
      ),
      mechanisms = list(
        confounding = TRUE,
        heterogeneity = FALSE,
        nonlinearity = FALSE
      )
    ),
    list(
      id = "S01_heterogeneous",
      name = "Heterogeneous Effects",
      description = "Treatment effect varies with covariates",
      params = list(
        ate_true = 2.0,
        effect_type = "heterogeneous",
        confounding_strength = 0.5
      ),
      mechanisms = list(
        confounding = TRUE,
        heterogeneity = TRUE,
        nonlinearity = FALSE
      )
    ),
    list(
      id = "S02_nonlinear",
      name = "Nonlinear DGP",
      description = "Nonlinear outcome model",
      params = list(
        ate_true = 2.0,
        effect_type = "constant",
        confounding_strength = 0.5
      ),
      mechanisms = list(
        confounding = TRUE,
        heterogeneity = FALSE,
        nonlinearity = TRUE
      )
    )
  )
}
```

### Step A6: シミュレーション実行

**ファイル**: `packages/cfomicsSim/R/simulate.R`

```r
#' Simulate data from a scenario
#'
#' @param scenario Character (scenario ID) or scenario specification list
#' @param n Integer. Number of samples
#' @param p Integer. Number of covariates
#' @param seed Integer or NULL. Random seed for reproducibility
#' @param ... Additional parameters to override scenario defaults
#'
#' @return A list with:
#'   - data: data.frame with simulated data
#'   - truth: list with true causal quantities (ITE, ATE, Y0, Y1)
#'   - scenario: the scenario specification used
#'   - seed: the random seed used
#'
#' @export
#'
#' @examples
#' # Simple simulation
#' sim <- cf_simulate("S00_linear", n = 100, seed = 42)
#' head(sim$data)
#' sim$truth$ate
#'
#' # Extract true values
#' cf_truth(sim)
cf_simulate <- function(scenario, n = 1000, p = 10, seed = NULL, ...) {
  # Get scenario spec
  if (is.character(scenario)) {
    spec <- cf_scenario_get(scenario)
  } else if (is.list(scenario)) {
    spec <- scenario
  } else {
    rlang::abort("scenario must be a character ID or a list specification")
  }

  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Override params
  dots <- list(...)
  for (nm in names(dots)) {
    if (nm %in% names(spec$params)) {
      spec$params[[nm]] <- dots[[nm]]
    }
  }

  # Generate data based on scenario type
  result <- switch(spec$id,
    "S00_linear" = .simulate_linear(n, p, spec),
    "S01_heterogeneous" = .simulate_heterogeneous(n, p, spec),
    "S02_nonlinear" = .simulate_nonlinear(n, p, spec),
    .simulate_linear(n, p, spec)  # default
  )

  # Add metadata
  result$scenario <- spec
  result$seed <- seed
  result$n <- n
  result$p <- p
  result$created_at <- Sys.time()

  structure(result, class = "cf_simulation")
}

#' Extract truth from simulation
#'
#' @param x A cf_simulation object
#'
#' @return A list with true causal quantities
#' @export
cf_truth <- function(x) {
  UseMethod("cf_truth")
}

#' @export
cf_truth.cf_simulation <- function(x) {
  x$truth
}

#' @export
cf_truth.default <- function(x) {
  if (is.list(x) && "truth" %in% names(x)) {
    return(x$truth)
  }

  rlang::abort("Cannot extract truth from this object")
}

#' @export
print.cf_simulation <- function(x, ...) {
  cli::cli_h1("cf_simulation object")
  cli::cli_text("Scenario: {.val {x$scenario$id}} - {x$scenario$name}")
  cli::cli_text("Samples: {.val {x$n}}, Covariates: {.val {x$p}}")
  cli::cli_text("True ATE: {.val {round(x$truth$ate, 4)}}")
  if (!is.null(x$seed)) {
    cli::cli_text("Seed: {.val {x$seed}}")
  }
  invisible(x)
}

# Internal simulation functions
.simulate_linear <- function(n, p, spec) {
  # Covariates
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  # Confounding through X1
  conf_strength <- spec$params$confounding_strength %||% 0.5
  propensity <- plogis(conf_strength * X[, 1])
  T <- rbinom(n, 1, propensity)

  # Potential outcomes
  ate <- spec$params$ate_true %||% 2.0
  Y0 <- 1 + 0.5 * X[, 1] + 0.3 * X[, 2] + rnorm(n)
  Y1 <- Y0 + ate

  # Observed outcome
  Y <- ifelse(T == 1, Y1, Y0)

  # Build data.frame
  data <- as.data.frame(X)
  data$T <- T
  data$Y <- Y

  # Truth
  truth <- list(
    ite = Y1 - Y0,
    ate = mean(Y1 - Y0),
    y0 = Y0,
    y1 = Y1
  )

  list(data = data, truth = truth)
}

.simulate_heterogeneous <- function(n, p, spec) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  conf_strength <- spec$params$confounding_strength %||% 0.5
  propensity <- plogis(conf_strength * X[, 1])
  T <- rbinom(n, 1, propensity)

  # Heterogeneous treatment effect
  base_ate <- spec$params$ate_true %||% 2.0
  ite_true <- base_ate + 0.5 * X[, 1] - 0.3 * X[, 2]

  Y0 <- 1 + 0.5 * X[, 1] + 0.3 * X[, 2] + rnorm(n)
  Y1 <- Y0 + ite_true

  Y <- ifelse(T == 1, Y1, Y0)

  data <- as.data.frame(X)
  data$T <- T
  data$Y <- Y

  truth <- list(
    ite = ite_true,
    ate = mean(ite_true),
    y0 = Y0,
    y1 = Y1
  )

  list(data = data, truth = truth)
}

.simulate_nonlinear <- function(n, p, spec) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))

  conf_strength <- spec$params$confounding_strength %||% 0.5
  propensity <- plogis(conf_strength * X[, 1] + 0.2 * X[, 1]^2)
  T <- rbinom(n, 1, propensity)

  ate <- spec$params$ate_true %||% 2.0

  # Nonlinear outcome
  Y0 <- 1 + sin(X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n)
  Y1 <- Y0 + ate

  Y <- ifelse(T == 1, Y1, Y0)

  data <- as.data.frame(X)
  data$T <- T
  data$Y <- Y

  truth <- list(
    ite = Y1 - Y0,
    ate = mean(Y1 - Y0),
    y0 = Y0,
    y1 = Y1
  )

  list(data = data, truth = truth)
}
```

### Step A7: テスト

**ファイル**: `packages/cfomicsSim/tests/testthat.R`

```r
library(testthat)
library(cfomicsSim)

test_check("cfomicsSim")
```

**ファイル**: `packages/cfomicsSim/tests/testthat/test-scenarios.R`

```r
test_that("cf_scenarios returns data.frame", {
  scenarios <- cf_scenarios()

  expect_s3_class(scenarios, "data.frame")
  expect_true("id" %in% names(scenarios))
  expect_true("name" %in% names(scenarios))
  expect_true(nrow(scenarios) >= 3)
})

test_that("cf_scenario_get returns valid scenario", {
  scenario <- cf_scenario_get("S00_linear")

  expect_type(scenario, "list")
  expect_equal(scenario$id, "S00_linear")
  expect_true("params" %in% names(scenario))
})

test_that("cf_scenario_get errors for unknown scenario", {
  expect_error(
    cf_scenario_get("nonexistent"),
    class = "cfomics_scenario_error"
  )
})
```

**ファイル**: `packages/cfomicsSim/tests/testthat/test-simulate.R`

```r
test_that("cf_simulate works with scenario ID", {
  sim <- cf_simulate("S00_linear", n = 100, seed = 42)

  expect_s3_class(sim, "cf_simulation")
  expect_equal(nrow(sim$data), 100)
  expect_true("Y" %in% names(sim$data))
  expect_true("T" %in% names(sim$data))
})

test_that("cf_simulate is reproducible with seed", {
  sim1 <- cf_simulate("S00_linear", n = 50, seed = 123)
  sim2 <- cf_simulate("S00_linear", n = 50, seed = 123)

  expect_equal(sim1$data$Y, sim2$data$Y)
  expect_equal(sim1$truth$ite, sim2$truth$ite)
})

test_that("cf_truth extracts truth correctly", {
  sim <- cf_simulate("S00_linear", n = 100, seed = 42)
  truth <- cf_truth(sim)

  expect_type(truth, "list")
  expect_true("ite" %in% names(truth))
  expect_true("ate" %in% names(truth))
  expect_true("y0" %in% names(truth))
  expect_true("y1" %in% names(truth))
  expect_equal(length(truth$ite), 100)
})

test_that("truth is self-consistent", {
  sim <- cf_simulate("S00_linear", n = 100, seed = 42)
  truth <- cf_truth(sim)

  # ITE = Y1 - Y0
  expect_equal(truth$ite, truth$y1 - truth$y0)

  # ATE = mean(ITE)
  expect_equal(truth$ate, mean(truth$ite))
})

test_that("heterogeneous scenario has varying ITE", {
  sim <- cf_simulate("S01_heterogeneous", n = 100, seed = 42)
  truth <- cf_truth(sim)

  # ITE should have variation
  expect_true(sd(truth$ite) > 0)
})
```

### Step A8: NEWS.md

**ファイル**: `packages/cfomicsSim/NEWS.md`

```markdown
# cfomicsSim 0.1.0

## Initial Release

* Scenario management: `cf_scenarios()`, `cf_scenario_get()`
* Simulation: `cf_simulate()` with seed reproducibility
* Truth extraction: `cf_truth()`
* Built-in scenarios:
  - S00_linear: Linear DGP with confounding
  - S01_heterogeneous: Heterogeneous treatment effects
  - S02_nonlinear: Nonlinear outcome model
```

---

## Part B: cfomicsBench

### Step B1: ディレクトリ構造の作成

```bash
mkdir -p packages/cfomicsBench/R
mkdir -p packages/cfomicsBench/tests/testthat
mkdir -p packages/cfomicsBench/man
```

### Step B2: DESCRIPTION

**ファイル**: `packages/cfomicsBench/DESCRIPTION`

```
Package: cfomicsBench
Title: Benchmarking Framework for cfomics
Version: 0.1.0
Authors@R:
    person("Yusuke", "Matsui", , "matsui@example.com", role = c("aut", "cre"),
           comment = c(ORCID = "0000-0000-0000-0000"))
Description: Provides a benchmarking framework for evaluating causal inference
    methods in the cfomics ecosystem. Supports reproducible experiments with
    multiple scenarios, methods, and replications.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.1
Depends:
    R (>= 4.2.0)
Imports:
    cfomics,
    cfomicsSim,
    rlang,
    cli
Suggests:
    testthat (>= 3.0.0),
    ggplot2
Config/testthat/edition: 3
```

### Step B3: LICENSE

**ファイル**: `packages/cfomicsBench/LICENSE`

```
YEAR: 2024
COPYRIGHT HOLDER: Yusuke Matsui
```

### Step B4: パッケージドキュメント

**ファイル**: `packages/cfomicsBench/R/cfomicsBench-package.R`

```r
#' @keywords internal
"_PACKAGE"

#' @import cfomics
#' @import cfomicsSim
#' @importFrom rlang abort warn
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
NULL
```

### Step B5: ベンチマーク実行

**ファイル**: `packages/cfomicsBench/R/benchmark_run.R`

```r
#' Run benchmark experiments
#'
#' @param scenarios Character vector of scenario IDs
#' @param methods Character vector of method names
#' @param n_reps Integer. Number of replications per scenario-method combination
#' @param n_samples Integer. Number of samples per simulation
#' @param seeds Integer vector or NULL. Seeds for reproducibility
#' @param verbose Logical. Show progress?
#' @param ... Additional arguments passed to cf_fit
#'
#' @return A cf_benchmark object
#' @export
#'
#' @examples
#' \donttest{
#' results <- cf_benchmark_run(
#'   scenarios = "S00_linear",
#'   methods = "gformula",
#'   n_reps = 2,
#'   n_samples = 100
#' )
#' }
cf_benchmark_run <- function(scenarios,
                              methods,
                              n_reps = 10,
                              n_samples = 500,
                              seeds = NULL,
                              verbose = TRUE,
                              ...) {
  # Validate inputs
  available_scenarios <- cfomicsSim::cf_scenarios()$id
  unknown_scenarios <- setdiff(scenarios, available_scenarios)
  if (length(unknown_scenarios) > 0) {
    rlang::abort(paste("Unknown scenarios:", paste(unknown_scenarios, collapse = ", ")))
  }

  available_methods <- cfomics::cf_methods()$id
  unknown_methods <- setdiff(methods, available_methods)
  if (length(unknown_methods) > 0) {
    rlang::warn(paste("Methods not available:", paste(unknown_methods, collapse = ", ")))
    methods <- intersect(methods, available_methods)
    if (length(methods) == 0) {
      rlang::abort("No available methods to benchmark")
    }
  }

  # Generate seeds if not provided
  if (is.null(seeds)) {
    seeds <- seq_len(n_reps)
  }

  # Setup results storage
  results <- list()
  total <- length(scenarios) * length(methods) * length(seeds)

  if (verbose) {
    cli::cli_progress_bar("Running benchmarks", total = total)
  }

  idx <- 0
  for (scenario_id in scenarios) {
    for (method in methods) {
      for (seed in seeds) {
        idx <- idx + 1
        if (verbose) cli::cli_progress_update()

        # Run single experiment
        result <- tryCatch({
          .run_single_experiment(
            scenario_id = scenario_id,
            method = method,
            n_samples = n_samples,
            seed = seed,
            ...
          )
        }, error = function(e) {
          list(
            scenario = scenario_id,
            method = method,
            seed = seed,
            error = conditionMessage(e),
            metrics = NULL
          )
        })

        results[[idx]] <- result
      }
    }
  }

  if (verbose) cli::cli_progress_done()

  # Build output
  output <- list(
    results = results,
    config = list(
      scenarios = scenarios,
      methods = methods,
      n_reps = n_reps,
      n_samples = n_samples,
      seeds = seeds
    ),
    run_at = Sys.time()
  )

  structure(output, class = "cf_benchmark")
}

#' Run single experiment
#' @keywords internal
#' @noRd
.run_single_experiment <- function(scenario_id, method, n_samples, seed, ...) {
  # Simulate data
  sim <- cfomicsSim::cf_simulate(scenario_id, n = n_samples, seed = seed)
  truth <- cfomicsSim::cf_truth(sim)

  # Fit model
  start_time <- Sys.time()

  fit <- cfomics::cf_fit(
    Y ~ T | .,
    data = sim$data,
    method = method,
    ...
  )

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Calculate metrics
  metrics <- .calculate_metrics(fit, truth)
  metrics$elapsed_time <- elapsed

  list(
    scenario = scenario_id,
    method = method,
    seed = seed,
    error = NULL,
    metrics = metrics
  )
}

#' Calculate evaluation metrics
#' @keywords internal
#' @noRd
.calculate_metrics <- function(fit, truth) {
  # Extract estimates
  est_ate <- fit$fit$res$ate
  est_ite <- fit$fit$res$ite

  true_ate <- truth$ate
  true_ite <- truth$ite

  list(
    ate_bias = est_ate - true_ate,
    ate_error = abs(est_ate - true_ate),
    ite_rmse = sqrt(mean((est_ite - true_ite)^2)),
    ite_mae = mean(abs(est_ite - true_ite)),
    ite_correlation = cor(est_ite, true_ite)
  )
}

#' @export
print.cf_benchmark <- function(x, ...) {
  cli::cli_h1("cf_benchmark results")
  cli::cli_text("Scenarios: {.val {paste(x$config$scenarios, collapse = ', ')}}")
  cli::cli_text("Methods: {.val {paste(x$config$methods, collapse = ', ')}}")
  cli::cli_text("Replications: {.val {x$config$n_reps}}")
  cli::cli_text("Total experiments: {.val {length(x$results)}}")

  # Count errors
  n_errors <- sum(vapply(x$results, function(r) !is.null(r$error), logical(1)))
  if (n_errors > 0) {
    cli::cli_alert_warning("{n_errors} experiments failed")
  }

  invisible(x)
}
```

### Step B6: 結果集計

**ファイル**: `packages/cfomicsBench/R/benchmark_summarize.R`

```r
#' Summarize benchmark results
#'
#' @param x A cf_benchmark object
#'
#' @return A data.frame with summarized metrics
#' @export
cf_benchmark_summarize <- function(x) {
  if (!inherits(x, "cf_benchmark")) {
    rlang::abort("x must be a cf_benchmark object")
  }

  # Filter successful results
  successful <- Filter(function(r) is.null(r$error), x$results)

  if (length(successful) == 0) {
    rlang::warn("No successful experiments to summarize")
    return(data.frame())
  }

  # Build raw data
  rows <- lapply(successful, function(r) {
    data.frame(
      scenario = r$scenario,
      method = r$method,
      seed = r$seed,
      ate_bias = r$metrics$ate_bias,
      ate_error = r$metrics$ate_error,
      ite_rmse = r$metrics$ite_rmse,
      ite_mae = r$metrics$ite_mae,
      ite_correlation = r$metrics$ite_correlation,
      elapsed_time = r$metrics$elapsed_time,
      stringsAsFactors = FALSE
    )
  })

  raw <- do.call(rbind, rows)

  # Summarize by scenario and method
  groups <- split(raw, list(raw$scenario, raw$method))

  summary_rows <- lapply(names(groups), function(nm) {
    g <- groups[[nm]]
    parts <- strsplit(nm, "\\.")[[1]]

    data.frame(
      scenario = parts[1],
      method = parts[2],
      n_reps = nrow(g),
      ate_bias_mean = mean(g$ate_bias),
      ate_bias_sd = sd(g$ate_bias),
      ate_error_mean = mean(g$ate_error),
      ite_rmse_mean = mean(g$ite_rmse),
      ite_rmse_sd = sd(g$ite_rmse),
      ite_correlation_mean = mean(g$ite_correlation, na.rm = TRUE),
      elapsed_time_mean = mean(g$elapsed_time),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_rows)
}
```

### Step B7: 可視化

**ファイル**: `packages/cfomicsBench/R/benchmark_plot.R`

```r
#' Plot benchmark results
#'
#' @param x A cf_benchmark object
#' @param metric Character. Metric to plot
#' @param ... Additional arguments
#'
#' @return A ggplot object (if ggplot2 available) or invisible NULL
#' @export
cf_benchmark_plot <- function(x, metric = "ate_error", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    rlang::warn("ggplot2 is required for plotting. Install with: install.packages('ggplot2')")
    return(invisible(NULL))
  }

  summary <- cf_benchmark_summarize(x)

  if (nrow(summary) == 0) {
    rlang::warn("No data to plot")
    return(invisible(NULL))
  }

  # Select metric column
  metric_col <- paste0(metric, "_mean")
  if (!metric_col %in% names(summary)) {
    available <- grep("_mean$", names(summary), value = TRUE)
    rlang::abort(paste0(
      "Metric '", metric, "' not found.\n",
      "Available: ", paste(gsub("_mean$", "", available), collapse = ", ")
    ))
  }

  p <- ggplot2::ggplot(summary, ggplot2::aes(
    x = .data$method,
    y = .data[[metric_col]],
    fill = .data$method
  )) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~scenario) +
    ggplot2::labs(
      title = paste("Benchmark Results:", metric),
      x = "Method",
      y = metric
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}
```

### Step B8: テスト

**ファイル**: `packages/cfomicsBench/tests/testthat.R`

```r
library(testthat)
library(cfomicsBench)

test_check("cfomicsBench")
```

**ファイル**: `packages/cfomicsBench/tests/testthat/test-benchmark.R`

```r
test_that("cf_benchmark_run works with minimal inputs", {
  skip_if_not_installed("cfomics")
  skip_if_not_installed("cfomicsSim")

  # Minimal benchmark
  results <- cf_benchmark_run(
    scenarios = "S00_linear",
    methods = "gformula",
    n_reps = 1,
    n_samples = 50,
    verbose = FALSE
  )

  expect_s3_class(results, "cf_benchmark")
  expect_equal(length(results$results), 1)
})

test_that("cf_benchmark_summarize returns data.frame", {
  skip_if_not_installed("cfomics")
  skip_if_not_installed("cfomicsSim")

  results <- cf_benchmark_run(
    scenarios = "S00_linear",
    methods = "gformula",
    n_reps = 2,
    n_samples = 50,
    verbose = FALSE
  )

  summary <- cf_benchmark_summarize(results)

  expect_s3_class(summary, "data.frame")
  expect_true("ate_bias_mean" %in% names(summary))
  expect_true("ite_rmse_mean" %in% names(summary))
})

test_that("benchmark handles errors gracefully", {
  skip_if_not_installed("cfomics")
  skip_if_not_installed("cfomicsSim")

  # Request unavailable method
  results <- suppressWarnings(cf_benchmark_run(
    scenarios = "S00_linear",
    methods = c("gformula", "nonexistent_method"),
    n_reps = 1,
    n_samples = 50,
    verbose = FALSE
  ))

  # Should still have results for gformula
  expect_true(length(results$results) >= 1)
})
```

### Step B9: NEWS.md

**ファイル**: `packages/cfomicsBench/NEWS.md`

```markdown
# cfomicsBench 0.1.0

## Initial Release

* `cf_benchmark_run()`: Run benchmark experiments
* `cf_benchmark_summarize()`: Summarize results
* `cf_benchmark_plot()`: Visualize results (requires ggplot2)
* Metrics: ATE bias/error, ITE RMSE/MAE/correlation, elapsed time
```

---

## 検証手順

### Step A: 各パッケージのドキュメント生成

```bash
cd packages/cfomicsSim && Rscript -e "devtools::document()"
cd packages/cfomicsBench && Rscript -e "devtools::document()"
```

### Step B: パッケージチェック

```bash
# cfomicsSim
cd packages/cfomicsSim
R CMD build .
R CMD check cfomicsSim_*.tar.gz --as-cran

# cfomicsBench
cd packages/cfomicsBench
R CMD build .
R CMD check cfomicsBench_*.tar.gz --as-cran
```

### Step C: 統合テスト

```r
devtools::install("packages/cfomics")
devtools::install("packages/cfomicsSim")
devtools::install("packages/cfomicsBench")

library(cfomicsSim)
library(cfomicsBench)

# Simulation test
sim <- cf_simulate("S00_linear", n = 100, seed = 42)
print(sim)

# Benchmark test
results <- cf_benchmark_run(
  scenarios = c("S00_linear", "S01_heterogeneous"),
  methods = c("gformula", "grf"),
  n_reps = 3,
  n_samples = 200
)

cf_benchmark_summarize(results)
```

---

## DoD チェックリスト

### cfomicsSim
- [ ] `packages/cfomicsSim/` ディレクトリが存在
- [ ] `R CMD check` が ERROR なしで通る
- [ ] `cf_scenarios()` がシナリオ一覧を返す
- [ ] `cf_simulate()` がデータと真値を生成する
- [ ] seed で再現可能

### cfomicsBench
- [ ] `packages/cfomicsBench/` ディレクトリが存在
- [ ] `R CMD check` が ERROR なしで通る
- [ ] `cf_benchmark_run()` が動作する
- [ ] `cf_benchmark_summarize()` が集計結果を返す

---

## コミット手順

```bash
# cfomicsSim
git add packages/cfomicsSim/
git commit -m "feat: add cfomicsSim package for simulation

- Add scenario management (cf_scenarios, cf_scenario_get)
- Add simulation execution (cf_simulate)
- Add truth extraction (cf_truth)
- Include 3 built-in scenarios (linear, heterogeneous, nonlinear)

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"

# cfomicsBench
git add packages/cfomicsBench/
git commit -m "feat: add cfomicsBench package for benchmarking

- Add benchmark execution (cf_benchmark_run)
- Add result summarization (cf_benchmark_summarize)
- Add visualization (cf_benchmark_plot)
- Include standard metrics (ATE bias, ITE RMSE, etc.)

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## PR 作成

```bash
git push -u origin feature/sim-bench-skeleton

gh pr create --title "feat: add cfomicsSim and cfomicsBench packages" --body "$(cat <<'EOF'
## Summary

Add two new packages for simulation and benchmarking:
- `cfomicsSim`: Data generating processes with known ground truth
- `cfomicsBench`: Benchmarking framework for method evaluation

## cfomicsSim

### APIs
- `cf_scenarios()` - List available scenarios
- `cf_scenario_get(id)` - Get scenario specification
- `cf_simulate(scenario, n, seed)` - Generate simulated data
- `cf_truth(sim)` - Extract true causal quantities

### Built-in Scenarios
- S00_linear: Linear DGP with confounding
- S01_heterogeneous: Heterogeneous treatment effects
- S02_nonlinear: Nonlinear outcome model

## cfomicsBench

### APIs
- `cf_benchmark_run(scenarios, methods, n_reps)` - Run experiments
- `cf_benchmark_summarize(results)` - Aggregate metrics
- `cf_benchmark_plot(results)` - Visualize (requires ggplot2)

### Metrics
- ATE bias, ATE error
- ITE RMSE, ITE MAE, ITE correlation
- Elapsed time

## Test plan

- [ ] `R CMD check packages/cfomicsSim` passes
- [ ] `R CMD check packages/cfomicsBench` passes
- [ ] Simulation is reproducible with seed
- [ ] Benchmark runs with R-native methods

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

*Last updated: 2026-01-05*
