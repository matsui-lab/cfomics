# S12 Nonlinear Confounding Redesign Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Modify `dgp_nonlinear_outcome` so PS depends on the same nonlinear function h(X) as the outcome, creating omitted-variable bias that breaks gformula.

**Architecture:** Replace the linear PS in `dgp_nonlinear_outcome` with a PS that depends on h(X). Update the test that assumes same T across nonlinearity levels. Re-run S12 benchmark scenarios and regenerate reports.

**Tech Stack:** R, cfomics package, ggplot2, rmarkdown

---

### Task 1: Update test for new PS behavior

**Files:**
- Modify: `packages/cfomics/tests/testthat/test-benchmark-dgp.R:340-349`

**Step 1: Write the updated test**

Replace the test at lines 340-349 with:

```r
test_that("dgp_nonlinear_outcome moderate and severe produce different outcomes", {
  set.seed(42)
  r_mod <- dgp_nonlinear_outcome(n = 1000, p = 10, nonlinearity = "moderate")
  set.seed(42)
  r_sev <- dgp_nonlinear_outcome(n = 1000, p = 10, nonlinearity = "severe")
  # Different nonlinearity levels produce different Y
  expect_false(all(r_mod$Y == r_sev$Y))
  # Both have valid structure
  expect_equal(r_mod$dgp_name, "nonlinear_outcome")
  expect_equal(r_sev$dgp_name, "nonlinear_outcome")
})
```

The old test expected `r_mod$T == r_sev$T` (same treatment), which won't hold because PS now depends on h(X) which differs between moderate/severe.

**Step 2: Run test to verify it still passes (old code)**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-benchmark-dgp.R")'`

Expected: PASS (the new test is less restrictive than the old one)

**Step 3: Commit**

```bash
git add packages/cfomics/tests/testthat/test-benchmark-dgp.R
git commit -m "test(benchmark): relax S12 test for nonlinear confounding redesign"
```

---

### Task 2: Add test for nonlinear confounding bias

**Files:**
- Modify: `packages/cfomics/tests/testthat/test-benchmark-dgp.R`

**Step 1: Write the failing test**

Add after the existing S12 tests:

```r
test_that("dgp_nonlinear_outcome creates confounding that biases OLS", {
  set.seed(123)
  r <- dgp_nonlinear_outcome(n = 5000, p = 10, nonlinearity = "moderate")
  # OLS should be biased because it omits h(X) which drives both PS and Y
  fit <- lm(r$Y ~ r$T + r$X)
  ate_ols <- coef(fit)["r$T"]
  bias <- abs(ate_ols - r$true_ate)
  # With nonlinear confounding, OLS bias should be substantial (> 0.3)
  expect_gt(bias, 0.3)
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-benchmark-dgp.R")'`

Expected: FAIL — current linear PS produces no OLS bias, so `bias < 0.3`.

**Step 3: Commit**

```bash
git add packages/cfomics/tests/testthat/test-benchmark-dgp.R
git commit -m "test(benchmark): add OLS bias test for S12 nonlinear confounding"
```

---

### Task 3: Modify dgp_nonlinear_outcome to use nonlinear PS

**Files:**
- Modify: `packages/cfomics/R/benchmark_dgp.R:864-912`

**Step 1: Replace the function implementation**

Replace lines 864-912 (the entire `dgp_nonlinear_outcome` function) with:

```r
#' Generate nonlinear outcome DGP (S12)
#'
#' Generates data with nonlinear confounding: both the outcome and propensity
#' score depend on the same nonlinear function h(X). This creates omitted-variable
#' bias for methods using linear outcome models (gformula).
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param nonlinearity Character, "moderate" or "severe"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, propensity_score, dgp_name, dgp_params
#' @export
dgp_nonlinear_outcome <- function(n = 500, p = 500,
                                   nonlinearity = "moderate",
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  tau <- 2.0

  # Nonlinear function h(X) drives both outcome and treatment
  if (nonlinearity == "moderate") {
    h <- X[, 1]^2 + X[, 2] * X[, 3] + sin(X[, 1])
  } else {
    h <- X[, 1]^3 / 3 + exp(X[, 2] / 2) + X[, 3] * X[, 4] * (X[, 1] > 0) +
      sin(2 * X[, 1]) * cos(X[, 2])
  }

  # PS depends on h(X) — creates nonlinear confounding
  ps_coef <- if (nonlinearity == "moderate") 0.8 else 0.5
  h_centered <- h - mean(h)
  ps <- stats::plogis(ps_coef * h_centered)
  ps <- pmin(pmax(ps, 0.05), 0.95)
  T <- stats::rbinom(n, 1, ps)

  # Outcome depends on same h(X)
  Y <- as.numeric(h + tau * T + stats::rnorm(n))

  list(
    X = X,
    T = T,
    Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    dgp_name = "nonlinear_outcome",
    dgp_params = list(n = n, p = p, nonlinearity = nonlinearity)
  )
}
```

Key changes from old version:
- PS now uses `plogis(ps_coef * h_centered)` instead of `plogis(X %*% beta_t)` (linear)
- h is centered before PS to maintain ~50/50 treatment balance
- PS clipped to [0.05, 0.95]
- Moderate uses ps_coef=0.8 (stronger confounding), severe uses 0.5 (weaker PS signal but harder outcome)

**Step 2: Run tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-benchmark-dgp.R")'`

Expected: ALL PASS including the new OLS bias test from Task 2.

**Step 3: Commit**

```bash
git add packages/cfomics/R/benchmark_dgp.R
git commit -m "feat(benchmark): redesign S12 with nonlinear confounding"
```

---

### Task 4: Run S12 smoke test

**Files:**
- No file changes

**Step 1: Run quick benchmark (2 reps) on S12 scenarios**

```bash
Rscript -e '
source("benchmarks/config.R")
source("benchmarks/R/runner.R")
devtools::load_all("packages/cfomics")

s12_ids <- c("S12_nonlinear_outcome_moderate", "S12_nonlinear_outcome_severe")
cfg <- BENCHMARK_CONFIG
cfg$scenarios <- cfg$scenarios[sapply(cfg$scenarios, function(s) s$id %in% s12_ids)]
cfg$n_reps <- 2L

results <- run_benchmark(cfg)
cat("Jobs completed:", nrow(results), "\n")
print(aggregate(rmse_ate ~ method + scenario_id, results, mean))
'
```

Expected: gformula should have noticeably higher RMSE than hdps.

**Step 2: If differentiation is weak, adjust ps_coef**

If gformula RMSE is not clearly higher than others, increase `ps_coef` from 0.8 to 1.2 (moderate) or from 0.5 to 0.8 (severe) and re-test.

---

### Task 5: Re-run S12 production benchmark (50 reps)

**Files:**
- Output: `benchmarks/results/raw/S12_*.rds` (overwritten)
- Output: `benchmarks/results/summary/aggregated_summary.csv` (updated)

**Step 1: Delete old S12 raw results**

```bash
rm -f benchmarks/results/raw/S12_nonlinear_outcome_*.rds
```

**Step 2: Run full benchmark (only S12 scenarios, 50 reps)**

```bash
Rscript -e '
source("benchmarks/config.R")
source("benchmarks/R/runner.R")
source("benchmarks/R/aggregation.R")
devtools::load_all("packages/cfomics")

s12_ids <- c("S12_nonlinear_outcome_moderate", "S12_nonlinear_outcome_severe")
cfg <- BENCHMARK_CONFIG
cfg$scenarios <- cfg$scenarios[sapply(cfg$scenarios, function(s) s$id %in% s12_ids)]

results <- run_benchmark(cfg)
cat("S12 jobs:", nrow(results), "\n")
print(aggregate(rmse_ate ~ method + scenario_id, results, mean))
'
```

**Step 3: Re-aggregate all results**

```bash
Rscript -e '
source("benchmarks/R/aggregation.R")
aggregate_results("benchmarks/results/raw", "benchmarks/results/summary")
cat("Done\n")
'
```

**Step 4: Commit**

```bash
git add benchmarks/results/
git commit -m "data: S12 benchmark results with nonlinear confounding"
```

---

### Task 6: Regenerate reports

**Files:**
- Output: `benchmarks/report/benchmark_report_en.html`
- Output: `benchmarks/report/benchmark_report_ja.html`

**Step 1: Render English report**

```bash
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", output_format = "html_document", output_file = "benchmark_report_en.html", params = list(lang = "en", results_dir = "../results"))'
```

**Step 2: Render Japanese report**

```bash
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", output_format = "html_document", output_file = "benchmark_report_ja.html", params = list(lang = "ja", results_dir = "../results"))'
```

**Step 3: Commit**

```bash
git add benchmarks/report/benchmark_report_en.html benchmarks/report/benchmark_report_ja.html
git commit -m "docs: regenerate reports with S12 nonlinear confounding results"
```
