# Benchmark Audit Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix all issues identified in the 2026-02-07 benchmark audit, including metric naming inconsistencies, overly strict dependency checks, and missing validation.

**Architecture:** Address issues in priority order (Critical → High → Medium). Each fix includes verification. The metric rename is a search-and-replace across multiple files, while other fixes are targeted edits.

**Tech Stack:** R (testthat), Python (PyTorch)

---

## Phase 1: Critical Fix (Metric Name Inconsistency)

### Task 1: Fix mse_ate references in run_nn_benchmark.R

**Files:**
- Modify: `benchmarks/external/nn/run_nn_benchmark.R:143,166`

**Step 1: Fix the success case reference**

Change line 143 from:
```r
mse_ate = metrics$mse_ate,
```

To:
```r
squared_error_ate = metrics$squared_error_ate,
```

**Step 2: Fix the error case reference**

Change line 166 from:
```r
mse_ate = NA_real_,
```

To:
```r
squared_error_ate = NA_real_,
```

**Step 3: Verify syntax**

```bash
Rscript -e "parse(file='benchmarks/external/nn/run_nn_benchmark.R')"
```

Expected: No errors

**Step 4: Commit**

```bash
git add benchmarks/external/nn/run_nn_benchmark.R
git commit -m "fix(benchmark): rename mse_ate to squared_error_ate in NN benchmark runner"
```

---

### Task 2: Fix mse_ate references in aggregate_results.R

**Files:**
- Modify: `benchmarks/aggregate_results.R:68`

**Step 1: Update required_cols**

Change line 68 from:
```r
required_cols <- c("status", "scenario_id", "method", "mse_ate", "n", "p", "ate_true")
```

To:
```r
required_cols <- c("status", "scenario_id", "method", "squared_error_ate", "n", "p", "ate_true")
```

**Step 2: Find and update all other mse_ate references in the file**

```bash
grep -n "mse_ate" benchmarks/aggregate_results.R
```

Update each occurrence (except rmse_ate which is correct).

**Step 3: Commit**

```bash
git add benchmarks/aggregate_results.R
git commit -m "fix(benchmark): rename mse_ate to squared_error_ate in aggregate_results"
```

---

### Task 3: Fix mse_ate references in fig2_simulation.R

**Files:**
- Modify: `benchmarks/figures/fig2_simulation.R:25,63-64,193`

**Step 1: Update required_cols (line 25)**

Change from:
```r
required_cols <- c("scenario_id", "method", "status", "mse_ate", "pehe")
```

To:
```r
required_cols <- c("scenario_id", "method", "status", "squared_error_ate", "pehe")
```

**Step 2: Update aggregation code**

Find all `mse_ate` references (excluding `rmse_ate`) and replace with `squared_error_ate`.

**Step 3: Commit**

```bash
git add benchmarks/figures/fig2_simulation.R
git commit -m "fix(benchmark): rename mse_ate to squared_error_ate in fig2_simulation"
```

---

### Task 4: Fix mse_ate references in fig3_highdim.R

**Files:**
- Modify: `benchmarks/figures/fig3_highdim.R:23,55-56,140-141,236`

**Step 1: Replace all mse_ate references**

Use replace_all to change `mse_ate` to `squared_error_ate` throughout the file.

**Step 2: Verify rmse_ate references are untouched**

```bash
grep "rmse_ate" benchmarks/figures/fig3_highdim.R | wc -l
```

Expected: Multiple lines (these should remain as rmse_ate)

**Step 3: Commit**

```bash
git add benchmarks/figures/fig3_highdim.R
git commit -m "fix(benchmark): rename mse_ate to squared_error_ate in fig3_highdim"
```

---

### Task 5: Fix mse_ate references in fig4_tcga.R

**Files:**
- Modify: `benchmarks/figures/fig4_tcga.R:23,41-42,114-115,200,260`

**Step 1: Replace all mse_ate references**

Use replace_all to change `mse_ate` to `squared_error_ate` throughout the file.

**Step 2: Commit**

```bash
git add benchmarks/figures/fig4_tcga.R
git commit -m "fix(benchmark): rename mse_ate to squared_error_ate in fig4_tcga"
```

---

## Phase 2: High Priority Fixes

### Task 6: Fix CAVAE dependency check (don't require CUDA)

**Files:**
- Modify: `packages/cfomics/R/benchmark_runner.R:206-213`

**Step 1: Read current implementation**

Lines 206-213 currently check `torch$cuda$is_available()` which fails on CPU-only systems.

**Step 2: Fix the check to only verify PyTorch is available**

Change from:
```r
"cavae" = {
  torch_ok <- tryCatch({
    torch <- reticulate::import("torch", convert = FALSE)
    torch$cuda$is_available()
  }, error = function(e) FALSE)
  if (!torch_ok) {
    return("PyTorch with CUDA not available")
  }
```

To:
```r
"cavae" = {
  torch_ok <- tryCatch({
    torch <- reticulate::import("torch", convert = FALSE)
    TRUE
  }, error = function(e) FALSE)
  if (!torch_ok) {
    return("PyTorch not available")
  }
```

**Step 3: Update documentation**

```bash
Rscript -e "devtools::document('packages/cfomics')"
```

**Step 4: Run tests**

```bash
Rscript -e "devtools::test('packages/cfomics', filter = 'cavae')"
```

**Step 5: Commit**

```bash
git add packages/cfomics/R/benchmark_runner.R
git commit -m "fix(benchmark): remove CUDA requirement from CAVAE dependency check

CAVAE can run on CPU; only PyTorch availability should be checked."
```

---

## Phase 3: Medium Priority Fixes

### Task 7: Add validation to missing data DGP (S15)

**Files:**
- Modify: `packages/cfomics/R/benchmark_dgp.R` (dgp_missing_data function)
- Test: `packages/cfomics/tests/testthat/test-benchmark-dgp.R`

**Step 1: Add validation to dgp_missing_data**

After generating missing data, add validation:

```r
# Validate missing rate
actual_missing_rate <- mean(is.na(X))
if (abs(actual_missing_rate - missing_rate) > 0.1) {
  warning(sprintf("Actual missing rate (%.2f) differs from requested (%.2f)",
                  actual_missing_rate, missing_rate))
}

# Validate dimensions match
stopifnot(
  identical(dim(X), dim(X_complete)),
  all(!is.na(X_complete))
)
```

**Step 2: Add test for missing data DGP**

Add to test-benchmark-dgp.R:

```r
test_that("dgp_missing_data generates correct missing pattern", {
  result <- dgp_missing_data(n = 200, p = 10, missing_rate = 0.2, seed = 123)


  # Check structure
expect_true("X_complete" %in% names(result))
  expect_true("X" %in% names(result))
  expect_equal(dim(result$X), dim(result$X_complete))

  # Check missing values exist in X but not in X_complete
  expect_true(any(is.na(result$X)))
  expect_false(any(is.na(result$X_complete)))

  # Check approximate missing rate (within tolerance)
  actual_rate <- mean(is.na(result$X))
  expect_true(abs(actual_rate - 0.2) < 0.1)
})
```

**Step 3: Run tests**

```bash
Rscript -e "devtools::test('packages/cfomics', filter = 'dgp')"
```

**Step 4: Commit**

```bash
git add packages/cfomics/R/benchmark_dgp.R packages/cfomics/tests/testthat/test-benchmark-dgp.R
git commit -m "feat(benchmark): add validation and tests for missing data DGP (S15)"
```

---

### Task 8: Add input validation to sweep functions

**Files:**
- Modify: `packages/cfomics/R/benchmark_sweep.R`

**Step 1: Add validation helper function**

Add at top of file:

```r
.validate_sweep_params <- function(n_values = NULL, p_values = NULL,
                                    n_reps = NULL, seed = NULL) {
  if (!is.null(n_values)) {
    stopifnot(
      is.numeric(n_values),
      all(n_values > 0),
      all(n_values == as.integer(n_values))
    )
  }
  if (!is.null(p_values)) {
    stopifnot(
      is.numeric(p_values),
      all(p_values > 0),
      all(p_values == as.integer(p_values))
    )
  }
  if (!is.null(n_reps)) {
    stopifnot(
      is.numeric(n_reps),
      length(n_reps) == 1,
      n_reps > 0,
      n_reps == as.integer(n_reps)
    )
  }
  if (!is.null(seed)) {
    stopifnot(
      is.numeric(seed),
      length(seed) == 1
    )
  }
  invisible(TRUE)
}
```

**Step 2: Add validation calls to sweep functions**

Add at start of `cf_benchmark_dimension_sweep()`:
```r
.validate_sweep_params(n_values = n_values, p_values = p_values,
                       n_reps = n_reps, seed = seed)
```

Add at start of `cf_benchmark_heterogeneity_sweep()`:
```r
.validate_sweep_params(n_reps = n_reps, seed = seed)
stopifnot(is.numeric(strength_values), all(strength_values >= 0))
```

**Step 3: Run tests**

```bash
Rscript -e "devtools::test('packages/cfomics', filter = 'sweep')"
```

**Step 4: Commit**

```bash
git add packages/cfomics/R/benchmark_sweep.R
git commit -m "feat(benchmark): add input validation to sweep functions"
```

---

## Summary

| Task | Phase | Priority | Description |
|------|-------|----------|-------------|
| 1 | 1 | Critical | Fix mse_ate in run_nn_benchmark.R |
| 2 | 1 | Critical | Fix mse_ate in aggregate_results.R |
| 3 | 1 | Critical | Fix mse_ate in fig2_simulation.R |
| 4 | 1 | Critical | Fix mse_ate in fig3_highdim.R |
| 5 | 1 | Critical | Fix mse_ate in fig4_tcga.R |
| 6 | 2 | High | Fix CAVAE CUDA requirement |
| 7 | 3 | Medium | Add missing data DGP validation |
| 8 | 3 | Medium | Add sweep input validation |

Total: **8 Tasks** (5 Critical, 1 High, 2 Medium)
