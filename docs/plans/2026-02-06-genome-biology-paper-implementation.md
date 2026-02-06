# cfomics Genome Biology Paper Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Prepare cfomics for Genome Biology paper submission with comprehensive benchmarking, external package comparison, TCGA validation, and publication figures.

**Architecture:**
- External package wrappers in `benchmarks/external/wrappers/` follow unified interface returning `list(ate, ite, y0_hat, y1_hat, ci_lower, ci_upper)`
- Wrappers are registered in `benchmarks/external/registry.R` and integrated with existing `benchmarks/R/runner.R`
- TCGA semi-synthetic validation uses real covariates with synthetic outcomes
- Neural network methods (TARNet, CFRNet, DragonNet) are comparison-only scripts, not part of cfomics package

**Tech Stack:** R (cfomics, MatchIt, WeightIt, hdm, DoubleML, tmle3, BART), Python (PyTorch for NN comparisons), future/furrr for parallel execution

---

## Phase 1: External Package Wrappers

### Task 1: Create wrapper directory structure and registry

**Files:**
- Create: `benchmarks/external/wrappers/.gitkeep`
- Create: `benchmarks/external/registry.R`

**Step 1: Create directory structure**

```bash
mkdir -p benchmarks/external/wrappers
touch benchmarks/external/wrappers/.gitkeep
```

**Step 2: Write registry.R with unified interface definition**

```r
# benchmarks/external/registry.R
# Registry for external package wrappers

#' Run external method wrapper
#'
#' Unified interface for external causal inference packages.
#' All wrappers return the same structure for fair comparison.
#'
#' @param method Character, external method name
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param ... Additional method-specific arguments
#' @return List with: ate, ite, y0_hat, y1_hat, ci_lower, ci_upper, time_sec
run_external_method <- function(method, X, T, Y, ...) {
  wrapper_fn <- get_external_wrapper(method)
  if (is.null(wrapper_fn)) {
    stop(sprintf("Unknown external method: %s", method))
  }

  start_time <- Sys.time()
  result <- wrapper_fn(X = X, T = T, Y = Y, ...)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Validate result structure
  required_fields <- c("ate", "ite")
  missing <- setdiff(required_fields, names(result))
  if (length(missing) > 0) {
    stop(sprintf("Wrapper for %s missing required fields: %s",
                 method, paste(missing, collapse = ", ")))
  }

  # Add defaults for optional fields
  result$y0_hat <- result$y0_hat %||% rep(NA_real_, length(Y))
  result$y1_hat <- result$y1_hat %||% rep(NA_real_, length(Y))
  result$ci_lower <- result$ci_lower %||% NA_real_
  result$ci_upper <- result$ci_upper %||% NA_real_

  result$time_sec <- elapsed


  result
}

#' Get wrapper function for external method
#' @param method Character, method name
#' @return Function or NULL if not found
get_external_wrapper <- function(method) {
  wrappers <- list(
    "matchit" = wrapper_matchit,
    "weightit" = wrapper_weightit,
    "hdm" = wrapper_hdm,
    "doubleml" = wrapper_doubleml,
    "tmle3" = wrapper_tmle3,
    "bart" = wrapper_bart,
    "superlearner" = wrapper_superlearner

  )
  wrappers[[method]]
}

#' List available external methods
#' @return Character vector of method names
list_external_methods <- function() {
  c("matchit", "weightit", "hdm", "doubleml", "tmle3", "bart", "superlearner")
}

#' Check if external method dependencies are available
#' @param method Character, method name
#' @return NULL if available, character string with reason if not
check_external_dependencies <- function(method) {
  deps <- list(
    matchit = c("MatchIt", "marginaleffects"),
    weightit = c("WeightIt", "marginaleffects"),
    hdm = c("hdm"),
    doubleml = c("DoubleML", "mlr3", "mlr3learners"),
    tmle3 = c("tmle3", "sl3"),
    bart = c("dbarts"),
    superlearner = c("SuperLearner")
  )

  required <- deps[[method]]
  if (is.null(required)) return(sprintf("Unknown method: %s", method))

  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    return(sprintf("Missing packages: %s", paste(missing, collapse = ", ")))
  }
  NULL
}

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
```

**Step 3: Commit**

```bash
git add benchmarks/external/
git commit -m "feat(benchmark): add external wrapper registry structure"
```

---

### Task 2: Implement MatchIt wrapper

**Files:**
- Create: `benchmarks/external/wrappers/wrapper_matchit.R`
- Create: `benchmarks/external/wrappers/tests/test_matchit.R`

**Step 1: Write the test**

```r
# benchmarks/external/wrappers/tests/test_matchit.R
# Test MatchIt wrapper

test_that("wrapper_matchit returns correct structure", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("marginaleffects")

  source(here::here("benchmarks/external/wrappers/wrapper_matchit.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_matchit(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(is.numeric(result$ate))
  expect_true(abs(result$ate - 2) < 1)  # Rough check
})
```

**Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_matchit.R')"
```

Expected: FAIL with "could not find function wrapper_matchit"

**Step 3: Write wrapper implementation**

```r
# benchmarks/external/wrappers/wrapper_matchit.R
# MatchIt wrapper for benchmark comparison

#' MatchIt wrapper
#'
#' Uses propensity score matching via MatchIt package.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param method Matching method (default: "nearest")
#' @param estimand Estimand (default: "ATE")
#' @return List with ate, ite, ci_lower, ci_upper
wrapper_matchit <- function(X, T, Y, method = "nearest", estimand = "ATE") {
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("Package 'MatchIt' required")
  }
  if (!requireNamespace("marginaleffects", quietly = TRUE)) {
    stop("Package 'marginaleffects' required")
  }

  # Prepare data
  df <- data.frame(Y = Y, T = T, X)
  cov_names <- colnames(X)
  if (is.null(cov_names)) cov_names <- paste0("X", seq_len(ncol(X)))
  colnames(df) <- c("Y", "T", cov_names)

  # Build formula
  ps_formula <- as.formula(paste("T ~", paste(cov_names, collapse = " + ")))

  # Run matching
  m_out <- MatchIt::matchit(
    ps_formula,
    data = df,
    method = method,
    estimand = estimand
  )

  # Get matched data
  m_data <- MatchIt::match.data(m_out)

  # Fit outcome model on matched data
  out_formula <- as.formula(paste("Y ~ T +", paste(cov_names, collapse = " + ")))
  fit <- lm(out_formula, data = m_data, weights = weights)

  # Get ATE with CI using marginaleffects
  ate_est <- marginaleffects::avg_comparisons(
    fit,
    variables = "T",
    vcov = "HC3",
    newdata = m_data,
    wts = "weights"
  )

  ate <- ate_est$estimate
  ci_lower <- ate_est$conf.low
  ci_upper <- ate_est$conf.high

  # ITE: predict Y(1) - Y(0) for each unit
  df0 <- df1 <- df

  df0$T <- 0
  df1$T <- 1
  y0_hat <- predict(fit, newdata = df0)
  y1_hat <- predict(fit, newdata = df1)
  ite <- y1_hat - y0_hat

  list(
    ate = ate,
    ite = as.numeric(ite),
    y0_hat = as.numeric(y0_hat),
    y1_hat = as.numeric(y1_hat),
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}
```

**Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_matchit.R')"
```

Expected: PASS

**Step 5: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_matchit.R benchmarks/external/wrappers/tests/
git commit -m "feat(benchmark): add MatchIt wrapper"
```

---

### Task 3: Implement WeightIt wrapper

**Files:**
- Create: `benchmarks/external/wrappers/wrapper_weightit.R`
- Create: `benchmarks/external/wrappers/tests/test_weightit.R`

**Step 1: Write the test**

```r
# benchmarks/external/wrappers/tests/test_weightit.R
test_that("wrapper_weightit returns correct structure", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("marginaleffects")

  source(here::here("benchmarks/external/wrappers/wrapper_weightit.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_weightit(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1)
})
```

**Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_weightit.R')"
```

**Step 3: Write wrapper implementation**

```r
# benchmarks/external/wrappers/wrapper_weightit.R
# WeightIt wrapper for benchmark comparison

#' WeightIt wrapper
#'
#' Uses inverse probability weighting via WeightIt package.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param method Weighting method (default: "glm")
#' @param estimand Estimand (default: "ATE")
#' @return List with ate, ite, ci_lower, ci_upper
wrapper_weightit <- function(X, T, Y, method = "glm", estimand = "ATE") {
  if (!requireNamespace("WeightIt", quietly = TRUE)) {
    stop("Package 'WeightIt' required")
  }
  if (!requireNamespace("marginaleffects", quietly = TRUE)) {
    stop("Package 'marginaleffects' required")
  }

  # Prepare data
  df <- data.frame(Y = Y, T = T, X)
  cov_names <- colnames(X)
  if (is.null(cov_names)) cov_names <- paste0("X", seq_len(ncol(X)))
  colnames(df) <- c("Y", "T", cov_names)

  # Build formula
  ps_formula <- as.formula(paste("T ~", paste(cov_names, collapse = " + ")))

  # Get weights
  w_out <- WeightIt::weightit(
    ps_formula,
    data = df,
    method = method,
    estimand = estimand
  )

  df$weights <- w_out$weights

  # Fit weighted outcome model
  out_formula <- as.formula(paste("Y ~ T +", paste(cov_names, collapse = " + ")))
  fit <- lm(out_formula, data = df, weights = weights)

  # Get ATE with CI
  ate_est <- marginaleffects::avg_comparisons(
    fit,
    variables = "T",
    vcov = "HC3",
    newdata = df,
    wts = "weights"
  )

  ate <- ate_est$estimate
  ci_lower <- ate_est$conf.low
  ci_upper <- ate_est$conf.high

  # ITE
  df0 <- df1 <- df
  df0$T <- 0
  df1$T <- 1
  y0_hat <- predict(fit, newdata = df0)
  y1_hat <- predict(fit, newdata = df1)
  ite <- y1_hat - y0_hat

  list(
    ate = ate,
    ite = as.numeric(ite),
    y0_hat = as.numeric(y0_hat),
    y1_hat = as.numeric(y1_hat),
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}
```

**Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_weightit.R')"
```

**Step 5: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_weightit.R benchmarks/external/wrappers/tests/test_weightit.R
git commit -m "feat(benchmark): add WeightIt wrapper"
```

---

### Task 4: Implement hdm wrapper

**Files:**
- Create: `benchmarks/external/wrappers/wrapper_hdm.R`
- Create: `benchmarks/external/wrappers/tests/test_hdm.R`

**Step 1: Write the test**

```r
# benchmarks/external/wrappers/tests/test_hdm.R
test_that("wrapper_hdm returns correct structure", {
  skip_if_not_installed("hdm")

  source(here::here("benchmarks/external/wrappers/wrapper_hdm.R"))

  set.seed(123)
  n <- 200
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  ps <- plogis(0.5 * X[,1] + 0.3 * X[,2])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_hdm(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1.5)
})
```

**Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_hdm.R')"
```

**Step 3: Write wrapper implementation**

```r
# benchmarks/external/wrappers/wrapper_hdm.R
# hdm (high-dimensional metrics) wrapper for benchmark comparison

#' hdm wrapper
#'
#' Uses double selection Lasso via hdm package for high-dimensional data.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @return List with ate, ite, ci_lower, ci_upper
wrapper_hdm <- function(X, T, Y) {
  if (!requireNamespace("hdm", quietly = TRUE)) {
    stop("Package 'hdm' required")
  }

  # hdm::rlassoEffect for treatment effect with double selection
  fit <- hdm::rlassoEffect(
    x = X,
    y = Y,
    d = T,
    method = "double selection"
  )

  # Extract ATE and CI
  summ <- summary(fit)
  ate <- summ$coefficients["d1", "Estimate"]
  se <- summ$coefficients["d1", "Std. Error"]
  ci_lower <- ate - 1.96 * se
  ci_upper <- ate + 1.96 * se

  # ITE: hdm doesn't provide heterogeneous effects, use constant
  n <- length(Y)
  ite <- rep(ate, n)

  list(
    ate = ate,
    ite = ite,
    y0_hat = rep(NA_real_, n),
    y1_hat = rep(NA_real_, n),
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}
```

**Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_hdm.R')"
```

**Step 5: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_hdm.R benchmarks/external/wrappers/tests/test_hdm.R
git commit -m "feat(benchmark): add hdm wrapper"
```

---

### Task 5: Implement DoubleML wrapper

**Files:**
- Create: `benchmarks/external/wrappers/wrapper_doubleml.R`
- Create: `benchmarks/external/wrappers/tests/test_doubleml.R`

**Step 1: Write the test**

```r
# benchmarks/external/wrappers/tests/test_doubleml.R
test_that("wrapper_doubleml returns correct structure", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3learners")

  source(here::here("benchmarks/external/wrappers/wrapper_doubleml.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_doubleml(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1)
})
```

**Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_doubleml.R')"
```

**Step 3: Write wrapper implementation**

```r
# benchmarks/external/wrappers/wrapper_doubleml.R
# DoubleML wrapper for benchmark comparison

#' DoubleML wrapper
#'
#' Uses Double Machine Learning via DoubleML package.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param ml_g Learner for outcome model (default: "regr.ranger")
#' @param ml_m Learner for propensity model (default: "classif.ranger")
#' @return List with ate, ite, ci_lower, ci_upper
wrapper_doubleml <- function(X, T, Y, ml_g = "regr.ranger", ml_m = "classif.ranger") {
  if (!requireNamespace("DoubleML", quietly = TRUE)) {
    stop("Package 'DoubleML' required")
  }
  if (!requireNamespace("mlr3", quietly = TRUE)) {
    stop("Package 'mlr3' required")
  }
  if (!requireNamespace("mlr3learners", quietly = TRUE)) {
    stop("Package 'mlr3learners' required")
  }

  # Prepare data
  df <- data.frame(Y = Y, T = T, X)
  cov_names <- colnames(X)
  if (is.null(cov_names)) cov_names <- paste0("X", seq_len(ncol(X)))
  colnames(df) <- c("Y", "T", cov_names)

  # Create DoubleML data object
  dml_data <- DoubleML::DoubleMLData$new(
    df,
    y_col = "Y",
    d_cols = "T",
    x_cols = cov_names
  )

  # Create learners
  lgr::get_logger("mlr3")$set_threshold("warn")
  ml_g_obj <- mlr3::lrn(ml_g)
  ml_m_obj <- mlr3::lrn(ml_m)

  # Fit PLR model
  dml_plr <- DoubleML::DoubleMLPLR$new(
    dml_data,
    ml_g = ml_g_obj,
    ml_m = ml_m_obj
  )
  dml_plr$fit()

  # Extract results
  ate <- dml_plr$coef
  se <- dml_plr$se
  ci_lower <- ate - 1.96 * se
  ci_upper <- ate + 1.96 * se

  # ITE: PLR gives constant effect
  n <- length(Y)
  ite <- rep(ate, n)

  list(
    ate = as.numeric(ate),
    ite = as.numeric(ite),
    y0_hat = rep(NA_real_, n),
    y1_hat = rep(NA_real_, n),
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper)
  )
}
```

**Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_doubleml.R')"
```

**Step 5: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_doubleml.R benchmarks/external/wrappers/tests/test_doubleml.R
git commit -m "feat(benchmark): add DoubleML wrapper"
```

---

### Task 6: Implement tmle3 wrapper

**Files:**
- Create: `benchmarks/external/wrappers/wrapper_tmle3.R`
- Create: `benchmarks/external/wrappers/tests/test_tmle3.R`

**Step 1: Write the test**

```r
# benchmarks/external/wrappers/tests/test_tmle3.R
test_that("wrapper_tmle3 returns correct structure", {
  skip_if_not_installed("tmle3")
  skip_if_not_installed("sl3")

  source(here::here("benchmarks/external/wrappers/wrapper_tmle3.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_tmle3(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1)
})
```

**Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_tmle3.R')"
```

**Step 3: Write wrapper implementation**

```r
# benchmarks/external/wrappers/wrapper_tmle3.R
# tmle3 (tlverse) wrapper for benchmark comparison

#' tmle3 wrapper
#'
#' Uses Targeted Learning via tmle3/tlverse packages.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @return List with ate, ite, ci_lower, ci_upper
wrapper_tmle3 <- function(X, T, Y) {
  if (!requireNamespace("tmle3", quietly = TRUE)) {
    stop("Package 'tmle3' required")
  }
  if (!requireNamespace("sl3", quietly = TRUE)) {
    stop("Package 'sl3' required")
  }

  # Prepare data
  df <- data.frame(Y = Y, A = T, X)
  cov_names <- colnames(X)
  if (is.null(cov_names)) cov_names <- paste0("X", seq_len(ncol(X)))
  colnames(df) <- c("Y", "A", cov_names)

  # Define node list
  node_list <- list(
    W = cov_names,
    A = "A",
    Y = "Y"
  )

  # Simple learner stack
  lrnr_glm <- sl3::Lrnr_glm$new()
  lrnr_mean <- sl3::Lrnr_mean$new()
  sl <- sl3::Lrnr_sl$new(learners = list(lrnr_glm, lrnr_mean))

  learner_list <- list(A = sl, Y = sl)

  # Define spec for ATE
  ate_spec <- tmle3::tmle_ATE(
    treatment_level = 1,
    control_level = 0
  )

  # Fit TMLE
  tmle_fit <- tmle3::tmle3(
    ate_spec,
    data = df,
    node_list = node_list,
    learner_list = learner_list
  )

  # Extract results
  summ <- tmle_fit$summary
  ate <- summ$psi
  ci_lower <- summ$lower
  ci_upper <- summ$upper

  # ITE: TMLE gives population average, not individual effects
  n <- length(Y)
  ite <- rep(ate, n)

  list(
    ate = ate,
    ite = ite,
    y0_hat = rep(NA_real_, n),
    y1_hat = rep(NA_real_, n),
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}
```

**Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_tmle3.R')"
```

**Step 5: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_tmle3.R benchmarks/external/wrappers/tests/test_tmle3.R
git commit -m "feat(benchmark): add tmle3 wrapper"
```

---

### Task 7: Implement BART wrapper

**Files:**
- Create: `benchmarks/external/wrappers/wrapper_bart.R`
- Create: `benchmarks/external/wrappers/tests/test_bart.R`

**Step 1: Write the test**

```r
# benchmarks/external/wrappers/tests/test_bart.R
test_that("wrapper_bart returns correct structure", {
  skip_if_not_installed("dbarts")

  source(here::here("benchmarks/external/wrappers/wrapper_bart.R"))

  set.seed(123)
  n <- 200
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("X", 1:5)
  ps <- plogis(0.5 * X[,1])
  T <- rbinom(n, 1, ps)
  Y <- X[,1] + 2 * T + rnorm(n)

  result <- wrapper_bart(X, T, Y)

  expect_type(result, "list")
  expect_true("ate" %in% names(result))
  expect_true("ite" %in% names(result))
  expect_length(result$ite, n)
  expect_true(abs(result$ate - 2) < 1)
})
```

**Step 2: Run test to verify it fails**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_bart.R')"
```

**Step 3: Write wrapper implementation**

```r
# benchmarks/external/wrappers/wrapper_bart.R
# BART (dbarts) wrapper for benchmark comparison

#' BART wrapper
#'
#' Uses Bayesian Additive Regression Trees via dbarts package.
#'
#' @param X Numeric matrix of covariates
#' @param T Integer vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param n_trees Number of trees (default: 200)
#' @param n_burn Burn-in iterations (default: 250)
#' @param n_iter Post-burn iterations (default: 1000)
#' @return List with ate, ite, y0_hat, y1_hat, ci_lower, ci_upper
wrapper_bart <- function(X, T, Y, n_trees = 200L, n_burn = 250L, n_iter = 1000L) {
  if (!requireNamespace("dbarts", quietly = TRUE)) {
    stop("Package 'dbarts' required")
  }

  n <- length(Y)

  # Combine X and T for BART
  X_with_T <- cbind(T = T, X)

  # Fit BART
  bart_fit <- dbarts::bart(
    x.train = X_with_T,
    y.train = Y,
    ntree = n_trees,
    nskip = n_burn,
    ndpost = n_iter,
    verbose = FALSE
  )

  # Counterfactual predictions
  X_t0 <- cbind(T = 0, X)
  X_t1 <- cbind(T = 1, X)

  # Posterior predictions (matrix: n_iter x n)
  pred_t0 <- predict(bart_fit, newdata = X_t0)
  pred_t1 <- predict(bart_fit, newdata = X_t1)

  # Point estimates (posterior mean)
  y0_hat <- colMeans(pred_t0)
  y1_hat <- colMeans(pred_t1)
  ite <- y1_hat - y0_hat
  ate <- mean(ite)

  # CI from posterior
  ite_samples <- pred_t1 - pred_t0
  ate_samples <- rowMeans(ite_samples)
  ci_lower <- quantile(ate_samples, 0.025)
  ci_upper <- quantile(ate_samples, 0.975)

  list(
    ate = ate,
    ite = as.numeric(ite),
    y0_hat = as.numeric(y0_hat),
    y1_hat = as.numeric(y1_hat),
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper)
  )
}
```

**Step 4: Run test to verify it passes**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_bart.R')"
```

**Step 5: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_bart.R benchmarks/external/wrappers/tests/test_bart.R
git commit -m "feat(benchmark): add BART wrapper"
```

---

### Task 8: Integrate external wrappers into benchmark runner

**Files:**
- Modify: `benchmarks/R/runner.R`
- Modify: `benchmarks/config.R`

**Step 1: Update runner.R to support external methods**

Add after line 91 in `benchmarks/R/runner.R`:

```r
# Run a single benchmark job for EXTERNAL method
run_single_job_external <- function(scenario, method, rep, cfg) {
  gen <- generate_scenario_data(scenario, rep, cfg$base_seed)

  # Source external wrappers
  external_dir <- file.path(dirname(cfg$scenarios_path), "..", "external")
  source(file.path(external_dir, "registry.R"), local = TRUE)
  source(file.path(external_dir, "wrappers", paste0("wrapper_", method, ".R")), local = TRUE)

  result <- tryCatch({
    ext_result <- run_external_method(
      method = method,
      X = gen$data[, grep("^X", names(gen$data))],
      T = gen$data$T,
      Y = gen$data$Y
    )

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
      mse_ate = metrics$mse_ate,
      pehe = metrics$pehe,
      coverage_ate = metrics$coverage_ate,
      ci_len_ate = metrics$ci_len_ate,
      time_sec = ext_result$time_sec,
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
```

**Step 2: Update config.R to include external methods**

Add to `benchmark_config()` in `benchmarks/config.R`:

```r
    # External methods to benchmark (prefix "ext_" added automatically)
    external_methods = c("matchit", "weightit", "hdm", "doubleml", "bart"),
```

**Step 3: Commit**

```bash
git add benchmarks/R/runner.R benchmarks/config.R
git commit -m "feat(benchmark): integrate external wrappers into runner"
```

---

## Phase 2: Neural Network Comparison Scripts

### Task 9: Create TARNet/CFRNet/DragonNet comparison scripts

**Files:**
- Create: `benchmarks/external/nn/requirements.txt`
- Create: `benchmarks/external/nn/models.py`
- Create: `benchmarks/external/nn/run_nn_benchmark.py`
- Create: `benchmarks/external/nn/run_nn_benchmark.R`

**Step 1: Create Python requirements**

```txt
# benchmarks/external/nn/requirements.txt
torch>=2.0.0
numpy>=1.24.0
pandas>=2.0.0
scikit-learn>=1.3.0
```

**Step 2: Create Python models implementation**

```python
# benchmarks/external/nn/models.py
"""Neural network models for causal inference: TARNet, CFRNet, DragonNet"""

import torch
import torch.nn as nn
import numpy as np
from typing import Tuple, Optional


class RepresentationNet(nn.Module):
    """Shared representation network."""

    def __init__(self, input_dim: int, hidden_dim: int = 200, n_layers: int = 3):
        super().__init__()
        layers = [nn.Linear(input_dim, hidden_dim), nn.ELU()]
        for _ in range(n_layers - 1):
            layers.extend([nn.Linear(hidden_dim, hidden_dim), nn.ELU()])
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class TARNet(nn.Module):
    """Treatment-Agnostic Representation Network."""

    def __init__(self, input_dim: int, hidden_dim: int = 200, n_layers: int = 3):
        super().__init__()
        self.repr = RepresentationNet(input_dim, hidden_dim, n_layers)
        self.head0 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.ELU(),
            nn.Linear(hidden_dim, 1)
        )
        self.head1 = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.ELU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x: torch.Tensor, t: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        phi = self.repr(x)
        y0 = self.head0(phi).squeeze(-1)
        y1 = self.head1(phi).squeeze(-1)
        y_pred = t * y1 + (1 - t) * y0
        return y_pred, (y0, y1)


class CFRNet(TARNet):
    """Counterfactual Regression Network with IPM regularization."""

    def __init__(self, input_dim: int, hidden_dim: int = 200, n_layers: int = 3, alpha: float = 1.0):
        super().__init__(input_dim, hidden_dim, n_layers)
        self.alpha = alpha

    def ipm_loss(self, phi: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        """Wasserstein-1 (MMD) approximation for IPM."""
        t = t.bool()
        phi0 = phi[~t]
        phi1 = phi[t]
        if len(phi0) == 0 or len(phi1) == 0:
            return torch.tensor(0.0, device=phi.device)
        mean0 = phi0.mean(dim=0)
        mean1 = phi1.mean(dim=0)
        return torch.norm(mean0 - mean1, p=2)


class DragonNet(nn.Module):
    """DragonNet: propensity score co-learning."""

    def __init__(self, input_dim: int, hidden_dim: int = 200, n_layers: int = 3):
        super().__init__()
        self.repr = RepresentationNet(input_dim, hidden_dim, n_layers)
        self.head0 = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.ELU(), nn.Linear(hidden_dim, 1))
        self.head1 = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.ELU(), nn.Linear(hidden_dim, 1))
        self.propensity_head = nn.Sequential(nn.Linear(hidden_dim, 1), nn.Sigmoid())

    def forward(self, x: torch.Tensor, t: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        phi = self.repr(x)
        y0 = self.head0(phi).squeeze(-1)
        y1 = self.head1(phi).squeeze(-1)
        ps = self.propensity_head(phi).squeeze(-1)
        y_pred = t * y1 + (1 - t) * y0
        return y_pred, (y0, y1), ps


def train_model(model, X, T, Y, n_epochs=300, lr=1e-3, batch_size=64, alpha=1.0):
    """Train a causal NN model."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    X_t = torch.tensor(X, dtype=torch.float32, device=device)
    T_t = torch.tensor(T, dtype=torch.float32, device=device)
    Y_t = torch.tensor(Y, dtype=torch.float32, device=device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    n = len(Y)

    for epoch in range(n_epochs):
        perm = torch.randperm(n)
        for i in range(0, n, batch_size):
            idx = perm[i:i+batch_size]
            x_b, t_b, y_b = X_t[idx], T_t[idx], Y_t[idx]

            optimizer.zero_grad()

            if isinstance(model, DragonNet):
                y_pred, (y0, y1), ps = model(x_b, t_b)
                loss_y = nn.MSELoss()(y_pred, y_b)
                loss_ps = nn.BCELoss()(ps, t_b)
                loss = loss_y + loss_ps
            elif isinstance(model, CFRNet):
                y_pred, (y0, y1) = model(x_b, t_b)
                loss_y = nn.MSELoss()(y_pred, y_b)
                phi = model.repr(x_b)
                loss_ipm = model.ipm_loss(phi, t_b)
                loss = loss_y + alpha * loss_ipm
            else:  # TARNet
                y_pred, (y0, y1) = model(x_b, t_b)
                loss = nn.MSELoss()(y_pred, y_b)

            loss.backward()
            optimizer.step()

    return model


def predict_ite(model, X):
    """Predict ITE using trained model."""
    device = next(model.parameters()).device
    model.eval()

    with torch.no_grad():
        X_t = torch.tensor(X, dtype=torch.float32, device=device)
        T0 = torch.zeros(len(X), device=device)
        T1 = torch.ones(len(X), device=device)

        if isinstance(model, DragonNet):
            _, (y0, _), _ = model(X_t, T0)
            _, (_, y1), _ = model(X_t, T1)
        else:
            _, (y0, y1) = model(X_t, T0)

        ite = (y1 - y0).cpu().numpy()
        y0_hat = y0.cpu().numpy()
        y1_hat = y1.cpu().numpy()

    return {
        "ite": ite,
        "ate": float(np.mean(ite)),
        "y0_hat": y0_hat,
        "y1_hat": y1_hat
    }
```

**Step 3: Create R wrapper for NN benchmark**

```r
# benchmarks/external/nn/run_nn_benchmark.R
# R interface to run NN benchmarks via reticulate

run_nn_benchmark <- function(scenario, method, rep, cfg, nn_config = list()) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' required for NN benchmarks")
  }

  # Source scenarios
  source(cfg$scenarios_path, local = TRUE)
  gen <- generate_scenario_data(scenario, rep, cfg$base_seed)

  # Extract data
  X <- as.matrix(gen$data[, grep("^X", names(gen$data))])
  T <- gen$data$T
  Y <- gen$data$Y

  # Run Python
  nn_dir <- file.path(dirname(cfg$scenarios_path), "..", "external", "nn")
  reticulate::source_python(file.path(nn_dir, "models.py"))

  model_class <- switch(method,
    "tarnet" = TARNet,
    "cfrnet" = CFRNet,
    "dragonnet" = DragonNet,
    stop("Unknown NN method: ", method)
  )

  input_dim <- ncol(X)
  hidden_dim <- nn_config$hidden_dim %||% 200L
  n_layers <- nn_config$n_layers %||% 3L
  n_epochs <- nn_config$n_epochs %||% 300L

  start_time <- Sys.time()

  model <- model_class(input_dim, hidden_dim, n_layers)
  model <- train_model(model, X, T, Y, n_epochs = n_epochs)
  result <- predict_ite(model, X)

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Compute metrics
  metrics <- cfomics:::cf_benchmark_compute_metrics(
    ate_hat = result$ate,
    ite_hat = result$ite,
    summary_hat = NULL,
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
    ate_hat = result$ate,
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
}

`%||%` <- function(x, y) if (is.null(x)) y else x
```

**Step 4: Commit**

```bash
git add benchmarks/external/nn/
git commit -m "feat(benchmark): add TARNet/CFRNet/DragonNet comparison scripts"
```

---

## Phase 3: TCGA Validation

### Task 10: Create TCGA data preparation pipeline

**Files:**
- Create: `benchmarks/tcga/download_tcga.R`
- Create: `benchmarks/tcga/preprocess_tcga.R`

**Step 1: Write TCGA download script**

```r
# benchmarks/tcga/download_tcga.R
# Download TCGA data using TCGAbiolinks

download_tcga_data <- function(project = "TCGA-BRCA", output_dir = "benchmarks/tcga/data") {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Package 'TCGAbiolinks' required. Install from Bioconductor.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' required.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Query RNA-seq data
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # Download
  TCGAbiolinks::GDCdownload(query, directory = file.path(output_dir, "GDCdata"))


  # Prepare data
  se <- TCGAbiolinks::GDCprepare(query, directory = file.path(output_dir, "GDCdata"))

  # Save
  saveRDS(se, file.path(output_dir, paste0(project, "_se.rds")))

  message(sprintf("Downloaded %s: %d samples, %d genes",
                  project, ncol(se), nrow(se)))

  invisible(se)
}

# Download multiple projects
if (sys.nframe() == 0) {
  projects <- c("TCGA-BRCA", "TCGA-LUAD", "TCGA-COAD")
  for (proj in projects) {
    message(sprintf("\n=== Downloading %s ===", proj))
    tryCatch(
      download_tcga_data(proj),
      error = function(e) message(sprintf("Error downloading %s: %s", proj, e$message))
    )
  }
}
```

**Step 2: Write preprocessing script**

```r
# benchmarks/tcga/preprocess_tcga.R
# Preprocess TCGA data for semi-synthetic benchmarks

preprocess_tcga <- function(se, n_genes = 1000, variance_filter = TRUE) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' required")
  }

  # Get expression matrix (log2 + 1 normalized)
  counts <- SummarizedExperiment::assay(se, "unstranded")
  expr <- log2(counts + 1)

  # Filter genes by variance
  if (variance_filter && n_genes < nrow(expr)) {
    gene_vars <- apply(expr, 1, var, na.rm = TRUE)
    top_genes <- order(gene_vars, decreasing = TRUE)[1:n_genes]
    expr <- expr[top_genes, ]
  }

  # Transpose: samples x genes
  X <- t(expr)

  # Get clinical data
  clinical <- as.data.frame(SummarizedExperiment::colData(se))

  # Extract useful variables
  result <- list(
    X = X,
    sample_ids = rownames(X),
    gene_ids = colnames(X),
    clinical = clinical
  )

  # Add TP53 mutation status if available
  if ("paper_TP53_mut_status" %in% names(clinical)) {
    result$tp53_mut <- as.integer(clinical$paper_TP53_mut_status == "Mut")
  }

  result
}

# Create semi-synthetic benchmark data
create_semi_synthetic_data <- function(tcga_data, dgp_params = list()) {
  X <- tcga_data$X
  n <- nrow(X)
  p <- ncol(X)

  # Use PCA for dimension reduction if needed
  if (p > 100) {
    pca <- prcomp(X, center = TRUE, scale. = TRUE)
    X_reduced <- pca$x[, 1:100]
  } else {
    X_reduced <- scale(X)
  }

  # Synthetic treatment assignment
  tau_true <- dgp_params$tau %||% 2.0
  ps_coef <- dgp_params$ps_coef %||% c(0.5, 0.3, -0.2, rep(0, ncol(X_reduced) - 3))
  ps <- plogis(X_reduced %*% ps_coef[1:ncol(X_reduced)])
  T <- rbinom(n, 1, ps)

  # Synthetic outcome
  outcome_coef <- dgp_params$outcome_coef %||% c(1.0, -0.5, 0.3, rep(0, ncol(X_reduced) - 3))
  mu0 <- X_reduced %*% outcome_coef[1:ncol(X_reduced)]

  # Heterogeneous treatment effect
  if (dgp_params$heterogeneous %||% FALSE) {
    tau <- tau_true * (1 + 0.3 * X_reduced[, 1])
  } else {
    tau <- rep(tau_true, n)
  }

  Y <- mu0 + T * tau + rnorm(n, sd = dgp_params$noise_sd %||% 1.0)

  list(
    data = data.frame(Y = as.numeric(Y), T = T, X_reduced),
    truth = list(
      ate_true = mean(tau),
      ite_true = as.numeric(tau)
    ),
    X_original = X,
    propensity_score = as.numeric(ps)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
```

**Step 3: Commit**

```bash
git add benchmarks/tcga/
git commit -m "feat(benchmark): add TCGA data preparation pipeline"
```

---

### Task 11: Create TCGA benchmark runner

**Files:**
- Create: `benchmarks/tcga/run_tcga_benchmark.R`

**Step 1: Write TCGA benchmark runner**

```r
# benchmarks/tcga/run_tcga_benchmark.R
# Run benchmarks on TCGA semi-synthetic data

source("benchmarks/tcga/preprocess_tcga.R")

run_tcga_benchmark <- function(
    project = "TCGA-BRCA",
    methods = c("gformula", "hdml", "tmle"),
    n_reps = 20,
    dgp_params = list(),
    output_dir = "benchmarks/tcga/results"
) {
  # Load TCGA data
  se_path <- file.path("benchmarks/tcga/data", paste0(project, "_se.rds"))
  if (!file.exists(se_path)) {
    stop("TCGA data not found. Run download_tcga.R first.")
  }
  se <- readRDS(se_path)

  # Preprocess
  tcga_data <- preprocess_tcga(se, n_genes = 500)

  message(sprintf("TCGA %s: %d samples, %d genes",
                  project, nrow(tcga_data$X), ncol(tcga_data$X)))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- list()
  idx <- 1

  for (rep in seq_len(n_reps)) {
    set.seed(rep)

    # Create semi-synthetic data
    sim_data <- create_semi_synthetic_data(tcga_data, dgp_params)

    for (method in methods) {
      message(sprintf("  Rep %d, Method %s", rep, method))

      result <- tryCatch({
        # Build formula
        cov_names <- grep("^PC|^X", names(sim_data$data), value = TRUE)
        k <- min(10, length(cov_names))
        fml <- as.formula(paste("Y ~ T |", paste(cov_names[1:k], collapse = " + ")))

        fit_time <- system.time({
          fit <- cfomics::cf_fit(fml, data = sim_data$data, method = method)
        })

        ate_hat <- predict(fit, type = "ate")
        ite_hat <- predict(fit, type = "ite")
        summary_hat <- tryCatch(predict(fit, type = "summary"), error = function(e) NULL)

        metrics <- cfomics:::cf_benchmark_compute_metrics(
          ate_hat = ate_hat,
          ite_hat = ite_hat,
          summary_hat = summary_hat,
          truth = sim_data$truth
        )

        data.frame(
          project = project,
          method = method,
          rep = rep,
          n = nrow(sim_data$data),
          p = length(cov_names),
          ate_true = sim_data$truth$ate_true,
          ate_hat = ate_hat,
          bias_ate = metrics$bias_ate,
          mse_ate = metrics$mse_ate,
          pehe = metrics$pehe,
          coverage_ate = metrics$coverage_ate,
          time_sec = as.numeric(fit_time["elapsed"]),
          status = "ok",
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(
          project = project,
          method = method,
          rep = rep,
          n = NA, p = NA,
          ate_true = NA, ate_hat = NA,
          bias_ate = NA, mse_ate = NA, pehe = NA, coverage_ate = NA,
          time_sec = NA,
          status = "error",
          stringsAsFactors = FALSE
        )
      })

      results[[idx]] <- result
      idx <- idx + 1
    }
  }

  results_df <- do.call(rbind, results)

  # Save
  output_file <- file.path(output_dir, paste0(project, "_results.rds"))
  saveRDS(results_df, output_file)

  message(sprintf("Results saved to %s", output_file))
  results_df
}

# Main execution
if (sys.nframe() == 0) {
  library(cfomics)

  projects <- c("TCGA-BRCA")
  methods <- c("gformula", "hdml", "hdps", "tmle")

  for (proj in projects) {
    message(sprintf("\n=== Running %s benchmark ===", proj))
    run_tcga_benchmark(proj, methods = methods, n_reps = 20)
  }
}
```

**Step 2: Commit**

```bash
git add benchmarks/tcga/run_tcga_benchmark.R
git commit -m "feat(benchmark): add TCGA semi-synthetic benchmark runner"
```

---

## Phase 4: New Metrics

### Task 12: Add computational metrics to benchmark

**Files:**
- Modify: `packages/cfomics/R/benchmark_metrics.R`
- Modify: `benchmarks/R/runner.R`

**Step 1: Add memory tracking to benchmark_metrics.R**

Add to `packages/cfomics/R/benchmark_metrics.R`:

```r
#' Compute computational metrics
#'
#' @param fit_fn Function to call for fitting
#' @param ... Arguments to pass to fit_fn
#' @return List with time_sec, peak_memory_mb
#' @export
cf_benchmark_compute_computational_metrics <- function(fit_fn, ...) {
  # Memory before
  gc(reset = TRUE)
  mem_before <- sum(gc()[, 2])

  # Time
  time_result <- system.time({
    result <- fit_fn(...)
  })

  # Memory after
  gc()
  mem_after <- sum(gc()[, 2])
  peak_memory_mb <- mem_after - mem_before

  list(
    result = result,
    time_sec = as.numeric(time_result["elapsed"]),
    peak_memory_mb = peak_memory_mb
  )
}
```

**Step 2: Update runner to include memory metrics**

In `benchmarks/R/runner.R`, update `run_single_job` to track memory:

```r
# Add after line 30, replace fit_time block:
    gc(reset = TRUE)
    mem_before <- sum(gc()[, 2])

    fit_time <- system.time({
      fit <- do.call(cfomics::cf_fit, c(list(formula = fml, data = gen$data, method = method), extra_args))
    })

    gc()
    mem_after <- sum(gc()[, 2])
    peak_memory_mb <- max(0, mem_after - mem_before)
```

And add `peak_memory_mb` to the result data.frame.

**Step 3: Commit**

```bash
git add packages/cfomics/R/benchmark_metrics.R benchmarks/R/runner.R
git commit -m "feat(benchmark): add computational metrics (time, memory)"
```

---

### Task 13: Add missing data scenario (S15)

**Files:**
- Modify: `packages/cfomics/R/benchmark_dgp.R`
- Modify: `benchmarks/config.R`

**Step 1: Add missing data DGP**

Add to `packages/cfomics/R/benchmark_dgp.R`:

```r
#' Generate missing data DGP (S15)
#'
#' Generates data with missing values under MCAR or MAR mechanisms.
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param missing_rate Numeric, proportion of missing values (0-1)
#' @param missing_type Character, "MCAR" or "MAR"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, missing_mask
#' @export
dgp_missing_data <- function(n = 500, p = 500,
                              missing_rate = 0.2,
                              missing_type = "MCAR",
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate complete data first
  X_complete <- matrix(stats::rnorm(n * p), n, p)
  colnames(X_complete) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- stats::plogis(X_complete %*% beta_t)
  T <- stats::rbinom(n, 1, ps)

  tau <- 2.0
  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X_complete %*% beta_y + tau * T + stats::rnorm(n)

  # Create missing mask
  if (missing_type == "MCAR") {
    missing_mask <- matrix(
      stats::rbinom(n * p, 1, missing_rate),
      n, p
    )
  } else {
    # MAR: missingness depends on observed values
    missing_prob <- stats::plogis(-1 + 0.5 * X_complete[, 1])
    missing_mask <- matrix(0, n, p)
    for (j in 2:p) {
      missing_mask[, j] <- stats::rbinom(n, 1, missing_prob * missing_rate * 2)
    }
  }

  # Apply missing values
  X <- X_complete
  X[missing_mask == 1] <- NA

  list(
    X = X,
    X_complete = X_complete,
    T = T,
    Y = as.numeric(Y),
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = as.numeric(ps),
    missing_mask = missing_mask,
    missing_rate_actual = mean(missing_mask),
    dgp_name = "missing_data",
    dgp_params = list(n = n, p = p, missing_rate = missing_rate,
                      missing_type = missing_type)
  )
}
```

**Step 2: Add to config.R**

```r
      # S15: Missing data
      list(id = "S15_mcar_10", dgp = "dgp_missing_data", n = 1000L, p = 100L,
           params = list(missing_rate = 0.1, missing_type = "MCAR")),
      list(id = "S15_mcar_30", dgp = "dgp_missing_data", n = 1000L, p = 100L,
           params = list(missing_rate = 0.3, missing_type = "MCAR")),
      list(id = "S15_mar_10", dgp = "dgp_missing_data", n = 1000L, p = 100L,
           params = list(missing_rate = 0.1, missing_type = "MAR")),
      list(id = "S15_mar_30", dgp = "dgp_missing_data", n = 1000L, p = 100L,
           params = list(missing_rate = 0.3, missing_type = "MAR"))
```

**Step 3: Update generate_benchmark_data switch**

**Step 4: Commit**

```bash
git add packages/cfomics/R/benchmark_dgp.R benchmarks/config.R
git commit -m "feat(benchmark): add missing data scenario (S15)"
```

---

## Phase 5: Publication Figures

### Task 14: Create figure generation scripts

**Files:**
- Create: `benchmarks/figures/generate_figures.R`
- Create: `benchmarks/figures/fig1_workflow.R`
- Create: `benchmarks/figures/fig2_simulation.R`
- Create: `benchmarks/figures/fig3_highdim.R`
- Create: `benchmarks/figures/fig4_tcga.R`
- Create: `benchmarks/figures/fig5_guidelines.R`

**Step 1: Create main figure generation script**

```r
# benchmarks/figures/generate_figures.R
# Generate all publication figures

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Set theme
theme_paper <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      legend.position = "bottom"
    )
}

# Source individual figure scripts
source("benchmarks/figures/fig2_simulation.R")
source("benchmarks/figures/fig3_highdim.R")
source("benchmarks/figures/fig4_tcga.R")

# Generate all figures
generate_all_figures <- function(results_dir = "benchmarks/results",
                                  output_dir = "benchmarks/figures/output") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Load results
  results <- load_benchmark_results(results_dir)

  # Fig 2: Main simulation results
  fig2 <- generate_fig2_simulation(results)
  ggsave(file.path(output_dir, "fig2_simulation.pdf"), fig2,
         width = 180, height = 120, units = "mm")

  # Fig 3: High-dimensional stability
  fig3 <- generate_fig3_highdim(results)
  ggsave(file.path(output_dir, "fig3_highdim.pdf"), fig3,
         width = 180, height = 100, units = "mm")

  # Fig 4: TCGA validation
  tcga_results <- load_tcga_results("benchmarks/tcga/results")
  fig4 <- generate_fig4_tcga(tcga_results)
  ggsave(file.path(output_dir, "fig4_tcga.pdf"), fig4,
         width = 180, height = 100, units = "mm")

  message("All figures generated in ", output_dir)
}

load_benchmark_results <- function(results_dir) {
  raw_dir <- file.path(results_dir, "raw")
  files <- list.files(raw_dir, pattern = "\\.rds$", full.names = TRUE)

  results_list <- lapply(files, readRDS)
  do.call(rbind, results_list)
}

load_tcga_results <- function(results_dir) {
  files <- list.files(results_dir, pattern = "_results\\.rds$", full.names = TRUE)
  results_list <- lapply(files, readRDS)
  do.call(rbind, results_list)
}
```

**Step 2: Create Fig 2 simulation script**

```r
# benchmarks/figures/fig2_simulation.R
# Figure 2: Main simulation benchmark results

generate_fig2_simulation <- function(results) {
  # Filter to key scenarios
  key_scenarios <- c("S1_n1000_p100", "S5_combined", "S7_weak", "S12_nonlinear_outcome_moderate")

  df <- results %>%
    filter(scenario_id %in% key_scenarios, status == "ok") %>%
    group_by(scenario_id, method) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      mean_pehe = mean(pehe, na.rm = TRUE),
      coverage = mean(coverage_ate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      scenario_label = case_when(
        grepl("S1", scenario_id) ~ "Baseline",
        grepl("S5", scenario_id) ~ "Nonlinear confounding",
        grepl("S7", scenario_id) ~ "Weak overlap",
        grepl("S12", scenario_id) ~ "Nonlinear outcome"
      )
    )

  # Panel A: RMSE
  p_rmse <- ggplot(df, aes(x = method, y = rmse_ate, fill = method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~scenario_label, scales = "free_y") +
    labs(x = NULL, y = "RMSE (ATE)") +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Panel B: PEHE
  p_pehe <- ggplot(df, aes(x = method, y = mean_pehe, fill = method)) +
    geom_bar(stat = "identity") +
    facet_wrap(~scenario_label, scales = "free_y") +
    labs(x = NULL, y = "PEHE") +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p_rmse / p_pehe +
    plot_annotation(tag_levels = "A")
}
```

**Step 3: Create remaining figure scripts (abbreviated)**

```r
# benchmarks/figures/fig3_highdim.R
generate_fig3_highdim <- function(results) {
  # Filter S2 dimension sweep scenarios
  df <- results %>%
    filter(grepl("^S2", scenario_id), status == "ok") %>%
    mutate(p_n_ratio = p / n) %>%
    group_by(method, p_n_ratio) %>%
    summarise(rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)), .groups = "drop")

  ggplot(df, aes(x = p_n_ratio, y = rmse_ate, color = method)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    labs(x = "p/n ratio", y = "RMSE (ATE)", color = "Method") +
    theme_paper()
}
```

```r
# benchmarks/figures/fig4_tcga.R
generate_fig4_tcga <- function(results) {
  df <- results %>%
    filter(status == "ok") %>%
    group_by(project, method) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      mean_pehe = mean(pehe, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(df, aes(x = method, y = rmse_ate, fill = project)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Method", y = "RMSE (ATE)", fill = "Dataset") +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
```

**Step 4: Commit**

```bash
git add benchmarks/figures/
git commit -m "feat(benchmark): add publication figure generation scripts"
```

---

## Phase 6: Run Full Benchmark

### Task 15: Execute full benchmark suite

**Files:**
- Existing: `benchmarks/run_full_benchmark.R`

**Step 1: Run simulation benchmarks**

```bash
cd /Users/bmhi/Documents/GitHub/cfomics
Rscript benchmarks/run_full_benchmark.R
```

Expected: Results saved to `benchmarks/results/raw/`

**Step 2: Run external package benchmarks**

```bash
Rscript -e "
source('benchmarks/config.R')
source('benchmarks/R/runner.R')
cfg <- benchmark_config()
cfg$methods <- c()
cfg$external_methods <- c('matchit', 'weightit', 'hdm', 'doubleml', 'bart')
# Run external benchmark (separate script needed)
"
```

**Step 3: Run TCGA benchmarks**

```bash
Rscript benchmarks/tcga/run_tcga_benchmark.R
```

**Step 4: Generate figures**

```bash
Rscript benchmarks/figures/generate_figures.R
```

**Step 5: Commit results summary**

```bash
git add benchmarks/results/summary/ benchmarks/figures/output/
git commit -m "results: add full benchmark results and figures"
```

---

## Summary

| Phase | Tasks | Key Deliverables |
|-------|-------|------------------|
| 1 | Tasks 1-8 | External package wrappers (MatchIt, WeightIt, hdm, DoubleML, tmle3, BART) |
| 2 | Task 9 | TARNet/CFRNet/DragonNet comparison scripts |
| 3 | Tasks 10-11 | TCGA data pipeline and semi-synthetic benchmark |
| 4 | Tasks 12-13 | Computational metrics, missing data scenario |
| 5 | Task 14 | Publication figure scripts |
| 6 | Task 15 | Execute full benchmark suite |

Total: **15 Tasks**, each with TDD approach (test → implement → verify → commit)
