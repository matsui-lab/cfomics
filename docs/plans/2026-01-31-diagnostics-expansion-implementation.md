# Diagnostics Expansion Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add sensitivity analysis (E-value), model diagnostics (PS/outcome checks), CATE validation (GATES/BLP), and a unified diagnostics report to cfomics.

**Architecture:** 4 new R source files + 4 test files. Each diagnostic is a standalone function returning an S3 class with print/plot methods. The report function composes them all. All functions take `cfomics_result` objects as input.

**Tech Stack:** R, ggplot2 (optional), stats (base R), rlang (errors)

**Key reference files:**
- `packages/cfomics/R/result_contract.R` — Result object structure (`fit$res$ate`, `fit$res$summary$ate_ci_lower/upper`, `fit$res$ite`, `fit$res$y0_hat/y1_hat`)
- `packages/cfomics/R/diagnostics.R` — Existing `cf_balance_check()`, `cf_overlap_check()` patterns

---

### Task 1: cf_sensitivity — E-value computation

**Files:**
- Create: `packages/cfomics/R/sensitivity.R`
- Create: `packages/cfomics/tests/testthat/test-sensitivity.R`

**Step 1: Write failing test**

Create `packages/cfomics/tests/testthat/test-sensitivity.R`:

```r
test_that("cf_sensitivity computes E-value for a cfomics_result", {
  # Create minimal cfomics_result
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 1.5,
      ite = rep(1.5, 100),
      y0_hat = rnorm(100, 5, 1),
      y1_hat = rnorm(100, 6.5, 1),
      summary = list(ate = 1.5, ate_ci_lower = 1.0, ate_ci_upper = 2.0)
    )),
    meta = list(formula = Y ~ T | X, n = 100L, p = 5L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = paste0("X", 1:5))
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  result <- cf_sensitivity(mock_result)

  expect_s3_class(result, "cf_sensitivity")
  expect_true(is.numeric(result$evalue_point))
  expect_true(is.numeric(result$evalue_ci))
  expect_true(result$evalue_point > 1)
  expect_true(result$evalue_ci > 1)
  # E-value for point should be >= E-value for CI lower bound
  expect_gte(result$evalue_point, result$evalue_ci)
  expect_true(is.character(result$interpretation))
})

test_that("cf_sensitivity returns E-value = 1 when ATE is 0", {
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 0,
      ite = rep(0, 100),
      y0_hat = rnorm(100, 5, 1),
      y1_hat = rnorm(100, 5, 1),
      summary = list(ate = 0, ate_ci_lower = -0.5, ate_ci_upper = 0.5)
    )),
    meta = list(formula = Y ~ T | X, n = 100L, p = 5L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = paste0("X", 1:5))
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  result <- cf_sensitivity(mock_result)
  expect_equal(result$evalue_point, 1)
})

test_that("cf_sensitivity print method works", {
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 1.5,
      ite = rep(1.5, 100),
      y0_hat = rnorm(100, 5, 1),
      y1_hat = rnorm(100, 6.5, 1),
      summary = list(ate = 1.5, ate_ci_lower = 1.0, ate_ci_upper = 2.0)
    )),
    meta = list(formula = Y ~ T | X, n = 100L, p = 5L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = paste0("X", 1:5))
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  sens <- cf_sensitivity(mock_result)
  expect_output(print(sens), "E-value")
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-sensitivity.R")'`
Expected: FAIL — `cf_sensitivity` not found.

**Step 3: Implement cf_sensitivity**

Create `packages/cfomics/R/sensitivity.R`:

```r
#' Sensitivity Analysis for Unmeasured Confounding
#'
#' Computes the E-value (VanderWeele & Ding, 2017) to quantify robustness
#' of causal effect estimates to unmeasured confounding.
#'
#' @param result A cfomics_result object from cf_fit()
#' @param type Character, type of sensitivity analysis. Currently only "evalue".
#' @param alpha Numeric, significance level for CI-based E-value (default 0.05)
#' @return A list of class "cf_sensitivity"
#' @export
#' @examples
#' \dontrun{
#' fit <- cf_fit(Y ~ T | X1 + X2, data = df, method = "gformula")
#' sens <- cf_sensitivity(fit)
#' print(sens)
#' plot(sens)
#' }
cf_sensitivity <- function(result, type = "evalue", alpha = 0.05) {
  if (!inherits(result, "cfomics_result")) {
    rlang::abort("Input must be a cfomics_result object",
                 class = "cfomics_invalid_input")
  }

  type <- match.arg(type, choices = "evalue")

  ate <- result$fit$res$ate
  ci_lower <- result$fit$res$summary$ate_ci_lower
  ci_upper <- result$fit$res$summary$ate_ci_upper

  # Compute SD of observed outcomes from counterfactuals
  y0 <- result$fit$res$y0_hat
  y1 <- result$fit$res$y1_hat
  sd_y <- stats::sd(c(y0, y1))
  if (is.na(sd_y) || sd_y == 0) sd_y <- 1

  # Convert ATE to approximate risk ratio
  rr_point <- exp(abs(ate) / sd_y)
  rr_ci <- exp(abs(ci_lower) / sd_y)
  # Use the CI bound closer to null for conservative E-value
  if (ci_lower > 0) {
    rr_ci <- exp(ci_lower / sd_y)
  } else if (ci_upper < 0) {
    rr_ci <- exp(abs(ci_upper) / sd_y)
  } else {
    # CI crosses null — E-value for CI is 1
    rr_ci <- 1
  }

  # E-value formula
  evalue_fn <- function(rr) {
    if (rr <= 1) return(1)
    rr + sqrt(rr * (rr - 1))
  }

  ev_point <- evalue_fn(rr_point)
  ev_ci <- evalue_fn(rr_ci)

  # Interpretation
  interp <- if (ev_point > 2) {
    "Robust: an unmeasured confounder would need strong associations with both treatment and outcome to explain away this effect."
  } else if (ev_point > 1.5) {
    "Moderate robustness: an unmeasured confounder with moderate associations could potentially explain this effect."
  } else {
    "Weak robustness: even a weak unmeasured confounder could explain away this effect."
  }

  result_out <- list(
    evalue_point = ev_point,
    evalue_ci = ev_ci,
    ate = ate,
    ci = c(lower = ci_lower, upper = ci_upper),
    rr_point = rr_point,
    sd_y = sd_y,
    interpretation = interp,
    method = result$method
  )

  class(result_out) <- c("cf_sensitivity", "list")
  result_out
}

#' @export
print.cf_sensitivity <- function(x, ...) {
  cat("\nSensitivity Analysis (E-value)\n")
  cat("==============================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("ATE: %.3f [%.3f, %.3f]\n", x$ate, x$ci["lower"], x$ci["upper"]))
  cat(sprintf("E-value (point): %.2f\n", x$evalue_point))
  cat(sprintf("E-value (CI):    %.2f\n", x$evalue_ci))
  cat(sprintf("\n%s\n", x$interpretation))
  invisible(x)
}

#' @export
plot.cf_sensitivity <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for sensitivity plot")
    return(invisible(NULL))
  }

  # Bias plot: range of confounder-treatment and confounder-outcome RRs
  rr_seq <- seq(1, max(x$evalue_point * 1.5, 3), length.out = 100)
  # For each RR_UD (confounder-outcome), compute minimum RR_EU (confounder-treatment)
  # that would explain away the effect
  bias_data <- data.frame(
    rr_outcome = rr_seq,
    rr_treatment_point = sapply(rr_seq, function(rr_uy) {
      # Solve: bias_factor >= rr_point requires RR_EU * RR_UY / (RR_EU + RR_UY - 1) >= rr_point
      # Simplified: for each RR_UY, find threshold
      if (rr_uy <= 1) return(Inf)
      x$rr_point * (rr_uy - 1 + x$rr_point) / (rr_uy * x$rr_point - rr_uy + 1)
    })
  )
  bias_data$rr_treatment_point[bias_data$rr_treatment_point < 1] <- NA
  bias_data$rr_treatment_point[bias_data$rr_treatment_point > max(rr_seq) * 2] <- NA

  p <- ggplot2::ggplot(bias_data, ggplot2::aes(x = rr_outcome, y = rr_treatment_point)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_point(data = data.frame(x = x$evalue_point, y = x$evalue_point),
                        ggplot2::aes(x = x, y = y), color = "red", size = 3) +
    ggplot2::annotate("text", x = x$evalue_point, y = x$evalue_point,
                      label = sprintf("E = %.1f", x$evalue_point),
                      vjust = -1, color = "red") +
    ggplot2::labs(x = "Confounder-Outcome RR",
                  y = "Confounder-Treatment RR",
                  title = "Sensitivity Analysis: E-value Bias Plot",
                  subtitle = sprintf("ATE = %.3f, E-value = %.2f", x$ate, x$evalue_point)) +
    ggplot2::theme_minimal()
  print(p)
  invisible(p)
}
```

**Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-sensitivity.R")'`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add packages/cfomics/R/sensitivity.R packages/cfomics/tests/testthat/test-sensitivity.R
git commit -m "feat: add cf_sensitivity for E-value sensitivity analysis"
```

---

### Task 2: cf_model_diagnostics — PS and outcome model checks

**Files:**
- Create: `packages/cfomics/R/model_diagnostics.R`
- Create: `packages/cfomics/tests/testthat/test-model-diagnostics.R`

**Step 1: Write failing test**

Create `packages/cfomics/tests/testthat/test-model-diagnostics.R`:

```r
test_that("cf_model_diagnostics returns valid diagnostics", {
  set.seed(42)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  ps <- plogis(0.5 * X$X1 - 0.3 * X$X2)
  T_var <- rbinom(n, 1, ps)
  Y <- 1.0 * X$X1 + 2.0 * T_var + rnorm(n)
  df <- data.frame(Y = Y, T = T_var, X)

  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 2.0,
      ite = rep(2.0, n),
      y0_hat = Y - 2.0 * T_var,
      y1_hat = Y - 2.0 * T_var + 2.0,
      summary = list(ate = 2.0, ate_ci_lower = 1.5, ate_ci_upper = 2.5)
    )),
    meta = list(formula = Y ~ T | X1 + X2, n = as.integer(n), p = 2L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = c("X1", "X2"))
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  diag <- cf_model_diagnostics(mock_result, df, Y ~ T | X1 + X2)

  expect_s3_class(diag, "cf_model_diagnostics")
  expect_true(is.numeric(diag$ps_cstat))
  expect_true(diag$ps_cstat > 0.5 && diag$ps_cstat < 1.0)
  expect_true(is.numeric(diag$extreme_weight_pct))
  expect_true(is.data.frame(diag$residual_summary))
})

test_that("cf_model_diagnostics handles methods without PS", {
  n <- 100
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 1.0, ite = rep(1.0, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 1),
      summary = list(ate = 1.0, ate_ci_lower = 0.5, ate_ci_upper = 1.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  diag <- cf_model_diagnostics(mock_result, df, Y ~ T | X1)

  # Should still work — PS is estimated internally for diagnostics
  expect_s3_class(diag, "cf_model_diagnostics")
})

test_that("cf_model_diagnostics print works", {
  n <- 100
  mock_result <- list(
    method = "hdml",
    fit = list(res = list(
      ate = 1.0, ite = rep(1.0, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 1),
      summary = list(ate = 1.0, ate_ci_lower = 0.5, ate_ci_upper = 1.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  diag <- cf_model_diagnostics(mock_result, df, Y ~ T | X1)
  expect_output(print(diag), "Model Diagnostics")
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-model-diagnostics.R")'`
Expected: FAIL

**Step 3: Implement cf_model_diagnostics**

Create `packages/cfomics/R/model_diagnostics.R`:

```r
#' Model Diagnostics for Causal Inference Results
#'
#' Evaluates the quality of propensity score and outcome models used in
#' causal inference estimation.
#'
#' @param result A cfomics_result object
#' @param data Data frame used for fitting
#' @param formula Formula of the form Y ~ T | X1 + X2 + ...
#' @return A list of class "cf_model_diagnostics"
#' @export
cf_model_diagnostics <- function(result, data, formula) {
  if (!inherits(result, "cfomics_result")) {
    rlang::abort("Input must be a cfomics_result object",
                 class = "cfomics_invalid_input")
  }

  # Parse formula
  parsed <- parse_cf_formula(formula)
  outcome_name <- parsed$outcome
  treatment_name <- parsed$treatment
  covariate_names <- parsed$covariates

  Y <- data[[outcome_name]]
  T_var <- data[[treatment_name]]
  X <- as.matrix(data[, covariate_names, drop = FALSE])
  n <- length(Y)

  # --- PS diagnostics ---
  # Fit logistic PS model for diagnostics
  ps_formula <- stats::as.formula(
    paste(treatment_name, "~", paste(covariate_names, collapse = " + "))
  )
  ps_model <- stats::glm(ps_formula, data = data, family = stats::binomial())
  ps_hat <- stats::predict(ps_model, type = "response")

  # C-statistic (concordance / AUC via Mann-Whitney)
  ps_t1 <- ps_hat[T_var == 1]
  ps_t0 <- ps_hat[T_var == 0]
  cstat <- mean(outer(ps_t1, ps_t0, ">")) + 0.5 * mean(outer(ps_t1, ps_t0, "=="))

  # Extreme weights
  weights <- ifelse(T_var == 1, 1 / ps_hat, 1 / (1 - ps_hat))
  extreme_threshold <- stats::quantile(weights, 0.99)
  extreme_pct <- mean(weights > extreme_threshold) * 100

  # --- Outcome model diagnostics ---
  y0_hat <- result$fit$res$y0_hat
  y1_hat <- result$fit$res$y1_hat
  y_hat <- ifelse(T_var == 1, y1_hat, y0_hat)
  residuals <- Y - y_hat

  # Residual summary by treatment group
  resid_summary <- data.frame(
    group = c("Control (T=0)", "Treated (T=1)", "Overall"),
    mean_resid = c(
      mean(residuals[T_var == 0]),
      mean(residuals[T_var == 1]),
      mean(residuals)
    ),
    sd_resid = c(
      stats::sd(residuals[T_var == 0]),
      stats::sd(residuals[T_var == 1]),
      stats::sd(residuals)
    ),
    stringsAsFactors = FALSE
  )

  # Normality test (on residuals, capped at 5000)
  resid_sample <- if (n > 5000) sample(residuals, 5000) else residuals
  normality_p <- tryCatch(
    stats::shapiro.test(resid_sample)$p.value,
    error = function(e) NA_real_
  )

  out <- list(
    ps_cstat = cstat,
    ps_hat = ps_hat,
    extreme_weight_pct = extreme_pct,
    weights = weights,
    residuals = residuals,
    residual_summary = resid_summary,
    normality_pvalue = normality_p,
    y_hat = y_hat,
    treatment = T_var,
    method = result$method
  )

  class(out) <- c("cf_model_diagnostics", "list")
  out
}

#' @export
print.cf_model_diagnostics <- function(x, ...) {
  cat("\nModel Diagnostics\n")
  cat("==================\n")
  cat(sprintf("Method: %s\n\n", x$method))

  cat("PS Model:\n")
  cat(sprintf("  C-statistic (AUC): %.3f", x$ps_cstat))
  if (x$ps_cstat > 0.95) cat(" [WARNING: near-deterministic]")
  else if (x$ps_cstat < 0.55) cat(" [WARNING: poor discrimination]")
  else cat(" [OK]")
  cat("\n")
  cat(sprintf("  Extreme weights (>99th pct): %.1f%%\n\n", x$extreme_weight_pct))

  cat("Outcome Model Residuals:\n")
  print(x$residual_summary, row.names = FALSE)

  if (!is.na(x$normality_pvalue)) {
    cat(sprintf("\n  Shapiro-Wilk p-value: %.4f", x$normality_pvalue))
    if (x$normality_pvalue < 0.05) cat(" [non-normal]")
    cat("\n")
  }

  invisible(x)
}

#' @export
plot.cf_model_diagnostics <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for diagnostic plots")
    return(invisible(NULL))
  }

  df_resid <- data.frame(
    fitted = x$y_hat,
    residuals = x$residuals,
    group = ifelse(x$treatment == 1, "Treated", "Control")
  )
  df_ps <- data.frame(
    ps = x$ps_hat,
    group = ifelse(x$treatment == 1, "Treated", "Control")
  )

  p1 <- ggplot2::ggplot(df_ps, ggplot2::aes(x = ps, fill = group)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(title = "PS Distribution", x = "Propensity Score", fill = "Group") +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(df_resid, ggplot2::aes(x = fitted, y = residuals, color = group)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") +
    ggplot2::theme_minimal()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(p1 + p2)
  } else {
    print(p1)
    print(p2)
  }
  invisible(NULL)
}
```

**Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-model-diagnostics.R")'`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add packages/cfomics/R/model_diagnostics.R packages/cfomics/tests/testthat/test-model-diagnostics.R
git commit -m "feat: add cf_model_diagnostics for PS and outcome model checks"
```

---

### Task 3: cf_cate_diagnostics — GATES and BLP

**Files:**
- Create: `packages/cfomics/R/cate_diagnostics.R`
- Create: `packages/cfomics/tests/testthat/test-cate-diagnostics.R`

**Step 1: Write failing test**

Create `packages/cfomics/tests/testthat/test-cate-diagnostics.R`:

```r
test_that("cf_cate_diagnostics works with heterogeneous ITE", {
  set.seed(42)
  n <- 500
  X1 <- rnorm(n)
  T_var <- rbinom(n, 1, 0.5)
  ite <- 1.0 + 2.0 * X1  # heterogeneous
  Y <- X1 + ite * T_var + rnorm(n)
  df <- data.frame(Y = Y, T = T_var, X1 = X1)

  mock_result <- list(
    method = "grf",
    fit = list(res = list(
      ate = mean(ite), ite = ite,
      y0_hat = X1, y1_hat = X1 + ite,
      summary = list(ate = mean(ite), ate_ci_lower = 0.5, ate_ci_upper = 1.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  cate <- cf_cate_diagnostics(mock_result, df, Y ~ T | X1)

  expect_s3_class(cate, "cf_cate_diagnostics")
  expect_true(!is.null(cate$gates))
  expect_equal(nrow(cate$gates), 4)  # quartiles
  expect_true(!is.null(cate$blp))
  expect_true(cate$heterogeneity_detected)
})

test_that("cf_cate_diagnostics skips GATES/BLP for constant ITE", {
  n <- 200
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 2.0, ite = rep(2.0, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 2),
      summary = list(ate = 2.0, ate_ci_lower = 1.5, ate_ci_upper = 2.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  cate <- cf_cate_diagnostics(mock_result, df, Y ~ T | X1)

  expect_s3_class(cate, "cf_cate_diagnostics")
  expect_true(is.null(cate$gates))
  expect_true(is.null(cate$blp))
  expect_false(cate$heterogeneity_detected)
  expect_true(cate$constant_ite)
})

test_that("cf_cate_diagnostics print works", {
  n <- 200
  mock_result <- list(
    method = "grf",
    fit = list(res = list(
      ate = 2.0, ite = rnorm(n, 2, 1),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 2),
      summary = list(ate = 2.0, ate_ci_lower = 1.5, ate_ci_upper = 2.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  cate <- cf_cate_diagnostics(mock_result, df, Y ~ T | X1)
  expect_output(print(cate), "CATE")
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-cate-diagnostics.R")'`
Expected: FAIL

**Step 3: Implement cf_cate_diagnostics**

Create `packages/cfomics/R/cate_diagnostics.R`:

```r
#' CATE Diagnostics: Heterogeneity Validation
#'
#' Tests for treatment effect heterogeneity using GATES (Sorted Group ATE)
#' and BLP (Best Linear Predictor) approaches from Chernozhukov et al. (2018).
#'
#' @param result A cfomics_result object
#' @param data Data frame used for fitting
#' @param formula Formula of the form Y ~ T | X1 + X2 + ...
#' @param n_groups Integer, number of ITE quantile groups for GATES (default 4)
#' @return A list of class "cf_cate_diagnostics"
#' @export
cf_cate_diagnostics <- function(result, data, formula, n_groups = 4L) {
  if (!inherits(result, "cfomics_result")) {
    rlang::abort("Input must be a cfomics_result object",
                 class = "cfomics_invalid_input")
  }

  parsed <- parse_cf_formula(formula)
  Y <- data[[parsed$outcome]]
  T_var <- data[[parsed$treatment]]
  ite <- result$fit$res$ite
  n <- length(ite)

  # Check for constant ITE
  ite_var <- stats::var(ite)
  is_constant <- is.na(ite_var) || ite_var < 1e-10

  if (is_constant) {
    cli::cli_inform("ITE values are constant. GATES/BLP analysis skipped.")
    out <- list(
      constant_ite = TRUE,
      heterogeneity_detected = FALSE,
      ite_variance = 0,
      gates = NULL,
      blp = NULL,
      method = result$method
    )
    class(out) <- c("cf_cate_diagnostics", "list")
    return(out)
  }

  # --- GATES: Sorted Group ATE ---
  ite_quantiles <- stats::quantile(ite, probs = seq(0, 1, length.out = n_groups + 1))
  group <- cut(ite, breaks = ite_quantiles, labels = FALSE, include.lowest = TRUE)

  gates_results <- lapply(seq_len(n_groups), function(g) {
    idx <- group == g
    if (sum(idx & T_var == 1) < 2 || sum(idx & T_var == 0) < 2) {
      return(data.frame(group = g, gate = NA_real_, se = NA_real_,
                        ci_lower = NA_real_, ci_upper = NA_real_,
                        n = sum(idx)))
    }
    fit_g <- stats::lm(Y[idx] ~ T_var[idx])
    coefs <- summary(fit_g)$coefficients
    gate <- coefs[2, 1]
    se <- coefs[2, 2]
    data.frame(group = g, gate = gate, se = se,
               ci_lower = gate - 1.96 * se, ci_upper = gate + 1.96 * se,
               n = sum(idx))
  })
  gates_df <- do.call(rbind, gates_results)

  # --- BLP: Best Linear Predictor ---
  ite_centered <- ite - mean(ite)
  blp_fit <- stats::lm(Y ~ T_var + T_var:ite_centered)
  blp_summary <- summary(blp_fit)$coefficients

  blp_result <- list(
    ate_coef = blp_summary["T_var", "Estimate"],
    ate_se = blp_summary["T_var", "Std. Error"],
    ate_pvalue = blp_summary["T_var", "Pr(>|t|)"],
    het_coef = blp_summary["T_var:ite_centered", "Estimate"],
    het_se = blp_summary["T_var:ite_centered", "Std. Error"],
    het_pvalue = blp_summary["T_var:ite_centered", "Pr(>|t|)"]
  )

  # Heterogeneity detection
  het_detected <- blp_result$het_pvalue < 0.05

  out <- list(
    constant_ite = FALSE,
    heterogeneity_detected = het_detected,
    ite_variance = ite_var,
    gates = gates_df,
    blp = blp_result,
    method = result$method
  )

  class(out) <- c("cf_cate_diagnostics", "list")
  out
}

#' @export
print.cf_cate_diagnostics <- function(x, ...) {
  cat("\nCATE Diagnostics\n")
  cat("=================\n")
  cat(sprintf("Method: %s\n", x$method))

  if (x$constant_ite) {
    cat("\nITE is constant (no heterogeneity). GATES/BLP skipped.\n")
    return(invisible(x))
  }

  cat(sprintf("ITE variance: %.4f\n\n", x$ite_variance))

  cat("GATES (Sorted Group ATE):\n")
  if (!is.null(x$gates)) {
    print(x$gates, row.names = FALSE, digits = 3)
  }

  if (!is.null(x$blp)) {
    cat("\nBLP (Best Linear Predictor):\n")
    cat(sprintf("  ATE coefficient: %.3f (p = %.4f)\n",
                x$blp$ate_coef, x$blp$ate_pvalue))
    cat(sprintf("  Heterogeneity:   %.3f (p = %.4f)\n",
                x$blp$het_coef, x$blp$het_pvalue))
    if (x$heterogeneity_detected) {
      cat("  -> Significant heterogeneity detected\n")
    } else {
      cat("  -> No significant heterogeneity\n")
    }
  }

  invisible(x)
}

#' @export
plot.cf_cate_diagnostics <- function(x, ...) {
  if (x$constant_ite || is.null(x$gates)) {
    message("No GATES data to plot (constant ITE)")
    return(invisible(NULL))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for CATE diagnostic plots")
    return(invisible(NULL))
  }

  gates <- x$gates
  gates$group <- factor(gates$group)

  p <- ggplot2::ggplot(gates, ggplot2::aes(x = group, y = gate)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                           width = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(x = "ITE Quartile Group", y = "Group ATE",
                  title = "GATES: Sorted Group Average Treatment Effects",
                  subtitle = sprintf("Heterogeneity %s (BLP p = %.4f)",
                                     ifelse(x$heterogeneity_detected, "detected", "not detected"),
                                     x$blp$het_pvalue)) +
    ggplot2::theme_minimal()
  print(p)
  invisible(p)
}
```

**Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-cate-diagnostics.R")'`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add packages/cfomics/R/cate_diagnostics.R packages/cfomics/tests/testthat/test-cate-diagnostics.R
git commit -m "feat: add cf_cate_diagnostics with GATES and BLP"
```

---

### Task 4: cf_diagnostics_report — Unified report

**Files:**
- Create: `packages/cfomics/R/diagnostics_report.R`
- Create: `packages/cfomics/tests/testthat/test-diagnostics-report.R`

**Step 1: Write failing test**

Create `packages/cfomics/tests/testthat/test-diagnostics-report.R`:

```r
test_that("cf_diagnostics_report runs all diagnostics", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n)
  T_var <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- X1 + 2.0 * T_var + rnorm(n)
  df <- data.frame(Y = Y, T = T_var, X1 = X1)

  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 2.0, ite = rep(2.0, n),
      y0_hat = X1, y1_hat = X1 + 2.0,
      summary = list(ate = 2.0, ate_ci_lower = 1.5, ate_ci_upper = 2.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  report <- cf_diagnostics_report(mock_result, df, Y ~ T | X1)

  expect_s3_class(report, "cf_diagnostics_report")
  expect_s3_class(report$sensitivity, "cf_sensitivity")
  expect_s3_class(report$model, "cf_model_diagnostics")
  expect_s3_class(report$cate, "cf_cate_diagnostics")
  expect_s3_class(report$balance, "cf_balance")
  expect_true(is.data.frame(report$summary))
})

test_that("cf_diagnostics_report respects include parameter", {
  n <- 100
  mock_result <- list(
    method = "gformula",
    fit = list(res = list(
      ate = 1.0, ite = rep(1.0, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 1),
      summary = list(ate = 1.0, ate_ci_lower = 0.5, ate_ci_upper = 1.5)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  report <- cf_diagnostics_report(mock_result, df, Y ~ T | X1,
                                   include = "sensitivity")

  expect_s3_class(report$sensitivity, "cf_sensitivity")
  expect_null(report$model)
  expect_null(report$cate)
})

test_that("cf_diagnostics_report print works", {
  n <- 100
  mock_result <- list(
    method = "hdml",
    fit = list(res = list(
      ate = 1.5, ite = rep(1.5, n),
      y0_hat = rnorm(n), y1_hat = rnorm(n, 1.5),
      summary = list(ate = 1.5, ate_ci_lower = 1.0, ate_ci_upper = 2.0)
    )),
    meta = list(formula = Y ~ T | X1, n = as.integer(n), p = 1L,
                outcome_name = "Y", treatment_name = "T",
                covariate_names = "X1")
  )
  class(mock_result) <- c("cf_model", "cfomics_result")

  df <- data.frame(Y = rnorm(n), T = rbinom(n, 1, 0.5), X1 = rnorm(n))
  report <- cf_diagnostics_report(mock_result, df, Y ~ T | X1)
  expect_output(print(report), "Diagnostics Report")
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-diagnostics-report.R")'`
Expected: FAIL

**Step 3: Implement cf_diagnostics_report**

Create `packages/cfomics/R/diagnostics_report.R`:

```r
#' Unified Diagnostics Report
#'
#' Runs all available diagnostic checks and produces a comprehensive report.
#'
#' @param result A cfomics_result object
#' @param data Data frame used for fitting
#' @param formula Formula of the form Y ~ T | X1 + X2 + ...
#' @param include Character vector of diagnostics to include.
#'   Options: "sensitivity", "model", "cate". Default: all three.
#' @return A list of class "cf_diagnostics_report"
#' @export
cf_diagnostics_report <- function(result, data, formula,
                                   include = c("sensitivity", "model", "cate")) {
  if (!inherits(result, "cfomics_result")) {
    rlang::abort("Input must be a cfomics_result object",
                 class = "cfomics_invalid_input")
  }

  include <- match.arg(include, c("sensitivity", "model", "cate"),
                       several.ok = TRUE)

  parsed <- parse_cf_formula(formula)

  out <- list(
    sensitivity = NULL,
    model = NULL,
    cate = NULL,
    balance = NULL,
    summary = NULL,
    method = result$method,
    n = result$meta$n,
    p = result$meta$p
  )

  # Sensitivity
  if ("sensitivity" %in% include) {
    out$sensitivity <- cf_sensitivity(result)
  }

  # Model diagnostics
  if ("model" %in% include) {
    out$model <- cf_model_diagnostics(result, data, formula)
  }

  # CATE diagnostics
  if ("cate" %in% include) {
    out$cate <- cf_cate_diagnostics(result, data, formula)
  }

  # Balance check (always run)
  out$balance <- cf_balance_check(data, parsed$treatment, parsed$covariates)

  # Build summary
  summary_rows <- list()

  if (!is.null(out$sensitivity)) {
    ev <- out$sensitivity$evalue_point
    status <- if (ev > 2) "Pass" else if (ev > 1.5) "Warning" else "Fail"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "Sensitivity (E-value)",
      value = sprintf("%.2f", ev),
      status = status, stringsAsFactors = FALSE
    )))
  }

  if (!is.null(out$model)) {
    cstat <- out$model$ps_cstat
    status <- if (cstat >= 0.6 && cstat <= 0.85) "Pass"
              else if (cstat < 0.5 || cstat > 0.95) "Fail"
              else "Warning"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "PS C-statistic",
      value = sprintf("%.3f", cstat),
      status = status, stringsAsFactors = FALSE
    )))

    ewp <- out$model$extreme_weight_pct
    status <- if (ewp < 1) "Pass" else if (ewp < 5) "Warning" else "Fail"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "Extreme weights",
      value = sprintf("%.1f%%", ewp),
      status = status, stringsAsFactors = FALSE
    )))
  }

  if (!is.null(out$balance)) {
    max_smd <- max(abs(out$balance$smd))
    status <- if (max_smd < 0.1) "Pass" else if (max_smd < 0.25) "Warning" else "Fail"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "Balance (max SMD)",
      value = sprintf("%.3f", max_smd),
      status = status, stringsAsFactors = FALSE
    )))
  }

  if (!is.null(out$cate) && !out$cate$constant_ite) {
    het <- if (out$cate$heterogeneity_detected) "Detected" else "Not detected"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "CATE heterogeneity",
      value = het,
      status = "Info", stringsAsFactors = FALSE
    )))
  }

  out$summary <- do.call(rbind, summary_rows)

  class(out) <- c("cf_diagnostics_report", "list")
  out
}

#' @export
print.cf_diagnostics_report <- function(x, ...) {
  cat("\n")
  cat(cli::rule(left = "cfomics Diagnostics Report"))
  cat("\n")
  cat(sprintf("Method: %s | n = %d | p = %d\n\n", x$method, x$n, x$p))

  if (!is.null(x$summary)) {
    for (i in seq_len(nrow(x$summary))) {
      row <- x$summary[i, ]
      icon <- switch(row$status,
        Pass = cli::col_green(cli::symbol$tick),
        Warning = cli::col_yellow("!"),
        Fail = cli::col_red(cli::symbol$cross),
        Info = cli::col_blue("-")
      )
      cat(sprintf("%s %s: %s\n", icon, row$diagnostic, row$value))
    }
  }

  invisible(x)
}

#' @export
plot.cf_diagnostics_report <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for report plots")
    return(invisible(NULL))
  }

  plots <- list()

  if (!is.null(x$sensitivity)) {
    plots$sensitivity <- plot(x$sensitivity)
  }
  if (!is.null(x$model)) {
    plot(x$model)
  }
  if (!is.null(x$cate) && !x$cate$constant_ite) {
    plots$cate <- plot(x$cate)
  }
  if (!is.null(x$balance)) {
    plots$balance <- plot(x$balance)
  }

  invisible(NULL)
}
```

**Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::load_all("packages/cfomics"); testthat::test_file("packages/cfomics/tests/testthat/test-diagnostics-report.R")'`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add packages/cfomics/R/diagnostics_report.R packages/cfomics/tests/testthat/test-diagnostics-report.R
git commit -m "feat: add cf_diagnostics_report unified diagnostic report"
```

---

### Task 5: Update NAMESPACE and run full tests

**Files:**
- Modify: `packages/cfomics/NAMESPACE` (via roxygen2)

**Step 1: Run roxygen2**

```bash
Rscript -e 'devtools::document("packages/cfomics")'
```

**Step 2: Verify exports**

Check that NAMESPACE contains: `export(cf_sensitivity)`, `export(cf_model_diagnostics)`, `export(cf_cate_diagnostics)`, `export(cf_diagnostics_report)` and all print/plot S3 methods.

**Step 3: Run full test suite**

```bash
Rscript -e 'devtools::test("packages/cfomics")'
```
Expected: ALL PASS

**Step 4: Commit**

```bash
git add packages/cfomics/NAMESPACE packages/cfomics/man/
git commit -m "docs: update NAMESPACE and man pages for diagnostics expansion"
```
