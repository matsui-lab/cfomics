# High-Dimensional Causal Inference Methods Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add 4 high-dimensional causal inference methods (hdml, bcf, tmle, hdps) and comprehensive benchmark DGPs to cfomics.

**Architecture:** Each method follows the existing pattern in `methods_traditional.R`: a `cf_fit_*()` function returning the standard `cfomics_result` structure with `fit$res` containing `ite`, `ate`, `y0_hat`, `y1_hat`, and `summary`. Methods are registered in `cf_fit()` switch statement. DGPs are added to `benchmark_dgp.R` following the existing `.dgp_*()` pattern.

**Tech Stack:** R, testthat, glmnet, hdm (or DoubleML), bcf, dbarts, tmle3, sl3, MASS

---

## Phase 1: HDML (High-Dimensional Machine Learning / Debiased Lasso)

### Task 1.1: Create methods_hdml.R skeleton

**Files:**
- Create: `packages/cfomics/R/methods_hdml.R`
- Test: `packages/cfomics/tests/testthat/test-hdml.R`

**Step 1: Write the failing test**

Create `packages/cfomics/tests/testthat/test-hdml.R`:

```r
test_that("cf_fit_hdml basic execution", {

skip_if_not_installed("glmnet")

set.seed(123)
n <- 200
p <- 50
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)
beta_t <- c(rep(0.3, 5), rep(0, p - 5))
ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)
Y <- 2.0 * T + 0.5 * X[,1] + 0.3 * X[,2] + rnorm(n)
data <- data.frame(Y = Y, T = T, X)

fit <- cf_fit(Y ~ T | ., data = data, method = "hdml")

expect_s3_class(fit, "cf_model")
expect_s3_class(fit, "cfomics_result")
expect_equal(fit$method, "hdml")

ate <- predict(fit, type = "ate")
expect_type(ate, "double")
expect_true(ate > 1.0 && ate < 3.0)

summary_stats <- predict(fit, type = "summary")
expect_true("ate" %in% names(summary_stats))
expect_true("ate_ci_lower" %in% names(summary_stats))
expect_true("ate_ci_upper" %in% names(summary_stats))
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-hdml.R")'`
Expected: FAIL with "Unknown method: hdml"

**Step 3: Write minimal implementation**

Create `packages/cfomics/R/methods_hdml.R`:

```r
#' Fit High-Dimensional Machine Learning (Debiased Lasso) model
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param penalty Character, penalty type: "lasso", "ridge", "elastic_net"
#' @param lambda Character or numeric, lambda selection: "cv", "1se", or numeric
#' @param alpha Numeric, elastic net mixing parameter (1=lasso, 0=ridge)
#' @param nfolds Integer, number of CV folds
#' @param ... Additional arguments
#' @return List with model fit and results
#' @keywords internal
cf_fit_hdml <- function(X, T, Y,
                       covariate_names = NULL,
                       penalty = "lasso",
                       lambda = "cv",
                       alpha = NULL,
                       nfolds = 10L,
                       ...) {
if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Package 'glmnet' is required for method='hdml'. Please install it.")
}

X_mat <- as.matrix(X)
W <- as.numeric(T)
Y_vec <- as.numeric(Y)
n <- length(Y_vec)

if (is.null(covariate_names)) {
  covariate_names <- paste0("X", seq_len(ncol(X_mat)))
}

# Set alpha based on penalty type
if (is.null(alpha)) {
  alpha <- switch(penalty,
    "lasso" = 1,
    "ridge" = 0,
    "elastic_net" = 0.5,
    1
  )
}

# Step 1: Estimate propensity score with regularized logistic regression
ps_cv <- glmnet::cv.glmnet(X_mat, W, family = "binomial",
                           alpha = alpha, nfolds = nfolds)
ps_lambda <- if (lambda == "1se") ps_cv$lambda.1se else ps_cv$lambda.min
ps_hat <- as.numeric(predict(ps_cv, newx = X_mat, s = ps_lambda, type = "response"))
ps_hat <- pmax(pmin(ps_hat, 0.99), 0.01)  # Trim extreme values

# Step 2: Estimate outcome model with regularized regression
outcome_cv <- glmnet::cv.glmnet(cbind(X_mat, W), Y_vec,
                                 alpha = alpha, nfolds = nfolds)
outcome_lambda <- if (lambda == "1se") outcome_cv$lambda.1se else outcome_cv$lambda.min

# Step 3: Predict counterfactual outcomes
X_with_T0 <- cbind(X_mat, 0)
X_with_T1 <- cbind(X_mat, 1)
y0_hat <- as.numeric(predict(outcome_cv, newx = X_with_T0, s = outcome_lambda))
y1_hat <- as.numeric(predict(outcome_cv, newx = X_with_T1, s = outcome_lambda))

# Step 4: Compute AIPW/doubly-robust ATE estimator
ite_dr <- y1_hat - y0_hat +
  W * (Y_vec - y1_hat) / ps_hat -
  (1 - W) * (Y_vec - y0_hat) / (1 - ps_hat)

ate <- mean(ite_dr)

# Step 5: Estimate variance using influence function
influence_fn <- ite_dr - ate
ate_se <- sqrt(var(influence_fn) / n)

ate_ci_lower <- ate - 1.96 * ate_se
ate_ci_upper <- ate + 1.96 * ate_se

# ITE is approximated by y1_hat - y0_hat (model-based)
ite <- y1_hat - y0_hat

summary_stats <- list(
  ate = ate,
  ate_ci_lower = ate_ci_lower,
  ate_ci_upper = ate_ci_upper,
  ite_mean = mean(ite),
  ite_std = sd(ite),
  ite_quantiles = list(
    q05 = quantile(ite, 0.05),
    q50 = quantile(ite, 0.50),
    q95 = quantile(ite, 0.95)
  )
)

list(
  model = list(
    ps_model = ps_cv,
    outcome_model = outcome_cv,
    ps_lambda = ps_lambda,
    outcome_lambda = outcome_lambda
  ),
  res = list(
    ite = ite,
    ate = ate,
    y0_hat = y0_hat,
    y1_hat = y1_hat,
    summary = summary_stats,
    samples = list(
      y0 = y0_hat,
      y1 = y1_hat,
      ite = ite
    ),
    propensity = ps_hat,
    metadata = list(
      penalty = penalty,
      alpha = alpha,
      lambda = lambda
    )
  )
)
}

#' Predict from HDML model
#'
#' @param object cf_model object with method="hdml"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_hdml <- function(object, newdata = NULL, type = "ite", ...) {
if (is.null(newdata)) {
  res <- object$fit$res
  return(switch(
    type,
    "ate" = res$ate,
    "ite" = res$ite,
    "y0" = res$y0_hat,
    "y1" = res$y1_hat,
    "summary" = res$summary,
    "samples" = res$samples
  ))
}

# Prediction on new data
outcome_model <- object$fit$model$outcome_model
outcome_lambda <- object$fit$model$outcome_lambda

X_new <- as.matrix(newdata)
X_with_T0 <- cbind(X_new, 0)
X_with_T1 <- cbind(X_new, 1)

y0_hat <- as.numeric(predict(outcome_model, newx = X_with_T0, s = outcome_lambda))
y1_hat <- as.numeric(predict(outcome_model, newx = X_with_T1, s = outcome_lambda))
ite <- y1_hat - y0_hat

switch(
  type,
  "ate" = mean(ite),
  "ite" = ite,
  "y0" = y0_hat,
  "y1" = y1_hat,
  stop(sprintf("type '%s' not supported for newdata prediction", type))
)
}
```

**Step 4: Register method in cf_fit.R**

Modify `packages/cfomics/R/cf_fit.R:182-221`, add hdml case in the switch statement:

```r
"hdml" = cf_fit_hdml(X = X, T = T, Y = Y,
                     covariate_names = parsed$covariate_names,
                     ...),
```

Also update the method argument documentation and match.arg call to include "hdml".

**Step 5: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-hdml.R")'`
Expected: PASS

**Step 6: Commit**

```bash
git add packages/cfomics/R/methods_hdml.R packages/cfomics/R/cf_fit.R packages/cfomics/tests/testthat/test-hdml.R
git commit -m "feat(hdml): add High-Dimensional Machine Learning method

- Implement cf_fit_hdml() using glmnet for regularized estimation
- AIPW/doubly-robust estimator for ATE
- Support lasso, ridge, and elastic_net penalties
- Register method in cf_fit() switch statement"
```

---

### Task 1.2: Add HDML to DESCRIPTION and select_best_method

**Files:**
- Modify: `packages/cfomics/DESCRIPTION:52-56`
- Modify: `packages/cfomics/R/cf_fit.R:293-321` (select_best_method)

**Step 1: Update DESCRIPTION**

Add glmnet to Suggests (it may already be there, verify first):

```
Suggests:
    ...
    glmnet,
    hdm,
    ...
```

**Step 2: Update select_best_method to consider HDML**

Add to `packages/cfomics/R/cf_fit.R` in select_best_method():

```r
# After checking grf, before Python methods
if (requireNamespace("glmnet", quietly = TRUE)) {
  # For high-dimensional data, prefer hdml
  # This is a simple heuristic; could be enhanced with data inspection
  return("hdml")
}
```

**Step 3: Run R CMD check**

Run: `Rscript -e 'devtools::check("packages/cfomics", args = "--no-manual")'`
Expected: 0 errors, 0 warnings (some notes OK)

**Step 4: Commit**

```bash
git add packages/cfomics/DESCRIPTION packages/cfomics/R/cf_fit.R
git commit -m "feat(hdml): register in DESCRIPTION and method selection"
```

---

## Phase 2: HDPS (High-Dimensional Propensity Score)

### Task 2.1: Create methods_hdps.R

**Files:**
- Create: `packages/cfomics/R/methods_hdps.R`
- Test: `packages/cfomics/tests/testthat/test-hdps.R`

**Step 1: Write the failing test**

Create `packages/cfomics/tests/testthat/test-hdps.R`:

```r
test_that("cf_fit_hdps basic execution", {
skip_if_not_installed("glmnet")

set.seed(123)
n <- 200
p <- 50
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)
beta_t <- c(rep(0.3, 5), rep(0, p - 5))
ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)
Y <- 2.0 * T + 0.5 * X[,1] + 0.3 * X[,2] + rnorm(n)
data <- data.frame(Y = Y, T = T, X)

fit <- cf_fit(Y ~ T | ., data = data, method = "hdps")

expect_s3_class(fit, "cf_model")
expect_s3_class(fit, "cfomics_result")
expect_equal(fit$method, "hdps")

ate <- predict(fit, type = "ate")
expect_type(ate, "double")
expect_true(ate > 1.0 && ate < 3.0)
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-hdps.R")'`
Expected: FAIL with "Unknown method: hdps"

**Step 3: Write minimal implementation**

Create `packages/cfomics/R/methods_hdps.R`:

```r
#' Fit High-Dimensional Propensity Score model
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param lambda Character or numeric, lambda selection
#' @param trim Numeric vector of length 2, propensity score trimming bounds
#' @param nfolds Integer, number of CV folds
#' @param ... Additional arguments
#' @return List with model fit and results
#' @keywords internal
cf_fit_hdps <- function(X, T, Y,
                        covariate_names = NULL,
                        lambda = "cv",
                        trim = c(0.01, 0.99),
                        nfolds = 10L,
                        ...) {
if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Package 'glmnet' is required for method='hdps'. Please install it.")
}

X_mat <- as.matrix(X)
W <- as.numeric(T)
Y_vec <- as.numeric(Y)
n <- length(Y_vec)

if (is.null(covariate_names)) {
  covariate_names <- paste0("X", seq_len(ncol(X_mat)))
}

# Step 1: Estimate propensity score with LASSO logistic regression
ps_cv <- glmnet::cv.glmnet(X_mat, W, family = "binomial",
                           alpha = 1, nfolds = nfolds)
ps_lambda <- if (lambda == "1se") ps_cv$lambda.1se else ps_cv$lambda.min
ps_hat <- as.numeric(predict(ps_cv, newx = X_mat, s = ps_lambda, type = "response"))

# Trim propensity scores
ps_hat <- pmax(pmin(ps_hat, trim[2]), trim[1])

# Step 2: Compute IPW weights
weights <- W / ps_hat + (1 - W) / (1 - ps_hat)

# Step 3: Weighted mean difference
y1_weighted <- sum(W * Y_vec / ps_hat) / sum(W / ps_hat)
y0_weighted <- sum((1 - W) * Y_vec / (1 - ps_hat)) / sum((1 - W) / (1 - ps_hat))
ate <- y1_weighted - y0_weighted

# Step 4: Variance estimation via influence function
mu1 <- W * Y_vec / ps_hat
mu0 <- (1 - W) * Y_vec / (1 - ps_hat)
influence_fn <- mu1 - mu0 - ate
ate_se <- sqrt(var(influence_fn) / n)

ate_ci_lower <- ate - 1.96 * ate_se
ate_ci_upper <- ate + 1.96 * ate_se

# IPW doesn't estimate ITE, use ATE for all
ite <- rep(ate, n)
y0_hat <- Y_vec - W * ate
y1_hat <- Y_vec + (1 - W) * ate

summary_stats <- list(
  ate = ate,
  ate_ci_lower = ate_ci_lower,
  ate_ci_upper = ate_ci_upper,
  ite_mean = ate,
  ite_std = 0,
  ite_quantiles = list(
    q05 = ate,
    q50 = ate,
    q95 = ate
  )
)

list(
  model = list(
    ps_model = ps_cv,
    ps_lambda = ps_lambda
  ),
  res = list(
    ite = ite,
    ate = ate,
    y0_hat = y0_hat,
    y1_hat = y1_hat,
    summary = summary_stats,
    samples = list(
      y0 = y0_hat,
      y1 = y1_hat,
      ite = ite
    ),
    propensity = ps_hat,
    weights = weights,
    metadata = list(
      lambda = lambda,
      trim = trim,
      n_trimmed = sum(ps_hat == trim[1] | ps_hat == trim[2])
    )
  )
)
}

#' Predict from HDPS model
#'
#' @param object cf_model object with method="hdps"
#' @param newdata Not supported for HDPS
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_hdps <- function(object, newdata = NULL, type = "ite", ...) {
if (!is.null(newdata)) {
  stop("HDPS does not support prediction on new data")
}

res <- object$fit$res
switch(
  type,
  "ate" = res$ate,
  "ite" = res$ite,
  "y0" = res$y0_hat,
  "y1" = res$y1_hat,
  "summary" = res$summary,
  "samples" = res$samples
)
}
```

**Step 4: Register method in cf_fit.R**

Add to switch statement in `packages/cfomics/R/cf_fit.R`:

```r
"hdps" = cf_fit_hdps(X = X, T = T, Y = Y,
                     covariate_names = parsed$covariate_names,
                     ...),
```

Update method argument to include "hdps".

**Step 5: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-hdps.R")'`
Expected: PASS

**Step 6: Commit**

```bash
git add packages/cfomics/R/methods_hdps.R packages/cfomics/R/cf_fit.R packages/cfomics/tests/testthat/test-hdps.R
git commit -m "feat(hdps): add High-Dimensional Propensity Score method

- Implement cf_fit_hdps() using glmnet for PS estimation
- IPW estimator with propensity score trimming
- Influence function-based variance estimation"
```

---

## Phase 3: BCF (Bayesian Causal Forests)

### Task 3.1: Create methods_bcf.R

**Files:**
- Create: `packages/cfomics/R/methods_bcf.R`
- Test: `packages/cfomics/tests/testthat/test-bcf.R`

**Step 1: Write the failing test**

Create `packages/cfomics/tests/testthat/test-bcf.R`:

```r
test_that("cf_fit_bcf basic execution", {
skip_if_not_installed("bcf")

set.seed(123)
n <- 200
p <- 10
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)
ps <- plogis(0.5 * X[,1])
T <- rbinom(n, 1, ps)
tau_true <- 2.0 + 0.5 * X[,1]
Y <- tau_true * T + 0.5 * X[,1] + 0.3 * X[,2] + rnorm(n)
data <- data.frame(Y = Y, T = T, X)

fit <- cf_fit(Y ~ T | X1 + X2 + X3, data = data, method = "bcf",
              n_burn = 100, n_iter = 200)

expect_s3_class(fit, "cf_model")
expect_s3_class(fit, "cfomics_result")
expect_equal(fit$method, "bcf")

ate <- predict(fit, type = "ate")
expect_type(ate, "double")
expect_true(ate > 1.5 && ate < 3.5)

ite <- predict(fit, type = "ite")
expect_equal(length(ite), n)
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-bcf.R")'`
Expected: FAIL with "Unknown method: bcf"

**Step 3: Write minimal implementation**

Create `packages/cfomics/R/methods_bcf.R`:

```r
#' Fit Bayesian Causal Forests model
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param n_burn Integer, number of burn-in iterations
#' @param n_iter Integer, number of MCMC iterations after burn-in
#' @param n_chains Integer, number of chains (currently only 1 supported)
#' @param pihat Optional propensity scores; if NULL, estimated internally
#' @param ... Additional arguments passed to bcf()
#' @return List with model fit and results
#' @keywords internal
cf_fit_bcf <- function(X, T, Y,
                       covariate_names = NULL,
                       n_burn = 1000L,
                       n_iter = 2000L,
                       n_chains = 1L,
                       pihat = NULL,
                       ...) {
if (!requireNamespace("bcf", quietly = TRUE)) {
  stop("Package 'bcf' is required for method='bcf'. Please install it.")
}

X_mat <- as.matrix(X)
W <- as.numeric(T)
Y_vec <- as.numeric(Y)
n <- length(Y_vec)

if (is.null(covariate_names)) {
  covariate_names <- paste0("X", seq_len(ncol(X_mat)))
}

# Estimate propensity score if not provided
if (is.null(pihat)) {
  if (requireNamespace("glmnet", quietly = TRUE)) {
    ps_cv <- glmnet::cv.glmnet(X_mat, W, family = "binomial", alpha = 1)
    pihat <- as.numeric(predict(ps_cv, newx = X_mat, s = "lambda.min",
                                 type = "response"))
  } else {
    # Simple logistic regression fallback
    df_ps <- data.frame(W = W, X_mat)
    ps_model <- glm(W ~ ., data = df_ps, family = binomial)
    pihat <- predict(ps_model, type = "response")
  }
  pihat <- pmax(pmin(pihat, 0.99), 0.01)
}

# Fit BCF
bcf_fit <- bcf::bcf(
  y = Y_vec,
  z = W,
  x_control = X_mat,
  x_moderate = X_mat,
  pihat = pihat,
  nburn = as.integer(n_burn),
  nsim = as.integer(n_iter),
  ...
)

# Extract treatment effects
# bcf returns tau samples of dimension (n_iter x n)
tau_samples <- bcf_fit$tau
ite <- colMeans(tau_samples)
ate <- mean(ite)

# Credible intervals from posterior
ate_samples <- rowMeans(tau_samples)
ate_ci_lower <- quantile(ate_samples, 0.025)
ate_ci_upper <- quantile(ate_samples, 0.975)

# Posterior predictive for counterfactuals
# mu_samples: prognostic function samples
mu_samples <- bcf_fit$yhat - W * tau_samples
y0_hat <- colMeans(mu_samples)
y1_hat <- y0_hat + ite

summary_stats <- list(
  ate = ate,
  ate_ci_lower = ate_ci_lower,
  ate_ci_upper = ate_ci_upper,
  ite_mean = mean(ite),
  ite_std = sd(ite),
  ite_quantiles = list(
    q05 = quantile(ite, 0.05),
    q50 = quantile(ite, 0.50),
    q95 = quantile(ite, 0.95)
  )
)

# ITE credible intervals
ite_ci_lower <- apply(tau_samples, 2, quantile, 0.025)
ite_ci_upper <- apply(tau_samples, 2, quantile, 0.975)

list(
  model = bcf_fit,
  res = list(
    ite = ite,
    ate = ate,
    y0_hat = y0_hat,
    y1_hat = y1_hat,
    summary = summary_stats,
    samples = list(
      y0 = y0_hat,
      y1 = y1_hat,
      ite = ite,
      tau_posterior = tau_samples,
      ate_posterior = ate_samples
    ),
    credible_intervals = list(
      ite_lower = ite_ci_lower,
      ite_upper = ite_ci_upper
    ),
    propensity = pihat,
    metadata = list(
      n_burn = n_burn,
      n_iter = n_iter
    )
  )
)
}

#' Predict from BCF model
#'
#' @param object cf_model object with method="bcf"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_bcf <- function(object, newdata = NULL, type = "ite", ...) {
if (is.null(newdata)) {
  res <- object$fit$res
  return(switch(
    type,
    "ate" = res$ate,
    "ite" = res$ite,
    "y0" = res$y0_hat,
    "y1" = res$y1_hat,
    "summary" = res$summary,
    "samples" = res$samples
  ))
}

# BCF prediction on new data requires re-running with predict.bcf
# This is computationally expensive; for now, not fully supported
stop("BCF prediction on new data not yet implemented")
}
```

**Step 4: Register method in cf_fit.R**

Add to switch statement:

```r
"bcf" = cf_fit_bcf(X = X, T = T, Y = Y,
                   covariate_names = parsed$covariate_names,
                   ...),
```

Update DESCRIPTION to add bcf to Suggests.

**Step 5: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-bcf.R")'`
Expected: PASS (or SKIP if bcf not installed)

**Step 6: Commit**

```bash
git add packages/cfomics/R/methods_bcf.R packages/cfomics/R/cf_fit.R packages/cfomics/tests/testthat/test-bcf.R packages/cfomics/DESCRIPTION
git commit -m "feat(bcf): add Bayesian Causal Forests method

- Implement cf_fit_bcf() using bcf package
- Full posterior inference for ITE and ATE
- Credible intervals from MCMC samples
- Automatic propensity score estimation"
```

---

## Phase 4: TMLE (Targeted Maximum Likelihood Estimation)

### Task 4.1: Create methods_tmle.R

**Files:**
- Create: `packages/cfomics/R/methods_tmle.R`
- Test: `packages/cfomics/tests/testthat/test-tmle.R`

**Step 1: Write the failing test**

Create `packages/cfomics/tests/testthat/test-tmle.R`:

```r
test_that("cf_fit_tmle basic execution", {
skip_if_not_installed("glmnet")

set.seed(123)
n <- 200
p <- 20
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)
ps <- plogis(0.5 * X[,1] + 0.3 * X[,2])
T <- rbinom(n, 1, ps)
Y <- 2.0 * T + 0.5 * X[,1] + 0.3 * X[,2] + rnorm(n)
data <- data.frame(Y = Y, T = T, X)

fit <- cf_fit(Y ~ T | ., data = data, method = "tmle")

expect_s3_class(fit, "cf_model")
expect_s3_class(fit, "cfomics_result")
expect_equal(fit$method, "tmle")

ate <- predict(fit, type = "ate")
expect_type(ate, "double")
expect_true(ate > 1.0 && ate < 3.0)
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-tmle.R")'`
Expected: FAIL with "Unknown method: tmle"

**Step 3: Write minimal implementation**

Create `packages/cfomics/R/methods_tmle.R`:

```r
#' Fit Targeted Maximum Likelihood Estimation model
#'
#' Implements a basic TMLE estimator for ATE using glmnet for nuisance
#' parameter estimation. This is a simplified implementation; for full
#' Super Learner integration, consider using the tmle3 package directly.
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector
#' @param covariate_names Character vector of covariate names
#' @param sl_lib Character vector of Super Learner libraries (currently unused,
#'   for future compatibility with tmle3/sl3)
#' @param nfolds Integer, number of CV folds for nuisance estimation
#' @param ... Additional arguments
#' @return List with model fit and results
#' @keywords internal
cf_fit_tmle <- function(X, T, Y,
                        covariate_names = NULL,
                        sl_lib = NULL,
                        nfolds = 10L,
                        ...) {
if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Package 'glmnet' is required for method='tmle'. Please install it.")
}

X_mat <- as.matrix(X)
W <- as.numeric(T)
Y_vec <- as.numeric(Y)
n <- length(Y_vec)

if (is.null(covariate_names)) {
  covariate_names <- paste0("X", seq_len(ncol(X_mat)))
}

# Step 1: Initial outcome model Q(A, W)
# Fit outcome model using cross-validation
Q_cv <- glmnet::cv.glmnet(cbind(X_mat, W), Y_vec, alpha = 0.5, nfolds = nfolds)
Q_lambda <- Q_cv$lambda.min

# Get initial predictions
Q1_init <- as.numeric(predict(Q_cv, newx = cbind(X_mat, 1), s = Q_lambda))
Q0_init <- as.numeric(predict(Q_cv, newx = cbind(X_mat, 0), s = Q_lambda))
QA_init <- W * Q1_init + (1 - W) * Q0_init

# Step 2: Propensity score model g(W)
g_cv <- glmnet::cv.glmnet(X_mat, W, family = "binomial", alpha = 1, nfolds = nfolds)
g_lambda <- g_cv$lambda.min
g_hat <- as.numeric(predict(g_cv, newx = X_mat, s = g_lambda, type = "response"))
g_hat <- pmax(pmin(g_hat, 0.99), 0.01)  # Bound propensity scores

# Step 3: Clever covariate H(A, W)
H1 <- 1 / g_hat
H0 <- -1 / (1 - g_hat)
HA <- W * H1 + (1 - W) * H0

# Step 4: Targeting step - fit epsilon
# Regress Y - Q_init on H with offset Q_init
# Using logistic regression for bounded outcomes would be better,
# but for simplicity we use linear here
epsilon_model <- lm(Y_vec ~ -1 + offset(QA_init) + HA)
epsilon <- coef(epsilon_model)["HA"]
if (is.na(epsilon)) epsilon <- 0

# Step 5: Update predictions
Q1_star <- Q1_init + epsilon * H1
Q0_star <- Q0_init + epsilon * H0

# ATE estimate
ate <- mean(Q1_star - Q0_star)

# Step 6: Influence curve based variance
IC <- (Q1_star - Q0_star - ate) +
      H1 * W * (Y_vec - Q1_star) +
      H0 * (1 - W) * (Y_vec - Q0_star)
ate_se <- sqrt(var(IC) / n)

ate_ci_lower <- ate - 1.96 * ate_se
ate_ci_upper <- ate + 1.96 * ate_se

# ITE approximation (model-based, not true CATE)
ite <- Q1_star - Q0_star

summary_stats <- list(
  ate = ate,
  ate_ci_lower = ate_ci_lower,
  ate_ci_upper = ate_ci_upper,
  ite_mean = mean(ite),
  ite_std = sd(ite),
  ite_quantiles = list(
    q05 = quantile(ite, 0.05),
    q50 = quantile(ite, 0.50),
    q95 = quantile(ite, 0.95)
  )
)

list(
  model = list(
    Q_model = Q_cv,
    g_model = g_cv,
    epsilon = epsilon
  ),
  res = list(
    ite = ite,
    ate = ate,
    y0_hat = Q0_star,
    y1_hat = Q1_star,
    summary = summary_stats,
    samples = list(
      y0 = Q0_star,
      y1 = Q1_star,
      ite = ite
    ),
    propensity = g_hat,
    influence_curve = IC,
    metadata = list(
      epsilon = epsilon,
      sl_lib = sl_lib
    )
  )
)
}

#' Predict from TMLE model
#'
#' @param object cf_model object with method="tmle"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_tmle <- function(object, newdata = NULL, type = "ite", ...) {
if (is.null(newdata)) {
  res <- object$fit$res
  return(switch(
    type,
    "ate" = res$ate,
    "ite" = res$ite,
    "y0" = res$y0_hat,
    "y1" = res$y1_hat,
    "summary" = res$summary,
    "samples" = res$samples
  ))
}

# Prediction on new data
Q_model <- object$fit$model$Q_model
g_model <- object$fit$model$g_model
epsilon <- object$fit$model$epsilon

X_new <- as.matrix(newdata)

Q1_init <- as.numeric(predict(Q_model, newx = cbind(X_new, 1), s = "lambda.min"))
Q0_init <- as.numeric(predict(Q_model, newx = cbind(X_new, 0), s = "lambda.min"))

g_hat <- as.numeric(predict(g_model, newx = X_new, s = "lambda.min", type = "response"))
g_hat <- pmax(pmin(g_hat, 0.99), 0.01)

H1 <- 1 / g_hat
H0 <- -1 / (1 - g_hat)

Q1_star <- Q1_init + epsilon * H1
Q0_star <- Q0_init + epsilon * H0
ite <- Q1_star - Q0_star

switch(
  type,
  "ate" = mean(ite),
  "ite" = ite,
  "y0" = Q0_star,
  "y1" = Q1_star,
  stop(sprintf("type '%s' not supported", type))
)
}
```

**Step 4: Register method in cf_fit.R**

Add to switch statement:

```r
"tmle" = cf_fit_tmle(X = X, T = T, Y = Y,
                     covariate_names = parsed$covariate_names,
                     ...),
```

**Step 5: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-tmle.R")'`
Expected: PASS

**Step 6: Commit**

```bash
git add packages/cfomics/R/methods_tmle.R packages/cfomics/R/cf_fit.R packages/cfomics/tests/testthat/test-tmle.R
git commit -m "feat(tmle): add Targeted Maximum Likelihood Estimation method

- Implement cf_fit_tmle() using glmnet for nuisance parameters
- Doubly-robust targeting step
- Influence curve-based variance estimation
- Foundation for future Super Learner integration"
```

---

## Phase 5: Benchmark DGP Extensions

### Task 5.1: Add baseline DGP function

**Files:**
- Modify: `packages/cfomics/R/benchmark_dgp.R`
- Test: `packages/cfomics/tests/testthat/test-benchmark-dgp.R`

**Step 1: Write the failing test**

Add to `packages/cfomics/tests/testthat/test-benchmark-dgp.R`:

```r
test_that("dgp_baseline generates valid data", {
dgp <- dgp_baseline(n = 100, p = 50)

expect_type(dgp, "list")
expect_equal(nrow(dgp$X), 100)
expect_equal(ncol(dgp$X), 50)
expect_equal(length(dgp$T), 100)
expect_equal(length(dgp$Y), 100)
expect_true("true_ate" %in% names(dgp))
expect_true("true_ite" %in% names(dgp))
expect_true("propensity_score" %in% names(dgp))
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-benchmark-dgp.R")'`
Expected: FAIL with "dgp_baseline not found"

**Step 3: Write minimal implementation**

Add to `packages/cfomics/R/benchmark_dgp.R`:

```r
#' Generate baseline DGP data (S1)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, propensity_score
#' @export
dgp_baseline <- function(n = 500, p = 50, seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

# Sparse confounding: first 10 variables
n_conf <- min(10, p)
beta_t <- c(rep(0.3, n_conf), rep(0, p - n_conf))
beta_y <- c(rep(0.2, n_conf), rep(0, p - n_conf))

ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)

tau <- 2.0
Y <- X %*% beta_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = tau,
  true_ite = rep(tau, n),
  propensity_score = as.numeric(ps),
  dgp_name = "baseline",
  dgp_params = list(n = n, p = p)
)
}
```

**Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("packages/cfomics/tests/testthat/test-benchmark-dgp.R")'`
Expected: PASS

**Step 5: Commit**

```bash
git add packages/cfomics/R/benchmark_dgp.R packages/cfomics/tests/testthat/test-benchmark-dgp.R
git commit -m "feat(dgp): add baseline DGP function (S1 scenario)"
```

---

### Task 5.2: Add dimension sweep DGP

**Files:**
- Modify: `packages/cfomics/R/benchmark_dgp.R`
- Test: `packages/cfomics/tests/testthat/test-benchmark-dgp.R`

**Step 1: Add test**

```r
test_that("dgp_dimension_sweep generates valid high-dim data", {
dgp <- dgp_dimension_sweep(n = 100, p = 500)

expect_equal(nrow(dgp$X), 100)
expect_equal(ncol(dgp$X), 500)
expect_true(dgp$dgp_params$n_p_ratio < 1)
})
```

**Step 2: Write implementation**

Add to `packages/cfomics/R/benchmark_dgp.R`:

```r
#' Generate high-dimensional DGP data (S2-S3 dimension sweep)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, propensity_score
#' @export
dgp_dimension_sweep <- function(n = 200, p = 1000, seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

n_conf <- min(10, p)
beta_t <- c(rep(0.5, n_conf), rep(0, p - n_conf))
beta_y <- c(rep(0.3, n_conf), rep(0, p - n_conf))

ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)

tau <- 2.0
Y <- X %*% beta_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = tau,
  true_ite = rep(tau, n),
  propensity_score = as.numeric(ps),
  dgp_name = "dimension_sweep",
  dgp_params = list(n = n, p = p, n_p_ratio = n / p)
)
}
```

**Step 3: Run test and commit**

```bash
git add packages/cfomics/R/benchmark_dgp.R packages/cfomics/tests/testthat/test-benchmark-dgp.R
git commit -m "feat(dgp): add dimension sweep DGP (S2-S3)"
```

---

### Task 5.3: Add heterogeneous treatment effect DGPs (S4)

**Files:**
- Modify: `packages/cfomics/R/benchmark_dgp.R`

**Step 1: Write tests**

```r
test_that("dgp_heterogeneous_linear generates valid HTE data", {
dgp <- dgp_heterogeneous_linear(n = 100, p = 50, strength = 1.0)

expect_equal(length(unique(dgp$true_ite)), 100)  # ITE varies
expect_true(sd(dgp$true_ite) > 0)
})

test_that("dgp_heterogeneous_subgroup generates subgroup data", {
dgp <- dgp_heterogeneous_subgroup(n = 300, p = 50)

expect_true("subgroup" %in% names(dgp))
expect_equal(length(dgp$subgroup), 300)
})
```

**Step 2: Write implementations**

Add to `packages/cfomics/R/benchmark_dgp.R`:

```r
#' Generate linear heterogeneous treatment effect DGP (S4a/e/f)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param strength Numeric, heterogeneity strength
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_heterogeneous_linear <- function(n = 500, p = 500, strength = 1.0, seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

beta_t <- c(rep(0.3, 10), rep(0, p - 10))
ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)

tau_base <- 2.0
tau <- tau_base + strength * X[,1] + 0.3 * strength * X[,2]

beta_y <- c(rep(0.2, 10), rep(0, p - 10))
Y <- X %*% beta_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = mean(tau),
  true_ite = tau,
  ite_sd = sd(tau),
  propensity_score = as.numeric(ps),
  dgp_name = "heterogeneous_linear",
  dgp_params = list(n = n, p = p, strength = strength)
)
}

#' Generate nonlinear heterogeneous treatment effect DGP (S4b)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_heterogeneous_nonlinear <- function(n = 500, p = 500, seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

beta_t <- c(rep(0.3, 10), rep(0, p - 10))
ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)

tau <- 2.0 + 1.0 * sin(pi * X[,1]) + 0.5 * X[,2]^2

beta_y <- c(rep(0.2, 10), rep(0, p - 10))
Y <- X %*% beta_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = mean(tau),
  true_ite = tau,
  propensity_score = as.numeric(ps),
  dgp_name = "heterogeneous_nonlinear",
  dgp_params = list(n = n, p = p)
)
}

#' Generate subgroup heterogeneous treatment effect DGP (S4c/g/h)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param subgroup_props Numeric vector, proportion in each subgroup
#' @param subgroup_effects Numeric vector, treatment effect in each subgroup
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, subgroup
#' @export
dgp_heterogeneous_subgroup <- function(n = 500, p = 500,
                                       subgroup_props = c(1/3, 1/3, 1/3),
                                       subgroup_effects = c(0.5, 2.0, 4.0),
                                       seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

beta_t <- c(rep(0.3, 10), rep(0, p - 10))
ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)

# Subgroup assignment based on X1
breaks <- qnorm(cumsum(subgroup_props)[1:(length(subgroup_props)-1)])
subgroup <- cut(X[,1], breaks = c(-Inf, breaks, Inf), labels = FALSE)

tau <- subgroup_effects[subgroup]

beta_y <- c(rep(0.2, 10), rep(0, p - 10))
Y <- X %*% beta_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = sum(subgroup_props * subgroup_effects),
  true_ite = tau,
  subgroup = subgroup,
  subgroup_n = table(subgroup),
  propensity_score = as.numeric(ps),
  dgp_name = "heterogeneous_subgroup",
  dgp_params = list(n = n, p = p, subgroup_props = subgroup_props,
                    subgroup_effects = subgroup_effects)
)
}

#' Generate qualitative interaction DGP (S4d)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, prop_harmed
#' @export
dgp_heterogeneous_qualitative <- function(n = 500, p = 500, seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

beta_t <- c(rep(0.3, 10), rep(0, p - 10))
ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)

tau <- 2.0 * X[,1]  # Negative for X1 < 0

beta_y <- c(rep(0.2, 10), rep(0, p - 10))
Y <- X %*% beta_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = mean(tau),
  true_ite = tau,
  prop_harmed = mean(tau < 0),
  prop_benefited = mean(tau > 0),
  propensity_score = as.numeric(ps),
  dgp_name = "heterogeneous_qualitative",
  dgp_params = list(n = n, p = p)
)
}
```

**Step 3: Run tests and commit**

```bash
git add packages/cfomics/R/benchmark_dgp.R packages/cfomics/tests/testthat/test-benchmark-dgp.R
git commit -m "feat(dgp): add heterogeneous treatment effect DGPs (S4a-h)"
```

---

### Task 5.4: Add nonlinear confounding DGP (S5)

Add to `packages/cfomics/R/benchmark_dgp.R`:

```r
#' Generate nonlinear confounding DGP (S5)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param nonlinear_type Character: "quadratic", "trigonometric", "interaction",
#'   "combined", "threshold"
#' @param strength Numeric, nonlinearity strength (0=linear, 1=full nonlinear)
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_nonlinear_confounding <- function(n = 500, p = 500,
                                       nonlinear_type = "combined",
                                       strength = 1.0,
                                       seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

linear_t <- 0.3 * X[,1] + 0.3 * X[,2]
linear_y <- 0.2 * X[,1] + 0.2 * X[,2]

nonlinear_t <- nonlinear_y <- 0

if (nonlinear_type %in% c("quadratic", "combined")) {
  nonlinear_t <- nonlinear_t + 0.5 * X[,1]^2
  nonlinear_y <- nonlinear_y + 0.3 * X[,1]^2
}

if (nonlinear_type %in% c("trigonometric", "combined")) {
  nonlinear_t <- nonlinear_t + 0.3 * sin(pi * X[,3])
  nonlinear_y <- nonlinear_y + 0.2 * sin(pi * X[,2])
}

if (nonlinear_type %in% c("interaction", "combined")) {
  nonlinear_t <- nonlinear_t + 0.2 * X[,1] * X[,2]
  nonlinear_y <- nonlinear_y + 0.1 * X[,1] * X[,3]
}

if (nonlinear_type == "threshold") {
  nonlinear_t <- 0.5 * as.numeric(X[,1] > 0)
  nonlinear_y <- 0.3 * as.numeric(X[,2] > 0.5)
}

ps <- plogis(linear_t + strength * nonlinear_t)
T <- rbinom(n, 1, ps)

tau <- 2.0
Y <- linear_y + strength * nonlinear_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = tau,
  true_ite = rep(tau, n),
  propensity_score = as.numeric(ps),
  dgp_name = "nonlinear_confounding",
  dgp_params = list(n = n, p = p, nonlinear_type = nonlinear_type, strength = strength)
)
}
```

**Commit:**

```bash
git add packages/cfomics/R/benchmark_dgp.R
git commit -m "feat(dgp): add nonlinear confounding DGP (S5)"
```

---

### Task 5.5: Add dense confounding DGP (S6)

```r
#' Generate dense confounding DGP (S6)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param n_confounders Integer, number of confounding variables
#' @param coef_scaling Character: "fixed", "sqrt", "linear"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite
#' @export
dgp_dense_confounding <- function(n = 500, p = 500,
                                   n_confounders = 100,
                                   coef_scaling = "fixed",
                                   seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

base_coef <- 0.1
coef_size <- switch(coef_scaling,
  "fixed"  = base_coef,
  "sqrt"   = base_coef / sqrt(n_confounders),
  "linear" = base_coef * 10 / n_confounders,
  base_coef
)

beta_t <- c(rep(coef_size, n_confounders), rep(0, p - n_confounders))
beta_y <- c(rep(coef_size, n_confounders), rep(0, p - n_confounders))

ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)

tau <- 2.0
Y <- X %*% beta_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = tau,
  true_ite = rep(tau, n),
  propensity_score = as.numeric(ps),
  confounding_info = list(
    n_confounders = n_confounders,
    coef_size = coef_size,
    coef_scaling = coef_scaling
  ),
  dgp_name = "dense_confounding",
  dgp_params = list(n = n, p = p, n_confounders = n_confounders,
                    coef_scaling = coef_scaling)
)
}
```

---

### Task 5.6: Add weak overlap DGP (S7)

```r
#' Generate weak overlap DGP (S7a)
#'
#' @param n Integer, number of observations
#' @param p Integer, number of covariates
#' @param overlap_strength Character: "good", "moderate", "weak", "extreme"
#' @param seed Integer, random seed
#' @return List with X, T, Y, true_ate, true_ite, ps_summary
#' @export
dgp_weak_overlap <- function(n = 500, p = 500,
                             overlap_strength = "weak",
                             seed = NULL) {
if (!is.null(seed)) set.seed(seed)

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("X", 1:p)

coef_scale <- switch(overlap_strength,
  "good"     = 0.3,
  "moderate" = 0.6,
  "weak"     = 1.2,
  "extreme"  = 2.0,
  1.2
)

beta_t <- c(rep(coef_scale, 10), rep(0, p - 10))
ps <- plogis(X %*% beta_t)
T <- rbinom(n, 1, ps)

ps_summary <- list(
  min = min(ps),
  max = max(ps),
  mean = mean(ps),
  prop_extreme = mean(ps < 0.05 | ps > 0.95)
)

tau <- 2.0
beta_y <- c(rep(0.3, 10), rep(0, p - 10))
Y <- X %*% beta_y + tau * T + rnorm(n)

list(
  X = X,
  T = T,
  Y = as.numeric(Y),
  true_ate = tau,
  true_ite = rep(tau, n),
  propensity_score = as.numeric(ps),
  ps_summary = ps_summary,
  overlap_strength = overlap_strength,
  dgp_name = "weak_overlap",
  dgp_params = list(n = n, p = p, overlap_strength = overlap_strength)
)
}
```

---

### Task 5.7-5.10: Add remaining DGPs (S8-S11)

Continue adding:
- `dgp_covariate_shift()` (S8)
- `dgp_correlated_confounding_block()` (S9a-c)
- `dgp_correlated_confounding_ar1()` (S9d)
- `dgp_correlated_confounding_factor()` (S9e)
- `dgp_unobserved_confounding()` (S10)
- `dgp_collider()` (S11)

Each follows the same pattern. See the design document for full implementations.

**Final commit for DGPs:**

```bash
git add packages/cfomics/R/benchmark_dgp.R packages/cfomics/tests/testthat/test-benchmark-dgp.R
git commit -m "feat(dgp): add complete benchmark DGP suite (S1-S11)

- S1: Baseline
- S2-S3: Dimension sweep
- S4a-h: Heterogeneous treatment effects
- S5a-e: Nonlinear confounding
- S6: Dense confounding
- S7a-c: Weak overlap
- S8a-d: Covariate shift
- S9a-e: Correlated confounding
- S10a-e: Unobserved confounding
- S11a-e: Collider bias"
```

---

## Phase 6: Benchmark Sweep Functions

### Task 6.1: Add sweep experiment functions

**Files:**
- Create: `packages/cfomics/R/benchmark_sweep.R`
- Test: `packages/cfomics/tests/testthat/test-benchmark-sweep.R`

Add functions:
- `cf_benchmark_dimension_sweep()`
- `cf_benchmark_heterogeneity_sweep()`
- `cf_benchmark_nonlinearity_sweep()`
- `cf_benchmark_density_sweep()`
- `cf_benchmark_overlap_sweep()`
- `cf_benchmark_covariate_shift_sweep()`
- `cf_benchmark_correlation_sweep()`
- `cf_benchmark_unmeasured_sweep()`
- `cf_benchmark_collider_sweep()`

See design document for implementations.

---

## Phase 7: Final Integration and Documentation

### Task 7.1: Update NAMESPACE

Run: `Rscript -e 'devtools::document("packages/cfomics")'`

### Task 7.2: Update DESCRIPTION

Ensure all new dependencies are in Suggests:
- glmnet
- bcf
- MASS
- hdm (optional)

### Task 7.3: Run full test suite

Run: `Rscript -e 'devtools::test("packages/cfomics")'`

### Task 7.4: Run R CMD check

Run: `Rscript -e 'devtools::check("packages/cfomics")'`

### Task 7.5: Final commit

```bash
git add .
git commit -m "feat: complete high-dimensional methods and benchmark suite

Added methods:
- hdml: High-Dimensional Machine Learning (Debiased Lasso)
- hdps: High-Dimensional Propensity Score
- bcf: Bayesian Causal Forests
- tmle: Targeted Maximum Likelihood Estimation

Added DGPs:
- 39 sub-scenarios covering S1-S11
- 9 sweep experiment functions

All tests passing, R CMD check clean."
```

---

## Summary

| Phase | Tasks | Est. Steps |
|-------|-------|------------|
| 1. HDML | 2 | ~12 |
| 2. HDPS | 1 | ~6 |
| 3. BCF | 1 | ~6 |
| 4. TMLE | 1 | ~6 |
| 5. DGPs | 10 | ~30 |
| 6. Sweeps | 1 | ~6 |
| 7. Integration | 5 | ~10 |
| **Total** | **21** | **~76** |

Each step is designed to be 2-5 minutes of focused work.
