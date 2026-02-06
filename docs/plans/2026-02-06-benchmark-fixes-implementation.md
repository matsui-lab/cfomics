# Benchmark Implementation Fixes Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix critical and important issues identified in the benchmark implementation verification.

**Architecture:** Address issues in order of severity (Critical → Important) across four areas: external wrappers, TCGA pipeline, NN models, and metrics. Each fix includes test verification where applicable.

**Tech Stack:** R (testthat), Python (PyTorch), ggplot2/dplyr for visualization

---

## Phase 1: Critical Fixes (External Wrappers)

### Task 1: Fix MatchIt wrapper weights reference

**Files:**
- Modify: `benchmarks/external/wrappers/wrapper_matchit.R:61`
- Test: `benchmarks/external/wrappers/tests/test_matchit.R`

**Step 1: Read current implementation**

```bash
cat benchmarks/external/wrappers/wrapper_matchit.R | head -70
```

**Step 2: Fix the weights reference bug**

Change line 61 from:
```r
fit <- lm(out_formula, data = m_data, weights = weights)
```

To:
```r
fit <- lm(out_formula, data = m_data, weights = m_data$weights)
```

**Step 3: Run tests to verify fix**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_matchit.R')"
```

Expected: All tests PASS

**Step 4: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_matchit.R
git commit -m "fix(benchmark): correct weights reference in MatchIt wrapper"
```

---

### Task 2: Fix WeightIt wrapper weights reference

**Files:**
- Modify: `benchmarks/external/wrappers/wrapper_weightit.R:43`
- Test: `benchmarks/external/wrappers/tests/test_weightit.R`

**Step 1: Read current implementation**

```bash
cat benchmarks/external/wrappers/wrapper_weightit.R | head -50
```

**Step 2: Fix the weights reference bug**

Change line 43 from:
```r
fit <- lm(out_formula, data = df, weights = weights)
```

To:
```r
fit <- lm(out_formula, data = df, weights = df$weights)
```

**Step 3: Run tests to verify fix**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_weightit.R')"
```

Expected: All tests PASS

**Step 4: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_weightit.R
git commit -m "fix(benchmark): correct weights reference in WeightIt wrapper"
```

---

### Task 3: Fix DragonNet missing targeted regularization (document as simplified)

**Files:**
- Modify: `benchmarks/external/nn/models.py:140-193` (DragonNet class)
- Modify: `benchmarks/external/nn/models.py:240-250` (train_model DragonNet branch)

**Step 1: Add documentation noting this is a simplified implementation**

Update DragonNet class docstring (around line 140):

```python
class DragonNet(nn.Module):
    """DragonNet: propensity score co-learning (Simplified Version).

    Based on Shi et al. (2019) "Adapting Neural Networks for Treatment Effects".

    NOTE: This is a simplified implementation that includes joint propensity
    learning but omits the targeted regularization loss from the original paper.
    For full DragonNet, see the official implementation at:
    https://github.com/claudiashi57/dragonnet

    The simplified version still provides benefits from representation sharing
    between outcome and propensity models, but may have slightly higher bias
    than the full targeted regularization approach.

    Args:
        input_dim: Number of input features
        hidden_dim: Hidden layer dimension (default: 200)
        n_layers: Number of representation layers (default: 3)
    """
```

**Step 2: Update training function comment**

Update the DragonNet branch in train_model (around line 240):

```python
            if isinstance(model, DragonNet):
                # Simplified DragonNet loss: outcome + propensity (no targeted regularization)
                y_pred, (y0, y1), ps = model(x_b, t_b)
                loss_y = nn.MSELoss()(y_pred, y_b)
                loss_ps = nn.BCELoss()(ps, t_b)
                loss = loss_y + loss_ps
```

**Step 3: Run Python syntax check**

```bash
python -m py_compile benchmarks/external/nn/models.py
```

Expected: No errors

**Step 4: Run quick test**

```bash
cd benchmarks/external/nn && python -c "from models import DragonNet; print('DragonNet loads OK')"
```

**Step 5: Commit**

```bash
git add benchmarks/external/nn/models.py
git commit -m "docs(benchmark): document DragonNet as simplified implementation without targeted regularization"
```

---

## Phase 2: Important Fixes (TCGA Pipeline)

### Task 4: Add PCA bounds checking in preprocess_tcga.R

**Files:**
- Modify: `benchmarks/tcga/preprocess_tcga.R:49-54`

**Step 1: Read current PCA code**

```bash
sed -n '45,60p' benchmarks/tcga/preprocess_tcga.R
```

**Step 2: Add bounds checking**

Replace lines 49-54:
```r
  # Use PCA for dimension reduction if needed
  if (p > 100) {
    pca <- prcomp(X, center = TRUE, scale. = TRUE)
    X_reduced <- pca$x[, 1:100]
  } else {
    X_reduced <- scale(X)
  }
```

With:
```r
  # Use PCA for dimension reduction if needed
  n_target_components <- 100
  if (p > n_target_components) {
    pca <- prcomp(X, center = TRUE, scale. = TRUE)
    # PCA returns min(n-1, p) components; ensure we don't exceed available
    n_available <- ncol(pca$x)
    n_components <- min(n_target_components, n_available)
    X_reduced <- pca$x[, 1:n_components, drop = FALSE]
    message(sprintf("PCA: reduced %d features to %d components", p, n_components))
  } else {
    X_reduced <- scale(X)
  }
```

**Step 3: Commit**

```bash
git add benchmarks/tcga/preprocess_tcga.R
git commit -m "fix(benchmark): add PCA bounds checking for small sample sizes"
```

---

### Task 5: Improve seed strategy in run_tcga_benchmark.R

**Files:**
- Modify: `benchmarks/tcga/run_tcga_benchmark.R:187`

**Step 1: Find seed line**

```bash
grep -n "set.seed" benchmarks/tcga/run_tcga_benchmark.R
```

**Step 2: Improve seed strategy**

Change from:
```r
set.seed(job$rep)
```

To:
```r
# Use robust seed: base_seed * multiplier + rep offset for better randomness
set.seed(20260206 * 1000 + job$rep)
```

**Step 3: Commit**

```bash
git add benchmarks/tcga/run_tcga_benchmark.R
git commit -m "fix(benchmark): improve seed strategy for TCGA benchmark reproducibility"
```

---

## Phase 3: Important Fixes (NN Models)

### Task 6: Fix IPM loss documentation in CFRNet

**Files:**
- Modify: `benchmarks/external/nn/models.py:117-137`

**Step 1: Update docstring to correctly describe the loss**

Change the ipm_loss method docstring from:
```python
    def ipm_loss(self, phi: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        """Wasserstein-1 (MMD) approximation for IPM."""
```

To:
```python
    def ipm_loss(self, phi: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        """Mean-matching IPM loss (L2 distance between group means).

        This is a simple linear IPM that penalizes differences in the mean
        representations between treated and control groups. It is computationally
        efficient but less powerful than kernel-based MMD or Wasserstein distance.

        For the original CFRNet paper's MMD implementation, see:
        https://github.com/clinicalml/cfrnet

        Args:
            phi: Representation tensor of shape (batch_size, hidden_dim)
            t: Treatment indicator tensor of shape (batch_size,)

        Returns:
            Scalar tensor with the L2 distance between group means
        """
```

**Step 2: Commit**

```bash
git add benchmarks/external/nn/models.py
git commit -m "docs(benchmark): correct IPM loss documentation in CFRNet (mean-matching, not MMD)"
```

---

### Task 7: Optimize predict_ite to avoid double forward pass

**Files:**
- Modify: `benchmarks/external/nn/models.py:273-310` (predict_ite function)

**Step 1: Read current predict_ite**

```bash
sed -n '273,320p' benchmarks/external/nn/models.py
```

**Step 2: Optimize to single forward pass**

Replace the prediction logic:
```python
def predict_ite(model, X):
    """Predict ITE using trained model.

    Args:
        model: Trained TARNet, CFRNet, or DragonNet model
        X: Covariate matrix (numpy array)

    Returns:
        Dictionary with 'ite', 'ate', 'y0_hat', 'y1_hat'
    """
    device = next(model.parameters()).device
    model.eval()

    with torch.no_grad():
        X_t = torch.tensor(X, dtype=torch.float32, device=device)
        # Use dummy treatment - potential outcomes are computed independently of t
        dummy_t = torch.zeros(len(X), device=device)

        if isinstance(model, DragonNet):
            _, (y0, y1), _ = model(X_t, dummy_t)
        else:
            _, (y0, y1) = model(X_t, dummy_t)

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

**Step 3: Run quick test**

```bash
cd benchmarks/external/nn && python -c "
from models import TARNet, predict_ite, train_model
import numpy as np
X = np.random.randn(100, 10)
T = np.random.binomial(1, 0.5, 100)
Y = X[:, 0] + 2 * T + np.random.randn(100)
model = TARNet(10)
model = train_model(model, X, T, Y, n_epochs=10)
result = predict_ite(model, X)
print('ITE shape:', result['ite'].shape)
print('ATE:', result['ate'])
"
```

Expected: Prints ITE shape and ATE without errors

**Step 4: Commit**

```bash
git add benchmarks/external/nn/models.py
git commit -m "perf(benchmark): optimize predict_ite to single forward pass"
```

---

## Phase 4: Important Fixes (Metrics)

### Task 8: Rename mse_ate to squared_error_ate for clarity

**Files:**
- Modify: `packages/cfomics/R/benchmark_metrics.R:43,66,68`
- Modify: `benchmarks/R/runner.R` (multiple lines referencing mse_ate)

**Step 1: Update benchmark_metrics.R**

Change:
```r
mse_ate <- bias_ate^2
```
to:
```r
# Note: This is the squared error for a single estimate, not MSE.
# True MSE requires averaging over replications: MSE = E[(estimate - true)^2]
# Aggregation to MSE happens at the analysis/reporting stage.
squared_error_ate <- bias_ate^2
```

And update the return list:
```r
list(
    bias_ate = as.numeric(bias_ate),
    abs_bias_ate = as.numeric(abs_bias_ate),
    squared_error_ate = as.numeric(squared_error_ate),  # Changed from mse_ate
    ...
)
```

**Step 2: Update runner.R references**

Find and replace all occurrences of `mse_ate` with `squared_error_ate` in:
- `benchmarks/R/runner.R`

```bash
grep -n "mse_ate" benchmarks/R/runner.R
```

Update each occurrence.

**Step 3: Update TCGA runner**

```bash
grep -n "mse_ate" benchmarks/tcga/run_tcga_benchmark.R
```

Update each occurrence.

**Step 4: Run package check**

```bash
Rscript -e "devtools::document('packages/cfomics'); devtools::check('packages/cfomics', args = '--no-tests')"
```

**Step 5: Commit**

```bash
git add packages/cfomics/R/benchmark_metrics.R benchmarks/R/runner.R benchmarks/tcga/run_tcga_benchmark.R
git commit -m "refactor(benchmark): rename mse_ate to squared_error_ate for clarity

The metric computed per-replication is squared error, not MSE.
True MSE = E[(estimate - true)^2] requires averaging across replications,
which happens at the analysis stage."
```

---

### Task 9: Fix dgp_heterogeneous_subgroup ATE computation

**Files:**
- Modify: `packages/cfomics/R/benchmark_dgp.R:405`

**Step 1: Find the heterogeneous subgroup DGP**

```bash
grep -n "dgp_heterogeneous_subgroup" packages/cfomics/R/benchmark_dgp.R
```

**Step 2: Change theoretical ATE to empirical ATE**

Change from:
```r
true_ate = sum(subgroup_props * subgroup_effects),
```

To:
```r
# Use empirical mean of ITEs for consistency with other heterogeneous DGPs
# (theoretical ATE may differ from sample ATE due to finite sample variation)
true_ate = mean(tau),
```

**Step 3: Run tests**

```bash
Rscript -e "devtools::test('packages/cfomics', filter = 'dgp')"
```

**Step 4: Update documentation**

```bash
Rscript -e "devtools::document('packages/cfomics')"
```

**Step 5: Commit**

```bash
git add packages/cfomics/R/benchmark_dgp.R
git commit -m "fix(benchmark): use empirical ATE in dgp_heterogeneous_subgroup

Changed from theoretical weighted average to mean(tau) for consistency
with other heterogeneous DGPs and to match the actual sample."
```

---

## Phase 5: Suggested Improvements (Optional)

### Task 10: Expand tmle3 Super Learner stack (Optional)

**Files:**
- Modify: `benchmarks/external/wrappers/wrapper_tmle3.R:42-44`

**Step 1: Expand learner stack**

Change from:
```r
  # Simple learner stack
  lrnr_glm <- sl3::Lrnr_glm$new()
  lrnr_mean <- sl3::Lrnr_mean$new()
  sl <- sl3::Lrnr_sl$new(learners = list(lrnr_glm, lrnr_mean))
```

To:
```r
  # Super Learner stack with parametric and non-parametric learners
  # Note: For production use, consider adding Lrnr_ranger, Lrnr_xgboost
  lrnr_glm <- sl3::Lrnr_glm$new()
  lrnr_mean <- sl3::Lrnr_mean$new()

  # Add ridge regression for regularization
  lrnr_ridge <- tryCatch(
    sl3::Lrnr_glmnet$new(alpha = 0),
    error = function(e) {
      warning("glmnet not available, using GLM-only stack")
      NULL
    }
  )

  learners <- list(lrnr_glm, lrnr_mean)
  if (!is.null(lrnr_ridge)) {
    learners <- c(learners, list(lrnr_ridge))
  }

  sl <- sl3::Lrnr_sl$new(learners = learners)
```

**Step 2: Run tests**

```bash
Rscript -e "testthat::test_file('benchmarks/external/wrappers/tests/test_tmle3.R')"
```

**Step 3: Commit**

```bash
git add benchmarks/external/wrappers/wrapper_tmle3.R
git commit -m "feat(benchmark): expand tmle3 Super Learner stack with ridge regression"
```

---

## Summary

| Task | Phase | Priority | Description |
|------|-------|----------|-------------|
| 1 | 1 | Critical | Fix MatchIt weights reference |
| 2 | 1 | Critical | Fix WeightIt weights reference |
| 3 | 1 | Critical | Document DragonNet as simplified |
| 4 | 2 | Important | Add PCA bounds checking |
| 5 | 2 | Important | Improve TCGA seed strategy |
| 6 | 3 | Important | Fix CFRNet IPM loss documentation |
| 7 | 3 | Important | Optimize predict_ite |
| 8 | 4 | Important | Rename mse_ate to squared_error_ate |
| 9 | 4 | Important | Fix dgp_heterogeneous_subgroup ATE |
| 10 | 5 | Optional | Expand tmle3 Super Learner |

Total: **10 Tasks** (3 Critical, 6 Important, 1 Optional)
