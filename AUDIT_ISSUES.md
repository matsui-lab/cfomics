# cfomics Audit - GitHub Issue Templates

This document contains ready-to-use GitHub issue templates for all P0 and P1 findings.

---

## P0 Issues (Blockers)

### Issue #1: Replace placeholder author information in DESCRIPTION
**Labels:** `P0-blocker`, `documentation`

**Title:** Replace placeholder author information in DESCRIPTION

**Body:**
```markdown
## Problem
DESCRIPTION contains placeholder author information:
- Email: `matsui-lab@example.com`
- ORCID: `0000-0000-0000-0000`

This will cause rejection from CRAN and Bioconductor.

## Location
`cfomics/DESCRIPTION` lines 6-9

## Fix Required
Replace with actual author contact information and valid ORCID identifiers.

## Acceptance Criteria
- [ ] Valid email addresses that will receive package-related communication
- [ ] Valid ORCID format (https://orcid.org/)
- [ ] R CMD check passes without author-related notes

## Related Finding
CFOMICS-REV-001
```

---

### Issue #2: IPW method returns misleading ITE values
**Labels:** `P0-blocker`, `bug`, `causal-semantics`

**Title:** IPW method returns constant ITE pretending to be heterogeneous

**Body:**
```markdown
## Problem
The IPW method returns `ite <- rep(ate, nrow(df))` (line 164 in methods_traditional.R), making it appear that individual treatment effects are estimated when IPW fundamentally cannot estimate ITEs.

This creates a serious risk of user misinterpretation.

## Location
`cfomics/R/methods_traditional.R:164`

## Current Behavior
```r
fit <- cf_fit(..., method = "ipw")
predict(fit, type = "ite")  # Returns vector of identical values = ATE
```

## Expected Behavior
Options:
1. Return NA with informative message explaining IPW estimates population ATE only
2. Return error with guidance
3. At minimum, clear warning that ITE values are imputed from ATE

## Acceptance Criteria
- [ ] `predict(fit, type = "ite")` for IPW does NOT return fake ITE values
- [ ] User receives clear message about IPW's limitations
- [ ] Documentation updated to explain estimand capabilities per method

## Related Finding
CFOMICS-REV-009
```

---

### Issue #3: Estimand not stored in model metadata
**Labels:** `P0-blocker`, `enhancement`, `causal-semantics`

**Title:** Store estimand (ATE/ATT/ATO) in model metadata

**Body:**
```markdown
## Problem
The `cf_model` object does not store what estimand was estimated. Users cannot determine if they got ATE, ATT, or ATO without reading source code.

Different methods target different estimands:
- GRF: CATE (individual level), averaged to ATE
- IPW with stabilized weights: May target ATT
- G-formula: ATE (marginal)

## Location
`cfomics/R/cf_fit.R:163-173` (meta structure definition)

## Current Structure
```r
meta = list(
  formula, n, p, outcome_name, treatment_name,
  covariate_names, dag_spec, bootstrap, random_state
  # estimand is missing!
)
```

## Required Changes
1. Add `estimand` field to meta
2. Each method must declare its estimand
3. Document estimand in summary() output

## Acceptance Criteria
- [ ] `fit$meta$estimand` returns character (e.g., "ATE", "ATT", "CATE")
- [ ] summary(fit) displays estimand
- [ ] Documentation explains estimand per method

## Related Finding
CFOMICS-REV-015
```

---

### Issue #4: IPW target population unclear
**Labels:** `P0-blocker`, `documentation`, `causal-semantics`

**Title:** Document IPW target population (stabilized weights affect estimand)

**Body:**
```markdown
## Problem
IPW implementation uses stabilized weights via `ipwpoint()` but does not document that this affects the target population. Stabilized weights can shift estimation from ATE to something closer to ATT.

## Location
`cfomics/R/methods_traditional.R:137-152`

## Required Changes
1. Document that stabilized weights are used
2. Clarify what estimand results (ATE vs weighted average)
3. Consider adding `stabilized = TRUE/FALSE` option with different estimand labels
4. Update summary() to show target population

## Acceptance Criteria
- [ ] Documentation explicitly states IPW target population
- [ ] `fit$meta$estimand` reflects actual estimand
- [ ] Vignette includes IPW-specific guidance

## Related Finding
CFOMICS-REV-016
```

---

### Issue #5: Remove Python initialization from .onLoad
**Labels:** `P0-blocker`, `bug`, `CRAN-compliance`

**Title:** Remove Python side effects from .onLoad

**Body:**
```markdown
## Problem
`.onLoad()` calls `ensure_python_env()` for interactive sessions, which violates CRAN policy against package load side effects.

```r
if (interactive()) {
  try(ensure_python_env(), silent = TRUE)
}
```

This can:
- Trigger Python initialization unexpectedly
- Cause slow package loads
- Fail in restricted environments

## Location
`cfomics/R/aaa-onload.R:14`

## Required Changes
1. Remove `ensure_python_env()` call from `.onLoad()`
2. Keep only S3 method registration
3. Document that users must explicitly set up Python

## Acceptance Criteria
- [ ] `library(cfomics)` completes without Python initialization
- [ ] No network access on package load
- [ ] R CMD check passes without side-effect notes

## Related Finding
CFOMICS-REV-021
```

---

## P1 Issues (High Priority)

### Issue #6: Remove duplicate NAMESPACE exports
**Labels:** `P1-high`, `bug`

**Title:** Remove duplicate exports in NAMESPACE

**Body:**
```markdown
## Problem
`cf_benchmark_run`, `cf_benchmark_summarize`, etc. are exported twice in NAMESPACE.

## Location
`cfomics/NAMESPACE` lines 37-43 and 53-57

## Fix
Remove duplicate block (lines 53-57).

## Acceptance Criteria
- [ ] No duplicate exports
- [ ] R CMD check shows no export warnings

## Related Finding
CFOMICS-REV-002
```

---

### Issue #7: Convert examples from \dontrun to runnable
**Labels:** `P1-high`, `documentation`, `CRAN-compliance`

**Title:** Make R-only method examples runnable

**Body:**
```markdown
## Problem
All @examples use `\dontrun{}`, which CRAN discourages. R-only methods (grf, ipw, gformula) should have runnable examples.

## Required Changes
1. GRF, IPW, G-formula examples: Use `\donttest{}` (runs locally, not on CRAN)
2. Python methods: Keep `\dontrun{}` with explanatory note
3. Ensure examples complete in < 5 seconds

## Acceptance Criteria
- [ ] R CMD check --as-cran runs examples without error
- [ ] At least cf_fit() with method="grf" has runnable example

## Related Finding
CFOMICS-REV-003
```

---

### Issue #8: Standardize y0_hat/y1_hat computation across methods
**Labels:** `P1-high`, `bug`, `consistency`

**Title:** Standardize counterfactual outcome computation

**Body:**
```markdown
## Problem
Different methods compute y0_hat/y1_hat inconsistently:
- GRF: `y0_hat <- Y_vec - W * ite; y1_hat <- Y_vec + (1 - W) * ite`
- DRLearner: `y0_hat <- Y - T * ite; y1_hat <- Y + (1 - T) * ite`

Variable naming differs (W vs T) but more importantly, the semantics should be identical.

## Required Changes
1. Standardize to consistent formula: `y0 = Y - T*tau; y1 = Y + (1-T)*tau`
2. Or compute as: `y0 = Y*(1-T) + (Y - tau)*T; y1 = Y*T + (Y + tau)*(1-T)`
3. Document the counterfactual definition used

## Acceptance Criteria
- [ ] All methods use identical counterfactual formula
- [ ] Documentation explains counterfactual computation

## Related Finding
CFOMICS-REV-011
```

---

### Issue #9: DRLearner and CAVAE need confidence intervals
**Labels:** `P1-high`, `enhancement`

**Title:** Implement CI computation for DRLearner and CAVAE

**Body:**
```markdown
## Problem
DRLearner and CAVAE return `ate_ci_lower = NA, ate_ci_upper = NA`, making them less useful for inference.

## Options
1. Bootstrap CI (consistent with G-formula approach)
2. For DRLearner: Use EconML's built-in inference methods
3. For CAVAE: Posterior credible intervals from Pyro

## Acceptance Criteria
- [ ] CIs available for DRLearner ATE
- [ ] CIs available for CAVAE ATE
- [ ] CI method documented (bootstrap vs asymptotic vs Bayesian)

## Related Finding
CFOMICS-REV-010
```

---

### Issue #10: Document causal assumptions per method
**Labels:** `P1-high`, `documentation`

**Title:** Add assumption documentation for each method

**Body:**
```markdown
## Problem
Users may misapply methods without understanding required assumptions:
- Exchangeability (no unmeasured confounding)
- Positivity (all units have positive probability of each treatment)
- SUTVA (no interference, consistent treatment)

## Required Changes
1. Add "Assumptions" section to each method's documentation
2. Create vignette section on assumption checking
3. Link cf_balance_check() and cf_overlap_check() to relevant methods

## Acceptance Criteria
- [ ] Each method's help page lists required assumptions
- [ ] Vignette covers assumption diagnostics
- [ ] cf_fit() help mentions where to find assumption info

## Related Finding
CFOMICS-REV-017
```

---

### Issue #11: Python error translation
**Labels:** `P1-high`, `enhancement`, `usability`

**Title:** Translate Python errors to user-friendly R messages

**Body:**
```markdown
## Problem
Python errors bubble up as raw tracebacks, confusing R users.

## Required Changes
1. Wrap Python calls in tryCatch
2. Extract relevant error message
3. Use rlang::abort() with:
   - Clear problem description
   - Suggested fix (e.g., "Run cf_check_env() to diagnose")
   - Original Python error in `parent` for debugging

## Example
```r
tryCatch({
  # Python call
}, error = function(e) {
  rlang::abort(
    "DoWhy model fitting failed",
    body = c("Check that all variables are in the DAG",
             "Run cf_check_env() to verify Python setup"),
    parent = e
  )
})
```

## Acceptance Criteria
- [ ] Python errors show clear R message
- [ ] Original error preserved for debugging
- [ ] Actionable suggestions provided

## Related Finding
CFOMICS-REV-022
```

---

### Issue #12: CAVAE should work on CPU
**Labels:** `P1-high`, `bug`

**Title:** CAVAE method incorrectly requires CUDA

**Body:**
```markdown
## Problem
Benchmark runner skips CAVAE unless CUDA is available, but PyTorch/Pyro work on CPU.

## Location
`cfomics/R/benchmark_runner.R:207-212`
```r
torch_ok <- tryCatch({
  torch <- reticulate::import("torch", convert = FALSE)
  torch$cuda$is_available()  # Should not require CUDA
}, error = function(e) FALSE)
```

## Fix
```r
torch_ok <- tryCatch({
  reticulate::import("torch", convert = FALSE)
  TRUE
}, error = function(e) FALSE)
```

## Acceptance Criteria
- [ ] CAVAE runs on CPU-only systems
- [ ] Tests pass without GPU

## Related Finding
CFOMICS-REV-024
```

---

### Issue #13: Feature selection reproducibility
**Labels:** `P1-high`, `bug`, `reproducibility`

**Title:** Make SE/MAE feature selection reproducible with seed

**Body:**
```markdown
## Problem
`select_features_internal()` uses PCA or variance ranking without seed control. Results are not reproducible.

## Location
`cfomics/R/bioc_integration.R:357-384`

## Required Changes
1. Add seed parameter to feature selection
2. Store selected feature names in meta
3. For PCA: store loadings matrix

## Acceptance Criteria
- [ ] Same seed produces same feature selection
- [ ] Selected features accessible from fit object
- [ ] PCA loadings available for interpretation

## Related Finding
CFOMICS-REV-028
```

---

### Issue #14: Add omics-representative DGP scenarios
**Labels:** `P1-high`, `enhancement`, `benchmarking`

**Title:** Add high-dimensional correlated DGP for omics benchmarking

**Body:**
```markdown
## Problem
Current DGPs are low-dimensional (p=20) and don't represent omics characteristics:
- High-dimensional (p >> n)
- Strong feature correlations (gene pathways)
- Latent confounders (batch effects)

## Required Changes
Add new DGP scenarios:
1. `high_dimensional`: p=500, n=200, sparse effects
2. `correlated_features`: Block correlation structure
3. `latent_confounder`: Unobserved batch effect

## Acceptance Criteria
- [ ] At least one DGP with p > n
- [ ] At least one DGP with correlated features
- [ ] Documentation explains omics relevance

## Related Finding
CFOMICS-REV-033
```

---

## Summary Checklist

### Before CRAN Submission
- [ ] Fix all 5 P0 issues
- [ ] Fix P1 issues #6, #7 (NAMESPACE, examples)
- [ ] R CMD check --as-cran: 0 errors, 0 warnings
- [ ] Remove biocViews from DESCRIPTION

### Before Bioconductor Submission
- [ ] Fix all P0 and P1 issues
- [ ] Add example data
- [ ] BiocCheck passes
- [ ] Deep vignette coverage

### Before JOSS Submission
- [ ] All P0 issues fixed
- [ ] Causal semantics clearly documented
- [ ] Reproducibility ensured (seeds, versions)
- [ ] Comparative benchmarks included
