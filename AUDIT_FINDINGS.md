# cfomics Package Audit Report

**Audit Date:** 2026-01-04
**Package Version:** 0.3.0
**Auditor:** Claude (AI-assisted review)
**Purpose:** Pre-public release quality assessment for CRAN/Bioconductor/JOSS

---

## Executive Summary

This audit identified **42 findings** across 8 categories:
- **P0 (Blocker):** 5 issues - Must fix before public release
- **P1 (High):** 12 issues - Significant problems affecting usability/correctness
- **P2 (Medium):** 17 issues - Technical debt and maintainability concerns
- **P3 (Low):** 8 issues - Minor improvements and style issues

The package has a solid foundation but requires work in three critical areas:
1. **Estimand/uncertainty semantics** - Methods report different estimands without clear documentation
2. **Python integration robustness** - Error handling and environment management need improvement
3. **Bioconductor integration depth** - SE/MAE support is present but underspecified for "omics-focused" claims

---

## Findings by Category

### A. Public Release Requirements (CRAN/Bioconductor/JOSS)

| ID | Severity | Issue | Evidence | Risk | Fix Proposal | Acceptance Criteria |
|---|---|---|---|---|---|---|
| CFOMICS-REV-001 | P0 | **Authors@R has placeholder email/ORCID** | DESCRIPTION:6-9 | Bioconductor/CRAN rejection | Replace with real author contact info and valid ORCID | Valid ORCID format, working email |
| CFOMICS-REV-002 | P1 | **NAMESPACE has duplicate exports** | NAMESPACE:37-43 vs 53-57 | R CMD check warning | Remove duplicate `cf_benchmark_*` exports | No duplicate export warnings |
| CFOMICS-REV-003 | P1 | **All @examples use \dontrun{}** | R/*.R (all exported functions) | CRAN policy violation - examples must run unless truly platform-specific | Convert R-only method examples to runnable, keep \donttest{} for Python methods | At least GRF/IPW/G-formula examples run in R CMD check |
| CFOMICS-REV-004 | P2 | **Missing @return documentation on some functions** | env_diagnostics.R, config.R | R CMD check NOTE | Add explicit @return tags | No missing value documentation notes |
| CFOMICS-REV-005 | P2 | **LazyData: true but no data/ directory** | DESCRIPTION:61 | R CMD check warning | Either remove LazyData or add example datasets | No LazyData warning |
| CFOMICS-REV-006 | P2 | **biocViews present but package targets CRAN-style submission** | DESCRIPTION:22-28 | Confusion about target repository | Clarify submission target; if CRAN, remove biocViews | Clear submission strategy |
| CFOMICS-REV-007 | P3 | **igraph in Imports but only needed for DoWhy-GCM** | DESCRIPTION:38 | Unnecessary dependency for R-only users | Move to Suggests, add runtime check | Reduced install footprint |
| CFOMICS-REV-008 | P3 | **gfoRmula in Suggests but not used** | DESCRIPTION:56 | Package not actually used in gformula method | Remove or implement usage | Clean Suggests list |

### B. Unified API Contract (cf_fit/predict/summary/plot)

| ID | Severity | Issue | Evidence | Risk | Fix Proposal | Acceptance Criteria |
|---|---|---|---|---|---|---|
| CFOMICS-REV-009 | P0 | **IPW method returns constant ITE (=ATE) pretending to be heterogeneous** | methods_traditional.R:164 `ite <- rep(ate, nrow(df))` | User misinterpretation - IPW cannot estimate ITE | Return NA for ITE or add clear warning; document that IPW estimates ATE only | `predict(type="ite")` returns NA or error for IPW with informative message |
| CFOMICS-REV-010 | P1 | **DRLearner/CAVAE don't compute CIs** | methods_drlearner.R:63-65, methods_cavae.R:96-98 `ate_ci_lower = NA` | Inconsistent uncertainty quantification across methods | Implement bootstrap CI or clearly document limitation | Either CI computation or explicit documentation |
| CFOMICS-REV-011 | P1 | **y0_hat/y1_hat computation differs by method** | GRF uses `Y - W*ite`, DRLearner uses `Y - T*ite` | Inconsistent counterfactual definitions | Standardize formula across all methods | Consistent counterfactual computation |
| CFOMICS-REV-012 | P1 | **predict.cf_model dispatches by method with switch** | cf_predict.R:19-28 | Violates open-closed principle; hard to add methods | Consider S3 method dispatch per backend class | Cleaner extensibility |
| CFOMICS-REV-013 | P2 | **summary.cf_model prints to console but returns invisibly** | cf_predict.R:64-99 | Inconsistent with R conventions (should return, optionally print) | Return summary object, use print method for display | `s <- summary(fit)` captures values |
| CFOMICS-REV-014 | P2 | **No S3 plot method in NAMESPACE** | plot.cfomics_result registered but delegates to cf_plot_all | Missing base R plotting for R-only usage | Implement R-native plot fallback when Python unavailable | plot() works without Python |

### C. Causal Inference Semantics (Estimand & Assumptions)

| ID | Severity | Issue | Evidence | Risk | Fix Proposal | Acceptance Criteria |
|---|---|---|---|---|---|---|
| CFOMICS-REV-015 | P0 | **Estimand not stored in meta or documented** | cf_fit.R:163-173 meta structure | Users cannot determine if ATE vs ATT vs ATO | Add `estimand` field to meta; document per-method default | `fit$meta$estimand` returns "ATE" etc. |
| CFOMICS-REV-016 | P0 | **IPW uses stabilized weights but doesn't clarify target population** | methods_traditional.R:151 `svydesign` | Stabilized weights may target ATT not ATE | Document estimand; add option for unstabilized weights | Clear documentation of what IPW estimates |
| CFOMICS-REV-017 | P1 | **No documentation of causal assumptions required per method** | All method docs | Users may misapply methods | Add vignette section on assumptions (exchangeability, positivity, SUTVA) | Assumptions documented per method |
| CFOMICS-REV-018 | P1 | **CI interpretation varies** | GRF uses asymptotic SE, G-formula uses bootstrap | Comparing CIs across methods is misleading | Document CI construction method in summary output | summary() shows CI method |
| CFOMICS-REV-019 | P2 | **No positivity warning integration** | cf_overlap_check exists but not called by cf_fit | Extreme weights can cause estimation problems | Add optional positivity check in cf_fit with warning | Option to check positivity automatically |
| CFOMICS-REV-020 | P2 | **Propensity score diagnostics not saved in model object** | methods_traditional.R (IPW) | Can't post-hoc diagnose propensity model | Store propensity scores in fit$res | Access ps scores from fitted model |

### D. Python Backend Integration (reticulate)

| ID | Severity | Issue | Evidence | Risk | Fix Proposal | Acceptance Criteria |
|---|---|---|---|---|---|---|
| CFOMICS-REV-021 | P0 | **.onLoad calls ensure_python_env() in interactive sessions** | aaa-onload.R:14 | Side effects on package load violate CRAN policy | Remove .onLoad Python init; require explicit setup | Package loads without Python side effects |
| CFOMICS-REV-022 | P1 | **Python error messages not translated to R** | methods_dowhy_gcm.R, methods_ganite.R | Python tracebacks confuse R users | Wrap Python calls in tryCatch with rlang::abort translation | Clear R error messages for Python failures |
| CFOMICS-REV-023 | P1 | **ensure_python_env() has two conflicting implementations** | python_env.R:223-248 and ensure_conda_env() at 258-295 | Inconsistent behavior | Consolidate into single entry point with clear documentation | One clear way to set up Python |
| CFOMICS-REV-024 | P1 | **CAVAE module check requires CUDA** | benchmark_runner.R:207-212 `torch$cuda$is_available()` | CAVAE skipped on CPU-only systems even though it works on CPU | Change to torch availability check, not CUDA | CAVAE works on CPU |
| CFOMICS-REV-025 | P2 | **dowhy version pinned to 0.11.1** | python_env.R:7 | Limits future compatibility | Use compatible version range or document pinning rationale | Documented versioning strategy |
| CFOMICS-REV-026 | P2 | **No test for Python module fallback behavior** | test-*.R files | Can't verify graceful degradation | Add tests that mock Python unavailability | Tests for no-Python scenarios |
| CFOMICS-REV-027 | P3 | **sklearn check but package is scikit-learn** | env_diagnostics.R:160-172 | Minor inconsistency in naming | Consistent naming in diagnostics | Clean diagnostic output |

### E. Bioconductor Data Integration (SE/MAE)

| ID | Severity | Issue | Evidence | Risk | Fix Proposal | Acceptance Criteria |
|---|---|---|---|---|---|---|
| CFOMICS-REV-028 | P1 | **Feature selection seed not saved or controllable** | bioc_integration.R:357-384 select_features_internal | Non-reproducible feature selection | Add seed parameter, store selected features in meta | Reproducible feature selection |
| CFOMICS-REV-029 | P1 | **PCA feature extraction loses biological interpretability** | bioc_integration.R:369-381 | PC loadings not stored; can't map back to genes | Store loadings matrix in output | Feature-to-gene mapping available |
| CFOMICS-REV-030 | P2 | **MAE integration requires manual outcome/treatment specification** | bioc_integration.R:236-239 | Formula interface doesn't work for multi-omics | Improve formula parsing for MAE or document limitation | Clear MAE usage documentation |
| CFOMICS-REV-031 | P2 | **No gene-wise/feature-wise causal analysis mode** | No loop over features | "Omics" use case often involves per-gene effects | Add cf_fit_features() or document workflow | Feature-wise analysis supported or documented |
| CFOMICS-REV-032 | P3 | **Missing imputation options for assay data** | bioc_integration.R uses complete.cases only | Drops rows with any NA | Add imputation strategy parameter | Explicit NA handling |

### F. Benchmark Infrastructure (DGP/metrics/report)

| ID | Severity | Issue | Evidence | Risk | Fix Proposal | Acceptance Criteria |
|---|---|---|---|---|---|---|
| CFOMICS-REV-033 | P1 | **DGPs don't represent omics characteristics** | benchmark_dgp.R all scenarios | Package claims omics focus but DGPs are low-dimensional | Add high-dimensional correlated DGP, latent factor DGP | At least one omics-like DGP |
| CFOMICS-REV-034 | P1 | **No positivity stress DGP** | benchmark_dgp.R - strong_confounding is moderate | Can't test method behavior at positivity boundary | Add near-violation/violation DGPs | DGP with extreme propensity scores |
| CFOMICS-REV-035 | P2 | **DAG created in DGP is fixed to 2 confounders** | benchmark_dgp.R:197-210 .create_benchmark_dag | Doesn't scale with p | Generate DAG matching DGP structure | DAG matches actual data structure |
| CFOMICS-REV-036 | P2 | **Missing policy value / decision quality metrics** | benchmark_metrics.R | PEHE alone insufficient for method selection | Add policy regret, targeting efficiency | Richer metric set |
| CFOMICS-REV-037 | P3 | **No parallel execution option** | benchmark_runner.R uses sequential loop | Slow for large benchmarks | Add parallel option with future/foreach | Optional parallel execution |

### G. Security/Privacy and Safe Defaults

| ID | Severity | Issue | Evidence | Risk | Fix Proposal | Acceptance Criteria |
|---|---|---|---|---|---|---|
| CFOMICS-REV-038 | P2 | **Config stored in ~/.cfomics with no permission check** | config.R:43-48 | May fail on shared systems or containers | Add fallback to tempdir, check write permissions | Graceful fallback for config storage |
| CFOMICS-REV-039 | P2 | **cf_export() creates summary file without user confirmation** | reporting.R:166-168 | May overwrite existing files | Add file existence check with warning | No silent overwrites |
| CFOMICS-REV-040 | P3 | **Python path may appear in logs** | env_diagnostics.R:288 prints full path | Minor privacy concern | Option to mask paths in output | Privacy-respecting output option |

### H. Code Quality and Maintainability

| ID | Severity | Issue | Evidence | Risk | Fix Proposal | Acceptance Criteria |
|---|---|---|---|---|---|---|
| CFOMICS-REV-041 | P2 | **Method dispatch uses switch statement in cf_fit()** | cf_fit.R:117-157 | Adding methods requires modifying core function | Registry pattern for method plugins | Method registration API |
| CFOMICS-REV-042 | P2 | **Bootstrap implementation duplicated** | G-formula has inline bootstrap, not reusable | Inconsistent CI computation | Extract common bootstrap function | Shared bootstrap infrastructure |

---

## Priority Summary for Remediation

### P0 (Must Fix Before Release)
1. **CFOMICS-REV-001:** Real author contact information in DESCRIPTION
2. **CFOMICS-REV-009:** IPW ITE handling (returns misleading values)
3. **CFOMICS-REV-015:** Estimand not stored in meta
4. **CFOMICS-REV-016:** IPW target population unclear
5. **CFOMICS-REV-021:** .onLoad Python side effects

### P1 (High Priority)
- Fix duplicate NAMESPACE exports
- Make examples runnable
- Standardize y0_hat/y1_hat computation
- Document assumptions per method
- Python error message translation
- Feature selection reproducibility
- Omics-representative DGPs

### P2 (Medium Priority)
- Return value consistency
- Plot method for R-only usage
- Propensity diagnostics storage
- Config fallback
- Bootstrap standardization

### P3 (Low Priority)
- Dependency optimization
- Minor naming consistency
- Privacy options

---

## Public Release Pathway Recommendations

### For CRAN Submission
1. Remove biocViews from DESCRIPTION
2. Fix all P0 issues
3. Address examples (P1)
4. Remove .onLoad Python initialization
5. Run R CMD check --as-cran with 0 errors, 0 warnings, minimal notes

### For Bioconductor Submission
1. Keep biocViews, add BiocManager::valid() checks
2. Enhance SE/MAE integration (P1/P2 issues)
3. Add example data in data/
4. Follow BiocCheck requirements
5. Deeper vignette coverage of omics workflows

### For JOSS Publication
1. Focus on correctness (estimand clarity, assumption documentation)
2. Add Statement of Need emphasizing omics application
3. Benchmark against existing packages (grf, DoubleML)
4. Ensure reproducibility (seed management, versioning)

---

## Appendix: Files Reviewed

### R Source Files (23 files)
- R/cf_fit.R - Main API
- R/cf_predict.R - Prediction methods
- R/methods_traditional.R - GRF, IPW, G-formula
- R/methods_dowhy_gcm.R - DoWhy wrapper
- R/methods_ganite.R - GANITE wrapper
- R/methods_drlearner.R - DRLearner wrapper
- R/methods_cavae.R - CAVAE wrapper
- R/bioc_integration.R - SE/MAE support
- R/diagnostics.R - Balance/overlap checks
- R/benchmark_*.R - Benchmark infrastructure
- R/env_*.R - Environment management
- R/config.R - Configuration
- R/visualization.R - Plotting
- R/reporting.R - Export functions

### Tests (20 test files)
- Reviewed test coverage for all major features
- Identified gaps in Python error handling tests
- Good coverage for R-native methods

### Configuration
- DESCRIPTION
- NAMESPACE
- 4 vignettes (not fully reviewed)

---

*End of Audit Report*
