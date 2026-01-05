# cfomics 0.4.0 (Development)

## Package Quality Improvements

This release focuses on package maturity, stability, and distribution readiness.

### Breaking Changes

* Removed automatic Python initialization from `.onLoad`
  - Python methods now require explicit environment setup before use
  - Use `cf_check_env()` to verify configuration
  - This improves CRAN/Bioconductor compliance and HPC compatibility

* Removed `ensure_python_env()` and `ensure_conda_env()` from exports
  - Use `cf_has_python()` and `cf_require_python()` instead

### New Features

* **Lazy Python initialization**: New `cf_has_python()` and `cf_require_python()`
  - Check Python availability without initializing: `cf_has_python("dowhy_gcm")`
  - Require Python with informative errors: `cf_require_python("dowhy_gcm")`
  - Consistent error messages with setup instructions

* **Result Contract specification**: Formal documentation of `cfomics_result` structure
  - Standardized `meta` object with software versions and timestamps
  - Validation function `validate_cfomics_result()`
  - Helper function `create_cf_meta()` for method implementers

* **Enhanced documentation**:
  - Estimand definitions (ATE, ATT, ATC, ITE, CATE) in `cf_fit()` help
  - Causal assumptions (SUTVA, Ignorability, Positivity) documented
  - Method selection guide in documentation

* **Example data** (requires R package build):
  - `cfomics_example_df`: Simple data.frame for basic examples
  - `cfomics_example_se`: SummarizedExperiment for omics examples
  - Data generation script in `data-raw/`

* **GitHub Actions CI**:
  - R CMD check on Ubuntu, macOS, Windows
  - R release and development versions
  - Optional Python integration tests

### Improvements

* Vignettes now build without Python
  - R-only examples execute during build
  - Python examples use conditional evaluation
  - Graceful fallback messages when Python unavailable

* Removed unused `gfoRmula` package from Suggests
  - G-formula method uses base R's `lm()`, no external package needed

* Fixed NAMESPACE duplicate exports for benchmark functions

* Enhanced `cf_check_env()` diagnostics output

* All Python methods now use unified `cf_require_python()` gate

### Documentation

* New `cfomics_result` help page documenting result structure
* Updated README with clearer R-only vs Python sections
* Improved vignette structure with working R-only examples

---

# cfomics 0.3.0

## Major Changes

### Bioconductor Data Structure Integration

* **SummarizedExperiment support**: Use gene expression data directly
  - Extract outcomes from colData or assay rownames (gene names)
  - Automatic covariate extraction from colData
  - Feature selection from assay: "variance" or "pca"

* **MultiAssayExperiment support**: Multi-omics integration
  - Combine features from multiple assays
  - Common sample matching across experiments
  - Configurable features per assay

* **New generic function**: `as_cf_data()` for flexible data conversion
  - Methods for data.frame, SummarizedExperiment, MultiAssayExperiment

### Environment Setup Improvements

* **One-command setup**: `cf_setup()` for complete environment configuration
  - Three levels: "minimal" (R-only), "standard", "full" (deep learning)
  - GPU support option
  - Automatic Miniconda installation

* **Environment diagnostics**: `cf_check_env()` comprehensive checks
  - R version, Python availability, conda status
  - Method availability verification
  - System resource reporting
  - Actionable recommendations

* **Configuration management**: `cf_config()` persistent settings
  - YAML-based configuration (~/.cfomics/config.yaml)
  - Dot notation for nested keys
  - `cf_config_reset()` to restore defaults

* **Unified Python environment**: New "unified" and "unified_full" options
  - Single environment for most methods
  - Reduced complexity from 4 separate environments

### Research Workflow Features

* **Covariate balance check**: `cf_balance_check()`
  - Standardized mean difference (SMD) calculation
  - Configurable imbalance threshold
  - print() and plot() methods

* **Overlap/positivity check**: `cf_overlap_check()`
  - Propensity score distribution analysis
  - Extreme propensity score detection
  - Overlap region visualization

* **Publication-ready reporting**: `cf_table()`
  - Multiple formats: data.frame, markdown, latex, html
  - Configurable decimal places
  - ATE, CI, and ITE statistics

* **Results export**: `cf_export()`
  - CSV, RDS, JSON formats
  - Automatic summary file generation

* **Session info**: `cf_session_info()` for reproducibility

### API Improvements

* **Automatic method selection**: `method = "auto"` in `cf_fit()`
  - Selects best available method based on installed packages
  - Prefers R-native methods (grf) for reliability

* **Extended cf_fit()** parameters
  - `assay_name`: Specify which assay to use
  - `feature_select`: Feature selection method
  - `n_features`: Number of features to select

### Bioconductor Compliance

* Explicit NAMESPACE exports (removed `exportPattern`)
* Added `biocViews` for proper categorization
* Added proper `Authors@R`, `URL`, and `BugReports` fields
* Created LICENSE file

## Bug Fixes

* Fixed duplicate function definition for `summary.cfomics_result`
* Removed hardcoded Python paths from tests
* Added portable Python environment detection
* Fixed missing exports for benchmark API functions in manual NAMESPACE management

## Documentation

* Comprehensive README.md with examples
* NEWS.md for version tracking
* Improved roxygen2 documentation

---

# cfomics 0.2.0

## New Features

* Added benchmark infrastructure for comparing causal inference methods
  - `cf_benchmark_run()`: Run comprehensive benchmarks
  - `cf_benchmark_generate_data()`: Generate synthetic data with known ground truth
  - `cf_benchmark_compute_metrics()`: Compute ATE bias, ITE RMSE, and PEHE

* Added four data generating processes (DGPs):
  - `linear_homogeneous`: Constant treatment effect
  - `nonlinear_outcome`: Nonlinear outcome functions
  - `heterogeneous_ite`: Effect varies with covariates
  - `strong_confounding`: Challenging propensity score scenarios

## Improvements

* Enhanced test coverage for all methods
* Improved error messages for Python environment issues

---

# cfomics 0.1.0

## Initial Release

* Unified API for causal inference with `cf_fit()` and `cf_dowhy()`
* Support for 7 causal inference methods:
  - **Python-based**: DoWhy-GCM, GANITE, DRLearner, CAVAE
  - **R-native**: GRF, IPW, G-formula
* Prediction interface with `predict()` for ATE, ITE, counterfactuals
* Visualization functions:
  - `cf_plot_ite()`: ITE distribution
  - `cf_plot_outcome_shift()`: Counterfactual comparison
  - `cf_plot_uplift()`: Uplift curves
  - `cf_plot_ate_bootstrap()`: Bootstrap confidence intervals
* Python environment management:
  - `cf_install_python_env()`: Create conda environments
  - `cf_use_python_env()`: Activate environments
  - `setup_python_env()`: One-step setup
