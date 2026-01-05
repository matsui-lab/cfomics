# cfomics: Counterfactual Causal Inference for Omics Data

<!-- badges: start -->
[![R-CMD-check](https://github.com/matsui-lab/counterfactual_agingR/workflows/R-CMD-check/badge.svg)](https://github.com/matsui-lab/counterfactual_agingR/actions)
<!-- badges: end -->

## Overview

`cfomics` provides a unified interface for counterfactual causal inference methods, specifically designed for omics data analysis. It supports Bioconductor data structures (SummarizedExperiment, MultiAssayExperiment) and offers both Python-based deep learning methods and R-native statistical approaches.

**Version**: 0.4.0

## Features

- **Unified API**: Single `cf_fit()` function for all causal inference methods
- **Multiple Backends**: 7 causal inference methods (4 Python-based, 3 R-native)
- **Bioconductor Integration**: Full support for SummarizedExperiment and MultiAssayExperiment
- **Flexible Predictions**: ATE, ITE, counterfactual outcomes
- **Visualization**: Built-in plotting functions for treatment effects
- **Diagnostics**: Covariate balance and overlap checking tools
- **Reporting**: Publication-ready tables and export functionality
- **Benchmarking**: Infrastructure for comparing methods on synthetic data
- **One-Command Setup**: `cf_setup()` for complete environment configuration

## Installation

### From GitHub (development version)

```r
# Install remotes if needed
install.packages("remotes")

# Install cfomics
remotes::install_github("matsui-lab/counterfactual_agingR", subdir = "cfomics")
```

### Python Environment Setup

For Python-based methods (DoWhy, GANITE, DRLearner, CAVAE):

```r
library(cfomics)

# Option 1: One-command setup (recommended)
cf_setup()                           # Standard setup (DoWhy, EconML)
cf_setup(level = "minimal")          # R-only, no Python required
cf_setup(level = "full", gpu = TRUE) # All methods with GPU support

# Option 2: Method-specific setup
cf_setup_method("dowhy")

# Option 3: Manual environment management
cf_install_python_env("unified")
cf_use_python_env("unified")

# Check environment status
cf_check_env()
```

## Quick Start

### Basic Usage with data.frame

```r
library(cfomics)
library(igraph)

# Create sample data
set.seed(123)
n <- 200
data <- data.frame(
  X1 = rnorm(n),
  X2 = rnorm(n),
  T = rbinom(n, 1, 0.5),
  Y = rnorm(n)
)
data$Y <- data$Y + 2 * data$T + 0.5 * data$X1  # True ATE = 2

# Fit model using GRF (R-native, no Python required)
fit <- cf_fit(
  Y ~ T | X1 + X2,
  data = data,
  method = "grf"
)

# Get results
summary(fit)
predict(fit, type = "ate")  # Average Treatment Effect
predict(fit, type = "ite")  # Individual Treatment Effects

# Visualize
plot(fit)
```

### Using Automatic Method Selection

```r
# Let cfomics choose the best available method
fit <- cf_fit(
  Y ~ T | X1 + X2,
  data = data,
  method = "auto"  # Selects best available method
)
```

### Using DoWhy-GCM (Python-based)

```r
# Setup Python environment (one-time)
cf_setup_method("dowhy")

# Define causal DAG
g <- igraph::graph_from_edgelist(rbind(
  c("X1", "T"), c("X1", "Y"),
  c("X2", "Y"), c("T", "Y")
))

# Fit model
fit <- cf_fit(
  Y ~ T | X1 + X2,
  data = data,
  method = "dowhy_gcm",
  graph = g,
  bootstrap = TRUE
)

# Get results with confidence intervals
summary(fit)
```

## Available Methods

| Method | Backend | Python Required | Best For |
|--------|---------|-----------------|----------|
| `grf` | R (grf) | No | Heterogeneous treatment effects |
| `ipw` | R (ipw) | No | Inverse probability weighting |
| `gformula` | R | No | G-computation |
| `dowhy_gcm` | Python (DoWhy) | Yes | DAG-based causal inference |
| `drlearner` | Python (EconML) | Yes | Double/debiased ML |
| `ganite` | Python (TensorFlow) | Yes | Deep learning ITE |
| `cavae` | Python (Pyro) | Yes | Latent confounding |

## Formula Syntax

cfomics uses a special formula syntax:

```
Y ~ T | X1 + X2 + X3
```

- `Y`: Outcome variable (numeric)
- `T`: Treatment variable (binary: 0/1)
- `X1, X2, X3`: Covariates (numeric)

## Prediction Types

```r
# Average Treatment Effect
ate <- predict(fit, type = "ate")

# Individual Treatment Effects
ite <- predict(fit, type = "ite")

# Counterfactual outcomes
y0 <- predict(fit, type = "y0")  # Outcome if untreated
y1 <- predict(fit, type = "y1")  # Outcome if treated

# Summary statistics with CIs
summary_stats <- predict(fit, type = "summary")
```

## Visualization

```r
# ITE distribution
cf_plot_ite(predict(fit, type = "ite"))

# Outcome comparison
cf_plot_outcome_shift(
  predict(fit, type = "y0"),
  predict(fit, type = "y1")
)

# Uplift curve
cf_plot_uplift(predict(fit, type = "ite"))

# All plots at once
cf_plot_all(fit, output_dir = "plots/")
```

## Benchmarking

```r
# Run benchmark on synthetic data
results <- cf_benchmark_run(
  scenarios = c("linear_homogeneous", "heterogeneous_ite"),
  methods = c("grf", "ipw", "gformula"),
  n = 500,
  p = 10,
  n_reps = 20
)

# Summarize results
cf_benchmark_summarize(results)
```

## Diagnostics

cfomics provides tools for checking causal inference assumptions:

```r
# Check covariate balance between treatment groups
balance <- cf_balance_check(data, treatment = "T", covariates = c("X1", "X2"))
print(balance)
plot(balance)

# Check overlap/positivity assumption
overlap <- cf_overlap_check(data, treatment = "T", covariates = c("X1", "X2"))
print(overlap)
plot(overlap)
```

## Reporting

Generate publication-ready results:
```r
# Create formatted results table
cf_table(fit)                         # data.frame
cf_table(fit, format = "markdown")    # Markdown
cf_table(fit, format = "latex")       # LaTeX

# Export results to file
cf_export(fit, "results.csv")         # CSV format
cf_export(fit, "results.rds")         # RDS format
cf_export(fit, "results.json")        # JSON format

# Session information for reproducibility
cf_session_info(fit)
```

## Bioconductor Integration

Use cfomics with SummarizedExperiment and MultiAssayExperiment:

```r
library(SummarizedExperiment)

# Convert Bioconductor object to cfomics format
cf_data <- as_cf_data(
  se,                              # SummarizedExperiment object
  outcome = "gene_expression",     # Outcome variable
  treatment = "treatment",         # Treatment variable
  feature_select = "variance",     # Feature selection: "variance", "pca", "none"
  n_features = 50                  # Number of features to include
)

# Fit model
fit <- cf_fit(outcome ~ treatment | ., data = cf_data, method = "grf")
```

## Configuration

Persistent configuration management:

```r
# View all settings
cf_config()

# Get/set specific values
cf_config("methods.default")           # Get
cf_config("methods.default", "grf")    # Set

# Reset to defaults
cf_config_reset()
```

## Environment Management

```r
# Comprehensive environment check
cf_check_env()

# One-command setup
cf_setup()                           # Standard
cf_setup(level = "full")             # All methods

# Install specific environment
cf_install_python_env("unified")     # Unified environment
cf_install_python_env("dowhy")       # Method-specific

# Activate environment
cf_use_python_env("unified")

# Reset all cfomics environments
cf_reset_env()
```

## Requirements

### R Dependencies

**Required:**
- R (>= 4.2.0)
- Formula, igraph, reticulate, rlang, cli, glue

**Suggested (for specific methods):**
- grf (for GRF method)
- ipw, survey (for IPW method)
- ggplot2, gridExtra (for enhanced visualizations)
- SummarizedExperiment, MultiAssayExperiment (for Bioconductor integration)
- yaml, jsonlite (for configuration and export)
- knitr, rmarkdown (for vignettes)

Note: The G-formula method uses base R (`lm()`) and does not require additional packages.

### Python Dependencies

For Python-based methods, you need:

- Python 3.9+ (3.11 recommended)
- numpy, pandas, scipy, scikit-learn (base dependencies)

**Method-specific packages:**
- dowhy (for DoWhy-GCM)
- econml (for DRLearner)
- tensorflow or tensorflow-cpu (for GANITE)
- torch, pyro-ppl (for CAVAE)

All Python dependencies can be installed automatically using `cf_setup()`.

## Citation

```bibtex
@software{cfomics,
  author = {Yusuke Matsui},
  title = {cfomics: Counterfactual Causal Inference for Omics Data},
  version = {0.4.0},
  year = {2025},
  url = {https://github.com/matsui-lab/counterfactual_agingR}
}
```

## Development

### Running Package Checks

```bash
# Quick check (no vignettes, faster)
Rscript tools/check_package.R --quick

# Full CRAN-style check
Rscript tools/check_package.R --full

# With BiocCheck
Rscript tools/check_package.R --bioc
```

Or manually:

```bash
R CMD build .
R CMD check --as-cran cfomics_*.tar.gz
```

### Running Tests

```r
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-traditional.R")
```

### Building Vignettes

```r
devtools::build_vignettes()
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Related Projects

- [DoWhy](https://github.com/py-why/dowhy): Python library for causal inference
- [EconML](https://github.com/microsoft/EconML): Machine learning for causal effects
- [grf](https://github.com/grf-labs/grf): Generalized Random Forests
