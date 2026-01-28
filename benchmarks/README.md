# cfomics Benchmark Suite

Comprehensive benchmark comparing 5 causal inference methods across 47 scenarios from 11 distinct data generating processes (DGPs). This benchmark suite evaluates method performance under various challenging conditions common in high-dimensional omics data analysis.

## Methods

| Method | Description | Type |
|--------|-------------|------|
| `gformula` | G-computation with regularized regression | Parametric |
| `hdml` | High-Dimensional Machine Learning (Double ML) | Semiparametric |
| `hdps` | High-Dimensional Propensity Score | Propensity-based |
| `bcf` | Bayesian Causal Forests | Nonparametric |
| `tmle` | Targeted Maximum Likelihood Estimation | Semiparametric |

## Scenarios

The benchmark covers 11 distinct data generating processes (DGPs):

| ID | Scenario | Description | Variations |
|----|----------|-------------|------------|
| S1 | Baseline | Standard confounding setup | n={500,1000,2000}, p={50,100} |
| S2 | Dimension Sweep | High-dimensional settings | n={500,1000}, p={100,500,1000} |
| S3 | Heterogeneous Linear | Linear treatment effect heterogeneity | strength={0,0.5,1,2} |
| S4 | Heterogeneous (Other) | Nonlinear, subgroup, qualitative heterogeneity | 3 types |
| S5 | Nonlinear Confounding | Complex confounding relationships | 5 nonlinearity types |
| S6 | Dense Confounding | Many true confounders | n_confounders={10,50,100,200} |
| S7 | Weak Overlap | Propensity score overlap violations | 4 severity levels |
| S8 | Covariate Shift | Distribution shift between groups | mean, variance shifts |
| S9 | Correlated Confounding | Correlated covariate structures | block, AR1, factor |
| S10 | Unobserved Confounding | Unmeasured confounders | strength={0,0.5,1,2} |
| S11 | Collider Bias | Collider conditioning | strength={0,0.5,1,2} |

## Directory Structure

```
benchmarks/
├── README.md                 # This file
├── config.R                  # Centralized benchmark configuration
├── run_full_benchmark.R      # Main entry point for running benchmarks
├── aggregate_results.R       # Collect and summarize raw results
├── .gitignore                # Ignores results/ and generated files
├── R/
│   ├── scenarios.R           # DGP dispatch functions
│   ├── runner.R              # Parallel execution engine
│   ├── metrics.R             # Statistical tests (Friedman, Nemenyi)
│   └── reporting.R           # Visualization functions
├── report/
│   ├── benchmark_report.Rmd  # Main report template
│   └── sections/             # Language-specific sections (EN/JA)
│       ├── 01_introduction_en.Rmd
│       ├── 01_introduction_ja.Rmd
│       ├── 02_methods_en.Rmd
│       ├── 02_methods_ja.Rmd
│       ├── 03_simulation_design_en.Rmd
│       ├── 03_simulation_design_ja.Rmd
│       ├── 04_results_en.Rmd
│       ├── 04_results_ja.Rmd
│       ├── 05_discussion_en.Rmd
│       ├── 05_discussion_ja.Rmd
│       ├── 06_appendix_en.Rmd
│       └── 06_appendix_ja.Rmd
└── results/                  # Output directory (gitignored)
    ├── raw/                  # Individual job RDS files
    ├── summary/              # Aggregated CSV files
    ├── plots/                # Generated figures
    └── tables/               # Generated tables
```

## Prerequisites

### R Version
- R >= 4.2.0

### Required Packages

**cfomics** (the main package):
```r
# Install cfomics from source
devtools::install("packages/cfomics")
```

**Benchmark Infrastructure**:
```r
install.packages(c("future", "furrr", "progressr"))
```

**Statistical Analysis**:
```r
install.packages("PMCMRplus")  # For Nemenyi post-hoc tests
```

**Visualization and Reporting**:
```r
install.packages(c("ggplot2", "scales", "patchwork", "kableExtra", "rmarkdown"))
```

### Optional: Install All Dependencies at Once

```r
# All required packages
pkgs <- c(
  "future", "furrr", "progressr",      # Parallel execution
  "PMCMRplus",                          # Statistical tests
  "ggplot2", "scales", "patchwork",     # Visualization
  "kableExtra", "rmarkdown"             # Report generation
)
install.packages(pkgs)
```

## Quick Start

### 1. Smoke Test (~5-10 minutes)

Run a quick validation with 3 scenarios and 2 replications:

```bash
# From repository root
Rscript benchmarks/run_full_benchmark.R --smoke
```

This runs:
- 3 representative scenarios (S1_n500_p50, S3_str0, S7_good)
- 5 methods each
- 2 replications per scenario-method combination
- Reduced BCF MCMC iterations (100 burn-in, 200 iterations)

### 2. Full Benchmark (~3-6 hours)

Run the complete benchmark suite:

```bash
Rscript benchmarks/run_full_benchmark.R
```

This runs:
- 47 scenarios
- 5 methods each
- 50 replications per scenario-method combination
- Total: 11,750 individual jobs

### 3. Aggregate Results

After running the benchmark, aggregate the raw results:

```bash
Rscript benchmarks/aggregate_results.R
```

This produces:
- `results/summary/raw_results.csv` - All individual results combined
- `results/summary/aggregated_summary.csv` - Metrics aggregated by scenario-method

### 4. Generate Reports

**English Report:**
```bash
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="en"))'
```

**Japanese Report:**
```bash
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="ja"))'
```

**Specify Output Directory:**
```bash
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="en"), output_dir="results/reports")'
```

## Command-Line Options

### run_full_benchmark.R

```
Usage: Rscript run_full_benchmark.R [options]

Options:
  --smoke   Run quick smoke test (3 scenarios, 2 reps)
  --retry   Retry previously failed jobs
  --help    Show help message and exit

Examples:
  Rscript benchmarks/run_full_benchmark.R           # Full production run
  Rscript benchmarks/run_full_benchmark.R --smoke   # Quick smoke test
  Rscript benchmarks/run_full_benchmark.R --retry   # Retry failed jobs
  Rscript benchmarks/run_full_benchmark.R --smoke --retry
```

## Expected Outputs

### Raw Results (`results/raw/`)

Individual RDS files for each job:
```
S1_n500_p50__gformula__rep001.rds
S1_n500_p50__gformula__rep002.rds
...
S11_str2__tmle__rep050.rds
```

Each file contains a data frame with:
- `scenario_id`, `method`, `rep`, `n`, `p`, `seed`
- `ate_true`, `ate_hat`, `bias_ate`, `abs_bias_ate`, `mse_ate`
- `pehe` (Precision in Estimating Heterogeneous Effects)
- `coverage_ate`, `ci_len_ate`
- `time_sec`, `status`, `error_msg`

### Aggregated Summary (`results/summary/`)

**raw_results.csv**: All individual results combined

**aggregated_summary.csv**: Per scenario-method statistics:
- `scenario_id`, `method`, `n_reps`, `n`, `p`, `ate_true`
- `mean_bias`, `sd_bias`, `mean_abs_bias`
- `rmse_ate` (Root Mean Squared Error)
- `mean_pehe`, `sd_pehe`
- `coverage` (proportion of CIs containing true ATE)
- `mean_ci_length`, `sd_ci_length`
- `mean_time`, `sd_time`, `total_time`

### Reports (`report/`)

PDF and HTML reports including:
- Executive summary with key findings
- Method descriptions and notation
- Simulation design details
- Comprehensive results with figures and tables
- Statistical tests (Friedman test, Nemenyi post-hoc)
- Critical Difference diagrams
- Discussion and recommendations

## Metrics Computed

| Metric | Description | Formula |
|--------|-------------|---------|
| Bias | Average difference from true ATE | E[ATE_hat - ATE_true] |
| RMSE | Root mean squared error | sqrt(E[(ATE_hat - ATE_true)^2]) |
| PEHE | Precision in Estimating Heterogeneous Effects | sqrt(E[(ITE_hat - ITE_true)^2]) |
| Coverage | Proportion of 95% CIs containing true ATE | mean(CI_lower <= ATE_true <= CI_upper) |
| CI Length | Average confidence interval width | E[CI_upper - CI_lower] |

## Statistical Analysis

The benchmark uses nonparametric statistical tests for method comparison:

1. **Friedman Test**: Omnibus test for differences across methods
2. **Nemenyi Post-hoc Test**: Pairwise comparisons with family-wise error control
3. **Critical Difference (CD) Diagrams**: Visual representation of statistical differences

## Configuration

The benchmark configuration is centralized in `config.R`:

```r
benchmark_config <- function() {
  list(
    methods = c("gformula", "hdml", "hdps", "bcf", "tmle"),
    n_reps = 50L,              # Replications per scenario
    base_seed = 20260128L,     # Random seed base
    n_workers = 0L,            # 0 = auto-detect
    bcf_n_burn = 1000L,        # BCF burn-in
    bcf_n_iter = 2000L,        # BCF iterations
    formula_k = 10L,           # Covariates in formula
    results_dir = "benchmarks/results",
    scenarios = list(...)      # Scenario definitions
  )
}
```

## Interrupt/Resume Support

The benchmark automatically saves results after each job completes. If interrupted:

1. Simply re-run the same command
2. Completed jobs will be skipped automatically
3. Only pending jobs will be executed

To retry failed jobs:
```bash
Rscript benchmarks/run_full_benchmark.R --retry
```

## Parallel Execution

The benchmark uses the `future` framework for parallel execution:

- By default, uses all available cores minus 1
- Workers are isolated (no shared state)
- Progress is displayed via `progressr`

To control parallelism:
```r
# In config.R, set n_workers:
cfg$n_workers <- 4L  # Use 4 workers
```

## Troubleshooting

### "Cannot determine script directory"

Run from the repository root directory:
```bash
cd /path/to/cfomics
Rscript benchmarks/run_full_benchmark.R
```

Or set the `CFOMICS_ROOT` environment variable:
```bash
export CFOMICS_ROOT=/path/to/cfomics
Rscript benchmarks/run_full_benchmark.R
```

### "cfomics package not found"

Install cfomics first:
```r
devtools::install("packages/cfomics")
```

### "PMCMRplus not available"

The benchmark will fall back to simplified CD-based comparisons if PMCMRplus is not installed:
```r
install.packages("PMCMRplus")
```

### "No RDS files found"

Run the benchmark before aggregation:
```bash
Rscript benchmarks/run_full_benchmark.R --smoke  # Quick test first
Rscript benchmarks/aggregate_results.R
```

### "Memory errors during parallel execution"

Reduce the number of workers:
```r
# In config.R
cfg$n_workers <- 2L
```

Or run sequentially:
```r
# In R/runner.R, change plan to:
plan(sequential)
```

### "BCF taking too long"

BCF uses MCMC and is computationally intensive. For testing, reduce iterations:
```r
# In config.R
cfg$bcf_n_burn <- 100L
cfg$bcf_n_iter <- 200L
```

### "Report generation fails"

Ensure all required packages are installed:
```r
install.packages(c("rmarkdown", "knitr", "kableExtra", "ggplot2", "patchwork"))
```

For PDF output, ensure a LaTeX distribution is installed (TinyTeX recommended):
```r
install.packages("tinytex")
tinytex::install_tinytex()
```

### "Japanese characters not rendering in PDF"

Ensure XeLaTeX is used and appropriate fonts are available:
```r
# The report uses xelatex engine by default
# Install fonts that support Japanese characters
```

## Extending the Benchmark

### Adding a New Method

1. Implement the method in cfomics package
2. Add method name to `config.R`:
   ```r
   methods = c("gformula", "hdml", "hdps", "bcf", "tmle", "newmethod")
   ```

### Adding a New Scenario

Add to the `scenarios` list in `config.R`:
```r
list(id = "S12_custom", dgp = "dgp_custom", n = 1000L, p = 100L, params = list(...))
```

### Customizing Visualizations

Modify functions in `R/reporting.R` or add new plotting functions following the existing patterns.

## Citation

If you use this benchmark suite in your research, please cite:

```bibtex
@software{cfomics_benchmark,
  title = {cfomics Benchmark Suite},
  author = {cfomics Development Team},
  year = {2024},
  url = {https://github.com/your-repo/cfomics}
}
```

## License

This benchmark suite is part of the cfomics package and is licensed under the MIT License.
