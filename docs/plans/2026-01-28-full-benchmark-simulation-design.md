# Full Benchmark Simulation Design

## Overview

Production-scale benchmark simulation comparing 5 causal inference methods (gformula, hdml, hdps, bcf, tmle) across 13 DGP scenarios (S1-S11) with 50 replications. Output includes a publication-quality report in English and Japanese with figures, LaTeX tables, statistical tests, and interpretation.

## Directory Structure

```
benchmarks/
├── run_full_benchmark.R       # Main entry point
├── aggregate_results.R        # Collect raw RDS → summary CSV
├── config.R                   # Centralized parameter configuration
├── R/
│   ├── scenarios.R            # S1-S11 scenario definitions and DGP dispatch
│   ├── runner.R               # Parallel execution engine (future + furrr)
│   ├── metrics.R              # Metric computation + statistical tests
│   └── reporting.R            # Visualization + LaTeX table generation
├── report/
│   ├── benchmark_report.Rmd   # Main report (parameterized: lang=en/ja)
│   ├── sections/
│   │   ├── 01_introduction_en.Rmd
│   │   ├── 01_introduction_ja.Rmd
│   │   ├── 02_methods_en.Rmd
│   │   ├── 02_methods_ja.Rmd
│   │   ├── 03_simulation_design_en.Rmd
│   │   ├── 03_simulation_design_ja.Rmd
│   │   ├── 04_results_en.Rmd
│   │   ├── 04_results_ja.Rmd
│   │   ├── 05_discussion_en.Rmd
│   │   ├── 05_discussion_ja.Rmd
│   │   ├── 06_appendix_en.Rmd
│   │   └── 06_appendix_ja.Rmd
│   ├── figures/               # High-resolution figures (300dpi PDF/PNG)
│   ├── tables/                # LaTeX .tex files
│   └── bibliography.bib       # References
├── results/                   # Output directory (.gitignore)
│   ├── raw/                   # Per-replicate RDS files
│   ├── summary/               # Aggregated CSV
│   ├── plots/                 # Generated plots
│   └── tables/                # Generated LaTeX tables
└── README.md                  # How to run
```

## Methods

5 R-native methods:

| Method | Estimator Type | Key Property |
|--------|---------------|-------------|
| gformula | Outcome modeling | Simple, fast, not doubly robust |
| hdml | Debiased LASSO + AIPW | Doubly robust, high-dimensional |
| hdps | LASSO propensity + IPW | Propensity-based, high-dimensional |
| bcf | Bayesian CART (MCMC) | Posterior inference, heterogeneity |
| tmle | Targeted MLE | Doubly robust, semiparametric efficient |

## Scenarios

| ID | DGP Function | Varying Parameters | n | p |
|----|-------------|-------------------|---|---|
| S1 | `dgp_baseline` | - | 500, 1000, 2000 | 50, 100 |
| S2 | `dgp_dimension_sweep` | - | 500, 1000 | 100, 500, 1000 |
| S3 | `dgp_heterogeneous_linear` | strength: 0, 0.5, 1.0, 2.0 | 1000 | 100 |
| S4a | `dgp_heterogeneous_nonlinear` | - | 1000 | 100 |
| S4b | `dgp_heterogeneous_subgroup` | - | 1000 | 100 |
| S4c | `dgp_heterogeneous_qualitative` | - | 1000 | 100 |
| S5 | `dgp_nonlinear_confounding` | type: quadratic, trigonometric, interaction, combined, threshold | 1000 | 100 |
| S6 | `dgp_dense_confounding` | n_conf: 10, 50, 100, 200 | 1000 | 500 |
| S7 | `dgp_weak_overlap` | strength: good, moderate, weak, extreme | 1000 | 100 |
| S8 | `dgp_covariate_shift` | type: mean, variance | 1000 | 100 |
| S9 | `dgp_correlated_confounding` | type: block, ar1, factor | 1000 | 100 |
| S10 | `dgp_unobserved_confounding` | strength: 0, 0.5, 1.0, 2.0 | 1000 | 100 |
| S11 | `dgp_collider` | strength: 0, 0.5, 1.0, 2.0 | 1000 | 100 |

**Total:** ~50 scenario settings x 5 methods x 50 reps = ~12,500 jobs

## Parallel Execution

- Engine: `future` + `furrr` with `plan(multisession)`
- Workers: `parallel::detectCores() - 1`
- Parallelization unit: 1 job = 1 scenario setting x 1 method x 1 rep
- Each job has an independent seed for reproducibility
- Future-compatible with HPC (OpenPBS) via `future.batchtools`:
  ```r
  # Local
  plan(multisession, workers = 7)
  # HPC (swap later if needed)
  library(future.batchtools)
  plan(batchtools_openpbs, resources = list(nodes = 1, ncpus = 4, walltime = "2:00:00"))
  ```

## Interrupt/Resume

- Each job saves result to `results/raw/{scenario}__{method}__n{n}_p{p}__rep{NNN}.rds`
- Before execution, scan `results/raw/` and skip existing files
- Re-run to process only incomplete jobs
- Error jobs saved with `status = "error"` (retry option available)

## Metrics

### Per-replicate metrics

| Metric | Definition | Target |
|--------|-----------|--------|
| Bias | ATE_hat - ATE_true | ATE accuracy |
| RMSE | sqrt(mean(bias^2)) | ATE accuracy |
| PEHE | sqrt(mean((ITE_hat - ITE_true)^2)) | ITE accuracy |
| Coverage | 95% CI contains ATE_true | Uncertainty quantification |
| CI Length | Upper - Lower CI | Uncertainty quantification |
| Time (sec) | fit + predict elapsed | Computational cost |

### Aggregated metrics (over 50 reps)

- Mean, SD, Median, IQR for each metric
- Coverage with Bernoulli 95% CI

### Statistical tests

- **Friedman test** per scenario: nonparametric rank comparison across methods
- **Nemenyi post-hoc**: pairwise comparisons when Friedman is significant
- **Critical Difference (CD) diagram**: rank visualization across all scenarios

## Figures

| # | Figure | Purpose |
|---|--------|---------|
| Fig 1 | DAG diagrams (igraph) | Causal structure of each scenario |
| Fig 2 | Method comparison heatmap | RMSE/PEHE across method x scenario |
| Fig 3 | Bias boxplots | Distribution shape comparison |
| Fig 4 | Coverage heatmap | Actual vs nominal 95% coverage |
| Fig 5 | CD diagram | Friedman/Nemenyi rank comparison |
| Fig 6 | Sweep curves (6 panels) | Parameter sensitivity (S3, S5, S6, S7, S10, S11) |
| Fig 7 | Bias-Variance decomposition | RMSE = Bias^2 + Variance stacked bars |
| Fig 8 | Computation time (log scale) | Scalability comparison |

All figures: 300dpi, PDF + PNG, publication-quality with ggplot2 + patchwork.

## Tables

| # | Table | Content |
|---|-------|---------|
| Tab 1 | Main results | RMSE, PEHE, Bias (bold=best, underline=2nd) |
| Tab 2 | Coverage + CI Length | Uncertainty quantification evaluation |
| Tab 3 | Computation time | Median seconds per method |
| Tab 4 | Friedman test + ranks | p-values and mean ranks |
| Tab 5 | Method recommendation matrix | Scenario characteristic -> recommended method |

All tables: `kableExtra` LaTeX output, directly `\input{}`-able in papers.

## Report Structure (R Markdown -> PDF)

Parameterized report with `params$lang` (en/ja). Each section has English and Japanese variants.

| Section | Content |
|---------|---------|
| §1 Introduction | High-dimensional causal inference challenges, limitations of existing methods, benchmark objectives |
| §2 Methods | Estimator definitions, theoretical properties (double robustness, regularity conditions, Bayesian posteriors), expected strengths/weaknesses |
| §3 Simulation Design | Why each S1-S11 scenario is needed, mathematical DGP definitions, evaluation metric definitions |
| §4 Results | Per-scenario results with captioned figures and tables |
| §5 Discussion | Cross-scenario interpretation, practical method selection guidelines, limitations, future work |
| §6 Appendix | Full parameter tables, supplementary figures, session info |

## Execution Flow

```bash
# Step 1: Run simulations (3-6 hours on 8-core local machine)
Rscript benchmarks/run_full_benchmark.R

# Step 2: Aggregate results
Rscript benchmarks/aggregate_results.R

# Step 3: Generate reports (English + Japanese)
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="en"))'
Rscript -e 'rmarkdown::render("benchmarks/report/benchmark_report.Rmd", params=list(lang="ja"))'
```

## Dependencies

| Package | Purpose |
|---------|---------|
| future + furrr | Parallel execution |
| kableExtra | LaTeX tables |
| ggplot2 + patchwork | Figure layout |
| scales | Axis formatting |
| PMCMRplus | Friedman/Nemenyi tests |
| rmarkdown + bookdown | Report PDF generation |
| igraph | DAG drawing (Fig 1) |
| tikzDevice (optional) | LaTeX-quality figures |

## Estimated Runtime

- gformula/hdml/hdps/tmle: ~1-5 sec/job
- bcf (MCMC, n_burn=1000, n_iter=2000): ~30-120 sec/job
- Total ~12,500 jobs on 8 cores: **3-6 hours**
