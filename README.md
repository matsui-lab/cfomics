# cfomics Ecosystem

Counterfactual Causal Inference for Omics Data - A modular R package ecosystem.

## Packages

| Package | Description | Status |
|---------|-------------|--------|
| [cfomics](packages/cfomics/) | Core API and R-native methods (GRF, IPW, G-formula) | Active |
| cfomicsPython | Python backends (DoWhy, EconML, GANITE, CAVAE) | Planned |
| cfomicsSim | Simulation and Data Generating Processes | Planned |
| cfomicsBench | Benchmarking framework | Planned |
| cfomicsAdapters | Bioconductor integration (SE, MAE, SCE) | Planned |
| cfomicsDiagnose | Diagnostics and assumption checking | Planned |

## Installation

### Core Package

```r
# Install from GitHub
remotes::install_github("matsui-lab/cfomics", subdir = "packages/cfomics")

# Or with devtools
devtools::install_github("matsui-lab/cfomics", subdir = "packages/cfomics")
```

### Development Installation

```r
# Clone the repository
# git clone https://github.com/matsui-lab/cfomics.git
# cd cfomics

# Install core package locally
devtools::install("packages/cfomics")
```

## Quick Start

```r
library(cfomics)

# Generate example data
set.seed(42)
n <- 500
data <- data.frame(
  X1 = rnorm(n),
  X2 = rnorm(n),
  T = rbinom(n, 1, 0.5),
  Y = rnorm(n)
)
data$Y <- data$Y + 2 * data$T + 0.5 * data$X1

# Fit causal model using formula syntax: Y ~ T | covariates
result <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")

# View results
summary(result)
plot(result)
```

## Available Methods

### R-native (cfomics core)

| Method | Description |
|--------|-------------|
| `grf` | Generalized Random Forest for heterogeneous treatment effects |
| `ipw` | Inverse Probability Weighting estimator |
| `gformula` | G-computation formula estimator |

### Python-based (cfomicsPython, coming soon)

| Method | Backend | Description |
|--------|---------|-------------|
| `dowhy_gcm` | DoWhy | Graphical Causal Model |
| `drlearner` | EconML | Double Robust Learner |
| `ganite` | TensorFlow | Generative Adversarial Nets for ITE |
| `cavae` | PyTorch/Pyro | Causal Autoencoder VAE |

## Documentation

- [Core Package Documentation](packages/cfomics/README.md)
- [Python Setup Guide](packages/cfomics/vignettes/python_setup.Rmd)
- [Method Comparison](packages/cfomics/vignettes/method_comparison.Rmd)

## Architecture

This repository uses a monorepo structure with multiple R packages under `packages/`.

See [ECOSYSTEM_SPEC_V2.md](ECOSYSTEM_SPEC_V2.md) for the ecosystem design specification.

```
cfomics/
├── packages/
│   ├── cfomics/          # Core package
│   ├── cfomicsPython/    # Python backends (planned)
│   ├── cfomicsSim/       # Simulation (planned)
│   └── cfomicsBench/     # Benchmarking (planned)
├── tools/                # Monorepo development tools
├── CLAUDE.md             # AI development guide
└── ECOSYSTEM_SPEC_V2.md  # Design specification
```

## License

MIT License - see [LICENSE.md](packages/cfomics/LICENSE.md)

## Contributing

See [CLAUDE.md](CLAUDE.md) for development guidelines and coding conventions.

## Citation

If you use cfomics in your research, please cite:

```
@software{cfomics,
  author = {Matsui, Yusuke},
  title = {cfomics: Counterfactual Causal Inference for Omics Data},
  url = {https://github.com/matsui-lab/cfomics}
}
```
