# cfomics Genome Biology Paper Design

**Date**: 2026-02-06
**Status**: Draft
**Target Journal**: Genome Biology (primary) / NAR Methods (secondary)

## Overview

### Paper Title (Draft)
"Comprehensive benchmarking of causal inference methods for high-dimensional omics data: practical guidelines and a Bioconductor package"

### Main Contributions
1. Comprehensive benchmark: R-native to Deep Learning methods
2. Omics-specific evaluation metrics
3. Practical method selection guidelines
4. Bioconductor package cfomics with diagnostics

### Differentiation from Vollenweider et al. (2024)
| Aspect | Vollenweider et al. | cfomics (this paper) |
|--------|---------------------|---------------------|
| Focus | Theoretical bias analysis | Practical guidelines |
| Methods | Meta-learners + NN | Traditional + ML + DL + external packages |
| Software | None | Bioconductor package |
| External comparison | No | MatchIt, WeightIt, hdm, tmle3, DoubleML |
| Diagnostics | No | E-value, GATES/BLP, PS diagnostics |

---

## Architecture

### Package vs Benchmark Separation

```
cfomics package (Bioconductor):
├── R-native methods: GRF, TMLE, BCF, HDML, HDPS, IPW, G-formula
├── Python methods (OPTIONAL): DoWhy, DRLearner, GANITE, CAVAE
├── Bioconductor integration: SummarizedExperiment, MAE
└── Diagnostics: E-value, GATES/BLP, PS/outcome model checks

Paper benchmark (comparison only):
├── All cfomics methods
├── External R: MatchIt, WeightIt, hdm, DoubleML, tmle3, SuperLearner
└── External Python: TARNet, CFRNet, DragonNet (not in package)
```

### Directory Structure for Benchmarks

```
benchmarks/
├── external/
│   ├── wrappers/
│   │   ├── wrapper_matchit.R
│   │   ├── wrapper_weightit.R
│   │   ├── wrapper_bart.R
│   │   ├── wrapper_tmle3.R
│   │   ├── wrapper_hdm.R
│   │   ├── wrapper_doubleml.R
│   │   ├── wrapper_superlearner.R
│   │   └── nn/
│   │       ├── tarnet.py
│   │       ├── cfrnet.py
│   │       └── dragonnet.py
│   └── registry.R
├── tcga/
│   ├── data_preparation.R
│   ├── semi_synthetic_dgp.R
│   └── run_tcga_benchmark.R
└── results/
```

---

## Methods to Compare

### Priority 1 (Must include)

| Method | Source | Type |
|--------|--------|------|
| GRF | cfomics | R-native |
| TMLE | cfomics | R-native |
| BCF | cfomics | R-native |
| HDML | cfomics | R-native |
| HDPS | cfomics | R-native |
| IPW | cfomics | R-native |
| G-formula | cfomics | R-native |
| MatchIt | External | R |
| WeightIt | External | R |
| hdm | External | R |
| DoubleML | External | R |
| tmle3 | External | R |
| TARNet | External | Python (comparison only) |
| CFRNet | External | Python (comparison only) |
| DragonNet | External | Python (comparison only) |

### Priority 2 (Recommended)

| Method | Source | Type |
|--------|--------|------|
| DoWhy-GCM | cfomics | Python (optional) |
| DRLearner | cfomics | Python (optional) |
| GANITE | cfomics | Python (optional) |
| SuperLearner | External | R |
| BART | External | R |

### Priority 3 (If feasible)

| Method | Source | Type |
|--------|--------|------|
| CAVAE | cfomics | Python (optional) |
| CausalML | External | Python |
| causallib | External | Python |

---

## Benchmark Scenarios

### Existing (S1-S14, 47 settings)

- S1: Baseline (n/p sweep)
- S2: Dimension sweep (p >> n)
- S3: Heterogeneous linear (strength sweep)
- S4a-d: Heterogeneity types (nonlinear, subgroup, qualitative)
- S5: Nonlinear confounding (5 types)
- S6: Dense confounding
- S7: Weak overlap
- S8: Covariate shift
- S9: Correlated confounding
- S10: Unobserved confounding
- S11: Collider
- S12: Nonlinear outcome
- S13: Nonlinear propensity
- S14: Double nonlinear

### New Metrics to Add

| Metric | Measures | Implementation |
|--------|----------|----------------|
| High-dim stability | Performance degradation vs p/n ratio | Compute from S2 results |
| Computational scalability | Time and memory vs n×p | Add to benchmark_runner |
| Missing data robustness | Performance under MCAR/MAR | New DGP (S15) |
| Sparse ITE detection | TPR for sparse true effects | New metric function |

---

## TCGA Validation

### Strategy: Semi-synthetic

1. Use real TCGA covariates (gene expression)
2. Synthetic treatment: TP53 mutation or predicted drug response
3. Synthetic outcome: Known DGP (true τ known)
4. Compare estimated vs true ATE/ITE

### Datasets

| Dataset | Samples | Use |
|---------|---------|-----|
| TCGA-BRCA | ~1,100 | Primary |
| TCGA-LUAD | ~500 | Secondary |
| TCGA-COAD | ~450 | Secondary |

---

## Paper Structure

### Main Text

1. **Abstract**: Problem, approach, key results, contributions
2. **Background**: Causal inference challenges in omics, existing limitations
3. **Results**:
   - Simulation benchmark results
   - High-dimensional stability analysis
   - TCGA semi-synthetic validation
   - Method selection guidelines
4. **Discussion**: Practical recommendations, limitations, future work
5. **Methods**: cfomics design, DGP specifications, evaluation metrics

### Figures

| Figure | Content |
|--------|---------|
| Fig 1 | cfomics workflow overview |
| Fig 2 | Simulation results (key scenarios) |
| Fig 3 | High-dimensional stability (p/n vs performance) |
| Fig 4 | TCGA semi-synthetic results |
| Fig 5 | Method selection flowchart/guidelines |

### Tables

| Table | Content |
|-------|---------|
| Table 1 | Methods compared (characteristics) |
| Table 2 | Scenario-specific recommendations |

### Supplementary

- Full results for all 47 settings
- Additional figures
- Code availability statement

---

## Task List

### Phase 1: Benchmark Infrastructure

| Task | Priority | Effort |
|------|----------|--------|
| Implement P1 external wrappers (MatchIt, WeightIt, hdm, DoubleML, tmle3) | Highest | Medium |
| Implement TARNet/CFRNet/DragonNet comparison scripts | Highest | Medium |
| Run full S1-S14 benchmark (50 reps × 47 settings × all methods) | Highest | Large |
| Add computational time/memory metrics | Medium | Small |

### Phase 2: TCGA Validation

| Task | Priority | Effort |
|------|----------|--------|
| TCGA data acquisition and preprocessing | High | Medium |
| Implement semi-synthetic DGP | High | Small |
| Run TCGA benchmark | High | Medium |

### Phase 3: New Metrics

| Task | Priority | Effort |
|------|----------|--------|
| Missing data scenario (S15) | Medium | Small |
| Sparse ITE detection metric | Medium | Small |

### Phase 4: Paper Writing

| Task | Priority | Effort |
|------|----------|--------|
| Generate publication figures | High | Medium |
| Write Methods section | High | Medium |
| Write Results + Discussion | High | Large |
| Prepare Supplementary materials | Medium | Medium |

---

## Dependencies

```
Phase 1 (Benchmark) ──┬──→ Phase 4 (Figures) ──→ Phase 4 (Results)
Phase 2 (TCGA) ───────┘
Phase 3 (Metrics) ────→ Phase 1 (re-run with new metrics)
```

---

## Timeline Considerations

- Phase 1 is blocking for paper figures
- TCGA and simulation can run in parallel
- Writing can start once Phase 1 results are available

---

## Open Questions

1. Should we include DepMap data for true-outcome validation (like Vollenweider)?
2. Exact scope of "practical guidelines" - flowchart vs decision tree vs prose?
3. Bioconductor submission timing - before or after paper acceptance?
