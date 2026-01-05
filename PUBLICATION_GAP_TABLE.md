# cfomics Publication Pathway Gap Analysis

This table summarizes the gaps between current package state and requirements for each publication target.

## Gap Overview by Target

| Requirement | CRAN | Bioconductor | JOSS | Current Status |
|-------------|------|--------------|------|----------------|
| **Author Info** | Required | Required | Required | Placeholder (P0) |
| **Examples Run** | Required | Required | Nice-to-have | All \dontrun (P1) |
| **No .onLoad Side Effects** | Required | Required | N/A | Violates (P0) |
| **biocViews** | Must remove | Required | N/A | Present |
| **LazyData** | Must match data/ | Must match data/ | N/A | Mismatch (P2) |
| **SE/MAE Integration** | Optional | Expected | Expected | Basic (P1/P2) |
| **Example Data** | Optional | Expected | Nice-to-have | None |
| **Vignettes** | Required | Required | Required | 4 present |
| **Unit Tests** | Expected | Required | Required | 20 files |
| **BiocCheck** | N/A | Required | N/A | Not tested |
| **R CMD check** | 0 errors/warnings | 0 errors/warnings | Expected | Unknown |
| **Reproducibility** | N/A | Expected | Required | Gaps (P1) |
| **Estimand Documentation** | N/A | Nice-to-have | Required | Missing (P0) |
| **Comparison Study** | N/A | Nice-to-have | Required | Basic benchmarks |

---

## CRAN Submission Requirements

### Must Fix (Blockers)
| Issue | Description | Effort |
|-------|-------------|--------|
| CFOMICS-REV-001 | Real author info | Low |
| CFOMICS-REV-021 | Remove .onLoad Python init | Low |
| Remove biocViews | Edit DESCRIPTION | Low |
| CFOMICS-REV-002 | Fix duplicate exports | Low |
| CFOMICS-REV-003 | Runnable examples | Medium |
| CFOMICS-REV-005 | Fix LazyData | Low |

### Should Fix
| Issue | Description | Effort |
|-------|-------------|--------|
| CFOMICS-REV-004 | Complete @return docs | Low |
| CFOMICS-REV-007 | Move igraph to Suggests | Low |
| CFOMICS-REV-008 | Remove unused gfoRmula | Low |

### CRAN Check Preparation
```bash
R CMD build .
R CMD check --as-cran cfomics_0.3.0.tar.gz
```

Expected issues to address:
- NOTE: No visible global function definition (fix with `@importFrom`)
- NOTE: Examples run > 5s (use `\donttest{}`)
- WARNING: Dependencies not available (Suggests must be handled)

---

## Bioconductor Submission Requirements

### Must Fix (All CRAN requirements plus)
| Issue | Description | Effort |
|-------|-------------|--------|
| Add biocViews | Already present (keep) | None |
| CFOMICS-REV-028/029 | Feature selection reproducibility | Medium |
| CFOMICS-REV-031 | Document omics workflow | Medium |
| Add example data | Create SummarizedExperiment example | Medium |
| BiocCheck compliance | Run BiocCheck, fix issues | Variable |

### Should Fix
| Issue | Description | Effort |
|-------|-------------|--------|
| CFOMICS-REV-030 | Improve MAE formula support | High |
| CFOMICS-REV-032 | NA handling options | Low |

### Bioconductor-Specific Actions
1. Use BiocManager for installation instructions
2. Ensure all SE/MAE operations use Bioconductor generics
3. Add `inst/NEWS.Rd` or `NEWS.md`
4. Follow naming conventions (camelCase for functions)
5. Add `BugReports:` to DESCRIPTION (already present)

### BiocCheck Preparation
```r
BiocManager::install("BiocCheck")
BiocCheck::BiocCheck("cfomics", `new-package` = TRUE)
```

---

## JOSS Submission Requirements

### Must Address (Scholarly)
| Requirement | Description | Current Gap | Effort |
|-------------|-------------|-------------|--------|
| Statement of Need | Why this package? | Not written | Medium |
| Target Audience | Who benefits? | Implicit | Low |
| State of Field | Existing alternatives | Not documented | Medium |
| CFOMICS-REV-015/016/017 | Estimand clarity | Missing | High |
| Reproducibility | Seeds, versions | Gaps | Medium |
| CFOMICS-REV-033 | Omics-relevant benchmarks | Limited | High |
| Comparison | vs grf, DoubleML, etc. | Not done | High |

### Paper Structure (paper.md)
```markdown
# Summary
- What cfomics does
- Why it's needed

# Statement of Need
- Gap in existing tools
- Target users (bioinformaticians, epidemiologists)

# Functionality
- Core API (cf_fit, predict, summary)
- Supported methods
- Bioconductor integration

# Illustrative Example
- Realistic omics use case

# Comparison with Related Work
- grf (R): Single method
- EconML (Python): Not R-integrated
- DoubleML: Different focus
- cfomics: Unified R interface for omics

# Acknowledgments
# References
```

### JOSS Checklist
- [ ] paper.md written
- [ ] Zenodo DOI for code release
- [ ] CITATION.cff file
- [ ] License file (MIT present)
- [ ] Installation instructions work
- [ ] Example reproduces in <10 min
- [ ] Tests pass
- [ ] Documentation complete

---

## Recommended Path

### Path A: CRAN First (Faster, Lower Bar)
**Timeline: 2-4 weeks**

1. Week 1: Fix P0 issues (author info, .onLoad, estimand)
2. Week 1: Remove biocViews, fix NAMESPACE
3. Week 2: Add runnable examples, fix LazyData
4. Week 2-3: R CMD check iterations
5. Week 3-4: Submit to CRAN

**Pros:** Faster, establishes package presence
**Cons:** May need to resubmit for Bioconductor later

### Path B: Bioconductor First (Higher Bar, Better Fit)
**Timeline: 4-8 weeks**

1. Week 1-2: Fix all CRAN requirements
2. Week 2-3: Enhance SE/MAE integration
3. Week 3-4: Add example data, improve vignettes
4. Week 4-6: BiocCheck iterations
5. Week 6-8: Review process

**Pros:** Better target for omics focus, Bioconductor visibility
**Cons:** Longer, stricter requirements

### Path C: JOSS with CRAN (Academic Impact)
**Timeline: 6-10 weeks**

1. Week 1-4: CRAN submission (Path A)
2. Week 5-6: Write paper.md
3. Week 6-8: Run comparative benchmarks
4. Week 8-10: JOSS review

**Pros:** Citable publication, dual presence
**Cons:** Most work, need benchmark results

---

## Effort Estimates

| Task Category | CRAN | Bioconductor | JOSS |
|---------------|------|--------------|------|
| P0 Fixes | 2-3 days | 2-3 days | 2-3 days |
| P1 Fixes | 3-5 days | 5-7 days | 5-7 days |
| Documentation | 2-3 days | 5-7 days | 7-10 days |
| Examples/Data | 1-2 days | 3-5 days | 3-5 days |
| Review Iteration | 1-2 weeks | 2-4 weeks | 2-4 weeks |
| **Total** | **2-4 weeks** | **4-8 weeks** | **6-10 weeks** |

---

## Quick Wins (Do First)

1. **Fix author info** (DESCRIPTION:6-9) - 5 minutes
2. **Remove .onLoad Python init** (aaa-onload.R:14) - 5 minutes
3. **Remove duplicate exports** (NAMESPACE:53-57) - 5 minutes
4. **Remove biocViews if targeting CRAN** - 5 minutes
5. **Add estimand field to meta** (cf_fit.R) - 30 minutes

These 5 changes address 3 P0 issues and 1 P1 issue in under an hour.
