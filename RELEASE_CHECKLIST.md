# cfomics Release Checklist

This checklist ensures the package meets CRAN/Bioconductor quality standards before release.

## Pre-Release Checks

### 1. Package Quality (DoD-1)

- [ ] `R CMD check --as-cran` produces 0 ERROR, 0 WARNING
- [ ] Run with `_R_CHECK_CRAN_INCOMING_=true` for CRAN-level checks
- [ ] All tests pass: `devtools::test()`
- [ ] No NOTEs related to code quality (some informational NOTEs are acceptable)

```bash
cd cfomics
R CMD build .
_R_CHECK_CRAN_INCOMING_=true R CMD check --as-cran cfomics_*.tar.gz
```

Or use the helper script:
```r
source("tools/check_package.R")
```

### 2. Documentation (DoD-2)

- [ ] All exported functions have `@return` documentation
- [ ] All exported functions have `@examples` that run (not wrapped in `\dontrun{}`)
- [ ] Examples complete in < 10 seconds each
- [ ] Vignettes build without Python
- [ ] `devtools::document()` runs without warnings
- [ ] `devtools::check_man()` passes

### 3. Reproducibility (DoD-3)

- [ ] `cf_session_info()` captures all relevant version information
- [ ] `random_state` parameter is documented and works
- [ ] Results are reproducible with the same seed
- [ ] `meta` object includes software versions

### 4. API Contract (DoD-4)

- [ ] All methods return consistent `cfomics_result` structure
- [ ] `validate_cfomics_result()` passes for all methods
- [ ] `predict()` works for all types: ate, ite, y0, y1, summary, samples
- [ ] `summary()` works for all methods
- [ ] R-native methods serialize correctly (saveRDS/readRDS)

### 5. Test Coverage (DoD-5)

- [ ] Python tests skip gracefully when Python unavailable
- [ ] CI passes on all platforms (Ubuntu, macOS, Windows)
- [ ] CI passes with R release and R devel
- [ ] No tests depend on specific file paths or system configuration

## Version Bump

1. Update `DESCRIPTION`:
   - [ ] Increment version number
   - [ ] Check `Authors@R` is correct
   - [ ] Check `Maintainer` email is current

2. Update `NEWS.md`:
   - [ ] Add section for new version
   - [ ] Document breaking changes
   - [ ] Document new features
   - [ ] Document bug fixes

## Final Steps

### For CRAN Submission

1. [ ] Run `rhub::check_for_cran()` for additional checks
2. [ ] Run `devtools::spell_check()` for typos
3. [ ] Run `urlchecker::url_check()` for broken URLs
4. [ ] Submit via https://cran.r-project.org/submit.html
5. [ ] Monitor CRAN-incoming queue

### For Bioconductor

1. [ ] Ensure `biocViews` are appropriate
2. [ ] Check BiocCheck: `BiocCheck::BiocCheck(".")`
3. [ ] Follow Bioconductor submission guidelines

## Post-Release

- [ ] Tag release in git: `git tag -a v0.x.0 -m "Release 0.x.0"`
- [ ] Push tags: `git push origin --tags`
- [ ] Create GitHub release with changelog
- [ ] Update website/documentation if applicable

## Quick Commands

```bash
# Full check workflow
cd cfomics
Rscript -e "devtools::document()"
R CMD build .
R CMD check --as-cran cfomics_*.tar.gz

# Test only
Rscript -e "devtools::test()"

# Check examples
Rscript -e "devtools::run_examples()"

# Build vignettes
Rscript -e "devtools::build_vignettes()"
```

## Known Issues / Exceptions

- Python-based methods (dowhy_gcm, ganite, drlearner, cavae) require external Python setup
- Some NOTEs about package size may appear (acceptable if < 5MB)
- Bioconductor dependencies are in Suggests, not Imports
