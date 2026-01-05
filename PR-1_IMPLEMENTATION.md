# PR-1: packages/ モノレポ構造への移行

## 目的

現在の単一パッケージ構造を `packages/cfomics/` 配下に移動し、将来のマルチパッケージ（cfomicsPython, cfomicsSim 等）追加に備えたモノレポ構造を確立する。

## 前提条件

- [ ] main ブランチが安定（`R CMD check` が通る状態）
- [ ] ベースラインタグ `cfomics-monolith-v0.4.0` を打つ（推奨）

## ブランチ

```bash
git checkout -b feature/monorepo-structure
```

---

## 作業手順

### Step 1: ベースラインタグの作成（推奨）

```bash
git tag -a cfomics-monolith-v0.4.0 -m "Pre-monorepo baseline: single package structure"
```

### Step 2: packages/ ディレクトリの作成

```bash
mkdir -p packages/cfomics
```

### Step 3: パッケージファイルの移動

```bash
# R source files
git mv R packages/cfomics/

# Documentation
git mv man packages/cfomics/
git mv vignettes packages/cfomics/

# Tests
git mv tests packages/cfomics/

# Data
git mv data-raw packages/cfomics/

# Inst (Python, benchmarks, scripts)
git mv inst packages/cfomics/

# Package metadata
git mv DESCRIPTION packages/cfomics/
git mv NAMESPACE packages/cfomics/
git mv LICENSE packages/cfomics/
git mv LICENSE.md packages/cfomics/
git mv NEWS.md packages/cfomics/

# Package README (現在のものをパッケージ用に移動)
git mv README.md packages/cfomics/
```

### Step 4: testthat.R の作成（欠落している場合）

現在 `tests/testthat.R` が存在しないため、作成が必要：

```bash
cat > packages/cfomics/tests/testthat.R << 'EOF'
# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(cfomics)

test_check("cfomics")
EOF
```

### Step 5: ルート README の作成

```bash
cat > README.md << 'EOF'
# cfomics Ecosystem

Counterfactual Causal Inference for Omics Data - A modular R package ecosystem.

## Packages

| Package | Description | Status |
|---------|-------------|--------|
| [cfomics](packages/cfomics/) | Core API and R-native methods | [![R-CMD-check](https://github.com/matsui-lab/cfomics/workflows/R-CMD-check/badge.svg)](https://github.com/matsui-lab/cfomics/actions) |
| cfomicsPython | Python backend (DoWhy, EconML, GANITE) | Planned |
| cfomicsSim | Simulation and DGP | Planned |
| cfomicsBench | Benchmarking harness | Planned |

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

# Fit causal model
result <- cf_fit(Y ~ T | X1 + X2, data = data, method = "grf")

# View results
summary(result)
plot(result)
```

## Documentation

- [Core Package Documentation](packages/cfomics/README.md)
- [Python Setup Guide](packages/cfomics/vignettes/python_setup.Rmd)
- [Method Comparison](packages/cfomics/vignettes/method_comparison.Rmd)

## Architecture

See [ECOSYSTEM_SPEC_V2.md](ECOSYSTEM_SPEC_V2.md) for the ecosystem design specification.

## License

MIT License - see [LICENSE.md](packages/cfomics/LICENSE.md)

## Contributing

See [CLAUDE.md](CLAUDE.md) for development guidelines.
EOF
```

### Step 6: ルート .Rbuildignore の更新

ルートレベルには R パッケージがないため、`.Rbuildignore` は `packages/cfomics/` に移動するか、ルートに残して CI/ツール用とする。

**Option A: パッケージ用の .Rbuildignore を作成**

```bash
cat > packages/cfomics/.Rbuildignore << 'EOF'
^.*\.Rproj$
^\.Rproj\.user$
^data-raw$
^LICENSE\.md$
^\.github$
^codecov\.yml$
^_pkgdown\.yml$
^docs$
^pkgdown$
^cran-comments\.md$
^renv\.lock$
^renv$
EOF
```

**Option B: 現在の .Rbuildignore を移動**

```bash
git mv .Rbuildignore packages/cfomics/.Rbuildignore
```

その後、パッケージ用に編集（v2 関連ドキュメントの除外を削除）：

```bash
cat > packages/cfomics/.Rbuildignore << 'EOF'
^.*\.Rproj$
^\.Rproj\.user$
^data-raw$
^LICENSE\.md$
^\.github$
^codecov\.yml$
^_pkgdown\.yml$
^docs$
^pkgdown$
^cran-comments\.md$
^renv\.lock$
^renv$
EOF
```

### Step 7: ルートレベルの .gitignore 更新

```bash
cat > .gitignore << 'EOF'
# R
.Rproj.user
.Rhistory
.Rdata
.httr-oauth
.DS_Store

# Package build artifacts
*.tar.gz
*.Rcheck/

# IDE
.vscode/
.idea/

# Python
__pycache__/
*.py[cod]
.Python
*.egg-info/
.eggs/

# Environment
.env
.venv/
venv/

# renv (if used at root level)
renv/library/
renv/staging/
renv/cellar/
EOF
```

### Step 8: ルート tools/ の作成

```bash
mkdir -p tools
```

**tools/check.R**

```bash
cat > tools/check.R << 'EOF'
#!/usr/bin/env Rscript
# Check all packages in packages/

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

packages <- list.dirs("packages", recursive = FALSE, full.names = TRUE)

if (length(packages) == 0) {
  message("No packages found in packages/")
  quit(status = 1)
}

failed <- character()

for (pkg in packages) {
  pkg_name <- basename(pkg)
  message("\n", strrep("=", 60))
  message("Checking: ", pkg_name)
  message(strrep("=", 60), "\n")

  tryCatch({
    if (quick) {
      # Quick check: no vignettes, no manual
      result <- devtools::check(
        pkg,
        document = FALSE,
        build_args = c("--no-build-vignettes", "--no-manual"),
        args = c("--no-vignettes", "--no-manual"),
        quiet = FALSE
      )
    } else {
      result <- devtools::check(pkg, quiet = FALSE)
    }

    if (length(result$errors) > 0) {
      failed <- c(failed, pkg_name)
    }
  }, error = function(e) {
    message("Error checking ", pkg_name, ": ", e$message)
    failed <<- c(failed, pkg_name)
  })
}

message("\n", strrep("=", 60))
if (length(failed) > 0) {
  message("FAILED packages: ", paste(failed, collapse = ", "))
  quit(status = 1)
} else {
  message("All packages passed!")
  quit(status = 0)
}
EOF
chmod +x tools/check.R
```

**tools/check_changed.R**

```bash
cat > tools/check_changed.R << 'EOF'
#!/usr/bin/env Rscript
# Check only packages with changes (based on git diff)

# Get changed files
changed_files <- system("git diff --name-only HEAD~1 2>/dev/null || git diff --name-only HEAD", intern = TRUE)

if (length(changed_files) == 0) {
  message("No changed files detected")
  quit(status = 0)
}

# Extract package names from changed paths
pkg_pattern <- "^packages/([^/]+)/"
matches <- regmatches(changed_files, regexec(pkg_pattern, changed_files))
changed_packages <- unique(unlist(lapply(matches, function(m) if (length(m) > 1) m[2] else NULL)))

if (length(changed_packages) == 0) {
  message("No package changes detected")
  quit(status = 0)
}

message("Changed packages: ", paste(changed_packages, collapse = ", "))

failed <- character()

for (pkg_name in changed_packages) {
  pkg_path <- file.path("packages", pkg_name)

  if (!dir.exists(pkg_path)) {
    message("Package directory not found: ", pkg_path)
    next
  }

  message("\n", strrep("=", 60))
  message("Checking: ", pkg_name)
  message(strrep("=", 60), "\n")

  tryCatch({
    result <- devtools::check(pkg_path, quiet = FALSE)
    if (length(result$errors) > 0) {
      failed <- c(failed, pkg_name)
    }
  }, error = function(e) {
    message("Error checking ", pkg_name, ": ", e$message)
    failed <<- c(failed, pkg_name)
  })
}

message("\n", strrep("=", 60))
if (length(failed) > 0) {
  message("FAILED packages: ", paste(failed, collapse = ", "))
  quit(status = 1)
} else {
  message("All changed packages passed!")
  quit(status = 0)
}
EOF
chmod +x tools/check_changed.R
```

### Step 9: CLAUDE.md の更新

CLAUDE.md のパス参照を更新：

```bash
# CLAUDE.md を編集して packages/cfomics/ 構造に対応させる
# 主な変更点:
# - Repository Structure セクションのパス更新
# - Running Tests のパス更新
# - Package Checks のパス更新
```

具体的な編集内容（主要部分）:

```markdown
## Repository Structure

```
cfomics/
├── packages/
│   └── cfomics/              # Core package
│       ├── R/                # R source files
│       ├── inst/
│       │   ├── python/       # Python backend scripts
│       │   ├── benchmarks/   # Benchmark scripts
│       │   └── scripts/      # Utility scripts
│       ├── tests/testthat/   # Unit tests
│       ├── man/              # Generated documentation
│       ├── vignettes/        # Package vignettes
│       ├── data-raw/         # Data generation scripts
│       ├── DESCRIPTION
│       └── NAMESPACE
├── tools/                    # Monorepo tools
│   ├── check.R
│   └── check_changed.R
├── README.md                 # Ecosystem overview
├── CLAUDE.md                 # AI assistant guide
├── ECOSYSTEM_SPEC_V2.md      # Design specification
└── IMPLEMENTATION_PLAN_V2.md # Implementation plan
```
```

### Step 10: .Rproj ファイルの処理

`cfomics.Rproj` はルートに残すか、`packages/cfomics/` に移動するか選択：

**Option A: ルートに残す（モノレポ開発用）**

```bash
# そのまま残す（ルートでRStudioを開く用）
```

**Option B: パッケージに移動**

```bash
git mv cfomics.Rproj packages/cfomics/
```

**推奨**: ルートに残し、パッケージ用に別途作成も可能

### Step 11: 残りのルートファイルの整理

以下のファイルはルートに残す（エコシステム管理用）：

```
.gitignore           # ルート用
CLAUDE.md            # AI開発ガイド
ECOSYSTEM_SPEC_V2.md # 仕様書
IMPLEMENTATION_PLAN_V2.md # 実装計画
AUDIT_FINDINGS.md    # 監査結果（ルートで管理）
AUDIT_ISSUES.md      # 監査課題
PUBLICATION_GAP_TABLE.md # 論文ギャップ
RELEASE_CHECKLIST.md # リリース手順
renv.lock            # ルートレベルの依存管理（任意）
cfomics.Rproj        # RStudio プロジェクト
```

---

## 検証手順

### Step A: ディレクトリ構造の確認

```bash
# 期待される構造
tree -L 3 packages/
# または
find packages -maxdepth 3 -type f | head -20
```

期待出力：
```
packages/
└── cfomics/
    ├── DESCRIPTION
    ├── LICENSE
    ├── LICENSE.md
    ├── NAMESPACE
    ├── NEWS.md
    ├── R/
    │   ├── aaa-onload.R
    │   ├── cf_fit.R
    │   └── ...
    ├── inst/
    │   ├── python/
    │   └── ...
    ├── man/
    ├── tests/
    │   ├── testthat/
    │   └── testthat.R
    └── vignettes/
```

### Step B: パッケージビルドの確認

```bash
cd packages/cfomics
R CMD build .
```

期待: `cfomics_0.4.0.tar.gz` が生成される

### Step C: パッケージチェックの確認

```bash
cd packages/cfomics
R CMD check cfomics_*.tar.gz --as-cran
```

期待: ERROR なし（WARNING/NOTE は許容）

### Step D: devtools での確認

```r
# R コンソールで
devtools::check("packages/cfomics")
```

### Step E: インストールテスト

```r
# ローカルインストール
devtools::install("packages/cfomics")

# 動作確認
library(cfomics)
cf_fit(Y ~ T | X1, data = data.frame(Y=1:10, T=rep(0:1,5), X1=rnorm(10)), method="gformula")
```

---

## DoD チェックリスト

- [ ] `packages/cfomics/` に全パッケージファイルが移動済み
- [ ] `packages/cfomics/tests/testthat.R` が存在する
- [ ] `R CMD build packages/cfomics` が成功する
- [ ] `R CMD check` が ERROR なしで通る
- [ ] ルート README にインストール手順が記載されている
- [ ] `tools/check.R` が動作する
- [ ] CLAUDE.md のパス参照が更新されている

---

## コミット手順

### Commit 1: ファイル移動

```bash
git add -A
git commit -m "refactor: move package files to packages/cfomics/

Prepare monorepo structure for future package separation.
All package files moved to packages/cfomics/ subdirectory.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 2: ルートファイル作成

```bash
git add README.md tools/ .gitignore
git commit -m "feat: add monorepo root files and tools

- Add ecosystem README with installation instructions
- Add tools/check.R for checking all packages
- Add tools/check_changed.R for CI optimization
- Update .gitignore for monorepo structure

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 3: testthat.R 追加（必要な場合）

```bash
git add packages/cfomics/tests/testthat.R
git commit -m "fix: add missing testthat.R

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 4: CLAUDE.md 更新

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md for monorepo structure

Update file paths and repository structure documentation
to reflect packages/ layout.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## PR 作成

```bash
git push -u origin feature/monorepo-structure

gh pr create --title "refactor: migrate to packages/ monorepo structure" --body "$(cat <<'EOF'
## Summary

- Move all package files to `packages/cfomics/` subdirectory
- Add ecosystem-level README with installation instructions
- Add `tools/check.R` and `tools/check_changed.R` for monorepo management
- Update CLAUDE.md for new structure

This prepares the repository for future package separation (cfomicsPython, cfomicsSim, cfomicsBench, etc.) while maintaining backward compatibility.

## Changes

- `packages/cfomics/` - Core package (all existing files)
- `README.md` - Ecosystem overview
- `tools/check.R` - Check all packages
- `tools/check_changed.R` - Check changed packages only

## Installation (after merge)

```r
remotes::install_github("matsui-lab/cfomics", subdir = "packages/cfomics")
```

## Test plan

- [ ] `R CMD check packages/cfomics` passes
- [ ] `devtools::install("packages/cfomics")` works
- [ ] `tools/check.R` runs successfully
- [ ] Package loads and basic functions work

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## トラブルシューティング

### 問題: NAMESPACE が見つからない

```bash
# packages/cfomics/ で実行
cd packages/cfomics
Rscript -e "devtools::document()"
```

### 問題: man/ ファイルが古い

```bash
cd packages/cfomics
Rscript -e "devtools::document()"
git add man/
git commit -m "docs: regenerate documentation"
```

### 問題: testthat が動かない

```bash
# testthat.R の存在確認
ls packages/cfomics/tests/testthat.R

# なければ作成（Step 4 参照）
```

### 問題: install_github が失敗

`subdir` パラメータの確認：
```r
remotes::install_github("matsui-lab/cfomics", subdir = "packages/cfomics")
```

---

## 次のステップ

PR-1 がマージされたら:

1. **PR-2**: Core registry + as_cf_data generic の実装
2. **CI 更新**: `.github/workflows/` のパス更新（別 PR でも可）

---

*Last updated: 2026-01-05*
