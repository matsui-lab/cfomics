# PR-7: CI Path Filter Setup

## 目的

GitHub Actions の CI を設定し、変更されたパッケージのみをチェックする path filter を導入する。Python パッケージは別 workflow で隔離する。

## 前提条件

- [ ] PR-1 がマージ済み（packages/ 構造）

## ブランチ

```bash
git checkout main
git pull origin main
git checkout -b feature/ci-path-filter
```

---

## 作業手順

### Step 1: ディレクトリ作成

```bash
mkdir -p .github/workflows
```

### Step 2: R-only パッケージ用 workflow

**ファイル**: `.github/workflows/r-cmd-check.yml`

```yaml
name: R-CMD-check

on:
  push:
    branches: [main, develop]
    paths:
      - 'packages/cfomics/**'
      - 'packages/cfomicsSim/**'
      - 'packages/cfomicsBench/**'
      - 'packages/cfomicsAdapters/**'
      - 'packages/cfomicsDiagnose/**'
      - '.github/workflows/r-cmd-check.yml'
  pull_request:
    branches: [main, develop]
    paths:
      - 'packages/cfomics/**'
      - 'packages/cfomicsSim/**'
      - 'packages/cfomicsBench/**'
      - 'packages/cfomicsAdapters/**'
      - 'packages/cfomicsDiagnose/**'
      - '.github/workflows/r-cmd-check.yml'

jobs:
  detect-changes:
    runs-on: ubuntu-latest
    outputs:
      cfomics: ${{ steps.filter.outputs.cfomics }}
      cfomicsSim: ${{ steps.filter.outputs.cfomicsSim }}
      cfomicsBench: ${{ steps.filter.outputs.cfomicsBench }}
      cfomicsAdapters: ${{ steps.filter.outputs.cfomicsAdapters }}
      cfomicsDiagnose: ${{ steps.filter.outputs.cfomicsDiagnose }}
    steps:
      - uses: actions/checkout@v4
      - uses: dorny/paths-filter@v3
        id: filter
        with:
          filters: |
            cfomics:
              - 'packages/cfomics/**'
            cfomicsSim:
              - 'packages/cfomicsSim/**'
            cfomicsBench:
              - 'packages/cfomicsBench/**'
            cfomicsAdapters:
              - 'packages/cfomicsAdapters/**'
            cfomicsDiagnose:
              - 'packages/cfomicsDiagnose/**'

  check-cfomics:
    needs: detect-changes
    if: needs.detect-changes.outputs.cfomics == 'true'
    runs-on: ${{ matrix.config.os }}
    name: cfomics (${{ matrix.config.os }}, R ${{ matrix.config.r }})
    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: ubuntu-latest, r: 'release'}
          - {os: ubuntu-latest, r: 'devel'}
          - {os: macos-latest, r: 'release'}
          - {os: windows-latest, r: 'release'}

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-pandoc@v2

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}
          use-public-rspm: true

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          working-directory: packages/cfomics
          extra-packages: any::rcmdcheck
          needs: check

      - uses: r-lib/actions/check-r-package@v2
        with:
          working-directory: packages/cfomics
          upload-snapshots: true

  check-cfomicsSim:
    needs: [detect-changes, check-cfomics]
    if: |
      always() &&
      needs.detect-changes.outputs.cfomicsSim == 'true' &&
      (needs.check-cfomics.result == 'success' || needs.check-cfomics.result == 'skipped')
    runs-on: ubuntu-latest
    name: cfomicsSim

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          working-directory: packages/cfomicsSim
          extra-packages: any::rcmdcheck
          needs: check

      # Install cfomics first
      - name: Install cfomics
        run: |
          install.packages("remotes")
          remotes::install_local("packages/cfomics", dependencies = TRUE)
        shell: Rscript {0}

      - uses: r-lib/actions/check-r-package@v2
        with:
          working-directory: packages/cfomicsSim

  check-cfomicsBench:
    needs: [detect-changes, check-cfomics, check-cfomicsSim]
    if: |
      always() &&
      needs.detect-changes.outputs.cfomicsBench == 'true' &&
      (needs.check-cfomics.result == 'success' || needs.check-cfomics.result == 'skipped') &&
      (needs.check-cfomicsSim.result == 'success' || needs.check-cfomicsSim.result == 'skipped')
    runs-on: ubuntu-latest
    name: cfomicsBench

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          working-directory: packages/cfomicsBench
          extra-packages: any::rcmdcheck
          needs: check

      # Install dependencies
      - name: Install dependencies
        run: |
          install.packages("remotes")
          remotes::install_local("packages/cfomics", dependencies = TRUE)
          remotes::install_local("packages/cfomicsSim", dependencies = TRUE)
        shell: Rscript {0}

      - uses: r-lib/actions/check-r-package@v2
        with:
          working-directory: packages/cfomicsBench

  check-cfomicsAdapters:
    needs: [detect-changes, check-cfomics]
    if: |
      always() &&
      needs.detect-changes.outputs.cfomicsAdapters == 'true' &&
      (needs.check-cfomics.result == 'success' || needs.check-cfomics.result == 'skipped')
    runs-on: ubuntu-latest
    name: cfomicsAdapters

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true

      # Install Bioconductor
      - name: Setup Bioconductor
        run: |
          install.packages("BiocManager")
          BiocManager::install(version = "3.18")
        shell: Rscript {0}

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          working-directory: packages/cfomicsAdapters
          extra-packages: any::rcmdcheck
          needs: check

      - name: Install cfomics
        run: |
          remotes::install_local("packages/cfomics", dependencies = TRUE)
        shell: Rscript {0}

      - uses: r-lib/actions/check-r-package@v2
        with:
          working-directory: packages/cfomicsAdapters

  check-cfomicsDiagnose:
    needs: [detect-changes, check-cfomics]
    if: |
      always() &&
      needs.detect-changes.outputs.cfomicsDiagnose == 'true' &&
      (needs.check-cfomics.result == 'success' || needs.check-cfomics.result == 'skipped')
    runs-on: ubuntu-latest
    name: cfomicsDiagnose

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          working-directory: packages/cfomicsDiagnose
          extra-packages: any::rcmdcheck
          needs: check

      - name: Install cfomics
        run: |
          remotes::install_local("packages/cfomics", dependencies = TRUE)
        shell: Rscript {0}

      - uses: r-lib/actions/check-r-package@v2
        with:
          working-directory: packages/cfomicsDiagnose
```

### Step 3: Python パッケージ用 workflow（隔離）

**ファイル**: `.github/workflows/r-cmd-check-python.yml`

```yaml
name: R-CMD-check-Python

on:
  push:
    branches: [main, develop]
    paths:
      - 'packages/cfomicsPython/**'
      - '.github/workflows/r-cmd-check-python.yml'
  pull_request:
    branches: [main, develop]
    paths:
      - 'packages/cfomicsPython/**'
      - '.github/workflows/r-cmd-check-python.yml'

jobs:
  check-cfomicsPython:
    runs-on: ${{ matrix.config.os }}
    name: cfomicsPython (${{ matrix.config.os }}, R ${{ matrix.config.r }}, Python ${{ matrix.config.python }})
    strategy:
      fail-fast: false
      matrix:
        config:
          - {os: ubuntu-latest, r: 'release', python: '3.10'}
          - {os: ubuntu-latest, r: 'release', python: '3.11'}
          - {os: macos-latest, r: 'release', python: '3.10'}

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes
      RETICULATE_PYTHON_ENV: cfomics-ci

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-pandoc@v2

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.config.r }}
          use-public-rspm: true

      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.config.python }}

      - name: Create Python virtual environment
        run: |
          python -m venv .venv
          source .venv/bin/activate
          pip install --upgrade pip
          pip install -r packages/cfomicsPython/inst/python/requirements.txt
        shell: bash

      - name: Set RETICULATE_PYTHON
        run: |
          echo "RETICULATE_PYTHON=$(pwd)/.venv/bin/python" >> $GITHUB_ENV
        shell: bash

      # Install cfomics first
      - name: Install cfomics
        run: |
          install.packages("remotes")
          remotes::install_local("packages/cfomics", dependencies = TRUE)
        shell: Rscript {0}

      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          working-directory: packages/cfomicsPython
          extra-packages: any::rcmdcheck
          needs: check

      - uses: r-lib/actions/check-r-package@v2
        with:
          working-directory: packages/cfomicsPython
          upload-snapshots: true

      - name: Test Python integration
        run: |
          library(cfomics)
          library(cfomicsPython)

          # Check environment
          env_check <- cf_check_env()
          print(env_check)

          # Check methods registered
          methods <- cf_methods(include_unavailable = TRUE)
          print(methods)

          # Verify Python methods are registered
          stopifnot("dowhy_gcm" %in% methods$id)
          stopifnot("drlearner" %in% methods$id)
        shell: Rscript {0}
```

### Step 4: Nightly 全パッケージチェック（オプション）

**ファイル**: `.github/workflows/r-cmd-check-nightly.yml`

```yaml
name: R-CMD-check-Nightly

on:
  schedule:
    - cron: '0 3 * * *'  # 毎日 3:00 UTC
  workflow_dispatch:  # 手動実行も可能

jobs:
  check-all-packages:
    runs-on: ubuntu-latest
    name: Full check (all packages)

    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}

    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: 'release'
          use-public-rspm: true

      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'

      - name: Setup Python
        run: |
          python -m venv .venv
          source .venv/bin/activate
          pip install --upgrade pip
          if [ -f packages/cfomicsPython/inst/python/requirements.txt ]; then
            pip install -r packages/cfomicsPython/inst/python/requirements.txt
          fi
          echo "RETICULATE_PYTHON=$(pwd)/.venv/bin/python" >> $GITHUB_ENV
        shell: bash

      - name: Install R dependencies
        run: |
          install.packages(c("remotes", "rcmdcheck", "devtools"))
        shell: Rscript {0}

      - name: Check all packages
        run: |
          source("tools/check.R")
        shell: Rscript {0}

      - name: Upload check results
        if: failure()
        uses: actions/upload-artifact@v4
        with:
          name: check-results
          path: packages/*/*.Rcheck/
```

### Step 5: PR 用のテスト workflow

**ファイル**: `.github/workflows/test-pr.yml`

```yaml
name: PR Tests

on:
  pull_request:
    branches: [main, develop]

jobs:
  lint:
    runs-on: ubuntu-latest
    name: Lint
    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-r@v2

      - name: Install lintr
        run: install.packages("lintr")
        shell: Rscript {0}

      - name: Lint changed R files
        run: |
          files <- system("git diff --name-only origin/main...HEAD | grep '\\.R$'", intern = TRUE)
          if (length(files) > 0) {
            results <- lintr::lint(files)
            print(results)
            if (length(results) > 0) {
              quit(status = 1)
            }
          }
        shell: Rscript {0}

  spell-check:
    runs-on: ubuntu-latest
    name: Spell Check
    steps:
      - uses: actions/checkout@v4

      - uses: r-lib/actions/setup-r@v2

      - name: Install spelling
        run: install.packages("spelling")
        shell: Rscript {0}

      - name: Check spelling
        run: |
          # Check each package that exists
          packages <- list.dirs("packages", recursive = FALSE, full.names = TRUE)
          for (pkg in packages) {
            if (file.exists(file.path(pkg, "DESCRIPTION"))) {
              message("Checking spelling in: ", basename(pkg))
              errors <- spelling::spell_check_package(pkg)
              if (nrow(errors) > 0) {
                print(errors)
              }
            }
          }
        shell: Rscript {0}
```

### Step 6: ステータスバッジの更新

ルート README にバッジを追加:

```markdown
# cfomics Ecosystem

[![R-CMD-check](https://github.com/matsui-lab/cfomics/workflows/R-CMD-check/badge.svg)](https://github.com/matsui-lab/cfomics/actions/workflows/r-cmd-check.yml)
[![R-CMD-check-Python](https://github.com/matsui-lab/cfomics/workflows/R-CMD-check-Python/badge.svg)](https://github.com/matsui-lab/cfomics/actions/workflows/r-cmd-check-python.yml)
```

---

## 検証手順

### Step A: workflow 構文チェック

```bash
# actionlint がインストールされている場合
actionlint .github/workflows/*.yml

# または YAML 構文チェック
python -c "import yaml; yaml.safe_load(open('.github/workflows/r-cmd-check.yml'))"
```

### Step B: ローカルでの path filter テスト

```bash
# cfomics のみ変更した場合
touch packages/cfomics/R/test.R
git diff --name-only

# どの job が実行されるか確認
# cfomics のパスにマッチ → check-cfomics job が実行される
```

### Step C: GitHub Actions での確認

1. PR を作成
2. Actions タブで実行を確認
3. 変更したパッケージのみが check されることを確認

---

## DoD チェックリスト

- [ ] `.github/workflows/r-cmd-check.yml` が作成されている
- [ ] `.github/workflows/r-cmd-check-python.yml` が作成されている
- [ ] path filter が正しく動作する
- [ ] cfomics の変更で cfomicsPython job が実行されない
- [ ] cfomicsPython の変更で Python workflow が実行される
- [ ] 依存関係の順序が正しい（cfomics → cfomicsSim → cfomicsBench）

---

## コミット手順

```bash
git add .github/workflows/
git commit -m "ci: add GitHub Actions with path filtering

- Add r-cmd-check.yml for R-only packages
- Add r-cmd-check-python.yml for Python packages (isolated)
- Add r-cmd-check-nightly.yml for scheduled full checks
- Add test-pr.yml for lint and spell check
- Use path filters to only check changed packages
- Respect package dependencies in job ordering

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## PR 作成

```bash
git push -u origin feature/ci-path-filter

gh pr create --title "ci: add GitHub Actions with path filtering" --body "$(cat <<'EOF'
## Summary

Add comprehensive GitHub Actions CI configuration with path filtering to only check changed packages.

## Workflows

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| r-cmd-check.yml | Push/PR to packages/* | Check R-only packages |
| r-cmd-check-python.yml | Push/PR to cfomicsPython | Check Python packages |
| r-cmd-check-nightly.yml | Daily 3:00 UTC | Full check all packages |
| test-pr.yml | PR | Lint and spell check |

## Features

- **Path filtering**: Only check packages with changes
- **Dependency ordering**: cfomics → cfomicsSim → cfomicsBench
- **Python isolation**: Separate workflow for Python to avoid flaky tests
- **Multi-platform**: Ubuntu, macOS, Windows for core package
- **Multi-R-version**: release and devel for core package

## Test plan

- [ ] Changing cfomics triggers cfomics check
- [ ] Changing cfomicsSim triggers cfomicsSim check (after cfomics)
- [ ] Changing cfomicsPython triggers only Python workflow
- [ ] Changing docs doesn't trigger any check

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## トラブルシューティング

### 問題: path filter が動作しない

```yaml
# dorny/paths-filter の代わりに on.push.paths を使用
on:
  push:
    paths:
      - 'packages/cfomics/**'
```

### 問題: 依存パッケージがインストールできない

```yaml
# モノレポ内のパッケージを先にインストール
- name: Install local dependencies
  run: |
    install.packages("remotes")
    remotes::install_local("packages/cfomics", dependencies = TRUE)
  shell: Rscript {0}
```

### 問題: Python 環境が認識されない

```yaml
# RETICULATE_PYTHON を明示的に設定
- name: Set RETICULATE_PYTHON
  run: |
    echo "RETICULATE_PYTHON=$(which python)" >> $GITHUB_ENV
```

---

*Last updated: 2026-01-05*
