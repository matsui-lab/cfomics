# PR-4: Python Methods Migration

## 目的

`packages/cfomics/` から Python 関連のコードを `packages/cfomicsPython/` に移植し、Core を Python 非依存にする。

## 前提条件

- [ ] PR-1 がマージ済み（packages/ 構造）
- [ ] PR-2 がマージ済み（registry システム）
- [ ] PR-3 がマージ済み（cfomicsPython skeleton）

## ブランチ

```bash
git checkout main
git pull origin main
git checkout -b feature/python-migration
```

---

## 作業手順

### Step 1: 移動対象ファイルの確認

**移動元（packages/cfomics/）**:
- `R/python_env.R`
- `R/methods_dowhy_gcm.R`
- `R/methods_drlearner.R`
- `R/methods_ganite.R`
- `R/methods_cavae.R`
- `R/env_setup.R`
- `R/env_diagnostics.R`
- `inst/python/gcm_dowhy.py`
- `inst/python/ganite.py`
- `inst/python/visualization.py`
- `inst/python/requirements.txt`
- `inst/python/test_gcm_wrapper.py`

**移動先（packages/cfomicsPython/）**:
- `R/` ディレクトリ
- `inst/python/` ディレクトリ

### Step 2: R ファイルの移動

```bash
# Python environment management
git mv packages/cfomics/R/python_env.R packages/cfomicsPython/R/python_env.R
git mv packages/cfomics/R/env_setup.R packages/cfomicsPython/R/env_setup.R
git mv packages/cfomics/R/env_diagnostics.R packages/cfomicsPython/R/env_diagnostics.R

# Method implementations
git mv packages/cfomics/R/methods_dowhy_gcm.R packages/cfomicsPython/R/methods_dowhy_gcm.R
git mv packages/cfomics/R/methods_drlearner.R packages/cfomicsPython/R/methods_drlearner.R
git mv packages/cfomics/R/methods_ganite.R packages/cfomicsPython/R/methods_ganite.R
git mv packages/cfomics/R/methods_cavae.R packages/cfomicsPython/R/methods_cavae.R
```

### Step 3: Python ファイルの移動

```bash
git mv packages/cfomics/inst/python/* packages/cfomicsPython/inst/python/
```

### Step 4: スタブファイルの削除

PR-3 で作成したスタブファイルを削除（実際の実装に置き換え）:

```bash
rm packages/cfomicsPython/R/methods_stub.R
```

### Step 5: cfomicsPython-package.R の更新

**ファイル**: `packages/cfomicsPython/R/cfomicsPython-package.R`

スタブの `.register_python_methods()` を実際の関数参照に更新:

```r
#' @keywords internal
"_PACKAGE"

#' @import cfomics
#' @importFrom reticulate py_available py_module_available use_condaenv
#' @importFrom reticulate py_config import source_python
#' @importFrom rlang abort warn
#' @importFrom cli cli_alert_success cli_alert_warning cli_alert_info
NULL

.onLoad <- function(libname, pkgname) {
  # Register Python methods with cfomics
  # NOTE: This only registers methods, does NOT initialize Python
  .register_python_methods()

  invisible()
}

#' Register Python-based methods with cfomics registry
#' @keywords internal
#' @noRd
.register_python_methods <- function() {
  # DoWhy GCM
  cfomics::cf_register_method(
    id = "dowhy_gcm",
    fit_fun = cf_fit_dowhy_gcm,
    predict_fun = NULL,
    requires_python = TRUE,
    package = "cfomicsPython",
    description = "DoWhy GCM (Graphical Causal Model)"
  )

  # DRLearner (EconML)
  cfomics::cf_register_method(
    id = "drlearner",
    fit_fun = cf_fit_drlearner,
    predict_fun = predict_cf_drlearner,
    requires_python = TRUE,
    package = "cfomicsPython",
    description = "Double Robust Learner (EconML)"
  )

  # GANITE
  cfomics::cf_register_method(
    id = "ganite",
    fit_fun = cf_fit_ganite,
    predict_fun = predict_cf_ganite,
    requires_python = TRUE,
    package = "cfomicsPython",
    description = "GANITE (Generative Adversarial Nets for ITE)"
  )

  # CAVAE
  cfomics::cf_register_method(
    id = "cavae",
    fit_fun = cf_fit_cavae,
    predict_fun = predict_cf_cavae,
    requires_python = TRUE,
    package = "cfomicsPython",
    description = "Causal Autoencoder VAE"
  )

  invisible()
}
```

### Step 6: Core から Python 依存を削除

**packages/cfomics/R/aaa-onload.R** から Python 関連コードを削除（必要な場合）

**packages/cfomics/DESCRIPTION** の Suggests から reticulate を削除（または維持して optional に）:

```
Suggests:
    testthat (>= 3.0.0),
    ...
    # reticulate は cfomicsPython に移動
```

### Step 7: Core のテストから Python テストを移動

**移動対象テスト**:
- `test-python-env.R`
- `test-python-methods-structure.R`
- `test-cavae.R`
- `test-drlearner.R`
- `test-ganite.R`
- `test-serialization-python-policy.R`
- `test-env-setup.R`
- `helper-reticulate.R`

```bash
git mv packages/cfomics/tests/testthat/test-python-env.R packages/cfomicsPython/tests/testthat/
git mv packages/cfomics/tests/testthat/test-python-methods-structure.R packages/cfomicsPython/tests/testthat/
git mv packages/cfomics/tests/testthat/test-cavae.R packages/cfomicsPython/tests/testthat/
git mv packages/cfomics/tests/testthat/test-drlearner.R packages/cfomicsPython/tests/testthat/
git mv packages/cfomics/tests/testthat/test-ganite.R packages/cfomicsPython/tests/testthat/
git mv packages/cfomics/tests/testthat/test-serialization-python-policy.R packages/cfomicsPython/tests/testthat/
git mv packages/cfomics/tests/testthat/test-env-setup.R packages/cfomicsPython/tests/testthat/
git mv packages/cfomics/tests/testthat/helper-reticulate.R packages/cfomicsPython/tests/testthat/
```

### Step 8: cfomicsPython DESCRIPTION の更新

**ファイル**: `packages/cfomicsPython/DESCRIPTION`

```
Package: cfomicsPython
Title: Python Backends for cfomics
Version: 0.1.0
Authors@R:
    person("Yusuke", "Matsui", , "matsui@example.com", role = c("aut", "cre"),
           comment = c(ORCID = "0000-0000-0000-0000"))
Description: Provides Python-based causal inference methods for the cfomics
    ecosystem. Includes backends for DoWhy, EconML, GANITE, and CAVAE.
    Requires Python and relevant Python packages to be installed.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.1
Depends:
    R (>= 4.2.0)
Imports:
    cfomics,
    reticulate (>= 1.28),
    rlang,
    cli,
    jsonlite
Suggests:
    testthat (>= 3.0.0),
    withr
Config/testthat/edition: 3
```

### Step 9: cfomicsPython NAMESPACE の更新

```bash
cd packages/cfomicsPython
Rscript -e "devtools::document()"
```

### Step 10: Core の互換性ラッパー（オプション）

Core に Python 関連関数の互換性ラッパーを残す場合:

**ファイル**: `packages/cfomics/R/python_compat.R`

```r
#' Check Python availability
#'
#' @param method Optional method name to check specific requirements
#' @return Logical
#' @export
cf_has_python <- function(method = NULL) {
  # Check if cfomicsPython is available
  if (!requireNamespace("cfomicsPython", quietly = TRUE)) {
    return(FALSE)
  }

  # Delegate to cfomicsPython
  cfomicsPython::cf_has_python(method)
}

#' Require Python for a method
#'
#' @param method Method name
#' @return Invisible TRUE, or error
#' @export
cf_require_python <- function(method) {
  if (!requireNamespace("cfomicsPython", quietly = TRUE)) {
    rlang::abort(
      paste0(
        "Python methods require the 'cfomicsPython' package.\n",
        "Install with: remotes::install_github('matsui-lab/cfomics', subdir='packages/cfomicsPython')"
      ),
      class = "cfomics_missing_package_error"
    )
  }

  cfomicsPython::cf_require_python(method)
}
```

### Step 11: man/ ファイルの移動

Python 関連のドキュメントを移動:

```bash
# Python 関連の man ファイルを確認
ls packages/cfomics/man/ | grep -E "(python|dowhy|ganite|drlearner|cavae|env)"

# 該当ファイルを移動
git mv packages/cfomics/man/cf_check_env.Rd packages/cfomicsPython/man/ 2>/dev/null || true
git mv packages/cfomics/man/cf_install_python_env.Rd packages/cfomicsPython/man/ 2>/dev/null || true
git mv packages/cfomics/man/cf_use_python_env.Rd packages/cfomicsPython/man/ 2>/dev/null || true
git mv packages/cfomics/man/cf_list_python_envs.Rd packages/cfomicsPython/man/ 2>/dev/null || true
# ... その他の Python 関連 man ファイル
```

### Step 12: ドキュメント再生成

```bash
cd packages/cfomics
Rscript -e "devtools::document()"

cd ../cfomicsPython
Rscript -e "devtools::document()"
```

---

## 検証手順

### Step A: Core の Python 非依存確認

```bash
cd packages/cfomics

# reticulate なしでビルド・チェック
R CMD build .
R CMD check cfomics_*.tar.gz --as-cran
```

### Step B: Core テストが Python なしで通ること

```r
# Python 環境変数をクリアしてテスト
Sys.unsetenv("RETICULATE_PYTHON")
devtools::test("packages/cfomics")
```

### Step C: cfomicsPython のビルド・チェック

```bash
cd packages/cfomicsPython
R CMD build .
R CMD check cfomicsPython_*.tar.gz --as-cran
```

### Step D: 統合テスト

```r
# 両パッケージをインストール
devtools::install("packages/cfomics")
devtools::install("packages/cfomicsPython")

library(cfomics)
library(cfomicsPython)

# Python methods が見えることを確認
cf_methods(include_unavailable = TRUE)

# Python があれば実行テスト
if (reticulate::py_available()) {
  df <- data.frame(Y = rnorm(100), T = rbinom(100, 1, 0.5), X1 = rnorm(100))
  # result <- cf_fit(Y ~ T | X1, data = df, method = "drlearner")
}
```

---

## DoD チェックリスト

- [ ] Python 関連の R ファイルが cfomicsPython に移動済み
- [ ] Python スクリプトが cfomicsPython/inst/python/ に移動済み
- [ ] Python テストが cfomicsPython に移動済み
- [ ] Core (cfomics) が Python 非依存で `R CMD check` 通過
- [ ] cfomicsPython が `R CMD check` 通過
- [ ] `cf_methods(include_unavailable = TRUE)` で Python methods が見える
- [ ] Core のみインストール時、Python method 呼び出しで適切なエラー

---

## コミット手順

### Commit 1: R ファイル移動

```bash
git add -A
git commit -m "refactor: move Python R files to cfomicsPython

Move Python environment management and method implementations
from cfomics to cfomicsPython package.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 2: Python スクリプト移動

```bash
git add -A
git commit -m "refactor: move Python scripts to cfomicsPython

Move inst/python/* to cfomicsPython package.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 3: テスト移動

```bash
git add -A
git commit -m "refactor: move Python tests to cfomicsPython

Move Python-dependent tests from cfomics to cfomicsPython.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 4: Core クリーンアップ

```bash
git add -A
git commit -m "refactor: remove Python dependencies from cfomics core

- Remove Python-specific imports
- Add compatibility wrappers for cf_has_python/cf_require_python
- Update NAMESPACE

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Commit 5: ドキュメント更新

```bash
git add -A
git commit -m "docs: regenerate documentation after Python migration

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## PR 作成

```bash
git push -u origin feature/python-migration

gh pr create --title "refactor: migrate Python code to cfomicsPython package" --body "$(cat <<'EOF'
## Summary

Move all Python-related code from `cfomics` to `cfomicsPython` package, making the core package Python-independent.

## Changes

### Moved to cfomicsPython
- `R/python_env.R`, `R/env_setup.R`, `R/env_diagnostics.R`
- `R/methods_dowhy_gcm.R`, `R/methods_drlearner.R`, `R/methods_ganite.R`, `R/methods_cavae.R`
- `inst/python/*`
- Python-related tests

### cfomics Core
- Now Python-independent
- Compatibility wrappers for `cf_has_python()`, `cf_require_python()`
- Helpful error messages when Python methods requested without cfomicsPython

## Benefits

- Core package is lighter and more stable for CRAN/Bioc
- Python dependencies isolated in separate package
- Easier CI (Core tests don't need Python setup)
- Users without Python can still use R-native methods

## Test plan

- [ ] `R CMD check packages/cfomics` passes without Python
- [ ] `R CMD check packages/cfomicsPython` passes
- [ ] `cf_methods()` shows R-native methods
- [ ] `cf_methods(include_unavailable = TRUE)` shows all methods
- [ ] Error message helpful when calling Python method without cfomicsPython

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## トラブルシューティング

### 問題: Core テストが Python を要求する

```r
# テストファイルを確認
grep -r "reticulate" packages/cfomics/tests/

# Python 依存テストを cfomicsPython に移動し忘れがないか確認
```

### 問題: NAMESPACE エラー

```bash
cd packages/cfomics
Rscript -e "devtools::document()"

cd ../cfomicsPython
Rscript -e "devtools::document()"
```

### 問題: 関数が見つからない

cfomicsPython で cfomics の関数を使う場合、`cfomics::` プレフィックスを確認:

```r
# Good
cfomics::cf_register_method(...)
cfomics::validate_cfomics_result(...)

# Bad (NAMESPACE に import されていない場合)
cf_register_method(...)
```

---

*Last updated: 2026-01-05*
