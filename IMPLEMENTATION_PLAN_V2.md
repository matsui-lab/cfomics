# cfomics Ecosystem v2.0 実装計画書

**Monorepo-first / Split-ready Implementation Plan**

---

## 0. 目的（本計画で達成する Claims）

本実装計画は、cfomics を以下の状態にする：

1. **Split-ready**: `packages/` 配下に独立パッケージとして配置し、将来の別リポジトリ化が「ディレクトリ単位の切り出し」で済む
2. **Contract-driven**: `cf_data`（入力）と `cfomics_result`（出力）が v2 として固定され、全パッケージがそれに従う
3. **Registry-based**: method 拡張は registry 登録で実現（Core に switch を増やさない）
4. **Python-independent Core**: `packages/cfomics` 単体で R-only methods が動作し、Python 非依存で CI が通る
5. **Reproducible Benchmarks**: Sim/Bench 分離により、論文用ベンチマークの再現性を担保

---

## 1. スコープ

### 1.1 対象（v2 で移行・実装するもの）

| パッケージ | 状態 | 優先度 |
|-----------|------|--------|
| cfomics (Core) | 既存を移行・リファクタ | 最高 |
| cfomicsPython | 既存コードを分離 | 高 |
| cfomicsSim | 新規作成（既存 benchmark から抽出） | 高 |
| cfomicsBench | 新規作成（既存 benchmark から抽出） | 高 |
| cfomicsAdapters | 既存 bioc_integration を分離 | 中 |
| cfomicsDiagnose | 既存 diagnostics を分離 | 中 |
| cfomicsSensitivity | 新規作成 | 低（後続） |
| cfomicsReport | 新規作成 | 低（後続） |

### 1.2 実装モード

- **Mode A（CI 最小）**: R-only tests のみ、Python/公開データ不要
- **Mode B（ローカル/論文再現）**: Python 環境 + シミュレーション + ベンチマーク完全実行

---

## 2. 前提条件（現状の把握）

### 2.1 現行ファイル構成

現 `cfomics` に同居しているコンポーネント：

```
R/
├── aaa-onload.R           # Load hook（副作用なし方針）
├── cf_fit.R               # Core: 統一 API
├── cf_predict.R           # Core: 予測
├── config.R               # Core: 設定管理
├── formula_utils.R        # Core: 数式パース
├── result_contract.R      # Core: 結果契約
├── methods_traditional.R  # Core: R-only methods (GRF/IPW/G-formula)
├── python_env.R           # Python: 環境管理
├── methods_dowhy_gcm.R    # Python: DoWhy
├── methods_drlearner.R    # Python: EconML
├── methods_ganite.R       # Python: GANITE
├── methods_cavae.R        # Python: CAVAE
├── bioc_integration.R     # Adapters: Bioc 変換
├── benchmark_*.R          # Bench: ベンチマーク
├── diagnostics.R          # Diagnose: 診断
└── visualization.R        # Viz: 可視化
```

### 2.2 既知の課題

- GitHub ライセンス表示が "Unknown, MIT licenses found"
- Python 依存コードが Core に混在
- method 追加に `cf_fit()` の switch 修正が必要

---

## 3. リポジトリ構成（v2 目標）

```
cfomics/
├── packages/
│   ├── cfomics/
│   │   ├── R/
│   │   │   ├── aaa-onload.R
│   │   │   ├── cf_fit.R
│   │   │   ├── cf_predict.R
│   │   │   ├── config.R
│   │   │   ├── formula_utils.R
│   │   │   ├── result_contract.R
│   │   │   ├── methods_traditional.R
│   │   │   ├── registry.R          # NEW: method registry
│   │   │   └── as_cf_data.R        # NEW: generic + data.frame method
│   │   ├── inst/
│   │   ├── man/
│   │   ├── tests/testthat/
│   │   ├── vignettes/
│   │   ├── DESCRIPTION
│   │   ├── NAMESPACE
│   │   ├── LICENSE
│   │   ├── NEWS.md
│   │   └── README.md
│   ├── cfomicsPython/
│   │   ├── R/
│   │   │   ├── python_env.R
│   │   │   ├── methods_dowhy_gcm.R
│   │   │   ├── methods_drlearner.R
│   │   │   ├── methods_ganite.R
│   │   │   ├── methods_cavae.R
│   │   │   └── register_methods.R  # registry 登録
│   │   ├── inst/python/
│   │   ├── tests/testthat/
│   │   ├── DESCRIPTION
│   │   └── NAMESPACE
│   ├── cfomicsSim/
│   │   ├── R/
│   │   │   ├── simulate.R
│   │   │   ├── scenarios.R
│   │   │   └── truth.R
│   │   ├── inst/extdata/scenarios/
│   │   ├── tests/testthat/
│   │   ├── DESCRIPTION
│   │   └── NAMESPACE
│   ├── cfomicsBench/
│   │   ├── R/
│   │   │   ├── benchmark_run.R
│   │   │   ├── benchmark_summarize.R
│   │   │   ├── benchmark_plot.R
│   │   │   └── metrics.R
│   │   ├── tests/testthat/
│   │   ├── DESCRIPTION
│   │   └── NAMESPACE
│   ├── cfomicsAdapters/
│   │   ├── R/
│   │   │   └── as_cf_data_bioc.R
│   │   ├── tests/testthat/
│   │   ├── DESCRIPTION
│   │   └── NAMESPACE
│   └── cfomicsDiagnose/
│       ├── R/
│       │   └── diagnose.R
│       ├── tests/testthat/
│       ├── DESCRIPTION
│       └── NAMESPACE
├── tools/
│   ├── check.R              # 全パッケージチェック
│   ├── check_changed.R      # 変更パッケージのみチェック
│   └── check_package.R      # 単一パッケージチェック
├── .github/workflows/
│   ├── r-cmd-check.yml      # R-only packages
│   └── r-cmd-check-python.yml  # Python packages（隔離）
├── README.md                # エコシステム概要
├── CLAUDE.md                # AI 開発ガイド（更新）
├── ECOSYSTEM_SPEC_V2.md     # 仕様書
├── IMPLEMENTATION_PLAN_V2.md # 本計画書
└── LICENSE.md               # MIT（GitHub 認識用）
```

---

## 4. 設計原則（実装上の制約）

### 4.1 境界侵食の禁止

```r
# 禁止: 他パッケージの内部関数呼び出し
cfomicsPython:::internal_function()  # NG

# 許可: export された関数のみ
cfomics::cf_register_method()        # OK
```

### 4.2 依存方向の固定

```
cfomics (Core) ← 他パッケージから依存される
    ↑
    ├── cfomicsPython (Imports: cfomics)
    ├── cfomicsAdapters (Imports: cfomics)
    ├── cfomicsDiagnose (Imports: cfomics)
    └── cfomicsBench (Imports: cfomics, cfomicsSim)
            ↑
            └── cfomicsSim (Imports: cfomics)
```

### 4.3 Python 初期化の遅延

```r
# .onLoad() での初期化は禁止
.onLoad <- function(libname, pkgname) {

  # Python 初期化 NG
  # reticulate::py_config()  # 禁止

  # registry 登録のみ OK
  cfomics::cf_register_method("dowhy_gcm", ...)
}

# 実行時に初期化
cf_fit_dowhy_gcm <- function(...) {
  cf_require_python("dowhy_gcm")  # ここで初期化
  ...
}
```

### 4.4 Contract 準拠の強制

```r
# 全 method は返却前に検証必須
cf_fit_grf <- function(...) {

  result <- ...

  validate_cfomics_result(result)  # 必須

  result
}
```

---

## 5. コア API 仕様（Phase 2 で実装）

### 5.1 Method Registry

**ファイル**: `packages/cfomics/R/registry.R`

```r
#' @title Method Registry
#' @description cfomics method の登録・取得を管理
#' @name registry

# 内部環境（registry storage）
.cfomics_registry <- new.env(parent = emptyenv())
.cfomics_registry$methods <- list()

#' Register a causal inference method
#' @param id Method identifier (e.g., "grf", "dowhy_gcm")
#' @param fit_fun Fitting function
#' @param predict_fun Prediction function (optional)
#' @param requires_python Logical, whether Python is required
#' @param package Package name that provides this method
#' @param description Human-readable description
#' @export
cf_register_method <- function(id,
                                fit_fun,
                                predict_fun = NULL,
                                requires_python = FALSE,
                                package = NULL,
                                description = NULL) {

  stopifnot(is.character(id), length(id) == 1)
  stopifnot(is.function(fit_fun))

  .cfomics_registry$methods[[id]] <- list(
    id = id,
    fit_fun = fit_fun,
    predict_fun = predict_fun,
    requires_python = requires_python,
    package = package %||% "unknown",
    description = description %||% ""
  )
  invisible(TRUE)
}

#' List available methods
#' @param include_unavailable Include methods that require unavailable dependencies
#' @return data.frame with method information
#' @export
cf_methods <- function(include_unavailable = FALSE) {


  methods <- .cfomics_registry$methods
  if (length(methods) == 0) {
    return(data.frame(
      id = character(),
      package = character(),
      requires_python = logical(),
      available = logical(),
      description = character(),
      stringsAsFactors = FALSE
    ))
  }

  df <- do.call(rbind, lapply(methods, function(m) {
    available <- if (m$requires_python) {
      cf_has_python(m$id)
    } else {
      TRUE
    }
    data.frame(
      id = m$id,
      package = m$package,
      requires_python = m$requires_python,
      available = available,
      description = m$description,
      stringsAsFactors = FALSE
    )
  }))

  if (!include_unavailable) {
    df <- df[df$available, ]
  }

  df
}

#' Get method by ID
#' @param id Method identifier
#' @return Method list or NULL
#' @keywords internal
cf_get_method <- function(id) {
  .cfomics_registry$methods[[id]]
}
```

**受入基準**:
- `cf_register_method()` で method 登録可能
- `cf_methods()` で一覧取得可能
- 未登録 method へのアクセスで適切なエラー

### 5.2 as_cf_data generic

**ファイル**: `packages/cfomics/R/as_cf_data.R`

```r
#' @title Convert to cf_data
#' @description Generic function to convert various data types to cf_data
#' @param x Input data
#' @param outcome_name Name of outcome variable
#' @param treatment_name Name of treatment variable
#' @param covariate_names Names of covariate variables
#' @param ... Additional arguments
#' @return cf_data object
#' @export
as_cf_data <- function(x, outcome_name, treatment_name, covariate_names, ...) {

  UseMethod("as_cf_data")
}

#' @rdname as_cf_data
#' @export
as_cf_data.data.frame <- function(x,
                                   outcome_name,
                                   treatment_name,
                                   covariate_names,
                                   ...) {
  # 検証
 stopifnot(outcome_name %in% names(x))
  stopifnot(treatment_name %in% names(x))
  stopifnot(all(covariate_names %in% names(x)))

  # 構築
  structure(
    list(
      data = x,
      outcome_name = outcome_name,
      treatment_name = treatment_name,
      covariate_names = covariate_names,
      meta = list(
        n = nrow(x),
        p = length(covariate_names),
        types = vapply(x[, c(outcome_name, treatment_name, covariate_names)],
                       class, character(1)),
        created_at = Sys.time()
      )
    ),
    class = "cf_data"
  )
}

#' @export
print.cf_data <- function(x, ...) {
  cli::cli_h1("cf_data object")
  cli::cli_text("Samples: {x$meta$n}")
  cli::cli_text("Covariates: {x$meta$p}")
  cli::cli_text("Outcome: {x$outcome_name}")
  cli::cli_text("Treatment: {x$treatment_name}")
  invisible(x)
}

#' Validate cf_data object
#' @param x Object to validate
#' @return TRUE if valid, error otherwise
#' @export
validate_cf_data <- function(x) {
  if (!inherits(x, "cf_data")) {
    rlang::abort("x must be a cf_data object", class = "cfomics_validation_error")
  }

  required <- c("data", "outcome_name", "treatment_name", "covariate_names", "meta")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    rlang::abort(
      paste("Missing required fields:", paste(missing, collapse = ", ")),
      class = "cfomics_validation_error"
    )
  }

  TRUE
}
```

**受入基準**:
- `as_cf_data.data.frame()` が動作
- `print.cf_data()` で情報表示
- `validate_cf_data()` で検証可能

---

## 6. テスト仕様

### 6.1 Core Tests（packages/cfomics/tests/）

| テストファイル | 内容 | 実行条件 |
|--------------|------|---------|
| test-registry.R | registry 登録・取得 | 常時 |
| test-cf_data.R | cf_data 生成・検証 | 常時 |
| test-result_contract.R | cfomics_result 検証 | 常時 |
| test-traditional.R | R-only methods | 常時 |
| test-cf_fit.R | 統一 API | 常時 |

### 6.2 Python Tests（packages/cfomicsPython/tests/）

| テストファイル | 内容 | 実行条件 |
|--------------|------|---------|
| test-python-env.R | 環境管理 | Python 利用可能時 |
| test-python-methods.R | Python methods | Python + modules 利用可能時 |
| test-register.R | registry 登録確認 | 常時（構造のみ） |

### 6.3 Skip ルール

```r
# Python テストのスキップパターン
skip_if_not(reticulate::py_available(initialize = FALSE))
skip_if_not(reticulate::py_module_available("dowhy"))
```

---

## 7. 作業計画（Phase 別・受入基準付き）

### Phase 0: 移行準備（ベースライン固定）

**目的**: 移行作業中に「壊したのか、元からか」を曖昧にしない

**タスク**

| # | タスク | 担当 | DoD |
|---|-------|------|-----|
| 0.1 | ベースラインタグ付け | Human | `cfomics-monolith-v0.4.0` タグ作成 |
| 0.2 | 現状 check 確認 | AI | `R CMD check` が通る |
| 0.3 | ライセンス方針確定 | Human | MIT で統一決定 |

**受入基準**
- [ ] main が安定（R-only で check が通る）
- [ ] ベースラインをタグで再現できる

---

### Phase 1: packages/ モノレポ化

**目的**: 以降の分割を「ディレクトリ単位の移動」で済ませる基盤

**タスク**

| # | タスク | 詳細 |
|---|-------|------|
| 1.1 | packages/ ディレクトリ作成 | `mkdir packages` |
| 1.2 | 既存パッケージ移動 | `git mv R inst man tests vignettes data-raw DESCRIPTION NAMESPACE LICENSE LICENSE.md NEWS.md packages/cfomics/` |
| 1.3 | ルート README 更新 | エコシステム概要に置き換え |
| 1.4 | パッケージ README 作成 | `packages/cfomics/README.md` |
| 1.5 | インストール導線更新 | `subdir="packages/cfomics"` |
| 1.6 | ルート tools 作成 | `tools/check.R`, `tools/check_changed.R` |
| 1.7 | .Rbuildignore 更新 | ルートレベルの除外設定 |
| 1.8 | CLAUDE.md 更新 | packages/ 構造に対応 |

**受入基準**
- [ ] `packages/cfomics` 単体で `R CMD check` が通る
- [ ] `remotes::install_github(..., subdir="packages/cfomics")` が動く
- [ ] ルート README にインストール手順記載

---

### Phase 2: Core 境界固定

**目的**: 越境ポイントを先に整備し、後の分割をファイル移動で済ませる

**タスク**

| # | タスク | 詳細 |
|---|-------|------|
| 2.1 | registry.R 作成 | `cf_register_method()`, `cf_methods()`, `cf_get_method()` |
| 2.2 | as_cf_data.R 作成 | generic + data.frame method |
| 2.3 | cf_fit.R リファクタ | registry-based dispatch |
| 2.4 | R-only methods 登録 | `.onLoad()` で grf/ipw/gformula 登録 |
| 2.5 | test-registry.R 作成 | registry テスト |
| 2.6 | test-cf_data.R 作成 | cf_data テスト |
| 2.7 | NAMESPACE 更新 | 新規 export 追加 |

**受入基準**
- [ ] `cf_fit(..., method="grf"/"ipw"/"gformula")` が従来どおり動作
- [ ] `cf_methods()` が R-only 3 種を返す
- [ ] Python なしで Core テスト全パス

---

### Phase 3: 新パッケージ骨組み作成

**目的**: AI 駆動でも壊れにくいよう、空パッケージを先に作って check を通す

**タスク（各パッケージ共通）**

| # | パッケージ | 最小 export |
|---|-----------|------------|
| 3.1 | cfomicsPython | `cf_install_python_env()`, `cf_use_python_env()` (stub) |
| 3.2 | cfomicsSim | `cf_simulate()`, `cf_truth()` (stub) |
| 3.3 | cfomicsBench | `cf_benchmark_run()` (stub) |
| 3.4 | cfomicsAdapters | `as_cf_data.SummarizedExperiment()` (stub) |
| 3.5 | cfomicsDiagnose | `cf_diagnose()` (stub) |

**各パッケージの必須ファイル**

```
packages/<pkg>/
├── R/
│   └── <pkg>-package.R  # package doc + stub functions
├── tests/
│   └── testthat/
│       ├── testthat.R
│       └── test-stub.R
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
└── NEWS.md
```

**受入基準**
- [ ] 各パッケージが単独で `R CMD check` を通る
- [ ] 依存が片方向（Core に依存、Core は依存しない）

---

### Phase 4: Python 分離

**目的**: Core の安定性を担保しつつ、Python 拡張を分離

**移植対象**

| 移動元 (packages/cfomics/) | 移動先 (packages/cfomicsPython/) |
|---------------------------|----------------------------------|
| R/python_env.R | R/python_env.R |
| R/methods_dowhy_gcm.R | R/methods_dowhy_gcm.R |
| R/methods_drlearner.R | R/methods_drlearner.R |
| R/methods_ganite.R | R/methods_ganite.R |
| R/methods_cavae.R | R/methods_cavae.R |
| inst/python/* | inst/python/* |

**タスク**

| # | タスク | 詳細 |
|---|-------|------|
| 4.1 | ファイル移動 | `git mv` で移動 |
| 4.2 | register_methods.R 作成 | `.onLoad()` で registry 登録 |
| 4.3 | DESCRIPTION 更新 | Imports: cfomics 追加 |
| 4.4 | NAMESPACE 更新 | export/import 調整 |
| 4.5 | テスト移動 | Python テストを cfomicsPython へ |
| 4.6 | Core テスト修正 | Python 依存テストを削除/スキップ化 |
| 4.7 | cf_has_python 調整 | Core 側は cfomicsPython への参照に |

**受入基準**
- [ ] `packages/cfomics` 単体で Python 不要、check 通過
- [ ] `packages/cfomicsPython` 追加で `cf_methods()` に Python methods 出現
- [ ] Python 未設定時に明確なエラーメッセージ

---

### Phase 5: Sim / Bench 分離

**目的**: B1 論文の核（DGP + ベンチハーネス）を Core から分離

#### Phase 5.1: cfomicsSim

**移植/新規作成**

| ファイル | 内容 |
|---------|------|
| R/simulate.R | `cf_simulate()` 実装 |
| R/scenarios.R | `cf_scenarios()`, `cf_scenario_get()` |
| R/truth.R | `cf_truth()` |
| inst/extdata/scenarios/*.yml | シナリオ定義 |

**最小シナリオ（v2）**

```yaml
# inst/extdata/scenarios/S00_linear.yml
id: S00_linear
name: Linear DGP
description: Simple linear treatment effect
params:
  n_default: 1000
  p_default: 10
  ate_true: 2.0
  treatment_effect_type: constant
mechanisms:
  confounding: true
  heterogeneity: false
```

**受入基準**
- [ ] `cf_simulate("S00_linear", n = 100)` が動作
- [ ] seed 固定で同一結果
- [ ] `cf_truth()` で ATE/ITE/y0/y1 取得可能

#### Phase 5.2: cfomicsBench

**移植/新規作成**

| ファイル | 内容 |
|---------|------|
| R/benchmark_run.R | `cf_benchmark_run()` |
| R/benchmark_summarize.R | `cf_benchmark_summarize()` |
| R/benchmark_plot.R | `cf_benchmark_plot()` |
| R/metrics.R | 評価指標（RMSE, MAE, coverage 等） |

**受入基準**
- [ ] `cf_benchmark_run(scenarios, methods, n_reps)` が動作
- [ ] 結果に provenance（seed, versions）含む
- [ ] スモークテストが通る（1 scenario × 1 method × 1 rep）

---

### Phase 6: Adapters / Diagnose 分離

**目的**: S1/A1 論文を加速する運用層を独立化

#### Phase 6.1: cfomicsAdapters

**移植対象**: `R/bioc_integration.R`

**実装**

```r
#' @rdname as_cf_data
#' @export
as_cf_data.SummarizedExperiment <- function(x,
                                             outcome_name,
                                             treatment_name,
                                             covariate_names,
                                             assay_name = 1,
                                             feature_select = "variance",
                                             n_features = 100,
                                             ...) {
  # Implementation
}
```

**受入基準**
- [ ] Adapters 導入で Bioc オブジェクト対応
- [ ] Adapters 未導入でも Core は動作

#### Phase 6.2: cfomicsDiagnose

**移植対象**: `R/diagnostics.R`

**実装**

```r
#' @title Diagnose causal inference assumptions
#' @param x cf_data or cfomics_result object
#' @return cf_diagnostics object
#' @export
cf_diagnose <- function(x, ...) {
  UseMethod("cf_diagnose")
}
```

**受入基準**
- [ ] `cf_diagnose()` が動作
- [ ] 診断結果を `fit$diagnostics` に付加可能

---

## 8. CI/開発ツール整備

### 8.1 GitHub Actions

**r-cmd-check.yml**（R-only packages）

```yaml
name: R-CMD-check

on:
  push:
    paths:
      - 'packages/cfomics/**'
      - 'packages/cfomicsSim/**'
      - 'packages/cfomicsBench/**'
      - 'packages/cfomicsAdapters/**'
      - 'packages/cfomicsDiagnose/**'
  pull_request:
    paths:
      - 'packages/cfomics/**'
      # ... same paths

jobs:
  check:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        package: [cfomics, cfomicsSim, cfomicsBench]
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
      - name: Check package
        run: |
          cd packages/${{ matrix.package }}
          R CMD build .
          R CMD check *.tar.gz --as-cran
```

**r-cmd-check-python.yml**（Python packages、隔離）

```yaml
name: R-CMD-check-Python

on:
  push:
    paths:
      - 'packages/cfomicsPython/**'
  pull_request:
    paths:
      - 'packages/cfomicsPython/**'

jobs:
  check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
      - uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - name: Install Python deps
        run: pip install -r packages/cfomicsPython/inst/python/requirements.txt
      # ... rest of check
```

### 8.2 ルート tools

**tools/check.R**

```r
#!/usr/bin/env Rscript
# Check all packages in packages/

packages <- list.dirs("packages", recursive = FALSE, full.names = TRUE)

for (pkg in packages) {
  message("Checking: ", basename(pkg))
  result <- devtools::check(pkg, quiet = TRUE)
  if (length(result$errors) > 0 || length(result$warnings) > 0) {
    stop("Check failed for: ", basename(pkg))
  }
}

message("All packages passed!")
```

**tools/check_changed.R**

```r
#!/usr/bin/env Rscript
# Check only changed packages (based on git diff)

changed_files <- system("git diff --name-only HEAD~1", intern = TRUE)
changed_packages <- unique(dirname(dirname(changed_files[grepl("^packages/", changed_files)])))

if (length(changed_packages) == 0) {
  message("No packages changed")
  quit(status = 0)
}

for (pkg in changed_packages) {
  if (dir.exists(pkg)) {
    message("Checking: ", basename(pkg))
    devtools::check(pkg)
  }
}
```

---

## 9. ライセンス整備

### 9.1 方針

- ルート: GitHub 判定用 `LICENSE.md`（MIT 本文）
- 各パッケージ: CRAN/Bioc 向け `LICENSE` + `LICENSE.md`

### 9.2 タスク

| # | タスク | 詳細 |
|---|-------|------|
| 9.1 | ルート LICENSE.md 確認 | MIT 本文のみ |
| 9.2 | 各パッケージ LICENSE 作成 | `YEAR: 2024`, `COPYRIGHT HOLDER: Yusuke Matsui` |
| 9.3 | DESCRIPTION 確認 | `License: MIT + file LICENSE` |

**受入基準**
- [ ] GitHub がリポジトリを MIT と認識
- [ ] 各パッケージの DESCRIPTION で MIT 明記

---

## 10. 受入基準（Definition of Done）

### v2.0 Phase 1-2（Core 安定化）

- [ ] `packages/cfomics` が単独で `R CMD check` 通過
- [ ] `cf_methods()` が利用可能 method を返す
- [ ] `as_cf_data()` が generic として動作
- [ ] R-only methods が従来どおり動作

### v2.0 Phase 3-4（Python 分離）

- [ ] `packages/cfomicsPython` が単独で check 通過
- [ ] Core が Python 非依存で check 通過
- [ ] Python methods が registry 経由で利用可能

### v2.0 Phase 5（Sim/Bench 分離）

- [ ] `packages/cfomicsSim` でシナリオ生成可能
- [ ] `packages/cfomicsBench` でベンチマーク実行可能
- [ ] スモークベンチが再現可能（seed 固定）

### v2.0 Phase 6（Adapters/Diagnose 分離）

- [ ] Bioc オブジェクトが `as_cf_data()` で変換可能
- [ ] `cf_diagnose()` が診断結果を返す

---

## 11. 実装チケット（PR 単位）

### PR-1: packages/ への移行（Phase 1）

**ブランチ**: `feature/monorepo-structure`

**変更内容**:
- `packages/cfomics/` へ全ファイル移動
- ルート README 更新
- インストール導線更新

**DoD**: `packages/cfomics` の check が通る

**推定作業量**: 小

---

### PR-2: Core registry + as_cf_data（Phase 2）

**ブランチ**: `feature/core-registry`

**変更内容**:
- `R/registry.R` 新規作成
- `R/as_cf_data.R` 新規作成
- `cf_fit.R` リファクタ
- テスト追加

**DoD**: R-only methods が registry 経由で動作

**推定作業量**: 中

---

### PR-3: cfomicsPython skeleton（Phase 3）

**ブランチ**: `feature/cfomics-python-skeleton`

**変更内容**:
- `packages/cfomicsPython/` 骨組み作成
- stub 関数
- 最小テスト

**DoD**: 単独 check 通過

**推定作業量**: 小

---

### PR-4: Python 移植（Phase 4）

**ブランチ**: `feature/python-separation`

**変更内容**:
- Python 関連ファイル移動
- registry 登録実装
- テスト移動

**DoD**: Core 単独 green, Python 追加で methods 可視化

**推定作業量**: 大

---

### PR-5: cfomicsSim/cfomicsBench skeleton（Phase 3）

**ブランチ**: `feature/sim-bench-skeleton`

**変更内容**:
- 両パッケージの骨組み作成
- stub 関数
- 最小テスト

**DoD**: 各パッケージ単独 check 通過

**推定作業量**: 小

---

### PR-6: Sim/Bench 実装（Phase 5）

**ブランチ**: `feature/sim-bench-implementation`

**変更内容**:
- シナリオ実装
- ベンチマーク実装
- 既存 benchmark コード移植

**DoD**: スモークベンチが走る

**推定作業量**: 大

---

### PR-7: CI path filter（CI 整備）

**ブランチ**: `feature/ci-path-filter`

**変更内容**:
- `.github/workflows/r-cmd-check.yml` 作成
- `.github/workflows/r-cmd-check-python.yml` 作成

**DoD**: 変更パッケージのみ check 実行

**推定作業量**: 小

---

### PR-8: Adapters/Diagnose 分離（Phase 6、任意）

**ブランチ**: `feature/adapters-diagnose`

**変更内容**:
- `packages/cfomicsAdapters/` 作成・移植
- `packages/cfomicsDiagnose/` 作成・移植

**DoD**: 各パッケージ単独 check 通過

**推定作業量**: 中

---

## 12. 実行順序（推奨）

```
Phase 0 (準備)
    ↓
Phase 1 (packages/ 化) ← PR-1
    ↓
Phase 2 (registry) ← PR-2
    ↓
Phase 3 (skeletons) ← PR-3, PR-5 [並行可]
    ↓
Phase 4 (Python 分離) ← PR-4
    ↓
Phase 5 (Sim/Bench 実装) ← PR-6
    ↓
CI 整備 ← PR-7
    ↓
Phase 6 (Adapters/Diagnose) ← PR-8 [必要時]
```

---

## 13. AI 駆動開発向け実装指示書テンプレート

各 PR に対して、以下の形式で実装指示書を作成可能：

```markdown
# PR-X: [タイトル]

## 目的
[1-2文で目的を記載]

## 前提条件
- [ ] PR-Y がマージ済み
- [ ] main が green

## 作業手順

### Step 1: [ステップ名]
```bash
[コマンド]
```

### Step 2: [ステップ名]
[ファイル作成/編集の詳細]

## 検証手順
```bash
cd packages/[pkg]
R CMD build .
R CMD check *.tar.gz
```

## DoD チェックリスト
- [ ] check が通る
- [ ] テストが通る
- [ ] [機能要件]

## コミットメッセージ
```
[type]: [description]

[body]
```
```

---

## 付録: 現行ファイル → v2 パッケージ対応表（確定版）

| 現行ファイル | 移行先パッケージ | 備考 |
|-------------|-----------------|------|
| R/aaa-onload.R | cfomics | 維持 |
| R/cf_fit.R | cfomics | registry dispatch 化 |
| R/cf_predict.R | cfomics | 維持 |
| R/config.R | cfomics | Python 設定は分離検討 |
| R/formula_utils.R | cfomics | 維持 |
| R/result_contract.R | cfomics | v2 contract として固定 |
| R/methods_traditional.R | cfomics | 維持 |
| R/python_env.R | cfomicsPython | 移動 |
| R/methods_dowhy_gcm.R | cfomicsPython | 移動 |
| R/methods_drlearner.R | cfomicsPython | 移動 |
| R/methods_ganite.R | cfomicsPython | 移動 |
| R/methods_cavae.R | cfomicsPython | 移動 |
| R/bioc_integration.R | cfomicsAdapters | 移動 |
| R/benchmark_*.R | cfomicsSim/cfomicsBench | 分解移植 |
| R/diagnostics.R | cfomicsDiagnose | 移動 |
| R/visualization.R | cfomics | 当面維持 |
| inst/python/* | cfomicsPython | 移動 |
| inst/benchmarks/* | cfomicsBench | 移動 |

---

*Last updated: 2026-01-05*
*Version: 2.0.0*
