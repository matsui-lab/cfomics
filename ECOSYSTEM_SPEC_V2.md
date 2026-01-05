# cfomics Ecosystem 詳細仕様設計 v2.0

**Monorepo-first / Split-ready Architecture**

---

## 0. 目的とスコープ

### 0.1 目的

AI 駆動開発（Claude Code / Devin）で高速に実装を進めつつ、後日以下に耐える**分割可能な設計**へ整流化する：

- 論文化の分割（S1: Software, B1: Benchmark, A1: Application）
- 配布（CRAN/Bioc）や DOI 付与
- 依存・CI の安定化

### 0.2 v2 の中核方針

- レポジトリは当面 **monorepo** を維持するが、**`packages/` 配下に複数Rパッケージ**として配置
- 将来の独立 repo 化は「ディレクトリ単位の切り出し」で済む状態を作る
- すべての連携は **Contract（契約）**で行い、**越境（他パッケージの内部関数・内部構造への依存）を禁止**

### 0.3 非スコープ（v2でやらない）

- 応用領域向けの大規模解析テンプレ（A1 用の研究別パイプラインの完成形）
- すべての causal estimand を完全網羅（まず ATE/ITE/y0/y1 を軸）

---

## 1. 全体アーキテクチャ（責務分離）

### 1.1 コンポーネント一覧

| パッケージ | 役割 | 説明 |
|-----------|------|------|
| **cfomics** (Core) | 統一API | 統一 API・結果契約・最小の R-only 推定器 |
| **cfomicsPython** (Extension) | Python連携 | Python バックエンド・環境管理（reticulate） |
| **cfomicsAdapters** (Extension) | Bioc統合 | Bioconductor（SE/MAE/SCE）入力統一 |
| **cfomicsSim** (DGP) | シミュレーション | 包括 DGP・シナリオ・真値保持 |
| **cfomicsBench** (Benchmark) | ベンチマーク | 実行・評価・集計・再現性ハーネス |
| **cfomicsDiagnose** (Validity) | 診断 | overlap/weights/balance 等 |
| **cfomicsSensitivity** (Robustness) | 感度分析 | 未測定交絡/欠測/誤差等 |
| **cfomicsReport** (Ops) | レポート | Quarto テンプレート生成 |

### 1.2 依存方向（循環禁止）

```
cfomics (Core) ← 依存される側（ecosystem内の他パッケージに依存しない）
    ↑
    ├── cfomicsPython
    ├── cfomicsAdapters
    ├── cfomicsDiagnose
    ├── cfomicsSensitivity
    ├── cfomicsReport
    └── cfomicsBench → cfomicsSim → cfomics
```

**原則**: Core は ecosystem 内の他パッケージに依存しない。依存は片方向に固定する。

---

## 2. Cross-package Contract（v2 の最重要）

### 2.1 cf_data contract（入力の統一表現）

**目的**: data.frame / Bioc オブジェクト差を吸収し、推定器・診断・ベンチが共通に扱う。

**最小仕様（v2）**

```r
# class: cf_data
list(
  data = data.frame(),           # サンプル行が基本
  outcome_name = character(1),   # 文字列
  treatment_name = character(1), # 文字列（v2 は二値 0/1 を原則）
  covariate_names = character(), # 文字列ベクトル
  meta = list(                   # メタ情報
    n = integer(1),
    p = integer(1),
    types = list(),              # 型情報
    missing = list(),            # 欠測情報
    batch = NULL,                # バッチ情報
    time = NULL,                 # 時系列情報
    cluster = NULL               # クラスター情報
  )
)
```

**生成方法**

- `cfomics::as_cf_data(x, ...)` を **generic** として Core に定義
- `as_cf_data.data.frame` は Core に実装
- `as_cf_data.SummarizedExperiment` 等は `cfomicsAdapters` が実装

### 2.2 cfomics_result contract（結果オブジェクト）

**目的**: 推定器が異なっても同一形式で downstream（診断・感度・ベンチ・レポート）が動く。

**必須フィールド（v2固定）**

```r
# class: c("cf_model", "cfomics_result")
list(
  method = character(1),         # method id
  fit = list(                    # method-specific fit object
    model = ...,                 # 基底モデル
    res = list(                  # 必須結果
      ite = numeric(n),          # Individual Treatment Effects
      ate = numeric(1),          # Average Treatment Effect
      y0_hat = numeric(n),       # Counterfactual: untreated
      y1_hat = numeric(n),       # Counterfactual: treated
      summary = list(            # 要約統計
        ate = numeric(1),
        ate_ci_lower = numeric(1),
        ate_ci_upper = numeric(1),
        ite_mean = numeric(1),
        ite_std = numeric(1),
        ite_quantiles = list(q05, q50, q95)
      )
    )
  ),
  meta = list(                   # メタデータ
    formula = formula,
    n = integer(1),
    p = integer(1),
    outcome_name = character(1),
    treatment_name = character(1),
    covariate_names = character(),
    estimand = character(1),     # "ATE"/"ITE"/"y0"/"y1"
    software_versions = list()   # 再現性情報
  )
)
```

**検証**: `cfomics::validate_cfomics_result()` が contract を検証し、全メソッドは返却直前に必ず呼ぶ。

### 2.3 Method registry contract（拡張の差し込み方式）

**目的**: Core に新規 backend をハードコードせず、拡張パッケージが method を登録できる。

**Core が提供する API**

```r
cf_methods()
# → 利用可能 method の一覧（メタ情報付き）

cf_register_method(
  id,                            # method identifier
  fit_fun,                       # fitting function
  predict_fun,                   # prediction function
  requires_python = FALSE,       # Python 依存フラグ
  ...                            # 追加メタ情報
)
```

**重要ルール**

- `cfomicsPython` は `.onLoad()` で method 登録してよいが、**Python 初期化は禁止**（実行時に `cf_require_python()`）

### 2.4 Simulation scenario contract（cfomicsSim）

**目的**: シナリオが再現可能で、真値が常に取れる。

**必須要素**

```r
list(
  scenario_id = character(1),    # 安定ID（例: "S00_linear"）
  params = list(),               # パラメータ（yaml/json 互換）
  truth = list(                  # 真値
    ite = numeric(n),
    ate = numeric(1),
    y0 = numeric(n),
    y1 = numeric(n)
  ),
  mechanisms = list(             # メカニズムフラグ
    batch = FALSE,
    lod = FALSE,
    zero_inflation = FALSE,
    mnar = FALSE,
    cluster = FALSE,
    longitudinal = FALSE,
    unmeasured_confounding = FALSE
  )
)
```

**出力**: data.frame または SummarizedExperiment/MAE（`metadata$truth` を必須化）

### 2.5 Diagnostics / Sensitivity attachment contract

**目的**: 診断・感度分析を「結果に同梱」して運用可能にする。

```r
# cfomics_result に追加
fit$diagnostics    # 診断結果（標準スキーマ）
fit$sensitivity    # 感度分析結果（標準スキーマ）
```

`cfomicsReport` はこれらがある場合、レポートに自動埋め込み可能。

---

## 3. Package specs（各パッケージの設計仕様）

### 3.1 packages/cfomics（Core）

**責務**

- 統一 API（`cf_fit`, `predict`, `summary`, `plot` の入口）
- contract（`cf_data`, `cfomics_result`）
- method registry
- 最小の R-only 方法（gformula/ipw/grf）

**主要 API（v2）**

| 関数 | 説明 |
|-----|------|
| `cf_fit(formula, data, method, ...)` | 統一フィッティング API |
| `as_cf_data(x, ...)` | generic 入力変換 |
| `validate_cfomics_result(x)` | 結果検証 |
| `cf_methods()` | 利用可能メソッド一覧 |
| `cf_register_method()` | メソッド登録 |

**依存方針**: Imports は最小（現在の igraph 等は Suggests）

**内部構成**

```
R/
├── cf_fit.R
├── cf_predict.R
├── result_contract.R
├── formula_utils.R
├── registry.R
├── as_cf_data.R
└── methods_traditional.R  # R-only methods
```

### 3.2 packages/cfomicsPython（Python backend / env）

**責務**

- Python 環境管理（install/use/check）
- DoWhy/EconML/GANITE/CAVAE の fit/predict 実装
- registry への method 登録

**主要 API**

| 関数 | 説明 |
|-----|------|
| `cf_install_python_env()` | Python 環境インストール |
| `cf_use_python_env()` | Python 環境アクティベート |
| `cf_check_env()` | 環境診断 |

**Method IDs**: `dowhy_gcm`, `drlearner`, `ganite`, `cavae`

**重要**: `.onLoad()` での Python 初期化は禁止

**内部構成**

```
R/
├── python_env.R
├── methods_dowhy_gcm.R
├── methods_drlearner.R
├── methods_ganite.R
├── methods_cavae.R
└── register_methods.R
inst/python/
├── gcm_dowhy.py
├── ganite.py
├── visualization.py
└── requirements.txt
```

### 3.3 packages/cfomicsAdapters（Bioconductor integration）

**責務**

- `as_cf_data.SummarizedExperiment`
- `as_cf_data.MultiAssayExperiment`
- `as_cf_data.SingleCellExperiment`
- feature selection（variance/pca/none）等の「入力整形」のみ

**主要 API**: `as_cf_data()` の S3 method 追加

### 3.4 packages/cfomicsSim（DGP）

**責務**

- scenario カタログ（S00..）
- 包括 DGP（batch/LOD/ZI/MNAR/cluster/longitudinal/unmeasured confounding）
- truth を必ず同梱

**主要 API**

| 関数 | 説明 |
|-----|------|
| `cf_simulate(scenario, n, ...)` | シミュレーション実行 |
| `cf_scenarios()` | 利用可能シナリオ一覧 |
| `cf_scenario_get(id)` | シナリオ取得 |
| `cf_truth(x)` | 真値抽出 |

### 3.5 packages/cfomicsBench（Benchmark harness）

**責務**

- 繰り返し実行（methods × scenarios × reps）
- 評価（推定誤差 + 診断性能 + 計算コスト）
- 集計・可視化
- 再現性（seed、環境、設定、出力アーカイブ）

**主要 API**

| 関数 | 説明 |
|-----|------|
| `cf_benchmark_run(scenarios, methods, n_reps, ...)` | ベンチマーク実行 |
| `cf_benchmark_summarize(results)` | 結果集計 |
| `cf_benchmark_plot(results)` | 可視化 |

### 3.6 packages/cfomicsDiagnose（診断）

**責務**

- overlap/positivity
- weights 診断
- バランス検証
- 影響点検出

**主要 API**: `cf_diagnose()` が `cf_data` または `cfomics_result` を受け取って標準結果を返す

### 3.7 packages/cfomicsSensitivity（感度分析）

**責務**

- 未測定交絡
- 欠測機構
- 測定誤差

**主要 API**: `cf_sensitivity()` を提供し、`fit$sensitivity` に格納可能

### 3.8 packages/cfomicsReport（レポート）

**責務**: Quarto/Markdown テンプレで `cfomics_result` を入力し、推定→診断→感度→結論を生成

---

## 4. Monorepo repository layout

```
cfomics/
├── packages/
│   ├── cfomics/
│   ├── cfomicsPython/
│   ├── cfomicsAdapters/
│   ├── cfomicsSim/
│   ├── cfomicsBench/
│   ├── cfomicsDiagnose/
│   ├── cfomicsSensitivity/
│   └── cfomicsReport/
├── tools/
│   ├── check.R                 # パッケージ単位チェック
│   └── check_changed.R         # 変更パッケージのみチェック
├── .github/workflows/
│   ├── r-cmd-check.yml         # path filter で package 別に実行
│   └── r-cmd-check-python.yml  # Python job 隔離
├── README.md                   # エコシステム概要
├── CLAUDE.md                   # モノレポ運用ルール（越境禁止・AI向け）
└── renv.lock                   # 開発用（任意）
```

**ルート README の役割**

- ecosystem の入口
- 各 package の目的
- インストール導線
- 論文化ロードマップ

---

## 5. CI / Test / Release

### 5.1 CI（必須）

- **変更された package のみ** `R CMD check`（path filter）
- `cfomicsPython` は別 workflow（Python 不安定性を隔離）
- キャッシュ:
  - R package cache（renv か pak）
  - Python env は可能なら cache（再現性優先で戦略を固定）

### 5.2 テスト方針

| パッケージ | テスト方針 |
|-----------|----------|
| cfomics | R-only で完走 |
| cfomicsPython | Python がない場合は skip（構造契約テストは R-only で通るよう分離） |
| cfomicsSim | truth の整合性テスト（ATE と ITE の整合、seed 再現等） |
| cfomicsBench | 最小シナリオ×最小 method のスモークテスト |

### 5.3 バージョン・タグ（必須）

- タグは package 単位で打つ（将来の切り出しに耐える）
  - 例: `cfomics-v0.5.0`
- 各 package に `NEWS.md` を置く（monorepo でも package 単位で管理）

---

## 6. 現状実装からの移行計画

### Phase A: packages/ 化（最優先）

**手順**

1. ルートに `packages/` を作成
2. 現在のパッケージ一式を `packages/cfomics/` に移動
3. ルート README を ecosystem 用に置き換え
4. インストール導線を `subdir="packages/cfomics"` に更新
5. CI のパスを修正

**受入基準**: `packages/cfomics` が単独で `R CMD check` を通る

### Phase B: Core contract 固定と registry 導入

**手順**

1. `as_cf_data()` を Core に generic として固定
2. registry（`cf_register_method`, `cf_methods`）を導入
3. `validate_cfomics_result` を v2 contract として固定

**受入基準**

- R-only（gformula/ipw/grf）が `cf_fit` で動く
- `cf_methods()` で一覧が取れる

### Phase C: cfomicsPython の新設＋移植

**移植対象**

- `R/python_env.R`
- `R/methods_dowhy_gcm.R`, `methods_drlearner.R`, `methods_ganite.R`, `methods_cavae.R`
- `inst/python/*`

**受入基準**

- `cfomics` 単体では Python 非依存で完走
- `cfomics + cfomicsPython` を入れると Python method が `method=` で見える

### Phase D: Sim / Bench の新設

**手順**

1. `packages/cfomicsSim` を作り、scenario schema と `cf_simulate()/cf_truth()` を定義
2. `packages/cfomicsBench` を作り、`cf_benchmark_run/summarize` を定義
3. 現在の `benchmark_*.R` を Sim と Bench に分解移植

**受入基準**: 最小シナリオ×最小 method でスモークベンチが走る（R-only）

### Phase E（任意・後続）: Adapters / Diagnose / Sensitivity / Report

必要になった順に分割し、contract に従って attach する。

---

## 7. 現行ファイル → v2 分割先の対応表

| 現行ファイル | 移行先 |
|-------------|-------|
| `R/cf_fit.R` | cfomics |
| `R/cf_predict.R` | cfomics |
| `R/formula_utils.R` | cfomics |
| `R/result_contract.R` | cfomics |
| `R/methods_traditional.R` | cfomics |
| `R/config.R` | cfomics（Python 固有設定は cfomicsPython へ） |
| `R/python_env.R` | cfomicsPython |
| `R/methods_dowhy_gcm.R` | cfomicsPython |
| `R/methods_ganite.R` | cfomicsPython |
| `R/methods_drlearner.R` | cfomicsPython |
| `R/methods_cavae.R` | cfomicsPython |
| `inst/python/*` | cfomicsPython |
| `R/bioc_integration.R` | cfomicsAdapters（早期移行推奨） |
| `R/benchmark_*.R` | cfomicsSim / cfomicsBench（分割） |
| `R/diagnostics.R` | cfomicsDiagnose |
| `R/visualization.R` | cfomics（後で cfomicsViz 作成時に移行） |

---

## 8. AI 駆動開発向けガードレール

monorepo＋AI 開発で最も壊れやすいのは「境界侵食」なので、以下をルール化し CI でも検知する：

### 禁止事項

- 他 package の `:::`（内部呼び出し）禁止
- 他 package の未 export 関数呼び出し禁止
- `tools/` からの `source()` で共通コード共有禁止
  - 共通化したいなら package 化して export

### 必須事項

- contract 変更は必ず `validate_*` とテストを同時更新
- 破壊的変更を CI で検知

---

## 付録: v2 で最初に確定すべき仕様（最優先）

### 1. Contract v2.0

- `cf_data`
- `cfomics_result`
- registry
- scenario/truth

### 2. 依存方向

- Core は ecosystem 内に依存しない
- Bench→Sim→Core の片方向

これが固まれば、実装は AI 駆動でも破綻しにくく、後日の分割（別 repo 化）も低コストになる。

---

## 次のアクション候補

1. **「packages/ モノレポ移行」実装チケット一式**
   - 手順、変更ファイル、受入基準、CI雛形まで

2. **「Contract v2.0」仕様書（厳密版）**
   - cfomics_result / cf_data / scenario/truth の JSON schema 相当
   - テスト設計

---

*Last updated: 2026-01-05*
*Version: 2.0.0*
