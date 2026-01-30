# 診断・仮定チェック拡充 設計ドキュメント

## 概要

cfomics パッケージの診断機能を拡充する。現在のバランスチェック（SMD）とオーバーラップチェック（PS分布）に加え、感度分析、機能的診断、CATE検証を追加し、統合診断レポート関数でまとめる。

## 新規関数

### 1. `cf_sensitivity(result, type = "evalue", alpha = 0.05)`

E-value（VanderWeele & Ding 2017）を計算。未観測交絡への頑健性を定量化。

- `evalue = RR + sqrt(RR * (RR - 1))`
- ATE → RR 変換: `RR = exp(ATE / SD(Y))`
- 戻り値: `cf_sensitivity` クラス（evalue_point, evalue_ci, interpretation）
- S3: `print.cf_sensitivity()`, `plot.cf_sensitivity()`（バイアスプロット）

### 2. `cf_model_diagnostics(result, data, formula)`

PS・結果モデルの適合度診断。

- PS診断: C-statistic (AUC), 極端な重み検出
- 結果モデル診断: 残差の正規性検定、残差 vs 予測値
- 既存 `cf_balance_check` の結果を統合
- 戻り値: `cf_model_diagnostics` クラス
- S3: `print.cf_model_diagnostics()`, `plot.cf_model_diagnostics()`（2x2パネル）

### 3. `cf_cate_diagnostics(result, data, formula)`

処置効果の異質性を検証。

- 異質性検定（Wald test）
- GATES（四分位グループATE）— Chernozhukov et al. (2018)
- BLP（Best Linear Predictor）
- 定数 ITE の場合は自動スキップ
- 戻り値: `cf_cate_diagnostics` クラス
- S3: `print.cf_cate_diagnostics()`, `plot.cf_cate_diagnostics()`（GATES棒グラフ）

### 4. `cf_diagnostics_report(result, data, formula, include = c("sensitivity", "model", "cate"))`

全診断を一括実行。

- 各診断の結果 + Pass/Warning/Fail サマリー
- 戻り値: `cf_diagnostics_report` クラス
- S3: `print.cf_diagnostics_report()`, `plot.cf_diagnostics_report()`（サマリー図）

## 判定基準

| 診断 | Pass | Warning | Fail |
|------|------|---------|------|
| E-value | > 2.0 | 1.5-2.0 | < 1.5 |
| PS C-statistic | 0.6-0.85 | 0.5-0.6 or 0.85-0.95 | < 0.5 or > 0.95 |
| 極端な重み | < 1% | 1-5% | > 5% |
| バランス (max SMD) | < 0.1 | 0.1-0.25 | > 0.25 |

## ファイル構成

- `packages/cfomics/R/sensitivity.R` — cf_sensitivity
- `packages/cfomics/R/model_diagnostics.R` — cf_model_diagnostics
- `packages/cfomics/R/cate_diagnostics.R` — cf_cate_diagnostics
- `packages/cfomics/R/diagnostics_report.R` — cf_diagnostics_report（統合）
- `packages/cfomics/tests/testthat/test-sensitivity.R`
- `packages/cfomics/tests/testthat/test-model-diagnostics.R`
- `packages/cfomics/tests/testthat/test-cate-diagnostics.R`
- `packages/cfomics/tests/testthat/test-diagnostics-report.R`
