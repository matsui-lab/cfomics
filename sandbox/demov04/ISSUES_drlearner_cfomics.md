# cfomics DR-Learner 実装の問題点

> 後日 GitHub Issue として登録予定
> 発見日: 2026-02-19 (demov04 ベンチマーク実施時)

---

## 概要

`packages/cfomics/R/methods_drlearner.R` の DR-Learner (EconML `DRLearner`) ラッパー実装に、
ITE/CATE 推定の安定性・精度に影響する複数の問題がある。
demov04 ベンチマークで τ̂ が ±50,000〜90,000 に爆発する事象が発生し、原因調査で判明した。

---

## Issue 1: モデル引数が EconML に渡されていない

**重要度**: High

`cf_fit_drlearner()` は `model_propensity` と `model_regression` を引数として受け取るが、
実際の EconML `DRLearner()` コンストラクタに渡していない。

```r
# 現状のコード（概要）
cf_fit_drlearner <- function(..., model_propensity = NULL, model_regression = NULL, ...) {
  # model_propensity, model_regression は使われずに捨てられる
  dr <- econml$dr$DRLearner(cv = cv)
  # ↑ model_propensity, model_regression が渡されていない
}
```

**影響**:
- EconML のデフォルト（LogisticRegressionCV / WeightedLassoCVWrapper）が常に使用される
- ユーザーが RandomForest 等を指定しても無視される
- 線形モデルが非線形な傾向スコアを捉えられず、DR 疑似アウトカムが爆発する

**修正案**:
```r
dr <- econml$dr$DRLearner(
  model_propensity = model_propensity,
  model_regression = model_regression,
  cv = cv
)
```

---

## Issue 2: 傾向スコアのクリッピングがない

**重要度**: High

DR-Learner の疑似アウトカムは以下の式で計算される:

```
Ỹ = μ̂₁(X) - μ̂₀(X) + T(Y - μ̂₁(X))/ê(X) - (1-T)(Y - μ̂₀(X))/(1-ê(X))
```

ê(X) が 0 や 1 に近いとき、`1/ê(X)` や `1/(1-ê(X))` が爆発する。
EconML 自体にはある程度の内部クリッピングがあるが、cfomics 側で明示的な制御がない。

**影響**:
- 小サンプル（n=100〜500）で τ̂ が数万〜数十万に発散
- 1 件の異常値が RMSE 全体を支配する
- R=1 でも reproducibility が不安定

**修正案**:
- `min_propensity` 引数を追加（デフォルト 0.01〜0.05）
- EconML の `min_propensity` パラメータを利用するか、後処理でクリップ

---

## Issue 3: model_final が指定されていない

**重要度**: Medium

EconML `DRLearner` の `model_final`（CATE を推定する最終段モデル）がデフォルトのまま。
デフォルトは `StatsModelsLinearRegression`（正則化なし OLS）。

**影響**:
- 共変量が多い場合に過学習
- DR 疑似アウトカムの外れ値に対して非ロバスト
- 外れ値 1 件が回帰係数を大きく歪める

**修正案**:
- `model_final` 引数を追加
- デフォルトを正則化付きモデル（Lasso/Ridge）に変更
- または `NonParamDML` / `ForestDRLearner` への切り替えオプション

---

## Issue 4: アウトカム Y の標準化がない

**重要度**: Medium

Y をそのまま Python に渡しているため、Y のスケールが大きい場合に
DR 疑似アウトカムのスケールも比例して大きくなる。

**影響**:
- Y のスケールに依存して数値不安定性が増幅
- 学習率やペナルティの適切な設定が困難

**修正案**:
- Y の z-score 標準化を前処理として追加
- 予測時に逆変換して元のスケールに戻す
- オプションとして `standardize_y = TRUE` (デフォルト) を追加

---

## Issue 5: 乱数シードの不完全な制御

**重要度**: Low

`random_state` 引数は NumPy のグローバルシードのみを設定しており、
scikit-learn の内部 CV 分割や EconML の内部処理のシードが制御されていない。

```r
# 現状
np <- reticulate::import("numpy")
np$random$seed(as.integer(random_state))
```

**影響**:
- 同一シードでも結果が完全に再現しない場合がある
- CV の分割がランダムに変わり、特定の分割で異常値が発生

**修正案**:
- `DRLearner(random_state=seed, ...)` として EconML に直接渡す
- `sklearn.model_selection` の `random_state` も揃える

---

## Issue 6: 信頼区間が未実装

**重要度**: Low

`ate_ci_lower` / `ate_ci_upper` が常に `NA` で返される。
EconML は `effect_inference()` で信頼区間を計算可能。

**修正案**:
```python
inference = model.effect_inference(X)
ci = inference.conf_int(alpha=0.05)
```

---

## Issue 7: newdata 予測で type="y0"/"y1" が未サポート

**重要度**: Low

`predict.cf_model` で `type = "y0"` や `type = "y1"` を指定しても
反事実予測が返されない（ITE のみ対応）。

**修正案**:
- EconML の `model.predict(Y, T, X)` や内部の `model_regression` を使って
  μ̂₀(X), μ̂₁(X) を取得し返す

---

## demov04 での暫定対処

ベンチマークスクリプト側で以下のクリップ処理を適用:

```r
tau_hat_raw <- -as.numeric(predict(fit, newdata = ..., type = "ite"))
tau_hat <- pmin(20, pmax(-20, tau_hat_raw))
```

これは根本修正ではなく、cfomics パッケージ側での対応が必要。

---

## 優先度まとめ

| # | Issue | 重要度 | 影響範囲 |
|---|-------|--------|----------|
| 1 | モデル引数が渡されない | High | 全ユーザー |
| 2 | 傾向スコアクリッピングなし | High | 小〜中サンプル |
| 3 | model_final 未指定 | Medium | 高次元データ |
| 4 | Y 標準化なし | Medium | スケール依存 |
| 5 | 乱数シード不完全 | Low | 再現性 |
| 6 | 信頼区間未実装 | Low | 推論時 |
| 7 | y0/y1 予測未サポート | Low | 反事実分析時 |
