# DR-Learner 実装の修正内容

## 問題の発生条件

**DR-Learner, n=500, linear, SexGene** で RMSE=19.94（100%クリップ失敗）

## 修正内容

### **FIX #1: モデル引数が EconML に渡されていない（HIGH priority）**

**修正前（`methods_drlearner.R:36`）**:
```r
dr <- DRLearner(cv = as.integer(cv))
# model_propensity, model_regression が無視されている！
```

**修正後**:
```r
dr <- DRLearner(
  model_propensity = model_prop_py,
  model_regression = model_reg_py,
  model_final = model_final_py,        # 新規追加
  cv = as.integer(cv),
  min_propensity = min_propensity,     # 新規追加
  random_state = as.integer(random_state)  # 新規追加
)
```

**効果**:
- ユーザー指定のモデルが実際に使用される
- デフォルトを正則化付きモデルに変更（Ridge, LassoCV）
- 外れ値に対する耐性が向上

---

### **FIX #2: 傾向スコアのクリッピング追加（HIGH priority）**

**修正前**: クリッピングなし → `1/ê(X)` が爆発

**修正後**:
```r
# 新規引数
min_propensity = 0.05  # デフォルト

# DRLearner に渡す
DRLearner(
  ...
  min_propensity = min_propensity
)
```

**効果**:
- 傾向スコア ê(X) が [0.05, 0.95] にクリップされる
- DR 疑似アウトカムの爆発を防ぐ
- 小サンプル（n=100〜500）での安定性が向上

---

### **FIX #3: model_final の明示的指定（MEDIUM priority）**

**修正前**: デフォルト（OLS、正則化なし）→ 外れ値に弱い

**修正後**:
```r
# デフォルトを LassoCV に変更
if (identical(model_final, "auto")) {
  LassoCV <- sklearn$linear_model$LassoCV
  model_final_py <- LassoCV(
    cv = as.integer(cv),
    random_state = as.integer(random_state)
  )
}
```

**効果**:
- DR 疑似アウトカムの外れ値に対してロバスト
- 高次元データでの過学習を防止
- 正則化による安定性向上

---

### **FIX #4: アウトカム Y の標準化（MEDIUM priority）**

**修正前**: Y をそのまま Python に渡す → スケール依存の不安定性

**修正後**:
```r
# 新規引数
standardize_y = TRUE  # デフォルト

# 標準化
if (standardize_y) {
  Y_mean <- mean(Y)
  Y_sd <- sd(Y)
  Y <- (Y - Y_mean) / Y_sd
}

# フィット後に元に戻す
if (standardize_y) {
  ite <- ite * Y_sd  # ITE は差分なのでスケールのみ
}
```

**効果**:
- Y のスケールに依存しない数値安定性
- 学習率やペナルティの調整が容易
- 異なるデータセット間での比較が容易

---

### **FIX #5: 乱数シード制御の改善（LOW priority）**

**修正前**:
```r
# NumPy のグローバルシードのみ
np$random$seed(as.integer(random_state))
```

**修正後**:
```r
# NumPy + Python random + EconML 本体
np$random$seed(as.integer(random_state))
random <- reticulate::import("random")
random$seed(as.integer(random_state))

# DRLearner に直接渡す
DRLearner(
  ...
  random_state = as.integer(random_state)
)
```

**効果**:
- 完全な再現性の保証
- CV fold 分割の固定
- デバッグが容易

---

## デフォルトモデルの変更

| 用途 | 旧実装 | 新実装 | 理由 |
|------|--------|--------|------|
| `model_propensity` | LogisticRegressionCV | LogisticRegressionCV (max_iter=1000) | 収束性向上 |
| `model_regression` | WeightedLassoCVWrapper | **RidgeCV** | より安定 |
| `model_final` | StatsModelsLinearRegression | **LassoCV** | 外れ値耐性 |

---

## 修正の優先度

1. **HIGH**: FIX #1 (モデル引数), FIX #2 (傾向スコアクリップ)
   - これらがないと n=500 失敗は解決しない
2. **MEDIUM**: FIX #3 (model_final), FIX #4 (Y標準化)
   - 安定性と汎用性の向上
3. **LOW**: FIX #5 (乱数シード)
   - 再現性の保証

---

## 検証方法

`test_fix_n500.R` スクリプトで失敗条件（n=500, linear, SexGene）を再実行：

```r
# 修正版を使用
source("methods_drlearner_FIXED.R")

# n=500, linear, SexGene で再検証
# 期待結果: RMSE < 1.0（クリップなし）
```

---

## 今後の対応

1. **即時対応**: `packages/cfomics/R/methods_drlearner.R` を修正版で置換
2. **テスト追加**: `tests/testthat/test-drlearner-stability.R` を作成
3. **ドキュメント更新**: `min_propensity`, `standardize_y` 引数の説明を追加
4. **GitHub Issue 登録**: `ISSUES_drlearner_cfomics.md` の内容を Issue として登録
5. **ベンチマーク再実行**: demov04 全体を再実行して改善を確認
