# LassoCV の詳細解説

## 基本概念

### Lasso とは？

**Lasso (Least Absolute Shrinkage and Selection Operator)**
- L1正則化付き線形回帰
- 多くの係数をゼロにする（スパース推定）
- 高次元データ (p >> n) で特徴選択と推定を同時実行

数式:
```
最小化: ||y - Xβ||² + α * ||β||₁
```

- `||β||₁ = |β₁| + |β₂| + ... + |βₚ|` (L1ノルム)
- `α`: 正則化パラメータ（大きいほどスパース）

### LassoCV = Lasso + 自動α選択

**問題**: α をどう決める？
- α が小さい → 過学習（係数が多く残る）
- α が大きい → 過度にスパース（重要な変数も消える）

**解決**: Cross-Validation で最適 α を自動選択

## Lasso vs LassoCV の違い

### Lasso（手動）

```python
from sklearn.linear_model import Lasso

# α を手動で指定
model = Lasso(alpha=0.1)  # ユーザーが決める
model.fit(X, y)

# 問題: 0.1 が最適かわからない
```

**手順**:
1. ユーザーが α を決定
2. その α でモデル学習
3. 精度が悪ければ α を手動調整して再実行

### LassoCV（自動）

```python
from sklearn.linear_model import LassoCV

# α は自動で選択
model = LassoCV(cv=5)  # CV fold 数のみ指定
model.fit(X, y)

# 自動で最適 α が選ばれる
print(f"最適 α: {model.alpha_}")
```

**手順**:
1. 複数の α 候補を自動生成（例: 0.001, 0.01, 0.1, 1, 10）
2. 各 α で 5-fold CV を実行
3. 最も CV スコアが良い α を選択
4. 全データでその α を使って最終モデル学習

## LassoCV の内部動作

### ステップ1: α 候補の生成

```python
# デフォルト（alphas=None の場合）
alphas = np.logspace(-4, 1, 100)
# 0.0001, 0.00013, 0.00017, ..., 5.62, 7.50, 10.0
# 100個の α を対数スケールで生成
```

### ステップ2: Cross-Validation

```python
# 疑似コード
best_score = -np.inf
best_alpha = None

for alpha in alphas:
    cv_scores = []

    for fold in range(5):  # 5-fold CV
        X_train_fold, X_val_fold = split(X_train, fold)
        y_train_fold, y_val_fold = split(y_train, fold)

        model = Lasso(alpha=alpha)
        model.fit(X_train_fold, y_train_fold)
        score = model.score(X_val_fold, y_val_fold)
        cv_scores.append(score)

    mean_score = np.mean(cv_scores)

    if mean_score > best_score:
        best_score = mean_score
        best_alpha = alpha

# 最適 α で全データを学習
final_model = Lasso(alpha=best_alpha)
final_model.fit(X_train, y_train)
```

### ステップ3: 最終モデル

```python
# ユーザー側
model = LassoCV(cv=5)
model.fit(X_train, y_train)

# 結果
print(f"最適 α: {model.alpha_}")
print(f"選択特徴数: {np.sum(model.coef_ != 0)}")
```

## パラメータ詳細

### 主要パラメータ

```python
from sklearn.linear_model import LassoCV

model = LassoCV(
    alphas=None,        # α候補（None なら自動生成）
    cv=5,               # CV fold 数
    max_iter=10000,     # 最適化の最大反復
    tol=1e-4,           # 収束判定の閾値
    n_jobs=-1,          # 並列化（-1 = 全CPU使用）
    random_state=42     # 再現性
)
```

### alphas の指定

#### パターン1: 自動（推奨）

```python
model = LassoCV(alphas=None, cv=5)
# 自動で 100個の α を生成
```

#### パターン2: 手動指定

```python
# より細かく探索
alphas = np.logspace(-5, 2, 200)
model = LassoCV(alphas=alphas, cv=5)
```

#### パターン3: 範囲指定

```python
# 特定の範囲のみ
alphas = [0.001, 0.01, 0.1, 1.0]
model = LassoCV(alphas=alphas, cv=5)
```

### cv の選択

| cv | 意味 | 計算時間 | 精度 | 推奨 |
|----|------|---------|------|------|
| 3 | 3-fold CV | 低 | 中 | 大規模データ |
| **5** | 5-fold CV | 中 | **高** | **デフォルト推奨** |
| 10 | 10-fold CV | 高 | 最高 | 小規模データ |

## 高次元データでの使用例

### RNA-seq (n=500, p=20,000)

```python
from sklearn.linear_model import LassoCV
import numpy as np

# データ準備
X_train.shape  # (400, 20000) - 20,000遺伝子
y_train.shape  # (400,)

# LassoCV で学習
model = LassoCV(
    cv=5,
    max_iter=10000,
    n_jobs=-1,
    random_state=42
)

model.fit(X_train, y_train)

# 結果確認
print(f"最適 α: {model.alpha_}")
# 出力例: 最適 α: 0.0234

print(f"選択遺伝子数: {np.sum(model.coef_ != 0)} / {X_train.shape[1]}")
# 出力例: 選択遺伝子数: 47 / 20000

# 重要遺伝子のインデックス
important_genes = np.where(model.coef_ != 0)[0]
print(f"重要遺伝子ID: {important_genes}")
# 出力例: 重要遺伝子ID: [23, 156, 489, 1024, ...]

# 予測
y_pred = model.predict(X_test)
rmse = np.sqrt(np.mean((y_pred - y_test)**2))
print(f"RMSE: {rmse:.3f}")
```

## 他の CV 付きモデルとの比較

### Lasso ファミリー

| モデル | 正則化 | 特徴選択 | 速度 | 推奨度 |
|--------|--------|---------|------|--------|
| **LassoCV** | L1 | ◎ | 中 | ★★★★★ |
| **ElasticNetCV** | L1 + L2 | ◎ | 中 | ★★★★★ |
| **LassoLarsCV** | L1 (LARS) | ◎ | 高 | ★★★★☆ |
| **RidgeCV** | L2 | × | 高 | ★★★☆☆ |

### ElasticNetCV（LassoCV の改良版）

```python
from sklearn.linear_model import ElasticNetCV

model = ElasticNetCV(
    l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0],  # L1 vs L2 比率
    cv=5,
    max_iter=10000
)
# l1_ratio=1.0 のとき LassoCV と同等
# l1_ratio=0.0 のとき RidgeCV と同等
```

**利点**:
- L1 (スパース) + L2 (安定性) の両立
- 高次元データでは **LassoCV より ElasticNetCV 推奨**

## DRLearner での使用例

### パターン1: model_final に LassoCV

```python
from econml.dr import DRLearner
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV, LassoCV

est = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=ElasticNetCV(cv=5),
    model_final=LassoCV(cv=5, max_iter=10000),  # 最終モデルに LassoCV
    cv=5
)
```

### パターン2: 全段階で CV 付き

```python
est = DRLearner(
    model_propensity=LogisticRegressionCV(  # Propensity も CV
        cv=5,
        penalty='l1',
        solver='saga',
        max_iter=10000
    ),
    model_regression=ElasticNetCV(  # Outcome も CV
        cv=5,
        l1_ratio=[0.5, 0.9],
        max_iter=10000
    ),
    model_final=LassoCV(  # Final も CV
        cv=5,
        max_iter=10000
    ),
    cv=5  # DR 全体の CV
)
```

**注意**: CV が多重にネストされる
- DR の CV (cv=5)
- 各モデルの内部 CV (cv=5)
- 計算時間は長くなるが精度は向上

## まとめ

### LassoCV の特徴

1. **α を自動選択**（手動調整不要）
2. **Cross-Validation で最適化**
3. **高次元データに最適**（p >> n OK）
4. **スパース推定**（特徴選択込み）

### いつ使うか？

**以下の場合に使用**:
- ✓ 高次元データ (p > n)
- ✓ 特徴選択が必要
- ✓ α を手動調整したくない
- ✓ 線形関係を仮定

**使わない場合**:
- × ElasticNetCV の方がより安定（推奨）
- × 非線形 → RandomForest, GradientBoosting
- × 特徴選択不要 → RidgeCV

### LassoCV vs ElasticNetCV

| 項目 | LassoCV | ElasticNetCV |
|------|---------|--------------|
| 正則化 | L1 のみ | L1 + L2 |
| スパース性 | 最大 | 中〜大 |
| 安定性 | 中 | **高** |
| 高次元耐性 | 高 | **最高** |
| **推奨度** | ★★★★☆ | **★★★★★** |

**結論**: 高次元データでは **ElasticNetCV** の方が推奨されることが多いです。
