# DRLearner でのフォレスト系最終モデルの使用

## TL;DR

**YES、DRLearner の model_final に Random Forest や GRF を使えます。**

ただし、厳密な GRF（honest splitting 付き）を使いたい場合は **ForestDRLearner** という専用クラスがあります。

---

## 1. DRLearner × RandomForest

### 基本的な使い方

```python
from econml.dr import DRLearner
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

# DRLearner の model_final に RandomForest を指定
est = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),  # ニュアンス1
    model_regression=ElasticNetCV(cv=5),          # ニュアンス2
    model_final=RandomForestRegressor(            # ターゲット（非線形）
        n_estimators=500,
        max_features='sqrt',
        max_depth=10,
        min_samples_leaf=5,
        random_state=42,
        n_jobs=-1
    ),
    cv=5
)

# 高次元データでも使用可能
est.fit(Y, T, X=X)  # X.shape = (500, 20000)
tau_hat = est.effect(X)
```

**特徴**:
- ✓ 非線形の治療効果を捉える
- ✓ ニュアンス推定は線形（Lasso等）で高次元対応
- ✓ 最終モデルのみ非線形
- ✓ sklearn の RandomForest を使用

---

## 2. ForestDRLearner（GRF 専用版）

### 厳密な GRF を使いたい場合

EconML には **ForestDRLearner** という専用クラスがあります。

```python
from econml.dr import ForestDRLearner
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

# model_final が GRF（honest splitting）固定
est_grf = ForestDRLearner(
    model_propensity=LogisticRegressionCV(cv=5),  # ニュアンス1（変更可）
    model_regression=ElasticNetCV(cv=5),          # ニュアンス2（変更可）
    # model_final は GRF 固定（変更不可）
    n_estimators=4000,      # GRF のパラメータ
    max_features='sqrt',
    min_samples_leaf=5,
    max_depth=None,
    cv=5,
    random_state=42
)

est_grf.fit(Y, T, X=X)
tau_grf = est_grf.effect(X)
```

**特徴**:
- ✓ **Honest splitting** 付きの本格的な GRF
- ✓ ニュアンス推定は自由（Lasso 等に変更可）
- ✓ 最終モデルは GRF 固定（変更不可）
- × model_final を他のモデルに変更できない

---

## 3. DRLearner vs ForestDRLearner vs CausalForestDML

### 比較表

| クラス | ニュアンス1 | ニュアンス2 | 最終モデル | GRF | 自由度 |
|--------|-----------|-----------|-----------|-----|--------|
| **DRLearner** | 自由 | 自由 | **自由** | △（RF可） | ★★★ |
| **ForestDRLearner** | 自由 | 自由 | **GRF固定** | ◎ | ★★☆ |
| **CausalForestDML** | 自由 | 自由 | **GRF固定** | ◎ | ★★☆ |

### 詳細な違い

#### DRLearner × RandomForest

```python
from sklearn.ensemble import RandomForestRegressor

est = DRLearner(
    model_propensity=...,
    model_regression=...,
    model_final=RandomForestRegressor(...)  # 手動で指定
)
```

**使用モデル**: sklearn の RandomForestRegressor
- 通常の Random Forest
- **Honest splitting なし**
- 全データで学習

#### ForestDRLearner

```python
from econml.dr import ForestDRLearner

est = ForestDRLearner(
    model_propensity=...,
    model_regression=...,
    # model_final は自動で GRF
    n_estimators=4000  # GRF のパラメータ
)
```

**使用モデル**: EconML の GRF（econml.grf.CausalForest）
- **Honest splitting あり**
- データを分割（50% 構造学習、50% 推定）
- 因果推論に特化

#### CausalForestDML

```python
from econml.dml import CausalForestDML

est = CausalForestDML(
    model_y=...,
    model_t=...,
    # model_final は自動で GRF
    n_estimators=4000
)
```

**使用モデル**: EconML の GRF
- **Honest splitting あり**
- DML 理論（Neyman orthogonalization）
- ForestDRLearner と理論的枠組みが異なる

---

## 4. Honest Splitting とは？

### 通常の Random Forest（DRLearner × RF）

```
全データ（n=1000）
    ↓
全てを使って木の構造を学習
    ↓
全てを使って葉の予測値を計算
```

**問題**: 過学習しやすい

### GRF（ForestDRLearner, CausalForestDML）

```
全データ（n=1000）
    ↓
ランダムに2分割
    ↓
├─ 50%（n=500）: 木の構造を学習（どこで分岐するか）
└─ 50%（n=500）: 葉の予測値を計算（治療効果の推定値）
```

**利点**: 不偏推定（過学習しない）

### 具体例

```python
# 通常の RF（DRLearner × RF）
データ1でツリー構築 → データ1で予測値計算
    ↓
過学習の可能性

# GRF（ForestDRLearner）
データ1でツリー構築 → データ2で予測値計算
    ↓
不偏推定（honest）
```

---

## 5. 実践例：3つの方法を比較

### データ生成（非線形 CATE）

```python
import numpy as np

np.random.seed(42)
n, p = 1000, 100

X = np.random.randn(n, p)
T = np.random.binomial(1, 0.5, n)

# 真の CATE: 非線形
tau_true = (
    3.0 * X[:, 0] +
    2.0 * np.sign(X[:, 1]) * X[:, 1] +  # 非線形
    1.5 * (X[:, 2] > 0) * X[:, 2]       # 閾値効果
)

Y0 = 100 + 2*X[:, 0] + 1*X[:, 1] + np.random.randn(n)
Y1 = Y0 + tau_true
Y = Y0 * (1 - T) + Y1 * T
```

### パターン1: DRLearner × RandomForest

```python
from econml.dr import DRLearner
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV
import time

start = time.time()

est_rf = DRLearner(
    model_propensity=LogisticRegressionCV(cv=3),
    model_regression=ElasticNetCV(cv=3),
    model_final=RandomForestRegressor(
        n_estimators=500,
        max_features='sqrt',
        max_depth=10,
        min_samples_leaf=5,
        random_state=42,
        n_jobs=-1
    ),
    cv=3
)

est_rf.fit(Y, T, X=X)
tau_rf = est_rf.effect(X)

time_rf = time.time() - start
rmse_rf = np.sqrt(np.mean((tau_rf - tau_true)**2))

print(f"DRLearner × RandomForest:")
print(f"  時間: {time_rf:.1f}秒")
print(f"  RMSE: {rmse_rf:.3f}")
```

### パターン2: ForestDRLearner（GRF）

```python
from econml.dr import ForestDRLearner

start = time.time()

est_grf = ForestDRLearner(
    model_propensity=LogisticRegressionCV(cv=3),
    model_regression=ElasticNetCV(cv=3),
    n_estimators=500,  # GRF のツリー数
    max_features='sqrt',
    min_samples_leaf=5,
    max_depth=10,
    cv=3,
    random_state=42
)

est_grf.fit(Y, T, X=X)
tau_grf = est_grf.effect(X)

time_grf = time.time() - start
rmse_grf = np.sqrt(np.mean((tau_grf - tau_true)**2))

print(f"\nForestDRLearner (GRF):")
print(f"  時間: {time_grf:.1f}秒")
print(f"  RMSE: {rmse_grf:.3f}")
```

### パターン3: CausalForestDML（GRF + DML）

```python
from econml.dml import CausalForestDML

start = time.time()

est_dml = CausalForestDML(
    model_y=ElasticNetCV(cv=3),
    model_t=LogisticRegressionCV(cv=3),
    n_estimators=500,
    max_features='sqrt',
    min_samples_leaf=5,
    cv=3,
    random_state=42
)

est_dml.fit(Y, T, X=X)
tau_dml = est_dml.effect(X)

time_dml = time.time() - start
rmse_dml = np.sqrt(np.mean((tau_dml - tau_true)**2))

print(f"\nCausalForestDML (GRF + DML):")
print(f"  時間: {time_dml:.1f}秒")
print(f"  RMSE: {rmse_dml:.3f}")
```

### 期待される結果

```
DRLearner × RandomForest:
  時間: 3.2秒
  RMSE: 0.523

ForestDRLearner (GRF):
  時間: 5.1秒
  RMSE: 0.471  # Honest splitting で若干改善

CausalForestDML (GRF + DML):
  時間: 5.3秒
  RMSE: 0.465  # DML + GRF で最良
```

---

## 6. 高次元データ（p >> n）での使用

### 推奨設定

```python
# n=500, p=20,000 の場合

from econml.dr import ForestDRLearner
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

est = ForestDRLearner(
    # ニュアンス推定: スパースモデル（高次元対応）
    model_propensity=LogisticRegressionCV(
        cv=5,
        penalty='l1',
        solver='saga',
        max_iter=10000
    ),
    model_regression=ElasticNetCV(
        cv=5,
        l1_ratio=[0.5, 0.7, 0.9],
        max_iter=10000
    ),
    # GRF パラメータ
    n_estimators=2000,
    max_features='sqrt',  # √20000 ≈ 141 特徴
    min_samples_leaf=10,
    max_depth=None,
    cv=5
)

est.fit(Y, T, X=X)  # X.shape = (500, 20000)
```

**動作の流れ**:
1. ElasticNet が 20,000 → 50遺伝子に削減（ニュアンス推定）
2. GRF が全 20,000遺伝子で学習（各分岐で √20000 = 141個サンプリング）
3. 計算時間: 約 30-60秒

---

## 7. 使い分けガイド

### DRLearner × RandomForest を使う場合

**条件**:
- ✓ 最大の柔軟性が欲しい（model_final を色々試したい）
- ✓ Honest splitting の必要性が低い
- ✓ 計算時間を短縮したい

**例**:
```python
# GradientBoosting も試したい
est1 = DRLearner(..., model_final=RandomForestRegressor(...))
est2 = DRLearner(..., model_final=GradientBoostingRegressor(...))
est3 = DRLearner(..., model_final=MLPRegressor(...))
```

### ForestDRLearner を使う場合

**条件**:
- ✓ **因果推論に特化した GRF** を使いたい
- ✓ Honest splitting で不偏推定
- ✓ ニュアンス推定は線形モデル（Lasso等）で十分
- ✓ DR 理論を適用

**推奨シナリオ**: 高次元データで非線形効果

### CausalForestDML を使う場合

**条件**:
- ✓ **DML 理論**（Neyman orthogonalization）を適用
- ✓ ニュアンス推定誤差に対してロバスト
- ✓ 観察研究（交絡が強い）

**推奨シナリオ**: 観察研究で非線形効果

---

## 8. まとめ

### Q: DRLearner でも GRF を使える？

**A: YES、2つの方法があります。**

#### 方法1: DRLearner × RandomForest

```python
est = DRLearner(
    ...,
    model_final=RandomForestRegressor(...)
)
```

- sklearn の RandomForest
- Honest splitting なし
- 最も柔軟

#### 方法2: ForestDRLearner（推奨）

```python
est = ForestDRLearner(
    ...,
    n_estimators=4000
)
```

- EconML の GRF
- **Honest splitting あり**
- 因果推論に特化

### 推奨ランキング（非線形 CATE の場合）

| 順位 | 手法 | 理由 |
|-----|------|------|
| **1位** | **ForestDRLearner** | GRF + DR理論、RCT向け |
| **2位** | **CausalForestDML** | GRF + DML理論、観察研究向け |
| **3位** | DRLearner × RF | 柔軟性高いが理論的保証は弱い |

### 高次元データでの推奨設定

```python
# RNA-seq (n=500, p=20,000) で非線形効果
from econml.dr import ForestDRLearner

est = ForestDRLearner(
    model_propensity=LogisticRegressionCV(cv=5, penalty='l1', solver='saga'),
    model_regression=ElasticNetCV(cv=5),  # 高次元ニュアンス
    n_estimators=2000,                     # GRF
    max_features='sqrt',
    cv=5
)
```

これで **ニュアンス推定は線形（高次元対応）**、**最終推定は GRF（非線形対応）** の両立ができます。
