# DML、DR、OrthoForestの詳細解説

## 重要な理解: 推定器は「組み合わせ」

これらの手法は**フレームワーク**であり、内部で使用する推定器（モデル）を**自由に選択・組み合わせ可能**です。

---

## 1. DML（Double Machine Learning）の詳細

### 構造: 3段階の推定器

DMLは**3つの異なる推定器**を組み合わせて動作します：

```
DML = ニュアンス推定器（2つ） + 最終推定器（1つ）
```

#### 推定器の役割

| 段階 | 推定器 | パラメータ名 | 推定対象 | 選択可能？ |
|------|--------|------------|---------|----------|
| **第1段階** | アウトカムモデル | `model_y` | E[Y\|X,W] | ✅ 自由 |
| **第1段階** | 処置モデル | `model_t` | E[T\|X,W] or P(T=1\|X,W) | ✅ 自由 |
| **第2段階** | 最終モデル | `model_final` | τ(X) | ✅ 自由 |

### アルゴリズムの詳細

#### ステップ1: ニュアンス推定（Cross-fitting）

```python
# K-fold cross-fittingループ
for fold in 1..K:
    # 訓練データでモデル学習
    model_y.fit(X_train, Y_train)
    model_t.fit(X_train, T_train)

    # 検証データで予測
    Y_hat[fold] = model_y.predict(X_val)
    T_hat[fold] = model_t.predict(X_val)
```

**重要**: 各foldで異なるデータを使用（過学習防止）

#### ステップ2: 残差計算（直交化）

```python
# アウトカムの残差
Y_res = Y - Y_hat

# 処置の残差
T_res = T - T_hat
```

**直交化の意味**: XとWの影響を「除去」

#### ステップ3: 最終推定

```python
# 残差同士で回帰
model_final.fit(X, Y_res / T_res)  # 簡略化した式
tau_hat = model_final.predict(X_test)
```

### DMLの具体的な組み合わせ例

#### 例1: 全て線形（LinearDML）

```python
from econml.dml import LinearDML
from sklearn.linear_model import LassoCV, RidgeCV

dml = LinearDML(
    model_y=RidgeCV(cv=5),              # アウトカム: Ridge
    model_t=LassoCV(cv=5),              # 処置: Lasso
    discrete_treatment=True,
    cv=5
)
# model_final は自動的に LinearRegression
```

**推定器の組み合わせ**:
- `model_y`: RidgeCV（正則化線形回帰）
- `model_t`: LassoCV（正則化線形回帰）
- `model_final`: LinearRegression（**固定**、変更不可）

#### 例2: ニュアンスは非線形、最終は線形（推奨）

```python
from econml.dml import LinearDML
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

dml = LinearDML(
    model_y=RandomForestRegressor(n_estimators=500),  # 非線形
    model_t=RandomForestClassifier(n_estimators=500), # 非線形
    discrete_treatment=True,
    cv=5
)
# model_final は LinearRegression（固定）
```

**なぜこの組み合わせ？**
- ニュアンス推定: 複雑な関数形を捉える（RF）
- 最終推定: 解釈可能な線形（Linear）
- **理論保証**: ニュアンスが非線形でも最終が線形ならOK

#### 例3: 全て非線形（NonParamDML）

```python
from econml.dml import NonParamDML
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier

dml = NonParamDML(
    model_y=GradientBoostingRegressor(),     # 非線形
    model_t=GradientBoostingClassifier(),    # 非線形
    model_final=GradientBoostingRegressor(), # 非線形 ←ここも選択可能！
    discrete_treatment=True,
    cv=5
)
```

**推定器の組み合わせ**:
- 全段階で**任意の**sklearn互換モデルを選択可能
- 最も柔軟だが、過学習リスク高

#### 例4: 高次元特化（SparseLinearDML）

```python
from econml.dml import SparseLinearDML
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

dml = SparseLinearDML(
    model_y=RandomForestRegressor(),     # p >> n でも使える
    model_t=RandomForestClassifier(),    # p >> n でも使える
    alpha='auto',                        # Lasso正則化を自動調整
    cv=5
)
# model_final は LassoCV（固定、α自動選択）
```

**推定器の組み合わせ**:
- `model_y`: 任意（RFなど）
- `model_t`: 任意（RFなど）
- `model_final`: LassoCV（**固定**、高次元用）

#### 例5: 因果フォレスト（CausalForestDML）

```python
from econml.dml import CausalForestDML

cf_dml = CausalForestDML(
    model_y=RandomForestRegressor(),     # ニュアンス推定
    model_t=RandomForestClassifier(),    # ニュアンス推定
    n_estimators=4000,                   # 最終モデルのパラメータ
    min_samples_leaf=5,
    cv=5
)
# model_final は GRF（一般化ランダムフォレスト、固定）
```

**推定器の組み合わせ**:
- `model_y`: 任意
- `model_t`: 任意
- `model_final`: GRF Causal Forest（**固定**、非線形）

### DML手法ごとの自由度

| DML手法 | model_y | model_t | model_final | 自由度 |
|---------|---------|---------|-------------|--------|
| **DML** | ✅ 自由 | ✅ 自由 | ✅ **自由** | 最大 |
| **LinearDML** | ✅ 自由 | ✅ 自由 | ❌ Linear固定 | 高 |
| **SparseLinearDML** | ✅ 自由 | ✅ 自由 | ❌ Lasso固定 | 高 |
| **NonParamDML** | ✅ 自由 | ✅ 自由 | ✅ **自由** | 最大 |
| **KernelDML** | ✅ 自由 | ✅ 自由 | ❌ Kernel固定 | 中 |
| **CausalForestDML** | ✅ 自由 | ✅ 自由 | ❌ GRF固定 | 中 |

---

## 2. DR（Doubly Robust）の詳細

### 構造: 3段階の推定器（DMLと同じ）

```
DR = ニュアンス推定器（2つ） + 最終推定器（1つ）
```

#### 推定器の役割

| 段階 | 推定器 | パラメータ名 | 推定対象 | 選択可能？ |
|------|--------|------------|---------|----------|
| **第1段階** | 傾向スコア | `model_propensity` | P(T=1\|X) | ✅ 自由 |
| **第1段階** | アウトカム回帰 | `model_regression` | E[Y\|T,X] | ✅ 自由 |
| **第2段階** | 最終モデル | `model_final` | τ(X) | ✅ 自由 |

### DMLとの違い

| 項目 | DML | DR |
|------|-----|-----|
| ニュアンス1 | E[Y\|X] | **P(T=1\|X)** (傾向スコア) |
| ニュアンス2 | E[T\|X] | **E[Y\|T,X]** (アウトカム) |
| スコア関数 | 直交化スコア | **DR スコア** |
| 頑健性 | ニュアンス誤差に | **モデル誤指定に** |

### アルゴリズムの詳細

#### ステップ1: ニュアンス推定

```python
# 傾向スコア
model_propensity.fit(X, T)
e_hat = model_propensity.predict_proba(X)[:, 1]

# アウトカム回帰（T=0とT=1で別々）
model_regression.fit(X[T==0], Y[T==0])  # μ₀
mu0_hat = model_regression.predict(X)

model_regression.fit(X[T==1], Y[T==1])  # μ₁
mu1_hat = model_regression.predict(X)
```

#### ステップ2: 疑似アウトカム（DR変換）

```python
# DR疑似アウトカム
pseudo_outcome = (
    T * (Y - mu1_hat) / e_hat
    - (1 - T) * (Y - mu0_hat) / (1 - e_hat)
    + mu1_hat - mu0_hat
)
```

**二重頑健性**: この式が、どちらかのモデルが正しければ不偏。

#### ステップ3: 最終推定

```python
model_final.fit(X, pseudo_outcome)
tau_hat = model_final.predict(X_test)
```

### DRの具体的な組み合わせ例

#### 例1: 汎用DRLearner（demov04使用）

```python
from econml.dr import DRLearner
from sklearn.linear_model import LogisticRegressionCV, RidgeCV, LassoCV

dr = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),  # 分類器
    model_regression=RidgeCV(cv=5),               # 回帰器
    model_final=LassoCV(cv=5),                    # 回帰器
    cv=5
)
```

**推定器の組み合わせ**:
- `model_propensity`: LogisticRegressionCV（**分類器**）
- `model_regression`: RidgeCV（回帰器）
- `model_final`: LassoCV（回帰器）
- **全て自由に選択可能**

#### 例2: LinearDRLearner

```python
from econml.dr import LinearDRLearner
from sklearn.linear_model import LogisticRegression

dr = LinearDRLearner(
    model_propensity=LogisticRegression(),   # 分類器
    model_regression='auto',                 # 自動選択（Ridge）
    cv=5
)
# model_final は LinearRegression（固定）
```

**推定器の組み合わせ**:
- `model_propensity`: 任意の分類器
- `model_regression`: 'auto' または任意の回帰器
- `model_final`: LinearRegression（**固定**）

#### 例3: SparseLinearDRLearner（高次元用）

```python
from econml.dr import SparseLinearDRLearner

dr = SparseLinearDRLearner(
    model_propensity='auto',    # LogisticRegressionCV
    model_regression='auto',    # LassoCV
    alpha='auto',               # Lasso正則化パラメータ
    cv=5
)
# model_final は LassoCV（固定）
```

**推定器の組み合わせ**:
- `model_propensity`: 'auto'（LogisticRegressionCV）
- `model_regression`: 'auto'（LassoCV）
- `model_final`: LassoCV（**固定**、α自動調整）

#### 例4: ForestDRLearner

```python
from econml.dr import ForestDRLearner
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestRegressor

dr = ForestDRLearner(
    model_propensity=LogisticRegression(),
    model_regression=RandomForestRegressor(),
    n_estimators=4000,        # 最終モデルのパラメータ
    min_samples_leaf=5,
    cv=5
)
# model_final は GRF Causal Forest（固定）
```

**推定器の組み合わせ**:
- `model_propensity`: 任意の分類器
- `model_regression`: 任意の回帰器
- `model_final`: GRF Causal Forest（**固定**）

#### 例5: 全て非線形（最大の柔軟性）

```python
from econml.dr import DRLearner
from sklearn.ensemble import (
    GradientBoostingClassifier,
    GradientBoostingRegressor
)

dr = DRLearner(
    model_propensity=GradientBoostingClassifier(n_estimators=500),
    model_regression=GradientBoostingRegressor(n_estimators=500),
    model_final=GradientBoostingRegressor(n_estimators=500),
    cv=5
)
```

**全て自由に選択！**

### DR手法ごとの自由度

| DR手法 | model_propensity | model_regression | model_final | 自由度 |
|--------|-----------------|-----------------|-------------|--------|
| **DRLearner** | ✅ 自由 | ✅ 自由 | ✅ **自由** | 最大 |
| **LinearDRLearner** | ✅ 自由 | ✅ 自由 | ❌ Linear固定 | 高 |
| **SparseLinearDRLearner** | ✅ 自由 | ✅ 自由 | ❌ Lasso固定 | 高 |
| **ForestDRLearner** | ✅ 自由 | ✅ 自由 | ❌ GRF固定 | 中 |

---

## 3. OrthoForest の詳細

### 構造: 局所ニュアンス推定

OrthoForestは**ニュアンス推定を局所的に行う**点が異なります。

```
OrthoForest = フォレスト分割 + 局所ニュアンス推定 + 局所CATE推定
```

#### 推定器の役割

| 段階 | 推定器 | パラメータ名 | 推定対象 | 選択可能？ |
|------|--------|------------|---------|----------|
| **分割** | フォレスト | n_trees, max_depth | X空間の分割 | パラメータのみ |
| **局所ニュアンス** | アウトカム/処置 | `model_Y`, `model_T` | 局所E[Y], E[T] | ✅ 自由 |
| **局所CATE** | - | - | 各葉での τ | 内部で自動 |

### アルゴリズムの詳細

#### ステップ1: データ駆動的な分割

```python
# フォレストでX空間を分割
for tree in forest:
    # 分割基準: 局所的なCATE推定の精度を最大化
    # （通常のRFとは異なる基準）
    split(X)
```

**重要**: 分割はCATEの異質性を捉えるように最適化

#### ステップ2: 各葉ノードで局所推定

```python
for leaf in tree:
    # この葉のデータのみ使用
    X_leaf = X[samples_in_leaf]
    Y_leaf = Y[samples_in_leaf]
    T_leaf = T[samples_in_leaf]

    # 局所ニュアンス推定
    model_Y.fit(X_leaf, Y_leaf)
    model_T.fit(X_leaf, T_leaf)

    # 局所CATE推定
    tau_leaf = estimate_local_cate(...)
```

**DML/DRとの違い**: ニュアンス推定が**葉ごとに異なる**

#### ステップ3: 予測

```python
# テストサンプルを適切な葉に割り当て
leaf_id = forest.get_leaf(X_test)

# その葉の局所CATEを使用
tau_hat = tau_leaf[leaf_id]
```

### OrthoForestの具体的な組み合わせ例

#### 例1: DMLOrthoForest

```python
from econml.orf import DMLOrthoForest
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

orf = DMLOrthoForest(
    n_trees=4000,              # フォレストのパラメータ
    min_leaf_size=10,
    max_depth=10,
    subsample_ratio=0.7,
    model_Y=RandomForestRegressor(n_estimators=100),  # 局所ニュアンス
    model_T=RandomForestClassifier(n_estimators=100), # 局所ニュアンス
    cv=5
)
orf.fit(Y, T, X=X, W=W)
```

**推定器の組み合わせ**:
- フォレスト構造: パラメータで制御
- `model_Y`: 任意の回帰器（局所で使用）
- `model_T`: 任意の分類/回帰器（局所で使用）
- 最終CATE: 自動（局所平均）

#### 例2: DROrthoForest

```python
from econml.orf import DROrthoForest
from sklearn.linear_model import LogisticRegression, Ridge

orf = DROrthoForest(
    n_trees=4000,
    min_leaf_size=10,
    propensity_model=LogisticRegression(),  # 局所傾向スコア
    model_Y=Ridge(),                        # 局所アウトカム
    cv=5
)
```

**推定器の組み合わせ**:
- `propensity_model`: 任意の分類器（局所で使用）
- `model_Y`: 任意の回帰器（局所で使用）
- DRスコアを局所的に計算

### グローバル vs 局所の比較

| 特性 | DML/DR（グローバル） | OrthoForest（局所） |
|------|-------------------|------------------|
| **ニュアンス推定** | 全データで1つのモデル | **葉ごとに異なるモデル** |
| **適応性** | 固定的 | **データ駆動的** |
| **パラメータ数** | 少ない | **多い**（葉の数だけ） |
| **計算コスト** | 低 | **高** |
| **過学習リスク** | 低 | やや高 |
| **用途** | 標準的CATE推定 | **局所的異質性探索** |

### 具体例: 局所適応の効果

```python
# データ: X₁ < 0.5 と X₁ >= 0.5 で全く異なるDGP

# グローバル（DML）
dml = LinearDML(...)
# → X₁ < 0.5 と X₁ >= 0.5 で同じニュアンスモデルを使用
# → 両方のパターンの「平均」を捉える

# 局所（OrthoForest）
orf = DMLOrthoForest(...)
# → X₁ < 0.5 と X₁ >= 0.5 で**異なる**ニュアンスモデル
# → 各領域に特化したモデルを学習
```

---

## 推定器選択のガイド

### 1. DML: どの組み合わせを選ぶ？

#### シナリオ別推奨

| シナリオ | 推奨DML手法 | model_y | model_t | model_final |
|---------|-----------|---------|---------|-------------|
| **標準（線形DGP）** | LinearDML | RidgeCV | LassoCV | Linear（固定） |
| **高次元（p>>n）** | SparseLinearDML | RF | RF | Lasso（固定） |
| **非線形DGP** | CausalForestDML | RF | RF | GRF（固定） |
| **探索的** | NonParamDML | GBM | GBM | GBM |
| **最大柔軟性** | DML | 任意 | 任意 | **任意** |

#### 実装例: 高次元RNA-seq

```python
from econml.dml import SparseLinearDML
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

dml = SparseLinearDML(
    model_y=RandomForestRegressor(
        n_estimators=500,
        max_features='sqrt'  # √p 個の特徴
    ),
    model_t=RandomForestClassifier(
        n_estimators=500,
        max_features='sqrt'
    ),
    alpha='auto',  # Lasso α を自動調整
    cv=5
)
# → 20,000遺伝子でも動作、重要遺伝子のみ選択
```

### 2. DR: どの組み合わせを選ぶ？

#### シナリオ別推奨

| シナリオ | 推奨DR手法 | model_propensity | model_regression | model_final |
|---------|----------|-----------------|-----------------|-------------|
| **標準（RCT）** | DRLearner | LogisticCV | RidgeCV | LassoCV |
| **高次元** | SparseLinearDRLearner | 'auto' | 'auto' | Lasso（固定） |
| **非線形** | ForestDRLearner | Logistic | RF | GRF（固定） |
| **最大柔軟性** | DRLearner | 任意 | 任意 | **任意** |

#### 実装例: RCT + 高次元

```python
from econml.dr import DRLearner
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

dr = DRLearner(
    model_propensity=LogisticRegressionCV(
        penalty='l1',
        solver='saga',
        cv=5
    ),
    model_regression=ElasticNetCV(
        l1_ratio=[0.5, 0.9],
        cv=5
    ),
    model_final=ElasticNetCV(
        l1_ratio=[0.5, 0.9],
        cv=5
    ),
    cv=5
)
# → 全段階でスパース推定、遺伝子選択
```

### 3. DML vs DR: どちらを選ぶ？

| 基準 | DML | DR |
|------|-----|-----|
| **RCTデータ** | ○ | ◎（傾向スコア既知で有利） |
| **観察研究** | ◎ | ◎ |
| **高次元** | ◎ | ◎ |
| **頑健性重視** | ○ | ◎ |
| **実装の簡潔さ** | 同等 | 同等 |

**結論**: RCTなら**DR**、観察研究なら**同等**

### 4. OrthoForest: いつ使う？

| 使用場面 | 理由 |
|---------|------|
| ✅ 局所的に異なるパターン疑い | 適応的に検出 |
| ✅ 探索的分析 | 異質性の発見 |
| ✅ 中規模データ（n=1000-10000） | 計算可能 |
| ❌ 標準的CATE推定 | DML/DRで十分 |
| ❌ 小サンプル（n<500） | 過学習リスク |
| ❌ 超大規模（n>100000） | 計算コスト高 |

---

## まとめ

### 推定器は「組み合わせ」

```
手法 = フレームワーク + 推定器の組み合わせ

例: DRLearner = DR理論 + (LogisticCV, RidgeCV, LassoCV)
```

### 選択の自由度

| 手法 | 第1段階 | 第2段階 | 最終段階 | 自由度 |
|------|---------|---------|---------|--------|
| **DML** | ✅ | ✅ | ✅ | ★★★ |
| **DRLearner** | ✅ | ✅ | ✅ | ★★★ |
| **LinearDML** | ✅ | ✅ | ❌ | ★★ |
| **LinearDRLearner** | ✅ | ✅ | ❌ | ★★ |
| **SparseLinearDML** | ✅ | ✅ | ❌ | ★★ |
| **SparseLinearDRLearner** | ✅ | ✅ | ❌ | ★★ |
| **CausalForestDML** | ✅ | ✅ | ❌ | ★★ |
| **ForestDRLearner** | ✅ | ✅ | ❌ | ★★ |
| **OrthoForest** | ✅ | ✅ | - | ★ |

### DemoV04の選択

```python
DRLearner(
    model_propensity = LogisticRegressionCV(cv=5),  # 選択
    model_regression = RidgeCV(cv=5),                # 選択
    model_final = LassoCV(cv=5),                     # 選択
    cv=5
)
```

**全て自由に選択可能！** 線形DGPに最適化した組み合わせ。

### DemoV05への推奨（高次元想定）

```python
# Option 1: SparseLinearDRLearner（シンプル）
SparseLinearDRLearner(
    model_propensity='auto',  # LogisticCV
    model_regression='auto',  # LassoCV
    alpha='auto'              # 自動調整
)

# Option 2: DRLearner（柔軟）
DRLearner(
    model_propensity=LogisticRegressionCV(penalty='l1'),
    model_regression=ElasticNetCV(l1_ratio=[0.5,0.9]),
    model_final=ElasticNetCV(l1_ratio=[0.5,0.9])
)
```

推定器は**自由に組み合わせて最適化可能**！
