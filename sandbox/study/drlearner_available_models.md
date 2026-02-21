# DR-Learnerで使用可能なモデル一覧

## 概要

EconMLのDR-Learnerは3つの段階でモデルを使用します：
1. **model_propensity**: 傾向スコア P(T=1|X) の推定
2. **model_regression**: アウトカム E[Y|T,X] の推定
3. **model_final**: 最終的なCATE τ(X) の推定

各段階で、scikit-learn互換の任意のモデルを指定できます。

---

## 1. model_propensity（分類モデル）

傾向スコア推定には**分類モデル**を使用（binary treatment の場合）。

### 線形分類モデル（sklearn.linear_model）

| モデル | 説明 | 正則化 | CV版 |
|--------|------|--------|------|
| **LogisticRegression** | ロジスティック回帰 | L2 (Ridge) | LogisticRegressionCV ✓ |
| **SGDClassifier** | 確率的勾配降下法 | L1/L2/ElasticNet | - |
| **RidgeClassifier** | Ridge分類 | L2 | RidgeClassifierCV ✓ |
| **PassiveAggressiveClassifier** | オンライン学習 | なし | - |
| **Perceptron** | パーセプトロン | なし | - |

### アンサンブル分類モデル（sklearn.ensemble）

| モデル | 説明 | 特徴 |
|--------|------|------|
| **RandomForestClassifier** | ランダムフォレスト | 非線形、高次元に強い |
| **ExtraTreesClassifier** | 極端ランダム化木 | RFより高速 |
| **GradientBoostingClassifier** | 勾配ブースティング | 高精度だが遅い |
| **HistGradientBoostingClassifier** | ヒストグラムGB | 大規模データに最適 |
| **AdaBoostClassifier** | AdaBoost | 弱学習器をブースト |
| **BaggingClassifier** | バギング | 分散削減 |
| **VotingClassifier** | 投票分類器 | 複数モデルの組み合わせ |
| **StackingClassifier** | スタッキング | メタ学習 |

### その他の分類モデル

| カテゴリ | モデル | 説明 |
|----------|--------|------|
| **SVM** (sklearn.svm) | SVC, LinearSVC, NuSVC | サポートベクターマシン |
| **NN** (sklearn.neural_network) | MLPClassifier | 多層パーセプトロン |
| **Tree** (sklearn.tree) | DecisionTreeClassifier | 決定木 |
| **KNN** (sklearn.neighbors) | KNeighborsClassifier, RadiusNeighborsClassifier | 近傍法 |
| **Naive Bayes** (sklearn.naive_bayes) | GaussianNB, MultinomialNB, BernoulliNB | ナイーブベイズ |

### 推奨設定

```python
# デフォルト（現在の実装）
model_propensity = LogisticRegressionCV(cv=5)

# 非線形効果が強い場合
model_propensity = RandomForestClassifier(n_estimators=100)

# 大規模データ
model_propensity = HistGradientBoostingClassifier()
```

---

## 2. model_regression（回帰モデル）

アウトカム推定には**回帰モデル**を使用（Y(0)とY(1)を別々に推定）。

### 線形回帰モデル（sklearn.linear_model）

#### 基本線形モデル

| モデル | 説明 | 正則化 | CV版 | 用途 |
|--------|------|--------|------|------|
| **LinearRegression** | 通常の最小二乗法 | なし | - | ベースライン |
| **Ridge** | Ridge回帰 | L2 | **RidgeCV** ✓ | 多重共線性対策 |
| **Lasso** | Lasso回帰 | L1 | **LassoCV** ✓ | 特徴選択 |
| **ElasticNet** | ElasticNet | L1+L2 | **ElasticNetCV** ✓ | バランス型 |
| **Lars** | Least Angle Regression | - | **LarsCV** ✓ | 高次元データ |
| **LassoLars** | Lasso + Lars | L1 | **LassoLarsCV** ✓ | 高速Lasso |

#### ロバスト回帰モデル

| モデル | 説明 | 特徴 |
|--------|------|------|
| **HuberRegressor** | Huber損失 | 外れ値にロバスト |
| **RANSACRegressor** | RANSAC | 外れ値除去 |
| **TheilSenRegressor** | Theil-Sen推定 | 非常にロバスト（遅い） |
| **QuantileRegressor** | 分位点回帰 | 条件付き分位点推定 |

#### ベイズ回帰モデル

| モデル | 説明 | 特徴 |
|--------|------|------|
| **BayesianRidge** | ベイジアンRidge | 事後分布を推定 |
| **ARDRegression** | 自動関連度決定 | スパース事前分布 |

#### 一般化線形モデル（GLM）

| モデル | 説明 | 分布 | 用途 |
|--------|------|------|------|
| **PoissonRegressor** | ポアソン回帰 | Poisson | カウントデータ |
| **GammaRegressor** | ガンマ回帰 | Gamma | 正の連続値 |
| **TweedieRegressor** | Tweedie回帰 | Tweedie | 複合ポアソン |

#### オンライン学習

| モデル | 説明 | 特徴 |
|--------|------|------|
| **SGDRegressor** | 確率的勾配降下 | 大規模データ、オンライン学習 |
| **PassiveAggressiveRegressor** | Passive-Aggressive | オンライン学習 |

### アンサンブル回帰モデル（sklearn.ensemble）

| モデル | 説明 | 速度 | 精度 | 解釈性 |
|--------|------|------|------|--------|
| **RandomForestRegressor** | ランダムフォレスト | 中 | 高 | 低 |
| **ExtraTreesRegressor** | 極端ランダム化木 | 高 | 高 | 低 |
| **GradientBoostingRegressor** | 勾配ブースティング | 低 | 最高 | 低 |
| **HistGradientBoostingRegressor** | ヒストグラムGB | **最高** | 最高 | 低 |
| **AdaBoostRegressor** | AdaBoost | 中 | 中 | 低 |
| **BaggingRegressor** | バギング | 中 | 中 | 低 |
| **VotingRegressor** | 投票回帰 | - | 高 | 低 |
| **StackingRegressor** | スタッキング | 低 | 最高 | 低 |

### その他の回帰モデル

| カテゴリ | モデル | 説明 | 用途 |
|----------|--------|------|------|
| **SVM** | SVR, LinearSVR, NuSVR | サポートベクター回帰 | 非線形回帰 |
| **NN** | MLPRegressor | ニューラルネットワーク | 深い非線形性 |
| **Tree** | DecisionTreeRegressor | 決定木 | 非線形、解釈可能 |
| **KNN** | KNeighborsRegressor | k近傍法 | 局所的パターン |

### 推奨設定

```python
# デフォルト（現在の実装）
model_regression = RidgeCV(cv=5)

# 非線形DGP
model_regression = RandomForestRegressor(n_estimators=100)

# 高精度重視
model_regression = HistGradientBoostingRegressor()

# スパース性重視
model_regression = LassoCV(cv=5)
```

---

## 3. model_final（最終CATE推定モデル）

最終的なCATE推定には**回帰モデル**を使用。
model_regressionと同じモデルが使用可能。

### 現在の実装（demov04）

```python
# デフォルト設定
model_final = LassoCV(cv=5)
```

### 推奨設定（DGPに応じて）

#### 線形DGP

```python
# L2正則化（Ridge）- 多重共線性対策
model_final = RidgeCV(cv=5)

# L1正則化（Lasso）- スパース推定
model_final = LassoCV(cv=5)  # 現在の設定

# ElasticNet - バランス型
model_final = ElasticNetCV(cv=5, l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99])
```

#### 非線形DGP

```python
# ランダムフォレスト
model_final = RandomForestRegressor(n_estimators=500, max_depth=10)

# 勾配ブースティング（高精度）
model_final = GradientBoostingRegressor(n_estimators=500, max_depth=5)

# ヒストグラムGB（大規模データ）
model_final = HistGradientBoostingRegressor(max_iter=500)
```

#### その他の特殊ケース

```python
# カウントデータ（RNA-seqカウント など）
model_final = PoissonRegressor(alpha=1.0)

# 外れ値が多い
model_final = HuberRegressor(epsilon=1.35)

# ニューラルネットワーク（深い非線形性）
from sklearn.neural_network import MLPRegressor
model_final = MLPRegressor(hidden_layer_sizes=(100, 50), max_iter=500)
```

---

## モデル選択のガイドライン

### 1. データサイズに応じて

| サンプルサイズ | 推奨モデル | 理由 |
|---------------|-----------|------|
| **小 (n < 200)** | RidgeCV, LassoCV | 正則化が重要 |
| **中 (200 ≤ n < 1000)** | RidgeCV, RandomForest | バランス型 |
| **大 (n ≥ 1000)** | HistGradientBoosting, RandomForest | 高精度モデル使用可能 |

### 2. 真のDGPの形状に応じて

| DGP | model_regression | model_final |
|-----|-----------------|-------------|
| **線形** | RidgeCV | RidgeCV or LassoCV |
| **非線形（滑らか）** | RandomForest | GradientBoosting |
| **非線形（不連続）** | RandomForest | RandomForest |
| **スパース** | LassoCV | LassoCV |

### 3. 計算速度重視

| 優先度 | モデル | 理由 |
|--------|--------|------|
| **最速** | LinearRegression, Ridge | 閉形式解 |
| **高速** | HistGradientBoosting | 最適化済み |
| **中速** | RandomForest | 並列化可能 |
| **低速** | GradientBoosting, SVM | 逐次処理 |

---

## 実装例

### 例1: 全て線形モデル（高速、解釈可能）

```python
from sklearn.linear_model import LogisticRegressionCV, RidgeCV, LassoCV

cf_fit_drlearner(
    X, T, Y,
    model_propensity = LogisticRegressionCV(cv=5),
    model_regression = RidgeCV(cv=5),
    model_final = LassoCV(cv=5)
)
```

### 例2: 全て非線形モデル（高精度）

```python
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

cf_fit_drlearner(
    X, T, Y,
    model_propensity = RandomForestClassifier(n_estimators=100),
    model_regression = RandomForestRegressor(n_estimators=100),
    model_final = RandomForestRegressor(n_estimators=500)
)
```

### 例3: ハイブリッド（バランス型）

```python
from sklearn.linear_model import LogisticRegressionCV, RidgeCV
from sklearn.ensemble import GradientBoostingRegressor

cf_fit_drlearner(
    X, T, Y,
    model_propensity = LogisticRegressionCV(cv=5),  # 線形
    model_regression = RidgeCV(cv=5),               # 線形
    model_final = GradientBoostingRegressor()      # 非線形
)
```

### 例4: 大規模データ用

```python
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor

cf_fit_drlearner(
    X, T, Y,
    model_propensity = HistGradientBoostingClassifier(),
    model_regression = HistGradientBoostingRegressor(),
    model_final = HistGradientBoostingRegressor()
)
```

---

## まとめ

### 利用可能なモデル数

- **線形モデル**: 22種類
- **アンサンブル**: 8種類
- **SVM**: 3種類
- **ニューラルネット**: 1種類
- **決定木**: 2種類
- **近傍法**: 2種類
- **その他**: 多数

**合計: 40種類以上の回帰モデルが使用可能**

### 現在の実装（demov04）

```python
model_propensity = LogisticRegressionCV(cv=5)  # L2正則化ロジスティック回帰
model_regression = RidgeCV(cv=5)                # L2正則化Ridge回帰
model_final = LassoCV(cv=5)                     # L1正則化Lasso回帰
```

この設定は**線形DGPに対して最適**で、n=1000でR²=0.992を達成。

### 参考文献

- [EconML DRLearner Documentation](https://www.pywhy.org/EconML/_autosummary/econml.dr.DRLearner.html)
- [Scikit-learn Linear Models](https://scikit-learn.org/stable/modules/linear_model.html)
- [Scikit-learn Ensemble Methods](https://scikit-learn.org/stable/modules/ensemble.html)
- Scikit-learn version: **1.6.1** (現在の環境)
- EconML version: **0.16.0**
