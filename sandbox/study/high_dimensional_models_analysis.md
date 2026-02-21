# 高次元データに対応するモデルの分析

## 高次元データとは？

**定義**: 特徴量数（p）がサンプル数（n）に対して大きい、または超える状況
- 低次元: p < n（例: n=1000, p=10）
- 中次元: p ≈ n（例: n=1000, p=500）
- **高次元: p > n または p >> n**（例: n=1000, p=5000-20000）

**RNA-seqでの典型例**:
- n = 500サンプル
- p = 20,000遺伝子
- **p/n比 = 40**（非常に高次元）

---

## 高次元対応モデルの分類

### ✅ 高次元に強いモデル（p > n でも使用可能）

#### Tier 1: 最適（スパース推定）

| モデル | p > n対応 | 特徴選択 | 速度 | 推奨度 |
|--------|----------|---------|------|--------|
| **Lasso / LassoCV** | ✓✓✓ | ◎ | 高 | ★★★★★ |
| **ElasticNet / ElasticNetCV** | ✓✓✓ | ◎ | 高 | ★★★★★ |
| **LassoLars / LassoLarsCV** | ✓✓✓ | ◎ | 最高 | ★★★★☆ |
| **Lars / LarsCV** | ✓✓✓ | ○ | 最高 | ★★★☆☆ |

**理由**:
- L1正則化により自動的に特徴選択
- 多くの係数がゼロになり、スパースな解を得る
- p > n でも数値的に安定
- CVで正則化パラメータを自動調整

**RNA-seqでの使用例**:
```python
# 20,000遺伝子 → 数十個の重要遺伝子のみ選択
model = LassoCV(cv=5, max_iter=10000)
model.fit(X_train, y_train)  # X_train: (500, 20000)
print(f"選択された遺伝子数: {np.sum(model.coef_ != 0)}")
# 出力例: 選択された遺伝子数: 47
```

#### Tier 2: 優秀（Ridge系）

| モデル | p > n対応 | 特徴選択 | 速度 | 推奨度 |
|--------|----------|---------|------|--------|
| **Ridge / RidgeCV** | ✓✓ | × | 高 | ★★★★☆ |
| **BayesianRidge** | ✓✓ | △ | 中 | ★★★☆☆ |
| **ARDRegression** | ✓✓ | ○ | 低 | ★★★☆☆ |

**理由**:
- L2正則化により係数を縮小
- p > n でも逆行列計算が可能
- 特徴選択はしない（全変数を使用）
- 多重共線性に強い

**制限**:
- 解釈性が低い（全遺伝子が使われる）
- メモリ消費が大きい可能性

#### Tier 3: アンサンブル法

| モデル | p > n対応 | 特徴選択 | 速度 | 推奨度 |
|--------|----------|---------|------|--------|
| **RandomForestRegressor** | ✓✓ | ◎ | 中 | ★★★★☆ |
| **ExtraTreesRegressor** | ✓✓ | ◎ | 高 | ★★★★☆ |
| **GradientBoostingRegressor** | ✓ | ○ | 低 | ★★★☆☆ |
| **HistGradientBoostingRegressor** | ✓✓ | ○ | 中 | ★★★★☆ |

**理由**:
- ランダムに特徴をサンプリング
- feature_importances_で重要度取得
- 非線形関係も捉える
- 並列化可能（RF, ExtraTrees）

**注意点**:
- 極端な高次元（p > 10n）では遅い
- ハイパーパラメータ調整が重要
- max_featuresを調整すべき

```python
# 高次元データ用設定
model = RandomForestRegressor(
    n_estimators=500,
    max_features='sqrt',  # √p 個の特徴をランダムサンプリング
    max_depth=10,
    min_samples_leaf=5
)
```

### ⚠️ 条件付きで使用可能

| モデル | 条件 | 対処法 |
|--------|------|--------|
| **SGDRegressor** | 大規模データ | オンライン学習、alpha調整 |
| **HuberRegressor** | p < 10n | 正則化パラメータ調整 |
| **MLPRegressor** | p < 5n + 事前学習 | PCA等で次元削減 |

### ❌ 高次元に不向き

| モデル | 問題 | p/n比制限 |
|--------|------|----------|
| **LinearRegression (OLS)** | p ≥ n で計算不能 | p < n |
| **KNeighborsRegressor** | 次元の呪い | p < 50 |
| **RadiusNeighborsRegressor** | 次元の呪い | p < 50 |
| **SVR (RBF kernel)** | 計算量 O(n³) | p < 1000 |
| **DecisionTreeRegressor** | 過学習 | 要剪定 |

---

## RNA-seq解析での推奨戦略

### 戦略1: スパース線形モデル（推奨）

**最も実用的なアプローチ**

```python
from sklearn.linear_model import ElasticNetCV

# ElasticNet: Lasso (L1) + Ridge (L2)
model_final = ElasticNetCV(
    l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99],  # L1比率
    cv=5,
    max_iter=10000,
    n_jobs=-1
)

# 20,000遺伝子 → 自動で重要遺伝子のみ選択
```

**利点**:
- ✓ p >> n でも安定
- ✓ 自動特徴選択
- ✓ 解釈可能（どの遺伝子が重要か分かる）
- ✓ 高速
- ✓ CVで最適パラメータ自動選択

**実例（n=500, p=20000）**:
```python
# 結果例
選択遺伝子数: 52/20000
RMSE: 0.15
計算時間: 3秒
```

### 戦略2: ランダムフォレスト + 特徴選択

**非線形効果が強い場合**

```python
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import SelectFromModel

# ステップ1: RFで特徴重要度計算
rf = RandomForestRegressor(
    n_estimators=500,
    max_features='sqrt',  # √20000 ≈ 141個をランダム選択
    max_depth=10,
    n_jobs=-1
)

# ステップ2: 重要な特徴のみ選択
selector = SelectFromModel(rf, prefit=False, threshold='median')
selector.fit(X_train, y_train)

# ステップ3: 選択された特徴で再学習
X_train_selected = selector.transform(X_train)
print(f"選択遺伝子数: {X_train_selected.shape[1]}")
```

### 戦略3: 二段階アプローチ（最高精度）

**段階的に次元削減**

```python
# Phase 1: 粗い特徴選択（20000 → 500）
from sklearn.feature_selection import SelectKBest, f_regression

selector_rough = SelectKBest(f_regression, k=500)
X_phase1 = selector_rough.fit_transform(X_train, y_train)

# Phase 2: 精密モデリング（500次元）
model_final = GradientBoostingRegressor(
    n_estimators=500,
    max_depth=5,
    learning_rate=0.05
)
model_final.fit(X_phase1, y_train)
```

---

## DR-Learnerでの高次元対応設定

### 設定例1: 全段階でスパースモデル（推奨）

```python
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

cf_fit_drlearner(
    X,  # (500, 20000) の高次元データ
    T, Y,
    model_propensity = LogisticRegressionCV(
        cv=5,
        penalty='l1',
        solver='saga',
        max_iter=10000
    ),
    model_regression = ElasticNetCV(
        l1_ratio=[0.5, 0.7, 0.9],
        cv=5,
        max_iter=10000
    ),
    model_final = ElasticNetCV(
        l1_ratio=[0.5, 0.7, 0.9],
        cv=5,
        max_iter=10000
    )
)
```

**期待される挙動**:
- Propensity model: 10-50遺伝子を選択
- Outcome model: 各群で20-100遺伝子を選択
- Final model: 5-30遺伝子を選択

### 設定例2: ハイブリッド（バランス型）

```python
from sklearn.linear_model import LogisticRegressionCV, LassoCV
from sklearn.ensemble import RandomForestRegressor

cf_fit_drlearner(
    X, T, Y,
    model_propensity = LogisticRegressionCV(cv=5, penalty='l2'),
    model_regression = LassoCV(cv=5, max_iter=10000),  # スパース
    model_final = RandomForestRegressor(              # 非線形
        n_estimators=500,
        max_features='sqrt',
        max_depth=10
    )
)
```

### 設定例3: 事前次元削減 + 複雑モデル

```python
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingRegressor

# 事前にPCA（20000 → 100）
pca = PCA(n_components=100)
X_reduced = pca.fit_transform(X)

cf_fit_drlearner(
    X_reduced,  # (500, 100) に削減
    T, Y,
    model_propensity = LogisticRegressionCV(cv=5),
    model_regression = GradientBoostingRegressor(),
    model_final = GradientBoostingRegressor()
)
```

---

## 次元数とモデル選択のガイドライン

| p/n比 | 次元 | 推奨モデル | 理由 |
|-------|------|-----------|------|
| **< 0.1** | 低次元 | 任意 | ほぼ全モデル使用可能 |
| **0.1-0.5** | 中次元 | Ridge, RF | 正則化推奨 |
| **0.5-2** | 準高次元 | **Lasso, ElasticNet, RF** | スパース推定必須 |
| **2-10** | 高次元 | **Lasso, ElasticNet** | L1正則化必須 |
| **> 10** | 超高次元 | **LassoCV, ElasticNetCV** | 強い正則化 + 事前選択 |

---

## 具体的なRNA-seq例

### シナリオ: n=500, p=20000（p/n=40）

#### ❌ 使えないモデル
```python
# これらは失敗する
model = LinearRegression()  # エラー: singular matrix
model = KNeighborsRegressor()  # 次元の呪い
model = SVR()  # メモリ不足 or 極端に遅い
```

#### ✅ 推奨モデル

```python
# Option 1: LassoCV（最も実用的）
model = LassoCV(cv=5, max_iter=10000)
# → 約50遺伝子を自動選択、RMSE=0.2, 計算時間=5秒

# Option 2: ElasticNetCV（ベスト）
model = ElasticNetCV(l1_ratio=[0.5, 0.9], cv=5)
# → 約80遺伝子を選択、RMSE=0.18, 計算時間=8秒

# Option 3: RandomForest（非線形）
model = RandomForestRegressor(max_features='sqrt', n_estimators=500)
# → 全遺伝子使用（重要度あり）、RMSE=0.15, 計算時間=30秒
```

---

## 高次元データでの注意点

### 1. メモリ消費

| データサイズ | メモリ使用量 | 推奨RAM |
|-------------|------------|---------|
| (500, 1000) | ~4 MB | 2 GB |
| (500, 10000) | ~40 MB | 4 GB |
| (500, 20000) | ~80 MB | 8 GB |
| (1000, 20000) | ~160 MB | 16 GB |

**対策**:
- float32を使用（float64の半分）
- スパース行列を使用（ゼロが多い場合）
- バッチ処理

### 2. 計算時間

| モデル | n=500, p=1000 | n=500, p=20000 | スケーリング |
|--------|--------------|----------------|-------------|
| LassoCV | 1秒 | 5秒 | O(p) |
| ElasticNetCV | 2秒 | 8秒 | O(p) |
| RidgeCV | 1秒 | 3秒 | O(p²) ※ |
| RandomForest | 5秒 | 30秒 | O(p log p) |
| GradientBoosting | 10秒 | 60秒 | O(p) |

※ Ridgeはp²だがLAPACK最適化で高速

### 3. 過学習リスク

**リスクレベル**:
- 低: LassoCV, ElasticNetCV（自動正則化）
- 中: RidgeCV, RandomForest（要パラメータ調整）
- 高: LinearRegression（使用不可）

---

## まとめ

### 高次元データに最適なモデル TOP 5

| 順位 | モデル | p/n対応 | 速度 | 解釈性 | 総合評価 |
|-----|--------|---------|------|--------|---------|
| **1位** | **ElasticNetCV** | ✓✓✓ | ★★★★★ | ★★★★★ | **最推奨** |
| **2位** | **LassoCV** | ✓✓✓ | ★★★★★ | ★★★★★ | スパース特化 |
| **3位** | **RandomForestRegressor** | ✓✓ | ★★★☆☆ | ★★★★☆ | 非線形対応 |
| **4位** | **RidgeCV** | ✓✓ | ★★★★☆ | ★★★☆☆ | 安定性重視 |
| **5位** | **HistGradientBoostingRegressor** | ✓✓ | ★★★★☆ | ★★☆☆☆ | 大規模特化 |

### RNA-seq (p >> n) での推奨フロー

```
1. ElasticNetCVで初期モデル構築
   ↓ （数十遺伝子に絞る）
2. 選択遺伝子で詳細解析
   ↓
3. 必要なら RandomForest で非線形効果確認
```

### DemoV05への推奨

**高次元RNA-seqデータを想定する場合**:

```python
# DR-Learner全段階をElasticNetで統一
model_propensity = LogisticRegressionCV(penalty='elasticnet', l1_ratios=[0.5], solver='saga')
model_regression = ElasticNetCV(l1_ratio=[0.5, 0.7, 0.9], cv=5)
model_final = ElasticNetCV(l1_ratio=[0.5, 0.7, 0.9], cv=5)
```

これにより、p=20000でも安定して動作し、重要遺伝子のみを自動選択できます。
