# 高次元データにおける因果推論手法の選択ガイド

## TL;DR - 第一選択

**YES、DRLearner × Lasso/ElasticNet が第一選択です（多くの場合）**

ただし、状況によって最適な選択は変わります。以下、詳細な意思決定フローを示します。

---

## 意思決定フローチャート

```
高次元データ (p >> n) で因果推論
    |
    ├─ データ構造は？
    |   ├─ 線形 or 準線形 → DRLearner/LinearDML × Lasso/ElasticNet ★第一選択★
    |   ├─ 非線形（滑らかな曲線） → CausalForestDML (GRF backend)
    |   └─ 非線形（局所的変化） → DMLOrthoForest/DROrthoForest
    |
    ├─ 解釈性の優先度は？
    |   ├─ 高い（どの遺伝子が重要か知りたい） → SparseLinearDRLearner ★解釈性特化★
    |   ├─ 中程度 → DRLearner × ElasticNet
    |   └─ 低い（予測精度のみ） → ForestDRLearner or DML × RandomForest
    |
    ├─ 計算リソースは？
    |   ├─ 限定的（速度優先） → LinearDML × LassoCV ★最速★
    |   ├─ 中程度 → DRLearner × ElasticNetCV
    |   └─ 豊富 → DMLOrthoForest × Lasso（局所適応）
    |
    └─ 研究デザインは？
        ├─ RCT（ランダム化実験） → DRLearner系 ★DR理論的に有利★
        ├─ 観察研究（交絡強い） → DML系（orthogonalization重視）
        └─ 操作変数あり → DMLIV, LinearDMLIV
```

---

## シナリオ別推奨手法

### ✅ シナリオ1: RNA-seq標準解析（最も一般的）

**状況**:
- n = 500-2000サンプル
- p = 20,000遺伝子
- 線形 or 準線形関係を仮定
- どの遺伝子が重要か知りたい
- RCT or 観察研究

**第一選択**: **DRLearner × ElasticNetCV**

```python
from econml.dr import DRLearner
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

est = DRLearner(
    model_propensity=LogisticRegressionCV(
        cv=5,
        penalty='l1',
        solver='saga',
        max_iter=10000
    ),
    model_regression=ElasticNetCV(
        l1_ratio=[0.5, 0.7, 0.9],
        cv=5,
        max_iter=10000
    ),
    model_final=ElasticNetCV(
        l1_ratio=[0.5, 0.7, 0.9],
        cv=5,
        max_iter=10000
    ),
    cv=5
)
```

**理由**:
- ✓ ElasticNet = Lasso (特徴選択) + Ridge (安定性)
- ✓ 自動で重要遺伝子を数十個に絞る
- ✓ p >> n でも数値的に安定
- ✓ CV で正則化パラメータ自動調整
- ✓ DR理論により propensity/outcome どちらかの推定が正しければ一致推定
- ✓ 計算時間: 数秒〜数十秒

**代替案**:
- さらに解釈性重視 → **SparseLinearDRLearner**（model_final が Lasso 固定）
- さらに速度重視 → **LinearDML × LassoCV**（DR より若干速い）

---

### ✅ シナリオ2: 最大限の解釈性（遺伝子選択が主目的）

**状況**:
- どの遺伝子が治療効果に寄与するか厳密に知りたい
- スパース性を最大化
- レポート/論文でモデルの係数を報告

**第一選択**: **SparseLinearDRLearner**

```python
from econml.dr import SparseLinearDRLearner

est = SparseLinearDRLearner(
    model_propensity=LogisticRegressionCV(cv=5, penalty='l1', solver='saga'),
    model_regression=ElasticNetCV(cv=5),
    # model_final は自動で Lasso（交差検証付き）
    cv=5,
    random_state=42
)

# 重要遺伝子の取得
est.fit(Y, T, X=X)
coef = est.coef_  # 治療効果の係数
important_genes = np.where(coef != 0)[0]
print(f"選択遺伝子数: {len(important_genes)} / {p}")
```

**理由**:
- ✓ **model_final が Lasso 固定** → 最大のスパース性
- ✓ 係数の解釈が明確（どの遺伝子が +/- 効果）
- ✓ 統計的推論（信頼区間）も可能
- ✓ DRLearner より制約が強いが解釈性は最高

**vs DRLearner の違い**:
- DRLearner: model_final を自由に選べる（ElasticNet も可）
- SparseLinearDRLearner: model_final = Lasso で固定（最もスパース）

---

### ✅ シナリオ3: 非線形効果が疑われる場合

**状況**:
- 遺伝子発現と治療効果の関係が非線形
- 閾値効果やU字型の効果
- 線形モデルでは RMSE が悪い

**第一選択**: **CausalForestDML** (GRF backend)

```python
from econml.dml import CausalForestDML
from sklearn.linear_model import LassoCV

est = CausalForestDML(
    model_y=LassoCV(cv=5, max_iter=10000),  # ニュアンス推定は Lasso
    model_t=LogisticRegressionCV(cv=5),
    # model_final は Random Forest（自動、変更不可）
    n_estimators=4000,
    min_samples_leaf=5,
    max_depth=None,
    max_features='sqrt',  # √p 個の特徴をランダムサンプリング
    cv=5
)
```

**理由**:
- ✓ ニュアンス推定は Lasso で高次元対応
- ✓ 最終モデルは Random Forest で非線形捕捉
- ✓ GRF の honest splitting で不偏性
- ✓ 特徴重要度も取得可能
- ✓ 滑らかな非線形関係に強い

**注意点**:
- × 解釈性は低い（どの遺伝子が重要かは分かるが、係数は無い）
- × 計算時間は長め（数分）
- × 局所的な急激な変化には弱い → その場合は OrthoForest

---

### ✅ シナリオ4: 局所的な治療効果の変化（サブグループ解析）

**状況**:
- 特定の遺伝子プロファイルで治療効果が急激に変わる
- サブグループごとに異なるモデルが必要
- 個別化医療の観点

**第一選択**: **DMLOrthoForest** or **DROrthoForest**

```python
from econml.orf import DMLOrthoForest
from sklearn.linear_model import LassoCV, LogisticRegressionCV

est = DMLOrthoForest(
    model_Y=LassoCV(cv=5),  # グローバルニュアンス推定
    model_T=LogisticRegressionCV(cv=5),
    # 各リーフで局所的にニュアンスを再推定
    n_trees=500,
    min_leaf_size=10,
    max_depth=10,
    subsample_ratio=0.5
)
```

**理由**:
- ✓ **各リーフで異なるニュアンスモデル** → 局所適応
- ✓ GRF よりも柔軟（global nuisance → local refinement）
- ✓ 急激な変化に強い

**代償**:
- × 計算時間が長い（GRF の 10-100倍）
- × 大規模データ（n > 5000）では現実的でない可能性
- × 解釈性は GRF よりさらに低い

**使い分け**:
- 滑らかな非線形 → **CausalForestDML (GRF)**
- 局所的な急変 → **DMLOrthoForest**

---

### ✅ シナリオ5: 極限の計算速度（大規模データ、反復実験）

**状況**:
- n = 10,000+ の大規模データ
- ベンチマーク実験で何百回も実行
- 速度が最優先

**第一選択**: **LinearDML × LassoCV**

```python
from econml.dml import LinearDML
from sklearn.linear_model import LassoCV, LogisticRegressionCV

est = LinearDML(
    model_y=LassoCV(cv=3, max_iter=5000, n_jobs=-1),  # CV を 3-fold に削減
    model_t=LogisticRegressionCV(cv=3, max_iter=5000, n_jobs=-1),
    # model_final は線形回帰（自動）
    cv=3,  # CV fold 数削減
    random_state=42
)
```

**理由**:
- ✓ **LinearDML は最速**（DR より若干速い）
- ✓ model_final が線形回帰固定 → 計算量最小
- ✓ n_jobs=-1 で並列化
- ✓ CV=3 で精度と速度のバランス

**vs DRLearner**:
- LinearDML: 約 20-30% 高速
- DRLearner: 若干遅いが DR 理論の利点

---

### ✅ シナリオ6: 観察研究で交絡が強い場合

**状況**:
- RCT ではなく観察研究
- 未観測交絡の懸念は低いが、観測交絡は強い
- propensity score の推定が難しい

**第一選択**: **DML系（LinearDML, DML）**

```python
from econml.dml import DML
from sklearn.linear_model import ElasticNetCV, LogisticRegressionCV
from sklearn.ensemble import RandomForestRegressor

est = DML(
    model_y=ElasticNetCV(cv=5),
    model_t=LogisticRegressionCV(cv=5, penalty='l1', solver='saga'),
    model_final=RandomForestRegressor(
        n_estimators=500,
        max_features='sqrt',
        max_depth=10
    ),
    cv=5
)
```

**理由**:
- ✓ **Neyman orthogonalization** が交絡バイアスに強い
- ✓ propensity/outcome 両方の推定誤差を1次で除去
- ✓ DR より理論的に安定（観察研究では）

**vs DRLearner**:
- RCT → **DRLearner** 有利（propensity score 既知）
- 観察研究 → **DML** 有利（orthogonalization）
- 実用上の差は小さい

---

## まとめ：高次元データでの推奨ランキング

### 総合ランキング（汎用性 × 性能）

| 順位 | 手法 | 使用モデル | 適用範囲 | 理由 |
|-----|------|----------|---------|------|
| **1位** | **DRLearner** | × ElasticNetCV | ★★★★★ | **最もバランスが良い** |
| **2位** | **LinearDML** | × LassoCV | ★★★★☆ | 速度重視、大規模データ |
| **3位** | **SparseLinearDRLearner** | × Lasso（固定） | ★★★★☆ | 解釈性最優先 |
| **4位** | **CausalForestDML** | × Lasso + RF | ★★★☆☆ | 非線形効果 |
| **5位** | **DMLOrthoForest** | × Lasso（局所） | ★★☆☆☆ | 局所的変化 |

### 特性別ランキング

#### 解釈性（重要遺伝子の同定）
1. **SparseLinearDRLearner** - Lasso 係数で明確
2. **DRLearner × Lasso** - 係数 + 柔軟性
3. **LinearDML × Lasso** - 係数のみ
4. CausalForestDML - 特徴重要度のみ
5. OrthoForest - 解釈困難

#### 計算速度
1. **LinearDML × LassoCV** - 最速
2. **DRLearner × LassoCV** - 高速
3. SparseLinearDRLearner - 高速
4. CausalForestDML - 中速
5. OrthoForest - 低速

#### 予測精度（非線形効果あり）
1. **DMLOrthoForest** - 局所適応最強
2. **CausalForestDML** - 滑らか非線形
3. DRLearner × RandomForest - 柔軟
4. DRLearner × ElasticNet - 線形最適
5. LinearDML - 線形固定

#### 高次元耐性（p/n 比）
1. **DRLearner × ElasticNetCV** - p/n > 100 OK
2. **SparseLinearDRLearner** - p/n > 100 OK
3. **LinearDML × LassoCV** - p/n > 100 OK
4. CausalForestDML × Lasso - p/n < 50 推奨
5. OrthoForest - p/n < 20 推奨

---

## 実践的な選択フロー（簡易版）

### Step 1: 線形性のチェック

```python
# 簡易チェック: 線形モデルで RMSE を確認
from econml.dr import LinearDRLearner
est_linear = LinearDRLearner(...)
est_linear.fit(Y, T, X=X)
tau_linear = est_linear.effect(X)
rmse_linear = np.sqrt(np.mean((tau_linear - tau_true)**2))
```

- **RMSE < 0.3** → 線形で十分 → **DRLearner × Lasso/ElasticNet**
- **RMSE > 0.5** → 非線形を試す → **CausalForestDML**

### Step 2: 解釈性の必要度

- **係数が必要** → **SparseLinearDRLearner**
- **重要度で十分** → **DRLearner × ElasticNet**
- **不要** → **CausalForestDML**

### Step 3: 計算時間の制約

- **厳しい** → **LinearDML × LassoCV**
- **普通** → **DRLearner × ElasticNetCV**
- **余裕あり** → **OrthoForest**（非線形強い場合のみ）

---

## DemoV05 への推奨設定

### 基本設定（ベースライン）

```python
from econml.dr import DRLearner
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

est = DRLearner(
    model_propensity=LogisticRegressionCV(
        cv=5,
        penalty='elasticnet',
        l1_ratios=[0.5],  # ElasticNet with L1 ratio=0.5
        solver='saga',
        max_iter=10000,
        random_state=42
    ),
    model_regression=ElasticNetCV(
        l1_ratio=[0.5, 0.7, 0.9],  # 複数試して最適化
        cv=5,
        max_iter=10000,
        random_state=42
    ),
    model_final=ElasticNetCV(
        l1_ratio=[0.5, 0.7, 0.9],
        cv=5,
        max_iter=10000,
        random_state=42
    ),
    cv=5,
    random_state=42
)
```

### 比較実験用のバリエーション

```python
# Variant 1: 最大スパース性
est_sparse = SparseLinearDRLearner(
    model_propensity=LogisticRegressionCV(cv=5, penalty='l1', solver='saga'),
    model_regression=ElasticNetCV(cv=5),
    cv=5
)

# Variant 2: 非線形対応
est_forest = CausalForestDML(
    model_y=LassoCV(cv=5),
    model_t=LogisticRegressionCV(cv=5),
    n_estimators=4000,
    max_features='sqrt'
)

# Variant 3: 速度優先
est_fast = LinearDML(
    model_y=LassoCV(cv=3, max_iter=5000),
    model_t=LogisticRegressionCV(cv=3, max_iter=5000),
    cv=3
)
```

---

## 結論

### Q: 高次元データで DRLearner × Lasso が第一選択か？

**A: YES（ただし ElasticNet がより推奨）**

理由:
1. ✓ **ElasticNet = Lasso + Ridge** → より安定
2. ✓ **p >> n でも数値的に安定**
3. ✓ **自動特徴選択** → 解釈可能
4. ✓ **DR 理論** → propensity/outcome の保険
5. ✓ **CV で自動調整** → ハイパーパラメータ不要

### 違う場合（代替案）:

| 状況 | 推奨手法 | 理由 |
|------|---------|------|
| **解釈性最優先** | **SparseLinearDRLearner** | Lasso 固定、最もスパース |
| **速度最優先** | **LinearDML × LassoCV** | DR より 20-30% 高速 |
| **非線形効果** | **CausalForestDML** | Random Forest で柔軟 |
| **局所的変化** | **DMLOrthoForest** | リーフごとに適応 |
| **観察研究** | **LinearDML/DML** | Orthogonalization 強い |

### 最終推奨（DemoV05）

**メイン**: `DRLearner × ElasticNetCV`
**比較1**: `SparseLinearDRLearner`（解釈性）
**比較2**: `CausalForestDML`（非線形）

この3つで **線形 vs 非線形**、**スパース vs 柔軟性** のトレードオフを評価できます。
