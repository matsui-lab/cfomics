# 最終モデルの変数選択メカニズム

## TL;DR

**重要**: 最終モデル（model_final）には **全ての X 変数が渡される**

- ニュアンスモデルで特徴選択されても、**model_final は全変数を受け取る**
- model_final 自身が変数選択を行う（Lasso なら自動、RF ならランダムサンプリング）
- ニュアンス選択 ≠ 最終モデル選択（独立）

---

## 1. 基本的な流れ

### DRLearner の例（n=500, p=20,000）

```python
from econml.dr import DRLearner
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

X.shape  # (500, 20000) - 20,000遺伝子

est = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),  # ニュアンス1
    model_regression=ElasticNetCV(cv=5),          # ニュアンス2
    model_final=ElasticNetCV(cv=5)                # ターゲット
)

est.fit(Y, T, X=X)
```

### 内部の動作（詳細）

#### Step 1: ニュアンス推定（Propensity）

```python
# model_propensity: P(T=1|X) を推定
model_prop = LogisticRegressionCV(cv=5)
model_prop.fit(X, T)  # 入力: 全 20,000 遺伝子

# Lasso により自動で特徴選択
print(f"Propensity で選択された遺伝子: {np.sum(model_prop.coef_ != 0)}")
# 出力例: Propensity で選択された遺伝子: 23

selected_genes_prop = np.where(model_prop.coef_[0] != 0)[0]
print(f"選択遺伝子ID: {selected_genes_prop}")
# 出力例: 選択遺伝子ID: [45, 234, 567, 1024, ...]
```

#### Step 2: ニュアンス推定（Regression）

```python
# model_regression: E[Y|T,X] を推定
model_reg = ElasticNetCV(cv=5)

# 治療群で学習
X_treated = X[T == 1]
Y_treated = Y[T == 1]
model_reg_1 = ElasticNetCV(cv=5)
model_reg_1.fit(X_treated, Y_treated)  # 入力: 全 20,000 遺伝子

print(f"Regression(T=1) で選択された遺伝子: {np.sum(model_reg_1.coef_ != 0)}")
# 出力例: Regression(T=1) で選択された遺伝子: 67

selected_genes_reg1 = np.where(model_reg_1.coef_ != 0)[0]

# 対照群で学習
X_control = X[T == 0]
Y_control = Y[T == 0]
model_reg_0 = ElasticNetCV(cv=5)
model_reg_0.fit(X_control, Y_control)  # 入力: 全 20,000 遺伝子

print(f"Regression(T=0) で選択された遺伝子: {np.sum(model_reg_0.coef_ != 0)}")
# 出力例: Regression(T=0) で選択された遺伝子: 58
```

**重要**: ここまでで、異なるモデルが異なる遺伝子を選択
- Propensity: 23遺伝子
- Regression(T=1): 67遺伝子
- Regression(T=0): 58遺伝子

#### Step 3: 疑似アウトカム計算

```python
# 疑似アウトカムの計算（全サンプルで）
p_hat = model_prop.predict_proba(X)[:, 1]  # 入力: 全 20,000 遺伝子
mu0_hat = model_reg_0.predict(X)            # 入力: 全 20,000 遺伝子
mu1_hat = model_reg_1.predict(X)            # 入力: 全 20,000 遺伝子

# DR 疑似アウトカム
pseudo_outcome = (
    (T - p_hat) / (p_hat * (1 - p_hat)) *
    (Y - (T * mu1_hat + (1 - T) * mu0_hat))
)
```

#### Step 4: 最終モデル学習（重要！）

```python
# model_final: τ(X) を推定
model_final = ElasticNetCV(cv=5)

# 疑似アウトカムを使って学習
model_final.fit(X, pseudo_outcome)  # 入力: 全 20,000 遺伝子！
                                     #       ニュアンスの選択遺伝子ではない

print(f"Final model で選択された遺伝子: {np.sum(model_final.coef_ != 0)}")
# 出力例: Final model で選択された遺伝子: 31

selected_genes_final = np.where(model_final.coef_ != 0)[0]
print(f"選択遺伝子ID: {selected_genes_final}")
# 出力例: 選択遺伝子ID: [0, 12, 45, 89, 234, ...]
```

**キーポイント**:
- model_final への入力は **X（全 20,000 遺伝子）**
- ニュアンスモデルで選ばれた遺伝子に限定されない
- model_final が独自に特徴選択（31遺伝子）

---

## 2. なぜ全変数を渡すのか？

### 理由1: ニュアンスと治療効果は別の現象

**ニュアンスで重要な変数** ≠ **治療効果で重要な変数**

**例**:
```python
# 真のモデル
Y = 100 + 5*gene1 + 3*gene2 + T * (2*gene3 + 1*gene4) + noise
                              ↑治療効果

# ニュアンス E[Y|X]
E[Y|X] = 100 + 5*gene1 + 3*gene2 + 0.5*T*(2*gene3 + 1*gene4)
         ↑gene1, gene2 が重要（主効果が大きい）

# 治療効果 τ(X)
τ(X) = 2*gene3 + 1*gene4
       ↑gene3, gene4 が重要（交互作用）
```

**結論**:
- ニュアンスモデルは gene1, gene2 を選択
- でも治療効果には gene3, gene4 が重要
- **だから model_final には全変数を渡すべき**

### 理由2: 理論的な要請（DML/DR）

**DML/DR 理論**:
- ニュアンス推定と治療効果推定は **分離されるべき**
- ニュアンスの選択に依存すると推定量の性質が変わる
- 治療効果は独立に推定すべき

### 理由3: 実用上の利点

**全変数を渡す利点**:
- ✓ ニュアンスで見逃した重要変数も拾える
- ✓ model_final の選択能力を信頼
- ✓ 柔軟性が高い

---

## 3. 各手法での変数選択

### DRLearner × Lasso/ElasticNet

```python
est = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=ElasticNetCV(cv=5),
    model_final=ElasticNetCV(cv=5)
)

est.fit(Y, T, X=X)  # X.shape = (500, 20000)
```

**model_final への入力**: 全 20,000 遺伝子

**model_final の選択**:
```python
# ElasticNet が自動で重要遺伝子を選択
coef = est.model_final.coef_
selected = np.where(coef != 0)[0]

print(f"最終モデルが使用した遺伝子数: {len(selected)}")
# 出力例: 最終モデルが使用した遺伝子数: 31

print(f"遺伝子リスト: {selected}")
# 出力例: 遺伝子リスト: [0, 12, 45, 89, 234, ...]
```

**特徴**:
- L1 正則化で自動的に係数をゼロに
- 重要遺伝子のみが非ゼロ係数
- 解釈可能（どの遺伝子が重要か明確）

### ForestDRLearner × GRF

```python
est = ForestDRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=ElasticNetCV(cv=5),
    n_estimators=2000,
    max_features='sqrt'  # √20000 ≈ 141
)

est.fit(Y, T, X=X)  # X.shape = (500, 20000)
```

**model_final への入力**: 全 20,000 遺伝子

**model_final の選択**:
- 各ノードの分岐で **ランダムに √20000 ≈ 141 遺伝子をサンプリング**
- 141遺伝子の中から最良の分岐変数を選択
- 全ツリーを通じて全 20,000 遺伝子が考慮される

```python
# 特徴重要度の取得
importance = est.feature_importances_

# Top 20 重要遺伝子
top_genes = np.argsort(importance)[-20:]
print(f"Top 20 重要遺伝子: {top_genes}")
# 出力例: Top 20 重要遺伝子: [0, 12, 23, 45, 67, ...]
```

**特徴**:
- ランダムサンプリングで全変数を探索
- 特徴重要度で寄与度を評価
- 非線形関係も捉える

### SparseLinearDRLearner

```python
est = SparseLinearDRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=ElasticNetCV(cv=5)
    # model_final は Lasso 固定
)

est.fit(Y, T, X=X)  # X.shape = (500, 20000)
```

**model_final への入力**: 全 20,000 遺伝子

**model_final の選択**:
```python
# Lasso が最もスパースに選択
coef = est.coef_
selected = np.where(coef != 0)[0]

print(f"最終モデルが使用した遺伝子数: {len(selected)}")
# 出力例: 最終モデルが使用した遺伝子数: 18（最もスパース）
```

---

## 4. 誤解: ニュアンスの選択変数のみ使う？

### ❌ 間違った理解

```python
# これは EconML の実装ではない

# Step 1: ニュアンスで特徴選択
model_reg = ElasticNetCV(cv=5)
model_reg.fit(X, Y)
selected_genes = np.where(model_reg.coef_ != 0)[0]  # 50遺伝子

# Step 2: 選択された特徴のみで最終モデル
X_selected = X[:, selected_genes]  # (500, 50)
model_final.fit(X_selected, pseudo_outcome)  # ❌ これは行われない
```

**問題**:
- ニュアンスの選択に依存してしまう
- 治療効果に重要だがニュアンスで選ばれなかった変数が消える
- 理論的保証が失われる

### ✅ 正しい実装（EconML）

```python
# EconML の実際の実装

# Step 1: ニュアンスで特徴選択（内部的に）
model_reg = ElasticNetCV(cv=5)
model_reg.fit(X, Y)  # 50遺伝子を選択（内部処理）

# Step 2: 疑似アウトカム計算（全変数で予測）
pseudo_outcome = ...

# Step 3: 最終モデルは全変数で学習
model_final.fit(X, pseudo_outcome)  # ✅ 全 20,000 遺伝子
```

---

## 5. 実験：変数選択の可視化

### データ生成

```python
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)
n, p = 1000, 100

X = np.random.randn(n, p)
T = np.random.binomial(1, 0.5, n)

# 真のモデル（変数を分離）
# ニュアンスに重要: gene 0-4
# 治療効果に重要: gene 5-9
Y0 = 100 + 10*X[:, 0] + 8*X[:, 1] + 6*X[:, 2] + 4*X[:, 3] + 2*X[:, 4]
tau_true = 5*X[:, 5] + 4*X[:, 6] + 3*X[:, 7] + 2*X[:, 8] + 1*X[:, 9]

Y1 = Y0 + tau_true
Y = Y0 * (1 - T) + Y1 * T + np.random.randn(n)
```

### DRLearner で推定

```python
from econml.dr import DRLearner
from sklearn.linear_model import LogisticRegressionCV, LassoCV

est = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=LassoCV(cv=5),
    model_final=LassoCV(cv=5)
)

est.fit(Y, T, X=X)
```

### 選択変数の確認

```python
# ニュアンス（regression）で選択された変数
# 注: 内部モデルへのアクセスは実装依存
# 以下は概念的な例

# 最終モデルで選択された変数
coef_final = est.model_final.coef_
selected_final = np.where(np.abs(coef_final) > 0.01)[0]

print("真に重要な変数（治療効果）: [5, 6, 7, 8, 9]")
print(f"最終モデルが選択した変数: {selected_final}")
# 出力例: 最終モデルが選択した変数: [5, 6, 7, 8, 9, 12, 23]

# 可視化
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.bar(range(p), np.abs(coef_final))
plt.axvline(4.5, color='r', linestyle='--', label='Nuisance/Effect boundary')
plt.xlabel('Gene ID')
plt.ylabel('|Coefficient|')
plt.title('Final Model Coefficients')
plt.legend()

plt.subplot(1, 2, 2)
true_coef = np.zeros(p)
true_coef[5:10] = [5, 4, 3, 2, 1]
plt.scatter(true_coef, np.abs(coef_final), alpha=0.6)
plt.plot([0, 5], [0, 5], 'r--', label='Perfect')
plt.xlabel('True Coefficient')
plt.ylabel('Estimated |Coefficient|')
plt.title('True vs Estimated')
plt.legend()

plt.tight_layout()
plt.savefig('final_model_variable_selection.png', dpi=150)
```

---

## 6. 特殊ケース: 事前に次元削減する場合

### パターン1: PCA で次元削減してから推定

```python
from sklearn.decomposition import PCA

# 事前に PCA（20,000 → 100）
pca = PCA(n_components=100)
X_reduced = pca.fit_transform(X)  # (500, 100)

# 削減されたデータで DRLearner
est = DRLearner(...)
est.fit(Y, T, X=X_reduced)
```

**この場合**:
- model_final への入力は 100 主成分
- 解釈性は低下（主成分は遺伝子の線形結合）

### パターン2: SelectKBest で事前選択

```python
from sklearn.feature_selection import SelectKBest, f_regression

# 事前に Top 500 遺伝子を選択
selector = SelectKBest(f_regression, k=500)
X_selected = selector.fit_transform(X, Y)  # (500, 500)

# 選択されたデータで DRLearner
est = DRLearner(...)
est.fit(Y, T, X=X_selected)
```

**この場合**:
- model_final への入力は 500 遺伝子
- Y との相関で選択（治療効果とは無関係）
- **推奨されない**（治療効果に重要な変数を見逃す可能性）

---

## 7. まとめ

### Q: 最終モデルに使う変数はどう選択される？

**A: model_final には全ての X 変数が渡され、model_final 自身が選択する**

### 変数選択のフロー

```
入力データ X (500, 20000)
    ↓
Step 1: ニュアンス推定
  - model_propensity: 全 20,000 を使用 → 23遺伝子を選択（内部）
  - model_regression: 全 20,000 を使用 → 67遺伝子を選択（内部）
    ↓
Step 2: 疑似アウトカム計算
  - 全 20,000 遺伝子で予測値を計算
    ↓
Step 3: 最終モデル学習
  - **入力: 全 20,000 遺伝子**
  - model_final が独自に選択 → 31遺伝子を選択
    ↓
結果: τ(X) の推定（31遺伝子を使用）
```

### 各モデルタイプでの選択方法

| model_final | 選択方法 | 選択数 |
|------------|---------|--------|
| **LassoCV** | L1正則化 | 最小（10-50） |
| **ElasticNetCV** | L1+L2正則化 | 小〜中（20-100） |
| **RandomForest/GRF** | ランダムサンプリング | 全変数（重要度で評価） |

### キーポイント

1. **全変数が model_final に渡される**（ニュアンスの選択に依存しない）
2. **model_final が独自に選択**（Lasso/RF の能力を活用）
3. **ニュアンス変数 ≠ 治療効果変数**（異なる現象）
4. **事前選択は推奨されない**（DML/DR 理論に反する）
