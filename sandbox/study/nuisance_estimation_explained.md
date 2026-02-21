# ニュアンス推定（Nuisance Estimation）の詳細解説

## TL;DR

**ニュアンス推定 = 本当に知りたいもの（治療効果）を推定するために必要な「邪魔なパラメータ」の推定**

- **本当に知りたいもの**: τ(X) = 治療効果（CATE）
- **ニュアンス（邪魔者）**: E[Y|X], E[T|X], P(T=1|X) など
- これらを正確に推定しないと、治療効果の推定が偏る

---

## 1. ニュアンス（Nuisance）とは？

### 統計学の定義

**Nuisance Parameter（ニュアンスパラメータ）**:
- 直接興味はないが、推定のために必要なパラメータ
- "nuisance" = 「邪魔者」「厄介なもの」

### 因果推論での具体例

**目的**: 治療効果 τ(X) = E[Y¹|X] - E[Y⁰|X] を推定したい

**問題**: Y¹ と Y⁰ を直接観測できない（counterfactual）

**解決**: 以下の補助的なパラメータを推定する（これが nuisance）:
1. **E[Y|X, T]** - アウトカムモデル
2. **E[T|X]** または **P(T=1|X)** - 傾向スコア

これらは「知りたいわけではないが、治療効果推定に必要」→ ニュアンス

---

## 2. 各手法でのニュアンス

### DRLearner（Doubly Robust）のニュアンス

**3段階推定**:

```python
from econml.dr import DRLearner

est = DRLearner(
    model_propensity=...,  # ニュアンス1: P(T=1|X)
    model_regression=...,  # ニュアンス2: E[Y|T,X]
    model_final=...        # ターゲット: τ(X)
)
```

| モデル | 推定対象 | 種類 | 目的 |
|--------|---------|------|------|
| **model_propensity** | P(T=1\|X) | **ニュアンス** | 傾向スコア |
| **model_regression** | E[Y\|T,X] | **ニュアンス** | アウトカム予測 |
| **model_final** | τ(X) | **ターゲット** | **治療効果（本命）** |

**なぜニュアンスが必要か？**

治療効果の推定式（疑似アウトカム）:
```
τ(X) を推定するには、まず以下を計算:

疑似アウトカム = [T - P(T=1|X)] / P(T=1|X) * [Y - E[Y|T,X]]
                 ↑ニュアンス1        ↑ニュアンス2

この疑似アウトカムを使って τ(X) を推定
```

### DML（Double Machine Learning）のニュアンス

**3段階推定**:

```python
from econml.dml import DML

est = DML(
    model_y=...,     # ニュアンス1: E[Y|X]
    model_t=...,     # ニュアンス2: E[T|X]
    model_final=...  # ターゲット: τ(X)
)
```

| モデル | 推定対象 | 種類 | 目的 |
|--------|---------|------|------|
| **model_y** | E[Y\|X] | **ニュアンス** | Y の期待値 |
| **model_t** | E[T\|X] | **ニュアンス** | T の期待値 |
| **model_final** | τ(X) | **ターゲット** | **治療効果（本命）** |

**DML の推定式**:
```
Step 1: ニュアンス推定
  Y_res = Y - E[Y|X]  ← model_y で推定
  T_res = T - E[T|X]  ← model_t で推定

Step 2: 治療効果推定
  τ(X) = E[Y_res | T_res, X]  ← model_final で推定
```

### CausalForestDML のニュアンス

```python
from econml.dml import CausalForestDML

est = CausalForestDML(
    model_y=LassoCV(cv=5),               # ニュアンス1
    model_t=LogisticRegressionCV(cv=5),  # ニュアンス2
    # model_final は GRF（自動、変更不可）
)
```

| モデル | 推定対象 | 種類 | カスタマイズ可能 |
|--------|---------|------|----------------|
| **model_y** | E[Y\|X] | **ニュアンス** | ✓（Lasso等） |
| **model_t** | E[T\|X] | **ニュアンス** | ✓（Logistic等） |
| **model_final** | τ(X) | **ターゲット** | ×（GRF固定） |

---

## 3. なぜ「ニュアンス」と呼ぶのか？

### 興味の対象ではない

**研究者の質問**:
- 「この薬は血圧を下げるか？」
- 「遺伝子発現が高い患者で効果は大きいか？」

**知りたいこと**: τ(X) = 治療効果

**知りたくないこと（でも必要）**:
- P(T=1|X) = 治療を受ける確率
- E[Y|X] = アウトカムの期待値

これらは「治療効果を推定するための補助的情報」→ ニュアンス

### 推定誤差の影響

**ニュアンス推定の誤差** → **治療効果推定の誤差**

```
ニュアンスが間違っている場合:

P(T=1|X) の推定が悪い
    ↓
疑似アウトカムが歪む
    ↓
τ(X) の推定も偏る
```

**だから厄介（nuisance）**: 興味はないが、正確に推定しないといけない

---

## 4. 高次元データでのニュアンス推定

### 問題: ニュアンス推定が難しい

**RNA-seq の例**:
- n = 500 サンプル
- p = 20,000 遺伝子
- ニュアンス E[Y|X] を推定 → 20,000次元の回帰問題

**通常の回帰では失敗**:
```python
# これは失敗する
from sklearn.linear_model import LinearRegression
model_y = LinearRegression()
model_y.fit(X, Y)  # X.shape = (500, 20000) → エラー
```

### 解決策: スパース推定

**Lasso/ElasticNet を使う**:

```python
from sklearn.linear_model import LassoCV

# 高次元ニュアンス推定
model_y = LassoCV(cv=5, max_iter=10000)
model_y.fit(X, Y)  # X.shape = (500, 20000) → OK

# 20,000 → 50遺伝子に自動削減
print(f"選択遺伝子数: {np.sum(model_y.coef_ != 0)}")
# 出力例: 選択遺伝子数: 52
```

**これが CausalForestDML で Lasso を使う理由**:

```python
cf_dml = CausalForestDML(
    model_y=LassoCV(cv=5),  # 高次元ニュアンス推定に最適
    model_t=LogisticRegressionCV(cv=5),
    # 最終モデル（GRF）は50遺伝子で推定すれば良い
)
```

---

## 5. ニュアンス推定の重要性

### DR（Doubly Robust）理論

**DRLearner の利点**:

```
治療効果が一致推定となる条件（どちらか一方でOK）:

条件1: P(T=1|X) が正しい OR
条件2: E[Y|T,X] が正しい

→ 両方間違っていなければ良い（保険がある）
```

これが「Doubly Robust（二重にロバスト）」の意味

### DML（Neyman Orthogonalization）理論

**DML の利点**:

```
ニュアンス推定の誤差の「1次の項」を除去

通常: τの推定誤差 ∝ ニュアンスの誤差
DML:  τの推定誤差 ∝ (ニュアンスの誤差)²

→ ニュアンスが少し間違っていても大丈夫
```

---

## 6. 実例：ニュアンス推定の可視化

### データ生成

```python
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)
n = 1000
p = 5

# 共変量
X = np.random.randn(n, p)

# ニュアンス1: E[Y|X]（線形）
E_Y_given_X = 100 + 2*X[:, 0] + 1.5*X[:, 1]

# ニュアンス2: P(T=1|X)（ロジスティック）
logit_p = 0.5*X[:, 0] - 0.3*X[:, 1]
P_T_given_X = 1 / (1 + np.exp(-logit_p))
T = np.random.binomial(1, P_T_given_X)

# 治療効果（ターゲット）
tau_true = 5 + 3*X[:, 0]  # 遺伝子0に依存

# アウトカム生成
Y0 = E_Y_given_X + np.random.randn(n)
Y1 = Y0 + tau_true
Y = Y0 * (1 - T) + Y1 * T
```

### ニュアンス推定の確認

```python
from sklearn.linear_model import LassoCV, LogisticRegressionCV

# ニュアンス1: E[Y|X] 推定
model_y = LassoCV(cv=5)
model_y.fit(X, Y)
E_Y_hat = model_y.predict(X)

# ニュアンス2: P(T=1|X) 推定
model_t = LogisticRegressionCV(cv=5)
model_t.fit(X, T)
P_T_hat = model_t.predict_proba(X)[:, 1]

# 可視化
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# ニュアンス1
axes[0].scatter(E_Y_given_X, E_Y_hat, alpha=0.5, s=10)
axes[0].plot([90, 110], [90, 110], 'r--', label='Perfect')
axes[0].set_xlabel('True E[Y|X]')
axes[0].set_ylabel('Estimated E[Y|X]')
axes[0].set_title('Nuisance 1: Outcome Model')
axes[0].legend()

# ニュアンス2
axes[1].scatter(P_T_given_X, P_T_hat, alpha=0.5, s=10)
axes[1].plot([0, 1], [0, 1], 'r--', label='Perfect')
axes[1].set_xlabel('True P(T=1|X)')
axes[1].set_ylabel('Estimated P(T=1|X)')
axes[1].set_title('Nuisance 2: Propensity Score')
axes[1].legend()

plt.tight_layout()
plt.savefig('nuisance_estimation_quality.png', dpi=150)
print("ニュアンス推定の精度を可視化しました")
```

### 治療効果推定（ニュアンス使用）

```python
from econml.dr import DRLearner

# DRLearner（ニュアンス使用）
est = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),  # ニュアンス2
    model_regression=LassoCV(cv=5),               # ニュアンス1
    model_final=LassoCV(cv=5)                     # ターゲット
)

est.fit(Y, T, X=X)
tau_hat = est.effect(X)

# 精度確認
rmse = np.sqrt(np.mean((tau_hat - tau_true)**2))
r2 = 1 - np.sum((tau_hat - tau_true)**2) / np.sum((tau_true - np.mean(tau_true))**2)

print(f"\n治療効果推定:")
print(f"  RMSE: {rmse:.3f}")
print(f"  R²: {r2:.3f}")

# 可視化
plt.figure(figsize=(6, 6))
plt.scatter(tau_true, tau_hat, alpha=0.5, s=10)
plt.plot([0, 15], [0, 15], 'r--', label='Perfect')
plt.xlabel('True CATE')
plt.ylabel('Estimated CATE')
plt.title(f'Treatment Effect Estimation (RMSE={rmse:.2f})')
plt.legend()
plt.savefig('cate_estimation_with_nuisance.png', dpi=150)
```

---

## 7. ニュアンス推定が失敗する例

### 悪いニュアンス推定

```python
# ニュアンスを Random で推定（最悪）
from sklearn.dummy import DummyRegressor, DummyClassifier

est_bad = DRLearner(
    model_propensity=DummyClassifier(strategy='prior'),  # 全員0.5
    model_regression=DummyRegressor(strategy='mean'),    # 全員平均値
    model_final=LassoCV(cv=5)
)

est_bad.fit(Y, T, X=X)
tau_bad = est_bad.effect(X)

rmse_bad = np.sqrt(np.mean((tau_bad - tau_true)**2))
print(f"\n悪いニュアンス推定:")
print(f"  RMSE: {rmse_bad:.3f}")  # 非常に悪い
```

**期待される出力**:
```
悪いニュアンス推定:
  RMSE: 4.821  # 良いニュアンスなら 0.3程度
```

### 良いニュアンス推定

```python
# ニュアンスを適切に推定
est_good = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=LassoCV(cv=5),
    model_final=LassoCV(cv=5)
)

est_good.fit(Y, T, X=X)
tau_good = est_good.effect(X)

rmse_good = np.sqrt(np.mean((tau_good - tau_true)**2))
print(f"\n良いニュアンス推定:")
print(f"  RMSE: {rmse_good:.3f}")  # 良好
```

**期待される出力**:
```
良いニュアンス推定:
  RMSE: 0.287
```

---

## 8. まとめ

### ニュアンス推定とは？

**定義**: 治療効果（ターゲット）を推定するために必要な補助的パラメータの推定

| 用語 | 意味 | 例 |
|------|------|-----|
| **Target Parameter** | 本当に知りたいパラメータ | τ(X) = 治療効果 |
| **Nuisance Parameter** | 補助的パラメータ | E[Y\|X], P(T=1\|X) |
| **Nuisance Estimation** | ニュアンスを推定すること | Lasso で E[Y\|X] 推定 |

### 各手法のニュアンス

| 手法 | ニュアンス1 | ニュアンス2 | ターゲット |
|------|-----------|-----------|-----------|
| **DRLearner** | P(T=1\|X) | E[Y\|T,X] | τ(X) |
| **DML** | E[Y\|X] | E[T\|X] | τ(X) |
| **CausalForestDML** | E[Y\|X] | E[T\|X] | τ(X) |

### 高次元データでの推奨

**ニュアンス推定に Lasso/ElasticNet を使う**:

```python
# 推奨設定
est = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),  # ニュアンス
    model_regression=ElasticNetCV(cv=5),          # ニュアンス
    model_final=ElasticNetCV(cv=5)                # ターゲット
)
```

理由:
1. ✓ p >> n でも推定可能
2. ✓ 自動特徴選択（20,000 → 50遺伝子）
3. ✓ ニュアンス推定の精度向上
4. ✓ 結果として治療効果推定も改善

### キーポイント

1. **ニュアンス = 邪魔者だが必要**: 直接興味はないが推定に不可欠
2. **ニュアンスの精度が重要**: 誤差が治療効果推定に伝播
3. **高次元ではスパース推定**: Lasso/ElasticNet で自動削減
4. **DR/DML 理論**: ニュアンス誤差に対してロバスト
