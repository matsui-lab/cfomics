# EconML 因果推論手法の分類と違い

## 誤解を正す

❌ **誤解**: EconMLの手法は全てMeta-learner
✅ **正解**: Meta-learnerは6つのカテゴリのうちの1つに過ぎない

---

## 手法の分類体系

EconMLの手法は、**推定の数理的フレームワーク**によって6つに分類されます：

```
EconML 因果推論手法
│
├─ 1. Meta-Learners（メタ学習器）
│   └─ 既存MLを「メタ的」に使用
│
├─ 2. DML（Double Machine Learning）
│   └─ 直交化による偏り除去
│
├─ 3. DR（Doubly Robust）
│   └─ 二重頑健性による保険
│
├─ 4. OrthoForest（直交フォレスト）
│   └─ 局所適応的推定
│
├─ 5. IV（Instrumental Variables）
│   └─ 操作変数法
│
└─ 6. Policy Learning（政策学習）
    └─ 直接最適化
```

---

## 1. Meta-Learners（4種類）

### 定義
既存の機械学習手法を「そのまま」使って、アルゴリズム的にCATEを計算する手法群。

### 特徴
- ✓ **数理的な特別な仮定なし**
- ✓ 任意のMLモデルを「プラグイン」可能
- ✓ 実装が最もシンプル
- ✗ 理論保証が弱い（一部を除く）

### 手法

| 手法 | アルゴリズム | 数式 |
|------|------------|------|
| **S-Learner** | 1つのモデルでY〜T+X | τ̂(x) = μ̂(x,1) - μ̂(x,0) |
| **T-Learner** | Y(0)とY(1)を別々に | τ̂(x) = μ̂₁(x) - μ̂₀(x) |
| **X-Learner** | T-Learner + 傾向スコア加重 | τ̂(x) = e(x)τ̂₀(x) + (1-e(x))τ̂₁(x) |
| **Domain Adaptation** | ドメイン適応 | 複雑 |

### 例（T-Learner）

```python
# Step 1: 対照群でモデル学習
model_0 = RandomForestRegressor()
model_0.fit(X[T==0], Y[T==0])

# Step 2: 治療群でモデル学習
model_1 = RandomForestRegressor()
model_1.fit(X[T==1], Y[T==1])

# Step 3: 差を取る
tau_hat = model_1.predict(X_test) - model_0.predict(X_test)
```

**ポイント**: 既存のMLモデルをそのまま使うだけ。特別な数理なし。

---

## 2. DML（Double Machine Learning）（6種類）

### 定義
ニュアンスパラメータ（μ, e）をMLで推定し、**直交化（Neyman orthogonalization）**により偽相関を除去する手法。

### 数理的フレームワーク

#### ステップ1: ニュアンス推定（Cross-fitting）
```
μ̂₀(X) = E[Y|T=0,X]  # アウトカムモデル
μ̂₁(X) = E[Y|T=1,X]
ê(X) = P(T=1|X)      # 傾向スコア
```

#### ステップ2: 疑似アウトカム（Pseudo-outcome）
```
ψ̂ᵢ = (Yᵢ - μ̂₀(Xᵢ)) - (Tᵢ - ê(Xᵢ))(μ̂₁(Xᵢ) - μ̂₀(Xᵢ)) / [ê(Xᵢ)(1-ê(Xᵢ))]
```

#### ステップ3: 最終推定
```
τ̂(X) = argmin E[(ψ̂ - θ(X))²]
```

### 特徴
- ✓ **理論保証**: √n-一致性、漸近正規性
- ✓ **高次元対応**: ニュアンス推定の誤差が許容される
- ✓ Cross-fittingで過学習防止
- ✗ 実装がやや複雑

### 手法

| 手法 | 最終モデル | 用途 |
|------|-----------|------|
| **DML** | 汎用 | 柔軟 |
| **LinearDML** | 線形 | 解釈性 |
| **SparseLinearDML** | Lasso | 高次元 |
| **NonParamDML** | 非パラメトリック | 探索 |
| **KernelDML** | カーネル | 滑らか |
| **CausalForestDML** | ランダムフォレスト | 非線形 |

### 例（LinearDML）

```python
from econml.dml import LinearDML

dml = LinearDML(
    model_y=RandomForestRegressor(),  # ニュアンス: Y
    model_t=RandomForestClassifier(),  # ニュアンス: T
    discrete_treatment=True,
    cv=5  # Cross-fitting
)
dml.fit(Y, T, X=X)
tau_hat = dml.effect(X_test)
```

**ポイント**:
- ニュアンス推定と最終推定が**分離**
- 直交化により偏り除去
- 理論的に保証された推定

---

## 3. DR（Doubly Robust）（4種類）

### 定義
傾向スコアとアウトカムモデルの**どちらか一方が正しければ一致推定量**となる性質を持つ。

### 数理的フレームワーク

#### DR推定量（Kennedy, 2023）
```
τ̂_DR(x) = 1/n Σᵢ [
    (Tᵢ(Yᵢ - μ̂₁(Xᵢ))) / ê(Xᵢ)
  - ((1-Tᵢ)(Yᵢ - μ̂₀(Xᵢ))) / (1-ê(Xᵢ))
  + μ̂₁(x) - μ̂₀(x)
]
```

### 二重頑健性の意味

| 状況 | 一致性 | 理由 |
|------|--------|------|
| μ̂ 正しい, ê 正しい | ✓ | 両方正しい |
| μ̂ 正しい, ê 間違い | ✓ | アウトカムモデルで救済 |
| μ̂ 間違い, ê 正しい | ✓ | 傾向スコアで救済 |
| μ̂ 間違い, ê 間違い | ✗ | 両方必要 |

### 特徴
- ✓ **頑健性最高**: 2つのモデルで保険
- ✓ RCTでは傾向スコアが既知（0.5）なので更に頑健
- ✓ DMLと似た理論保証
- ✗ 計算コストやや高

### 手法

| 手法 | 最終モデル | 用途 |
|------|-----------|------|
| **DRLearner** | 汎用回帰 | 標準（**demov04使用**） |
| **LinearDRLearner** | 線形 | シンプル |
| **SparseLinearDRLearner** | Lasso | 高次元 |
| **ForestDRLearner** | ランダムフォレスト | 非線形 |

### 例（DRLearner）

```python
from econml.dr import DRLearner

dr = DRLearner(
    model_propensity=LogisticRegressionCV(),  # 傾向スコア
    model_regression=RidgeCV(),               # アウトカム
    model_final=LassoCV(),                    # 最終CATE
    cv=5
)
dr.fit(Y, T, X=X)
tau_hat = dr.effect(X_test)
```

**ポイント**: DMLと似ているが、異なるスコア関数を使用。

---

## 4. OrthoForest（2種類）

### 定義
ニュアンスパラメータを**局所的に**推定し、データ駆動的にXの空間を分割。

### 通常のDML/DRとの違い

| 特性 | DML/DR | OrthoForest |
|------|--------|------------|
| ニュアンス推定 | **グローバル** | **局所** |
| 適応性 | 固定 | データ駆動 |
| 計算コスト | 低 | 高 |

### アルゴリズム

```
1. データをフォレストで分割
2. 各葉ノード（局所領域）で：
   - その領域のデータのみ使用
   - ニュアンスパラメータを推定
   - CATEを推定
3. 局所推定を統合
```

### 手法

| 手法 | ベース | 特徴 |
|------|-------|------|
| **DMLOrthoForest** | DML | 局所DML |
| **DROrthoForest** | DR | 局所DR |

### 例

```python
from econml.orf import DMLOrthoForest

orf = DMLOrthoForest(
    n_trees=4000,
    min_leaf_size=10,
    model_Y=RandomForestRegressor(),
    model_T=RandomForestClassifier()
)
```

**ポイント**: 局所的に異なるCATEパターンを自動検出。

---

## 5. IV（Instrumental Variables）

### 定義
未観測交絡がある場合に、**操作変数（Z）**を用いて因果効果を推定。

### 問題設定

```
Z → T → Y
     ↑
     X (観測)
     U (未観測交絡)
```

### 操作変数の条件

1. **関連性**: Z は T に影響する
2. **除外制約**: Z は Y に T を通じてのみ影響
3. **独立性**: Z は U と独立

### 他の手法との違い

| 手法 | 前提 | 必要なもの |
|------|------|----------|
| Meta/DML/DR | **無交絡** | X のみ |
| **IV** | **交絡OK** | Z（操作変数） |

### 例

```python
from econml.iv.dml import DMLIV

dmliv = DMLIV(
    model_Y_X=RandomForestRegressor(),
    model_T_X=RandomForestRegressor(),
    model_T_XZ=RandomForestClassifier(),
    model_final=LassoCV()
)
# Z: 操作変数（例: ランダム化の推奨）
dmliv.fit(Y, T, Z=Z, X=X)
```

**ポイント**: **完全に異なる問題設定**。観察研究で必須。

---

## 6. Policy Learning（4種類）

### 定義
CATEを推定するのではなく、**最適な治療割り当てルール**を直接学習。

### 目的関数

```
π*(x) = argmax E[Y | do(T=t), X=x]
         t∈{0,1}
```

### 他の手法との違い

| 手法 | 出力 | 目的 |
|------|------|------|
| Meta/DML/DR | **τ̂(x)** (連続値) | 効果の推定 |
| **Policy Learning** | **π(x) ∈ {0,1}** (離散) | 意思決定 |

### 手法

| 手法 | アルゴリズム | 特徴 |
|------|------------|------|
| **PolicyTree** | 決定木 | 解釈可能 |
| **PolicyForest** | ランダムフォレスト | 高精度 |
| **DRPolicyTree** | DR + 決定木 | 頑健 |
| **DRPolicyForest** | DR + フォレスト | 頑健 + 高精度 |

### 例

```python
from econml.policy import DRPolicyTree

policy = DRPolicyTree(max_depth=3)
policy.fit(Y, T, X=X)

# 出力: 治療すべきか否か（0 or 1）
treatment_decision = policy.predict(X_test)
```

**ポイント**: CATEではなく、直接アクションを決定。

---

## 分類の比較表

### 数理的フレームワーク

| カテゴリ | 理論的基礎 | 頑健性 | 高次元 | 計算コスト |
|---------|----------|--------|--------|-----------|
| **Meta-Learners** | なし（アルゴリズム的） | △ | △ | 低 |
| **DML** | 直交化 | ○ | ◎ | 中 |
| **DR** | 二重頑健性 | ◎ | ◎ | 中 |
| **OrthoForest** | 局所直交化 | ○ | ○ | 高 |
| **IV** | 操作変数法 | ○ | △ | 中 |
| **Policy** | 直接最適化 | - | - | 中 |

### 使用場面

| カテゴリ | 使用場面 |
|---------|---------|
| **Meta-Learners** | クイック分析、ベースライン |
| **DML** | 標準的CATE推定、高次元データ |
| **DR** | 頑健性重視、RCT、高次元 |
| **OrthoForest** | 局所的異質性探索 |
| **IV** | 観察研究、未観測交絡 |
| **Policy** | 意思決定、個別化医療 |

---

## 主要な違いのまとめ

### 1. Meta-Learners vs DML/DR

| 項目 | Meta-Learners | DML/DR |
|------|--------------|--------|
| **数理的フレームワーク** | なし | あり（直交化/DR） |
| **理論保証** | 弱い | 強い（√n-一致性） |
| **ニュアンス推定** | 直接使用 | 影響を除去 |
| **実装** | シンプル | やや複雑 |
| **高次元対応** | 限定的 | 優秀 |

**結論**: DML/DRの方が理論的に洗練されている。

### 2. DML vs DR

| 項目 | DML | DR |
|------|-----|-----|
| **理論的基礎** | 直交化 | 二重頑健性 |
| **頑健性** | ニュアンス推定誤差に強い | モデル誤指定に強い |
| **RCTでの利点** | 標準 | 傾向スコア既知で更に頑健 |
| **実装** | やや複雑 | やや複雑 |

**結論**: RCTではDRがわずかに有利、観察研究では同等。

### 3. グローバル vs 局所

| 項目 | DML/DR | OrthoForest |
|------|--------|------------|
| **推定範囲** | グローバル | 局所 |
| **適応性** | 固定 | データ駆動 |
| **計算コスト** | 低 | 高 |
| **用途** | 標準 | 異質性探索 |

**結論**: 標準はDML/DR、探索的にはOrthoForest。

---

## DemoV04での選択

### 現在の実装: DRLearner

**理由**:
1. ✓ 線形DGPに最適（LinearモデルをFinalに使用）
2. ✓ 二重頑健性でRCTに最適
3. ✓ 高次元対応（ElasticNetCV使用可能）
4. ✓ 理論保証あり

### なぜMeta-Learnerではないのか？

- Meta-Learner（例: T-Learner）でも推定可能
- しかし**DRLearnerの方が理論的に優れている**：
  - 直交化により偏り除去
  - 二重頑健性で保険
  - 高次元データへの理論保証

### なぜDMLではなくDRなのか？

- DMLでも同等の性能
- DRは**RCTで傾向スコアが既知（0.5）**なので更に頑健
- 実装上の違いは小さい

---

## まとめ

### 誤解を正す

❌ **全てMeta-learner**
✅ **6つの異なる数理的フレームワーク**

### 分類体系

1. **Meta-Learners**: アルゴリズム的アプローチ
2. **DML**: 直交化による偏り除去
3. **DR**: 二重頑健性
4. **OrthoForest**: 局所適応
5. **IV**: 操作変数法（完全に別問題）
6. **Policy**: 直接最適化（出力が異なる）

### 実用的には

- **標準**: DRLearner または LinearDML
- **高次元**: SparseLinearDRLearner
- **非線形**: CausalForestDML
- **観察研究**: DMLIV
- **意思決定**: DRPolicyTree

DemoV04での**DRLearner**選択は理論的に適切。
