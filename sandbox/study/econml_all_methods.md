# EconML 因果推論手法 完全一覧

**EconML version: 0.16.0**

EconMLは、機械学習と計量経済学を組み合わせた因果推論のためのPythonパッケージです。
異質な治療効果（Heterogeneous Treatment Effects / CATE）の推定に特化しています。

---

## 手法の分類

| カテゴリ | 手法数 | 特徴 |
|---------|-------|------|
| **Meta-Learners** | 4種類 | シンプル、解釈しやすい |
| **DML (Double ML)** | 6種類 | 高次元対応、理論保証 |
| **DR (Doubly Robust)** | 4種類 | 頑健性が高い |
| **OrthoForest** | 2種類 | 局所適応的 |
| **IV (操作変数法)** | 複数 | 内生性に対応 |
| **Policy Learning** | 4種類 | 意思決定最適化 |

**合計: 20種類以上**

---

## 1. Meta-Learners（メタ学習器）

### 概要
既存のML手法を「メタ的」に使用してCATEを推定する手法群。
最もシンプルで実装しやすい。

### 利用可能な手法

| 手法 | クラス名 | 説明 | 使用場面 |
|------|---------|------|---------|
| **S-Learner** | `SLearner` | 単一モデルでY(0)とY(1)を同時推定 | 最もシンプル、ベースライン |
| **T-Learner** | `TLearner` | Y(0)とY(1)を別々のモデルで推定 | 中サンプルサイズ |
| **X-Learner** | `XLearner` | T-Learnerを改良、傾向スコア加重 | 小サンプル、不均衡データ |
| **Domain Adaptation Learner** | `DomainAdaptationLearner` | ドメイン適応を利用 | 処置群/対照群の分布が異なる |

#### 詳細

**S-Learner**
```python
from econml.metalearners import SLearner
from sklearn.ensemble import RandomForestRegressor

slearner = SLearner(overall_model=RandomForestRegressor())
slearner.fit(Y, T, X=X)
cate = slearner.effect(X_test)
```
- **メリット**: 最もシンプル、全データ使用
- **デメリット**: 処置効果が小さいと検出困難

**T-Learner**
```python
from econml.metalearners import TLearner

tlearner = TLearner(
    models=[RandomForestRegressor(), RandomForestRegressor()]
)
tlearner.fit(Y, T, X=X)
cate = tlearner.effect(X_test)
```
- **メリット**: 柔軟、各群で異なるモデル
- **デメリット**: サンプルサイズが小さいと分散大

**X-Learner**
```python
from econml.metalearners import XLearner

xlearner = XLearner(
    models=[RandomForestRegressor(), RandomForestRegressor()],
    propensity_model=LogisticRegression()
)
xlearner.fit(Y, T, X=X)
cate = xlearner.effect(X_test)
```
- **メリット**: 小サンプルに強い、不均衡データ対応
- **デメリット**: 傾向スコア推定が必要

**Domain Adaptation Learner**
```python
from econml.metalearners import DomainAdaptationLearner

dalearner = DomainAdaptationLearner(
    models=[RandomForestRegressor(), RandomForestRegressor()],
    final_models=[RandomForestRegressor(), RandomForestRegressor()]
)
```
- **メリット**: 群間分布差に対応
- **デメリット**: 複雑、計算コスト高

---

## 2. DML（Double Machine Learning / 二重機械学習）

### 概要
Chernozhukov et al. (2018)による手法。
ニュアンスパラメータ（傾向スコア、アウトカム）をMLで推定し、直交化により偽相関を除去。

### 利用可能な手法

| 手法 | クラス名 | 特徴 | 高次元対応 | 推奨度 |
|------|---------|------|-----------|--------|
| **DML (汎用)** | `DML` | 最も柔軟な基底クラス | ◎ | ★★★★☆ |
| **Linear DML** | `LinearDML` | 線形最終モデル | ◎ | ★★★★★ |
| **Sparse Linear DML** | `SparseLinearDML` | Lasso正則化 | ◎◎ | ★★★★★ |
| **Non-Parametric DML** | `NonParamDML` | 非パラメトリック | ○ | ★★★☆☆ |
| **Kernel DML** | `KernelDML` | カーネル法 | ○ | ★★★☆☆ |
| **Causal Forest DML** | `CausalForestDML` | 因果フォレスト | ◎ | ★★★★★ |

#### 詳細

**LinearDML**（最推奨）
```python
from econml.dml import LinearDML

dml = LinearDML(
    model_y=RandomForestRegressor(),
    model_t=RandomForestClassifier(),
    discrete_treatment=True,
    cv=5
)
dml.fit(Y, T, X=X, W=W)
cate = dml.effect(X_test)
```
- **メリット**: 理論保証、高次元対応、解釈可能
- **デメリット**: 最終モデルは線形のみ
- **用途**: 標準的なCATE推定

**SparseLinearDML**（高次元データ用）
```python
from econml.dml import SparseLinearDML

sparse_dml = SparseLinearDML(
    model_y=RandomForestRegressor(),
    model_t=RandomForestClassifier(),
    cv=5,
    alpha='auto'  # Lasso正則化パラメータ
)
```
- **メリット**: p >> n でも使用可能、特徴選択
- **デメリット**: 線形性の仮定
- **用途**: RNA-seqなど高次元遺伝子データ

**CausalForestDML**（非線形 + 高精度）
```python
from econml.dml import CausalForestDML

cf_dml = CausalForestDML(
    model_y=RandomForestRegressor(),
    model_t=RandomForestClassifier(),
    n_estimators=4000,
    min_samples_leaf=5,
    cv=5
)
```
- **メリット**: 非線形CATE、信頼区間
- **デメリット**: 計算コスト高
- **用途**: 複雑な非線形効果

**NonParamDML**（完全非パラメトリック）
```python
from econml.dml import NonParamDML

nparam_dml = NonParamDML(
    model_y=GradientBoostingRegressor(),
    model_t=GradientBoostingClassifier(),
    model_final=GradientBoostingRegressor(),
    cv=5
)
```
- **メリット**: 柔軟性最大
- **デメリット**: 過学習リスク、解釈困難
- **用途**: 探索的分析

**KernelDML**
```python
from econml.dml import KernelDML

kernel_dml = KernelDML(
    model_y=RandomForestRegressor(),
    model_t=RandomForestClassifier(),
    cv=5
)
```
- **メリット**: 滑らかなCATE推定
- **デメリット**: カーネル選択が難しい

---

## 3. DR（Doubly Robust / 二重頑健推定）

### 概要
傾向スコアとアウトカムモデルのどちらかが正しければ一致推定量となる。
DMLと類似だが、理論的に頑健性が高い。

### 利用可能な手法

| 手法 | クラス名 | 特徴 | 高次元対応 | 推奨度 |
|------|---------|------|-----------|--------|
| **DRLearner** | `DRLearner` | 汎用DR推定 | ◎ | ★★★★★ |
| **Linear DRLearner** | `LinearDRLearner` | 線形最終モデル | ◎ | ★★★★★ |
| **Sparse Linear DRLearner** | `SparseLinearDRLearner` | Lasso正則化 | ◎◎ | ★★★★★ |
| **Forest DRLearner** | `ForestDRLearner` | 因果フォレスト | ◎ | ★★★★☆ |

#### 詳細

**DRLearner**（現在の実装）
```python
from econml.dr import DRLearner
from sklearn.linear_model import LassoCV, RidgeCV, LogisticRegressionCV

dr = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=RidgeCV(cv=5),
    model_final=LassoCV(cv=5),
    cv=5
)
dr.fit(Y, T, X=X)
cate = dr.effect(X_test)
```
- **メリット**: 頑健性最高、線形DGPに最適
- **デメリット**: 最終モデルは回帰のみ
- **用途**: 標準的なCATE推定（demov04で使用）

**LinearDRLearner**
```python
from econml.dr import LinearDRLearner

linear_dr = LinearDRLearner(
    model_propensity=LogisticRegression(),
    model_regression='auto',
    cv=5
)
```
- **メリット**: シンプル、高速
- **用途**: 線形効果の場合

**SparseLinearDRLearner**（高次元用）
```python
from econml.dr import SparseLinearDRLearner

sparse_dr = SparseLinearDRLearner(
    model_propensity='auto',
    model_regression='auto',
    alpha='auto',
    cv=5
)
```
- **メリット**: p >> n 対応、自動特徴選択
- **用途**: 遺伝子データ、高次元

**ForestDRLearner**
```python
from econml.dr import ForestDRLearner

forest_dr = ForestDRLearner(
    model_propensity=LogisticRegression(),
    model_regression=RandomForestRegressor(),
    n_estimators=4000,
    min_samples_leaf=5
)
```
- **メリット**: 非線形CATE、信頼区間
- **用途**: 複雑な非線形効果

---

## 4. OrthoForest（直交フォレスト）

### 概要
局所的にニュアンスパラメータを推定し、適応的な推論を行う。
データ駆動的にXの空間を学習。

### 利用可能な手法

| 手法 | クラス名 | 基底手法 | 特徴 |
|------|---------|---------|------|
| **DML OrthoForest** | `DMLOrthoForest` | DML | 局所DML |
| **DR OrthoForest** | `DROrthoForest` | DR | 局所DR |

#### 詳細

**DMLOrthoForest**
```python
from econml.orf import DMLOrthoForest

orf_dml = DMLOrthoForest(
    n_trees=4000,
    min_leaf_size=10,
    max_depth=10,
    subsample_ratio=0.7,
    model_Y=RandomForestRegressor(),
    model_T=RandomForestClassifier()
)
orf_dml.fit(Y, T, X=X, W=W)
```
- **メリット**: 局所適応、異質性検出
- **デメリット**: 計算コスト高、複雑
- **用途**: 局所的に異なるCATEパターン

**DROrthoForest**
```python
from econml.orf import DROrthoForest

orf_dr = DROrthoForest(
    n_trees=4000,
    min_leaf_size=10,
    propensity_model=LogisticRegression(),
    model_Y=RandomForestRegressor()
)
```
- **メリット**: DR + 局所適応
- **用途**: 頑健性と局所性の両立

---

## 5. IV（Instrumental Variables / 操作変数法）

### 概要
未観測交絡がある場合に、操作変数を用いて因果効果を推定。
内生性（endogeneity）に対応。

### 利用可能な手法

EconMLのIVモジュールには複数の手法が含まれていますが、
主要なものは以下の通り：

- **DMLIV**: Double ML with IV
- **DRIV**: Doubly Robust with IV
- **IntentToTreatDRIV**: ITT推定（コンプライアンス不完全）
- **LinearDMLIV**: 線形モデル + IV

#### 使用例

```python
from econml.iv.dml import DMLIV

dmliv = DMLIV(
    model_Y_X=RandomForestRegressor(),
    model_T_X=RandomForestRegressor(),
    model_T_XZ=RandomForestClassifier(),
    model_final=LassoCV()
)
# Z: 操作変数
dmliv.fit(Y, T, Z=Z, X=X)
```

**用途**:
- 無作為化されていない観察研究
- コンプライアンスが不完全なRCT
- 内生性がある場合

---

## 6. Policy Learning（政策学習）

### 概要
最適な治療割り当てルール（Policy）を学習。
CATEではなく、直接最適化を行う。

### 利用可能な手法

| 手法 | クラス名 | 説明 | 出力 |
|------|---------|------|------|
| **Policy Tree** | `PolicyTree` | 決定木ベース | 解釈可能なルール |
| **Policy Forest** | `PolicyForest` | ランダムフォレスト | 高精度 |
| **DR Policy Tree** | `DRPolicyTree` | DR + Tree | 頑健 + 解釈可能 |
| **DR Policy Forest** | `DRPolicyForest` | DR + Forest | 頑健 + 高精度 |

#### 使用例

```python
from econml.policy import DRPolicyTree

policy_tree = DRPolicyTree(
    max_depth=2,
    min_samples_leaf=10,
    model_regression=RandomForestRegressor(),
    model_propensity=LogisticRegression()
)
policy_tree.fit(Y, T, X=X)

# 最適な治療割り当て
treatment_recommendation = policy_tree.predict(X_test)
```

**用途**:
- 個別化医療（誰に治療すべきか）
- マーケティング（誰にオファーすべきか）
- 意思決定ルールの構築

---

## 7. その他のツール

### Validation / Scoring

| ツール | 説明 |
|-------|------|
| **RScorer** | CATE推定のR²スコア計算 |
| **EnsembleCateEstimator** | 複数手法のアンサンブル |

### Inference（推論）

- **Bootstrap**: ブートストラップ信頼区間
- **Normal Inference**: 正規近似信頼区間

### Utilities

- **SHAP統合**: SHAP値でCATEの要因分解
- **DoWhy統合**: 因果グラフとの統合

---

## 手法選択ガイド

### データ特性による選択

| データ特性 | 推奨手法 | 理由 |
|-----------|---------|------|
| **低次元 (p < n)** | TLearner, DRLearner | シンプルで十分 |
| **高次元 (p >> n)** | SparseLinearDML/DR, ElasticNet | L1正則化必須 |
| **非線形効果** | CausalForestDML, ForestDRLearner | ツリーベース |
| **小サンプル** | XLearner, DRLearner | 頑健性重視 |
| **大サンプル** | LinearDML, NonParamDML | 柔軟なモデル |
| **内生性あり** | DMLIV, DRIV | 操作変数使用 |

### 目的による選択

| 目的 | 推奨手法 |
|------|---------|
| **CATE推定（標準）** | DRLearner, LinearDML |
| **CATE推定（高次元）** | SparseLinearDRLearner, SparseLinearDML |
| **CATE推定（非線形）** | CausalForestDML, NonParamDML |
| **政策決定** | DRPolicyTree, DRPolicyForest |
| **頑健性重視** | DRLearner, XLearner |
| **解釈性重視** | LinearDML, PolicyTree |

### RNA-seq解析での推奨

| シナリオ | 推奨手法 | 設定 |
|---------|---------|------|
| **標準（n=500, p=20000）** | SparseLinearDRLearner | alpha='auto' |
| **探索的** | DRLearner + ElasticNetCV | l1_ratio=[0.5, 0.9] |
| **非線形効果疑い** | CausalForestDML | n_estimators=4000 |
| **政策決定** | DRPolicyTree | max_depth=3 |

---

## まとめ

### 総数
- **20種類以上**の因果推論手法が利用可能
- Meta-Learners（4）+ DML（6）+ DR（4）+ OrthoForest（2）+ IV（複数）+ Policy（4）

### 最も実用的な手法 TOP 5

| 順位 | 手法 | 用途 | 理由 |
|-----|------|------|------|
| **1** | **DRLearner** | 標準的CATE推定 | 頑健、高次元対応、理論保証 |
| **2** | **LinearDML** | 線形CATE推定 | DMLの理論、解釈可能 |
| **3** | **SparseLinearDRLearner** | 高次元CATE推定 | p >> n対応、自動特徴選択 |
| **4** | **CausalForestDML** | 非線形CATE推定 | 高精度、信頼区間 |
| **5** | **XLearner** | 小サンプルCATE推定 | シンプル、不均衡データ対応 |

### DemoV04での使用

現在の実装: **DRLearner**
- 線形DGPに最適
- n=1000でR²=0.992達成
- 高次元対応可能

### DemoV05への推奨

**高次元RNA-seqを想定する場合**:
```python
from econml.dr import SparseLinearDRLearner

sparse_dr = SparseLinearDRLearner(
    model_propensity='auto',
    model_regression='auto',
    alpha='auto',
    cv=5
)
```

または

```python
from econml.dr import DRLearner
from sklearn.linear_model import ElasticNetCV

dr = DRLearner(
    model_regression=ElasticNetCV(l1_ratio=[0.5, 0.9]),
    model_final=ElasticNetCV(l1_ratio=[0.5, 0.9]),
    cv=5
)
```

---

## 参考文献

- [EconML Official Documentation](https://www.pywhy.org/EconML/)
- [Meta-Learners Documentation](https://www.pywhy.org/EconML/spec/estimation/metalearners.html)
- [Forest-Based Methods Documentation](https://www.pywhy.org/EconML/spec/estimation/forest.html)
- [GitHub Repository](https://github.com/py-why/EconML)
- EconML version: **0.16.0**

---

**保存日時**: 2026-02-20
**保存場所**: `/home/rstudio/work/cfomics/sandbox/econml_all_methods.md`
