# SparseLinearDRLearner の詳細解説

## DRLearner ファミリーの比較

| クラス | model_propensity | model_regression | model_final | 自由度 | 用途 |
|--------|-----------------|------------------|-------------|-------|------|
| **DRLearner** | 自由に選択可能 | 自由に選択可能 | **自由に選択可能** | ★★★ | 最も柔軟 |
| **LinearDRLearner** | 自由に選択可能 | 自由に選択可能 | **線形回帰（固定）** | ★★☆ | 速度重視 |
| **SparseLinearDRLearner** | 自由に選択可能 | 自由に選択可能 | **Lasso（固定）** | ★★☆ | 解釈性重視 |
| **ForestDRLearner** | 自由に選択可能 | 自由に選択可能 | **RandomForest（固定）** | ★★☆ | 非線形 |

## SparseLinearDRLearner の特徴

### 内部実装（固定されている部分）

```python
# SparseLinearDRLearner の内部（ユーザーは model_final を指定できない）
class SparseLinearDRLearner:
    def __init__(self,
                 model_propensity=...,  # ユーザーが選択
                 model_regression=...,  # ユーザーが選択
                 # model_final は指定できない（内部で自動設定）
                 ):
        # 内部で自動的に Lasso with CV を設定
        self.model_final = MultiTaskLassoCV(
            cv=...,
            alphas=None,  # 自動で範囲を決定
            max_iter=10000
        )
```

### 使用例

```python
from econml.dr import SparseLinearDRLearner
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV

# ユーザーが指定できるのは propensity と regression のみ
est = SparseLinearDRLearner(
    model_propensity=LogisticRegressionCV(cv=5),  # 自由
    model_regression=ElasticNetCV(cv=5),          # 自由
    # model_final は自動で Lasso（指定不可）
    cv=5,
    random_state=42
)

est.fit(Y, T, X=X)

# 係数を取得（スパース）
coef = est.coef_  # shape: (n_features,)
selected_features = np.where(coef != 0)[0]
print(f"選択特徴数: {len(selected_features)} / {X.shape[1]}")
# 出力例: 選択特徴数: 23 / 20000
```

## なぜ Lasso 固定なのか？

### 目的: **最大のスパース性 + 統計的推論**

1. **スパース性の保証**
   - Lasso (L1正則化) は多くの係数をゼロにする
   - ElasticNet や Ridge ではゼロにならない係数も残る
   - 論文・レポートで「この遺伝子のみが重要」と明言したい場合に最適

2. **統計的推論（信頼区間）**
   - `SparseLinearDRLearner` は係数の信頼区間も計算可能
   - どの遺伝子が統計的に有意か検定できる

   ```python
   # 信頼区間の取得
   est.fit(Y, T, X=X, inference='auto')
   summary = est.summary()
   print(summary)
   # 出力: 各遺伝子の係数、標準誤差、p値、信頼区間
   ```

3. **解釈性の最大化**
   - 係数が明確（β1 * gene1 + β2 * gene2 + ...）
   - ゼロ係数は「効果なし」と解釈可能
   - 非ゼロ係数のみ報告すれば良い

## DRLearner × Lasso との違い

### DRLearner × Lasso（手動設定）

```python
from econml.dr import DRLearner
from sklearn.linear_model import LassoCV

est = DRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=ElasticNetCV(cv=5),
    model_final=LassoCV(cv=5),  # 手動で Lasso 指定
    cv=5
)
```

**特徴**:
- ✓ 柔軟性高い（ElasticNet に変更も可能）
- × 統計的推論は手動で実装が必要
- × インターフェースが若干複雑

### SparseLinearDRLearner

```python
from econml.dr import SparseLinearDRLearner

est = SparseLinearDRLearner(
    model_propensity=LogisticRegressionCV(cv=5),
    model_regression=ElasticNetCV(cv=5),
    cv=5
)
```

**特徴**:
- ✓ **統計的推論が組み込み済み**（.summary(), .coef_stderr()）
- ✓ インターフェースがシンプル
- ✓ Lasso に特化した最適化
- × model_final は変更不可

## 実践的な使い分け

| 状況 | 推奨 | 理由 |
|------|------|------|
| **論文で係数を報告** | **SparseLinearDRLearner** | 信頼区間、p値が必要 |
| **重要遺伝子の同定** | **SparseLinearDRLearner** | 最もスパース |
| **柔軟性が欲しい** | **DRLearner × LassoCV** | ElasticNet に変更可能 |
| **予測精度重視** | **DRLearner × ElasticNetCV** | スパース性より安定性 |
| **非線形効果** | **DRLearner × RandomForest** | Lasso では捉えられない |

## 内部アルゴリズム

### Lasso の正則化パラメータ選択

SparseLinearDRLearner は内部で以下を実行：

```python
# 疑似コード
alphas = np.logspace(-4, 1, 100)  # 0.0001 から 10 まで
cv_scores = []

for alpha in alphas:
    model = Lasso(alpha=alpha)
    score = cross_val_score(model, X_pseudo, tau_pseudo, cv=5)
    cv_scores.append(score.mean())

best_alpha = alphas[np.argmax(cv_scores)]
final_model = Lasso(alpha=best_alpha)
```

### DRLearner × LassoCV も同様だが...

- **DRLearner × LassoCV**: ユーザーが alpha 範囲を指定可能
- **SparseLinearDRLearner**: 範囲は自動決定（最適化済み）

## まとめ

### SparseLinearDRLearner の特徴

1. **model_final = Lasso（CV付き）で固定**
2. **統計的推論が組み込み**（信頼区間、p値）
3. **最大のスパース性**（最も多くの係数がゼロ）
4. **解釈性特化**（論文・レポート向け）

### いつ使うか？

**以下のいずれかに該当する場合**:
- ✓ 論文で「遺伝子Aの効果は β=0.5 (95%CI: 0.3-0.7, p<0.01)」と報告したい
- ✓ 重要遺伝子を10-50個に厳密に絞りたい
- ✓ モデルの解釈性が最優先
- ✓ 線形効果を仮定できる

**使わない場合**:
- × 非線形効果が強い → **DRLearner × RandomForest**
- × 予測精度最優先 → **DRLearner × ElasticNet**
- × 柔軟に試したい → **DRLearner（自由度最大）**
