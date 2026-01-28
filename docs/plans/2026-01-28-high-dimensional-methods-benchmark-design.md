# 高次元因果推論手法とベンチマーク設計

**作成日**: 2026-01-28
**ステータス**: ドラフト
**対象パッケージ**: cfomics (R-native手法), cfomicsPython (Python手法)

---

## 1. 概要

### 1.1 目的

cfomicsエコシステムに高次元オミクスデータ向けの因果推論手法を追加し、包括的なベンチマークフレームワークを構築する。

### 1.2 設計目標

1. **統一API** - 全手法が `cf_fit()` から利用可能
2. **包括的ベンチマーク** - 手法間の公平な比較
3. **実用ガイダンス** - データ特性に基づく手法選択支援
4. **R-native優先** - Python依存を最小化

### 1.3 パッケージ分離方針

| パッケージ | 内容 | Python依存 |
|-----------|------|-----------|
| cfomics | R-native手法（既存 + 新規高次元手法） | なし |
| cfomicsPython | Python依存手法（DoWhy, GANITE, CAVAE等） | 必須 |

分離理由：
- ローカル環境（Apple Silicon）ではGPU不十分
- HPC環境でのみPython/GPUベース手法を実行
- 依存関係の明確化

---

## 2. 手法設計

### 2.1 手法一覧

#### 既存手法（低〜中次元向け）

| メソッド名 | 説明 | 依存パッケージ |
|-----------|------|---------------|
| `grf` | Generalized Random Forest | grf |
| `ipw` | Inverse Probability Weighting | ipw, survey |
| `gformula` | G-computation | base R (lm) |

#### 新規手法（高次元向け）

| メソッド名 | 説明 | 依存パッケージ | 主な用途 |
|-----------|------|---------------|---------|
| `hdml` | High-Dimensional Machine Learning (Debiased Lasso) | hdm または DoubleML | 高次元ATE推定 |
| `bcf` | Bayesian Causal Forests | bcf, dbarts | 異質性効果 + 不確実性定量化 |
| `tmle` | Targeted Maximum Likelihood Estimation | tmle3, sl3 | 二重ロバスト推定 |
| `hdps` | High-Dimensional Propensity Score | glmnet | 正則化傾向スコア |

### 2.2 統一インターフェース

```r
# 全手法で同一構文
result <- cf_fit(Y ~ T | X1 + X2 + ... + Xp,
                 data = omics_data,
                 method = "hdml")

# 結果も同一構造
result$fit$res$ate      # 平均処置効果
result$fit$res$ite      # 個別処置効果（対応手法のみ）
result$fit$res$y0_hat   # 反事実 Y(0)
result$fit$res$y1_hat   # 反事実 Y(1)
```

### 2.3 手法固有パラメータ

```r
# hdml: 正則化の種類
cf_fit(..., method = "hdml",
       penalty = "lasso",           # "lasso", "ridge", "elastic_net"
       lambda = "cv")               # "cv", "1se", または数値

# bcf: MCMC設定
cf_fit(..., method = "bcf",
       n_burn = 1000,
       n_iter = 2000,
       n_chains = 4)

# tmle: Super Learnerのライブラリ
cf_fit(..., method = "tmle",
       sl_lib = c("SL.glmnet", "SL.ranger", "SL.xgboost"))

# hdps: 傾向スコア推定設定
cf_fit(..., method = "hdps",
       lambda = "cv",
       trim = c(0.01, 0.99))        # PSのトリミング範囲
```

### 2.4 ITE出力対応

| 手法 | ITE出力 | 備考 |
|------|---------|------|
| GRF | ✓ | 個別予測を直接返す |
| BCF | ✓ | 個別予測 + 事後分布 |
| G-formula | ✓ | Y(1) - Y(0) を計算可能 |
| HDML | △ | 実装による（ATE中心の場合が多い） |
| TMLE | △ | 主にATE、拡張でCATE可能 |
| IPW | × | ATEのみ |
| HDPS | × | ATEのみ |

ITE非対応手法では `fit$res$ite <- NA` を返す。

### 2.5 ファイル構成

```
packages/cfomics/R/
├── methods_hdml.R      # Debiased Lasso / Double ML
├── methods_bcf.R       # Bayesian Causal Forests
├── methods_tmle.R      # TMLE + Super Learner
├── methods_hdps.R      # HD-Propensity Score
├── benchmark_dgp.R     # DGP定義
├── benchmark_run.R     # ベンチマーク実行
└── benchmark_viz.R     # ベンチマーク可視化
```

---

## 3. ベンチマーク設計

### 3.1 シナリオ一覧

#### 基本シナリオ

| ID | 名称 | n | p | 交絡構造 | 目的 | 実応用例 |
|----|------|---|---|---------|------|---------|
| S1 | 基準 | 500 | 50 | 線形・疎・独立 | ベースライン | 臨床試験補助解析、少数バイオマーカー研究 |
#### S5: 非線形交絡（サブシナリオ）

| ID | 非線形タイプ | 説明 | 目的 | 実応用例 |
|----|-------------|------|------|---------|
| S5a | 二乗項のみ | Y, T ∝ X² | 単純な非線形 | 用量反応の二次効果 |
| S5b | 三角関数 | Y, T ∝ sin(πX) | 周期的パターン | 周期的生物リズム |
| S5c | 交互作用 | Y, T ∝ X₁·X₂ | 変数間相互作用 | パスウェイ相互作用 |
| S5d | 複合非線形 | 上記の組み合わせ | 複雑なパターン | 現実的な生物学的関係 |
| S5e | 閾値効果 | Y, T ∝ I(X > c) | ステップ関数的 | バイオマーカー閾値 |

非線形強度スイープ: `nonlinearity_strength = c(0, 0.25, 0.5, 0.75, 1.0)`
（0 = 完全線形, 1.0 = 完全非線形）

#### S6: 密交絡（スイープ実験）

| パラメータ | スイープ値 |
|-----------|-----------|
| 交絡変数数 | 10, 25, 50, 100, 200, 500 |
| 係数スケーリング | fixed, sqrt, linear |

係数スケーリング:
- `fixed`: 係数固定（0.1）→ 交絡総効果が増加
- `sqrt`: 係数 ∝ 1/√n_conf → 分散一定
- `linear`: 係数 ∝ 1/n_conf → 効果一定

実応用例: 環境疫学、バッチ効果、多数の生活習慣因子

#### S2-S3: 次元数スイープ実験

個別シナリオではなく、次元数を段階的に変化させるスイープ実験として設計。

| パラメータ | 値 |
|-----------|-----|
| n（固定） | 200 |
| p（変動） | 50, 100, 200, 500, 1000, 2000, 5000, 10000 |
| n/p比 | 4.0, 2.0, 1.0, 0.4, 0.2, 0.1, 0.04, 0.02 |
| 交絡構造 | 線形・疎・独立（固定） |
| 目的 | 各手法の次元限界点を特定 |
| 実応用例 | RNA-seq、プロテオミクス、マルチオミクス統合 |

#### S4: 異質な処置効果（サブシナリオ）

| ID | 異質性タイプ | 強さ | サンプル偏り | 目的 | 実応用例 |
|----|-------------|------|-------------|------|---------|
| S4a | 線形 | 中 | なし | 基本的なHTE | 薬剤反応性の遺伝的背景 |
| S4b | 非線形 | 中 | なし | 非線形パターンの検出 | 複雑な効果修飾 |
| S4c | サブグループ | 中 | なし | 離散的異質性 | がんサブタイプ別効果 |
| S4d | 質的交互作用 | 中 | なし | 有害サブグループ検出 | 一部患者で有害な治療 |
| S4e | 線形 | 弱 | なし | 検出力の評価 | 微弱な効果修飾 |
| S4f | 線形 | 強 | なし | 強い異質性での精度 | 明確なレスポンダー/非レスポンダー |
| S4g | サブグループ | 中 | 中程度 | 偏りへのロバスト性 | 一般的なコホート構成 |
| S4h | サブグループ | 中 | 極端 | 希少サブグループ推定 | 希少疾患、希少遺伝子型 |

#### S7: 弱overlap（サブシナリオ）

| ID | タイプ | PS範囲 | 目的 | 実応用例 |
|----|--------|--------|------|---------|
| S7a-good | 係数スケール | [0.1, 0.9] | ベースライン（良好なoverlap） | 理想的観察研究 |
| S7a-moderate | 係数スケール | [0.05, 0.95] | 軽度違反 | 一般的な観察研究 |
| S7a-weak | 係数スケール | [0.01, 0.99] | 実質的違反 | 強い適応基準 |
| S7a-extreme | 係数スケール | [0.001, 0.999] | 重度違反 | ほぼ決定的な処置割り当て |
| S7b | 構造的違反 | 一部領域でPS≈0 | 完全な違反への対処 | 禁忌条件、希少疾患 |
| S7c | 非対称 | 片側のみ極端 | 群間不均衡の影響 | 処置群過多/過少 |

#### S8: 共変量シフト（サブシナリオ）

| ID | シフトタイプ | パラメータ | 目的 | 実応用例 |
|----|-------------|-----------|------|---------|
| S8a | 平均シフトのみ | mean_shift=0.5,1.0,2.0 | 位置のずれ | 異なる母集団 |
| S8b | 分散シフトのみ | var_ratio=1.0,1.5,2.0,3.0 | スケールのずれ | 測定精度の違い |
| S8c | 相関シフトのみ | cor_shift=0.1,0.2,0.3,0.4 | 構造のずれ | 技術的変動 |
| S8d | 複合シフト | 上記の組み合わせ | 現実的状況 | コホート間比較 |

シフト強度スイープ: `shift_magnitude = c("mild", "moderate", "severe")`

実応用例: 異なるコホート間比較、バッチ効果、前向きvs後ろ向き研究

#### S9: 相関交絡（サブシナリオ）

| ID | 相関構造 | パラメータ | 目的 | 実応用例 |
|----|---------|-----------|------|---------|
| S9a | ブロック相関（弱） | within_cor=0.3 | 弱い共発現 | 遠縁遺伝子 |
| S9b | ブロック相関（中） | within_cor=0.5 | 中程度共発現 | 同一パスウェイ |
| S9c | ブロック相関（強） | within_cor=0.7 | 強い共発現 | 遺伝子クラスター |
| S9d | AR(1)構造 | rho=0.3,0.5,0.7 | 順序依存 | ゲノム位置依存 |
| S9e | 因子モデル | n_factors=5,10,20 | 潜在因子構造 | 転写因子調節 |

相関強度スイープ: `within_cor = c(0.1, 0.3, 0.5, 0.7, 0.9)`

実応用例: 遺伝子共発現、パスウェイ内遺伝子群、転写因子ネットワーク

#### S10: 未観測交絡（サブシナリオ）

| ID | 設定 | パラメータ | 目的 | 実応用例 |
|----|------|-----------|------|---------|
| S10a | 弱い未観測交絡 | strength=0.2 | 軽微なバイアス | 軽度のバッチ効果 |
| S10b | 中程度 | strength=0.5 | 中程度バイアス | 一般的な未測定交絡 |
| S10c | 強い未観測交絡 | strength=0.8 | 深刻なバイアス | 重大な欠落変数 |
| S10d | 多数の未観測 | n_unobserved=1,5,10,20 | 未観測変数の数 | 複合的バイアス |
| S10e | 観測割合 | obs_ratio=0.9,0.7,0.5 | 交絡の観測割合 | 部分的測定 |

未観測強度スイープ: `strength = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0)`

追加評価指標:
- 実バイアス: 推定ATE - 真ATE
- 理論バイアス比: 実バイアス / 理論予測バイアス
- 感度パラメータ: どの程度の未観測交絡で推定が逆転するか

実応用例: バッチ効果、測定されない生活習慣、技術的変動

#### S11: 合流点（サブシナリオ）

| ID | 設定 | パラメータ | 目的 | 実応用例 |
|----|------|-----------|------|---------|
| S11a | 弱いcollider | effect=0.2 | 軽微な選択バイアス | 軽度の選択 |
| S11b | 中程度 | effect=0.5 | 中程度バイアス | 一般的な選択バイアス |
| S11c | 強いcollider | effect=0.8 | 深刻なバイアス | 強い選択圧 |
| S11d | collider数 | n_colliders=1,3,5,10 | 複数のcollider | 多段階選択 |
| S11e | 混合構造 | 交絡+collider混在 | 変数選択の難しさ | 現実的状況 |

collider強度スイープ: `effect = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0)`

追加評価指標:
- Collider調整バイアス: colliderを含めた場合 vs 除外した場合
- 変数選択正確度: 交絡のみを正しく選択できたか
- Oracle比較: 真の交絡のみで調整した場合との差

実応用例: 選択バイアス、生存者バイアス、入院バイアス

### 3.2 DGP実装

#### S1: 基準

```r
dgp_baseline <- function(n = 500, p = 50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 疎な交絡: 最初の10変数のみが関与
  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  beta_y <- c(rep(0.2, 10), rep(0, p - 10))

  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    dgp_name = "baseline",
    dgp_params = list(n = n, p = p)
  )
}
```

#### S2-S3: 次元数スイープ

```r
dgp_dimension_sweep <- function(n = 200, p = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 疎な交絡: 最初の10変数のみ（pに関わらず固定）
  n_conf <- min(10, p)
  beta_t <- c(rep(0.5, n_conf), rep(0, p - n_conf))
  beta_y <- c(rep(0.3, n_conf), rep(0, p - n_conf))

  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    dgp_name = "dimension_sweep",
    dgp_params = list(n = n, p = p, n_p_ratio = n / p)
  )
}
```

#### S4a: 線形異質性

```r
dgp_heterogeneous_linear <- function(n = 500, p = 500,
                                      strength = 1.0,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  # 異質な処置効果: X1とX2が効果修飾因子
  tau_base <- 2.0
  tau <- tau_base + strength * X[,1] + 0.3 * strength * X[,2]

  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = mean(tau),
    true_ite = tau,
    ite_sd = sd(tau),
    propensity_score = ps,
    dgp_name = "heterogeneous_linear",
    dgp_params = list(n = n, p = p, strength = strength)
  )
}
```

#### S4b: 非線形異質性

```r
dgp_heterogeneous_nonlinear <- function(n = 500, p = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  # 非線形な処置効果
  tau <- 2.0 + 1.0 * sin(pi * X[,1]) + 0.5 * X[,2]^2

  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = mean(tau),
    true_ite = tau,
    propensity_score = ps,
    dgp_name = "heterogeneous_nonlinear",
    dgp_params = list(n = n, p = p)
  )
}
```

#### S4c: サブグループ異質性

```r
dgp_heterogeneous_subgroup <- function(n = 500, p = 500,
                                        subgroup_props = c(1/3, 1/3, 1/3),
                                        subgroup_effects = c(0.5, 2.0, 4.0),
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  # サブグループ割り当て（X1に基づく、または確率的）
  if (length(subgroup_props) == 3) {
    breaks <- qnorm(cumsum(subgroup_props)[1:2])
    subgroup <- cut(X[,1], breaks = c(-Inf, breaks, Inf), labels = 1:3)
    subgroup <- as.integer(subgroup)
  }

  tau <- subgroup_effects[subgroup]

  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = sum(subgroup_props * subgroup_effects),
    true_ite = tau,
    subgroup = subgroup,
    subgroup_n = table(subgroup),
    propensity_score = ps,
    dgp_name = "heterogeneous_subgroup",
    dgp_params = list(n = n, p = p,
                      subgroup_props = subgroup_props,
                      subgroup_effects = subgroup_effects)
  )
}
```

#### S4d: 質的交互作用

```r
dgp_heterogeneous_qualitative <- function(n = 500, p = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.3, 10), rep(0, p - 10))
  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  # 一部のサブグループでは処置が有害
  tau <- 2.0 * X[,1]  # X1 < 0 なら負の効果

  beta_y <- c(rep(0.2, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = mean(tau),
    true_ite = tau,
    prop_harmed = mean(tau < 0),
    prop_benefited = mean(tau > 0),
    propensity_score = ps,
    dgp_name = "heterogeneous_qualitative",
    dgp_params = list(n = n, p = p)
  )
}
```

#### S5: 非線形交絡（拡張版）

```r
# S5: 非線形交絡の統一DGP（タイプと強度を指定可能）
dgp_nonlinear_confounding <- function(n = 500, p = 500,
                                       nonlinear_type = "combined",
                                       strength = 1.0,
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 線形成分（常に含む）
  linear_t <- 0.3 * X[,1] + 0.3 * X[,2]
  linear_y <- 0.2 * X[,1] + 0.2 * X[,2]

  # 非線形成分（タイプにより異なる）
  nonlinear_t <- nonlinear_y <- 0

  if (nonlinear_type == "quadratic" || nonlinear_type == "combined") {
    # S5a: 二乗項
    nonlinear_t <- nonlinear_t + 0.5 * X[,1]^2
    nonlinear_y <- nonlinear_y + 0.3 * X[,1]^2
  }

  if (nonlinear_type == "trigonometric" || nonlinear_type == "combined") {
    # S5b: 三角関数
    nonlinear_t <- nonlinear_t + 0.3 * sin(pi * X[,3])
    nonlinear_y <- nonlinear_y + 0.2 * sin(pi * X[,2])
  }

  if (nonlinear_type == "interaction" || nonlinear_type == "combined") {
    # S5c: 交互作用
    nonlinear_t <- nonlinear_t + 0.2 * X[,1] * X[,2]
    nonlinear_y <- nonlinear_y + 0.1 * X[,1] * X[,3]
  }

  if (nonlinear_type == "threshold") {
    # S5e: 閾値効果
    nonlinear_t <- nonlinear_t + 0.5 * as.numeric(X[,1] > 0)
    nonlinear_y <- nonlinear_y + 0.3 * as.numeric(X[,2] > 0.5)
  }

  # 強度で非線形成分をスケーリング
  ps <- plogis(linear_t + strength * nonlinear_t)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- linear_y + strength * nonlinear_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    dgp_name = "nonlinear_confounding",
    dgp_params = list(n = n, p = p,
                      nonlinear_type = nonlinear_type,
                      strength = strength)
  )
}
```

#### S6: 密交絡（拡張版）

```r
dgp_dense_confounding <- function(n = 500, p = 500,
                                   n_confounders = 100,
                                   coef_scaling = "fixed",
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 係数スケーリング方式
  base_coef <- 0.1
  coef_size <- switch(coef_scaling,
    "fixed"  = base_coef,                        # 固定: 総効果増加
    "sqrt"   = base_coef / sqrt(n_confounders),  # √n: 分散一定
    "linear" = base_coef * 10 / n_confounders    # n: 効果一定
  )

  beta_t <- c(rep(coef_size, n_confounders), rep(0, p - n_confounders))
  beta_y <- c(rep(coef_size, n_confounders), rep(0, p - n_confounders))

  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + rnorm(n)

  # 診断情報
  total_confounding_t <- sum(abs(beta_t))
  total_confounding_y <- sum(abs(beta_y))

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    confounding_info = list(
      n_confounders = n_confounders,
      coef_size = coef_size,
      coef_scaling = coef_scaling,
      total_confounding_t = total_confounding_t,
      total_confounding_y = total_confounding_y
    ),
    dgp_name = "dense_confounding",
    dgp_params = list(n = n, p = p, n_confounders = n_confounders,
                      coef_scaling = coef_scaling)
  )
}
```

#### S7a: 弱overlap（係数スケール）

```r
dgp_weak_overlap <- function(n = 500, p = 500,
                              overlap_strength = "weak",
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # overlap強度に応じた係数
  coef_scale <- switch(overlap_strength,
    "good"     = 0.3,
    "moderate" = 0.6,
    "weak"     = 1.2,
    "extreme"  = 2.0
  )

  beta_t <- c(rep(coef_scale, 10), rep(0, p - 10))
  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  ps_summary <- list(
    min = min(ps),
    max = max(ps),
    mean = mean(ps),
    prop_extreme = mean(ps < 0.05 | ps > 0.95)
  )

  tau <- 2.0
  beta_y <- c(rep(0.3, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    ps_summary = ps_summary,
    overlap_strength = overlap_strength,
    dgp_name = "weak_overlap",
    dgp_params = list(n = n, p = p, overlap_strength = overlap_strength)
  )
}
```

#### S7b: 構造的positivity違反

```r
dgp_structural_violation <- function(n = 500, p = 500,
                                      violation_region = 0.2,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.5, 10), rep(0, p - 10))
  ps <- plogis(X %*% beta_t)

  # X1が極端に低い領域では処置を受けない
  threshold <- qnorm(violation_region)
  ps[X[,1] < threshold] <- 0.01

  T <- rbinom(n, 1, ps)

  violation_info <- list(
    n_violation_region = sum(X[,1] < threshold),
    n_treated_in_violation = sum(T[X[,1] < threshold] == 1),
    threshold = threshold
  )

  tau <- 2.0
  beta_y <- c(rep(0.3, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    violation_info = violation_info,
    dgp_name = "structural_violation",
    dgp_params = list(n = n, p = p, violation_region = violation_region)
  )
}
```

#### S7c: 非対称overlap

```r
dgp_asymmetric_overlap <- function(n = 500, p = 500,
                                    shift = 1.0,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  beta_t <- c(rep(0.8, 10), rep(0, p - 10))
  linear_pred <- X %*% beta_t

  # 正の領域を強調（処置群に偏る）
  ps <- plogis(linear_pred + shift)
  T <- rbinom(n, 1, ps)

  treatment_rate <- mean(T)

  tau <- 2.0
  beta_y <- c(rep(0.3, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    treatment_rate = treatment_rate,
    dgp_name = "asymmetric_overlap",
    dgp_params = list(n = n, p = p, shift = shift)
  )
}
```

#### S8: 共変量シフト

```r
dgp_covariate_shift <- function(n = 500, p = 500,
                                 mean_shift = 1.0,
                                 var_ratio = 1.5,
                                 cor_shift = 0.2,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_t <- n / 2
  n_c <- n / 2
  conf_idx <- 1:10

  # 対照群の相関構造
  Sigma_c <- diag(p)
  Sigma_c[conf_idx, conf_idx] <- 0.3
  diag(Sigma_c) <- 1

  # 処置群の相関構造（より強い相関、異なる分散）
  Sigma_t <- diag(p)
  Sigma_t[conf_idx, conf_idx] <- 0.3 + cor_shift
  diag(Sigma_t) <- var_ratio

  # 各群からサンプリング
  mu_c <- rep(0, p)
  mu_t <- c(rep(mean_shift, 10), rep(0, p - 10))

  X_control <- MASS::mvrnorm(n_c, mu = mu_c, Sigma = Sigma_c)
  X_treated <- MASS::mvrnorm(n_t, mu = mu_t, Sigma = Sigma_t)

  X <- rbind(X_control, X_treated)
  colnames(X) <- paste0("X", 1:p)
  T <- c(rep(0, n_c), rep(1, n_t))

  # シャッフル
  idx <- sample(n)
  X <- X[idx, ]
  T <- T[idx]

  tau <- 2.0
  beta_y <- c(rep(0.3, 10), rep(0, p - 10))
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    shift_params = list(
      mean_shift = mean_shift,
      var_ratio = var_ratio,
      cor_shift = cor_shift
    ),
    dgp_name = "covariate_shift",
    dgp_params = list(n = n, p = p, mean_shift = mean_shift,
                      var_ratio = var_ratio, cor_shift = cor_shift)
  )
}
```

#### S9: 相関交絡（拡張版）

```r
# S9a-c: ブロック相関構造
dgp_correlated_confounding_block <- function(n = 500, p = 500,
                                              n_blocks = 10,
                                              within_cor = 0.7,
                                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # ブロック相関構造（遺伝子共発現を模倣）
  block_size <- p / n_blocks
  Sigma <- diag(p)

  for (b in 1:n_blocks) {
    idx <- ((b-1)*block_size + 1):(b*block_size)
    Sigma[idx, idx] <- within_cor
    diag(Sigma)[idx] <- 1
  }

  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)

  # 各ブロックの代表変数が交絡
  beta_t <- beta_y <- rep(0, p)
  representative_idx <- seq(1, p, by = block_size)
  beta_t[representative_idx] <- 0.5
  beta_y[representative_idx] <- 0.3

  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    correlation_structure = list(
      type = "block",
      n_blocks = n_blocks,
      block_size = block_size,
      within_cor = within_cor
    ),
    dgp_name = "correlated_confounding_block",
    dgp_params = list(n = n, p = p, n_blocks = n_blocks, within_cor = within_cor)
  )
}

# S9d: AR(1)構造
dgp_correlated_confounding_ar1 <- function(n = 500, p = 500,
                                            rho = 0.5,
                                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # AR(1)相関構造: Σ[i,j] = rho^|i-j|
  Sigma <- rho^abs(outer(1:p, 1:p, "-"))

  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)

  # 等間隔で交絡変数を選択
  conf_idx <- seq(1, p, length.out = 10)
  beta_t <- beta_y <- rep(0, p)
  beta_t[conf_idx] <- 0.5
  beta_y[conf_idx] <- 0.3

  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    correlation_structure = list(
      type = "ar1",
      rho = rho
    ),
    dgp_name = "correlated_confounding_ar1",
    dgp_params = list(n = n, p = p, rho = rho)
  )
}

# S9e: 因子モデル
dgp_correlated_confounding_factor <- function(n = 500, p = 500,
                                               n_factors = 10,
                                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # X = L·F + ε, where F is n_factors latent factors
  L <- matrix(rnorm(p * n_factors, sd = 0.5), p, n_factors)
  F <- matrix(rnorm(n * n_factors), n, n_factors)
  noise <- matrix(rnorm(n * p, sd = 0.5), n, p)
  X <- F %*% t(L) + noise
  colnames(X) <- paste0("X", 1:p)

  # 潜在因子が交絡（観測変数を通じて）
  # 最初の5因子がTとYに影響
  gamma_t <- c(rep(0.5, 5), rep(0, n_factors - 5))
  gamma_y <- c(rep(0.3, 5), rep(0, n_factors - 5))

  ps <- plogis(F %*% gamma_t)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- F %*% gamma_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    correlation_structure = list(
      type = "factor",
      n_factors = n_factors,
      n_confounding_factors = 5
    ),
    dgp_name = "correlated_confounding_factor",
    dgp_params = list(n = n, p = p, n_factors = n_factors)
  )
}
```

#### S10: 未観測交絡（拡張版）

```r
dgp_unobserved_confounding <- function(n = 500, p = 500,
                                        n_unobserved = 5,
                                        unobserved_strength = 0.5,
                                        obs_ratio = NULL,
                                        seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # obs_ratioが指定された場合、n_unobservedを計算
  n_total_conf <- 10  # 真の交絡変数の総数
  if (!is.null(obs_ratio)) {
    n_observed_conf <- round(n_total_conf * obs_ratio)
    n_unobserved <- n_total_conf - n_observed_conf
  }

  # 観測される共変量
  X_obs <- matrix(rnorm(n * p), n, p)
  colnames(X_obs) <- paste0("X", 1:p)

  # 未観測の交絡因子
  U <- matrix(rnorm(n * n_unobserved), n, n_unobserved)

  # 観測交絡の効果
  n_observed_conf <- n_total_conf - n_unobserved
  beta_t_obs <- c(rep(0.3, n_observed_conf), rep(0, p - n_observed_conf))
  beta_y_obs <- c(rep(0.2, n_observed_conf), rep(0, p - n_observed_conf))

  # 処置と結果の両方に未観測交絡が影響
  ps <- plogis(X_obs %*% beta_t_obs +
               U %*% rep(unobserved_strength, n_unobserved))
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X_obs %*% beta_y_obs +
       U %*% rep(unobserved_strength, n_unobserved) +
       tau * T + rnorm(n)

  # バイアスの理論値（近似）
  # OVB ≈ γ_t * γ_y * var(U) / var(T) の合計
  theoretical_bias <- unobserved_strength^2 * n_unobserved

  list(
    X = X_obs, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    unobserved_info = list(
      n_unobserved = n_unobserved,
      n_observed_conf = n_observed_conf,
      obs_ratio = n_observed_conf / n_total_conf,
      unobserved_strength = unobserved_strength,
      theoretical_bias = theoretical_bias
    ),
    dgp_name = "unobserved_confounding",
    dgp_params = list(n = n, p = p, n_unobserved = n_unobserved,
                      unobserved_strength = unobserved_strength)
  )
}
```

#### S11: 合流点（拡張版）

```r
dgp_collider <- function(n = 500, p = 500,
                          n_confounders = 5,
                          n_colliders = 3,
                          collider_strength = 0.5,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 真の交絡: X[1:n_confounders]
  beta_t_conf <- c(rep(0.3, n_confounders), rep(0, p - n_confounders))
  beta_y_conf <- c(rep(0.2, n_confounders), rep(0, p - n_confounders))

  ps <- plogis(X %*% beta_t_conf)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y_conf + tau * T + rnorm(n)

  # 合流点の生成（TとYの共通結果）
  collider_start <- n_confounders + 1
  collider_end <- n_confounders + n_colliders
  for (j in 1:n_colliders) {
    col_idx <- n_confounders + j
    # 各colliderはTとYの両方から影響を受ける
    X[, col_idx] <- collider_strength * T +
                    collider_strength * Y +
                    rnorm(n, sd = 0.5)
  }

  # 変数の役割を記録
  variable_roles <- list(
    confounders = 1:n_confounders,
    colliders = collider_start:collider_end,
    neutral = (collider_end + 1):p
  )

  # collider調整時のバイアス推定（理論値）
  # colliderを調整すると、TとYの間に擬似的な関連が生じる
  collider_bias_if_adjusted <- collider_strength^2 * n_colliders

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    collider_info = list(
      n_confounders = n_confounders,
      n_colliders = n_colliders,
      collider_strength = collider_strength,
      variable_roles = variable_roles,
      collider_bias_if_adjusted = collider_bias_if_adjusted
    ),
    dgp_name = "collider",
    dgp_params = list(n = n, p = p, n_confounders = n_confounders,
                      n_colliders = n_colliders,
                      collider_strength = collider_strength)
  )
}

# S11e: 混合構造（交絡とcolliderが混在、識別困難）
dgp_collider_mixed <- function(n = 500, p = 500,
                                n_confounders = 5,
                                n_colliders = 5,
                                n_mediators = 3,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 真の交絡
  beta_t_conf <- c(rep(0.3, n_confounders), rep(0, p - n_confounders))
  beta_y_conf <- c(rep(0.2, n_confounders), rep(0, p - n_confounders))

  ps <- plogis(X %*% beta_t_conf)
  T <- rbinom(n, 1, ps)

  # 媒介変数（調整すべきでない）: T -> M -> Y
  mediator_start <- n_confounders + 1
  for (j in 1:n_mediators) {
    m_idx <- n_confounders + j
    X[, m_idx] <- 0.5 * T + rnorm(n, sd = 0.5)
  }

  # 媒介変数を通じた効果
  mediator_idx <- mediator_start:(n_confounders + n_mediators)
  Y_through_mediator <- 0.3 * rowSums(X[, mediator_idx, drop = FALSE])

  tau <- 2.0
  Y <- X %*% beta_y_conf + tau * T + Y_through_mediator + rnorm(n)

  # 合流点
  collider_start <- n_confounders + n_mediators + 1
  for (j in 1:n_colliders) {
    col_idx <- collider_start + j - 1
    X[, col_idx] <- 0.4 * T + 0.4 * Y + rnorm(n, sd = 0.5)
  }

  variable_roles <- list(
    confounders = 1:n_confounders,
    mediators = mediator_start:(n_confounders + n_mediators),
    colliders = collider_start:(collider_start + n_colliders - 1),
    neutral = (collider_start + n_colliders):p
  )

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,  # 直接効果のみ
    true_total_effect = tau + 0.5 * 0.3 * n_mediators,  # 総効果
    propensity_score = ps,
    variable_roles = variable_roles,
    dgp_name = "collider_mixed",
    dgp_params = list(n = n, p = p, n_confounders = n_confounders,
                      n_colliders = n_colliders, n_mediators = n_mediators)
  )
}
```

### 3.3 評価指標

#### ATE推定指標

| 指標 | 定義 | 説明 |
|------|------|------|
| Bias | `ATE_hat - ATE_true` | 推定バイアス |
| RMSE | `sqrt(mean((ATE_hat - ATE_true)^2))` | 二乗平均平方根誤差 |
| Coverage | `mean(CI_lower <= ATE_true & ATE_true <= CI_upper)` | 95%信頼区間被覆率 |
| CI Width | `mean(CI_upper - CI_lower)` | 信頼区間幅 |

#### ITE推定指標（対応手法のみ）

| 指標 | 定義 | 説明 |
|------|------|------|
| PEHE | `sqrt(mean((ITE_hat - ITE_true)^2))` | 異質性推定精度 |
| ITE Correlation | `cor(ITE_hat, ITE_true)` | 個別効果の相関 |
| Subgroup PEHE | サブグループ別PEHE | 層別推定精度 |

#### 異質性検出指標（S4シナリオ用）

| 指標 | 定義 | 説明 |
|------|------|------|
| 異質性検出検定力 | H₀: τᵢ = ATE の棄却率 | 異質性の有無を検出できるか |
| 質的交互作用検出率 | 有害サブグループの検出率 | 危険な患者を見つけられるか |
| 希少群カバレッジ | 希少サブグループの95%CI被覆率 | 少数派の推定精度 |

#### 計算性能指標

| 指標 | 定義 | 説明 |
|------|------|------|
| Runtime | 実行時間（秒） | 計算効率 |
| Memory | メモリ使用量（MB） | リソース要件 |
| Failure Rate | 推定失敗・エラーの割合 | 安定性 |

### 3.4 ベンチマーク実行フレームワーク

```r
# 基本的なベンチマーク実行
cf_benchmark <- function(
  dgp = c("baseline", "nonlinear", "dense_confounding",
          "correlated_confounding", "unobserved_confounding", "collider"),
  methods = c("grf", "ipw", "gformula", "hdml", "bcf", "tmle", "hdps"),
  n_rep = 100,
  metrics = c("ate_bias", "ate_rmse", "coverage", "pehe", "runtime"),
  parallel = TRUE,
  n_cores = NULL,
  seed = 42
) {

  # 実装
}

# 次元数スイープ実験
cf_benchmark_dimension_sweep <- function(
  n = 200,
  p_seq = c(50, 100, 200, 500, 1000, 2000, 5000, 10000),
  methods = c("grf", "ipw", "gformula", "hdml", "bcf", "tmle", "hdps"),
  n_rep = 50,
  ...
) {
  # 各次元数でベンチマークを実行
  # 結果を統合して限界点を特定
}

# 異質性強度スイープ実験
cf_benchmark_heterogeneity_sweep <- function(
  strength_seq = c(0, 0.3, 0.5, 1.0, 1.5, 2.0),
  methods = c("grf", "bcf", "hdml"),
  n_rep = 50,
  ...
) {
  # 各強度でのPEHEを計算
}

# overlap強度スイープ実験
cf_benchmark_overlap_sweep <- function(
  overlap_seq = c("good", "moderate", "weak", "extreme"),
  methods = c("grf", "ipw", "gformula", "hdml", "bcf", "tmle", "hdps"),
  n_rep = 50,
  ...
) {
  # 各overlap強度での性能を評価
}

# 非線形強度スイープ実験
cf_benchmark_nonlinearity_sweep <- function(
  nonlinear_type = "combined",
  strength_seq = c(0, 0.25, 0.5, 0.75, 1.0),
  methods = c("grf", "bcf", "hdml", "gformula", "tmle"),
  n_rep = 50,
  ...
) {
  results <- list()
  for (strength in strength_seq) {
    res <- cf_benchmark(
      dgp = function() dgp_nonlinear_confounding(
        nonlinear_type = nonlinear_type,
        strength = strength
      ),
      methods = methods,
      n_rep = n_rep,
      ...
    )
    res$strength <- strength
    results[[as.character(strength)]] <- res
  }
  structure(results, class = "cf_nonlinearity_sweep")
}

# 密度スイープ実験
cf_benchmark_density_sweep <- function(
  n_confounders_seq = c(10, 25, 50, 100, 200, 500),
  coef_scaling = "sqrt",
  methods = c("grf", "ipw", "hdml", "tmle", "bcf"),
  n_rep = 50,
  ...
) {
  results <- list()
  for (n_conf in n_confounders_seq) {
    res <- cf_benchmark(
      dgp = function() dgp_dense_confounding(
        n_confounders = n_conf,
        coef_scaling = coef_scaling
      ),
      methods = methods,
      n_rep = n_rep,
      ...
    )
    res$n_confounders <- n_conf
    results[[as.character(n_conf)]] <- res
  }
  structure(results, class = "cf_density_sweep")
}

# 共変量シフト強度スイープ実験
cf_benchmark_covariate_shift_sweep <- function(
  shift_type = "combined",
  magnitude_seq = c("mild", "moderate", "severe"),
  methods = c("grf", "ipw", "hdml", "tmle", "bcf"),
  n_rep = 50,
  ...
) {
  # 強度設定
  params <- list(
    mild = list(mean_shift = 0.5, var_ratio = 1.2, cor_shift = 0.1),
    moderate = list(mean_shift = 1.0, var_ratio = 1.5, cor_shift = 0.2),
    severe = list(mean_shift = 2.0, var_ratio = 2.0, cor_shift = 0.4)
  )

  results <- list()
  for (mag in magnitude_seq) {
    p <- params[[mag]]
    res <- cf_benchmark(
      dgp = function() dgp_covariate_shift(
        mean_shift = p$mean_shift,
        var_ratio = p$var_ratio,
        cor_shift = p$cor_shift
      ),
      methods = methods,
      n_rep = n_rep,
      ...
    )
    res$magnitude <- mag
    results[[mag]] <- res
  }
  structure(results, class = "cf_shift_sweep")
}

# 相関強度スイープ実験
cf_benchmark_correlation_sweep <- function(
  correlation_type = "block",
  cor_seq = c(0.1, 0.3, 0.5, 0.7, 0.9),
  methods = c("grf", "hdml", "tmle", "bcf"),
  n_rep = 50,
  ...
) {
  results <- list()
  for (cor_val in cor_seq) {
    dgp_fn <- switch(correlation_type,
      "block" = function() dgp_correlated_confounding_block(within_cor = cor_val),
      "ar1"   = function() dgp_correlated_confounding_ar1(rho = cor_val)
    )
    res <- cf_benchmark(
      dgp = dgp_fn,
      methods = methods,
      n_rep = n_rep,
      ...
    )
    res$correlation <- cor_val
    results[[as.character(cor_val)]] <- res
  }
  structure(results, class = "cf_correlation_sweep")
}

# 未観測交絡強度スイープ実験
cf_benchmark_unmeasured_sweep <- function(
  strength_seq = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0),
  n_unobserved = 5,
  methods = c("grf", "ipw", "hdml", "tmle", "bcf"),
  n_rep = 50,
  ...
) {
  results <- list()
  for (strength in strength_seq) {
    res <- cf_benchmark(
      dgp = function() dgp_unobserved_confounding(
        n_unobserved = n_unobserved,
        unobserved_strength = strength
      ),
      methods = methods,
      n_rep = n_rep,
      ...
    )
    res$strength <- strength
    res$expected_bias <- strength^2 * n_unobserved
    results[[as.character(strength)]] <- res
  }
  structure(results, class = "cf_unmeasured_sweep")
}

# Collider強度スイープ実験
cf_benchmark_collider_sweep <- function(
  strength_seq = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0),
  n_colliders = 3,
  methods = c("grf", "hdml", "tmle", "bcf"),
  n_rep = 50,
  compare_adjustment = TRUE,
  ...
) {
  results <- list()
  for (strength in strength_seq) {
    res <- cf_benchmark(
      dgp = function() dgp_collider(
        n_colliders = n_colliders,
        collider_strength = strength
      ),
      methods = methods,
      n_rep = n_rep,
      ...
    )
    res$strength <- strength
    res$expected_bias_if_adjusted <- strength^2 * n_colliders
    results[[as.character(strength)]] <- res
  }
  structure(results, class = "cf_collider_sweep")
}
```

### 3.5 結果の可視化

```r
# ヒートマップ（シナリオ × 手法）
cf_benchmark_plot(result, type = "heatmap", metric = "ate_rmse")

# ボックスプロット（手法間比較）
cf_benchmark_plot(result, type = "boxplot", metric = "ate_bias")

# 次元数スイープ曲線
plot.cf_dimension_sweep <- function(x, metric = "ate_rmse") {
  # X軸: p（対数スケール）
  # Y軸: 指定指標
  # 各手法の性能曲線と限界点を表示
}

# 限界点サマリー
cf_find_breaking_point <- function(sweep_result, threshold = 0.5) {
  # 各手法が閾値を超えるpを特定
}
```

---

## 4. 手法の前提条件と選択ガイドライン

### 4.1 前提条件比較表

| 手法 | 線形性 | 疎性 | overlap | n要件 | p要件 | 分布仮定 |
|------|--------|------|---------|-------|-------|---------|
| GRF | 不要 | 不要 | 中程度 | 中〜大 | 低〜中 | なし |
| IPW | 傾向のみ | 不要 | 厳格 | 中 | 低 | なし |
| G-formula | 必要 | 不要 | 緩い | 小〜中 | 低 | 正規 |
| HDML | 必要 | 必要 | 中程度 | 小〜中 | 高 | 正規 |
| BCF | 不要 | 不要 | 中程度 | 中 | 中〜高 | なし |
| TMLE | 不要 | 不要 | 緩い | 中〜大 | 中〜高 | なし |
| HDPS | 傾向のみ | 必要 | 厳格 | 中 | 高 | なし |

### 4.2 計算効率と適用スケール

| 手法 | 時間計算量 | 空間計算量 | 推奨n | 推奨p | 実行時間目安* |
|------|-----------|-----------|-------|-------|--------------|
| GRF | O(n log n × p) | O(n × p) | 100-10,000 | 10-500 | 秒〜分 |
| IPW | O(n × p) | O(n × p) | 50-5,000 | 10-100 | 秒 |
| G-formula | O(n × p²) | O(p²) | 50-1,000 | 10-50 | 秒 |
| HDML | O(n × p) | O(n × p) | 50-1,000 | 100-10,000 | 秒〜分 |
| BCF | O(n × p × iter) | O(n × p) | 100-2,000 | 50-1,000 | 分〜時間 |
| TMLE | O(n × p × |lib|) | O(n × p) | 200-5,000 | 50-2,000 | 分 |
| HDPS | O(n × p) | O(n × p) | 100-5,000 | 100-10,000 | 秒〜分 |

*n=500, p=500での目安

### 4.3 手法選択フローチャート

```
データを取得
    │
    ▼
┌─────────────────┐
│ p > n か？       │
└────────┬────────┘
         │
    ┌────┴────┐
    Yes       No
    │         │
    ▼         ▼
高次元手法    標準手法候補
(hdml,bcf,   (grf,ipw,
 tmle,hdps)   gformula)
    │         │
    ▼         ▼
┌─────────────────┐
│ ITE推定が必要？  │
└────────┬────────┘
    ┌────┴────┐
    Yes       No
    │         │
    ▼         ▼
 bcf,grf    hdml,tmle,
            ipw,hdps
    │         │
    ▼         ▼
┌─────────────────┐
│ 計算時間の制約？ │
└────────┬────────┘
    ┌────┴────┐
   厳しい    余裕あり
    │         │
    ▼         ▼
 hdml,hdps   bcf,tmle
 grf,ipw
```

### 4.4 診断・推奨機能

```r
# データ特性の診断
cf_diagnose_data <- function(data, formula) {
  # n, p, n/p比
  # 疎性の推定
  # overlap診断
  # 非線形性の検出
}

# 手法の推奨
cf_recommend_method <- function(data, formula) {
  diagnosis <- cf_diagnose_data(data, formula)

  # データ特性に基づく推奨
  # 推奨理由の説明
  # 非推奨手法と理由
  # 想定実行時間
}

# インタラクティブガイド
cf_guide <- function(data, formula) {
  # cf_diagnose_data + cf_recommend_method + 可視化
}
```

#### 出力例

```
=== cfomics 手法選択ガイド ===

データ特性:
  サンプル数 (n): 200
  変数数 (p): 1500
  n/p 比: 0.13 (高次元)
  overlap: 良好 (PS range: 0.15-0.85)

推奨手法 (優先順):
  1. hdml  - 高次元に最適化、高速
  2. bcf   - ITE推定も必要な場合（計算時間: 長）
  3. tmle  - 理論的保証重視の場合

非推奨:
  - ipw, gformula: p >> n で推定不安定
  - grf: 超高次元では精度低下の可能性

想定実行時間:
  hdml: ~30秒, bcf: ~15分, tmle: ~3分
```

### 4.5 overlap診断

```r
cf_diagnose_overlap <- function(data, formula) {
  ps <- estimate_propensity(data, formula)

  # 傾向スコア分布プロット
  # overlap領域のサマリー
  # 推奨事項（トリミング等）
}
```

---

## 5. 実装計画

### 5.1 フェーズ1: 新規手法の追加

1. `methods_hdml.R` - Debiased Lasso / Double ML
2. `methods_bcf.R` - Bayesian Causal Forests
3. `methods_tmle.R` - TMLE + Super Learner
4. `methods_hdps.R` - HD-Propensity Score
5. 各手法のテスト

### 5.2 フェーズ2: ベンチマークフレームワーク

1. `benchmark_dgp.R` - 全DGPの実装
2. `benchmark_run.R` - ベンチマーク実行エンジン
3. `benchmark_viz.R` - 可視化機能
4. 統合テスト

### 5.3 フェーズ3: ガイダンス機能

1. `cf_diagnose_data()` - データ診断
2. `cf_recommend_method()` - 手法推奨
3. `cf_guide()` - インタラクティブガイド
4. ドキュメンテーション

### 5.4 フェーズ4: ベンチマーク実行と論文化

1. 全シナリオでのベンチマーク実行
2. 結果の分析
3. 論文/ビネット作成

---

## 6. 依存パッケージ

### 6.1 新規依存（Suggests）

```
hdm,
DoubleML,
bcf,
dbarts,
tmle3,
sl3,
MASS
```

### 6.2 既存依存の確認

```
grf,
ipw,
survey,
glmnet,
ggplot2
```

---

## 付録A: シナリオ一覧（完全版）

### 基本シナリオ

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S1 | 基準 | n=500, p=50, 線形・疎 | ベースライン |

### S2-S3: 次元数スイープ

| 設計 | 目的 |
|------|------|
| n=200, p=50-10000（8点） | 手法の次元限界を特定 |

### S4: 異質な処置効果（8サブシナリオ）

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S4a | 線形異質性 | 強度=中 | 基本的HTE |
| S4b | 非線形異質性 | sin, 二乗項 | 非線形パターン |
| S4c | サブグループ異質性 | 3グループ | 離散的異質性 |
| S4d | 質的交互作用 | 一部で害 | 有害検出 |
| S4e | 線形異質性（弱） | 強度=弱 | 検出力 |
| S4f | 線形異質性（強） | 強度=強 | 強異質性 |
| S4g | サブグループ（偏り中） | 60/30/10% | 中程度偏り |
| S4h | サブグループ（偏り大） | 80/15/5% | 希少サブグループ |

### S5: 非線形交絡（5サブシナリオ + スイープ）

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S5a | 二乗項 | Y, T ∝ X² | 単純な非線形 |
| S5b | 三角関数 | Y, T ∝ sin(πX) | 周期的パターン |
| S5c | 交互作用 | Y, T ∝ X₁·X₂ | 変数間相互作用 |
| S5d | 複合非線形 | 上記の組み合わせ | 複雑なパターン |
| S5e | 閾値効果 | Y, T ∝ I(X > c) | ステップ関数的 |
| スイープ | 非線形強度 | strength=0-1.0（5点） | 強度の影響 |

### S6: 密交絡（スイープ）

| 設計 | 目的 |
|------|------|
| n_confounders=10-500（6点） | 交絡変数数の限界 |
| coef_scaling=fixed/sqrt/linear | スケーリング方式の影響 |

### S7: 弱overlap（6サブシナリオ）

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S7a-good | 係数スケール | PS∈[0.1, 0.9] | ベースライン |
| S7a-moderate | 係数スケール | PS∈[0.05, 0.95] | 軽度違反 |
| S7a-weak | 係数スケール | PS∈[0.01, 0.99] | 実質的違反 |
| S7a-extreme | 係数スケール | PS∈[0.001, 0.999] | 重度違反 |
| S7b | 構造的違反 | 一部でPS≈0 | 完全違反 |
| S7c | 非対称overlap | 片側極端 | 群間不均衡 |

### S8: 共変量シフト（4サブシナリオ + スイープ）

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S8a | 平均シフトのみ | mean_shift変動 | 位置のずれ |
| S8b | 分散シフトのみ | var_ratio変動 | スケールのずれ |
| S8c | 相関シフトのみ | cor_shift変動 | 構造のずれ |
| S8d | 複合シフト | 全パラメータ変動 | 現実的状況 |
| スイープ | シフト強度 | mild/moderate/severe | 強度の影響 |

### S9: 相関交絡（5サブシナリオ + スイープ）

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S9a | ブロック相関（弱） | within_cor=0.3 | 弱い共発現 |
| S9b | ブロック相関（中） | within_cor=0.5 | 中程度共発現 |
| S9c | ブロック相関（強） | within_cor=0.7 | 強い共発現 |
| S9d | AR(1)構造 | rho=0.3-0.7 | 順序依存 |
| S9e | 因子モデル | n_factors=5-20 | 潜在因子構造 |
| スイープ | 相関強度 | cor=0.1-0.9（5点） | 強度の影響 |

### S10: 未観測交絡（5サブシナリオ + スイープ）

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S10a | 弱い未観測交絡 | strength=0.2 | 軽微なバイアス |
| S10b | 中程度 | strength=0.5 | 中程度バイアス |
| S10c | 強い未観測交絡 | strength=0.8 | 深刻なバイアス |
| S10d | 多数の未観測 | n_unobserved=1-20 | 未観測変数の数 |
| S10e | 観測割合 | obs_ratio=0.5-0.9 | 交絡の観測割合 |
| スイープ | 未観測強度 | strength=0.1-1.0（6点） | 強度の影響 |

### S11: 合流点（5サブシナリオ + スイープ）

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S11a | 弱いcollider | effect=0.2 | 軽微な選択バイアス |
| S11b | 中程度 | effect=0.5 | 中程度バイアス |
| S11c | 強いcollider | effect=0.8 | 深刻なバイアス |
| S11d | collider数 | n_colliders=1-10 | 複数のcollider |
| S11e | 混合構造 | 交絡+媒介+collider | 変数選択の難しさ |
| スイープ | collider強度 | effect=0.1-1.0（6点） | 強度の影響 |

### シナリオ数サマリー

| カテゴリ | サブシナリオ数 | スイープ実験 |
|---------|--------------|-------------|
| S1 | 1 | - |
| S2-S3 | - | 次元数（8点） |
| S4 | 8 | 異質性強度 |
| S5 | 5 | 非線形強度（5点） |
| S6 | - | 密度（6点）× スケーリング（3種） |
| S7 | 6 | overlap（4点） |
| S8 | 4 | シフト強度（3点） |
| S9 | 5 | 相関強度（5点） |
| S10 | 5 | 未観測強度（6点） |
| S11 | 5 | collider強度（6点） |
| **合計** | **39** | **9種のスイープ** |

---

## 付録B: 実応用シナリオ対応表

| シナリオ | 実応用例 |
|---------|---------|
| S1 | 臨床試験補助解析、少数バイオマーカー研究 |
| S2-S3 | RNA-seq、プロテオミクス、マルチオミクス統合 |
| S4a-h | 個別化医療、薬剤反応性、サブタイプ別効果、希少レスポンダー |
| S5a | 用量反応の二次効果 |
| S5b | 概日リズム、周期的生物プロセス |
| S5c | パスウェイ相互作用、シグナル伝達クロストーク |
| S5d | 複雑な生物学的ネットワーク |
| S5e | バイオマーカー閾値、臨床的カットオフ |
| S6 | 環境疫学、多数の生活習慣因子、バッチ効果 |
| S7a-c | 希少疾患、禁忌条件、強い適応基準 |
| S8a | 異なる母集団間の比較 |
| S8b | 測定プラットフォームの違い |
| S8c | 技術的変動、バッチ効果 |
| S8d | 異なるコホート統合、メタ解析 |
| S9a-c | 遺伝子共発現ネットワーク（弱〜強） |
| S9d | ゲノム位置依存の相関、連鎖不平衡 |
| S9e | 転写因子調節ネットワーク、潜在クラスター |
| S10a-c | バッチ効果（軽度〜重度） |
| S10d | 複合的な未測定変数 |
| S10e | 部分的オミクス測定、コスト制約 |
| S11a-c | 入院バイアス、生存者バイアス（軽度〜重度） |
| S11d | 多段階の選択プロセス |
| S11e | 複雑な因果構造、観察研究の典型例 |

## 付録C: スイープ実験一覧

| スイープ名 | パラメータ | 値の範囲 | 目的 |
|-----------|-----------|---------|------|
| dimension_sweep | p | 50-10000（8点） | 次元限界の特定 |
| heterogeneity_sweep | strength | 0-2.0（6点） | 異質性検出力 |
| nonlinearity_sweep | strength | 0-1.0（5点） | 非線形対応力 |
| density_sweep | n_confounders | 10-500（6点） | 密交絡対応力 |
| overlap_sweep | overlap_strength | good-extreme（4点） | positivity違反耐性 |
| shift_sweep | magnitude | mild-severe（3点） | 分布シフト耐性 |
| correlation_sweep | within_cor | 0.1-0.9（5点） | 相関構造対応力 |
| unmeasured_sweep | strength | 0.1-1.0（6点） | 未観測交絡のバイアス |
| collider_sweep | effect | 0.1-1.0（6点） | collider調整の影響 |
