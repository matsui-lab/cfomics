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
| S5 | 非線形 | 500 | 500 | 非線形交絡 | モデル柔軟性 | 複雑なパスウェイ相互作用、用量反応関係 |
| S6 | 密交絡 | 500 | 500 | 線形・密・独立 | 多数の交絡 | 環境疫学、バッチ効果 |
| S9 | 相関交絡 | 500 | 500 | ブロック相関構造 | 共発現パターン | 遺伝子共発現、パスウェイ内遺伝子群 |
| S10 | 未観測交絡 | 500 | 500 | 一部が観測不能 | 感度分析評価 | バッチ効果、測定されない生活習慣 |
| S11 | 合流点 | 500 | 500 | collider含む | 不適切な調整の検出 | 選択バイアス、生存者バイアス |

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

#### S8: 共変量シフト

| パラメータ | 値 |
|-----------|-----|
| n | 500 |
| p | 500 |
| 平均シフト | 1.0（処置群で一部変数が+1.0） |
| 分散比 | 1.5（処置群で分散が1.5倍） |
| 相関シフト | 0.2（処置群で相関が+0.2） |
| 目的 | 分布の不一致への対処 |
| 実応用例 | 異なるコホート間比較、バッチ効果、前向きvs後ろ向き研究 |

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

#### S5: 非線形交絡

```r
dgp_nonlinear_confounding <- function(n = 500, p = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 非線形な傾向スコアモデル
  ps_linear <- 0.3 * X[,1] + 0.3 * X[,2]
  ps_nonlinear <- 0.5 * X[,1]^2 + 0.3 * sin(pi * X[,3]) + 0.2 * X[,1] * X[,2]
  ps <- plogis(ps_linear + ps_nonlinear)
  T <- rbinom(n, 1, ps)

  # 非線形な結果モデル
  tau <- 2.0
  Y_baseline <- 0.3 * X[,1]^2 + 0.2 * sin(pi * X[,2]) + 0.1 * X[,1] * X[,3]
  Y <- Y_baseline + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    dgp_name = "nonlinear_confounding",
    dgp_params = list(n = n, p = p)
  )
}
```

#### S6: 密交絡

```r
dgp_dense_confounding <- function(n = 500, p = 500,
                                   n_confounders = 100,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 多数の交絡変数（係数は小さいが多い）
  coef_size <- 0.1
  beta_t <- c(rep(coef_size, n_confounders), rep(0, p - n_confounders))
  beta_y <- c(rep(coef_size, n_confounders), rep(0, p - n_confounders))

  ps <- plogis(X %*% beta_t)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y + tau * T + rnorm(n)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    n_confounders = n_confounders,
    dgp_name = "dense_confounding",
    dgp_params = list(n = n, p = p, n_confounders = n_confounders)
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

#### S9: 相関交絡

```r
dgp_correlated_confounding <- function(n = 500, p = 500,
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
    block_structure = list(
      n_blocks = n_blocks,
      block_size = block_size,
      within_cor = within_cor
    ),
    dgp_name = "correlated_confounding",
    dgp_params = list(n = n, p = p, n_blocks = n_blocks, within_cor = within_cor)
  )
}
```

#### S10: 未観測交絡

```r
dgp_unobserved_confounding <- function(n = 500, p = 500,
                                        n_unobserved = 5,
                                        unobserved_strength = 0.5,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 観測される共変量
  X_obs <- matrix(rnorm(n * p), n, p)
  colnames(X_obs) <- paste0("X", 1:p)

  # 未観測の交絡因子
  U <- matrix(rnorm(n * n_unobserved), n, n_unobserved)

  # 観測交絡の効果
  beta_t_obs <- c(rep(0.3, 10), rep(0, p - 10))
  beta_y_obs <- c(rep(0.2, 10), rep(0, p - 10))

  # 処置と結果の両方に未観測交絡が影響
  ps <- plogis(X_obs %*% beta_t_obs +
               U %*% rep(unobserved_strength, n_unobserved))
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X_obs %*% beta_y_obs +
       U %*% rep(unobserved_strength, n_unobserved) +
       tau * T + rnorm(n)

  # バイアスの理論値（近似）
  theoretical_bias <- unobserved_strength^2 * n_unobserved

  list(
    X = X_obs, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    unobserved_info = list(
      n_unobserved = n_unobserved,
      unobserved_strength = unobserved_strength,
      theoretical_bias = theoretical_bias
    ),
    dgp_name = "unobserved_confounding",
    dgp_params = list(n = n, p = p, n_unobserved = n_unobserved,
                      unobserved_strength = unobserved_strength)
  )
}
```

#### S11: 合流点（collider）

```r
dgp_collider <- function(n = 500, p = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)

  # 真の交絡: X1-X5
  beta_t_conf <- c(rep(0.3, 5), rep(0, p - 5))
  beta_y_conf <- c(rep(0.2, 5), rep(0, p - 5))

  ps <- plogis(X %*% beta_t_conf)
  T <- rbinom(n, 1, ps)

  tau <- 2.0
  Y <- X %*% beta_y_conf + tau * T + rnorm(n)

  # 合流点: X6-X10はTとYの共通の結果（調整すべきでない）
  # X6 = f(T, Y) + noise
  X[, 6] <- 0.5 * T + 0.3 * Y + rnorm(n, sd = 0.5)
  X[, 7] <- 0.4 * T + 0.4 * Y + rnorm(n, sd = 0.5)
  X[, 8] <- 0.3 * T + 0.5 * Y + rnorm(n, sd = 0.5)

  list(
    X = X, T = T, Y = Y,
    true_ate = tau,
    true_ite = rep(tau, n),
    propensity_score = ps,
    collider_info = list(
      confounder_idx = 1:5,
      collider_idx = 6:8,
      neutral_idx = 9:p
    ),
    dgp_name = "collider",
    dgp_params = list(n = n, p = p)
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

| ID | 名称 | 設計 | 目的 |
|----|------|------|------|
| S1 | 基準 | n=500, p=50, 線形・疎 | ベースライン |
| S2-S3 | 次元数スイープ | n=200, p=50-10000 | 手法の次元限界 |
| S4a | 線形異質性 | 強度=中 | 基本的HTE |
| S4b | 非線形異質性 | sin, 二乗項 | 非線形パターン |
| S4c | サブグループ異質性 | 3グループ | 離散的異質性 |
| S4d | 質的交互作用 | 一部で害 | 有害検出 |
| S4e | 線形異質性（弱） | 強度=弱 | 検出力 |
| S4f | 線形異質性（強） | 強度=強 | 強異質性 |
| S4g | サブグループ（偏り中） | 60/30/10% | 中程度偏り |
| S4h | サブグループ（偏り大） | 80/15/5% | 希少サブグループ |
| S5 | 非線形交絡 | sin, 二乗, 交互作用 | モデル柔軟性 |
| S6 | 密交絡 | 100変数が交絡 | 多数交絡 |
| S7a | 弱overlap（係数） | good/moderate/weak/extreme | 段階的違反 |
| S7b | 構造的違反 | 一部でPS≈0 | 完全違反 |
| S7c | 非対称overlap | 片側極端 | 群間不均衡 |
| S8 | 共変量シフト | 平均・分散・相関シフト | 分布不一致 |
| S9 | 相関交絡 | ブロック相関 | 共発現 |
| S10 | 未観測交絡 | 一部観測不能 | 感度分析 |
| S11 | 合流点 | collider含む | 不適切調整 |

---

## 付録B: 実応用シナリオ対応表

| シナリオ | 実応用例 |
|---------|---------|
| S1 | 臨床試験補助解析、少数バイオマーカー研究 |
| S2-S3 | RNA-seq、プロテオミクス、マルチオミクス統合 |
| S4a-h | 個別化医療、薬剤反応性、サブタイプ別効果 |
| S5 | パスウェイ相互作用、用量反応関係 |
| S6 | 環境疫学、バッチ効果 |
| S7a-c | 希少疾患、禁忌条件、強い適応基準 |
| S8 | 異なるコホート比較、前向きvs後ろ向き研究 |
| S9 | 遺伝子共発現、パスウェイ内遺伝子群 |
| S10 | バッチ効果、測定されない生活習慣 |
| S11 | 選択バイアス、生存者バイアス |
