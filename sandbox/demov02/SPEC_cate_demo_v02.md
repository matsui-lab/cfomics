# CATE Demo v0.2 — 仕様書

> **スクリプト**: `sandbox/run_cate_demo_v02.R`
> **実行日**: 2026-02-17
> **目的**: 高次元オミクス模擬データ（RNA-seq 50遺伝子）上でCATE推定手法4種のベンチマーク比較

---

## 1. 実験条件

| 設定項目 | 値 |
|---|---|
| サンプル数 n | 100, 500, 1000 |
| 効果タイプ | linear, threshold3 |
| 手法 | causal_forest, drlearner, ganite, cevae |
| 入力前処理 | log1p_z（z-score of log1p counts） |
| 分割比（train:test） | 80:20（層化ランダム分割 by T） |
| seed_base | 20260215 |
| seed_condition | `seed_base + 10000*(idx_n−1) + 100*(idx_effect−1)` |

---

## 2. データ生成過程（DGP）

### 2.1 概要

二層構造：潜在ログ発現量 → NB カウント観測（RNA-seq模擬）

```
潜在ログ発現: E[i,j] = b0_j + bA_j×A_i + Σ_k(Λ_jk×G_ik) + ε_ij
観測カウント: C[i,j] ~ NegBinom(size=θ_NB, mu=exp(E[i,j])×s_i)
```

### 2.2 パラメータ

| パラメータ | 記号 | 値 | 説明 |
|---|---|---|---|
| 遺伝子数 | p | 50 | |
| 潜在因子数 | K | 5 | 遺伝子間相関の源 |
| 因子ローディングSD | σ_Λ | 0.3 | |
| 発現ノイズSD | σ_E | 0.3 | |
| NB過分散パラメータ | θ_NB | 10 | size パラメータ（大きいほどPoisson寄り） |
| ライブラリサイズSD | σ_s | 0.2 | log正規分布 |
| 予測遺伝子（効果修飾） | M_genes | gene 1–5 | τ を決める遺伝子群 |
| 予後遺伝子 | P_genes | gene 6–10 | Y(0)を決める遺伝子群 |
| ノイズ遺伝子 | — | gene 11–50 | CATEにもY(0)にも無関係 |

### 2.3 共変量

| 変数 | 定義 |
|---|---|
| A_i | 年齢（標準正規分布） |
| T_i | 治療割り付け（Bernoulli(0.5)、RCT） |
| G_i | n×K 潜在因子（標準正規） |

### 2.4 効果修飾スコア u

```
u_i = Σ_{j∈M_genes} E[i,j]   （予測遺伝子の潜在ログ発現の合計）
```

※ モデルへの入力はNBカウントの前処理値であり、u は直接観測不可能

### 2.5 真のCATE（τ）

**Linear effect**

```
τ_i = (u_i − mean(u)) × (s_tau / SD(u_centered))
```

**Threshold3 effect**

```
τ_i = -s_tau   if u_i ≤ Q_{0.33}(u)
      0         if Q_{0.33} < u_i ≤ Q_{0.66}(u)
      +s_tau    if u_i > Q_{0.66}(u)
```

- SD正規化によりいずれも `SD(τ) = s_tau = 5` に固定

### 2.6 予後スコア（η）と観測アウトカム（Y）

```
η_i = (β_A×A_i + β_P×vc_i) × (s_eta / SD(η_raw))

Y(0)_i = μ_0 + η_i + ε_i    （ε_i ~ N(0, σ_y²)）
Y(1)_i = Y(0)_i + τ_i
Y_i    = Y(0)_i + T_i × τ_i  （観測値）
```

| パラメータ | 値 |
|---|---|
| μ_0（ベースライン） | 120 |
| β_A（年齢効果） | 3 |
| β_P（予後遺伝子効果） | 1 |
| s_eta（予後SDスケール） | 10 |
| σ_y（残差SD） | 5 |

---

## 3. 前処理（log1p_z）

1. 生カウントC を log1p 変換 → L
2. **訓練セットの統計量**（遺伝子ごとの mean・SD）を計算
3. 全サンプル（train + test）にその統計量でz-score正規化（情報漏洩防止）
4. モデル入力 X = `[A, z-score(L_1), ..., z-score(L_50)]`（次元: n × 51）

---

## 4. 手法と設定

### 4.1 Causal Forest（GRF）

| パラメータ | 値 |
|---|---|
| num.trees | 2000 |
| random_state | seed_condition |
| 実装 | cfomics::cf_fit(method="grf") |

### 4.2 DR-Learner（EconML）

| パラメータ | 値 |
|---|---|
| CV fold数 | 3（n_train<200）/ 5（n_train≥200） |
| random_state | seed_condition |
| 実装 | cfomics::cf_fit(method="drlearner") |

### 4.3 GANITE（TensorFlow）

| パラメータ | 値 |
|---|---|
| iterations | 500（n_train<200）/ 1000（n_train≥200） |
| batch_size | max(16, min(n_train÷4, 64)) |
| h_dim | 64 |
| alpha | 1.0 |
| 実装 | cfomics Python backend（GANITEModel 直接操作） |
| 備考 | cfomics ラッパーは予測未対応のため GANITEModel クラスを直接使用 |

### 4.4 CEVAE（PyTorch / Pyro）

| パラメータ | 値 |
|---|---|
| num_epochs | 20 |
| batch_size | min(n_train, 64) |
| latent_dim | 5 |
| hidden_dim | 64 |
| num_layers | 2 |
| learning_rate | 1e-3 |
| learning_rate_decay | 0.1 |
| weight_decay | 1e-4 |
| num_samples | 10 |
| 実装 | pyro.contrib.cevae.CEVAE 直接使用 |

> **注**: 全手法とも CPU 実行を前提とした保守的な設定

---

## 5. 評価指標

```
RMSE_τ = sqrt( mean( (τ̂_i − τ_i)² ) )    （テストセットのみ）
```

- τ の SD = 5（s_tau）が自然な基準。RMSE < 5 でτのSDを下回る精度

---

## 6. 結果

### RMSE_τ（テストセット）

**Linear effect**

| 手法 | n=100 | n=500 | n=1000 |
|---|---|---|---|
| Causal Forest | NA ※ | 3.87 | **2.94** |
| DR-Learner | 11.55 | 5.14 | 3.27 |
| GANITE | 10.59 | 13.09 | 10.32 |
| CEVAE | 6.03 | 4.91 | 4.99 |

**Threshold3 effect**

| 手法 | n=100 | n=500 | n=1000 |
|---|---|---|---|
| Causal Forest | 4.29 | 4.30 | **4.08** |
| DR-Learner | 17.08 | 5.94 | 4.99 |
| GANITE | 13.64 | 10.03 | 11.98 |
| CEVAE | 5.40 | 4.80 | 4.98 |

※ Causal Forest n=100 linear: honesty fraction 問題により NA

### 所見

- **Causal Forest**: linear では n 増加とともに最良の精度。threshold3 でも安定
- **DR-Learner**: n=100 で不安定（特に threshold3）、n≥500 で急改善
- **GANITE**: 全条件で RMSE > 10。CPU・少イテレーション設定での限界
- **CEVAE**: n=100 から比較的安定（RMSE ≈ 5–6）。linear/threshold3 いずれも均一

---

## 7. 出力ファイル

```
sandbox/benchmark_results/cate_demo_v02/
├── results_demo_omics50_v02.csv         # 24行：手法×条件ごとのRMSE
├── predictions_demo_omics50_v02.csv     # 2560行：テスト個体ごとのτ_true/τ̂
├── session_info.txt                     # R/Python環境情報
├── dgp_viz_01_age_treatment.png         # 年齢・治療の分布
├── dgp_viz_02_count_distributions.png  # 代表遺伝子のカウント分布
├── dgp_viz_03_all_genes_boxplot.png    # 全50遺伝子のlog1p箱ひげ図
├── dgp_viz_04_zscore_input.png         # z-score後の分布（モデル入力）
├── dgp_viz_05_outcome_tau.png          # 深度因子・τ・Y の分布
├── dgp_viz_06_gene_correlation.png     # 遺伝子間Spearman相関ヒートマップ
├── dgp_viz_07_tau_distribution.png     # τ分布（linear vs threshold3）
├── dgp_viz_08_outcome_by_treatment.png # Y分布（治療群別）
├── dgp_viz_09_tau_vs_Y0.png           # τ vs Y(0) 散布図
├── dgp_viz_10_combined_tau_outcome.png # τ・Y 統合4パネル図
├── dgp_viz_11_residuals_*.png          # 残差プロット（手法別、4枚）
└── dgp_viz_12_u_vs_cate_*.png         # u vs τ/τ̂ 散布図（手法別、4枚）
```

---

## 8. 実行方法

```bash
# リポジトリルートから実行
Rscript sandbox/run_cate_demo_v02.R

# 可視化（個別実行）
Rscript sandbox/visualize_dgp.R
Rscript sandbox/visualize_cate_outcome.R
Rscript sandbox/visualize_residuals.R
Rscript sandbox/visualize_u_vs_cate.R
```

**前提条件**:
- `sandbox/.venv` に Python 仮想環境が存在すること（`sandbox/setup_env.R` で構築）
- R パッケージ: `cfomics`, `ggplot2`, `patchwork`
- Python パッケージ: `torch`, `pyro-ppl`, `econml`, `tensorflow`
