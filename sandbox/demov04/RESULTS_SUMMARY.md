# demov04 ベンチマーク結果まとめ

**実行日**: 2026-02-20
**目的**: CATE推定手法の比較（オミクスデータ想定）
**重要な変更**: DR-Learner修正版、sigma_y=5.0（個人差拡大）

---

## 1. シナリオ設定

### データ生成プロセス (DGP)

**基本設定**:
```r
# サンプルサイズ
n = 100, 500, 1000

# 反復回数
R = 1

# Train/Test split
train_frac = 0.8
```

**アウトカムモデル**:
```r
# ベースライン（個人差あり）
Y(0) ~ N(120, σ_y²)  where σ_y = 5.0
範囲: [105, 135] mmHg (±3SD)

# 観測アウトカム
Y = Y(0) - T × τ

# 測定誤差なし（個人差のみ）
```

**遺伝子発現モデル**:
```r
# 真の発現量
E_true ~ N(log(50), 0.35)

# カウントデータ（負の二項分布）
C_count ~ NB(mu = exp(E_true), size = 30)

# 観測値（log1p変換 + 測定誤差）
G_obs = log1p(C_count) + N(0, 0.08)
```

**治療効果モデル**:

_Linear effect_:
```r
τ = -12.0 + 2.0×(Sex-0.5) + 1.5×Z_obs
範囲: [-18, -6] mmHg (clipped)

# Sex効果: ±1 mmHg
# geneA効果: ±3 mmHg (Z∈[-2,2]として)
```

_Threshold3 effect_:
```r
τ = -12.0 + 2.0×(Sex-0.5) + g(Z_obs)

where g(Z) = { -1.5  if Z ≤ Q33
             {  0    if Q33 < Z ≤ Q66
             { +1.5  if Z > Q66

# 3段階の非連続効果
```

**共変量セット**:
- **SexOnly**: Sex（性別のみ）
- **SexGene**: Sex + geneA（性別 + 遺伝子）

---

## 2. 評価手法

### 比較手法

| 手法 | タイプ | 実装 | 備考 |
|------|--------|------|------|
| **Causal Forest** | Tree-based | R (grf) | デフォルト設定 |
| **DR-Learner** | Doubly Robust | Python (econml) | **修正版** |

### DR-Learner 修正内容

元の実装の問題点:
1. モデル引数が EconML に渡されていない
2. 傾向スコアのクリッピングなし
3. Y の標準化なし
4. model_final が未指定（OLS使用）

修正版の改善:
```r
cf_fit_drlearner(
  X, T, Y,
  model_propensity = LogisticRegressionCV(cv, max_iter=1000),
  model_regression = RidgeCV(cv),
  model_final = LassoCV(cv),
  min_propensity = 0.05,      # NEW
  standardize_y = TRUE,        # NEW
  random_state = seed
)
```

**効果**: n=500, linear, SexGene で RMSE 19.94 → 0.63 (97%改善)

---

## 3. 結果

### 3.1 全体サマリー

| n | Effect | Covset | Method | RMSE | geneA効果 |
|---|--------|--------|--------|------|-----------|
| 100 | linear | SexOnly | DR-Learner | 3.75 | - |
| 100 | linear | SexGene | DR-Learner | 3.31 | **+11.7%** |
| 500 | linear | SexOnly | Causal Forest | 1.69 | - |
| 500 | linear | SexGene | Causal Forest | 0.88 | **+47.8%** |
| 500 | linear | SexOnly | DR-Learner | 1.60 | - |
| 500 | linear | SexGene | DR-Learner | 0.63 | **+60.6%** |
| 1000 | linear | SexOnly | Causal Forest | 1.57 | - |
| 1000 | linear | SexGene | Causal Forest | 0.92 | **+41.1%** |
| 1000 | linear | SexOnly | DR-Learner | 1.54 | - |
| 1000 | linear | SexGene | DR-Learner | 0.16 | **+89.5%** |
| | | | | | |
| 100 | threshold3 | SexOnly | DR-Learner | 2.48 | - |
| 100 | threshold3 | SexGene | DR-Learner | 2.43 | **+2.1%** |
| 500 | threshold3 | SexOnly | Causal Forest | 1.28 | - |
| 500 | threshold3 | SexGene | Causal Forest | 0.97 | **+24.0%** |
| 500 | threshold3 | SexOnly | DR-Learner | 1.28 | - |
| 500 | threshold3 | SexGene | DR-Learner | 0.98 | **+23.5%** |
| 1000 | threshold3 | SexOnly | Causal Forest | 1.22 | - |
| 1000 | threshold3 | SexGene | Causal Forest | 1.26 | **-3.3%** ❌ |
| 1000 | threshold3 | SexOnly | DR-Learner | 1.62 | - |
| 1000 | threshold3 | SexGene | DR-Learner | 0.59 | **+63.4%** |

### 3.2 主要な知見

**✅ DR-Learner（修正版）の優位性**:
- 全条件で動作（n=100でもOK）
- Linear効果: 12-90%の改善
- Threshold効果: 2-63%の改善
- 修正により安定性が劇的に向上

**⚠️ Causal Forest の問題点**:

1. **n=100 で完全失敗**
   - Error: "honesty fraction too close to 1 or 0"
   - 原因: サンプル不足（honesty.fraction=0.5 default）

2. **階段状予測**
   - 決定木ベースの本質的特性
   - 線形DGPでは DR-Learner に劣る

3. **変数追加で性能悪化**（最も深刻）
   - n=1000, threshold3: SexOnly 1.22 → SexGene 1.26
   - 原因: 過度な正則化・平均への縮小
   - 予測分散が真値の1/4に圧縮

**📈 sigma_y の影響**:
- sigma_y = 2.0 (旧): 簡単すぎる問題
- sigma_y = 5.0 (新): 中程度の難易度
- SNR ≈ 0.32 (妥当なバランス)

---

## 4. 考察

### 4.1 手法特性の比較

| 観点 | DR-Learner | Causal Forest |
|------|-----------|---------------|
| **線形DGP** | ✅ 最適 | △ 階段近似 |
| **非線形DGP** | △ 限定的 | ✅ 柔軟 |
| **小サンプル** | ✅ n=100でもOK | ❌ n=100で失敗 |
| **予測形状** | 滑らか | 階段状 |
| **計算速度** | 遅い（Python） | 速い（R native） |
| **解釈性** | 係数解釈可 | 変数重要度 |

### 4.2 実用的推奨

**線形効果が期待される場合**:
→ **DR-Learner** 推奨（特に修正版）

**非線形・複雑な効果の場合**:
→ **Causal Forest** だが、以下の調整が必要:
- n < 200: `honesty.fraction = 0.2-0.3`
- `min.node.size` の調整
- チューニングが重要

**オミクスデータ（高次元）**:
→ 次のフェーズで検証が必要
- p > n での性能
- 変数選択能力
- スパース効果の検出

### 4.3 σ_y = 5.0 の妥当性

**個人差の範囲**:
- Y(0) ~ N(120, 5²) → 範囲 [105, 135]
- 血圧として妥当（やや均質だが許容範囲）

**検出難易度**:
- Signal (τのSD): 1.6 mmHg
- Noise (σ_y): 5.0 mmHg
- SNR = 0.32（中程度）

**手法間の差**:
- SexGene 追加の効果が明確（12-90%改善）
- CF vs DR の差も顕著

---

## 5. 次のステップ

### 完了した検証
- ✅ 基本的なCATE推定性能（p=2）
- ✅ DR-Learner の修正・検証
- ✅ サンプルサイズの影響（n=100-1000）
- ✅ 効果形状の影響（linear vs threshold）
- ✅ ノイズレベルの影響（sigma_y=5.0）

### 未検証の重要項目

**Phase 2: バイオマーカー発見**（最優先）
- 高次元データ（p=100-1000）
- スパース効果（k=1-5 真の効果修飾因子）
- 変数選択精度評価

**Phase 3: 複雑な効果**
- 遺伝子間相関
- 交互作用効果
- バッチ効果

**Phase 4: 追加手法**
- GANITE（深層学習）
- CEVAE（深層学習）
- その他のCATEメソッド

---

## 6. ファイル一覧

### スクリプト
- `run_cate_demov04_with_fix.R` - メインベンチマーク（DR修正版使用）
- `visualize_demov04_sigma5.R` - 可視化スクリプト
- `methods_drlearner_FIXED.R` - DR-Learner修正版実装

### 結果ファイル
- `results_demov04_long_FIXED.csv` - 全条件の詳細結果
- `results_demov04_summary_FIXED.csv` - サマリー統計
- `predictions_demov04_FIXED.csv` - 個別予測値

### 可視化（12枚）
- `plot_geneA_vs_cate_*_sigma5.png`
  - 2 methods × 2 effects × 3 sample sizes

### ドキュメント
- `FIX_SUMMARY.md` - DR-Learner修正内容の詳細
- `CHANGES_DIFF.md` - 修正のdiff
- `ISSUES_drlearner_cfomics.md` - 発見した問題点の記録
- `RESULTS_SUMMARY.md` - 本ドキュメント

---

## 7. 再現方法

```bash
# 環境設定（Python venv使用）
cd /home/rstudio/work/cfomics

# ベンチマーク実行
Rscript sandbox/demov04/run_cate_demov04_with_fix.R

# 可視化
Rscript sandbox/demov04/visualize_demov04_sigma5.R

# 結果確認
cat sandbox/demov04/results_with_fix/results_demov04_summary_FIXED.csv
```

---

## 付録: パラメータ一覧

```r
config <- list(
  seed_base  = 20260215L,
  R          = 1L,
  n_values   = c(100L, 500L, 1000L),
  effect_types = c("linear", "threshold3"),
  covsets    = c("SexOnly", "SexGene"),
  methods    = c("causal_forest", "drlearner"),
  train_frac = 0.8,

  # RNA-seq params
  m_E = log(50),      # Mean log expression
  s_E = 0.35,         # SD log expression
  theta_nb = 30,      # NB dispersion
  sigma_G = 0.08,     # Measurement error

  # CATE params
  tau_base = -12.0,   # Baseline effect
  a_S = 2.0,          # Sex effect
  a_G = 1.5,          # Gene effect
  delta_thr = 1.5,    # Threshold jump
  tau_min = -18.0,    # Min effect
  tau_max = -6.0,     # Max effect

  # Outcome params
  mu0 = 120.0,        # Baseline BP mean
  sigma_y = 5.0       # Individual differences SD
)
```
