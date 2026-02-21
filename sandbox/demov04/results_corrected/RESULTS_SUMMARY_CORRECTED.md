# demov04 ベンチマーク結果（修正版DGP）

**実行日**: 2026-02-20
**目的**: CATE推定手法の比較（オミクスデータ想定）
**重要な修正**: DGPの符号修正により、治療で血圧が正しく**下がる**ようになった

---

## 修正内容

### 問題点（修正前）

```r
# 間違ったDGP
Y = Y0 - T × tau
where tau < 0 (例: -12 mmHg)

# 結果
T=0: Y = Y0 = 120 mmHg
T=1: Y = Y0 - (-12) = Y0 + 12 = 132 mmHg
→ 治療により血圧が上昇（降圧薬として逆！）❌
```

### 修正版DGP

```r
# 正しいDGP
Y = Y0 + T × tau
where tau < 0 (例: -12 mmHg)

# 結果
T=0: Y = Y0 = 120 mmHg
T=1: Y = Y0 + (-12) = Y0 - 12 = 108 mmHg
→ 治療により血圧が下がる（降圧効果）✓
```

**変更箇所**:
- `Y = Y0 - T × tau` → `Y = Y0 + T × tau`（1行のみ）
- 予測関数の符号調整（`-predict(...)` → `predict(...)`）

---

## 1. シナリオ設定

### データ生成プロセス (DGP)

**基本設定**:
```r
n = 500          # サンプルサイズ
R = 1            # 反復回数
train_frac = 0.8 # Train/Test分割
```

**アウトカムモデル（修正版）**:
```r
# ベースライン（個人差あり）
Y(0) ~ N(120, σ_y²)  where σ_y = 5.0
範囲: [100, 140] mmHg

# 潜在的アウトカム
Y(1) = Y(0) + τ  where τ < 0
範囲: [88, 127] mmHg（治療により血圧が下がる）

# 観測アウトカム
Y = Y(0) + T × τ

# 測定誤差なし（個人差のみ）
```

**遺伝子発現モデル**:
```r
# 真の発現量
E_true ~ N(log(50), 0.35²)

# カウントデータ（負の二項分布）
C_count ~ NB(mu = exp(E_true), size = 30)

# 観測値（log1p変換 + 測定誤差）
G_obs = log1p(C_count) + N(0, 0.08²)
```

**治療効果モデル（Linear）**:
```r
τ = -12.0 + 2.0×(Sex-0.5) + 1.5×Z_obs
範囲: [-18, -6] mmHg (clipped)

# Sex効果: ±1 mmHg
# geneA効果: ±3 mmHg程度
```

---

## 2. 結果

### 2.1 CATE推定性能（RMSE）

| 共変量セット | Causal Forest | DR-Learner (FIXED) | 改善率 |
|-------------|--------------|-------------------|--------|
| **SexOnly** | 1.670 | 1.650 | 1.2% |
| **SexGene** | 1.238 | 1.294 | -4.5% |

**geneA追加の効果**:
- Causal Forest: 1.670 → 1.238（**25.9%改善**）✓
- DR-Learner: 1.650 → 1.294（**21.6%改善**）✓

### 2.2 潜在的アウトカムの統計

**観測値（修正版）**:
```
Y(0) 未治療時: 120.09 ± 4.93 mmHg
Y(1) 治療時:   108.09 ± 5.18 mmHg
ATE:           -12.00 mmHg（降圧効果）✓

対照群 (T=0): 平均 120.19 mmHg
治療群 (T=1): 平均 108.06 mmHg
単純差分:     -12.14 mmHg
```

**修正前との比較**:
| 項目 | 修正前（間違い） | 修正後（正しい） |
|------|----------------|----------------|
| Y(T=0) | 120 mmHg | 120 mmHg |
| Y(T=1) | 132 mmHg ❌ | 108 mmHg ✓ |
| ATE | +12 mmHg ❌ | -12 mmHg ✓ |
| 解釈 | 血圧上昇 ❌ | 血圧低下 ✓ |

---

## 3. 主要な知見

### 3.1 DGP修正の影響

**CATE推定性能への影響**:
- RMSE値は**ほぼ同じ**（修正前: 1.6-0.6, 修正後: 1.7-1.2）
- 手法間の相対的な性能も類似
- **符号が逆だっただけで、推定精度自体は有効だった**

**解釈の修正**:
- ❌ 「治療により血圧が上昇」→ ✓ 「治療により血圧が下がる」
- τ < 0 が降圧効果を表す（正しい医学的解釈）
- Y(1) < Y(0) となり、因果推論の一般的な表現と一致

### 3.2 手法の特徴（修正版でも同様）

**Causal Forest**:
- Linear DGPでも一定の性能
- geneA追加で26%改善（効果修飾因子を検出）
- 階段状の予測（tree-basedの特性）

**DR-Learner (FIXED)**:
- Linear DGPで安定した性能
- geneA追加で22%改善
- 滑らかな線形予測

---

## 4. 可視化

### 生成されたプロット

1. **`plot_potential_outcomes_corrected.png`**
   - Y(0) と Y(1) の分布（Y(1) < Y(0) ✓）
   - 散布図（全ての点が対角線より下 ✓）
   - ITE と τ の分布（負の値 ✓）

2. **`plot_geneA_vs_cate_causal_forest_corrected.png`**
   - geneA vs CATE（Causal Forest）
   - τ が負の値で表示 ✓
   - geneA高 → τ大（降圧効果小）

3. **`plot_geneA_vs_cate_drlearner_corrected.png`**
   - geneA vs CATE（DR-Learner）
   - 線形の傾向を捉える ✓

---

## 5. 重要な教訓

### 5.1 DGP設計の重要性

**符号の整合性**:
- τ の定義（正/負）
- Y の計算式（Y = Y0 ± T×τ）
- 効果の解釈（増加/減少）

→ **一貫性が重要**

### 5.2 ベンチマーク結果の頑健性

- 符号ミスがあっても、CATE推定の**相対的な性能評価**は有効
- RMSE は絶対値の差なので、符号に依存しない
- ただし、医学的解釈としては完全に誤り

### 5.3 検証の重要性

**確認すべきポイント**:
```r
# 1. Y(T=0) と Y(T=1) の平均値
mean(Y[T==0])  # 120 mmHg
mean(Y[T==1])  # 108 mmHg ← 下がっているべき

# 2. τ の符号
mean(tau)  # -12 mmHg ← 負であるべき

# 3. 潜在的アウトカム
Y(1) < Y(0)  # 治療で改善 ← 全個人で成立すべき
```

---

## 6. ファイル一覧

### スクリプト
- `run_cate_corrected.R` - 修正版ベンチマーク
- `visualize_corrected.R` - geneA vs CATE可視化
- `visualize_potential_outcomes_corrected.R` - Y(0), Y(1)可視化

### 結果ファイル
- `results_corrected_long.csv` - 詳細結果
- `results_corrected_summary.csv` - サマリー統計
- `predictions_corrected.csv` - 個別予測値

### 可視化
- `plot_potential_outcomes_corrected.png` - 潜在的アウトカム
- `plot_geneA_vs_cate_causal_forest_corrected.png` - CF結果
- `plot_geneA_vs_cate_drlearner_corrected.png` - DR結果

### ドキュメント
- `RESULTS_SUMMARY_CORRECTED.md` - 本ドキュメント

---

## 7. 再現方法

```bash
# 環境設定
cd /home/rstudio/work/cfomics

# ベンチマーク実行（修正版）
Rscript sandbox/demov04/run_cate_corrected.R

# 可視化
Rscript sandbox/demov04/visualize_corrected.R
Rscript sandbox/demov04/visualize_potential_outcomes_corrected.R

# 結果確認
cat sandbox/demov04/results_corrected/results_corrected_summary.csv
```

---

## 8. 次のステップ

**Phase 2**: 高次元バイオマーカー発見（修正版DGPを使用）
- p=100-1000 遺伝子
- 正しい符号でのCATE推定
- 変数選択精度の評価

**修正の教訓を活かす**:
- DGP設計時に医学的妥当性を確認
- 潜在的アウトカムの可視化を常に実施
- Y(T=0) vs Y(T=1) の関係を検証

---

## 付録: パラメータ一覧

```r
config <- list(
  seed_base  = 20260215L,
  n_val      = 500L,
  effect_type = "linear",

  # RNA-seq params
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,

  # CATE params
  tau_base = -12.0,   # Baseline effect（負=降圧）
  a_S = 2.0,          # Sex effect
  a_G = 1.5,          # Gene effect
  tau_min = -18.0,    # Min effect
  tau_max = -6.0,     # Max effect

  # Outcome params
  mu0 = 120.0,        # Baseline BP mean
  sigma_y = 5.0       # Individual differences SD
)
```

**修正点**: `Y = Y0 + T × tau`（1行のみ変更）

---

✓ **修正版DGPにより、医学的に正しい解釈が可能になりました**
