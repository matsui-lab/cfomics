# cfomicsSim DGP Suite 設計 (v1.0)

cfomicsBench / cfomicsSim の中核となる「包括的シミュレーションモデル（DGP：Data Generating Process）スイート」の設計仕様書です。

単なる統計的おもちゃではなく、オミックス解析で実際に起きる破綻要因（バッチ、LOD、ゼロ過剰、MNAR、階層、時変交絡、非線形・交互作用、潜在交絡、マルチオミックス）を「モジュールとして組み合わせ可能」にし、**真値（ATE/CATE/ITE/反実仮想）を常に保持**することを目指します。

---

## 0. 設計原則（ベンチマークとしての要件）

1. **真値が定義できる**：Y^(a)（潜在アウトカム）をシミュレータ内部で生成し、ATE/CATE/ITE、y0, y1 を保存する。

2. **omics-native**：単一アウトカムだけでなく、**高次元特徴（p=10^3–10^5）**と、**MultiAssayExperiment/SE/SCE 互換出力**を前提にする。

3. **現実的破綻要因を"独立モジュール"で表現**：バッチ・欠測・LOD・ゼロ過剰・階層・縦断・潜在交絡などを、オン/オフまたは強度パラメータで制御する。

4. **シナリオ名＝再現可能な研究単位**：`scenario_id` が同じなら同じ機構（分布・因果構造・ノイズ）を必ず再現できる。

5. **推定対象（estimand）を先に固定**：ベンチで評価する主ターゲットは
   - ATE / CATE / ITE（臨床アウトカムY）
   - 反実仮想予測（y0, y1）

   必要なら追加タスクとして「特徴量ごとの因果効果（omics outcome）」を別枠で定義する。

---

## 1. 変数定義（最小の統一因果モデル）

サンプル i=1..n、オミックス層 ℓ=1..L（例：RNA/Protein/Metabolite）、潜在因子次元 K。

| 変数 | 説明 |
|------|------|
| X_i ∈ ℝ^{p_x} | 観測共変量（連続＋カテゴリ混在可） |
| U_i ∈ ℝ^{p_u} | 潜在交絡（観測されない） |
| S_i | クラスター：施設/ロット等 |
| B_i | バッチ |
| A_i | 介入（Treatment）：基本は二値、拡張で連続/多値 |
| Z_i ∈ ℝ^K | 生物学的潜在状態（層間で共有） |
| M̃_{iℓj} | 層ごとの真のオミックス量（測定前）：特徴 j=1..p_ℓ |
| M_{iℓj} | 観測オミックス（バッチ、LOD、欠測等を反映） |
| Y_i | アウトカム（連続/二値/生存等。まずは連続＋ロジスティックを標準） |

---

## 2. 構造的因果モデル（SCM：Structural Causal Model）

### 2.1 共変量生成（相関・非線形を含む）

- **連続共変量**：X_i ~ N(0, Σ)（AR(1) などで相関制御）
- **カテゴリ**：X^{cat} を多項分布で生成し、ダミー化して使用可能
- **非線形特徴**：h(X)（例：X₁², sin(X₂), I(X₃>0)）を内部的に用意し、治療・アウトカムの非線形性に利用

### 2.2 潜在交絡

- U_i ~ N(0, I) または混合正規（heavy tail 用）
- "潜在交絡シナリオ"では U を **A と Y と Z（ひいては M）に同時に入れる**

### 2.3 介入割付（Treatment assignment）

基本：二値介入

```
A_i ~ Bernoulli(σ(α₀ + α_X^⊤ X_i + α_U^⊤ U_i + α_S(S_i) + α_{BX}(B_i, X_i)))
```

- **positivity/overlap 制御**：α のスケール、または "ほぼ決定的割付" で重なりを崩す（評価上とても重要）
- **非線形割付**：α_X^⊤ h(X_i) を入れる（線形仮定の破綻）
- **バッチと治療の相関**：α_{BX} で B → A を意図的に作る（現実の悪条件）

### 2.4 生物学的潜在状態（因子）

```
Z_i = β₀ + β_X X_i + β_A A_i + β_U U_i + η_i,  η_i ~ N(0, I_K)
```

- ここが **層間相関（マルチオミックス）**と **生物学的共有変動**の源泉

### 2.5 真のオミックス量（測定前の log-scale を標準）

層 ℓ、特徴 j：

```
M̃_{iℓj} = μ_{ℓj} + λ_{ℓj}^⊤ Z_i + γ_{ℓj} A_i + φ_{ℓj}^⊤ X_i + ψ_{ℓj}^⊤ U_i + δ^{(batch)}_{ℓ,B_i,j} + ε_{iℓj}
```

| パラメータ | 説明 |
|------------|------|
| λ_{ℓj} | **因子負荷量**（疎な構造にして生物学モジュールっぽくする） |
| γ_{ℓj} | 治療の真の効果（因果シグナル）。多くは 0、少数が非ゼロ（スパース） |
| φ, ψ | 観測/潜在交絡による影響（"交絡で差が出る特徴"を作る） |
| δ^{(batch)} | 特徴ごとのバッチ効果（後述の生成機構） |
| ε | 層依存ノイズ（分散も特徴依存にできる） |

### 2.6 観測オミックス（分布族＋観測過程）

層ごとに「分布族」を切り替える：

#### RNA-seq（count）

```
C_{iℓj} ~ NegBin(mean = exp(M̃_{iℓj}) · s_i, disp = θ_{ℓj})
```

その後、logCPM などへ変換して M としても良い（ただし"前処理が因果仮定を壊す"評価用に、raw と normalized の両方を出せる設計が望ましい）

#### Proteomics/Metabolomics（intensity）

```
I_{iℓj} = exp(M̃_{iℓj}) · ε,  log I = M̃ + log ε,  log ε ~ N(0, σ²)
```

#### Compositional（microbiome など）（オプション）

Dirichlet-multinomial/log-ratio 系（後述シナリオで追加）

---

## 3. 現実的破綻要因モジュール（オン/オフ・強度制御）

### 3.1 バッチ効果（Batch）

**目的**：単なるノイズではなく、**治療や共変量と相関して"疑似効果"を作る**。

生成例：

- バッチ数：n_{batch}
- バッチ効果：
  ```
  δ^{(batch)}_{ℓ,b,j} ~ N(0, σ²_{ℓ,batch})  （特徴ごと）
  ```
- さらに現実性を上げるなら、バッチ効果に低ランク成分を混ぜる：
  ```
  δ_{ℓ,b,·} = L_ℓ q_b + r_{ℓ,b},  q_b ~ N(0, I)
  ```
- "バッチが治療群に偏る" を作る：P(B|A) を偏らせる、または A の割付式に B を入れる

### 3.2 LOD（検出限界）と左側打ち切り

- 特徴ごとの LOD：LOD_{ℓj}
- 観測：
  - left-censor：M_{obs} = max(M, LOD) か
  - "未検出＝0 or NA" として返す（メソッドによって扱いが変わる点が評価対象）

### 3.3 ゼロ過剰（Zero inflation）

ゼロ化確率を強度化：

```
ZI_{iℓj} ~ Bernoulli(π_{iℓj}),  π_{iℓj} = σ(κ₀ + κ₁ M̃_{iℓj} + κ₂ B_i)
```

ZI=1 のとき M=0（または NA）に置換

### 3.4 欠測（MCAR / MAR / MNAR）

| タイプ | 定義 |
|--------|------|
| MCAR | 一定確率で欠測 |
| MAR | P(R=1\|X, B) |
| MNAR | P(R=1\|M̃, X, B)（最も重要。観測研究オミックスの難所） |

```
R_{iℓj} ~ Bernoulli(σ(ρ₀ + ρ₁ M̃_{iℓj} + ρ₂ A_i + ρ₃ B_i))
```

### 3.5 潜在交絡（Unmeasured confounding）

- U を A, Y, Z/M に同時に入れる（上の式の α_U, β_U, ψ をオン）
- 強度パラメータ：||α_U||, ||β_U||, ||ψ|| を調整

### 3.6 階層・クラスター（site/family/repeated）

施設 S_i のランダム効果：

```
u^{(site)}_{S_i} ~ N(0, σ²_{site})
```

を A と Y と M̃ に入れられる

反復測定（縦断）を導入する場合は後述の 3.7

### 3.7 縦断・時変交絡（Longitudinal / time-varying）

時点 t=0..T：

- 時変共変量 X_{it}、介入 A_{it}、潜在状態 Z_{it}、アウトカム Y_{it}、オミックス M_{it}
- 典型的な時変交絡（過去の治療が将来の共変量を変える）：
  ```
  X_{i,t+1} = g(X_{it}, A_{it}, U_i) + ξ
  A_{it} ~ Bernoulli(σ(... + α_X^⊤ X_{it} + α_A A_{i,t-1}))
  ```
- 評価ターゲット：最終時点 Y_{iT} の ATE、または治療レジーム間比較

---

## 4. アウトカムモデル（推定対象の真値を保証する設計）

アウトカム Y は「真の効果が分かる」形で設計する。

### 4.1 連続アウトカム（標準）

```
Y_i = f(X_i) + τ(X_i) · A_i + ω^⊤ U_i + ν^⊤ Z_i + ε_i
```

| 成分 | 説明 |
|------|------|
| f(X) | 非線形可（ベースライン難度） |
| τ(X) | 効果不均一性（CATE/ITE の真値）。例：τ(X) = τ₀ + τ₁ X₁ + τ₂ I(X₂>0) |
| U, Z | 潜在交絡・媒介を入れることで難度調整 |

### 4.2 二値アウトカム（オプション）

```
P(Y_i = 1) = σ(f(X_i) + τ(X_i) A_i + ...)
```

（ATE の定義はリスク差・オッズ比などで揺れるので、まずは連続が無難。二値は別タスク化が安全）

---

## 5. シナリオ設計（"包括的"を破綻させない最小カタログ）

「包括的」は"軸設計"で担保し、シナリオ数は制御します。推奨は **12–16 シナリオ**。

### 5.1 軸（factor）定義

#### 因果難度：{低, 中, 高}

- 交絡強度（観測/潜在）
- 非線形性
- 効果不均一性
- positivity 崩壊

#### 測定難度：{低, 中, 高}

- バッチ（治療と相関）
- LOD / ゼロ過剰
- MAR/MNAR 欠測

#### 構造難度：{単層, マルチオミックス, 縦断, クラスター}

### 5.2 推奨シナリオ（v1.0 最小セット）

| シナリオID | 名称 | 説明 |
|------------|------|------|
| S00 | linear_clean | 線形・交絡弱・ノイズ低 |
| S01 | hetero_nonlinear | τ(X)不均一＋非線形 |
| S02 | strong_confounding | 観測交絡強（X→A,Y強） |
| S03 | positivity_violation | overlap 崩壊（重み暴れ） |
| S04 | hidden_confounder | U あり（未測定交絡） |
| S05 | batch_correlated | バッチが治療と相関（疑似効果） |
| S06 | lod_zero_inflated | LOD＋ゼロ過剰（MS/メタボ寄り） |
| S07 | mnar_missing | MNAR 欠測（強） |
| S08 | highdim_p>>n | 高次元共変量（omics を confounder 側に置くモード） |
| S09 | multiomics_shared_factors | 3層（共有Z＋層特異因子） |
| S10 | longitudinal_tvc | 時変交絡（縦断） |
| S11 | clustered_sites | 施設クラスター＋ランダム効果 |
| S12 | compositional_microbiome | （オプション）Compositional データ |

---

## 6. 真値（Ground truth）の保持方法

ベンチマークとして最重要です。

### 6.1 ITE / CATE / ATE

個体 ITE（真値）：

```
ITE_i = E[Y_i^{(1)} - Y_i^{(0)} | X_i]
```

ここで期待は U, η, ε に対して取る（観測不能でもシミュレータは知っている）。

- **線形なら解析的に計算**できる設計にしておく（推奨）
- 非線形の場合は **固定Xで Monte Carlo**（m=200–1000）で近似し、truth として保存
- ATE：(1/n) Σ_i ITE_i（サンプルATE）も併記

### 6.2 反実仮想アウトカム

- y0_i = E[Y_i^{(0)} | X_i]、y1_i を保存
- 観測アウトカムは Y_i = Y_i^{(A_i)}

### 6.3 参考真値（診断評価用）

- **真の propensity**：e(X, U) = P(A=1 | X, U, ...)（推定側には出さないが、overlap 破綻度の真値評価に使える）
- **真の「効果あり特徴」集合**：{j : γ_{ℓj} ≠ 0}（omics outcome タスク用）

---

## 7. 出力仕様（Bioconductor 互換＋ベンチ向けメタ）

### 7.1 推奨出力

| データタイプ | 出力形式 |
|--------------|----------|
| 単層 | `SummarizedExperiment` |
| マルチオミックス | `MultiAssayExperiment` |
| 縦断 | `MultiAssayExperiment`（assay に時点を持たせる）または `SummarizedExperiment` + `colData$time` |

### 7.2 フィールド設計

| フィールド | 内容 |
|------------|------|
| `colData` | A、観測共変量 X、batch/site/time など（推定に渡す） |
| `rowData` | 特徴タイプ、LOD、真の効果フラグ（bench 用メタとして `metadata` 側に隠してもよい） |
| `metadata$truth` | ATE/CATE/ITE、y0, y1、propensity（隠し）、生成パラメータ、seed |
| `metadata$mechanism` | オンにしたモジュール（batch/LOD/MNAR等）と強度 |

---

## 8. シナリオ記述（コンフィグ化の推奨）

実装では「シナリオ＝YAML/JSON（またはR list）で宣言」し、DGP を組み立てるのが保守性最強です。

### 構成ブロック例

```yaml
scenario:
  id: S01_hetero_nonlinear

  core:
    n: 1000
    px: 10
    K: 3
    layers:
      - name: RNA
        p: 5000
        distribution: negbin
      - name: Protein
        p: 1000
        distribution: lognormal

  assignment:
    model: logistic
    nonlinear: true
    positivity_level: high  # high/medium/low

  outcome:
    baseline: nonlinear  # f(X)
    heterogeneity: true  # τ(X)
    noise_sd: 1.0

  omics:
    factor_sparsity: 0.3
    effect_sparsity: 0.05
    cross_layer_sharing: 0.7

  nuisance:
    batch:
      enabled: true
      n_batches: 5
      treatment_correlated: false
    lod:
      enabled: false
    zero_inflation:
      enabled: false
    missingness:
      type: none  # none/mcar/mar/mnar
    cluster:
      enabled: false
    longitudinal:
      enabled: false
```

---

## 9. 実装上の要点（HPC/大規模 p に耐えるため）

- p が大きい層は **チャンク生成**（1e4ずつ）し、必要ならオンディスク（HDF5/Arrow）へ
- 乱数は **シナリオ seed ＋ rep seed ＋ chunk seed** で分岐（並列でも再現性維持）
- truth の計算（MC）は n 全部で重いので、
  - 解析解が出る設計を優先
  - 非線形はサブサンプルで truth を補助的に出す（ただしベンチの評価指標と整合させる）

---

## 10. 「包括的」と呼べる理由（論文化の要点）

1. **オミックス特有ノイズ**（batch/LOD/zero/MNAR）を明示的にモデル化
2. **因果の難しさ**（positivity、潜在交絡、時変交絡、階層）を軸で制御
3. **真値**（ATE/CATE/ITE/反実仮想）を常に保持し、診断性能も評価可能

---

## 次のステップ

1. cfomicsSim の「シナリオ仕様（schema）」の実装
2. v1.0 の 12 シナリオのデフォルトパラメータ表の作成
3. cfomicsBench 側の評価指標（推定誤差＋overlap 診断検出＋計算量）との統合
4. Paper 2 の Methods セクション草案への展開
