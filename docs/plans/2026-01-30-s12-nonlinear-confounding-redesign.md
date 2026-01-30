# S12 DGP 再設計: 非線形交絡

## 問題

現在の S12 `dgp_nonlinear_outcome` は結果モデルのみ非線形で、PS は線形。
処置割り付けが X の線形関数のため、OLS (`Y ~ T + X`) でも T の係数にバイアスが入らず、gformula を壊せない。

50rep ベンチマーク結果:

| シナリオ | gformula | hdml | hdps | tmle |
|---|---|---|---|---|
| S12 moderate | 0.123 | 0.131 | 0.136 | 0.131 |
| S12 severe | 0.106 | 0.111 | 0.126 | 0.112 |

手法間の差が小さく、gformula が最良になっている。

## 原因分析

定数処置効果 `tau = 2.0` かつ PS が X の線形関数の場合:
- OLS は `Y = X*beta + tau*T + e` をフィット
- 非線形項 (X1², sin(X1) 等) は残差に入るが、T とほぼ無相関（PS が線形なので）
- T の係数推定にバイアスが入らない

**必要条件**: 非線形項が PS にも影響する（＝非線形交絡）ことで、OLS の脱落変数バイアスが発生する。

## 設計

### 修正内容

PS と Y が同じ非線形関数 `h(X)` に依存するように変更。

**Moderate:**
```r
h <- X[, 1]^2 + X[, 2] * X[, 3] + sin(X[, 1])
logit_ps <- 0.8 * h
Y <- h + tau * T + noise
```

**Severe:**
```r
h <- X[, 1]^3 / 3 + exp(X[, 2] / 2) + X[, 3] * X[, 4] * (X[, 1] > 0) +
     sin(2 * X[, 1]) * cos(X[, 2])
h_centered <- h - mean(h)
logit_ps <- 0.5 * h_centered
Y <- h + tau * T + noise
```

PS クリッピング: `[0.05, 0.95]`。Severe では `h` を中心化して処置群バランスを維持。

### 期待される手法間差異

| 手法 | 影響 | 理由 |
|---|---|---|
| gformula | 高バイアス | OLS が h(X) を線形近似 → 脱落変数バイアス |
| hdml | 中バイアス | AIPW だが両方のモデルが誤特定 → DR 保護なし |
| hdps | 中バイアス | IPW のみ。PS 線形モデルが h(X) を捉えられない |
| tmle | 中バイアス | ターゲティングで部分補正あり |

### 変更ファイル

- `packages/cfomics/R/benchmark_dgp.R`: `dgp_nonlinear_outcome` の PS 生成を修正
- `packages/cfomics/tests/testthat/test-benchmark-dgp.R`: テスト更新（T の分布が変わるため）
- ベンチマーク再実行 + レポート再生成

### テスト更新

既存テスト「moderate vs severe share same T」は PS が変わるため成立しなくなる。
修正: 「moderate vs severe produce different Y」のみ検証に変更。
