#!/usr/bin/env Rscript
# explain_internal_splits.R
# 各手法がtrainデータをどう使っているか説明

cat("=== 各手法のtrainデータの使い方 ===\n\n")

cat("【1】外側の分割（共通）\n")
cat("────────────────────────────────\n")
cat("全データ (n) → Train (80%) + Test (20%)\n")
cat("  n=500:  Train=400, Test=100\n")
cat("  n=1000: Train=800, Test=200\n")
cat("\n→ この Train データで各手法を学習\n")
cat("→ Test データで性能評価（RMSE計算）\n\n")

cat("【2】Causal Forest (GRF) の内部処理\n")
cat("────────────────────────────────\n")
cat("GRFは「Honest Trees」を使用:\n\n")

cat("各決定木ごとに Train データをさらに分割:\n")
cat("  - Building subsample (50%): 木構造の構築に使用\n")
cat("  - Estimation subsample (50%): 葉ノードでの効果推定に使用\n\n")

cat("例: n=1000の場合\n")
cat("  Train = 800\n")
cat("  → 各木で Building用 ≈ 400, Estimation用 ≈ 400\n")
cat("  → ランダムフォレスト（デフォルト2000本）なので、\n")
cat("     各サンプルは複数の木で使用される\n\n")

cat("honesty = TRUE (デフォルト)の効果:\n")
cat("  - 過学習を防ぐ\n")
cat("  - 推定値のバイアスを減らす\n")
cat("  - 信頼区間の構築が可能\n\n")

cat("【3】DR-Learner の内部処理\n")
cat("────────────────────────────────\n")
cat("3つのステップで異なるモデルを学習、各モデルでCVを使用:\n\n")

cat("ステップ1: Propensity Score モデル (LogisticRegressionCV)\n")
cat("  - Train データ全体を使用\n")
cat("  - 5-fold CV で正則化パラメータを選択\n")
cat("  - n_train=800なら、各foldで640 train / 160 validation\n\n")

cat("ステップ2: Outcome Regression モデル (RidgeCV)\n")
cat("  - Y(T=0)とY(T=1)を別々に推定\n")
cat("  - 各モデルで 5-fold CV\n")
cat("  - 各群 n≈400なら、各foldで320 train / 80 validation\n\n")

cat("ステップ3: Final CATE モデル (LassoCV)\n")
cat("  - 疑似アウトカム（DR変換後）を使用\n")
cat("  - 5-fold CV で正則化パラメータを選択\n")
cat("  - n_train=800なら、各foldで640 train / 160 validation\n\n")

cat("CV設定:\n")
cat("  - n_train >= 200: cv = 5 (5-fold)\n")
cat("  - n_train < 200:  cv = 3 (3-fold)\n\n")

cat("【4】実際の使用データ量の比較\n")
cat("────────────────────────────────\n")
cat("\nn=1000 (Train=800) の場合:\n\n")

cat("Causal Forest:\n")
cat("  - 各木で Building用 ≈400 + Estimation用 ≈400\n")
cat("  - 2000本の木で反復使用\n")
cat("  - 実効的にはTrain全体（800）を活用\n\n")

cat("DR-Learner:\n")
cat("  - Propensity model: 800全体、5-fold CV\n")
cat("  - Outcome models: 各群≈400、5-fold CV\n")
cat("  - Final model: 800全体、5-fold CV\n")
cat("  - CVにより各モデルで 1/5 (20%) は常にvalidationに使用\n")
cat("  - 最終的には最適なハイパーパラメータで全データを再学習\n\n")

cat("【5】なぜDR-Learnerの方がよくfitするのか？\n")
cat("────────────────────────────────\n")
cat("\n1. **線形DGPに対する適合性**\n")
cat("   - DR-Learner: 線形モデル（Ridge, Lasso）使用\n")
cat("   - GRF: 非線形な木ベース\n")
cat("   → 真のDGPが線形なら線形モデルが有利\n\n")

cat("2. **CVによる適切な正則化**\n")
cat("   - 各モデルで最適な正則化パラメータを選択\n")
cat("   - 過学習と過小学習のバランスが良い\n\n")

cat("3. **サンプルサイズの効果**\n")
cat("   - n=500: DR RMSE=0.610\n")
cat("   - n=1000: DR RMSE=0.168 (72%改善)\n")
cat("   → 線形モデルはサンプルサイズ増加の恩恵を受けやすい\n\n")

cat("4. **Doubly Robust性質**\n")
cat("   - Propensity modelかOutcome modelのどちらかが正しければ一致推定量\n")
cat("   - RCTなのでPropensity scoreは既知（0.5）だが、推定しても頑健\n\n")

cat("【まとめ】\n")
cat("────────────────────────────────\n")
cat("\n✓ 両手法ともTrainデータを効率的に活用\n")
cat("✓ Causal Forest: Honesty splitting で過学習防止\n")
cat("✓ DR-Learner: 5-fold CV で各モデルの正則化を最適化\n")
cat("✓ 線形DGPでは線形モデルベースのDR-Learnerが優位\n")
cat("✓ n=1000で DR-Learner は準完璧なフィット（R²=0.992）を達成\n")
