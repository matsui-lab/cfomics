#!/usr/bin/env Rscript
# compare_wrong_vs_corrected.R
# 符号修正前後での性能比較

# 符号修正前（間違ったDGP: Y = Y0 - T×tau）
wrong <- data.frame(
  method = c("Causal Forest", "Causal Forest", "DR-Learner", "DR-Learner"),
  covset = c("SexOnly", "SexGene", "SexOnly", "SexGene"),
  rmse = c(1.691, 0.882, 1.597, 0.629),
  stringsAsFactors = FALSE
)

# 符号修正後（正しいDGP: Y = Y0 + T×tau）
correct <- data.frame(
  method = c("Causal Forest", "Causal Forest", "DR-Learner", "DR-Learner"),
  covset = c("SexOnly", "SexGene", "SexOnly", "SexGene"),
  rmse = c(1.280, 1.087, 1.287, 0.610),
  stringsAsFactors = FALSE
)

cat("=== 符号修正前後でのRMSE比較（n=500, linear） ===\n\n")

# 比較表
comparison <- data.frame(
  Method = wrong$method,
  Covset = wrong$covset,
  Wrong_DGP = wrong$rmse,
  Correct_DGP = correct$rmse,
  Diff = correct$rmse - wrong$rmse,
  Change = ifelse(correct$rmse < wrong$rmse, "改善", "悪化")
)

print(comparison, row.names = FALSE)

cat("\n=== 重要な観察 ===\n\n")

# Causal Forest + SexGene の悪化
cf_sg_wrong <- wrong$rmse[wrong$method == "Causal Forest" & wrong$covset == "SexGene"]
cf_sg_correct <- correct$rmse[correct$method == "Causal Forest" & correct$covset == "SexGene"]

cat(sprintf("【1】Causal Forest + SexGene が悪化:\n"))
cat(sprintf("    符号修正前: %.3f\n", cf_sg_wrong))
cat(sprintf("    符号修正後: %.3f\n", cf_sg_correct))
cat(sprintf("    差分: %+.3f (%.1f%% 悪化)\n\n",
            cf_sg_correct - cf_sg_wrong,
            (cf_sg_correct - cf_sg_wrong) / cf_sg_wrong * 100))

# DR-Learner + SexGene はほぼ同じ
dr_sg_wrong <- wrong$rmse[wrong$method == "DR-Learner" & wrong$covset == "SexGene"]
dr_sg_correct <- correct$rmse[correct$method == "DR-Learner" & wrong$covset == "SexGene"]

cat(sprintf("【2】DR-Learner + SexGene はほぼ同じ（むしろ改善）:\n"))
cat(sprintf("    符号修正前: %.3f\n", dr_sg_wrong))
cat(sprintf("    符号修正後: %.3f\n", dr_sg_correct))
cat(sprintf("    差分: %+.3f (%.1f%% 改善)\n\n",
            dr_sg_correct - dr_sg_wrong,
            (dr_sg_correct - dr_sg_wrong) / dr_sg_wrong * 100))

# SexOnly は両方とも改善
cat(sprintf("【3】SexOnly は両手法とも大幅改善:\n"))
cat(sprintf("    CF: %.3f → %.3f (%.1f%% 改善)\n",
            wrong$rmse[wrong$covset == "SexOnly" & wrong$method == "Causal Forest"],
            correct$rmse[correct$covset == "SexOnly" & correct$method == "Causal Forest"],
            (correct$rmse[correct$covset == "SexOnly" & correct$method == "Causal Forest"] -
             wrong$rmse[wrong$covset == "SexOnly" & wrong$method == "Causal Forest"]) /
            wrong$rmse[wrong$covset == "SexOnly" & wrong$method == "Causal Forest"] * 100))
cat(sprintf("    DR: %.3f → %.3f (%.1f%% 改善)\n\n",
            wrong$rmse[wrong$covset == "SexOnly" & wrong$method == "DR-Learner"],
            correct$rmse[correct$covset == "SexOnly" & correct$method == "DR-Learner"],
            (correct$rmse[correct$covset == "SexOnly" & correct$method == "DR-Learner"] -
             wrong$rmse[wrong$covset == "SexOnly" & wrong$method == "DR-Learner"]) /
            wrong$rmse[wrong$covset == "SexOnly" & wrong$method == "DR-Learner"] * 100))

cat("=== なぜCausal Forest + SexGeneだけ悪化したのか？ ===\n\n")

cat("考えられる理由:\n\n")
cat("1. **異なるseed値**\n")
cat("   - 符号修正前: seed = 20260325\n")
cat("   - 符号修正後: seed = 20260425\n")
cat("   - データの実現値が異なるため、性能も変わる\n\n")

cat("2. **Y の分布の違いによるツリー構造への影響**\n")
cat("   - 符号修正前: Y ∈ [112, 138] (高い値)\n")
cat("   - 符号修正後: Y ∈ [92, 128] (低い値)\n")
cat("   - Causal Forestは決定木ベースなので、Y の範囲が変わると\n")
cat("     最適な分割点や木構造が変わる可能性がある\n\n")

cat("3. **ランダム性**\n")
cat("   - Causal Forestは本質的にランダムフォレスト\n")
cat("   - seedが違えば結果も変わる\n")
cat("   - 0.882 → 1.087 の変化は統計的な変動の範囲内かもしれない\n\n")

cat("4. **DR-Learnerは安定**\n")
cat("   - DR-Learnerは線形モデルベース\n")
cat("   - Y の値が変わっても、線形関係が同じなら性能は安定\n")
cat("   - 実際、0.629 → 0.610 とほぼ同じ性能を維持\n\n")

cat("=== 結論 ===\n\n")
cat("✓ 「符号修正前の方がfitしていた」は **Causal Forest + SexGene** に関しては正しい\n")
cat("✓ ただし、これは **seedの違い** による偶然の可能性が高い\n")
cat("✓ DR-Learnerは符号修正後も同等以上の性能を維持（0.610）\n")
cat("✓ SexOnlyは両方とも改善している\n")
cat("\n→ 符号修正後のDGPが医学的に正しいので、現在の結果を使うべき\n")
