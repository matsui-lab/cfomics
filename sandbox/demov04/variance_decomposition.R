#!/usr/bin/env Rscript
# variance_decomposition.R
# 性別とgeneAがτのばらつきをどれくらい説明しているか

config <- list(
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0,
  mu0 = 120.0, sigma_y = 5.0
)

seed_cond <- 20260425
set.seed(seed_cond)

n <- 500
S <- rbinom(n, 1, 0.5)
Sc <- S - 0.5

# geneA生成
E_true <- rnorm(n, config$m_E, config$s_E)
mu_nb <- exp(E_true)
C_count <- rnbinom(n, size = config$theta_nb, mu = mu_nb)
G_raw <- log1p(C_count)
delta <- rnorm(n, 0, config$sigma_G)
G_obs <- G_raw + delta
Z_obs <- (G_obs - mean(G_obs)) / sd(G_obs)

# τの計算
tau_sex <- config$a_S * Sc              # 性別による寄与
tau_gene <- config$a_G * Z_obs          # geneAによる寄与
tau_indiv <- tau_sex + tau_gene         # 合計の個人差
tau_raw <- config$tau_base + tau_indiv  # ベース効果追加
tau <- pmin(config$tau_max, pmax(config$tau_min, tau_raw))  # クリッピング

cat("=== τのばらつきの分散分解 ===\n\n")

# 分散の計算
var_sex <- var(tau_sex)
var_gene <- var(tau_gene)
var_indiv <- var(tau_indiv)
var_tau <- var(tau)

cat("各成分の分散:\n")
cat(sprintf("  性別の寄与 (a_S × Sc):        Var = %.4f\n", var_sex))
cat(sprintf("  geneAの寄与 (a_G × Z_obs):    Var = %.4f\n", var_gene))
cat(sprintf("  個人差合計 (tau_indiv):       Var = %.4f\n", var_indiv))
cat(sprintf("  最終的なτ (クリッピング後):  Var = %.4f\n", var_tau))

# 標準偏差（mmHg単位）
sd_sex <- sd(tau_sex)
sd_gene <- sd(tau_gene)
sd_indiv <- sd(tau_indiv)
sd_tau <- sd(tau)

cat("\n各成分の標準偏差 (mmHg):\n")
cat(sprintf("  性別の寄与:    SD = %.3f mmHg\n", sd_sex))
cat(sprintf("  geneAの寄与:   SD = %.3f mmHg\n", sd_gene))
cat(sprintf("  個人差合計:    SD = %.3f mmHg\n", sd_indiv))
cat(sprintf("  最終的なτ:     SD = %.3f mmHg\n", sd_tau))

# 分散の説明率
prop_sex <- var_sex / var_indiv * 100
prop_gene <- var_gene / var_indiv * 100

cat("\nクリッピング前の個人差に対する説明率:\n")
cat(sprintf("  性別:   %.2f%%\n", prop_sex))
cat(sprintf("  geneA:  %.2f%%\n", prop_gene))
cat(sprintf("  合計:   %.2f%% (独立なので和が100%%)\n", prop_sex + prop_gene))

# 理論値との比較
cat("\n理論的な分散（独立性を仮定）:\n")
cat(sprintf("  Var(Sc) = %.4f (二値変数 p=0.5)\n", var(Sc)))
cat(sprintf("  Var(Z_obs) = %.4f (標準化変数なので ≈ 1)\n", var(Z_obs)))
cat(sprintf("  理論 Var(a_S × Sc) = a_S² × Var(Sc) = %.2f² × %.4f = %.4f\n",
            config$a_S, var(Sc), config$a_S^2 * var(Sc)))
cat(sprintf("  理論 Var(a_G × Z_obs) = a_G² × Var(Z_obs) = %.2f² × %.4f = %.4f\n",
            config$a_G, var(Z_obs), config$a_G^2 * var(Z_obs)))

# 範囲の比較
cat("\nτの範囲への寄与:\n")
cat(sprintf("  性別による変動幅:   [%.2f, %.2f] = %.2f mmHg\n",
            min(tau_sex), max(tau_sex), max(tau_sex) - min(tau_sex)))
cat(sprintf("  geneAによる変動幅:  [%.2f, %.2f] = %.2f mmHg\n",
            min(tau_gene), max(tau_gene), max(tau_gene) - min(tau_gene)))
cat(sprintf("  個人差合計の範囲:   [%.2f, %.2f] = %.2f mmHg\n",
            min(tau_indiv), max(tau_indiv), max(tau_indiv) - min(tau_indiv)))
cat(sprintf("  最終的なτの範囲:    [%.2f, %.2f] = %.2f mmHg\n",
            min(tau), max(tau), max(tau) - min(tau)))

# 性別ごとの統計
cat("\n性別ごとのτの統計:\n")
tau_female <- tau[S == 0]
tau_male <- tau[S == 1]
cat(sprintf("  Female (n=%d): 平均 = %.2f, SD = %.2f, 範囲 = [%.2f, %.2f]\n",
            length(tau_female), mean(tau_female), sd(tau_female),
            min(tau_female), max(tau_female)))
cat(sprintf("  Male (n=%d):   平均 = %.2f, SD = %.2f, 範囲 = [%.2f, %.2f]\n",
            length(tau_male), mean(tau_male), sd(tau_male),
            min(tau_male), max(tau_male)))
cat(sprintf("  性別間の平均の差: %.2f mmHg\n", mean(tau_male) - mean(tau_female)))

cat("\n=== まとめ ===\n")
cat(sprintf("性別は τ のばらつきの %.1f%% を説明\n", prop_sex))
cat(sprintf("geneA は τ のばらつきの %.1f%% を説明\n", prop_gene))
cat(sprintf("\n血圧の幅でいうと:\n"))
cat(sprintf("  性別による違い: 約 %.1f mmHg の変動\n", sd_sex * 2))
cat(sprintf("  geneAによる違い: 約 %.1f mmHg の変動\n", sd_gene * 2))
cat(sprintf("  合計の個人差:   約 %.1f mmHg の変動\n", sd_tau * 2))
