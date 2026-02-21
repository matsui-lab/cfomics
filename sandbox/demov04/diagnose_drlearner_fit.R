#!/usr/bin/env Rscript
# diagnose_drlearner_fit.R
# なぜ線形DGPで完璧にフィットしないのか調査

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

config <- list(
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0,
  mu0 = 120.0, sigma_y = 5.0
)

# データ生成
seed_cond <- 20260425
set.seed(seed_cond)

n <- 500
S <- rbinom(n, 1, 0.5)
Sc <- S - 0.5
E_true <- rnorm(n, config$m_E, config$s_E)
mu_nb <- exp(E_true)
C_count <- rnbinom(n, size = config$theta_nb, mu = mu_nb)
G_raw <- log1p(C_count)
delta <- rnorm(n, 0, config$sigma_G)
G_obs <- G_raw + delta
T_vec <- rbinom(n, 1, 0.5)
Z_obs <- (G_obs - mean(G_obs)) / sd(G_obs)

tau_indiv <- config$a_S * Sc + config$a_G * Z_obs
tau_raw <- config$tau_base + tau_indiv
tau <- pmin(config$tau_max, pmax(config$tau_min, tau_raw))

Y0 <- config$mu0 + rnorm(n, 0, config$sigma_y)
Y <- Y0 + T_vec * tau

# Train/test分割
set.seed(seed_cond)
idx0 <- which(T_vec == 0); idx1 <- which(T_vec == 1)
n0_tr <- round(length(idx0) * 0.8)
n1_tr <- round(length(idx1) * 0.8)
tr_idx <- sort(c(sample(idx0, n0_tr), sample(idx1, n1_tr)))
te_idx <- setdiff(1:n, tr_idx)

cat("=== 線形DGPでの完璧なフィット可能性の検証 ===\n\n")
cat(sprintf("Train: %d, Test: %d\n\n", length(tr_idx), length(te_idx)))

# ============================================================
# ベースライン1: 真のパラメータを使った予測（理論的最良）
# ============================================================
cat("【1】真のパラメータでの予測（理論的上限）\n")
tau_true_test <- tau[te_idx]

# 真のパラメータ: tau = -12.0 + 2.0*Sc + 1.5*Z_obs
Sc_test <- S[te_idx] - 0.5
Z_test <- Z_obs[te_idx]
tau_oracle <- config$tau_base + config$a_S * Sc_test + config$a_G * Z_test
tau_oracle_clip <- pmin(config$tau_max, pmax(config$tau_min, tau_oracle))

rmse_oracle <- sqrt(mean((tau_true_test - tau_oracle_clip)^2))
cat(sprintf("  RMSE = %.6f (ほぼ0になるはず)\n", rmse_oracle))
cat(sprintf("  → クリッピングの影響のみ\n\n"))

# ============================================================
# ベースライン2: 単純なOLS回帰（正則化なし）
# ============================================================
cat("【2】単純なOLS線形回帰（正則化なし）\n")

# Y(0)とY(1)を別々に推定
df_train <- data.frame(
  Y = Y[tr_idx],
  S = S[tr_idx],
  geneA = G_obs[tr_idx],
  T = T_vec[tr_idx]
)

# 対照群のアウトカムモデル
df_train0 <- df_train[df_train$T == 0, ]
lm0 <- lm(Y ~ S + geneA, data = df_train0)

# 治療群のアウトカムモデル
df_train1 <- df_train[df_train$T == 1, ]
lm1 <- lm(Y ~ S + geneA, data = df_train1)

X_test <- data.frame(S = S[te_idx], geneA = G_obs[te_idx])
mu0_hat <- predict(lm0, newdata = X_test)
mu1_hat <- predict(lm1, newdata = X_test)
tau_ols <- mu1_hat - mu0_hat

rmse_ols <- sqrt(mean((tau_true_test - tau_ols)^2))
cat(sprintf("  RMSE = %.3f\n", rmse_ols))
cat(sprintf("  R² = %.3f\n", 1 - (rmse_ols^2 / var(tau_true_test))))

# 係数確認
cat("\n  対照群モデル係数:\n")
print(coef(lm0))
cat("\n  治療群モデル係数:\n")
print(coef(lm1))
cat("\n  推定された効果（係数の差）:\n")
cat(sprintf("    切片の差 (ベース効果): %.3f (真値: %.1f)\n",
            coef(lm1)[1] - coef(lm0)[1], config$tau_base))
cat(sprintf("    性別の差: %.3f (真値: %.1f)\n",
            coef(lm1)[2] - coef(lm0)[2], config$a_S))
cat(sprintf("    geneAの差: %.3f (真値: %.1f)\n\n",
            coef(lm1)[3] - coef(lm0)[3], config$a_G))

# ============================================================
# ベースライン3: 直接τを予測（チート的方法）
# ============================================================
cat("【3】τを直接OLS回帰（非現実的だが参考）\n")

# Trainデータでτを直接観測できる仮定
df_train_cheat <- data.frame(
  tau = tau[tr_idx],
  S = S[tr_idx],
  geneA = G_obs[tr_idx]
)

lm_cheat <- lm(tau ~ S + geneA, data = df_train_cheat)
tau_cheat <- predict(lm_cheat, newdata = X_test)

rmse_cheat <- sqrt(mean((tau_true_test - tau_cheat)^2))
cat(sprintf("  RMSE = %.3f\n", rmse_cheat))
cat(sprintf("  R² = %.3f\n", 1 - (rmse_cheat^2 / var(tau_true_test))))
cat("\n  推定された係数:\n")
print(coef(lm_cheat))
cat(sprintf("\n    真のパラメータ: 切片=%.1f, S=%.1f, geneA=%.1f\n\n",
            config$tau_base, config$a_S, config$a_G))

# ============================================================
# 実際のDR-Learnerの結果
# ============================================================
cat("【4】DR-Learner (RidgeCV + LassoCV, 正則化あり)\n")
df_pred <- read.csv("sandbox/demov04/results_corrected/predictions_corrected.csv")
sub <- subset(df_pred, method == "drlearner" & covset == "SexGene")
rmse_dr <- sqrt(mean((sub$tau_true - sub$tau_hat)^2, na.rm = TRUE))
cat(sprintf("  RMSE = %.3f\n", rmse_dr))
cat(sprintf("  R² = %.3f\n\n", 1 - (rmse_dr^2 / var(sub$tau_true))))

# ============================================================
# まとめ
# ============================================================
cat("=== まとめ ===\n\n")
cat("| 手法 | RMSE | 説明 |\n")
cat("|------|------|------|\n")
cat(sprintf("| 真のパラメータ | %.6f | クリッピングのみの影響 |\n", rmse_oracle))
cat(sprintf("| τ直接回帰(OLS) | %.3f | 非現実的だが最良の線形近似 |\n", rmse_cheat))
cat(sprintf("| Y(0),Y(1)別々のOLS | %.3f | 正則化なしのDR的アプローチ |\n", rmse_ols))
cat(sprintf("| DR-Learner (正則化) | %.3f | 実際の実装（Ridge+Lasso） |\n\n", rmse_dr))

cat("【結論】\n")
cat("1. 完璧なフィットは理論的には可能（RMSE ≈ 0）\n")
cat("2. 単純なOLSでも RMSE ≈ ", sprintf("%.3f", rmse_ols), "\n")
cat("3. DR-Learnerは正則化により ", sprintf("%.3f", rmse_dr), " に留まる\n")
cat("4. 正則化の影響: RMSE差 ≈ ", sprintf("%.3f", rmse_dr - rmse_ols), " mmHg\n\n")

if (rmse_ols < 1.0 && rmse_dr > rmse_ols) {
  cat("→ 正則化（Ridge/Lasso）がバイアスを導入している可能性が高い\n")
  cat("  有限サンプルでの過学習防止のため、やや保守的な推定になる\n")
} else if (rmse_ols > 1.0) {
  cat("→ サンプルサイズ不足または分割の影響が大きい\n")
}

# 可視化
out_dir <- "sandbox/demov04/results_corrected"

df_compare <- data.frame(
  tau_true = rep(tau_true_test, 3),
  tau_hat = c(tau_ols, tau_cheat, sub$tau_hat),
  Method = rep(c("OLS (Y0,Y1別)", "OLS (τ直接)", "DR-Learner"),
               each = length(tau_true_test))
)

p <- ggplot(df_compare, aes(x = tau_true, y = tau_hat)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.4, size = 2) +
  facet_wrap(~Method) +
  labs(title = "線形DGPでの予測精度比較",
       subtitle = sprintf("OLS: %.3f | τ直接: %.3f | DR-Learner: %.3f",
                         rmse_ols, rmse_cheat, rmse_dr),
       x = "真値 τ", y = "予測値 τ_hat") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "linear_fit_comparison.png"),
       p, width = 12, height = 4, dpi = 150)

cat("\n✓ 図を保存: linear_fit_comparison.png\n")
