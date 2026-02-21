#!/usr/bin/env Rscript
# DGPの違いを確認

config <- list(
  seed_base = 20260215L,
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0, mu0 = 120.0, sigma_y = 5.0
)

seed <- 20260425
set.seed(seed)

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

# 元のDGP
Y_old <- Y0 - T_vec * tau

# 修正版DGP
Y_new <- Y0 + T_vec * tau

cat("=== DGP比較（同じseed） ===\n\n")

cat("τ（治療効果）:\n")
cat(sprintf("  平均: %.3f\n", mean(tau)))
cat(sprintf("  SD:   %.3f\n", sd(tau)))
cat(sprintf("  範囲: [%.2f, %.2f]\n", min(tau), max(tau)))

cat("\nY0（ベースライン）:\n")
cat(sprintf("  平均: %.3f\n", mean(Y0)))
cat(sprintf("  SD:   %.3f\n", sd(Y0)))

cat("\n元のDGP (Y = Y0 - T×τ):\n")
cat(sprintf("  Y(T=0) 平均: %.3f\n", mean(Y_old[T_vec==0])))
cat(sprintf("  Y(T=1) 平均: %.3f\n", mean(Y_old[T_vec==1])))
cat(sprintf("  Y 全体 SD:   %.3f\n", sd(Y_old)))
cat(sprintf("  差分: %.3f\n", mean(Y_old[T_vec==1]) - mean(Y_old[T_vec==0])))

cat("\n修正版DGP (Y = Y0 + T×τ):\n")
cat(sprintf("  Y(T=0) 平均: %.3f\n", mean(Y_new[T_vec==0])))
cat(sprintf("  Y(T=1) 平均: %.3f\n", mean(Y_new[T_vec==1])))
cat(sprintf("  Y 全体 SD:   %.3f\n", sd(Y_new)))
cat(sprintf("  差分: %.3f\n", mean(Y_new[T_vec==1]) - mean(Y_new[T_vec==0])))

cat("\n重要な観察:\n")
cat("  - τ は両方で同じ（同じseed、同じ計算式）\n")
cat("  - Y0 も同じ\n")
cat("  - でも Y の値は完全に異なる\n")
cat("  - Y の分散も異なる可能性がある\n")
cat("\n  → アウトカムモデルが学習するYが違うので、\n")
cat("    CATE推定の性能も変わる可能性がある\n")
