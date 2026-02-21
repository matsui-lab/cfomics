#!/usr/bin/env Rscript
# analyze_failure.R - n=500 linear SexGene 失敗原因の調査

config <- list(
  seed_base  = 20260215L,
  train_frac = 0.8,
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0, mu0 = 120.0, sigma_y = 2.0
)

generate_data_v04 <- function(n, effect_type, seed, cfg) {
  set.seed(seed)
  S  <- rbinom(n, 1, 0.5)
  Sc <- S - 0.5
  E_true  <- rnorm(n, cfg$m_E, cfg$s_E)
  mu_nb   <- exp(E_true)
  C_count <- rnbinom(n, size = cfg$theta_nb, mu = mu_nb)
  G_raw   <- log1p(C_count)
  delta   <- rnorm(n, 0, cfg$sigma_G)
  G_obs   <- G_raw + delta
  T_vec   <- rbinom(n, 1, 0.5)

  Z_obs <- (G_obs - mean(G_obs)) / sd(G_obs)
  tau_indiv <- cfg$a_S * Sc + cfg$a_G * Z_obs
  tau_raw <- cfg$tau_base + tau_indiv
  tau <- pmin(cfg$tau_max, pmax(cfg$tau_min, tau_raw))

  # Y を生成
  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y <- Y0 - T_vec * tau

  list(S = S, G_obs = G_obs, T_vec = T_vec, tau = tau, Y = Y)
}

stratified_split <- function(T_vec, train_frac, seed) {
  set.seed(seed)
  idx0 <- which(T_vec == 0); idx1 <- which(T_vec == 1)
  n0_tr <- round(length(idx0) * train_frac)
  n1_tr <- round(length(idx1) * train_frac)
  tr_idx <- sort(c(sample(idx0)[seq_len(n0_tr)], sample(idx1)[seq_len(n1_tr)]))
  te_idx <- sort(setdiff(seq_along(T_vec), tr_idx))
  list(train = tr_idx, test = te_idx)
}

# 失敗条件を再生成
seed_cond <- 20260215
dat <- generate_data_v04(500, "linear", seed_cond, config)
spl <- stratified_split(dat$T_vec, config$train_frac, seed_cond)

tr <- spl$train
te <- spl$test

cat("=== n=500, linear, SexGene 失敗条件の分析 ===\n\n")

cat("1. サンプルサイズ:\n")
cat("   Train:", length(tr), "/ Test:", length(te), "\n\n")

cat("2. Treatment 分布 (train):\n")
print(table(T = dat$T_vec[tr]))
cat("\n")

cat("3. 共変量の統計 (train):\n")
cat("   Sex (S):    ", sprintf("%.3f ± %.3f", mean(dat$S[tr]), sd(dat$S[tr])), "\n")
cat("   geneA:      ", sprintf("%.3f ± %.3f", mean(dat$G_obs[tr]), sd(dat$G_obs[tr])), "\n")
cat("   geneA range:", sprintf("[%.2f, %.2f]", min(dat$G_obs[tr]), max(dat$G_obs[tr])), "\n\n")

cat("4. アウトカム Y の統計 (train):\n")
cat("   Overall:    ", sprintf("%.2f ± %.2f", mean(dat$Y[tr]), sd(dat$Y[tr])), "\n")
cat("   T=0 群:     ", sprintf("%.2f ± %.2f (n=%d)",
                              mean(dat$Y[tr][dat$T_vec[tr]==0]),
                              sd(dat$Y[tr][dat$T_vec[tr]==0]),
                              sum(dat$T_vec[tr]==0)), "\n")
cat("   T=1 群:     ", sprintf("%.2f ± %.2f (n=%d)",
                              mean(dat$Y[tr][dat$T_vec[tr]==1]),
                              sd(dat$Y[tr][dat$T_vec[tr]==1]),
                              sum(dat$T_vec[tr]==1)), "\n")
cat("   Y range:    ", sprintf("[%.2f, %.2f]", min(dat$Y[tr]), max(dat$Y[tr])), "\n\n")

cat("5. 真の CATE (τ) の分布 (train):\n")
cat("   Mean:       ", sprintf("%.3f", mean(dat$tau[tr])), "\n")
cat("   SD:         ", sprintf("%.3f", sd(dat$tau[tr])), "\n")
cat("   Range:      ", sprintf("[%.2f, %.2f]", min(dat$tau[tr]), max(dat$tau[tr])), "\n\n")

# 傾向スコアの推定（簡易版）
cat("6. 傾向スコアの分布（ロジスティック回帰）:\n")

# SexOnly モデル
ps_model_sexonly <- glm(T_vec ~ S, data = data.frame(T_vec = dat$T_vec[tr], S = dat$S[tr]),
                        family = binomial())
ps_sexonly <- predict(ps_model_sexonly, type = "response")
cat("   SexOnly モデル:\n")
cat("     e(X) range: ", sprintf("[%.4f, %.4f]", min(ps_sexonly), max(ps_sexonly)), "\n")
cat("     極端値 (<0.05 or >0.95): ", sum(ps_sexonly < 0.05 | ps_sexonly > 0.95),
    "/", length(ps_sexonly), "\n")

# SexGene モデル
ps_model_sexgene <- glm(T_vec ~ S + G_obs,
                        data = data.frame(T_vec = dat$T_vec[tr], S = dat$S[tr], G_obs = dat$G_obs[tr]),
                        family = binomial())
ps_sexgene <- predict(ps_model_sexgene, type = "response")
cat("   SexGene モデル:\n")
cat("     e(X) range: ", sprintf("[%.4f, %.4f]", min(ps_sexgene), max(ps_sexgene)), "\n")
cat("     極端値 (<0.05 or >0.95): ", sum(ps_sexgene < 0.05 | ps_sexgene > 0.95),
    "/", length(ps_sexgene), "\n\n")

cat("7. 共変量バランス (standardized difference):\n")
t0_idx <- tr[dat$T_vec[tr] == 0]
t1_idx <- tr[dat$T_vec[tr] == 1]

smd_S <- abs(mean(dat$S[t1_idx]) - mean(dat$S[t0_idx])) /
         sqrt((var(dat$S[t1_idx]) + var(dat$S[t0_idx])) / 2)
smd_G <- abs(mean(dat$G_obs[t1_idx]) - mean(dat$G_obs[t0_idx])) /
         sqrt((var(dat$G_obs[t1_idx]) + var(dat$G_obs[t0_idx])) / 2)

cat("   Sex:        ", sprintf("%.3f", smd_S),
    ifelse(smd_S > 0.1, " ⚠️ 不均衡", " ✓ バランス良好"), "\n")
cat("   geneA:      ", sprintf("%.3f", smd_G),
    ifelse(smd_G > 0.1, " ⚠️ 不均衡", " ✓ バランス良好"), "\n")
