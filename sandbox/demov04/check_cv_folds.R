#!/usr/bin/env Rscript
# CV fold 内での傾向スコア・アウトカムモデルを検証

# Manual CV fold creation

config <- list(
  seed_base = 20260215L, train_frac = 0.8,
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
  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y <- Y0 - T_vec * tau
  list(S = S, G_obs = G_obs, T_vec = T_vec, tau = tau, Y = Y, Y0 = Y0)
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

cat("=== CV fold 内での極端な傾向スコアの発生 ===\n\n")

seed_cond <- 20260215
dat <- generate_data_v04(500, "linear", seed_cond, config)
spl <- stratified_split(dat$T_vec, config$train_frac, seed_cond)
tr <- spl$train

df_train <- data.frame(
  Y = dat$Y[tr],
  T = dat$T_vec[tr],
  S = dat$S[tr],
  G = dat$G_obs[tr]
)

# 手動でCV fold を作成（stratified by T）
set.seed(seed_cond)
n_folds <- 5

# Stratified CV folds
idx_t0 <- which(df_train$T == 0)
idx_t1 <- which(df_train$T == 1)

folds_t0 <- split(sample(idx_t0), rep(1:n_folds, length.out = length(idx_t0)))
folds_t1 <- split(sample(idx_t1), rep(1:n_folds, length.out = length(idx_t1)))

folds <- lapply(1:n_folds, function(i) sort(c(folds_t0[[i]], folds_t1[[i]])))

cat("Train サイズ:", nrow(df_train), "\n")
cat("CV folds:", n_folds, "\n\n")

all_extreme_ps <- list()

for (fold_idx in seq_along(folds)) {
  test_fold <- folds[[fold_idx]]
  train_fold <- setdiff(seq_len(nrow(df_train)), test_fold)

  df_cv_train <- df_train[train_fold, ]
  df_cv_test <- df_train[test_fold, ]

  cat(sprintf("--- Fold %d ---\n", fold_idx))
  cat("  Train:", nrow(df_cv_train), "/ Test:", nrow(df_cv_test), "\n")

  # 傾向スコアモデル
  ps_model <- glm(T ~ S + G, data = df_cv_train, family = binomial())
  ps_train <- predict(ps_model, newdata = df_cv_train, type = "response")
  ps_test <- predict(ps_model, newdata = df_cv_test, type = "response")

  cat("  傾向スコア (train):\n")
  cat("    Range:", sprintf("[%.4f, %.4f]", min(ps_train), max(ps_train)), "\n")
  extreme_train <- sum(ps_train < 0.05 | ps_train > 0.95)
  cat("    極端値(<0.05 or >0.95):", extreme_train, "/", length(ps_train), "\n")

  cat("  傾向スコア (test, CV評価用):\n")
  cat("    Range:", sprintf("[%.4f, %.4f]", min(ps_test), max(ps_test)), "\n")
  extreme_test <- sum(ps_test < 0.05 | ps_test > 0.95)
  cat("    極端値(<0.05 or >0.95):", extreme_test, "/", length(ps_test), "\n")

  if (extreme_test > 0) {
    idx_extreme <- which(ps_test < 0.05 | ps_test > 0.95)
    all_extreme_ps[[fold_idx]] <- list(
      fold = fold_idx,
      ps = ps_test[idx_extreme],
      T = df_cv_test$T[idx_extreme]
    )
  }

  # アウトカムモデル
  df_cv_t0 <- subset(df_cv_train, T == 0)
  df_cv_t1 <- subset(df_cv_train, T == 1)

  model_y0 <- lm(Y ~ S + G, data = df_cv_t0)
  model_y1 <- lm(Y ~ S + G, data = df_cv_t1)

  # test での予測と残差
  mu0_test <- predict(model_y0, newdata = df_cv_test)
  mu1_test <- predict(model_y1, newdata = df_cv_test)

  residual_t0_test <- df_cv_test$Y[df_cv_test$T == 0] -
                      mu0_test[df_cv_test$T == 0]
  residual_t1_test <- df_cv_test$Y[df_cv_test$T == 1] -
                      mu1_test[df_cv_test$T == 1]

  cat("  残差 SD (test):\n")
  if (length(residual_t0_test) > 0) {
    cat("    T=0:", sprintf("%.2f", sd(residual_t0_test)), "\n")
  }
  if (length(residual_t1_test) > 0) {
    cat("    T=1:", sprintf("%.2f", sd(residual_t1_test)), "\n")
  }

  cat("\n")
}

if (length(all_extreme_ps) > 0) {
  cat("\n=== 極端な傾向スコアが発生した fold の詳細 ===\n")
  for (item in all_extreme_ps) {
    cat(sprintf("\nFold %d:\n", item$fold))
    for (i in seq_along(item$ps)) {
      cat(sprintf("  ê(X) = %.4f, T = %d\n", item$ps[i], item$T[i]))
    }
  }
} else {
  cat("\n✓ 全てのfoldで極端な傾向スコアは発生しませんでした\n")
}

cat("\n結論:\n")
cat("CV fold内でも極端な傾向スコアは発生していない可能性が高い。\n")
cat("では、なぜn=500で失敗したのか？\n\n")
cat("可能性:\n")
cat("1. EconML内部のCV実装が異なる（stratified ではない？）\n")
cat("2. DR疑似アウトカムの計算時に数値誤差が蓄積\n")
cat("3. Y標準化なし + 中程度の残差 → 累積的な不安定性\n")
cat("4. model_final (OLS) が外れ値に過剰反応\n")
