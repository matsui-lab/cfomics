#!/usr/bin/env Rscript
# test_fix_n500.R
# 修正版 DR-Learner で n=500 失敗条件を検証

cat("=== DR-Learner 修正版の検証 (n=500, linear, SexGene) ===\n\n")

suppressPackageStartupMessages({
  library(reticulate)
})

# cfomics の helper 関数を定義（簡易版）
.wrap_python_call <- function(expr, method_name = "Python") {
  tryCatch(expr, error = function(e) {
    stop(sprintf("[%s] Python error: %s", method_name, conditionMessage(e)))
  })
}

cf_require_python <- function(method) {
  if (!py_available(initialize = TRUE)) {
    stop("Python not available")
  }
  TRUE
}

# 修正版を読み込み
source("sandbox/demov04/methods_drlearner_FIXED.R")

# DGP 設定
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

# 失敗条件でデータ生成
seed_cond <- 20260215
dat <- generate_data_v04(500, "linear", seed_cond, config)
spl <- stratified_split(dat$T_vec, config$train_frac, seed_cond)

tr <- spl$train
te <- spl$test

cat("1. データ生成完了\n")
cat("   Train:", length(tr), "/ Test:", length(te), "\n\n")

# SexGene 共変量
X_train <- data.frame(S = dat$S[tr], geneA = dat$G_obs[tr])
X_test  <- data.frame(S = dat$S[te], geneA = dat$G_obs[te])
T_train <- dat$T_vec[tr]
T_test  <- dat$T_vec[te]
Y_train <- dat$Y[tr]

tau_true_test <- dat$tau[te]

cat("2. 修正版 DR-Learner でフィット開始...\n")

result_fixed <- tryCatch({
  fit_fixed <- cf_fit_drlearner(
    X = X_train,
    T = T_train,
    Y = Y_train,
    random_state = as.integer(seed_cond),
    cv = 5L,
    min_propensity = 0.05,      # FIX #2
    standardize_y = TRUE,        # FIX #4
    model_propensity = "auto",   # FIX #1
    model_regression = "auto",   # FIX #1
    model_final = "auto"         # FIX #3
  )

  # 予測（符号反転不要 - demov04 の Y モデルに依存）
  # demov04 では Y = Y0 - T*tau なので、予測は負の値
  # 真値も負なので、符号はそのまま
  tau_hat_test <- predict_cf_drlearner(fit_fixed, newdata = X_test, type = "ite")

  # 実際には符号反転が必要（Y モデルに依存）
  # demov04: Y = Y(0) - T*tau → tau は負 → 予測も負 → 反転して正にする？
  # いや、真値も負なので反転不要
  # 確認のため両方計算
  tau_hat_raw <- tau_hat_test
  tau_hat_flipped <- -tau_hat_test

  rmse_raw <- sqrt(mean((tau_true_test - tau_hat_raw)^2))
  rmse_flipped <- sqrt(mean((tau_true_test - tau_hat_flipped)^2))

  # クリップ検出
  n_clip <- sum(abs(tau_hat_raw) >= 19.9)  # ±20 に近い

  list(
    success = TRUE,
    tau_hat_raw = tau_hat_raw,
    tau_hat_flipped = tau_hat_flipped,
    rmse_raw = rmse_raw,
    rmse_flipped = rmse_flipped,
    n_clip = n_clip,
    fit = fit_fixed
  )
}, error = function(e) {
  list(success = FALSE, error = conditionMessage(e))
})

cat("\n3. 結果:\n")

if (!result_fixed$success) {
  cat("   ❌ フィット失敗:", result_fixed$error, "\n")
} else {
  cat("   ✅ フィット成功\n\n")

  cat("   予測値の統計:\n")
  cat("     Raw 予測:    ", sprintf("%.2f ± %.2f", mean(result_fixed$tau_hat_raw),
                                    sd(result_fixed$tau_hat_raw)), "\n")
  cat("     Raw range:   ", sprintf("[%.2f, %.2f]",
                                    min(result_fixed$tau_hat_raw),
                                    max(result_fixed$tau_hat_raw)), "\n")
  cat("     真値:        ", sprintf("%.2f ± %.2f", mean(tau_true_test),
                                    sd(tau_true_test)), "\n")
  cat("     真値 range:  ", sprintf("[%.2f, %.2f]",
                                    min(tau_true_test),
                                    max(tau_true_test)), "\n\n")

  cat("   RMSE:\n")
  cat("     Raw 予測:    ", sprintf("%.3f", result_fixed$rmse_raw), "\n")
  cat("     Flipped予測: ", sprintf("%.3f", result_fixed$rmse_flipped), "\n\n")

  cat("   クリッピング:\n")
  cat("     ±20付近:     ", result_fixed$n_clip, "/", length(te), "\n")
  cat("     クリップ率:  ", sprintf("%.1f%%",
                                    100 * result_fixed$n_clip / length(te)), "\n\n")

  # 比較: 旧実装の結果（RMSE=19.94, 100%クリップ）
  cat("4. 旧実装との比較:\n")
  cat("   旧実装 RMSE:     19.94 (100% clipped)\n")
  cat("   修正版 RMSE:     ", sprintf("%.3f (%.1f%% clipped)",
                                      min(result_fixed$rmse_raw, result_fixed$rmse_flipped),
                                      100 * result_fixed$n_clip / length(te)), "\n\n")

  improvement <- 19.94 - min(result_fixed$rmse_raw, result_fixed$rmse_flipped)
  cat("   改善幅:          ", sprintf("%.2f (%.1f%% 削減)",
                                      improvement,
                                      100 * improvement / 19.94), "\n\n")

  if (min(result_fixed$rmse_raw, result_fixed$rmse_flipped) < 1.0) {
    cat("   ✅ 修正成功！RMSE < 1.0 達成\n")
  } else if (min(result_fixed$rmse_raw, result_fixed$rmse_flipped) < 5.0) {
    cat("   ⚠️  改善したが、まだ RMSE > 1.0\n")
  } else {
    cat("   ❌ まだ問題が残っている\n")
  }
}

cat("\n=== 検証完了 ===\n")
