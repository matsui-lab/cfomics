#!/usr/bin/env Rscript
# run_cate_demov04_with_fix.R
# DR-Learner のみ修正版を使用

# ============================================================
# 0. ENVIRONMENT SETUP (元のスクリプトと同じ)
# ============================================================
suppressPackageStartupMessages({
  venv_python <- file.path(getwd(), "sandbox", ".venv", "bin", "python3")
  Sys.setenv(RETICULATE_PYTHON = venv_python)
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0")
  devtools::load_all("packages/cfomics", quiet = TRUE)
})

# 修正版 DR-Learner を読み込み（上書き）
cat("Loading FIXED DR-Learner implementation...\n")
source("sandbox/demov04/methods_drlearner_FIXED.R")

# 元のベンチマークスクリプトを読み込み
cat("Loading benchmark script...\n")

# ============================================================
# CONFIG
# ============================================================
config <- list(
  seed_base  = 20260215L,
  R          = 1L,
  n_values   = c(100L, 500L, 1000L),
  effect_types = c("linear", "threshold3"),
  covsets    = c("SexOnly", "SexGene"),
  methods    = c("causal_forest", "drlearner"),
  train_frac = 0.8,
  out_dir    = "sandbox/demov04/results_with_fix",

  # DGP params (案A - 現実的設定, sigma_y=5.0)
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5, delta_thr = 1.5,
  tau_min = -18.0, tau_max = -6.0,
  mu0 = 120.0, sigma_y = 5.0
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# DGP
# ============================================================
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
  Z_true  <- (E_true - cfg$m_E) / cfg$s_E
  T_vec   <- rbinom(n, 1, 0.5)

  Z_obs <- (G_obs - mean(G_obs)) / sd(G_obs)

  if (effect_type == "linear") {
    tau_indiv <- cfg$a_S * Sc + cfg$a_G * Z_obs
  } else if (effect_type == "threshold3") {
    q1 <- quantile(Z_obs, 0.33)
    q2 <- quantile(Z_obs, 0.66)
    g_Z <- ifelse(Z_obs <= q1, -cfg$delta_thr,
                  ifelse(Z_obs <= q2, 0, cfg$delta_thr))
    tau_indiv <- cfg$a_S * Sc + g_Z
  }

  tau_raw <- cfg$tau_base + tau_indiv
  tau <- pmin(cfg$tau_max, pmax(cfg$tau_min, tau_raw))

  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y <- Y0 - T_vec * tau

  list(S = S, Age = 50, geneA = G_obs, Z_true = Z_true, Z_obs = Z_obs,
       T_vec = T_vec, tau = tau, Y = Y)
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

# ============================================================
# FIT FUNCTIONS
# ============================================================
fit_predict_causal_forest <- function(X_train, T_train, Y_train, X_test, seed, n_train) {
  tryCatch({
    df_train <- data.frame(Y = Y_train, Trt = T_train, X_train, check.names = FALSE)
    fml <- as.formula(paste("Y ~ Trt |", paste(colnames(X_train), collapse = " + ")))
    fit <- cf_fit(fml, data = df_train, method = "grf", seed = as.integer(seed))
    tau_hat_raw <- -as.numeric(predict(fit, newdata = as.data.frame(X_test), type = "ite"))
    tau_hat <- pmin(20, pmax(-20, tau_hat_raw))
    n_clip <- sum(tau_hat != tau_hat_raw)
    clip_note <- if (n_clip > 0) sprintf(", clipped=%d", n_clip) else ""
    list(tau_hat = tau_hat, notes = sprintf("OK%s", clip_note))
  }, error = function(e) list(tau_hat = NA_real_, notes = conditionMessage(e)))
}

fit_predict_drlearner <- function(X_train, T_train, Y_train, X_test, seed, n_train) {
  tryCatch({
    cv <- if (n_train < 200L) 3L else 5L

    # 修正版を使用（cf_fit経由ではなく直接呼び出し）
    cat("    Using FIXED DR-Learner (standardize_y=TRUE, min_propensity=0.05)\n")

    fit <- cf_fit_drlearner(
      X = as.matrix(X_train),
      T = T_train,
      Y = Y_train,
      random_state = as.integer(seed),
      cv = cv,
      min_propensity = 0.05,      # FIX #2
      standardize_y = TRUE,        # FIX #4
      model_propensity = "auto",   # FIX #1
      model_regression = "auto",   # FIX #1
      model_final = "auto"         # FIX #3
    )

    # 予測（修正版のpredict関数を使用）
    tau_hat_raw <- -predict_cf_drlearner(
      list(fit = fit),
      newdata = as.data.frame(X_test),
      type = "ite"
    )

    # クリップ（念のため）
    tau_hat <- pmin(20, pmax(-20, tau_hat_raw))
    n_clip <- sum(tau_hat != tau_hat_raw)
    clip_note <- if (n_clip > 0) sprintf(", clipped=%d", n_clip) else ""

    list(tau_hat = tau_hat, notes = sprintf("OK (cv=%d, FIXED%s)", cv, clip_note))
  }, error = function(e) {
    list(tau_hat = NA_real_, notes = paste0("ERROR: ", conditionMessage(e)))
  })
}

# ============================================================
# MAIN LOOP
# ============================================================
results_list <- list()
pred_list <- list()

cat("\n=== Starting benchmark with FIXED DR-Learner ===\n\n")

for (rep in seq_len(config$R)) {
  for (idx_n in seq_along(config$n_values)) {
    n <- config$n_values[idx_n]

    for (idx_effect in seq_along(config$effect_types)) {
      effect_type <- config$effect_types[idx_effect]
      seed_condition <- config$seed_base + (rep - 1L) * 1000L + idx_n * 100L + idx_effect * 10L

      cat(sprintf("[Rep %d/%d] n=%d, effect=%s (seed=%d)\n",
                  rep, config$R, n, effect_type, seed_condition))

      dat <- generate_data_v04(n, effect_type, seed_condition, config)
      spl <- stratified_split(dat$T_vec, config$train_frac, seed_condition)

      tr_idx <- spl$train
      te_idx <- spl$test
      n_train <- length(tr_idx)

      for (covset_name in config$covsets) {
        if (covset_name == "SexOnly") {
          X_cols <- "S"
        } else {
          X_cols <- c("S", "geneA")
        }

        if (covset_name == "SexOnly") {
          X_train <- data.frame(S = dat$S[tr_idx])
          X_test  <- data.frame(S = dat$S[te_idx])
        } else {
          X_train <- data.frame(S = dat$S[tr_idx], geneA = dat$geneA[tr_idx])
          X_test  <- data.frame(S = dat$S[te_idx], geneA = dat$geneA[te_idx])
        }
        T_train <- dat$T_vec[tr_idx]
        Y_train <- dat$Y[tr_idx]
        tau_true_test <- dat$tau[te_idx]

        for (method_name in config$methods) {
          cat(sprintf("  %s + %s ... ", covset_name, method_name))

          result <- switch(
            method_name,
            causal_forest = fit_predict_causal_forest(X_train, T_train, Y_train, X_test, seed_condition, n_train),
            drlearner = fit_predict_drlearner(X_train, T_train, Y_train, X_test, seed_condition, n_train)
          )

          tau_hat <- result$tau_hat
          notes <- result$notes

          if (all(is.na(tau_hat))) {
            rmse_tau <- NA_real_
            cat(sprintf("FAILED: %s\n", notes))
          } else {
            rmse_tau <- sqrt(mean((tau_true_test - tau_hat)^2, na.rm = TRUE))
            cat(sprintf("RMSE=%.3f (%s)\n", rmse_tau, notes))
          }

          results_list[[length(results_list) + 1]] <- data.frame(
            rep = rep, n = n, effect_type = effect_type, covset = covset_name,
            method = method_name, seed_condition = seed_condition,
            rmse_tau = rmse_tau, notes = notes, stringsAsFactors = FALSE
          )

          for (i in seq_along(te_idx)) {
            pred_list[[length(pred_list) + 1]] <- data.frame(
              rep = rep, n = n, effect_type = effect_type, covset = covset_name,
              method = method_name, seed_condition = seed_condition,
              id = te_idx[i], tau_true = tau_true_test[i],
              tau_hat = if (all(is.na(tau_hat))) NA_real_ else tau_hat[i],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
}

# ============================================================
# SAVE RESULTS
# ============================================================
df_long <- do.call(rbind, results_list)
df_pred <- do.call(rbind, pred_list)

write.csv(df_long, file.path(config$out_dir, "results_demov04_long_FIXED.csv"),
          row.names = FALSE)
write.csv(df_pred, file.path(config$out_dir, "predictions_demov04_FIXED.csv"),
          row.names = FALSE)

# Summary
df_summary <- aggregate(rmse_tau ~ n + effect_type + covset + method,
                       data = df_long, FUN = function(x) {
                         c(mean = mean(x, na.rm = TRUE),
                           sd = sd(x, na.rm = TRUE),
                           n = sum(!is.na(x)))
                       })
df_summary <- do.call(data.frame, df_summary)
names(df_summary)[5:7] <- c("mean_rmse", "sd_rmse", "n_success")

write.csv(df_summary, file.path(config$out_dir, "results_demov04_summary_FIXED.csv"),
          row.names = FALSE)

cat("\n=== Results saved ===\n")
cat("Directory:", config$out_dir, "\n")
print(df_summary)

cat("\n=== Comparison: n=500, linear, SexGene, drlearner ===\n")
row_target <- subset(df_summary, n == 500 & effect_type == "linear" &
                     covset == "SexGene" & method == "drlearner")
if (nrow(row_target) > 0) {
  cat(sprintf("FIXED version RMSE: %.3f\n", row_target$mean_rmse))
  cat("Original version RMSE: 19.940\n")
  improvement <- 19.940 - row_target$mean_rmse
  cat(sprintf("Improvement: %.3f (%.1f%% reduction)\n",
              improvement, 100 * improvement / 19.940))
} else {
  cat("Target condition not found in results.\n")
}
