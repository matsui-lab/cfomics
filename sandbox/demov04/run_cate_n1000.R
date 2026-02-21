#!/usr/bin/env Rscript
# run_cate_n1000.R
# n=1000での解析

suppressPackageStartupMessages({
  venv_python <- file.path(getwd(), "sandbox", ".venv", "bin", "python3")
  Sys.setenv(RETICULATE_PYTHON = venv_python)
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0")
  devtools::load_all("packages/cfomics", quiet = TRUE)
})

cat("Loading FIXED DR-Learner implementation...\n")
source("sandbox/demov04/methods_drlearner_FIXED.R")

# ============================================================
# CONFIG
# ============================================================
config <- list(
  seed_base  = 20260215L,
  R          = 1L,
  n_values   = c(100L, 500L, 1000L),
  effect_types = c("linear"),
  covsets    = c("SexOnly", "SexGene"),
  methods    = c("causal_forest", "drlearner"),
  train_frac = 0.8,
  out_dir    = "sandbox/demov04/results_n1000",

  # DGP params
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0,
  mu0 = 120.0, sigma_y = 5.0
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# DGP
# ============================================================
generate_data_corrected <- function(n, effect_type, seed, cfg) {
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

  if (effect_type == "linear") {
    tau_indiv <- cfg$a_S * Sc + cfg$a_G * Z_obs
  }

  tau_raw <- cfg$tau_base + tau_indiv
  tau <- pmin(cfg$tau_max, pmax(cfg$tau_min, tau_raw))

  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y <- Y0 + T_vec * tau

  list(S = S, Age = 50, geneA = G_obs, Z_obs = Z_obs,
       T_vec = T_vec, tau = tau, Y = Y, Y0 = Y0)
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
    tau_hat_raw <- as.numeric(predict(fit, newdata = as.data.frame(X_test), type = "ite"))
    tau_hat <- pmin(20, pmax(-20, tau_hat_raw))
    n_clip <- sum(tau_hat != tau_hat_raw)
    clip_note <- if (n_clip > 0) sprintf(", clipped=%d", n_clip) else ""
    list(tau_hat = tau_hat, notes = sprintf("OK%s", clip_note))
  }, error = function(e) list(tau_hat = NA_real_, notes = conditionMessage(e)))
}

fit_predict_drlearner <- function(X_train, T_train, Y_train, X_test, seed, n_train) {
  tryCatch({
    cv <- if (n_train < 200L) 3L else 5L

    fit <- cf_fit_drlearner(
      X = as.matrix(X_train),
      T = T_train,
      Y = Y_train,
      random_state = as.integer(seed),
      cv = cv,
      min_propensity = 0.05,
      standardize_y = TRUE,
      model_propensity = "auto",
      model_regression = "auto",
      model_final = "auto"
    )

    tau_hat_raw <- predict_cf_drlearner(
      list(fit = fit),
      newdata = as.data.frame(X_test),
      type = "ite"
    )

    tau_hat <- pmin(20, pmax(-20, tau_hat_raw))
    n_clip <- sum(tau_hat != tau_hat_raw)
    clip_note <- if (n_clip > 0) sprintf(", clipped=%d", n_clip) else ""

    list(tau_hat = tau_hat, notes = sprintf("OK (cv=%d%s)", cv, clip_note))
  }, error = function(e) {
    list(tau_hat = NA_real_, notes = paste0("ERROR: ", conditionMessage(e)))
  })
}

# ============================================================
# MAIN LOOP
# ============================================================
results_list <- list()
pred_list <- list()

cat("\n=== n=1000での解析 ===\n\n")

for (rep in seq_len(config$R)) {
  for (idx_n in seq_along(config$n_values)) {
    n <- config$n_values[idx_n]

    # n=1000のみ実行
    if (n != 1000L) next

    for (idx_effect in seq_along(config$effect_types)) {
      effect_type <- config$effect_types[idx_effect]
      seed_condition <- config$seed_base + (rep - 1L) * 1000L + idx_n * 100L + idx_effect * 10L

      cat(sprintf("[Rep %d/%d] n=%d, effect=%s (seed=%d)\n",
                  rep, config$R, n, effect_type, seed_condition))

      dat <- generate_data_corrected(n, effect_type, seed_condition, config)
      spl <- stratified_split(dat$T_vec, config$train_frac, seed_condition)

      tr_idx <- spl$train
      te_idx <- spl$test
      n_train <- length(tr_idx)

      cat(sprintf("  Train: %d, Test: %d\n", n_train, length(te_idx)))

      for (covset_name in config$covsets) {
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

      # 全データも保存
      saveRDS(dat, file.path(config$out_dir, "full_data_n1000.rds"))
      saveRDS(spl, file.path(config$out_dir, "split_n1000.rds"))
    }
  }
}

# ============================================================
# SAVE RESULTS
# ============================================================
df_long <- do.call(rbind, results_list)
df_pred <- do.call(rbind, pred_list)

write.csv(df_long, file.path(config$out_dir, "results_n1000_long.csv"),
          row.names = FALSE)
write.csv(df_pred, file.path(config$out_dir, "predictions_n1000.csv"),
          row.names = FALSE)

df_summary <- aggregate(rmse_tau ~ n + effect_type + covset + method,
                       data = df_long, FUN = function(x) {
                         c(mean = mean(x, na.rm = TRUE),
                           sd = sd(x, na.rm = TRUE),
                           n = sum(!is.na(x)))
                       })
df_summary <- do.call(data.frame, df_summary)
names(df_summary)[5:7] <- c("mean_rmse", "sd_rmse", "n_success")

write.csv(df_summary, file.path(config$out_dir, "results_n1000_summary.csv"),
          row.names = FALSE)

cat("\n=== Results saved (n=1000) ===\n")
print(df_summary)

cat("\n✓ n=1000での解析完了\n")
