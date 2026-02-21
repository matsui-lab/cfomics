#!/usr/bin/env Rscript
# ============================================================
# CATE Demo v04
# 性別 + geneA(NBカウント) による CATE 推定
# covset SexOnly vs SexGene × effect_type linear/threshold3
#
# Script: sandbox/demov04/run_cate_demov04.R
# Run from repo root: Rscript sandbox/demov04/run_cate_demov04.R
# ============================================================

# ============================================================
# 0. ENVIRONMENT SETUP
# ============================================================
suppressPackageStartupMessages({
  venv_python <- file.path(getwd(), "sandbox", ".venv", "bin", "python3")
  Sys.setenv(RETICULATE_PYTHON = venv_python)
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0")
  devtools::load_all("packages/cfomics", quiet = TRUE)
})

cat("=== CATE Demo v04 ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

# ============================================================
# 1. CONFIG
# ============================================================
config <- list(
  seed_base     = 20260215L,
  R             = 1L,
  n_values      = c(100L, 500L, 1000L),
  effect_types  = c("linear", "threshold3"),
  covsets       = c("SexOnly", "SexGene"),
  methods       = c("causal_forest", "drlearner"),
  train_frac    = 0.8,
  out_dir       = "sandbox/demov04/results",

  # geneA (NB count) params - 案1: 高品質RNA-seq
  m_E       = log(50),   # 潜在ログ平均
  s_E       = 0.35,      # 潜在ログSD（生物学的ばらつき）
  theta_nb  = 30,        # NB 過分散パラメータ（サンプリングノイズ小）
  sigma_G   = 0.08,      # 測定誤差 SD（技術的再現性高）

  # effect params - 案A: 現実的設定
  tau_base  = -12.0,     # 全員に共通の平均降圧効果 (mmHg)
  a_S       = 2.0,       # τ への性別効果: 2 mmHg 差
  a_G       = 1.5,       # τ への geneA 効果: 6 mmHg 差（Z範囲±2として）
  delta_thr = 1.5,       # threshold3 用の段階幅
  tau_min   = -18.0,     # τ クリップ下限
  tau_max   = -6.0,      # τ クリップ上限

  # outcome params
  mu0       = 120.0,     # ベースライン血圧
  sigma_y   = 2.0        # アウトカムノイズ SD
)

dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. DGP
# ============================================================
generate_data_v04 <- function(n, effect_type, seed, cfg) {
  set.seed(seed)

  # --- 2.1 性別 ---
  S  <- rbinom(n, 1, 0.5)
  Sc <- S - 0.5   # centered

  # --- 2.2 geneA (NB count → log1p + 測定誤差) ---
  E_true   <- rnorm(n, cfg$m_E, cfg$s_E)            # 潜在ログ平均
  mu_nb    <- exp(E_true)                             # NB mean
  C_count  <- rnbinom(n, size = cfg$theta_nb, mu = mu_nb)  # NBカウント
  G_raw    <- log1p(C_count)                          # log1p 変換
  delta    <- rnorm(n, 0, cfg$sigma_G)                # 測定誤差
  G_obs    <- G_raw + delta                           # 観測geneA

  # 潜在スコア（参考用、τ生成には使わない）
  Z_true <- (E_true - cfg$m_E) / cfg$s_E

  # --- 2.3 治療割付 (RCT) ---
  T_vec <- rbinom(n, 1, 0.5)

  # --- 2.4 真のCATE τ（観測値 G_obs から生成）---
  # IMPORTANT: τ は観測可能な G_obs から生成される
  # これにより測定ノイズの影響がなくなり、SexGene の改善が明確になる

  # G_obs を標準化（平均0、SD=1）
  Z_obs <- (G_obs - mean(G_obs)) / sd(G_obs)

  if (effect_type == "linear") {
    # 基本効果 + 性別 + geneA による個人差
    tau_indiv <- cfg$a_S * Sc + cfg$a_G * Z_obs

  } else if (effect_type == "threshold3") {
    # 基本効果 + 性別 + geneA 3段階
    q1 <- quantile(Z_obs, 0.33)
    q2 <- quantile(Z_obs, 0.66)
    g_Z <- ifelse(Z_obs <= q1, -cfg$delta_thr,
                  ifelse(Z_obs <= q2, 0, cfg$delta_thr))
    tau_indiv <- cfg$a_S * Sc + g_Z

  } else {
    stop(sprintf("Unknown effect_type: %s", effect_type))
  }

  # tau = 基本降圧効果 + 個人差
  tau_raw <- cfg$tau_base + tau_indiv
  sd_tau_raw <- sd(tau_raw)

  # クリップ（ほぼ全員マイナス）
  tau <- pmin(cfg$tau_max, pmax(cfg$tau_min, tau_raw))
  sd_tau_clipped <- sd(tau)
  n_clipped <- sum(tau_raw != tau)
  notes_dgp <- ""

  # --- 2.5 アウトカム ---
  eps_y <- rnorm(n, 0, cfg$sigma_y)
  Y0    <- cfg$mu0 + eps_y
  Y1    <- Y0 - tau       # 治療あり: 血圧低下
  Y     <- Y0 - T_vec * tau

  # ログ
  cat(sprintf("    [DGP] tau: mean=%.2f, SD=%.2f, range=[%.2f, %.2f]\n",
              mean(tau), sd_tau_clipped, min(tau), max(tau)))
  cat(sprintf("    [DGP] clipped: %d/%d samples | T=1: %d, T=0: %d\n",
              n_clipped, n, sum(T_vec), sum(T_vec == 0)))
  if (effect_type == "threshold3") {
    tbl <- table(cut(Z_true, breaks = c(-Inf, q1, q2, Inf),
                     labels = c("low", "mid", "high")))
    cat(sprintf("    [DGP] threshold3 counts: %s\n",
                paste(names(tbl), tbl, sep = "=", collapse = ", ")))
  }

  list(
    S = S, Sc = Sc, E_true = E_true, Z_true = Z_true, Z_obs = Z_obs,
    C_count = C_count, G_obs = G_obs,
    T_vec = T_vec, tau = tau, Y0 = Y0, Y1 = Y1, Y = Y,
    notes_dgp        = notes_dgp,
    sd_tau_clipped   = sd_tau_clipped,
    tau_mean         = mean(tau),
    tau_min          = min(tau),
    tau_max          = max(tau),
    n_clipped        = n_clipped
  )
}

# ============================================================
# 3. STRATIFIED SPLIT
# ============================================================
stratified_split <- function(T_vec, train_frac, seed) {
  set.seed(seed)
  idx0  <- which(T_vec == 0); idx1 <- which(T_vec == 1)
  n0_tr <- round(length(idx0) * train_frac)
  n1_tr <- round(length(idx1) * train_frac)
  tr_idx <- sort(c(sample(idx0)[seq_len(n0_tr)], sample(idx1)[seq_len(n1_tr)]))
  te_idx <- sort(setdiff(seq_along(T_vec), tr_idx))
  list(train = tr_idx, test = te_idx)
}

# ============================================================
# 4. PREPROCESSING
# ============================================================
make_X <- function(dat, train_idx, covset) {
  S_bin <- dat$S

  if (covset == "SexOnly") {
    X <- matrix(S_bin, ncol = 1)
    colnames(X) <- "S"
  } else {
    # G_obs: z-score (train 統計量でリーク防止)
    mu_G <- mean(dat$G_obs[train_idx])
    sd_G <- max(sd(dat$G_obs[train_idx]), 1e-6)
    G_z  <- (dat$G_obs - mu_G) / sd_G
    X    <- cbind(S = S_bin, G = G_z)
  }

  list(
    X_train = X[train_idx,  , drop = FALSE],
    X_test  = X[-train_idx, , drop = FALSE]
  )
}

# ============================================================
# 5. RMSE
# ============================================================
rmse_fn <- function(pred, true) sqrt(mean((pred - true)^2))

# ============================================================
# 6. METHOD WRAPPERS
# ============================================================

# --- Causal Forest (GRF) ---
fit_predict_grf <- function(X_train, T_train, Y_train, X_test, seed) {
  tryCatch({
    df_train <- data.frame(Y = Y_train, Trt = T_train, X_train, check.names = FALSE)
    fml <- as.formula(paste("Y ~ Trt |", paste(colnames(X_train), collapse = " + ")))
    fit <- cf_fit(fml, data = df_train, method = "grf",
                  num.trees = 2000L, random_state = as.integer(seed))
    # 符号反転: モデルは E[Y(1)-Y(0)|X] = -τ を推定するが、
    # DGP は Y = Y(0) - T*τ (τ>0 = 血圧低下) なので negate する
    tau_hat <- -as.numeric(predict(fit, newdata = as.data.frame(X_test), type = "ite"))
    list(tau_hat = tau_hat, notes = "OK")
  }, error = function(e) list(tau_hat = NA_real_, notes = conditionMessage(e)))
}

# --- DR-Learner ---
fit_predict_drlearner <- function(X_train, T_train, Y_train, X_test, seed, n_train) {
  tryCatch({
    cv       <- if (n_train < 200L) 3L else 5L
    df_train <- data.frame(Y = Y_train, Trt = T_train, X_train, check.names = FALSE)
    fml <- as.formula(paste("Y ~ Trt |", paste(colnames(X_train), collapse = " + ")))
    fit <- cf_fit(fml, data = df_train, method = "drlearner",
                  random_state = as.integer(seed), cv = cv)
    # 同上: 符号反転 + クリップ（DR疑似スコア爆発対策）
    tau_hat_raw <- -as.numeric(predict(fit, newdata = as.data.frame(X_test), type = "ite"))
    tau_hat <- pmin(20, pmax(-20, tau_hat_raw))
    n_clip  <- sum(tau_hat != tau_hat_raw)
    clip_note <- if (n_clip > 0) sprintf(", clipped=%d", n_clip) else ""
    list(tau_hat = tau_hat, notes = sprintf("OK (cv=%d%s)", cv, clip_note))
  }, error = function(e) list(tau_hat = NA_real_, notes = conditionMessage(e)))
}

# ============================================================
# 7. MAIN LOOP
# ============================================================
results_list <- list()
pred_list    <- list()

for (rep in seq_len(config$R)) {
  for (idx_n in seq_along(config$n_values)) {
    n <- config$n_values[idx_n]

    for (idx_effect in seq_along(config$effect_types)) {
      effect_type <- config$effect_types[idx_effect]

      # seed: (rep, n, effect_type) ごとにユニーク
      seed_cond <- config$seed_base +
        (rep - 1L) * 100000L +
        (idx_n - 1L) * 10000L +
        (idx_effect - 1L) * 100L

      cat(sprintf("\n=== rep=%d, n=%d, effect=%s, seed=%d ===\n",
                  rep, n, effect_type, seed_cond))
      cat("  Generating data...\n")

      dat    <- generate_data_v04(n, effect_type, seed_cond, config)
      spl    <- stratified_split(dat$T_vec, config$train_frac, seed_cond)
      n_train <- length(spl$train)
      n_test  <- length(spl$test)
      cat(sprintf("  n_train=%d, n_test=%d\n", n_train, n_test))

      T_train       <- dat$T_vec[spl$train]
      Y_train       <- dat$Y[spl$train]
      tau_test_true <- dat$tau[spl$test]

      for (covset in config$covsets) {
        pre     <- make_X(dat, spl$train, covset)
        X_train <- pre$X_train
        X_test  <- pre$X_test

        method_fns <- list(
          causal_forest = function()
            fit_predict_grf(X_train, T_train, Y_train, X_test, seed_cond),
          drlearner = function()
            fit_predict_drlearner(X_train, T_train, Y_train, X_test,
                                  seed_cond, n_train)
        )

        for (method in config$methods) {
          cat(sprintf("  [%s | %s | %s] fitting...\n", effect_type, covset, method))
          t0    <- proc.time()["elapsed"]
          res_m <- method_fns[[method]]()
          elapsed <- proc.time()["elapsed"] - t0

          tau_hat_vec <- res_m$tau_hat
          if (length(tau_hat_vec) == n_test && !anyNA(tau_hat_vec)) {
            rmse_val <- rmse_fn(tau_hat_vec, tau_test_true)
            cat(sprintf("  [%s | %s | %s] %.1fs | %s | RMSE=%.4f\n",
                        effect_type, covset, method, elapsed, res_m$notes, rmse_val))
          } else {
            tau_hat_vec <- rep(NA_real_, n_test)
            rmse_val    <- NA_real_
            cat(sprintf("  [%s | %s | %s] %.1fs | FAILED: %s\n",
                        effect_type, covset, method, elapsed, res_m$notes))
          }

          notes_full <- paste0(
            if (nchar(dat$notes_dgp) > 0) paste0(dat$notes_dgp, "; ") else "",
            res_m$notes
          )

          results_list[[length(results_list) + 1]] <- data.frame(
            seed_base        = config$seed_base,
            seed_condition   = seed_cond,
            rep              = rep,
            n                = n,
            effect_type      = effect_type,
            covset           = covset,
            method           = method,
            rmse_tau         = rmse_val,
            n_train          = n_train,
            n_test           = n_test,
            notes            = notes_full,
            sd_tau_clipped   = dat$sd_tau_clipped,
            tau_mean         = dat$tau_mean,
            tau_min_realized = dat$tau_min,
            tau_max_realized = dat$tau_max,
            n_clipped        = dat$n_clipped,
            stringsAsFactors = FALSE
          )

          # 予測値保存
          pred_list[[length(pred_list) + 1]] <- data.frame(
            seed_condition = seed_cond,
            rep            = rep,
            n              = n,
            effect_type    = effect_type,
            covset         = covset,
            method         = method,
            id             = spl$test,
            split          = "test",
            tau_true       = tau_test_true,
            tau_hat        = tau_hat_vec,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

# ============================================================
# 8. SAVE OUTPUTS
# ============================================================
cat("\n=== Saving outputs ===\n")

df_long <- do.call(rbind, results_list)
path_long <- file.path(config$out_dir, "results_demov04_long.csv")
write.csv(df_long, path_long, row.names = FALSE)
cat(sprintf("  Saved: results_demov04_long.csv (%d rows)\n", nrow(df_long)))

df_pred <- do.call(rbind, pred_list)
path_pred <- file.path(config$out_dir, "predictions_demov04.csv")
write.csv(df_pred, path_pred, row.names = FALSE)
cat(sprintf("  Saved: predictions_demov04.csv (%d rows)\n", nrow(df_pred)))

# サマリー
df_valid   <- subset(df_long, !is.na(rmse_tau))
df_summary <- do.call(rbind, lapply(
  split(df_valid,
        list(df_valid$n, df_valid$effect_type, df_valid$covset, df_valid$method),
        drop = TRUE),
  function(g) {
    data.frame(
      n           = g$n[1],
      effect_type = g$effect_type[1],
      covset      = g$covset[1],
      method      = g$method[1],
      mean_rmse   = mean(g$rmse_tau),
      sd_rmse     = if (nrow(g) > 1) sd(g$rmse_tau) else NA_real_,
      se_rmse     = if (nrow(g) > 1) sd(g$rmse_tau) / sqrt(nrow(g)) else NA_real_,
      n_success   = nrow(g),
      stringsAsFactors = FALSE
    )
  }
))
df_summary <- df_summary[order(df_summary$effect_type, df_summary$covset,
                                df_summary$n, df_summary$method), ]
path_summary <- file.path(config$out_dir, "results_demov04_summary.csv")
write.csv(df_summary, path_summary, row.names = FALSE)
cat(sprintf("  Saved: results_demov04_summary.csv (%d rows)\n", nrow(df_summary)))

# session info
sink(file.path(config$out_dir, "session_info.txt"))
cat("=== R Session Info ===\n"); print(sessionInfo())
cat("\n=== Config ===\n"); print(config)
cat(sprintf("\n=== cfomics git hash ===\n"))
cat(system("git rev-parse HEAD", intern = TRUE), "\n")
sink()
cat("  Saved: session_info.txt\n")

# ============================================================
# 9. VISUALIZATION
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

method_labels <- c(causal_forest = "Causal Forest", drlearner = "DR-Learner")
method_colors <- c("Causal Forest" = "#E41A1C", "DR-Learner" = "#377EB8")

# サマリーにラベル付き factor
df_summary$method_label <- factor(df_summary$method,
                                   levels = names(method_labels),
                                   labels = method_labels)
df_summary$covset <- factor(df_summary$covset, levels = c("SexOnly", "SexGene"))

# --- 折れ線 (effect_type 別) ---
for (et in config$effect_types) {
  d <- subset(df_summary, effect_type == et)

  p <- ggplot(d, aes(x = n, y = mean_rmse, color = method_label,
                      linetype = covset, group = interaction(method_label, covset))) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    scale_color_manual(values = method_colors, name = "Method") +
    scale_linetype_manual(values = c("SexOnly" = "dashed", "SexGene" = "solid"),
                          name = "Covariates") +
    scale_x_continuous(breaks = c(100, 500, 1000)) +
    labs(
      title    = sprintf("CATE RMSE — %s effect  [demov04 realistic, R=%d]", et, config$R),
      subtitle = "Solid: SexGene (S+G)  |  Dashed: SexOnly (S)  |  τ ≈ -8 mmHg, SD ≈ 0.5-0.7 mmHg",
      x        = "Sample size (n)",
      y        = expression(RMSE[tau] ~ "(mmHg)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = "right"
    )

  fname <- sprintf("plot_rmse_%s.png", et)
  ggsave(file.path(config$out_dir, fname), p, width = 9, height = 5, dpi = 150)
  cat(sprintf("  Saved: %s\n", fname))
}

# --- ヒートマップ (effect_type × covset) ---
df_heat <- df_summary
df_heat$label <- ifelse(is.na(df_heat$mean_rmse), "NA",
                        sprintf("%.2f", df_heat$mean_rmse))
df_heat$cell_val <- pmin(df_heat$mean_rmse, 20)
df_heat$panel <- paste0(df_heat$effect_type, " / ", df_heat$covset)

fig_heat <- ggplot(df_heat, aes(x = factor(n), y = method_label, fill = cell_val)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.5, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "#FFFFBF", high = "#D6604D",
                       midpoint = 8, limits = c(0, 20),
                       na.value = "grey80",
                       name = expression(RMSE[tau])) +
  facet_wrap(~panel, ncol = 2) +
  labs(
    title    = sprintf("CATE RMSE ヒートマップ  [demov04, R=%d]", config$R),
    subtitle = "青=良好 / 赤=不良  |  上限20キャップ",
    x = "n", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    strip.text    = element_text(size = 10, face = "bold")
  )

ggsave(file.path(config$out_dir, "plot_heatmap_rmse.png"),
       fig_heat, width = 11, height = 5, dpi = 150)
cat("  Saved: plot_heatmap_rmse.png\n")

# ============================================================
# 10. PRINT SUMMARY
# ============================================================
cat("\n=================================================================\n")
cat("                     RESULTS SUMMARY\n")
cat("=================================================================\n")

for (et in config$effect_types) {
  for (cs in config$covsets) {
    cat(sprintf("\n--- effect=%s, covset=%s ---\n", et, cs))
    sub <- subset(df_summary, effect_type == et & covset == cs)
    n_vals <- sort(unique(sub$n))

    header <- sprintf("  %15s", "method")
    for (nv in n_vals) header <- paste0(header, sprintf("  %8s", paste0("n=", nv)))
    cat(header, "\n")

    for (m in config$methods) {
      row_str <- sprintf("  %15s", m)
      for (nv in n_vals) {
        val <- sub$mean_rmse[sub$method == m & sub$n == nv]
        row_str <- paste0(row_str, sprintf("  %8s",
          if (length(val) > 0 && !is.na(val)) sprintf("%.4f", val) else "NA"))
      }
      cat(row_str, "\n")
    }
  }
}

cat(sprintf("\n=== Finished: %s ===\n", Sys.time()))
cat("=================================================================\n")
