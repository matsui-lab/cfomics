#!/usr/bin/env Rscript
# ============================================================
# CATE Demo v03
# Age・Sex・geneA が治療効果（血圧低下）を修飾するCATE推定
# covset AS（geneA なし）vs ASG（geneA あり）の比較
#
# Script: sandbox/demov03/run_cate_demov03.R
# Run from repo root: Rscript sandbox/demov03/run_cate_demov03.R
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

cat("=== CATE Demo v03 ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

# ============================================================
# 1. CONFIG
# ============================================================
config <- list(
  seed_base    = 20260215L,
  R            = 1L,            # 反復回数（本番は20）
  n_values     = c(100L, 500L, 1000L),
  covsets      = c("AS", "ASG"),
  methods      = c("causal_forest", "drlearner", "ganite", "cevae"),
  train_frac   = 0.8,
  out_dir      = "sandbox/demov03/results",

  # DGP: geneA 構造
  rho_A        = 0.4,   # geneA と Age の相関
  rho_S        = 0.2,   # geneA と Sex の相関
  sigma_G      = 0.2,   # geneA 測定誤差 SD

  # τ の重み
  w_A          = 1.0,
  w_S          = 1.0,
  w_G          = 1.0,
  s_tau        = 5.0,   # τ の SD スケール（mmHg）

  # 予後の重み
  beta_A       = 3.0,
  beta_S       = 2.0,
  s_eta        = 10.0,  # η の SD スケール（mmHg）
  mu0          = 120.0, # ベースライン血圧
  sigma_y      = 5.0    # アウトカムノイズ SD
)

dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. DGP
# ============================================================
generate_data_v03 <- function(n, seed, cfg) {
  set.seed(seed)

  # --- 共変量 ---
  A  <- rnorm(n, 0, 1)
  S  <- rbinom(n, 1, 0.5)
  Sc <- S - 0.5   # centered

  # --- geneA（真値 + 測定誤差）---
  coeff_Z <- sqrt(max(1 - cfg$rho_A^2 - cfg$rho_S^2, 0))
  Z       <- rnorm(n, 0, 1)
  G_true  <- cfg$rho_A * A + cfg$rho_S * Sc + coeff_Z * Z
  e_G     <- rnorm(n, 0, cfg$sigma_G)
  G_obs   <- G_true + e_G   # モデルへの入力（観測値）

  # --- 治療割付（RCT）---
  T_vec <- rbinom(n, 1, 0.5)

  # --- 真のCATE τ（血圧低下量）---
  tau_raw  <- cfg$w_A * A + cfg$w_S * Sc + cfg$w_G * G_true
  sd_tau   <- sd(tau_raw)
  notes_dgp <- ""
  if (sd_tau < 1e-8) {
    tau <- rep(0, n)
    notes_dgp <- "SD(tau_raw)<1e-8; tau=0"
    warning("SD(tau_raw) < 1e-8; setting tau = 0")
  } else {
    tau <- tau_raw * (cfg$s_tau / sd_tau)
  }

  # --- 予後 η（ベースライン血圧に影響、geneA 非依存）---
  eta_raw <- cfg$beta_A * A + cfg$beta_S * Sc
  sd_eta  <- sd(eta_raw)
  eta     <- if (sd_eta < 1e-8) rep(0, n) else eta_raw * (cfg$s_eta / sd_eta)

  # --- アウトカム（Y = Y(0) - T*τ：低下量がプラス）---
  eps_y <- rnorm(n, 0, cfg$sigma_y)
  Y0    <- cfg$mu0 + eta + eps_y
  Y1    <- Y0 - tau      # 治療あり：血圧低下
  Y     <- Y0 - T_vec * tau

  cat(sprintf("    [DGP] tau: mean=%.2f, SD=%.2f (target=%.1f), range=[%.2f, %.2f]\n",
              mean(tau), sd(tau), cfg$s_tau, min(tau), max(tau)))
  cat(sprintf("    [DGP] eta: SD=%.2f (target=%.1f) | T=1: %d, T=0: %d\n",
              sd(eta), cfg$s_eta, sum(T_vec), sum(T_vec == 0)))

  list(
    A = A, S = S, Sc = Sc, G_true = G_true, G_obs = G_obs,
    T_vec = T_vec, tau = tau, Y0 = Y0, Y1 = Y1, Y = Y,
    notes_dgp       = notes_dgp,
    sd_tau_realized = sd(tau),
    sd_eta_realized = sd(eta),
    tau_mean        = mean(tau),
    tau_min         = min(tau),
    tau_max         = max(tau)
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
# 4. PREPROCESSING（z-score, リーク防止）
# ============================================================
make_X <- function(dat, train_idx, covset) {
  # A: z-score（train統計量）
  mu_A <- mean(dat$A[train_idx]); sd_A <- max(sd(dat$A[train_idx]), 1e-6)
  A_z  <- (dat$A - mu_A) / sd_A

  # S: 0/1のまま（バイナリ）
  S_bin <- dat$S

  if (covset == "AS") {
    X <- cbind(A = A_z, S = S_bin)
  } else {
    # G: z-score（train統計量）
    mu_G <- mean(dat$G_obs[train_idx]); sd_G <- max(sd(dat$G_obs[train_idx]), 1e-6)
    G_z  <- (dat$G_obs - mu_G) / sd_G
    X    <- cbind(A = A_z, S = S_bin, G = G_z)
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

# --- 6.1 Causal Forest (GRF) ---
fit_predict_grf <- function(X_train, T_train, Y_train, X_test, seed) {
  tryCatch({
    df_train <- data.frame(Y = Y_train, Trt = T_train, X_train, check.names = FALSE)
    fml <- as.formula(paste("Y ~ Trt |", paste(colnames(X_train), collapse = " + ")))
    fit <- cf_fit(fml, data = df_train, method = "grf",
                  num.trees = 2000L, random_state = as.integer(seed))
    tau_hat <- as.numeric(predict(fit, newdata = as.data.frame(X_test), type = "ite"))
    list(tau_hat = tau_hat, notes = "OK")
  }, error = function(e) list(tau_hat = NA_real_, notes = conditionMessage(e)))
}

# --- 6.2 DR-Learner ---
fit_predict_drlearner <- function(X_train, T_train, Y_train, X_test, seed, n_train) {
  tryCatch({
    cv       <- if (n_train < 200L) 3L else 5L
    df_train <- data.frame(Y = Y_train, Trt = T_train, X_train, check.names = FALSE)
    fml <- as.formula(paste("Y ~ Trt |", paste(colnames(X_train), collapse = " + ")))
    fit <- cf_fit(fml, data = df_train, method = "drlearner",
                  random_state = as.integer(seed), cv = cv)
    tau_hat <- as.numeric(predict(fit, newdata = as.data.frame(X_test), type = "ite"))
    list(tau_hat = tau_hat, notes = sprintf("OK (cv=%d)", cv))
  }, error = function(e) list(tau_hat = NA_real_, notes = conditionMessage(e)))
}

# --- 6.3 GANITE ---
fit_predict_ganite <- function(X_train, T_train, Y_train, X_test, seed, n_train) {
  tryCatch({
    ganite_mod <- reticulate::import_from_path(
      "ganite", path = system.file("python", package = "cfomics")
    )
    np <- reticulate::import("numpy", convert = FALSE)

    X_tr_np <- np$array(X_train,             dtype = "float32")
    T_tr_np <- np$array(as.numeric(T_train),  dtype = "float32")
    Y_tr_np <- np$array(as.numeric(Y_train),  dtype = "float32")
    X_te_np <- np$array(X_test,              dtype = "float32")

    batch_sz <- as.integer(max(16L, min(n_train %/% 4L, 64L)))
    iters    <- if (n_train < 200L) 1000L else 2000L

    model <- ganite_mod$GANITEModel(
      input_dim    = as.integer(ncol(X_train)),
      h_dim        = 64L,
      alpha        = 1.0,
      random_state = as.integer(seed)
    )
    model$fit(X_tr_np, T_tr_np, Y_tr_np,
              iterations = iters, batch_size = batch_sz, verbose = FALSE)

    y_hat   <- reticulate::py_to_r(model$predict(X_te_np))
    tau_hat <- y_hat[, 2] - y_hat[, 1]

    list(tau_hat = as.numeric(tau_hat),
         notes = sprintf("OK (iters=%d, batch=%d)", iters, batch_sz))
  }, error = function(e) list(tau_hat = NA_real_, notes = conditionMessage(e)))
}

# --- 6.4 CEVAE ---
fit_predict_cevae <- function(X_train, T_train, Y_train, X_test, seed, n_train) {
  tryCatch({
    torch     <- reticulate::import("torch",  convert = FALSE)
    pyro      <- reticulate::import("pyro",   convert = FALSE)
    cevae_mod <- reticulate::import("pyro.contrib.cevae", convert = FALSE)
    CEVAE     <- cevae_mod$CEVAE

    pyro$set_rng_seed(as.integer(seed))

    num_epochs <- 200L
    batch_sz   <- as.integer(min(n_train, 64L))

    x_tr <- torch$tensor(X_train,             dtype = torch$float32)
    t_tr <- torch$tensor(as.numeric(T_train),  dtype = torch$float32)
    y_tr <- torch$tensor(as.numeric(Y_train),  dtype = torch$float32)
    x_te <- torch$tensor(X_test,              dtype = torch$float32)

    cevae <- CEVAE(
      feature_dim  = as.integer(ncol(X_train)),
      outcome_dist = "normal",
      latent_dim   = 5L,
      hidden_dim   = 64L,
      num_layers   = 2L,
      num_samples  = 10L
    )
    cevae$fit(x_tr, t_tr, y_tr,
              num_epochs          = num_epochs,
              batch_size          = batch_sz,
              learning_rate       = 1e-3,
              learning_rate_decay = 0.1,
              weight_decay        = 1e-4)

    tau_hat <- as.numeric(cevae$ite(x_te)$detach()$cpu()$numpy())

    list(tau_hat = tau_hat,
         notes = sprintf("OK (epochs=%d, batch=%d)", num_epochs, batch_sz))
  }, error = function(e) list(tau_hat = NA_real_, notes = conditionMessage(e)))
}

# ============================================================
# 7. MAIN LOOP
# ============================================================
results_list <- list()

for (rep in seq_len(config$R)) {
  for (idx_n in seq_along(config$n_values)) {
    n <- config$n_values[idx_n]

    # seed: (rep, n) ごとにユニーク
    seed_rep <- config$seed_base + (rep - 1L) * 10000L + (idx_n - 1L) * 100L

    cat(sprintf("\n=== rep=%d/%d, n=%d, seed=%d ===\n",
                rep, config$R, n, seed_rep))

    # --- データ生成（covset/method 共通）---
    dat    <- generate_data_v03(n, seed_rep, config)
    spl    <- stratified_split(dat$T_vec, config$train_frac, seed_rep)
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
          fit_predict_grf(X_train, T_train, Y_train, X_test, seed_rep),
        drlearner = function()
          fit_predict_drlearner(X_train, T_train, Y_train, X_test, seed_rep, n_train),
        ganite = function()
          fit_predict_ganite(X_train, T_train, Y_train, X_test, seed_rep, n_train),
        cevae = function()
          fit_predict_cevae(X_train, T_train, Y_train, X_test, seed_rep, n_train)
      )

      for (method in config$methods) {
        cat(sprintf("  [%s | %s] fitting...\n", covset, method))
        t0    <- proc.time()["elapsed"]
        res_m <- method_fns[[method]]()
        elapsed <- proc.time()["elapsed"] - t0

        tau_hat_vec <- res_m$tau_hat
        if (length(tau_hat_vec) == n_test && !anyNA(tau_hat_vec)) {
          rmse_val <- rmse_fn(tau_hat_vec, tau_test_true)
          cat(sprintf("  [%s | %s] %.1fs | %s | RMSE=%.4f\n",
                      covset, method, elapsed, res_m$notes, rmse_val))
        } else {
          tau_hat_vec <- rep(NA_real_, n_test)
          rmse_val    <- NA_real_
          cat(sprintf("  [%s | %s] %.1fs | FAILED: %s\n",
                      covset, method, elapsed, res_m$notes))
        }

        notes_full <- paste0(
          if (nchar(dat$notes_dgp) > 0) paste0(dat$notes_dgp, "; ") else "",
          res_m$notes
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          seed_base        = config$seed_base,
          seed_rep         = seed_rep,
          rep              = rep,
          n                = n,
          covset           = covset,
          method           = method,
          rmse_tau         = rmse_val,
          n_train          = n_train,
          n_test           = n_test,
          notes            = notes_full,
          sd_tau_realized  = dat$sd_tau_realized,
          sd_eta_realized  = dat$sd_eta_realized,
          tau_mean         = dat$tau_mean,
          tau_min          = dat$tau_min,
          tau_max          = dat$tau_max,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# ============================================================
# 8. SAVE OUTPUTS
# ============================================================
cat("\n=== Saving outputs ===\n")

df_long <- do.call(rbind, results_list)
path_long <- file.path(config$out_dir, "results_demov03_long.csv")
write.csv(df_long, path_long, row.names = FALSE)
cat(sprintf("  Saved: results_demov03_long.csv (%d rows)\n", nrow(df_long)))

# サマリー（R=1 では sd/se は NA）
df_valid <- subset(df_long, !is.na(rmse_tau))
df_summary <- do.call(rbind, lapply(
  split(df_valid, list(df_valid$n, df_valid$covset, df_valid$method), drop = TRUE),
  function(g) {
    data.frame(
      n         = g$n[1],
      covset    = g$covset[1],
      method    = g$method[1],
      mean_rmse = mean(g$rmse_tau),
      sd_rmse   = if (nrow(g) > 1) sd(g$rmse_tau) else NA_real_,
      se_rmse   = if (nrow(g) > 1) sd(g$rmse_tau) / sqrt(nrow(g)) else NA_real_,
      n_success = nrow(g),
      stringsAsFactors = FALSE
    )
  }
))
df_summary <- df_summary[order(df_summary$covset, df_summary$n, df_summary$method), ]
path_summary <- file.path(config$out_dir, "results_demov03_summary.csv")
write.csv(df_summary, path_summary, row.names = FALSE)
cat(sprintf("  Saved: results_demov03_summary.csv (%d rows)\n", nrow(df_summary)))

# session info
sink(file.path(config$out_dir, "session_info.txt"))
cat("=== R Session Info ===\n"); print(sessionInfo())
cat("\n=== Config ===\n"); print(config)
cat(sprintf("\n=== cfomics git hash ===\n"))
cat(system("git rev-parse HEAD", intern = TRUE), "\n")
sink()
cat("  Saved: session_info.txt\n")

# ============================================================
# 9. PRINT SUMMARY TABLE
# ============================================================
cat("\n=================================================================\n")
cat("                     RESULTS SUMMARY\n")
cat("=================================================================\n")
for (cs in config$covsets) {
  cat(sprintf("\n--- covset: %s ---\n", cs))
  sub    <- subset(df_summary, covset == cs)
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

cat(sprintf("\n=== Finished: %s ===\n", Sys.time()))
cat("=================================================================\n")
