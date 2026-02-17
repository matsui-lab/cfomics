#!/usr/bin/env Rscript
# ============================================================
# CATE Demo v0.2
# omics50 / two-layer RNA / log-mechanism / log-input only
#
# Script: sandbox/run_cate_demo_v02.R
# Run from repo root: Rscript sandbox/run_cate_demo_v02.R
# ============================================================

# ============================================================
# 0. ENVIRONMENT SETUP
# ============================================================
suppressPackageStartupMessages({
  # Use file.path(getwd(), ...) to avoid normalizePath() following symlinks
  # which would resolve to the system Python instead of the venv
  venv_python <- file.path(getwd(), "sandbox", ".venv", "bin", "python3")
  Sys.setenv(RETICULATE_PYTHON = venv_python)
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0")
  devtools::load_all("packages/cfomics", quiet = TRUE)
})

cat("=== CATE Demo v0.2 ===\n")
cat(sprintf("Started: %s\n", Sys.time()))

# ============================================================
# 1. CONFIG
# ============================================================
config <- list(
  seed_base    = 20260215L,
  n_values     = c(100L, 500L, 1000L),
  effect_types = c("linear", "threshold3"),
  methods      = c("causal_forest", "drlearner", "ganite", "cevae"),
  train_frac   = 0.8,
  out_dir      = "sandbox/benchmark_results/cate_demo_v02",

  # DGP
  p_genes      = 50L,
  K_latent     = 5L,
  sigma_lambda = 0.3,
  sigma_e      = 0.3,
  theta_nb     = 10,
  libsize_sd   = 0.2,

  # Gene sets
  M_genes = 1:5,   # predictive (effect modifying)
  P_genes = 6:10,  # prognostic

  # Effect params
  s_tau  = 5,

  # Prognosis params
  s_eta  = 10,
  mu0    = 120,
  beta_A = 3,
  beta_P = 1,
  sigma_y = 5,

  # Preprocessing
  x_mode = "log1p_z"
)

dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. DATA GENERATION FUNCTION
# ============================================================
generate_data <- function(n, effect_type, seed, cfg) {
  set.seed(seed)

  p <- cfg$p_genes
  K <- cfg$K_latent

  # --- 2.1 Latent factors (gene correlation source) ---
  A      <- rnorm(n, 0, 1)                              # Age
  G      <- matrix(rnorm(n * K), n, K)                  # n x K latent factors
  Lambda <- matrix(rnorm(p * K, 0, cfg$sigma_lambda), p, K)  # p x K loadings

  # Gene baselines and age effects
  b0 <- rnorm(p, log(50), 0.2)           # per-gene baseline
  bA <- rnorm(p, 0.15, 0.05)             # per-gene age effect
  eps_E <- matrix(rnorm(n * p, 0, cfg$sigma_e), n, p)

  # --- 2.2 Latent log expression: n x p ---
  # E_ij = b0_j + bA_j * A_i + sum_k(Lambda_jk * G_ik) + eps_E_ij
  E <- outer(rep(1, n), b0) +
       outer(A, bA) +
       G %*% t(Lambda) +
       eps_E

  # --- 2.3 NB counts (RNA-seq observation) ---
  s   <- exp(rnorm(n, 0, cfg$libsize_sd))            # per-sample depth factor
  mu_mat <- exp(E) * s                                # NB mean: n x p
  C   <- matrix(
    rnbinom(n * p, size = cfg$theta_nb, mu = as.vector(mu_mat)),
    n, p
  )
  colnames(C) <- paste0("gene", seq_len(p))

  # --- 2.4 Observed log ---
  L <- log1p(C)
  colnames(L) <- paste0("gene", seq_len(p))

  # --- 3. Treatment assignment (RCT) ---
  T_vec <- rbinom(n, 1, 0.5)

  # --- 4. True CATE (log mechanism) ---
  M  <- cfg$M_genes
  u  <- rowSums(E[, M, drop = FALSE])
  uc <- u - mean(u)

  if (effect_type == "linear") {
    tau_raw <- uc

  } else if (effect_type == "threshold3") {
    q1 <- quantile(u, 0.33)
    q2 <- quantile(u, 0.66)
    tau_raw <- ifelse(u <= q1, -1, ifelse(u <= q2, 0, 1))

  } else {
    stop(sprintf("Unknown effect_type: %s", effect_type))
  }

  # SD normalization (guard against degenerate case)
  sd_tau <- sd(tau_raw)
  notes_dgp <- ""
  if (sd_tau < 1e-8) {
    warning("SD(tau_raw) < 1e-8; setting tau = 0")
    tau <- rep(0, n)
    notes_dgp <- "tau_raw SD<1e-8; tau=0"
  } else {
    tau <- tau_raw * (cfg$s_tau / sd_tau)
  }

  # --- 5. Prognosis (log mechanism) ---
  P  <- cfg$P_genes
  v  <- rowSums(E[, P, drop = FALSE])
  vc <- v - mean(v)

  eta_raw <- cfg$beta_A * A + cfg$beta_P * vc
  sd_eta  <- sd(eta_raw)
  if (sd_eta < 1e-8) {
    eta <- rep(0, n)
  } else {
    eta <- eta_raw * (cfg$s_eta / sd_eta)
  }

  eps_y <- rnorm(n, 0, cfg$sigma_y)
  Y0    <- cfg$mu0 + eta + eps_y
  Y1    <- Y0 + tau
  Y     <- Y0 + T_vec * tau

  # --- Verification logs ---
  cat(sprintf("    [DGP] tau: mean=%.2f, SD=%.2f (target SD=%.1f)\n",
              mean(tau), sd(tau), cfg$s_tau))
  if (effect_type == "threshold3") {
    tbl <- table(cut(u, breaks = c(-Inf, quantile(u, 0.33),
                                   quantile(u, 0.66), Inf),
                     labels = c("low", "mid", "high")))
    cat(sprintf("    [DGP] threshold3 counts: %s\n",
                paste(names(tbl), tbl, sep = "=", collapse = ", ")))
  }

  list(
    A = A, E = E, C = C, L = L,
    T_vec = T_vec, tau = tau,
    Y0 = Y0, Y1 = Y1, Y = Y,
    notes_dgp = notes_dgp
  )
}

# ============================================================
# 3. STRATIFIED SPLIT
# ============================================================
stratified_split <- function(T_vec, train_frac, seed) {
  set.seed(seed)
  idx0 <- which(T_vec == 0)
  idx1 <- which(T_vec == 1)

  n0_tr <- round(length(idx0) * train_frac)
  n1_tr <- round(length(idx1) * train_frac)

  tr_idx <- sort(c(
    sample(idx0)[seq_len(n0_tr)],
    sample(idx1)[seq_len(n1_tr)]
  ))
  te_idx <- sort(setdiff(seq_along(T_vec), tr_idx))

  list(train = tr_idx, test = te_idx)
}

# ============================================================
# 4. PREPROCESSING (log1p_z, train stats -> train/test)
# ============================================================
preprocess_log1p_z <- function(L, A, train_idx, p_genes = 50L) {
  # Gene-wise z-score: compute stats on train, apply to all
  mu_tr <- colMeans(L[train_idx, , drop = FALSE])
  sd_tr <- apply(L[train_idx, , drop = FALSE], 2, sd)
  sd_tr <- pmax(sd_tr, 1e-6)

  L_scaled <- sweep(L, 2, mu_tr, "-")
  L_scaled <- sweep(L_scaled, 2, sd_tr, "/")
  colnames(L_scaled) <- paste0("gene", seq_len(p_genes))

  # Full X = [A, scaled_genes]: age is always included
  X_full <- cbind(A = A, L_scaled)

  list(
    X_train      = X_full[train_idx,  , drop = FALSE],
    X_test       = X_full[-train_idx, , drop = FALSE],
    covar_names  = colnames(X_full)
  )
}

# ============================================================
# 5. RMSE
# ============================================================
rmse_fn <- function(pred, true) sqrt(mean((pred - true)^2))

# ============================================================
# 6. METHOD WRAPPERS
# Each returns list(tau_hat = numeric_vec_or_NA, notes = string)
# ============================================================

# --- 6.1 Causal Forest (GRF) ---
fit_predict_grf <- function(X_train, T_train, Y_train, X_test, seed) {
  tryCatch({
    covar_names <- colnames(X_train)
    df_train <- data.frame(Y = Y_train, Trt = T_train, X_train,
                           check.names = FALSE)
    fml <- as.formula(
      paste("Y ~ Trt |", paste(covar_names, collapse = " + "))
    )
    fit <- cf_fit(fml, data = df_train, method = "grf",
                  num.trees = 2000L, random_state = as.integer(seed))
    df_test <- as.data.frame(X_test)
    tau_hat <- predict(fit, newdata = df_test, type = "ite")

    list(tau_hat = as.numeric(tau_hat), notes = "OK")
  }, error = function(e) {
    list(tau_hat = NA_real_, notes = conditionMessage(e))
  })
}

# --- 6.2 DR-Learner ---
fit_predict_drlearner <- function(X_train, T_train, Y_train, X_test,
                                  seed, n_train) {
  tryCatch({
    cv <- if (n_train < 200L) 3L else 5L

    covar_names <- colnames(X_train)
    df_train <- data.frame(Y = Y_train, Trt = T_train, X_train,
                           check.names = FALSE)
    fml <- as.formula(
      paste("Y ~ Trt |", paste(covar_names, collapse = " + "))
    )
    fit <- cf_fit(fml, data = df_train, method = "drlearner",
                  random_state = as.integer(seed), cv = cv)
    df_test <- as.data.frame(X_test)
    tau_hat <- predict(fit, newdata = df_test, type = "ite")

    list(tau_hat = as.numeric(tau_hat),
         notes = sprintf("OK (cv=%d)", cv))
  }, error = function(e) {
    list(tau_hat = NA_real_, notes = conditionMessage(e))
  })
}

# --- 6.3 GANITE (direct Python GANITEModel) ---
fit_predict_ganite <- function(X_train, T_train, Y_train, X_test,
                               seed, n_train) {
  tryCatch({
    # Import Python module directly (bypassing run_ganite wrapper
    # so we can call model.predict(X_test))
    ganite_mod <- reticulate::import_from_path(
      "ganite",
      path = system.file("python", package = "cfomics")
    )
    np <- reticulate::import("numpy", convert = FALSE)

    # Convert to float32 numpy arrays
    X_tr_np  <- np$array(X_train,            dtype = "float32")
    T_tr_np  <- np$array(as.numeric(T_train), dtype = "float32")
    Y_tr_np  <- np$array(as.numeric(Y_train), dtype = "float32")
    X_te_np  <- np$array(X_test,             dtype = "float32")

    # Conservative hyperparameters for CPU
    h_dim    <- 64L
    batch_sz <- as.integer(max(16L, min(n_train %/% 4L, 64L)))
    iters    <- if (n_train < 200L) 500L else 1000L

    model <- ganite_mod$GANITEModel(
      input_dim    = as.integer(ncol(X_train)),
      h_dim        = h_dim,
      alpha        = 1.0,
      random_state = as.integer(seed)
    )
    model$fit(
      X_tr_np, T_tr_np, Y_tr_np,
      iterations = iters,
      batch_size = batch_sz,
      verbose    = FALSE
    )

    # predict() returns (n_test, 2): col0=Y0, col1=Y1
    y_hat    <- reticulate::py_to_r(model$predict(X_te_np))
    tau_hat  <- y_hat[, 2] - y_hat[, 1]  # Y1 - Y0

    list(tau_hat = as.numeric(tau_hat),
         notes = sprintf("OK (iters=%d, batch=%d, h_dim=%d)",
                         iters, batch_sz, h_dim))
  }, error = function(e) {
    list(tau_hat = NA_real_, notes = conditionMessage(e))
  })
}

# --- 6.4 CEVAE (direct pyro.contrib.cevae) ---
fit_predict_cevae <- function(X_train, T_train, Y_train, X_test,
                              seed, n_train) {
  tryCatch({
    torch     <- reticulate::import("torch",  convert = FALSE)
    pyro      <- reticulate::import("pyro",   convert = FALSE)
    cevae_mod <- reticulate::import("pyro.contrib.cevae", convert = FALSE)
    CEVAE     <- cevae_mod$CEVAE

    pyro$set_rng_seed(as.integer(seed))

    # Conservative hyperparameters for CPU
    latent_dim <- 5L
    hidden_dim <- 64L
    num_layers <- 2L
    num_epochs <- 20L
    batch_sz   <- as.integer(min(n_train, 64L))

    x_tr <- torch$tensor(X_train,            dtype = torch$float32)
    t_tr <- torch$tensor(as.numeric(T_train), dtype = torch$float32)
    y_tr <- torch$tensor(as.numeric(Y_train), dtype = torch$float32)
    x_te <- torch$tensor(X_test,             dtype = torch$float32)

    cevae <- CEVAE(
      feature_dim  = as.integer(ncol(X_train)),
      outcome_dist = "normal",
      latent_dim   = latent_dim,
      hidden_dim   = hidden_dim,
      num_layers   = num_layers,
      num_samples  = 10L
    )
    cevae$fit(
      x_tr, t_tr, y_tr,
      num_epochs          = num_epochs,
      batch_size          = batch_sz,
      learning_rate       = 1e-3,
      learning_rate_decay = 0.1,
      weight_decay        = 1e-4
    )

    ite_py  <- cevae$ite(x_te)
    ite_np  <- ite_py$detach()$cpu()$numpy()
    tau_hat <- as.numeric(ite_np)

    list(tau_hat = tau_hat,
         notes = sprintf("OK (epochs=%d, batch=%d, latent=%d)",
                         num_epochs, batch_sz, latent_dim))
  }, error = function(e) {
    list(tau_hat = NA_real_, notes = conditionMessage(e))
  })
}

# ============================================================
# 7. MAIN EXPERIMENT LOOP
# ============================================================
results_list <- list()
pred_list    <- list()

for (idx_n in seq_along(config$n_values)) {
  n <- config$n_values[idx_n]

  for (idx_effect in seq_along(config$effect_types)) {
    effect_type <- config$effect_types[idx_effect]

    seed_cond <- config$seed_base +
      10000L * (idx_n - 1L) +
      100L   * (idx_effect - 1L)

    cat(sprintf("\n=== n=%d, effect=%s, seed=%d ===\n",
                n, effect_type, seed_cond))

    # --- Generate data (shared across all methods) ---
    cat("  Generating data...\n")
    dat <- generate_data(n, effect_type, seed_cond, config)

    # --- Stratified split (shared across all methods) ---
    spl <- stratified_split(dat$T_vec, config$train_frac, seed_cond)

    # --- Preprocess (log1p_z, train stats → train/test) ---
    pre <- preprocess_log1p_z(dat$L, dat$A, spl$train, config$p_genes)

    X_train       <- pre$X_train
    X_test        <- pre$X_test
    T_train       <- dat$T_vec[spl$train]
    Y_train       <- dat$Y[spl$train]
    tau_test_true <- dat$tau[spl$test]
    n_train       <- length(spl$train)
    n_test        <- length(spl$test)

    cat(sprintf("  n_train=%d, n_test=%d\n", n_train, n_test))

    # --- Dispatch table ---
    method_fns <- list(
      causal_forest = function()
        fit_predict_grf(X_train, T_train, Y_train, X_test, seed_cond),
      drlearner     = function()
        fit_predict_drlearner(X_train, T_train, Y_train, X_test,
                              seed_cond, n_train),
      ganite        = function()
        fit_predict_ganite(X_train, T_train, Y_train, X_test,
                           seed_cond, n_train),
      cevae         = function()
        fit_predict_cevae(X_train, T_train, Y_train, X_test,
                          seed_cond, n_train)
    )

    for (method in config$methods) {
      cat(sprintf("  [%s] fitting...\n", method))
      t0 <- proc.time()["elapsed"]

      res_m <- method_fns[[method]]()

      elapsed <- proc.time()["elapsed"] - t0
      cat(sprintf("  [%s] %.1fs | %s\n", method, elapsed, res_m$notes))

      # RMSE (on test only)
      tau_hat_vec <- res_m$tau_hat
      if (length(tau_hat_vec) == n_test && !anyNA(tau_hat_vec)) {
        rmse_val <- rmse_fn(tau_hat_vec, tau_test_true)
        cat(sprintf("  [%s] RMSE_tau = %.4f\n", method, rmse_val))
      } else {
        tau_hat_vec <- rep(NA_real_, n_test)
        rmse_val    <- NA_real_
      }

      notes_full <- paste0(
        if (nchar(dat$notes_dgp) > 0) paste0(dat$notes_dgp, "; ") else "",
        res_m$notes
      )

      # Append summary row
      results_list[[length(results_list) + 1]] <- data.frame(
        seed_base      = config$seed_base,
        seed_condition = seed_cond,
        n              = n,
        effect_type    = effect_type,
        x_mode         = config$x_mode,
        method         = method,
        rmse_tau       = rmse_val,
        n_train        = n_train,
        n_test         = n_test,
        notes          = notes_full,
        stringsAsFactors = FALSE
      )

      # Append prediction rows
      pred_list[[length(pred_list) + 1]] <- data.frame(
        seed_base      = config$seed_base,
        seed_condition = seed_cond,
        n              = n,
        effect_type    = effect_type,
        x_mode         = config$x_mode,
        method         = method,
        id             = spl$test,
        split          = "test",
        tau_true       = tau_test_true,
        tau_hat        = tau_hat_vec,
        stringsAsFactors = FALSE
      )
    }  # method loop
  }  # effect_type loop
}  # n loop

# ============================================================
# 8. OUTPUT
# ============================================================
cat("\n=== Saving outputs ===\n")

# --- 8.1 Summary CSV ---
results_df <- do.call(rbind, results_list)
csv_summary <- file.path(config$out_dir, "results_demo_omics50_v02.csv")
write.csv(results_df, csv_summary, row.names = FALSE)
cat(sprintf("  Saved: %s (%d rows)\n", csv_summary, nrow(results_df)))

# --- 8.2 Predictions CSV ---
pred_df <- do.call(rbind, pred_list)
csv_pred <- file.path(config$out_dir, "predictions_demo_omics50_v02.csv")
write.csv(pred_df, csv_pred, row.names = FALSE)
cat(sprintf("  Saved: %s (%d rows)\n", csv_pred, nrow(pred_df)))

# --- 8.3 Plots ---
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  method_colors <- c(
    causal_forest = "#E41A1C",
    drlearner     = "#377EB8",
    ganite        = "#4DAF4A",
    cevae         = "#FF7F00"
  )

  for (et in config$effect_types) {
    sub <- subset(results_df, effect_type == et & !is.na(rmse_tau))
    sub$n      <- factor(sub$n, levels = config$n_values)
    sub$method <- factor(sub$method, levels = config$methods)

    p <- ggplot(sub, aes(x = n, y = rmse_tau,
                         color = method, group = method)) +
      geom_line(linewidth = 1.0, na.rm = TRUE) +
      geom_point(size = 2.5, na.rm = TRUE) +
      scale_color_manual(values = method_colors) +
      labs(
        title = sprintf("CATE RMSE by sample size (%s)", et),
        subtitle = sprintf(
          "s_tau=%.0f, K=%d, theta=%d, seed=%d",
          config$s_tau, config$K_latent,
          config$theta_nb, config$seed_base
        ),
        x = "n (sample size)",
        y = "CATE RMSE (test set)",
        color = "Method"
      ) +
      theme_bw(base_size = 13) +
      theme(legend.position = "bottom")

    # Add NA-method annotations if any were skipped
    na_methods <- setdiff(config$methods,
                          unique(sub$method))
    if (length(na_methods) > 0) {
      p <- p + labs(caption = paste("NA:", paste(na_methods, collapse = ", ")))
    }

    png_path <- file.path(
      config$out_dir,
      sprintf("plot_rmse_by_n_%s.png", et)
    )
    ggsave(png_path, p, width = 7, height = 5, dpi = 150)
    cat(sprintf("  Saved: %s\n", png_path))
  }
} else {
  cat("  ggplot2 not available; skipping plots\n")
}

# --- 8.4 Session info ---
session_path <- file.path(config$out_dir, "session_info.txt")
sink(session_path)
cat("=== R Session Info ===\n")
print(sessionInfo())
cat("\n=== Python Environment ===\n")
cat(sprintf("RETICULATE_PYTHON: %s\n", Sys.getenv("RETICULATE_PYTHON")))
tryCatch({
  py_cfg <- reticulate::py_config()
  cat(sprintf("Python: %s\n", py_cfg$version))
  cat(sprintf("numpy: %s\n", reticulate::py_module_available("numpy")))
  cat(sprintf("torch: %s\n", reticulate::py_module_available("torch")))
  cat(sprintf("pyro:  %s\n", reticulate::py_module_available("pyro")))
  cat(sprintf("econml: %s\n", reticulate::py_module_available("econml")))
  cat(sprintf("tensorflow: %s\n", reticulate::py_module_available("tensorflow")))
}, error = function(e) cat("Python config error:", conditionMessage(e), "\n"))
cat("\n=== cfomics git hash ===\n")
tryCatch({
  hash <- system("git -C packages/cfomics rev-parse HEAD 2>/dev/null ||
                   git rev-parse HEAD 2>/dev/null", intern = TRUE)
  cat(hash[1], "\n")
}, error = function(e) cat("git hash unavailable\n"))
sink()
cat(sprintf("  Saved: %s\n", session_path))

# --- 8.5 Print summary tables ---
cat("\n")
cat("=================================================================\n")
cat("                     RESULTS SUMMARY\n")
cat("=================================================================\n")

for (et in config$effect_types) {
  cat(sprintf("\n--- effect_type: %s ---\n", et))
  sub <- subset(results_df, effect_type == et)

  # Wide table: rows = method, cols = n
  tbl <- reshape(
    sub[, c("method", "n", "rmse_tau")],
    idvar    = "method",
    timevar  = "n",
    direction = "wide"
  )
  colnames(tbl) <- gsub("rmse_tau\\.", "n=", colnames(tbl))
  row.names(tbl) <- NULL

  # Round for display
  for (col in grep("^n=", colnames(tbl), value = TRUE)) {
    tbl[[col]] <- round(tbl[[col]], 4)
  }

  print(tbl, row.names = FALSE)
}

cat(sprintf("\n=== Finished: %s ===\n", Sys.time()))
cat("=================================================================\n")
