#!/usr/bin/env Rscript
#
# CATE Demo: Single Biomarker Effect Modification (v0.1)
#
# Evaluates 5 causal inference methods on CATE estimation where the true
# effect modifier is a single omics feature (gene 1).
#
# DGP: 30 NB-count genes + age, RCT design, continuous outcome (blood pressure).
# Effect shapes: linear and threshold (responder/non-responder).
# Methods: Causal Forest (GRF), BCF, DR-Learner, GANITE, CAVAE.
# Evaluation: RMSE of CATE on test set.
#
# Usage:
#   cd /path/to/cfomics
#   Rscript packages/cfomics/inst/benchmarks/run_cate_demo_v01.R
#
# Output:
#   benchmark_results/cate_demo_v01/results_demo_single_biomarker_v0.1.csv
#   benchmark_results/cate_demo_v01/predictions_demo_single_biomarker_v0.1.csv
#   benchmark_results/cate_demo_v01/session_info.txt
#

# ==============================================================================
# Section 0: Experiment Configuration
# ==============================================================================

CONFIG <- list(
  seed        = 20260215L,
  n_values    = c(100L, 500L),
  effect_types = c("linear", "threshold"),
  x_modes     = c("count_raw", "log1p_z"),
  methods     = c("grf", "bcf", "drlearner", "ganite", "cavae"),
  train_frac  = 0.8,
  out_dir     = file.path("benchmark_results", "cate_demo_v01")
)

# Global: path to ganite.py (resolved at runtime in main())
GANITE_PY_PATH <- NULL

# DGP fixed parameters (from spec v0.1)
DGP_PARAMS <- list(
  n_genes   = 30L,
  theta     = 10,       # NB dispersion
  b0_mean   = log(50),  # intercept mean for log(mu)
  b0_sd     = 0.2,      # intercept SD
  bA_mean   = 0.15,     # age->omics coefficient mean
  bA_sd     = 0.05,     # age->omics coefficient SD
  m         = 4.0,      # centering for log1p (log1p(50) ~ 3.93)
  s         = 1.0,      # scaling denominator
  s_tau     = 5,        # target SD of tau
  s_eta     = 10,       # target SD of prognostic score
  beta_A    = 3,        # age coefficient for prognostic
  sigma     = 5,        # noise SD
  mu_0      = 120,      # baseline blood pressure
  effect_j  = 1L,       # effect modifier gene index (j*)
  prog_idx  = 11:20,    # prognostic gene indices
  noise_idx = 21:30,    # noise gene indices
  w_linear  = 1,        # weight for linear effect
  c_thresh  = 0         # threshold cutoff for u
)

# ==============================================================================
# Section 1: DGP — Data Generating Process
# ==============================================================================

#' Generate CATE demo data following spec v0.1
#'
#' @param n Integer, sample size
#' @param effect_type Character, "linear" or "threshold"
#' @param seed Integer, random seed
#' @return List with A, C, L, T_vec, Y, Y0, Y1, tau_true, and metadata
generate_cate_demo_data <- function(n, effect_type = c("linear", "threshold"),
                                    seed = CONFIG$seed) {
  effect_type <- match.arg(effect_type)
  P <- DGP_PARAMS

  set.seed(seed)

  # Step 1: Age
  A <- rnorm(n, mean = 0, sd = 1)

  # Step 2: NB counts (30 genes)
  # Gene-specific parameters (deterministic given seed)
  b0 <- rnorm(P$n_genes, mean = P$b0_mean, sd = P$b0_sd)
  bA <- rnorm(P$n_genes, mean = P$bA_mean, sd = P$bA_sd)

  C <- matrix(NA_integer_, nrow = n, ncol = P$n_genes)
  for (j in seq_len(P$n_genes)) {
    log_mu <- b0[j] + bA[j] * A
    mu <- exp(log_mu)
    C[, j] <- rnbinom(n, size = P$theta, mu = mu)
  }
  colnames(C) <- paste0("gene", seq_len(P$n_genes))

  # Step 3: Internal log1p representation (used for effect/prognostic definition)
  L <- log1p(C)

  # Step 4: True CATE — effect modification by gene j* = 1 only
  u <- (L[, P$effect_j] - P$m) / P$s

  if (effect_type == "linear") {
    tau_raw <- P$w_linear * u
  } else {
    tau_raw <- as.numeric(u > P$c_thresh)
  }

  # Scale to target SD
  sd_tau_raw <- sd(tau_raw)
  if (sd_tau_raw < 1e-10) {
    # All-zero or constant: no heterogeneity
    tau <- rep(0, n)
    warning("tau_raw has zero variance; setting tau = 0")
  } else {
    tau <- tau_raw * (P$s_tau / sd_tau_raw)
  }

  # Step 5: Baseline Y(0) — prognostic from age + genes 11-20
  Z_P <- (L[, P$prog_idx] - P$m) / P$s

  # Prognostic coefficients: +1/-1 half each (deterministic given seed)
  set.seed(seed + 1L)
  beta_P <- sample(rep(c(1, -1), length.out = length(P$prog_idx)))

  eta_raw <- P$beta_A * A + as.numeric(Z_P %*% beta_P)

  sd_eta_raw <- sd(eta_raw)
  if (sd_eta_raw < 1e-10) {
    eta <- rep(0, n)
  } else {
    eta <- eta_raw * (P$s_eta / sd_eta_raw)
  }

  epsilon <- rnorm(n, mean = 0, sd = P$sigma)
  Y0 <- P$mu_0 + eta + epsilon

  # Step 6: Treatment assignment (RCT) and observed outcome
  set.seed(seed + 2L)
  T_vec <- rbinom(n, size = 1, prob = 0.5)

  Y1 <- Y0 + tau
  Y <- T_vec * Y1 + (1 - T_vec) * Y0

  list(
    A        = A,
    C        = C,
    L        = L,
    T_vec    = T_vec,
    Y        = Y,
    Y0       = Y0,
    Y1       = Y1,
    tau_true = tau,
    n        = n,
    effect_type = effect_type,
    seed     = seed,
    b0       = b0,
    bA       = bA,
    beta_P   = beta_P,
    ate_true = mean(tau)
  )
}

# ==============================================================================
# Section 2: Preprocessing
# ==============================================================================

#' Preprocess covariates according to x_mode
#'
#' @param A Numeric vector of age
#' @param C Integer matrix of counts (n x 30)
#' @param x_mode Character, "count_raw" or "log1p_z"
#' @param train_idx Integer vector of training indices (for computing stats)
#' @return List with X_full (n x 31 matrix) and stats
preprocess_covariates <- function(A, C, x_mode = c("count_raw", "log1p_z"),
                                  train_idx) {
  x_mode <- match.arg(x_mode)
  n <- nrow(C)
  p <- ncol(C)

  if (x_mode == "count_raw") {
    # Clip each gene at 99th percentile (computed from train)
    q99 <- apply(C[train_idx, , drop = FALSE], 2, quantile, probs = 0.99)
    C_clip <- C
    for (j in seq_len(p)) {
      C_clip[, j] <- pmin(C[, j], q99[j])
    }
    X_omics <- as.matrix(C_clip)
  } else {
    # log1p + z-score standardization (fit on train)
    L <- log1p(C)
    L_train <- L[train_idx, , drop = FALSE]
    means <- colMeans(L_train)
    sds <- apply(L_train, 2, sd)
    sds[sds < 1e-10] <- 1  # guard against zero variance

    X_omics <- sweep(L, 2, means, "-")
    X_omics <- sweep(X_omics, 2, sds, "/")
  }

  X_full <- cbind(age = A, X_omics)
  colnames(X_full) <- c("age", paste0("gene", seq_len(p)))
  X_full
}

# ==============================================================================
# Section 3: Stratified Train/Test Split
# ==============================================================================

#' Stratified train/test split by treatment
#'
#' @param T_vec Binary treatment vector
#' @param ratio Numeric, fraction for training
#' @param seed Integer, random seed for splitting
#' @return List with train and test index vectors
stratified_split <- function(T_vec, ratio = 0.8, seed = CONFIG$seed) {
  set.seed(seed + 100L)

  idx_t0 <- which(T_vec == 0)
  idx_t1 <- which(T_vec == 1)

  n_train_t0 <- round(length(idx_t0) * ratio)
  n_train_t1 <- round(length(idx_t1) * ratio)

  train_t0 <- sample(idx_t0, n_train_t0)
  train_t1 <- sample(idx_t1, n_train_t1)

  train_idx <- sort(c(train_t0, train_t1))
  test_idx  <- setdiff(seq_along(T_vec), train_idx)

  list(train = train_idx, test = test_idx)
}

# ==============================================================================
# Section 4: Method Availability Check
# ==============================================================================

#' Check if a method's dependencies are available
#'
#' @param method Character, method name
#' @return NULL if available, character error message if not
check_method_available <- function(method) {
  switch(method,
    "grf" = {
      if (!requireNamespace("grf", quietly = TRUE))
        return("grf package not installed")
      NULL
    },
    "bcf" = {
      if (!requireNamespace("bcf", quietly = TRUE))
        return("bcf package not installed")
      NULL
    },
    "drlearner" = {
      if (!reticulate::py_available(initialize = FALSE))
        return("Python not available")
      tryCatch({
        reticulate::py_available(initialize = TRUE)
        if (!reticulate::py_module_available("econml"))
          return("econml Python module not available")
        NULL
      }, error = function(e) conditionMessage(e))
    },
    "ganite" = {
      if (!reticulate::py_available(initialize = FALSE))
        return("Python not available")
      tryCatch({
        reticulate::py_available(initialize = TRUE)
        if (!reticulate::py_module_available("tensorflow"))
          return("tensorflow Python module not available")
        NULL
      }, error = function(e) conditionMessage(e))
    },
    "cavae" = {
      if (!reticulate::py_available(initialize = FALSE))
        return("Python not available")
      tryCatch({
        reticulate::py_available(initialize = TRUE)
        if (!reticulate::py_module_available("torch"))
          return("torch Python module not available")
        if (!reticulate::py_module_available("pyro"))
          return("pyro Python module not available")
        NULL
      }, error = function(e) conditionMessage(e))
    },
    sprintf("Unknown method: %s", method)
  )
}

# ==============================================================================
# Section 5: Method-Specific Fit + Predict Wrappers
# ==============================================================================

# Each function: fit on train, predict tau_hat on test
# Returns list(tau_hat_test = numeric, notes = character)

# --- 5a: GRF (Causal Forest) -------------------------------------------------

fit_predict_grf <- function(X_train, T_train, Y_train, X_test, seed) {
  forest <- grf::causal_forest(
    X           = X_train,
    Y           = Y_train,
    W           = T_train,
    W.hat       = rep(0.5, length(T_train)),  # RCT
    num.trees   = 2000L,
    seed        = seed
  )

  tau_hat_test <- as.numeric(predict(forest, newdata = X_test)$predictions)

  list(tau_hat_test = tau_hat_test, notes = NA_character_)
}

# --- 5b: BCF (Bayesian Causal Forests) ---------------------------------------
# BCF does not support out-of-sample prediction.
# Fallback: fit on full data, extract test indices.

fit_predict_bcf <- function(X_full, T_full, Y_full, test_idx, seed) {
  if (!is.null(seed)) set.seed(seed)

  pihat <- rep(0.5, length(T_full))  # RCT

  bcf_fit <- bcf::bcf(
    y          = Y_full,
    z          = T_full,
    x_control  = X_full,
    x_moderate = X_full,
    pihat      = pihat,
    nburn      = 1000L,
    nsim       = 2000L
  )

  tau_samples <- bcf_fit$tau
  ite_all     <- colMeans(tau_samples)
  tau_hat_test <- as.numeric(ite_all[test_idx])

  list(
    tau_hat_test = tau_hat_test,
    notes = "BCF: fitted on full data (no OOS prediction); test tau is in-sample"
  )
}

# --- 5c: DR-Learner (EconML) -------------------------------------------------

fit_predict_drlearner <- function(X_train, T_train, Y_train, X_test, seed) {
  np <- reticulate::import("numpy")
  DRLearner <- reticulate::import("econml.dr")$DRLearner

  if (!is.null(seed)) np$random$seed(as.integer(seed))

  X_train_py <- reticulate::r_to_py(X_train)
  T_train_py <- reticulate::r_to_py(as.numeric(T_train))
  Y_train_py <- reticulate::r_to_py(as.numeric(Y_train))
  X_test_py  <- reticulate::r_to_py(X_test)

  n_train <- nrow(X_train)
  cv_folds <- if (n_train < 100) 3L else 5L

  dr <- DRLearner(cv = cv_folds)
  dr$fit(Y_train_py, T_train_py, X = X_train_py)

  tau_hat_py   <- dr$effect(X_test_py, T0 = 0L, T1 = 1L)
  tau_hat_test <- as.numeric(reticulate::py_to_r(tau_hat_py))

  list(tau_hat_test = tau_hat_test, notes = NA_character_)
}

# --- 5d: GANITE (GAN-based ITE) ----------------------------------------------
# Use GANITEModel class directly for proper train/test.

fit_predict_ganite <- function(X_train, T_train, Y_train, X_test, seed) {
  # Resolve ganite.py path (works without cfomics loaded)
  ganite_path <- GANITE_PY_PATH
  if (is.null(ganite_path) || !dir.exists(ganite_path)) {
    stop("Cannot find ganite.py. Set GANITE_PY_PATH before running.")
  }
  ganite_mod <- reticulate::import_from_path(
    module = "ganite",
    path   = ganite_path
  )
  np <- reticulate::import("numpy")

  X_train_np <- np$array(X_train, dtype = "float64")
  T_train_np <- np$array(as.numeric(T_train), dtype = "float64")
  Y_train_np <- np$array(as.numeric(Y_train), dtype = "float64")
  X_test_np  <- np$array(X_test, dtype = "float64")

  n_train    <- nrow(X_train)
  batch_size <- min(256L, as.integer(n_train))

  model <- ganite_mod$GANITEModel(
    input_dim    = as.integer(ncol(X_train)),
    h_dim        = 100L,
    alpha        = 1.0,
    random_state = as.integer(seed)
  )

  model$fit(
    X_train_np, T_train_np, Y_train_np,
    iterations = 5000L,
    batch_size = batch_size,
    verbose    = FALSE
  )

  # predict returns (n_test, 2) with [:, 0]=y0, [:, 1]=y1
  y_hat_test <- reticulate::py_to_r(model$predict(X_test_np))
  tau_hat_test <- as.numeric(y_hat_test[, 2] - y_hat_test[, 1])

  list(tau_hat_test = tau_hat_test, notes = NA_character_)
}

# --- 5e: CAVAE (Causal Effect VAE) -------------------------------------------

fit_predict_cavae <- function(X_train, T_train, Y_train, X_test, seed) {
  torch     <- reticulate::import("torch", convert = FALSE)
  pyro      <- reticulate::import("pyro", convert = FALSE)
  cevae_mod <- reticulate::import("pyro.contrib.cevae", convert = FALSE)
  CEVAE     <- cevae_mod$CEVAE

  pyro$set_rng_seed(as.integer(seed))

  n_train     <- nrow(X_train)
  feature_dim <- ncol(X_train)
  batch_size  <- min(100L, as.integer(n_train))

  x_train_t <- torch$tensor(X_train, dtype = torch$float32)
  t_train_t <- torch$tensor(as.numeric(T_train), dtype = torch$float32)
  y_train_t <- torch$tensor(as.numeric(Y_train), dtype = torch$float32)
  x_test_t  <- torch$tensor(X_test, dtype = torch$float32)

  cevae <- CEVAE(
    feature_dim  = as.integer(feature_dim),
    outcome_dist = "normal",  # continuous outcome
    latent_dim   = 20L,
    hidden_dim   = 200L,
    num_layers   = 3L,
    num_samples  = 10L
  )

  cevae$fit(
    x_train_t, t_train_t, y_train_t,
    num_epochs          = 30L,
    batch_size          = batch_size,
    learning_rate       = 1e-3,
    learning_rate_decay = 0.1,
    weight_decay        = 1e-4
  )

  ite_py <- cevae$ite(x_test_t)
  ite_np <- ite_py$detach()$cpu()$numpy()
  tau_hat_test <- as.numeric(ite_np)

  list(tau_hat_test = tau_hat_test, notes = NA_character_)
}

# ==============================================================================
# Section 6: Evaluation
# ==============================================================================

compute_rmse <- function(tau_hat, tau_true) {
  sqrt(mean((tau_hat - tau_true)^2))
}

# ==============================================================================
# Section 7: Main Experiment Loop
# ==============================================================================

run_experiment <- function(config = CONFIG) {
  results_list     <- list()
  predictions_list <- list()
  run_idx <- 0L

  # Pre-check method availability
  method_status <- setNames(
    lapply(config$methods, check_method_available),
    config$methods
  )

  available_methods <- config$methods[sapply(method_status, is.null)]
  skipped_methods   <- config$methods[!sapply(method_status, is.null)]

  if (length(skipped_methods) > 0) {
    message("\n--- Skipped methods (missing dependencies) ---")
    for (m in skipped_methods) {
      message(sprintf("  %s: %s", m, method_status[[m]]))
    }
    message("")
  }

  total_runs <- length(config$n_values) * length(config$effect_types) *
    length(config$x_modes) * length(config$methods)
  message(sprintf("Total planned runs: %d (available: %d, skipped: %d per condition)",
                  total_runs,
                  total_runs - length(skipped_methods) *
                    length(config$n_values) * length(config$effect_types) *
                    length(config$x_modes),
                  length(skipped_methods)))

  for (n_val in config$n_values) {
    for (eff in config$effect_types) {
      # Generate data once per (n, effect_type) — shared across x_modes
      message(sprintf("\n=== Generating data: n=%d, effect=%s ===", n_val, eff))
      dgp <- generate_cate_demo_data(n = n_val, effect_type = eff,
                                     seed = config$seed)

      # Split once per (n, effect_type) — shared across x_modes and methods
      split <- stratified_split(dgp$T_vec, ratio = config$train_frac,
                                seed = config$seed)

      for (xm in config$x_modes) {
        # Preprocess covariates
        X_full <- preprocess_covariates(dgp$A, dgp$C, x_mode = xm,
                                        train_idx = split$train)
        X_train <- X_full[split$train, , drop = FALSE]
        X_test  <- X_full[split$test, , drop = FALSE]
        T_train <- dgp$T_vec[split$train]
        Y_train <- dgp$Y[split$train]
        tau_test <- dgp$tau_true[split$test]

        for (method in config$methods) {
          run_idx <- run_idx + 1L

          # Check availability
          if (method %in% skipped_methods) {
            message(sprintf("  [%d/%d] n=%d, %s, %s, %s -> SKIPPED",
                            run_idx, total_runs, n_val, eff, xm, method))

            results_list[[run_idx]] <- data.frame(
              seed        = config$seed,
              n           = n_val,
              effect_type = eff,
              x_mode      = xm,
              method      = method,
              rmse_tau    = NA_real_,
              n_train     = length(split$train),
              n_test      = length(split$test),
              notes       = paste("SKIPPED:", method_status[[method]]),
              stringsAsFactors = FALSE
            )

            predictions_list[[run_idx]] <- data.frame(
              seed        = config$seed,
              n           = n_val,
              effect_type = eff,
              x_mode      = xm,
              method      = method,
              id          = split$test,
              split       = "test",
              tau_true    = tau_test,
              tau_hat     = NA_real_,
              stringsAsFactors = FALSE
            )
            next
          }

          message(sprintf("  [%d/%d] n=%d, %s, %s, %s ...",
                          run_idx, total_runs, n_val, eff, xm, method))

          run_result <- tryCatch({
            t0 <- proc.time()

            res <- switch(method,
              "grf" = fit_predict_grf(
                X_train, T_train, Y_train, X_test, seed = config$seed
              ),
              "bcf" = fit_predict_bcf(
                X_full, dgp$T_vec, dgp$Y, split$test, seed = config$seed
              ),
              "drlearner" = fit_predict_drlearner(
                X_train, T_train, Y_train, X_test, seed = config$seed
              ),
              "ganite" = fit_predict_ganite(
                X_train, T_train, Y_train, X_test, seed = config$seed
              ),
              "cavae" = fit_predict_cavae(
                X_train, T_train, Y_train, X_test, seed = config$seed
              )
            )

            elapsed <- (proc.time() - t0)["elapsed"]
            rmse <- compute_rmse(res$tau_hat_test, tau_test)

            message(sprintf("    -> RMSE = %.4f (%.1f sec)", rmse, elapsed))

            list(
              status       = "ok",
              tau_hat_test = res$tau_hat_test,
              rmse_tau     = rmse,
              time_sec     = as.numeric(elapsed),
              notes        = res$notes
            )
          }, error = function(e) {
            message(sprintf("    -> ERROR: %s", conditionMessage(e)))
            list(
              status       = "error",
              tau_hat_test = rep(NA_real_, length(split$test)),
              rmse_tau     = NA_real_,
              time_sec     = NA_real_,
              notes        = paste("ERROR:", conditionMessage(e))
            )
          })

          # Append summary row
          results_list[[run_idx]] <- data.frame(
            seed        = config$seed,
            n           = n_val,
            effect_type = eff,
            x_mode      = xm,
            method      = method,
            rmse_tau    = run_result$rmse_tau,
            n_train     = length(split$train),
            n_test      = length(split$test),
            notes       = ifelse(is.na(run_result$notes), "", run_result$notes),
            stringsAsFactors = FALSE
          )

          # Append per-sample predictions
          predictions_list[[run_idx]] <- data.frame(
            seed        = config$seed,
            n           = n_val,
            effect_type = eff,
            x_mode      = xm,
            method      = method,
            id          = split$test,
            split       = "test",
            tau_true    = tau_test,
            tau_hat     = run_result$tau_hat_test,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  list(
    summary     = do.call(rbind, results_list),
    predictions = do.call(rbind, predictions_list)
  )
}

# ==============================================================================
# Section 8: Entry Point
# ==============================================================================

main <- function() {
  message("============================================================")
  message("  cfomics CATE Demo v0.1 — Single Biomarker Effect Modifier")
  message("============================================================")
  message(sprintf("Date: %s", Sys.time()))
  message(sprintf("Seed: %d", CONFIG$seed))
  message(sprintf("Output: %s", CONFIG$out_dir))
  message("")

  # --- Resolve ganite.py path (no cfomics package load needed) ---
  py_path <- NULL
  candidates <- c(
    file.path(getwd(), "packages", "cfomics", "inst", "python"),
    file.path(getwd(), "cfomics", "inst", "python"),
    file.path(dirname(dirname(getwd())), "cfomics", "inst", "python")
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "ganite.py"))) {
      py_path <- cand
      break
    }
  }
  GANITE_PY_PATH <<- py_path
  if (!is.null(py_path)) {
    message(sprintf("Python scripts: %s", py_path))
  } else {
    message("WARNING: ganite.py not found — GANITE method will fail")
  }

  # --- Create output directory ---
  dir.create(CONFIG$out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Run experiment ---
  results <- run_experiment(CONFIG)

  # --- Write CSVs ---
  summary_path <- file.path(CONFIG$out_dir,
                            "results_demo_single_biomarker_v0.1.csv")
  write.csv(results$summary, summary_path, row.names = FALSE)
  message(sprintf("\nSummary saved: %s", summary_path))

  pred_path <- file.path(CONFIG$out_dir,
                         "predictions_demo_single_biomarker_v0.1.csv")
  write.csv(results$predictions, pred_path, row.names = FALSE)
  message(sprintf("Predictions saved: %s", pred_path))

  # --- Print summary table ---
  message("\n===== Results Summary =====")
  summary_df <- results$summary[, c("n", "effect_type", "x_mode",
                                     "method", "rmse_tau", "n_train",
                                     "n_test", "notes")]
  print(summary_df, row.names = FALSE)

  # --- Save session info ---
  si_path <- file.path(CONFIG$out_dir, "session_info.txt")
  sink(si_path)
  cat("CATE Demo v0.1 — Session Info\n")
  cat(sprintf("Timestamp: %s\n", Sys.time()))
  cat(sprintf("Seed: %d\n\n", CONFIG$seed))
  print(sessionInfo())
  sink()
  message(sprintf("\nSession info saved: %s", si_path))

  message("\nDone.")
  invisible(results)
}

# Run if executed as a script (not sourced)
if (!interactive() || identical(Sys.getenv("RUN_CATE_DEMO"), "1")) {
  main()
}
