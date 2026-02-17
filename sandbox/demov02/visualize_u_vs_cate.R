#!/usr/bin/env Rscript
# visualize_u_vs_cate.R
# 潜在効果修飾スコア u vs CATE（真値・予測値）の散布図
# X軸: u = rowSums(E[, M_genes])  (DGP内部の効果修飾変数)
# Y軸: τ_true（青）、τ̂（赤）
# モデルごとに1 figure（6条件パネル）
# Run from repo root: Rscript sandbox/visualize_u_vs_cate.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/benchmark_results/cate_demo_v02"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. DGPパラメータ（run_cate_demo_v02.R と完全同一）
# ============================================================
cfg <- list(
  seed_base    = 20260215L,
  n_values     = c(100L, 500L, 1000L),
  effect_types = c("linear", "threshold3"),
  p_genes      = 50L,
  K_latent     = 5L,
  sigma_lambda = 0.3,
  sigma_e      = 0.3,
  theta_nb     = 10,
  libsize_sd   = 0.2,
  M_genes      = 1:5,
  P_genes      = 6:10,
  s_tau        = 5,
  s_eta        = 10,
  mu0          = 120,
  beta_A       = 3,
  beta_P       = 1,
  sigma_y      = 5,
  train_frac   = 0.8
)

# ============================================================
# 2. DGP再生成（u と test_idx のみ取得）
# ============================================================
regenerate_u <- function(n, effect_type, seed, cfg) {
  set.seed(seed)
  p <- cfg$p_genes; K <- cfg$K_latent

  A      <- rnorm(n, 0, 1)
  G      <- matrix(rnorm(n * K), n, K)
  Lambda <- matrix(rnorm(p * K, 0, cfg$sigma_lambda), p, K)
  b0     <- rnorm(p, log(50), 0.2)
  bA     <- rnorm(p, 0.15, 0.05)
  eps_E  <- matrix(rnorm(n * p, 0, cfg$sigma_e), n, p)
  E      <- outer(rep(1, n), b0) + outer(A, bA) + G %*% t(Lambda) + eps_E

  # NB counts (seed消費順を合わせる)
  s      <- exp(rnorm(n, 0, cfg$libsize_sd))
  mu_mat <- exp(E) * s
  C      <- matrix(rnbinom(n * p, size = cfg$theta_nb, mu = as.vector(mu_mat)), n, p)

  T_vec  <- rbinom(n, 1, 0.5)

  # Effect modifier score
  u <- rowSums(E[, cfg$M_genes, drop = FALSE])

  list(u = u, T_vec = T_vec)
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
# 3. 各条件の u（testセット）を生成
# ============================================================
conditions <- expand.grid(
  idx_n      = seq_along(cfg$n_values),
  idx_effect = seq_along(cfg$effect_types),
  stringsAsFactors = FALSE
)

u_data_list <- vector("list", nrow(conditions))

for (i in seq_len(nrow(conditions))) {
  idx_n      <- conditions$idx_n[i]
  idx_effect <- conditions$idx_effect[i]
  n          <- cfg$n_values[idx_n]
  effect_type <- cfg$effect_types[idx_effect]
  seed_cond  <- cfg$seed_base + 10000L * (idx_n - 1L) + 100L * (idx_effect - 1L)

  dgp <- regenerate_u(n, effect_type, seed_cond, cfg)
  spl <- stratified_split(dgp$T_vec, cfg$train_frac, seed_cond)

  u_data_list[[i]] <- data.frame(
    n           = n,
    effect_type = effect_type,
    seed_condition = seed_cond,
    id          = spl$test,
    u           = dgp$u[spl$test]
  )
}

u_df <- do.call(rbind, u_data_list)
cat(sprintf("u_df: %d rows\n", nrow(u_df)))

# ============================================================
# 4. predictions CSV と結合
# ============================================================
pred <- read.csv(file.path(out_dir, "predictions_demo_omics50_v02.csv"),
                 stringsAsFactors = FALSE)
pred_test <- subset(pred, split == "test")

merged <- merge(
  pred_test,
  u_df,
  by = c("n", "effect_type", "seed_condition", "id"),
  all.x = TRUE
)
cat(sprintf("merged: %d rows (NA u: %d)\n", nrow(merged),
            sum(is.na(merged$u))))

# ============================================================
# 5. プロット
# ============================================================
methods_all <- c("causal_forest", "drlearner", "ganite", "cevae")
method_labels <- c(
  causal_forest = "Causal Forest (GRF)",
  drlearner     = "DR-Learner",
  ganite        = "GANITE",
  cevae         = "CEVAE"
)

# 条件ラベル（factor順）
merged$condition <- factor(
  paste0("n=", merged$n, "\n", merged$effect_type),
  levels = c(
    "n=100\nlinear",    "n=100\nthreshold3",
    "n=500\nlinear",    "n=500\nthreshold3",
    "n=1000\nlinear",   "n=1000\nthreshold3"
  )
)

make_figure <- function(method_id) {
  dm <- subset(merged, method == method_id)
  if (nrow(dm) == 0) return(NULL)

  panels <- lapply(levels(dm$condition), function(cond) {
    dc <- subset(dm, condition == cond)
    if (nrow(dc) == 0) return(NULL)

    # 有効な予測がある行のみ（τ̂がNA）
    dc_valid <- subset(dc, !is.na(tau_hat))
    has_pred  <- nrow(dc_valid) > 0

    # 軸範囲
    all_tau <- c(dc$tau_true, if (has_pred) dc_valid$tau_hat else NULL)
    y_range <- range(all_tau, na.rm = TRUE)
    y_pad   <- diff(y_range) * 0.05
    x_range <- range(dc$u, na.rm = TRUE)

    # RMSE annotation
    rmse_label <- if (has_pred) {
      sprintf("RMSE=%.2f", sqrt(mean((dc_valid$tau_hat - dc_valid$tau_true)^2)))
    } else "RMSE=NA"

    p <- ggplot() +
      # 真値（灰色）
      geom_point(data = dc,
                 aes(x = u, y = tau_true),
                 color = "#2166AC", alpha = 0.45, size = 1.0) +
      geom_smooth(data = dc,
                  aes(x = u, y = tau_true),
                  method = "loess", se = FALSE,
                  color = "#2166AC", linewidth = 0.9, span = 0.6) +
      labs(x = "Effect modifier score (u)",
           y = expression(tau),
           title = cond) +
      coord_cartesian(xlim = x_range,
                      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)) +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
               label = rmse_label, size = 3.0, color = "grey20") +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(size = 9, hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 7)
      )

    # 予測値（赤）
    if (has_pred) {
      p <- p +
        geom_point(data = dc_valid,
                   aes(x = u, y = tau_hat),
                   color = "#D6604D", alpha = 0.55, size = 1.0) +
        geom_smooth(data = dc_valid,
                    aes(x = u, y = tau_hat),
                    method = "loess", se = FALSE,
                    color = "#D6604D", linewidth = 0.9, span = 0.6)
    }
    p
  })

  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0) return(NULL)

  # 凡例パネル（ダミー）
  legend_data <- data.frame(
    u   = c(0, 0),
    tau = c(0, 0),
    type = factor(c("True τ", "Predicted τ̂"), levels = c("True τ", "Predicted τ̂"))
  )
  legend_p <- ggplot(legend_data, aes(x = u, y = tau, color = type)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("True τ" = "#2166AC", "Predicted τ̂" = "#D6604D"),
                       name = NULL) +
    theme_void() +
    theme(
      legend.position  = "right",
      legend.key.size  = unit(0.8, "lines"),
      legend.text      = element_text(size = 9)
    )
  legend_grob <- cowplot::get_legend(legend_p)

  wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title    = paste0(method_labels[method_id],
                        "  —  Effect modifier score (u) vs CATE (τ)"),
      subtitle = "Blue: true τ  |  Red: predicted τ̂  |  Lines: loess smooth  |  Test set only",
      theme    = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    )
}

# cowplot があるか確認
has_cowplot <- requireNamespace("cowplot", quietly = TRUE)

make_figure_simple <- function(method_id) {
  dm <- subset(merged, method == method_id)
  if (nrow(dm) == 0) return(NULL)

  # 真値データ（全サンプルのテスト）と予測値データを縦に結合して ggplot
  d_true <- dm[, c("condition", "u", "tau_true")]
  d_true$type    <- "True τ"
  d_true$tau_val <- d_true$tau_true

  dm_valid <- subset(dm, !is.na(tau_hat))
  d_pred <- dm_valid[, c("condition", "u", "tau_hat")]
  d_pred$type    <- "Predicted τ̂"
  d_pred$tau_val <- d_pred$tau_hat

  d_combined <- rbind(
    d_true[, c("condition", "u", "type", "tau_val")],
    d_pred[, c("condition", "u", "type", "tau_val")]
  )
  d_combined$type <- factor(d_combined$type,
                             levels = c("True τ", "Predicted τ̂"))

  # RMSEアノテーション
  rmse_labels <- aggregate(
    cbind(sq_err = (tau_hat - tau_true)^2) ~ condition,
    data = dm_valid,
    FUN  = function(x) sprintf("RMSE=%.2f", sqrt(mean(x)))
  )
  if (nrow(rmse_labels) == 0) {
    rmse_labels <- data.frame(condition = levels(dm$condition),
                               sq_err = "RMSE=NA",
                               stringsAsFactors = FALSE)
  }

  ggplot(d_combined, aes(x = u, y = tau_val, color = type)) +
    geom_point(alpha = 0.45, size = 0.9) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, span = 0.6,
                aes(group = type)) +
    scale_color_manual(
      values = c("True τ" = "#2166AC", "Predicted τ̂" = "#D6604D"),
      name   = NULL
    ) +
    geom_text(data = rmse_labels,
              aes(x = Inf, y = Inf, label = sq_err),
              hjust = 1.1, vjust = 1.5,
              size = 3.0, color = "grey20",
              inherit.aes = FALSE) +
    facet_wrap(~condition, ncol = 2, scales = "free") +
    labs(
      title    = paste0(method_labels[method_id],
                        "  —  Effect modifier score (u) vs CATE (τ)"),
      subtitle = "Blue: true τ  |  Red: predicted τ̂  |  loess smooth  |  Test set only",
      x        = "Effect modifier score (u)",
      y        = expression(CATE ~ (tau))
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 10, color = "grey40"),
      legend.position  = "bottom",
      strip.text       = element_text(size = 9),
      legend.text      = element_text(size = 10)
    )
}

for (m in methods_all) {
  fig <- make_figure_simple(m)
  if (is.null(fig)) {
    cat(sprintf("  [%s] No data — skipped\n", m))
    next
  }
  fname <- sprintf("dgp_viz_12_u_vs_cate_%s.png", m)
  ggsave(file.path(out_dir, fname), fig,
         width = 10, height = 10, dpi = 150)
  cat(sprintf("Saved: %s\n", fname))
}

cat("\n全図の保存完了:", out_dir, "\n")
