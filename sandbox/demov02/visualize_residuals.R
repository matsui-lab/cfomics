#!/usr/bin/env Rscript
# visualize_residuals.R
# 各モデルの残差プロット（tau_hat - tau_true vs tau_true）
# 条件（n × effect_type）別パネル、モデルごとに1 figure
# Run from repo root: Rscript sandbox/visualize_residuals.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/benchmark_results/cate_demo_v02"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred <- read.csv(file.path(out_dir, "predictions_demo_omics50_v02.csv"),
                 stringsAsFactors = FALSE)

# テストセットのみ使用、tau_hat が NA でない行のみ
pred_test <- subset(pred, split == "test" & !is.na(tau_hat))
pred_test$residual <- pred_test$tau_hat - pred_test$tau_true

# 条件ラベル
pred_test$condition <- factor(
  paste0("n=", pred_test$n, "\n", pred_test$effect_type),
  levels = c(
    "n=100\nlinear",    "n=100\nthreshold3",
    "n=500\nlinear",    "n=500\nthreshold3",
    "n=1000\nlinear",   "n=1000\nthreshold3"
  )
)

methods_all <- c("causal_forest", "drlearner", "ganite", "cevae")
method_labels <- c(
  causal_forest = "Causal Forest (GRF)",
  drlearner     = "DR-Learner",
  ganite        = "GANITE",
  cevae         = "CEVAE"
)

# 方法ごとにRMSEを計算
rmse_df <- aggregate(
  cbind(residual2 = residual^2) ~ method + n + effect_type,
  data  = pred_test,
  FUN   = mean
)
rmse_df$rmse <- sqrt(rmse_df$residual2)
rmse_df$condition <- factor(
  paste0("n=", rmse_df$n, "\n", rmse_df$effect_type),
  levels = levels(pred_test$condition)
)

# 色設定（effect_type別）
eff_colors <- c(linear = "#E41A1C", threshold3 = "#377EB8")

make_figure <- function(method_id) {
  d <- subset(pred_test, method == method_id)
  if (nrow(d) == 0) return(NULL)

  r <- subset(rmse_df, method == method_id)

  panels <- lapply(levels(d$condition), function(cond) {
    dc <- subset(d, condition == cond)
    rc <- subset(r, condition == cond)
    rmse_val <- if (nrow(rc) > 0) rc$rmse[1] else NA

    # effect_type を cond から取り出す
    eff <- if (grepl("linear", cond)) "linear" else "threshold3"
    pt_color <- eff_colors[eff]

    # 軸範囲（対称）
    max_abs_res <- max(abs(dc$residual), na.rm = TRUE)
    max_abs_res <- ceiling(max_abs_res / 5) * 5   # 5刻みに丸め

    ggplot(dc, aes(x = tau_true, y = residual)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
      geom_point(alpha = 0.5, size = 1.2, color = pt_color) +
      geom_smooth(method = "loess", se = TRUE,
                  color = "black", fill = "grey80",
                  linewidth = 0.7, span = 0.8) +
      coord_cartesian(ylim = c(-max_abs_res, max_abs_res)) +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
               label = if (!is.na(rmse_val)) sprintf("RMSE=%.2f", rmse_val) else "RMSE=NA",
               size = 3.2, color = "grey20") +
      labs(
        title = cond,
        x     = expression("True CATE (" * tau * ")"),
        y     = expression(hat(tau) - tau)
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title   = element_text(size = 9, hjust = 0.5),
        axis.title   = element_text(size = 8),
        axis.text    = element_text(size = 7)
      )
  })

  # 有効なパネルのみ
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) == 0) return(NULL)

  # 2列配置
  n_panels <- length(panels)
  n_cols   <- 2
  n_rows   <- ceiling(n_panels / n_cols)

  wrap_plots(panels, ncol = n_cols, nrow = n_rows) +
    plot_annotation(
      title    = paste0(method_labels[method_id], " — Residuals (τ̂ − τ) vs True CATE (τ)"),
      subtitle = "Test set only  |  loess smooth ± 95% CI  |  dashed line: zero residual",
      theme    = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    )
}

# 出力
for (m in methods_all) {
  fig <- make_figure(m)
  if (is.null(fig)) {
    cat(sprintf("  [%s] No data — skipped\n", m))
    next
  }
  fname <- sprintf("dgp_viz_11_residuals_%s.png", m)
  ggsave(file.path(out_dir, fname), fig,
         width = 10, height = 7, dpi = 150)
  cat(sprintf("Saved: %s\n", fname))
}

cat("\n全残差プロットの保存完了:", out_dir, "\n")
