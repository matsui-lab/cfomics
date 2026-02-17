#!/usr/bin/env Rscript
# visualize_demov03.R
# demov03 結果の可視化
# Run from repo root: Rscript sandbox/demov03/visualize_demov03.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/demov03/results"

# ============================================================
# データ読み込み
# ============================================================
df_long    <- read.csv(file.path(out_dir, "results_demov03_long.csv"),
                       stringsAsFactors = FALSE)
df_summary <- read.csv(file.path(out_dir, "results_demov03_summary.csv"),
                       stringsAsFactors = FALSE)

# 異常値除外（RMSE > 1000 は数値的失敗とみなす）
df_long_clean    <- subset(df_long,    is.na(rmse_tau) | rmse_tau < 1000)
df_summary_clean <- subset(df_summary, is.na(mean_rmse) | mean_rmse < 1000)

# 手法・covset の factor 化
method_levels <- c("causal_forest", "drlearner", "ganite", "cevae")
method_labels <- c("Causal Forest", "DR-Learner", "GANITE", "CEVAE")
method_colors <- c("#E41A1C", "#377EB8", "#FF7F00", "#4DAF4A")
names(method_colors) <- method_labels   # factor後のラベルに合わせる

df_summary_clean$method <- factor(df_summary_clean$method,
                                   levels = method_levels,
                                   labels = method_labels)
df_summary_clean$covset <- factor(df_summary_clean$covset,
                                   levels = c("AS", "ASG"))
df_summary_clean$n <- as.integer(df_summary_clean$n)

# ============================================================
# 図1: 折れ線グラフ（covset = AS）
# ============================================================
make_rmse_plot <- function(cs, title_str, tau_sd = 5) {
  d <- subset(df_summary_clean, covset == cs)

  ggplot(d, aes(x = n, y = mean_rmse, color = method, group = method)) +
    geom_hline(yintercept = tau_sd, linetype = "dashed",
               color = "grey50", linewidth = 0.7) +
    annotate("text", x = max(d$n, na.rm = TRUE) * 0.95, y = tau_sd * 1.05,
             label = sprintf("τ SD = %.0f", tau_sd),
             hjust = 1, size = 3.2, color = "grey40") +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    scale_color_manual(values = method_colors,
                       labels = method_labels,
                       name   = "Method") +
    scale_x_continuous(breaks = c(100, 500, 1000)) +
    labs(
      title    = title_str,
      subtitle = sprintf("covset = %s  |  R=%d  |  RMSE < 1000 のみ表示",
                         cs, max(df_long$rep, na.rm = TRUE)),
      x        = "Sample size (n)",
      y        = expression(RMSE[tau])
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = "right"
    )
}

p_AS  <- make_rmse_plot("AS",  "Covariates: Age + Sex only")
p_ASG <- make_rmse_plot("ASG", "Covariates: Age + Sex + geneA")

ggsave(file.path(out_dir, "plot_rmse_by_n_covset_AS.png"),
       p_AS,  width = 8, height = 5, dpi = 150)
cat("Saved: plot_rmse_by_n_covset_AS.png\n")

ggsave(file.path(out_dir, "plot_rmse_by_n_covset_ASG.png"),
       p_ASG, width = 8, height = 5, dpi = 150)
cat("Saved: plot_rmse_by_n_covset_ASG.png\n")

# ============================================================
# 図2: AS vs ASG 並列比較（1枚）
# ============================================================
fig2 <- (p_AS | p_ASG) +
  plot_annotation(
    title    = "CATE RMSE: geneA なし（AS）vs あり（ASG）  [demov03, R=1]",
    subtitle = "Dashed line: τ SD = 5 (oracle baseline)  |  Outliers (RMSE>1000) excluded",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(file.path(out_dir, "plot_rmse_AS_vs_ASG.png"),
       fig2, width = 14, height = 5.5, dpi = 150)
cat("Saved: plot_rmse_AS_vs_ASG.png\n")

# ============================================================
# 図3: AS vs ASG の RMSE 差（ASG - AS、負なら geneA で改善）
# ============================================================
df_wide <- merge(
  subset(df_summary_clean, covset == "AS",  select = c("n", "method", "mean_rmse")),
  subset(df_summary_clean, covset == "ASG", select = c("n", "method", "mean_rmse")),
  by = c("n", "method"), suffixes = c("_AS", "_ASG")
)
df_wide$delta_rmse <- df_wide$mean_rmse_ASG - df_wide$mean_rmse_AS

fig3 <- ggplot(df_wide, aes(x = n, y = delta_rmse, color = method, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  annotate("text", x = 50, y = 0.3, label = "← geneA 追加で改善",
           hjust = 0, size = 3.2, color = "grey40") +
  annotate("text", x = 50, y = -0.3, label = "→ geneA 追加で悪化",
           hjust = 0, size = 3.2, color = "grey40") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors, name = "Method") +
  scale_x_continuous(breaks = c(100, 500, 1000)) +
  labs(
    title    = "ΔRMSE = RMSE(ASG) − RMSE(AS)  [geneA 追加の効果]",
    subtitle = "負の値 = geneA を加えることで RMSE が改善  |  demov03, R=1",
    x        = "Sample size (n)",
    y        = expression(Delta * RMSE[tau] ~ "(ASG − AS)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "plot_delta_rmse_geneA_effect.png"),
       fig3, width = 8, height = 5, dpi = 150)
cat("Saved: plot_delta_rmse_geneA_effect.png\n")

# ============================================================
# 図4: 全結果のヒートマップ（method × n × covset）
# ============================================================
df_heat <- df_summary_clean
df_heat$label <- ifelse(is.na(df_heat$mean_rmse), "NA",
                        sprintf("%.2f", df_heat$mean_rmse))
df_heat$cell_val <- pmin(df_heat$mean_rmse, 20)   # 20でキャップ（外れ値対策）

fig4 <- ggplot(df_heat, aes(x = factor(n), y = method, fill = cell_val)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.2, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "#FFFFBF", high = "#D6604D",
                       midpoint = 10, limits = c(0, 20),
                       na.value = "grey80",
                       name = expression(RMSE[tau])) +
  facet_wrap(~covset, labeller = label_both) +
  labs(
    title    = "CATE RMSE ヒートマップ  [demov03, R=1]",
    subtitle = "色: 青=良好 / 赤=不良  |  上限20でキャップ  |  外れ値(>1000)は除外",
    x        = "Sample size (n)",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    strip.text    = element_text(size = 11, face = "bold")
  )

ggsave(file.path(out_dir, "plot_heatmap_rmse.png"),
       fig4, width = 10, height = 5, dpi = 150)
cat("Saved: plot_heatmap_rmse.png\n")

cat("\n全図の保存完了:", out_dir, "\n")
