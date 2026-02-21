#!/usr/bin/env Rscript
# show_raw_counts_distribution.R
# geneAの対数変換前（カウントデータ）の分布を表示

suppressPackageStartupMessages({
  library(ggplot2)
})

config <- list(
  seed_base = 20260215L,
  m_E = log(50), s_E = 0.35, theta_nb = 30,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0,
  mu0 = 120.0, sigma_y = 5.0, sigma_G = 0.08
)

seed_cond <- 20260425
set.seed(seed_cond)

n <- 500
S <- rbinom(n, 1, 0.5)
Sc <- S - 0.5

# RNA-seqカウントデータの生成
E_true <- rnorm(n, config$m_E, config$s_E)
mu_nb <- exp(E_true)
C_count <- rnbinom(n, size = config$theta_nb, mu = mu_nb)

# log1p変換後のデータ（実際に解析で使用）
geneA_log <- log1p(C_count)
delta <- rnorm(n, 0, config$sigma_G)
geneA_obs <- geneA_log + delta

# 統計量の出力
cat("=== geneA 対数変換前（カウントデータ）の分布 ===\n\n")
cat(sprintf("サンプルサイズ: %d\n", n))
cat(sprintf("シード: %d\n\n", seed_cond))

cat("カウントデータ (C_count):\n")
cat(sprintf("  範囲:     [%d, %d]\n", min(C_count), max(C_count)))
cat(sprintf("  平均:     %.2f\n", mean(C_count)))
cat(sprintf("  中央値:   %.0f\n", median(C_count)))
cat(sprintf("  SD:       %.2f\n", sd(C_count)))
cat(sprintf("  四分位点: Q1=%.0f, Q3=%.0f\n",
            quantile(C_count, 0.25), quantile(C_count, 0.75)))

cat("\nlog1p変換後 (geneA_log = log1p(C_count)):\n")
cat(sprintf("  範囲:     [%.3f, %.3f]\n", min(geneA_log), max(geneA_log)))
cat(sprintf("  平均:     %.3f\n", mean(geneA_log)))
cat(sprintf("  中央値:   %.3f\n", median(geneA_log)))
cat(sprintf("  SD:       %.3f\n", sd(geneA_log)))

cat("\nノイズ追加後 (geneA_obs = geneA_log + ε):\n")
cat(sprintf("  範囲:     [%.3f, %.3f]\n", min(geneA_obs), max(geneA_obs)))
cat(sprintf("  平均:     %.3f\n", mean(geneA_obs)))
cat(sprintf("  中央値:   %.3f\n", median(geneA_obs)))
cat(sprintf("  SD:       %.3f\n", sd(geneA_obs)))

# 分布の形状（簡易統計）
cat("\n分布の形状:\n")
skew_simple <- mean((C_count - mean(C_count))^3) / sd(C_count)^3
cat(sprintf("  歪度(簡易):     %.3f (カウントデータ)\n", skew_simple))
cat(sprintf("  変動係数:       %.3f (SD/平均)\n", sd(C_count)/mean(C_count)))

# 可視化
out_dir <- "sandbox/demov04/results_corrected"

# カウントデータのヒストグラム
p1 <- ggplot(data.frame(count = C_count), aes(x = count)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "#4DAF4A", alpha = 0.7, color = "black") +
  geom_density(color = "blue", linewidth = 1.2) +
  labs(title = "geneA カウントデータの分布",
       subtitle = sprintf("Negative Binomial (μ=exp(E), size=%.0f) | n=%d",
                         config$theta_nb, n),
       x = "RNA-seq カウント", y = "密度") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# log1p変換後のヒストグラム
p2 <- ggplot(data.frame(log_count = geneA_log), aes(x = log_count)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "#377EB8", alpha = 0.7, color = "black") +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "log1p変換後の分布",
       subtitle = "geneA_log = log1p(count)",
       x = "log1p(カウント)", y = "密度") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ノイズ追加後のヒストグラム
p3 <- ggplot(data.frame(obs = geneA_obs), aes(x = obs)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "#E41A1C", alpha = 0.7, color = "black") +
  geom_density(color = "darkgreen", linewidth = 1.2) +
  labs(title = "測定ノイズ追加後の分布（解析に使用）",
       subtitle = sprintf("geneA_obs = log1p(count) + N(0, %.3f²)",
                         config$sigma_G),
       x = "観測値 geneA", y = "密度") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# 変換プロセスの可視化
library(patchwork)
layout <- p1 / p2 / p3
layout <- layout + plot_annotation(
  title = "geneA の変換プロセス（n=500）",
  subtitle = "カウントデータ → log1p変換 → 測定ノイズ追加",
  theme = theme(plot.title = element_text(face = "bold", size = 14))
)

ggsave(file.path(out_dir, "geneA_transformation_process.png"),
       layout, width = 10, height = 12, dpi = 150)

cat("\n✓ 図を保存しました: geneA_transformation_process.png\n")
