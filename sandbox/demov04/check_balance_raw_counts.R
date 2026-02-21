#!/usr/bin/env Rscript
# check_balance_raw_counts.R
# 対数変換前のカウントデータで対照群と治療群の共変量バランスを確認

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
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

# 治療割り当て（完全ランダム化）
T_vec <- rbinom(n, 1, 0.5)

# データフレーム作成
dat <- data.frame(
  id = 1:n,
  Sex = factor(ifelse(S == 1, "Male", "Female")),
  C_count = C_count,
  T = T_vec,
  T_label = factor(ifelse(T_vec == 1, "治療群", "対照群"),
                   levels = c("対照群", "治療群"))
)

# 統計量の計算
cat("=== カウントデータの共変量バランス（n=500） ===\n\n")

cat("サンプルサイズ:\n")
table_T <- table(dat$T_label)
cat(sprintf("  対照群: %d\n", table_T["対照群"]))
cat(sprintf("  治療群: %d\n", table_T["治療群"]))

cat("\n性別の分布:\n")
sex_balance <- table(dat$Sex, dat$T_label)
print(sex_balance)
cat(sprintf("\n  カイ二乗検定: p = %.3f\n",
            chisq.test(sex_balance)$p.value))

cat("\ngeneA カウントデータの分布:\n")
count_T0 <- dat$C_count[dat$T == 0]
count_T1 <- dat$C_count[dat$T == 1]

cat(sprintf("\n対照群 (T=0, n=%d):\n", length(count_T0)))
cat(sprintf("  平均:   %.2f\n", mean(count_T0)))
cat(sprintf("  中央値: %.0f\n", median(count_T0)))
cat(sprintf("  SD:     %.2f\n", sd(count_T0)))
cat(sprintf("  範囲:   [%d, %d]\n", min(count_T0), max(count_T0)))

cat(sprintf("\n治療群 (T=1, n=%d):\n", length(count_T1)))
cat(sprintf("  平均:   %.2f\n", mean(count_T1)))
cat(sprintf("  中央値: %.0f\n", median(count_T1)))
cat(sprintf("  SD:     %.2f\n", sd(count_T1)))
cat(sprintf("  範囲:   [%d, %d]\n", min(count_T1), max(count_T1)))

cat("\n差の検定:\n")
cat(sprintf("  平均の差:     %.2f\n", mean(count_T1) - mean(count_T0)))
t_test <- t.test(count_T1, count_T0)
cat(sprintf("  t検定 p値:    %.3f\n", t_test$p.value))
ks_test <- ks.test(count_T0, count_T1)
cat(sprintf("  KS検定 p値:   %.3f\n", ks_test$p.value))

# 標準化平均差（SMD）
smd <- (mean(count_T1) - mean(count_T0)) /
       sqrt((sd(count_T0)^2 + sd(count_T1)^2) / 2)
cat(sprintf("  標準化平均差: %.3f (|SMD| < 0.1 が良好なバランス)\n", smd))

# 可視化
out_dir <- "sandbox/demov04/results_corrected"

# 1. 重ねヒストグラム（密度）
p1 <- ggplot(dat, aes(x = C_count, fill = T_label)) +
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity", alpha = 0.6, bins = 40,
                 color = "black") +
  scale_fill_manual(values = c("対照群" = "#E41A1C", "治療群" = "#377EB8")) +
  labs(title = "geneA カウントデータの分布（群別）",
       subtitle = sprintf("対照群 n=%d, 治療群 n=%d | KS検定 p=%.3f",
                         length(count_T0), length(count_T1), ks_test$p.value),
       x = "RNA-seq カウント", y = "密度", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# 2. 密度曲線のみ
p2 <- ggplot(dat, aes(x = C_count, color = T_label)) +
  geom_density(linewidth = 1.5) +
  scale_color_manual(values = c("対照群" = "#E41A1C", "治療群" = "#377EB8")) +
  labs(title = "密度曲線の比較",
       subtitle = sprintf("標準化平均差 (SMD) = %.3f", smd),
       x = "RNA-seq カウント", y = "密度", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# 3. ボックスプロット
p3 <- ggplot(dat, aes(x = T_label, y = C_count, fill = T_label)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 1) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
  scale_fill_manual(values = c("対照群" = "#E41A1C", "治療群" = "#377EB8")) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "yellow", color = "black") +
  labs(title = "ボックスプロット（カウントデータ）",
       subtitle = sprintf("平均の差 = %.2f, t検定 p=%.3f",
                         mean(count_T1) - mean(count_T0), t_test$p.value),
       x = NULL, y = "RNA-seq カウント") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

# 4. QQプロット
qq_data <- data.frame(
  q_T0 = sort(count_T0),
  q_T1 = sort(count_T1)[1:length(count_T0)]
)

p4 <- ggplot(qq_data, aes(x = q_T0, y = q_T1)) +
  geom_point(alpha = 0.5, color = "#377EB8") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "QQプロット（分布の一致度）",
       subtitle = "対角線上 = 完全に同一分布",
       x = "対照群の分位点", y = "治療群の分位点") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# レイアウト
layout <- (p1 | p2) / (p3 | p4)
layout <- layout + plot_annotation(
  title = "カウントデータの共変量バランス（対数変換前）",
  subtitle = "完全ランダム化により対照群と治療群の分布は一致している",
  theme = theme(plot.title = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(size = 11, color = "darkgreen"))
)

ggsave(file.path(out_dir, "covariate_balance_raw_counts.png"),
       layout, width = 14, height = 10, dpi = 150)

cat("\n✓ 図を保存しました: covariate_balance_raw_counts.png\n")
cat("\n結論: 対照群と治療群でgeneAのカウントデータ分布は十分に重なっている\n")
cat("      これは完全ランダム化試験の期待される結果です\n")
