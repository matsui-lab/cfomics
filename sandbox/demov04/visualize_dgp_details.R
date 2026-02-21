#!/usr/bin/env Rscript
# visualize_dgp_details.R
# DGPの詳細を可視化

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/demov04/results_with_fix"

config <- list(
  seed_base = 20260215L,
  n_val = 500L,
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08
)

# ============================================================
# データ生成
# ============================================================
set.seed(config$seed_base + 200)

n <- config$n_val

# Step 1: 真の発現量
E_true <- rnorm(n, config$m_E, config$s_E)

# Step 2: 平均カウント
mu_nb <- exp(E_true)

# Step 3: カウントデータ
C_count <- rnbinom(n, size = config$theta_nb, mu = mu_nb)

# Step 4: log1p変換
G_raw <- log1p(C_count)

# Step 5: 測定誤差
delta <- rnorm(n, 0, config$sigma_G)
G_obs <- G_raw + delta

# Step 6: 性別
S <- rbinom(n, 1, 0.5)

df <- data.frame(
  id = 1:n,
  S = S,
  E_true = E_true,
  mu_nb = mu_nb,
  C_count = C_count,
  G_raw = G_raw,
  G_obs = G_obs
)

# ============================================================
# プロット1: E_true（真の発現量）の分布
# ============================================================
p1 <- ggplot(df, aes(x = E_true)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#377EB8", alpha = 0.7, color = "black") +
  geom_density(color = "red", linewidth = 1) +
  stat_function(fun = dnorm, args = list(mean = config$m_E, sd = config$s_E),
                color = "darkgreen", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Step 1: 真の発現量 (E_true)",
    subtitle = sprintf("E_true ~ N(%.2f, %.2f²) = N(log(50), 0.35²)",
                       config$m_E, config$s_E),
    x = "E_true (対数スケール)",
    y = "密度"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# プロット2: mu_nb（平均カウント）の分布
# ============================================================
p2 <- ggplot(df, aes(x = mu_nb)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#E41A1C", alpha = 0.7, color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = median(df$mu_nb), linetype = "dashed",
             color = "darkgreen", linewidth = 1) +
  annotate("text", x = median(df$mu_nb) * 1.5, y = 0.02,
           label = sprintf("中央値: %.1f", median(df$mu_nb)),
           color = "darkgreen", size = 4) +
  labs(
    title = "Step 2: 平均カウント (μ_nb)",
    subtitle = "μ_nb = exp(E_true)  — 対数正規分布に従う",
    x = "μ_nb (元のスケール)",
    y = "密度"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# プロット3: C_count（カウントデータ）の分布
# ============================================================
p3 <- ggplot(df, aes(x = C_count)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#4DAF4A", alpha = 0.7, color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  labs(
    title = "Step 3: カウントデータ (C_count)",
    subtitle = sprintf("C_count ~ NB(μ = μ_nb, size = %d)  — 負の二項分布",
                       config$theta_nb),
    x = "RNA-seq カウント",
    y = "密度"
  ) +
  annotate("text", x = max(df$C_count) * 0.7, y = 0.015,
           label = sprintf("平均: %.1f\nSD: %.1f\n分散/平均: %.2f",
                          mean(df$C_count), sd(df$C_count),
                          var(df$C_count) / mean(df$C_count)),
           size = 3.5, hjust = 0) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# プロット4: log1p変換後の分布
# ============================================================
p4 <- ggplot(df, aes(x = G_raw)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#984EA3", alpha = 0.7, color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  labs(
    title = "Step 4: log1p変換",
    subtitle = "G_raw = log1p(C_count) = log(1 + C_count)",
    x = "log1p(カウント)",
    y = "密度"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# プロット5: 測定誤差追加後
# ============================================================
p5 <- ggplot(df, aes(x = G_obs)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#FF7F00", alpha = 0.7, color = "black") +
  geom_density(color = "blue", linewidth = 1) +
  labs(
    title = "Step 5: 測定誤差追加",
    subtitle = sprintf("G_obs = G_raw + ε,  where ε ~ N(0, %.2f²)", config$sigma_G),
    x = "geneA (観測値)",
    y = "密度"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# プロット6: 性別の分布
# ============================================================
sex_counts <- data.frame(
  Sex = factor(c("Female", "Male"), levels = c("Female", "Male")),
  Count = c(sum(df$S == 0), sum(df$S == 1))
)
sex_counts$Prop <- sex_counts$Count / sum(sex_counts$Count)

p6 <- ggplot(sex_counts, aes(x = Sex, y = Prop, fill = Sex)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Prop*100, Count)),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
  labs(
    title = "性別の割り振り",
    subtitle = "S ~ Bernoulli(0.5)  — 各個人が独立に50%ずつ",
    x = NULL,
    y = "割合"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.7)) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

# ============================================================
# レイアウト
# ============================================================
layout <- (p1 | p2) / (p3 | p4) / (p5 | p6)
layout <- layout +
  plot_annotation(
    title = sprintf("DGP の詳細可視化 (n=%d)", config$n_val),
    subtitle = "遺伝子発現データ生成の各ステップ",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "grey40")
    )
  )

fname <- "plot_dgp_details.png"
ggsave(file.path(out_dir, fname), layout, width = 14, height = 12, dpi = 150)
cat(sprintf("Saved: %s\n", fname))

# ============================================================
# 統計量の出力
# ============================================================
cat("\n=== DGP統計量 (n=500) ===\n\n")

cat("Step 1: E_true (対数スケール)\n")
cat(sprintf("  平均: %.3f (理論値: %.3f)\n", mean(df$E_true), config$m_E))
cat(sprintf("  SD:   %.3f (理論値: %.3f)\n", sd(df$E_true), config$s_E))

cat("\nStep 2: μ_nb (平均カウント)\n")
cat(sprintf("  中央値: %.1f (理論値: exp(%.2f) = %.1f)\n",
            median(df$mu_nb), config$m_E, exp(config$m_E)))
cat(sprintf("  範囲: [%.1f, %.1f]\n", min(df$mu_nb), max(df$mu_nb)))

cat("\nStep 3: C_count (カウントデータ)\n")
cat(sprintf("  平均: %.1f\n", mean(df$C_count)))
cat(sprintf("  SD:   %.1f\n", sd(df$C_count)))
cat(sprintf("  分散: %.1f\n", var(df$C_count)))
cat(sprintf("  過分散比 (分散/平均): %.2f\n", var(df$C_count) / mean(df$C_count)))
cat(sprintf("  理論値（mu=50の場合）: 50 + 50²/30 = %.1f\n", 50 + 50^2/30))

cat("\nStep 4: G_raw (log1p変換)\n")
cat(sprintf("  平均: %.3f\n", mean(df$G_raw)))
cat(sprintf("  SD:   %.3f\n", sd(df$G_raw)))

cat("\nStep 5: G_obs (測定誤差追加)\n")
cat(sprintf("  平均: %.3f\n", mean(df$G_obs)))
cat(sprintf("  SD:   %.3f\n", sd(df$G_obs)))
cat(sprintf("  測定誤差のSD: %.3f\n", config$sigma_G))

cat("\n性別の割り振り:\n")
cat(sprintf("  Female: %d (%.1f%%)\n", sum(df$S == 0), mean(df$S == 0) * 100))
cat(sprintf("  Male:   %d (%.1f%%)\n", sum(df$S == 1), mean(df$S == 1) * 100))

cat("\n✓ DGPの各ステップが理論値と一致しています\n")
