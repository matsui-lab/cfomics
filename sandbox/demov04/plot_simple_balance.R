#!/usr/bin/env Rscript
# plot_simple_balance.R
# シンプルな分布比較

suppressPackageStartupMessages({
  library(ggplot2)
})

config <- list(
  m_E = log(50), s_E = 0.35, theta_nb = 30
)

seed_cond <- 20260425
set.seed(seed_cond)

n <- 500
S <- rbinom(n, 1, 0.5)

# RNA-seqカウントデータ
E_true <- rnorm(n, config$m_E, config$s_E)
mu_nb <- exp(E_true)
C_count <- rnbinom(n, size = config$theta_nb, mu = mu_nb)

# 治療割り当て
T_vec <- rbinom(n, 1, 0.5)

# データフレーム（明確なラベル）
dat <- data.frame(
  C_count = C_count,
  Treatment = factor(ifelse(T_vec == 1, "治療あり (T=1)", "治療なし (T=0)"),
                     levels = c("治療なし (T=0)", "治療あり (T=1)"))
)

# 統計量
stats <- aggregate(C_count ~ Treatment, dat, function(x) {
  c(n = length(x), mean = mean(x), sd = sd(x))
})
stats_df <- do.call(data.frame, stats)

cat("=== geneA カウントデータの分布 ===\n\n")
cat("治療なし (T=0):\n")
cat(sprintf("  n = %d, 平均 = %.2f, SD = %.2f\n",
            stats_df$C_count.n[1], stats_df$C_count.mean[1], stats_df$C_count.sd[1]))
cat("\n治療あり (T=1):\n")
cat(sprintf("  n = %d, 平均 = %.2f, SD = %.2f\n",
            stats_df$C_count.n[2], stats_df$C_count.mean[2], stats_df$C_count.sd[2]))

out_dir <- "sandbox/demov04/results_corrected"

# シンプルな重ね合わせヒストグラム
p <- ggplot(dat, aes(x = C_count, fill = Treatment)) +
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity", alpha = 0.5, bins = 40,
                 color = "black", linewidth = 0.3) +
  geom_density(aes(color = Treatment), linewidth = 1.5, fill = NA) +
  scale_fill_manual(values = c("治療なし (T=0)" = "#E41A1C",
                                "治療あり (T=1)" = "#377EB8")) +
  scale_color_manual(values = c("治療なし (T=0)" = "#8B0000",
                                 "治療あり (T=1)" = "#00008B")) +
  labs(title = "geneA カウントデータの分布（対数変換前）",
       subtitle = sprintf("治療なし: n=%d, 平均=%.1f | 治療あり: n=%d, 平均=%.1f",
                         stats_df$C_count.n[1], stats_df$C_count.mean[1],
                         stats_df$C_count.n[2], stats_df$C_count.mean[2]),
       x = "RNA-seq カウント", y = "密度",
       fill = NULL, color = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11))

ggsave(file.path(out_dir, "balance_simple.png"),
       p, width = 10, height = 7, dpi = 150)

cat("\n✓ 保存: balance_simple.png\n")
cat("\n結論: 両群の分布はほぼ完全に重なっている\n")
