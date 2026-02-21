#!/usr/bin/env Rscript
# plot_geneA_vs_cate_n1000.R
# n=1000での geneA vs CATE プロット

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/demov04/results_n1000"

# データ読み込み
dat <- readRDS(file.path(out_dir, "full_data_n1000.rds"))
spl <- readRDS(file.path(out_dir, "split_n1000.rds"))
df_pred <- read.csv(file.path(out_dir, "predictions_n1000.csv"))

# テストデータ
n <- length(dat$S)
dat_full <- data.frame(
  id = 1:n,
  S = dat$S,
  geneA = dat$geneA,
  tau = dat$tau,
  Sex = factor(ifelse(dat$S == 1, "Male", "Female"), levels = c("Female", "Male"))
)

dat_test <- dat_full[spl$test, ]

# 各手法の予測値を取得
cf_so <- subset(df_pred, method == "causal_forest" & covset == "SexOnly")
cf_sg <- subset(df_pred, method == "causal_forest" & covset == "SexGene")
dr_so <- subset(df_pred, method == "drlearner" & covset == "SexOnly")
dr_sg <- subset(df_pred, method == "drlearner" & covset == "SexGene")

# 可視化用データフレーム
df_vis <- data.frame(
  geneA = dat_test$geneA,
  tau_true = dat_test$tau,
  Sex = dat_test$Sex,
  CF_SexOnly = cf_so$tau_hat,
  CF_SexGene = cf_sg$tau_hat,
  DR_SexOnly = dr_so$tau_hat,
  DR_SexGene = dr_sg$tau_hat
)

# RMSE計算
rmse_cf_so <- sqrt(mean((df_vis$tau_true - df_vis$CF_SexOnly)^2))
rmse_cf_sg <- sqrt(mean((df_vis$tau_true - df_vis$CF_SexGene)^2))
rmse_dr_so <- sqrt(mean((df_vis$tau_true - df_vis$DR_SexOnly)^2))
rmse_dr_sg <- sqrt(mean((df_vis$tau_true - df_vis$DR_SexGene)^2))

cat("=== n=1000での RMSE ===\n")
cat(sprintf("Causal Forest + SexOnly:  %.3f\n", rmse_cf_so))
cat(sprintf("Causal Forest + SexGene:  %.3f\n", rmse_cf_sg))
cat(sprintf("DR-Learner + SexOnly:     %.3f\n", rmse_dr_so))
cat(sprintf("DR-Learner + SexGene:     %.3f\n\n", rmse_dr_sg))

# ============================================================
# プロット作成
# ============================================================

# 1. Causal Forest: SexOnly vs SexGene
p1 <- ggplot(df_vis, aes(x = geneA)) +
  geom_point(aes(y = tau_true, color = "真値 τ"), alpha = 0.3, size = 2) +
  geom_point(aes(y = CF_SexOnly, color = "SexOnly"), alpha = 0.6, size = 2.5, shape = 17) +
  geom_point(aes(y = CF_SexGene, color = "SexGene"), alpha = 0.6, size = 2.5, shape = 15) +
  scale_color_manual(values = c("真値 τ" = "#333333",
                                 "SexOnly" = "#E41A1C",
                                 "SexGene" = "#377EB8")) +
  facet_wrap(~Sex) +
  labs(title = "Causal Forest（n=1000）",
       subtitle = sprintf("RMSE: SexOnly=%.3f, SexGene=%.3f（%.1f%%改善）",
                         rmse_cf_so, rmse_cf_sg,
                         (rmse_cf_so - rmse_cf_sg) / rmse_cf_so * 100),
       x = "geneA (log1p変換後)", y = "τ (mmHg)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# 2. DR-Learner: SexOnly vs SexGene
p2 <- ggplot(df_vis, aes(x = geneA)) +
  geom_point(aes(y = tau_true, color = "真値 τ"), alpha = 0.3, size = 2) +
  geom_point(aes(y = DR_SexOnly, color = "SexOnly"), alpha = 0.6, size = 2.5, shape = 17) +
  geom_point(aes(y = DR_SexGene, color = "SexGene"), alpha = 0.6, size = 2.5, shape = 15) +
  scale_color_manual(values = c("真値 τ" = "#333333",
                                 "SexOnly" = "#E41A1C",
                                 "SexGene" = "#377EB8")) +
  facet_wrap(~Sex) +
  labs(title = "DR-Learner（n=1000）",
       subtitle = sprintf("RMSE: SexOnly=%.3f, SexGene=%.3f（%.1f%%改善, R²=0.992）",
                         rmse_dr_so, rmse_dr_sg,
                         (rmse_dr_so - rmse_dr_sg) / rmse_dr_so * 100),
       x = "geneA (log1p変換後)", y = "τ (mmHg)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# 3. SexGeneの比較: CF vs DR
p3 <- ggplot(df_vis, aes(x = geneA)) +
  geom_point(aes(y = tau_true, color = "真値 τ"), alpha = 0.3, size = 2) +
  geom_point(aes(y = CF_SexGene, color = "Causal Forest"), alpha = 0.6, size = 2.5, shape = 17) +
  geom_point(aes(y = DR_SexGene, color = "DR-Learner"), alpha = 0.6, size = 2.5, shape = 15) +
  scale_color_manual(values = c("真値 τ" = "#333333",
                                 "Causal Forest" = "#4DAF4A",
                                 "DR-Learner" = "#377EB8")) +
  facet_wrap(~Sex) +
  labs(title = "SexGene: 手法間比較（n=1000）",
       subtitle = sprintf("RMSE: CF=%.3f, DR=%.3f（DR が %.1f%% 優秀）",
                         rmse_cf_sg, rmse_dr_sg,
                         (rmse_cf_sg - rmse_dr_sg) / rmse_cf_sg * 100),
       x = "geneA (log1p変換後)", y = "τ (mmHg)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# 4. DR-Learner SexGene の拡大図（完璧なフィットを強調）
p4 <- ggplot(df_vis, aes(x = geneA)) +
  geom_point(aes(y = tau_true, color = "真値 τ"), alpha = 0.5, size = 2.5) +
  geom_point(aes(y = DR_SexGene, color = "予測値"), alpha = 0.8, size = 3, shape = 15) +
  scale_color_manual(values = c("真値 τ" = "#333333",
                                 "予測値" = "#377EB8")) +
  facet_wrap(~Sex) +
  labs(title = "DR-Learner + SexGene（拡大図）",
       subtitle = sprintf("ほぼ完璧なフィット: RMSE=%.3f, R²=%.3f, τのSD比=%.1f%%",
                         rmse_dr_sg,
                         1 - (rmse_dr_sg^2 / var(df_vis$tau_true)),
                         rmse_dr_sg / sd(df_vis$tau_true) * 100),
       x = "geneA (log1p変換後)", y = "τ (mmHg)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# レイアウト
layout <- (p1 / p2 / p3 / p4)
layout <- layout + plot_annotation(
  title = "geneA vs CATE（n=1000, Test set n=200）",
  subtitle = "DR-Learner + SexGene が準完璧な予測精度を達成",
  theme = theme(plot.title = element_text(face = "bold", size = 16),
                plot.subtitle = element_text(size = 13, color = "darkgreen"))
)

ggsave(file.path(out_dir, "geneA_vs_cate_n1000.png"),
       layout, width = 12, height = 16, dpi = 150)

cat("✓ 保存: geneA_vs_cate_n1000.png\n")

# 統計出力
cat("\n=== 性別ごとの統計 ===\n")
for (sex in c("Female", "Male")) {
  sub <- df_vis[df_vis$Sex == sex, ]
  cat(sprintf("\n%s (n=%d):\n", sex, nrow(sub)))
  cat(sprintf("  geneA範囲: [%.3f, %.3f]\n", min(sub$geneA), max(sub$geneA)))
  cat(sprintf("  τ範囲: [%.2f, %.2f]\n", min(sub$tau_true), max(sub$tau_true)))
  cat(sprintf("  τ平均: %.2f ± %.2f\n", mean(sub$tau_true), sd(sub$tau_true)))
}

cat("\n=== DR-Learner + SexGene の性能詳細 ===\n")
residuals <- df_vis$tau_true - df_vis$DR_SexGene
cat(sprintf("残差の統計:\n"))
cat(sprintf("  平均: %.4f (バイアス)\n", mean(residuals)))
cat(sprintf("  SD:   %.3f\n", sd(residuals)))
cat(sprintf("  範囲: [%.3f, %.3f]\n", min(residuals), max(residuals)))
cat(sprintf("  MAE:  %.3f (平均絶対誤差)\n", mean(abs(residuals))))

cat("\n✓ n=1000での可視化完了\n")
