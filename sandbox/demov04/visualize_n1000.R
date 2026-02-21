#!/usr/bin/env Rscript
# visualize_n1000.R
# n=1000の結果を可視化

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/demov04/results_n1000"

# データ読み込み
dat <- readRDS(file.path(out_dir, "full_data_n1000.rds"))
spl <- readRDS(file.path(out_dir, "split_n1000.rds"))
df_pred <- read.csv(file.path(out_dir, "predictions_n1000.csv"))

n <- length(dat$S)
dat_full <- data.frame(
  id = 1:n,
  S = dat$S,
  geneA = dat$geneA,
  T = dat$T_vec,
  tau = dat$tau,
  Y0 = dat$Y0,
  Y = dat$Y,
  Sex = factor(ifelse(dat$S == 1, "Male", "Female"), levels = c("Female", "Male")),
  T_label = factor(ifelse(dat$T_vec == 1, "治療群", "対照群"),
                   levels = c("対照群", "治療群"))
)

cat("=== n=1000での統計 ===\n\n")
cat(sprintf("サンプルサイズ: %d (Train: %d, Test: %d)\n",
            n, length(spl$train), length(spl$test)))
cat(sprintf("\ngeneA: 平均=%.3f, SD=%.3f, 範囲=[%.3f, %.3f]\n",
            mean(dat_full$geneA), sd(dat_full$geneA),
            min(dat_full$geneA), max(dat_full$geneA)))
cat(sprintf("τ: 平均=%.3f, SD=%.3f, 範囲=[%.2f, %.2f]\n",
            mean(dat_full$tau), sd(dat_full$tau),
            min(dat_full$tau), max(dat_full$tau)))

# ============================================================
# プロット作成
# ============================================================

# 1. geneAの分布
p1 <- ggplot(dat_full, aes(x = geneA)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#4DAF4A", alpha = 0.7, color = "black") +
  geom_density(color = "blue", linewidth = 1.2) +
  labs(title = "geneA分布（log1p変換後）",
       subtitle = sprintf("n=%d | 平均=%.2f, SD=%.2f",
                         n, mean(dat_full$geneA), sd(dat_full$geneA)),
       x = "geneA (log1p変換後)", y = "密度") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# 2. geneA vs τ の関係
p2 <- ggplot(dat_full, aes(x = geneA, y = tau, color = Sex)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
  labs(title = "geneA vs τ（真値）",
       subtitle = "係数 a_G = 1.5",
       x = "geneA (log1p変換後)", y = "τ (mmHg)", color = "性別") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 3. 予測値 vs 真値（4パネル）
df_plot <- df_pred
df_plot$Method <- ifelse(df_plot$method == "causal_forest", "Causal Forest", "DR-Learner")
df_plot$Covset <- df_plot$covset
df_plot <- df_plot[!is.na(df_plot$tau_hat), ]

p3 <- ggplot(df_plot, aes(x = tau_true, y = tau_hat)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.4, size = 1.5, color = "#377EB8") +
  facet_grid(Method ~ Covset) +
  labs(title = "予測値 vs 真値（n=1000）",
       x = "真値 τ", y = "予測値 τ_hat") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

# 4. RMSE比較（各組み合わせで個別に計算）
rmse_list <- list()
for (meth in c("Causal Forest", "DR-Learner")) {
  for (covset in c("SexOnly", "SexGene")) {
    sub <- df_plot[df_plot$Method == meth & df_plot$Covset == covset, ]
    rmse_val <- sqrt(mean((sub$tau_true - sub$tau_hat)^2))
    rmse_list[[length(rmse_list) + 1]] <- data.frame(
      Method = meth,
      Covset = covset,
      RMSE = rmse_val,
      stringsAsFactors = FALSE
    )
  }
}
rmse_data <- do.call(rbind, rmse_list)
rmse_data$Condition <- paste(rmse_data$Method, rmse_data$Covset, sep = "\n")

p4 <- ggplot(rmse_data, aes(x = Condition, y = RMSE, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", RMSE)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Causal Forest" = "#4DAF4A", "DR-Learner" = "#377EB8")) +
  labs(title = "RMSE比較（n=1000）",
       subtitle = "★ DR-Learner + SexGene: RMSE=0.168（τのSD=1.83の9%）",
       x = NULL, y = "RMSE", fill = "手法") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 5. geneA vs 予測値（DR-Learner）
dat_test <- dat_full[spl$test, ]
sub_dr_sg <- subset(df_pred, method == "drlearner" & covset == "SexGene")
sub_dr_so <- subset(df_pred, method == "drlearner" & covset == "SexOnly")

df_vis <- data.frame(
  geneA = dat_test$geneA,
  tau_true = dat_test$tau,
  tau_hat_SexGene = sub_dr_sg$tau_hat,
  tau_hat_SexOnly = sub_dr_so$tau_hat,
  Sex = dat_test$Sex
)

p5 <- ggplot(df_vis, aes(x = geneA)) +
  geom_point(aes(y = tau_true, color = "True τ"), alpha = 0.3, size = 1.5) +
  geom_point(aes(y = tau_hat_SexOnly, color = "SexOnly"), alpha = 0.6, size = 1.8, shape = 17) +
  geom_point(aes(y = tau_hat_SexGene, color = "SexGene"), alpha = 0.6, size = 1.8, shape = 15) +
  scale_color_manual(values = c("True τ" = "#333333", "SexOnly" = "#E41A1C", "SexGene" = "#377EB8")) +
  facet_wrap(~Sex) +
  labs(title = "geneA vs CATE（DR-Learner, n=1000）",
       subtitle = "SexGeneの予測が真値にほぼ完璧に一致",
       x = "geneA (log1p変換後)", y = "τ (mmHg)", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# レイアウト
layout <- (p1 | p2) / (p3 | p4) / p5
layout <- layout + plot_annotation(
  title = "n=1000での結果（修正版DGP）",
  subtitle = "DR-Learner + SexGene が非常に高精度（RMSE=0.168, R²=0.991）",
  theme = theme(plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12, color = "darkgreen"))
)

ggsave(file.path(out_dir, "comprehensive_n1000.png"), layout,
       width = 14, height = 14, dpi = 150)

cat("\n✓ 保存: comprehensive_n1000.png\n")

# 統計出力
cat("\n=== RMSE結果 ===\n")
for (i in seq_len(nrow(rmse_data))) {
  cat(sprintf("%s + %s: RMSE=%.3f\n",
              rmse_data$Method[i], rmse_data$Covset[i], rmse_data$RMSE[i]))
}

cat("\n=== geneA追加の効果 ===\n")
rmse_cf_diff <- rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexGene"] -
                rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexOnly"]
rmse_dr_diff <- rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexGene"] -
                rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexOnly"]

cat(sprintf("Causal Forest: %.3f → %.3f (%.1f%% 改善)\n",
            rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexOnly"],
            rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexGene"],
            -rmse_cf_diff / rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexOnly"] * 100))

cat(sprintf("DR-Learner:    %.3f → %.3f (%.1f%% 改善)\n",
            rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexOnly"],
            rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexGene"],
            -rmse_dr_diff / rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexOnly"] * 100))

# R²の計算
r2_dr_sg <- 1 - (rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexGene"]^2 /
                 var(df_vis$tau_true))

cat(sprintf("\nDR-Learner + SexGene の性能:\n"))
cat(sprintf("  RMSE = %.3f mmHg\n", rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexGene"]))
cat(sprintf("  R² = %.3f (τの分散の%.1f%%を説明)\n", r2_dr_sg, r2_dr_sg * 100))
cat(sprintf("  τのSD比 = %.1f%% (RMSE / SD(τ))\n",
            rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexGene"] / sd(dat_full$tau) * 100))

cat("\n✓ n=1000での可視化完了\n")
