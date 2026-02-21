#!/usr/bin/env Rscript
# visualize_comprehensive.R
# 修正版DGPの包括的可視化

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/demov04/results_corrected"

# ============================================================
# CONFIG & DGP
# ============================================================
config <- list(
  seed_base = 20260215L,
  n_val = 500L,
  effect_type = "linear",
  train_frac = 0.8,
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0, mu0 = 120.0, sigma_y = 5.0
)

generate_full_data_corrected <- function(n, effect_type, seed, cfg) {
  set.seed(seed)
  S <- rbinom(n, 1, 0.5)
  Sc <- S - 0.5
  E_true <- rnorm(n, cfg$m_E, cfg$s_E)
  mu_nb <- exp(E_true)
  C_count <- rnbinom(n, size = cfg$theta_nb, mu = mu_nb)
  G_raw <- log1p(C_count)
  delta <- rnorm(n, 0, cfg$sigma_G)
  G_obs <- G_raw + delta
  T_vec <- rbinom(n, 1, 0.5)
  Z_obs <- (G_obs - mean(G_obs)) / sd(G_obs)

  tau_indiv <- cfg$a_S * Sc + cfg$a_G * Z_obs
  tau_raw <- cfg$tau_base + tau_indiv
  tau <- pmin(cfg$tau_max, pmax(cfg$tau_min, tau_raw))

  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y1 <- Y0 + tau
  Y_obs <- ifelse(T_vec == 1, Y1, Y0)

  data.frame(
    id = 1:n, S = S, geneA = G_obs, T = T_vec,
    tau = tau, Y0 = Y0, Y1 = Y1, Y_obs = Y_obs
  )
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
# データ生成
# ============================================================
seed_cond <- 20260425  # 元と同じseed
dat <- generate_full_data_corrected(config$n_val, config$effect_type, seed_cond, config)
spl <- stratified_split(dat$T, config$train_frac, seed_cond)

dat$Sex <- factor(ifelse(dat$S == 1, "Male", "Female"), levels = c("Female", "Male"))
dat$T_label <- factor(ifelse(dat$T == 1, "治療群 (T=1)", "対照群 (T=0)"),
                      levels = c("対照群 (T=0)", "治療群 (T=1)"))
dat$Split <- ifelse(seq_len(nrow(dat)) %in% spl$train, "Train", "Test")

# 予測結果読み込み
df_pred <- read.csv(file.path(out_dir, "predictions_corrected.csv"))

# ============================================================
# 1. アウトカムの分布
# ============================================================

# 1a. 観測アウトカム Y の分布（T別）
p1a <- ggplot(dat, aes(x = Y_obs, fill = T_label)) +
  geom_density(alpha = 0.6) +
  geom_vline(data = aggregate(Y_obs ~ T_label, dat, mean),
             aes(xintercept = Y_obs, color = T_label),
             linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("対照群 (T=0)" = "#4DAF4A", "治療群 (T=1)" = "#E41A1C")) +
  scale_color_manual(values = c("対照群 (T=0)" = "#4DAF4A", "治療群 (T=1)" = "#E41A1C"),
                     guide = "none") +
  labs(title = "観測アウトカム Y の分布",
       subtitle = sprintf("対照群: %.1f mmHg | 治療群: %.1f mmHg | 差分: %.1f mmHg",
                         mean(dat$Y_obs[dat$T==0]), mean(dat$Y_obs[dat$T==1]),
                         mean(dat$Y_obs[dat$T==1]) - mean(dat$Y_obs[dat$T==0])),
       x = "血圧 (mmHg)", y = "密度", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 1b. Y0 と Y1 の分布
df_pot <- rbind(
  data.frame(Y = dat$Y0, Type = "Y(0) 未治療", stringsAsFactors = FALSE),
  data.frame(Y = dat$Y1, Type = "Y(1) 治療", stringsAsFactors = FALSE)
)
df_pot$Type <- factor(df_pot$Type, levels = c("Y(0) 未治療", "Y(1) 治療"))

p1b <- ggplot(df_pot, aes(x = Y, fill = Type)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Y(0) 未治療" = "#4DAF4A", "Y(1) 治療" = "#E41A1C")) +
  labs(title = "潜在的アウトカムの分布",
       subtitle = sprintf("Y(0): %.1f mmHg | Y(1): %.1f mmHg | ATE: %.1f mmHg",
                         mean(dat$Y0), mean(dat$Y1), mean(dat$Y1 - dat$Y0)),
       x = "血圧 (mmHg)", y = "密度", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# ============================================================
# 2. 説明変数の分布
# ============================================================

# 2a. Sex の分布（T別）
sex_counts <- as.data.frame(table(dat$T_label, dat$Sex))
names(sex_counts) <- c("Treatment", "Sex", "Count")
sex_counts$Prop <- ave(sex_counts$Count, sex_counts$Treatment, FUN = function(x) x/sum(x))

p2a <- ggplot(sex_counts, aes(x = Treatment, y = Prop, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Prop*100)),
            position = position_dodge(0.9), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.7)) +
  labs(title = "性別の分布（治療群別）", x = NULL, y = "割合", fill = "性別") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 2b. geneA の分布（T別）
p2b <- ggplot(dat, aes(x = geneA, fill = T_label)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("対照群 (T=0)" = "#4DAF4A", "治療群 (T=1)" = "#E41A1C")) +
  labs(title = "geneA 発現量の分布（治療群別）",
       subtitle = sprintf("対照群: %.2f±%.2f | 治療群: %.2f±%.2f",
                         mean(dat$geneA[dat$T==0]), sd(dat$geneA[dat$T==0]),
                         mean(dat$geneA[dat$T==1]), sd(dat$geneA[dat$T==1])),
       x = "geneA (log1p count)", y = "密度", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# ============================================================
# 3. CATE τ の分布
# ============================================================

# 3a. τ の全体分布
p3a <- ggplot(dat, aes(x = tau)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#984EA3", alpha = 0.7, color = "black") +
  geom_density(color = "blue", linewidth = 1.2) +
  geom_vline(xintercept = mean(dat$tau), color = "red",
             linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = mean(dat$tau) - 0.5, y = 0.2,
           label = sprintf("平均 = %.2f\nSD = %.2f", mean(dat$tau), sd(dat$tau)),
           hjust = 1, size = 4) +
  labs(title = "治療効果 τ の分布",
       subtitle = "τ < 0 は降圧効果（血圧が下がる）",
       x = "τ (mmHg)", y = "密度") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# 3b. τ の分布（Sex × geneA別）
dat$geneA_group <- cut(dat$geneA,
                       breaks = quantile(dat$geneA, c(0, 0.33, 0.66, 1)),
                       labels = c("Low", "Medium", "High"),
                       include.lowest = TRUE)

p3b <- ggplot(dat, aes(x = geneA_group, y = tau, fill = Sex)) +
  geom_boxplot(alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
               position = position_dodge(0.75), color = "red") +
  scale_fill_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
  labs(title = "τ の分布（性別 × 遺伝子発現別）",
       subtitle = "赤いダイヤ: 平均値",
       x = "geneA 発現レベル", y = "τ (mmHg)", fill = "性別") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# ============================================================
# 4. 予測値 vs 真値
# ============================================================

# Testデータのみ
dat_test <- dat[spl$test, ]

# 4つの条件の予測値を結合
plot_data_list <- list()

for (meth in c("causal_forest", "drlearner")) {
  for (covset in c("SexOnly", "SexGene")) {
    sub <- subset(df_pred, method == meth & covset == covset)

    plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
      tau_true = sub$tau_true,
      tau_hat = sub$tau_hat,
      Method = ifelse(meth == "causal_forest", "Causal Forest", "DR-Learner"),
      Covset = covset,
      stringsAsFactors = FALSE
    )
  }
}

df_plot <- do.call(rbind, plot_data_list)
df_plot <- df_plot[!is.na(df_plot$tau_hat), ]
df_plot$Condition <- paste(df_plot$Method, df_plot$Covset, sep = " + ")

# 4a. 予測値 vs 真値の散布図（4パネル）
p4a <- ggplot(df_plot, aes(x = tau_true, y = tau_hat)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(alpha = 0.4, size = 1.5, color = "#377EB8") +
  facet_grid(Method ~ Covset) +
  labs(title = "予測値 vs 真値",
       subtitle = "赤線: 完全予測 (τ_hat = τ_true)",
       x = "真値 τ (mmHg)", y = "予測値 τ_hat (mmHg)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

# 4b. 残差プロット
df_plot$residual <- df_plot$tau_hat - df_plot$tau_true

p4b <- ggplot(df_plot, aes(x = tau_true, y = residual)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(alpha = 0.4, size = 1.5, color = "#E41A1C") +
  facet_grid(Method ~ Covset) +
  labs(title = "残差プロット",
       subtitle = "残差 = τ_hat - τ_true",
       x = "真値 τ (mmHg)", y = "残差 (mmHg)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

# ============================================================
# 5. RMSE比較
# ============================================================

# RMSEを計算
rmse_data <- aggregate(residual ~ Method + Covset, df_plot,
                      function(x) sqrt(mean(x^2)))
names(rmse_data)[3] <- "RMSE"
rmse_data$Condition <- paste(rmse_data$Method, rmse_data$Covset, sep = "\n")

p5 <- ggplot(rmse_data, aes(x = Condition, y = RMSE, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", RMSE)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Causal Forest" = "#4DAF4A", "DR-Learner" = "#377EB8")) +
  labs(title = "CATE推定精度（RMSE）",
       subtitle = "低いほど良い",
       x = NULL, y = "RMSE (mmHg)", fill = "手法") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# ============================================================
# 6. geneA vs Y の関係
# ============================================================

p6 <- ggplot(dat, aes(x = geneA, y = Y_obs, color = T_label, shape = Sex)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("対照群 (T=0)" = "#4DAF4A", "治療群 (T=1)" = "#E41A1C")) +
  labs(title = "geneA vs 観測アウトカム Y",
       subtitle = "治療群の方が血圧が低い（降圧効果）",
       x = "geneA (log1p count)", y = "血圧 Y (mmHg)",
       color = NULL, shape = "性別") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# ============================================================
# レイアウト
# ============================================================

# ページ1: アウトカムと説明変数
layout1 <- (p1a | p1b) / (p2a | p2b)
layout1 <- layout1 + plot_annotation(
  title = "1. アウトカムと説明変数の分布（修正版DGP, n=500）",
  theme = theme(plot.title = element_text(face = "bold", size = 15))
)

# ページ2: CATE分布と予測性能
layout2 <- (p3a | p3b) / (p4a)
layout2 <- layout2 + plot_annotation(
  title = "2. CATE分布と予測値 vs 真値",
  theme = theme(plot.title = element_text(face = "bold", size = 15))
)

# ページ3: 残差とRMSE
layout3 <- p4b / (p5 | p6)
layout3 <- layout3 + plot_annotation(
  title = "3. 残差分析とRMSE",
  theme = theme(plot.title = element_text(face = "bold", size = 15))
)

# 保存
ggsave(file.path(out_dir, "comprehensive_page1.png"), layout1,
       width = 14, height = 10, dpi = 150)
ggsave(file.path(out_dir, "comprehensive_page2.png"), layout2,
       width = 14, height = 10, dpi = 150)
ggsave(file.path(out_dir, "comprehensive_page3.png"), layout3,
       width = 14, height = 10, dpi = 150)

cat("Saved: comprehensive_page1.png (アウトカムと説明変数)\n")
cat("Saved: comprehensive_page2.png (CATE分布と予測値)\n")
cat("Saved: comprehensive_page3.png (残差とRMSE)\n")

# ============================================================
# 統計サマリー出力
# ============================================================

cat("\n=== データ統計サマリー（修正版DGP, n=500, seed=20260425） ===\n\n")

cat("【アウトカム】\n")
cat(sprintf("  Y(T=0): %.2f ± %.2f mmHg (n=%d)\n",
            mean(dat$Y_obs[dat$T==0]), sd(dat$Y_obs[dat$T==0]), sum(dat$T==0)))
cat(sprintf("  Y(T=1): %.2f ± %.2f mmHg (n=%d)\n",
            mean(dat$Y_obs[dat$T==1]), sd(dat$Y_obs[dat$T==1]), sum(dat$T==1)))
cat(sprintf("  差分: %.2f mmHg（治療により血圧が下がる）✓\n",
            mean(dat$Y_obs[dat$T==1]) - mean(dat$Y_obs[dat$T==0])))

cat("\n【潜在的アウトカム】\n")
cat(sprintf("  Y(0): %.2f ± %.2f mmHg\n", mean(dat$Y0), sd(dat$Y0)))
cat(sprintf("  Y(1): %.2f ± %.2f mmHg\n", mean(dat$Y1), sd(dat$Y1)))
cat(sprintf("  ATE: %.2f mmHg\n", mean(dat$Y1 - dat$Y0)))

cat("\n【説明変数】\n")
cat(sprintf("  Sex (Male比率): T=0で%.1f%%, T=1で%.1f%%\n",
            mean(dat$S[dat$T==0])*100, mean(dat$S[dat$T==1])*100))
cat(sprintf("  geneA: T=0で%.3f±%.3f, T=1で%.3f±%.3f\n",
            mean(dat$geneA[dat$T==0]), sd(dat$geneA[dat$T==0]),
            mean(dat$geneA[dat$T==1]), sd(dat$geneA[dat$T==1])))

cat("\n【治療効果 τ】\n")
cat(sprintf("  平均: %.2f mmHg\n", mean(dat$tau)))
cat(sprintf("  SD:   %.2f mmHg\n", sd(dat$tau)))
cat(sprintf("  範囲: [%.2f, %.2f] mmHg\n", min(dat$tau), max(dat$tau)))

cat("\n【CATE推定性能（RMSE）】\n")
for (i in seq_len(nrow(rmse_data))) {
  cat(sprintf("  %s + %s: %.3f\n",
              rmse_data$Method[i], rmse_data$Covset[i], rmse_data$RMSE[i]))
}

cat("\n【geneA追加の効果】\n")
rmse_cf_so <- rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexOnly"]
rmse_cf_sg <- rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexGene"]
rmse_dr_so <- rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexOnly"]
rmse_dr_sg <- rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexGene"]

cat(sprintf("  Causal Forest: %.3f → %.3f (%.1f%%改善)\n",
            rmse_cf_so, rmse_cf_sg, 100*(rmse_cf_so - rmse_cf_sg)/rmse_cf_so))
cat(sprintf("  DR-Learner:    %.3f → %.3f (%.1f%%改善)\n",
            rmse_dr_so, rmse_dr_sg, 100*(rmse_dr_so - rmse_dr_sg)/rmse_dr_so))

cat("\n✓ 包括的可視化が完了しました\n")
