#!/usr/bin/env Rscript
# visualize_raw_counts.R
# カウントデータ版の包括的可視化

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/demov04/results_raw_counts"

config <- list(
  seed_base = 20260215L,
  n_val = 500L,
  train_frac = 0.8,
  m_E = log(50), s_E = 0.35, theta_nb = 30,
  tau_base = -12.0, a_S = 2.0, a_G = 0.05,
  tau_min = -18.0, tau_max = -6.0,
  mu0 = 120.0, sigma_y = 5.0
)

generate_full_data_raw_counts <- function(n, seed, cfg) {
  set.seed(seed)
  S <- rbinom(n, 1, 0.5)
  Sc <- S - 0.5
  E_true <- rnorm(n, cfg$m_E, cfg$s_E)
  mu_nb <- exp(E_true)
  C_count <- rnbinom(n, size = cfg$theta_nb, mu = mu_nb)
  geneA <- C_count  # カウントデータそのまま
  T_vec <- rbinom(n, 1, 0.5)

  Z_obs <- (geneA - mean(geneA)) / sd(geneA)
  tau_indiv <- cfg$a_S * Sc + cfg$a_G * Z_obs
  tau_raw <- cfg$tau_base + tau_indiv
  tau <- pmin(cfg$tau_max, pmax(cfg$tau_min, tau_raw))

  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y1 <- Y0 + tau
  Y_obs <- ifelse(T_vec == 1, Y1, Y0)

  data.frame(
    id = 1:n, S = S, geneA = geneA, T = T_vec,
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

seed_cond <- 20260425
dat <- generate_full_data_raw_counts(config$n_val, seed_cond, config)
spl <- stratified_split(dat$T, config$train_frac, seed_cond)

dat$Sex <- factor(ifelse(dat$S == 1, "Male", "Female"), levels = c("Female", "Male"))
dat$T_label <- factor(ifelse(dat$T == 1, "治療群", "対照群"),
                      levels = c("対照群", "治療群"))

df_pred <- read.csv(file.path(out_dir, "predictions_raw_counts.csv"))

# ============================================================
# プロット作成
# ============================================================

# 1. geneA（カウント）の分布
p1 <- ggplot(dat, aes(x = geneA)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#4DAF4A", alpha = 0.7, color = "black") +
  geom_density(color = "blue", linewidth = 1.2) +
  labs(title = "geneA カウントデータの分布",
       subtitle = sprintf("範囲: [%d, %d] カウント | 平均: %.1f | SD: %.1f",
                         min(dat$geneA), max(dat$geneA),
                         mean(dat$geneA), sd(dat$geneA)),
       x = "geneA (RNA-seq カウント)", y = "密度") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# 2. geneA vs τ の関係
p2 <- ggplot(dat, aes(x = geneA, y = tau, color = Sex)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
  labs(title = "geneA vs τ（真値）",
       subtitle = sprintf("係数 a_G = %.3f（カウントに対する効果）", config$a_G),
       x = "geneA (カウント)", y = "τ (mmHg)", color = "性別") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 3. 予測値 vs 真値（4パネル）
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

p3 <- ggplot(df_plot, aes(x = tau_true, y = tau_hat)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.4, size = 1.5, color = "#377EB8") +
  facet_grid(Method ~ Covset) +
  labs(title = "予測値 vs 真値（カウントデータ）",
       x = "真値 τ", y = "予測値 τ_hat") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

# 4. RMSE比較
df_plot$residual <- df_plot$tau_hat - df_plot$tau_true
rmse_data <- aggregate(residual ~ Method + Covset, df_plot,
                      function(x) sqrt(mean(x^2)))
names(rmse_data)[3] <- "RMSE"
rmse_data$Condition <- paste(rmse_data$Method, rmse_data$Covset, sep = "\n")

p4 <- ggplot(rmse_data, aes(x = Condition, y = RMSE, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", RMSE)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Causal Forest" = "#4DAF4A", "DR-Learner" = "#377EB8")) +
  labs(title = "RMSE比較（カウントデータ）",
       subtitle = "★ SexGene追加で性能が悪化！",
       x = NULL, y = "RMSE", fill = "手法") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 5. geneA vs 予測値（DR-Learner）
dat_test <- dat[spl$test, ]
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
  labs(title = "geneA vs CATE（DR-Learner）",
       subtitle = "カウントデータではSexGeneの予測がばらつく",
       x = "geneA (カウント)", y = "τ (mmHg)", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# レイアウト
layout <- (p1 | p2) / (p3 | p4) / p5
layout <- layout + plot_annotation(
  title = "カウントデータ（log変換なし）での結果（n=500, a_G=0.05）",
  subtitle = "問題: geneA追加で性能が悪化している",
  theme = theme(plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12, color = "red"))
)

ggsave(file.path(out_dir, "comprehensive_raw_counts.png"), layout,
       width = 14, height = 14, dpi = 150)

cat("Saved: comprehensive_raw_counts.png\n")

# 統計出力
cat("\n=== カウントデータ版の統計 ===\n\n")
cat(sprintf("geneA (カウント):\n"))
cat(sprintf("  範囲: [%d, %d]\n", min(dat$geneA), max(dat$geneA)))
cat(sprintf("  平均: %.1f ± %.1f\n", mean(dat$geneA), sd(dat$geneA)))
cat(sprintf("  中央値: %.0f\n", median(dat$geneA)))

cat(sprintf("\nτ (治療効果):\n"))
cat(sprintf("  平均: %.2f mmHg\n", mean(dat$tau)))
cat(sprintf("  SD: %.2f mmHg\n", sd(dat$tau)))
cat(sprintf("  範囲: [%.2f, %.2f]\n", min(dat$tau), max(dat$tau)))

cat(sprintf("\nRMSE:\n"))
for (i in seq_len(nrow(rmse_data))) {
  cat(sprintf("  %s + %s: %.3f\n", rmse_data$Method[i], rmse_data$Covset[i], rmse_data$RMSE[i]))
}

cat(sprintf("\ngeneA追加の影響:\n"))
rmse_cf_diff <- rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexGene"] -
                rmse_data$RMSE[rmse_data$Method == "Causal Forest" & rmse_data$Covset == "SexOnly"]
rmse_dr_diff <- rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexGene"] -
                rmse_data$RMSE[rmse_data$Method == "DR-Learner" & rmse_data$Covset == "SexOnly"]

cat(sprintf("  Causal Forest: %+.3f（%s）\n", rmse_cf_diff,
            ifelse(rmse_cf_diff > 0, "悪化", "改善")))
cat(sprintf("  DR-Learner:    %+.3f（%s）\n", rmse_dr_diff,
            ifelse(rmse_dr_diff > 0, "悪化", "改善")))

cat("\n✓ カウントデータ版の可視化が完了しました\n")
cat("\n注意: a_G=0.05 では効果が小さすぎる可能性があります\n")
cat("     geneAのSDが約25なので、τへの寄与は 0.05 × 25 = 1.25 mmHg程度\n")
