#!/usr/bin/env Rscript
# check_covariate_balance.R
# 治療群と対照群での共変量分布の確認

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/demov04/results_with_fix"

# ============================================================
# CONFIG
# ============================================================
config <- list(
  seed_base  = 20260215L,
  n_values   = c(100L, 500L, 1000L),
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5, delta_thr = 1.5,
  tau_min = -18.0, tau_max = -6.0, mu0 = 120.0, sigma_y = 5.0
)

# ============================================================
# DGP
# ============================================================
generate_data_v04 <- function(n, effect_type, seed, cfg) {
  set.seed(seed)
  S  <- rbinom(n, 1, 0.5)
  Sc <- S - 0.5
  E_true  <- rnorm(n, cfg$m_E, cfg$s_E)
  mu_nb   <- exp(E_true)
  C_count <- rnbinom(n, size = cfg$theta_nb, mu = mu_nb)
  G_raw   <- log1p(C_count)
  delta   <- rnorm(n, 0, cfg$sigma_G)
  G_obs   <- G_raw + delta
  T_vec   <- rbinom(n, 1, 0.5)

  Z_obs <- (G_obs - mean(G_obs)) / sd(G_obs)

  if (effect_type == "linear") {
    tau_indiv <- cfg$a_S * Sc + cfg$a_G * Z_obs
  } else {
    q1 <- quantile(Z_obs, 0.33)
    q2 <- quantile(Z_obs, 0.66)
    g_Z <- ifelse(Z_obs <= q1, -cfg$delta_thr,
                  ifelse(Z_obs <= q2, 0, cfg$delta_thr))
    tau_indiv <- cfg$a_S * Sc + g_Z
  }

  tau_raw <- cfg$tau_base + tau_indiv
  tau <- pmin(cfg$tau_max, pmax(cfg$tau_min, tau_raw))
  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y <- Y0 - T_vec * tau

  data.frame(
    id = 1:n,
    S = S,
    geneA = G_obs,
    T = T_vec,
    tau = tau,
    Y = Y
  )
}

# ============================================================
# プロット作成
# ============================================================
plot_list <- list()

for (n_val in config$n_values) {
  seed_cond <- config$seed_base + which(config$n_values == n_val) * 100

  # linear効果で代表データ生成
  dat <- generate_data_v04(n_val, "linear", seed_cond, config)

  dat$T_label <- factor(ifelse(dat$T == 1, "治療群 (T=1)", "対照群 (T=0)"),
                        levels = c("対照群 (T=0)", "治療群 (T=1)"))
  dat$Sex_label <- factor(ifelse(dat$S == 1, "Male", "Female"),
                          levels = c("Female", "Male"))

  # --- Sex の分布（棒グラフ） ---
  sex_counts <- as.data.frame(table(dat$T_label, dat$Sex_label))
  names(sex_counts) <- c("Treatment", "Sex", "Count")

  # 割合計算
  sex_counts$Prop <- ave(sex_counts$Count, sex_counts$Treatment,
                         FUN = function(x) x / sum(x))

  p_sex <- ggplot(sex_counts, aes(x = Treatment, y = Prop, fill = Sex)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Prop*100, Count)),
              position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    labs(
      title = sprintf("性別の分布 (n=%d)", n_val),
      x = NULL,
      y = "割合",
      fill = "性別"
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 0.7)) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )

  # --- geneA の分布（密度プロット + 箱ひげ図） ---

  # 密度プロット
  p_geneA_dens <- ggplot(dat, aes(x = geneA, fill = T_label)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("対照群 (T=0)" = "#4DAF4A",
                                  "治療群 (T=1)" = "#984EA3")) +
    labs(
      title = sprintf("geneA の分布 (n=%d)", n_val),
      x = "geneA (log1p count + noise)",
      y = "密度",
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )

  # 箱ひげ図
  p_geneA_box <- ggplot(dat, aes(x = T_label, y = geneA, fill = T_label)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "red", color = "red") +
    scale_fill_manual(values = c("対照群 (T=0)" = "#4DAF4A",
                                  "治療群 (T=1)" = "#984EA3")) +
    labs(
      title = sprintf("geneA の分布 (箱ひげ図, n=%d)", n_val),
      x = NULL,
      y = "geneA",
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )

  # --- 統計量の計算 ---
  t0 <- subset(dat, T == 0)
  t1 <- subset(dat, T == 1)

  # Sex の比率
  sex_t0 <- mean(t0$S)
  sex_t1 <- mean(t1$S)
  sex_diff <- sex_t1 - sex_t0

  # geneA の平均・SD
  geneA_t0_mean <- mean(t0$geneA)
  geneA_t0_sd <- sd(t0$geneA)
  geneA_t1_mean <- mean(t1$geneA)
  geneA_t1_sd <- sd(t1$geneA)
  geneA_diff <- geneA_t1_mean - geneA_t0_mean

  # t検定
  sex_test <- prop.test(c(sum(t1$S), sum(t0$S)), c(nrow(t1), nrow(t0)))
  geneA_test <- t.test(t1$geneA, t0$geneA)

  stats_text <- sprintf(
    "性別 (Male比率):\n  T=0: %.1f%% | T=1: %.1f%% | 差: %.1f%% (p=%.3f)\n\ngeneA (平均±SD):\n  T=0: %.3f±%.3f | T=1: %.3f±%.3f\n  差: %.3f (p=%.3f)",
    sex_t0*100, sex_t1*100, sex_diff*100, sex_test$p.value,
    geneA_t0_mean, geneA_t0_sd, geneA_t1_mean, geneA_t1_sd,
    geneA_diff, geneA_test$p.value
  )

  p_stats <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = stats_text,
             hjust = 0.5, vjust = 0.5, size = 3.5, family = "mono") +
    theme_void() +
    labs(title = sprintf("統計量 (n=%d)", n_val)) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

  # --- レイアウト ---
  layout <- (p_sex | p_geneA_dens) / (p_geneA_box | p_stats)
  layout <- layout +
    plot_annotation(
      title = sprintf("共変量のバランス確認 (n=%d, seed=%d)", n_val, seed_cond),
      subtitle = "ランダム化により治療群と対照群で共変量分布は同じはず",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey40")
      )
    )

  # 保存
  fname <- sprintf("plot_covariate_balance_n%d.png", n_val)
  ggsave(file.path(out_dir, fname), layout, width = 12, height = 8, dpi = 150)
  cat(sprintf("Saved: %s\n", fname))

  # 統計情報を出力
  cat(sprintf("\n=== n=%d ===\n", n_val))
  cat(sprintf("対照群: n=%d | 治療群: n=%d\n", nrow(t0), nrow(t1)))
  cat(stats_text, "\n\n")
}

cat("\nAll covariate balance plots saved.\n")
