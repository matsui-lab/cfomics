#!/usr/bin/env Rscript
# visualize_potential_outcomes.R
# Y(0), Y(1) の潜在的アウトカムを可視化

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
  n_val      = 500L,
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5, delta_thr = 1.5,
  tau_min = -18.0, tau_max = -6.0, mu0 = 120.0, sigma_y = 5.0
)

# ============================================================
# DGP
# ============================================================
generate_full_data <- function(n, effect_type, seed, cfg) {
  set.seed(seed)

  # 共変量
  S  <- rbinom(n, 1, 0.5)
  Sc <- S - 0.5

  # 遺伝子発現
  E_true  <- rnorm(n, cfg$m_E, cfg$s_E)
  mu_nb   <- exp(E_true)
  C_count <- rnbinom(n, size = cfg$theta_nb, mu = mu_nb)
  G_raw   <- log1p(C_count)
  delta   <- rnorm(n, 0, cfg$sigma_G)
  G_obs   <- G_raw + delta

  # 標準化
  Z_obs <- (G_obs - mean(G_obs)) / sd(G_obs)

  # 治療効果
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

  # 潜在的アウトカム
  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y1 <- Y0 - tau  # tau < 0 なので、治療により減少

  # 治療割付
  T_vec <- rbinom(n, 1, 0.5)

  # 観測アウトカム
  Y_obs <- ifelse(T_vec == 1, Y1, Y0)

  data.frame(
    id = 1:n,
    S = S,
    geneA = G_obs,
    T = T_vec,
    tau = tau,
    Y0 = Y0,
    Y1 = Y1,
    Y_obs = Y_obs,
    ITE = Y1 - Y0  # = -tau
  )
}

# ============================================================
# プロット関数
# ============================================================
plot_potential_outcomes <- function(dat, effect_type, seed) {

  dat$Sex <- factor(ifelse(dat$S == 1, "Male", "Female"),
                    levels = c("Female", "Male"))

  # --- プロット1: Y(0) と Y(1) の分布 ---
  df_y0 <- data.frame(Y = dat$Y0, Type = "Y(0) 未治療時", Sex = dat$Sex)
  df_y1 <- data.frame(Y = dat$Y1, Type = "Y(1) 治療時", Sex = dat$Sex)
  df_both <- rbind(df_y0, df_y1)
  df_both$Type <- factor(df_both$Type,
                         levels = c("Y(0) 未治療時", "Y(1) 治療時"))

  p1 <- ggplot(df_both, aes(x = Y, fill = Type)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Y(0) 未治療時" = "#4DAF4A",
                                  "Y(1) 治療時" = "#E41A1C")) +
    geom_vline(xintercept = mean(dat$Y0), color = "#4DAF4A",
               linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = mean(dat$Y1), color = "#E41A1C",
               linetype = "dashed", linewidth = 1) +
    labs(
      title = "潜在的アウトカムの分布",
      subtitle = sprintf("平均: Y(0)=%.1f, Y(1)=%.1f, ATE=%.1f",
                        mean(dat$Y0), mean(dat$Y1), mean(dat$Y1 - dat$Y0)),
      x = "血圧 (mmHg)",
      y = "密度",
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  # --- プロット2: Y(0) と Y(1) の散布図 ---
  p2 <- ggplot(dat, aes(x = Y0, y = Y1, color = Sex)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray50", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = mean(dat$ITE),
                color = "red", linewidth = 1) +
    scale_color_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    labs(
      title = "潜在的アウトカムの関係",
      subtitle = sprintf("点線: Y(1)=Y(0) | 赤線: Y(1)=Y(0)+ATE (ATE=%.1f)",
                        mean(dat$ITE)),
      x = "Y(0) 未治療時 (mmHg)",
      y = "Y(1) 治療時 (mmHg)",
      color = "性別"
    ) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  # --- プロット3: ITE (Individual Treatment Effect) の分布 ---
  p3 <- ggplot(dat, aes(x = ITE)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "#984EA3", alpha = 0.7, color = "black") +
    geom_density(color = "blue", linewidth = 1.2) +
    geom_vline(xintercept = mean(dat$ITE), color = "red",
               linetype = "dashed", linewidth = 1.2) +
    annotate("text", x = mean(dat$ITE) - 1, y = 0.2,
             label = sprintf("ATE = %.2f", mean(dat$ITE)),
             color = "red", size = 4.5, hjust = 1) +
    labs(
      title = "個人治療効果 (ITE) の分布",
      subtitle = sprintf("ITE = Y(1) - Y(0) = -τ | SD = %.2f", sd(dat$ITE)),
      x = "ITE (mmHg)",
      y = "密度"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # --- プロット4: 観測 vs 潜在 ---
  df_obs_t0 <- data.frame(
    Y = dat$Y_obs[dat$T == 0],
    Type = sprintf("観測 (T=0, n=%d)", sum(dat$T == 0)),
    stringsAsFactors = FALSE
  )
  df_obs_t1 <- data.frame(
    Y = dat$Y_obs[dat$T == 1],
    Type = sprintf("観測 (T=1, n=%d)", sum(dat$T == 1)),
    stringsAsFactors = FALSE
  )
  df_pot_y0 <- data.frame(
    Y = dat$Y0,
    Type = sprintf("潜在 Y(0) (n=%d)", nrow(dat)),
    stringsAsFactors = FALSE
  )
  df_pot_y1 <- data.frame(
    Y = dat$Y1,
    Type = sprintf("潜在 Y(1) (n=%d)", nrow(dat)),
    stringsAsFactors = FALSE
  )

  df_compare <- rbind(df_obs_t0, df_obs_t1, df_pot_y0, df_pot_y1)
  df_compare$Type <- factor(df_compare$Type,
                            levels = c(df_obs_t0$Type[1], df_obs_t1$Type[1],
                                      df_pot_y0$Type[1], df_pot_y1$Type[1]))

  p4 <- ggplot(df_compare, aes(x = Y, fill = Type)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = c(
      setNames("#4DAF4A", df_obs_t0$Type[1]),
      setNames("#E41A1C", df_obs_t1$Type[1]),
      setNames("#377EB8", df_pot_y0$Type[1]),
      setNames("#FF7F00", df_pot_y1$Type[1])
    )) +
    labs(
      title = "観測 vs 潜在的アウトカム",
      subtitle = "観測できるのは各個人でY(0)かY(1)のどちらか一方のみ",
      x = "血圧 (mmHg)",
      y = "密度",
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  # --- プロット5: τ（CATE）の分布 ---
  p5 <- ggplot(dat, aes(x = tau)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "#FF7F00", alpha = 0.7, color = "black") +
    geom_density(color = "blue", linewidth = 1.2) +
    geom_vline(xintercept = mean(dat$tau), color = "red",
               linetype = "dashed", linewidth = 1.2) +
    annotate("text", x = mean(dat$tau) + 0.5, y = 0.2,
             label = sprintf("平均τ = %.2f\nSD = %.2f",
                           mean(dat$tau), sd(dat$tau)),
             color = "red", size = 4, hjust = 0) +
    labs(
      title = "治療効果 τ (CATE) の分布",
      subtitle = "τ = -(Y(1) - Y(0)) = -ITE  (降圧効果なので負)",
      x = "τ (mmHg)",
      y = "密度"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # --- プロット6: Sex × geneA 別の τ ---
  dat$geneA_group <- cut(dat$geneA,
                         breaks = quantile(dat$geneA, c(0, 0.33, 0.66, 1)),
                         labels = c("Low", "Medium", "High"),
                         include.lowest = TRUE)

  p6 <- ggplot(dat, aes(x = geneA_group, y = tau, fill = Sex)) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23,
                 size = 3, position = position_dodge(0.75),
                 aes(fill = Sex), color = "red") +
    scale_fill_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8")) +
    labs(
      title = "τ の分布（性別 × 遺伝子発現別）",
      subtitle = "赤いダイヤ: 平均値",
      x = "geneA発現レベル",
      y = "τ (mmHg)",
      fill = "性別"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  # --- レイアウト ---
  layout <- (p1 | p2) / (p3 | p4) / (p5 | p6)
  layout <- layout +
    plot_annotation(
      title = sprintf("潜在的アウトカムの可視化 (n=%d, %s, seed=%d)",
                     nrow(dat), effect_type, seed),
      subtitle = sprintf("σ_y=%.1f (個人差のSD) | 全員のY(0)とY(1)を表示（因果推論では一方のみ観測可能）",
                        config$sigma_y),
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 11, color = "grey40")
      )
    )

  list(layout = layout, data = dat)
}

# ============================================================
# メイン実行
# ============================================================
for (effect_type in c("linear", "threshold3")) {
  seed_cond <- config$seed_base + 2 * 100 + ifelse(effect_type == "linear", 10, 20)

  cat(sprintf("\n=== %s effect (n=%d, seed=%d) ===\n",
              effect_type, config$n_val, seed_cond))

  dat <- generate_full_data(config$n_val, effect_type, seed_cond, config)
  result <- plot_potential_outcomes(dat, effect_type, seed_cond)

  # 保存
  fname <- sprintf("plot_potential_outcomes_n500_%s.png", effect_type)
  ggsave(file.path(out_dir, fname), result$layout,
         width = 14, height = 12, dpi = 150)
  cat(sprintf("Saved: %s\n", fname))

  # 統計量
  cat(sprintf("\nY(0) 未治療時:\n"))
  cat(sprintf("  平均: %.2f mmHg (理論値: %.0f)\n",
              mean(dat$Y0), config$mu0))
  cat(sprintf("  SD:   %.2f mmHg (理論値: %.1f)\n",
              sd(dat$Y0), config$sigma_y))
  cat(sprintf("  範囲: [%.1f, %.1f]\n", min(dat$Y0), max(dat$Y0)))

  cat(sprintf("\nY(1) 治療時:\n"))
  cat(sprintf("  平均: %.2f mmHg\n", mean(dat$Y1)))
  cat(sprintf("  SD:   %.2f mmHg\n", sd(dat$Y1)))
  cat(sprintf("  範囲: [%.1f, %.1f]\n", min(dat$Y1), max(dat$Y1)))

  cat(sprintf("\nITE = Y(1) - Y(0):\n"))
  cat(sprintf("  平均 (ATE): %.2f mmHg\n", mean(dat$ITE)))
  cat(sprintf("  SD:         %.2f mmHg\n", sd(dat$ITE)))
  cat(sprintf("  範囲: [%.2f, %.2f]\n", min(dat$ITE), max(dat$ITE)))

  cat(sprintf("\nτ = -ITE (降圧効果):\n"))
  cat(sprintf("  平均: %.2f mmHg (理論値: %.1f)\n",
              mean(dat$tau), config$tau_base))
  cat(sprintf("  SD:   %.2f mmHg\n", sd(dat$tau)))
  cat(sprintf("  範囲: [%.2f, %.2f] (設定: [%.1f, %.1f])\n",
              min(dat$tau), max(dat$tau), config$tau_min, config$tau_max))

  cat(sprintf("\n観測データ:\n"))
  cat(sprintf("  対照群 (T=0): n=%d, 平均Y=%.2f\n",
              sum(dat$T == 0), mean(dat$Y_obs[dat$T == 0])))
  cat(sprintf("  治療群 (T=1): n=%d, 平均Y=%.2f\n",
              sum(dat$T == 1), mean(dat$Y_obs[dat$T == 1])))
  cat(sprintf("  単純差分: %.2f mmHg\n",
              mean(dat$Y_obs[dat$T == 1]) - mean(dat$Y_obs[dat$T == 0])))
}

cat("\n✓ 潜在的アウトカムの可視化が完了しました\n")
