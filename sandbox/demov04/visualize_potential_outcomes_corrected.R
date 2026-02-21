#!/usr/bin/env Rscript
# visualize_potential_outcomes_corrected.R
# 修正版DGPでのY(0), Y(1)可視化

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/demov04/results_corrected"

# ============================================================
# CONFIG
# ============================================================
config <- list(
  seed_base  = 20260215L,
  n_val      = 500L,
  effect_type = "linear",
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0, mu0 = 120.0, sigma_y = 5.0
)

# ============================================================
# DGP (修正版)
# ============================================================
generate_full_data_corrected <- function(n, effect_type, seed, cfg) {
  set.seed(seed)

  S  <- rbinom(n, 1, 0.5)
  Sc <- S - 0.5

  E_true  <- rnorm(n, cfg$m_E, cfg$s_E)
  mu_nb   <- exp(E_true)
  C_count <- rnbinom(n, size = cfg$theta_nb, mu = mu_nb)
  G_raw   <- log1p(C_count)
  delta   <- rnorm(n, 0, cfg$sigma_G)
  G_obs   <- G_raw + delta

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

  # 潜在的アウトカム (修正版)
  Y0 <- cfg$mu0 + rnorm(n, 0, cfg$sigma_y)
  Y1 <- Y0 + tau  # tau < 0 なので Y1 < Y0（血圧が下がる）

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
    ITE = Y1 - Y0  # = tau
  )
}

# ============================================================
# データ生成
# ============================================================
seed_cond <- config$seed_base + 2 * 100 + 10  # n=500, linear

dat <- generate_full_data_corrected(config$n_val, config$effect_type, seed_cond, config)

dat$Sex <- factor(ifelse(dat$S == 1, "Male", "Female"),
                  levels = c("Female", "Male"))

# ============================================================
# プロット作成
# ============================================================

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
  annotate("text", x = mean(dat$Y0), y = 0.07,
           label = sprintf("%.1f", mean(dat$Y0)),
           color = "#4DAF4A", size = 4, vjust = -0.5) +
  annotate("text", x = mean(dat$Y1), y = 0.07,
           label = sprintf("%.1f", mean(dat$Y1)),
           color = "#E41A1C", size = 4, vjust = -0.5) +
  labs(
    title = "潜在的アウトカムの分布（修正版）",
    subtitle = sprintf("平均: Y(0)=%.1f, Y(1)=%.1f, ATE=%.1f mmHg（負=降圧効果）",
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

# --- プロット3: ITE の分布 ---
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
    subtitle = sprintf("ITE = Y(1) - Y(0) = τ (負=降圧効果) | SD = %.2f", sd(dat$ITE)),
    x = "ITE (mmHg)",
    y = "密度"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# --- プロット4: τ の分布 ---
p4 <- ggplot(dat, aes(x = tau)) +
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
    subtitle = "τ = Y(1) - Y(0) = ITE  (負=降圧効果)",
    x = "τ (mmHg)",
    y = "密度"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# --- レイアウト ---
layout <- (p1 | p2) / (p3 | p4)
layout <- layout +
  plot_annotation(
    title = sprintf("潜在的アウトカムの可視化（修正版DGP, n=%d, %s)",
                   nrow(dat), config$effect_type),
    subtitle = "Y = Y0 + T×τ で、τ < 0 なので治療により血圧が下がる",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )
  )

fname <- "plot_potential_outcomes_corrected.png"
ggsave(file.path(out_dir, fname), layout, width = 14, height = 10, dpi = 150)
cat(sprintf("Saved: %s\n", fname))

# ============================================================
# 統計量
# ============================================================
cat(sprintf("\n=== 潜在的アウトカム統計（修正版, n=%d） ===\n", config$n_val))

cat(sprintf("\nY(0) 未治療時:\n"))
cat(sprintf("  平均: %.2f mmHg (理論値: %.0f)\n", mean(dat$Y0), config$mu0))
cat(sprintf("  SD:   %.2f mmHg (理論値: %.1f)\n", sd(dat$Y0), config$sigma_y))
cat(sprintf("  範囲: [%.1f, %.1f]\n", min(dat$Y0), max(dat$Y0)))

cat(sprintf("\nY(1) 治療時:\n"))
cat(sprintf("  平均: %.2f mmHg\n", mean(dat$Y1)))
cat(sprintf("  SD:   %.2f mmHg\n", sd(dat$Y1)))
cat(sprintf("  範囲: [%.1f, %.1f]\n", min(dat$Y1), max(dat$Y1)))

cat(sprintf("\nITE = τ = Y(1) - Y(0) (降圧効果):\n"))
cat(sprintf("  平均 (ATE): %.2f mmHg\n", mean(dat$ITE)))
cat(sprintf("  SD:         %.2f mmHg\n", sd(dat$ITE)))
cat(sprintf("  範囲: [%.2f, %.2f]\n", min(dat$ITE), max(dat$ITE)))

cat(sprintf("\n観測データ:\n"))
cat(sprintf("  対照群 (T=0): n=%d, 平均Y=%.2f mmHg\n",
            sum(dat$T == 0), mean(dat$Y_obs[dat$T == 0])))
cat(sprintf("  治療群 (T=1): n=%d, 平均Y=%.2f mmHg\n",
            sum(dat$T == 1), mean(dat$Y_obs[dat$T == 1])))
cat(sprintf("  単純差分: %.2f mmHg (負=降圧効果)\n",
            mean(dat$Y_obs[dat$T == 1]) - mean(dat$Y_obs[dat$T == 0])))

cat("\n✓ 修正版DGPでは、治療により血圧が約12 mmHg下がります\n")
