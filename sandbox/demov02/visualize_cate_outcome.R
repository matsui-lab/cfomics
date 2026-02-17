#!/usr/bin/env Rscript
# visualize_cate_outcome.R
# 真のCATE（tau）とアウトカム（Y）の分布を linear / threshold3 両方で可視化
# Run from repo root: Rscript sandbox/visualize_cate_outcome.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/benchmark_results/cate_demo_v02"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- DGP パラメータ（run_cate_demo_v02.R と同一） ----
n            <- 500
p            <- 50
K            <- 5
sigma_lambda <- 0.3
sigma_e      <- 0.3
theta_nb     <- 10
libsize_sd   <- 0.2
M_genes      <- 1:5
P_genes      <- 6:10
s_tau        <- 5
s_eta        <- 10
beta_A       <- 3
beta_P       <- 1
sigma_y      <- 5
mu0          <- 120

# ---- 共通DGP（潜在構造・観測まで、seed固定） ----
generate_base <- function(n, seed) {
  set.seed(seed)
  A      <- rnorm(n, 0, 1)
  G      <- matrix(rnorm(n * K), n, K)
  Lambda <- matrix(rnorm(p * K, 0, sigma_lambda), p, K)
  b0     <- rnorm(p, log(50), 0.2)
  bA     <- rnorm(p, 0.15, 0.05)
  eps_E  <- matrix(rnorm(n * p, 0, sigma_e), n, p)
  E      <- outer(rep(1, n), b0) + outer(A, bA) + G %*% t(Lambda) + eps_E

  s      <- exp(rnorm(n, 0, libsize_sd))
  mu_mat <- exp(E) * s
  C      <- matrix(rnbinom(n * p, size = theta_nb, mu = as.vector(mu_mat)), n, p)

  T_vec  <- rbinom(n, 1, 0.5)

  # 予後（effect_typeに依存しない）
  v   <- rowSums(E[, P_genes])
  vc  <- v - mean(v)
  eta_raw <- beta_A * A + beta_P * vc
  eta <- eta_raw * (s_eta / sd(eta_raw))
  eps_y <- rnorm(n, 0, sigma_y)

  # 修飾スコア（effect_typeに依存）
  u  <- rowSums(E[, M_genes])
  uc <- u - mean(u)

  list(A = A, E = E, u = u, uc = uc, eta = eta, eps_y = eps_y,
       T_vec = T_vec, mu0 = mu0)
}

make_tau <- function(base, effect_type) {
  if (effect_type == "linear") {
    tau_raw <- base$uc
  } else {
    q1 <- quantile(base$u, 0.33)
    q2 <- quantile(base$u, 0.66)
    tau_raw <- ifelse(base$u <= q1, -1, ifelse(base$u <= q2, 0, 1))
  }
  sd_tau <- sd(tau_raw)
  if (sd_tau < 1e-8) rep(0, length(tau_raw)) else tau_raw * (s_tau / sd_tau)
}

make_outcomes <- function(base, tau) {
  Y0 <- base$mu0 + base$eta + base$eps_y
  Y1 <- Y0 + tau
  Y  <- Y0 + base$T_vec * tau
  list(Y0 = Y0, Y1 = Y1, Y = Y)
}

# n=500 の seed（linear: 20270215, threshold3: 20270315）を使用
# 同一潜在構造で比較するため linear の seed を共通で使用
base <- generate_base(n = 500, seed = 20270215)

results <- lapply(c("linear", "threshold3"), function(et) {
  tau <- make_tau(base, et)
  out <- make_outcomes(base, tau)
  list(
    effect_type = et,
    tau   = tau,
    Y     = out$Y,
    Y0    = out$Y0,
    Y1    = out$Y1,
    T_vec = base$T_vec
  )
})

# ================================================================
# 図1: tau の分布（linear vs threshold3）
# ================================================================
df_tau <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    tau         = r$tau,
    effect_type = r$effect_type
  )
}))
df_tau$effect_type <- factor(df_tau$effect_type,
                              levels = c("linear", "threshold3"),
                              labels = c("Linear", "Threshold3"))

tau_stats <- aggregate(tau ~ effect_type, data = df_tau,
                       FUN = function(x) c(mean = mean(x), sd = sd(x),
                                           q05 = quantile(x, 0.05),
                                           q95 = quantile(x, 0.95)))

p_tau_lin <- ggplot(subset(df_tau, effect_type == "Linear"), aes(x = tau)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "#FC8D59", color = "white", alpha = 0.85) +
  geom_density(color = "#D73027", linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("mean=%.2f\nSD=%.2f",
                           mean(subset(df_tau, effect_type=="Linear")$tau),
                           sd(subset(df_tau, effect_type=="Linear")$tau)),
           size = 3.5, color = "#D73027") +
  labs(title = "Linear effect", x = "True CATE (τ)", y = "Density") +
  theme_bw(base_size = 12)

p_tau_thr <- ggplot(subset(df_tau, effect_type == "Threshold3"), aes(x = tau)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20,
                 fill = "#91BFDB", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("mean=%.2f\nSD=%.2f\n3 mass points",
                           mean(subset(df_tau, effect_type=="Threshold3")$tau),
                           sd(subset(df_tau, effect_type=="Threshold3")$tau)),
           size = 3.5, color = "#4575B4") +
  labs(title = "Threshold3 effect", x = "True CATE (τ)", y = "Density") +
  theme_bw(base_size = 12)

fig1 <- (p_tau_lin | p_tau_thr) +
  plot_annotation(
    title    = "True CATE (τ) distribution  [n=500, s_tau=5]",
    subtitle = "Same latent structure, different effect mechanism",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 11, color = "grey40"))
  )

ggsave(file.path(out_dir, "dgp_viz_07_tau_distribution.png"),
       fig1, width = 10, height = 4.5, dpi = 150)
cat("Saved: dgp_viz_07_tau_distribution.png\n")


# ================================================================
# 図2: アウトカム Y の分布（治療群別、linear vs threshold3）
# ================================================================
df_Y <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    Y           = r$Y,
    T           = factor(r$T_vec, labels = c("Control (T=0)", "Treated (T=1)")),
    effect_type = r$effect_type
  )
}))
df_Y$effect_type <- factor(df_Y$effect_type,
                            levels = c("linear", "threshold3"),
                            labels = c("Linear", "Threshold3"))

p_Y_lin <- ggplot(subset(df_Y, effect_type == "Linear"),
                  aes(x = Y, fill = T, color = T)) +
  geom_density(alpha = 0.35, linewidth = 0.9) +
  scale_fill_manual(values  = c("#74C476", "#FD8D3C")) +
  scale_color_manual(values = c("#238B45", "#D94701")) +
  labs(title = "Linear effect", x = "Y (blood pressure)", y = "Density",
       fill = NULL, color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.78, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.7)))

p_Y_thr <- ggplot(subset(df_Y, effect_type == "Threshold3"),
                  aes(x = Y, fill = T, color = T)) +
  geom_density(alpha = 0.35, linewidth = 0.9) +
  scale_fill_manual(values  = c("#74C476", "#FD8D3C")) +
  scale_color_manual(values = c("#238B45", "#D94701")) +
  labs(title = "Threshold3 effect", x = "Y (blood pressure)", y = "Density",
       fill = NULL, color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.78, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.7)))

fig2 <- (p_Y_lin | p_Y_thr) +
  plot_annotation(
    title    = "Observed outcome Y by treatment group  [n=500]",
    subtitle = sprintf("mu0=%.0f, sigma_y=%.0f, s_eta=%.0f", mu0, sigma_y, s_eta),
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 11, color = "grey40"))
  )

ggsave(file.path(out_dir, "dgp_viz_08_outcome_by_treatment.png"),
       fig2, width = 10, height = 4.5, dpi = 150)
cat("Saved: dgp_viz_08_outcome_by_treatment.png\n")


# ================================================================
# 図3: tau と Y の散布図（τ vs Y0、τとの関係確認）
# ================================================================
df_scatter <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    tau         = r$tau,
    Y0          = r$Y0,
    T           = factor(r$T_vec, labels = c("Control", "Treated")),
    effect_type = r$effect_type
  )
}))
df_scatter$effect_type <- factor(df_scatter$effect_type,
                                  levels = c("linear", "threshold3"),
                                  labels = c("Linear", "Threshold3"))

fig3 <- ggplot(df_scatter, aes(x = tau, y = Y0, color = T)) +
  geom_point(alpha = 0.3, size = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("#238B45", "#D94701")) +
  facet_wrap(~effect_type, scales = "free_x") +
  labs(
    title    = "True CATE (τ) vs. baseline outcome Y(0)  [n=500]",
    subtitle = "Dashed line: OLS. Effect modifier and prognostic factors are independent by design.",
    x = "True CATE (τ)", y = "Y(0) — baseline outcome",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "dgp_viz_09_tau_vs_Y0.png"),
       fig3, width = 10, height = 5, dpi = 150)
cat("Saved: dgp_viz_09_tau_vs_Y0.png\n")


# ================================================================
# 図4: 4パネル統合図（tau + Y をまとめて比較）
# ================================================================
fig4 <- (p_tau_lin | p_tau_thr) / (p_Y_lin | p_Y_thr) +
  plot_annotation(
    title    = "True CATE and Outcome distributions — Linear vs Threshold3  [n=500]",
    tag_levels = "A",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave(file.path(out_dir, "dgp_viz_10_combined_tau_outcome.png"),
       fig4, width = 11, height = 9, dpi = 150)
cat("Saved: dgp_viz_10_combined_tau_outcome.png\n")

cat("\n全図の保存完了:", out_dir, "\n")
