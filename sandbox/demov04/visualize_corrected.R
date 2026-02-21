#!/usr/bin/env Rscript
# visualize_corrected.R
# 修正版DGPでのgeneA vs CATE可視化

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
  train_frac = 0.8,
  m_E = log(50), s_E = 0.35, theta_nb = 30, sigma_G = 0.08,
  tau_base = -12.0, a_S = 2.0, a_G = 1.5,
  tau_min = -18.0, tau_max = -6.0, mu0 = 120.0, sigma_y = 5.0
)

# ============================================================
# DGP (修正版)
# ============================================================
generate_data_corrected <- function(n, effect_type, seed, cfg) {
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

  list(S = S, G_obs = G_obs, T_vec = T_vec, tau = tau)
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
# データ読み込み
# ============================================================
df_pred <- read.csv(file.path(out_dir, "predictions_corrected.csv"),
                    stringsAsFactors = FALSE)

seed_cond <- unique(df_pred$seed_condition)
n_val <- unique(df_pred$n)
effect_type <- unique(df_pred$effect_type)

dat <- generate_data_corrected(n_val, effect_type, seed_cond, config)
spl <- stratified_split(dat$T_vec, config$train_frac, seed_cond)
te <- spl$test

method_names <- c(causal_forest = "Causal Forest", drlearner = "DR-Learner (FIXED)")

for (meth in c("causal_forest", "drlearner")) {
  method_label <- method_names[meth]

  sub_so <- subset(df_pred,
                   covset == "SexOnly" & method == meth)
  sub_sg <- subset(df_pred,
                   covset == "SexGene" & method == meth)

  df_m <- merge(
    sub_so[, c("id", "tau_true", "tau_hat")],
    sub_sg[, c("id", "tau_hat")],
    by = "id", suffixes = c("_so", "_sg")
  )
  names(df_m)[names(df_m) == "tau_hat_so"] <- "SexOnly"
  names(df_m)[names(df_m) == "tau_hat_sg"] <- "SexGene"

  df_m$geneA <- dat$G_obs[df_m$id]
  df_m$Sex   <- factor(ifelse(dat$S[df_m$id] == 1, "Male", "Female"),
                       levels = c("Female", "Male"))

  # 真値（test のみ）
  df_true <- data.frame(
    geneA = dat$G_obs[te],
    Sex   = factor(ifelse(dat$S[te] == 1, "Male", "Female"),
                   levels = c("Female", "Male")),
    CATE  = dat$tau[te],
    Type  = "True τ",
    stringsAsFactors = FALSE
  )

  df_so <- data.frame(geneA = df_m$geneA, Sex = df_m$Sex,
                      CATE = df_m$SexOnly, Type = "SexOnly",
                      stringsAsFactors = FALSE)
  df_sg <- data.frame(geneA = df_m$geneA, Sex = df_m$Sex,
                      CATE = df_m$SexGene, Type = "SexGene",
                      stringsAsFactors = FALSE)

  df_plot_pred <- rbind(df_so, df_sg)
  df_plot_pred <- subset(df_plot_pred, !is.na(CATE))
  df_plot_pred$Type <- factor(df_plot_pred$Type,
                              levels = c("SexOnly", "SexGene"))

  df_true$Type <- factor(df_true$Type, levels = "True τ")

  type_colors <- c("True τ" = "#333333",
                   "SexOnly" = "#E41A1C",
                   "SexGene" = "#377EB8")
  type_shapes <- c("True τ" = 16, "SexOnly" = 17, "SexGene" = 15)
  type_sizes  <- c("True τ" = 1.2, "SexOnly" = 1.8, "SexGene" = 1.8)
  type_alphas <- c("True τ" = 0.35, "SexOnly" = 0.6, "SexGene" = 0.6)

  y_range <- range(dat$tau[te], na.rm = TRUE)
  y_center <- mean(y_range)
  y_half_width <- max(diff(y_range) * 0.8, 2.5)
  y_lim <- c(y_center - y_half_width, y_center + y_half_width)

  p <- ggplot() +
    geom_point(data = df_true,
               aes(x = geneA, y = CATE, color = Type, shape = Type),
               alpha = type_alphas["True τ"], size = type_sizes["True τ"]) +
    geom_point(data = df_plot_pred,
               aes(x = geneA, y = CATE, color = Type, shape = Type,
                   alpha = Type, size = Type)) +
    scale_alpha_manual(values = type_alphas, guide = "none") +
    scale_size_manual(values = type_sizes, guide = "none") +
    scale_color_manual(values = type_colors, name = NULL) +
    scale_shape_manual(values = type_shapes, name = NULL) +
    coord_cartesian(ylim = y_lim) +
    facet_wrap(~Sex) +
    labs(
      title = sprintf("%s — %s, n=%d  [CORRECTED DGP]",
                     method_label, effect_type, n_val),
      subtitle = sprintf("Black: true τ (負=降圧効果) | Red: SexOnly | Blue: SexGene | Test: %d",
                        length(te)),
      x = "geneA (log1p count + noise)",
      y = expression(CATE ~ (tau ~ "[mmHg]"))
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "bottom",
      strip.text    = element_text(size = 11, face = "bold")
    )

  fname <- sprintf("plot_geneA_vs_cate_%s_corrected.png", meth)
  ggsave(file.path(out_dir, fname), p, width = 10, height = 5.5, dpi = 150)
  cat(sprintf("Saved: %s\n", fname))
}

cat("\n✓ geneA vs CATE プロットを保存しました（修正版DGP）\n")
