#!/usr/bin/env Rscript
# visualize_dgp.R
# DGPで生成したオミクス・年齢データの分布可視化
# Run from repo root: Rscript sandbox/visualize_dgp.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "sandbox/benchmark_results/cate_demo_v02"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- DGP (run_cate_demo_v02.R と同一パラメータ) ----
set.seed(20270215)   # n=500, linear の条件と同一seed
n         <- 500
p         <- 50
K         <- 5
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

# 潜在変数
A      <- rnorm(n, 0, 1)
G      <- matrix(rnorm(n * K), n, K)
Lambda <- matrix(rnorm(p * K, 0, sigma_lambda), p, K)
b0     <- rnorm(p, log(50), 0.2)
bA     <- rnorm(p, 0.15, 0.05)
eps_E  <- matrix(rnorm(n * p, 0, sigma_e), n, p)

# 潜在ログ発現
E <- outer(rep(1,n), b0) + outer(A, bA) + G %*% t(Lambda) + eps_E

# RNA-seq観測
s      <- exp(rnorm(n, 0, libsize_sd))
mu_mat <- exp(E) * s
C      <- matrix(rnbinom(n * p, size = theta_nb, mu = as.vector(mu_mat)), n, p)
L      <- log1p(C)
colnames(C) <- colnames(L) <- paste0("gene", 1:p)

T_vec <- rbinom(n, 1, 0.5)

# CATE (linear)
u      <- rowSums(E[, M_genes])
uc     <- u - mean(u)
tau    <- uc * (s_tau / sd(uc))

# 予後
v      <- rowSums(E[, P_genes])
vc     <- v - mean(v)
eta    <- (beta_A * A + beta_P * vc)
eta    <- eta * (s_eta / sd(eta))
Y0     <- mu0 + eta + rnorm(n, 0, sigma_y)
Y      <- Y0 + T_vec * tau

# z-score（trainのみで計算→全体に適用、デモ用に全体で計算）
train_idx <- sort(sample(n, round(n * 0.8)))
mu_tr  <- colMeans(L[train_idx, ])
sd_tr  <- pmax(apply(L[train_idx, ], 2, sd), 1e-6)
L_z    <- sweep(sweep(L, 2, mu_tr, "-"), 2, sd_tr, "/")

gene_role <- ifelse(1:p %in% M_genes, "predictive",
             ifelse(1:p %in% P_genes, "prognostic", "noise"))


# ================================================================
# 図1: 年齢・治療の分布
# ================================================================
df_base <- data.frame(A = A, T = factor(T_vec, labels = c("Control", "Treated")))

p_age <- ggplot(df_base, aes(x = A)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "#4292C6", color = "white", alpha = 0.8) +
  geom_density(color = "#08306B", linewidth = 0.8) +
  labs(title = "Age (A)", x = "A", y = "Density") +
  theme_bw(base_size = 12)

p_treat <- ggplot(df_base, aes(x = T, fill = T)) +
  geom_bar(width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#74C476", "#FD8D3C")) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "Treatment (T)", x = "", y = "Count") +
  theme_bw(base_size = 12)

fig1 <- p_age + p_treat +
  plot_annotation(title = "Covariates: Age & Treatment  [n=500, seed=20270215]",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(out_dir, "dgp_viz_01_age_treatment.png"),
       fig1, width = 8, height = 4, dpi = 150)
cat("Saved: dgp_viz_01_age_treatment.png\n")


# ================================================================
# 図2: カウント分布（生・log1p）遺伝子タイプ別
# ================================================================
# 代表遺伝子を各タイプから3つずつサンプリング
set.seed(1)
genes_show <- c(
  sample(M_genes, 3),
  sample(P_genes, 3),
  sample(which(!1:p %in% c(M_genes, P_genes)), 3)
)
role_labels <- c(rep("predictive", 3), rep("prognostic", 3), rep("noise", 3))

df_count <- do.call(rbind, lapply(seq_along(genes_show), function(i) {
  g <- genes_show[i]
  data.frame(
    count = C[, g],
    log1p = L[, g],
    gene  = paste0("gene", g),
    role  = role_labels[i]
  )
}))
df_count$role <- factor(df_count$role, levels = c("predictive", "prognostic", "noise"))

role_colors <- c(predictive = "#E41A1C", prognostic = "#377EB8", noise = "#999999")

p_raw <- ggplot(df_count, aes(x = count, fill = role, color = role)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  facet_wrap(~gene, scales = "free", ncol = 3) +
  scale_fill_manual(values = role_colors) +
  scale_color_manual(values = role_colors) +
  labs(title = "Raw counts (C)", x = "Count", y = "Frequency") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))

p_log <- ggplot(df_count, aes(x = log1p, fill = role, color = role)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  facet_wrap(~gene, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = role_colors) +
  scale_color_manual(values = role_colors) +
  labs(title = "Log-transformed counts (L = log1p(C))", x = "log1p(count)", y = "Frequency") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))

fig2 <- p_raw / p_log +
  plot_annotation(
    title = "Gene expression distributions by role  [9 representative genes]",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(out_dir, "dgp_viz_02_count_distributions.png"),
       fig2, width = 10, height = 10, dpi = 150)
cat("Saved: dgp_viz_02_count_distributions.png\n")


# ================================================================
# 図3: 全50遺伝子のlog1p分布まとめ（箱ひげ図）
# ================================================================
df_all <- data.frame(
  log1p = as.vector(L),
  gene  = rep(paste0("gene", 1:p), each = n),
  role  = rep(ifelse(1:p %in% M_genes, "predictive",
              ifelse(1:p %in% P_genes, "prognostic", "noise")), each = n),
  gene_num = rep(1:p, each = n)
)
df_all$role <- factor(df_all$role, levels = c("predictive", "prognostic", "noise"))
df_all$gene <- factor(df_all$gene, levels = paste0("gene", 1:p))

fig3 <- ggplot(df_all, aes(x = gene, y = log1p, fill = role, color = role)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.4, alpha = 0.6) +
  scale_fill_manual(values = role_colors) +
  scale_color_manual(values = role_colors) +
  labs(
    title = "log1p(count) distribution — all 50 genes  [n=500]",
    subtitle = "red=predictive (j=1-5), blue=prognostic (j=6-10), grey=noise (j=11-50)",
    x = "Gene", y = "log1p(count)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.position = "top"
  )

ggsave(file.path(out_dir, "dgp_viz_03_all_genes_boxplot.png"),
       fig3, width = 14, height = 5, dpi = 150)
cat("Saved: dgp_viz_03_all_genes_boxplot.png\n")


# ================================================================
# 図4: z-score後の分布（前処理確認）
# ================================================================
df_z <- data.frame(
  z     = as.vector(L_z[, c(M_genes, P_genes, sample(11:50, 5))]),
  role  = rep(c(rep("predictive", 5), rep("prognostic", 5), rep("noise", 5)), each = n),
  gene  = rep(c(paste0("gene", M_genes),
                paste0("gene", P_genes),
                paste0("gene", sample(11:50, 5))), each = n)
)
df_z$role <- factor(df_z$role, levels = c("predictive", "prognostic", "noise"))

fig4 <- ggplot(df_z, aes(x = z, fill = role, color = role)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35,
                 alpha = 0.5, position = "identity") +
  geom_density(linewidth = 0.7, alpha = 0) +
  facet_wrap(~gene, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = role_colors) +
  scale_color_manual(values = role_colors) +
  labs(
    title = "z-scored log1p — model input X (15 genes)",
    x = "z-score", y = "Density"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))

ggsave(file.path(out_dir, "dgp_viz_04_zscore_input.png"),
       fig4, width = 14, height = 7, dpi = 150)
cat("Saved: dgp_viz_04_zscore_input.png\n")


# ================================================================
# 図5: 深度因子・アウトカム・tau の分布
# ================================================================
df_outcome <- data.frame(
  s    = s,
  tau  = tau,
  Y    = Y,
  T    = factor(T_vec, labels = c("Control", "Treated"))
)

p_s <- ggplot(df_outcome, aes(x = s)) +
  geom_histogram(bins = 30, fill = "#9E9AC8", color = "white", alpha = 0.8) +
  labs(title = "Library size factor (s)", x = "s", y = "Count") +
  theme_bw(base_size = 11)

p_tau <- ggplot(df_outcome, aes(x = tau)) +
  geom_histogram(bins = 30, fill = "#FB6A4A", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = sprintf("True CATE (τ)  SD=%.2f", sd(tau)), x = "τ", y = "Count") +
  theme_bw(base_size = 11)

p_Y <- ggplot(df_outcome, aes(x = Y, fill = T, color = T)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values = c("#74C476", "#FD8D3C")) +
  scale_color_manual(values = c("#238B45", "#D94701")) +
  labs(title = "Observed outcome Y by treatment", x = "Y (blood pressure)", y = "Density") +
  theme_bw(base_size = 11) +
  theme(legend.position = "inside", legend.position.inside = c(0.15, 0.8))

fig5 <- (p_s | p_tau | p_Y) +
  plot_annotation(
    title = "DGP outputs: depth factor, true CATE, outcome  [n=500, linear]",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(out_dir, "dgp_viz_05_outcome_tau.png"),
       fig5, width = 13, height = 4, dpi = 150)
cat("Saved: dgp_viz_05_outcome_tau.png\n")


# ================================================================
# 図6: 遺伝子間相関ヒートマップ（ランク相関）
# ================================================================
genes_heatmap <- c(M_genes, P_genes, sample(11:50, 10, replace = FALSE))
cor_mat <- cor(L[, genes_heatmap], method = "spearman")
rownames(cor_mat) <- colnames(cor_mat) <- paste0(
  "gene", genes_heatmap,
  ifelse(genes_heatmap %in% M_genes, "(P)", ifelse(genes_heatmap %in% P_genes, "(G)", ""))
)

df_cor <- as.data.frame(as.table(cor_mat))
colnames(df_cor) <- c("Gene1", "Gene2", "Correlation")

role_order <- c(
  paste0("gene", M_genes, "(P)"),
  paste0("gene", P_genes, "(G)"),
  paste0("gene", sort(genes_heatmap[!genes_heatmap %in% c(M_genes, P_genes)]))
)
df_cor$Gene1 <- factor(df_cor$Gene1, levels = role_order)
df_cor$Gene2 <- factor(df_cor$Gene2, levels = role_order)

fig6 <- ggplot(df_cor, aes(x = Gene1, y = Gene2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(
    title = "Gene–gene Spearman correlation (log1p counts)",
    subtitle = "(P)=predictive, (G)=prognostic, others=noise",
    x = "", y = ""
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

ggsave(file.path(out_dir, "dgp_viz_06_gene_correlation.png"),
       fig6, width = 9, height = 8, dpi = 150)
cat("Saved: dgp_viz_06_gene_correlation.png\n")

cat("\n全図の保存完了: ", out_dir, "\n")
