#!/usr/bin/env Rscript
# =============================================================================
# Scenario C1: Treatment Selection Bias — CATE Estimation Accuracy
#
# Compare CATE estimation accuracy across 5 methods under confounding by
# severity.  True CATE = -10 (constant for all subjects).
# =============================================================================

devtools::load_all("packages/cfomics")

# -- 0. Python Environment Setup ----------------------------------------------
reticulate::use_condaenv("cfomics_unified", required = TRUE)

# -- 1. Data Generation -------------------------------------------------------

set.seed(42)
n <- 100

severity  <- sample(1:10, n, replace = TRUE)
age       <- rnorm(n, 50, 10)

# Propensity depends on severity: sev=1 ~17%, sev=5 ~50%, sev=10 ~88%
treatment <- rbinom(n, 1, plogis(-2 + 0.4 * severity))

# True CATE = -10 (constant)
true_cate <- rep(-10, n)
true_ate  <- -10

blood_pressure <- 130 + 3 * severity + 0.3 * (age - 50) + rnorm(n, sd = 8) +
  (-10) * treatment

dat <- data.frame(blood_pressure, treatment, severity, age)

cat(sprintf("Ground truth: ATE = %.2f, CATE = -10 (constant)\n", true_ate))
cat(sprintf("Sample: n = %d, treated = %d (%.0f%%)\n\n",
            n, sum(treatment), 100 * mean(treatment)))

# -- 2. Treatment Distribution Visualisation ----------------------------------

out_dir <- "packages/cfomics/sandbox/c1"
png(file.path(out_dir, "treatment_distribution.png"),
    width = 900, height = 400, res = 120)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Left: treatment rate by severity
trt_rate <- tapply(treatment, severity, mean)
barplot(trt_rate, names.arg = names(trt_rate),
        xlab = "Severity", ylab = "Treatment Rate",
        main = "Treatment Rate by Severity",
        col = "steelblue", ylim = c(0, 1))
abline(h = 0.5, lty = 2, col = "grey50")

# Right: stacked counts (treated / untreated)
cnt_trt   <- tapply(treatment, severity, sum)
cnt_untrt <- tapply(1 - treatment, severity, sum)
counts    <- rbind(Untreated = cnt_untrt, Treated = cnt_trt)
barplot(counts, names.arg = colnames(counts),
        xlab = "Severity", ylab = "Count",
        main = "Sample Distribution by Severity",
        col = c("grey70", "steelblue"), legend.text = TRUE,
        args.legend = list(x = "topleft", bty = "n"))

dev.off()
cat("Saved: treatment_distribution.png\n\n")

# -- 3. Fit 5 Methods ---------------------------------------------------------

fml     <- blood_pressure ~ treatment | severity + age
methods <- c("ipw", "grf", "drlearner", "dowhy_gcm", "cavae")

# DAG for dowhy_gcm: Python backend uses "T" for treatment and "Y" for outcome
dag <- igraph::graph_from_edgelist(rbind(
  c("severity", "T"),
  c("T", "Y"),
  c("severity", "Y"),
  c("age", "Y")
))

results <- list()

for (m in methods) {
  cat(sprintf("--- Fitting: %s ---\n", m))
  extra_args <- if (m == "dowhy_gcm") list(graph = dag) else list()
  if (m == "cavae") extra_args$outcome_dist <- "normal"
  fit <- tryCatch(
    do.call(cf_fit, c(list(formula = fml, data = dat, method = m), extra_args)),
    error = function(e) {
      cat(sprintf("  SKIP (%s)\n\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fit)) {
    results[[m]] <- list(status = "skip")
    next
  }

  ite_hat <- predict(fit, type = "ite")
  ate_hat <- predict(fit, type = "ate")

  bias <- ate_hat - true_ate
  pehe <- sqrt(mean((ite_hat - true_cate)^2))

  # Subgroup CATE: low severity (<=5) vs high severity (>5)
  cate_low  <- mean(ite_hat[severity <= 5])
  cate_high <- mean(ite_hat[severity > 5])

  results[[m]] <- list(
    status    = "ok",
    ate       = ate_hat,
    bias      = bias,
    pehe      = pehe,
    cate_low  = cate_low,
    cate_high = cate_high
  )
  cat(sprintf("  ATE = %.2f, Bias = %.2f, PEHE = %.2f\n",
              ate_hat, bias, pehe))
  cat(sprintf("  CATE_low = %.2f, CATE_high = %.2f\n\n",
              cate_low, cate_high))
}

# -- 4. Summary Table ---------------------------------------------------------

cat("============================================================\n")
cat("Scenario C1 Results — Treatment Selection Bias\n")
cat("============================================================\n")
cat(sprintf("%-10s | %7s | %7s | %6s | %10s | %11s | %s\n",
            "method", "ATE", "Bias", "PEHE", "CATE_low", "CATE_high", "status"))
cat(paste0(rep("-", 72), collapse = ""), "\n")
cat(sprintf("%-10s | %7.2f | %7.2f | %6s | %10.2f | %11.2f | %s\n",
            "truth", true_ate, 0, "---", -10, -10, "---"))

for (m in methods) {
  r <- results[[m]]
  if (r$status == "skip") {
    cat(sprintf("%-10s | %7s | %7s | %6s | %10s | %11s | %s\n",
                m, "---", "---", "---", "---", "---", "skip"))
  } else {
    cat(sprintf("%-10s | %7.2f | %7.2f | %6.2f | %10.2f | %11.2f | %s\n",
                m, r$ate, r$bias, r$pehe, r$cate_low, r$cate_high, r$status))
  }
}
cat("============================================================\n")
