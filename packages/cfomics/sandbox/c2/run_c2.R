#!/usr/bin/env Rscript
# =============================================================================
# Scenario C2: Blood Pressure CATE Benchmark
#
# Compare CATE estimation accuracy across 3 methods (grf, drlearner, ganite).
# True CATE differs by gender: male = -10, female = -5.
# =============================================================================

devtools::load_all("packages/cfomics")

# -- 0. Python Environment Setup ----------------------------------------------
# Use the unified conda env (econml + tensorflow)
reticulate::use_condaenv("cfomics_unified", required = TRUE)

# -- 1. Data Generation -------------------------------------------------------

set.seed(42)
n <- 500
gender <- rbinom(n, 1, 0.5)
age <- rnorm(n, mean = 50, sd = 10)

# Propensity: depends on gender and age
propensity <- plogis(-0.5 + 0.5 * gender + 0.02 * (age - 50))
drug_a <- rbinom(n, 1, propensity)

# True CATE: male = -10, female = -5
true_cate <- ifelse(gender == 1, -10, -5)

# Blood pressure outcome
baseline_bp <- 130 + 5 * gender + 0.3 * (age - 50) + rnorm(n, sd = 8)
blood_pressure <- baseline_bp + true_cate * drug_a

dat <- data.frame(blood_pressure, drug_a, gender, age)

# Ground truth
true_ate <- mean(true_cate)
cat(sprintf("Ground truth: ATE = %.2f, CATE_male = -10, CATE_female = -5\n\n", true_ate))

# -- 2. Fit 3 Methods ---------------------------------------------------------

fml <- blood_pressure ~ drug_a | gender + age
methods <- c("grf", "drlearner", "ganite")

results <- list()

for (m in methods) {
  cat(sprintf("--- Fitting: %s ---\n", m))
  fit <- tryCatch(
    cf_fit(fml, data = dat, method = m),
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

  cate_male   <- mean(ite_hat[gender == 1])
  cate_female <- mean(ite_hat[gender == 0])
  pehe <- sqrt(mean((ite_hat - true_cate)^2))

  results[[m]] <- list(
    status      = "ok",
    ate         = ate_hat,
    cate_male   = cate_male,
    cate_female = cate_female,
    pehe        = pehe
  )
  cat(sprintf("  ATE = %.2f, CATE_male = %.2f, CATE_female = %.2f, PEHE = %.2f\n\n",
              ate_hat, cate_male, cate_female, pehe))
}

# -- 3. Summary Table ---------------------------------------------------------

cat("============================================================\n")
cat("Scenario C2 Results\n")
cat("============================================================\n")
cat(sprintf("%-10s | %7s | %10s | %12s | %6s | %s\n",
            "method", "ATE", "CATE_male", "CATE_female", "PEHE", "status"))
cat(paste0(rep("-", 66), collapse = ""), "\n")
cat(sprintf("%-10s | %7.2f | %10.2f | %12.2f | %6s | %s\n",
            "truth", true_ate, -10, -5, "---", "---"))

for (m in methods) {
  r <- results[[m]]
  if (r$status == "skip") {
    cat(sprintf("%-10s | %7s | %10s | %12s | %6s | %s\n",
                m, "---", "---", "---", "---", "skip"))
  } else {
    cat(sprintf("%-10s | %7.2f | %10.2f | %12.2f | %6.2f | %s\n",
                m, r$ate, r$cate_male, r$cate_female, r$pehe, r$status))
  }
}
cat("============================================================\n")
