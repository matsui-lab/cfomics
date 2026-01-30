#' Model Diagnostics for Causal Inference Results
#'
#' Evaluates the quality of propensity score and outcome models used in
#' causal inference estimation.
#'
#' @param result A cfomics_result object
#' @param data Data frame used for fitting
#' @param formula Formula of the form Y ~ T | X1 + X2 + ...
#' @return A list of class "cf_model_diagnostics"
#' @export
cf_model_diagnostics <- function(result, data, formula) {
  if (!inherits(result, "cfomics_result")) {
    rlang::abort("Input must be a cfomics_result object",
                 class = "cfomics_invalid_input")
  }

  # Parse formula
  parsed <- parse_cf_formula(formula, data)
  outcome_name <- parsed$outcome_name
  treatment_name <- parsed$treatment_name
  covariate_names <- parsed$covariate_names

  Y <- data[[outcome_name]]
  T_var <- data[[treatment_name]]
  n <- length(Y)

  # --- PS diagnostics ---
  ps_formula <- stats::as.formula(
    paste(treatment_name, "~", paste(covariate_names, collapse = " + "))
  )
  ps_model <- stats::glm(ps_formula, data = data, family = stats::binomial())
  ps_hat <- stats::predict(ps_model, type = "response")

  # C-statistic (concordance / AUC via Mann-Whitney)
  ps_t1 <- ps_hat[T_var == 1]
  ps_t0 <- ps_hat[T_var == 0]
  cstat <- mean(outer(ps_t1, ps_t0, ">")) + 0.5 * mean(outer(ps_t1, ps_t0, "=="))

  # Extreme weights
  weights <- ifelse(T_var == 1, 1 / ps_hat, 1 / (1 - ps_hat))
  extreme_threshold <- stats::quantile(weights, 0.99)
  extreme_pct <- mean(weights > extreme_threshold) * 100

  # --- Outcome model diagnostics ---
  y0_hat <- result$fit$res$y0_hat
  y1_hat <- result$fit$res$y1_hat
  y_hat <- ifelse(T_var == 1, y1_hat, y0_hat)
  residuals <- Y - y_hat

  resid_summary <- data.frame(
    group = c("Control (T=0)", "Treated (T=1)", "Overall"),
    mean_resid = c(
      mean(residuals[T_var == 0]),
      mean(residuals[T_var == 1]),
      mean(residuals)
    ),
    sd_resid = c(
      stats::sd(residuals[T_var == 0]),
      stats::sd(residuals[T_var == 1]),
      stats::sd(residuals)
    ),
    stringsAsFactors = FALSE
  )

  # Normality test
  resid_sample <- if (n > 5000) sample(residuals, 5000) else residuals
  normality_p <- tryCatch(
    stats::shapiro.test(resid_sample)$p.value,
    error = function(e) NA_real_
  )

  out <- list(
    ps_cstat = cstat,
    ps_hat = ps_hat,
    extreme_weight_pct = extreme_pct,
    weights = weights,
    residuals = residuals,
    residual_summary = resid_summary,
    normality_pvalue = normality_p,
    y_hat = y_hat,
    treatment = T_var,
    method = result$method
  )

  class(out) <- c("cf_model_diagnostics", "list")
  out
}

#' @export
print.cf_model_diagnostics <- function(x, ...) {
  cat("\nModel Diagnostics\n")
  cat("==================\n")
  cat(sprintf("Method: %s\n\n", x$method))

  cat("PS Model:\n")
  cat(sprintf("  C-statistic (AUC): %.3f", x$ps_cstat))
  if (x$ps_cstat > 0.95) cat(" [WARNING: near-deterministic]")
  else if (x$ps_cstat < 0.55) cat(" [WARNING: poor discrimination]")
  else cat(" [OK]")
  cat("\n")
  cat(sprintf("  Extreme weights (>99th pct): %.1f%%\n\n", x$extreme_weight_pct))

  cat("Outcome Model Residuals:\n")
  print(x$residual_summary, row.names = FALSE)

  if (!is.na(x$normality_pvalue)) {
    cat(sprintf("\n  Shapiro-Wilk p-value: %.4f", x$normality_pvalue))
    if (x$normality_pvalue < 0.05) cat(" [non-normal]")
    cat("\n")
  }

  invisible(x)
}

#' @export
plot.cf_model_diagnostics <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for diagnostic plots")
    return(invisible(NULL))
  }

  df_resid <- data.frame(
    fitted = x$y_hat,
    residuals = x$residuals,
    group = ifelse(x$treatment == 1, "Treated", "Control")
  )
  df_ps <- data.frame(
    ps = x$ps_hat,
    group = ifelse(x$treatment == 1, "Treated", "Control")
  )

  p1 <- ggplot2::ggplot(df_ps, ggplot2::aes(x = .data$ps, fill = .data$group)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(title = "PS Distribution", x = "Propensity Score", fill = "Group") +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(df_resid, ggplot2::aes(x = .data$fitted, y = .data$residuals, color = .data$group)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") +
    ggplot2::theme_minimal()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(p1 + p2)
  } else {
    print(p1)
    print(p2)
  }
  invisible(NULL)
}
