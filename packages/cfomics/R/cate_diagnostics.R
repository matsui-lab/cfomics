#' CATE Diagnostics: Heterogeneity Validation
#'
#' Tests for treatment effect heterogeneity using GATES (Sorted Group ATE)
#' and BLP (Best Linear Predictor) approaches from Chernozhukov et al. (2018).
#'
#' @param result A cfomics_result object
#' @param data Data frame used for fitting
#' @param formula Formula of the form Y ~ T | X1 + X2 + ...
#' @param n_groups Integer, number of ITE quantile groups for GATES (default 4)
#' @return A list of class "cf_cate_diagnostics"
#' @export
cf_cate_diagnostics <- function(result, data, formula, n_groups = 4L) {
  if (!inherits(result, "cfomics_result")) {
    rlang::abort("Input must be a cfomics_result object",
                 class = "cfomics_invalid_input")
  }

  parsed <- parse_cf_formula(formula, data)
  Y <- parsed$Y
  T_var <- parsed$T
  ite <- result$fit$res$ite
  n <- length(ite)

  # Check for constant ITE
  ite_var <- stats::var(ite)
  is_constant <- is.na(ite_var) || ite_var < 1e-10

  if (is_constant) {
    cli::cli_inform("ITE values are constant. GATES/BLP analysis skipped.")
    out <- list(
      constant_ite = TRUE,
      heterogeneity_detected = FALSE,
      ite_variance = 0,
      gates = NULL,
      blp = NULL,
      method = result$method
    )
    class(out) <- c("cf_cate_diagnostics", "list")
    return(out)
  }

  # --- GATES: Sorted Group ATE ---
  ite_quantiles <- stats::quantile(ite, probs = seq(0, 1, length.out = n_groups + 1))
  group <- cut(ite, breaks = ite_quantiles, labels = FALSE, include.lowest = TRUE)

  gates_results <- lapply(seq_len(n_groups), function(g) {
    idx <- group == g
    if (sum(idx & T_var == 1) < 2 || sum(idx & T_var == 0) < 2) {
      return(data.frame(group = g, gate = NA_real_, se = NA_real_,
                        ci_lower = NA_real_, ci_upper = NA_real_,
                        n = sum(idx)))
    }
    fit_g <- stats::lm(Y[idx] ~ T_var[idx])
    coefs <- summary(fit_g)$coefficients
    gate <- coefs[2, 1]
    se <- coefs[2, 2]
    data.frame(group = g, gate = gate, se = se,
               ci_lower = gate - 1.96 * se, ci_upper = gate + 1.96 * se,
               n = sum(idx))
  })
  gates_df <- do.call(rbind, gates_results)

  # --- BLP: Best Linear Predictor ---
  ite_centered <- ite - mean(ite)
  blp_fit <- stats::lm(Y ~ T_var + T_var:ite_centered)
  blp_summary <- summary(blp_fit)$coefficients

  blp_result <- list(
    ate_coef = blp_summary["T_var", "Estimate"],
    ate_se = blp_summary["T_var", "Std. Error"],
    ate_pvalue = blp_summary["T_var", "Pr(>|t|)"],
    het_coef = blp_summary["T_var:ite_centered", "Estimate"],
    het_se = blp_summary["T_var:ite_centered", "Std. Error"],
    het_pvalue = blp_summary["T_var:ite_centered", "Pr(>|t|)"]
  )

  # Heterogeneity detection
  het_detected <- blp_result$het_pvalue < 0.05

  out <- list(
    constant_ite = FALSE,
    heterogeneity_detected = het_detected,
    ite_variance = ite_var,
    gates = gates_df,
    blp = blp_result,
    method = result$method
  )

  class(out) <- c("cf_cate_diagnostics", "list")
  out
}

#' Print CATE diagnostics
#' @param x A cf_cate_diagnostics object
#' @param ... Additional arguments (ignored)
#' @export
print.cf_cate_diagnostics <- function(x, ...) {
  cat("\nCATE Diagnostics\n")
  cat("=================\n")
  cat(sprintf("Method: %s\n", x$method))

  if (x$constant_ite) {
    cat("\nITE is constant (no heterogeneity). GATES/BLP skipped.\n")
    return(invisible(x))
  }

  cat(sprintf("ITE variance: %.4f\n\n", x$ite_variance))

  cat("GATES (Sorted Group ATE):\n")
  if (!is.null(x$gates)) {
    print(x$gates, row.names = FALSE, digits = 3)
  }

  if (!is.null(x$blp)) {
    cat("\nBLP (Best Linear Predictor):\n")
    cat(sprintf("  ATE coefficient: %.3f (p = %.4f)\n",
                x$blp$ate_coef, x$blp$ate_pvalue))
    cat(sprintf("  Heterogeneity:   %.3f (p = %.4f)\n",
                x$blp$het_coef, x$blp$het_pvalue))
    if (x$heterogeneity_detected) {
      cat("  -> Significant heterogeneity detected\n")
    } else {
      cat("  -> No significant heterogeneity\n")
    }
  }

  invisible(x)
}

#' Plot CATE diagnostics
#' @param x A cf_cate_diagnostics object
#' @param ... Additional arguments (ignored)
#' @export
plot.cf_cate_diagnostics <- function(x, ...) {
  if (x$constant_ite || is.null(x$gates)) {
    message("No GATES data to plot (constant ITE)")
    return(invisible(NULL))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for CATE diagnostic plots")
    return(invisible(NULL))
  }

  gates <- x$gates
  gates$group <- factor(gates$group)

  p <- ggplot2::ggplot(gates, ggplot2::aes(x = group, y = gate)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                           width = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(x = "ITE Quartile Group", y = "Group ATE",
                  title = "GATES: Sorted Group Average Treatment Effects",
                  subtitle = sprintf("Heterogeneity %s (BLP p = %.4f)",
                                     ifelse(x$heterogeneity_detected, "detected", "not detected"),
                                     x$blp$het_pvalue)) +
    ggplot2::theme_minimal()
  print(p)
  invisible(p)
}
