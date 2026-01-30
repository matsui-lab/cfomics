#' Unified Diagnostics Report
#'
#' Runs all available diagnostic checks and produces a comprehensive report.
#'
#' @param result A cfomics_result object
#' @param data Data frame used for fitting
#' @param formula Formula of the form Y ~ T | X1 + X2 + ...
#' @param include Character vector of diagnostics to include.
#'   Options: "sensitivity", "model", "cate". Default: all three.
#' @return A list of class "cf_diagnostics_report"
#' @export
cf_diagnostics_report <- function(result, data, formula,
                                   include = c("sensitivity", "model", "cate")) {
  if (!inherits(result, "cfomics_result")) {
    rlang::abort("Input must be a cfomics_result object",
                 class = "cfomics_invalid_input")
  }

  include <- match.arg(include, c("sensitivity", "model", "cate"),
                       several.ok = TRUE)

  parsed <- parse_cf_formula(formula, data)

  out <- list(
    sensitivity = NULL,
    model = NULL,
    cate = NULL,
    balance = NULL,
    summary = NULL,
    method = result$method,
    n = result$meta$n,
    p = result$meta$p
  )

  # Sensitivity
  if ("sensitivity" %in% include) {
    out$sensitivity <- cf_sensitivity(result)
  }

  # Model diagnostics
  if ("model" %in% include) {
    out$model <- cf_model_diagnostics(result, data, formula)
  }

  # CATE diagnostics
  if ("cate" %in% include) {
    out$cate <- cf_cate_diagnostics(result, data, formula)
  }

  # Balance check (always run)
  out$balance <- cf_balance_check(data, parsed$treatment_name, parsed$covariate_names)

  # Build summary
  summary_rows <- list()

  if (!is.null(out$sensitivity)) {
    ev <- out$sensitivity$evalue_point
    status <- if (ev > 2) "Pass" else if (ev > 1.5) "Warning" else "Fail"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "Sensitivity (E-value)",
      value = sprintf("%.2f", ev),
      status = status, stringsAsFactors = FALSE
    )))
  }

  if (!is.null(out$model)) {
    cstat <- out$model$ps_cstat
    status <- if (cstat >= 0.6 && cstat <= 0.85) "Pass"
              else if (cstat < 0.5 || cstat > 0.95) "Fail"
              else "Warning"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "PS C-statistic",
      value = sprintf("%.3f", cstat),
      status = status, stringsAsFactors = FALSE
    )))

    ewp <- out$model$extreme_weight_pct
    status <- if (ewp < 1) "Pass" else if (ewp < 5) "Warning" else "Fail"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "Extreme weights",
      value = sprintf("%.1f%%", ewp),
      status = status, stringsAsFactors = FALSE
    )))
  }

  if (!is.null(out$balance)) {
    max_smd <- max(abs(out$balance$smd))
    status <- if (max_smd < 0.1) "Pass" else if (max_smd < 0.25) "Warning" else "Fail"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "Balance (max SMD)",
      value = sprintf("%.3f", max_smd),
      status = status, stringsAsFactors = FALSE
    )))
  }

  if (!is.null(out$cate) && !out$cate$constant_ite) {
    het <- if (out$cate$heterogeneity_detected) "Detected" else "Not detected"
    summary_rows <- c(summary_rows, list(data.frame(
      diagnostic = "CATE heterogeneity",
      value = het,
      status = "Info", stringsAsFactors = FALSE
    )))
  }

  out$summary <- do.call(rbind, summary_rows)

  class(out) <- c("cf_diagnostics_report", "list")
  out
}

#' @export
print.cf_diagnostics_report <- function(x, ...) {
  cat("\n")
  cat(cli::rule(left = "cfomics Diagnostics Report"))
  cat("\n")
  cat(sprintf("Method: %s | n = %d | p = %d\n\n", x$method, x$n, x$p))

  if (!is.null(x$summary)) {
    for (i in seq_len(nrow(x$summary))) {
      row <- x$summary[i, ]
      icon <- switch(row$status,
        Pass = cli::col_green(cli::symbol$tick),
        Warning = cli::col_yellow("!"),
        Fail = cli::col_red(cli::symbol$cross),
        Info = cli::col_blue("-")
      )
      cat(sprintf("%s %s: %s\n", icon, row$diagnostic, row$value))
    }
  }

  invisible(x)
}

#' @export
plot.cf_diagnostics_report <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for report plots")
    return(invisible(NULL))
  }

  plots <- list()

  if (!is.null(x$sensitivity)) {
    plots$sensitivity <- plot(x$sensitivity)
  }
  if (!is.null(x$model)) {
    plot(x$model)
  }
  if (!is.null(x$cate) && !x$cate$constant_ite) {
    plots$cate <- plot(x$cate)
  }
  if (!is.null(x$balance)) {
    plots$balance <- plot(x$balance)
  }

  invisible(NULL)
}
