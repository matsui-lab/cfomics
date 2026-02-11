#' Generate Publication-Ready Results Table
#'
#' Creates a formatted table of causal inference results suitable for
#' publication or reporting.
#'
#' @param fit A cf_model object from cf_fit()
#' @param format Output format: "data.frame" (default), "markdown", "latex", "html"
#' @param digits Number of decimal places (default 3)
#' @param include_ite Logical, whether to include ITE statistics (default TRUE)
#' @return Formatted table (data.frame or character string)
#' @export
#' @examples
#' # Create sample data and fit model
#' set.seed(123)
#' demo_data <- data.frame(
#'   X1 = rnorm(50),
#'   X2 = rnorm(50),
#'   T = rbinom(50, 1, 0.5)
#' )
#' demo_data$Y <- 1.5 * demo_data$T + 0.3 * demo_data$X1 + rnorm(50)
#' fit <- cf_fit(Y ~ T | X1 + X2, data = demo_data, method = "gformula")
#'
#' # Generate results table
#' cf_table(fit)
#'
#' \donttest{
#' # Markdown format (for reports)
#' cf_table(fit, format = "markdown")
#' }
cf_table <- function(fit, format = c("data.frame", "markdown", "latex", "html"),
                     digits = 3, include_ite = TRUE) {
  format <- match.arg(format)

  if (!inherits(fit, "cf_model")) {
    stop("fit must be a cf_model object from cf_fit()", call. = FALSE)
  }

  # Extract results
  res <- fit$fit$res
  meta <- fit$meta

  # Build summary table
  rows <- list(
    c("Method", fit$method),
    c("Sample size (N)", meta$n),
    c("Covariates (p)", meta$p),
    c("Outcome", meta$outcome_name),
    c("Treatment", meta$treatment_name)
  )

  # ATE
  ate_val <- round(res$ate, digits)
  rows <- c(rows, list(c("ATE", ate_val)))

  # Confidence interval if available
  if (!is.null(res$summary) && !is.null(res$summary$ate_ci_lower)) {
    ci_str <- sprintf("[%.*f, %.*f]",
                      digits, res$summary$ate_ci_lower,
                      digits, res$summary$ate_ci_upper)
    rows <- c(rows, list(c("ATE 95% CI", ci_str)))
  }

  # ITE statistics
  if (include_ite && !is.null(res$summary)) {
    if (!is.null(res$summary$ite_mean)) {
      rows <- c(rows, list(
        c("ITE Mean", round(res$summary$ite_mean, digits)),
        c("ITE SD", round(res$summary$ite_std, digits))
      ))
    }

    if (!is.null(res$summary$ite_quantiles)) {
      q <- res$summary$ite_quantiles
      rows <- c(rows, list(
        c("ITE 5th percentile", round(q$q05, digits)),
        c("ITE Median", round(q$q50, digits)),
        c("ITE 95th percentile", round(q$q95, digits))
      ))
    }
  }

  # Convert to data.frame
  summary_df <- data.frame(
    Metric = vapply(rows, `[`, character(1), 1),
    Value = vapply(rows, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )

  if (format == "data.frame") {
    return(summary_df)
  }

  # Format as table
  if (!requireNamespace("knitr", quietly = TRUE)) {
    message("knitr package not available. Returning data.frame.")
    return(summary_df)
  }

  knitr::kable(summary_df, format = format, row.names = FALSE)
}

#' Export Results to File
#'
#' Exports causal inference results to various file formats.
#'
#' @param fit A cf_model object from cf_fit()
#' @param file Output file path
#' @param format File format: auto-detected from extension, or specify
#'   "csv", "rds", "json"
#' @param include_ite Logical, whether to include individual treatment effects
#' @param overwrite Logical, whether to overwrite existing files. If FALSE
#'   (default), an error is raised when the file already exists.
#' @return Invisibly returns the file path
#' @export
#' @examples
#' # Fit a model
#' set.seed(123)
#' demo_data <- data.frame(
#'   X = rnorm(30),
#'   T = rbinom(30, 1, 0.5)
#' )
#' demo_data$Y <- 2 * demo_data$T + rnorm(30)
#' fit <- cf_fit(Y ~ T | X, data = demo_data, method = "gformula")
#'
#' # Export to temporary file (RDS format - recommended for R)
#' tmp_file <- tempfile(fileext = ".rds")
#' cf_export(fit, tmp_file, format = "rds")
#' unlink(tmp_file)  # Clean up
#'
#' \donttest{
#' # Export to CSV (creates two files: data and summary)
#' tmp_csv <- tempfile(fileext = ".csv")
#' cf_export(fit, tmp_csv, format = "csv")
#' unlink(tmp_csv)
#' unlink(sub("\\.csv$", "_summary.csv", tmp_csv))
#' }
cf_export <- function(fit, file, format = NULL, include_ite = TRUE,
                      overwrite = FALSE) {
  if (!inherits(fit, "cf_model")) {
    stop("fit must be a cf_model object from cf_fit()", call. = FALSE)
  }

  # Auto-detect format from extension
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(file))
    format <- switch(ext,
      csv = "csv",
      rds = "rds",
      json = "json",
      "csv"  # default
    )
  }

  # Check for existing files (including summary file for CSV)
  files_to_check <- file
  if (format == "csv") {
    summary_file <- sub("\\.csv$", "_summary.csv", file)
    files_to_check <- c(file, summary_file)
  }

  existing_files <- files_to_check[file.exists(files_to_check)]
  if (length(existing_files) > 0 && !overwrite) {
    stop(
      sprintf(
        "File(s) already exist: %s\nUse overwrite = TRUE to replace existing files.",
        paste(existing_files, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Prepare results
  res <- fit$fit$res

  # Convert meta for JSON serialization (formula -> character)
  meta_json <- fit$meta
  if (!is.null(meta_json$formula)) {
    meta_json$formula <- deparse(meta_json$formula)
  }

  results <- list(
    method = fit$method,
    ate = res$ate,
    meta = meta_json
  )

  # Add confidence intervals if available
  if (!is.null(res$summary)) {
    results$ate_ci <- c(
      lower = res$summary$ate_ci_lower,
      upper = res$summary$ate_ci_upper
    )
    results$ite_summary <- list(
      mean = res$summary$ite_mean,
      sd = res$summary$ite_std,
      quantiles = res$summary$ite_quantiles
    )
  }

  if (include_ite && !is.null(res$ite)) {
    results$ite <- res$ite
    results$y0_hat <- res$y0_hat
    results$y1_hat <- res$y1_hat
  }

  # Export based on format
  if (format == "csv") {
    # Create a flat data frame for CSV
    df <- data.frame(
      sample_id = seq_along(res$ite),
      ite = res$ite
    )
    if (!is.null(res$y0_hat)) df$y0_hat <- res$y0_hat
    if (!is.null(res$y1_hat)) df$y1_hat <- res$y1_hat

    utils::write.csv(df, file, row.names = FALSE)

    # Also save summary to separate file
    summary_file <- sub("\\.csv$", "_summary.csv", file)
    summary_df <- cf_table(fit, format = "data.frame")
    utils::write.csv(summary_df, summary_file, row.names = FALSE)

    message(sprintf("Results exported to:\n  - %s (ITE)\n  - %s (summary)",
                    file, summary_file))
  } else if (format == "rds") {
    saveRDS(results, file)
    message(sprintf("Results exported to: %s", file))
  } else if (format == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("jsonlite package required for JSON export", call. = FALSE)
    }
    jsonlite::write_json(results, file, auto_unbox = TRUE, pretty = TRUE)
    message(sprintf("Results exported to: %s", file))
  }

  invisible(file)
}

#' Generate Session Information for Reproducibility
#'
#' Creates a summary of the R session, including package versions,
#' for reproducibility documentation.
#'
#' @param fit Optional cf_model object to include fit-specific info
#' @return A list with session information
#' @export
#' @examples
#' # Get current session info
#' session <- cf_session_info()
#' print(session)
#'
#' \donttest{
#' # Include model fit information
#' set.seed(123)
#' demo_data <- data.frame(X = rnorm(30), T = rbinom(30, 1, 0.5))
#' demo_data$Y <- 2 * demo_data$T + rnorm(30)
#' fit <- cf_fit(Y ~ T | X, data = demo_data, method = "gformula")
#' session_with_fit <- cf_session_info(fit)
#' }
cf_session_info <- function(fit = NULL) {
  info <- list(
    timestamp = Sys.time(),
    cfomics_version = tryCatch(
      as.character(utils::packageVersion("cfomics")),
      error = function(e) "unknown"
    ),
    r_version = R.version.string,
    platform = R.version$platform,
    packages = list(
      cfomics = tryCatch(utils::packageVersion("cfomics"), error = function(e) NA),
      reticulate = tryCatch(utils::packageVersion("reticulate"), error = function(e) NA),
      grf = tryCatch(utils::packageVersion("grf"), error = function(e) NA)
    )
  )

  # Python info if available
  if (requireNamespace("reticulate", quietly = TRUE)) {
    info$python <- tryCatch({
      if (reticulate::py_available(initialize = FALSE)) {
        config <- reticulate::py_config()
        list(
          version = config$version,
          path = config$python
        )
      } else {
        list(available = FALSE)
      }
    }, error = function(e) list(available = FALSE))
  }

  # Fit-specific info
  if (!is.null(fit) && inherits(fit, "cf_model")) {
    info$fit <- list(
      method = fit$method,
      n = fit$meta$n,
      p = fit$meta$p,
      random_state = fit$meta$random_state
    )
  }

  class(info) <- c("cf_session_info", "list")
  info
}

#' Print method for cf_session_info
#' @param x cf_session_info object
#' @param ... Additional arguments (unused)
#' @export
print.cf_session_info <- function(x, ...) {
  cat("\n")
  cat("cfomics Session Information\n")
  cat("============================\n")
  cat(sprintf("Timestamp: %s\n", x$timestamp))
  cat(sprintf("cfomics version: %s\n", x$cfomics_version))
  cat(sprintf("R version: %s\n", x$r_version))
  cat(sprintf("Platform: %s\n", x$platform))

  if (!is.null(x$python) && isTRUE(x$python$available) || !is.null(x$python$version)) {
    cat(sprintf("Python: %s\n", x$python$version))
  }

  if (!is.null(x$fit)) {
    cat("\n")
    cat("Fit Information:\n")
    cat(sprintf("  Method: %s\n", x$fit$method))
    cat(sprintf("  N: %d, p: %d\n", x$fit$n, x$fit$p))
    cat(sprintf("  Random state: %d\n", x$fit$random_state))
  }

  invisible(x)
}
