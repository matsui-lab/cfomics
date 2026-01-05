#' Print method for cfomics_result
#' @param x A cfomics_result object
#' @param ... Additional arguments (unused)
#' @export
print.cfomics_result <- function(x, ...) {
  # Validate input
if (!inherits(x, c("cf_model", "cfomics_result"))) {
    rlang::abort(
      message = "x must be a cf_model from cf_fit()",
      class = "cfomics_invalid_input"
    )
  }

  cat("cfomics model\n")
  cat("  method   :", x$method, "\n")
  cat("  outcome  :", x$meta$outcome_name, "\n")
  cat("  treatment:", x$meta$treatment_name, "\n")
  cat("  n, p     :", x$meta$n, ",", x$meta$p, "\n")
  a <- tryCatch(predict(x, type = "ate"), error = function(e) NA)
  if (!is.na(a)) cat("  ATE      :", round(a, 4), "\n")
  invisible(x)
}

#' Plot method for cfomics_result
#'
#' Creates a histogram of individual treatment effects (ITE).
#'
#' @param x A cfomics_result object from cf_fit()
#' @param ... Additional arguments passed to plotting functions
#' @return A ggplot2 object if ggplot2 is available, otherwise NULL invisibly
#' @export
#' @examples
#' # Fit a model
#' set.seed(42)
#' demo_data <- data.frame(
#'   X = rnorm(50),
#'   T = rbinom(50, 1, 0.5)
#' )
#' demo_data$Y <- 1.5 * demo_data$T + 0.3 * demo_data$X + rnorm(50)
#' fit <- cf_fit(Y ~ T | X, data = demo_data, method = "gformula")
#'
#' # Plot ITE distribution
#' plot(fit)
plot.cfomics_result <- function(x, ...) {
  # Validate input
  if (!inherits(x, c("cf_model", "cfomics_result"))) {
    rlang::abort(
      message = "x must be a cf_model from cf_fit()",
      class = "cfomics_invalid_input"
    )
  }

  if (is.null(x$method) || is.null(x$fit)) {
    rlang::abort(
      message = "Invalid cf_model structure: missing method or fit",
      class = "cfomics_invalid_result"
    )
  }

  ite <- tryCatch(predict(x, type = "ite"), error = function(e) NULL)

  if (is.null(ite)) {
    message("No ITE available for plotting.")
    return(invisible(NULL))
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(ITE = ite)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = ITE)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Distribution of Individual Treatment Effects (ITE)",
                    x = "ITE", y = "Count")
    return(p)
  } else {
    hist(ite, main = "Distribution of ITE", xlab = "ITE", col = "steelblue", border = "white")
  }

  invisible(NULL)
}
