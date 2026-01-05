#' Check Covariate Balance Between Treatment Groups
#'
#' Computes standardized mean differences (SMD) to assess balance of covariates
#' between treatment and control groups. Useful for checking the plausibility
#' of causal inference assumptions.
#'
#' @param data Data frame or SummarizedExperiment object
#' @param treatment Name of the treatment variable
#' @param covariates Character vector of covariate names to check. If NULL,
#'   all numeric columns (except treatment) are used.
#' @param threshold SMD threshold for imbalance warning (default 0.1)
#' @return A data.frame of class "cf_balance" with balance statistics
#' @export
#' @examples
#' \dontrun{
#' # Check balance
#' balance <- cf_balance_check(my_data, "treatment", c("age", "sex"))
#' print(balance)
#' plot(balance)
#' }
cf_balance_check <- function(data, treatment, covariates = NULL, threshold = 0.1) {
  # Handle Bioconductor objects
  if (is_se(data)) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package required", call. = FALSE)
    }
    data <- as.data.frame(SummarizedExperiment::colData(data))
  }

  if (!(treatment %in% names(data))) {
    stop("Treatment variable '", treatment, "' not found in data", call. = FALSE)
  }

  # Get covariates
  if (is.null(covariates)) {
    numeric_cols <- sapply(data, is.numeric)
    covariates <- names(data)[numeric_cols]
    covariates <- setdiff(covariates, treatment)
  }

  if (length(covariates) == 0) {
    stop("No numeric covariates found to check", call. = FALSE)
  }

  T_vec <- data[[treatment]]

  # Compute SMD for each covariate
  results <- lapply(covariates, function(var) {
    x <- data[[var]]
    if (!is.numeric(x)) return(NULL)

    # Remove NAs
    valid <- !is.na(x) & !is.na(T_vec)
    x <- x[valid]
    t <- T_vec[valid]

    # Compute means and SDs
    mean_t1 <- mean(x[t == 1], na.rm = TRUE)
    mean_t0 <- mean(x[t == 0], na.rm = TRUE)
    var_t1 <- stats::var(x[t == 1], na.rm = TRUE)
    var_t0 <- stats::var(x[t == 0], na.rm = TRUE)

    # Pooled SD
    sd_pooled <- sqrt((var_t1 + var_t0) / 2)

    # Standardized mean difference
    if (sd_pooled > 0) {
      smd <- (mean_t1 - mean_t0) / sd_pooled
    } else {
      smd <- 0
    }

    data.frame(
      variable = var,
      mean_treated = round(mean_t1, 4),
      mean_control = round(mean_t0, 4),
      diff = round(mean_t1 - mean_t0, 4),
      smd = round(smd, 4),
      balanced = abs(smd) < threshold,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, Filter(Negate(is.null), results))
  rownames(result) <- NULL

  class(result) <- c("cf_balance", "data.frame")
  attr(result, "threshold") <- threshold
  attr(result, "treatment") <- treatment

  result
}

#' Print method for cf_balance
#' @param x cf_balance object
#' @param ... Additional arguments (unused)
#' @export
print.cf_balance <- function(x, ...) {
  threshold <- attr(x, "threshold")
  treatment <- attr(x, "treatment")
  n_imbalanced <- sum(!x$balanced)

  cat("\n")
  cat("Covariate Balance Check\n")
  cat("========================\n")
  cat(sprintf("Treatment: %s\n", treatment))
  cat(sprintf("Threshold: |SMD| < %.2f\n", threshold))
  cat(sprintf("Variables: %d total, %d imbalanced\n\n",
              nrow(x), n_imbalanced))

  # Print table
  print.data.frame(x, row.names = FALSE)

  if (n_imbalanced > 0) {
    cat("\n")
    cat("Warning: The following variables show imbalance:\n")
    imbalanced <- x$variable[!x$balanced]
    cat(paste("  -", imbalanced, collapse = "\n"), "\n")
    cat("\nConsider adjusting for these variables or using matching methods.\n")
  } else {
    cat("\nAll covariates are balanced.\n")
  }

  invisible(x)
}

#' Plot Covariate Balance
#' @param x cf_balance object
#' @param ... Additional arguments passed to plot
#' @export
plot.cf_balance <- function(x, ...) {
  threshold <- attr(x, "threshold")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(x, ggplot2::aes(
      x = smd,
      y = stats::reorder(variable, abs(smd))
    )) +
      ggplot2::geom_point(ggplot2::aes(color = balanced), size = 3) +
      ggplot2::geom_vline(xintercept = c(-threshold, threshold),
                          linetype = "dashed", color = "red", alpha = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", alpha = 0.5) +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "forestgreen", "FALSE" = "red"),
        labels = c("TRUE" = "Balanced", "FALSE" = "Imbalanced")
      ) +
      ggplot2::labs(
        x = "Standardized Mean Difference (SMD)",
        y = "Covariate",
        title = "Covariate Balance Plot",
        color = "Status"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

    return(p)
  } else {
    # Base R plot
    old_par <- graphics::par(mar = c(5, 10, 4, 2))
    on.exit(graphics::par(old_par))

    colors <- ifelse(x$balanced, "forestgreen", "red")
    graphics::barplot(
      x$smd,
      names.arg = x$variable,
      horiz = TRUE,
      col = colors,
      main = "Covariate Balance",
      xlab = "Standardized Mean Difference",
      las = 1
    )
    graphics::abline(v = c(-threshold, threshold), lty = 2, col = "red")
    graphics::abline(v = 0)
  }

  invisible(NULL)
}

#' Check Overlap/Positivity Assumption
#'
#' Examines the distribution of propensity scores to assess whether the
#' positivity assumption holds (sufficient overlap between treatment groups).
#'
#' @param data Data frame or SummarizedExperiment
#' @param treatment Name of treatment variable
#' @param covariates Names of covariates for propensity score model
#' @param method Propensity score estimation method: "logistic" (default)
#' @return A list of class "cf_overlap" with overlap diagnostics
#' @export
#' @examples
#' \dontrun{
#' overlap <- cf_overlap_check(my_data, "treatment", c("age", "sex"))
#' print(overlap)
#' plot(overlap)
#' }
cf_overlap_check <- function(data, treatment, covariates,
                             method = c("logistic")) {
  method <- match.arg(method)

  # Handle Bioconductor objects
  if (is_se(data)) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment package required", call. = FALSE)
    }
    data <- as.data.frame(SummarizedExperiment::colData(data))
  }

  # Fit propensity score model
  formula <- stats::as.formula(
    paste(treatment, "~", paste(covariates, collapse = " + "))
  )

  ps_model <- stats::glm(formula, data = data, family = stats::binomial())
  ps <- stats::predict(ps_model, type = "response")

  T_vec <- data[[treatment]]

  # Summary statistics by group
  ps_t1 <- ps[T_vec == 1]
  ps_t0 <- ps[T_vec == 0]

  result <- list(
    propensity_scores = ps,
    treatment = T_vec,
    summary = data.frame(
      group = c("Control (T=0)", "Treated (T=1)"),
      n = c(sum(T_vec == 0), sum(T_vec == 1)),
      min_ps = c(min(ps_t0), min(ps_t1)),
      max_ps = c(max(ps_t0), max(ps_t1)),
      mean_ps = c(mean(ps_t0), mean(ps_t1)),
      median_ps = c(stats::median(ps_t0), stats::median(ps_t1))
    ),
    overlap_region = c(
      lower = max(min(ps_t0), min(ps_t1)),
      upper = min(max(ps_t0), max(ps_t1))
    ),
    model = ps_model
  )

  # Check for violations
  result$violations <- list(
    no_overlap = result$overlap_region["upper"] < result$overlap_region["lower"],
    extreme_ps = any(ps < 0.01 | ps > 0.99)
  )

  class(result) <- c("cf_overlap", "list")
  result
}

#' Print method for cf_overlap
#' @param x cf_overlap object
#' @param ... Additional arguments (unused)
#' @export
print.cf_overlap <- function(x, ...) {
  cat("\n")
  cat("Overlap/Positivity Check\n")
  cat("=========================\n\n")

  cat("Propensity Score Summary by Group:\n")
  print(x$summary, row.names = FALSE)

  cat("\n")
  cat(sprintf("Overlap region: [%.4f, %.4f]\n",
              x$overlap_region["lower"], x$overlap_region["upper"]))

  if (x$violations$no_overlap) {
    cat("\n[WARNING] No overlap between treatment groups!\n")
    cat("Causal inference may not be valid.\n")
  }

  if (x$violations$extreme_ps) {
    n_extreme <- sum(x$propensity_scores < 0.01 | x$propensity_scores > 0.99)
    cat(sprintf("\n[WARNING] %d observations have extreme propensity scores (< 0.01 or > 0.99)\n",
                n_extreme))
    cat("Consider trimming or other adjustments.\n")
  }

  if (!x$violations$no_overlap && !x$violations$extreme_ps) {
    cat("\nPositivity assumption appears satisfied.\n")
  }

  invisible(x)
}

#' Plot Overlap Diagnostics
#' @param x cf_overlap object
#' @param ... Additional arguments
#' @export
plot.cf_overlap <- function(x, ...) {
  df <- data.frame(
    ps = x$propensity_scores,
    group = ifelse(x$treatment == 1, "Treated", "Control")
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = ps, fill = group)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::geom_vline(xintercept = x$overlap_region,
                          linetype = "dashed", alpha = 0.7) +
      ggplot2::scale_fill_manual(values = c("Control" = "blue", "Treated" = "red")) +
      ggplot2::labs(
        x = "Propensity Score",
        y = "Density",
        title = "Propensity Score Distribution",
        subtitle = sprintf("Overlap: [%.3f, %.3f]",
                           x$overlap_region["lower"], x$overlap_region["upper"]),
        fill = "Group"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")

    return(p)
  } else {
    # Base R plot
    graphics::hist(x$propensity_scores[x$treatment == 0],
                   col = rgb(0, 0, 1, 0.5), main = "Propensity Score Distribution",
                   xlab = "Propensity Score", xlim = c(0, 1), freq = FALSE)
    graphics::hist(x$propensity_scores[x$treatment == 1],
                   col = rgb(1, 0, 0, 0.5), add = TRUE, freq = FALSE)
    graphics::legend("topright", c("Control", "Treated"),
                     fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))
    graphics::abline(v = x$overlap_region, lty = 2)
  }

  invisible(NULL)
}
