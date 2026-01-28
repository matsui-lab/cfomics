# benchmarks/R/metrics.R
# Statistical tests for benchmark method comparison
# Uses Friedman test with Nemenyi post-hoc for nonparametric comparisons
#
# Dependencies: PMCMRplus package for Nemenyi test

#' Compute ranks per scenario
#'
#' For each scenario, ranks the methods based on the specified metric.
#' This is the first step for Friedman/Nemenyi nonparametric comparison.
#'
#' @param df Data frame with columns: grouping column(s), method, and metric column
#' @param metric Character, name of metric column to rank
#' @param groups Character vector of column names to group by (default "scenario_id")
#' @param lower_is_better Logical, TRUE if lower metric values are better (default TRUE for RMSE)
#' @return Data frame with grouping column(s), method, rank columns
#' @export
compute_ranks <- function(df, metric, groups = "scenario_id", lower_is_better = TRUE) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame")
  }

  # Validate groups parameter
  if (!all(groups %in% names(df))) {
    stop("'groups' columns not found in data frame: ",
         paste(setdiff(groups, names(df)), collapse = ", "))
  }

  if (!"method" %in% names(df)) {
    stop("'df' must contain 'method' column")
  }
  if (!metric %in% names(df)) {
    stop("Metric '", metric, "' not found in data frame")
  }

  # Use first column of groups for scenario-level grouping
  group_col <- groups[1]
  group_values <- unique(df[[group_col]])
  methods <- unique(df$method)

  ranks_list <- lapply(group_values, function(g) {
    sub_df <- df[df[[group_col]] == g, ]

    # Handle missing methods for this group
    if (nrow(sub_df) == 0) {
      return(NULL)
    }

    values <- sub_df[[metric]]

    # Handle NA values by assigning worst rank
    na_mask <- is.na(values)
    if (all(na_mask)) {
      # All NA: assign average rank
      r <- rep((length(values) + 1) / 2, length(values))
    } else {
      # Rank non-NA values, NA gets worst rank
      if (lower_is_better) {
        # Lower is better: smaller values get lower (better) ranks
        r <- rank(values, ties.method = "average", na.last = "keep")
        # Assign NA values the worst rank (highest)
        r[na_mask] <- length(values)  # Worst rank
      } else {
        # Higher is better: negate values so higher gets lower rank
        r <- rank(-values, ties.method = "average", na.last = "keep")
        r[na_mask] <- length(values)  # Worst rank
      }
    }

    # Build result data frame with the group column name
    result <- data.frame(
      group_value = g,
      method = sub_df$method,
      rank = r,
      stringsAsFactors = FALSE
    )
    names(result)[1] <- group_col
    result
  })

  # Remove NULL entries and combine
  ranks_list <- Filter(Negate(is.null), ranks_list)
  do.call(rbind, ranks_list)
}

#' Friedman test for comparing methods across scenarios
#'
#' Performs the Friedman test for comparing multiple methods across multiple
#' scenarios. This is a nonparametric alternative to repeated measures ANOVA.
#'
#' @param ranks_df Data frame from compute_ranks() with grouping column, method, rank columns
#' @param group_col Character, name of the grouping column (default "scenario_id")
#' @return List with statistic, df, p_value, mean_ranks, n_scenarios, n_methods
#' @export
friedman_test <- function(ranks_df, group_col = "scenario_id") {
  if (!is.data.frame(ranks_df)) {
    stop("'ranks_df' must be a data.frame")
  }
  if (!group_col %in% names(ranks_df)) {
    stop("'", group_col, "' column not found in ranks_df")
  }
  if (!all(c("method", "rank") %in% names(ranks_df))) {
    stop("'ranks_df' must contain 'method' and 'rank' columns")
  }

  methods <- unique(ranks_df$method)
  scenarios <- unique(ranks_df[[group_col]])

  n <- length(scenarios)
  k <- length(methods)

  if (n < 2) {
    stop("Friedman test requires at least 2 scenarios (blocks)")
  }
  if (k < 2) {
    stop("Friedman test requires at least 2 methods (treatments)")
  }

  # Build rank matrix: rows = scenarios, cols = methods
  rank_matrix <- matrix(NA_real_, nrow = n, ncol = k)
  colnames(rank_matrix) <- methods
  rownames(rank_matrix) <- scenarios

  for (i in seq_along(scenarios)) {
    for (j in seq_along(methods)) {
      sub <- ranks_df[ranks_df[[group_col]] == scenarios[i] &
                      ranks_df$method == methods[j], ]
      if (nrow(sub) == 1) {
        rank_matrix[i, j] <- sub$rank
      }
    }
  }

  # Check for complete data
  n_complete <- sum(complete.cases(rank_matrix))
  if (n_complete < n) {
    warning(sprintf("Only %d of %d scenarios have complete data for all methods",
                    n_complete, n))
  }

  # Compute mean ranks per method
  mean_ranks <- colMeans(rank_matrix, na.rm = TRUE)

  # Friedman chi-squared statistic
  # Formula: chi^2 = (12*n / (k*(k+1))) * sum((R_j - (k+1)/2)^2)
  chi_sq <- (12 * n / (k * (k + 1))) * sum((mean_ranks - (k + 1) / 2)^2)

  # Degrees of freedom
  df <- k - 1

  # P-value from chi-squared distribution
  p_value <- stats::pchisq(chi_sq, df = df, lower.tail = FALSE)

  list(
    statistic = chi_sq,
    df = df,
    p_value = p_value,
    mean_ranks = mean_ranks,
    n_scenarios = n,
    n_methods = k,
    rank_matrix = rank_matrix
  )
}

#' Compute critical difference for CD diagrams
#'
#' Calculates the critical difference value used in Critical Difference (CD)
#' diagrams for visualizing method comparisons. Two methods are significantly
#' different if their average rank difference exceeds the CD.
#'
#' @param n_methods Number of methods being compared
#' @param n_scenarios Number of scenarios (blocks)
#' @param alpha Significance level (default 0.05)
#' @return Critical difference value
#' @export
compute_critical_difference <- function(n_methods, n_scenarios, alpha = 0.05) {
  if (n_methods < 2) {
    stop("n_methods must be at least 2")
  }
  if (n_scenarios < 1) {
    stop("n_scenarios must be at least 1")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }

  # Critical value from Studentized Range distribution
  # For Nemenyi test: q_alpha from Tukey distribution
  q_alpha <- stats::qtukey(1 - alpha, n_methods, Inf) / sqrt(2)

  # Critical difference formula
  cd <- q_alpha * sqrt(n_methods * (n_methods + 1) / (6 * n_scenarios))

  cd
}

#' Nemenyi post-hoc test for pairwise comparisons
#'
#' Performs the Nemenyi post-hoc test after a significant Friedman test.
#' This provides pairwise comparisons between all methods with family-wise
#' error rate control.
#'
#' @param ranks_df Data frame from compute_ranks() with grouping column, method, rank columns
#' @param alpha Significance level (default 0.05)
#' @param group_col Character, name of the grouping column (default "scenario_id")
#' @return List with pvalues matrix, mean_ranks, critical_difference, and significant_pairs
#' @export
nemenyi_posthoc <- function(ranks_df, alpha = 0.05, group_col = "scenario_id") {
  if (!is.data.frame(ranks_df)) {
    stop("'ranks_df' must be a data.frame")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
  if (!group_col %in% names(ranks_df)) {
    stop("'", group_col, "' column not found in ranks_df")
  }
  if (!all(c("method", "rank") %in% names(ranks_df))) {
    stop("'ranks_df' must contain 'method' and 'rank' columns")
  }

  methods <- unique(ranks_df$method)
  scenarios <- unique(ranks_df[[group_col]])

  k <- length(methods)
  n <- length(scenarios)

  if (k < 2) {
    stop("Nemenyi test requires at least 2 methods")
  }

  # Build rank matrix: rows = scenarios, cols = methods
  rank_matrix <- matrix(NA_real_, nrow = n, ncol = k)
  colnames(rank_matrix) <- methods
  rownames(rank_matrix) <- scenarios

  for (i in seq_along(scenarios)) {
    for (j in seq_along(methods)) {
      sub <- ranks_df[ranks_df[[group_col]] == scenarios[i] &
                      ranks_df$method == methods[j], ]
      if (nrow(sub) == 1) {
        rank_matrix[i, j] <- sub$rank
      }
    }
  }

  # Compute mean ranks
  mean_ranks <- colMeans(rank_matrix, na.rm = TRUE)

  # Compute critical difference
  cd <- compute_critical_difference(k, n, alpha)

  # Check if PMCMRplus is available for full Nemenyi test
  has_pmcmrplus <- requireNamespace("PMCMRplus", quietly = TRUE)

  if (has_pmcmrplus) {
    # Use PMCMRplus for Nemenyi test
    nemenyi_result <- PMCMRplus::frdAllPairsNemenyiTest(rank_matrix)

    # Extract p-values matrix (lower triangular, NA on diagonal and upper)
    pvalues <- nemenyi_result$p.value

    # Identify significant pairs from lower triangular matrix
    # PMCMRplus returns a lower triangular matrix with NA elsewhere
    sig_pairs <- which(pvalues < alpha & !is.na(pvalues), arr.ind = TRUE)
    if (nrow(sig_pairs) > 0) {
      # In lower triangular: row index > col index
      # rownames = second method (from index 2 onward), colnames = first method (from index 1 onward)
      significant_pairs <- data.frame(
        method1 = colnames(pvalues)[sig_pairs[, 2]],
        method2 = rownames(pvalues)[sig_pairs[, 1]],
        p_value = pvalues[sig_pairs],
        stringsAsFactors = FALSE
      )
    } else {
      significant_pairs <- data.frame(
        method1 = character(0),
        method2 = character(0),
        p_value = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    list(
      pvalues = pvalues,
      mean_ranks = mean_ranks,
      critical_difference = cd,
      significant_pairs = significant_pairs,
      n_scenarios = n,
      n_methods = k,
      alpha = alpha
    )
  } else {
    # Fallback: compute based on mean rank differences and CD
    warning("PMCMRplus package not available. Using simplified comparison based on critical difference.")

    # Create pairwise comparison matrix based on mean rank differences
    rank_diffs <- outer(mean_ranks, mean_ranks, function(x, y) abs(x - y))
    rownames(rank_diffs) <- methods
    colnames(rank_diffs) <- methods

    # Determine significant pairs based on CD
    sig_matrix <- rank_diffs > cd
    diag(sig_matrix) <- FALSE

    sig_pairs <- which(sig_matrix, arr.ind = TRUE)
    if (nrow(sig_pairs) > 0) {
      significant_pairs <- data.frame(
        method1 = methods[sig_pairs[, 1]],
        method2 = methods[sig_pairs[, 2]],
        rank_diff = rank_diffs[sig_pairs],
        stringsAsFactors = FALSE
      )
      # Remove duplicates
      significant_pairs <- significant_pairs[significant_pairs$method1 < significant_pairs$method2, ]
    } else {
      significant_pairs <- data.frame(
        method1 = character(0),
        method2 = character(0),
        rank_diff = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    list(
      rank_differences = rank_diffs,
      mean_ranks = mean_ranks,
      critical_difference = cd,
      significant_pairs = significant_pairs,
      n_scenarios = n,
      n_methods = k,
      alpha = alpha
    )
  }
}

#' Summarize statistical test results
#'
#' Prints a formatted summary of Friedman test and Nemenyi post-hoc results.
#'
#' @param friedman_result Result from friedman_test()
#' @param nemenyi_result Result from nemenyi_posthoc() (optional)
#' @return Invisibly returns NULL
#' @export
summarize_statistical_tests <- function(friedman_result, nemenyi_result = NULL) {
  message("\n=== Friedman Test Results ===")
  message(sprintf("Chi-squared statistic: %.4f", friedman_result$statistic))
  message(sprintf("Degrees of freedom: %d", friedman_result$df))
  message(sprintf("P-value: %.6f", friedman_result$p_value))
  message(sprintf("Number of scenarios: %d", friedman_result$n_scenarios))
  message(sprintf("Number of methods: %d", friedman_result$n_methods))

  # Interpretation
  if (friedman_result$p_value < 0.05) {
    message("\nConclusion: Significant differences exist between methods (p < 0.05)")
  } else {
    message("\nConclusion: No significant differences between methods (p >= 0.05)")
  }

  # Mean ranks (sorted)
  message("\n--- Mean Ranks (lower is better) ---")
  sorted_ranks <- sort(friedman_result$mean_ranks)
  for (i in seq_along(sorted_ranks)) {
    message(sprintf("  %d. %s: %.3f", i, names(sorted_ranks)[i], sorted_ranks[i]))
  }

  if (!is.null(nemenyi_result)) {
    message("\n=== Nemenyi Post-hoc Test Results ===")
    message(sprintf("Critical Difference (alpha=%.2f): %.4f",
                    nemenyi_result$alpha, nemenyi_result$critical_difference))

    if (nrow(nemenyi_result$significant_pairs) > 0) {
      message("\n--- Significantly Different Pairs ---")
      for (i in seq_len(nrow(nemenyi_result$significant_pairs))) {
        pair <- nemenyi_result$significant_pairs[i, ]
        if ("p_value" %in% names(pair)) {
          message(sprintf("  %s vs %s (p = %.4f)", pair$method1, pair$method2, pair$p_value))
        } else {
          message(sprintf("  %s vs %s (rank diff = %.4f)", pair$method1, pair$method2, pair$rank_diff))
        }
      }
    } else {
      message("\nNo significantly different pairs found at alpha = ", nemenyi_result$alpha)
    }
  }

  invisible(NULL)
}
