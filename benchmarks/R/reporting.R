# benchmarks/R/reporting.R
# Visualization and LaTeX table generation for benchmark reports
#
# Dependencies: ggplot2, patchwork, scales, kableExtra
#
# Functions:
#   - plot_method_heatmap()    : Heatmap of method x scenario performance
#   - plot_bias_boxplot()      : Boxplots comparing bias distributions
#   - plot_coverage_heatmap()  : Coverage heatmap (actual vs nominal 95%)
#   - plot_cd_diagram()        : Critical Difference diagram
#   - plot_sweep_curves()      : Parameter sensitivity curves
#   - plot_bias_variance()     : Bias-variance decomposition (stacked bars)
#   - plot_computation_time()  : Log-scale computation time comparison
#   - generate_latex_table()   : kableExtra LaTeX output

library(ggplot2)
library(scales)

# Color palette for consistent method coloring
METHOD_COLORS <- c(

  "gformula" = "#E69F00",
  "hdml"     = "#56B4E9",
  "hdps"     = "#009E73",

  "bcf"      = "#F0E442",
  "tmle"     = "#0072B2",
  "grf"      = "#D55E00",
  "ipw"      = "#CC79A7"
)

#' Method comparison heatmap
#'
#' Creates a heatmap showing method performance across scenarios for a given metric.
#' Values are displayed in each cell with appropriate color scaling.
#'
#' @param df Aggregated results data frame with columns: method, scenario_id, and metric column
#' @param metric Character, column name to visualize (e.g., "rmse_ate", "mean_pehe")
#' @param lower_is_better Logical, use reverse color scale if TRUE (default TRUE)
#' @param title Optional plot title (default uses metric name)
#' @return ggplot object
#' @export
plot_method_heatmap <- function(df, metric, lower_is_better = TRUE, title = NULL) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame")
  }
  if (nrow(df) == 0) {
    stop("'df' contains no rows")
  }
  if (!metric %in% names(df)) {
    stop("Metric '", metric, "' not found in data frame")
  }
  if (!all(c("method", "scenario_id") %in% names(df))) {
    stop("'df' must contain 'method' and 'scenario_id' columns")
  }

  # Handle NA values for midpoint calculation
  metric_vals <- df[[metric]]
  if (all(is.na(metric_vals))) {
    warning("All values for metric '", metric, "' are NA")
    midpoint <- 0
  } else {
    midpoint <- median(metric_vals, na.rm = TRUE)
  }

  # Create plot title
  if (is.null(title)) {
    title <- paste("Method Performance:", metric)
  }

  # Define color scale based on lower_is_better
  if (lower_is_better) {
    low_color <- "#2166AC"   # Blue (good)
    high_color <- "#B2182B"  # Red (bad)
  } else {
    low_color <- "#B2182B"   # Red (bad)
    high_color <- "#2166AC"  # Blue (good)
  }

  p <- ggplot(df, aes(x = method, y = scenario_id, fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", .data[[metric]])),
              size = 2.5, color = "black") +
    scale_fill_gradient2(
      low = low_color,
      mid = "white",
      high = high_color,
      midpoint = midpoint,
      na.value = "grey80"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      x = "Method",
      y = "Scenario",
      fill = metric,
      title = title
    )

  p
}

#' Bias distribution boxplots
#'
#' Creates boxplots comparing the distribution of bias across methods.
#' A dashed red line at zero indicates unbiased estimation.
#'
#' @param df Raw results data frame with columns: method, bias_ate (or specified bias_col)
#' @param bias_col Character, name of the bias column (default "bias_ate")
#' @param title Optional plot title
#' @return ggplot object
#' @export
plot_bias_boxplot <- function(df, bias_col = "bias_ate", title = NULL) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame")
  }
  if (nrow(df) == 0) {
    stop("'df' contains no rows")
  }
  if (!bias_col %in% names(df)) {
    stop("Bias column '", bias_col, "' not found in data frame")
  }
  if (!"method" %in% names(df)) {
    stop("'df' must contain 'method' column")
  }

  if (is.null(title)) {
    title <- "Bias Distribution by Method"
  }

  # Get available methods that exist in our color palette
  methods_in_data <- unique(df$method)
  colors_to_use <- METHOD_COLORS[intersect(methods_in_data, names(METHOD_COLORS))]
  # For methods not in palette, use grey
  missing_methods <- setdiff(methods_in_data, names(METHOD_COLORS))
  if (length(missing_methods) > 0) {
    extra_colors <- setNames(rep("grey50", length(missing_methods)), missing_methods)
    colors_to_use <- c(colors_to_use, extra_colors)
  }

  p <- ggplot(df, aes(x = method, y = .data[[bias_col]], fill = method)) +
    geom_boxplot(alpha = 0.8, outlier.size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_manual(values = colors_to_use, guide = "none") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Method",
      y = "Bias",
      title = title
    )

  p
}

#' Coverage heatmap
#'
#' Creates a heatmap showing actual coverage rates across methods and scenarios.
#' Highlights deviation from nominal 95% coverage with diverging color scale.
#'
#' @param df Aggregated results data frame with columns: method, scenario_id, coverage
#' @param coverage_col Character, name of the coverage column (default "coverage")
#' @param nominal Numeric, nominal coverage level (default 0.95)
#' @param title Optional plot title
#' @return ggplot object
#' @export
plot_coverage_heatmap <- function(df, coverage_col = "coverage", nominal = 0.95, title = NULL) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame")
  }
  if (nrow(df) == 0) {
    stop("'df' contains no rows")
  }
  if (!coverage_col %in% names(df)) {
    stop("Coverage column '", coverage_col, "' not found in data frame")
  }
  if (!all(c("method", "scenario_id") %in% names(df))) {
    stop("'df' must contain 'method' and 'scenario_id' columns")
  }

  if (is.null(title)) {
    title <- sprintf("Coverage (nominal = %.0f%%)", nominal * 100)
  }

  p <- ggplot(df, aes(x = method, y = scenario_id, fill = .data[[coverage_col]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", .data[[coverage_col]])),
              size = 2.5, color = "black") +
    scale_fill_gradient2(
      low = "#B2182B",     # Red for under-coverage
      mid = "#2166AC",     # Blue for nominal
      high = "#1B7837",    # Green for over-coverage
      midpoint = nominal,
      limits = c(0, 1),
      na.value = "grey80"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      x = "Method",
      y = "Scenario",
      fill = "Coverage",
      title = title
    )

  p
}

#' Critical Difference diagram
#'
#' Creates a Critical Difference (CD) diagram for visualizing statistical
#' comparisons between methods. Methods connected by a horizontal bar are
#' not significantly different. Based on Demsar (2006) methodology.
#'
#' @param friedman_result List output from friedman_test() containing mean_ranks
#' @param nemenyi_result List output from nemenyi_posthoc() containing critical_difference
#' @param title Optional plot title
#' @return ggplot object
#' @export
plot_cd_diagram <- function(friedman_result, nemenyi_result, title = NULL) {
  if (!is.list(friedman_result) || !"mean_ranks" %in% names(friedman_result)) {
    stop("'friedman_result' must be output from friedman_test()")
  }
  if (!is.list(nemenyi_result) || !"critical_difference" %in% names(nemenyi_result)) {
    stop("'nemenyi_result' must be output from nemenyi_posthoc()")
  }

  mean_ranks <- friedman_result$mean_ranks
  cd <- nemenyi_result$critical_difference

  if (is.null(names(mean_ranks))) {
    stop("'mean_ranks' in friedman_result must be a named numeric vector")
  }

  # Sort methods by rank (best = lowest rank first)
  sorted_idx <- order(mean_ranks)
  methods <- names(mean_ranks)[sorted_idx]
  ranks <- mean_ranks[sorted_idx]
  n <- length(methods)

  if (is.null(title)) {
    title <- sprintf("Critical Difference Diagram (CD = %.2f)", cd)
  }

  # Create data for plotting
  df <- data.frame(
    method = factor(methods, levels = rev(methods)),  # Reverse for bottom-to-top
    rank = ranks,
    y = seq_len(n),
    stringsAsFactors = FALSE
  )

  # Find groups of methods that are not significantly different
  # Methods are in the same group if their rank difference <= CD
  groups <- list()
  for (i in seq_len(n)) {
    for (j in seq(i, n)) {
      if (abs(ranks[i] - ranks[j]) <= cd) {
        # Find or create group containing method i
        found <- FALSE
        for (g in seq_along(groups)) {
          if (i %in% groups[[g]]) {
            groups[[g]] <- union(groups[[g]], j)
            found <- TRUE
            break
          }
        }
        if (!found) {
          groups[[length(groups) + 1]] <- c(i, j)
        }
      }
    }
  }

  # Merge overlapping groups
  merged <- TRUE
  while (merged) {
    merged <- FALSE
    for (i in seq_along(groups)) {
      if (i > length(groups)) break
      if (i < length(groups)) {
        for (j in (i + 1):length(groups)) {
          if (j > length(groups)) break
          if (length(intersect(groups[[i]], groups[[j]])) > 0) {
            groups[[i]] <- union(groups[[i]], groups[[j]])
            groups[[j]] <- NULL
            groups <- Filter(Negate(is.null), groups)
            merged <- TRUE
            break
          }
        }
      }
      if (merged) break
    }
  }

  # Build segments for groups (horizontal bars connecting methods)
  group_segments <- data.frame(
    x = numeric(0), xend = numeric(0),
    y = numeric(0), yend = numeric(0)
  )

  if (length(groups) > 0) {
    y_offset <- 0.3
    for (g in seq_along(groups)) {
      grp <- groups[[g]]
      if (length(grp) > 1) {
        grp_ranks <- ranks[grp]
        group_segments <- rbind(group_segments, data.frame(
          x = min(grp_ranks),
          xend = max(grp_ranks),
          y = n + y_offset + (g - 1) * 0.4,
          yend = n + y_offset + (g - 1) * 0.4
        ))
      }
    }
  }

  # Create plot
  p <- ggplot(df, aes(x = rank, y = y)) +
    # Method labels and rank markers
    geom_segment(aes(x = rank, xend = rank, y = y - 0.2, yend = y + 0.2),
                 linewidth = 1) +
    geom_text(aes(label = method), hjust = 1.1, vjust = 0.5, size = 3) +
    # CD indicator at top
    geom_segment(
      data = data.frame(x = 1, xend = 1 + cd, y = n + 0.8, yend = n + 0.8),
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 2, color = "darkred",
      inherit.aes = FALSE
    ) +
    annotate("text", x = 1 + cd / 2, y = n + 1.1,
             label = sprintf("CD = %.2f", cd), size = 3, fontface = "bold")

  # Add group connection bars if present

  if (nrow(group_segments) > 0) {
    p <- p + geom_segment(
      data = group_segments,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 2, color = "grey40",
      inherit.aes = FALSE
    )
  }

  p <- p +
    scale_x_continuous(
      limits = c(0, n + 0.5),
      breaks = seq(1, n),
      name = "Average Rank"
    ) +
    scale_y_continuous(limits = c(0.5, n + 1.5)) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    labs(title = title)

  p
}

#' Parameter sensitivity curves
#'
#' Creates line plots showing how metrics change across a parameter sweep.
#' Useful for visualizing sensitivity analysis results.
#'
#' @param df Data frame with columns: x_var (parameter), method, and metric columns
#' @param x_var Character, name of the x-axis variable (parameter being swept)
#' @param metrics Character vector of metric column names to plot
#' @param title Optional plot title
#' @return ggplot object (patchwork combined if multiple metrics)
#' @export
plot_sweep_curves <- function(df, x_var, metrics, title = NULL) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame")
  }
  if (nrow(df) == 0) {
    stop("'df' contains no rows")
  }
  if (!x_var %in% names(df)) {
    stop("x_var '", x_var, "' not found in data frame")
  }
  if (!all(metrics %in% names(df))) {
    missing <- setdiff(metrics, names(df))
    stop("Metrics not found in data frame: ", paste(missing, collapse = ", "))
  }
  if (!"method" %in% names(df)) {
    stop("'df' must contain 'method' column")
  }

  # Get available methods for coloring
  methods_in_data <- unique(df$method)
  colors_to_use <- METHOD_COLORS[intersect(methods_in_data, names(METHOD_COLORS))]
  missing_methods <- setdiff(methods_in_data, names(METHOD_COLORS))
  if (length(missing_methods) > 0) {
    extra_colors <- setNames(
      scales::hue_pal()(length(missing_methods)),
      missing_methods
    )
    colors_to_use <- c(colors_to_use, extra_colors)
  }

  # Create individual plots for each metric
  plots <- lapply(metrics, function(m) {
    ggplot(df, aes(x = .data[[x_var]], y = .data[[m]],
                   color = method, group = method)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      scale_color_manual(values = colors_to_use) +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      ) +
      labs(x = x_var, y = m)
  })

  # Combine plots using patchwork if available, otherwise return list
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = 1) +
      patchwork::plot_layout(guides = "collect") &
      theme(legend.position = "bottom")

    if (!is.null(title)) {
      combined <- combined + patchwork::plot_annotation(title = title)
    }
    return(combined)
  } else {
    warning("patchwork package not available. Returning list of plots.")
    return(plots)
  }
}

#' Bias-variance decomposition plot
#'
#' Creates stacked bar plots showing the decomposition of MSE into
#' bias squared and variance components for each method.
#'
#' @param df Aggregated results data frame with columns: method, mean_bias, sd_bias
#'           (or custom column names specified)
#' @param bias_col Character, name of mean bias column (default "mean_bias")
#' @param var_col Character, name of variance/sd column (default "sd_bias").
#'        If this contains standard deviation, it will be squared.
#' @param is_sd Logical, if TRUE (default), var_col is treated as SD and squared
#' @param title Optional plot title
#' @return ggplot object
#' @export
plot_bias_variance <- function(df, bias_col = "mean_bias", var_col = "sd_bias",
                               is_sd = TRUE, title = NULL) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame")
  }
  if (nrow(df) == 0) {
    stop("'df' contains no rows")
  }
  if (!bias_col %in% names(df)) {
    stop("Bias column '", bias_col, "' not found in data frame")
  }
  if (!var_col %in% names(df)) {
    stop("Variance column '", var_col, "' not found in data frame")
  }
  if (!"method" %in% names(df)) {
    stop("'df' must contain 'method' column")
  }

  if (is.null(title)) {
    title <- "Bias-Variance Decomposition"
  }

  # Compute bias^2 and variance
  bias_sq <- df[[bias_col]]^2
  variance <- if (is_sd) df[[var_col]]^2 else df[[var_col]]

  # Handle case where data has multiple scenarios - aggregate by method
  if ("scenario_id" %in% names(df)) {
    agg_df <- aggregate(
      cbind(bias_sq = bias_sq, variance = variance),
      by = list(method = df$method),
      FUN = function(x) mean(x, na.rm = TRUE)
    )
  } else {
    agg_df <- data.frame(
      method = df$method,
      bias_sq = bias_sq,
      variance = variance
    )
  }

  # Reshape to long format for stacking
  long_df <- rbind(
    data.frame(
      method = agg_df$method,
      component = "Bias^2",
      value = agg_df$bias_sq
    ),
    data.frame(
      method = agg_df$method,
      component = "Variance",
      value = agg_df$variance
    )
  )

  # Order methods by total MSE
  mse_order <- agg_df$bias_sq + agg_df$variance
  method_order <- agg_df$method[order(mse_order)]
  long_df$method <- factor(long_df$method, levels = method_order)

  p <- ggplot(long_df, aes(x = method, y = value, fill = component)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.9) +
    scale_fill_manual(
      values = c("Bias^2" = "#E69F00", "Variance" = "#56B4E9")
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.title = element_blank()
    ) +
    labs(
      x = "Method",
      y = "Component Value",
      title = title
    )

  p
}

#' Computation time comparison plot
#'
#' Creates a bar plot comparing computation times across methods
#' with optional log scale for wide-ranging times.
#'
#' @param df Aggregated results data frame with columns: method, mean_time
#'           (or custom column name specified)
#' @param time_col Character, name of time column (default "mean_time")
#' @param log_scale Logical, use log10 scale for y-axis (default TRUE)
#' @param show_error_bars Logical, show error bars if sd_time column exists (default TRUE)
#' @param title Optional plot title
#' @return ggplot object
#' @export
plot_computation_time <- function(df, time_col = "mean_time", log_scale = TRUE,
                                  show_error_bars = TRUE, title = NULL) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame")
  }
  if (nrow(df) == 0) {
    stop("'df' contains no rows")
  }
  if (!time_col %in% names(df)) {
    stop("Time column '", time_col, "' not found in data frame")
  }
  if (!"method" %in% names(df)) {
    stop("'df' must contain 'method' column")
  }

  if (is.null(title)) {
    title <- "Computation Time by Method"
  }

  # Aggregate by method if data has multiple scenarios
  if ("scenario_id" %in% names(df)) {
    agg_df <- aggregate(
      df[[time_col]],
      by = list(method = df$method),
      FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
    )
    # Handle edge case where aggregate returns vector instead of matrix for single method
    if (!is.matrix(agg_df$x)) {
      agg_df$x <- matrix(agg_df$x, nrow = 1, dimnames = list(NULL, c("mean", "sd")))
    }
    agg_df <- data.frame(
      method = agg_df$method,
      mean_time = agg_df$x[, "mean"],
      sd_time = agg_df$x[, "sd"]
    )
  } else {
    agg_df <- df
    if (!"sd_time" %in% names(agg_df)) {
      agg_df$sd_time <- NA_real_
    }
  }

  # Order methods by time
  agg_df <- agg_df[order(agg_df$mean_time), ]
  agg_df$method <- factor(agg_df$method, levels = agg_df$method)

  # Get colors
  methods_in_data <- as.character(unique(agg_df$method))
  colors_to_use <- METHOD_COLORS[intersect(methods_in_data, names(METHOD_COLORS))]
  missing_methods <- setdiff(methods_in_data, names(METHOD_COLORS))
  if (length(missing_methods) > 0) {
    extra_colors <- setNames(rep("grey50", length(missing_methods)), missing_methods)
    colors_to_use <- c(colors_to_use, extra_colors)
  }

  p <- ggplot(agg_df, aes(x = method, y = mean_time, fill = method)) +
    geom_bar(stat = "identity", alpha = 0.9) +
    scale_fill_manual(values = colors_to_use, guide = "none")

  # Add error bars if requested and data available
  if (show_error_bars && "sd_time" %in% names(agg_df) && !all(is.na(agg_df$sd_time))) {
    p <- p + geom_errorbar(
      aes(ymin = pmax(0.001, mean_time - sd_time),
          ymax = mean_time + sd_time),
      width = 0.3, linewidth = 0.5
    )
  }

  # Apply log scale if requested
  if (log_scale) {
    if (any(agg_df$mean_time <= 0, na.rm = TRUE)) {
      warning("Zero or negative time values detected; log scale may produce warnings")
    }
    p <- p + scale_y_log10(
      labels = scales::label_number(accuracy = 0.01)
    ) +
      labs(y = "Time (seconds, log scale)")
  } else {
    p <- p + labs(y = "Time (seconds)")
  }

  p <- p +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Method",
      title = title
    )

  p
}

#' Generate LaTeX table with kableExtra
#'
#' Converts a data frame to a formatted LaTeX table suitable for publication.
#' Uses booktabs style and provides options for styling.
#'
#' @param df Data frame to convert
#' @param caption Table caption
#' @param label LaTeX label for cross-referencing
#' @param digits Number of decimal places for numeric columns (default 3)
#' @param font_size Font size (NULL for default, or e.g., "small", "footnotesize")
#' @return Character string with LaTeX code
#' @export
generate_latex_table <- function(df, caption, label, digits = 3,
                                 font_size = NULL) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame")
  }
  if (nrow(df) == 0) {
    stop("'df' contains no rows")
  }
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("kableExtra package required for LaTeX table generation. ",
         "Install with: install.packages('kableExtra')")
  }

  # Round numeric columns
  numeric_cols <- sapply(df, is.numeric)
  df[numeric_cols] <- lapply(df[numeric_cols], function(x) round(x, digits))

  # Create base table
  tbl <- kableExtra::kbl(
    df,
    format = "latex",
    digits = digits,
    caption = caption,
    label = label,
    booktabs = TRUE,
    escape = TRUE
  )

  # Apply styling
  latex_opts <- c("striped", "hold_position")
  if (!is.null(font_size)) {
    tbl <- kableExtra::kable_styling(tbl, latex_options = latex_opts,
                                     font_size = NULL)
    tbl <- kableExtra::kable_styling(tbl, font_size = font_size)
  } else {
    tbl <- kableExtra::kable_styling(tbl, latex_options = latex_opts)
  }

  # Return as character
  as.character(tbl)
}

#' Create comprehensive benchmark report figure
#'
#' Combines multiple visualization plots into a single publication-ready figure
#' using patchwork for layout.
#'
#' @param agg_df Aggregated results data frame
#' @param raw_df Raw results data frame (optional, for bias boxplot)
#' @param friedman_result Output from friedman_test() (optional, for CD diagram)
#' @param nemenyi_result Output from nemenyi_posthoc() (optional, for CD diagram)
#' @param title Overall figure title
#' @return Combined ggplot/patchwork object
#' @export
create_benchmark_report_figure <- function(agg_df, raw_df = NULL,
                                           friedman_result = NULL,
                                           nemenyi_result = NULL,
                                           title = "Benchmark Results Summary") {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package required for combined figures. ",
         "Install with: install.packages('patchwork')")
  }

  plots <- list()

  # Panel A: RMSE heatmap
  if ("rmse_ate" %in% names(agg_df)) {
    plots$A <- plot_method_heatmap(agg_df, "rmse_ate",
                                   lower_is_better = TRUE,
                                   title = "A: RMSE by Method x Scenario")
  }

  # Panel B: Coverage heatmap
  if ("coverage" %in% names(agg_df)) {
    plots$B <- plot_coverage_heatmap(agg_df, "coverage",
                                     title = "B: Coverage (95% nominal)")
  }

  # Panel C: Bias boxplot (requires raw data)
  if (!is.null(raw_df) && "bias_ate" %in% names(raw_df)) {
    plots$C <- plot_bias_boxplot(raw_df, title = "C: Bias Distribution")
  }

  # Panel D: Bias-variance decomposition
  if (all(c("mean_bias", "sd_bias") %in% names(agg_df))) {
    plots$D <- plot_bias_variance(agg_df, title = "D: Bias-Variance Decomposition")
  }

  # Panel E: Computation time
  if ("mean_time" %in% names(agg_df)) {
    plots$E <- plot_computation_time(agg_df, title = "E: Computation Time")
  }

  # Panel F: CD diagram (requires statistical test results)
  if (!is.null(friedman_result) && !is.null(nemenyi_result)) {
    plots$F <- plot_cd_diagram(friedman_result, nemenyi_result,
                               title = "F: Critical Difference Diagram")
  }

  # Combine available plots
  if (length(plots) == 0) {
    stop("No valid plots could be created from the provided data")
  }

  combined <- patchwork::wrap_plots(plots) +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  combined
}

#' Save plot with publication-quality settings
#'
#' Wrapper around ggsave with sensible defaults for publication figures.
#'
#' @param plot ggplot object to save
#' @param filename Output filename (extension determines format)
#' @param width Width in inches (default 8)
#' @param height Height in inches (default 6)
#' @param dpi Resolution in dots per inch (default 300)
#' @param ... Additional arguments passed to ggsave
#' @return Invisibly returns the filename
#' @export
save_benchmark_plot <- function(plot, filename, width = 8, height = 6,
                                dpi = 300, ...) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
  invisible(filename)
}
