# benchmarks/figures/fig3_highdim.R
# Figure 3: High-dimensional stability analysis
#
# This figure demonstrates method stability as the dimensionality (p/n ratio)
# increases. It uses the S2 dimension sweep scenarios where n is fixed and
# p varies from 100 to 1000.
#
# The figure shows:
#   - Line plot of RMSE vs p/n ratio by method
#   - Highlights which methods remain stable in high dimensions
#   - Key for selecting appropriate methods for omics data (p >> n)

#' Generate Figure 3: High-dimensional stability
#'
#' @param results Data frame with benchmark results from load_benchmark_results()
#' @return A ggplot object showing RMSE vs p/n ratio
generate_fig3_highdim <- function(results) {
  # Validate input
  if (!is.data.frame(results)) {
    stop("'results' must be a data.frame")
  }

  required_cols <- c("scenario_id", "method", "status", "mse_ate", "n", "p")
  missing_cols <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Filter to S2 dimension sweep scenarios
  df_s2 <- results %>%
    filter(
      grepl("^S2", scenario_id),
      status == "ok"
    )

  if (nrow(df_s2) == 0) {
    # Fall back to any scenarios with varying p
    warning("S2 scenarios not found. Attempting to use scenarios with varying dimensions.")
    df_s2 <- results %>%
      filter(status == "ok") %>%
      group_by(method) %>%
      filter(n_distinct(p) > 1) %>%
      ungroup()

    if (nrow(df_s2) == 0) {
      stop("No scenarios with dimension variation found")
    }
  }

  # Calculate p/n ratio and aggregate
  df_agg <- df_s2 %>%
    mutate(p_n_ratio = p / n) %>%
    group_by(method, p_n_ratio, n, p) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      rmse_se = sqrt(var(mse_ate, na.rm = TRUE) / n()),
      mean_pehe = mean(pehe, na.rm = TRUE),
      n_reps = n(),
      .groups = "drop"
    )

  # Get method colors
  methods_in_data <- unique(df_agg$method)
  colors_to_use <- METHOD_COLORS[intersect(methods_in_data, names(METHOD_COLORS))]
  missing_methods <- setdiff(methods_in_data, names(METHOD_COLORS))
  if (length(missing_methods) > 0) {
    extra_colors <- setNames(
      scales::hue_pal()(length(missing_methods)),
      missing_methods
    )
    colors_to_use <- c(colors_to_use, extra_colors)
  }

  # Main plot: RMSE vs p/n ratio
  p <- ggplot(df_agg, aes(x = p_n_ratio, y = rmse_ate, color = method, group = method)) +
    # Add ribbon for standard error
    geom_ribbon(
      aes(ymin = pmax(0, rmse_ate - rmse_se),
          ymax = rmse_ate + rmse_se,
          fill = method),
      alpha = 0.15, color = NA
    ) +
    # Add lines and points
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    # Log scale for x-axis (p/n ratio spans orders of magnitude)
    scale_x_log10(
      breaks = c(0.1, 0.2, 0.5, 1.0, 2.0),
      labels = c("0.1", "0.2", "0.5", "1.0", "2.0")
    ) +
    # Color scales
    scale_color_manual(values = colors_to_use) +
    scale_fill_manual(values = colors_to_use) +
    # Annotations
    labs(
      x = expression(paste("Dimension ratio (", italic("p/n"), ")")),
      y = "RMSE (ATE)",
      color = "Method",
      fill = "Method"
    ) +
    # Add vertical line at p = n (ratio = 1)
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    annotate(
      "text", x = 1.05, y = max(df_agg$rmse_ate, na.rm = TRUE) * 0.95,
      label = "p = n", hjust = 0, size = 3, color = "grey40"
    ) +
    # Theme
    theme_paper() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor.x = element_blank()
    )

  p
}


#' Generate extended Figure 3 with multiple sample sizes
#'
#' Creates a faceted version showing high-dimensional behavior
#' across different sample sizes (n = 500, 1000).
#'
#' @param results Data frame with benchmark results
#' @return A patchwork object with faceted panels
generate_fig3_extended <- function(results) {
  # Filter to S2 scenarios
  df <- results %>%
    filter(
      grepl("^S2", scenario_id),
      status == "ok"
    ) %>%
    mutate(
      p_n_ratio = p / n,
      n_label = paste0("n = ", n)
    ) %>%
    group_by(method, n, p, p_n_ratio, n_label) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      rmse_se = sqrt(var(mse_ate, na.rm = TRUE) / n()),
      mean_pehe = mean(pehe, na.rm = TRUE),
      pehe_se = sd(pehe, na.rm = TRUE) / sqrt(n()),
      mean_time = mean(time_sec, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(df) == 0) {
    stop("No S2 scenarios found in results")
  }

  # Get colors
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

  # Panel A: RMSE across dimensions
  p_rmse <- ggplot(df, aes(x = p, y = rmse_ate, color = method, group = method)) +
    geom_ribbon(
      aes(ymin = pmax(0, rmse_ate - rmse_se),
          ymax = rmse_ate + rmse_se, fill = method),
      alpha = 0.1, color = NA
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~n_label, scales = "free_y") +
    scale_x_log10() +
    scale_color_manual(values = colors_to_use) +
    scale_fill_manual(values = colors_to_use) +
    labs(x = "Number of covariates (p)", y = "RMSE (ATE)") +
    theme_paper() +
    theme(legend.position = "bottom")

  # Panel B: PEHE across dimensions
  p_pehe <- ggplot(df, aes(x = p, y = mean_pehe, color = method, group = method)) +
    geom_ribbon(
      aes(ymin = pmax(0, mean_pehe - pehe_se),
          ymax = mean_pehe + pehe_se, fill = method),
      alpha = 0.1, color = NA
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~n_label, scales = "free_y") +
    scale_x_log10() +
    scale_color_manual(values = colors_to_use) +
    scale_fill_manual(values = colors_to_use) +
    labs(x = "Number of covariates (p)", y = "PEHE") +
    theme_paper() +
    theme(legend.position = "bottom")

  # Panel C: Computation time (log scale)
  p_time <- ggplot(df, aes(x = p, y = mean_time, color = method, group = method)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~n_label, scales = "free_y") +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = colors_to_use) +
    labs(x = "Number of covariates (p)", y = "Time (seconds, log scale)") +
    theme_paper() +
    theme(legend.position = "bottom")

  # Combine
  combined <- (p_rmse / p_pehe / p_time) +
    plot_layout(guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(size = 12, face = "bold"))
    ) &
    theme(legend.position = "bottom")

  combined
}


#' Generate stability heatmap
#'
#' Creates a heatmap showing relative performance degradation
#' as dimensions increase.
#'
#' @param results Data frame with benchmark results
#' @return A ggplot object with stability heatmap
generate_stability_heatmap <- function(results) {
  # Calculate relative degradation compared to lowest dimension
  df <- results %>%
    filter(grepl("^S2", scenario_id), status == "ok") %>%
    group_by(method, n, p) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    group_by(method, n) %>%
    mutate(
      # Reference: lowest p for this n
      rmse_ref = rmse_ate[which.min(p)],
      relative_change = (rmse_ate - rmse_ref) / rmse_ref * 100
    ) %>%
    ungroup()

  if (nrow(df) == 0) {
    stop("No S2 scenarios found")
  }

  ggplot(df, aes(x = factor(p), y = method, fill = relative_change)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%+.0f%%", relative_change)),
      size = 2.5
    ) +
    facet_wrap(~paste("n =", n), scales = "free_x") +
    scale_fill_gradient2(
      low = "#2166AC",   # Blue (improvement)
      mid = "white",
      high = "#B2182B",  # Red (degradation)
      midpoint = 0,
      name = "% change\nvs. baseline"
    ) +
    labs(
      x = "Number of covariates (p)",
      y = "Method",
      title = "Performance Stability Across Dimensions"
    ) +
    theme_paper() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0)
    )
}
