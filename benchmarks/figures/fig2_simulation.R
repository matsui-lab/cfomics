# benchmarks/figures/fig2_simulation.R
# Figure 2: Main simulation benchmark results
#
# This figure presents the core benchmark results comparing cfomics methods
# across key simulation scenarios. It consists of two panels:
#   A: RMSE (ATE) comparison across scenarios
#   B: PEHE (heterogeneous effects) comparison across scenarios
#
# Key scenarios highlighted:
#   - S1_n1000_p100: Baseline linear confounding
#   - S5_combined: Nonlinear confounding
#   - S7_weak: Weak propensity score overlap
#   - S12_nonlinear_outcome_moderate: Nonlinear outcome model

#' Generate Figure 2: Main simulation benchmark results
#'
#' @param results Data frame with benchmark results from load_benchmark_results()
#' @return A patchwork object with two panels (RMSE and PEHE)
generate_fig2_simulation <- function(results) {
  # Validate input
  if (!is.data.frame(results)) {
    stop("'results' must be a data.frame")
  }

  required_cols <- c("scenario_id", "method", "status", "squared_error_ate", "pehe")
  missing_cols <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Key scenarios to highlight in the main figure
  # These represent diverse challenges in causal inference
  key_scenarios <- c(
    "S1_n1000_p100",                    # Baseline: linear confounding, moderate dimension
    "S5_combined",                       # Nonlinear confounding (combined effects)
    "S7_weak",                           # Weak propensity score overlap
    "S12_nonlinear_outcome_moderate"     # Nonlinear outcome model
  )

  # Filter to successful results in key scenarios
  df_filtered <- results %>%
    filter(
      scenario_id %in% key_scenarios,
      status == "ok"
    )

  if (nrow(df_filtered) == 0) {
    # Fall back to any available scenarios
    warning("Key scenarios not found. Using all available scenarios.")
    available_scenarios <- unique(results$scenario_id[results$status == "ok"])
    if (length(available_scenarios) == 0) {
      stop("No successful results found")
    }
    key_scenarios <- head(available_scenarios, 4)
    df_filtered <- results %>%
      filter(scenario_id %in% key_scenarios, status == "ok")
  }

  # Aggregate by scenario and method
  df_agg <- df_filtered %>%
    group_by(scenario_id, method) %>%
    summarise(
      rmse_ate = sqrt(mean(squared_error_ate, na.rm = TRUE)),
      rmse_ate_se = sqrt(var(squared_error_ate, na.rm = TRUE) / n()),
      mean_pehe = mean(pehe, na.rm = TRUE),
      pehe_se = sd(pehe, na.rm = TRUE) / sqrt(n()),
      coverage = mean(coverage_ate, na.rm = TRUE),
      n_reps = n(),
      .groups = "drop"
    ) %>%
    mutate(
      # Create human-readable scenario labels
      scenario_label = case_when(
        grepl("^S1", scenario_id) ~ "Baseline\n(linear confounding)",
        grepl("S5_combined", scenario_id) ~ "Nonlinear\nconfounding",
        grepl("S7_weak", scenario_id) ~ "Weak\noverlap",
        grepl("S12", scenario_id) ~ "Nonlinear\noutcome",
        TRUE ~ scenario_id
      )
    )

  # Order scenarios logically
  scenario_order <- c(
    "Baseline\n(linear confounding)",
    "Nonlinear\nconfounding",
    "Weak\noverlap",
    "Nonlinear\noutcome"
  )
  df_agg$scenario_label <- factor(
    df_agg$scenario_label,
    levels = scenario_order
  )

  # Get method colors (handle missing methods gracefully)
  methods_in_data <- unique(df_agg$method)
  colors_to_use <- METHOD_COLORS[intersect(methods_in_data, names(METHOD_COLORS))]

  # Generate colors for methods not in our palette
  missing_methods <- setdiff(methods_in_data, names(METHOD_COLORS))
  if (length(missing_methods) > 0) {
    extra_colors <- setNames(
      scales::hue_pal()(length(missing_methods)),
      missing_methods
    )
    colors_to_use <- c(colors_to_use, extra_colors)
  }

  # Panel A: RMSE (ATE) by scenario and method
  p_rmse <- ggplot(df_agg, aes(x = method, y = rmse_ate, fill = method)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
    geom_errorbar(
      aes(ymin = pmax(0, rmse_ate - rmse_ate_se),
          ymax = rmse_ate + rmse_ate_se),
      width = 0.25, linewidth = 0.4
    ) +
    facet_wrap(~scenario_label, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = colors_to_use) +
    labs(
      x = NULL,
      y = "RMSE (ATE)",
      title = NULL
    ) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "none",
      strip.text = element_text(size = 9)
    )

  # Panel B: PEHE by scenario and method
  p_pehe <- ggplot(df_agg, aes(x = method, y = mean_pehe, fill = method)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
    geom_errorbar(
      aes(ymin = pmax(0, mean_pehe - pehe_se),
          ymax = mean_pehe + pehe_se),
      width = 0.25, linewidth = 0.4
    ) +
    facet_wrap(~scenario_label, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = colors_to_use) +
    labs(
      x = NULL,
      y = "PEHE",
      title = NULL
    ) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "none",
      strip.text = element_text(size = 9)
    )

  # Combine panels with patchwork
  combined <- p_rmse / p_pehe +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(size = 12, face = "bold")
      )
    )

  combined
}


#' Generate extended Figure 2 with additional scenarios
#'
#' Creates a more comprehensive figure including heterogeneity and
#' high-dimensional scenarios.
#'
#' @param results Data frame with benchmark results
#' @return A patchwork object with three panels
generate_fig2_extended <- function(results) {
  # Extended scenario set
  extended_scenarios <- c(
    # Baseline variations
    "S1_n500_p50", "S1_n1000_p100", "S1_n2000_p100",
    # Nonlinear confounding
    "S5_quadratic", "S5_combined", "S5_threshold",
    # Overlap issues
    "S7_good", "S7_moderate", "S7_weak",
    # Heterogeneity
    "S3_str0", "S3_str1", "S3_str2",
    # High-dimensional
    "S2_n1000_p100", "S2_n1000_p500", "S2_n1000_p1000"
  )

  # Filter and aggregate
  df <- results %>%
    filter(status == "ok") %>%
    filter(scenario_id %in% extended_scenarios | any(scenario_id %in% extended_scenarios)) %>%
    group_by(scenario_id, method) %>%
    summarise(
      rmse_ate = sqrt(mean(squared_error_ate, na.rm = TRUE)),
      mean_pehe = mean(pehe, na.rm = TRUE),
      coverage = mean(coverage_ate, na.rm = TRUE),
      mean_time = mean(time_sec, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      scenario_group = case_when(
        grepl("^S1", scenario_id) ~ "Baseline",
        grepl("^S2", scenario_id) ~ "High-dimensional",
        grepl("^S3", scenario_id) ~ "Heterogeneity",
        grepl("^S5", scenario_id) ~ "Nonlinear confounding",
        grepl("^S7", scenario_id) ~ "Overlap",
        TRUE ~ "Other"
      )
    )

  if (nrow(df) == 0) {
    stop("No matching scenarios found in results")
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

  # Panel A: RMSE by scenario group
  p_rmse <- ggplot(df, aes(x = method, y = rmse_ate, fill = method)) +
    geom_boxplot(alpha = 0.8, outlier.size = 1) +
    facet_wrap(~scenario_group, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = colors_to_use) +
    labs(x = NULL, y = "RMSE (ATE)") +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Panel B: PEHE by scenario group
  p_pehe <- ggplot(df, aes(x = method, y = mean_pehe, fill = method)) +
    geom_boxplot(alpha = 0.8, outlier.size = 1) +
    facet_wrap(~scenario_group, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = colors_to_use) +
    labs(x = NULL, y = "PEHE") +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Panel C: Coverage (should be ~0.95)
  p_cov <- ggplot(df, aes(x = method, y = coverage, fill = method)) +
    geom_boxplot(alpha = 0.8, outlier.size = 1) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.5) +
    facet_wrap(~scenario_group, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = colors_to_use) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Coverage (95% nominal)") +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Combine
  combined <- p_rmse / p_pehe / p_cov +
    plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(size = 12, face = "bold"))
    )

  combined
}
