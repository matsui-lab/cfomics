# benchmarks/figures/fig4_tcga.R
# Figure 4: TCGA semi-synthetic validation results
#
# This figure presents the real-world validation using TCGA gene expression
# data with semi-synthetic treatment effects. It demonstrates that cfomics
# methods perform well on realistic biological data structures.
#
# The figure shows:
#   - RMSE comparison across TCGA cancer types (BRCA, LUAD, COAD)
#   - Method performance on real covariate distributions
#   - Validation that simulation results transfer to real data

#' Generate Figure 4: TCGA validation results
#'
#' @param results Data frame with TCGA benchmark results from load_tcga_results()
#' @return A ggplot object comparing methods across TCGA projects
generate_fig4_tcga <- function(results) {
  # Validate input
  if (!is.data.frame(results)) {
    stop("'results' must be a data.frame")
  }

  required_cols <- c("project", "method", "status", "mse_ate")
  missing_cols <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Filter to successful results
  df_ok <- results %>%
    filter(status == "ok")

  if (nrow(df_ok) == 0) {
    stop("No successful TCGA results found")
  }

  # Aggregate by project and method
  df_agg <- df_ok %>%
    group_by(project, method) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      rmse_se = sqrt(var(mse_ate, na.rm = TRUE) / n()),
      mean_pehe = mean(pehe, na.rm = TRUE),
      pehe_se = sd(pehe, na.rm = TRUE) / sqrt(n()),
      coverage = mean(coverage_ate, na.rm = TRUE),
      mean_time = mean(time_sec, na.rm = TRUE),
      n_reps = n(),
      .groups = "drop"
    ) %>%
    mutate(
      # Simplify project names for display
      dataset = gsub("^TCGA-", "", project)
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

  # Main plot: RMSE by method and dataset
  p <- ggplot(df_agg, aes(x = method, y = rmse_ate, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, alpha = 0.85) +
    geom_errorbar(
      aes(ymin = pmax(0, rmse_ate - rmse_se),
          ymax = rmse_ate + rmse_se),
      position = position_dodge(width = 0.8),
      width = 0.25, linewidth = 0.4
    ) +
    # Dataset-based coloring (different from method colors)
    scale_fill_manual(
      values = c(
        "BRCA" = "#E41A1C",  # Red
        "LUAD" = "#377EB8",  # Blue
        "COAD" = "#4DAF4A"   # Green
      ),
      name = "Cancer Type"
    ) +
    labs(
      x = "Method",
      y = "RMSE (ATE)",
      title = NULL
    ) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  p
}


#' Generate extended Figure 4 with multiple panels
#'
#' Creates a comprehensive TCGA validation figure with RMSE, PEHE,
#' and coverage panels.
#'
#' @param results Data frame with TCGA benchmark results
#' @return A patchwork object with three panels
generate_fig4_extended <- function(results) {
  # Filter and aggregate
  df <- results %>%
    filter(status == "ok") %>%
    group_by(project, method) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      rmse_se = sqrt(var(mse_ate, na.rm = TRUE) / n()),
      mean_pehe = mean(pehe, na.rm = TRUE),
      pehe_se = sd(pehe, na.rm = TRUE) / sqrt(n()),
      coverage = mean(coverage_ate, na.rm = TRUE),
      mean_bias = mean(bias_ate, na.rm = TRUE),
      mean_time = mean(time_sec, na.rm = TRUE),
      n_reps = n(),
      .groups = "drop"
    ) %>%
    mutate(dataset = gsub("^TCGA-", "", project))

  if (nrow(df) == 0) {
    stop("No successful TCGA results found")
  }

  # Color palette for datasets
  dataset_colors <- c(
    "BRCA" = "#E41A1C",
    "LUAD" = "#377EB8",
    "COAD" = "#4DAF4A",
    "OV"   = "#984EA3",
    "KIRC" = "#FF7F00"
  )

  # Panel A: RMSE by method and dataset
  p_rmse <- ggplot(df, aes(x = method, y = rmse_ate, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = pmax(0, rmse_ate - rmse_se), ymax = rmse_ate + rmse_se),
      position = position_dodge(0.8), width = 0.2
    ) +
    scale_fill_manual(values = dataset_colors, name = "Dataset") +
    labs(x = NULL, y = "RMSE (ATE)") +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Panel B: PEHE by method and dataset
  p_pehe <- ggplot(df, aes(x = method, y = mean_pehe, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = pmax(0, mean_pehe - pehe_se), ymax = mean_pehe + pehe_se),
      position = position_dodge(0.8), width = 0.2
    ) +
    scale_fill_manual(values = dataset_colors, name = "Dataset") +
    labs(x = NULL, y = "PEHE") +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Panel C: Coverage (should be ~0.95)
  p_cov <- ggplot(df, aes(x = method, y = coverage, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.5) +
    scale_fill_manual(values = dataset_colors, name = "Dataset") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Coverage (95% nominal)") +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Combine with shared legend
  combined <- (p_rmse + p_pehe + p_cov) +
    plot_layout(guides = "collect", ncol = 3) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(size = 12, face = "bold"))
    ) &
    theme(legend.position = "bottom")

  combined
}


#' Generate TCGA method comparison heatmap
#'
#' Creates a heatmap showing relative performance across cancer types.
#'
#' @param results Data frame with TCGA benchmark results
#' @param metric Which metric to display ("rmse" or "pehe")
#' @return A ggplot heatmap
generate_tcga_heatmap <- function(results, metric = c("rmse", "pehe")) {
  metric <- match.arg(metric)

  df <- results %>%
    filter(status == "ok") %>%
    group_by(project, method) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      mean_pehe = mean(pehe, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(dataset = gsub("^TCGA-", "", project))

  if (nrow(df) == 0) {
    stop("No successful TCGA results found")
  }

  # Select metric
  if (metric == "rmse") {
    df$value <- df$rmse_ate
    fill_label <- "RMSE (ATE)"
  } else {
    df$value <- df$mean_pehe
    fill_label <- "PEHE"
  }

  # Calculate midpoint for color scale
  midpoint <- median(df$value, na.rm = TRUE)

  ggplot(df, aes(x = dataset, y = method, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%.3f", value)),
      size = 3, color = "black"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",   # Blue (good)
      mid = "white",
      high = "#B2182B",  # Red (bad)
      midpoint = midpoint,
      name = fill_label
    ) +
    labs(
      x = "Cancer Type",
      y = "Method",
      title = paste("TCGA Validation:", fill_label)
    ) +
    theme_paper() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0)
    )
}


#' Generate TCGA ranking summary
#'
#' Creates a summary plot showing method rankings across cancer types.
#'
#' @param results Data frame with TCGA benchmark results
#' @return A ggplot object with mean ranks
generate_tcga_ranking <- function(results) {
  # Calculate ranks per dataset
  df <- results %>%
    filter(status == "ok") %>%
    group_by(project, method) %>%
    summarise(
      rmse_ate = sqrt(mean(mse_ate, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    group_by(project) %>%
    mutate(rank = rank(rmse_ate, ties.method = "average")) %>%
    ungroup() %>%
    mutate(dataset = gsub("^TCGA-", "", project))

  if (nrow(df) == 0) {
    stop("No successful TCGA results found")
  }

  # Calculate mean rank across datasets
  mean_ranks <- df %>%
    group_by(method) %>%
    summarise(
      mean_rank = mean(rank),
      sd_rank = sd(rank),
      n_datasets = n(),
      .groups = "drop"
    ) %>%
    arrange(mean_rank)

  # Get colors
  methods_in_data <- mean_ranks$method
  colors_to_use <- METHOD_COLORS[intersect(methods_in_data, names(METHOD_COLORS))]
  missing_methods <- setdiff(methods_in_data, names(METHOD_COLORS))
  if (length(missing_methods) > 0) {
    extra_colors <- setNames(
      scales::hue_pal()(length(missing_methods)),
      missing_methods
    )
    colors_to_use <- c(colors_to_use, extra_colors)
  }

  # Order methods by mean rank
  mean_ranks$method <- factor(mean_ranks$method, levels = mean_ranks$method)

  ggplot(mean_ranks, aes(x = method, y = mean_rank, fill = method)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_errorbar(
      aes(ymin = pmax(1, mean_rank - sd_rank),
          ymax = mean_rank + sd_rank),
      width = 0.25
    ) +
    scale_fill_manual(values = colors_to_use, guide = "none") +
    scale_y_reverse() +  # Lower rank = better
    labs(
      x = "Method",
      y = "Mean Rank (lower is better)",
      title = "Method Rankings Across TCGA Datasets"
    ) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
