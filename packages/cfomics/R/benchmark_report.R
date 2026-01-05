#' Plot benchmark results
#'
#' This function creates visualization plots for benchmark results.
#' Requires ggplot2 to be installed; if not available, a message is printed
#' and the function returns invisibly without error.
#'
#' @param summary_df data.frame from cf_benchmark_summarize()
#' @param out_dir Character, optional directory path to save PDF plots.
#'   If NULL, plots are displayed but not saved.
#'
#' @return Invisibly returns a list of ggplot objects (if ggplot2 is available),
#'   or NULL if ggplot2 is not available.
#'
#' @export
cf_benchmark_plot <- function(summary_df, out_dir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not available; skipping plots")
    return(invisible(NULL))
  }
  
  plots <- list()
  
  plot_data <- summary_df[!is.na(summary_df$rmse_ate), ]
  
  if (nrow(plot_data) == 0) {
    message("No valid data for plotting (all RMSE values are NA)")
    return(invisible(NULL))
  }
  
  p_rmse <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$method, y = .data$rmse_ate, fill = .data$method)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::facet_wrap(~ scenario_id, scales = "free_y") +
    ggplot2::labs(
      title = "RMSE of ATE by Method and Scenario",
      x = "Method",
      y = "RMSE(ATE)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  plots$rmse_ate <- p_rmse
  
  pehe_data <- plot_data[!is.na(plot_data$mean_pehe), ]
  
  if (nrow(pehe_data) > 0) {
    p_pehe <- ggplot2::ggplot(
      pehe_data,
      ggplot2::aes(x = .data$method, y = .data$mean_pehe, fill = .data$method)
    ) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::facet_wrap(~ scenario_id, scales = "free_y") +
      ggplot2::labs(
        title = "Mean PEHE by Method and Scenario",
        x = "Method",
        y = "Mean PEHE"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    plots$pehe <- p_pehe
  }
  
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    ggplot2::ggsave(
      filename = file.path(out_dir, "rmse_ate.pdf"),
      plot = plots$rmse_ate,
      width = 10,
      height = 6
    )
    message(sprintf("Saved: %s", file.path(out_dir, "rmse_ate.pdf")))
    
    if (!is.null(plots$pehe)) {
      ggplot2::ggsave(
        filename = file.path(out_dir, "pehe.pdf"),
        plot = plots$pehe,
        width = 10,
        height = 6
      )
      message(sprintf("Saved: %s", file.path(out_dir, "pehe.pdf")))
    }
  }
  
  invisible(plots)
}
