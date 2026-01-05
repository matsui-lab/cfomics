#' @title Visualization Functions for Causal Inference Results
#' @description R wrappers for Python visualization functions that create
#'   plots for ITE distributions, outcome shifts, uplift curves, and
#'   bootstrap ATE distributions.

.get_viz_module <- function() {
  reticulate::import_from_path(
    module = "visualization",
    path   = system.file("python", package = "cfomics")
  )
}

#' Plot ITE Distribution
#'
#' Creates a histogram of individual treatment effects with mean, median,
#' and 90% confidence interval annotations.
#'
#' @param ite Numeric vector of individual treatment effects.
#' @param save_path Optional path to save the plot as an image file.
#' @param figsize Numeric vector of length 2 specifying figure size (width, height).
#' @return Invisibly returns a list with matplotlib figure and axes objects.
#' @export
#' @examples
#' \dontrun{
#' fit <- cf_dowhy(Y ~ T | X1 + X2, data = mydata, graph = g)
#' ite <- predict(fit, type = "ite")
#' cf_plot_ite(ite, save_path = "ite_distribution.png")
#' }
cf_plot_ite <- function(ite, save_path = NULL, figsize = c(10, 6)) {
  ensure_python_env()
  viz <- .get_viz_module()
  
  result <- viz$plot_ite_distribution(
    ite = as.numeric(ite),
    save_path = save_path,
    figsize = reticulate::tuple(as.integer(figsize[1]), as.integer(figsize[2]))
  )
  
  invisible(result)
}

#' Plot Outcome Shift
#'
#' Creates a comparison plot showing the distribution of outcomes under
#' control (T=0) vs treatment (T=1), including both histograms and a
#' scatter plot of individual outcome shifts.
#'
#' @param y0 Numeric vector of counterfactual outcomes under T=0.
#' @param y1 Numeric vector of counterfactual outcomes under T=1.
#' @param save_path Optional path to save the plot as an image file.
#' @param figsize Numeric vector of length 2 specifying figure size (width, height).
#' @return Invisibly returns a list with matplotlib figure and axes objects.
#' @export
#' @examples
#' \dontrun{
#' fit <- cf_dowhy(Y ~ T | X1 + X2, data = mydata, graph = g)
#' y0 <- predict(fit, type = "y0")
#' y1 <- predict(fit, type = "y1")
#' cf_plot_outcome_shift(y0, y1, save_path = "outcome_shift.png")
#' }
cf_plot_outcome_shift <- function(y0, y1, save_path = NULL, figsize = c(14, 6)) {
  ensure_python_env()
  viz <- .get_viz_module()
  
  result <- viz$plot_outcome_shift(
    y0 = as.numeric(y0),
    y1 = as.numeric(y1),
    save_path = save_path,
    figsize = reticulate::tuple(as.integer(figsize[1]), as.integer(figsize[2]))
  )
  
  invisible(result)
}

#' Plot Uplift Curve
#'
#' Creates an uplift curve showing the cumulative treatment effect when
#' targeting individuals in order of their predicted treatment effect
#' (highest first). Compares against random targeting baseline.
#'
#' @param ite Numeric vector of individual treatment effects.
#' @param save_path Optional path to save the plot as an image file.
#' @param figsize Numeric vector of length 2 specifying figure size (width, height).
#' @return Invisibly returns a list with matplotlib figure and axes objects.
#' @export
#' @examples
#' \dontrun{
#' fit <- cf_dowhy(Y ~ T | X1 + X2, data = mydata, graph = g)
#' ite <- predict(fit, type = "ite")
#' cf_plot_uplift(ite, save_path = "uplift_curve.png")
#' }
cf_plot_uplift <- function(ite, save_path = NULL, figsize = c(10, 6)) {
  ensure_python_env()
  viz <- .get_viz_module()
  
  result <- viz$plot_uplift_curve(
    ite = as.numeric(ite),
    save_path = save_path,
    figsize = reticulate::tuple(as.integer(figsize[1]), as.integer(figsize[2]))
  )
  
  invisible(result)
}

#' Plot Bootstrap ATE Distribution
#'
#' Creates a histogram of bootstrap ATE samples with 95% confidence
#' interval and point estimate annotations.
#'
#' @param bootstrap_samples Numeric vector of bootstrap ATE samples.
#' @param ate_point Optional point estimate of ATE to highlight.
#' @param save_path Optional path to save the plot as an image file.
#' @param figsize Numeric vector of length 2 specifying figure size (width, height).
#' @return Invisibly returns a list with matplotlib figure and axes objects.
#' @export
#' @examples
#' \dontrun{
#' fit <- cf_dowhy(Y ~ T | X1 + X2, data = mydata, graph = g, bootstrap = TRUE)
#' samples <- fit$fit$res$samples$ate_bootstrap
#' ate <- predict(fit, type = "ate")
#' cf_plot_ate_bootstrap(samples, ate_point = ate, save_path = "ate_bootstrap.png")
#' }
cf_plot_ate_bootstrap <- function(bootstrap_samples, ate_point = NULL, 
                                   save_path = NULL, figsize = c(10, 6)) {
  ensure_python_env()
  viz <- .get_viz_module()
  
  result <- viz$plot_ate_bootstrap_distribution(
    bootstrap_samples = as.numeric(bootstrap_samples),
    ate_point = if (!is.null(ate_point)) as.numeric(ate_point) else NULL,
    save_path = save_path,
    figsize = reticulate::tuple(as.integer(figsize[1]), as.integer(figsize[2]))
  )
  
  invisible(result)
}

#' Plot All Causal Inference Results
#'
#' Convenience function to generate all standard plots from a cf_model object.
#'
#' @param object A cf_model object from cf_fit or cf_dowhy.
#' @param output_dir Directory to save plots. If NULL, plots are displayed but not saved.
#' @param prefix Prefix for output file names.
#' @return Invisibly returns a list of plot results.
#' @export
#' @examples
#' \dontrun{
#' fit <- cf_dowhy(Y ~ T | X1 + X2, data = mydata, graph = g, bootstrap = TRUE)
#' cf_plot_all(fit, output_dir = "plots/", prefix = "analysis_")
#' }
cf_plot_all <- function(object, output_dir = NULL, prefix = "") {
  if (!inherits(object, "cf_model")) {
    stop("object must be a cf_model from cf_fit() or cf_dowhy()")
  }
  
  results <- list()
  
  # Get predictions
  ite <- predict(object, type = "ite")
  y0 <- predict(object, type = "y0")
  y1 <- predict(object, type = "y1")
  ate <- predict(object, type = "ate")
  
  # Generate file paths if output_dir is provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    ite_path <- file.path(output_dir, paste0(prefix, "ite_distribution.png"))
    outcome_path <- file.path(output_dir, paste0(prefix, "outcome_shift.png"))
    uplift_path <- file.path(output_dir, paste0(prefix, "uplift_curve.png"))
  } else {
    ite_path <- NULL
    outcome_path <- NULL
    uplift_path <- NULL
  }
  
  # Generate plots
  results$ite <- cf_plot_ite(ite, save_path = ite_path)
  results$outcome_shift <- cf_plot_outcome_shift(y0, y1, save_path = outcome_path)
  results$uplift <- cf_plot_uplift(ite, save_path = uplift_path)
  
  # Bootstrap plot if available
  if (!is.null(object$fit$res$samples$ate_bootstrap)) {
    if (!is.null(output_dir)) {
      bootstrap_path <- file.path(output_dir, paste0(prefix, "ate_bootstrap.png"))
    } else {
      bootstrap_path <- NULL
    }
    results$bootstrap <- cf_plot_ate_bootstrap(
      object$fit$res$samples$ate_bootstrap,
      ate_point = ate,
      save_path = bootstrap_path
    )
  }
  
  invisible(results)
}
