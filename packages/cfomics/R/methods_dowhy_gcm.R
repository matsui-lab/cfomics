#' @title DoWhy GCM Method Implementation
#' @description Internal functions for the DoWhy GCM causal inference method.

.get_dowhy_module <- function() {
  reticulate::import_from_path(
    module = "gcm_dowhy",
    path   = system.file("python", package = "cfomics")
  )
}

#' Fit DoWhy GCM Model
#'
#' Internal function to fit a causal model using DoWhy's Graphical Causal Models.
#'
#' @param X Numeric matrix of covariates.
#' @param T Binary treatment vector.
#' @param Y Numeric outcome vector.
#' @param dag_spec List of edge pairs defining the DAG.
#' @param random_state Integer random seed for reproducibility.
#' @param covariate_names Character vector of covariate names.
#' @param bootstrap Logical, whether to compute bootstrap confidence intervals.
#' @param n_bootstrap Integer, number of bootstrap iterations.
#' @param config_path Path to JSON configuration file.
#' @param return_metadata Logical, whether to return metadata.
#' @param ... Additional arguments passed to Python.
#' @return A cf_model_dowhy_gcm object.
#' @keywords internal
cf_fit_dowhy_gcm <- function(X, T, Y, dag_spec,
                             random_state = 0L,
                             covariate_names = NULL,
                             bootstrap = FALSE,
                             n_bootstrap = 100L,
                             config_path = NULL,
                             return_metadata = FALSE,
                             ...) {
  # Ensure Python environment is properly configured
  cf_require_python("dowhy_gcm")
  mod <- .get_dowhy_module()

  # R -> Python conversion
  # We pass covariate_names to fix the column name issue in Python DataFrame

  res <- mod$run_dowhy_gcm(
    X              = X,          # matrix -> numpy
    T              = as.numeric(T),
    Y              = as.numeric(Y),
    edges_list     = dag_spec,
    variable_names = covariate_names,
    random_state   = as.integer(random_state),
    config_path    = config_path,
    return_metadata = return_metadata,
    bootstrap      = bootstrap,
    n_bootstrap    = as.integer(n_bootstrap)
  )

  structure(
    list(
      backend = "dowhy_gcm",
      res     = res
    ),
    class = "cf_model_dowhy_gcm"
  )
}
