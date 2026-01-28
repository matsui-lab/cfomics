#' Fit a Causal Inference Model
#'
#' @description
#' Fits a causal inference model using various backends with a unified API.
#' Supports data.frame, SummarizedExperiment, and MultiAssayExperiment objects.
#'
#' @details
#' \strong{Estimand:}
#'
#' All methods estimate the \strong{Average Treatment Effect (ATE)}:
#' \deqn{ATE = E[Y(1) - Y(0)]}
#'
#' where \eqn{Y(1)} and \eqn{Y(0)} are potential outcomes under treatment and
#' control, respectively. Some methods (grf, drlearner, ganite, cavae) also
#' estimate \strong{Individual Treatment Effects (ITE)}:
#' \deqn{ITE_i = Y_i(1) - Y_i(0)}
#'
#' \strong{Causal Assumptions:}
#'
#' Valid causal inference requires:
#' \enumerate{
#'   \item \strong{SUTVA}: No interference between units
#'   \item \strong{Ignorability}: \eqn{(Y(0), Y(1)) \perp T | X} (no unmeasured confounding)
#'   \item \strong{Positivity}: \eqn{0 < P(T=1|X) < 1} for all X
#' }
#'
#' Use \code{\link{cf_balance_check}} and \code{\link{cf_overlap_check}} to
#' assess these assumptions before analysis.
#'
#' \strong{Method Selection Guide:}
#' \itemize{
#'   \item \strong{grf}: Best for heterogeneous effects, robust, R-native
#'   \item \strong{ipw}: Classical approach, easy to interpret, requires correct PS model
#'   \item \strong{gformula}: Outcome modeling, straightforward implementation
#'   \item \strong{dowhy_gcm}: DAG-based, explicit causal assumptions
#'   \item \strong{drlearner}: Double robustness, requires Python
#'   \item \strong{ganite}: Deep learning ITE, requires TensorFlow
#'   \item \strong{cavae}: Handles latent confounding, requires Pyro
#' }
#'
#' @param formula Formula specifying the model: \code{Y ~ T | X1 + X2 + ...}
#'   \itemize{
#'     \item Y: Outcome variable (numeric)
#'     \item T: Treatment variable (binary: 0/1)
#'     \item X1, X2, ...: Covariates (numeric)
#'   }
#'   For SummarizedExperiment, Y can be a gene name from the assay or a
#'   column from colData.
#' @param data Data source: data.frame, SummarizedExperiment, or
#'   MultiAssayExperiment object.
#' @param method Character, the causal inference method to use:
#'   \describe{
#'     \item{"grf"}{Generalized Random Forest (R-native, recommended for most cases)}
#'     \item{"ipw"}{Inverse Probability Weighting (R-native)}
#'     \item{"gformula"}{G-computation/parametric g-formula (R-native, no extra packages)}
#'     \item{"hdml"}{High-Dimensional Machine Learning/Debiased Lasso (R-native, requires glmnet)}
#'     \item{"hdps"}{High-Dimensional Propensity Score/IPW (R-native, requires glmnet)}
#'     \item{"dowhy_gcm"}{DoWhy Graphical Causal Model (Python, requires DAG)}
#'     \item{"drlearner"}{Double/Debiased Machine Learning (Python)}
#'     \item{"ganite"}{GAN-based ITE estimation (Python/TensorFlow)}
#'     \item{"cavae"}{Causal Effect VAE (Python/Pyro)}
#'     \item{"auto"}{Automatically select best available method}
#'   }
#' @param graph igraph object representing the causal DAG. Required for
#'   "dowhy_gcm" method. Nodes should match variable names in formula.
#' @param bootstrap Logical, whether to compute bootstrap confidence intervals
#'   for ATE (default FALSE). Currently only for "dowhy_gcm" method.
#' @param n_bootstrap Integer, number of bootstrap iterations (default 100).
#' @param return_metadata Logical, whether to include extended metadata in
#'   results (default FALSE).
#' @param random_state Integer, random seed for reproducibility (default 0).
#' @param config_path Optional path to JSON configuration file for advanced
#'   settings (e.g., dimensionality reduction, hyperparameters).
#' @param assay_name For SummarizedExperiment: name or index of assay to use
#'   (default 1).
#' @param feature_select For SummarizedExperiment/MultiAssayExperiment:
#'   feature selection method ("none", "variance", "pca").
#' @param n_features Number of features to select from assay (default 100).
#' @param ... Additional method-specific options passed to the backend.
#'
#' @return Object of class \code{c("cf_model", "cfomics_result")} with components:
#'   \describe{
#'     \item{method}{Character string identifying the method used.}
#'     \item{fit}{Method-specific fitted model with \code{res} containing results.}
#'     \item{meta}{Metadata including formula, dimensions, and software versions.}
#'   }
#'   See \code{\link{cfomics_result}} for the full structure specification.
#'
#' @seealso
#'   \code{\link{cfomics_result}} for result object structure,
#'   \code{\link{cf_balance_check}} and \code{\link{cf_overlap_check}} for assumption checks,
#'   \code{\link{cf_check_env}} for environment diagnostics
#' @export
#' @examples
#' # Basic example with G-formula (no external packages required)
#' set.seed(123)
#' n <- 100
#' demo_data <- data.frame(
#'   X1 = rnorm(n),
#'   X2 = rnorm(n),
#'   T = rbinom(n, 1, 0.5)
#' )
#' demo_data$Y <- 2 * demo_data$T + 0.5 * demo_data$X1 + rnorm(n)
#'
#' # Fit model using G-formula
#' fit <- cf_fit(Y ~ T | X1 + X2, data = demo_data, method = "gformula")
#'
#' # Get average treatment effect
#' ate <- predict(fit, type = "ate")
#' print(ate)
#'
#' # Get individual treatment effects
#' ite <- predict(fit, type = "ite")
#' head(ite)
#'
#' # Model summary
#' summary(fit)
#'
#' \donttest{
#' # Using GRF (requires 'grf' package)
#' if (requireNamespace("grf", quietly = TRUE)) {
#'   fit_grf <- cf_fit(Y ~ T | X1 + X2, data = demo_data, method = "grf")
#'   predict(fit_grf, type = "ate")
#' }
#' }
#'
#' \dontrun{
#' # Using DoWhy-GCM (requires Python and dowhy package)
#' library(igraph)
#' g <- graph_from_edgelist(rbind(
#'   c("X1", "T"), c("T", "Y"), c("X2", "Y")
#' ))
#' fit_dowhy <- cf_fit(Y ~ T | X1 + X2, data = demo_data,
#'                     method = "dowhy_gcm", graph = g)
#'
#' # With SummarizedExperiment
#' library(SummarizedExperiment)
#' fit_se <- cf_fit(Gene1 ~ treatment | age, data = se, method = "grf")
#' }
cf_fit <- function(formula, data,
                   method = c("grf", "ipw", "gformula", "hdml", "hdps",
                              "dowhy_gcm", "ganite", "drlearner", "cavae", "auto"),
                   graph = NULL,
                   bootstrap = FALSE, n_bootstrap = 100L,
                   return_metadata = FALSE, random_state = 0L,
                   config_path = NULL,
                   assay_name = 1,
                   feature_select = "none",
                   n_features = 100,
                   ...) {
  method <- match.arg(method)

  # Handle "auto" method selection
  if (method == "auto") {
    method <- select_best_method(data)
    message(sprintf("Auto-selected method: %s", method))
  }

  # Parse data based on input type
  if (is_se(data)) {
    parsed <- as_cf_data.SummarizedExperiment(
      data, formula = formula,
      assay_name = assay_name,
      feature_select = feature_select,
      n_features = n_features
    )
  } else if (is_mae(data)) {
    parsed <- as_cf_data.MultiAssayExperiment(
      data, formula = formula,
      feature_select = feature_select,
      n_features_per_assay = n_features,
      ...
    )
  } else {
    # Standard data.frame path
    parsed <- parse_cf_formula(formula, data)
  }
  Y <- parsed$Y
  T <- parsed$T
  X <- parsed$X

  dag_spec <- if (!is.null(graph)) dag_to_edges_list(graph) else NULL

  fit <- switch(
    method,
    "dowhy_gcm" = {
      if (is.null(graph)) {
        stop("graph argument is required for method 'dowhy_gcm'")
      }
      cf_fit_dowhy_gcm(X = X, T = T, Y = Y,
                       dag_spec = dag_spec,
                       covariate_names = parsed$covariate_names,
                       bootstrap = bootstrap,
                       n_bootstrap = as.integer(n_bootstrap),
                       return_metadata = return_metadata,
                       random_state = as.integer(random_state),
                       config_path = config_path,
                       ...)
    },
    "ganite" = cf_fit_ganite(X = X, T = T, Y = Y,
                             dag_spec = dag_spec,
                             covariate_names = parsed$covariate_names,
                             random_state = as.integer(random_state),
                             ...),
    "drlearner" = cf_fit_drlearner(X = X, T = T, Y = Y,
                                   covariate_names = parsed$covariate_names,
                                   random_state = as.integer(random_state),
                                   ...),
    "cavae" = cf_fit_cavae(X = X, T = T, Y = Y,
                           covariate_names = parsed$covariate_names,
                           random_state = as.integer(random_state),
                           ...),
    "grf" = cf_fit_grf(X = X, T = T, Y = Y,
                       covariate_names = parsed$covariate_names,
                       random_state = as.integer(random_state),
                       ...),
    "ipw" = cf_fit_ipw(X = X, T = T, Y = Y,
                       covariate_names = parsed$covariate_names,
                       ...),
    "gformula" = cf_fit_gformula(X = X, T = T, Y = Y,
                                 covariate_names = parsed$covariate_names,
                                 ...),
    "hdml" = cf_fit_hdml(X = X, T = T, Y = Y,
                         covariate_names = parsed$covariate_names,
                         ...),
    "hdps" = cf_fit_hdps(X = X, T = T, Y = Y,
                         covariate_names = parsed$covariate_names,
                         ...),
    stop(sprintf("Unknown method: %s", method))
  )

  # Create standardized meta object
  meta <- create_cf_meta(
    formula = formula,
    n = length(Y),
    p = ncol(X),
    outcome_name = parsed$outcome_name,
    treatment_name = parsed$treatment_name,
    covariate_names = parsed$covariate_names,
    method = method,
    estimand = "ATE",  # All methods estimate ATE
    target_population = "all",
    dag_spec = dag_spec,
    bootstrap = bootstrap,
    random_state = random_state,
    preprocess = parsed$preprocess  # From feature selection if applicable
  )

  result <- structure(
    list(
      method = method,
      fit    = fit,  # method specific object
      meta   = meta
    ),
    class = c("cf_model", "cfomics_result")
  )

  # Validate result structure (contract enforcement)
  tryCatch({
    validate_cfomics_result(result)
  }, error = function(e) {
    rlang::warn(
      message = sprintf("Result validation warning for method '%s': %s", method, e$message),
      class = "cfomics_validation_warning"
    )
  })

  result
}

#' Convenience wrapper for DoWhy-GCM
#'
#' @param formula Formula Y ~ T | X...
#' @param data data.frame
#' @param graph igraph DAG
#' @param bootstrap Logical, whether to compute bootstrap confidence intervals.
#' @param n_bootstrap Integer, number of bootstrap iterations.
#' @param return_metadata Logical, whether to include metadata.
#' @param random_state Integer, random seed.
#' @param config_path Optional path to JSON configuration file.
#' @param ... Additional arguments passed to cf_fit
#' @export
cf_dowhy <- function(formula, data, graph, bootstrap = FALSE,
                     n_bootstrap = 100L, return_metadata = FALSE,
                     random_state = 0L, config_path = NULL, ...) {
  cf_fit(formula = formula, data = data,
         method = "dowhy_gcm", graph = graph,
         bootstrap = bootstrap, n_bootstrap = n_bootstrap,
         return_metadata = return_metadata, random_state = random_state,
         config_path = config_path, ...)
}

#' Select Best Available Method
#'
#' Automatically selects the best available causal inference method
#' based on installed packages and Python availability.
#'
#' @param data Data object (used for future data-dependent selection)
#' @return Character string with method name
#' @noRd
select_best_method <- function(data = NULL) {
  # Prefer R-only methods for reliability
  if (requireNamespace("grf", quietly = TRUE)) {
    return("grf")
  }

  # HDML is good for high-dimensional data with glmnet
  if (requireNamespace("glmnet", quietly = TRUE)) {
    return("hdml")
  }

  # Check Python availability
  py_available <- tryCatch({
    reticulate::py_available(initialize = FALSE)
  }, error = function(e) FALSE)

  if (py_available) {
    # Check for Python packages
    if (tryCatch(reticulate::py_module_available("dowhy"), error = function(e) FALSE)) {
      return("dowhy_gcm")
    }
    if (tryCatch(reticulate::py_module_available("econml"), error = function(e) FALSE)) {
      return("drlearner")
    }
  }

  # R-only fallbacks
  if (requireNamespace("ipw", quietly = TRUE)) {
    return("ipw")
  }

  # G-formula is always available (built-in)
  return("gformula")
}
