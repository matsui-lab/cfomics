#' @title Result Contract: cfomics_result Object Structure
#'
#' @description
#' This file documents the formal structure (contract) for \code{cfomics_result}
#' objects returned by \code{\link{cf_fit}}. All causal inference methods in
#' cfomics MUST return objects conforming to this structure.
#'
#' @section Object Structure:
#'
#' A \code{cfomics_result} object is a list with the following components:
#'
#' \describe{
#'   \item{method}{Character string identifying the method used (e.g., "grf",
#'     "dowhy_gcm").}
#'   \item{fit}{Method-specific fitted model object. Structure varies by method.}
#'   \item{meta}{Metadata list (see Meta Structure below).}
#' }
#'
#' @section Meta Structure:
#'
#' The \code{meta} component is a list containing:
#'
#' \describe{
#'   \item{formula}{The formula used for fitting.}
#'   \item{n}{Integer, number of observations.}
#'   \item{p}{Integer, number of covariates.}
#'   \item{outcome_name}{Character, name of the outcome variable.}
#'   \item{treatment_name}{Character, name of the treatment variable.}
#'   \item{covariate_names}{Character vector of covariate names.}
#'   \item{estimand}{Character, the causal estimand being estimated (see Estimands).}
#'   \item{target_population}{Character, the target population ("all", "treated",
#'     or "control").}
#'   \item{dag_spec}{List of edge pairs if a DAG was provided, NULL otherwise.}
#'   \item{bootstrap}{Logical, whether bootstrap was used.}
#'   \item{random_state}{Integer, the random seed used.}
#'   \item{preprocess}{List describing preprocessing steps applied.}
#'   \item{software_versions}{List of R and Python package versions used.}
#' }
#'
#' @section Fit Structure (res component):
#'
#' The \code{fit$res} component contains the actual results and MUST include:
#'
#' \describe{
#'   \item{ite}{Numeric vector of Individual Treatment Effects. Length equals n.}
#'   \item{ate}{Numeric scalar, Average Treatment Effect.}
#'   \item{y0_hat}{Numeric vector of counterfactual outcomes under no treatment.}
#'   \item{y1_hat}{Numeric vector of counterfactual outcomes under treatment.}
#'   \item{summary}{List containing summary statistics (see Summary Structure).}
#' }
#'
#' @section Summary Structure:
#'
#' The \code{fit$res$summary} list MUST contain:
#'
#' \describe{
#'   \item{ate}{Numeric, point estimate of ATE.}
#'   \item{ate_ci_lower}{Numeric, lower bound of 95\% CI for ATE.}
#'   \item{ate_ci_upper}{Numeric, upper bound of 95\% CI for ATE.}
#'   \item{ite_mean}{Numeric, mean of ITEs.}
#'   \item{ite_std}{Numeric, standard deviation of ITEs.}
#'   \item{ite_quantiles}{List with q05, q50, q95 quantiles of ITEs.}
#' }
#'
#' @section Estimands:
#'
#' cfomics supports the following causal estimands:
#'
#' \describe{
#'   \item{ATE (Average Treatment Effect)}{
#'     \eqn{E[Y(1) - Y(0)]}
#'     The expected difference in outcomes between treatment and control
#'     across the entire population.
#'   }
#'   \item{ATT (Average Treatment Effect on the Treated)}{
#'     \eqn{E[Y(1) - Y(0) | T=1]}
#'     The expected effect for those who actually received treatment.
#'   }
#'   \item{ATC (Average Treatment Effect on the Control)}{
#'     \eqn{E[Y(1) - Y(0) | T=0]}
#'     The expected effect for those who did not receive treatment.
#'   }
#'   \item{ITE (Individual Treatment Effect)}{
#'     \eqn{Y_i(1) - Y_i(0)}
#'     The treatment effect for individual i.
#'   }
#'   \item{CATE (Conditional Average Treatment Effect)}{
#'     \eqn{E[Y(1) - Y(0) | X=x]}
#'     The expected effect conditional on covariates.
#'   }
#' }
#'
#' @section Causal Assumptions:
#'
#' All methods in cfomics rely on the following fundamental assumptions:
#'
#' \describe{
#'   \item{SUTVA (Stable Unit Treatment Value Assumption)}{
#'     No interference between units and no hidden variations of treatment.
#'   }
#'   \item{Ignorability (Conditional Independence)}{
#'     \eqn{(Y(0), Y(1)) \perp T | X}
#'     Treatment assignment is independent of potential outcomes given covariates.
#'   }
#'   \item{Positivity (Overlap)}{
#'     \eqn{0 < P(T=1|X) < 1}
#'     Every unit has a positive probability of receiving either treatment.
#'   }
#' }
#'
#' Use \code{\link{cf_balance_check}} and \code{\link{cf_overlap_check}} to
#' assess these assumptions in your data.
#'
#' @section S3 Methods:
#'
#' Objects of class \code{cfomics_result} support the following S3 methods:
#'
#' \describe{
#'   \item{\code{predict(object, type, ...)}}{Extract predictions (ate, ite, y0, y1).}
#'   \item{\code{summary(object, ...)}}{Display summary with CIs.}
#'   \item{\code{print(object, ...)}}{Print overview of results.}
#'   \item{\code{plot(object, ...)}}{Visualize treatment effects.}
#' }
#'
#' @section Class Hierarchy:
#'
#' Objects have dual classes: \code{c("cf_model", "cfomics_result")}.
#' This allows for both generic and specific method dispatch.
#'
#' @name cfomics_result
#' @aliases cf_model
#' @seealso \code{\link{cf_fit}}, \code{\link{predict.cf_model}},
#'   \code{\link{summary.cf_model}}
NULL

#' Validate cfomics_result Structure
#'
#' Internal function to validate that a result object conforms to the
#' cfomics_result contract.
#'
#' @param x Object to validate
#' @param strict Logical, whether to enforce all optional fields (default FALSE)
#' @return TRUE if valid, throws error otherwise
#' @keywords internal
validate_cfomics_result <- function(x, strict = FALSE) {
  # Check class

  if (!inherits(x, "cfomics_result")) {
    rlang::abort(
      "Object must have class 'cfomics_result'",
      class = "cfomics_invalid_result"
    )
  }

  # Check required top-level components
  required_toplevel <- c("method", "fit", "meta")
  missing_toplevel <- setdiff(required_toplevel, names(x))
  if (length(missing_toplevel) > 0) {
    rlang::abort(
      sprintf("Missing required components: %s", paste(missing_toplevel, collapse = ", ")),
      class = "cfomics_invalid_result"
    )
  }

  # Check method
  if (!is.character(x$method) || length(x$method) != 1) {
    rlang::abort(
      "'method' must be a single character string",
      class = "cfomics_invalid_result"
    )
  }

  # Check meta
  required_meta <- c("formula", "n", "p", "outcome_name", "treatment_name",
                     "covariate_names")
  missing_meta <- setdiff(required_meta, names(x$meta))
  if (length(missing_meta) > 0) {
    rlang::abort(
      sprintf("Missing required meta fields: %s", paste(missing_meta, collapse = ", ")),
      class = "cfomics_invalid_result"
    )
  }

  # Check fit$res structure (if res exists)
  if (!is.null(x$fit$res)) {
    res <- x$fit$res
    required_res <- c("ite", "ate", "y0_hat", "y1_hat", "summary")
    missing_res <- setdiff(required_res, names(res))
    if (length(missing_res) > 0) {
      rlang::abort(
        sprintf("Missing required fit$res fields: %s", paste(missing_res, collapse = ", ")),
        class = "cfomics_invalid_result"
      )
    }

    # Validate summary structure
    if (!is.null(res$summary)) {
      required_summary <- c("ate", "ate_ci_lower", "ate_ci_upper")
      missing_summary <- setdiff(required_summary, names(res$summary))
      if (length(missing_summary) > 0) {
        rlang::abort(
          sprintf("Missing required summary fields: %s", paste(missing_summary, collapse = ", ")),
          class = "cfomics_invalid_result"
        )
      }
    }
  }

  invisible(TRUE)
}

#' Create Standard Meta Object
#'
#' Internal helper to create a standardized meta object for cfomics_result.
#'
#' @param formula The formula used for fitting
#' @param n Number of observations
#' @param p Number of covariates
#' @param outcome_name Name of outcome variable
#' @param treatment_name Name of treatment variable
#' @param covariate_names Character vector of covariate names
#' @param method Character string of method name
#' @param estimand Character, the estimand being estimated (default "ATE")
#' @param target_population Character, target population (default "all")
#' @param dag_spec DAG specification if provided
#' @param bootstrap Logical, whether bootstrap was used
#' @param random_state Random seed used
#' @param preprocess List of preprocessing steps
#' @return A list with standardized meta structure
#' @keywords internal
create_cf_meta <- function(formula,
                           n,
                           p,
                           outcome_name,
                           treatment_name,
                           covariate_names,
                           method,
                           estimand = "ATE",
                           target_population = "all",
                           dag_spec = NULL,
                           bootstrap = FALSE,
                           random_state = NULL,
                           preprocess = NULL) {

  # Capture software versions
  r_version <- paste(R.version$major, R.version$minor, sep = ".")
  pkg_versions <- list(
    cfomics = tryCatch(as.character(utils::packageVersion("cfomics")),
                       error = function(e) NA_character_),
    reticulate = tryCatch(as.character(utils::packageVersion("reticulate")),
                          error = function(e) NA_character_)
  )

  # Method-specific package versions
  method_pkg <- switch(method,
    grf = "grf",
    ipw = "ipw",
    gformula = NULL,
    NULL
  )

  if (!is.null(method_pkg)) {
    pkg_versions[[method_pkg]] <- tryCatch(
      as.character(utils::packageVersion(method_pkg)),
      error = function(e) NA_character_
    )
  }

  list(
    formula = formula,
    n = as.integer(n),
    p = as.integer(p),
    outcome_name = outcome_name,
    treatment_name = treatment_name,
    covariate_names = covariate_names,
    estimand = estimand,
    target_population = target_population,
    dag_spec = dag_spec,
    bootstrap = bootstrap,
    random_state = random_state,
    preprocess = preprocess,
    software_versions = list(
      R = r_version,
      cfomics = pkg_versions$cfomics,
      timestamp = Sys.time(),
      packages = pkg_versions
    )
  )
}
