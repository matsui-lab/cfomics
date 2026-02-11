#' Example Dataset for cfomics
#'
#' A simulated dataset for demonstrating causal inference methods.
#' Contains 200 observations with a binary treatment, continuous outcome,
#' and two confounding covariates.
#'
#' @format A data frame with 200 rows and 4 variables:
#' \describe{
#'   \item{Y}{Numeric. Continuous outcome variable.}
#'   \item{T}{Integer. Binary treatment indicator (0 or 1).}
#'   \item{X1}{Numeric. First confounding covariate.}
#'   \item{X2}{Numeric. Second confounding covariate.}
#' }
#'
#' @details
#' Data generating process:
#' \itemize{
#'   \item Propensity: logit(P(T=1)) = 0.5*X1 - 0.3*X2
#'   \item Outcome: Y = 2*T + 0.8*X1 + 0.4*X2 + noise
#'   \item True ATE = 2
#' }
#'
#' @source Simulated data (seed = 42)
#'
#' @examples
#' data(cfomics_example)
#' head(cfomics_example)
#' attr(cfomics_example, "true_ate")  # 2
"cfomics_example"
