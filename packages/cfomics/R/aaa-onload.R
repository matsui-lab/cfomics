# Declare global variables used in ggplot2 NSE to avoid R CMD check NOTEs
utils::globalVariables(c(
  "ITE", "balanced", "ci_lower", "ci_upper", "gate", "group", "ps", "smd", "variable"
))

#' @import methods
#' @importFrom stats model.frame as.formula var glm binomial predict rnorm rbinom setNames median reorder prcomp complete.cases plogis sd quantile lm coef confint
#' @importFrom Formula Formula
#' @importFrom rlang abort warn inform .data
#' @importFrom cli cli_alert_success cli_alert_warning cli_alert_info cli_abort cli_h1 cli_h2 cli_text cli_progress_step col_green col_red col_yellow
#' @importFrom glue glue
#' @importFrom utils askYesNo install.packages packageVersion
#' @importFrom graphics hist legend lines par plot abline
#' @importFrom grDevices rgb
#' @importFrom parallel detectCores
#' @keywords internal
"_PACKAGE"

# NOTE: S3 methods are registered in NAMESPACE via roxygen2 @export tags.
# No .onLoad needed for S3 method registration.
# Python initialization is deliberately NOT performed at package load.
# Python environments should be configured by the user before calling
# Python-based methods. Use cf_check_env() to verify environment status.
