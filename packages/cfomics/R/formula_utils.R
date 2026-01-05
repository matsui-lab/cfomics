
parse_cf_formula <- function(formula, data) {
  # Early NA check on raw data before model.frame processes it
  all_v <- all.vars(formula)
  for (v in all_v) {
    if (v %in% names(data) && any(is.na(data[[v]]))) {
      stop(sprintf("NA values found in variable '%s'. Please remove or impute missing values.", v), call. = FALSE)
    }
  }

  f <- Formula::Formula(formula)
  mf <- stats::model.frame(f, data = data)

  # Outcome Y
  outcome_name <- all.vars(formula)[1]
  Y <- mf[[outcome_name]]
  if (!is.numeric(Y)) {
    stop("Outcome Y must be numeric.", call. = FALSE)
  }

  # Treatment T (assumed to be the first RHS variable or explicitly identified if we had syntax,
  # but spec says Y ~ T | X1 + ... implies first var on RHS is T?
  # Wait, Formula `Y ~ T | X1 + X2` implies 2 parts on RHS.
  # Part 1: T, Part 2: X.
  # Let's check the spec "Y ~ T | X1 + X2".
  # Formula uses `|` to separate parts.
  #
  # Spec 4.1 Implementation Image:
  # `rhs_vars <- all.vars(formula)[-1]`
  # `treatment_name <- rhs_vars[1]`
  # `covariate_names <- setdiff(rhs_vars, treatment_name)`
  #
  # This implementation assumes `all.vars` order.
  # `all.vars(Y ~ T | X1)` returns `c("Y", "T", "X1")`.
  # This works for simple formulas.
  #
  # However, strictly speaking using Formula package allows separating parts.
  # `f <- Formula(Y ~ T | X)`
  # `model.part(f, data = mf, rhs = 1)` -> T
  # `model.part(f, data = mf, rhs = 2)` -> X
  #
  # The spec implementation uses `all.vars`. I will follow the spec logic as it handles the names directly.
  #
  # `all.vars` extracts variables in order of appearance.
  # Y ~ T | X1 + X2 -> Y, T, X1, X2.

  all_v <- all.vars(formula)
  outcome_name <- all_v[1]

  # Check if T is separated by |
  # If user writes `Y ~ T + X1`, then T is just the first var?
  # The spec example is `Y ~ T | X1 + X2`.
  # But the code provided in 4.1 uses `all.vars` which flattens everything.
  # So it implicitly assumes T is the *first* variable on the RHS.

  rhs_vars <- all_v[-1]
  treatment_name <- rhs_vars[1]
  T <- mf[[treatment_name]]

  # T check
  if (is.logical(T)) T <- as.integer(T)
  uT <- sort(unique(T))
  if (!all(uT %in% c(0, 1))) {
    stop("Treatment must be binary (0/1) or logical.", call. = FALSE)
  }

  # Covariates X
  covariate_names <- setdiff(rhs_vars, treatment_name)
  X <- mf[, covariate_names, drop = FALSE]

  # v1 constraint: X must be numeric
  is_numeric_X <- vapply(X, is.numeric, logical(1))
  if (!all(is_numeric_X)) {
    bad <- covariate_names[!is_numeric_X]
    stop(
      paste0(
        "v1 of cfomics supports only numeric covariates.\n",
        "Please encode factors/characters manually before cf_fit().\n",
        "Non-numeric columns: ", paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Variable name constraint
  all_names <- c(outcome_name, treatment_name, covariate_names)
  if (any(grepl("[^A-Za-z0-9_]", all_names))) {
    stop(
      paste0(
        "v1 of cfomics requires variable names to use only letters, digits, and underscore.\n",
        "Please rename variables such as 'x.y' -> 'x_y' before using cf_fit()."
      ),
      call. = FALSE
    )
  }

  # Sample size check
  n <- length(Y)
  if (n == 0) {
    stop("Data has no observations.", call. = FALSE)
  }
  if (n < 10) {
    stop("Insufficient sample size. At least 10 observations required.", call. = FALSE)
  }

  # Treatment group check
  n_treated <- sum(T == 1)
  n_control <- sum(T == 0)
  if (n_treated == 0 || n_control == 0) {
    stop("Both treatment and control groups must have observations.", call. = FALSE)
  }

  list(
    Y = as.numeric(Y),
    T = as.integer(T),
    X = as.matrix(X),
    outcome_name    = outcome_name,
    treatment_name  = treatment_name,
    covariate_names = covariate_names
  )
}
