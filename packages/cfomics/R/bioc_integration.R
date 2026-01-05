#' Convert Data to cfomics Format
#'
#' Generic function to convert various data structures to a format
#' suitable for cfomics analysis. Supports data.frame, SummarizedExperiment,
#' and MultiAssayExperiment objects.
#'
#' @param data Input data object
#' @param ... Additional arguments passed to methods
#' @return A list with components:
#'   \itemize{
#'     \item Y: numeric vector of outcomes
#'     \item T: integer vector of treatment assignments (0/1)
#'     \item X: numeric matrix of covariates
#'     \item outcome_name: character string
#'     \item treatment_name: character string
#'     \item covariate_names: character vector
#'   }
#' @export
#' @examples
#' \dontrun{
#' # With data.frame
#' result <- as_cf_data(my_df, Y ~ T | X1 + X2)
#'
#' # With SummarizedExperiment
#' result <- as_cf_data(se, Gene1 ~ treatment | age + sex)
#' }
as_cf_data <- function(data, ...) {
  UseMethod("as_cf_data")
}

#' @rdname as_cf_data
#' @param formula Formula specifying Y ~ T | X1 + X2 + ...
#' @export
as_cf_data.data.frame <- function(data, formula, ...) {
  # Delegate to existing parse_cf_formula
  parse_cf_formula(formula, data)
}

#' @rdname as_cf_data
#' @export
as_cf_data.default <- function(data, ...) {
  stop("as_cf_data() does not support objects of class: ",
       paste(class(data), collapse = ", "),
       "\nSupported classes: data.frame, SummarizedExperiment, MultiAssayExperiment",
       call. = FALSE)
}

#' Convert SummarizedExperiment to cfomics format
#'
#' @rdname as_cf_data
#' @param outcome Character, name of outcome variable in colData or assay
#' @param treatment Character, name of treatment variable in colData
#' @param covariates Character vector, names of covariate variables in colData
#' @param assay_name Name or index of assay to use (default 1)
#' @param feature_select Method for feature selection from assay:
#'   "none" (default), "variance" (top variance), or "pca"
#' @param n_features Number of features to select if feature_select != "none"
#' @export
as_cf_data.SummarizedExperiment <- function(data,
                                            formula = NULL,
                                            outcome = NULL,
                                            treatment = NULL,
                                            covariates = NULL,
                                            assay_name = 1,
                                            feature_select = c("none", "variance", "pca"),
                                            n_features = 100,
                                            ...) {
  # Check for SummarizedExperiment package
 if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required. Install with:\n",
         "  BiocManager::install('SummarizedExperiment')",
         call. = FALSE)
  }

  feature_select <- match.arg(feature_select)

  # Extract colData as data.frame
  col_data <- as.data.frame(SummarizedExperiment::colData(data))
  rownames(col_data) <- colnames(data)

  # If formula is provided, use formula-based extraction
  if (!is.null(formula)) {
    return(extract_from_formula_se(data, formula, col_data, assay_name,
                                   feature_select, n_features))
  }

  # Manual specification requires outcome and treatment
  if (is.null(outcome) || is.null(treatment)) {
    stop("Either 'formula' or both 'outcome' and 'treatment' must be specified",
         call. = FALSE)
  }

  extract_manual_se(data, outcome, treatment, covariates, col_data,
                    assay_name, feature_select, n_features)
}

#' Extract data from SummarizedExperiment using formula
#' @noRd
extract_from_formula_se <- function(data, formula, col_data, assay_name,
                                    feature_select, n_features) {
  vars <- all.vars(formula)
  outcome_name <- vars[1]

  # Check if outcome is in colData or assay
  if (outcome_name %in% colnames(col_data)) {
    # Outcome is a phenotype variable - use standard formula parsing
    return(parse_cf_formula(formula, col_data))
  }

  # Check if outcome is a feature (gene) name in assay
  assay_data <- SummarizedExperiment::assay(data, assay_name)

  if (outcome_name %in% rownames(assay_data)) {
    # Add feature value to colData and parse
    col_data[[outcome_name]] <- as.numeric(assay_data[outcome_name, ])
    return(parse_cf_formula(formula, col_data))
  }

  stop("Outcome '", outcome_name, "' not found in colData or assay rownames",
       call. = FALSE)
}

#' Extract data from SummarizedExperiment manually
#' @noRd
extract_manual_se <- function(data, outcome, treatment, covariates, col_data,
                              assay_name, feature_select, n_features) {
  assay_data <- SummarizedExperiment::assay(data, assay_name)

  # Extract outcome (Y)
  if (outcome %in% rownames(assay_data)) {
    Y <- as.numeric(assay_data[outcome, ])
  } else if (outcome %in% colnames(col_data)) {
    Y <- as.numeric(col_data[[outcome]])
  } else {
    stop("Outcome '", outcome, "' not found in assay rownames or colData",
         call. = FALSE)
  }

  # Extract treatment (T)
  if (!(treatment %in% colnames(col_data))) {
    stop("Treatment '", treatment, "' must be in colData", call. = FALSE)
  }
  T_vec <- col_data[[treatment]]

  # Convert to binary
  if (is.logical(T_vec)) T_vec <- as.integer(T_vec)
  if (is.factor(T_vec)) T_vec <- as.integer(T_vec) - 1L
  T_vec <- as.integer(T_vec)

  if (!all(unique(T_vec) %in% c(0L, 1L))) {
    stop("Treatment must be binary (0/1)", call. = FALSE)
  }

  # Extract covariates (X) from colData
  if (is.null(covariates)) {
    # Use all numeric columns except outcome and treatment
    numeric_cols <- sapply(col_data, is.numeric)
    covariates <- names(col_data)[numeric_cols]
    covariates <- setdiff(covariates, c(outcome, treatment))
  }

  # Filter to existing columns
  covariates <- intersect(covariates, colnames(col_data))

  # Build covariate matrix
  if (length(covariates) > 0) {
    X_coldata <- as.matrix(col_data[, covariates, drop = FALSE])
  } else {
    X_coldata <- matrix(nrow = ncol(data), ncol = 0)
  }

  # Feature selection from assay
  if (feature_select != "none" && n_features > 0) {
    feature_matrix <- select_features_internal(
      assay_data,
      method = feature_select,
      n = n_features
    )
    X <- cbind(X_coldata, t(feature_matrix))
    all_covariate_names <- c(covariates, rownames(feature_matrix))
  } else {
    X <- X_coldata
    all_covariate_names <- covariates
  }

  # Ensure X is numeric matrix
  X <- apply(X, 2, as.numeric)
  if (!is.matrix(X)) X <- matrix(X, ncol = 1)

  list(
    Y = Y,
    T = T_vec,
    X = X,
    outcome_name = outcome,
    treatment_name = treatment,
    covariate_names = all_covariate_names,
    source = "SummarizedExperiment",
    n_samples = length(Y),
    n_features_from_assay = if (feature_select != "none") n_features else 0L
  )
}

#' Convert MultiAssayExperiment to cfomics format
#'
#' @rdname as_cf_data
#' @param experiments Character vector of experiment names to include
#' @param n_features_per_assay Number of features to select from each assay
#' @export
as_cf_data.MultiAssayExperiment <- function(data,
                                            formula = NULL,
                                            outcome = NULL,
                                            treatment = NULL,
                                            covariates = NULL,
                                            experiments = NULL,
                                            feature_select = c("variance", "pca", "none"),
                                            n_features_per_assay = 50,
                                            ...) {
  # Check for MultiAssayExperiment package
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("MultiAssayExperiment package required. Install with:\n",
         "  BiocManager::install('MultiAssayExperiment')",
         call. = FALSE)
  }

  feature_select <- match.arg(feature_select)

  # Extract colData
  col_data <- as.data.frame(MultiAssayExperiment::colData(data))

  # If formula provided with colData variables only
  if (!is.null(formula)) {
    vars <- all.vars(formula)
    if (all(vars %in% colnames(col_data))) {
      return(parse_cf_formula(formula, col_data))
    }
  }

  # For multi-omics integration, need outcome and treatment
  if (is.null(outcome) || is.null(treatment)) {
    stop("For MultiAssayExperiment, 'outcome' and 'treatment' must be specified",
         call. = FALSE)
  }

  # Get experiments to use
  if (is.null(experiments)) {
    experiments <- names(MultiAssayExperiment::experiments(data))
  }

  # Extract and integrate features from each experiment
  integrated_features <- list()

  for (exp_name in experiments) {
    exp_data <- MultiAssayExperiment::experiments(data)[[exp_name]]

    # Get assay matrix depending on object type
    if (inherits(exp_data, "SummarizedExperiment")) {
      assay_mat <- SummarizedExperiment::assay(exp_data, 1)
    } else if (inherits(exp_data, "matrix")) {
      assay_mat <- exp_data
    } else {
      warning("Skipping experiment '", exp_name, "': unsupported type ",
              class(exp_data)[1])
      next
    }

    if (nrow(assay_mat) == 0) next

    # Select top features
    if (feature_select != "none" && n_features_per_assay > 0) {
      selected <- select_features_internal(
        assay_mat,
        method = feature_select,
        n = min(n_features_per_assay, nrow(assay_mat))
      )
    } else {
      selected <- assay_mat
    }

    # Prefix feature names with experiment name
    rownames(selected) <- paste0(exp_name, "_", rownames(selected))
    integrated_features[[exp_name]] <- selected
  }

  if (length(integrated_features) == 0) {
    stop("No features extracted from experiments", call. = FALSE)
  }

  # Find common samples across all experiments
  sample_lists <- lapply(integrated_features, colnames)
  common_samples <- Reduce(intersect, sample_lists)

  if (length(common_samples) == 0) {
    stop("No common samples across experiments", call. = FALSE)
  }

  # Subset to common samples and combine
  feature_matrix <- do.call(rbind, lapply(integrated_features, function(x) {
    x[, common_samples, drop = FALSE]
  }))

  # Subset colData to common samples
  col_data <- col_data[common_samples, , drop = FALSE]

  # Extract Y, T
  if (outcome %in% colnames(col_data)) {
    Y <- as.numeric(col_data[[outcome]])
  } else {
    stop("Outcome '", outcome, "' not found in colData", call. = FALSE)
  }

  T_vec <- col_data[[treatment]]
  if (is.logical(T_vec)) T_vec <- as.integer(T_vec)
  if (is.factor(T_vec)) T_vec <- as.integer(T_vec) - 1L
  T_vec <- as.integer(T_vec)

  if (!all(unique(T_vec) %in% c(0L, 1L))) {
    stop("Treatment must be binary (0/1)", call. = FALSE)
  }

  # Covariates from colData
  if (is.null(covariates)) {
    covariates <- character(0)
  }
  covariates <- intersect(covariates, colnames(col_data))

  if (length(covariates) > 0) {
    covar_data <- as.matrix(col_data[, covariates, drop = FALSE])
    covar_data <- apply(covar_data, 2, as.numeric)
  } else {
    covar_data <- matrix(nrow = nrow(col_data), ncol = 0)
  }

  # Combine covariates and features
  X <- cbind(covar_data, t(feature_matrix))
  all_covariate_names <- c(covariates, rownames(feature_matrix))

  list(
    Y = Y,
    T = T_vec,
    X = X,
    outcome_name = outcome,
    treatment_name = treatment,
    covariate_names = all_covariate_names,
    source = "MultiAssayExperiment",
    n_samples = length(Y),
    experiments_used = names(integrated_features),
    n_features_per_assay = n_features_per_assay
  )
}

#' Select Features from Expression Matrix
#'
#' Internal function to select top features from a matrix by variance or PCA.
#'
#' @param mat Numeric matrix (features x samples)
#' @param method Selection method: "variance" or "pca"
#' @param n Number of features/components to select
#' @return Numeric matrix with selected features/components
#' @noRd
select_features_internal <- function(mat, method = "variance", n = 100) {
  if (n >= nrow(mat)) {
    return(mat)
  }

  if (method == "variance") {
    # Select top n features by variance
    vars <- apply(mat, 1, stats::var, na.rm = TRUE)
    top_idx <- order(vars, decreasing = TRUE)[seq_len(n)]
    return(mat[top_idx, , drop = FALSE])
  }

  if (method == "pca") {
    # PCA - transpose so samples are rows
    mat_clean <- mat[complete.cases(mat), , drop = FALSE]
    if (nrow(mat_clean) < n) {
      n <- nrow(mat_clean)
    }

    pca <- stats::prcomp(t(mat_clean), center = TRUE, scale. = TRUE)
    n_pcs <- min(n, ncol(pca$x))
    pc_features <- t(pca$x[, seq_len(n_pcs), drop = FALSE])
    rownames(pc_features) <- paste0("PC", seq_len(n_pcs))
    return(pc_features)
  }

  mat
}

#' Check if object is SummarizedExperiment
#' @noRd
is_se <- function(x) {
  inherits(x, "SummarizedExperiment")
}

#' Check if object is MultiAssayExperiment
#' @noRd
is_mae <- function(x) {
  inherits(x, "MultiAssayExperiment")
}
