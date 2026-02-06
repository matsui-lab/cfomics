# benchmarks/tcga/preprocess_tcga.R
# Preprocess TCGA data for semi-synthetic benchmarks

preprocess_tcga <- function(se, n_genes = 1000, variance_filter = TRUE) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' required")
  }

  # Get expression matrix (log2 + 1 normalized)
  counts <- SummarizedExperiment::assay(se, "unstranded")
  expr <- log2(counts + 1)

  # Filter genes by variance
  if (variance_filter && n_genes < nrow(expr)) {
    gene_vars <- apply(expr, 1, var, na.rm = TRUE)
    top_genes <- order(gene_vars, decreasing = TRUE)[1:n_genes]
    expr <- expr[top_genes, ]
  }

  # Transpose: samples x genes
  X <- t(expr)

  # Get clinical data
  clinical <- as.data.frame(SummarizedExperiment::colData(se))

  # Extract useful variables
  result <- list(
    X = X,
    sample_ids = rownames(X),
    gene_ids = colnames(X),
    clinical = clinical
  )

  # Add TP53 mutation status if available
  if ("paper_TP53_mut_status" %in% names(clinical)) {
    result$tp53_mut <- as.integer(clinical$paper_TP53_mut_status == "Mut")
  }

  result
}

# Create semi-synthetic benchmark data
create_semi_synthetic_data <- function(tcga_data, dgp_params = list()) {
  X <- tcga_data$X
  n <- nrow(X)
  p <- ncol(X)

  # Use PCA for dimension reduction if needed
  n_target_components <- 100
  if (p > n_target_components) {
    pca <- prcomp(X, center = TRUE, scale. = TRUE)
    # PCA returns min(n-1, p) components; ensure we don't exceed available
    n_available <- ncol(pca$x)
    n_components <- min(n_target_components, n_available)
    X_reduced <- pca$x[, 1:n_components, drop = FALSE]
    message(sprintf("PCA: reduced %d features to %d components", p, n_components))
  } else {
    X_reduced <- scale(X)
  }

  # Synthetic treatment assignment
  tau_true <- dgp_params$tau %||% 2.0
  ps_coef <- dgp_params$ps_coef %||% c(0.5, 0.3, -0.2, rep(0, ncol(X_reduced) - 3))
  ps <- plogis(X_reduced %*% ps_coef[1:ncol(X_reduced)])
  T <- rbinom(n, 1, ps)

  # Synthetic outcome
  outcome_coef <- dgp_params$outcome_coef %||% c(1.0, -0.5, 0.3, rep(0, ncol(X_reduced) - 3))
  mu0 <- X_reduced %*% outcome_coef[1:ncol(X_reduced)]

  # Heterogeneous treatment effect
  if (dgp_params$heterogeneous %||% FALSE) {
    tau <- tau_true * (1 + 0.3 * X_reduced[, 1])
  } else {
    tau <- rep(tau_true, n)
  }

  Y <- mu0 + T * tau + rnorm(n, sd = dgp_params$noise_sd %||% 1.0)

  list(
    data = data.frame(Y = as.numeric(Y), T = T, X_reduced),
    truth = list(
      ate_true = mean(tau),
      ite_true = as.numeric(tau)
    ),
    X_original = X,
    propensity_score = as.numeric(ps)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
