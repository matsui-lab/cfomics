# Tests for Bioconductor integration (SummarizedExperiment, MultiAssayExperiment)

# Helper function to create a mock SummarizedExperiment
create_test_se <- function(n_samples = 50, n_genes = 100) {
  skip_if_not_installed("SummarizedExperiment")

  set.seed(42)

  # Gene expression matrix (genes x samples)
  counts <- matrix(
    rnorm(n_genes * n_samples, mean = 10, sd = 3),
    nrow = n_genes, ncol = n_samples
  )
  rownames(counts) <- paste0("gene", 1:n_genes)
  colnames(counts) <- paste0("sample", 1:n_samples)

  # Sample metadata
  col_data <- data.frame(
    treatment = rbinom(n_samples, 1, 0.5),
    age = rnorm(n_samples, 50, 10),
    outcome = rnorm(n_samples),
    row.names = colnames(counts)
  )
  # Make outcome correlated with treatment and age
  col_data$outcome <- col_data$treatment * 1.5 + col_data$age * 0.05 + rnorm(n_samples)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )
}

# Helper function to create a mock MultiAssayExperiment
create_test_mae <- function(n_samples = 30, n_genes = 50, n_proteins = 40) {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  set.seed(42)

  # Shared sample names
  sample_names <- paste0("sample", 1:n_samples)

  # RNAseq experiment
  rnaseq <- matrix(
    rnorm(n_genes * n_samples, mean = 10, sd = 3),
    nrow = n_genes, ncol = n_samples
  )
  rownames(rnaseq) <- paste0("gene", 1:n_genes)
  colnames(rnaseq) <- sample_names

  # Proteomics experiment
  proteomics <- matrix(
    rnorm(n_proteins * n_samples, mean = 5, sd = 2),
    nrow = n_proteins, ncol = n_samples
  )
  rownames(proteomics) <- paste0("protein", 1:n_proteins)
  colnames(proteomics) <- sample_names

  # Sample metadata
  col_data <- data.frame(
    treatment = rbinom(n_samples, 1, 0.5),
    age = rnorm(n_samples, 50, 10),
    outcome = rnorm(n_samples),
    row.names = sample_names
  )
  col_data$outcome <- col_data$treatment * 1.5 + col_data$age * 0.05 + rnorm(n_samples)

  # Create experiment list
  exp_list <- list(
    rnaseq = SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = rnaseq)
    ),
    proteomics = SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = proteomics)
    )
  )

  # Sample map
  sample_map <- data.frame(
    assay = rep(c("rnaseq", "proteomics"), each = n_samples),
    primary = rep(sample_names, 2),
    colname = rep(sample_names, 2)
  )

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = exp_list,
    colData = col_data,
    sampleMap = sample_map
  )
}

# --- SummarizedExperiment tests ---

test_that("is_se correctly identifies SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 10, n_genes = 20)
  df <- data.frame(a = 1:10)


  expect_true(is_se(se))
  expect_false(is_se(df))
  expect_false(is_se(NULL))
  expect_false(is_se(list()))
})

test_that("as_cf_data.SummarizedExperiment extracts colData variables", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(se, formula = outcome ~ treatment | age)

  expect_equal(result$outcome_name, "outcome")
  expect_equal(result$treatment_name, "treatment")
  expect_equal(result$covariate_names, "age")
  expect_equal(length(result$Y), 50)
  expect_equal(length(result$T), 50)
  expect_true(all(result$T %in% c(0, 1)))
})

test_that("as_cf_data.SummarizedExperiment extracts gene as outcome", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(se, formula = gene1 ~ treatment | age)

  expect_equal(result$outcome_name, "gene1")
  expect_equal(result$treatment_name, "treatment")
  expect_equal(length(result$Y), 50)
})

test_that("as_cf_data.SummarizedExperiment fails for non-existent variable", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  expect_error(
    as_cf_data(se, formula = nonexistent ~ treatment | age),
    "not found"
  )
})

test_that("as_cf_data.SummarizedExperiment requires formula or outcome/treatment", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  expect_error(
    as_cf_data(se),
    "formula|outcome|treatment"
  )
})

test_that("as_cf_data.SummarizedExperiment with manual specification", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(
    se,
    outcome = "outcome",
    treatment = "treatment",
    covariates = "age"
  )

  expect_equal(result$outcome_name, "outcome")
  expect_equal(result$treatment_name, "treatment")
  expect_equal(result$covariate_names, "age")
})

test_that("as_cf_data.SummarizedExperiment with gene as outcome (manual)", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(
    se,
    outcome = "gene1",
    treatment = "treatment",
    covariates = "age"
  )

  expect_equal(result$outcome_name, "gene1")
})

test_that("as_cf_data.SummarizedExperiment with variance feature selection", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(
    se,
    outcome = "outcome",
    treatment = "treatment",
    covariates = "age",
    feature_select = "variance",
    n_features = 10
  )

  # Should have age + 10 features from assay
  expect_equal(length(result$covariate_names), 11)
  expect_true("age" %in% result$covariate_names)
  expect_equal(ncol(result$X), 11)
})

test_that("as_cf_data.SummarizedExperiment with PCA feature selection", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  result <- as_cf_data(
    se,
    outcome = "outcome",
    treatment = "treatment",
    covariates = "age",
    feature_select = "pca",
    n_features = 5
  )

  # Should have age + 5 PC components
  expect_equal(length(result$covariate_names), 6)
  expect_true("age" %in% result$covariate_names)
  expect_true(any(grepl("^PC", result$covariate_names)))
})

test_that("cf_fit works with SummarizedExperiment", {
  skip_if_not_installed("SummarizedExperiment")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  fit <- cf_fit(
    outcome ~ treatment | age,
    data = se,
    method = "gformula",
    nsimul = 20L
  )

  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "gformula")
  expect_equal(fit$meta$n, 50)
})

test_that("cf_fit works with SummarizedExperiment and feature selection", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("grf")

  se <- create_test_se(n_samples = 50, n_genes = 100)

  fit <- cf_fit(
    outcome ~ treatment | age,
    data = se,
    method = "grf",
    feature_select = "variance",
    n_features = 10
  )

  expect_s3_class(fit, "cf_model")
  # Should include age + 10 variance-selected features
  expect_equal(fit$meta$p, 11)
})

# --- MultiAssayExperiment tests ---

test_that("is_mae correctly identifies MultiAssayExperiment", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 10, n_genes = 20, n_proteins = 15)
  df <- data.frame(a = 1:10)

  expect_true(is_mae(mae))
  expect_false(is_mae(df))
  expect_false(is_mae(NULL))
})

test_that("as_cf_data.MultiAssayExperiment extracts colData variables", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 30, n_genes = 50, n_proteins = 40)

  # When all variables are in colData, parse like data.frame
  result <- as_cf_data(mae, formula = outcome ~ treatment | age)

  expect_equal(result$outcome_name, "outcome")
  expect_equal(result$treatment_name, "treatment")
  expect_equal(result$covariate_names, "age")
})

test_that("as_cf_data.MultiAssayExperiment with multi-omics features", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 30, n_genes = 50, n_proteins = 40)

  result <- as_cf_data(
    mae,
    outcome = "outcome",
    treatment = "treatment",
    covariates = "age",
    feature_select = "variance",
    n_features_per_assay = 5
  )

  expect_equal(result$outcome_name, "outcome")
  expect_equal(result$treatment_name, "treatment")
  expect_true("age" %in% result$covariate_names)
  # Should have features from both rnaseq and proteomics
  expect_true(any(grepl("^rnaseq_", result$covariate_names)))
  expect_true(any(grepl("^proteomics_", result$covariate_names)))
})

test_that("as_cf_data.MultiAssayExperiment requires outcome and treatment", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 30, n_genes = 50, n_proteins = 40)

  expect_error(
    as_cf_data(mae, outcome = "outcome"),  # missing treatment
    "treatment"
  )
})

test_that("cf_fit works with MultiAssayExperiment", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")

  mae <- create_test_mae(n_samples = 30, n_genes = 50, n_proteins = 40)

  fit <- cf_fit(
    outcome ~ treatment | age,
    data = mae,
    method = "gformula",
    nsimul = 20L
  )

  expect_s3_class(fit, "cf_model")
  expect_equal(fit$method, "gformula")
})

# --- Feature selection tests ---

test_that("select_features_internal with variance method", {
  set.seed(42)
  mat <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  rownames(mat) <- paste0("feature", 1:100)

  # Add some high-variance features
  mat[1, ] <- mat[1, ] * 10
  mat[2, ] <- mat[2, ] * 8

  result <- select_features_internal(mat, method = "variance", n = 10)

  expect_equal(nrow(result), 10)
  expect_true("feature1" %in% rownames(result))
  expect_true("feature2" %in% rownames(result))
})

test_that("select_features_internal with PCA method", {
  set.seed(42)
  mat <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  rownames(mat) <- paste0("feature", 1:100)

  result <- select_features_internal(mat, method = "pca", n = 5)

  expect_equal(nrow(result), 5)
  expect_true(all(grepl("^PC", rownames(result))))
})

test_that("select_features_internal returns all if n >= nrow", {
  set.seed(42)
  mat <- matrix(rnorm(10 * 50), nrow = 10, ncol = 50)
  rownames(mat) <- paste0("feature", 1:10)

  result <- select_features_internal(mat, method = "variance", n = 20)

  expect_equal(nrow(result), 10)  # All features returned
})
