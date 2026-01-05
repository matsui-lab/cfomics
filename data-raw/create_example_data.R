# Script to create example datasets for cfomics package
# Run this script to regenerate the example data in data/

# ==============================================================================
# Example 1: Simple data.frame for basic examples
# ==============================================================================

set.seed(42)
n <- 200

# Generate confounders
age <- rnorm(n, mean = 50, sd = 10)
batch <- sample(c("A", "B", "C"), n, replace = TRUE)
baseline_score <- rnorm(n)

# Treatment assignment (depends on confounders)
propensity <- plogis(-1 + 0.02 * age + 0.3 * baseline_score)
treatment <- rbinom(n, 1, propensity)

# Outcome (treatment effect = 2, with heterogeneity by age)
# True ITE = 2 + 0.05 * (age - 50)
true_ite <- 2 + 0.05 * (age - 50)
noise <- rnorm(n, sd = 1)
outcome <- 10 + 0.1 * age + 0.5 * baseline_score + treatment * true_ite + noise

# Create data.frame
cfomics_example_df <- data.frame(
  outcome = outcome,
  treatment = treatment,
  age = age,
  batch = batch,
  baseline_score = baseline_score,
  true_ite = true_ite
)

# Save as RData
usethis::use_data(cfomics_example_df, overwrite = TRUE)

# ==============================================================================
# Example 2: SummarizedExperiment (if Bioconductor available)
# ==============================================================================

if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  library(SummarizedExperiment)

  set.seed(42)
  n_samples <- 100
  n_genes <- 500

  # Generate gene expression matrix (genes x samples)
  # Some genes are affected by treatment, others are noise
  expression_matrix <- matrix(
    rnorm(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(expression_matrix) <- paste0("Gene", seq_len(n_genes))
  colnames(expression_matrix) <- paste0("Sample", seq_len(n_samples))

  # Sample metadata
  sample_age <- rnorm(n_samples, mean = 55, sd = 12)
  sample_treatment <- rbinom(n_samples, 1, 0.5)

  # Add treatment effect to first 50 genes
  treatment_effect_genes <- 1:50
  for (g in treatment_effect_genes) {
    # Gene-specific treatment effect
    gene_effect <- 0.5 + 0.3 * (g / 50)
    expression_matrix[g, sample_treatment == 1] <-
      expression_matrix[g, sample_treatment == 1] + gene_effect
  }

  # Create outcome as first gene expression (affected by treatment)
  outcome_gene <- expression_matrix[1, ]

  col_data <- DataFrame(
    sample_id = colnames(expression_matrix),
    treatment = sample_treatment,
    age = sample_age,
    outcome = outcome_gene,
    batch = sample(c("Batch1", "Batch2"), n_samples, replace = TRUE)
  )
  rownames(col_data) <- colnames(expression_matrix)

  # Row data (gene annotations)
  row_data <- DataFrame(
    gene_id = rownames(expression_matrix),
    is_treatment_affected = seq_len(n_genes) %in% treatment_effect_genes
  )
  rownames(row_data) <- rownames(expression_matrix)

  # Create SummarizedExperiment
  cfomics_example_se <- SummarizedExperiment(
    assays = list(counts = expression_matrix),
    colData = col_data,
    rowData = row_data
  )

  usethis::use_data(cfomics_example_se, overwrite = TRUE)
} else {
  message("SummarizedExperiment not available. Skipping SE example data.")
}

cat("Example data created successfully!\n")
cat("- cfomics_example_df: Basic data.frame example\n")
if (exists("cfomics_example_se")) {
  cat("- cfomics_example_se: SummarizedExperiment example\n")
}
