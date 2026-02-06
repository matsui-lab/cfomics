# benchmarks/tcga/download_tcga.R
# Download TCGA data using TCGAbiolinks

download_tcga_data <- function(project = "TCGA-BRCA", output_dir = "benchmarks/tcga/data") {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Package 'TCGAbiolinks' required. Install from Bioconductor.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' required.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Query RNA-seq data
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # Download
  TCGAbiolinks::GDCdownload(query, directory = file.path(output_dir, "GDCdata"))

  # Prepare data
  se <- TCGAbiolinks::GDCprepare(query, directory = file.path(output_dir, "GDCdata"))

  # Save
  saveRDS(se, file.path(output_dir, paste0(project, "_se.rds")))

  message(sprintf("Downloaded %s: %d samples, %d genes",
                  project, ncol(se), nrow(se)))

  invisible(se)
}

# Download multiple projects
if (sys.nframe() == 0) {
  projects <- c("TCGA-BRCA", "TCGA-LUAD", "TCGA-COAD")
  for (proj in projects) {
    message(sprintf("\n=== Downloading %s ===", proj))
    tryCatch(
      download_tcga_data(proj),
      error = function(e) message(sprintf("Error downloading %s: %s", proj, e$message))
    )
  }
}
