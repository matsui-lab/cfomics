#!/usr/bin/env Rscript
# Generate documentation for cfomics package
#
# Run this script from the cfomics directory:
#   Rscript inst/scripts/generate_docs.R
#
# Or from R:
#   source("inst/scripts/generate_docs.R")

# Ensure we're in the right directory
if (!file.exists("DESCRIPTION")) {
  stop("Please run this script from the cfomics package directory")
}

# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

# Install roxygen2 if needed
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  install.packages("roxygen2", repos = "https://cloud.r-project.org")
}

# Generate documentation
message("Generating documentation with roxygen2...")
devtools::document()

# Check if man/ directory was created
if (dir.exists("man") && length(list.files("man", pattern = "\\.Rd$")) > 0) {
  message("Documentation generated successfully!")
  message(sprintf("Created %d Rd files in man/",
                  length(list.files("man", pattern = "\\.Rd$"))))
} else {
  warning("No Rd files were generated. Check roxygen2 comments in R files.")
}

# Optionally run R CMD check
if (interactive()) {
  run_check <- readline("Run R CMD check? (y/n): ")
  if (tolower(run_check) == "y") {
    devtools::check()
  }
}
