# Helper to configure Python environment for tests
# This file is sourced before tests run

if (requireNamespace("reticulate", quietly = TRUE)) {
  # Priority 1: Environment variable

  cfomics_python_path <- Sys.getenv("CFOMICS_PYTHON_PATH", unset = "")

  if (nzchar(cfomics_python_path) && file.exists(cfomics_python_path)) {
    message("Using Python from CFOMICS_PYTHON_PATH: ", cfomics_python_path)
    reticulate::use_python(cfomics_python_path, required = FALSE)
  } else {
    # Priority 2: Try to find cfomics conda environments
    tryCatch({
      envs <- reticulate::conda_list()
      cfomics_envs <- envs[grepl("cfomics", envs$name, ignore.case = TRUE), ]

      if (nrow(cfomics_envs) > 0) {
        # Prefer unified environment, then dowhy, then any cfomics env
        preferred_order <- c("cfomics_unified", "cfomics_dowhy", "cfomics_econml")
        for (env_name in preferred_order) {
          if (env_name %in% cfomics_envs$name) {
            message("Using conda environment: ", env_name)
            reticulate::use_condaenv(env_name, required = FALSE)
            break
          }
        }
        # If none of preferred found, use first available
        if (!reticulate::py_available(initialize = FALSE)) {
          message("Using conda environment: ", cfomics_envs$name[1])
          reticulate::use_condaenv(cfomics_envs$name[1], required = FALSE)
        }
      }
    }, error = function(e) {
      # Conda not available or other error - continue without Python
      message("Conda not available. R-only tests will run.")
    })
  }

  # Try to initialize Python (silently fail if not available)
  tryCatch({
    reticulate::py_available(initialize = TRUE)
  }, error = function(e) {
    message("Python initialization failed. R-only tests will run.")
  })
}

# Helper function to check if Python method is available
skip_if_no_python <- function() {
  if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python not available")
  }
}

# Helper function to check if specific Python package is available
skip_if_no_python_package <- function(package) {
  skip_if_no_python()
  if (!reticulate::py_module_available(package)) {
    testthat::skip(paste0("Python package '", package, "' not available"))
  }
}
