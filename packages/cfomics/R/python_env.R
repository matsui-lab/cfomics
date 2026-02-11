# Environment specifications for each method
.cfomics_env_specs <- list(
  # Unified environment with all dependencies (recommended)
  unified = list(
    name = "cfomics_unified",
    packages = c("numpy", "pandas", "scikit-learn", "scipy", "networkx", "matplotlib"),
    pip_packages = c("dowhy==0.11.1", "econml")
  ),

  # Unified environment with deep learning (full)
  unified_full = list(
    name = "cfomics_unified_full",
    packages = c("numpy", "pandas", "scikit-learn", "scipy", "networkx", "matplotlib"),
    pip_packages = c("dowhy==0.11.1", "econml", "tensorflow", "torch", "pyro-ppl")
  ),

  # Method-specific environments (legacy)
  dowhy = list(
    name = "cfomics_dowhy",
    packages = c("numpy", "pandas", "scikit-learn", "networkx", "matplotlib"),
    pip_packages = c("dowhy==0.11.1")
  ),

  ganite = list(
    name = "cfomics_ganite",
    packages = c("numpy", "pandas"),
    pip_packages = c("tensorflow")
  ),

  econml = list(
    name = "cfomics_econml",
    packages = c("numpy", "pandas", "scikit-learn"),
    pip_packages = c("econml")
  ),

  pyro = list(
    name = "cfomics_pyro",
    packages = c("numpy", "pandas"),
    pip_packages = c("torch", "pyro-ppl")
  )
)

#' Install Python environment for cfomics
#'
#' Creates conda environments with the required Python packages for
#' causal inference methods. The "unified" environment is recommended
#' as it contains dependencies for most methods.
#'
#' @param method The environment to install:
#'   \itemize{
#'     \item "unified": Recommended. Contains DoWhy and EconML (no deep learning)
#'     \item "unified_full": Full environment with TensorFlow, PyTorch, Pyro
#'     \item "dowhy", "ganite", "econml", "pyro": Method-specific environments
#'     \item "all": Install all method-specific environments (legacy)
#'   }
#' @param python_version The Python version to use (default "3.11")
#' @param force Logical, whether to recreate the environment if it exists
#' @return Invisibly returns TRUE on success
#' @export
#' @examples
#' \dontrun{
#' cf_install_python_env("unified")       # Recommended
#' cf_install_python_env("unified_full")  # With deep learning
#' cf_install_python_env("dowhy")         # DoWhy only
#' }
cf_install_python_env <- function(method = c("unified", "unified_full", "dowhy",
                                             "ganite", "econml", "pyro", "all"),
                                  python_version = "3.11",
                                  force = FALSE) {
  method <- match.arg(method)
  
  tryCatch({
    reticulate::install_miniconda()
  }, error = function(e) {
    message("Miniconda installation skipped or failed: ", e$message)
  })
  
  methods_to_install <- if (method == "all") {
    c("dowhy", "ganite", "econml", "pyro")
  } else {
    method
  }
  
  for (m in methods_to_install) {
    spec <- .cfomics_env_specs[[m]]
    envname <- spec$name
    
    existing_envs <- tryCatch(
      reticulate::conda_list()$name,
      error = function(e) character(0)
    )
    
    if (envname %in% existing_envs && !force) {
      message(sprintf("Environment '%s' already exists. Use force=TRUE to recreate.", envname))
      next
    }
    
    if (envname %in% existing_envs && force) {
      message(sprintf("Removing existing environment '%s'...", envname))
      reticulate::conda_remove(envname)
    }
    
    message(sprintf("Creating environment '%s' for method '%s'...", envname, m))
    
    reticulate::conda_create(
      envname,
      packages = c(sprintf("python=%s", python_version), spec$packages)
    )
    
    if (length(spec$pip_packages) > 0) {
      reticulate::conda_install(
        envname,
        packages = spec$pip_packages,
        pip = TRUE
      )
    }
    
    message(sprintf("Environment '%s' created successfully.", envname))
  }
  
  invisible(TRUE)
}

#' Use the cfomics Python environment
#'
#' Activates a cfomics conda environment for the current R session.
#'
#' @param method The environment to use:
#'   \itemize{
#'     \item "unified": Recommended unified environment
#'     \item "unified_full": Full environment with deep learning
#'     \item "dowhy", "ganite", "econml", "pyro": Method-specific environments
#'   }
#' @return Invisibly returns TRUE on success
#' @export
#' @examples
#' \dontrun{
#' cf_use_python_env("unified")
#' cf_use_python_env("dowhy")
#' }
cf_use_python_env <- function(method = c("unified", "unified_full", "dowhy",
                                         "ganite", "econml", "pyro")) {
  method <- match.arg(method)
  spec <- .cfomics_env_specs[[method]]
  envname <- spec$name
  reticulate::use_condaenv(envname, required = TRUE)
  invisible(TRUE)
}

#' Setup Python environment for a specific method
#'
#' This is a convenience function that checks if the environment exists,
#' installs it if needed, and activates it.
#'
#' @param method The environment to setup:
#'   \itemize{
#'     \item "unified": Recommended unified environment (default)
#'     \item "unified_full": Full environment with deep learning
#'     \item "dowhy", "ganite", "econml", "pyro": Method-specific environments
#'   }
#' @param python_version The Python version to use if installing
#' @return Invisibly returns TRUE on success
#' @export
#' @examples
#' \dontrun{
#' setup_python_env("unified")
#' setup_python_env("dowhy")
#' }
setup_python_env <- function(method = c("unified", "unified_full", "dowhy",
                                        "ganite", "econml", "pyro"),
                             python_version = "3.11") {
  method <- match.arg(method)
  spec <- .cfomics_env_specs[[method]]
  envname <- spec$name
  
  existing_envs <- tryCatch(
    reticulate::conda_list()$name,
    error = function(e) character(0)
  )
  
  if (!(envname %in% existing_envs)) {
    message(sprintf("Environment '%s' not found. Installing...", envname))
    cf_install_python_env(method = method, python_version = python_version)
  }
  
  cf_use_python_env(method = method)
  invisible(TRUE)
}

#' List available cfomics Python environments
#'
#' @return A data.frame with environment names and their status
#' @export
cf_list_python_envs <- function() {
  existing_envs <- tryCatch(
    reticulate::conda_list()$name,
    error = function(e) character(0)
  )
  
  env_status <- data.frame(
    method = names(.cfomics_env_specs),
    env_name = vapply(.cfomics_env_specs, function(x) x$name, character(1)),
    installed = vapply(.cfomics_env_specs, function(x) x$name %in% existing_envs, logical(1)),
    stringsAsFactors = FALSE
  )
  
  env_status
}

#' Ensure Python environment is available
#'
#' This function checks that a Python environment is available for use.
#' cfomics does NOT manage Python environments - users must prepare their own
#' conda/mamba/virtualenv environment with the required dependencies.
#'
#' The function respects the following priority:
#' 1. If Python is already initialized (e.g., via use_condaenv()), use it
#' 2. If RETICULATE_PYTHON is set, use that Python
#' 3. Otherwise, error with instructions
#'
#' @return Invisible TRUE on success
#' @export
ensure_python_env <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("reticulate is required")
  }


  if (reticulate::py_available(initialize = FALSE)) {
    return(invisible(TRUE))
  }


  if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
    return(invisible(TRUE))
  }

  stop(
    "No Python environment configured.\n",
    "cfomics does not manage Python environments.\n",
    "Please prepare a conda/mamba environment with required dependencies and set:\n",
    "  export RETICULATE_PYTHON=/path/to/python\n",
    "Or use reticulate::use_condaenv() / reticulate::use_python() before calling cfomics functions.\n",
    "See the package README for details on required Python packages.",
    call. = FALSE
  )
}

#' Ensure conda-based Python environment is available (for specific methods)
#'
#' This function uses conda environments for method-specific dependencies.
#' Use this when you need specific method backends (dowhy, ganite, econml, pyro).
#'
#' @param method The method to setup: "dowhy", "ganite", "econml", or "pyro"
#' @return Invisible TRUE on success
#' @export
ensure_conda_env <- function(method = "dowhy") {
  if (reticulate::py_available(initialize = FALSE)) return(invisible(TRUE))
  
  if (!(method %in% names(.cfomics_env_specs))) {
    method <- "dowhy"
  }
  
  spec <- .cfomics_env_specs[[method]]
  envname <- spec$name
  
  ok <- tryCatch({
    reticulate::use_condaenv(envname, required = FALSE)
    reticulate::py_available(initialize = TRUE)
  }, error = function(e) FALSE)
  
  if (!ok) {
    msg <- paste0(
      "Python environment '", envname, "' for method '", method, "' is not available.\n",
      " - Run cf_install_python_env(method = '", method, "') to create it,\n",
      " - Or run setup_python_env(method = '", method, "') to install and activate.\n"
    )
    
    if (interactive()) {
      ans <- utils::askYesNo(
        paste0(msg, "\nDo you want to install it now?")
      )
      if (isTRUE(ans)) {
        cf_install_python_env(method = method)
        reticulate::use_condaenv(envname, required = TRUE)
        return(invisible(TRUE))
      }
    }
    
    stop(msg, call. = FALSE)
  }

  invisible(TRUE)
}

# =============================================================================
# Lazy Python Initialization (Phase 3 API)
# =============================================================================

#' Check if Python is Available for a Method
#'
#' Checks whether Python and required packages are available for a specific
#' cfomics method WITHOUT initializing Python. This function is safe to call
#' at any time and will not trigger Python initialization or conda environment
#' creation.
#'
#' @param method Character string specifying the method to check. One of:
#'   \itemize{
#'     \item \code{"any"}: Check if any Python is available (default)
#'     \item \code{"dowhy_gcm"}: Check for DoWhy package
#'     \item \code{"drlearner"}: Check for EconML package
#'     \item \code{"ganite"}: Check for TensorFlow package
#'     \item \code{"cavae"}: Check for Pyro package
#'   }
#' @return Logical indicating whether Python (and required packages for the
#'   specified method) is available.
#'
#' @details
#' This function is designed to be called before attempting to use Python-based
#' methods. It allows code to gracefully degrade or provide informative messages
#' when Python is not available.
#'
#' For method-specific checks, Python must be already initialized (e.g., via
#' a prior call to \code{reticulate::use_condaenv()} or environment variable).
#' If Python is not initialized, only the \code{"any"} check will attempt to
#' detect Python availability.
#'
#' @export
#' @examples
#' # Check if any Python is available
#' if (cf_has_python()) {
#'   message("Python is available")
#' }
#'
#' # Check for specific method
#' if (cf_has_python("dowhy_gcm")) {
#'   message("DoWhy is available")
#' }
#'
#' # Use in conditional code
#' \dontrun{
#' if (cf_has_python("grf")) {
#'   fit <- cf_fit(Y ~ T | X, data = df, method = "grf")
#' } else {
#'   fit <- cf_fit(Y ~ T | X, data = df, method = "gformula")
#' }
#' }
cf_has_python <- function(method = c("any", "dowhy_gcm", "drlearner", "ganite", "cavae")) {
  method <- match.arg(method)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }

  # Check if Python is available (without forcing initialization)
  py_available <- tryCatch({
    reticulate::py_available(initialize = FALSE)
  }, error = function(e) FALSE)

  # For "any", just check Python availability
  if (method == "any") {
    # If not initialized, try to detect if Python could be available
    if (!py_available) {
      # Check if RETICULATE_PYTHON is set
      if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
        return(TRUE)
      }
      # Check if any cfomics conda env exists
      cfomics_envs <- tryCatch({
        envs <- reticulate::conda_list()
        any(grepl("cfomics", envs$name, ignore.case = TRUE))
      }, error = function(e) FALSE)
      return(cfomics_envs)
    }
    return(py_available)
  }

  # For method-specific checks, Python must be initialized
  if (!py_available) {
    return(FALSE)
  }

  # Map method to required Python package
  required_pkg <- switch(method,
    dowhy_gcm = "dowhy",
    drlearner = "econml",
    ganite = "tensorflow",
    cavae = "pyro",
    NULL
  )

  if (is.null(required_pkg)) {
    return(FALSE)
  }

  # Check if the required package is available
  tryCatch({
    reticulate::py_module_available(required_pkg)
  }, error = function(e) FALSE)
}

#' Require Python for a Method
#'
#' Ensures that Python and required packages are available for a specific
#' cfomics method. If not available, raises an informative error with
#' instructions for setup.
#'
#' Unlike \code{\link{cf_has_python}}, this function may initialize Python
#' if it hasn't been initialized yet and an appropriate environment can be
#' found.
#'
#' @param method Character string specifying the method that requires Python.
#'   One of: \code{"dowhy_gcm"}, \code{"drlearner"}, \code{"ganite"},
#'   \code{"cavae"}.
#' @param initialize Logical indicating whether to attempt Python initialization
#'   if not already initialized. Default is \code{TRUE}.
#'
#' @return Invisibly returns \code{TRUE} if Python is available. Throws an
#'   error if Python or required packages are not available.
#'
#' @details
#' This function is called internally by Python-based method implementations
#' (e.g., \code{cf_fit_dowhy_gcm}) to ensure the environment is properly
#' configured before attempting to run Python code.
#'
#' The function checks:
#' \enumerate{
#'   \item Whether Python is available
#'   \item Whether the required Python package for the method is installed
#' }
#'
#' If checks fail, an error is raised with:
#' \itemize{
#'   \item A clear description of what's missing
#'   \item Instructions on how to set up the environment
#'   \item A suggestion to run \code{cf_check_env()} for diagnostics
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Called internally by cf_fit when method = "dowhy_gcm"
#' cf_require_python("dowhy_gcm")
#'
#' # After this call, Python is guaranteed to be available
#' # and dowhy package can be imported
#' }
cf_require_python <- function(method = c("dowhy_gcm", "drlearner", "ganite", "cavae"),
                               initialize = TRUE) {
  method <- match.arg(method)

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    rlang::abort(
      message = "Package 'reticulate' is required for Python methods.",
      class = "cfomics_missing_reticulate",
      i = "Install it with: install.packages('reticulate')"
    )
  }

  # Map method to required Python package and environment
  pkg_info <- switch(method,
    dowhy_gcm = list(pkg = "dowhy", env = "dowhy", display = "DoWhy"),
    drlearner = list(pkg = "econml", env = "econml", display = "EconML"),
    ganite = list(pkg = "tensorflow", env = "ganite", display = "TensorFlow"),
    cavae = list(pkg = "pyro", env = "pyro", display = "Pyro")
  )

  # Check/initialize Python
  py_available <- tryCatch({
    reticulate::py_available(initialize = FALSE)
  }, error = function(e) FALSE)

  if (!py_available && initialize) {
    # Try to find and use a cfomics environment
    py_available <- tryCatch({
      # First try unified environment
      envs <- reticulate::conda_list()
      if ("cfomics_unified" %in% envs$name) {
        reticulate::use_condaenv("cfomics_unified", required = TRUE)
        TRUE
      } else if ("cfomics_unified_full" %in% envs$name) {
        reticulate::use_condaenv("cfomics_unified_full", required = TRUE)
        TRUE
      } else {
        # Try method-specific environment
        spec <- .cfomics_env_specs[[pkg_info$env]]
        if (!is.null(spec) && spec$name %in% envs$name) {
          reticulate::use_condaenv(spec$name, required = TRUE)
          TRUE
        } else {
          FALSE
        }
      }
    }, error = function(e) FALSE)
  }

  if (!py_available) {
    rlang::abort(
      message = sprintf("Python is not configured for method '%s'.", method),
      class = "cfomics_no_python",
      i = c(
        "cfomics does not auto-configure Python environments.",
        sprintf("To use %s, you must set up Python manually:", pkg_info$display),
        "",
        "Option 1: Use cfomics helper functions",
        sprintf("  cf_install_python_env('%s')", pkg_info$env),
        sprintf("  cf_use_python_env('%s')", pkg_info$env),
        "",
        "Option 2: Use your own conda/virtual environment",
        "  reticulate::use_condaenv('your_env', required = TRUE)",
        "",
        "Option 3: Set RETICULATE_PYTHON environment variable",
        "  Sys.setenv(RETICULATE_PYTHON = '/path/to/python')",
        "",
        "Run cf_check_env() to diagnose your environment."
      )
    )
  }

  # Check for required package
  pkg_available <- tryCatch({
    reticulate::py_module_available(pkg_info$pkg)
  }, error = function(e) FALSE)

  if (!pkg_available) {
    rlang::abort(
      message = sprintf(
        "Python package '%s' is not installed (required for method '%s').",
        pkg_info$pkg, method
      ),
      class = "cfomics_missing_python_package",
      x = sprintf("The '%s' Python package was not found.", pkg_info$pkg),
      i = c(
        sprintf("To install %s, run:", pkg_info$display),
        sprintf("  cf_install_python_env('%s')", pkg_info$env),
        "",
        "Or install manually in your Python environment:",
        sprintf("  pip install %s", pkg_info$pkg),
        "",
        "Run cf_check_env() to see all available packages."
      )
    )
  }

  invisible(TRUE)
}
