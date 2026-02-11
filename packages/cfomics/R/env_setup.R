#' One-Command Environment Setup for cfomics
#'
#' Sets up the complete Python environment for cfomics with a single command.
#' Automatically installs Miniconda if needed and creates the appropriate
#' conda environment with all dependencies.
#'
#' @param level Setup level: "minimal" (R-only), "standard" (common methods),
#'   or "full" (all methods including deep learning)
#' @param gpu Logical, whether to install GPU-enabled packages (TensorFlow, PyTorch)
#' @param force Logical, whether to recreate environments if they exist
#' @param verbose Logical, whether to show progress messages
#' @return Invisibly returns TRUE on success
#' @export
#' @examples
#' \dontrun{
#' cf_setup()                      # Standard setup
#' cf_setup(level = "minimal")     # R-only, no Python
#' cf_setup(level = "full", gpu = TRUE)  # Everything with GPU
#' }
cf_setup <- function(level = c("standard", "minimal", "full"),
                     gpu = FALSE,
                     force = FALSE,
                     verbose = TRUE) {
  level <- match.arg(level)

  if (verbose) {
    cat("\n")
    cat("=======================================================\n")
    cat("              cfomics Setup Wizard                     \n")
    cat("=======================================================\n")
    cat(sprintf("Setup level: %s\n", level))
    cat(sprintf("GPU support: %s\n", if (gpu) "Yes" else "No"))
    cat("\n")
  }

  # Step 1: Check R packages
  if (verbose) cat("[1/4] Checking R dependencies...\n")
  install_r_dependencies(level, verbose)

  # Step 2: For minimal, we're done
  if (level == "minimal") {
    if (verbose) {
      cat("\n")
      cat("[OK] Minimal setup complete!\n")
      cat("Available methods: grf, ipw, gformula\n")
      cat("\nRun cf_check_env() to verify installation.\n")
    }
    return(invisible(TRUE))
  }

  # Step 3: Ensure Miniconda
  if (verbose) cat("[2/4] Checking Python/Conda installation...\n")
  ensure_miniconda_installed(verbose)

  # Step 4: Create unified environment
  if (verbose) cat("[3/4] Creating cfomics Python environment...\n")
  create_unified_environment(level, gpu, force, verbose)

  # Step 5: Activate and verify
  if (verbose) cat("[4/4] Verifying installation...\n")
  activate_cfomics_env(verbose)

  # Final check
  if (verbose) {
    check <- cf_check_env(verbose = FALSE)
    n_available <- sum(vapply(check$methods, `[[`, logical(1), "available"))

    cat("\n")
    cat("=======================================================\n")
    cat("[OK] Setup complete!\n")
    cat(sprintf("     %d of %d methods available\n",
                n_available, length(check$methods)))
    cat("\nRun cf_check_env() to see details.\n")
    cat("=======================================================\n")
  }

  invisible(TRUE)
}

#' Install R dependencies
#' @noRd
install_r_dependencies <- function(level, verbose = TRUE) {
  required <- c("Formula", "igraph", "reticulate", "rlang", "cli", "glue")

  suggested <- switch(level,
    minimal = c("grf", "ipw"),
    standard = c("grf", "ipw", "ggplot2"),
    full = c("grf", "ipw", "survey", "ggplot2", "gridExtra")
  )

  all_pkgs <- c(required, suggested)

  for (pkg in all_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (verbose) cat(sprintf("     Installing R package: %s\n", pkg))
      tryCatch({
        utils::install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      }, error = function(e) {
        if (verbose) cat(sprintf("     [!!] Failed to install %s: %s\n", pkg, e$message))
      })
    }
  }

  if (verbose) cat("     R dependencies OK\n")
}

#' Ensure Miniconda is installed
#' @noRd
ensure_miniconda_installed <- function(verbose = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    if (verbose) cat("     Installing reticulate...\n")
    utils::install.packages("reticulate", repos = "https://cloud.r-project.org", quiet = TRUE)
  }

  conda <- tryCatch({
    reticulate::conda_binary()
  }, error = function(e) NULL)

  if (is.null(conda)) {
    if (verbose) cat("     Installing Miniconda...\n")
    tryCatch({
      reticulate::install_miniconda()
      if (verbose) cat("     Miniconda installed successfully\n")
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("     [!!] Failed to install Miniconda: %s\n", e$message))
        cat("     Please install Miniconda manually from:\n")
        cat("     https://docs.conda.io/en/latest/miniconda.html\n")
      }
      stop("Miniconda installation failed", call. = FALSE)
    })
  } else {
    if (verbose) cat(sprintf("     Conda found: %s\n", conda))
  }
}

#' Create unified Python environment
#' @noRd
create_unified_environment <- function(level, gpu, force, verbose = TRUE) {
  envname <- "cfomics_unified"

  # Check if exists
  existing <- tryCatch({
    reticulate::conda_list()$name
  }, error = function(e) character(0))

  if (envname %in% existing) {
    if (force) {
      if (verbose) cat(sprintf("     Removing existing environment: %s\n", envname))
      tryCatch({
        reticulate::conda_remove(envname)
      }, error = function(e) {
        if (verbose) cat(sprintf("     [!!] Failed to remove: %s\n", e$message))
      })
    } else {
      if (verbose) {
        cat(sprintf("     Environment '%s' already exists\n", envname))
        cat("     Use force=TRUE to recreate\n")
      }
      return(invisible(TRUE))
    }
  }

  # Base packages
  conda_packages <- c("python=3.11", "numpy", "pandas", "scipy", "scikit-learn")
  pip_packages <- c()

  # Level-specific packages
  if (level %in% c("standard", "full")) {
    pip_packages <- c(pip_packages, "dowhy==0.11.1", "econml")
  }

  if (level == "full") {
    if (gpu) {
      pip_packages <- c(pip_packages, "tensorflow", "torch", "pyro-ppl")
      if (verbose) cat("     Note: GPU packages will be installed\n")
    } else {
      pip_packages <- c(pip_packages, "tensorflow-cpu", "torch", "pyro-ppl")
    }
  }

  # Create environment
  if (verbose) {
    cat(sprintf("     Creating environment: %s\n", envname))
    cat(sprintf("     Conda packages: %s\n", paste(conda_packages, collapse = ", ")))
    if (length(pip_packages) > 0) {
      cat(sprintf("     Pip packages: %s\n", paste(pip_packages, collapse = ", ")))
    }
  }

  tryCatch({
    reticulate::conda_create(envname, packages = conda_packages)

    if (length(pip_packages) > 0) {
      reticulate::conda_install(envname, packages = pip_packages, pip = TRUE)
    }

    if (verbose) cat("     Environment created successfully\n")
  }, error = function(e) {
    if (verbose) cat(sprintf("     [!!] Failed to create environment: %s\n", e$message))
    stop("Environment creation failed", call. = FALSE)
  })

  invisible(TRUE)
}

#' Activate cfomics environment
#' @noRd
activate_cfomics_env <- function(verbose = TRUE) {
  envname <- "cfomics_unified"

  tryCatch({
    reticulate::use_condaenv(envname, required = TRUE)
    if (verbose) cat(sprintf("     Activated environment: %s\n", envname))
  }, error = function(e) {
    # Try method-specific env as fallback
    fallback_envs <- c("cfomics_dowhy", "cfomics_econml")
    for (env in fallback_envs) {
      result <- tryCatch({
        reticulate::use_condaenv(env, required = TRUE)
        if (verbose) cat(sprintf("     Activated fallback environment: %s\n", env))
        TRUE
      }, error = function(e) FALSE)
      if (result) return(invisible(TRUE))
    }
    if (verbose) cat("     [!!] Could not activate cfomics environment\n")
  })

  invisible(TRUE)
}

#' Reset cfomics Environment
#'
#' Removes all cfomics conda environments for a clean reinstall.
#'
#' @param confirm Logical, whether to ask for confirmation (default TRUE)
#' @return Invisibly returns TRUE on success
#' @export
#' @examples
#' \dontrun{
#' cf_reset_env()
#' cf_reset_env(confirm = FALSE)  # Skip confirmation
#' }
cf_reset_env <- function(confirm = TRUE) {
  envs <- c("cfomics_unified", "cfomics_dowhy", "cfomics_ganite",
            "cfomics_econml", "cfomics_pyro")

  existing <- tryCatch({
    reticulate::conda_list()$name
  }, error = function(e) character(0))

  to_remove <- intersect(envs, existing)

  if (length(to_remove) == 0) {
    message("No cfomics environments found.")
    return(invisible(TRUE))
  }

  cat(sprintf("Found %d cfomics environment(s): %s\n",
              length(to_remove), paste(to_remove, collapse = ", ")))

  if (confirm && interactive()) {
    cat("This will remove all cfomics Python environments.\n")
    ans <- utils::askYesNo("Proceed?")
    if (!isTRUE(ans)) {
      message("Cancelled.")
      return(invisible(FALSE))
    }
  }

  for (env in to_remove) {
    cat(sprintf("Removing %s...\n", env))
    tryCatch({
      reticulate::conda_remove(env)
    }, error = function(e) {
      cat(sprintf("[!!] Failed to remove %s: %s\n", env, e$message))
    })
  }

  message("Environment reset complete.")
  invisible(TRUE)
}

#' Quick Setup for Specific Method
#'
#' Sets up the environment for a specific causal inference method.
#'
#' @param method The method to set up: "dowhy", "ganite", "econml", or "pyro"
#' @param force Logical, whether to recreate if exists
#' @return Invisibly returns TRUE on success
#' @export
#' @examples
#' \dontrun{
#' cf_setup_method("dowhy")
#' }
cf_setup_method <- function(method = c("dowhy", "ganite", "econml", "pyro"),
                            force = FALSE) {
  method <- match.arg(method)

  cat(sprintf("Setting up environment for method: %s\n", method))

  # Use existing function
  cf_install_python_env(method = method, force = force)
  cf_use_python_env(method = method)

  cat("Setup complete. Run cf_check_env() to verify.\n")
  invisible(TRUE)
}
