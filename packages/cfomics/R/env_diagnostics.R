#' Check cfomics Environment Status
#'
#' Performs comprehensive diagnostics of the cfomics environment including
#' R version, Python availability, required packages, and system resources.
#'
#' @param verbose Logical, whether to print detailed report (default TRUE)
#' @param methods Character vector of methods to check, or "all" (default)
#' @return Invisibly returns a list with diagnostic results
#' @export
#' @examples
#' \dontrun{
#' cf_check_env()
#' cf_check_env(methods = c("grf", "dowhy_gcm"))
#' }
cf_check_env <- function(verbose = TRUE, methods = "all") {
  checks <- list(
    timestamp = Sys.time(),
    r_version = check_r_version(),
    os = check_os_info(),
    python = check_python_status(),
    conda = check_conda_status(),
    packages = list(
      r = check_r_packages(),
      python = check_python_packages()
    ),
    methods = check_method_availability(methods),
    resources = check_system_resources()
  )

  class(checks) <- c("cf_env_check", "list")

  if (verbose) {
    print(checks)
  }

  invisible(checks)
}

#' Check R version
#' @noRd
check_r_version <- function() {
  v <- getRversion()
  list(
    version = as.character(v),
    ok = v >= "4.2.0",
    message = if (v >= "4.2.0") "OK" else "R >= 4.2.0 required"
  )
}

#' Check OS information
#' @noRd
check_os_info <- function() {
  info <- Sys.info()
  list(
    sysname = info["sysname"],
    release = info["release"],
    machine = info["machine"]
  )
}

#' Check Python status
#' @noRd
check_python_status <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(list(available = FALSE, message = "reticulate not installed"))
  }

  available <- tryCatch({
    reticulate::py_available(initialize = TRUE)
  }, error = function(e) FALSE)

  if (!available) {
    return(list(available = FALSE, message = "Python not configured"))
  }

  config <- tryCatch({
    reticulate::py_config()
  }, error = function(e) NULL)

  if (is.null(config)) {
    return(list(available = TRUE, message = "Python available but config failed"))
  }

  list(
    available = TRUE,
    version = config$version,
    path = config$python,
    message = "OK"
  )
}

#' Check Conda status
#' @noRd
check_conda_status <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(list(available = FALSE, message = "reticulate not installed"))
  }

  conda <- tryCatch({
    reticulate::conda_binary()
  }, error = function(e) NULL)

  if (is.null(conda)) {
    return(list(available = FALSE, message = "Conda not found"))
  }

  envs <- tryCatch({
    reticulate::conda_list()
  }, error = function(e) data.frame(name = character(0), python = character(0)))

  cfomics_envs <- envs[grepl("cfomics", envs$name, ignore.case = TRUE), ]

  list(
    available = TRUE,
    binary = conda,
    total_environments = nrow(envs),
    cfomics_envs = if (nrow(cfomics_envs) > 0) cfomics_envs$name else character(0),
    message = "OK"
  )
}

#' Check R packages
#' @noRd
check_r_packages <- function() {
  required <- c("Formula", "igraph", "reticulate", "rlang", "cli", "glue")
  suggested <- c("grf", "ipw", "survey", "ggplot2",
                 "SummarizedExperiment", "MultiAssayExperiment")

  check_pkg <- function(pkg) {
    installed <- requireNamespace(pkg, quietly = TRUE)
    version <- if (installed) {
      tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA)
    } else {
      NA
    }
    list(installed = installed, version = version)
  }

  list(
    required = lapply(stats::setNames(required, required), check_pkg),
    suggested = lapply(stats::setNames(suggested, suggested), check_pkg)
  )
}

#' Check Python packages
#' @noRd
check_python_packages <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(list())
  }

  py_avail <- tryCatch({
    reticulate::py_available(initialize = FALSE)
  }, error = function(e) FALSE)

  if (!py_avail) {
    return(list())
  }

  packages <- c("numpy", "pandas", "scipy", "sklearn",
                "dowhy", "econml", "tensorflow", "torch", "pyro")

  check_py_pkg <- function(pkg) {
    tryCatch({
      reticulate::py_module_available(pkg)
    }, error = function(e) FALSE)
  }

  result <- lapply(stats::setNames(packages, packages), check_py_pkg)
  # sklearn is imported as sklearn but package is scikit-learn
  names(result)[names(result) == "sklearn"] <- "scikit-learn"
  result
}

#' Check method availability
#' @noRd
check_method_availability <- function(methods = "all") {
  all_methods <- c("grf", "ipw", "gformula", "dowhy_gcm",
                   "drlearner", "ganite", "cavae")

  if (identical(methods, "all")) {
    methods <- all_methods
  }

  r_only <- c("grf", "ipw", "gformula")

  check_method <- function(m) {
    if (m %in% r_only) {
      # R-only method
      pkg <- switch(m,
        grf = "grf",
        ipw = "ipw",
        gformula = NULL  # built-in
      )
      available <- is.null(pkg) || requireNamespace(pkg, quietly = TRUE)
      return(list(
        available = available,
        type = "R",
        message = if (available) "OK" else paste(pkg, "not installed")
      ))
    } else {
      # Python method
      py_avail <- tryCatch({
        reticulate::py_available(initialize = FALSE)
      }, error = function(e) FALSE)

      if (!py_avail) {
        return(list(
          available = FALSE,
          type = "Python",
          message = "Python not configured"
        ))
      }

      required_pkg <- switch(m,
        dowhy_gcm = "dowhy",
        drlearner = "econml",
        ganite = "tensorflow",
        cavae = "pyro"
      )

      available <- tryCatch({
        reticulate::py_module_available(required_pkg)
      }, error = function(e) FALSE)

      list(
        available = available,
        type = "Python",
        message = if (available) "OK" else paste(required_pkg, "not installed")
      )
    }
  }

  lapply(stats::setNames(methods, methods), check_method)
}

#' Check system resources
#' @noRd
check_system_resources <- function() {
  # Memory info (platform-dependent)
  mem <- tryCatch({
    if (Sys.info()["sysname"] == "Linux") {
      meminfo <- readLines("/proc/meminfo", n = 3)
      total <- as.numeric(gsub("[^0-9]", "", meminfo[1])) / 1024 / 1024
      available <- as.numeric(gsub("[^0-9]", "", meminfo[3])) / 1024 / 1024
      list(total_gb = round(total, 1), available_gb = round(available, 1))
    } else {
      list(total_gb = NA, available_gb = NA)
    }
  }, error = function(e) list(total_gb = NA, available_gb = NA))

  cores <- tryCatch({
    parallel::detectCores()
  }, error = function(e) NA)

  list(
    memory = mem,
    cores = cores
  )
}

#' Print method for cf_env_check
#' @param x cf_env_check object
#' @param ... Additional arguments (unused)
#' @export
print.cf_env_check <- function(x, ...) {
  # Header
  cat("\n")
  cat("=======================================================\n")
  cat("           cfomics Environment Diagnostic Report       \n")
  cat("=======================================================\n")
  cat(sprintf("Timestamp: %s\n", x$timestamp))
  cat("\n")

  # R Version
  r_icon <- if (x$r_version$ok) "[OK]" else "[!!]"
  cat(sprintf("%s R version: %s (%s)\n",
              r_icon, x$r_version$version, x$r_version$message))

  # OS
  cat(sprintf("[--] OS: %s %s (%s)\n",
              x$os$sysname, x$os$release, x$os$machine))

  # Python
  py_icon <- if (isTRUE(x$python$available)) "[OK]" else "[!!]"
  if (isTRUE(x$python$available)) {
    cat(sprintf("%s Python: %s (%s)\n",
                py_icon, x$python$version, x$python$path))
  } else {
    cat(sprintf("%s Python: %s\n", py_icon, x$python$message))
  }

  # Conda
  conda_icon <- if (isTRUE(x$conda$available)) "[OK]" else "[--]"
  if (isTRUE(x$conda$available)) {
    n_cfomics <- length(x$conda$cfomics_envs)
    cat(sprintf("%s Conda: %d cfomics environment(s) found\n",
                conda_icon, n_cfomics))
    if (n_cfomics > 0) {
      cat(sprintf("        Environments: %s\n",
                  paste(x$conda$cfomics_envs, collapse = ", ")))
    }
  } else {
    cat(sprintf("%s Conda: %s\n", conda_icon, x$conda$message))
  }

  cat("\n")
  cat("-------------------------------------------------------\n")
  cat("                    Method Availability                \n")
  cat("-------------------------------------------------------\n")

  for (m in names(x$methods)) {
    info <- x$methods[[m]]
    icon <- if (info$available) "[OK]" else "[!!]"
    cat(sprintf("%s %-12s [%s]: %s\n", icon, m, info$type, info$message))
  }

  cat("\n")
  cat("-------------------------------------------------------\n")
  cat("                    System Resources                   \n")
  cat("-------------------------------------------------------\n")

  if (!is.na(x$resources$memory$total_gb)) {
    cat(sprintf("Memory: %.1f GB available / %.1f GB total\n",
                x$resources$memory$available_gb,
                x$resources$memory$total_gb))
  } else {
    cat("Memory: Unable to determine\n")
  }

  if (!is.na(x$resources$cores)) {
    cat(sprintf("CPU cores: %d\n", x$resources$cores))
  }

  # Recommendations
  unavailable <- names(x$methods)[!vapply(x$methods, `[[`, logical(1), "available")]
  if (length(unavailable) > 0) {
    cat("\n")
    cat("-------------------------------------------------------\n")
    cat("                    Recommendations                    \n")
    cat("-------------------------------------------------------\n")
    cat("To enable unavailable methods, run:\n")
    cat("  cf_setup(level = 'full')\n")
    cat("\nOr install specific environments:\
")
    for (m in unavailable) {
      if (x$methods[[m]]$type == "Python") {
        env <- switch(m,
          dowhy_gcm = "dowhy",
          drlearner = "econml",
          ganite = "ganite",
          cavae = "pyro"
        )
        cat(sprintf("  cf_install_python_env('%s')
", env))
      } else {
        pkg <- switch(m, grf = "grf", ipw = "ipw", NULL)
        if (!is.null(pkg)) {
          cat(sprintf("  install.packages('%s')\n", pkg))
        }
      }
    }
  }

  cat("\n")
  cat("=======================================================\n")

  invisible(x)
}
