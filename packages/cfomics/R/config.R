#' Get or Set cfomics Configuration
#'
#' Access and modify cfomics configuration settings stored in
#' ~/.cfomics/config.yaml. Configuration persists across R sessions.
#'
#' @param key Configuration key using dot notation (e.g., "python.path",
#'   "methods.default"). If NULL, returns all configuration.
#' @param value Value to set. If NULL, returns current value for the key.
#' @return Current configuration value or invisibly TRUE after setting
#' @export
#' @examples
#' \dontrun{
#' cf_config()                              # Show all config
#' cf_config("methods.default")             # Get specific value
#' cf_config("methods.default", "grf")      # Set value
#' }
cf_config <- function(key = NULL, value = NULL) {
  config_path <- get_config_path()
  config <- load_config(config_path)

  if (is.null(key) && is.null(value)) {
    # Show all config
    print_config(config)
    return(invisible(config))
  }

  if (is.null(value)) {
    # Get value
    val <- get_nested_value(config, key)
    return(val)
  }

  # Set value
  config <- set_nested_value(config, key, value)
  save_config(config, config_path)
  message(sprintf("Set %s = %s", key, as.character(value)))
  invisible(TRUE)
}

#' Get configuration file path
#' @noRd
get_config_path <- function() {
  config_dir <- file.path(Sys.getenv("HOME"), ".cfomics")
  if (!dir.exists(config_dir)) {
    dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
  }
  file.path(config_dir, "config.yaml")
}

#' Get default configuration
#' @noRd
get_default_config <- function() {
  list(
    python = list(
      path = "auto",
      env_name = "cfomics_unified"
    ),
    methods = list(
      default = "grf",
      fallback_to_r = TRUE
    ),
    preferences = list(
      auto_install = FALSE,
      verbose = TRUE
    )
  )
}

#' Load configuration from file
#' @noRd
load_config <- function(path) {
  if (!file.exists(path)) {
    return(get_default_config())
  }

  if (!requireNamespace("yaml", quietly = TRUE)) {
    message("yaml package not installed. Using default configuration.")
    return(get_default_config())
  }

  tryCatch({
    yaml::read_yaml(path)
  }, error = function(e) {
    message(sprintf("Failed to read config: %s", e$message))
    get_default_config()
  })
}

#' Save configuration to file
#' @noRd
save_config <- function(config, path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    message("yaml package not installed. Cannot save configuration.")
    message("Install with: install.packages('yaml')")
    return(invisible(FALSE))
  }

  # Ensure directory exists
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  tryCatch({
    yaml::write_yaml(config, path)
    invisible(TRUE)
  }, error = function(e) {
    message(sprintf("Failed to save config: %s", e$message))
    invisible(FALSE)
  })
}

#' Get nested value from list
#' @noRd
get_nested_value <- function(lst, key) {
  keys <- strsplit(key, "\\.")[[1]]
  val <- lst
  for (k in keys) {
    if (is.null(val) || !is.list(val)) return(NULL)
    val <- val[[k]]
  }
  val
}

#' Set nested value in list
#' @noRd
set_nested_value <- function(lst, key, value) {
  keys <- strsplit(key, "\\.")[[1]]

  if (length(keys) == 1) {
    lst[[keys[1]]] <- value
  } else {
    if (is.null(lst[[keys[1]]])) lst[[keys[1]]] <- list()
    lst[[keys[1]]] <- set_nested_value(lst[[keys[1]]],
                                        paste(keys[-1], collapse = "."),
                                        value)
  }
  lst
}

#' Print configuration
#' @noRd
print_config <- function(config) {
  cat("\n")
  cat("=======================================================\n")
  cat("               cfomics Configuration                   \n")
  cat("=======================================================\n")
  cat(sprintf("Config file: %s\n", get_config_path()))
  cat("\n")

  print_nested <- function(lst, indent = 0) {
    prefix <- paste(rep("  ", indent), collapse = "")
    for (name in names(lst)) {
      val <- lst[[name]]
      if (is.list(val)) {
        cat(sprintf("%s%s:\n", prefix, name))
        print_nested(val, indent + 1)
      } else {
        cat(sprintf("%s%s: %s\n", prefix, name, as.character(val)))
      }
    }
  }

  print_nested(config)
  cat("\n")
  cat("Use cf_config('key', value) to modify settings.\n")
  cat("Use cf_config_reset() to restore defaults.\n")
  cat("=======================================================\n")
}

#' Reset Configuration to Defaults
#'
#' Removes the configuration file and restores default settings.
#'
#' @param confirm Logical, whether to ask for confirmation (default TRUE)
#' @return Invisibly returns TRUE on success
#' @export
#' @examples
#' \dontrun{
#' cf_config_reset()
#' }
cf_config_reset <- function(confirm = TRUE) {
  config_path <- get_config_path()

  if (!file.exists(config_path)) {
    message("No configuration file found. Already using defaults.")
    return(invisible(TRUE))
  }

  if (confirm && interactive()) {
    cat(sprintf("This will delete: %s\n", config_path))
    ans <- utils::askYesNo("Reset configuration to defaults?")
    if (!isTRUE(ans)) {
      message("Cancelled.")
      return(invisible(FALSE))
    }
  }

  tryCatch({
    file.remove(config_path)
    message("Configuration reset to defaults.")
    invisible(TRUE)
  }, error = function(e) {
    message(sprintf("Failed to reset config: %s", e$message))
    invisible(FALSE)
  })
}

#' Get Default Method
#'
#' Returns the default causal inference method from configuration,
#' or "grf" if not configured.
#'
#' @return Character string with the default method name
#' @noRd
get_default_method <- function() {
  method <- cf_config("methods.default")
  if (is.null(method)) "grf" else method
}

#' Check if R Fallback is Enabled
#'
#' Returns whether to fall back to R methods when Python is unavailable.
#'
#' @return Logical
#' @noRd
use_r_fallback <- function() {
  fallback <- cf_config("methods.fallback_to_r")
  if (is.null(fallback)) TRUE else as.logical(fallback)
}
