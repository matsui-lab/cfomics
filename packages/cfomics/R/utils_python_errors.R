#' @title Python Error Handling Utilities
#' @description Internal functions for wrapping Python calls with user-friendly
#' error messages.

#' Wrap Python calls with user-friendly error handling
#'
#' Executes a Python expression and translates any Python errors into
#' user-friendly R error messages with class "cfomics_python_error".
#'
#' @param expr Expression to evaluate (typically a Python call via reticulate)
#' @param method_name Character string with the method name for error messages
#'   (e.g., "DoWhy-GCM", "GANITE", "DRLearner", "CAVAE")
#'
#' @return The result of evaluating `expr` if successful
#'
#' @details
#' This function wraps Python calls to provide consistent, informative error

#' messages when Python code fails. It:
#' \itemize{
#'   \item Catches Python errors from reticulate
#'   \item Extracts the underlying Python error message
#'   \item Provides suggestions for common issues
#'   \item Preserves the original error as the parent condition
#' }
#'
#' @keywords internal
#' @noRd
.wrap_python_call <- function(expr, method_name) {
  tryCatch(
    expr,
    error = function(e) {
      msg <- conditionMessage(e)
      error_class <- class(e)

      # Check if it's a Python error (from reticulate)
      is_python_error <- any(grepl("python", tolower(error_class))) ||
        grepl("Error in py_call_impl", msg, fixed = TRUE) ||
        grepl("reticulate", msg, ignore.case = TRUE) ||
        inherits(e, "python.builtin.Exception")

      if (is_python_error) {
        # Extract clean Python error message
        clean_msg <- .extract_python_error_message(msg)

        # Build user-friendly message with suggestions
        rlang::abort(
          c(
            paste0(method_name, " fitting failed"),
            "x" = paste("Python error:", clean_msg),
            "i" = "Check that your data is numeric and properly formatted",
            "i" = "Ensure treatment is binary (0/1) and outcome is numeric",
            "i" = "Run cf_check_env() to verify your Python environment"
          ),
          class = "cfomics_python_error",
          parent = e,
          method = method_name
        )
      } else {
        # Re-throw non-Python errors with proper parent
        rlang::abort(
          msg,
          class = "cfomics_error",
          parent = e
        )
      }
    }
  )
}

#' Extract clean Python error message
#'
#' Parses a reticulate error message to extract the underlying Python
#' error message, removing R call stack noise.
#'
#' @param msg Character string with the full error message
#'
#' @return Character string with the cleaned Python error message
#'
#' @keywords internal
#' @noRd
.extract_python_error_message <- function(msg) {
  # Remove common R/reticulate prefixes
  clean <- gsub("Error in py_call_impl\\([^)]*\\)\\s*:\\s*", "", msg)

  # Remove Python traceback file paths if present
  clean <- gsub("File \"[^\"]+\", line \\d+, in [^\n]+\n?", "", clean)

  # Extract the actual error type and message if present
  # Pattern: ErrorType: message
  if (grepl("^[A-Z][a-zA-Z]+Error:", clean)) {
    # Keep the Python error type and message
    clean <- sub("^.*?([A-Z][a-zA-Z]+Error:.*)$", "\\1", clean)
  }

  # Trim whitespace and newlines

  clean <- trimws(clean)

  # Truncate very long messages

  if (nchar(clean) > 500) {
    clean <- paste0(substr(clean, 1, 497), "...")
  }

  # If we still have the original message, provide a generic fallback

  if (clean == msg || nchar(clean) < 5) {
    clean <- "An error occurred in the Python backend"
  }

  clean
}

#' Wrap Python module import with user-friendly error handling
#'
#' @param module_name Character string with the Python module name
#' @param method_name Character string with the cfomics method name
#'
#' @return The imported Python module
#'
#' @keywords internal
#' @noRd
.wrap_python_import <- function(module_name, method_name) {
  tryCatch(
    reticulate::import(module_name),
    error = function(e) {
      rlang::abort(
        c(
          paste0("Failed to import Python module '", module_name, "' for ", method_name),
          "x" = conditionMessage(e),
          "i" = paste0("Install the required package: pip install ", module_name),
          "i" = "Run cf_check_env() to verify your Python environment"
        ),
        class = "cfomics_python_import_error",
        parent = e,
        method = method_name,
        module = module_name
      )
    }
  )
}
