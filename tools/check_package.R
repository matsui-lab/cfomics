#!/usr/bin/env Rscript
#' Package Quality Check Script
#'
#' Run this script to perform comprehensive package checks locally.
#' Usage: Rscript tools/check_package.R [--quick|--full|--bioc]
#'
#' Options:
#'   --quick  Quick check (no vignettes, no manual)
#'   --full   Full CRAN-style check (default)
#'   --bioc   Run BiocCheck in addition to R CMD check

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[1] else "--full"

cat("=======================================================\n")
cat("           cfomics Package Quality Check               \n")
cat("=======================================================\n")
cat(sprintf("Mode: %s\n", mode))
cat(sprintf("Time: %s\n", Sys.time()))
cat("\n")

# Ensure we're in the package directory
pkg_dir <- getwd()
if (!file.exists("DESCRIPTION")) {
  # Try parent directory
  if (file.exists("cfomics/DESCRIPTION")) {
    pkg_dir <- file.path(getwd(), "cfomics")
  } else {
    stop("Cannot find DESCRIPTION file. Run from package root or repository root.")
  }
}

cat(sprintf("Package directory: %s\n\n", pkg_dir))

# Step 1: Document
cat("[1/4] Generating documentation...\n")
tryCatch({
  devtools::document(pkg_dir)
  cat("      Documentation generated successfully.\n")
}, error = function(e) {
  cat(sprintf("      Warning: %s\n", e$message))
})

# Step 2: Build
cat("\n[2/4] Building package...\n")
build_args <- character(0)
if (mode == "--quick") {
  build_args <- c("--no-build-vignettes", "--no-manual")
}

pkg_file <- tryCatch({
  devtools::build(pkg_dir, args = build_args)
}, error = function(e) {
  stop(sprintf("Build failed: %s", e$message))
})
cat(sprintf("      Built: %s\n", basename(pkg_file)))

# Step 3: Check
cat("\n[3/4] Running R CMD check...\n")
check_args <- c("--as-cran")
if (mode == "--quick") {
  check_args <- c(check_args, "--no-vignettes", "--no-manual", "--no-examples")
}

check_result <- tryCatch({
  rcmdcheck::rcmdcheck(
    pkg_file,
    args = check_args,
    error_on = "never"
  )
}, error = function(e) {
  stop(sprintf("Check failed: %s", e$message))
})

# Print summary
cat("\n")
cat("-------------------------------------------------------\n")
cat("                    Check Summary                      \n")
cat("-------------------------------------------------------\n")

n_errors <- length(check_result$errors)
n_warnings <- length(check_result$warnings)
n_notes <- length(check_result$notes)

cat(sprintf("Errors:   %d\n", n_errors))
cat(sprintf("Warnings: %d\n", n_warnings))
cat(sprintf("Notes:    %d\n", n_notes))

if (n_errors > 0) {
  cat("\nERRORS:\n")
  for (e in check_result$errors) {
    cat(sprintf("  - %s\n", substr(e, 1, 200)))
  }
}

if (n_warnings > 0) {
  cat("\nWARNINGS:\n")
  for (w in check_result$warnings) {
    cat(sprintf("  - %s\n", substr(w, 1, 200)))
  }
}

if (n_notes > 0) {
  cat("\nNOTES:\n")
  for (n in check_result$notes) {
    cat(sprintf("  - %s\n", substr(n, 1, 200)))
  }
}

# Step 4: BiocCheck (optional)
if (mode == "--bioc") {
  cat("\n[4/4] Running BiocCheck...\n")
  if (requireNamespace("BiocCheck", quietly = TRUE)) {
    tryCatch({
      BiocCheck::BiocCheck(pkg_dir, `new-package` = TRUE)
    }, error = function(e) {
      cat(sprintf("      BiocCheck error: %s\n", e$message))
    })
  } else {
    cat("      BiocCheck not installed. Install with:\n")
    cat("      BiocManager::install('BiocCheck')\n")
  }
} else {
  cat("\n[4/4] Skipping BiocCheck (use --bioc to enable)\n")
}

# Final status
cat("\n")
cat("=======================================================\n")
if (n_errors == 0 && n_warnings == 0) {
  cat("[OK] Package check passed!\n")
  if (n_notes > 0) {
    cat(sprintf("     (%d NOTE(s) - review if submitting to CRAN)\n", n_notes))
  }
} else {
  cat("[!!] Package check found issues.\n")
  cat("     Please fix errors and warnings before release.\n")
}
cat("=======================================================\n")

# Return exit code
quit(status = if (n_errors > 0) 1 else 0)
