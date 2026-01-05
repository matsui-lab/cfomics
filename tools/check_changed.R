#!/usr/bin/env Rscript
# Check only packages with changes (based on git diff)

# Get changed files compared to main branch or last commit
changed_files <- tryCatch({
  # Try comparing to main branch first
  files <- system("git diff --name-only origin/main...HEAD 2>/dev/null", intern = TRUE)
  if (length(files) == 0) {
    # Fall back to comparing with last commit
    files <- system("git diff --name-only HEAD~1 2>/dev/null", intern = TRUE)
  }
  if (length(files) == 0) {
    # Fall back to staged/unstaged changes
    files <- system("git diff --name-only HEAD", intern = TRUE)
  }
  files
}, error = function(e) {
  character(0)
})

if (length(changed_files) == 0) {
  message("No changed files detected")
  quit(status = 0)
}

# Extract package names from changed paths
pkg_pattern <- "^packages/([^/]+)/"
matches <- regmatches(changed_files, regexec(pkg_pattern, changed_files))
changed_packages <- unique(unlist(lapply(matches, function(m) if (length(m) > 1) m[2] else NULL)))

if (length(changed_packages) == 0) {
  message("No package changes detected")
  message("Changed files: ", paste(changed_files, collapse = ", "))
  quit(status = 0)
}

message("Changed packages: ", paste(changed_packages, collapse = ", "))

failed <- character()

for (pkg_name in changed_packages) {
  pkg_path <- file.path("packages", pkg_name)

  if (!dir.exists(pkg_path)) {
    message("Package directory not found: ", pkg_path)
    next
  }

  if (!file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    message("Skipping ", pkg_name, " (no DESCRIPTION)")
    next
  }

  message("\n", strrep("=", 60))
  message("Checking: ", pkg_name)
  message(strrep("=", 60), "\n")

  tryCatch({
    result <- devtools::check(pkg_path, quiet = FALSE)
    if (length(result$errors) > 0) {
      failed <- c(failed, pkg_name)
    }
  }, error = function(e) {
    message("Error checking ", pkg_name, ": ", e$message)
    failed <<- c(failed, pkg_name)
  })
}

message("\n", strrep("=", 60))
if (length(failed) > 0) {
  message("FAILED packages: ", paste(failed, collapse = ", "))
  quit(status = 1)
} else {
  message("All changed packages passed!")
  quit(status = 0)
}
