#!/usr/bin/env Rscript
# Check all packages in packages/

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

packages <- list.dirs("packages", recursive = FALSE, full.names = TRUE)

if (length(packages) == 0) {
  message("No packages found in packages/")
  quit(status = 1)
}
failed <- character()

for (pkg in packages) {
  pkg_name <- basename(pkg)

  # Skip if not a valid R package
  if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
    message("Skipping ", pkg_name, " (no DESCRIPTION)")
    next
  }

  message("\n", strrep("=", 60))
  message("Checking: ", pkg_name)
  message(strrep("=", 60), "\n")

  tryCatch({
    if (quick) {
      # Quick check: no vignettes, no manual
      result <- devtools::check(
        pkg,
        document = FALSE,
        build_args = c("--no-build-vignettes", "--no-manual"),
        args = c("--no-vignettes", "--no-manual"),
        quiet = FALSE
      )
    } else {
      result <- devtools::check(pkg, quiet = FALSE)
    }

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
  message("All packages passed!")
  quit(status = 0)
}
