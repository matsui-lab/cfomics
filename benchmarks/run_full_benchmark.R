#!/usr/bin/env Rscript

# benchmarks/run_full_benchmark.R
# Main entry point for running cfomics benchmark simulations
#
# Usage:
#   Rscript benchmarks/run_full_benchmark.R           # Full production run
#   Rscript benchmarks/run_full_benchmark.R --smoke   # Quick smoke test
#   Rscript benchmarks/run_full_benchmark.R --retry   # Retry failed jobs
#   Rscript benchmarks/run_full_benchmark.R --smoke --retry

# Determine script directory robustly
get_script_dir <- function() {
  # Method 1: commandArgs (works with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }

  # Method 2: Check if benchmarks/config.R exists relative to cwd
  if (file.exists("benchmarks/config.R")) {
    return("benchmarks")
  }

  # Method 3: Check if config.R exists in current directory (running from benchmarks/)

  if (file.exists("config.R")) {
    return(".")
  }

  stop("Cannot determine script directory. Run from project root or benchmarks/ directory.")
}

script_dir <- get_script_dir()

# Load benchmark infrastructure
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "R", "scenarios.R"))
source(file.path(script_dir, "R", "runner.R"))

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Handle --help flag
if ("--help" %in% args) {
  cat("Usage: Rscript run_full_benchmark.R [options]\n")
  cat("\n")
  cat("Options:\n")
  cat("  --smoke   Run quick smoke test (3 scenarios, 2 reps)\n")
  cat("  --retry   Retry previously failed jobs\n")
  cat("  --help    Show this help message and exit\n")
  cat("\n")
  cat("Examples:\n")
  cat("  Rscript benchmarks/run_full_benchmark.R           # Full production run\n")
  cat("  Rscript benchmarks/run_full_benchmark.R --smoke   # Quick smoke test\n")
  cat("  Rscript benchmarks/run_full_benchmark.R --retry   # Retry failed jobs\n")
  cat("  Rscript benchmarks/run_full_benchmark.R --smoke --retry\n")
  quit(status = 0, save = "no")
}

# Validate arguments - warn on unknown flags
known_args <- c("--smoke", "--retry", "--help")
unknown <- setdiff(args, known_args)
if (length(unknown) > 0) {
  warning("Unknown arguments ignored: ", paste(unknown, collapse = ", "),
          "\nUse --help for usage information.")
}

smoke <- "--smoke" %in% args
retry <- "--retry" %in% args

# Select configuration
if (smoke) {
  message("Running smoke test...")
  cfg <- benchmark_config_smoke()
} else {
  message("Running full production benchmark...")
  cfg <- benchmark_config()
}

# Display configuration summary
message(sprintf("  Methods: %s", paste(cfg$methods, collapse = ", ")))
message(sprintf("  Scenarios: %d", length(cfg$scenarios)))
message(sprintf("  Replications: %d", cfg$n_reps))
message(sprintf("  Total jobs: %d", length(cfg$scenarios) * length(cfg$methods) * cfg$n_reps))
message(sprintf("  Results dir: %s", cfg$results_dir))
if (retry) {
  message("  Retry mode: ON (will re-run failed jobs)")
}
message("")

# Run benchmark
results <- run_benchmark(cfg, retry_errors = retry)

# Print summary
message("")
message("=== Benchmark Complete ===")
message(sprintf("Total jobs: %d", results$total))
message(sprintf("Previously completed: %d", results$completed))
message(sprintf("Executed this run: %d", results$executed))

# Count successes and failures from this run
if (results$executed > 0) {
  statuses <- vapply(results$results, function(r) r$status, character(1))
  n_ok <- sum(statuses == "ok")
  n_err <- sum(statuses == "error")
  message(sprintf("  Succeeded: %d", n_ok))
  message(sprintf("  Failed: %d", n_err))

  if (n_err > 0) {
    message("\nFailed jobs:")
    for (r in results$results) {
      if (r$status == "error") {
        message(sprintf("  - %s / %s / rep%d: %s",
                        r$scenario_id, r$method, r$rep, r$error_msg))
      }
    }
    message("\nRe-run with --retry to attempt failed jobs again.")
  }
}

message("")
message(sprintf("Results saved to: %s", file.path(cfg$results_dir, "raw")))

# Exit with non-zero status if any jobs failed
if (results$executed > 0 && n_err > 0) {
  quit(status = 1, save = "no")
}
