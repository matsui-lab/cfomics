#!/usr/bin/env Rscript
#
# Benchmark Runner for cfomics Causal Inference Methods
#
# This script runs systematic benchmarks comparing multiple causal inference
# methods across different data generating scenarios.
#
# Usage:
#   Rscript cfomics/inst/benchmarks/run_benchmark.R
#
# Output:
#   - benchmark_results/raw_results.rds
#   - benchmark_results/summary.csv
#   - benchmark_results/plots/*.pdf (if ggplot2 is available)
#

# ==============================================================================
# Configuration (edit these settings as needed)
# ==============================================================================

scenarios <- c(
  "linear_homogeneous",
  "nonlinear_outcome",
  "heterogeneous_ite",
  "strong_confounding"
)

methods <- c(
  "gformula",
  "ipw",
  "grf",
  "dowhy_gcm",
  "drlearner",
  "cavae"
)

n <- 300L
p <- 20L
n_reps <- 10L
base_seed <- 42L
effect_size <- 1.0
noise_sd <- 1.0
covariates_k <- 10L

allow_python_env_install <- FALSE

out_dir <- file.path(getwd(), "benchmark_results")

# ==============================================================================
# Setup
# ==============================================================================

message("cfomics Benchmark Runner")
message("========================")
message(sprintf("Output directory: %s", out_dir))
message(sprintf("Scenarios: %s", paste(scenarios, collapse = ", ")))
message(sprintf("Methods: %s", paste(methods, collapse = ", ")))
message(sprintf("n=%d, p=%d, n_reps=%d, base_seed=%d", n, p, n_reps, base_seed))
message(sprintf("allow_python_env_install: %s", allow_python_env_install))
message("")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

use_devtools <- FALSE

pkg_path <- file.path(getwd(), "cfomics")
if (!dir.exists(pkg_path)) {
  pkg_path <- file.path(dirname(dirname(dirname(getwd()))), "cfomics")
}

if (dir.exists(pkg_path) && file.exists(file.path(pkg_path, "R", "benchmark_runner.R"))) {
  message(sprintf("Loading cfomics from source: %s", pkg_path))
  devtools::load_all(pkg_path, export_all = TRUE)
  use_devtools <- TRUE
} else if (requireNamespace("cfomics", quietly = TRUE)) {
  library(cfomics)
} else {
  stop("cfomics package not found. Please install it or run from the repo root.")
}

# ==============================================================================
# Run Benchmark
# ==============================================================================

message("Starting benchmark runs...")
start_time <- Sys.time()

get_benchmark_fn <- function(fn_name) {
  if (use_devtools && exists(fn_name, envir = globalenv())) {
    return(get(fn_name, envir = globalenv()))
  }
  tryCatch({
    return(get(fn_name, envir = asNamespace("cfomics")))
  }, error = function(e) {
    stop(sprintf("Function '%s' not found. Make sure cfomics is properly installed or run from repo root.", fn_name))
  })
}

results_df <- get_benchmark_fn("cf_benchmark_run")(
  scenarios = scenarios,
  methods = methods,
  n = n,
  p = p,
  n_reps = n_reps,
  base_seed = base_seed,
  effect_size = effect_size,
  noise_sd = noise_sd,
  covariates_k = covariates_k,
  allow_python_env_install = allow_python_env_install
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
message(sprintf("Benchmark completed in %.2f minutes", as.numeric(elapsed)))

# ==============================================================================
# Save Raw Results
# ==============================================================================

raw_results_path <- file.path(out_dir, "raw_results.rds")
saveRDS(results_df, raw_results_path)
message(sprintf("Saved raw results: %s", raw_results_path))

# ==============================================================================
# Summarize Results
# ==============================================================================

message("Summarizing results...")
summary_df <- get_benchmark_fn("cf_benchmark_summarize")(results_df)

summary_path <- file.path(out_dir, "summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
message(sprintf("Saved summary: %s", summary_path))

message("")
message("Summary Statistics:")
message("-------------------")
print(summary_df)

# ==============================================================================
# Generate Plots (optional)
# ==============================================================================

message("")
message("Generating plots...")
plots_dir <- file.path(out_dir, "plots")
get_benchmark_fn("cf_benchmark_plot")(summary_df, out_dir = plots_dir)

# ==============================================================================
# Session Info
# ==============================================================================

session_info_path <- file.path(out_dir, "session_info.txt")
sink(session_info_path)
cat("Benchmark Session Information\n")
cat("=============================\n\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Platform: %s\n", R.version$platform))
cat("\n")
cat("Configuration:\n")
cat(sprintf("  scenarios: %s\n", paste(scenarios, collapse = ", ")))
cat(sprintf("  methods: %s\n", paste(methods, collapse = ", ")))
cat(sprintf("  n: %d\n", n))
cat(sprintf("  p: %d\n", p))
cat(sprintf("  n_reps: %d\n", n_reps))
cat(sprintf("  base_seed: %d\n", base_seed))
cat(sprintf("  effect_size: %.2f\n", effect_size))
cat(sprintf("  noise_sd: %.2f\n", noise_sd))
cat(sprintf("  covariates_k: %d\n", covariates_k))
cat(sprintf("  allow_python_env_install: %s\n", allow_python_env_install))
cat("\n")
cat("Elapsed time: ")
cat(sprintf("%.2f minutes\n", as.numeric(elapsed)))
cat("\n")
cat("Session Info:\n")
print(sessionInfo())
sink()
message(sprintf("Saved session info: %s", session_info_path))

message("")
message("Benchmark complete!")
message(sprintf("Results saved to: %s", out_dir))
