#!/usr/bin/env Rscript

# benchmarks/aggregate_results.R
# Collect raw RDS results and aggregate into summary CSV files
#
# Usage:
#   Rscript benchmarks/aggregate_results.R
#
# Outputs:
#   results/summary/raw_results.csv       - All individual results combined
#   results/summary/aggregated_summary.csv - Metrics aggregated by scenario-method

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

# Directories
raw_dir <- file.path(script_dir, "results", "raw")
summary_dir <- file.path(script_dir, "results", "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# Validate raw directory exists
if (!dir.exists(raw_dir)) {
  stop("Raw results directory does not exist: ", raw_dir,
       "\nRun the benchmark first with: Rscript benchmarks/run_full_benchmark.R")
}

# Read all RDS files
rds_files <- list.files(raw_dir, pattern = "\\.rds$", full.names = TRUE)
message(sprintf("Found %d result files in %s", length(rds_files), raw_dir))

if (length(rds_files) == 0) {
  stop("No RDS files found in ", raw_dir,
       "\nRun the benchmark first with: Rscript benchmarks/run_full_benchmark.R")
}

# Load and combine results
message("Loading result files...")
results_list <- lapply(rds_files, function(f) {
  tryCatch({
    res <- readRDS(f)
    # Ensure it's a data.frame with required columns
    if (!is.data.frame(res) || !("status" %in% names(res))) {
      warning("Skipping malformed file: ", basename(f))
      return(NULL)
    }
    res
  }, error = function(e) {
    warning("Failed to read file: ", basename(f), " - ", conditionMessage(e))
    NULL
  })
})

# Remove NULL entries (corrupted files)
results_list <- Filter(Negate(is.null), results_list)

if (length(results_list) == 0) {
  stop("No valid result files found. All files may be corrupted.")
}

# Combine into single data.frame
raw_df <- do.call(rbind, results_list)
message(sprintf("Combined %d total results (%d successful, %d errors)",
                nrow(raw_df),
                sum(raw_df$status == "ok"),
                sum(raw_df$status == "error")))

# Save raw combined results
raw_csv_path <- file.path(summary_dir, "raw_results.csv")
write.csv(raw_df, raw_csv_path, row.names = FALSE)
message(sprintf("Saved raw results to: %s", raw_csv_path))

# Filter to successful results for aggregation
ok_df <- raw_df[raw_df$status == "ok", ]

if (nrow(ok_df) == 0) {
  warning("No successful results to aggregate. All jobs failed.")
  message("\n=== Aggregation Complete (no successful results) ===")
  message(sprintf("Raw results: %s", raw_csv_path))
  quit(status = 0, save = "no")
}

# Aggregate by scenario_id and method
message("\nAggregating metrics by scenario and method...")

# Helper to compute safe statistics
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_sd <- function(x) if (all(is.na(x)) || sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
safe_sum <- function(x) if (all(is.na(x))) NA_integer_ else sum(!is.na(x))

# Group by scenario_id and method
groups <- split(ok_df, list(ok_df$scenario_id, ok_df$method), drop = TRUE)

agg_list <- lapply(names(groups), function(grp_name) {
  grp <- groups[[grp_name]]
  if (nrow(grp) == 0) return(NULL)

  # Extract scenario_id and method from group name
  # Format is "scenario_id.method"
  parts <- strsplit(grp_name, "\\.")[[1]]
  scenario_id <- parts[1]
  method <- paste(parts[-1], collapse = ".")  # Handle methods with dots if any

  # Compute RMSE from individual MSE values
  # RMSE = sqrt(mean(MSE)) since MSE = (estimate - true)^2
  rmse_ate <- sqrt(safe_mean(grp$mse_ate))

  data.frame(
    scenario_id = scenario_id,
    method = method,
    n_reps = nrow(grp),
    n = grp$n[1],
    p = grp$p[1],
    ate_true = grp$ate_true[1],

    # Bias metrics
    mean_bias = safe_mean(grp$bias_ate),
    sd_bias = safe_sd(grp$bias_ate),
    mean_abs_bias = safe_mean(grp$abs_bias_ate),

    # RMSE
    rmse_ate = rmse_ate,

    # PEHE (Precision in Estimating Heterogeneous Effects)
    mean_pehe = safe_mean(grp$pehe),
    sd_pehe = safe_sd(grp$pehe),

    # Coverage (proportion of CIs containing true ATE)
    coverage = safe_mean(grp$coverage_ate),

    # CI length
    mean_ci_length = safe_mean(grp$ci_len_ate),
    sd_ci_length = safe_sd(grp$ci_len_ate),

    # Computation time
    mean_time = safe_mean(grp$time_sec),
    sd_time = safe_sd(grp$time_sec),
    total_time = sum(grp$time_sec, na.rm = TRUE),

    stringsAsFactors = FALSE
  )
})

# Combine aggregated results
agg_df <- do.call(rbind, Filter(Negate(is.null), agg_list))

# Sort by scenario_id then method
agg_df <- agg_df[order(agg_df$scenario_id, agg_df$method), ]
rownames(agg_df) <- NULL

# Save aggregated summary
agg_csv_path <- file.path(summary_dir, "aggregated_summary.csv")
write.csv(agg_df, agg_csv_path, row.names = FALSE)
message(sprintf("Saved aggregated summary to: %s", agg_csv_path))

# Print summary statistics
message("\n=== Aggregation Summary ===")
message(sprintf("Total result files: %d", length(rds_files)))
message(sprintf("Valid results: %d", nrow(raw_df)))
message(sprintf("Successful results: %d", nrow(ok_df)))
message(sprintf("Failed results: %d", sum(raw_df$status == "error")))
message(sprintf("Unique scenarios: %d", length(unique(ok_df$scenario_id))))
message(sprintf("Unique methods: %d", length(unique(ok_df$method))))
message(sprintf("Scenario-method combinations: %d", nrow(agg_df)))

# Print method-level summary
message("\n--- Per-Method Summary ---")
method_summary <- aggregate(
  cbind(n_reps, rmse_ate, coverage, mean_time) ~ method,
  data = agg_df,
  FUN = function(x) round(mean(x, na.rm = TRUE), 4)
)
names(method_summary) <- c("Method", "Avg Reps", "Avg RMSE", "Avg Coverage", "Avg Time(s)")
print(method_summary, row.names = FALSE)

# Print scenario-level summary (abbreviated)
message("\n--- Top 10 Scenarios by RMSE (best methods) ---")
best_per_scenario <- do.call(rbind, lapply(split(agg_df, agg_df$scenario_id), function(x) {
  x[which.min(x$rmse_ate), c("scenario_id", "method", "rmse_ate", "coverage", "mean_time")]
}))
best_per_scenario <- best_per_scenario[order(best_per_scenario$rmse_ate), ]
print(head(best_per_scenario, 10), row.names = FALSE)

message("\n=== Aggregation Complete ===")
message(sprintf("Raw results: %s", raw_csv_path))
message(sprintf("Aggregated summary: %s", agg_csv_path))
