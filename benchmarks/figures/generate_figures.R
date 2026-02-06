# benchmarks/figures/generate_figures.R
# Generate all publication figures for Genome Biology paper
#
# Usage:
#   Rscript benchmarks/figures/generate_figures.R
#   Rscript benchmarks/figures/generate_figures.R --output-dir path/to/output
#
# Outputs:
#   - fig2_simulation.pdf: Main simulation benchmark results
#   - fig3_highdim.pdf: High-dimensional stability analysis
#   - fig4_tcga.pdf: TCGA semi-synthetic validation

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Determine script directory robustly
get_script_dir <- function() {
  # Method 1: commandArgs (works with Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }

  # Method 2: Check if benchmarks/figures exists relative to cwd
  if (dir.exists("benchmarks/figures")) {
    return("benchmarks/figures")
  }

  # Method 3: Check if we're in the figures directory
  if (file.exists("generate_figures.R")) {
    return(".")
  }

  stop("Cannot determine script directory. Run from project root or benchmarks/figures/ directory.")
}

script_dir <- get_script_dir()

# Publication-quality theme for Genome Biology
# Journal guidelines: 180mm width for full-page figures
theme_paper <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      # Remove minor grid lines for cleaner look
      panel.grid.minor = element_blank(),
      # Subtle background for facet strips
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = base_size),
      # Legend at bottom to save horizontal space
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 0, b = 0),
      # Consistent text sizing
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(size = base_size, color = "grey40"),
      # Panel spacing for faceted plots
      panel.spacing = unit(0.5, "lines")
    )
}

# Consistent method colors across all figures
METHOD_COLORS <- c(
  # cfomics native methods
  "gformula" = "#E69F00",  # Orange
  "hdml"     = "#56B4E9",  # Light blue
  "hdps"     = "#009E73",  # Green
  "tmle"     = "#0072B2",  # Dark blue
  "grf"      = "#D55E00",  # Red-orange
  "ipw"      = "#CC79A7",  # Pink

  # External methods (prefixed with ext_)
  "ext_matchit"  = "#F0E442",  # Yellow

  "ext_weightit" = "#999999",  # Grey
  "ext_hdm"      = "#882255",  # Wine
  "ext_doubleml" = "#44AA99",  # Teal
  "ext_bart"     = "#332288",  # Indigo
  "ext_tmle3"    = "#117733",  # Forest green

  # Neural network methods
  "tarnet"    = "#AA4499",  # Purple
  "cfrnet"    = "#88CCEE",  # Cyan
  "dragonnet" = "#DDCC77"   # Tan
)

# Source individual figure scripts
source(file.path(script_dir, "fig2_simulation.R"))
source(file.path(script_dir, "fig3_highdim.R"))
source(file.path(script_dir, "fig4_tcga.R"))


#' Load benchmark results from raw directory
#'
#' @param results_dir Path to benchmarks/results directory
#' @return Data frame with all benchmark results
load_benchmark_results <- function(results_dir) {
  raw_dir <- file.path(results_dir, "raw")

  if (!dir.exists(raw_dir)) {
    stop("Raw results directory not found: ", raw_dir,
         "\nRun the benchmark first with: Rscript benchmarks/run_full_benchmark.R")
  }

  files <- list.files(raw_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No result files found in ", raw_dir)
  }

  message(sprintf("Loading %d result files from %s...", length(files), raw_dir))

  results_list <- lapply(files, function(f) {
    tryCatch({
      res <- readRDS(f)
      if (!is.data.frame(res)) {
        warning("Skipping malformed file: ", basename(f))
        return(NULL)
      }
      res
    }, error = function(e) {
      warning("Failed to read file: ", basename(f), " - ", e$message)
      NULL
    })
  })

  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0) {
    stop("No valid result files found")
  }

  results <- do.call(rbind, results_list)
  message(sprintf("Loaded %d results (%d successful, %d errors)",
                  nrow(results),
                  sum(results$status == "ok", na.rm = TRUE),
                  sum(results$status == "error", na.rm = TRUE)))

  results
}


#' Load TCGA benchmark results
#'
#' @param results_dir Path to benchmarks/tcga/results directory
#' @return Data frame with TCGA benchmark results
load_tcga_results <- function(results_dir) {
  if (!dir.exists(results_dir)) {
    stop("TCGA results directory not found: ", results_dir,
         "\nRun the TCGA benchmark first with: Rscript benchmarks/tcga/run_tcga_benchmark.R")
  }

  files <- list.files(results_dir, pattern = "_results\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No TCGA result files found in ", results_dir)
  }

  message(sprintf("Loading %d TCGA result files from %s...", length(files), results_dir))

  results_list <- lapply(files, function(f) {
    tryCatch({
      res <- readRDS(f)
      if (!is.data.frame(res)) {
        warning("Skipping malformed file: ", basename(f))
        return(NULL)
      }
      res
    }, error = function(e) {
      warning("Failed to read file: ", basename(f), " - ", e$message)
      NULL
    })
  })

  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0) {
    stop("No valid TCGA result files found")
  }

  results <- do.call(rbind, results_list)
  message(sprintf("Loaded %d TCGA results", nrow(results)))

  results
}


#' Generate all publication figures
#'
#' @param results_dir Path to benchmarks/results directory
#' @param tcga_dir Path to benchmarks/tcga/results directory
#' @param output_dir Path to save generated figures
#' @param format Output format ("pdf" or "png")
#' @param dpi Resolution for PNG output (default 300)
#' @return Invisibly returns list of generated figure paths
generate_all_figures <- function(
    results_dir = file.path(dirname(script_dir), "results"),
    tcga_dir = file.path(dirname(script_dir), "tcga", "results"),
    output_dir = file.path(script_dir, "output"),
    format = "pdf",
    dpi = 300
) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Validate format
  format <- match.arg(format, c("pdf", "png"))
  ext <- paste0(".", format)

  generated_files <- list()

  # Load simulation results
  message("\n=== Loading Simulation Results ===")
  results <- tryCatch({
    load_benchmark_results(results_dir)
  }, error = function(e) {
    warning("Could not load simulation results: ", e$message)
    NULL
  })

  # Figure 2: Main simulation results
  if (!is.null(results)) {
    message("\n=== Generating Figure 2: Simulation Results ===")
    fig2 <- tryCatch({
      generate_fig2_simulation(results)
    }, error = function(e) {
      warning("Failed to generate Figure 2: ", e$message)
      NULL
    })

    if (!is.null(fig2)) {
      fig2_path <- file.path(output_dir, paste0("fig2_simulation", ext))
      ggsave(fig2_path, fig2,
             width = 180, height = 120, units = "mm", dpi = dpi)
      message(sprintf("Saved Figure 2 to: %s", fig2_path))
      generated_files$fig2 <- fig2_path
    }

    # Figure 3: High-dimensional stability
    message("\n=== Generating Figure 3: High-dimensional Stability ===")
    fig3 <- tryCatch({
      generate_fig3_highdim(results)
    }, error = function(e) {
      warning("Failed to generate Figure 3: ", e$message)
      NULL
    })

    if (!is.null(fig3)) {
      fig3_path <- file.path(output_dir, paste0("fig3_highdim", ext))
      ggsave(fig3_path, fig3,
             width = 180, height = 100, units = "mm", dpi = dpi)
      message(sprintf("Saved Figure 3 to: %s", fig3_path))
      generated_files$fig3 <- fig3_path
    }
  }

  # Figure 4: TCGA validation
  message("\n=== Loading TCGA Results ===")
  tcga_results <- tryCatch({
    load_tcga_results(tcga_dir)
  }, error = function(e) {
    warning("Could not load TCGA results: ", e$message)
    NULL
  })

  if (!is.null(tcga_results)) {
    message("\n=== Generating Figure 4: TCGA Validation ===")
    fig4 <- tryCatch({
      generate_fig4_tcga(tcga_results)
    }, error = function(e) {
      warning("Failed to generate Figure 4: ", e$message)
      NULL
    })

    if (!is.null(fig4)) {
      fig4_path <- file.path(output_dir, paste0("fig4_tcga", ext))
      ggsave(fig4_path, fig4,
             width = 180, height = 100, units = "mm", dpi = dpi)
      message(sprintf("Saved Figure 4 to: %s", fig4_path))
      generated_files$fig4 <- fig4_path
    }
  }

  # Summary
  message("\n=== Figure Generation Complete ===")
  if (length(generated_files) > 0) {
    message(sprintf("Generated %d figures in: %s", length(generated_files), output_dir))
    for (name in names(generated_files)) {
      message(sprintf("  - %s: %s", name, basename(generated_files[[name]])))
    }
  } else {
    warning("No figures were generated. Check that benchmark results exist.")
  }

  invisible(generated_files)
}


# Main execution
if (sys.nframe() == 0) {
  # Parse command-line arguments
  args <- commandArgs(trailingOnly = TRUE)

  # Handle --help flag
  if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: Rscript generate_figures.R [options]\n")
    cat("\n")
    cat("Options:\n")
    cat("  --output-dir DIR   Output directory for figures (default: benchmarks/figures/output)\n")
    cat("  --format FORMAT    Output format: pdf or png (default: pdf)\n")
    cat("  --dpi DPI          DPI for PNG output (default: 300)\n")
    cat("  --help, -h         Show this help message and exit\n")
    cat("\n")
    cat("Examples:\n")
    cat("  Rscript benchmarks/figures/generate_figures.R\n")
    cat("  Rscript benchmarks/figures/generate_figures.R --format png --dpi 600\n")
    cat("  Rscript benchmarks/figures/generate_figures.R --output-dir manuscript/figures\n")
    quit(status = 0, save = "no")
  }

  # Parse --output-dir
  output_dir <- NULL
  output_idx <- which(args == "--output-dir")
  if (length(output_idx) > 0 && output_idx < length(args)) {
    output_dir <- args[output_idx + 1]
  }

  # Parse --format
  format <- "pdf"
  format_idx <- which(args == "--format")
  if (length(format_idx) > 0 && format_idx < length(args)) {
    format <- args[format_idx + 1]
  }

  # Parse --dpi
  dpi <- 300
  dpi_idx <- which(args == "--dpi")
  if (length(dpi_idx) > 0 && dpi_idx < length(args)) {
    dpi <- as.integer(args[dpi_idx + 1])
  }

  # Run generation
  if (is.null(output_dir)) {
    generate_all_figures(format = format, dpi = dpi)
  } else {
    generate_all_figures(output_dir = output_dir, format = format, dpi = dpi)
  }
}
