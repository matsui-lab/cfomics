# packages/cfomics/data-raw/example_data.R
# Script to generate example datasets for cfomics

set.seed(42)

# Simple example data.frame
n <- 200
X1 <- rnorm(n)
X2 <- rnorm(n)
ps <- plogis(0.5 * X1 - 0.3 * X2)
T <- rbinom(n, 1, ps)
Y <- 2 * T + 0.8 * X1 + 0.4 * X2 + rnorm(n)

cfomics_example <- data.frame(
  Y = Y,
  T = T,
  X1 = X1,
  X2 = X2
)

# True ATE = 2
attr(cfomics_example, "true_ate") <- 2

usethis::use_data(cfomics_example, overwrite = TRUE)
