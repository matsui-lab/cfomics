test_that("cf_benchmark_generate_data returns correct structure", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 100L,
    p = 10L,
    seed = 123L
  )
  
  expect_type(result, "list")
  expect_true("data" %in% names(result))
  expect_true("truth" %in% names(result))
  expect_true("graph" %in% names(result))
  expect_true("meta" %in% names(result))
  
  expect_s3_class(result$data, "data.frame")
  expect_equal(nrow(result$data), 100)
  
  expect_true("Y" %in% names(result$data))
  expect_true("T" %in% names(result$data))
  expect_true("X1" %in% names(result$data))
  expect_true("X10" %in% names(result$data))
})

test_that("cf_benchmark_generate_data produces reproducible results with same seed", {
  result1 <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 50L,
    p = 5L,
    seed = 42L
  )
  
  result2 <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 50L,
    p = 5L,
    seed = 42L
  )
  
  expect_equal(result1$data$Y, result2$data$Y)
  expect_equal(result1$data$T, result2$data$T)
  expect_equal(result1$data$X1, result2$data$X1)
  expect_equal(result1$truth$ate_true, result2$truth$ate_true)
  expect_equal(result1$truth$ite_true, result2$truth$ite_true)
})

test_that("cf_benchmark_generate_data treatment is binary 0/1 integer", {
  for (scenario in c("linear_homogeneous", "nonlinear_outcome", 
                     "heterogeneous_ite", "strong_confounding")) {
    result <- cfomics:::cf_benchmark_generate_data(
      scenario = scenario,
      n = 100L,
      p = 5L,
      seed = 1L
    )
    
    T_vals <- result$data$T
    expect_type(T_vals, "integer")
    expect_true(all(T_vals %in% c(0L, 1L)))
  }
})

test_that("cf_benchmark_generate_data truth values have correct types", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "heterogeneous_ite",
    n = 100L,
    p = 5L,
    seed = 1L
  )
  
  expect_type(result$truth$ate_true, "double")
  expect_length(result$truth$ate_true, 1)
  
  expect_type(result$truth$ite_true, "double")
  expect_length(result$truth$ite_true, 100)
  
  expect_type(result$truth$mu0_true, "double")
  expect_length(result$truth$mu0_true, 100)
  
  expect_type(result$truth$mu1_true, "double")
  expect_length(result$truth$mu1_true, 100)
  
  expect_type(result$truth$propensity_true, "double")
  expect_length(result$truth$propensity_true, 100)
})

test_that("cf_benchmark_generate_data linear_homogeneous has constant ITE", {
  effect_size <- 2.5
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 100L,
    p = 5L,
    seed = 1L,
    effect_size = effect_size
  )
  
  expect_equal(result$truth$ate_true, effect_size)
  expect_true(all(result$truth$ite_true == effect_size))
})

test_that("cf_benchmark_generate_data heterogeneous_ite has varying ITE", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "heterogeneous_ite",
    n = 100L,
    p = 5L,
    seed = 1L,
    effect_size = 1.0
  )
  
  ite_sd <- sd(result$truth$ite_true)
  expect_true(ite_sd > 0)
  
  expect_equal(result$truth$ate_true, mean(result$truth$ite_true))
})

test_that("cf_benchmark_generate_data returns igraph when return_graph=TRUE", {
  skip_if_not_installed("igraph")
  
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 50L,
    p = 5L,
    seed = 1L,
    return_graph = TRUE
  )
  
  expect_s3_class(result$graph, "igraph")
  
  edges <- igraph::as_data_frame(result$graph, what = "edges")
  expect_true(any(edges$from == "X1" & edges$to == "T"))
  expect_true(any(edges$from == "T" & edges$to == "Y"))
  expect_true(any(edges$from == "X1" & edges$to == "Y"))
})

test_that("cf_benchmark_generate_data returns NULL graph when return_graph=FALSE", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "linear_homogeneous",
    n = 50L,
    p = 5L,
    seed = 1L,
    return_graph = FALSE
  )
  
  expect_null(result$graph)
})

test_that("cf_benchmark_generate_data meta contains correct values", {
  result <- cfomics:::cf_benchmark_generate_data(
    scenario = "strong_confounding",
    n = 200L,
    p = 15L,
    seed = 99L,
    effect_size = 3.0,
    noise_sd = 0.5
  )
  
  expect_equal(result$meta$scenario, "strong_confounding")
  expect_equal(result$meta$n, 200L)
  expect_equal(result$meta$p, 15L)
  expect_equal(result$meta$seed, 99L)
  expect_equal(result$meta$effect_size, 3.0)
  expect_equal(result$meta$noise_sd, 0.5)
})

test_that("cf_benchmark_generate_data all scenarios work", {
  scenarios <- c("linear_homogeneous", "nonlinear_outcome", 
                 "heterogeneous_ite", "strong_confounding")
  
  for (scenario in scenarios) {
    result <- cfomics:::cf_benchmark_generate_data(
      scenario = scenario,
      n = 50L,
      p = 5L,
      seed = 1L
    )
    
    expect_equal(nrow(result$data), 50)
    expect_true(is.numeric(result$truth$ate_true))
    expect_length(result$truth$ite_true, 50)
  }
})
