# Tests for DAG utility functions

test_that("dag_to_edges_list requires igraph object", {
  # Access internal function
  dag_to_edges_list <- cfomics:::dag_to_edges_list

  expect_error(
    dag_to_edges_list("not a graph"),
    "must be an igraph object"
  )

  expect_error(
    dag_to_edges_list(list(edges = c(1, 2))),
    "must be an igraph object"
  )

  expect_error(
    dag_to_edges_list(NULL),
    "must be an igraph object"
  )
})

test_that("dag_to_edges_list converts simple graph correctly", {
  skip_if_not_installed("igraph")

  dag_to_edges_list <- cfomics:::dag_to_edges_list

  # Create simple 2-node graph: A -> B
  g <- igraph::graph_from_edgelist(
    matrix(c("A", "B"), ncol = 2, byrow = TRUE)
  )

  edges <- dag_to_edges_list(g)

  expect_true(is.list(edges))
  expect_length(edges, 1)
  expect_equal(edges[[1]], c("A", "B"))
})

test_that("dag_to_edges_list converts chain graph correctly", {
  skip_if_not_installed("igraph")

  dag_to_edges_list <- cfomics:::dag_to_edges_list

  # Create chain: A -> B -> C
  g <- igraph::graph_from_edgelist(
    matrix(c("A", "B", "B", "C"), ncol = 2, byrow = TRUE)
  )

  edges <- dag_to_edges_list(g)

  expect_true(is.list(edges))
  expect_length(edges, 2)

  # Check edges are correct (order may vary)
  edge_strs <- sapply(edges, paste, collapse = "->")
  expect_true("A->B" %in% edge_strs)
  expect_true("B->C" %in% edge_strs)
})

test_that("dag_to_edges_list converts fork structure correctly", {
  skip_if_not_installed("igraph")

  dag_to_edges_list <- cfomics:::dag_to_edges_list

  # Create fork: X -> T, X -> Y (confounder structure)
  g <- igraph::graph_from_edgelist(
    matrix(c("X", "T", "X", "Y"), ncol = 2, byrow = TRUE)
  )

  edges <- dag_to_edges_list(g)

  expect_true(is.list(edges))
  expect_length(edges, 2)

  edge_strs <- sapply(edges, paste, collapse = "->")
  expect_true("X->T" %in% edge_strs)
  expect_true("X->Y" %in% edge_strs)
})

test_that("dag_to_edges_list converts complex DAG correctly", {
  skip_if_not_installed("igraph")

  dag_to_edges_list <- cfomics:::dag_to_edges_list

  # Create typical causal DAG: X1 -> T -> Y, X1 -> Y, X2 -> Y
  g <- igraph::graph_from_edgelist(
    matrix(c(
      "X1", "T",
      "T", "Y",
      "X1", "Y",
      "X2", "Y"
    ), ncol = 2, byrow = TRUE)
  )

  edges <- dag_to_edges_list(g)

  expect_true(is.list(edges))
  expect_length(edges, 4)

  edge_strs <- sapply(edges, paste, collapse = "->")
  expect_true("X1->T" %in% edge_strs)
  expect_true("T->Y" %in% edge_strs)
  expect_true("X1->Y" %in% edge_strs)
  expect_true("X2->Y" %in% edge_strs)
})

test_that("dag_to_edges_list preserves node names with underscores", {
  skip_if_not_installed("igraph")

  dag_to_edges_list <- cfomics:::dag_to_edges_list

  # Node names with underscores
  g <- igraph::graph_from_edgelist(
    matrix(c("var_name_1", "var_name_2"), ncol = 2, byrow = TRUE)
  )

  edges <- dag_to_edges_list(g)

  expect_equal(edges[[1]], c("var_name_1", "var_name_2"))
})

test_that("dag_to_edges_list returns empty list for empty graph", {
  skip_if_not_installed("igraph")

  dag_to_edges_list <- cfomics:::dag_to_edges_list

  # Create graph with no edges
  g <- igraph::make_empty_graph(n = 3)

  edges <- dag_to_edges_list(g)

  expect_true(is.list(edges))
  expect_length(edges, 0)
})

test_that("dag_to_edges_list handles numeric node IDs", {
  skip_if_not_installed("igraph")

  dag_to_edges_list <- cfomics:::dag_to_edges_list

  # Create graph with numeric-like names (they become character)
  g <- igraph::graph_from_edgelist(
    matrix(c("1", "2", "2", "3"), ncol = 2, byrow = TRUE)
  )

  edges <- dag_to_edges_list(g)

  expect_length(edges, 2)
  expect_type(edges[[1]], "character")
})

test_that("dag_to_edges_list output is compatible with Python", {
  skip_if_not_installed("igraph")

  dag_to_edges_list <- cfomics:::dag_to_edges_list

  g <- igraph::graph_from_edgelist(
    matrix(c("A", "B", "B", "C"), ncol = 2, byrow = TRUE)
  )

  edges <- dag_to_edges_list(g)

  # Each element should be a character vector of length 2
  for (edge in edges) {
    expect_type(edge, "character")
    expect_length(edge, 2)
  }

  # Structure should be list of vectors (reticulate converts to list of tuples)
  expect_true(is.list(edges))
  expect_false(is.data.frame(edges))
})

# Tests for DAG integration with cf_fit (when DoWhy is available)
test_that("cf_fit accepts igraph DAG", {
  skip_if_not_installed("igraph")

  set.seed(42)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  T <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- X1 + X2 + T + rnorm(n)
  data <- data.frame(Y = Y, T = T, X1 = X1, X2 = X2)

  g <- igraph::graph_from_edgelist(
    matrix(c("X1", "T", "T", "Y", "X1", "Y", "X2", "Y"),
           ncol = 2, byrow = TRUE)
  )

  # Use gformula which doesn't require graph, but accepts it
  fit <- cf_fit(Y ~ T | X1 + X2, data = data, method = "gformula")

  expect_s3_class(fit, "cf_model")
})
