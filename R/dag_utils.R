
dag_to_edges_list <- function(graph) {
  if (!inherits(graph, "igraph")) {
      stop("graph must be an igraph object")
  }
  edges_df <- igraph::as_data_frame(graph, what = "edges")[, 1:2]

  # Convert to list of lists (or list of vectors) for reticulate
  # Reticulate converts list(c("A", "B"), c("B", "C")) to [("A", "B"), ("B", "C")] in Python
  # which works for networkx.

  # apply returns a matrix if result length is constant, or list if not.
  # Here result is length 2 char vector. apply might return matrix 2 x N.
  # We want a list of vectors.

  # Using simple loop or split to be safe with reticulate conversion.
  # Actually `apply(..., 1, function(x) x)` returns a matrix (cols=samples) if simplified.
  # We want a list.

  edges_list <- lapply(seq_len(nrow(edges_df)), function(i) {
      as.character(edges_df[i, ])
  })

  edges_list
}
