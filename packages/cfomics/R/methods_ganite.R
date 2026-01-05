#' @title GANITE Method Implementation
#' @description Internal functions for the GANITE causal inference method.
#'
#' GANITE (Generative Adversarial Nets for Inference of Individualized Treatment Effects)
#' uses GANs to estimate individual treatment effects.
#'
#' Reference: Jinsung Yoon, James Jordon, Mihaela van der Schaar,
#' "GANITE: Estimation of Individualized Treatment Effects using Generative Adversarial Nets",
#' International Conference on Learning Representations (ICLR), 2018.

#' Get GANITE Python module
#' @keywords internal
.get_ganite_module <- function() {
  reticulate::import_from_path(
    module = "ganite",
    path = system.file("python", package = "cfomics")
  )
}

#' Fit GANITE model
#'
#' Internal function to fit a causal model using GANITE (GAN-based ITE estimation).
#'
#' @param X Numeric matrix of covariates.
#' @param T Binary treatment vector.
#' @param Y Numeric outcome vector.
#' @param dag_spec List of edge pairs defining the DAG (kept for API consistency).
#' @param covariate_names Character vector of covariate names.
#' @param random_state Integer random seed for reproducibility.
#' @param h_dim Integer, hidden layer dimension (default: 100).
#' @param alpha Numeric, GAN loss weight (default: 1.0).
#' @param iterations Integer, number of training iterations (default: 10000).
#' @param batch_size Integer, mini-batch size (default: 256).
#' @param verbose Logical, whether to print training progress (default: FALSE).
#' @param ... Additional arguments passed to Python.
#' @return A cf_model_ganite object.
#' @keywords internal
cf_fit_ganite <- function(X, T, Y, dag_spec = NULL,
                          covariate_names = NULL,
                          random_state = 0L,
                          h_dim = 100L,
                          alpha = 1.0,
                          iterations = 10000L,
                          batch_size = 256L,
                          verbose = FALSE,
                          ...) {
  # Ensure Python environment is properly configured
  cf_require_python("ganite")

  mod <- .get_ganite_module()

  res <- mod$run_ganite(
    X = X,
    T = as.numeric(T),
    Y = as.numeric(Y),
    edges_list = dag_spec,
    variable_names = covariate_names,
    random_state = as.integer(random_state),
    h_dim = as.integer(h_dim),
    alpha = as.numeric(alpha),
    iterations = as.integer(iterations),
    batch_size = as.integer(batch_size),
    verbose = verbose,
    ...
  )

  structure(
    list(
      backend = "ganite",
      res = res
    ),
    class = "cf_model_ganite"
  )
}

#' Predict from GANITE model
#'
#' @param object cf_model object with method="ganite"
#' @param type Type of prediction: "ate", "ite", "y0", "y1", "summary", "samples"
#' @param ... Additional arguments (unused)
#' @return Predicted values
#' @keywords internal
predict_cf_ganite <- function(object,
                              type = c("ate", "ite", "y0", "y1", "summary", "samples"),
                              ...) {
  type <- match.arg(type)
  res <- object$fit$res

  switch(
    type,
    "ate" = res$ate,
    "ite" = res$ite,
    "y0" = res$y0_hat,
    "y1" = res$y1_hat,
    "summary" = res$summary,
    "samples" = res$samples
  )
}
