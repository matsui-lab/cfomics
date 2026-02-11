#' Fit CAVAE (Causal Effect Variational Autoencoder) model using Pyro
#'
#' @param X Covariate matrix
#' @param T Treatment vector (binary)
#' @param Y Outcome vector (binary for bernoulli, continuous for normal)
#' @param covariate_names Character vector of covariate names
#' @param random_state Integer random seed
#' @param outcome_dist Outcome distribution: "bernoulli" or "normal" (default: "bernoulli")
#' @param latent_dim Integer, latent dimension (default: 20)
#' @param hidden_dim Integer, hidden layer dimension (default: 200)
#' @param num_layers Integer, number of hidden layers (default: 3)
#' @param num_samples Integer, number of samples for ITE estimation (default: 10)
#' @param num_epochs Integer, number of training epochs (default: 30)
#' @param batch_size Integer, batch size (default: 100)
#' @param learning_rate Numeric, learning rate (default: 1e-3)
#' @param ... Additional arguments passed to CEVAE
#' @return List with model fit and results
#' @keywords internal
cf_fit_cavae <- function(X, T, Y,
                         covariate_names = NULL,
                         random_state = 0L,
                         outcome_dist = "bernoulli",
                         latent_dim = 20L,
                         hidden_dim = 200L,
                         num_layers = 3L,
                         num_samples = 10L,
                         num_epochs = 30L,
                         batch_size = 100L,
                         learning_rate = 1e-3,
                         ...) {
  # Ensure Python environment is properly configured
  cf_require_python("cavae")
  
  torch <- reticulate::import("torch", convert = FALSE)
  pyro <- reticulate::import("pyro", convert = FALSE)
  cevae_mod <- reticulate::import("pyro.contrib.cevae", convert = FALSE)
  CEVAE <- cevae_mod$CEVAE
  
  if (!is.null(random_state)) {
    pyro$set_rng_seed(as.integer(random_state))
  }
  
  X_mat <- as.matrix(X)
  feature_dim <- ncol(X_mat)
  n <- nrow(X_mat)
  
  x_tensor <- torch$tensor(X_mat, dtype = torch$float32)
  t_tensor <- torch$tensor(as.numeric(T), dtype = torch$float32)
  y_tensor <- torch$tensor(as.numeric(Y), dtype = torch$float32)
  
  cevae <- CEVAE(
    feature_dim = as.integer(feature_dim),
    outcome_dist = outcome_dist,
    latent_dim = as.integer(latent_dim),
    hidden_dim = as.integer(hidden_dim),
    num_layers = as.integer(num_layers),
    num_samples = as.integer(num_samples)
  )

  # Wrap Python fitting call with user-friendly error handling
  .wrap_python_call(
    cevae$fit(
      x_tensor,
      t_tensor,
      y_tensor,
      num_epochs = as.integer(num_epochs),
      batch_size = as.integer(batch_size),
      learning_rate = learning_rate,
      learning_rate_decay = 0.1,
      weight_decay = 1e-4
    ),
    method_name = "CAVAE"
  )

  # Wrap Python ITE estimation with user-friendly error handling
  ite_py <- .wrap_python_call(
    cevae$ite(x_tensor),
    method_name = "CAVAE"
  )
  
  # Ensure proper type conversion from torch tensor to R numeric vector
  # First detach from GPU, convert to numpy, then to R
  np <- reticulate::import("numpy", convert = TRUE)
  ite_np <- ite_py$detach()$cpu()$numpy()
  ite <- as.numeric(ite_np)
  
  # Ensure T and Y are also R numeric vectors for arithmetic operations
  T_num <- as.numeric(T)
  Y_num <- as.numeric(Y)
  
  ate <- mean(ite)
  
  y0_hat <- Y_num - T_num * ite
  y1_hat <- Y_num + (1 - T_num) * ite
  
  ite_summary <- list(
    mean = mean(ite),
    std = sd(ite),
    quantile_05 = quantile(ite, 0.05),
    quantile_50 = quantile(ite, 0.50),
    quantile_95 = quantile(ite, 0.95)
  )
  
  summary_stats <- list(
    ate = ate,
    ate_ci_lower = NA,
    ate_ci_upper = NA,
    ite_mean = ite_summary$mean,
    ite_std = ite_summary$std,
    ite_quantiles = list(
      q05 = ite_summary$quantile_05,
      q50 = ite_summary$quantile_50,
      q95 = ite_summary$quantile_95
    )
  )
  
  list(
    model = cevae,
    res = list(
      ite = ite,
      ate = ate,
      y0_hat = y0_hat,
      y1_hat = y1_hat,
      summary = summary_stats,
      samples = list(
        y0 = y0_hat,
        y1 = y1_hat,
        ite = ite
      ),
      metadata = NULL
    )
  )
}

#' Predict from CAVAE model
#'
#' @param object cf_model object with method="cavae"
#' @param newdata Optional new data for prediction
#' @param type Type of prediction: "ate", "ite", "y0", "y1"
#' @param ... Additional arguments
#' @return Predicted values
#' @keywords internal
predict_cf_cavae <- function(object, newdata = NULL, type = "ite", ...) {
  if (is.null(newdata)) {
    res <- object$fit$res
    return(switch(
      type,
      "ate" = res$ate,
      "ite" = res$ite,
      "y0" = res$y0_hat,
      "y1" = res$y1_hat,
      "summary" = res$summary,
      "samples" = res$samples
    ))
  }
  
  torch <- reticulate::import("torch", convert = FALSE)
  X_new <- torch$tensor(as.matrix(newdata), dtype = torch$float32)
  
  cevae <- object$fit$model
  ite_py <- cevae$ite(X_new)
  
  # Ensure proper type conversion from torch tensor to R numeric vector
  ite_np <- ite_py$detach()$cpu()$numpy()
  ite <- as.numeric(ite_np)
  
  switch(
    type,
    "ate" = mean(ite),
    "ite" = ite,
    stop(sprintf("type '%s' not supported for newdata prediction", type))
  )
}
