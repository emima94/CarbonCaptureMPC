get_state_covariance <- function(model, pred_state_i) {
  # Get state names from the model's systems (not from pred_state_i directly)
  state_names <- names(model$getSystems())
  
  # Number of states
  n_states <- length(state_names)
  
  # Initialize the state covariance matrix
  P <- matrix(NA_real_, nrow = n_states, ncol = n_states,
              dimnames = list(state_names, state_names))
  
  # Build the covariance matrix using the pred_state_i (variance and covariance)
  for (i in 1:n_states) {
    for (j in 1:n_states) {
      cov_name <- paste0("cov.", state_names[i], ".", state_names[j])
      var_name <- paste0("var.", state_names[i])
      
      if (i == j) {
        # Diagonal: variance of state i
        P[i, j] <- pred_state_i[[var_name]]
      } else {
        # Off-diagonal: covariance between states i and j
        P[i, j] <- pred_state_i[[cov_name]]
      }
    }
  }
  
  return(P)
}
