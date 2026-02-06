compute_obs_variance_covariance <- function(model, df, pred, par_free = NULL) {
  # Get number of timesteps in the prediction
  n_timesteps <- nrow(df)
  
  # Get the observation names directly from the model
  obs_names <- names(model$getObservations())
  
  # Loop over each timestep and compute observation variances and covariances
  for (i in 1:n_timesteps) {
    
    # Evaluate the Jacobian for the observations at timestep i
    C <- evaluate_obs_jacobian(model, df[i, ], pred$states[i, ], par_free)
    
    # Get the state covariance matrix for timestep i
    P <- get_state_covariance(model, pred$states[i, ])
    
    # Get the measurement variances at timestep i
    S <- get_meas_variance(model, df[i, ], pred$states[i, ], par_free)
    
    # Compute the full observation covariance matrix R
    R <- C %*% P %*% t(C) + S
    
    # Loop over the observation names and assign the variance and covariance to the right spots
    for (j in 1:length(obs_names)) {
      obs_name <- obs_names[j]
      
      # Assign variances
      pred$observations[[paste0("var.", obs_name)]][i] <- R[j, j]
      
      # Assign covariances
      if (j + 1 <= length(obs_names)) {
        for (k in (j + 1):length(obs_names)) {
          other_obs_name <- obs_names[k]
          
          # Assign covariance terms
          pred$observations[[paste0("cov.", obs_name, ".", other_obs_name)]][i] <- R[j, k]
          pred$observations[[paste0("cov.", other_obs_name, ".", obs_name)]][i] <- R[k, j]
        }
      }
    }
  }
  
  # Return updated pred with the added variances and covariances
  return(pred)
}
