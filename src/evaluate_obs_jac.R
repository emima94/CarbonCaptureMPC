evaluate_obs_jacobian <- function(model, input_i, pred_state_i, par_free = NULL) {
  # Extract the Jacobian expressions for the observations
  jac_expr_list <- model$.private()$diff.terms.obs
  
  # Get observation names from the model
  obs_names <- names(jac_expr_list)
  
  # Get state names from the model's systems (not from pred_state_i directly)
  state_names <- names(model$getSystems())
  
  # Number of observations and states
  n_obs <- length(obs_names)
  n_states <- length(state_names)
  
  # Initialize the Jacobian matrix C
  C <- matrix(NA_real_, nrow = n_obs, ncol = n_states,
              dimnames = list(obs_names, state_names))
  
  # Build the evaluation context: inputs, states, and parameters
  context <- c(as.list(input_i), as.list(pred_state_i))
  par_df <- model$getParameters()
  
  # Check if we need to use initial or estimated parameters
  if (is.null(par_free)) {
    use_initial <- any(is.na(par_df$estimate))
    par_values <- if (use_initial) par_df$initial else par_df$estimate
  }
  else {
    # Check that the length of par_free matches the number of 'free' parameters in par_df
    n_free <- sum(par_df$type == "free")
    
    if (length(par_free) != n_free) {
      stop(sprintf(
        "Length of par_free (%d) does not match the number of free parameters in model$getParameters() (%d).",
        length(par_free), n_free
      ))
    }
    par_df[names(par_free),"estimate"] <- par_free
    par_values <- par_df$estimate
  }
    
  names(par_values) <- rownames(par_df)
  
  # Add parameters to the context
  context <- c(context, as.list(par_values))
  env <- list2env(context)
  
  # Loop through the observations and compute the Jacobian
  for (obs in obs_names) {
    for (state in state_names) {
      expr <- jac_expr_list[[obs]][[state]]
      C[obs, state] <- tryCatch(eval(expr, envir = env), error = function(e) NA)
    }
  }
  
  return(C)
}
