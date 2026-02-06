#' Get Measurement Variance
#'
#' This function calculates the measurement variance matrix for a given timestep
#' based on the model, inputs, and state predictions.
#'
#' @param model A ctsmTMB object representing the model.
#' @param df_i A data frame representing the inputs at timestep i.
#' @param pred_states_i A data frame representing the predicted states at timestep i.
#' 
#' @return A square matrix of measurement variances for the given timestep.
get_meas_variance <- function(model, df_i, pred_states_i, par_free = NULL) {
  # Check if the model is a ctsmTMB object
  if (!inherits(model, "ctsmTMB")) {
    stop("Error: model must be a 'ctsmTMB' object.")
  }
  
  # Check if df_i and pred_states_i are data frames
  if (!is.data.frame(df_i)) {
    stop("Error: df_i must be a data frame.")
  }
  
  if (!is.data.frame(pred_states_i)) {
    stop("Error: pred_states_i must be a data frame.")
  }
  
  # Grab those variances from the model, ready to go!
  var_expr_list <- model$getVariances()
  
  # Time to grab the parameters and make sure we're using the right ones (initial or estimated)
  par_df <- model$getParameters()
  
  # Are we using the initial ones or estimated ones? Let’s make sure!
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
  
  # Get ready to set up our square matrix where the magic happens!
  R_meas_i_matrix <- matrix(0, nrow = length(var_expr_list), ncol = length(var_expr_list))
  rownames(R_meas_i_matrix) <- names(var_expr_list)
  colnames(R_meas_i_matrix) <- names(var_expr_list)
  
  # Loop through each observation and get those variances
  for (obs_name in names(var_expr_list)) {
    # Let's get that variance expression for this observation
    var_expr <- var_expr_list[[obs_name]]
    
    # Creating the context to evaluate the expression (states, inputs, and parameters)
    context <- c(as.list(df_i), as.list(pred_states_i), as.list(par_values))
    env <- list2env(context)
    
    # Replace the variables in the expression with the actual values and evaluate it
    rhs_expr <- rlang::f_rhs(var_expr)
    
    # A little safety net for any errors during evaluation – we handle it with grace
    R_meas_i_value <- tryCatch(
      eval(rhs_expr, envir = env),
      error = function(e) {
        warning(sprintf("Uh-oh! Couldn't evaluate measurement variance for %s at timestep %d: %s", 
                        obs_name, df_i$t[1], e$message))
        NA
      }
    )
    
    # Finally, setting the variance value on the diagonal
    R_meas_i_matrix[obs_name, obs_name] <- R_meas_i_value
  }
  
  # Voila! The matrix is ready and looking sharp!
  return(R_meas_i_matrix)
}
