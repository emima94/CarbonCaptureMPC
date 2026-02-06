evaluate_algebraics <- function(model, inputs, pred, par_free = NULL) {
  
  
  # Extract model parameters (either "estimate" or "initial" depending on NA values in estimate)
  par <- model$getParameters()
  
  if (is.null(par_free)) {
    # Check if any parameter in 'estimate' is NA
    if (any(is.na(par[,"estimate"]))) {
      # If any 'estimate' is NA, use 'initial' values for all parameters
      par_values <- as.list(par[,"initial"])
    } else {
      # Otherwise, use 'estimate' values
      par_values <- as.list(par[,"estimate"])
    } 
  }
  else {
    # Check that the length of par_free matches the number of 'free' parameters in par
    n_free <- sum(par$type == "free")
    
    if (length(par_free) != n_free) {
      stop(sprintf(
        "Length of par_free (%d) does not match the number of free parameters in model$getParameters() (%d).",
        length(par_free), n_free
      ))
    }
    par[names(par_free),"estimate"] <- par_free
    par_values <- as.list(par[,"estimate"])
    
  }
  
  # Assign parameter names
  names(par_values) <- rownames(par)

  print(names(par_values))
  print(par_values)
  
  
  algebraics <- model$getAlgebraics()
  n_steps <- nrow(inputs)
  
  if (nrow(pred$states) != n_steps) {
    stop("Inputs and pred$states must have the same number of rows")
  }
  
  alg_results <- vector("list", n_steps)
  
  for (i in seq_len(n_steps)) {
    # Combine inputs, states, and parameters for this timestep
    context <- c(
      setNames(as.list(inputs[i, , drop = TRUE]), names(inputs)),
      setNames(as.list(pred$states[i, , drop = TRUE]), names(pred$states)),
      setNames(par_values, names(par_values))  # Add parameters here
    )
    
    env <- list2env(context)
    alg_row <- list()
    
    for (alg_name in names(algebraics)) {
      rhs_expr <- rlang::f_rhs(algebraics[[alg_name]])
      
      value <- tryCatch(
        eval(rhs_expr, envir = env),
        error = function(e) {
          warning(sprintf("Could not evaluate algebraic '%s' at row %d: %s", 
                          alg_name, i, e$message))
          NA
        }
      )
      
      env[[alg_name]] <- value
      alg_row[[alg_name]] <- value
    }
    
    alg_results[[i]] <- alg_row
  }
  
  pred$algebraics <- do.call(rbind, lapply(alg_results, as.data.frame))
  return(pred)
}
