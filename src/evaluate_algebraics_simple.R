evaluate_algebraics_simple <- function(model, data, par_free = NULL) {
  
  # 1. Extract model parameters 
  par <- model$getParameters()
  
  if (is.null(par_free)) {
    par_values <- if (any(is.na(par[,"estimate"]))) par[,"initial"] else par[,"estimate"]
  } else {
    # Check and update free parameters if provided
    n_free <- sum(par$type == "free")
    if (length(par_free) != n_free) stop("par_free length mismatch.")
    par[par$type == "free", "estimate"] <- par_free
    par_values <- par[,"estimate"]
  }
  names(par_values) <- rownames(par)
  
  # 2. Build the environment
  # Start with parameters, then add the specific 'data' (states/inputs)
  env <- list2env(as.list(par_values))
  
  # Overlay the user-provided data (states and inputs)
  # This keeps them logically separate in the call, even if they merge in the env
  list2env(as.list(data), envir = env)
  
  # 3. Sequential evaluation
  algebraics <- model$getAlgebraics()
  results <- list()
  
  for (alg_name in names(algebraics)) {
    rhs_expr <- rlang::f_rhs(algebraics[[alg_name]])
    
    # Eval and update env so algebraics can depend on each other
    val <- eval(rhs_expr, envir = env)
    env[[alg_name]] <- val
    results[[alg_name]] <- val
  }
  
  return(as.data.frame(results))
}