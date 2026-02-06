compute_CI <- function(model, pred) {
  # CI multiplier for ~95% confidence
  z <- 1.96
  
  # ---- STATES ----
  state_names <- names(model$getSystems())
  for (state in state_names) {
    var_name <- paste0("var.", state)
    if (is.null(pred$states[[var_name]])) {
      stop(paste0("Variance column '", var_name, "' is missing in pred$states."))
    }
    pred$states[[paste0(state, "_CImin")]] <- pred$states[[state]] - z * sqrt(pred$states[[var_name]])
    pred$states[[paste0(state, "_CImax")]] <- pred$states[[state]] + z * sqrt(pred$states[[var_name]])
  }
  
  # ---- OBSERVATIONS ----
  obs_names <- names(model$getObservations())
  for (obs in obs_names) {
    var_name <- paste0("var.", obs)
    if (is.null(pred$observations[[var_name]])) {
      stop(paste0("Variance column '", var_name, "' is missing in pred$observations."))
    }
    pred$observations[[paste0(obs, "_CImin")]] <- pred$observations[[obs]] - z * sqrt(pred$observations[[var_name]])
    pred$observations[[paste0(obs, "_CImax")]] <- pred$observations[[obs]] + z * sqrt(pred$observations[[var_name]])
  }
  
  
  
  return(pred)
}
