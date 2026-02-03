# This function adds algebraics, observation variance and state and observation CI
post_process_pred <- function(model, data, pred, par_free = NULL){

# Add algebraics to the prediction
pred <- evaluate_algebraics(model, data, pred, par_free)

# Add observation variance to the prediction
pred <- compute_obs_variance_covariance(model, data, pred, par_free)

# Add CI to states and observations
pred <- compute_CI(model,pred)

return(pred)
}