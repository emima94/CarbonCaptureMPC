# Define constraint function for steady state
g_ss_fn <- function(par) {
    alg_ss <- evaluate_algebraics_simple(model, df_ss_inputs, par_free = par)         

    # Compute sum of mean squared errors
    residuals <- c(alg_ss$Na - df_ss_outputs$Na,
        alg_ss$Nd - df_ss_outputs$Na
        )    

    g <- mean(residuals^2)
    return(g)
}