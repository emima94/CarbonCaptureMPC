compile_nll <- function(df, train_idxs, x0, P0, 
                        nn_settings,
                        method = "laplace",
                        ode_solver = "rk4",
                        ode_stepsize = 10 / 3600,
                        incl_intermediate_storage = TRUE) {
    ### Compile NLL functions ###
    nll_fun_list <- list()
    for (i in seq_along(train_idxs)) {
        message(sprintf("Compiling NLL function for training series %d... (%d/%d)", train_idxs[i], i, length(train_idxs)))

        df_train <- df[[train_idxs[i]]]

        model_i <- create_ctsm_model(x0[[i]], P0,   
                                    incl_intermediate_storage = incl_intermediate_storage,
                                    nn_settings)

        nll_fn_i <- model_i$likelihood(df_train, 
                method = method,
                ode.solver = ode_solver, 
                ode.stepsize = ode_stepsize,
                silent = TRUE
        )
        nll_fun_list[[i]] <- nll_fn_i
    }

    return(nll_fun_list)    
}
