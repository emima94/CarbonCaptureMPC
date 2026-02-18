# Create the Neural Network-embedded CTSM model

create_ctsm_model <- function(x0, P0,
                                incl_intermediate_storage = FALSE,
                                nn_settings) {

    model <- ctsmTMB$new()    



    if (incl_intermediate_storage) {
        model$addSystem(
            dca ~ 1/Va * (F * (ct - ca) + Na * mol2kmol) * dt + sa * dw1,
            dcd ~ 1/Vd * (F * (ca - cd) - Nd * mol2kmol) * dt + sd * dw2,
            dct ~ 1/Vt * (F * (cd - ct)) * dt + st * dw3
        )                               

        model$setParameter(
            Vt = 0.17,
            st = 0.05
        )

    } else {
        model$addSystem(
            dca ~ 1/Va * (F * (cd - ca) + Na * mol2kmol) * dt + sa * dw1,
            dcd ~ 1/Vd * (F * (ca - cd) - Nd * mol2kmol) * dt + sd * dw2
        )
    }

    model$addObs(
        yga ~ cga / cga_tot * 100
    )

    model$setVariance(
        yga ~ s_yga^2
    )

    # Required inputs
    # F, cgina, Fgina
    for (input_name in input_names) {
        do.call(model$addInput, list(as.name(input_name)))
    }

    # Constants
    model$setAlgebraics(
        mol2kmol ~ 1e-3,
        R ~ 8.314,
        Ta ~ 313.15,
        Pa ~ 101325
        #Na ~ ka * (ca_star - ca)
    )
        
    # Insert neural network for Na
    model <- create_ctsm_NN(
        model = model,
        par_name = "a",
        nn_settings = nn_settings$Na
    )

    # Insert neural network for Nd
    model <- create_ctsm_NN(
        model = model,
        par_name = "d",
        nn_settings = nn_settings$Nd
    )

    # Algebraic equations for cga and cga_tot
    model$setAlgebraics(
        cga_tot ~ Pa / (R * Ta),
        cga ~ cgina - Na / Fgina
    )

    model$setParameter(
        #ka = c(init = 2208, lb = 1e-5, ub = 1e5),
        #ca_star = c(init = 2.904, lb = 0, ub = 10),
        #ka = 2208,
        #ca_star = 2.904,
        Va = 0.17,
        Vd = 0.17,
        sa = c(init = 0.05, lower = 0.00001, upper = 10),
        sd = 0.05,
        #sd = c(init = 0.05, lower = 0.00001, upper = 1),
        s_yga = 0.2
    )

    # Initial states
    model$setInitialState(list(x0, P0))

    return(model)


}