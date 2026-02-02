# Create the Neural Network-embedded CTSM model

create_ctsm_model <- function(x0, P0,
                                input_names,
                                incl_intermediate_storage = FALSE,
                                nn_settings) {

    model <- ctsmTMB$new()    

    if (incl_intermediate_storage) {
        model$addSystem(
            dca ~ 1/Va * (F * (ct - ca) + Na) * dt + sa * dw1,
            dcd ~ 1/Vd * (F * (ca - cd) - Nd) * dt + sd * dw2,
            dct ~ 1/Vt * (F * (cd - ct)) * dt + st * dw3
        )                               

        model$setParameter(
            Vt = 0.17,
            st = 0.01
        )

    } else {
        model$addSystem(
            dca ~ 1/Va * (F * (cd - ca) + Na) * dt + sa * dw1,
            dcd ~ 1/Vd * (F * (ca - cd) - Nd) * dt + sd * dw2
        )
    }

    model$addObs(
        yga ~ cga / cga_tot * 100
    )

    model$setVariance(
        yga ~ s_yga^2
    )

    # Required inputs
    # F, cgin_a, Fgin_a
    for (input_name in input_names) {
        do.call(model$addInput, list(as.name(input_name)))
    }

    model$setAlgebraics(
        cga_tot ~ Pa / (R * Ta),
        cga ~ c_gina - Na / F_gina
    )
        
    # Insert neural network for Na
    model <- create_ctsm_NN(
        model = model,
        output_name = "Na",
        par_name = "a",
        nn_settings = nn_settings$Na
    )
    # Insert neural network for Nd
    model <- create_ctsm_NN(
        model = model,
        output_name = "Nd",
        par_name = "d",
        nn_settings = nn_settings$Nd
    )

    model$setParameter(
        Va = 0.16,
        Vd = 0.12,
        R = 8.314,
        Ta = 313.15,
        Pa = 101325,
        sa = 0.01,
        sd = 0.01,
        s_yga = 0.05
    )

    # Initial states
    model$setInitialState(list(x0, P0))

    return(model)


}