# Functions to create a neural network in a ctsmTMB model

act_fun_string <- function(act_fun, input_expr) {
    if (act_fun == "tanh") {
        return(
            #sprintf("(exp(2*%s)-1)/(exp(2*%s)+1)", input_expr, input_expr)
            sprintf("tanh(%s)", input_expr)
            )
    } else if (act_fun == "sigmoid") {
        return(sprintf("1/(1 + exp(-%s))", input_expr))
    } else {
        stop("Unsupported activation function. Choose 'tanh' or 'sigmoid'.")
    }
}

create_ctsm_NN <- function(model, output_name, par_name, nn_settings) {
    
    # Check if nn_settings is valid
    if (is.null(nn_settings$input)) {
        stop("nn_settings must contain 'input' field (a vector of input variable names as strings).")
    }
    if (is.null(nn_settings$no_neurons)) {
        stop("nn_settings must contain 'no_neurons' field (a vector specifying number of neurons in each hidden layer).")
    }
    if (is.null(nn_settings$act_fun)) {
        stop("nn_settings must contain 'act_fun' field (activation function as string).")
    }

    # Extract neural network settings
    input_names <- nn_settings$input
    no_neurons <- nn_settings$no_neurons
    no_hidden_layers <- if (is.na(no_neurons)) 0 else length(no_neurons)
    act_fun <- nn_settings$act_fun

    message(sprintf("Creating NN for %s with %d hidden layers.", output_name, no_hidden_layers))

    
    # Inputs
    for (i in seq_along(input_names)) {
        model$setAlgebraics(
            as.formula(paste0("h", par_name, "0_", i, " ~ ", input_names[i])
                            )
        )
    }

    # Initialise number of neurons in previous layer
    no_neurons_prev <- length(input_names)

    # Hidden layers
    i <- 0
    if (no_hidden_layers > 0) {
        for (i in seq_len(no_hidden_layers)) {
            # Iterate over neurons in the layer
            for (j in 1:no_neurons[i]) {
                # Linear transformation
                model$setAlgebraics(
                    as.formula(
                        paste0(
                            "z", par_name, i, "_", j, " ~ ",
                            paste(
                                sprintf("W%s%d_%d%d * h%s%d_%d", par_name, i, j, 1:no_neurons_prev, par_name, i - 1, 1:no_neurons_prev),
                                collapse = " + "
                            ),
                            " + b", par_name, i, "_", j
                        )
                    )
                )
                # Activation function
                model$setAlgebraics(
                    as.formula(
                        paste0(
                            "h", par_name, i, "_", j, " ~ ", act_fun_string(act_fun, sprintf("z%s%d_%d", par_name, i, j))
                        )
                    )
                )
                # Add parameters for weights and biases in the hidden layer i
                for (k in 1:no_neurons_prev) {
                    pname <- paste0("W", par_name, i, "_", j, k)
                    vec <- c(initial = runif(1, -1, 1), lower = -10, upper = 10)
                    # Use do.call to build the named argument
                    args <- list()
                    args[[pname]] <- vec
                    do.call(model$setParameter, args)
                }
                pname <- paste0("b", par_name, i, "_", j)
                vec <- c(initial = 0.0, lower = -10, upper = 10)
                args <- list()
                args[[pname]] <- vec
                do.call(model$setParameter, args)
            }
            # Update number of neurons in previous layer
            no_neurons_prev <- no_neurons[i]
        }
    }

    # Output layer (we assume single output neuron)
    model$setAlgebraics(
        as.formula(
            paste0(
                "z", par_name, i + 1, "_1 ~ ",
                paste(
                    sprintf("W%s%d_1%d * h%s%d_%d",
                        par_name, i + 1, 1:no_neurons_prev, par_name, i, 1:no_neurons_prev), collapse = " + "), 
                " + b", par_name, i + 1, "_1"
            )
        )
    )
    # Add parameters for weights and biases in output layer
    for (k in 1:no_neurons_prev) {
        pname <- paste0("W", par_name, i + 1, "_1", k)
        vec <- c(initial = runif(1, -1, 1), lower = -10, upper = 10)
        # Use do.call to build the named argument
        args <- list()
        args[[pname]] <- vec
        do.call(model$setParameter, args)
    }
    pname <- paste0("b", par_name, i + 1, "_1")
    vec <- c(initial = 0.0, lower = -10, upper = 10)
    args <- list()
    args[[pname]] <- vec
    do.call(model$setParameter, args)

    ## Final output assignment
    model$setAlgebraics(
        as.formula(
            paste0(
                output_name, " ~ z", par_name, i + 1, "_1"
            )
        )
    )

    return(model)
}