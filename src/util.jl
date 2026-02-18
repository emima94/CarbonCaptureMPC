function extract_NN_params(par_estimates, name, n_neurons, nn_settings)
    n_layers = length(n_neurons) - 1

    mu_input = Float64.(nn_settings["input"]["means"])
    sigma_input = Float64.(nn_settings["input"]["sds"])

    mu_output = Float64.(nn_settings["output"]["mean"])
    sigma_output = Float64.(nn_settings["output"]["sd"])

    nt = (input = (mu = mu_input, sigma = sigma_input),
            output = (mu = mu_output, sigma = sigma_output))

    for i in 1:n_layers
        W_idx = findall(row -> startswith(row.Column1, "W$(name)$(i)"), eachrow(par_estimates))
        W = reshape(par_estimates[W_idx,2], (n_neurons[i], n_neurons[i+1]))'
        
        b_idx = findall(row -> startswith(row.Column1, "b$(name)$(i)"), eachrow(par_estimates))
        b = par_estimates[b_idx,2]

        layer_i = (W = W, b = b)

        nt = (; nt..., Symbol("layer_", i) => layer_i)
        
    end         
    return nt
end



# Evaluate NN
function evaluate_NN(NN_params, z0)
    mu_input = NN_params.input.mu
    sigma_input = NN_params.input.sigma

    mu_output = NN_params.output.mu
    sigma_output = NN_params.output.sigma

    z0_normalized = (z0 .- mu_input) ./ sigma_input

    layer = 0
    z = z0_normalized
    for i in 1:(length(NN_params)-3) # Exclude input and output pars and output layer
        layer += 1
        W = NN_params[Symbol("layer_", layer)].W
        b = NN_params[Symbol("layer_", layer)].b
        z = tanh.(W * z + b)
    end

    # Output layer
    W = NN_params[Symbol("layer_", layer+1)].W
    b = NN_params[Symbol("layer_", layer+1)].b
    z = W * z + b

    # Re-scale output
    z_output = z .* sigma_output .+ mu_output

    return z_output
end


function get_inputs_and_disturbances(cgina_val, Fgina_val, F_val, Q_val, N)

    cgina = repeat(cgina_val, inner=Int(N/length(cgina_val)))
    Fgina = repeat(Fgina_val, inner=Int(N/length(Fgina_val)))
    F = repeat(F_val, inner=Int(N/length(F_val)))
    Q = repeat(Q_val, inner=Int(N/length(Q_val)))

    D = hcat(cgina, Fgina)'
    U = hcat(F, Q)'
    return U, D
end


# Callback
function callback(state, loss_val)
    if state.iter % 10 == 0  # Print every 10 iterations
        println("Iteration: $(state.iter), Loss: $(loss_val)")
    end

    if state.iter >= 500
        println("Maximum iterations reached. Stopping optimization.")
        return true  # Return true to stop optimization
    else 
        return false
    end
end

function PID_controller(setpoint, measurement, K_PID, dt, integral_prev, error_prev)

    Kp, Ki, Kd = K_PID

    error = setpoint - measurement
    integral = integral_prev + error * dt
    derivative = (error - error_prev) / dt

    output = Kp * error + Ki * integral + Kd * derivative

    return output, integral, error
end

