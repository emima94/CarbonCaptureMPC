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

    z = z0_normalized
    for i in 1:(length(NN_params)-2) # Exclude input and output
        W = NN_params[Symbol("layer_", i)].W
        b = NN_params[Symbol("layer_", i)].b
        z = tanh.(W * z + b)
    end

    z_output = z .* sigma_output .+ mu_output
    return z_output
end