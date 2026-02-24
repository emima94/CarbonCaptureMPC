
# ============================================================
# Generate ODE problems for each dataset 
# ============================================================

function make_ode_problem(ffun, dataset, x0, p, p_fixed)
    function ode!(dx, x, _p, t)
        u = ConstantInterpolation(dataset.u, dataset.t; extrapolation = ExtrapolationType.Constant)
        d = ConstantInterpolation(dataset.dist, dataset.t; extrapolation = ExtrapolationType.Constant)
        dx .= ffun(x, u, d, _p, t, p_fixed)
    end

    return ODEProblem(ode!, x0, (dataset.t[1], dataset.t[end]), p)
end

# ============================================================
# Simulate one dataset, return predicted Y (n_y x N)
# ============================================================


function simulate_dataset(prob, u, d, x0, p, t_vec, sensealg, p_fixed)

    _prob = remake(prob, u0=x0, p=p)
    X  = Array(solve(_prob, Tsit5();
        tstops  = t_vec,
        saveat  = t_vec,
        abstol  = 1e-6,
        reltol  = 1e-6,
        sensealg = sensealg))

    Y_pred = hcat([hfun(X[:,j], u[:,j], d[:,j], p, p_fixed) for j in 1:length(t_vec)]...)
    return Y_pred
end

# ============================================================
# Loss function
# ============================================================
function loss(theta, datasets, prob_list, sensealg, p_fixed)
    p   = theta.p
    x0 = theta.x0

    mse = zero(eltype(p))

    for (i, ds) in enumerate(datasets)
        
        Y_pred = simulate_dataset(prob_list[i], ds.u, ds.dist, x0[i], p, ds.t, sensealg, p_fixed)

        error = Y_pred .- ds.Y
        mse += mean(error .^ 2)
    end
    
    mse_tot = mse / length(datasets)

    return sqrt(mse_tot)
end
    07
# callback
function callback_opt(state, l)
    push!(iter_log, (state.iter, l))
    mod(state.iter, 5) == 0 && println("  iter $(state.iter): loss = $(round(l, digits=8))")
    false
end

# Load data 
function load_data_sys_id(; every_k = 6)
    df = Vector{DataFrame}(undef, 10)
    t = Vector{Vector{Float64}}(undef, 10)
    Y = Vector{Matrix}(undef, 10)
    U = Vector{Matrix}(undef, 10)
    D = Vector{Matrix}(undef, 10)
    for i in 1:10
        df[i] = CSV.read("data/Tiller_edit/series_$(i).csv", DataFrame)
        t[i] = df[i].t[1:every_k:end]
        Y[i] = hcat(df[i].ygA[1:every_k:end])'
        U[i] =  hcat(df[i].F_co_A2D_vol, df[i].P_reb)[1:every_k:end,:]'#, df[i].cginA, df[i].FginA))
        D[i] = hcat(df[i].cginA, df[i].FginA)[1:every_k:end,:]'
    end

    datasets = [ (t = t[i], Y = Y[i], u = U[i], dist = D[i], N = length(t[i])) for i in 1:10 ]
    return datasets
end

# Model:

function ffun(x, u, d, p, t, p_fixed)

    # Extract states
    ca, cd, ct = x
    
    Va = 0.17
    Vd = 0.17
    Vt = 0.17

    # Conversion Constants
    M_CO2 = 44.01 # g/mol
    g2kg = 1e-3 # g to kg
    mol2kmol = 1e-3

    # Inputs #
    F, Q = u(t)

    # Disturbances #
    #cgina, Fgina = d(t)

    # Evaluate Na
    input_Na = [ca]
    Na = evaluate_NN(p.Na_params, p_fixed.Na, input_Na)[1] # Output is scalar

    # Evaluate Nd
    input_Nd = [ca, cd, F, Q]
    Nd = evaluate_NN(p.Nd_params, p_fixed.Nd, input_Nd)[1] # Output is scalar
    # ODEs:
    dx = [(F * (ct - ca) + Na * mol2kmol) / Va,
        (F * (ca - cd) - Nd * mol2kmol) / Vd,
        F * (cd - ct) / Vt
    ] 
    
    return dx

end

function ffun_Jac(x, u, d, p, t, p_fixed)

    return Zygote.jacobian(_x -> ffun(_x, u, d, p, t, p_fixed), x)[1]

end



# Create ODE problem for 1. and 2. moment equations for each dataset
moment_equations = function(x, u, dist, p, t, p_fixed)

    # Number of states
    nx = p_fixed.nx
    # Extract mean and covariance from x
    x_mean = x[1:nx]
    P = reshape(x[nx+1:end], (nx, nx))
    # Extract noise parameters from p
    sx = p.sx

    f = ffun(x_mean, u, dist, p, t, p_fixed)

    # Linearize f around current state to get A matrix (df/dx)
    #A = ForwardDiff.jacobian(_x -> ffun(_x, u, dist, p, t, p_fixed), x_mean)
    A = ffun_Jac(x_mean, u, dist, p, t, p_fixed)
    dPdt = A * P + P * A' + sx.^2 .* I(nx) # Process noise added to covariance dynamics

    return vcat(f, reshape(dPdt, :))

end

function hfun(x, u, d, p, p_fixed)
    # Extract states
    ca = x[1]

    # Conversion Constants
    R = 8.314 # J/(mol*K)
    Ta = 313.15 # K
    Pa = 101325 # Pa

    # Inputs #
    F, Q = u

    # Disturbances #
    cgina, Fgina = d

    # Evaluate Na
    input_Na = [ca]
    Na = evaluate_NN(p.Na_params, p_fixed.Na, input_Na)[1] # Output is scalar

    # Mass balance, CO2, gas phase, absorber, steady state:
    cga = cgina - Na / Fgina
    
    # Total mol concentration in gas phase (assuming ideal gas behavior and constant T and P)
    cga_tot = Pa / (R * Ta)

    # Output CO2 concentration [vol%]:
    yga = cga / cga_tot * 100

    return [yga]
end

function hfun_Jac(x, u, d, p, p_fixed)
    
    return Zygote.jacobian(_x -> hfun(_x, u, d, p, p_fixed), x)[1]

end


# Evaluate NN
function evaluate_NN(NN_params, NN_normalization, z0)
    mu_input        = NN_normalization.input.mu
    sigma_input     = NN_normalization.input.sigma
    mu_output       = NN_normalization.output.mu
    sigma_output    = NN_normalization.output.sigma

    z0_normalized   = (z0 .- mu_input) ./ sigma_input

    n_layers        = length(keys(NN_params))
    layer = 0
    z = z0_normalized

    for i in 1:(n_layers-1) # Exclude input and output pars and output layer
        layer += 1
        W = getproperty(NN_params, Symbol("layer_", layer)).W
        b = getproperty(NN_params, Symbol("layer_", layer)).b
        z = tanh.(W * z + b)
    end

    # Output layer
    W = getproperty(NN_params, Symbol("layer_", layer + 1)).W
    b = getproperty(NN_params, Symbol("layer_", layer + 1)).b
    z = W * z + b

    # Re-scale output
    z_output = z .* sigma_output .+ mu_output

    return z_output
end




function get_NN_settings()
    # NN settings:
    # Na NN:
    input_names_Na = ["ca"]
    output_names_Na = ["Na"]
    # Specify number of neurons in each layer:
    n_neurons_Na = Float64[] # (No hidden layers)

    # Normalization parameters for input and output
    mu_input_Na = [2.60]
    sigma_input_Na = [0.05]
    mu_output_Na = [648]
    sigma_output_Na = [48]

    Na = (input = (mu = mu_input_Na, sigma = sigma_input_Na),
            output = (mu = mu_output_Na, sigma = sigma_output_Na), n_neurons = n_neurons_Na)

    # Nd NN:
    input_names_Nd = ["ca", "cd", "F", "Q"]
    output_names_Nd = ["Nd"]
    # Specify number of neurons in each layer:
    n_neurons_Nd = [8] # 1 hidden layer with 8 neurons

    # Normalization parameters for input and output, Nd
    mu_input_Nd = [2.60, 0.3, 0.3, 28]
    sigma_input_Nd = [0.05, 0.22, 0.10, 3.95]
    mu_output_Nd = [648]
    sigma_output_Nd = [48]

    Nd = (input = (mu = mu_input_Nd, sigma = sigma_input_Nd),
            output = (mu = mu_output_Nd, sigma = sigma_output_Nd), n_neurons = n_neurons_Nd)

    NN_settings = (Na = Na, Nd = Nd)

    return NN_settings
end

function initialize_NN_params_from_settings(nn_setting)
    n_neurons = nn_setting.n_neurons
    n_input   = length(nn_setting.input.mu)
    n_output  = length(nn_setting.output.mu)

    layer_sizes = [n_input; Int.(n_neurons); n_output]
    n_layers    = length(layer_sizes) - 1

    layer_params = []
    for i in 1:n_layers
        W = randn(layer_sizes[i+1], layer_sizes[i]) .* 0.1
        b = zeros(layer_sizes[i+1])
        push!(layer_params, (W = W, b = b))
    end

    layer_nt = NamedTuple{Tuple(Symbol("layer_", i) for i in 1:n_layers)}(
        Tuple(layer_params)
    )

    return layer_nt  # pure NamedTuple — safe for ComponentArray
end

function get_NN_normalization(NN_settings)
    return (
        Na = (input  = NN_settings.Na.input,  output = NN_settings.Na.output),
        Nd = (input  = NN_settings.Nd.input,  output = NN_settings.Nd.output)
    )
end


## Prediction step for dataset 1 from time 1 to time 2
function ekf_predict(prob, x0, p, t_span, sensealg, p_fixed)

    nx = p_fixed.nx

    # _prob = remake(prob, u0=x0, p=p, tspan=t_span)
    # sol = Array(solve(_prob, Tsit5();
    #     tstops  = t_span,
    #     saveat  = t_span,
    #     abstol  = 1e-6,
    #     reltol  = 1e-6,
    #     sensealg = sensealg))

    # x_prior = sol[1:nx, end]
    # P_prior = reshape(sol[nx+1:end, end], (nx, nx))
    _prob = remake(prob, u0=x0, p=p, tspan=t_span)
    sol = solve(_prob, Tsit5();
        tstops  = [t_span[2]],
        saveat  = [t_span[2]],
        abstol  = 1e-6,
        reltol  = 1e-6,
        sensealg = sensealg)

    x_prior = sol[1:nx, end]
    P_prior = reshape(sol[nx+1:end, end], (nx, nx))    

    return x_prior, P_prior
end

function ekf_update(x_prior, P_prior, y_meas, u, dist, p, p_fixed)

    nx = p_fixed.nx
    ny = p_fixed.ny
    sy = p_fixed.sy

    h = hfun(x_prior, u, dist, p, p_fixed) # Predicted measurement
    e = y_meas - h # Innovation error

    # dhdx
    #H = ForwardDiff.jacobian(x -> hfun(x, u, dist, p, p_fixed), x_prior)
    H = hfun_Jac(x_prior, u, dist, p, p_fixed)

    Rv = Symmetric(sy.^2 .* I(ny)) # Measurement noise covariance

    R = Symmetric(H * P_prior * H' + Rv) # Innovation covariance
    R_chol = cholesky(R)
    K = (R_chol \ (H * P_prior))' # Kalman gain using Cholesky factorization for numerical stability

    IKH = I(nx) - K * H
    x_post = x_prior + K * e # Updated state estimate
    P_post = Symmetric(IKH * P_prior * IKH' + K * Rv * K') # Updated covariance estimate using Joseph form for numerical stability

    return x_post, P_post, R_chol, e

end

function loss_NLL(p, x0, ds, prob, sensealg, p_fixed)

    NLL = zero(eltype(p))

    x = copy(x0)
    P = copy(p_fixed.P0)

    nx = p_fixed.nx

    x0_P0 = similar(x0, nx + nx^2) # Allocate memory for combined state and covariance vector

    # Step through each time step in the dataset
    for k in 1:(length(ds.t)-1)

            # Prediction step
            x0_P0[1:nx] .= x
            x0_P0[nx+1:end] .= reshape(P, nx^2)

            x_prior, P_prior = ekf_predict(prob, x0_P0, p, (ds.t[k], ds.t[k+1]), sensealg, p_fixed)
    
            # Update step with measurement at time k+1
            y_meas = ds.Y[:,k+1]
            u = ds.u[:,k+1]
            dist = ds.dist[:,k+1]
    
            x_post, P_post, R_chol, e = ekf_update(x_prior, P_prior, y_meas, u, dist, p, p_fixed)
    
            # Compute NLL contribution from this time step
            NLL += 0.5 * (2 * sum(log.(diag(R_chol.L))) + sum(abs2, R_chol.L \ e))
    
            # Update state and covariance for next iteration
            x = x_post
            P = P_post
        end
    return NLL
end

function NLL_aggregate(θ, datasets, prob_list, sensealg, p_fixed)
    p = θ.p
    x0_vec = θ.x0

    #total_NLL = zero(eltype(p))

    NLL_parts = zeros(eltype(p), length(datasets))

    Threads.@threads for i in eachindex(datasets)
        ds = datasets[i]
        NLL_parts[i] = loss_NLL(p, x0_vec[i], ds, prob_list[i], sensealg, p_fixed)
    end

    return sum(NLL_parts) / length(datasets)
end


callback_opt_02 = (θ, loss) -> begin
    push!(loss_history, loss)
    
    # Print every 10 iterations
    if length(loss_history) % 2 == 0
        println("Iter $(length(loss_history)): loss = $loss")
    end
    
    # Early stopping: stop if loss hasn't improved in 50 iters
    if length(loss_history) > 50
        recent_improvement = loss_history[end-50] - loss
        if recent_improvement < 1e-4
            println("Early stopping at iter $(length(loss_history))")
            return true  # stop optimizer
        end
    end
    
    false
end