function nn_forward(input, NN, params)
    n_hidden = length(NN.n_neurons)

    # Normalize input
    x = (input .- NN.input.mu) ./ NN.input.sigma

    # Hidden layers
    for i in 1:n_hidden
        W = params[Symbol("W_$(i)")]
        b = params[Symbol("b_$(i)")]
        x = tanh.(W * x .+ b)
    end

    # Output layer + denormalize
    y_norm = params[Symbol("W_$(n_hidden + 1)")] * x .+ params[Symbol("b_$(n_hidden + 1)")]
    y = y_norm .* NN.output.sigma .+ NN.output.mu

    return y
end


function system_equations_sys_id(x, u, d, p)

    # Extract states
    ca, cd, ct = x

    # Parameters #
    NN, params, s, V = p
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Conversion Constants
    M_CO2 = 44.01 # g/mol
    g2kg = 1e-3 # g to kg
    mol2kmol = 1e-3

    # Inputs #
    F, Q = u

    # Disturbances #
    cgina, Fgina = d

    # Evaluate Na
    input_Na = [ca]
    Na = nn_forward(input_Na, NN[:Na], params[:Na])[1] # Output is scalar

    # Evaluate Nd
    input_Nd = [ca, cd, F, Q]
    Nd = nn_forward(input_Nd, NN[:Nd], params[:Nd])[1] # Output is scalar
    
    # ODEs:
    dx = [(F * (ct - ca) + Na * mol2kmol) / Va,
        (F * (ca - cd) - Nd * mol2kmol) / Vd,
        F * (cd - ct) / Vt
    ] 
    
    return dx

end

function system_equations_sys_id_w_cum_CO2(x, u, d, p)

    # Extract states
    ca, cd, ct = x

    # Parameters #
    NN, params, s, V = p
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Conversion Constants
    M_CO2 = 44.01 # g/mol
    g2kg = 1e-3 # g to kg
    mol2kmol = 1e-3

    # Inputs #
    F, Q = u

    # Disturbances #
    cgina, Fgina = d

    # Evaluate Na
    input_Na = [ca]
    Na = nn_forward(input_Na, NN[:Na], params[:Na])[1] # Output is scalar

    # Evaluate Nd
    input_Nd = [ca, cd, F, Q]
    Nd = nn_forward(input_Nd, NN[:Nd], params[:Nd])[1] # Output is scalar
    
    # ODEs:
    dx = [(F * (ct - ca) + Na * mol2kmol) / Va,
        (F * (ca - cd) - Nd * mol2kmol) / Vd,
        F * (cd - ct) / Vt,
        Na * M_CO2 * g2kg, # Cumulative CO2 captured
        Fgina * cgina * M_CO2 * g2kg # Cumulative CO2 in
    ] 
    
    return dx
end

function output_equations_sys_id(x, u, d, p)
    # Extract states
    ca = x[1]
    cd = x[2]

    # Parameters #
    NN, params, s, V = p
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Conversion Constants
    R = 8.314 # J/(mol*K)
    Ta = 313.15 # K
    Pa = 101325 # Pa
    M_CO2 = 44.01 # g/mol
    g2kg = 1e-3 # g to kg
    mol2kmol = 1e-3

    # Inputs #
    F, Q = u

    # Disturbances #
    cgina, Fgina = d

    # Evaluate Na
    input_Na = [ca]
    Na = nn_forward(input_Na, NN[:Na], params[:Na])[1] # Output is scalar

    # Evaluate Nd
    #input_Nd = [ca, cd, F, Q]
    #Nd = nn_forward(input_Nd, NN[:Nd], params[:Nd])[1] # Output is scalar

    # Mass balance, CO2, gas phase, absorber, steady state:
    cga = cgina - Na / Fgina
    
    # Total mol concentration in gas phase (assuming ideal gas behavior and constant T and P)
    cga_tot = Pa / (R * Ta)

    # Inlet CO2 concentration [vol%]:
    #ygina = cgina / cga_tot * 100

    # Output CO2 concentration [vol%]:
    yga = cga / cga_tot * 100

    # Capture efficiency:
    #cap_eff = (ygina - yga) / ygina * 100

    # Capture rate [kg/h]
    #ma = Na * M_CO2 * g2kg

    z = yga

    return z

end

function output_equations_sys_id_detailed(x, u, d, p)
    # Extract states
    ca = x[1]
    cd = x[2]

    # Parameters #
    NN, params, s, V = p
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Conversion Constants
    R = 8.314 # J/(mol*K)
    Ta = 313.15 # K
    Pa = 101325 # Pa
    M_CO2 = 44.01 # g/mol
    g2kg = 1e-3 # g to kg
    mol2kmol = 1e-3

    # Inputs #
    F, Q = u

    # Disturbances #
    cgina, Fgina = d

    # Evaluate Na
    input_Na = [ca]
    Na = nn_forward(input_Na, NN[:Na], params[:Na])[1] # Output is scalar

    # Evaluate Nd
    input_Nd = [ca, cd, F, Q]
    Nd = nn_forward(input_Nd, NN[:Nd], params[:Nd])[1] # Output is scalar

    # Mass balance, CO2, gas phase, absorber, steady state:
    cga = cgina - Na / Fgina
    
    # Total mol concentration in gas phase (assuming ideal gas behavior and constant T and P)
    cga_tot = Pa / (R * Ta)

    # Inlet CO2 concentration [vol%]:
    ygina = cgina / cga_tot * 100

    # Output CO2 concentration [vol%]:
    yga = cga / cga_tot * 100

    # Capture efficiency:
    cap_eff = (ygina - yga) / ygina * 100

    # Capture rate [kg/h]
    ma = Na * M_CO2 * g2kg

    z = [yga, ygina, cap_eff, ma, Na, Nd]

    return z

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

# Load data 
function load_data_sys_id(; every_k = 6)
    df = Vector{DataFrame}(undef, 10)
    t = Vector{Vector{Float64}}(undef, 10)
    Y = Vector{Vector{Float64}}(undef, 10)
    U = Vector{Matrix}(undef, 10)
    D = Vector{Matrix}(undef, 10)
    for i in 1:10
        df[i] = CSV.read("data/Tiller_edit/series_$(i).csv", DataFrame)
        t[i] = df[i].t[1:every_k:end]
        Y[i] = df[i].ygA[1:every_k:end]
        U[i] =  hcat(df[i].F_co_A2D_vol, df[i].P_reb)[1:every_k:end,:]'#, df[i].cginA, df[i].FginA))
        D[i] = hcat(df[i].cginA, df[i].FginA)[1:every_k:end,:]'
    end

    datasets = [ (t = t[i], Y = Y[i], u = U[i], dist = D[i], N = length(t[i])) for i in 1:10 ]
    return datasets
end
