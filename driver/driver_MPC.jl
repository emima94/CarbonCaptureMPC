# Load pkgs
using CSV, DataFrames, JSON


# Read parameter estimates from R
par_estimates = CSV.read("results/parameter_estimates_Tiller.csv", DataFrame)
par_fixed = CSV.read("results/parameter_fixed_Tiller.csv", DataFrame)

par_est_nt = (; [Symbol(row.Column1) => row.x for row in eachrow(par_estimates)]...)
par_fixed_nt = (; [Symbol(row.Column1) => row.x for row in eachrow(par_fixed)]...)

par_nt = (; par_fixed_nt..., par_est_nt...)

# Estimated parameters
par_est_nt = (; [Symbol(row.Column1) => row.x for row in eachrow(par_estimates)]...)
# Read JSON file
nn_settings = JSON.parsefile("results/nn_settings_Tiller.json")

# Access values
nn_settings["Na"]["input"]["names"]  # ["ca"]
nn_settings["Nd"]["input"]["means"]  # [2.6042461, 0.3707650, ...]

# Define system equations
function system_equations!(dx, x, u, p, nn_settings, include_intermediate_storage, t)
    
    # Extract states
    ca = x[1]  
    cd = x[2]
    if include_intermediate_storage
        ct = x[3]  # Intermediate storage state
    end

    # Extract inputs
    F = u[1]  # Inlet flow rate
    P_reb = u[2]  # Reboiler duty
    cgina = u[3]  # Inlet gas CO2 concentration
    Fgina = u[4]  # Inlet gas flow rate

    # Extract parameters


    
    # Define equations 
    if include_intermediate_storage
        dx[1] = (F * (ct - ca) + Na) / Va
        dx[2] = (F * (ca - cd) - Nd) / Vd
        dx[3] = F * (cd - ct) / Vt
    else
        dx[1] = (F * (cd - ca) + Na) / Va
        dx[2] = (F * (ca - cd) - Nd) / Vd
    end


end