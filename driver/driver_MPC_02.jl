using Pkg
Pkg.activate(".")

using JuMP, Ipopt
using CSV, DataFrames, JSON, Revise, JLD2
using NLsolve
using LinearAlgebra

using Plots
using LaTeXStrings
using Measures
using Printf
default(fontfamily="times", titlefontsize = 9, guidefontsize = 9)

# Set seed for reproducibility
using Random
seed = 1234
Random.seed!(seed)

# tickfontsize, legendfontsize

includet("../src/util.jl")
#includet("../src/system_dyn_JuMP.jl")
includet("../src/sys_id_util.jl")
includet("../src/models_JuMP.jl")
includet("../src/MPC_scenarios_JuMP.jl")
includet("../src/plots_JuMP.jl")
includet("../src/simulate_JuMP.jl")

# ----------- Load parameters --------
@load "results/sys_id_results.jld2" p_opt x0_list

p_opt
x0_list

NN, params, s, V = p_opt

## Reference control
scen_funs = [define_scenario01A, define_scenario01B, define_scenario01C]

for scen_fun in scen_funs

    p_scen = scen_fun()
    model = generate_JuMP_model_ref(p_scen, p)
    optimize!(model)
    println(solution_summary(model))
    plt = plot_ref_MPC(model, p_scen, p)
    savefig(plt, "figures/MPC_" * p_scen.name * ".pdf")
    plt    
end

## Scenario parameters
# Scenario 02A - EMPC with optimized F and Q
p02A = define_scenario02A()

#Generate JuMP model for EMPC scenario 02A
model_2A = generate_JuMP_model_EMPC(p02A, p)
optimize!(model_2A)
println(solution_summary(model_2A))
plt = plot_EMPC(model_2A, p02A, p)
savefig(plt, "figures/MPC_scenario02A.pdf")

# p02B = (; p02A..., lambda = [1e4, 1e2, 1e6]) # Keep control input and slack penalties, but remove price from objective
model_2B = generate_JuMP_model_EMPC_no_price(p02A, p)
optimize!(model_2B)
println(solution_summary(model_2B))
plt_2B = plot_EMPC(model_2B, p02A, p)
savefig(plt_2B, "figures/MPC_scenario02B.pdf")

p02C = (; p02A..., lambda = [1e5, 1e1, 1e6]) # Increase minimum capture efficiency requirement
model_2C = generate_JuMP_model_EMPC_no_price(p02C, p)
optimize!(model_2C)
println(solution_summary(model_2C))
plt_2C = plot_EMPC(model_2C, p02C, p)
savefig(plt_2C, "figures/MPC_scenario02C.pdf")  

p02D = (; p02A..., eta_min = 0.90) 
model_2D = generate_JuMP_model_EMPC(p02D, p)
optimize!(model_2D)
println(solution_summary(model_2D))
plt_2D = plot_EMPC(model_2D, p02D, p)
savefig(plt_2D, "figures/MPC_scenario02D.pdf")

p02E = (; p02A..., eta_min = 0.90, lambda = [1e5, 1e1, 1e6]) 
model_2E = generate_JuMP_model_EMPC_no_price(p02E, p)
optimize!(model_2E)
println(solution_summary(model_2E))
plt_2E = plot_EMPC(model_2E, p02E, p)
savefig(plt_2E, "figures/MPC_scenario02E.pdf")

## Plot profit as function of min capture efficiency requirement for a price sensitive and non-pricesensitive model
eta_min_vals = range(0.70, 0.98, length = 10)
profit_price_sensitive = zeros(length(eta_min_vals))
profit_non_price_sensitive = zeros(length(eta_min_vals))

for (i, eta_min) in enumerate(eta_min_vals)
    println("Simulating for eta_min = $(round(eta_min, digits=2))")
    p_temp = (; p02A..., eta_min = eta_min)
    model_temp_price = generate_JuMP_model_EMPC(p_temp, p)
    optimize!(model_temp_price)
    _, profit_temp = plot_EMPC(model_temp_price, p_temp, p)
    profit_price_sensitive[i] = profit_temp

    p_temp_no_price = (; p02A..., eta_min = eta_min, lambda = [1e5, 1e1, 1e6])
    model_temp_no_price = generate_JuMP_model_EMPC_no_price(p_temp_no_price, p)
    optimize!(model_temp_no_price)
    _, profit_temp = plot_EMPC(model_temp_no_price, p_temp_no_price, p)
    profit_non_price_sensitive[i] = profit_temp
end

# Plot profit vs eta_min for price sensitive and non-price sensitive models
plt_profit_eta = plot(eta_min_vals*100, profit_price_sensitive, label="Price sensitive", xlabel="Minimum average capture efficiency requirement [%]", ylabel="Total profit [EUR]")
plot!(plt_profit_eta, eta_min_vals*100, profit_non_price_sensitive, label="Non-price sensitive")
savefig(plt_profit_eta, "figures/profit_vs_eta_min.pdf")


## Simulate pre-specified control trajectories

# p_sim = define_scenario_sim()
# model_sim = generate_JuMP_model_EMPC(p_sim, p)
# plt_sim = plot_EMPC(model_sim, p_sim, p, p_sim.U_sim)


## Plot Nd vs. Q and F

F_val = range(0.22, 0.32, length=250)
Q_val = range(22,35, length=250)

# Simulate with constant control inputs and plot Ma[tf] and total profit as function of F and Q
Ma_tf = zeros(length(F_val), length(Q_val))
profit = zeros(length(F_val), length(Q_val))
cap_eff = zeros(length(F_val), length(Q_val))

x0 = [2.63; 0.3; 0.3; 0.0; 0.0]

ffun = system_equations_w_cum_CO2
hfun = output_equations
t0 = 0
tf = 12
dt = 0.25
N = Int((tf - t0) / dt)  # number of intervals
t = range(t0, tf, length=N+1)

Nx_per_interval = 10
Nx = N * Nx_per_interval + 1   # total number of state nodes
dt_sim = dt / Nx_per_interval
t_sim = range(t0, tf, length=Nx)

cgina = 4.4
Fgina = 180
Dsim = hcat(repeat([cgina], inner=Nx-1), repeat([Fgina], inner=Nx-1))'
p_el_val = [0.2, 0.2, 1.0, 1.0, 0.2, 0.2] # EUR/kWh
p_el = repeat(p_el_val, inner=Int((Nx-1)/length(p_el_val)))
p_CO2 = 0.07 * 10  #EUR/kg

for (i,F) in enumerate(F_val)
    for (j,Q) in enumerate(Q_val)
        println("Simulating for F = $F, Q = $Q")

        F_sim = repeat([F], inner = Nx-1)
        Q_sim = repeat([Q], inner = Nx-1)

        U_sim = hcat(F_sim, Q_sim)'

        X, Z = euler_simulation(ffun, hfun, x0, U_sim, Dsim, p, t_sim, dt_sim)

        Ma_tf[i,j] = X[4,end] # Cumulative CO2 captured at final time step
        profit[i,j] = -sum(p_el .* Q_sim) * dt_sim + p_CO2 * Ma_tf[i,j]
        cap_eff[i,j] = X[4,end] / X[5,end]
    end
end

# Plot Ma[tf] and profit as function of F and Q
plt_Ma = contour(F_val, Q_val, Ma_tf', xlabel="F [m3/h]", ylabel="Q [kW]", title="Cumulative CO2 captured at tf", fill=true, colorbar_title="Ma[tf]", levels=20)
plt_profit = contour(F_val, Q_val, profit', xlabel="F [m3/h]", ylabel="Q [kW]", title="Total profit [EUR]", fill=true, colorbar_title="Profit [EUR]", levels=20)
plt_cap_eff = contour(F_val, Q_val, cap_eff'*100, xlabel="F [m3/h]", ylabel="Q [kW]", title="Average capture efficiency", fill=true, colorbar_title="Capture efficiency [%]", levels=20)

savefig(plt_Ma, "figures/Ma_tf_vs_F_Q.pdf")
savefig(plt_profit, "figures/profit_vs_F_Q.pdf")
savefig(plt_cap_eff, "figures/cap_eff_vs_F_Q.pdf")

## Closed-loop MPC simulations

seed = 1234
Random.seed!(seed)

# Simulation horizon for closed-loop simulation
t0 = 0
tf = 12
dt = 0.25
N_sim = Int((tf - t0) / dt)  # number of intervals

# Define scenario:
# Initial state (prior):
x0 = [2.63; 0.3; 0.3; 0.0; 0.0]
P0 = Diagonal([0.04^2, 0.02^2, 0.02^2, 0.001^2, 0.001^2]) # Initial state covariance matrix (uncertainty in initial state)

# Initial controls:
F_prev = 0.3
Q_prev = 28

# Electricity and CO2 prices
p_el_val = [0.1, 0.12, 0.12, 0.24, 0.26, 0.31, 0.25, 0.15, 0.02, 0.02, 0.05, 0.03, 0.16, 0.44, 0.65, 0.56, 0.65, 0.7, 0.10, 0.04, 0.02, 0.10, 0.09, 0.08] # EUR/kWh
p_el = repeat(p_el_val, inner=Int(N_sim/length(p_el_val)))
p_CO2 = 0.07 * 10  #EUR/kg

# Disturbances
cgina_val = [4.4, 4.5, 4.3, 4.4, 4.6, 4.5]
Fgina_val = [180]
cgina = repeat(cgina_val, inner=Int(N_sim/length(cgina_val))) # Expexted disturbance (known to controller)
Fgina = repeat(Fgina_val, inner=Int(N_sim/length(Fgina_val))) # Expexted disturbance (known to controller)

D = hcat(cgina, Fgina)'

# Process and measurement noise  
sigma_x = [0.04, 0.02, 0.02, 0.0, 0.0] # Standard deviation of process noise for each state variable
sigma_y = 0.3 # Standard deviation of measurement noise for output variable yga

# MPC tuning
# MPC tuning Parameters
lambda_F = 1e3
lambda_Q = 1e-1
lambda_s = 1e6
lambda = [lambda_F, lambda_Q, lambda_s]

# Required average capture efficiency at final time step
eta_min = 0.90

# ODE rhs:
ffun = system_equations_w_cum_CO2
# Output equations
hfun = output_equations

# Control horizon for MPC
t_horizon = 4.0
N = Int(t_horizon / dt)  # number of intervals in control horizon

Nx_per_interval = 10
Nx = N * Nx_per_interval + 1   # total number of state nodes in horizon
dt_sim = dt / Nx_per_interval

interval_of = repeat(1:N, inner=Nx_per_interval)
push!(interval_of, N)  # last state belongs to last interval

# Bounds 
F_min, F_max = 0.22, 0.32
Q_min, Q_max = 22.0, 35.0
s_min = 0.0

lb = (F_min = F_min, Q_min = Q_min, s_min = s_min)
ub = (F_max = F_max, Q_max = Q_max, s_max = Inf)

# Create parameter tuple for scenario and update as horizon moves forward in closed-loop simulation
p_CL = (ffun = ffun, hfun = hfun, x0 = x0, D = D, t_sim = t_sim, t = t, dt = dt, dt_sim = dt_sim, Q_prev = Q_prev, F_prev = F_prev, p_el = p_el, p_CO2 = p_CO2, lambda = lambda, eta_min = eta_min, N = N, Nx = Nx, interval_of = interval_of, lb = lb, ub = ub)

# Allocate arrays to store closed-loop optimal trajectories and solver information
F_opt_CL = Vector{Vector{Float64}}(undef, N_sim)
Q_opt_CL = Vector{Vector{Float64}}(undef, N_sim)


# step 1:
# Solve OCP
x0_horizon = x0 # In first iteration, use true initial state as initial condition for MPC. In subsequent iterations, this will be updated with state estimates from EKF.
i = 1 # index for current time step in closed-loop simulation
idx_horizon = i:(i+N-1)
D_horizon = D[:, idx_horizon]
p_el_horizon = p_el[idx_horizon]

# update parameter tuple for MPC horizon
p_horizon = (; p_CL..., x0 = x0_horizon, D = D_horizon, p_el = p_el_horizon)

# Generate JuMP model
model_horizon = generate_JuMP_model_EMPC(p_horizon, p)
# Solve OCP
optimize!(model_horizon)

# Extract entire optimal control trajectory
F_opt = value.(model_horizon[:F])
Q_opt = value.(model_horizon[:Q])

# Store optimal control trajectory in closed-loop arrays
F_opt_CL[i] = F_opt
Q_opt_CL[i] = Q_opt



# Simulate with constant control inputs and plot Ma[tf] and total profit as function of F and Q



