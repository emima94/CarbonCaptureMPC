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

## Evaluate on training data

datasets = load_data_sys_id(;every_k = 6)
train_idxs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
nd_train = length(train_idxs)
datasets_train = datasets[train_idxs]

ffun = system_equations_sys_id

function hfun(x, u, dist, p)
    return [output_equations_sys_id(x, u, dist, p)]
end

plt_list = []

for d in 1:nd_train
    x0_d = x0_list[d]
    u_d = datasets_train[d].u
    dist_d = datasets_train[d].dist
    t_vec = datasets_train[d].t
    dt_sim = 180/3600

    # Find steady state for initial state
    #sol_ss = nlsolve(x -> system_equations_sys_id(x, u_d[:,1], dist_d[:,1], p_opt), x0_d)
    #x0_d = sol_ss.zero

    X_sim, Z_sim = euler_simulation(ffun, hfun, x0_d, u_d, dist_d, p_opt, t_vec, dt_sim)

    # Plot with data as dots
    plt = scatter(t_vec, datasets_train[d].Y[:,1], 
        ms = 2.5, 
        label="Data", 
        color = :red, markerstrokecolor = :red,
        title = "Series $(d)",
        ylims = (0, 4.0))

    plot!(plt, t_vec[1:end-1], Z_sim[1,:], label="Prediction", color = :black)
    
    push!(plt_list, plt)
end

for i in 6:10
    plot!(plt_list[i], xlabel = "Time [h]")
end
for i in [1,6]
    plot!(plt_list[i], ylabel = "CO2 conc. gas, " * L"y^g_A" * " [mol/m3]")
end
plt = plot(plt_list..., layout=(2,5), size=(1200,600),left_margin = 6mm, bottom_margin = 6mm)

savefig(plt, "figures/sys_id_fit.pdf")

## Simulate new trajectories
seed = 987
Random.seed!(seed)
# Get min and max values of Inputs and disturbances in training data for simulation
U_min = minimum(hcat([minimum(dataset.u, dims = 2) for dataset in datasets_train]...), dims = 2)
U_max = maximum(hcat([maximum(dataset.u, dims = 2) for dataset in datasets_train]...), dims = 2)

D_min = minimum(hcat([minimum(dataset.dist, dims = 2) for dataset in datasets_train]...), dims = 2)
D_max = maximum(hcat([maximum(dataset.dist, dims = 2) for dataset in datasets_train]...), dims = 2)

# Override with specific values for better visualization
D_min[1] = 4.0


t0 = 0
tf = 18
dt = 10/3600
Nx = Int((tf - t0) / dt) + 1
t_vec = range(t0, tf, length=Nx)

# Number of control intervals
N = 3
Nx_per_interval = Int((Nx-1) / N)

# Simulate with random inputs and disturbances within training data range
# Simulate 5 different trajectories
# Number of simulations 
N_sim = 5

for i in 1:N_sim
    U_val = rand(N,2)' .* (U_max - U_min) .+ U_min
    D_val = rand(N,2)' .* (D_max - D_min) .+ D_min
    U_sim = repeat(U_val, inner=(1,Nx_per_interval))
    D_sim = repeat(D_val, inner=(1,Nx_per_interval))

    # Find steady state for initial state
    sol_ss = nlsolve(x -> system_equations_sys_id(x, U_val[:,1], D_val[:,1], p_opt), x0_list[1])
    x0_sim = sol_ss.zero

    X_sim, Z_sim = euler_simulation(ffun, output_equations_sys_id_detailed, x0_sim, U_sim, D_sim, p_opt, t_vec, dt)
    plt = plot_simulation(X_sim, Z_sim, t_vec, U_sim,D_sim, p_opt, dt, U_min, U_max, D_min, D_max)
    savefig(plt, "figures/simulation_$(i).pdf")
end

## Get (F,Q)-grid for steady state values of x
N_grid = 50
rx = range(U_min[1], U_max[1], length = N_grid)
ry = range(U_min[2], U_max[2]+0, length = N_grid)

cap_eff = zeros(N_grid, N_grid)
SRD = zeros(N_grid, N_grid)
profit = zeros(N_grid, N_grid)

D_sim = [4.5, 180]

p_el = 90 # EUR/MWh
p_CO2 = 100  #EUR/ton

x0_guess = [2.5, 0.5, 0.5]

s2h = 1/3600
kg2ton = 1/1000
kW2MW = 1/1000

for (i, F) in enumerate(rx)
    for (j, Q) in enumerate(ry)
        U_ij = [F; Q]
        D_test = D_sim

        sol = nlsolve(x -> system_equations_sys_id(x, U_ij, D_sim, p_opt), x0_guess)
        x_ss = sol.zero

        # Compute output at steady state
        z_ss = output_equations_sys_id_detailed(x_ss, U_ij, D_sim, p_opt)
        ma = z_ss[4] # CO2 mass flow rate
        cap_eff[i,j] = z_ss[3] # Capture efficiency
        SRD[i,j] = (Q / s2h) / ma * 1e-3 # kJ/h per kg/h of CO2 captured
        profit[i,j] = -Q * kW2MW * p_el + ma * kg2ton * p_CO2

    end
end

plt_cap_eff = contour(rx, ry, cap_eff', xlabel="F [m3/h]", ylabel="Q [kW]", title="Steady state capture efficiency", fill=true, colorbar_title="Capture efficiency [%]", levels=20)
plt_SRD = contour(rx, ry, SRD', xlabel="F [m3/h]", ylabel="Q [kW]", title="Steady state specific reboiler duty", fill=true, colorbar_title="SRD [kJ/kg]", levels=20)
plt_profit = contour(rx, ry, profit', xlabel="F [m3/h]", ylabel="Q [kW]", title="Steady state profit, p_el = $(p_el) EUR/MWh, p_CO2 = $(p_CO2) EUR/ton", fill=true, colorbar_title="Profit [EUR/h]", levels=20)

# Save figures
savefig(plt_cap_eff, "figures/cap_eff_grid.pdf")
savefig(plt_SRD, "figures/SRD_grid.pdf")
savefig(plt_profit, "figures/profit_grid.pdf")
## Compute max proft as a function of CO2 price for a fixed electricity price
p_CO2_vals = range(50, 250, length=10)
profit_max = zeros(length(p_CO2_vals))

for (i, p_CO2_val) in enumerate(p_CO2_vals)
    profit_val = zeros(N_grid, N_grid)
    for (j, F) in enumerate(rx)
        for (k, Q) in enumerate(ry)
            U_ij = [F; Q]
            D_test = D_sim

            sol = nlsolve(x -> system_equations_sys_id(x, U_ij, D_sim, p_opt), x0_guess)
            x_ss = sol.zero

            # Compute output at steady state
            z_ss = output_equations_sys_id_detailed(x_ss, U_ij, D_sim, p_opt)
            ma = z_ss[4] # CO2 mass flow rate
            profit_val[j,k] = -Q * kW2MW * p_el + ma * kg2ton * p_CO2_val
        end
    end
    profit_max[i] = maximum(profit_val)
end

plt = plot(p_CO2_vals, profit_max, xlabel="CO2 price [EUR/ton]", ylabel="Maximum profit [EUR]", title="Maximum SS profit vs CO2 price, p_el = $(p_el) EUR/MWh")
savefig(plt, "figures/profit_vs_CO2_price.pdf")

## Compute the SRD vs. L/G ratio U-curve (steady state) with eta = 90 % capture eff.
function U_curve_wrap(X, F, D, p, eta)
    Q = X[1]
    x = X[2:end]
    f_eval = system_equations_sys_id(x, [F; Q], D, p)
    z_eval = output_equations_sys_id_detailed(x, [F; Q], D, p)
    cap_eff = z_eval[3]

    return [f_eval; cap_eff - eta]
end

cap_eff = [70, 80, 90, 95]

SRD = zeros(N_grid, length(cap_eff))
F_val = range(0.22, 0.50, length=N_grid)


for (i, F) in enumerate(rx)
    for (j, eta) in enumerate(cap_eff)
        sol = nlsolve(X -> U_curve_wrap(X, F, D_sim, p_opt, eta), [30.0; x0_guess])
        Q_sol = sol.zero[1]
        x_sol = sol.zero[2:end]
        z_sol = output_equations_sys_id_detailed(x_sol, [F; Q_sol], D_sim, p_opt)
        ma_sol = z_sol[4] # CO2 mass flow rate
        SRD[i,j] = (Q_sol / s2h) / ma_sol * 1e-3 # kJ/h per kg/h of CO2 captured
    end
end
plt = plot(xlabel = "L/G ratio [-]", ylabel = "SRD [MJ/kg]", title = "Steady state SRD vs L/G ratio for different capture efficiencies")
for j in eachindex(cap_eff)
    plot!(plt, F_val, SRD[:,j], lw = 2, label = L"\eta = " * "$(cap_eff[j]) %")
end
savefig(plt, "figures/SRD_vs_LG.pdf")



