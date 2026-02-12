# Load pkgs
using Pkg

Pkg.activate(".")

using CSV, DataFrames, JSON, Revise
using DifferentialEquations, DataInterpolations
using Optimization, OptimizationOptimJL
using Plots

using ForwardDiff, ADTypes

includet("../src/util.jl")
includet("../src/model.jl")
includet("../src/simulate.jl")
includet("../src/loss_functions.jl")
includet("../src/MPC.jl")
includet("../src/MPC_scenarios.jl")

# Read parameter estimates from R
par_estimates = CSV.read("results/parameter_estimates_Tiller.csv", DataFrame)
par_fixed = CSV.read("results/parameter_fixed_Tiller.csv", DataFrame)
# Read JSON file
nn_settings = JSON.parsefile("results/nn_settings_Tiller.json")

# Convert to NamedTuple for easier access
par_estimates_nt = (; [Symbol(row.Column1) => row.x for row in eachrow(par_estimates)]...)
par_fixed_nt = (; [Symbol(row.Column1) => row.x for row in eachrow(par_fixed)]...)

# Extract NN parameters
n_neurons_A = [1, 1]
n_neurons_D = [4, 8, 1]

NA_params = extract_NN_params(par_estimates, "a", n_neurons_A, nn_settings["Na"])
ND_params = extract_NN_params(par_estimates, "d", n_neurons_D, nn_settings["Nd"])

# Extract other parameters
s = (sa = par_estimates_nt.sa, sd = par_fixed_nt.sd, st = par_fixed_nt.st, s_yga = par_fixed_nt.s_yga)
V = (Va = par_fixed_nt.Va, Vd = par_fixed_nt.Vd, Vt = par_fixed_nt.Vt)
out_allocate = (u_out = zeros(2), d_out = zeros(2)) # Placeholder for output allocation parameters, if needed

# Define parameter tuple
theta = (NA_params = NA_params, ND_params = ND_params, s = s, V = V, out_allocate = out_allocate)

# Define p tuple for ODE problem (will be updated in simulate function)
p = (u = (F = NaN, Q = NaN), d = (cgina = NaN, Fgina = NaN), theta = theta)


# Simulation Time
t0 = 0.0
tf = 6.0
dt = 15/60   # Lenght of each control interval (hr)
N = Int((tf - t0) / dt) # Number of control intervals
t_vec = range(t0, tf, length=N+1)


### Scenario 1 ###
obj_scenario01, F0, lb, ub, sd1 = define_scenario01(t0, tf, p, N);

# Test objective function
obj_val = obj_scenario01(F0, p)
println("Initial objective value: ", obj_val)

# Solve OCP
sol = solve_OCP(obj_scenario01, F0, lb, ub)

# Simulate with optimized F values

Usim = hcat(sol.u, sd1.U0[2,:])' # Combine optimized F with original Q
Dsim = sd1.D0

Xsim, Zsim, t_fine_sim = simulate_fine_grid(Usim, Dsim, sd1.prob, sd1.t_vec, 50)

# Plot results
function plot_results()
    plt1 = plot(t_vec, vcat(sd1.cap_eff_ref[1], sd1.cap_eff_ref), seriestype=:steppre, linestyle=:dash, label="Reference", title="Capture efficiency", xlabel="Time (hr)", ylabel="Capture Efficiency (%)", ylims=(50,100))
    plot!(plt1, t_fine_sim[2:end], Zsim[3,:], label="Cap eff.", title="Capture efficiency", xlabel="Time (hr)", ylabel="Capture Efficiency (%)", ylims=(50,100))

    plt2 = plot(t_fine_sim[2:end], Zsim[1:2,:]', labels=["yga" "ygina"], title="Output Trajectories with Optimized F", xlabel="Time (hr)", ylabel="Concentration (vol%)", ylims = (-2,14))
    plt3 = plot(t_fine_sim[2:end], Zsim[4,:], labels=["Na"],  title = "Capture rate", xlabel="Time (hr)", ylabel="Capture Rate (mol/hr)", ylims=(500,900))
    plt4 = plot(t_fine_sim, Xsim', labels=["ca" "cd" "ct"], title="State Trajectories with Optimized F", xlabel="Time (hr)", ylabel="Concentration", ylims = (-1,5))
    plt5 = plot(t_vec, vcat(Usim[1,1], Usim[1,:]), seriestype = :steppre, label="F", title="Optimized Control Input", xlabel="Time (hr)", ylabel="Flow Rate", ylims=(0.0,0.5))
    plot!(plt5, t_vec, vcat(lb[1], lb), seriestype=:steppre, linestyle=:dash, label="Min F", color=:red)
    plot!(plt5, t_vec, vcat(ub[1], ub), seriestype=:steppre, linestyle=:dash, label="Max F", color=:red)
    plt6 = plot(t_vec, vcat(Usim[2,1], Usim[2,:]), seriestype = :steppre, label="Q", title="Optimized Control Input", xlabel="Time (hr)", ylabel="Flow Rate", ylims = (20,35))
    plt7 = plot(t_vec, vcat(Dsim[1,1], Dsim[1,:]), seriestype = :steppre, label="cgina", title="Disturbance", xlabel="Time (hr)", ylabel="Concentration", ylims=(4,5))
    plt8 = plot(t_vec, vcat(Dsim[2,1], Dsim[2,:]), seriestype = :steppre, label="Fgina", title="Disturbance", xlabel="Time (hr)", ylabel="Flow Rate", ylims=(140,220))

    plot(plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, 
        layout = (3,4),
        size = (1600, 1200)
    )
end


plot_results()


# To do...
# Clean up
# Add Nd to plots
# Save plots for different scenarios