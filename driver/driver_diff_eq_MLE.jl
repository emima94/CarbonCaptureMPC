# System identification by minimizing NLL function
using Pkg
Pkg.activate(".")


using DifferentialEquations
import SciMLSensitivity as SMS
using Optimization, OptimizationOptimisers
using OptimizationOptimJL
using LineSearches
using ForwardDiff, Zygote
using DataInterpolations
using ComponentArrays
using Statistics
using LinearAlgebra
using Plots
using Printf
using CSV, DataFrames, Revise, JLD2
using BenchmarkTools
using LaTeXStrings
using Dates
using Random 

seed = 1234
Random.seed!(seed)

default(fontfamily="times", titlefontsize = 9, guidefontsize = 9)

includet("../src/util_diffeq_sys_id.jl")

# Load data
println("Loading data...")
datasets = load_data_sys_id(;every_k = 6)

# Training dataset idxs:
#train_idxs = [2, 3, 4, 5, 6, 7, 8, 9, 10]
train_idxs = [3]
nd_train = length(train_idxs)
datasets_train = datasets[train_idxs]

# NN settings:
NN_settings             = get_NN_settings()
NN_normalization        = get_NN_normalization(NN_settings)

# NN parameters to be estimated:
Na_params = initialize_NN_params_from_settings(NN_settings.Na)
Nd_params = initialize_NN_params_from_settings(NN_settings.Nd)

# Parameters for estimation
p0 = (Na_params = Na_params, Nd_params = Nd_params, sx = [0.1, 0.2, 0.1])

# Initial guess for initial states (to be estimated)
x0_guess = [[2.63; 0.3; 0.3] for _ in 1:length(datasets_train)]

# All parameters for estimation as ComponentArray:
θ0 = ComponentArray((p = p0, x0 = x0_guess))
    
# Fixed parameters:
P0_fixed = 0.001 * I(3)
sy = [0.25]
nx = length(x0_guess[1])
ny = size(datasets_train[1].Y, 1)
p_fixed = (NN_normalization..., P0 = P0_fixed, sy = sy, nx = nx, ny = ny)

# Make ODE problem for each dataset
prob_list = [make_ode_problem(moment_equations, ds, x0_guess[i], θ0.p, p_fixed) for (i, ds) in enumerate(datasets_train)];

### AD type ###
# Forward mode
sensealg_forward = SMS.ForwardSensitivity()
ad_type_forward = Optimization.AutoForwardDiff()
# Reverse mode
sensealg_reverse  = SMS.GaussAdjoint(autojacvec = SMS.ZygoteVJP())
#sensealg_reverse = SMS.InterpolatingAdjoint(autojacvec = SMS.ZygoteVJP())
#ad_type_reverse = Optimization.AutoReverseDiff()
ad_type_reverse = Optimization.AutoZygote() # Zygote is faster than ReverseDiff for this problem

# Choose AD type for optimization
sensealg = sensealg_forward
ad_type = ad_type_forward

# Initial evaluation of NLL
NLL_initial = NLL_aggregate(θ0, datasets_train, prob_list, sensealg, p_fixed)
println("Initial NLL: ", NLL_initial)

# Optimize parameters to minimize NLL
loss_history = Float64[]

# Optimization function
optf = Optimization.OptimizationFunction((θ, _) -> NLL_aggregate(θ, datasets_train, prob_list, sensealg, p_fixed), ad_type)

# Step 1: Adam optimizer
optimizer_adam = OptimizationOptimisers.ADAM(0.02)
prob_adam = Optimization.OptimizationProblem(optf, θ0)

println("Step 1: ADAM...")
result_adam = Optimization.solve(prob_adam, optimizer_adam;
    maxiters = 500,
    callback = callback_opt_02)

iter_adam = length(loss_history)

# Step 2: LBFGS optimizer
optimizer_lbfgs = OptimizationOptimJL.LBFGS(
    linesearch = LineSearches.BackTracking()
)
prob_lbfgs = Optimization.OptimizationProblem(optf, result_adam.u)

println("\nStep 2: LBFGS...")
result_lbfgs = Optimization.solve(prob_lbfgs, optimizer_lbfgs;
    maxiters = 100,
    callback = callback_opt_02)

sol = result_lbfgs.u
p_opt = NamedTuple(sol.p)
x0_opt = copy.(sol.x0)

# Save results
time_str = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
@save "results/sys_id_diffeq_MLE_results_$time_str.jld2" p_opt x0_opt

